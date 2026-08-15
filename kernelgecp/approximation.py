"""Explicit-rank bounds and numerical separated approximations."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.interpolate import BarycentricInterpolator

from .grids import chebyshev_lobatto
from .kernels import FermionicKernel
from .types import FloatArray


def exp_family_rank_bound(max_rate: float, tolerance: float) -> int:
    r"""Return the explicit Lemma 4.4 exponential-family node count.

    The count is
    ``(floor(log2(M)) + 1) * (floor(log4(1/tolerance)) + 1)``.
    This function reports the published arithmetic bound; it does not claim
    that the repository has yet formalized the paper's construction in Lean.
    """

    if not np.isfinite(max_rate) or max_rate < 1:
        raise ValueError("max_rate must be finite and at least one")
    if not np.isfinite(tolerance) or not (0 < tolerance < 1):
        raise ValueError("tolerance must lie strictly between zero and one")
    scales = math.floor(math.log2(max_rate)) + 1
    orders = math.floor(math.log(1.0 / tolerance, 4.0)) + 1
    return scales * orders


def dlr_rank_bound(cutoff: float, tolerance: float) -> int:
    """A conservative explicit count for positive, negative, and central bands."""

    band = exp_family_rank_bound(cutoff, tolerance / 3.0)
    central = math.ceil(math.log(3.0 / tolerance)) + 1
    return 2 * band + central


def _lagrange_values(
    nodes: FloatArray, weights: FloatArray, x: FloatArray
) -> FloatArray:
    differences = x[:, None] - nodes[None, :]
    exact = differences == 0
    safe = np.where(exact, 1.0, differences)
    numerators = weights[None, :] / safe
    values = numerators / np.sum(numerators, axis=1, keepdims=True)
    for row, matches in enumerate(exact):
        if np.any(matches):
            values[row, :] = 0.0
            values[row, int(np.flatnonzero(matches)[0])] = 1.0
    return np.asarray(values, dtype=np.float64)


@dataclass(slots=True)
class SeparatedApproximation:
    """Frequency-interpolation representation ``Σ K(t,ω_j) L_j(ω)``."""

    kernel: FermionicKernel
    omega_nodes: FloatArray
    barycentric_weights: FloatArray
    band_slices: list[tuple[int, int]]
    band_bounds: list[tuple[float, float]]

    @property
    def rank(self) -> int:
        return int(self.omega_nodes.size)

    def evaluate(self, t: npt.ArrayLike, omega: npt.ArrayLike) -> FloatArray:
        t_values = np.atleast_1d(np.asarray(t, dtype=np.float64))
        omega_values = np.atleast_1d(np.asarray(omega, dtype=np.float64))
        coefficients = np.asarray(
            self.kernel(t_values[:, None], self.omega_nodes[None, :]),
            dtype=np.float64,
        )
        basis = np.zeros((omega_values.size, self.rank), dtype=np.float64)
        assigned = np.zeros(omega_values.size, dtype=bool)
        for (start, stop), (left, right) in zip(
            self.band_slices, self.band_bounds, strict=True
        ):
            mask = (~assigned) & (omega_values >= left) & (omega_values <= right)
            if np.any(mask):
                basis[mask, start:stop] = _lagrange_values(
                    self.omega_nodes[start:stop],
                    self.barycentric_weights[start:stop],
                    omega_values[mask],
                )
                assigned[mask] = True
        if not np.all(assigned):
            raise ValueError(
                "evaluation frequency lies outside the approximation domain"
            )
        return np.asarray(coefficients @ basis.T, dtype=np.float64)


def fermionic_separated_approximation(
    cutoff: float, order: int = 12
) -> SeparatedApproximation:
    """Build a deterministic composite-Chebyshev frequency interpolant."""

    kernel = FermionicKernel(cutoff)
    positive_bounds: list[tuple[float, float]] = [(0.0, 1.0)]
    left = 1.0
    while left < cutoff:
        right = min(cutoff, 2.0 * left)
        positive_bounds.append((left, right))
        left = right
    negative_bounds = [(-right, -left) for left, right in reversed(positive_bounds)]
    bounds = negative_bounds + positive_bounds
    node_parts = [chebyshev_lobatto(left, right, order) for left, right in bounds]
    weight_parts = [
        np.asarray(
            BarycentricInterpolator(nodes, np.zeros_like(nodes)).wi,
            dtype=np.float64,
        )
        for nodes in node_parts
    ]
    slices: list[tuple[int, int]] = []
    offset = 0
    for nodes in node_parts:
        slices.append((offset, offset + nodes.size))
        offset += nodes.size
    return SeparatedApproximation(
        kernel,
        np.concatenate(node_parts),
        np.concatenate(weight_parts),
        slices,
        bounds,
    )
