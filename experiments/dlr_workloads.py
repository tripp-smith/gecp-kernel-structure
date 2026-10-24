"""Independent finite-pole Green functions and common physical time meshes.

These workloads are synthetic physical forward evaluations, not measurements or
continuous-domain certificates. Their truth values do not use a fitted DLR or
GECP representation. All methods receive the same training mesh; validation
midpoints are held out from fitting, rank selection and stopping.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real
from typing import Any

import numpy as np
import numpy.typing as npt

from kernelgecp.grids import chebyshev_lobatto

FloatArray = npt.NDArray[np.float64]


def _positive_scalar(name: str, value: float) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a real scalar")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _physical_parameters(beta: float, wmax: float) -> tuple[float, float, float]:
    beta = _positive_scalar("beta", beta)
    wmax = _positive_scalar("wmax", wmax)
    cutoff = beta * wmax
    if not math.isfinite(cutoff) or cutoff == 0:
        raise ValueError("beta * wmax must be finite and positive in float64")
    return beta, wmax, cutoff


def _real_vector(name: str, values: Any) -> FloatArray:
    raw = np.asarray(values)
    if raw.ndim != 1 or not raw.size:
        raise ValueError(f"{name} must be a nonempty one-dimensional real array")
    if raw.dtype.kind not in "fiu":
        raise ValueError(f"{name} must contain real numeric values")
    array = np.asarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class PoleWorkload:
    """A positive normalized real spectral measure with distinct finite poles.

    Inputs are copied into immutable tuples. Normalization is checked to an
    absolute tolerance of ``1e-14`` to accommodate binary64 addition; weights
    are never silently renormalized.
    """

    name: str
    poles: tuple[float, ...]
    weights: tuple[float, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a nonempty string")
        poles = _real_vector("poles", self.poles)
        weights = _real_vector("weights", self.weights)
        if poles.shape != weights.shape:
            raise ValueError("poles and weights must have the same shape")
        if np.unique(poles).size != poles.size:
            raise ValueError("poles must be distinct")
        if np.any(weights <= 0):
            raise ValueError("weights must be positive")
        if np.any(weights > 1.0 + 1e-14):
            raise ValueError("weights must be normalized to one")
        if not math.isclose(math.fsum(weights), 1.0, rel_tol=0, abs_tol=1e-14):
            raise ValueError("weights must be normalized to one")
        object.__setattr__(self, "poles", tuple(float(pole) for pole in poles))
        object.__setattr__(self, "weights", tuple(float(weight) for weight in weights))


def make_pole_workloads(*, beta: float, wmax: float) -> tuple[PoleWorkload, ...]:
    """Return three fixed physical fixtures, scaled only by ``wmax``.

    The pole locations and weights are truth data used to form forward values.
    They are not candidate atoms, sampling nodes or fitting inputs.
    """

    _, wmax, _ = _physical_parameters(beta, wmax)
    return (
        PoleWorkload("central_pole", (0.0,), (1.0,)),
        PoleWorkload(
            "asymmetric_satellites",
            tuple(wmax * fraction for fraction in (-0.75, 0.08, 0.55)),
            (0.2, 0.5, 0.3),
        ),
        PoleWorkload(
            "gapped",
            tuple(wmax * fraction for fraction in (-0.8, -0.3, 0.35, 0.9)),
            (0.15, 0.35, 0.3, 0.2),
        ),
    )


def green_function(
    workload: PoleWorkload,
    tau: Any,
    *,
    beta: float,
    wmax: float,
) -> FloatArray:
    r"""Evaluate ``G(tau) = -sum w exp(-tau*p)/(1+exp(-beta*p))``.

    For negative poles the algebraically equivalent exponent
    ``exp((beta-tau)*p)/(1+exp(beta*p))`` avoids overflow. Every exponential
    argument is nonpositive. Inputs and output are binary64; underflow can
    occur and is not a claim that the exact forward value is zero.
    """

    beta, wmax, _ = _physical_parameters(beta, wmax)
    times = _real_vector("tau", tau)
    if np.any(times < 0) or np.any(times > beta):
        raise ValueError("tau must lie in the physical interval [0, beta]")
    poles = np.asarray(workload.poles, dtype=np.float64)
    if np.any(np.abs(poles) > wmax):
        raise ValueError("poles must lie in the physical interval [-wmax, wmax]")
    weights = np.asarray(workload.weights, dtype=np.float64)
    terms = np.empty((times.size, poles.size), dtype=np.float64)
    positive = poles >= 0
    # Local underflow is permitted even when callers enable NumPy traps.
    with np.errstate(under="ignore"):
        terms[:, positive] = np.exp(-times[:, None] * poles[positive]) / (
            1.0 + np.exp(-beta * poles[positive])
        )
        terms[:, ~positive] = np.exp((beta - times[:, None]) * poles[~positive]) / (
            1.0 + np.exp(beta * poles[~positive])
        )
    return np.asarray(-(terms @ weights), dtype=np.float64)


def validate_tau_meshes(training_tau: Any, validation_tau: Any, *, beta: float) -> None:
    """Reject nonphysical, unordered, duplicate or overlapping sample meshes."""

    beta = _positive_scalar("beta", beta)
    training = _real_vector("training_tau", training_tau)
    validation = _real_vector("validation_tau", validation_tau)
    for name, nodes in (("training_tau", training), ("validation_tau", validation)):
        if np.any(nodes < 0) or np.any(nodes > beta):
            raise ValueError(f"{name} must lie in [0, beta]")
        if np.any(np.diff(nodes) <= 0):
            raise ValueError(f"{name} must be strictly increasing without duplicates")
    if np.intersect1d(training, validation).size:
        raise ValueError("training and validation meshes must be disjoint")


def common_tau_meshes(
    *, beta: float, wmax: float, order: int = 16
) -> tuple[FloatArray, FloatArray]:
    """Return common composite-Chebyshev training and midpoint validation times.

    The canonical dimensionless grid resolves endpoint layers using only
    ``max(1, beta*wmax)``. Scaling maps it to physical time ``[0, beta]``.
    Validation has one interior midpoint per training interval, with no forward
    values consulted. Coordinate-resolution failures are rejected rather than
    silently discarding repeated or overlapping samples.
    """

    beta, _, cutoff = _physical_parameters(beta, wmax)
    if isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer)):
        raise ValueError("order must be an integer of at least two")
    if order < 2:
        raise ValueError("order must be an integer of at least two")
    first_right = min(0.5, 1.0 / max(1.0, cutoff))
    panels: list[FloatArray] = []
    left, right = 0.0, first_right
    while left < 0.5:
        nodes = chebyshev_lobatto(left, right, int(order))
        # Pin common panel endpoints before joining. Ordinary midpoint
        # arithmetic can round nominally equal Chebyshev endpoints to adjacent
        # floats; their reflection would introduce duplicate training times.
        nodes[0], nodes[-1] = left, right
        panels.append(nodes if not panels else nodes[1:])
        left, right = right, min(0.5, 2.0 * right)
    half = np.concatenate(panels)
    dimensionless = np.concatenate((half, 1.0 - half[-2::-1]))
    training = np.asarray(dimensionless * beta)
    training[0], training[-1] = 0.0, beta
    validation = training[:-1] + (training[1:] - training[:-1]) / 2.0
    validate_tau_meshes(training, validation, beta=beta)
    return training, validation
