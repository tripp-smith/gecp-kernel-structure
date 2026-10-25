"""A posteriori pivot quality for the true selected cross on a finite grid.

The exact physical kernel/core are interval enclosed independently of the
facade's rounded normalization, LU and exported cardinal coefficients. This
does not certify a continuous domain or an exact complete-pivot trajectory.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Any

import numpy as np
import numpy.typing as npt

from .basis import FermionicBasis, _tau_array
from .fitting import _numeric_array
from .kernels import FermionicKernel
from .precision import _IntervalBasis, _upper_abs
from .types import NumericalBreakdownError


@dataclass(frozen=True, slots=True)
class GridPrecisionCertificate:
    """Content-bound finite-grid discrepancy and multiplicative pivot bound.

    ``None`` bounds mean no enclosure was completed. A valid enclosure can
    still fail to establish positive binary64 pivot quality. ``eta`` bounds
    the candidate's exact residual magnitude divided by the exact grid maximum
    when that maximum is positive; an exactly zero grid has quality one.
    It says nothing about a maximum between these particular grid points.
    """

    certified: bool
    termination_reason: str
    content_sha256: str
    grid_shape: tuple[int, int]
    rank: int
    precision_bits: int
    estimated_work: int
    candidate_index: tuple[int, int] | None = None
    candidate_tau: float | None = None
    candidate_omega: float | None = None
    rounded_candidate_magnitude: float | None = None
    rounded_global_upper: float | None = None
    discrepancy_upper: float | None = None
    true_candidate_lower: float | None = None
    true_global_upper: float | None = None
    eta: float | None = None
    exact_tie_count: int = 0
    near_tie_count: int = 0
    core_pivot_lower_bounds: tuple[float, ...] = ()
    conditioning_scope: str = "ideal_core_pivot_regularity_not_a_condition_number"
    tie_rtol: float = 0.0
    tie_scope: str = "exact_binary64_grid_residual_magnitudes"
    objective_scope: str = "exact_physical_selected_cross_on_supplied_finite_grid"
    arithmetic: str = "owned_mpmath_interval_core_lu_and_outward_exports"
    coordinate_precision_bits: int = 53
    bound_precision_bits: int = 53
    continuous_certificate: bool = False


def _bounded_integer(value: int, name: str, lower: int, upper: int) -> int:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or not lower <= value <= upper
    ):
        raise ValueError(f"{name} must be an integer in [{lower}, {upper}]")
    return int(value)


def _directed(value: Fraction, upward: bool) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ArithmeticError("finite bound exceeds binary64 export range")
    needs_adjustment = Fraction(result) < value if upward else Fraction(result) > value
    if needs_adjustment:
        result = math.nextafter(result, math.inf if upward else -math.inf)
    if not math.isfinite(result):
        raise ArithmeticError("finite bound exceeds binary64 export range")
    return result


def certify_grid_pivot_precision(
    basis: FermionicBasis,
    tau: npt.ArrayLike,
    omega: npt.ArrayLike,
    *,
    precision_bits: int = 256,
    max_rank: int = 64,
    max_work: int = 2**22,
) -> GridPrecisionCertificate:
    """Enclose true/rounded selected-cross discrepancy, then transfer quality.

    Times are physical; omega values are dimensionless. The candidate is the
    first row-major maximum of the actual binary64 residual. Exact dyadic input
    values define the ideal target; every division by beta is enclosed.
    Budgets and zero-containing interval pivots yield an uncertified report.
    No explicit inverse or a priori BLAS rounding model is used.
    """
    if not isinstance(basis, FermionicBasis):
        raise TypeError("basis must be a FermionicBasis")
    bits = _bounded_integer(precision_bits, "precision_bits", 53, 4096)
    rank_budget = _bounded_integer(max_rank, "max_rank", 1, 256)
    work_budget = _bounded_integer(max_work, "max_work", 1, 2**30)
    times = np.array(_tau_array(tau, basis.beta), copy=True)
    frequencies = np.array(_numeric_array(omega, "omega"), copy=True)
    if frequencies.ndim != 1 or np.iscomplexobj(frequencies):
        raise ValueError("omega must be a one-dimensional real array")
    if not times.size or not frequencies.size:
        raise ValueError("both finite query grids must be nonempty")
    if np.any(np.abs(frequencies) > basis.cutoff):
        raise ValueError("omega must lie within the basis dimensionless cutoff")
    # The estimate counts precision-weighted dense scalar work, including the
    # Cartesian residual dots; rejecting precedes interval allocation.
    rank = basis.rank
    work = (
        rank**3
        + times.size * rank**2
        + times.size * frequencies.size * (rank + 1)
        + rank * frequencies.size
    ) * ((bits + 63) // 64)
    digest = hashlib.sha256(b"kernelgecp-grid-precision-v1")
    for array in (
        np.array([basis.beta, basis.wmax], dtype=np.float64),
        basis.tau_nodes,
        basis.omega_nodes,
        times,
        frequencies,
    ):
        digest.update(str(array.shape).encode())
        digest.update(array.dtype.str.encode())
        digest.update(array.tobytes())
    metadata: dict[str, Any] = {
        "content_sha256": digest.hexdigest(),
        "grid_shape": (int(times.size), int(frequencies.size)),
        "rank": rank,
        "precision_bits": bits,
        "estimated_work": work,
    }

    def incomplete(reason: str) -> GridPrecisionCertificate:
        return GridPrecisionCertificate(False, reason, **metadata)

    if rank > rank_budget:
        return incomplete("rank_budget")
    if work > work_budget:
        return incomplete("work_budget")
    try:
        exact = _IntervalBasis(basis, bits)
        metadata["core_pivot_lower_bounds"] = tuple(
            exact.transpose_lu.pivot_lower_bounds
        )
        rows = exact.tau_rows(times)
        columns = [
            [exact.kernel(float(node), float(w)) for w in frequencies]
            for node in basis.tau_nodes
        ]
        with np.errstate(all="ignore"):
            rounded_design = basis.design_tau(times)
            kernel = FermionicKernel(basis.cutoff)
            rounded = kernel(
                times[:, None] / basis.beta, frequencies[None, :]
            ) - rounded_design @ kernel(
                basis.tau_nodes[:, None] / basis.beta, frequencies[None, :]
            )
        if not np.all(np.isfinite(rounded)):
            return incomplete("nonfinite_rounded_residual")
        # Include actual rounded values in the report identity. Query/basis
        # inputs alone cannot identify an implementation-dependent BLAS result.
        digest.update(rounded.dtype.str.encode())
        digest.update(rounded.tobytes())
        metadata["content_sha256"] = digest.hexdigest()
        discrepancy = 0.0
        for i, time in enumerate(times):
            for j, frequency in enumerate(frequencies):
                ideal = exact.kernel(float(time), float(frequency)) - sum(
                    (rows[i][k] * columns[k][j] for k in range(rank)),
                    exact.ctx.mpf(0),
                )
                discrepancy = max(discrepancy, _upper_abs(ideal - float(rounded[i, j])))
        magnitudes = np.abs(rounded)
        location = np.unravel_index(int(np.argmax(magnitudes)), magnitudes.shape)
        i, j = int(location[0]), int(location[1])
        magnitude = float(magnitudes[i, j])
        ties = int(np.count_nonzero(magnitudes == magnitude))
        lower = max(Fraction(0), Fraction(magnitude) - Fraction(discrepancy))
        upper = Fraction(magnitude) + Fraction(discrepancy)
        quality = _directed(lower / upper, False) if upper else 1.0
        positive = quality > 0
        return GridPrecisionCertificate(
            certified=positive,
            termination_reason=(
                "exact_zero_grid_residual"
                if not upper
                else (
                    "positive_finite_grid_quality"
                    if positive
                    else "discrepancy_dominates_candidate"
                )
            ),
            candidate_index=(i, j),
            candidate_tau=float(times[i]),
            candidate_omega=float(frequencies[j]),
            rounded_candidate_magnitude=magnitude,
            rounded_global_upper=magnitude,
            discrepancy_upper=discrepancy,
            true_candidate_lower=_directed(lower, False),
            true_global_upper=_directed(upper, True),
            eta=quality,
            exact_tie_count=ties,
            near_tie_count=ties,
            **metadata,
        )
    except (ArithmeticError, ValueError, NumericalBreakdownError) as error:
        return incomplete(f"interval_inconclusive:{type(error).__name__}")
