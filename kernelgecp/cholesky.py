"""Canonical diagonal-pivoted Cholesky for positive-semidefinite matrices."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .diagnostics import conditioning_diagnostic, conditioning_summary
from .types import (
    DEFAULT_GECP_CONFIG,
    CholeskyResult,
    GECPConfig,
    NumericalBreakdownError,
    NumericalMetadata,
)


def pivoted_cholesky(
    matrix: npt.ArrayLike, *, config: GECPConfig = DEFAULT_GECP_CONFIG
) -> CholeskyResult:
    """Factor a symmetric PSD matrix with lexicographic diagonal pivoting."""

    original = np.asarray(matrix, dtype=np.float64)
    if original.ndim != 2 or original.shape[0] != original.shape[1]:
        raise ValueError("matrix must be square")
    if not np.all(np.isfinite(original)):
        raise ValueError("matrix must contain only finite values")
    scale = max(1.0, float(np.max(np.abs(original), initial=0.0)))
    if not np.allclose(original, original.T, rtol=1e-12, atol=1e-14 * scale):
        raise ValueError("matrix must be symmetric")

    size = original.shape[0]
    residual = original / 2 + original.T / 2
    factor = np.zeros((size, min(size, config.max_rank)), dtype=np.float64)
    indices: list[int] = []
    pivots: list[float] = []
    history: list[float] = []
    conditioning = []
    tie_counts: list[int] = []
    tie_rtol = (
        config.tie_rtol if config.tie_rtol is not None else 32 * np.finfo(float).eps
    )
    converged = False
    stop_reason = "max_rank"
    rank_limit = min(size, config.max_rank)

    for rank in range(rank_limit + 1):
        diagonal = np.diag(residual)
        if not np.all(np.isfinite(residual)):
            raise NumericalBreakdownError("nonfinite Cholesky residual")
        if float(np.min(diagonal, initial=0.0)) < -100 * np.finfo(float).eps * scale:
            raise NumericalBreakdownError("residual has a materially negative diagonal")
        maximum = float(np.max(diagonal, initial=0.0))
        history.append(maximum)
        if maximum <= config.tol:
            if (
                float(np.max(np.abs(residual), initial=0.0))
                > config.tol + 100 * np.finfo(float).eps * scale
            ):
                raise NumericalBreakdownError(
                    "small diagonal does not bound residual; PSD assumption failed"
                )
            converged = True
            stop_reason = "tolerance"
            break
        if rank == rank_limit:
            break
        tolerance = tie_rtol * max(1.0, maximum)
        eligible = (diagonal > 0) & (np.abs(diagonal - maximum) <= tolerance)
        eligible[indices] = False
        candidates = np.flatnonzero(eligible)
        if not candidates.size:
            raise NumericalBreakdownError(
                "no usable unselected positive diagonal pivot"
            )
        tie_counts.append(int(candidates.size))
        pivot_index = int(candidates[0])
        pivot = float(residual[pivot_index, pivot_index])
        if pivot <= 0:
            raise NumericalBreakdownError("positive pivot required for Cholesky")
        column = residual[:, pivot_index] / np.sqrt(pivot)
        factor[:, rank] = column
        residual -= np.outer(column, column)
        residual = residual / 2 + residual.T / 2
        indices.append(pivot_index)
        pivots.append(pivot)
        conditioning.append(
            conditioning_diagnostic(original[np.ix_(indices, indices)], "selected_core")
        )

    returned_factor = factor[:, : len(indices)]
    residual = original - returned_factor @ returned_factor.T
    if not np.all(np.isfinite(residual)):
        raise NumericalBreakdownError("nonfinite returned Cholesky residual")
    returned_error = float(np.max(np.abs(residual), initial=0.0))
    if returned_error <= config.tol:
        converged = True
        stop_reason = "tolerance"
    elif converged:
        converged = False
        stop_reason = "roundoff_or_symmetry_limit"

    return CholeskyResult(
        indices=indices,
        pivots=pivots,
        residual_history=history,
        rank=len(indices),
        converged=converged,
        stop_reason=stop_reason,
        factor=returned_factor,
        residual_matrix=residual,
        metadata=NumericalMetadata(
            algorithm="pivoted_cholesky_grid",
            arithmetic="float64",
            precision_bits=53,
            requested_precision_bits=config.precision_bits,
            precision_status="as_requested"
            if config.precision_bits == 53
            else "float64_only",
            conditioning_status=conditioning_summary(conditioning),
            conditioning_method="float64_selected_core_svd",
            conditioning_history=conditioning,
            tie_status="recorded" if indices else "not_applicable",
            tie_policy="first_unselected_positive_index_within_tolerance",
            tie_rtol=tie_rtol,
            tie_tolerance_rule="tie_rtol*max(1,maximum_diagonal)",
            tie_counts=tie_counts,
            certification_scope="sampled_matrix_no_psd_or_rounding_certificate",
            convergence_scope="original_sampled_matrix_reconstruction_infinity_error",
            converged=converged,
            stop_reason=stop_reason,
            configuration={
                "requested_config": config.normalized(),
                "actual_pivot": "diagonal_grid",
                "rank_limit": rank_limit,
                "roundoff_slack": 100 * np.finfo(float).eps * scale,
                "returned_error": returned_error,
                "symmetry_error": float(
                    np.max(np.abs(original - original.T), initial=0.0)
                ),
            },
        ),
    )
