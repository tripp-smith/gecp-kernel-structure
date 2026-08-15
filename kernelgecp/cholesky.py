"""Canonical diagonal-pivoted Cholesky for positive-semidefinite matrices."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .types import (
    DEFAULT_GECP_CONFIG,
    CholeskyResult,
    GECPConfig,
    NumericalBreakdownError,
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
    scale = max(1.0, float(np.max(np.abs(original))))
    if not np.allclose(original, original.T, rtol=1e-12, atol=1e-14 * scale):
        raise ValueError("matrix must be symmetric")

    size = original.shape[0]
    residual = original.copy()
    factor = np.zeros((size, min(size, config.max_rank)), dtype=np.float64)
    indices: list[int] = []
    pivots: list[float] = []
    history: list[float] = []
    converged = False
    stop_reason = "max_rank"
    rank_limit = min(size, config.max_rank)

    for rank in range(rank_limit):
        diagonal = np.diag(residual)
        maximum = float(np.max(diagonal))
        history.append(maximum)
        if maximum <= config.tol:
            if maximum < -100 * np.finfo(np.float64).eps * scale:
                raise NumericalBreakdownError(
                    "residual has a materially negative diagonal"
                )
            converged = True
            stop_reason = "tolerance"
            break
        tolerance = (config.tie_rtol or 32 * np.finfo(np.float64).eps) * max(
            1.0, maximum
        )
        candidates = np.flatnonzero(np.abs(diagonal - maximum) <= tolerance)
        pivot_index = int(candidates[0])
        pivot = float(residual[pivot_index, pivot_index])
        if pivot <= 0:
            raise NumericalBreakdownError("positive pivot required for Cholesky")
        column = residual[:, pivot_index] / np.sqrt(pivot)
        factor[:, rank] = column
        residual -= np.outer(column, column)
        residual = (residual + residual.T) / 2.0
        indices.append(pivot_index)
        pivots.append(pivot)
    else:
        final_maximum = float(np.max(np.diag(residual))) if size else 0.0
        history.append(final_maximum)
        if final_maximum <= config.tol:
            converged = True
            stop_reason = "tolerance"

    return CholeskyResult(
        indices=indices,
        pivots=pivots,
        residual_history=history,
        rank=len(indices),
        converged=converged,
        stop_reason=stop_reason,
        factor=factor[:, : len(indices)],
        residual_matrix=residual,
    )
