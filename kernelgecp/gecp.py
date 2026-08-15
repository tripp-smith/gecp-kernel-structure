"""Deterministic finite-grid Gaussian elimination with complete pivoting."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt

from .grids import composite_chebyshev_grids
from .kernels import FermionicKernel
from .types import (
    DEFAULT_GECP_CONFIG,
    FloatArray,
    GECPConfig,
    GECPResult,
    NumericalBreakdownError,
)

KernelCallable = Callable[[npt.ArrayLike, npt.ArrayLike], npt.ArrayLike]


def _kernel_matrix(
    kernel: KernelCallable, t_grid: FloatArray, omega_grid: FloatArray
) -> FloatArray:
    values = np.asarray(kernel(t_grid[:, None], omega_grid[None, :]), dtype=np.float64)
    if values.shape != (t_grid.size, omega_grid.size):
        values = np.asarray(
            [[kernel(float(t), float(w)) for w in omega_grid] for t in t_grid],
            dtype=np.float64,
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("kernel evaluation produced a non-finite value")
    return values


def _tie_tolerance(config: GECPConfig) -> float:
    if config.tie_rtol is not None:
        return config.tie_rtol
    return 32.0 * np.finfo(np.float64).eps


def gecp_matrix(
    matrix: npt.ArrayLike,
    *,
    config: GECPConfig = DEFAULT_GECP_CONFIG,
    t_grid: npt.ArrayLike | None = None,
    omega_grid: npt.ArrayLike | None = None,
) -> GECPResult:
    """Run complete pivoting on a finite matrix using rank-one residual updates."""

    original = np.asarray(matrix, dtype=np.float64)
    if original.ndim != 2 or min(original.shape, default=0) == 0:
        raise ValueError("matrix must be nonempty and two-dimensional")
    if not np.all(np.isfinite(original)):
        raise ValueError("matrix must contain only finite values")
    t_nodes = (
        np.arange(original.shape[0], dtype=np.float64)
        if t_grid is None
        else np.asarray(t_grid, dtype=np.float64)
    )
    omega_nodes = (
        np.arange(original.shape[1], dtype=np.float64)
        if omega_grid is None
        else np.asarray(omega_grid, dtype=np.float64)
    )
    if t_nodes.shape != (original.shape[0],) or omega_nodes.shape != (
        original.shape[1],
    ):
        raise ValueError("grid dimensions must match the matrix")

    residual = original.copy()
    rows: list[int] = []
    columns: list[int] = []
    pivots: list[float] = []
    residual_history: list[float] = []
    det_history: list[float] = []
    sigma_history: list[float] = []
    tie_counts: list[int] = []
    converged = False
    stop_reason = "max_rank"
    rank_limit = min(config.max_rank, min(original.shape))
    tie_rtol = _tie_tolerance(config)

    for _ in range(rank_limit):
        magnitudes = np.abs(residual)
        maximum = float(np.max(magnitudes))
        residual_history.append(maximum)
        if maximum <= config.tol:
            converged = True
            stop_reason = "tolerance"
            break
        threshold = tie_rtol * max(1.0, maximum)
        near_ties = np.flatnonzero(np.abs(magnitudes.ravel() - maximum) <= threshold)
        flat_index = int(near_ties[0])
        row, column = np.unravel_index(flat_index, residual.shape)
        pivot = float(residual[row, column])
        if not np.isfinite(pivot) or abs(pivot) <= np.finfo(np.float64).tiny:
            raise NumericalBreakdownError(
                f"unusable pivot {pivot!r} while residual maximum is {maximum!r}"
            )

        rows.append(int(row))
        columns.append(int(column))
        pivots.append(pivot)
        tie_counts.append(int(near_ties.size))
        residual -= np.outer(residual[:, column], residual[row, :]) / pivot

        core = original[np.ix_(rows, columns)]
        det_history.append(float(np.linalg.det(core)))
        sigma_history.append(float(np.linalg.svd(core, compute_uv=False)[-1]))
    else:
        final_maximum = float(np.max(np.abs(residual)))
        residual_history.append(final_maximum)
        if final_maximum <= config.tol:
            converged = True
            stop_reason = "tolerance"

    return GECPResult(
        t_nodes=t_nodes[rows],
        omega_nodes=omega_nodes[columns],
        row_indices=rows,
        column_indices=columns,
        pivots=pivots,
        residual_history=residual_history,
        core_det_history=det_history,
        core_sigma_min_history=sigma_history,
        tie_counts=tie_counts,
        rank=len(pivots),
        converged=converged,
        stop_reason=stop_reason,
        config=config.normalized(),
        precision_bits=53,
        residual_matrix=residual,
    )


def gecp(
    kernel: KernelCallable,
    *,
    config: GECPConfig = DEFAULT_GECP_CONFIG,
    t_grid: npt.ArrayLike | None = None,
    omega_grid: npt.ArrayLike | None = None,
) -> GECPResult:
    """Run GECP on supplied grids or the canonical fermionic grids."""

    if t_grid is None or omega_grid is None:
        if not isinstance(kernel, FermionicKernel):
            raise ValueError(
                "explicit t_grid and omega_grid are required for this kernel"
            )
        default_t, default_omega = composite_chebyshev_grids(
            kernel.cutoff, config.grid_order
        )
        t_grid = default_t if t_grid is None else t_grid
        omega_grid = default_omega if omega_grid is None else omega_grid
    t_values = np.asarray(t_grid, dtype=np.float64)
    omega_values = np.asarray(omega_grid, dtype=np.float64)
    if t_values.ndim != 1 or omega_values.ndim != 1:
        raise ValueError("grids must be one-dimensional")
    matrix = _kernel_matrix(kernel, t_values, omega_values)
    return gecp_matrix(matrix, config=config, t_grid=t_values, omega_grid=omega_values)


def cross_approximation(
    kernel: KernelCallable,
    t_nodes: npt.ArrayLike,
    omega_nodes: npt.ArrayLike,
    t: npt.ArrayLike,
    omega: npt.ArrayLike,
) -> FloatArray:
    """Evaluate a selected cross by solving with its core matrix."""

    selected_t = np.asarray(t_nodes, dtype=np.float64)
    selected_omega = np.asarray(omega_nodes, dtype=np.float64)
    targets_t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    targets_omega = np.atleast_1d(np.asarray(omega, dtype=np.float64))
    if selected_t.size != selected_omega.size:
        raise ValueError("the selected row and column counts must agree")
    if selected_t.size == 0:
        return np.zeros((targets_t.size, targets_omega.size), dtype=np.float64)
    core = _kernel_matrix(kernel, selected_t, selected_omega)
    columns = _kernel_matrix(kernel, targets_t, selected_omega)
    rows = _kernel_matrix(kernel, selected_t, targets_omega)
    solved = np.linalg.solve(core, rows)
    return np.asarray(columns @ solved, dtype=np.float64)


def evaluate_residual(
    kernel: KernelCallable,
    result: GECPResult,
    t: npt.ArrayLike,
    omega: npt.ArrayLike,
) -> FloatArray:
    """Evaluate the original kernel minus the direct selected cross."""

    targets_t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    targets_omega = np.atleast_1d(np.asarray(omega, dtype=np.float64))
    original = _kernel_matrix(kernel, targets_t, targets_omega)
    return original - cross_approximation(
        kernel, result.t_nodes, result.omega_nodes, targets_t, targets_omega
    )
