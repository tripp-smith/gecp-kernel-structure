"""Deterministic finite-grid Gaussian elimination with complete pivoting."""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import mpmath as mp
import numpy as np
import numpy.typing as npt

from .certified import interval_certified_pivot
from .grids import composite_chebyshev_grids
from .kernels import FermionicKernel
from .types import (
    DEFAULT_GECP_CONFIG,
    FloatArray,
    GECPConfig,
    GECPResult,
    HighPrecisionData,
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


def _mp_value(value: object) -> mp.mpf:
    """Convert numeric input through its deterministic decimal representation."""

    if isinstance(value, mp.mpf):
        return value
    if isinstance(value, (float, np.floating)):
        return mp.mpf(repr(float(value)))
    return mp.mpf(str(value))


def _mp_decimal(value: mp.mpf, precision_bits: int) -> str:
    digits = math.ceil(precision_bits * math.log10(2)) + 3
    return str(mp.nstr(value, n=digits, strip_zeros=False))


class _FermionicResidualEvaluator:
    """Point and interval evaluation of successive fermionic GECP residuals."""

    def __init__(self, cutoff: float) -> None:
        self.cutoff = cutoff
        self.t_nodes: list[mp.mpf] = []
        self.omega_nodes: list[mp.mpf] = []
        self.pivots: list[mp.mpf] = []
        self.pivot_columns: list[list[mp.mpf]] = []
        self.pivot_rows: list[list[mp.mpf]] = []

    @staticmethod
    def _point_kernel(t: mp.mpf, omega: mp.mpf) -> mp.mpf:
        if omega >= 0:
            return mp.exp(-t * omega) / (1 + mp.exp(-omega))
        return mp.exp((1 - t) * omega) / (1 + mp.exp(omega))

    @staticmethod
    def _point_kernel_with_gradient(
        t: mp.mpf, omega: mp.mpf
    ) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        value = _FermionicResidualEvaluator._point_kernel(t, omega)
        if omega >= 0:
            exp_negative = mp.exp(-omega)
            logistic = exp_negative / (1 + exp_negative)
        else:
            logistic = 1 / (1 + mp.exp(omega))
        return value, -omega * value, (logistic - t) * value

    @staticmethod
    def _interval_kernel(t: Any, omega: Any) -> Any:
        if hasattr(omega, "a") and hasattr(omega, "b"):
            lower = float(omega.a)
            upper = float(omega.b)
        else:
            lower = upper = float(omega)
        t_interval = t if hasattr(t, "a") else mp.iv.mpf(t)
        omega_interval = omega if hasattr(omega, "a") else mp.iv.mpf(omega)
        if lower >= 0:
            return mp.iv.exp(-t_interval * omega_interval) / (
                1 + mp.iv.exp(-omega_interval)
            )
        if upper <= 0:
            return mp.iv.exp((1 - t_interval) * omega_interval) / (
                1 + mp.iv.exp(omega_interval)
            )
        raise ValueError("interval frequency cell must not cross zero")

    @staticmethod
    def _interval_kernel_with_gradient(t: Any, omega: Any) -> tuple[Any, Any, Any]:
        value = _FermionicResidualEvaluator._interval_kernel(t, omega)
        t_interval = t if hasattr(t, "a") else mp.iv.mpf(t)
        omega_interval = omega if hasattr(omega, "a") else mp.iv.mpf(omega)
        lower = float(omega_interval.a)
        upper = float(omega_interval.b)
        if lower >= 0:
            exp_negative = mp.iv.exp(-omega_interval)
            logistic = exp_negative / (1 + exp_negative)
        elif upper <= 0:
            logistic = 1 / (1 + mp.iv.exp(omega_interval))
        else:
            raise ValueError("interval frequency cell must not cross zero")
        return (
            value,
            -omega_interval * value,
            (logistic - t_interval) * value,
        )

    def _evaluate(self, t: Any, omega: Any, *, interval: bool) -> Any:
        kernel = self._interval_kernel if interval else self._point_kernel
        value = kernel(t, omega)
        columns: list[Any] = []
        rows: list[Any] = []
        for j, pivot in enumerate(self.pivots):
            column = kernel(t, self.omega_nodes[j])
            row = kernel(self.t_nodes[j], omega)
            for earlier in range(j):
                column -= (
                    columns[earlier]
                    * self.pivot_rows[j][earlier]
                    / self.pivots[earlier]
                )
                row -= (
                    rows[earlier]
                    * self.pivot_columns[j][earlier]
                    / self.pivots[earlier]
                )
            columns.append(column)
            rows.append(row)
            value -= column * row / pivot
        return value

    def _evaluate_with_gradient(
        self, t: Any, omega: Any, *, interval: bool
    ) -> tuple[Any, Any, Any]:
        kernel = (
            self._interval_kernel_with_gradient
            if interval
            else self._point_kernel_with_gradient
        )
        value, derivative_t, derivative_omega = kernel(t, omega)
        columns: list[Any] = []
        column_derivatives: list[Any] = []
        rows: list[Any] = []
        row_derivatives: list[Any] = []
        for j, pivot in enumerate(self.pivots):
            column, column_dt, _ = kernel(t, self.omega_nodes[j])
            row, _, row_dw = kernel(self.t_nodes[j], omega)
            for earlier in range(j):
                row_coefficient = self.pivot_rows[j][earlier]
                column_coefficient = self.pivot_columns[j][earlier]
                column -= columns[earlier] * row_coefficient / self.pivots[earlier]
                column_dt -= (
                    column_derivatives[earlier] * row_coefficient / self.pivots[earlier]
                )
                row -= rows[earlier] * column_coefficient / self.pivots[earlier]
                row_dw -= (
                    row_derivatives[earlier] * column_coefficient / self.pivots[earlier]
                )
            columns.append(column)
            column_derivatives.append(column_dt)
            rows.append(row)
            row_derivatives.append(row_dw)
            value -= column * row / pivot
            derivative_t -= column_dt * row / pivot
            derivative_omega -= column * row_dw / pivot
        return value, derivative_t, derivative_omega

    def point(self, t: mp.mpf, omega: mp.mpf) -> mp.mpf:
        return self._evaluate(t, omega, interval=False)

    def interval(self, t: Any, omega: Any) -> Any:
        return self._evaluate(t, omega, interval=True)

    def interval_gradient(self, t: Any, omega: Any) -> tuple[Any, Any]:
        _, derivative_t, derivative_omega = self._evaluate_with_gradient(
            t, omega, interval=True
        )
        return derivative_t, derivative_omega

    def add_pivot(self, t: mp.mpf, omega: mp.mpf, pivot: mp.mpf) -> None:
        columns: list[mp.mpf] = []
        rows: list[mp.mpf] = []
        for j in range(len(self.pivots)):
            column = self._point_kernel(t, self.omega_nodes[j])
            row = self._point_kernel(self.t_nodes[j], omega)
            for earlier in range(j):
                column -= (
                    columns[earlier]
                    * self.pivot_rows[j][earlier]
                    / self.pivots[earlier]
                )
                row -= (
                    self.pivot_columns[j][earlier]
                    * rows[earlier]
                    / self.pivots[earlier]
                )
            columns.append(column)
            rows.append(row)
        self.t_nodes.append(+t)
        self.omega_nodes.append(+omega)
        self.pivots.append(+pivot)
        self.pivot_columns.append(columns)
        self.pivot_rows.append(rows)


def _gecp_adaptive_fermionic(
    kernel: FermionicKernel, *, config: GECPConfig
) -> GECPResult:
    """Run continuous fermionic GECP with an interval certificate per pivot."""

    evaluator = _FermionicResidualEvaluator(kernel.cutoff)
    pivots_mp: list[mp.mpf] = []
    residual_history_mp: list[mp.mpf] = []
    det_history_mp: list[mp.mpf] = []
    sigma_history_mp: list[mp.mpf] = []
    certificates = []
    converged = False
    stop_reason = "max_rank"

    with mp.workprec(config.precision_bits):
        for step in range(config.max_rank + 1):
            certificate = interval_certified_pivot(
                evaluator.point,
                evaluator.interval,
                omega_bounds=(-kernel.cutoff, kernel.cutoff),
                precision_bits=config.precision_bits,
                abs_tol=config.certificate_abs_tol,
                rel_tol=config.certificate_rel_tol,
                max_cells=config.max_cells,
                gradient_interval_function=evaluator.interval_gradient,
            )
            residual_history_mp.append(_mp_value(certificate.upper_bound))
            if not certificate.certified:
                stop_reason = "uncertified_pivot"
                break
            if certificate.upper_bound <= config.tol:
                converged = True
                stop_reason = "certified_tolerance"
                break
            if step == config.max_rank:
                stop_reason = "max_rank"
                break

            t = _mp_value(certificate.t)
            omega = _mp_value(certificate.omega)
            pivot = evaluator.point(t, omega)
            if not mp.isfinite(pivot) or pivot == 0:
                raise NumericalBreakdownError(
                    "unusable interval-certified pivot while residual is nonzero"
                )
            evaluator.add_pivot(t, omega, pivot)
            pivots_mp.append(+pivot)
            certificates.append(certificate)

            core = mp.matrix(
                [
                    [
                        evaluator._point_kernel(row_t, column_omega)
                        for column_omega in evaluator.omega_nodes
                    ]
                    for row_t in evaluator.t_nodes
                ]
            )
            det_history_mp.append(+mp.det(core))
            singular_values = mp.svd(core, compute_uv=False)
            sigma_history_mp.append(+singular_values[singular_values.rows - 1])

        def decimal(value: mp.mpf) -> str:
            return _mp_decimal(value, config.precision_bits)

        high_precision = HighPrecisionData(
            t_nodes=[decimal(value) for value in evaluator.t_nodes],
            omega_nodes=[decimal(value) for value in evaluator.omega_nodes],
            pivots=[decimal(value) for value in pivots_mp],
            residual_history=[decimal(value) for value in residual_history_mp],
            core_det_history=[decimal(value) for value in det_history_mp],
            core_sigma_min_history=[decimal(value) for value in sigma_history_mp],
        )
        return GECPResult(
            t_nodes=np.asarray(
                [float(value) for value in evaluator.t_nodes], dtype=np.float64
            ),
            omega_nodes=np.asarray(
                [float(value) for value in evaluator.omega_nodes], dtype=np.float64
            ),
            row_indices=[],
            column_indices=[],
            pivots=[float(value) for value in pivots_mp],
            residual_history=[float(value) for value in residual_history_mp],
            core_det_history=[float(value) for value in det_history_mp],
            core_sigma_min_history=[float(value) for value in sigma_history_mp],
            tie_counts=[1] * len(pivots_mp),
            pivot_certificates=certificates,
            rank=len(pivots_mp),
            converged=converged,
            stop_reason=stop_reason,
            config=config.normalized(),
            precision_bits=config.precision_bits,
            residual_matrix=None,
            high_precision=high_precision,
        )


def _gecp_mpmath_matrix(
    matrix: list[list[mp.mpf]],
    *,
    config: GECPConfig,
    t_grid: list[mp.mpf],
    omega_grid: list[mp.mpf],
) -> GECPResult:
    """Run deterministic finite-grid GECP entirely in mpmath arithmetic."""

    with mp.workprec(config.precision_bits):
        row_count = len(matrix)
        column_count = len(matrix[0]) if matrix else 0
        if row_count == 0 or column_count == 0:
            raise ValueError("matrix must be nonempty and two-dimensional")
        if any(len(row) != column_count for row in matrix):
            raise ValueError("matrix rows must have equal length")
        if len(t_grid) != row_count or len(omega_grid) != column_count:
            raise ValueError("grid dimensions must match the matrix")
        if any(not mp.isfinite(value) for row in matrix for value in row):
            raise ValueError("matrix must contain only finite values")

        original = [[+value for value in row] for row in matrix]
        residual = [[+value for value in row] for row in matrix]
        rows: list[int] = []
        columns: list[int] = []
        pivots_mp: list[mp.mpf] = []
        residual_history_mp: list[mp.mpf] = []
        det_history_mp: list[mp.mpf] = []
        sigma_history_mp: list[mp.mpf] = []
        tie_counts: list[int] = []
        converged = False
        stop_reason = "max_rank"
        rank_limit = min(config.max_rank, row_count, column_count)
        tolerance = _mp_value(config.tol)
        tie_rtol = (
            _mp_value(config.tie_rtol)
            if config.tie_rtol is not None
            else mp.power(2, 5 - config.precision_bits)
        )

        for _ in range(rank_limit):
            maximum = mp.mpf("-1")
            pivot_row = 0
            pivot_column = 0
            for i, residual_row in enumerate(residual):
                for j, value in enumerate(residual_row):
                    magnitude = abs(value)
                    if magnitude > maximum:
                        maximum = magnitude
                        pivot_row = i
                        pivot_column = j
            residual_history_mp.append(+maximum)
            if maximum <= tolerance:
                converged = True
                stop_reason = "tolerance"
                break

            threshold = tie_rtol * max(mp.mpf(1), maximum)
            near_ties = 0
            chosen = False
            for i, residual_row in enumerate(residual):
                for j, value in enumerate(residual_row):
                    if abs(abs(value) - maximum) <= threshold:
                        near_ties += 1
                        if not chosen:
                            pivot_row = i
                            pivot_column = j
                            chosen = True
            pivot = residual[pivot_row][pivot_column]
            if not mp.isfinite(pivot) or pivot == 0:
                raise NumericalBreakdownError(
                    "unusable arbitrary-precision pivot while residual is nonzero"
                )

            rows.append(pivot_row)
            columns.append(pivot_column)
            pivots_mp.append(+pivot)
            tie_counts.append(near_ties)
            pivot_column_values = [residual[i][pivot_column] for i in range(row_count)]
            pivot_row_values = list(residual[pivot_row])
            for i in range(row_count):
                for j in range(column_count):
                    residual[i][j] -= (
                        pivot_column_values[i] * pivot_row_values[j] / pivot
                    )

            core = mp.matrix([[original[i][j] for j in columns] for i in rows])
            det_history_mp.append(+mp.det(core))
            singular_values = mp.svd(core, compute_uv=False)
            sigma_history_mp.append(+singular_values[singular_values.rows - 1])
        else:
            final_maximum = max(abs(value) for row in residual for value in row)
            residual_history_mp.append(+final_maximum)
            if final_maximum <= tolerance:
                converged = True
                stop_reason = "tolerance"

        def decimal(value: mp.mpf) -> str:
            return _mp_decimal(value, config.precision_bits)

        high_precision = HighPrecisionData(
            t_nodes=[decimal(t_grid[i]) for i in rows],
            omega_nodes=[decimal(omega_grid[j]) for j in columns],
            pivots=[decimal(value) for value in pivots_mp],
            residual_history=[decimal(value) for value in residual_history_mp],
            core_det_history=[decimal(value) for value in det_history_mp],
            core_sigma_min_history=[decimal(value) for value in sigma_history_mp],
        )
        return GECPResult(
            t_nodes=np.asarray([float(t_grid[i]) for i in rows], dtype=np.float64),
            omega_nodes=np.asarray(
                [float(omega_grid[j]) for j in columns], dtype=np.float64
            ),
            row_indices=rows,
            column_indices=columns,
            pivots=[float(value) for value in pivots_mp],
            residual_history=[float(value) for value in residual_history_mp],
            core_det_history=[float(value) for value in det_history_mp],
            core_sigma_min_history=[float(value) for value in sigma_history_mp],
            tie_counts=tie_counts,
            rank=len(pivots_mp),
            converged=converged,
            stop_reason=stop_reason,
            config=config.normalized(),
            precision_bits=config.precision_bits,
            residual_matrix=np.asarray(
                [[float(value) for value in row] for row in residual],
                dtype=np.float64,
            ),
            high_precision=high_precision,
        )


def gecp_matrix(
    matrix: npt.ArrayLike,
    *,
    config: GECPConfig = DEFAULT_GECP_CONFIG,
    t_grid: npt.ArrayLike | None = None,
    omega_grid: npt.ArrayLike | None = None,
) -> GECPResult:
    """Run complete pivoting on a finite matrix using rank-one residual updates."""

    if config.precision_bits > 53:
        raw = np.asarray(matrix, dtype=object)
        if raw.ndim != 2 or min(raw.shape, default=0) == 0:
            raise ValueError("matrix must be nonempty and two-dimensional")
        raw_t = (
            np.arange(raw.shape[0], dtype=np.float64)
            if t_grid is None
            else np.asarray(t_grid, dtype=object)
        )
        raw_omega = (
            np.arange(raw.shape[1], dtype=np.float64)
            if omega_grid is None
            else np.asarray(omega_grid, dtype=object)
        )
        if raw_t.shape != (raw.shape[0],) or raw_omega.shape != (raw.shape[1],):
            raise ValueError("grid dimensions must match the matrix")
        with mp.workprec(config.precision_bits):
            matrix_mp = [
                [_mp_value(raw[i, j]) for j in range(raw.shape[1])]
                for i in range(raw.shape[0])
            ]
            return _gecp_mpmath_matrix(
                matrix_mp,
                config=config,
                t_grid=[_mp_value(value) for value in raw_t],
                omega_grid=[_mp_value(value) for value in raw_omega],
            )

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

    if config.pivot == "adaptive":
        if not isinstance(kernel, FermionicKernel):
            raise ValueError("adaptive GECP currently requires FermionicKernel")
        if t_grid is not None or omega_grid is not None:
            raise ValueError("adaptive GECP does not accept explicit grids")
        return _gecp_adaptive_fermionic(kernel, config=config)

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
    if config.precision_bits > 53:
        with mp.workprec(config.precision_bits):
            t_mp = [_mp_value(value) for value in t_values]
            omega_mp = [_mp_value(value) for value in omega_values]
            if isinstance(kernel, FermionicKernel):
                matrix_mp = [
                    [
                        _FermionicResidualEvaluator._point_kernel(t, omega)
                        for omega in omega_mp
                    ]
                    for t in t_mp
                ]
            else:
                matrix_mp = [
                    [_mp_value(kernel(t, omega)) for omega in omega_mp] for t in t_mp
                ]
            return _gecp_mpmath_matrix(
                matrix_mp,
                config=config,
                t_grid=t_mp,
                omega_grid=omega_mp,
            )
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
