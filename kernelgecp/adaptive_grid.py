"""Residual-driven nested finite grids, with separate final sampled validation.

All maxima and hierarchical surpluses here are numerical sample diagnostics.
They are neither continuous supremum bounds nor interval certificates.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from dataclasses import asdict, dataclass
from decimal import (
    MAX_EMAX,
    MIN_EMIN,
    ROUND_HALF_EVEN,
    Context,
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    localcontext,
)
from typing import Any

import mpmath as mp
import numpy as np

from .adaptive_grids import (
    Panel,
    initial_panels,
    mark_panels,
    probe_nodes,
    refine_panels,
    training_nodes,
    validation_nodes,
)
from .gecp import gecp
from .kernels import FermionicKernel
from .types import FloatArray, GECPConfig, GECPResult, NumericalBreakdownError


def _decimal_context(bits: int = 192) -> Context:
    return Context(
        prec=max(80, math.ceil(bits * math.log10(2)) + 20),
        rounding=ROUND_HALF_EVEN,
        Emin=MIN_EMIN,
        Emax=MAX_EMAX,
        traps=[InvalidOperation, DivisionByZero, Overflow],
    )


@dataclass(frozen=True, slots=True)
class AdaptiveGridConfig:
    """Frozen budgets for a sampled refinement study, not a certificate request."""

    tolerances: tuple[float, ...] = (1e-6, 1e-8, 1e-10, 1e-12)
    precision_bits: int = 128
    check_precision_bits: int = 192
    generation_margin: float = 0.25
    max_rank: int = 256
    max_stages: int = 5
    max_nodes: int = 512
    max_marked: int = 4

    def __post_init__(self) -> None:
        for name, minimum in (
            ("precision_bits", 54),
            ("check_precision_bits", 55),
            ("max_rank", 0),
            ("max_stages", 1),
            ("max_nodes", 2),
            ("max_marked", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if self.check_precision_bits <= self.precision_bits:
            raise ValueError("check precision must exceed elimination precision")
        if (
            not isinstance(self.tolerances, tuple)
            or not self.tolerances
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not np.isfinite(value)
                or not 0 < value < 1
                for value in self.tolerances
            )
            or len(set(self.tolerances)) != len(self.tolerances)
            or tuple(sorted(self.tolerances, reverse=True)) != self.tolerances
        ):
            raise ValueError("tolerances must be distinct decreasing values in (0,1)")
        if (
            isinstance(self.generation_margin, bool)
            or not isinstance(self.generation_margin, (int, float))
            or not np.isfinite(self.generation_margin)
            or not 0 < self.generation_margin < 1
            or self.generation_tolerance <= 0
        ):
            raise ValueError(
                "generation_margin must give a positive stricter tolerance"
            )

    @property
    def generation_tolerance(self) -> float:
        with localcontext(_decimal_context()):
            return float(
                Decimal(str(min(self.tolerances)))
                * Decimal(str(self.generation_margin))
            )


DEFAULT_ADAPTIVE_GRID_CONFIG = AdaptiveGridConfig()


@dataclass(slots=True)
class AdaptiveGridCase:
    """One actual refinement trajectory and its final tolerance-prefix checks."""

    cutoff: float
    config: AdaptiveGridConfig
    stages: list[dict[str, Any]]
    final_result: GECPResult
    t_grid: FloatArray
    omega_grid: FloatArray
    t_panels: list[Panel]
    omega_panels: list[Panel]
    prefixes: list[dict[str, Any]]
    converged: bool
    stop_reason: str


def _coordinates(values: FloatArray) -> list[str]:
    return [repr(float(value)) for value in values]


def _encode(value: Any, bits: int) -> str:
    if not mp.isfinite(value):
        raise NumericalBreakdownError("nonfinite adaptive-grid diagnostic")
    return str(mp.nstr(value, math.ceil(bits * math.log10(2)) + 3))


def _kernel(t: Any, omega: Any) -> Any:
    """Independent sign-stable expression; never call the elimination evaluator."""
    if omega >= 0:
        return mp.exp(-t * omega) / (1 + mp.exp(-omega))
    return mp.exp((1 - t) * omega) / (1 + mp.exp(omega))


def _lu_solve(core: Any, right: Any) -> Any:
    """One LU factorization and batched forward/back substitution, no inverse."""
    permutation, lower, upper = mp.lu(core)
    size = core.rows
    forward = mp.matrix(size, right.cols)
    solution = mp.matrix(size, right.cols)
    ordering = [next(j for j in range(size) if permutation[i, j]) for i in range(size)]
    for i in range(size):
        for column in range(right.cols):
            forward[i, column] = (
                right[ordering[i], column]
                - mp.fsum(lower[i, j] * forward[j, column] for j in range(i))
            ) / lower[i, i]
    for i in range(size - 1, -1, -1):
        for column in range(right.cols):
            solution[i, column] = (
                forward[i, column]
                - mp.fsum(upper[i, j] * solution[j, column] for j in range(i + 1, size))
            ) / upper[i, i]
    return solution


def _residual_matrix(
    result: GECPResult,
    rank: int,
    times: FloatArray,
    frequencies: FloatArray,
    bits: int,
) -> tuple[Any, dict[str, Any]]:
    """Reconstruct the selected cross at higher precision on explicit axes.

    Coordinates are interpreted through repr(float), as in grid GECP. Selected
    nodes use its canonical decimal strings. This is point arithmetic, not an
    enclosure of rounding error or an interval residual model.
    """
    high = result.high_precision
    if high is None or not 0 <= rank <= result.rank:
        raise ValueError("a valid arbitrary-precision trajectory prefix is required")
    if bits <= result.precision_bits:
        raise ValueError("cross-check precision must exceed source precision")
    with mp.workprec(bits):
        t_values = [mp.mpf(value) for value in _coordinates(times)]
        omega_values = [mp.mpf(value) for value in _coordinates(frequencies)]
        if rank:
            selected_t = [mp.mpf(value) for value in high.t_nodes[:rank]]
            selected_w = [mp.mpf(value) for value in high.omega_nodes[:rank]]
            core = mp.matrix([[_kernel(t, w) for w in selected_w] for t in selected_t])
            right = mp.matrix(
                [[_kernel(t, w) for w in omega_values] for t in selected_t]
            )
            try:
                solution = _lu_solve(core, right)
            except (ZeroDivisionError, ValueError) as exc:
                raise NumericalBreakdownError(
                    "adaptive-grid selected-core solve failed"
                ) from exc
            left = mp.matrix([[_kernel(t, w) for w in selected_w] for t in t_values])
            approximation = left * solution
            # Bound the diagnostic cost: this is an observed solve residual on
            # at most three RHS columns, explicitly not an all-column bound.
            columns = sorted({0, right.cols // 2, right.cols - 1})
            numerator = max(
                abs(
                    mp.fsum(core[i, j] * solution[j, c] for j in range(rank))
                    - right[i, c]
                )
                for i in range(rank)
                for c in columns
            )
            denominator = max(abs(right[i, c]) for i in range(rank) for c in columns)
            solve_residual = (
                _encode(numerator / denominator, bits) if denominator else "0.0"
            )
        else:
            approximation = mp.matrix(len(t_values), len(omega_values))
            columns = []
            solve_residual = None
        residual = mp.matrix(
            [
                [
                    _kernel(t, w) - approximation[i, j]
                    for j, w in enumerate(omega_values)
                ]
                for i, t in enumerate(t_values)
            ]
        )
        if any(not mp.isfinite(value) for value in residual):
            raise NumericalBreakdownError("nonfinite adaptive-grid cross residual")
    return residual, {
        "solve_method": "one_mpmath_lu_batched_triangular_solves"
        if rank
        else "rank_zero",
        "factorizations": 1 if rank else 0,
        "solve_residual": solve_residual,
        "solve_residual_scope": "relative_max_entry_on_at_most_three_rhs_columns",
        "solve_residual_rhs_indices": columns,
    }


def _training_max(result: GECPResult, rank: int) -> str:
    if result.high_precision is None:
        raise ValueError("canonical high-precision history is required")
    return result.high_precision.residual_history[rank]


def _probe_metrics(
    result: GECPResult,
    rank: int,
    t_grid: FloatArray,
    omega_grid: FloatArray,
    t_panels: list[Panel],
    omega_panels: list[Panel],
    bits: int,
) -> dict[str, Any]:
    """Measure probes and interpolation surpluses once for this actual prefix."""
    probe_t, probe_w = probe_nodes(t_panels), probe_nodes(omega_panels)
    if np.intersect1d(probe_t, t_grid).size or np.intersect1d(probe_w, omega_grid).size:
        raise NumericalBreakdownError(
            "marking probes coincide with training coordinates"
        )
    union_t, union_w = np.union1d(t_grid, probe_t), np.union1d(omega_grid, probe_w)
    residual, solve = _residual_matrix(result, rank, union_t, union_w, bits)
    index_t = {float(value): i for i, value in enumerate(union_t)}
    index_w = {float(value): i for i, value in enumerate(union_w)}
    metrics: dict[str, Any] = {
        "rank": rank,
        "precision_bits": bits,
        "t_nodes": _coordinates(probe_t),
        "omega_nodes": _coordinates(probe_w),
        "coordinate_precision_bits": 53,
        "coordinate_semantics": "repr_float64_decimal",
        "certification_status": "not_certified",
        "scope": "marking_and_stopping_samples_not_independent_final_validation",
        **solve,
    }
    with mp.workprec(bits):
        probe_max = max(
            abs(residual[index_t[float(t)], index_w[float(w)]])
            for t in probe_t
            for w in probe_w
        )
        surplus_max = mp.mpf(0)
        for axis, panels, probes, grid, index, other, other_index in (
            ("t", t_panels, probe_t, t_grid, index_t, probe_w, index_w),
            ("omega", omega_panels, probe_w, omega_grid, index_w, probe_t, index_t),
        ):
            panel_metrics = {}
            for panel in panels:
                peak, surplus = mp.mpf(0), mp.mpf(0)
                for point in probes:
                    if not panel.left < point < panel.right:
                        continue
                    position = bisect_left(grid, point)
                    low, high = float(grid[position - 1]), float(grid[position])
                    fraction = (mp.mpf(repr(float(point))) - mp.mpf(repr(low))) / (
                        mp.mpf(repr(high)) - mp.mpf(repr(low))
                    )
                    p, lo, hi = index[float(point)], index[low], index[high]
                    for coordinate in other:
                        j = other_index[float(coordinate)]
                        if axis == "t":
                            value = residual[p, j]
                            interpolated = (1 - fraction) * residual[
                                lo, j
                            ] + fraction * residual[hi, j]
                        else:
                            value = residual[j, p]
                            interpolated = (1 - fraction) * residual[
                                j, lo
                            ] + fraction * residual[j, hi]
                        peak = max(peak, abs(value))
                        surplus = max(surplus, abs(value - interpolated))
                panel_metrics[panel.identifier] = {
                    "probe_max": _encode(peak, bits),
                    "surplus": _encode(surplus, bits),
                }
                surplus_max = max(surplus_max, surplus)
            metrics[f"{axis}_scores"] = panel_metrics
        metrics["probe_max"] = _encode(probe_max, bits)
        metrics["surplus_max"] = _encode(surplus_max, bits)
    return metrics


def _classify_probe(
    metrics: dict[str, Any], tolerance: float, training_max: str
) -> dict[str, Any]:
    threshold = Decimal(str(tolerance))
    classified = dict(metrics)
    with localcontext(_decimal_context(metrics["precision_bits"])):
        for axis in ("t", "omega"):
            classified[f"{axis}_scores"] = {
                identifier: {
                    **values,
                    "score": str(
                        max(
                            Decimal(values["probe_max"]) / threshold,
                            4 * Decimal(values["surplus"]) / threshold,
                        )
                    ),
                }
                for identifier, values in metrics[f"{axis}_scores"].items()
            }
        training_passed = Decimal(training_max) <= threshold
        probe_passed = Decimal(metrics["probe_max"]) <= threshold
        surplus_passed = Decimal(metrics["surplus_max"]) <= threshold / 4
    classified.update(
        tolerance=str(tolerance),
        training_max=training_max,
        training_passed=training_passed,
        probe_passed=probe_passed,
        surplus_passed=surplus_passed,
        passed=training_passed and probe_passed and surplus_passed,
    )
    return classified


def _validation_metrics(
    result: GECPResult, rank: int, times: FloatArray, frequencies: FloatArray, bits: int
) -> dict[str, Any]:
    residual, solve = _residual_matrix(result, rank, times, frequencies, bits)
    with mp.workprec(bits):
        maximum = _encode(max(abs(value) for value in residual), bits)
    return {
        "rank": rank,
        "precision_bits": bits,
        "t_nodes": _coordinates(times),
        "omega_nodes": _coordinates(frequencies),
        "max_residual": maximum,
        "coordinate_precision_bits": 53,
        "coordinate_semantics": "repr_float64_decimal",
        "scope": "untouched_final_seventh_fraction_samples_after_refinement_stops",
        "certification_status": "not_certified",
        **solve,
    }


def _final_prefixes(
    result: GECPResult,
    t_grid: FloatArray,
    omega_grid: FloatArray,
    t_panels: list[Panel],
    omega_panels: list[Panel],
    config: AdaptiveGridConfig,
    cache: dict[int, dict[str, Any]],
    controller_passed: bool,
    controller_stop_reason: str,
) -> list[dict[str, Any]]:
    prefixes = []
    for tolerance in config.tolerances:
        threshold = Decimal(str(tolerance))
        training_rank = next(
            (
                rank
                for rank in range(result.rank + 1)
                if Decimal(_training_max(result, rank)) <= threshold
            ),
            None,
        )
        tested = []
        chosen = result.rank
        check: dict[str, Any] | None = None
        for rank in range(
            training_rank if training_rank is not None else result.rank, result.rank + 1
        ):
            if rank not in cache:
                cache[rank] = _probe_metrics(
                    result,
                    rank,
                    t_grid,
                    omega_grid,
                    t_panels,
                    omega_panels,
                    config.check_precision_bits,
                )
            check = _classify_probe(cache[rank], tolerance, _training_max(result, rank))
            tested.append(rank)
            chosen = rank
            if check["passed"]:
                break
        assert check is not None
        prefixes.append(
            {
                "tolerance": str(tolerance),
                "training_rank": training_rank,
                "chosen_rank": chosen,
                "tested_ranks": tested,
                "training_max": check["training_max"],
                "probe_max": check["probe_max"],
                "surplus_max": check["surplus_max"],
                "training_passed": check["training_passed"],
                "probe_passed": check["probe_passed"],
                "surplus_passed": check["surplus_passed"],
                "joint_passed": check["passed"],
                "probe": check,
                "fallback_to_deepest": chosen == result.rank
                and chosen != training_rank,
                "prefix_budget_exhausted": not check["passed"],
                "rank_scope": "first_jointly_passing_tested_prefix_not_global_minimum",
            }
        )
    # This call and every access to validation coordinates occur after all
    # controller decisions and all prefix choices are complete.
    final_t, final_w = validation_nodes(t_panels), validation_nodes(omega_panels)
    for validation_axis, training, marking in (
        (final_t, t_grid, probe_nodes(t_panels)),
        (final_w, omega_grid, probe_nodes(omega_panels)),
    ):
        if (
            np.intersect1d(validation_axis, training).size
            or np.intersect1d(validation_axis, marking).size
        ):
            raise NumericalBreakdownError(
                "final validation coordinates are not disjoint"
            )
    validation_cache = {}
    for prefix in prefixes:
        rank = prefix["chosen_rank"]
        if rank not in validation_cache:
            validation_cache[rank] = _validation_metrics(
                result, rank, final_t, final_w, config.check_precision_bits
            )
        validation = {
            **validation_cache[rank],
            "passed": Decimal(validation_cache[rank]["max_residual"])
            <= Decimal(prefix["tolerance"]),
        }
        fit_passed = bool(prefix["joint_passed"] and validation["passed"])
        converged = controller_passed and fit_passed
        prefix.update(
            validation=validation,
            fit_passed=fit_passed,
            controller_passed=controller_passed,
            converged=converged,
            stop_reason=(
                controller_stop_reason
                if not controller_passed
                else "sampled_tolerance"
                if converged
                else "validation_failed"
                if prefix["joint_passed"]
                else "prefix_budget"
            ),
        )
    return prefixes


def run_adaptive_grid(
    cutoff: float, *, config: AdaptiveGridConfig = DEFAULT_ADAPTIVE_GRID_CONFIG
) -> AdaptiveGridCase:
    """Run genuine nested sampled refinement; preserve budget/validation failures."""
    if isinstance(cutoff, bool) or not np.isfinite(cutoff) or cutoff < 1:
        raise ValueError("cutoff must be finite and at least one")
    t_panels, omega_panels = initial_panels(float(cutoff))
    t_grid, omega_grid = training_nodes(t_panels), training_nodes(omega_panels)
    if max(len(t_grid), len(omega_grid)) > config.max_nodes:
        raise ValueError("initial grid exceeds max_nodes; no legal stage can be run")
    strictest = min(config.tolerances)
    stages = []
    consecutive = 0
    controller_passed = False
    stop_reason = "max_stages"
    cache: dict[int, dict[str, Any]] = {}
    previous_t: FloatArray | None = None
    previous_w: FloatArray | None = None
    for index in range(config.max_stages):
        result = gecp(
            FermionicKernel(float(cutoff)),
            config=GECPConfig(
                tol=config.generation_tolerance,
                max_rank=config.max_rank,
                pivot="grid",
                precision_bits=config.precision_bits,
            ),
            t_grid=t_grid,
            omega_grid=omega_grid,
        )
        cache = {
            result.rank: _probe_metrics(
                result,
                result.rank,
                t_grid,
                omega_grid,
                t_panels,
                omega_panels,
                config.check_precision_bits,
            )
        }
        check = _classify_probe(
            cache[result.rank], strictest, _training_max(result, result.rank)
        )
        stage_passed = bool(result.converged and check["passed"])
        consecutive = consecutive + 1 if stage_passed else 0
        stage = {
            "index": index,
            "t_grid": _coordinates(t_grid),
            "omega_grid": _coordinates(omega_grid),
            "t_panels": [asdict(panel) for panel in t_panels],
            "omega_panels": [asdict(panel) for panel in omega_panels],
            "grid_shape": [len(t_grid), len(omega_grid)],
            "coordinate_precision_bits": 53,
            "coordinate_semantics": "repr_float64_decimal",
            "nested": previous_t is None
            or (
                previous_w is not None
                and set(previous_t).issubset(t_grid)
                and set(previous_w).issubset(omega_grid)
            ),
            "gecp": {
                "rank": result.rank,
                "converged": result.converged,
                "stop_reason": result.stop_reason,
                "high_precision": asdict(result.high_precision)
                if result.high_precision is not None
                else None,
                "tie_counts": result.tie_counts,
                "metadata": result.metadata.as_dict(),
            },
            "probe": check,
            "passed": stage_passed,
            "consecutive_passes": consecutive,
            "marked_t": [],
            "marked_omega": [],
            "refinement_applied": False,
            "budgets": asdict(config),
            "generation_tolerance": str(config.generation_tolerance),
            "certification_status": "not_certified",
            "stop_reason": "refine",
        }
        stages.append(stage)
        if index >= 1 and consecutive >= 2:
            controller_passed, stop_reason = True, "strictest_sampled_stages_passed"
            break
        if not result.converged:
            stop_reason = result.stop_reason
            break
        if index + 1 == config.max_stages:
            stop_reason = "max_stages"
            break
        marked_t = mark_panels(
            t_panels,
            {
                identifier: Decimal(values["score"])
                for identifier, values in check["t_scores"].items()
            },
            axis="t",
            max_marked=config.max_marked,
        )
        marked_w = mark_panels(
            omega_panels,
            {
                identifier: Decimal(values["score"])
                for identifier, values in check["omega_scores"].items()
            },
            axis="omega",
            max_marked=config.max_marked,
        )
        stage.update(marked_t=marked_t, marked_omega=marked_w)
        if not marked_t and not marked_w:
            stop_reason = "no_marked_panels"
            break
        try:
            next_t_panels = refine_panels(t_panels, marked_t)
            next_w_panels = refine_panels(omega_panels, marked_w)
            next_t = training_nodes(next_t_panels, previous=t_grid)
            next_w = training_nodes(next_w_panels, previous=omega_grid)
        except ValueError as exc:
            stop_reason = "coordinate_resolution"
            stage["failure"] = str(exc)
            break
        stage["proposed_grid_shape"] = [len(next_t), len(next_w)]
        if max(len(next_t), len(next_w)) > config.max_nodes:
            stop_reason = "max_nodes"
            break
        if len(next_t) == len(t_grid) and len(next_w) == len(omega_grid):
            stop_reason = "coordinate_resolution"
            break
        stage["refinement_applied"] = True
        previous_t, previous_w = t_grid, omega_grid
        t_grid, omega_grid = next_t, next_w
        t_panels, omega_panels = next_t_panels, next_w_panels
    stages[-1]["stop_reason"] = stop_reason
    prefixes = _final_prefixes(
        result,
        t_grid,
        omega_grid,
        t_panels,
        omega_panels,
        config,
        cache,
        controller_passed,
        stop_reason,
    )
    converged = controller_passed and all(prefix["converged"] for prefix in prefixes)
    if controller_passed:
        stop_reason = "sampled_tolerances" if converged else "validation_failed"
    return AdaptiveGridCase(
        float(cutoff),
        config,
        stages,
        result,
        t_grid,
        omega_grid,
        t_panels,
        omega_panels,
        prefixes,
        converged,
        stop_reason,
    )
