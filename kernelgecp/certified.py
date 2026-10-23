"""Deterministic branch-and-bound certificates for approximate pivots."""

from __future__ import annotations

import heapq
import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import mpmath as mp
import numpy as np

from .types import PivotCertificate

ScalarKernel = Callable[[float, float], float]
LipschitzBound = (
    float | tuple[float, float] | Callable[[float, float], tuple[float, float]]
)
MPPointFunction = Callable[[mp.mpf, mp.mpf], mp.mpf]
MPIntervalFunction = Callable[[Any, Any], Any]
MPGradientIntervalFunction = Callable[[Any, Any], tuple[Any, Any]]


@dataclass(order=True, slots=True)
class _Cell:
    priority: float
    serial: int
    t0: float = field(compare=False)
    t1: float = field(compare=False)
    w0: float = field(compare=False)
    w1: float = field(compare=False)
    value: float = field(compare=False)
    upper: float = field(compare=False)


@dataclass(order=True, slots=True)
class _IntervalCell:
    priority: float
    serial: int
    t0: float = field(compare=False)
    t1: float = field(compare=False)
    w0: float = field(compare=False)
    w1: float = field(compare=False)
    point_value: mp.mpf = field(compare=False)
    point_lower: float = field(compare=False)
    upper: float = field(compare=False)


def _lipschitz_at(bound: LipschitzBound, t: float, omega: float) -> tuple[float, float]:
    if callable(bound):
        lt, lw = bound(t, omega)
    elif isinstance(bound, tuple):
        lt, lw = bound
    else:
        lt = lw = bound
    if not np.isfinite(lt) or not np.isfinite(lw) or lt < 0 or lw < 0:
        raise ValueError("Lipschitz bounds must be finite and nonnegative")
    return float(lt), float(lw)


def _make_cell(
    function: ScalarKernel,
    bound: LipschitzBound,
    coordinates: tuple[float, float, float, float],
    serial: int,
) -> _Cell:
    t0, t1, w0, w1 = coordinates
    t = (t0 + t1) / 2.0
    omega = (w0 + w1) / 2.0
    value = float(function(t, omega))
    if not np.isfinite(value):
        raise ValueError("pivot objective returned a non-finite value")
    lt, lw = _lipschitz_at(bound, t, omega)
    upper = abs(value) + lt * (t1 - t0) / 2.0 + lw * (w1 - w0) / 2.0
    return _Cell(-upper, serial, t0, t1, w0, w1, value, upper)


def certified_pivot(
    function: ScalarKernel,
    *,
    t_bounds: tuple[float, float] = (0.0, 1.0),
    omega_bounds: tuple[float, float],
    lipschitz: LipschitzBound,
    abs_tol: float = 1e-12,
    rel_tol: float = 1e-8,
    max_cells: int = 100_000,
) -> PivotCertificate:
    """Certify an approximate absolute maximizer from a valid Lipschitz bound.

    The function assumes the supplied coordinatewise Lipschitz constants are
    rigorous. Its certificate is conditional on that analytic bound.
    """

    t0, t1 = t_bounds
    w0, w1 = omega_bounds
    if not all(np.isfinite((t0, t1, w0, w1))) or t0 >= t1 or w0 >= w1:
        raise ValueError("bounds must be finite, ordered intervals")
    if abs_tol < 0 or rel_tol < 0 or max_cells < 1:
        raise ValueError("invalid certificate tolerance or budget")

    serial = 0
    first = _make_cell(function, lipschitz, (t0, t1, w0, w1), serial)
    queue = [first]
    evaluated = 1
    best_value = abs(first.value)
    best_point = ((t0 + t1) / 2.0, (w0 + w1) / 2.0)

    while queue:
        global_upper = -queue[0].priority
        tolerance = abs_tol + rel_tol * best_value
        if global_upper - best_value <= tolerance:
            eta = 1.0 if global_upper == 0 else min(1.0, best_value / global_upper)
            return PivotCertificate(
                t=best_point[0],
                omega=best_point[1],
                value=float(function(*best_point)),
                lower_bound=best_value,
                upper_bound=global_upper,
                eta=eta,
                cells_evaluated=evaluated,
                certified=True,
                termination_reason="gap_tolerance",
            )
        if evaluated >= max_cells:
            break

        if evaluated + 2 > max_cells:
            break

        cell = heapq.heappop(queue)
        t_mid = (cell.t0 + cell.t1) / 2.0
        w_mid = (cell.w0 + cell.w1) / 2.0
        split_t = (cell.t1 - cell.t0) >= (cell.w1 - cell.w0)
        if split_t:
            children = (
                (cell.t0, t_mid, cell.w0, cell.w1),
                (t_mid, cell.t1, cell.w0, cell.w1),
            )
        else:
            children = (
                (cell.t0, cell.t1, cell.w0, w_mid),
                (cell.t0, cell.t1, w_mid, cell.w1),
            )
        for coordinates in children:
            serial += 1
            child = _make_cell(function, lipschitz, coordinates, serial)
            evaluated += 1
            heapq.heappush(queue, child)
            child_abs = abs(child.value)
            point = (
                (child.t0 + child.t1) / 2.0,
                (child.w0 + child.w1) / 2.0,
            )
            if child_abs > best_value or (
                child_abs == best_value and point < best_point
            ):
                best_value = child_abs
                best_point = point

    global_upper = max((-cell.priority for cell in queue), default=best_value)
    eta = 1.0 if global_upper == 0 else min(1.0, best_value / global_upper)
    return PivotCertificate(
        t=best_point[0],
        omega=best_point[1],
        value=float(function(*best_point)),
        lower_bound=best_value,
        upper_bound=global_upper,
        eta=eta,
        cells_evaluated=evaluated,
        certified=False,
        termination_reason="cell_budget",
    )


def _outward_upper(value: float) -> float:
    if value == 0:
        return 0.0
    return float(np.nextafter(value, math.inf))


def _outward_lower(value: float) -> float:
    return max(0.0, float(np.nextafter(value, -math.inf)))


def _midpoint(lower: float, upper: float) -> float:
    """Avoid overflow in the midpoint; enclosure radii use the actual result."""
    return lower / 2.0 + upper / 2.0


def _interval_absolute_bounds(enclosure: Any) -> tuple[float, float]:
    """Convert a valid interval to conservative float64 absolute bounds."""
    lower_endpoint = float(enclosure.a)
    upper_endpoint = float(enclosure.b)
    if (
        not np.isfinite(lower_endpoint)
        or not np.isfinite(upper_endpoint)
        or lower_endpoint > upper_endpoint
    ):
        raise ValueError("interval pivot objective produced a non-finite enclosure")
    lower = (
        0.0
        if lower_endpoint <= 0 <= upper_endpoint
        else _outward_lower(min(abs(lower_endpoint), abs(upper_endpoint)))
    )
    # A nonzero interval endpoint may underflow when converted to float. Only
    # an exactly zero interval may have a zero upper bound; otherwise one
    # outward float step also covers values below the minimum subnormal.
    upper = (
        0.0
        if enclosure.a == 0 and enclosure.b == 0
        else float(
            np.nextafter(max(abs(lower_endpoint), abs(upper_endpoint)), math.inf)
        )
    )
    if not np.isfinite(upper):
        raise ValueError("interval pivot bound exceeds the finite float64 range")
    return lower, upper


def _interval_point(
    point_function: MPPointFunction,
    interval_function: MPIntervalFunction,
    t: float,
    omega: float,
    precision_bits: int,
) -> tuple[mp.mpf, float, Any]:
    """An approximate point value is never used as a certified lower bound."""
    with mp.workprec(precision_bits):
        point_value = point_function(mp.mpf(t), mp.mpf(omega))
    if not mp.isfinite(point_value) or not np.isfinite(float(point_value)):
        raise ValueError("pivot objective returned a non-finite point value")
    enclosure = interval_function(mp.iv.mpf(t), mp.iv.mpf(omega))
    lower, _ = _interval_absolute_bounds(enclosure)
    return point_value, lower, enclosure


def _make_interval_cell(
    point_function: MPPointFunction,
    interval_function: MPIntervalFunction,
    coordinates: tuple[float, float, float, float],
    serial: int,
    precision_bits: int,
    gradient_interval_function: MPGradientIntervalFunction | None,
) -> _IntervalCell:
    t0, t1, w0, w1 = coordinates
    t, omega = _midpoint(t0, t1), _midpoint(w0, w1)
    point_value, point_lower, center_enclosure = _interval_point(
        point_function, interval_function, t, omega, precision_bits
    )
    t_box, omega_box = mp.iv.mpf([t0, t1]), mp.iv.mpf([w0, w1])
    enclosure = interval_function(t_box, omega_box)
    _, upper = _interval_absolute_bounds(enclosure)
    if gradient_interval_function is not None:
        dt_enclosure, dw_enclosure = gradient_interval_function(t_box, omega_box)
        _interval_absolute_bounds(dt_enclosure)
        _interval_absolute_bounds(dw_enclosure)
        # The rounded midpoint need not be the exact center of the box. All
        # radius, product and sum operations remain in interval arithmetic.
        mean_enclosure = (
            abs(center_enclosure)
            + abs(dt_enclosure) * abs(t_box - mp.iv.mpf(t))
            + abs(dw_enclosure) * abs(omega_box - mp.iv.mpf(omega))
        )
        _, mean_upper = _interval_absolute_bounds(mean_enclosure)
        upper = min(upper, mean_upper)
    if point_lower > upper:
        raise ValueError("singleton lower bound exceeds the cell upper bound")
    return _IntervalCell(
        -upper,
        serial,
        t0,
        t1,
        w0,
        w1,
        point_value,
        point_lower,
        upper,
    )


def interval_certified_pivot(
    point_function: MPPointFunction,
    interval_function: MPIntervalFunction,
    *,
    t_bounds: tuple[float, float] = (0.0, 1.0),
    omega_bounds: tuple[float, float],
    precision_bits: int = 128,
    abs_tol: float = 1e-12,
    rel_tol: float = 1e-8,
    max_cells: int = 100_000,
    gradient_interval_function: MPGradientIntervalFunction | None = None,
    tie_rtol: float | None = None,
    objective_scope: str = "caller_supplied_interval_function",
) -> PivotCertificate:
    """Bound the supplied interval objective and select a sampled lower bound.

    ``interval_function`` and the optional gradient must enclose the same fixed
    mathematical function over each box. ``point_function`` supplies only the
    approximate compatibility value. In adaptive GECP the fixed function is K
    minus the approximant defined by frozen rounded elimination coefficients;
    this does not certify an exact ideal GECP residual or a fresh float64 cross.
    Tie counts describe distinct observed singleton lower bounds only.
    """

    t0, t1 = t_bounds
    w0, w1 = omega_bounds
    if not all(np.isfinite((t0, t1, w0, w1))) or t0 >= t1 or w0 >= w1:
        raise ValueError("bounds must be finite, ordered intervals")
    if (
        precision_bits < 53
        or not np.isfinite(abs_tol)
        or not np.isfinite(rel_tol)
        or abs_tol < 0
        or rel_tol < 0
        or max_cells < 1
    ):
        raise ValueError("invalid interval certificate precision, tolerance, or budget")
    if tie_rtol is not None and (not np.isfinite(tie_rtol) or tie_rtol < 0):
        raise ValueError("tie_rtol must be finite and nonnegative")
    effective_tie_rtol = float(
        32 * np.finfo(float).eps if tie_rtol is None else tie_rtol
    )

    old_dps = mp.iv.dps
    mp.iv.dps = math.ceil(precision_bits * math.log10(2)) + 5
    try:
        initial_boxes = (
            ((t0, t1, w0, 0.0), (t0, t1, 0.0, w1))
            if w0 < 0 < w1
            else ((t0, t1, w0, w1),)
        )
        queue: list[_IntervalCell] = []
        evaluated = 0
        serial = 0
        best_value = -1.0
        best_point = (t0, w0)
        best_signed_value = mp.mpf(0)
        candidates: dict[tuple[float, float], float] = {}

        def observe(point: tuple[float, float], value: mp.mpf, lower: float) -> None:
            nonlocal best_value, best_point, best_signed_value
            candidates[point] = lower
            if lower > best_value or (lower == best_value and point < best_point):
                best_value, best_point, best_signed_value = lower, point, value

        def global_bound() -> float:
            upper = -queue[0].priority
            if not (0 <= best_value <= upper and np.isfinite(upper)):
                raise ValueError("inconsistent singleton and global interval bounds")
            return upper

        def result(certified: bool, reason: str) -> PivotCertificate:
            upper = global_bound()
            threshold = effective_tie_rtol * max(1.0, best_value)
            return PivotCertificate(
                t=best_point[0],
                omega=best_point[1],
                value=float(best_signed_value),
                lower_bound=best_value,
                upper_bound=upper,
                eta=1.0 if upper == 0 else _outward_lower(best_value / upper),
                cells_evaluated=evaluated,
                certified=certified,
                termination_reason=reason,
                candidate_count=len(candidates),
                exact_tie_count=int(
                    sum(value == best_value for value in candidates.values())
                ),
                near_tie_count=int(
                    sum(
                        best_value - value <= threshold for value in candidates.values()
                    )
                ),
                tie_rtol=effective_tie_rtol,
                tie_scope="distinct_sampled_singleton_lower_bounds_only",
                objective_scope=objective_scope,
            )

        for coordinates in initial_boxes:
            cell = _make_interval_cell(
                point_function,
                interval_function,
                coordinates,
                serial,
                precision_bits,
                gradient_interval_function,
            )
            serial += 1
            evaluated += 1
            heapq.heappush(queue, cell)
            point = (
                _midpoint(cell.t0, cell.t1),
                _midpoint(cell.w0, cell.w1),
            )
            observe(point, cell.point_value, cell.point_lower)

        seed_t = (t0, _midpoint(t0, t1), t1)
        seed_omega = (w0, _midpoint(w0, w1), w1)
        for t in seed_t:
            for omega in seed_omega:
                point = (t, omega)
                if point not in candidates:
                    signed_value, lower, _ = _interval_point(
                        point_function, interval_function, t, omega, precision_bits
                    )
                    observe(point, signed_value, lower)

        t_scale = t1 - t0
        omega_scale = w1 - w0
        while queue:
            if evaluated > max_cells:
                return result(False, "cell_budget")
            upper = global_bound()
            gap = _outward_upper(upper - best_value)
            tolerance = _outward_lower(abs_tol + _outward_lower(rel_tol * best_value))
            if gap <= tolerance:
                return result(True, "interval_gap_tolerance")
            if evaluated >= max_cells:
                break

            if evaluated + 2 > max_cells:
                break

            cell = heapq.heappop(queue)
            t_mid = _midpoint(cell.t0, cell.t1)
            w_mid = _midpoint(cell.w0, cell.w1)
            can_split_t = cell.t0 < t_mid < cell.t1
            can_split_w = cell.w0 < w_mid < cell.w1
            if not can_split_t and not can_split_w:
                heapq.heappush(queue, cell)
                return result(False, "coordinate_resolution")
            split_t = can_split_t and (
                not can_split_w
                or ((cell.t1 - cell.t0) / t_scale)
                >= ((cell.w1 - cell.w0) / omega_scale)
            )
            children = (
                (
                    (cell.t0, t_mid, cell.w0, cell.w1),
                    (t_mid, cell.t1, cell.w0, cell.w1),
                )
                if split_t
                else (
                    (cell.t0, cell.t1, cell.w0, w_mid),
                    (cell.t0, cell.t1, w_mid, cell.w1),
                )
            )
            for coordinates in children:
                child = _make_interval_cell(
                    point_function,
                    interval_function,
                    coordinates,
                    serial,
                    precision_bits,
                    gradient_interval_function,
                )
                serial += 1
                evaluated += 1
                heapq.heappush(queue, child)
                point = (
                    _midpoint(child.t0, child.t1),
                    _midpoint(child.w0, child.w1),
                )
                observe(point, child.point_value, child.point_lower)

        return result(False, "cell_budget")
    finally:
        mp.iv.dps = old_dps
