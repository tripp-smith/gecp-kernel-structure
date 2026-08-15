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
    return float(np.nextafter(value, math.inf))


def _outward_lower(value: float) -> float:
    return max(0.0, float(np.nextafter(value, -math.inf)))


def _make_interval_cell(
    point_function: MPPointFunction,
    interval_function: MPIntervalFunction,
    coordinates: tuple[float, float, float, float],
    serial: int,
    precision_bits: int,
    gradient_interval_function: MPGradientIntervalFunction | None,
) -> _IntervalCell:
    t0, t1, w0, w1 = coordinates
    t = (t0 + t1) / 2.0
    omega = (w0 + w1) / 2.0
    with mp.workprec(precision_bits):
        point_value = point_function(mp.mpf(repr(t)), mp.mpf(repr(omega)))
    enclosure = interval_function(mp.iv.mpf([t0, t1]), mp.iv.mpf([w0, w1]))
    lower_endpoint = float(enclosure.a)
    upper_endpoint = float(enclosure.b)
    upper = _outward_upper(max(abs(lower_endpoint), abs(upper_endpoint)))
    if gradient_interval_function is not None:
        center_enclosure = interval_function(mp.iv.mpf(t), mp.iv.mpf(omega))
        center_upper = _outward_upper(
            max(abs(float(center_enclosure.a)), abs(float(center_enclosure.b)))
        )
        dt_enclosure, dw_enclosure = gradient_interval_function(
            mp.iv.mpf([t0, t1]), mp.iv.mpf([w0, w1])
        )
        dt_upper = _outward_upper(
            max(abs(float(dt_enclosure.a)), abs(float(dt_enclosure.b)))
        )
        dw_upper = _outward_upper(
            max(abs(float(dw_enclosure.a)), abs(float(dw_enclosure.b)))
        )
        dt_term = _outward_upper(dt_upper * ((t1 - t0) / 2.0))
        dw_term = _outward_upper(dw_upper * ((w1 - w0) / 2.0))
        mean_upper = _outward_upper(
            _outward_upper(center_upper + dt_term) + dw_term
        )
        upper = min(upper, mean_upper)
    if not mp.isfinite(point_value) or not np.isfinite(upper):
        raise ValueError("interval pivot objective produced a non-finite enclosure")
    return _IntervalCell(
        -upper,
        serial,
        t0,
        t1,
        w0,
        w1,
        point_value,
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
) -> PivotCertificate:
    """Certify an absolute maximizer using outward-rounded interval enclosures."""

    t0, t1 = t_bounds
    w0, w1 = omega_bounds
    if not all(np.isfinite((t0, t1, w0, w1))) or t0 >= t1 or w0 >= w1:
        raise ValueError("bounds must be finite, ordered intervals")
    if precision_bits < 53 or abs_tol < 0 or rel_tol < 0 or max_cells < 1:
        raise ValueError("invalid interval certificate precision, tolerance, or budget")

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
                (cell.t0 + cell.t1) / 2.0,
                (cell.w0 + cell.w1) / 2.0,
            )
            magnitude = _outward_lower(abs(float(cell.point_value)))
            if magnitude > best_value or (
                magnitude == best_value and point < best_point
            ):
                best_value = magnitude
                best_point = point
                best_signed_value = cell.point_value

        seed_t = (t0, (t0 + t1) / 2.0, t1)
        seed_omega = (w0, (w0 + w1) / 2.0, w1)
        for t in seed_t:
            for omega in seed_omega:
                with mp.workprec(precision_bits):
                    signed_value = point_function(
                        mp.mpf(repr(t)), mp.mpf(repr(omega))
                    )
                magnitude = _outward_lower(abs(float(signed_value)))
                point = (t, omega)
                if magnitude > best_value or (
                    magnitude == best_value and point < best_point
                ):
                    best_value = magnitude
                    best_point = point
                    best_signed_value = signed_value

        t_scale = t1 - t0
        omega_scale = w1 - w0
        while queue:
            global_upper = -queue[0].priority
            tolerance = abs_tol + rel_tol * max(0.0, best_value)
            if global_upper - best_value <= tolerance:
                eta = (
                    1.0
                    if global_upper == 0
                    else min(1.0, best_value / global_upper)
                )
                return PivotCertificate(
                    t=best_point[0],
                    omega=best_point[1],
                    value=float(best_signed_value),
                    lower_bound=best_value,
                    upper_bound=global_upper,
                    eta=eta,
                    cells_evaluated=evaluated,
                    certified=True,
                    termination_reason="interval_gap_tolerance",
                )
            if evaluated >= max_cells:
                break

            if evaluated + 2 > max_cells:
                break

            cell = heapq.heappop(queue)
            t_mid = (cell.t0 + cell.t1) / 2.0
            w_mid = (cell.w0 + cell.w1) / 2.0
            split_t = ((cell.t1 - cell.t0) / t_scale) >= (
                (cell.w1 - cell.w0) / omega_scale
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
                    (child.t0 + child.t1) / 2.0,
                    (child.w0 + child.w1) / 2.0,
                )
                magnitude = _outward_lower(abs(float(child.point_value)))
                if magnitude > best_value or (
                    magnitude == best_value and point < best_point
                ):
                    best_value = magnitude
                    best_point = point
                    best_signed_value = child.point_value

        global_upper = max((-cell.priority for cell in queue), default=best_value)
        eta = (
            1.0 if global_upper == 0 else min(1.0, best_value / global_upper)
        )
        return PivotCertificate(
            t=best_point[0],
            omega=best_point[1],
            value=float(best_signed_value),
            lower_bound=best_value,
            upper_bound=global_upper,
            eta=eta,
            cells_evaluated=evaluated,
            certified=False,
            termination_reason="cell_budget",
        )
    finally:
        mp.iv.dps = old_dps
