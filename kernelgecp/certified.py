"""Deterministic branch-and-bound certificates for approximate pivots."""

from __future__ import annotations

import heapq
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from .types import PivotCertificate

ScalarKernel = Callable[[float, float], float]
LipschitzBound = (
    float | tuple[float, float] | Callable[[float, float], tuple[float, float]]
)


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
            if evaluated >= max_cells:
                break
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
