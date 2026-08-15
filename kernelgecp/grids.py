"""Deterministic composite Chebyshev grids for GECP experiments."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt


def chebyshev_lobatto(a: float, b: float, order: int) -> npt.NDArray[np.float64]:
    """Return ascending Chebyshev--Lobatto nodes on ``[a, b]``."""

    if not np.isfinite(a) or not np.isfinite(b) or a >= b:
        raise ValueError("expected finite endpoints with a < b")
    if order < 2:
        raise ValueError("order must be at least two")
    angles = np.linspace(0.0, np.pi, order)
    nodes = (a + b) / 2.0 + (b - a) * np.cos(angles) / 2.0
    return np.asarray(nodes[::-1], dtype=np.float64)


def _merge_nodes(parts: list[npt.NDArray[np.float64]]) -> npt.NDArray[np.float64]:
    rounded = np.concatenate(parts)
    return np.asarray(np.unique(rounded), dtype=np.float64)


def time_grid(order: int = 24) -> npt.NDArray[np.float64]:
    """Composite grid refined toward both time endpoints."""

    return _merge_nodes(
        [chebyshev_lobatto(0.0, 0.5, order), chebyshev_lobatto(0.5, 1.0, order)]
    )


def frequency_grid(cutoff: float, order: int = 24) -> npt.NDArray[np.float64]:
    """Dyadically partitioned frequency grid refined near zero and band edges."""

    if not np.isfinite(cutoff) or cutoff < 1:
        raise ValueError("cutoff must be finite and at least one")
    positive: list[npt.NDArray[np.float64]] = [chebyshev_lobatto(0.0, 1.0, order)]
    left = 1.0
    bands = max(0, math.ceil(math.log2(cutoff)))
    for _ in range(bands):
        right = min(cutoff, 2.0 * left)
        if right > left:
            positive.append(chebyshev_lobatto(left, right, order))
        left = right
        if left >= cutoff:
            break
    positive_nodes = _merge_nodes(positive)
    return np.concatenate((-positive_nodes[:0:-1], positive_nodes))


def composite_chebyshev_grids(
    cutoff: float, order: int = 24
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return the canonical time and frequency grids."""

    return time_grid(order), frequency_grid(cutoff, order)
