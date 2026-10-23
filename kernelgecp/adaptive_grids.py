"""Physical panels and deterministic float64 lattices for sampled refinement.

Coordinates are stored binary64 samples, even when kernel evaluation later uses
higher precision.  These helpers reject coordinate collapse rather than suggest
that an unrepresentable lattice supplies additional independent samples.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import MAX_EMAX, MIN_EMIN, Decimal, localcontext
from numbers import Real
from typing import Literal

import numpy as np

from kernelgecp.types import FloatArray

_TRAINING_FRACTIONS = (0.0, 0.25, 0.5, 0.75, 1.0)
_PROBE_FRACTIONS = (0.125, 0.375, 0.625, 0.875)
_VALIDATION_FRACTIONS = tuple(j / 7.0 for j in range(1, 7))


def _finite_float(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real number") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


@dataclass(frozen=True)
class Panel:
    """One physical leaf interval with a stable refinement-tree identifier."""

    identifier: str
    left: float
    right: float
    depth: int = 0
    parent: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.identifier, str) or not self.identifier:
            raise ValueError("panel identifier must be a nonempty string")
        left = _finite_float(self.left, "panel left endpoint")
        right = _finite_float(self.right, "panel right endpoint")
        if left >= right:
            raise ValueError("panel endpoints must satisfy left < right")
        if isinstance(self.depth, bool) or not isinstance(self.depth, int):
            raise ValueError("panel depth must be a nonnegative integer")
        if self.depth < 0:
            raise ValueError("panel depth must be a nonnegative integer")
        if self.parent is not None and (
            not isinstance(self.parent, str)
            or not self.parent
            or self.parent == self.identifier
        ):
            raise ValueError("panel parent must be a distinct nonempty identifier")
        object.__setattr__(self, "left", left)
        object.__setattr__(self, "right", right)


def _validated_panels(panels: Sequence[Panel]) -> tuple[Panel, ...]:
    result = tuple(panels)
    if not result or any(not isinstance(panel, Panel) for panel in result):
        raise ValueError("expected a nonempty sequence of panels")
    if len({panel.identifier for panel in result}) != len(result):
        raise ValueError("panel identifiers must be unique")
    if any(
        left.right != right.left
        for left, right in zip(result, result[1:], strict=False)
    ):
        raise ValueError("panels must form an ascending contiguous partition")
    return result


def _from_boundaries(boundaries: list[float], prefix: str) -> list[Panel]:
    return [
        Panel(f"{prefix}{index:04d}", left, right)
        for index, (left, right) in enumerate(
            zip(boundaries, boundaries[1:], strict=False)
        )
    ]


def initial_panels(cutoff: float) -> tuple[list[Panel], list[Panel]]:
    """Return endpoint-resolving time and dyadic-frequency physical panels.

    The supplied cutoff is retained exactly as the outer frequency endpoint.
    A cutoff whose time endpoint layer collapses at binary64 precision fails
    explicitly, instead of emitting duplicate panels or out-of-domain nodes.
    """

    cutoff = _finite_float(cutoff, "cutoff")
    if cutoff < 1.0:
        raise ValueError("cutoff must be finite and at least one")
    left_time = [0.0, min(0.5, 1.0 / cutoff)]
    while left_time[-1] < 0.5:
        left_time.append(min(0.5, 2.0 * left_time[-1]))
    time_boundaries = left_time + [1.0 - value for value in left_time[-2::-1]]
    if any(a >= b for a, b in zip(time_boundaries, time_boundaries[1:], strict=False)):
        raise ValueError("time endpoint panels collapse at binary64 precision")

    positive_frequency = [0.0, 1.0]
    while positive_frequency[-1] < cutoff:
        left = positive_frequency[-1]
        # Avoid overflowing an intermediate doubling for very large cutoffs.
        positive_frequency.append(cutoff if left >= cutoff / 2.0 else 2.0 * left)
    frequency_boundaries = [
        -value for value in positive_frequency[:0:-1]
    ] + positive_frequency
    time_panels = _from_boundaries(time_boundaries, "t")
    frequency_panels = _from_boundaries(frequency_boundaries, "omega")
    # An initial panel must already support the promised five-node lattice.
    training_nodes(time_panels)
    training_nodes(frequency_panels)
    return time_panels, frequency_panels


def _panel_nodes(panel: Panel, fractions: tuple[float, ...]) -> FloatArray:
    nodes = np.asarray(
        [
            panel.left
            if fraction == 0.0
            else panel.right
            if fraction == 1.0
            else (1.0 - fraction) * panel.left + fraction * panel.right
            for fraction in fractions
        ],
        dtype=np.float64,
    )
    if (
        not np.all(np.isfinite(nodes))
        or not np.all(np.diff(nodes) > 0.0)
        or nodes[0] < panel.left
        or nodes[-1] > panel.right
        or (fractions[0] > 0.0 and nodes[0] == panel.left)
        or (fractions[-1] < 1.0 and nodes[-1] == panel.right)
    ):
        raise ValueError(
            f"panel {panel.identifier!r} lattice collapses at binary64 precision"
        )
    return nodes


def _sample_nodes(panels: Sequence[Panel], fractions: tuple[float, ...]) -> FloatArray:
    leaves = _validated_panels(panels)
    return np.asarray(
        np.unique(np.concatenate([_panel_nodes(panel, fractions) for panel in leaves])),
        dtype=np.float64,
    )


def training_nodes(
    panels: Sequence[Panel], previous: FloatArray | None = None
) -> FloatArray:
    """Return five quarter-lattice nodes per panel, retaining old floats exactly.

    ``previous`` must be a finite ascending axis on the same physical domain.
    Unioning stored coordinates avoids loss from a different rounding path
    after bisection; it does not promote their coordinate precision beyond 53
    binary significant bits.
    """

    leaves = _validated_panels(panels)
    current = _sample_nodes(leaves, _TRAINING_FRACTIONS)
    if previous is None:
        return current
    old = np.asarray(previous, dtype=np.float64)
    if (
        old.ndim != 1
        or old.size < 2
        or not np.all(np.isfinite(old))
        or not np.all(np.diff(old) > 0.0)
        or old[0] != leaves[0].left
        or old[-1] != leaves[-1].right
    ):
        raise ValueError(
            "previous nodes must be a finite ascending axis on this domain"
        )
    return np.asarray(np.union1d(current, old), dtype=np.float64)


def probe_nodes(panels: Sequence[Panel]) -> FloatArray:
    """Return the four odd-eighth samples inside each leaf panel."""

    leaves = _validated_panels(panels)
    probes = _sample_nodes(leaves, _PROBE_FRACTIONS)
    if np.intersect1d(probes, _sample_nodes(leaves, _TRAINING_FRACTIONS)).size:
        raise ValueError("probe and training coordinates collide at binary64 precision")
    return probes


def validation_nodes(panels: Sequence[Panel]) -> FloatArray:
    """Return six seventh-fraction samples independent of the leaf lattices.

    These coordinates are reserved for final checking.  The function rejects
    any binary64 collision with the training or marking/stopping probe lattice.
    """

    leaves = _validated_panels(panels)
    validation = _sample_nodes(leaves, _VALIDATION_FRACTIONS)
    used = np.union1d(
        _sample_nodes(leaves, _TRAINING_FRACTIONS),
        _sample_nodes(leaves, _PROBE_FRACTIONS),
    )
    if np.intersect1d(validation, used).size:
        raise ValueError("validation coordinates collide at binary64 precision")
    return validation


def _reflection_mates(
    panels: tuple[Panel, ...], axis: Literal["t", "omega"]
) -> dict[str, str]:
    if axis == "t":
        if panels[0].left != 0.0 or panels[-1].right != 1.0:
            raise ValueError("time panels must partition [0, 1]")
        shift, scale = 1.0, 1.0
    elif axis == "omega":
        if panels[0].left != -panels[-1].right:
            raise ValueError("frequency panels must have symmetric physical endpoints")
        shift, scale = 0.0, panels[-1].right
    else:
        raise ValueError("axis must be 't' or 'omega'")
    result: dict[str, str] = {}
    for panel in panels:
        reflected_left = shift - panel.right
        reflected_right = shift - panel.left
        tolerance = 4.0 * max(
            math.ulp(scale),
            math.ulp(reflected_left),
            math.ulp(reflected_right),
        )
        matches = [
            candidate.identifier
            for candidate in panels
            if abs(candidate.left - reflected_left) <= tolerance
            and abs(candidate.right - reflected_right) <= tolerance
        ]
        if len(matches) != 1:
            raise ValueError(
                "panels lack unambiguous reflection mates at binary64 precision"
            )
        result[panel.identifier] = matches[0]
    if any(result[result[identifier]] != identifier for identifier in result):
        raise ValueError("panel reflection must be an involution")
    return result


def mark_panels(
    panels: Sequence[Panel],
    scores: Mapping[str, Decimal],
    axis: Literal["t", "omega"],
    max_marked: int = 4,
    relative_threshold: Decimal = Decimal("0.5"),
) -> list[str]:
    """Mark influential positive-score panels and their reflection mates.

    Select scores at least ``relative_threshold`` times the positive maximum,
    sort by descending score then identifier, and apply ``max_marked`` before
    adding mirror mates.  Missing scores are zero; invalid or unknown scores
    fail explicitly.  Returned identifiers follow physical panel order.
    """

    leaves = _validated_panels(panels)
    if (
        isinstance(max_marked, bool)
        or not isinstance(max_marked, int)
        or max_marked < 1
    ):
        raise ValueError("max_marked must be a positive integer")
    if (
        not isinstance(relative_threshold, Decimal)
        or not relative_threshold.is_finite()
        or not Decimal(0) < relative_threshold <= Decimal(1)
    ):
        raise ValueError("relative_threshold must be a finite Decimal in (0, 1]")
    identifiers = {panel.identifier for panel in leaves}
    for identifier, score in scores.items():
        if identifier not in identifiers:
            raise ValueError(f"unknown panel score identifier {identifier!r}")
        if not isinstance(score, Decimal) or not score.is_finite() or score < 0:
            raise ValueError("panel scores must be finite nonnegative Decimals")
    mates = _reflection_mates(leaves, axis)
    positive = [
        (identifier, score) for identifier, score in scores.items() if score > 0
    ]
    if not positive:
        return []
    maximum = max(score for _, score in positive)
    # Decimal multiplication and unary minus otherwise inherit the caller's
    # precision.  In particular, unary minus can erase a small score difference
    # and accidentally turn it into an identifier tie.
    with localcontext() as context:
        context.prec = len(maximum.as_tuple().digits) + len(
            relative_threshold.as_tuple().digits
        )
        context.Emax = MAX_EMAX
        context.Emin = MIN_EMIN
        threshold = relative_threshold * maximum
    candidates = sorted(
        (item for item in positive if item[1] >= threshold), key=lambda item: item[0]
    )
    candidates.sort(key=lambda item: item[1], reverse=True)
    selected = {identifier for identifier, _ in candidates[:max_marked]}
    selected.update(mates[identifier] for identifier in tuple(selected))
    return [panel.identifier for panel in leaves if panel.identifier in selected]


def refine_panels(panels: Sequence[Panel], marked_ids: Sequence[str]) -> list[Panel]:
    """Bisect exactly the marked panels and preserve all other leaf objects."""

    leaves = _validated_panels(panels)
    marked = set(marked_ids)
    if len(marked) != len(marked_ids):
        raise ValueError("marked panel identifiers must be unique")
    if not marked <= {panel.identifier for panel in leaves}:
        raise ValueError("marked identifiers must name current leaf panels")
    refined: list[Panel] = []
    for panel in leaves:
        if panel.identifier not in marked:
            refined.append(panel)
            continue
        midpoint = 0.5 * panel.left + 0.5 * panel.right
        if not panel.left < midpoint < panel.right:
            raise ValueError("panel bisection collapses at binary64 precision")
        children = [
            Panel(
                panel.identifier + "L",
                panel.left,
                midpoint,
                panel.depth + 1,
                panel.identifier,
            ),
            Panel(
                panel.identifier + "R",
                midpoint,
                panel.right,
                panel.depth + 1,
                panel.identifier,
            ),
        ]
        for child in children:
            _panel_nodes(child, _TRAINING_FRACTIONS)
        refined.extend(children)
    return list(_validated_panels(refined))
