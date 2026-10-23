"""Nested sampling, physical geometry, and deterministic refinement checks."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from decimal import Decimal, localcontext

import numpy as np
import pytest

from kernelgecp.adaptive_grids import (
    Panel,
    initial_panels,
    mark_panels,
    probe_nodes,
    refine_panels,
    training_nodes,
    validation_nodes,
)


@pytest.mark.parametrize("cutoff", [1, 1.2, 1.7, 10, 100, 1000, 10000, 100000, 1e6])
def test_initial_physical_panels_and_independent_lattices(cutoff: float) -> None:
    time, frequency = initial_panels(cutoff)
    for panels, left, right, reflected in (
        (time, 0.0, 1.0, 1.0),
        (frequency, -cutoff, cutoff, 0.0),
    ):
        train = training_nodes(panels)
        probes = probe_nodes(panels)
        validation = validation_nodes(panels)
        assert train[0] == left
        assert train[-1] == right
        assert len(train) == 4 * len(panels) + 1
        assert len(probes) == 4 * len(panels)
        assert len(validation) == 6 * len(panels)
        assert np.all(np.diff(train) > 0)
        assert np.all(probes > left) and np.all(probes < right)
        assert np.all(validation > left) and np.all(validation < right)
        assert not np.intersect1d(train, probes).size
        assert not np.intersect1d(validation, np.union1d(train, probes)).size
        np.testing.assert_allclose(train, reflected - train[::-1], rtol=0, atol=2e-16)
        assert all(a.right == b.left for a, b in zip(panels, panels[1:], strict=False))


def test_initial_panels_retain_existing_endpoint_and_frequency_breakpoints() -> None:
    time, frequency = initial_panels(10)
    np.testing.assert_array_equal(
        [panel.left for panel in time] + [time[-1].right],
        [0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1],
    )
    np.testing.assert_array_equal(
        [panel.left for panel in frequency] + [frequency[-1].right],
        [-10, -8, -4, -2, -1, 0, 1, 2, 4, 8, 10],
    )
    assert len(initial_panels(1)[0]) == 2
    assert len(initial_panels(1)[1]) == 2


def test_selective_refinement_retains_every_previous_coordinate_exactly() -> None:
    panels, _ = initial_panels(10)
    previous = training_nodes(panels)
    old_identifiers = [panel.identifier for panel in panels]
    for _ in range(3):
        scores = {panels[0].identifier: Decimal(10)}
        marked = mark_panels(panels, scores, "t", max_marked=1)
        assert marked == [panels[0].identifier, panels[-1].identifier]
        refined = refine_panels(panels, marked)
        assert len(refined) == len(panels) + 2
        assert refined[2:-2] == panels[1:-1]
        current = training_nodes(refined, previous)
        assert np.array_equal(np.intersect1d(previous, current), previous)
        assert len(current) > len(previous)
        assert refined[0].identifier == panels[0].identifier + "L"
        assert refined[0].parent == panels[0].identifier
        assert refined[0].depth == panels[0].depth + 1
        panels, previous = refined, current
    assert old_identifiers[1:-1] == [panel.identifier for panel in panels[4:-4]]


def test_previous_axis_union_preserves_rounding_path_and_deterministic_bytes() -> None:
    panels, _ = initial_panels(1.7)
    previous = training_nodes(panels)
    previous[1] = np.nextafter(previous[1], previous[2])
    refined = refine_panels(panels, [panels[0].identifier, panels[-1].identifier])
    current = training_nodes(refined, previous)
    assert previous[1] in current
    assert current.tobytes() == training_nodes(refined, previous).tobytes()
    assert validation_nodes(refined).tobytes() == validation_nodes(refined).tobytes()


@pytest.mark.parametrize("axis", ["t", "omega"])
def test_score_threshold_cap_reflection_and_tie_breaking(axis: str) -> None:
    time, frequency = initial_panels(10)
    panels = time if axis == "t" else frequency
    # Reverse insertion order deliberately; identifiers, rather than mapping
    # iteration order, decide which equally influential panel gets the cap.
    scores = {panel.identifier: Decimal(10) for panel in reversed(panels)}
    marked = mark_panels(panels, scores, axis, max_marked=1)  # type: ignore[arg-type]
    assert marked == [panels[0].identifier, panels[-1].identifier]
    scores = {
        panels[0].identifier: Decimal(10),
        panels[1].identifier: Decimal(5),
        panels[2].identifier: Decimal("4.999"),
    }
    marked = mark_panels(panels, scores, axis, max_marked=2)  # type: ignore[arg-type]
    assert marked == [
        panels[0].identifier,
        panels[1].identifier,
        panels[-2].identifier,
        panels[-1].identifier,
    ]
    assert mark_panels(panels, {}, axis) == []  # type: ignore[arg-type]
    assert mark_panels(panels, {panels[0].identifier: Decimal(0)}, axis) == []  # type: ignore[arg-type]


def test_center_panel_can_be_its_own_reflection() -> None:
    panels = [Panel("a", -1, -0.25), Panel("b", -0.25, 0.25), Panel("c", 0.25, 1)]
    assert mark_panels(panels, {"b": Decimal(1)}, "omega") == ["b"]


def test_score_decisions_do_not_inherit_decimal_rounding_precision() -> None:
    panels, _ = initial_panels(10)
    scores = {
        panels[0].identifier: Decimal("10"),
        panels[2].identifier: Decimal("10.000000000000000000000000000001"),
    }
    with localcontext() as context:
        context.prec = 3
        assert mark_panels(panels, scores, "t", max_marked=1) == [
            panels[2].identifier,
            panels[-3].identifier,
        ]
        scores[panels[2].identifier] = Decimal("20.000000000000000000000000000002")
        assert mark_panels(panels, scores, "t", max_marked=4) == [
            panels[2].identifier,
            panels[-3].identifier,
        ]


@pytest.mark.parametrize("cutoff", [0, 0.9, -1, float("nan"), float("inf"), True])
def test_invalid_cutoffs_fail(cutoff: float) -> None:
    with pytest.raises(ValueError):
        initial_panels(cutoff)


def test_unrepresentable_endpoint_layer_fails() -> None:
    with pytest.raises(ValueError, match="binary64"):
        initial_panels(1e20)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"identifier": ""},
        {"left": float("nan")},
        {"right": float("inf")},
        {"right": 0},
        {"depth": -1},
        {"depth": True},
        {"parent": "a"},
    ],
)
def test_malformed_panels_fail(kwargs: dict[str, object]) -> None:
    arguments: dict[str, object] = {"identifier": "a", "left": 0.0, "right": 1.0}
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        Panel(**arguments)  # type: ignore[arg-type]


def test_panels_are_immutable() -> None:
    panel = Panel("a", 0, 1)
    with pytest.raises(FrozenInstanceError):
        panel.left = 0.25  # type: ignore[misc]


@pytest.mark.parametrize(
    "panels",
    [
        [],
        [Panel("a", 0, 0.5), Panel("a", 0.5, 1)],
        [Panel("a", 0, 0.4), Panel("b", 0.5, 1)],
        [Panel("a", 0.5, 1), Panel("b", 0, 0.5)],
    ],
)
def test_malformed_partitions_fail(panels: list[Panel]) -> None:
    for sample in (training_nodes, probe_nodes, validation_nodes):
        with pytest.raises(ValueError):
            sample(panels)
    with pytest.raises(ValueError):
        refine_panels(panels, [])


@pytest.mark.parametrize(
    "previous",
    [
        np.array([]),
        np.array([[0.0, 1.0]]),
        np.array([0.0, 0.5, 0.5, 1.0]),
        np.array([0.0, 0.7, 0.5, 1.0]),
        np.array([0.0, float("nan"), 1.0]),
        np.array([-0.1, 0.5, 1.0]),
        np.array([0.0, 0.5, 1.1]),
    ],
)
def test_invalid_previous_coordinates_fail(previous: np.ndarray) -> None:  # type: ignore[type-arg]
    panels, _ = initial_panels(1)
    with pytest.raises(ValueError, match="previous nodes"):
        training_nodes(panels, previous)


def test_coordinate_budget_exhaustion_fails_instead_of_reusing_samples() -> None:
    tiny = [Panel("tiny", 1.0, np.nextafter(1.0, 2.0))]
    for sample in (training_nodes, probe_nodes, validation_nodes):
        with pytest.raises(ValueError, match="binary64"):
            sample(tiny)
    with pytest.raises(ValueError, match="binary64"):
        refine_panels(tiny, ["tiny"])


@pytest.mark.parametrize("score", [Decimal(-1), Decimal("NaN"), Decimal("Infinity")])
def test_invalid_scores_fail(score: Decimal) -> None:
    panels, _ = initial_panels(1)
    with pytest.raises(ValueError, match="scores"):
        mark_panels(panels, {panels[0].identifier: score}, "t")


@pytest.mark.parametrize("budget", [0, -1, True, 1.5])
def test_invalid_marking_budget_fails(budget: int) -> None:
    panels, _ = initial_panels(1)
    with pytest.raises(ValueError, match="max_marked"):
        mark_panels(panels, {}, "t", max_marked=budget)


@pytest.mark.parametrize("threshold", [Decimal(0), Decimal("1.1"), Decimal("NaN")])
def test_invalid_marking_threshold_fails(threshold: Decimal) -> None:
    panels, _ = initial_panels(1)
    with pytest.raises(ValueError, match="relative_threshold"):
        mark_panels(panels, {}, "t", relative_threshold=threshold)


def test_unknown_marking_and_refinement_identifiers_fail() -> None:
    panels, _ = initial_panels(1)
    with pytest.raises(ValueError, match="unknown"):
        mark_panels(panels, {"missing": Decimal(1)}, "t")
    with pytest.raises(ValueError, match="current leaf"):
        refine_panels(panels, ["missing"])
    with pytest.raises(ValueError, match="unique"):
        refine_panels(panels, [panels[0].identifier] * 2)


def test_asymmetric_partition_cannot_claim_reflection() -> None:
    panels = [Panel("a", 0, 0.25), Panel("b", 0.25, 1)]
    with pytest.raises(ValueError, match="reflection"):
        mark_panels(panels, {"a": Decimal(1)}, "t")
    with pytest.raises(ValueError, match="partition"):
        mark_panels([Panel("a", 0, 2)], {"a": Decimal(1)}, "t")
