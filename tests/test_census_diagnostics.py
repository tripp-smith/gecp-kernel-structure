"""Independent exact checks for schema-2 census geometry and prefixes."""

from __future__ import annotations

from dataclasses import asdict, replace
from decimal import Decimal
from fractions import Fraction
from typing import Any

import numpy as np
import pytest

from kernelgecp.census import _pivot_geometry, _records_from_result
from kernelgecp.types import GECPConfig, GECPResult, HighPrecisionData


def _fraction(value: str) -> Fraction:
    return Fraction(Decimal(value))


def test_pivot_geometry_matches_independent_exact_boundary_and_reflection() -> None:
    times = ["0", "1", "0.25", "0.75", "0.5"]
    frequencies = ["8", "-8", "2", "-3", "0"]
    geometry = _pivot_geometry(times, frequencies, "8", 128, "0.001")
    expected_partners = [1, 0, 3, 2, 4]
    expected_distances = [Fraction(0), Fraction(0), Fraction(1, 8), Fraction(1, 8), 0]
    expected_scales = [3, 3, 1, 2, 0]
    expected_panels = ["positive", "negative", "positive", "negative", "central"]
    for index, item in enumerate(geometry):
        time = Fraction(times[index])
        frequency = Fraction(frequencies[index])
        assert tuple(map(_fraction, item.t_boundary_distances)) == (time, 1 - time)
        assert tuple(map(_fraction, item.omega_boundary_distances)) == (
            frequency + 8,
            8 - frequency,
        )
        assert _fraction(item.reflection_t) == 1 - time
        assert _fraction(item.reflection_omega) == -frequency
        # Compute the minimizer independently, including lexicographic ties.
        distances = [
            max(abs(Fraction(t) - (1 - time)), abs(Fraction(w) + frequency) / 8)
            for t, w in zip(times, frequencies, strict=True)
        ]
        assert item.reflection_nearest_index == distances.index(min(distances))
        assert item.reflection_nearest_index == expected_partners[index]
        assert _fraction(item.reflection_distance) == expected_distances[index]
        assert item.reflection_matched is (index in (0, 1, 4))
        assert item.dyadic_scale == expected_scales[index]
        assert item.dyadic_panel == expected_panels[index]

    inclusive = _pivot_geometry(times, frequencies, "8", 128, "0.125")
    assert all(item.reflection_matched for item in inclusive)


@pytest.mark.parametrize(
    ("frequency", "scale", "panel"),
    [
        ("-8", 3, "negative"),
        ("-4", 2, "negative"),
        ("-2", 1, "negative"),
        ("-1", 0, "central"),
        ("0", 0, "central"),
        ("1", 0, "central"),
        ("1.25", 1, "positive"),
        ("2", 1, "positive"),
        ("2.25", 2, "positive"),
        ("4", 2, "positive"),
        ("8", 3, "positive"),
    ],
)
def test_dyadic_scale_uses_least_closed_upper_power_of_two(
    frequency: str, scale: int, panel: str
) -> None:
    item = _pivot_geometry(["0.5"], [frequency], "8", 128, "0")[0]
    assert item.dyadic_scale == scale
    assert item.dyadic_panel == panel
    magnitude = abs(Fraction(frequency))
    assert magnitude <= 2**scale
    if scale:
        assert 2 ** (scale - 1) < magnitude


def test_reflection_nearest_ties_choose_first_prefix_index() -> None:
    geometry = _pivot_geometry(["0", "1", "1"], ["0", "-1", "1"], "2", 128, "0")
    assert geometry[0].reflection_nearest_index == 1
    assert _fraction(geometry[0].reflection_distance) == Fraction(1, 2)
    assert not geometry[0].reflection_matched
    # The later exact partner is unavailable to the rank-one prefix.
    full = _pivot_geometry(["0", "1"], ["2", "-2"], "2", 128, "0")
    prefix = _pivot_geometry(["0"], ["2"], "2", 128, "0")
    assert full[0].reflection_nearest_index == 1
    assert full[0].reflection_matched
    assert prefix[0].reflection_nearest_index == 0
    assert _fraction(prefix[0].reflection_distance) == 2
    assert not prefix[0].reflection_matched


def _raw(
    tolerances: list[float], *, max_rank: int, pivot: str = "grid"
) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "cutoffs": [8],
        "tolerances": tolerances,
        "precision_bits": 128,
        "max_rank": max_rank,
        "grid_order": 4,
        "pivot": pivot,
        "reflection_match_tolerance": "1e-12",
    }


def _manual_result(*, rank_zero: bool = False, adaptive: bool = False) -> GECPResult:
    precision = HighPrecisionData(
        t_nodes=[] if rank_zero else ["0", "0.75"],
        omega_nodes=[] if rank_zero else ["8", "-4"],
        pivots=[] if rank_zero else ["1", "-0.25"],
        residual_history=["0.03125"] if rank_zero else ["1", "0.25", "0.0625"],
        core_det_history=[] if rank_zero else ["1", "-0.25"],
        core_sigma_min_history=[] if rank_zero else ["1", "0.2"],
    )
    configuration = GECPConfig(
        tol=0.01,
        max_rank=0 if rank_zero else 2,
        grid_order=4,
        pivot="adaptive" if adaptive else "grid",
        precision_bits=128,
    ).normalized()
    return GECPResult(
        t_nodes=np.asarray([float(value) for value in precision.t_nodes]),
        omega_nodes=np.asarray([float(value) for value in precision.omega_nodes]),
        row_indices=[] if rank_zero else [0, 1],
        column_indices=[] if rank_zero else [0, 1],
        pivots=[float(value) for value in precision.pivots],
        residual_history=[float(value) for value in precision.residual_history],
        core_det_history=[float(value) for value in precision.core_det_history],
        core_sigma_min_history=[
            float(value) for value in precision.core_sigma_min_history
        ],
        tie_counts=[] if rank_zero else [2, 1],
        rank=0 if rank_zero else 2,
        converged=False,
        stop_reason="uncertified_pivot" if adaptive else "max_rank",
        config=configuration,
        precision_bits=128,
        residual_matrix=None if adaptive else np.asarray([[0.03125]]),
        high_precision=precision,
    )


def test_census_prefixes_keep_budget_and_complete_history_alignment() -> None:
    result = _manual_result()
    records = _records_from_result(
        _raw([0.3, 0.1, 0.01], max_rank=2),
        8.0,
        "revision",
        "configuration-hash",
        result,
        {"git_dirty": False, "environment": {}},
    )
    assert [record.rank for record in records] == [1, 2, 2]
    assert [record.converged for record in records] == [True, True, False]
    assert [record.stop_reason for record in records] == [
        "tolerance",
        "tolerance",
        "max_rank",
    ]
    for record, tolerance in zip(records, [0.3, 0.1, 0.01], strict=True):
        data = asdict(record)
        rank = record.rank
        assert data["schema_version"] == 2
        assert data["git_commit"] == "revision"
        assert data["config_hash"] == "configuration-hash"
        assert Decimal(data["tolerance"]) == Decimal(str(tolerance))
        assert data["residual_history"] == ["1", "0.25", "0.0625"][: rank + 1]
        assert data["residual"] == data["residual_history"][-1]
        assert data["pivots"] == ["1", "-0.25"][:rank]
        assert list(map(Decimal, data["pivot_magnitudes"])) == list(
            map(Decimal, ["1", "0.25"][:rank])
        )
        assert data["core_det_history"] == ["1", "-0.25"][:rank]
        assert data["core_sigma_min"] == ["1", "0.2"][:rank]
        assert len(data["pivot_geometry"]) == rank
        assert data["tie_counts"] == [2, 1][:rank]
        metadata = data["metadata"]
        assert metadata["converged"] is record.converged
        assert metadata["stop_reason"] == record.stop_reason
        assert metadata["configuration"]["tol"] == tolerance
        assert len(metadata["conditioning_history"]) == rank
        assert metadata["certification_status"] == "not_certified"
    assert result.config["tol"] == 0.01
    assert result.stop_reason == "max_rank"


def test_rank_zero_census_has_one_residual_and_no_pivot_diagnostics() -> None:
    records = _records_from_result(
        _raw([0.125, 0.01], max_rank=0),
        8.0,
        "revision",
        "configuration-hash",
        _manual_result(rank_zero=True),
        {"git_dirty": False, "environment": {}},
    )
    assert [record.rank for record in records] == [0, 0]
    assert [record.converged for record in records] == [True, False]
    for record in records:
        data = asdict(record)
        assert data["residual_history"] == ["0.03125"]
        for field in (
            "pivots",
            "pivot_magnitudes",
            "core_det_history",
            "core_sigma_min",
            "t_nodes",
            "omega_nodes",
            "tie_counts",
            "pivot_geometry",
        ):
            assert data[field] == []
        assert data["metadata"]["conditioning_history"] == []


def test_uncertified_adaptive_terminal_bound_does_not_imply_convergence() -> None:
    result = _manual_result(adaptive=True)
    record = _records_from_result(
        _raw([0.1], max_rank=2, pivot="adaptive"),
        8.0,
        "revision",
        "configuration-hash",
        result,
        {"git_dirty": False, "environment": {}},
    )[0]
    assert Decimal(record.residual) <= Decimal("0.1")
    assert record.rank == 2
    assert not record.converged
    assert record.stop_reason == "uncertified_pivot"
    metadata = asdict(record)["metadata"]
    assert metadata["converged"] is False
    assert metadata["certification_status"] == "not_certified"
    assert metadata["stop_reason"] == "uncertified_pivot"


@pytest.mark.parametrize(
    "field",
    [
        "pivots",
        "t_nodes",
        "omega_nodes",
        "residual_history",
        "core_det_history",
        "core_sigma_min_history",
    ],
)
def test_misaligned_census_history_fails_closed(field: str) -> None:
    result = _manual_result()
    assert result.high_precision is not None
    result.high_precision = replace(result.high_precision, **{field: []})
    with pytest.raises(ValueError, match="misaligned"):
        _records_from_result(
            _raw([0.1, 0.01], max_rank=2),
            8.0,
            "revision",
            "configuration-hash",
            result,
            {"git_dirty": False, "environment": {}},
        )


@pytest.mark.parametrize("terminal", ["NaN", "Infinity", "-Infinity", "-0.0625"])
def test_nonfinite_or_negative_census_residual_is_rejected(terminal: str) -> None:
    result = _manual_result()
    assert result.high_precision is not None
    result.high_precision = replace(
        result.high_precision, residual_history=["1", "0.25", terminal]
    )
    with pytest.raises(ValueError):
        _records_from_result(
            _raw([0.1, 0.01], max_rank=2),
            8.0,
            "revision",
            "configuration-hash",
            result,
            {"git_dirty": False, "environment": {}},
        )
