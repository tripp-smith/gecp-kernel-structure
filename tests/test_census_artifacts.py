"""Audit committed census coverage and evidence scope without rerunning it."""

import hashlib
import json
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "experiments/data"
CUTOFFS = {Decimal(10) ** power for power in range(7)}
TOLERANCES = {Decimal(f"1e-{power}") for power in (6, 8, 10, 12)}


def _records(name: str) -> list[dict]:
    return [json.loads(line) for line in (DATA / name).read_text().splitlines()]


def test_committed_schema_two_census_covers_all_twenty_eight_cases() -> None:
    records = _records("gecp_census_v2.jsonl")
    assert len(records) == 28
    assert {(Decimal(r["cutoff"]), Decimal(r["tolerance"])) for r in records} == {
        (cutoff, tolerance) for cutoff in CUTOFFS for tolerance in TOLERANCES
    }
    assert {r["git_commit"] for r in records} == {
        "c0e8c5093134f369b93781c413750b7f1de46036"
    }
    for record in records:
        assert record["schema_version"] == 2
        assert record["converged"] and not record["git_dirty"]
        assert record["precision_bits"] == 128
        assert record["grid_coordinate_precision_bits"] == 53
        assert record["algorithm"] == "grid"
        rank = record["rank"]
        for field in (
            "pivots",
            "pivot_magnitudes",
            "t_nodes",
            "omega_nodes",
            "core_det_history",
            "core_sigma_min",
            "tie_counts",
            "pivot_geometry",
        ):
            assert len(record[field]) == rank
        assert len(record["residual_history"]) == rank + 1
        assert record["residual_history"][-1] == record["residual"]
        assert Decimal(record["residual"]) <= Decimal(record["tolerance"])
        assert all(Decimal(value).is_finite() for value in record["core_det_history"])
        assert all(Decimal(value) > 0 for value in record["core_sigma_min"])
        for pivot, magnitude in zip(
            record["pivots"], record["pivot_magnitudes"], strict=True
        ):
            assert Decimal(pivot).copy_abs() == Decimal(magnitude)
        for geometry in record["pivot_geometry"]:
            assert 0 <= geometry["reflection_nearest_index"] < rank
            assert Decimal(geometry["reflection_distance"]) >= 0
            assert all(Decimal(d) >= 0 for d in geometry["t_boundary_distances"])
            assert all(Decimal(d) >= 0 for d in geometry["omega_boundary_distances"])
        metadata = record["metadata"]
        assert metadata["converged"]
        assert metadata["convergence_scope"] == "sampled_grid"
        assert metadata["certification_status"] == "not_certified"
        assert len(metadata["conditioning_history"]) == rank


def test_resolution_ladder_contains_actual_denser_gecp_runs() -> None:
    canonical = {
        Decimal(r["cutoff"]): r
        for r in _records("gecp_census_v2.jsonl")
        if Decimal(r["tolerance"]) == Decimal("1e-12")
    }
    denser = _records("gecp_census_resolution.jsonl")
    assert len(denser) == 3
    assert {r["git_commit"] for r in denser} == {
        "ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61"
    }
    assert {Decimal(r["cutoff"]) for r in denser} == {
        Decimal(1),
        Decimal(1000),
        Decimal(1000000),
    }
    for record in denser:
        source = canonical[Decimal(record["cutoff"])]
        assert record["converged"] and not record["git_dirty"]
        assert record["metadata"]["configuration"]["grid_order"] == 32
        assert source["metadata"]["configuration"]["grid_order"] == 24
        assert all(
            a > b
            for a, b in zip(record["grid_shape"], source["grid_shape"], strict=True)
        )
        assert record["precision_bits"] == source["precision_bits"] == 128
        assert Decimal(record["residual"]) <= Decimal("1e-12")


def test_check_report_binds_actual_data_and_preserves_noncertification() -> None:
    report = json.loads((DATA / "gecp_census_checks.json").read_text())
    assert (
        report["input_sha256"]
        == hashlib.sha256((DATA / "gecp_census_v2.jsonl").read_bytes()).hexdigest()
    )
    assert report["precision_bits"] == 192
    assert report["source_record_count"] == 28
    assert report["strictest_trajectory_count"] == 7
    assert report["status"] == "passed_numerical_core_checks"
    assert report["certification_status"] == "not_certified"
    assert len(report["held_out_cross_checks"]) == 3
    ladder = report["resolution_ladder"]
    assert (
        ladder["input_sha256"]
        == hashlib.sha256(
            (DATA / "gecp_census_resolution.jsonl").read_bytes()
        ).hexdigest()
    )
    assert ladder["source_record_count"] == 3
    assert len(ladder["comparisons"]) == 3
    assert all(
        case["controls_compatible_except_grid_order"] for case in ladder["comparisons"]
    )
    assert all(case["status"] == "passed" for case in ladder["selected_core_checks"])
    for case in report["selected_core_checks"]:
        assert case["status"] == "passed"
    for case in report["held_out_cross_checks"]:
        assert case["certification_status"] == "not_certified"
        assert "denser_GECP_run" in case["not_claimed"]


def test_historical_census_artifacts_remain_unchanged() -> None:
    expected = {
        "gecp_census.jsonl": (
            "a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179"
        ),
        "gecp_census.png": (
            "566469437d18fead51f3e1851046805a7b31526f1bb0c0f54bbe97ca0a956c2d"
        ),
    }
    for name, digest in expected.items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest


def test_shared_tolerance_traces_preserve_historical_pivots_and_ranks() -> None:
    current = {
        (float(record["cutoff"]), float(record["tolerance"])): record
        for record in _records("gecp_census_v2.jsonl")
    }
    historical = _records("gecp_census.jsonl")
    assert len(historical) == 21
    for record in historical:
        shared = current[(float(record["cutoff"]), float(record["tolerance"]))]
        assert shared["rank"] == record["rank"]
        assert shared["pivots"] == record["pivots"]
