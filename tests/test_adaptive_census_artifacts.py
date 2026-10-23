"""Audit committed DP coverage/provenance without rerunning GECP or cross solves."""

import hashlib
import json
import math
import re
import struct
import subprocess
import tomllib
from decimal import Decimal, localcontext
from pathlib import Path

import pytest

from experiments.check_adaptive_census import validate_payload
from experiments.plot_adaptive_census import normalized_curves

ROOT = Path(__file__).parents[1]
DATA = ROOT / "experiments/data"
GENERATION_REVISION = "3cee4186fddf990abdd33dbe24c46ff908b1923f"
# Immutable historical inputs. Future source/dependency changes do not rewrite
# the experiment's generation manifest; shallow/archive exports may lack blobs.
GENERATION_SOURCES = {
    "kernelgecp/adaptive_grid.py": (
        "9a728caf111c271f03f5fa302c855786d362c1f410b05d64946747b116bb6ee5"
    ),
    "kernelgecp/adaptive_grids.py": (
        "238b342b957b9fba7c12db8ee4c04b06484e9d64a11e7e786da32ccc62237d1c"
    ),
    "kernelgecp/gecp.py": (
        "13a04be8395adc6db6c2bf89072880f8188e1d98e015895869198524df7f4e43"
    ),
    "kernelgecp/census.py": (
        "61ccd6f90b07fd5e971fb40c02586027857c0ec241374a5ae09c0ae156e5ea22"
    ),
    "kernelgecp/kernels.py": (
        "2150cc70834d2a20978f44cfb52be627d2e0d26503b029751790d137271b6a2b"
    ),
    "kernelgecp/types.py": (
        "4eb7ae7d3158a1ca5d2160a39486ebb23d4b6905a807687b8b63e78f593cda2e"
    ),
    "experiments/adaptive_census.py": (
        "620207b40e0e149532ba6207e9d97f7a0ba7c1bdafc0ab7a845764a5eaf70a46"
    ),
}
GENERATION_LOCK = "ec9204011e82b7d4bc7e1f626d4449ee7ed86d0c3e3e6375b90eea7fc026969c"
CHECKER_SOURCE = "711ac118e26cdab256ed43a3918f876fc7d00e3b0245885dfe994fb761af1316"
DATA_SHA256 = "730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c"
REPORT_SHA256 = "4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03"
PLOT_SHA256 = "dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6"
CUTOFFS = {Decimal(10) ** exponent for exponent in range(7)}
TOLERANCES = {Decimal(f"1e-{exponent}") for exponent in (6, 8, 10, 12)}


def _revision_blob(path: str) -> bytes | None:
    try:
        result = subprocess.run(
            ["git", "show", f"{GENERATION_REVISION}:{path}"],
            cwd=ROOT,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    return result.stdout if result.returncode == 0 else None


def _canonical_json(path: Path) -> dict:
    contents = path.read_bytes()
    payload = json.loads(contents)
    expected = (
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
        + b"\n"
    )
    assert contents == expected
    return payload


def _decimals(values: list[str]) -> list[Decimal]:
    assert all(isinstance(value, str) for value in values)
    numbers = [Decimal(value) for value in values]
    assert all(number.is_finite() for number in numbers)
    return numbers


@pytest.fixture(scope="module")
def payload() -> dict:
    assert hashlib.sha256((DATA / "adaptive_census.json").read_bytes()).hexdigest() == (
        DATA_SHA256
    )
    return _canonical_json(DATA / "adaptive_census.json")


@pytest.fixture(scope="module")
def report() -> dict:
    assert (
        hashlib.sha256((DATA / "adaptive_census_checks.json").read_bytes()).hexdigest()
        == REPORT_SHA256
    )
    return _canonical_json(DATA / "adaptive_census_checks.json")


def test_committed_adaptive_census_has_all_twenty_eight_converged_records(
    payload: dict,
) -> None:
    assert payload["schema_version"] == 3
    assert payload["kind"] == "residual_adaptive_sampled_grid"
    assert payload["git_commit"] == GENERATION_REVISION
    assert payload["git_dirty"] is False
    assert payload["certification_status"] == "not_certified"
    assert payload["error_scope"] == (
        "finite_training_probe_and_independent_validation_sets_only"
    )
    cases = payload["cases"]
    assert len(cases) == 7
    assert {Decimal(case["cutoff"]) for case in cases} == CUTOFFS
    pairs = set()
    for case in cases:
        assert case["converged"] is True
        assert case["stop_reason"] == "sampled_tolerances"
        assert len(case["records"]) == len(case["prefixes"]) == 4
        assert {Decimal(p["tolerance"]) for p in case["prefixes"]} == TOLERANCES
        for record, prefix in zip(case["records"], case["prefixes"], strict=True):
            pairs.add((Decimal(record["cutoff"]), Decimal(record["tolerance"])))
            assert record["converged"] is prefix["converged"] is True
            assert prefix["controller_passed"] is prefix["fit_passed"] is True
            assert prefix["joint_passed"] is prefix["validation"]["passed"] is True
            assert record["rank"] == prefix["chosen_rank"]
            assert record["tolerance"] == prefix["tolerance"]
            assert prefix["training_rank"] <= prefix["chosen_rank"]
            assert prefix["tested_ranks"] == list(
                range(prefix["training_rank"], prefix["chosen_rank"] + 1)
            )
            assert record["stop_reason"] == prefix["stop_reason"] == "sampled_tolerance"
    assert pairs == {
        (cutoff, tolerance) for cutoff in CUTOFFS for tolerance in TOLERANCES
    }


def test_generation_hashes_bind_recorded_revision_not_future_source(
    payload: dict,
) -> None:
    assert re.fullmatch(r"[0-9a-f]{40}", GENERATION_REVISION)
    assert payload["algorithm_source_sha256"] == GENERATION_SOURCES
    assert payload["lock_sha256"] == GENERATION_LOCK
    for path, digest in {**GENERATION_SOURCES, "uv.lock": GENERATION_LOCK}.items():
        blob = _revision_blob(path)
        if blob is not None:
            assert hashlib.sha256(blob).hexdigest() == digest
    config_bytes = json.dumps(
        payload["configuration"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert hashlib.sha256(config_bytes).hexdigest() == payload["config_hash"]
    config_blob = _revision_blob("experiments/adaptive_census.toml")
    if config_blob is not None:
        assert tomllib.loads(config_blob.decode()) == payload["configuration"]
    assert all(
        isinstance(payload["environment"][name], str) and payload["environment"][name]
        for name in ("python", "numpy", "scipy", "mpmath", "mpmath_backend")
    )


def test_each_case_records_genuine_nested_refinement_and_axis_retention(
    payload: dict,
) -> None:
    for case in payload["cases"]:
        stages = case["stages"]
        assert 2 <= len(stages) <= payload["configuration"]["max_stages"]
        assert stages[-1]["consecutive_passes"] >= 2
        assert stages[-1]["passed"] is stages[-2]["passed"] is True
        assert stages[-1]["stop_reason"] == "strictest_sampled_stages_passed"
        for index, stage in enumerate(stages):
            assert stage["index"] == index and stage["nested"] is True
            assert stage["certification_status"] == "not_certified"
            assert stage["coordinate_precision_bits"] == 53
            assert stage["coordinate_semantics"] == "repr_float64_decimal"
            assert stage["grid_shape"] == [
                len(stage["t_grid"]),
                len(stage["omega_grid"]),
            ]
            rank = stage["gecp"]["rank"]
            high = stage["gecp"]["high_precision"]
            for field in (
                "pivots",
                "t_nodes",
                "omega_nodes",
                "core_det_history",
                "core_sigma_min_history",
                "residual_history",
            ):
                assert len(high[field]) == rank + (field == "residual_history")
                _decimals(high[field])
            assert stage["gecp"]["metadata"]["precision_bits"] == 128
            assert stage["gecp"]["metadata"]["certification_status"] == "not_certified"
            assert stage["probe"]["precision_bits"] == 192
            _decimals(
                [
                    stage["probe"][name]
                    for name in ("training_max", "probe_max", "surplus_max")
                ]
            )
            for axis, name in (("t", "t_grid"), ("omega", "omega_grid")):
                nodes = _decimals(stage[name])
                assert nodes == sorted(set(nodes))
                assert len(nodes) <= payload["configuration"]["max_nodes"]
                assert set(stage["probe"][f"{axis}_scores"]) == {
                    panel["identifier"] for panel in stage[f"{axis}_panels"]
                }
                for score in stage["probe"][f"{axis}_scores"].values():
                    assert all(value >= 0 for value in _decimals(list(score.values())))
                if index:
                    assert set(stages[index - 1][name]).issubset(stage[name])
            if index:
                previous = stages[index - 1]
                assert previous["refinement_applied"] is True
                assert previous["marked_t"] or previous["marked_omega"]
                assert (
                    stage["t_grid"] != previous["t_grid"]
                    or stage["omega_grid"] != previous["omega_grid"]
                )
        assert case["final_axes"] == {
            "t": stages[-1]["t_grid"],
            "omega": stages[-1]["omega_grid"],
        }
        assert case["final_panels"] == {
            "t": stages[-1]["t_panels"],
            "omega": stages[-1]["omega_panels"],
        }


def test_final_coordinates_are_untouched_by_every_training_and_probe_stage(
    payload: dict,
) -> None:
    for case in payload["cases"]:
        for prefix in case["prefixes"]:
            validation = prefix["validation"]
            assert validation["precision_bits"] == 192
            assert validation["certification_status"] == "not_certified"
            assert validation["scope"] == (
                "untouched_final_seventh_fraction_samples_after_refinement_stops"
            )
            assert Decimal(validation["max_residual"]) <= Decimal(prefix["tolerance"])
            for axis, name in (("t", "t_grid"), ("omega", "omega_grid")):
                values = _decimals(validation[f"{axis}_nodes"])
                assert len(values) == 6 * len(case["final_panels"][axis])
                assert values == sorted(set(values))
                validation_bits = set(map(float, validation[f"{axis}_nodes"]))
                for stage in case["stages"]:
                    assert validation_bits.isdisjoint(map(float, stage[name]))
                    assert validation_bits.isdisjoint(
                        map(float, stage["probe"][f"{axis}_nodes"])
                    )


def test_full_prefix_histories_geometry_and_numerical_scope_are_serialized(
    payload: dict,
) -> None:
    for case in payload["cases"]:
        for record, prefix in zip(case["records"], case["prefixes"], strict=True):
            rank = record["rank"]
            assert record["schema_version"] == 3
            assert record["algorithm"] == "residual_adaptive_grid"
            assert record["git_commit"] == payload["git_commit"]
            assert record["config_hash"] == payload["config_hash"]
            assert record["git_dirty"] is False
            assert record["environment"] == payload["environment"]
            assert record["precision_bits"] == 128
            assert record["grid_coordinate_precision_bits"] == 53
            assert record["grid_generation"] == (
                "float64_nested_dyadic_panel_lattice_residual_refinement"
            )
            assert record["grid_shape"] == [
                len(case["final_axes"][axis]) for axis in ("t", "omega")
            ]
            for field in (
                "pivots",
                "pivot_magnitudes",
                "t_nodes",
                "omega_nodes",
                "core_det_history",
                "core_sigma_min",
            ):
                assert len(record[field]) == rank
                values = _decimals(record[field])
                if field in ("pivots", "core_det_history"):
                    assert all(value != 0 for value in values)
                elif field == "core_sigma_min":
                    assert all(value > 0 for value in values)
            history = _decimals(record["residual_history"])
            assert len(history) == rank + 1 and all(value >= 0 for value in history)
            assert (
                record["residual"]
                == record["residual_history"][-1]
                == prefix["training_max"]
            )
            assert Decimal(record["residual"]) <= Decimal(record["tolerance"])
            assert [value.copy_abs() for value in _decimals(record["pivots"])] == (
                _decimals(record["pivot_magnitudes"])
            )
            assert len(record["tie_counts"]) == len(record["pivot_geometry"]) == rank
            assert all(
                type(count) is int and count >= 1 for count in record["tie_counts"]
            )
            for geometry in record["pivot_geometry"]:
                for field in ("t_boundary_distances", "omega_boundary_distances"):
                    assert len(geometry[field]) == 2
                    assert all(value >= 0 for value in _decimals(geometry[field]))
                assert 0 <= geometry["reflection_nearest_index"] < rank
                _decimals(
                    [
                        geometry[name]
                        for name in (
                            "reflection_t",
                            "reflection_omega",
                            "reflection_distance",
                        )
                    ]
                )
            metadata = record["metadata"]
            assert metadata["converged"] is True
            assert metadata["arithmetic"] == "mpmath"
            assert (
                metadata["precision_bits"]
                == metadata["requested_precision_bits"]
                == 128
            )
            assert metadata["certification_status"] == "not_certified"
            assert metadata["certification_scope"] == "sampled_only_no_continuous_bound"
            assert (
                metadata["convergence_scope"]
                == "sampled_training_probe_and_independent_validation"
            )
            assert metadata["tie_counts"] == record["tie_counts"]
            assert metadata["tie_policy"] == "lexicographic_near_maximum"
            assert len(metadata["conditioning_history"]) == rank
            assert (
                metadata["conditioning_method"]
                == "selected_core_smallest_singular_value_only"
            )
            for diagnostic in metadata["conditioning_history"]:
                _decimals([diagnostic["sigma_min_decimal"]])
            assert metadata["configuration"]["check_precision_bits"] == 192
            assert (
                metadata["configuration"]["final_validation_used_for_refinement"]
                is False
            )


def test_checker_report_binds_data_and_preserves_all_twenty_eight_cross_checks(
    payload: dict,
    report: dict,
) -> None:
    assert report["schema_version"] == 1
    assert report["kind"] == "independent_residual_adaptive_sampled_grid_checks"
    assert (
        report["input_sha256"]
        == hashlib.sha256((DATA / "adaptive_census.json").read_bytes()).hexdigest()
    )
    assert report["checker_source_sha256"] == CHECKER_SOURCE
    blob = _revision_blob("experiments/check_adaptive_census.py")
    if blob is not None:
        assert hashlib.sha256(blob).hexdigest() == CHECKER_SOURCE
    assert report["precision_bits"] == 192
    assert report["status"] == "passed"
    assert report["certification_status"] == "not_certified"
    structure = report["structure"]
    assert (
        structure["status"] == "passed" and structure["full_requested_coverage"] is True
    )
    assert structure["case_count"] == 7 and structure["record_count"] == 28
    assert structure["all_cases_converged"] is True
    assert structure["provenance"]["git_commit"] == GENERATION_REVISION
    assert structure["provenance"]["config_hash"] == payload["config_hash"]
    assert structure["provenance"]["source_hashes"] == {
        **GENERATION_SOURCES,
        "uv.lock": GENERATION_LOCK,
    }
    assert structure["provenance"]["recorded_environment"] == payload["environment"]
    cases = {Decimal(case["cutoff"]): case for case in payload["cases"]}
    assert len(report["selected_cores"]) == 7
    assert {Decimal(case["cutoff"]) for case in report["selected_cores"]} == CUTOFFS
    for core in report["selected_cores"]:
        case = cases[Decimal(core["cutoff"])]
        deepest = max(case["records"], key=lambda record: record["rank"])
        rank = deepest["rank"]
        assert core["source_rank"] == rank
        assert [check["rank"] for check in core["checks"]] == sorted(
            {1, max(1, rank // 2), rank}
        )
        assert core["status"] == "passed"
        assert core["certification_status"] == "not_certified"
        assert core["scope"] == "original_kernel_selected_cores_first_middle_final"
        for check in core["checks"]:
            assert check["status"] == "passed"
            determinant, sigma, condition = _decimals(
                [
                    check[name]
                    for name in ("determinant", "sigma_min", "condition_number_2")
                ]
            )
            assert determinant != 0 and sigma > 0 and condition >= 1
            for name, bound in (
                ("stored_determinant_relative_error", "1e-18"),
                ("signed_pivot_product_relative_error", "1e-18"),
                ("stored_sigma_relative_error", "1e-15"),
            ):
                (error,) = _decimals([check[name]])
                assert 0 <= error <= Decimal(bound)
    crosses = report["final_crosses"]
    assert len(crosses) == 28
    assert {
        (Decimal(check["cutoff"]), Decimal(check["tolerance"])) for check in crosses
    } == {(cutoff, tolerance) for cutoff in CUTOFFS for tolerance in TOLERANCES}
    for cross in crosses:
        case = cases[Decimal(cross["cutoff"])]
        record, prefix = next(
            (record, prefix)
            for record, prefix in zip(case["records"], case["prefixes"], strict=True)
            if Decimal(record["tolerance"]) == Decimal(cross["tolerance"])
        )
        validation = prefix["validation"]
        assert cross["rank"] == record["rank"] == prefix["chosen_rank"]
        assert cross["precision_bits"] == 192
        assert cross["sample_shape"] == [
            len(validation[f"{axis}_nodes"]) for axis in ("t", "omega")
        ]
        assert cross["sample_shape"] == [
            6 * len(case["final_panels"][axis]) for axis in ("t", "omega")
        ]
        assert cross["stored_max_abs_residual"] == validation["max_residual"]
        assert cross["status"] == "passed"
        assert cross["within_requested_tolerance_on_sample"] is True
        assert cross["used_for_refinement_or_stopping"] is False
        assert cross["certification_status"] == "not_certified"
        assert cross["scope"] == (
            "stored_selected_cross_on_all_final_leaf_seventh_fraction_points"
        )
        maximum, stored, disagreement = _decimals(
            [
                cross[name]
                for name in (
                    "max_abs_residual",
                    "stored_max_abs_residual",
                    "stored_absolute_disagreement",
                )
            ]
        )
        assert 0 <= maximum <= Decimal(cross["tolerance"])
        assert stored >= 0 and disagreement >= 0
        with localcontext() as context:
            context.prec = 100
            bound = max(Decimal("1e-45"), maximum * Decimal("1e-25"))
            assert disagreement <= bound
            assert abs(stored - maximum) <= bound
        assert len(cross["worst_coordinates"]) == 2
        assert cross["worst_coordinates"][0] in validation["t_nodes"]
        assert cross["worst_coordinates"][1] in validation["omega_nodes"]
        assert cross["solve_method"] == (
            "mpmath_lu_once_batched_forward_backward_substitution"
        )
        (solve_residual,) = _decimals([cross["solve_relative_max_entry_residual"]])
        assert solve_residual >= 0


def test_full_structural_audit_does_not_rerun_numerical_experiments(
    payload: dict,
) -> None:
    # This validates recorded geometry/decisions only, never GECP, SVD or cross
    # evaluation; historical-blob availability is handled separately above.
    summary = validate_payload(payload, check_revision=False)
    assert summary["status"] == "passed" and summary["all_cases_converged"]
    assert all(case["real_refinement_count"] >= 1 for case in summary["cases"])


def test_committed_plot_has_seven_true_normalized_curves_and_full_png(
    payload: dict,
) -> None:
    curves = normalized_curves(payload)
    assert len(curves) == 7
    assert {Decimal(str(curve.cutoff)) for curve in curves} == CUTOFFS
    for curve in curves:
        assert curve.converged and all(curve.validation_passed)
        assert curve.stage_count >= 2
        assert curve.ranks == tuple(range(len(curve.residuals)))
        assert curve.normalized_ranks == tuple(
            rank / math.log1p(curve.cutoff) for rank in curve.ranks
        )
        assert len(curve.validation_errors) == len(curve.validation_ranks) == 4
    png = (DATA / "adaptive_census.png").read_bytes()
    assert hashlib.sha256(png).hexdigest() == PLOT_SHA256
    assert png[:8] == b"\x89PNG\r\n\x1a\n"
    assert png[12:16] == b"IHDR"
    assert struct.unpack(">II", png[16:24]) == (2240, 1440)
    assert len(png) > 10_000
