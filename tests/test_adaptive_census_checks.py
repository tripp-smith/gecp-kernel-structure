"""Independent adaptive-census checks detect altered samples and decisions."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
from decimal import ROUND_UP, Decimal, Inexact, localcontext
from pathlib import Path
from typing import Any

import mpmath as mp
import pytest

from experiments.adaptive_census import SOURCES, _run_case
from experiments.check_adaptive_census import (
    _batched_lu_solve,
    _mirrored_marks,
    _provenance,
    check_adaptive_census,
    check_final_cross,
    check_selected_cores,
    evaluate_cross,
    validate_payload,
)
from kernelgecp.adaptive_grid import AdaptiveGridConfig


def _record(rank: int = 2) -> dict[str, Any]:
    """Analytic two-node original-kernel fixture, independent of the controller."""
    with mp.workprec(192):
        a = 1 / (1 + mp.e)
        b = mp.e / (1 + mp.e)
        determinant = a * a - b * b
        return {
            "schema_version": 3,
            "algorithm": "residual_adaptive_grid",
            "cutoff": "1",
            "tolerance": "0.1",
            "precision_bits": 128,
            "rank": rank,
            "residual": ["1", "0.8", "0.01"][rank],
            "residual_history": ["1", "0.8", "0.01"][: rank + 1],
            "t_nodes": ["0", "1"][:rank],
            "omega_nodes": ["-1", "1"][:rank],
            "pivots": [mp.nstr(a, 55), mp.nstr(determinant / a, 55)][:rank],
            "core_det_history": [mp.nstr(a, 55), mp.nstr(determinant, 55)][:rank],
            "core_sigma_min": [mp.nstr(a, 55), mp.nstr(abs(a - b), 55)][:rank],
            "tie_counts": [1] * rank,
            "converged": False,
            "stop_reason": "pilot",
        }


@pytest.fixture(scope="module")
def pilot_payload() -> dict[str, Any]:
    config = AdaptiveGridConfig(
        tolerances=(0.1,), max_rank=8, max_stages=3, max_nodes=100, max_marked=1
    )
    configuration = {
        "schema_version": 3,
        "cutoffs": [1],
        "tolerances": [0.1],
        "precision_bits": 128,
        "check_precision_bits": 192,
        "generation_margin": 0.25,
        "max_rank": 8,
        "max_stages": 3,
        "max_nodes": 100,
        "max_marked": 1,
        "workers": 1,
        "require_clean_revision": False,
    }
    config_hash = hashlib.sha256(
        json.dumps(configuration, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    revision = "1" * 40
    environment = {"python": "test", "numpy": "test", "mpmath": "test", "scipy": "test"}
    provenance = {"git_dirty": True, "environment": environment}
    case = _run_case(1.0, config, revision, config_hash, provenance)
    payload = {
        "schema_version": 3,
        "kind": "residual_adaptive_sampled_grid",
        "configuration": configuration,
        "config_hash": config_hash,
        "git_commit": revision,
        **provenance,
        "lock_sha256": "2" * 64,
        "algorithm_source_sha256": {source: "3" * 64 for source in SOURCES},
        "cases": [case],
    }
    return json.loads(json.dumps(payload))


def test_independent_cores_recompute_signed_products_and_singular_values() -> None:
    report = check_selected_cores(_record())
    assert report["status"] == "passed"
    assert [check["rank"] for check in report["checks"]] == [1, 2]
    assert mp.mpf(report["checks"][-1]["determinant"]) < 0
    bad = _record()
    bad["core_det_history"][-1] = "-0.1"
    assert check_selected_cores(bad)["status"] == "disagreement"


def test_independent_factorization_handles_noninvolutive_permutation() -> None:
    with mp.workprec(192):
        matrix = mp.matrix([[0, 0, 2], [4, 0, 3], [1, 5, 6]])
        rhs = mp.matrix([[1, 2], [3, 4], [5, 6]])
        solved = _batched_lu_solve(matrix, rhs)
        assert mp.norm(matrix * solved - rhs) < mp.mpf("1e-55")


def test_rank_zero_check_preserves_failed_sample_and_has_no_solve() -> None:
    report = evaluate_cross(_record(0), ["0.2"], ["0.3"])
    assert report["solve_method"] == "not_applicable_rank_zero"
    assert not report["within_requested_tolerance_on_sample"]
    assert report["certification_status"] == "not_certified"
    with pytest.raises(ValueError, match="outside physical"):
        evaluate_cross(_record(), ["1.1"], ["0"])


def test_pilot_structure_proves_real_nested_stages_without_full_scope_claim(
    pilot_payload: dict[str, Any],
) -> None:
    report = validate_payload(pilot_payload, require_full=False, check_revision=False)
    assert report["status"] == "passed"
    assert not report["full_requested_coverage"]
    assert report["cases"][0]["real_refinement_count"] >= 1
    assert report["all_cases_converged"]
    with pytest.raises(ValueError, match="all seven"):
        validate_payload(pilot_payload, check_revision=False)


@pytest.mark.parametrize(
    "alteration",
    [
        "nested",
        "axis",
        "score",
        "marks",
        "validation",
        "provenance",
        "pass",
        "budget",
        "geometry",
        "stop",
    ],
)
def test_structure_rejects_tampered_decisions_and_sample_leakage(
    pilot_payload: dict[str, Any],
    alteration: str,
) -> None:
    payload = copy.deepcopy(pilot_payload)
    case = payload["cases"][0]
    stage = case["stages"][0]
    if alteration == "nested":
        stage["nested"] = False
    elif alteration == "axis":
        stage["t_grid"].pop(1)
    elif alteration == "score":
        next(iter(stage["probe"]["t_scores"].values()))["score"] = "1e99"
    elif alteration == "marks":
        stage["marked_t"] = []
    elif alteration == "validation":
        case["prefixes"][0]["validation"]["t_nodes"][0] = stage["t_grid"][0]
    elif alteration == "provenance":
        payload["configuration"]["max_rank"] += 1
    elif alteration == "pass":
        stage["probe"]["probe_passed"] = not stage["probe"]["probe_passed"]
    elif alteration == "budget":
        stage["budgets"]["max_rank"] += 1
    elif alteration == "geometry":
        case["records"][0]["pivot_geometry"][0]["dyadic_scale"] += 1
    else:
        case["stages"][-1]["stop_reason"] = "max_stages"
    with pytest.raises(ValueError):
        validate_payload(payload, require_full=False, check_revision=False)


def test_independent_final_cross_uses_exact_sevenths_and_rejects_false_pass(
    pilot_payload: dict[str, Any],
) -> None:
    case = pilot_payload["cases"][0]
    report = check_final_cross(
        case["records"][0], case["final_panels"], case["prefixes"][0]["validation"]
    )
    assert report["status"] == "passed"
    assert report["used_for_refinement_or_stopping"] is False
    bad = copy.deepcopy(case["prefixes"][0]["validation"])
    bad["passed"] = not bad["passed"]
    with pytest.raises(ValueError, match="pass flag"):
        check_final_cross(case["records"][0], case["final_panels"], bad)


def test_mark_ties_preserve_decimal_order_beyond_default_context() -> None:
    panels = [
        {
            "identifier": f"p{i}",
            "left": i / 4,
            "right": (i + 1) / 4,
            "depth": 0,
            "parent": None,
        }
        for i in range(4)
    ]
    scores = {
        "p0": {"score": "1.00000000000000000000000000001"},
        "p1": {"score": "1.00000000000000000000000000002"},
        "p2": {"score": "0"},
        "p3": {"score": "0"},
    }
    assert _mirrored_marks(panels, scores, "t", 1, Decimal("0.5")) == ["p1", "p2"]


def test_checker_decisions_ignore_ambient_decimal_precision_rounding_and_traps(
    pilot_payload: dict[str, Any],
) -> None:
    expected = validate_payload(pilot_payload, require_full=False, check_revision=False)
    with localcontext() as context:
        context.prec = 2
        context.rounding = ROUND_UP
        context.traps[Inexact] = True
        observed = validate_payload(
            pilot_payload, require_full=False, check_revision=False
        )
        assert observed == expected
        assert check_selected_cores(_record())["status"] == "passed"


@pytest.mark.parametrize("budget", ["max_stages", "max_nodes", "max_rank"])
def test_budget_exhaustion_remains_a_protocol_failure(
    pilot_payload: dict[str, Any],
    budget: str,
) -> None:
    payload = copy.deepcopy(pilot_payload)
    configuration = payload["configuration"]
    configuration[budget] = {"max_stages": 1, "max_nodes": 9, "max_rank": 0}[budget]
    payload["config_hash"] = hashlib.sha256(
        json.dumps(configuration, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    controls = {
        key: value
        for key, value in configuration.items()
        if key not in ("schema_version", "cutoffs", "workers", "require_clean_revision")
    }
    controls["tolerances"] = tuple(controls["tolerances"])
    case = _run_case(
        1.0,
        AdaptiveGridConfig(**controls),
        payload["git_commit"],
        payload["config_hash"],
        {"git_dirty": True, "environment": payload["environment"]},
    )
    payload["cases"] = [case]
    payload = json.loads(json.dumps(payload))
    report = validate_payload(payload, require_full=False, check_revision=False)
    assert not report["all_cases_converged"]
    assert report["cases"][0]["stop_reason"] == budget
    prefix = payload["cases"][0]["prefixes"][0]
    assert not prefix["controller_passed"]
    assert not prefix["converged"]
    assert not payload["cases"][0]["records"][0]["converged"]
    if budget != "max_rank":
        assert prefix["fit_passed"]


def test_source_hashes_verify_recorded_git_blobs_not_current_files(
    pilot_payload: dict[str, Any],
) -> None:
    payload = copy.deepcopy(pilot_payload)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    hashes = {}
    for path in ("kernelgecp/kernels.py", "uv.lock"):
        blob = subprocess.check_output(["git", "show", f"{revision}:{path}"])
        hashes[path] = hashlib.sha256(blob).hexdigest()
    payload.update(
        git_commit=revision,
        git_dirty=False,
        lock_sha256=hashes["uv.lock"],
        algorithm_source_sha256={
            "kernelgecp/kernels.py": hashes["kernelgecp/kernels.py"]
        },
    )
    assert _provenance(payload, check_revision=True)["revision_blobs_checked"]
    payload["algorithm_source_sha256"]["kernelgecp/kernels.py"] = "0" * 64
    with pytest.raises(ValueError, match="hash disagrees with revision"):
        _provenance(payload, check_revision=True)


def test_checker_repeat_bytes_and_cli_help(
    tmp_path: Path,
    pilot_payload: dict[str, Any],
) -> None:
    source = tmp_path / "pilot.json"
    source.write_text(json.dumps(pilot_payload), encoding="utf-8")
    a = check_adaptive_census(source, require_full=False, check_revision=False)
    b = check_adaptive_census(source, require_full=False, check_revision=False)
    assert a == b
    assert a["status"] == "passed"
    completed = subprocess.run(
        [sys.executable, "-m", "experiments.check_adaptive_census", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--allow-partial" in completed.stdout
