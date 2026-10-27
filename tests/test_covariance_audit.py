"""False-pass attacks on scientific EI receipt consumption."""

from __future__ import annotations

import copy
import subprocess

import numpy as np
import pytest

from experiments import covariance_benchmark as benchmark
from scripts.audit_covariance import (
    CovarianceAuditError,
    artifact,
    audit_case,
    bind_source_to_git,
    finite_number,
    read_json,
)


def _case():
    rng = np.random.default_rng(735)
    data = {
        "values": rng.normal(size=(benchmark.TEST_END, 6)),
        "timestamps": [str(j) for j in range(benchmark.TEST_END)],
        "ordinals": np.arange(benchmark.TEST_END, dtype=np.int64),
    }
    split = benchmark.split_data(data, 512, "test")
    noise = 0.8
    root = np.ones((512, 1))
    # The huge declared lengthscale rounds every kernel entry to exactly one.
    solution = (split["y"] - np.sum(split["y"]) / (512 + noise)) / noise
    mean = np.full(512, np.sum(solution))
    std = np.full(512, np.sqrt(noise + noise / (512 + noise)))
    residual = np.sum(solution) + noise * solution - split["y"]
    error = float(np.linalg.norm(residual) / np.linalg.norm(split["y"]))
    provenance = {
        "source_revision": "a" * 40,
        "source_clean_and_tracked": True,
        "source_hashes": {"example": "hash"},
        "versions": {"numpy": "2"},
        "threads": {key: "1" for key in benchmark.THREADS},
    }
    frozen = {
        "standard_deviation_definition": benchmark.STD_DEFINITION,
        "provenance": provenance,
        "dataset": {"identity": "unit-test"},
        "selected_ranks": {"512": 32},
        "lengthscale": 1e100,
        "noise": noise,
    }
    import json

    receipt = {
        "standard_deviation_definition": benchmark.STD_DEFINITION,
        "schema": "kernelgecp.covariance-case.v1",
        "status": "passed",
        "canonical": True,
        "size": 512,
        "repeat": 0,
        "actual_rank": 1,
        "rank_budget": 32,
        "provenance": provenance,
        "dataset": frozen["dataset"],
        "frozen_calibration_sha256": benchmark.sha256(
            json.dumps(frozen, sort_keys=True).encode()
        ),
        "checks": {
            key: True
            for key in (
                "dense_accuracy",
                "approximate_solve",
                "original_kernel_cg",
                "retained_state",
                "rss",
                "source_clean",
            )
        },
        "normalization": split["normalization"],
        "split": {
            key: split[key]
            for key in (
                "train_ordinal_sha256",
                "query_ordinal_sha256",
                "train_first",
                "train_last",
                "query_first",
                "query_last",
            )
        },
        "solve_relative_residual": error,
        "cg_relative_residual": error,
        "cg_info": 0,
        "accuracy": {
            key: 0.0 for key in ("mean_rms", "mean_max", "std_rms", "std_max")
        },
        "retained_numerical_bytes": 10000,
        "retained_components_bytes": {"state": 10000},
        "peak_rss_bytes": 100000000,
        "timing_seconds": {
            "construction": 0.1,
            "prediction": 0.1,
            "dense_reference": 0.2,
            "original_kernel_cg": 0.1,
            "full_workflow": 0.5,
        },
    }
    arrays = {
        "root": root,
        "mean": mean,
        "std": std,
        "expected_mean": mean.copy(),
        "expected_std": std.copy(),
        "query_y": split["query_y"],
        "cg_solution": solution,
        "cg_residual": residual,
        "cg_rhs": split["y"],
        "pivots": np.array([0]),
        "approximate_solution": solution.copy(),
    }
    return receipt, arrays, frozen, data


def test_independent_receipt_audit_accepts_a_valid_finite_rank_model():
    receipt, arrays, frozen, data = _case()
    result = audit_case(receipt, arrays, frozen, data, "a" * 40)
    assert result["solve_residual"] < 1e-10


@pytest.mark.parametrize(
    "attack",
    [
        "wrong_source",
        "pilot",
        "fake_flag",
        "normalization",
        "wrong_rank",
        "inaccurate_solve",
        "false_cg",
        "prediction",
        "wrong_data",
        "rss",
        "nan_timing",
    ],
)
def test_receipt_audit_rejects_stale_source_and_numerical_false_passes(attack):
    receipt, arrays, frozen, data = _case()
    if attack == "wrong_source":
        receipt["provenance"]["source_revision"] = "b" * 40
    elif attack == "pilot":
        receipt["canonical"] = False
    elif attack == "fake_flag":
        receipt["checks"]["approximate_solve"] = 1
    elif attack == "normalization":
        receipt["normalization"]["fit_rows"] = 16384
    elif attack == "wrong_rank":
        receipt["rank_budget"] = 64
    elif attack == "inaccurate_solve":
        arrays["approximate_solution"][:] = 0
    elif attack == "false_cg":
        arrays["cg_solution"][:] = 0
        arrays["cg_residual"][:] = 0
        receipt["cg_relative_residual"] = 0
    elif attack == "prediction":
        arrays["mean"] += 0.01
    elif attack == "wrong_data":
        arrays["query_y"] = np.zeros(512)
    elif attack == "rss":
        receipt["peak_rss_bytes"] = 5 * 1024**3
    elif attack == "nan_timing":
        receipt["timing_seconds"]["construction"] = float("nan")
    with pytest.raises(CovarianceAuditError):
        audit_case(receipt, arrays, frozen, data, "a" * 40)


def test_later_repeat_must_reproduce_actual_arrays_not_only_pass_flags():
    receipt, arrays, frozen, data = _case()
    previous = copy.deepcopy(arrays)
    arrays["mean"][0] += 1e-5
    with pytest.raises(CovarianceAuditError, match="fresh repeats"):
        audit_case(receipt, arrays, frozen, data, "a" * 40, previous=previous)


def test_json_rejects_duplicate_keys_nonfinite_constants_and_nonobjects(tmp_path):
    path = tmp_path / "receipt.json"
    for invalid in ('{"status":"failed","status":"passed"}', '{"error":NaN}', "[]"):
        path.write_text(invalid)
        with pytest.raises(CovarianceAuditError):
            read_json(path)


def test_artifacts_require_owning_paths_and_actual_content_hash(tmp_path):
    path = tmp_path / "data.npz"
    path.write_bytes(b"actual data")
    for name, digest in (
        ("../data.npz", "a" * 64),
        (path.name, "b" * 64),
        (path.name, "not-a-hash"),
    ):
        with pytest.raises(CovarianceAuditError):
            artifact(tmp_path, name, digest)
    assert artifact(tmp_path, path.name, benchmark.sha256(path.read_bytes())) == path
    for value in (True, np.inf, np.nan, "1"):
        with pytest.raises(CovarianceAuditError):
            finite_number(value, "test")


def test_source_binding_rejects_invented_revision_hashes_and_missing_keys(monkeypatch):
    from scripts import audit_covariance as audit

    contents = {name: f"frozen {name}".encode() for name in benchmark.SOURCE_FILES}
    record = {
        "source_revision": "a" * 40,
        "source_hashes": {
            name: benchmark.sha256(data) for name, data in contents.items()
        },
    }

    def git_output(command, **kwargs):
        if command[1] == "rev-parse":
            return "a" * 40
        return contents[command[-1].split(":", 1)[1]]

    monkeypatch.setattr(audit.subprocess, "check_output", git_output)
    bind_source_to_git(record, "a" * 40)
    changed = copy.deepcopy(record)
    changed["source_hashes"][benchmark.SOURCE_FILES[0]] = "b" * 64
    with pytest.raises(CovarianceAuditError, match="Git blob"):
        bind_source_to_git(changed)
    missing = copy.deepcopy(record)
    missing["source_hashes"].pop(benchmark.SOURCE_FILES[0])
    with pytest.raises(CovarianceAuditError, match="keyset"):
        bind_source_to_git(missing)
    with pytest.raises(CovarianceAuditError, match="revision mismatch"):
        bind_source_to_git(record, "b" * 40)


def test_source_binding_reads_actual_git_blobs_instead_of_worktree_claims(
    tmp_path, monkeypatch
):
    for name in benchmark.SOURCE_FILES:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"committed {name}".encode())
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "add", "."], cwd=tmp_path, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=Audit fixture",
            "-c",
            "user.email=audit@example.invalid",
            "commit",
            "-qm",
            "frozen blobs",
        ],
        cwd=tmp_path,
        check=True,
    )
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=tmp_path, text=True
    ).strip()
    record = {
        "source_revision": revision,
        "source_hashes": {
            name: benchmark.sha256((tmp_path / name).read_bytes())
            for name in benchmark.SOURCE_FILES
        },
    }
    monkeypatch.setattr(benchmark, "ROOT", tmp_path)
    bind_source_to_git(record, revision)
    changed = benchmark.SOURCE_FILES[0]
    (tmp_path / changed).write_bytes(b"uncommitted replacement")
    bind_source_to_git(record, revision)
    record["source_hashes"][changed] = benchmark.sha256(
        (tmp_path / changed).read_bytes()
    )
    with pytest.raises(CovarianceAuditError, match="Git blob"):
        bind_source_to_git(record, revision)
