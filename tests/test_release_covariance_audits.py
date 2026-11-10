"""Synthetic gate fixtures; no scientific claims or provider runs."""

import copy
import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

from scripts import release_covariance_audits as gate

P = "a" * 40


def case(n, r, frozen):
    x = dict.fromkeys(gate.CASE_KEYS)
    x.update(
        schema="kernelgecp.covariance-case.v1",
        status="passed",
        canonical=True,
        size=n,
        repeat=r,
        actual_rank=1,
        rank_budget=1,
        split=dict.fromkeys(gate.SPLIT, "x"),
        normalization={
            "features_mean": [0] * 5,
            "features_std": [1] * 5,
            "target_mean": 0,
            "target_std": 1,
            "fit_rows": n,
        },
        provenance={
            "source_revision": P,
            "source_hashes": frozen["provenance"]["source_hashes"],
        },
        dataset=frozen["dataset"],
        checks=dict.fromkeys(
            (
                "dense_accuracy",
                "approximate_solve",
                "original_kernel_cg",
                "retained_state",
                "rss",
                "source_clean",
            ),
            True,
        ),
        artifact={"path": f"n{n}-repeat{r}.npz", "sha256": "x"},
        timing_seconds={
            k: (None if k == "dense_reference" and n > 2048 else 1.0)
            for k in (
                "application_workflow",
                "construction",
                "prediction",
                "dense_reference",
                "original_kernel_cg",
                "full_workflow",
            )
        },
        peak_rss_bytes=1,
        retained_numerical_bytes=1,
        retained_components_bytes=dict.fromkeys(
            (
                "shifted_factor",
                "source_values",
                "source_ordinals",
                "normalized_split",
                "prediction_arrays",
                "solve_arrays",
            ),
            1,
        ),
    )
    x["provenance"] = {
        **dict.fromkeys(
            (
                "source_clean_and_tracked",
                "source_binding_scope",
                "development_lock_sha256_context_only",
                "versions",
                "python",
                "python_executable",
                "scientific_runtime_lock",
                "platform",
                "threads",
            )
        ),
        **x["provenance"],
    }
    x["factor_metadata"] = dict.fromkeys(
        (
            "algorithm",
            "arithmetic",
            "certification_scope",
            "certification_status",
            "conditioning_history",
            "conditioning_method",
            "conditioning_status",
            "configuration",
            "converged",
            "convergence_scope",
            "precision_bits",
            "precision_status",
            "refinement",
            "requested_precision_bits",
            "schema_version",
            "stop_reason",
            "tie_counts",
            "tie_policy",
            "tie_rtol",
            "tie_status",
            "tie_tolerance_rule",
        )
    )
    x["factor_metadata"]["configuration"] = dict.fromkeys(
        (
            "column_evaluations",
            "retained_numerical_bytes",
            "returned_diagonal_error",
            "roundoff_slack",
            "threshold",
        )
    )
    x["factor_metadata"]["refinement"] = dict.fromkeys(
        (
            "accepted",
            "evaluations",
            "max_evaluations",
            "message",
            "optimizer_status",
            "optimizer_success",
            "status",
        )
    )
    x["factor_metadata"]["conditioning_history"] = []
    x["accuracy"] = (
        dict.fromkeys(("mean_rms", "mean_max", "std_rms", "std_max"), 0.0)
        if n <= 2048
        else None
    )
    x["retained_numerical_bytes"] = 6
    return x


def fixture(tmp_path, monkeypatch):
    root = tmp_path / "payload"
    full = root / "full"
    full.mkdir(parents=True)
    frozen = {
        "selected_ranks": {str(n): 1 for n in gate.SIZES},
        "provenance": {"source_hashes": {"x": "y"}},
        "dataset": {"station": "fixture"},
    }
    items = []
    for n in gate.SIZES:
        for r in range(3):
            value = case(n, r, frozen)
            array = full / f"n{n}-repeat{r}.npz"
            array.write_bytes(b"synthetic")
            value["artifact"]["sha256"] = gate.digest(array)
            path = array.with_suffix(".json")
            path.write_text(json.dumps(value))
            path.with_suffix(".log").write_text("")
            items.append(
                {
                    "path": path.name,
                    "exit_code": 0,
                    "process_wall_seconds": 1.0,
                    "sha256": gate.digest(path),
                }
            )
    (full / "campaign.json").write_text(
        json.dumps(
            {
                "schema": "kernelgecp.covariance-campaign.v1",
                "provenance": {"source_revision": P, "source_hashes": {"x": "y"}},
                "frozen_path": "/historical",
                "cases": items,
            }
        )
    )
    (root / "n512-repeat0.json").write_bytes((full / "n512-repeat0.json").read_bytes())
    (root / "n512-repeat0.npz").write_bytes(b"synthetic")
    (root / "provider.npz").write_bytes(b"synthetic")
    provider = dict.fromkeys(gate.PROVIDER_KEYS)
    provider.update(
        schema="kernelgecp.covariance-provider.v1",
        status="passed",
        provider="gpytorch",
        size=512,
        dataset=frozen["dataset"],
        provenance={"source_revision": P, "source_hashes": {"x": "y"}},
        artifact={"path": "provider.npz", "sha256": gate.digest(root / "provider.npz")},
        checks=dict.fromkeys(("cpu_float64", "frozen_factor", "numerical"), True),
        metrics=dict.fromkeys(
            (
                "mean_rms",
                "mean_max",
                "std_rms",
                "std_max",
                "solve_max",
                "log_probability_error",
            ),
            0.0,
        ),
    )
    (root / "provider.json").write_text(json.dumps(provider))
    smallrow = {"size": 512, "repeat": 0, "solve_residual": 0.0, "cg_residual": 0.0}
    (root / "small-audit.json").write_text(
        json.dumps(
            {
                "schema": "kernelgecp.covariance-small-audit.v1",
                "status": "passed",
                "source_revision": P,
                "scope": "small_real_CI_case_not_full_EI_campaign_acceptance",
                "case": smallrow,
            }
        )
    )
    fullrows = [
        {**smallrow, "size": n, "repeat": r} for n in gate.SIZES for r in range(3)
    ]
    (root / "full-audit.json").write_text(
        json.dumps(
            {
                "schema": "kernelgecp.covariance-audit.v1",
                "status": "passed",
                "source_revision": P,
                "cases": fullrows,
                "real_data": True,
                "actual_provider": True,
                "independent_full_kernel_cg_recomputed_once_per_size": True,
            }
        )
    )
    monkeypatch.setattr(gate, "npz_headers", lambda *a, **kw: None)
    return root, frozen


def test_full_structural_fixture_and_partial_rejection(tmp_path, monkeypatch):
    root, frozen = fixture(tmp_path, monkeypatch)
    assert len(gate.payload_inventory(root)) == 61
    gate.validate_payload(root, frozen, P)
    path = root / "full/campaign.json"
    manifest = json.loads(path.read_text())
    manifest["cases"][-1] = copy.deepcopy(manifest["cases"][0])
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="duplicate"):
        gate.validate_payload(root, frozen, P)


def test_split_and_provider_checks_fail_closed(tmp_path, monkeypatch):
    root, frozen = fixture(tmp_path, monkeypatch)
    x = case(512, 0, frozen)
    x["split"].pop("query_last")
    with pytest.raises(ValueError, match="split keyset"):
        gate.validate_case(x, 512, 0, frozen, P)
    path = root / "provider.json"
    x = json.loads(path.read_text())
    x["checks"]["numerical"] = False
    path.write_text(json.dumps(x))
    with pytest.raises(ValueError, match="provider checks"):
        gate.validate_payload(root, frozen, P)


def test_strict_json_and_exact_layout(tmp_path, monkeypatch):
    root, _ = fixture(tmp_path, monkeypatch)
    path = root / "provider.json"
    for raw in ('{"a":1,"a":2}', '{"a":1e999}'):
        path.write_text(raw)
        with pytest.raises(ValueError):
            gate.obj(path)
    (root / "extra").mkdir()
    with pytest.raises(ValueError, match="inventory"):
        gate.payload_inventory(root)


def test_npz_preflight_bounds_keys_shapes_and_hardlinks(tmp_path):
    path = tmp_path / "a.npz"
    np.savez_compressed(
        path,
        root=np.zeros((512, 1)),
        mean=np.zeros(512),
        std=np.zeros(512),
        expected_mean=np.zeros(512),
        expected_std=np.zeros(512),
        solution=np.zeros(512),
        rhs=np.zeros(512),
    )
    gate.npz_headers(path, gate.PROVIDER_ARRAYS, 512, 1, provider=True)
    with pytest.raises(ValueError, match="keyset"):
        gate.npz_headers(path, gate.ARRAYS, 512, 1)
    with pytest.raises(ValueError, match="dtype/shape"):
        gate.npz_headers(path, gate.PROVIDER_ARRAYS, 512, 2, provider=True)
    import os

    os.link(path, tmp_path / "link")
    with pytest.raises(ValueError, match="unsafe"):
        gate.npz_headers(path, gate.PROVIDER_ARRAYS, 512, 1, provider=True)


def test_missing_data_never_launches_child(tmp_path, monkeypatch):
    root, _ = fixture(tmp_path, monkeypatch)

    def fail(*args, **kwargs):
        raise AssertionError("child must not run")

    result = gate.audit(
        root,
        P,
        "/env/bin/python",
        tmp_path / "missing.zip",
        tmp_path / "output",
        runner=fail,
    )
    assert result["status"] == "failed"
    assert (tmp_path / "output/receipt.json").is_file()


def test_bounded_runner_kills_timeout_group():
    import os

    with pytest.raises(subprocess.TimeoutExpired):
        gate.bounded_runner(
            ["python", "-c", "import time; time.sleep(5)"],
            cwd=Path.cwd(),
            env=dict(os.environ),
            text=True,
            capture_output=True,
            check=False,
            timeout=0.05,
        )


def test_controller_rejects_child_failure_even_with_valid_json(tmp_path, monkeypatch):
    from contextlib import contextmanager
    from types import SimpleNamespace

    root, _ = fixture(tmp_path, monkeypatch)
    data = tmp_path / "dataset.zip"
    data.write_bytes(b"synthetic")
    monkeypatch.setattr(gate, "dataset_check", lambda _: "synthetic-hash")
    monkeypatch.setattr(gate, "validate_payload", lambda *args: None)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=gate.ROOT, text=True
    ).strip()

    @contextmanager
    def snapshot(*args):
        yield gate.ROOT

    def failed(command, **kwargs):
        assert command[1] == "-I" and "scripts.audit_covariance" in command[3]
        assert "urlopen=no_network" in command[3]
        Path(command[-1]).write_text(
            json.dumps({"schema": "kernelgecp.covariance-audit.v1", "status": "passed"})
        )
        return SimpleNamespace(returncode=1, stdout="test", stderr="failed")

    result = gate.audit(
        root,
        revision,
        "/env/bin/python",
        data,
        tmp_path / "out",
        runner=failed,
        snapshot=snapshot,
        identity_probe=lambda _: {"revision": "b" * 40},
    )
    assert result["status"] == "failed"
    assert "auditor failed" in result["error"]
    assert (tmp_path / "out/science.log").is_file()


def test_runtime_accepts_venv_interpreter_symlink(tmp_path):
    prefix = tmp_path / "venv"
    (prefix / "bin").mkdir(parents=True)
    base = tmp_path / "base-python"
    base.write_text("stub")
    executable = prefix / "bin/python"
    executable.symlink_to(base)
    info = {
        "platform": "linux",
        "machine": "x86_64",
        "python": "3.12.8",
        "versions": {
            "numpy": "2.2.6",
            "scipy": "1.18.1",
            "torch": "2.14.1+cpu",
            "gpytorch": "1.15.2",
            "linear-operator": "0.6.1",
        },
        "threads": dict.fromkeys(gate.THREADS, "1"),
        "prefix": str(prefix),
        "executable": str(executable),
        "modules": {
            name: str(prefix / "lib" / (name + ".py"))
            for name in ("numpy", "scipy", "torch", "gpytorch", "linear_operator")
        },
    }
    gate.validate_runtime(info, executable)
    info["executable"] = str(base)
    with pytest.raises(ValueError, match="runtime origin"):
        gate.validate_runtime(info, executable)


def test_factor_and_resource_mutations_rejected():
    frozen = {
        "selected_ranks": {"512": 1},
        "provenance": {"source_hashes": {}},
        "dataset": {},
    }
    x = case(512, 0, frozen)
    x["factor_metadata"].pop("precision_bits")
    with pytest.raises(ValueError, match="factor metadata"):
        gate.validate_case(x, 512, 0, frozen, P)
    x = case(512, 0, frozen)
    x["retained_components_bytes"]["source_values"] = -1
    with pytest.raises(ValueError, match="retained resource values"):
        gate.validate_case(x, 512, 0, frozen, P)


def test_controller_delegates_both_replays(tmp_path, monkeypatch):
    from contextlib import contextmanager
    from types import SimpleNamespace

    root, _ = fixture(tmp_path, monkeypatch)
    data = tmp_path / "dataset.zip"
    data.write_bytes(b"synthetic")
    monkeypatch.setattr(gate, "dataset_check", lambda p: gate.digest(p))
    monkeypatch.setattr(gate, "validate_payload", lambda *args: None)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=gate.ROOT, text=True
    ).strip()
    runtime = {
        "platform": "linux",
        "machine": "x86_64",
        "python": "3.12.9",
        "executable": "/env/bin/python",
        "prefix": "/env",
        "versions": {
            "numpy": "2.2.6",
            "scipy": "1.18.1",
            "torch": "2.14.1+cpu",
            "gpytorch": "1.15.2",
            "linear-operator": "0.6.1",
        },
        "threads": dict.fromkeys(gate.THREADS, "1"),
        "modules": {
            n: f"/env/lib/{n}.py"
            for n in ("numpy", "scipy", "torch", "gpytorch", "linear_operator")
        },
    }
    seen = []

    @contextmanager
    def snapshot(*args):
        yield gate.ROOT

    def runner(cmd, **kwargs):
        kind = cmd[6]
        seen.append(kind)
        assert cmd[1] == "-I" and kwargs["timeout"] <= gate.TIMEOUT
        assert "urlopen=no_network" in cmd[3]
        Path(cmd[4]).write_text(json.dumps(runtime))
        if kind == "science":
            output = json.loads((root / "full-audit.json").read_text())
            output["source_revision"] = revision
        else:
            paths = {
                "provider_json": root / "provider.json",
                "provider_npz": root / "provider.npz",
                "case_json": root / "full/n512-repeat0.json",
                "case_npz": root / "full/n512-repeat0.npz",
                "manifest": root / "full/campaign.json",
                "calibration": gate.ROOT / "data/covariance/calibration.json",
                "benchmark": gate.ROOT / "experiments/covariance_benchmark.py",
                "data": data,
            }
            output = {
                "schema": "kernelgecp.covariance-provider-replay.v1",
                "status": "passed",
                "scope": "synthetic fixture",
                "historical_reported_log_probability_error": 0.0,
                "historical_log_probability_value_available": False,
                "scipy_reference_log_probability": 1.0,
                "input_sha256": {k: gate.digest(v) for k, v in paths.items()},
                "actual_replay": {
                    "errors": dict.fromkeys(
                        (
                            "solution_max_abs_error",
                            "mean_max_abs_error",
                            "std_max_abs_error",
                            "log_probability_abs_error",
                        ),
                        0.0,
                    ),
                    "log_probability": 1.0,
                    "device": "cpu",
                    "dtype": "torch.float64",
                    "finite_counts": dict.fromkeys(("solution", "mean", "std"), 512),
                    "versions": {
                        **{
                            k: v
                            for k, v in runtime["versions"].items()
                            if k != "linear-operator"
                        },
                        "linear_operator": "0.6.1",
                    },
                    "module_paths": runtime["modules"],
                },
            }
        Path(cmd[-1]).write_text(json.dumps(output))
        return SimpleNamespace(returncode=0, stdout="fixture", stderr="")

    result = gate.audit(
        root,
        revision,
        "/env/bin/python",
        data,
        tmp_path / "out",
        runner=runner,
        snapshot=snapshot,
        identity_probe=lambda _: {"revision": "b" * 40},
    )
    assert result["status"] == "passed", result.get("error")
    assert seen == ["science", "provider"]
