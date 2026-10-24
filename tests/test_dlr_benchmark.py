"""Matched forward-observation, cached fitting and frozen-worker regressions."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from scipy.linalg import lstsq

import experiments.dlr_benchmark as benchmark
from experiments.dlr_adapters import AdapterUnavailableError
from experiments.dlr_benchmark import (
    BenchmarkCase,
    BenchmarkConfig,
    CachedFitPlan,
    fingerprint,
    load_config,
    run_benchmark,
    run_case,
    source_snapshot,
)
from experiments.dlr_workloads import (
    common_tau_meshes,
    green_function,
    make_pole_workloads,
)


def _config(**changes: Any) -> BenchmarkConfig:
    return BenchmarkConfig(
        cases=(BenchmarkCase(beta=2.0, wmax=2.0, tolerance=1e-5),),
        grid_order=6,
        repeats=1,
        warmups=0,
        **changes,
    )


def test_cached_fit_factors_once_and_matches_independent_complex_qr_lstsq(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    matrix = np.column_stack((np.ones(5), np.linspace(-2, 2, 5)))
    real_rhs = np.array([0.1, 0.4, 0.8, -0.3, 1.2])
    complex_rhs = np.column_stack((real_rhs + 2j, 3 * real_rhs - 0.7j))
    svd = np.linalg.svd
    factorizations = []

    def count_svd(design: np.ndarray, **kwargs: Any) -> Any:
        factorizations.append(design.copy())
        return svd(design, **kwargs)

    monkeypatch.setattr(benchmark.np.linalg, "svd", count_svd)
    plan = CachedFitPlan(matrix, rcond=1e-14)
    # The factorization owns its matrix instead of observing later caller writes.
    matrix[:] = 999
    for rhs in (real_rhs, complex_rhs, real_rhs):
        expected = lstsq(factorizations[0], rhs, cond=1e-14, lapack_driver="gelsy")[0]
        np.testing.assert_allclose(plan.fit(rhs), expected, rtol=5e-15, atol=5e-16)
    assert len(factorizations) == 1
    assert plan.diagnostics["effective_rank"] == 2
    assert plan.diagnostics["precision_bits"] == 53
    assert plan.diagnostics["coefficient_coordinates"].endswith("not_basis_invariant")


def test_cached_fit_reports_rank_truncation_and_zero_design() -> None:
    plan = CachedFitPlan(np.diag([1.0, 1e-16]), rcond=1e-14)
    np.testing.assert_array_equal(plan.fit(np.ones(2)), [1.0, 0.0])
    assert plan.diagnostics["effective_rank"] == 1
    zero = CachedFitPlan(np.zeros((3, 2)), rcond=1e-14)
    np.testing.assert_array_equal(zero.fit(np.ones(3)), [0.0, 0.0])
    assert zero.diagnostics["effective_rank"] == 0
    assert zero.diagnostics["condition_number"] is None


def test_condition_overflow_is_json_safe_and_nonfinite_fit_is_rejected() -> None:
    plan = CachedFitPlan(np.diag([1.0, 1e-320]), rcond=1e-323)
    assert plan.diagnostics["condition_number"] is None
    assert plan.diagnostics["condition_status"] == "singular_or_outside_binary64_range"
    json.dumps(plan.diagnostics, allow_nan=False)
    with (
        np.errstate(over="ignore", invalid="ignore"),
        pytest.raises(ValueError, match="nonfinite fitted coefficients"),
    ):
        plan.fit(np.ones(2))


@pytest.mark.parametrize(
    "matrix",
    [
        1.0,
        [1.0, 2.0],
        np.empty((0, 2)),
        np.empty((2, 0)),
        np.ones((2, 2, 1)),
        [[1.0, float("nan")]],
        [[1.0, float("inf")]],
        [[1.0 + 1j]],
    ],
)
def test_cached_fit_rejects_invalid_design(matrix: Any) -> None:
    with pytest.raises(ValueError):
        CachedFitPlan(matrix, rcond=1e-14)


@pytest.mark.parametrize("rcond", [0, 1, -1, float("nan"), float("inf"), True])
def test_cached_fit_rejects_invalid_cutoff(rcond: float) -> None:
    with pytest.raises(ValueError):
        CachedFitPlan(np.eye(2), rcond=rcond)


@pytest.mark.parametrize(
    "rhs",
    [
        1.0,
        np.ones(3),
        np.ones((2, 1, 1)),
        np.array([1, float("nan")]),
        np.array([1, float("inf")]),
        np.array([1, 1 + complex(0, float("inf"))]),
    ],
)
def test_cached_fit_rejects_invalid_observation_axis_or_values(rhs: Any) -> None:
    plan = CachedFitPlan(np.eye(2), rcond=1e-14)
    with pytest.raises(ValueError, match="RHS"):
        plan.fit(rhs)


@pytest.mark.parametrize(
    "field,value",
    [
        ("beta", 0),
        ("beta", True),
        ("beta", np.bool_(True)),
        ("beta", [1.0]),
        ("beta", np.array(1.0)),
        ("beta", float("nan")),
        ("wmax", -1),
        ("wmax", float("inf")),
        ("tolerance", 0),
        ("tolerance", 1),
        ("tolerance", True),
    ],
)
def test_case_rejects_invalid_physical_scalars(field: str, value: Any) -> None:
    with pytest.raises(ValueError):
        replace(BenchmarkCase(1, 1), **{field: value})


@pytest.mark.parametrize("beta,wmax", [(0.5, 1), (1e300, 1e300)])
def test_case_rejects_out_of_protocol_dimensionless_cutoff(
    beta: float, wmax: float
) -> None:
    with pytest.raises(ValueError, match="beta\\*wmax"):
        BenchmarkCase(beta, wmax)


@pytest.mark.parametrize(
    "field,value",
    [
        ("cases", ()),
        ("cases", ("invalid",)),
        ("methods", ()),
        ("methods", ("kernelgecp", "kernelgecp")),
        ("methods", ("unknown",)),
        ("grid_order", 2),
        ("grid_order", True),
        ("max_rank", 0),
        ("max_rank", 2.5),
        ("repeats", 0),
        ("warmups", -1),
        ("warmups", True),
        ("fit_rcond", 0),
        ("fit_rcond", 1),
        ("fit_rcond", float("nan")),
        ("timeout_seconds", 0),
        ("timeout_seconds", float("inf")),
        ("timeout_seconds", True),
        ("timeout_seconds", np.bool_(True)),
        ("timeout_seconds", [1.0]),
    ],
)
def test_config_rejects_invalid_budget_method_and_scalar(
    field: str, value: Any
) -> None:
    with pytest.raises(ValueError):
        replace(_config(), **{field: value})


def test_config_loading_rejects_missing_and_unknown_fields(tmp_path: Path) -> None:
    config = tmp_path / "benchmark.toml"
    for text in (
        "schema_version = 2\n",
        "schema_version = 1\n",
        "schema_version = 1\nunknown = 1\n[[cases]]\nbeta = 1\nwmax = 1\n",
    ):
        config.write_text(text)
        with pytest.raises(ValueError):
            load_config(config)


def test_fingerprint_is_shape_sensitive_and_normalizes_real_storage() -> None:
    values = np.arange(6, dtype=np.float64).reshape(2, 3)
    assert fingerprint(values) == fingerprint(np.asfortranarray(values))
    assert fingerprint(values) == fingerprint(values.astype(">f8"))
    assert fingerprint(values) != fingerprint(values.ravel())


def test_small_gecp_case_has_physical_accuracy_and_native_conditioning() -> None:
    config = _config(methods=("kernelgecp",), max_rank=32)
    case = config.cases[0]
    record = run_case("kernelgecp", case, config)
    assert record["status"] == "accuracy_passed"
    assert 1 <= record["rank"] <= config.max_rank
    assert record["basis_metadata"]["native_converged"]
    assert record["basis_metadata"]["certification_status"] == "not_certified"
    assert record["certification_status"] == "not_certified_sampled_validation"
    conditioning = record["fit_diagnostics"]
    assert conditioning["effective_rank"] == record["rank"]
    assert np.isfinite(conditioning["condition_number"])
    assert conditioning["condition_number"] >= 1
    for workload in record["workloads"]:
        assert np.isfinite(workload["training_absolute_error"])
        assert workload["training_absolute_error"] <= case.tolerance
        assert workload["absolute_error"] <= case.tolerance
        assert workload["accuracy_passed"]
    assert record["measurements"]["rss_scope"] == "in_process_diagnostic_not_isolated"
    json.dumps(record, allow_nan=False)


def test_matched_methods_share_only_observations_and_freeze_before_holdout_truth(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    case = config.cases[0]
    training, validation = common_tau_meshes(
        beta=case.beta, wmax=case.wmax, order=config.grid_order
    )
    workload_count = len(make_pole_workloads(beta=case.beta, wmax=case.wmax))
    events: list[str] = []
    method_design_times: dict[str, list[np.ndarray]] = {}

    class ObservationOnlyBasis:
        rank = 3
        tau_nodes = np.array([0.0, case.beta / 2, case.beta])
        metadata = {"test_basis": "polynomial_without_truth_supports"}

        def __init__(self, method: str) -> None:
            self.method = method

        def design_matrix(self, times: np.ndarray) -> np.ndarray:
            method_design_times[self.method].append(times.copy())
            normalized = times / case.beta
            # Native coordinate scaling varies without changing the span.
            scale = 1 + benchmark.METHODS.index(self.method)
            return np.column_stack((np.ones_like(times), normalized, normalized**2)) * (
                np.array([1.0, scale, scale**2])
            )

    def build(
        method: str, seen_case: BenchmarkCase, seen_config: BenchmarkConfig
    ) -> Any:
        # The construction boundary receives physical parameters and budgets,
        # with neither workload measures nor forward/holdout values.
        assert seen_case is case
        assert seen_config is config
        method_design_times[method] = []
        events.append("build")
        return ObservationOnlyBasis(method)

    direct_green = benchmark.green_function
    direct_fit = CachedFitPlan.fit

    def observe_truth(workload: Any, times: np.ndarray, **kwargs: Any) -> np.ndarray:
        if np.array_equal(times, validation):
            assert events.count("initial_fit") == workload_count
            events.append("holdout_truth")
        else:
            np.testing.assert_array_equal(times, training)
            assert "build" not in events
            events.append("training_truth")
        return direct_green(workload, times, **kwargs)

    def observe_fit(plan: CachedFitPlan, values: np.ndarray) -> np.ndarray:
        if "holdout_truth" not in events:
            events.append("initial_fit")
        return direct_fit(plan, values)

    monkeypatch.setattr(benchmark, "_build_basis", build)
    monkeypatch.setattr(benchmark, "_preload_provider", lambda method: None)
    monkeypatch.setattr(benchmark, "green_function", observe_truth)
    monkeypatch.setattr(CachedFitPlan, "fit", observe_fit)
    records = []
    for method in benchmark.METHODS:
        events.clear()
        record = run_case(method, case, config)
        assert record["status"] in ("accuracy_passed", "accuracy_failed")
        assert events[:workload_count] == ["training_truth"] * workload_count
        assert events.count("initial_fit") == workload_count
        assert events.count("holdout_truth") == workload_count
        records.append(record)
    assert all(record["inputs"] == records[0]["inputs"] for record in records)
    for times in method_design_times.values():
        np.testing.assert_array_equal(times[0], training)
        np.testing.assert_array_equal(times[1], validation)
    for index in range(workload_count):
        errors = [record["workloads"][index]["absolute_error"] for record in records]
        np.testing.assert_allclose(errors, errors[0], rtol=1e-12, atol=2e-15)


@pytest.mark.parametrize(
    "error,status",
    [
        (AdapterUnavailableError("optional dependency absent"), "import_unavailable"),
        (RuntimeError("native construction failed"), "runtime_failed"),
    ],
)
def test_optional_import_and_runtime_failures_are_retained(
    monkeypatch: pytest.MonkeyPatch, error: Exception, status: str
) -> None:
    def fail(*args: Any) -> Any:
        raise error

    monkeypatch.setattr(benchmark, "_build_basis", fail)
    monkeypatch.setattr(benchmark, "_preload_provider", lambda method: None)
    config = _config()
    result = run_case("pydlr", config.cases[0], config, isolated=True)
    assert result["status"] == status
    assert str(error) in result["error"]
    assert result["inputs"]["observation_hashes"]
    assert "workloads" not in result
    assert "fresh_worker_process" in result["measurements"]["rss_scope"]
    json.dumps(result, allow_nan=False)


def test_nonfinite_holdout_design_fails_without_success_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class BadBasis:
        def design_matrix(self, times: np.ndarray) -> np.ndarray:
            # Training contains the endpoints and holdout contains neither.
            values = np.ones((len(times), 1))
            if times[0] != 0:
                values[:] = np.nan
            return values

    monkeypatch.setattr(benchmark, "_build_basis", lambda *args: BadBasis())
    monkeypatch.setattr(benchmark, "_preload_provider", lambda method: None)
    config = _config()
    result = run_case("pydlr", config.cases[0], config)
    assert result["status"] == "runtime_failed"
    assert "nonfinite validation design" in result["error"]
    assert "workloads" not in result


def _minimal_source_tree(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "source"
    root.mkdir()
    (root / "kernelgecp").mkdir()
    (root / "experiments").mkdir()
    (root / "kernelgecp" / "example.py").write_text("value = 1\n")
    (root / "experiments" / "example.py").write_text("value = 2\n")
    (root / "pyproject.toml").write_text("[project]\nname='test'\n")
    (root / "uv.lock").write_text("version = 1\n")
    config = root / "benchmark.toml"
    config.write_text("schema_version = 1\n")
    return root, config


def _mock_git(monkeypatch: pytest.MonkeyPatch, *, dirty: str = "") -> None:
    def git(command: list[str], **kwargs: Any) -> Any:
        assert command[0] == "git"
        if command[1] == "status":
            return SimpleNamespace(stdout=dirty)
        if command[1] == "ls-files":
            return SimpleNamespace(stdout="")
        return SimpleNamespace(
            stdout="tree\n" if "HEAD^{tree}" in command else "head\n"
        )

    monkeypatch.setattr(benchmark.subprocess, "run", git)


def test_source_snapshot_hashes_sources_config_lock_and_revision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, config = _minimal_source_tree(tmp_path)
    _mock_git(monkeypatch)
    snapshot = source_snapshot(root, config)
    assert snapshot["git_commit"] == "head"
    assert snapshot["git_tree"] == "tree"
    assert snapshot["config_sha256"] == hashlib.sha256(config.read_bytes()).hexdigest()
    assert set(snapshot["source_sha256"]) == {
        "kernelgecp/example.py",
        "experiments/example.py",
        "pyproject.toml",
        "uv.lock",
    }
    for name, digest in snapshot["source_sha256"].items():
        assert digest == hashlib.sha256((root / name).read_bytes()).hexdigest()


def test_source_snapshot_rejects_dirty_tracked_revision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, config = _minimal_source_tree(tmp_path)
    _mock_git(monkeypatch, dirty=" M experiments/example.py\n")
    with pytest.raises(ValueError, match="clean tracked"):
        source_snapshot(root, config)


def test_source_snapshot_rejects_untracked_experiment_code(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, config = _minimal_source_tree(tmp_path)

    def git(command: list[str], **kwargs: Any) -> Any:
        return SimpleNamespace(
            stdout="experiments/new.py\n" if command[1] == "ls-files" else ""
        )

    monkeypatch.setattr(benchmark.subprocess, "run", git)
    with pytest.raises(ValueError, match="untracked"):
        source_snapshot(root, config)


def test_workers_retain_timeout_exit_and_missing_output_with_pinned_threads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(timeout_seconds=0.1)
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: {"head": "frozen"})
    seen = []

    def worker(command: list[str], **kwargs: Any) -> Any:
        method = command[command.index("--worker") + 1]
        seen.append(method)
        assert kwargs["timeout"] == 0.1
        assert kwargs["cwd"] == tmp_path
        assert all(kwargs["env"][name] == "1" for name in benchmark.THREAD_VARIABLES)
        if method == "kernelgecp":
            raise subprocess.TimeoutExpired(command, 0.1)
        if method == "pydlr":
            return SimpleNamespace(returncode=7, stderr="x" * 5000)
        return SimpleNamespace(returncode=0, stderr="no output was written")

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    result = run_benchmark(tmp_path / "config.toml", tmp_path)
    assert seen == list(benchmark.METHODS)
    assert [record["status"] for record in result["records"]] == [
        "timeout",
        "worker_failed",
        "worker_failed",
    ]
    assert len(result["records"][1]["stderr"]) == 4000
    assert result["claim_level"] == "sampled_forward_pilot_not_speedup_or_certification"


def test_source_change_during_workers_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(methods=("kernelgecp",))
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    snapshots = iter(({"head": "first"}, {"head": "first"}, {"head": "changed"}))
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: next(snapshots))
    monkeypatch.setattr(
        benchmark.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=1, stderr="failure"),
    )
    with pytest.raises(RuntimeError, match="changed during benchmark"):
        run_benchmark(tmp_path / "config.toml", tmp_path)


def test_invalid_worker_json_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(methods=("kernelgecp",))
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: {"head": "frozen"})

    def worker(command: list[str], **kwargs: Any) -> Any:
        output = Path(command[command.index("--output") + 1])
        output.write_text("{broken")
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    result = run_benchmark(tmp_path / "config.toml", tmp_path)
    assert result["records"][0]["status"] == "worker_failed"
    assert "invalid worker payload" in result["records"][0]["error"]


def _worker_payload(
    config: BenchmarkConfig, snapshot: dict[str, Any]
) -> dict[str, Any]:
    case = config.cases[0]
    training, validation = common_tau_meshes(
        beta=case.beta, wmax=case.wmax, order=config.grid_order
    )
    return {
        "method": "kernelgecp",
        "case": asdict(case),
        "protocol": "shared_observations_cached_svd_v1",
        "status": "import_unavailable",
        "grid_order": config.grid_order,
        "worker_provenance": snapshot,
        "inputs": {
            "training_hash": fingerprint(training),
            "validation_hash": fingerprint(validation),
            "observation_hashes": [
                fingerprint(green_function(w, training, beta=case.beta, wmax=case.wmax))
                for w in make_pole_workloads(beta=case.beta, wmax=case.wmax)
            ],
        },
    }


def test_valid_frozen_worker_payload_is_accepted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(methods=("kernelgecp",))
    snapshot = {"head": "frozen"}
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: snapshot)
    # Use the real run_case output so any missing required worker metadata is
    # detected at the parent/worker integration boundary.
    payload = run_case("kernelgecp", config.cases[0], config, isolated=True)
    payload["worker_provenance"] = snapshot

    def worker(command: list[str], **kwargs: Any) -> Any:
        output = Path(command[command.index("--output") + 1])
        assert command[command.index("--expected-provenance") + 1] == (
            benchmark.provenance_digest(snapshot)
        )
        output.write_text(json.dumps(payload, allow_nan=False))
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    result = run_benchmark(tmp_path / "config.toml", tmp_path)
    assert result["records"] == [payload]


@pytest.mark.parametrize(
    "field,value",
    [
        ("method", "pydlr"),
        ("case", {"beta": 99, "wmax": 2, "tolerance": 1e-5}),
        ("protocol", "wrong"),
        ("status", "certified"),
        ("grid_order", 3),
        ("inputs", {"training_hash": "wrong", "validation_hash": "wrong"}),
        ("nonfinite_field", float("nan")),
    ],
)
def test_worker_identity_protocol_and_requested_mesh_cannot_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str, value: Any
) -> None:
    config = _config(methods=("kernelgecp",))
    snapshot = {"head": "frozen"}
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: snapshot)
    payload = _worker_payload(config, snapshot)
    payload[field] = value
    if field == "grid_order":
        # Even self-consistent hashes at an unrequested resolution are invalid.
        training, validation = common_tau_meshes(
            beta=config.cases[0].beta, wmax=config.cases[0].wmax, order=value
        )
        payload["inputs"] = {
            "training_hash": fingerprint(training),
            "validation_hash": fingerprint(validation),
        }

    def worker(command: list[str], **kwargs: Any) -> Any:
        Path(command[command.index("--output") + 1]).write_text(json.dumps(payload))
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    result = run_benchmark(tmp_path / "config.toml", tmp_path)
    assert result["records"][0]["status"] == "worker_failed"


def test_worker_provenance_mismatch_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(methods=("kernelgecp",))
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: {"head": "frozen"})
    payload = _worker_payload(config, {"head": "other"})

    def worker(command: list[str], **kwargs: Any) -> Any:
        Path(command[command.index("--output") + 1]).write_text(json.dumps(payload))
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    with pytest.raises(RuntimeError, match="does not match frozen parent"):
        run_benchmark(tmp_path / "config.toml", tmp_path)


def test_source_change_while_loading_config_fails_before_launch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    snapshots = iter(({"head": "first"}, {"head": "changed"}))
    monkeypatch.setattr(benchmark, "load_config", lambda path: _config())
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: next(snapshots))

    def unexpected_worker(*args: Any, **kwargs: Any) -> Any:
        pytest.fail("workers must not launch after configuration drift")

    monkeypatch.setattr(benchmark.subprocess, "run", unexpected_worker)
    with pytest.raises(RuntimeError, match="changed while loading configuration"):
        run_benchmark(tmp_path / "config.toml", tmp_path)


def test_main_never_replaces_existing_output_with_nonfinite_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "result.json"
    output.write_text('{"old":true}\n')
    monkeypatch.setattr(benchmark, "load_config", lambda path: _config())
    monkeypatch.setattr(
        benchmark, "run_benchmark", lambda *args: {"invalid": float("nan")}
    )
    monkeypatch.setattr(sys, "argv", ["dlr_benchmark", "--output", str(output)])
    with pytest.raises(ValueError, match="Out of range float"):
        benchmark.main()
    assert output.read_text() == '{"old":true}\n'
    assert not output.with_suffix(".json.tmp").exists()


@pytest.mark.parametrize("changed", [False, True])
def test_worker_main_stamps_frozen_provenance_and_blocks_changed_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, changed: bool
) -> None:
    config = _config(methods=("kernelgecp",))
    output = tmp_path / "worker.json"
    output.write_text('{"old":true}\n')
    before = {"head": "frozen"}
    after = {"head": "changed"} if changed else before
    snapshots = iter((before, after))
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: next(snapshots))
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    calls = []

    def run(
        method: str, case: BenchmarkCase, seen: BenchmarkConfig, **kwargs: Any
    ) -> Any:
        calls.append((method, case, seen, kwargs))
        return {"status": "import_unavailable"}

    monkeypatch.setattr(benchmark, "run_case", run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "dlr_benchmark",
            "--worker",
            "kernelgecp",
            "--output",
            str(output),
            "--expected-provenance",
            benchmark.provenance_digest(before),
        ],
    )
    if changed:
        with pytest.raises(RuntimeError, match="changed during execution"):
            benchmark.main()
        assert output.read_text() == '{"old":true}\n'
    else:
        benchmark.main()
        assert json.loads(output.read_text())["worker_provenance"] == before
    assert calls == [("kernelgecp", config.cases[0], config, {"isolated": True})]
