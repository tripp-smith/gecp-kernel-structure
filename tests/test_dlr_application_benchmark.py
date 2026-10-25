"""Isolation, matched calibration and fail-closed physical benchmark checks."""

from __future__ import annotations

import copy
import json
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from scipy.linalg import lstsq

import experiments.dlr_application_benchmark as benchmark
from experiments.dlr_adapters import AdapterUnavailableError
from experiments.dlr_application_benchmark import (
    ApplicationCase,
    ApplicationConfig,
    array_payload,
    fingerprint,
    run_calibration,
    run_case,
)
from experiments.dlr_workloads import common_tau_meshes
from kernelgecp import FermionicBasis, LeastSquaresPlan


@pytest.fixture(scope="module", autouse=True)
def synthetic_test_meshes() -> Any:
    """Exercise the harness without inspecting the frozen canonical holdouts."""
    patch = pytest.MonkeyPatch()

    def meshes(case: ApplicationCase, config: ApplicationConfig) -> tuple[Any, ...]:
        training, calibration = common_tau_meshes(
            beta=case.beta, wmax=case.wmax, order=config.grid_order
        )
        width = np.diff(training)
        final = np.sort(
            np.concatenate((training[:-1] + width / 3, training[:-1] + 2 * width / 3))
        )
        return (
            training,
            calibration,
            final,
            np.array([-64, -32, -16, -8, -4, -1, 0, 3, 7, 15, 31, 63]),
            np.array(
                [
                    -513,
                    -257,
                    -193,
                    -129,
                    -12,
                    -6,
                    -3,
                    -2,
                    1,
                    2,
                    5,
                    11,
                    128,
                    192,
                    256,
                    512,
                ]
            ),
        )

    patch.setattr(benchmark, "_meshes", meshes)
    patch.setattr(benchmark, "solver_index_mesh", lambda: np.arange(-16, 16))
    patch.setattr(
        benchmark, "solver_validation_index_mesh", lambda: np.array([-129, 128])
    )
    yield
    patch.undo()


def _config(**kwargs: Any) -> ApplicationConfig:
    return ApplicationConfig(
        cases=(ApplicationCase(2.0, 2.0, target=1e-5),),
        methods=("kernelgecp",),
        native_eps_candidates=(1e-6, 1e-8),
        grid_order=6,
        max_rank=32,
        repeats=1,
        warmups=0,
        **kwargs,
    )


@pytest.mark.parametrize("version", ["true", "1.0", "2"])
def test_config_schema_is_an_exact_supported_integer(tmp_path, version):
    path = tmp_path / "bad.toml"
    path.write_text(f"schema_version = {version}\n")
    with pytest.raises(ValueError, match="schema"):
        benchmark.load_config(path)


@pytest.fixture(scope="module")
def calibrated() -> tuple[ApplicationConfig, dict[str, Any], dict[str, Any]]:
    config = _config()
    calibration = run_calibration("kernelgecp", config.cases[0], config)
    measurement = run_case(
        "kernelgecp", config.cases[0], config, calibration=calibration
    )
    return config, calibration, measurement


def test_finite_complex_payload_hash_preserves_shape_dtype_and_imaginary_part() -> None:
    values = np.array([[1 + 2j, 3 - 4j]], dtype=np.complex128)
    assert fingerprint(values) == fingerprint(np.asfortranarray(values))
    assert fingerprint(values) == fingerprint(values.astype(">c16"))
    assert fingerprint(values) != fingerprint(values.conj())
    assert fingerprint(values) != fingerprint(values.ravel())
    assert fingerprint(values) != fingerprint(values.real)
    assert fingerprint(values) != fingerprint(values.astype(np.complex64))
    np.testing.assert_array_equal(
        benchmark._decode_array(array_payload(values)), values
    )
    corrupt = array_payload(values)
    corrupt["imag"][0][0] = 42
    with pytest.raises(ValueError, match="hash/content"):
        benchmark._decode_array(corrupt)
    for bad in ([np.inf], [complex(1, np.nan)], [True], ["a"]):
        with pytest.raises(ValueError, match="finite numeric"):
            fingerprint(bad)


@pytest.mark.parametrize(
    "field,value",
    [
        ("beta", True),
        ("wmax", 0),
        ("target", np.nan),
        ("target", 1.0),
        ("solver", 1),
    ],
)
def test_case_rejects_invalid_inputs(field: str, value: Any) -> None:
    with pytest.raises(ValueError):
        replace(ApplicationCase(1.0, 1.0), **{field: value})


@pytest.mark.parametrize(
    "field,value",
    [
        ("cases", ()),
        ("methods", ("kernelgecp", "kernelgecp")),
        ("methods", ("unknown",)),
        ("native_eps_candidates", ()),
        ("native_eps_candidates", (1e-8, 1e-6)),
        ("native_eps_candidates", (True,)),
        ("grid_order", True),
        ("max_rank", 0),
        ("precision_bits", 53),
        ("repeats", 0),
        ("warmups", -1),
        ("fit_rcond", 1),
        ("timeout_seconds", np.inf),
        ("solver_mixing", 0),
        ("solver_max_iterations", 0),
    ],
)
def test_config_rejects_invalid_inputs(field: str, value: Any) -> None:
    with pytest.raises(ValueError):
        replace(_config(), **{field: value})


def test_canonical_config_declares_four_cases_and_fixed_controls() -> None:
    config = benchmark.load_config(Path("experiments/dlr_application_benchmark.toml"))
    assert [(case.beta, case.wmax, case.solver) for case in config.cases] == [
        (10, 1, False),
        (10, 10, False),
        (20, 50, False),
        (10, 2, True),
    ]
    assert config.native_eps_candidates == (1e-6, 1e-8, 1e-10, 1e-12)
    assert config.solver_mixing == 0.5
    assert config.solver_tolerance == 1e-9
    assert config.solver_max_iterations == 200


def test_small_public_gecp_pass_and_numerical_fields_repeat(
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, measurement = calibrated
    assert calibration["status"] == "calibration_passed"
    assert measurement["status"] == "accuracy_passed"
    assert measurement["basis_metadata"]["api"] == "public_FermionicBasis"
    assert measurement["fit_diagnostics"]["factorization"] == "cached_thin_svd"
    assert measurement["fit_diagnostics"]["effective_rank"] == measurement["rank"]
    assert measurement["fit_diagnostics"]["singular_values"]
    assert "not_certified" in measurement["certification_status"]
    assert "in_process" in measurement["measurements"]["rss_scope"]
    repeat = run_case("kernelgecp", config.cases[0], config, calibration=calibration)
    assert {k: v for k, v in measurement.items() if k != "measurements"} == {
        k: v for k, v in repeat.items() if k != "measurements"
    }
    json.dumps(measurement, allow_nan=False)


def test_complex_svd_fit_has_independent_qr_oracle() -> None:
    x = np.linspace(-1, 1, 9)
    design = np.column_stack((np.ones(9), x + 0.2j, x**2 - 0.3j))
    values = np.stack((x + 1j, 2 * x - 0.7j, x**3 + 0.1j, -x), axis=1).reshape(9, 2, 2)
    plan = LeastSquaresPlan(design, rcond=1e-14)
    expected = lstsq(design, values.reshape(9, -1), cond=1e-14, lapack_driver="gelsy")[
        0
    ]
    np.testing.assert_allclose(plan.fit(values).reshape(3, -1), expected, atol=2e-15)
    assert plan.diagnostics.design_dtype == "complex128"


def test_calibration_retains_failed_candidates_and_never_reads_final_truth(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    case = config.cases[0]
    _, _, final, _, final_n = benchmark._meshes(case, config)
    original_time, original_frequency = (
        benchmark.green_function,
        benchmark.matsubara_green_function,
    )
    builds = []

    def time_truth(workload: Any, tau: Any, **kwargs: Any) -> Any:
        assert not np.array_equal(tau, final)
        return original_time(workload, tau, **kwargs)

    def frequency_truth(workload: Any, indices: Any, **kwargs: Any) -> Any:
        assert not np.array_equal(indices, final_n)
        return original_frequency(workload, indices, **kwargs)

    class InaccurateBasis:
        rank = 1
        tau_nodes = np.array([0.0])
        metadata = {"fixture": "inaccurate"}

        def design_tau(self, tau: np.ndarray) -> np.ndarray:
            return np.ones((len(tau), 1))

        def design_matsubara(self, indices: np.ndarray) -> np.ndarray:
            return np.zeros((len(indices), 1), dtype=complex)

    def construct(
        method: str, case: ApplicationCase, config: ApplicationConfig, eps: float
    ) -> Any:
        builds.append(eps)
        return (
            InaccurateBasis()
            if eps == 1e-6
            else FermionicBasis(
                case.beta,
                case.wmax,
                eps,
                grid_order=config.grid_order,
                max_rank=config.max_rank,
            )
        )

    monkeypatch.setattr(benchmark, "green_function", time_truth)
    monkeypatch.setattr(benchmark, "matsubara_green_function", frequency_truth)
    monkeypatch.setattr(benchmark, "_build_basis", construct)
    calibration = run_calibration("kernelgecp", case, config)
    assert builds == [1e-6, 1e-8]
    assert [candidate["status"] for candidate in calibration["candidates"]] == [
        "accuracy_failed",
        "accuracy_passed",
    ]
    assert benchmark._selection(calibration, config) == 1e-8
    # Same wrappers remain in place, so validation also cannot read final truth.
    assert benchmark._selection(calibration, config) == 1e-8


def test_import_and_arithmetic_failures_are_retained(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    monkeypatch.setattr(benchmark, "_preload_provider", lambda method: None)

    def fail(*args: Any) -> Any:
        raise ArithmeticError("intentional arithmetic failure")

    monkeypatch.setattr(benchmark, "_build_basis", fail)
    result = run_calibration("kernelgecp", config.cases[0], config)
    assert result["status"] == "calibration_failed"
    assert [c["status"] for c in result["candidates"]] == ["runtime_failed"] * 2
    assert benchmark._selection(result, config) is None

    def unavailable(method: str) -> Any:
        raise AdapterUnavailableError("not installed")

    monkeypatch.setattr(benchmark, "_preload_provider", unavailable)
    result = run_calibration("pydlr", config.cases[0], config)
    assert result["status"] == "import_unavailable"
    assert result["selected_native_eps"] is None


def test_worker_validation_rejects_identity_selection_and_physical_forgery(
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, measurement = calibrated
    snapshot = {"frozen": "revision"}
    calibration = {**calibration, "worker_provenance": snapshot}
    measurement = {
        **measurement,
        "worker_provenance": snapshot,
        "selection_digest": benchmark.selection_digest(calibration),
    }
    case = config.cases[0]
    benchmark._validated_worker(
        calibration, "kernelgecp", case, config, snapshot, "calibration"
    )
    benchmark._validated_worker(
        measurement, "kernelgecp", case, config, snapshot, "measurement", calibration
    )
    for stage, record, mutate in (
        ("calibration", calibration, lambda r: r.update(method="pydlr")),
        ("calibration", calibration, lambda r: r.update(status="unrecognized")),
        ("calibration", calibration, lambda r: r.update(selected_native_eps=1e-12)),
        (
            "calibration",
            calibration,
            lambda r: r["candidates"][0]["workloads"][0]["time"].update(
                absolute_error=0.0
            ),
        ),
        (
            "measurement",
            measurement,
            lambda r: r["inputs"]["matrix_amplitude"]["imag"][0].__setitem__(1, 0.0),
        ),
        ("measurement", measurement, lambda r: r.update(selected_native_eps=1e-12)),
        ("measurement", measurement, lambda r: r.update(selection_digest="forged")),
        (
            "measurement",
            measurement,
            lambda r: r["workloads"][0]["matsubara"]["truth"]["real"].__setitem__(
                0, [[0, 0], [0, 0]]
            ),
        ),
    ):
        broken = copy.deepcopy(record)
        mutate(broken)
        with pytest.raises((ValueError, RuntimeError)):
            benchmark._validated_worker(
                broken, "kernelgecp", case, config, snapshot, stage, calibration
            )
    with pytest.raises(RuntimeError, match="frozen parent"):
        benchmark._validated_worker(
            calibration, "kernelgecp", case, config, {}, "calibration"
        )


def test_runner_uses_two_workers_and_excludes_calibration_rss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, measurement = calibrated
    snapshot = {"frozen": "revision"}
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: snapshot)
    monkeypatch.setattr(benchmark, "load_config", lambda path: config)
    monkeypatch.setattr(benchmark.platform, "platform", lambda: "test-platform")
    calls = []

    def subprocess_run(command: list[str], **kwargs: Any) -> Any:
        calls.append((command, kwargs))
        stage = command[command.index("--stage") + 1]
        result = copy.deepcopy(calibration if stage == "calibration" else measurement)
        result["worker_provenance"] = snapshot
        if stage == "measurement":
            saved = json.loads(
                Path(command[command.index("--selection") + 1]).read_text()
            )
            result["selection_digest"] = saved["selection_digest"]
            assert "inputs" not in saved
            assert all(
                "workloads" not in candidate for candidate in saved["candidates"]
            )
            benchmark._validate_receipt(saved, "kernelgecp", config.cases[0], config)
            result["measurements"]["peak_process_rss_bytes"] = 123
        else:
            result["measurements"]["peak_process_rss_bytes"] = 456
        Path(command[command.index("--output") + 1]).write_text(json.dumps(result))
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(benchmark.subprocess, "run", subprocess_run)
    result = benchmark.run_benchmark(tmp_path / "config.toml", tmp_path)
    assert len(calls) == 2
    assert calls[0][0] != calls[1][0]
    assert all(
        kwargs["env"][name] == "1"
        for _, kwargs in calls
        for name in benchmark.THREAD_VARIABLES
    )
    record = result["records"][0]
    assert record["status"] == "accuracy_passed"
    assert record["calibration"]["measurements"]["peak_process_rss_bytes"] == 456
    assert record["measurement"]["measurements"]["peak_process_rss_bytes"] == 123


@pytest.mark.parametrize("failure", ["timeout", "returncode", "malformed"])
def test_runner_retains_worker_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    monkeypatch.setattr(benchmark, "source_snapshot", lambda *args: {})
    monkeypatch.setattr(benchmark, "load_config", lambda path: _config())
    monkeypatch.setattr(benchmark.platform, "platform", lambda: "test-platform")

    def subprocess_run(command: list[str], **kwargs: Any) -> Any:
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, 600)
        if failure == "malformed":
            Path(command[command.index("--output") + 1]).write_text("{truncated")
        return SimpleNamespace(
            returncode=1 if failure == "returncode" else 0, stderr="failed"
        )

    monkeypatch.setattr(benchmark.subprocess, "run", subprocess_run)
    result = benchmark.run_benchmark(tmp_path / "config.toml", tmp_path)
    assert result["records"][0]["status"] == (
        "timeout" if failure == "timeout" else "worker_failed"
    )
    assert "measurement" not in result["records"][0]


def test_solver_nonconvergence_cannot_report_overall_success() -> None:
    config = replace(
        _config(),
        cases=(ApplicationCase(2.0, 2.0, target=1e-5, solver=True),),
        solver_max_iterations=1,
    )
    case = config.cases[0]
    calibration = run_calibration("kernelgecp", case, config)
    record = run_case("kernelgecp", case, config, calibration=calibration)
    assert record["forward_status"] == "accuracy_passed"
    assert record["status"] == "accuracy_failed"
    assert record["overall_success"] is False
    assert record["solver"]["converged"] is False
    assert record["solver"]["status"] == "not_converged"
    assert record["solver"]["history"]
    record["worker_provenance"] = {}
    benchmark._validated_worker(
        record, "kernelgecp", case, config, {}, "measurement", calibration
    )
    record["overall_success"] = True
    with pytest.raises(ValueError, match="overall application"):
        benchmark._validated_worker(
            record, "kernelgecp", case, config, {}, "measurement", calibration
        )


def test_selected_measurement_uses_compact_receipt_without_calibration_arrays(
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config, calibration, measurement = calibrated
    receipt = benchmark._selection_receipt(
        {**calibration, "worker_provenance": {"source": "frozen"}}
    )
    assert "inputs" not in receipt
    assert "measurements" not in receipt
    assert all("workloads" not in candidate for candidate in receipt["candidates"])
    assert (
        benchmark._validate_receipt(receipt, "kernelgecp", config.cases[0], config)
        == calibration["selected_native_eps"]
    )

    def forbidden(*args: Any) -> Any:
        raise AssertionError("selected worker must not load/run calibration")

    monkeypatch.setattr(benchmark, "run_calibration", forbidden)
    monkeypatch.setattr(benchmark, "_selection", forbidden)
    selected = run_case("kernelgecp", config.cases[0], config, selection=receipt)
    assert selected["overall_success"]
    assert selected["workloads"] == measurement["workloads"]
    for mutate in (
        lambda r: r.update(selected_native_eps=1e-12),
        lambda r: r.update(selection_digest="invalid"),
        lambda r: r["candidates"][0]["max_absolute_errors"].update(time=1.0),
        lambda r: r["candidates"].append(copy.deepcopy(r["candidates"][0])),
    ):
        broken = copy.deepcopy(receipt)
        mutate(broken)
        with pytest.raises(ValueError):
            benchmark._validate_receipt(broken, "kernelgecp", config.cases[0], config)


@pytest.fixture(scope="module")
def converged_solver_record() -> tuple[
    ApplicationConfig, dict[str, Any], dict[str, Any]
]:
    config = replace(
        _config(), cases=(ApplicationCase(2.0, 2.0, target=1e-5, solver=True),)
    )
    calibration = run_calibration("kernelgecp", config.cases[0], config)
    record = run_case("kernelgecp", config.cases[0], config, calibration=calibration)
    assert record["solver"]["converged"]
    record["worker_provenance"] = {}
    benchmark._validated_worker(
        record, "kernelgecp", config.cases[0], config, {}, "measurement", calibration
    )
    return config, calibration, record


@pytest.mark.parametrize(
    "tamper",
    [
        "converged",
        "defect",
        "history_missing",
        "history_order",
        "history_continues",
        "negative_norm",
        "boolean_norm",
        "stop_reason",
        "tolerance",
        "mixing",
        "hopping",
        "beta",
        "max_iterations",
        "iterations",
        "history_metadata",
        "fit_cutoff",
        "values",
        "arithmetic",
    ],
)
def test_worker_rejects_forged_solver_convergence_evidence(
    converged_solver_record: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
    tamper: str,
) -> None:
    config, calibration, original = converged_solver_record
    record = copy.deepcopy(original)
    solver = record["solver"]
    if tamper == "converged":
        solver["converged"] = False
    elif tamper == "defect":
        solver["history"][-1]["defect_norm"] = 0.0
    elif tamper == "history_missing":
        solver["history"] = []
    elif tamper == "history_order":
        solver["history"][-1]["iteration"] = 1
    elif tamper == "history_continues":
        solver["history"][0]["defect_norm"] = 0.0
    elif tamper == "negative_norm":
        solver["history"][0]["step_norm"] = -1.0
    elif tamper == "boolean_norm":
        solver["history"][0]["projection_residual_norm"] = False
    elif tamper == "stop_reason":
        solver["stop_reason"] = "iteration_budget_exhausted"
    elif tamper in {"tolerance", "mixing", "hopping", "beta"}:
        solver["metadata"][tamper] *= 0.5
    elif tamper == "max_iterations" or tamper == "iterations":
        solver["metadata"][tamper] += 1
    elif tamper == "fit_cutoff":
        solver["fit_diagnostics"]["rcond"] = 1e-4
        solver["metadata"]["fit_diagnostics"] = copy.deepcopy(solver["fit_diagnostics"])
    elif tamper == "values":
        solver["frequency_values"] = array_payload(np.zeros(32, dtype=complex))
    elif tamper == "arithmetic":
        solver["metadata"]["arithmetic"] = "unknown"
    # Duplicate metadata cannot legitimize forged convergence/history evidence.
    solver["metadata"]["history"] = copy.deepcopy(solver["history"])
    solver["metadata"]["converged"] = solver["converged"]
    solver["metadata"]["stop_reason"] = solver["stop_reason"]
    if tamper != "iterations":
        solver["metadata"]["iterations"] = len(solver["history"])
    if tamper == "history_metadata":
        solver["metadata"]["history"][0]["step_norm"] += 1.0
    with pytest.raises(ValueError, match="solver"):
        benchmark._validated_worker(
            record,
            "kernelgecp",
            config.cases[0],
            config,
            {},
            "measurement",
            calibration,
        )


def test_actual_defect_threshold_is_not_relaxed_by_roundoff_comparison_allowance(
    converged_solver_record: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, original = converged_solver_record
    record = copy.deepcopy(original)
    solver = record["solver"]
    values = benchmark._decode_array(solver["frequency_values"])
    indices = benchmark._decode_array(solver["indices"])
    case = config.cases[0]
    z = 1j * (2.0 * indices + 1.0) * (np.pi / case.beta)

    def defect(perturbation: float) -> float:
        candidate = values + perturbation
        return float(
            np.max(np.abs(candidate - 1.0 / (z - (case.wmax / 2) ** 2 * candidate)))
        )

    low, high = 0.0, 1e-8
    boundary = config.solver_tolerance + 1e-15
    for _ in range(60):
        middle = (low + high) / 2
        if defect(middle) < boundary:
            low = middle
        else:
            high = middle
    actual = defect(high)
    assert config.solver_tolerance < actual < config.solver_tolerance + 2e-15
    solver["frequency_values"] = array_payload(values + high)
    solver["history"][-1]["defect_norm"] = config.solver_tolerance
    solver["metadata"]["history"] = copy.deepcopy(solver["history"])
    with pytest.raises(ValueError, match="actual defect"):
        benchmark._validated_worker(
            record, "kernelgecp", case, config, {}, "measurement", calibration
        )


def test_solver_rejects_history_beyond_configured_iteration_budget(
    converged_solver_record: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, original = converged_solver_record
    record = copy.deepcopy(original)
    iteration = record["solver"]["history"][0]
    record["solver"]["history"] = [
        {**iteration, "iteration": index}
        for index in range(1, config.solver_max_iterations + 2)
    ]
    with pytest.raises(ValueError, match="invalid length"):
        benchmark._validated_worker(
            record,
            "kernelgecp",
            config.cases[0],
            config,
            {},
            "measurement",
            calibration,
        )


def test_solver_budget_exhaustion_requires_complete_history() -> None:
    config = replace(
        _config(),
        cases=(ApplicationCase(2.0, 2.0, target=1e-5, solver=True),),
        solver_max_iterations=1,
    )
    calibration = run_calibration("kernelgecp", config.cases[0], config)
    record = run_case("kernelgecp", config.cases[0], config, calibration=calibration)
    record["worker_provenance"] = {}
    record["solver"]["history"] = []
    record["solver"]["metadata"].update(history=[], iterations=0)
    with pytest.raises(ValueError, match="incomplete history"):
        benchmark._validated_worker(
            record,
            "kernelgecp",
            config.cases[0],
            config,
            {},
            "measurement",
            calibration,
        )


@pytest.mark.parametrize(
    "fraction", [0.0, -0.1, 1.01, np.inf, np.nan, True, np.bool_(True), [0.1]]
)
def test_calibration_fraction_rejects_invalid_bounds_and_types(fraction: Any) -> None:
    with pytest.raises(ValueError, match="calibration_fraction"):
        replace(_config(), calibration_fraction=fraction)


def test_calibration_fraction_allows_one_and_rejects_target_product_underflow() -> None:
    assert replace(_config(), calibration_fraction=1.0).calibration_fraction == 1.0
    tiny_case = ApplicationCase(1.0, 1.0, target=float(np.nextafter(0.0, 1.0)))
    with pytest.raises(ValueError, match="calibration_absolute_target"):
        replace(_config(), cases=(tiny_case,), calibration_fraction=0.1)
    assert replace(_config(), cases=(tiny_case,), calibration_fraction=1.0)


def test_calibration_margin_is_explicit_and_final_target_stays_unchanged(
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, measurement = calibrated
    case = config.cases[0]
    assert case.target == 1e-5
    assert calibration["calibration_absolute_target"] == case.target * 0.1
    assert measurement["calibration_absolute_target"] == case.target * 0.1
    assert (
        benchmark.load_config(
            Path("experiments/dlr_application_benchmark.toml")
        ).calibration_fraction
        == 0.1
    )
    assert all(
        canonical_case.target == 1e-8
        for canonical_case in benchmark.load_config(
            Path("experiments/dlr_application_benchmark.toml")
        ).cases
    )
    # A residual between the two targets must pass final validation but fail
    # calibration, even if a worker supplies the looser final-target flag.
    original = calibration["candidates"][0]["workloads"][0]["time"]
    truth = benchmark._decode_array(original["truth"])
    falsely_passing = benchmark._error(truth + case.target / 2, truth, case.target)
    assert falsely_passing["accuracy_passed"]
    assert (
        falsely_passing["absolute_error"] > calibration["calibration_absolute_target"]
    )
    bad_calibration = copy.deepcopy(calibration)
    bad_calibration["candidates"][0]["workloads"][0]["time"] = falsely_passing
    with pytest.raises(ValueError, match="physical error"):
        benchmark._selection(bad_calibration, config)
    final_records = copy.deepcopy(measurement["workloads"])
    final_truth = benchmark._decode_array(final_records[0]["time"]["truth"])
    final_records[0]["time"] = benchmark._error(
        final_truth + case.target / 2, final_truth, case.target
    )
    benchmark._validate_workloads(final_records, case, config, "measurement")


def test_selection_receipt_cannot_replace_strict_calibration_target_with_final_target(
    calibrated: tuple[ApplicationConfig, dict[str, Any], dict[str, Any]],
) -> None:
    config, calibration, _ = calibrated
    case = config.cases[0]
    receipt = benchmark._selection_receipt({**calibration, "worker_provenance": {}})
    assert (
        receipt["calibration_absolute_target"]
        == case.target * config.calibration_fraction
    )
    inflated = copy.deepcopy(receipt)
    inflated["calibration_absolute_target"] = case.target
    with pytest.raises(ValueError, match="identity"):
        benchmark._validate_receipt(inflated, "kernelgecp", case, config)
    forged_pass = copy.deepcopy(receipt)
    forged_pass["candidates"][0]["max_absolute_errors"]["matsubara"] = case.target / 2
    with pytest.raises(ValueError, match="accuracy"):
        benchmark._validate_receipt(forged_pass, "kernelgecp", case, config)
    inflated_calibration = copy.deepcopy(calibration)
    inflated_calibration["calibration_absolute_target"] = case.target
    with pytest.raises(ValueError, match="absolute target"):
        benchmark._selection(inflated_calibration, config)
