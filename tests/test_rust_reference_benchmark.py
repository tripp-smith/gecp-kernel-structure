"""Fresh-worker Rust experiment contracts and fail-closed evidence checks."""

from __future__ import annotations

import copy
import subprocess
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from experiments import rust_reference_benchmark as bench


@pytest.fixture(scope="module")
def config() -> bench.RustBenchmarkConfig:
    return bench.RustBenchmarkConfig(
        beta=3.7,
        wmax=2.4,
        tolerance=1e-10,
        construction_precision="reference",
        grid_order=6,
        max_rank=32,
        cache_size=2,
        trials=3,
        query_count=3,
        rows_per_query=8,
        small_query_count=1,
        warm_hits=2,
    )


@pytest.fixture(scope="module")
def source() -> dict[str, Any]:
    return {"git_commit": "f" * 40, "native_binary": None}


@pytest.fixture(scope="module")
def baseline(
    config: bench.RustBenchmarkConfig, source: dict[str, Any]
) -> dict[str, Any]:
    # Only source-freeze plumbing is mocked. A real small public reference fit,
    # actual cache misses/hits and independent 512-bit QR checks all execute.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(bench, "frozen_snapshot", lambda *args: copy.deepcopy(source))
        for name in bench.THREAD_VARIABLES:
            patch.setenv(name, "1")
        record = bench.run_worker(
            {
                "config": asdict(config),
                "provider": "python",
                "trial": 0,
                "source": source,
                "root": ".",
                "config_path": "unused",
            }
        )
    assert record["status"] == "completed", record.get("error")
    return record


def test_real_small_public_worker_has_actual_misses_and_independent_checks(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
) -> None:
    validated = bench.validate_worker(baseline, config, "python", 0, source)
    assert validated["calls"]["stream"] == config.query_count
    assert validated["calls"]["warm"] == 0
    assert validated["oracle"]["time"]["max_error"] <= config.parity_absolute_target
    assert (
        validated["oracle"]["frequency"]["physical_max_error"]
        <= config.physical_absolute_target
    )
    assert (
        validated["basis_diagnostics"]["fit_backend"]
        == "numpy_svd_float64_or_complex128"
    )
    assert validated["measurements"]["native_setup_including_load_seconds"] == 0


def test_complete_complex_decoder_keeps_signed_zero_bits() -> None:
    values = np.zeros((2, 2), dtype=np.complex128)
    values.real[0, 0] = -0.0
    values.imag[1, 1] = -0.0
    assert (
        bench._decode_array(bench.array_payload(values)).tobytes() == values.tobytes()
    )
    payload = bench.array_payload(values)
    payload["real"][0][0] = False
    with pytest.raises(ValueError, match="boolean"):
        bench._decode_array(payload)


@pytest.mark.parametrize(
    "change",
    [
        {"precision_bits": 53},
        {"precision_bits": True},
        {"trials": True},
        {"trials": 0},
        {"beta": float("nan")},
        {"wmax": float("inf")},
        {"retention_fraction": 1.0},
        {"query_count": 189},
        {"small_query_count": 81},
        {"oracle_precision_bits": 128},
        {"construction_precision": "float64"},
        {"timeout_seconds": 0},
    ],
)
def test_invalid_protocol_controls_are_rejected(change: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        bench.RustBenchmarkConfig(**change)


def test_config_loader_requires_integer_schema_and_declares_default_protocol(
    tmp_path: Path,
) -> None:
    path = tmp_path / "bad.toml"
    path.write_text("schema_version = true\n")
    with pytest.raises(ValueError):
        bench.load_config(path)
    canonical = bench.load_config(
        Path(__file__).resolve().parents[1]
        / "experiments/rust_reference_benchmark.toml"
    )
    assert canonical.trials == 3 and canonical.query_count == 80
    assert canonical.rows_per_query == 64 and canonical.cache_size == 8
    assert canonical.small_query_count == 8 and canonical.warm_hits == 100
    grids = bench.query_grids(canonical)
    assert len({bench.fingerprint(row) for row in grids}) == 80
    assert np.all((grids > 0) & (grids < canonical.beta))
    assert not np.intersect1d(canonical.beta * bench.ORACLE_FRACTIONS, grids).size


@pytest.mark.parametrize(
    "field,value",
    [
        ("provider", "rust"),
        ("trial", False),
        ("source", {"git_commit": "0" * 40}),
        ("protocol", "faster_unchecked"),
        ("status", "complete"),
        ("accuracy_passed", False),
        ("protocol_passed", False),
    ],
)
def test_false_flags_and_forged_worker_identity_fail(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
    field: str,
    value: Any,
) -> None:
    record = copy.deepcopy(baseline)
    record[field] = value
    with pytest.raises(ValueError):
        bench.validate_worker(record, config, "python", 0, source)


@pytest.mark.parametrize(
    "mutation",
    ["hash", "imaginary", "shape", "oracle", "timing", "threads", "cache", "backend"],
)
def test_bad_complete_arrays_or_computation_metadata_are_rejected(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
    mutation: str,
) -> None:
    record = copy.deepcopy(baseline)
    if mutation == "hash":
        record["arrays"]["outputs"]["sha256"] = "0" * 64
    elif mutation == "imaginary":
        record["arrays"]["outputs"]["imag"][0][0][0][1] += 0.01
    elif mutation == "shape":
        record["arrays"]["coefficients"]["shape"] = [1]
    elif mutation == "oracle":
        record["oracle"]["time"]["max_error"] = 0.0
        record["oracle"]["time"]["canonical_max_error"] = "0.0"
    elif mutation == "timing":
        record["measurements"]["stream_seconds"] = float("inf")
    elif mutation == "threads":
        record["thread_controls"][bench.THREAD_VARIABLES[0]] = "4"
    elif mutation == "cache":
        record["calls"]["warm"] = 1
    else:
        record["basis_diagnostics"]["frequency_evaluation_backend"] = "rust"
    with pytest.raises(ValueError):
        bench.validate_worker(record, config, "python", 0, source)


def test_unavailable_native_backend_is_an_explicit_failed_record(
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(bench, "frozen_snapshot", lambda *args: source)
    monkeypatch.setattr(bench.platform, "platform", lambda: "test platform")
    record = bench.run_worker(
        {
            "config": asdict(config),
            "provider": "rust",
            "trial": 0,
            "source": source,
            "root": ".",
            "config_path": "unused",
        }
    )
    assert record["status"] == "native_unavailable"
    assert not record["accuracy_passed"] and not record["protocol_passed"]
    bench.validate_worker(record, config, "rust", 0, source)
    record["accuracy_passed"] = True
    with pytest.raises(ValueError, match="failed worker"):
        bench.validate_worker(record, config, "rust", 0, source)


def paired_records(
    baseline: dict[str, Any], config: bench.RustBenchmarkConfig
) -> list[dict[str, Any]]:
    records = []
    for trial in range(config.trials):
        for provider in bench.BACKENDS:
            record = copy.deepcopy(baseline)
            record.update(provider=provider, trial=trial)
            # Decision-function fixture: worker validation is tested separately.
            record["measurements"].update(
                application_seconds=40.0 if provider == "python" else 20.0,
                full_process_wall_seconds=42.0 if provider == "python" else 22.0,
                small_prefix_application_seconds=20.0 if provider == "python" else 19.0,
                warm_cached_seconds=0.001,
            )
            records.append(record)
    return records


def test_retention_requires_measured_full_process_gain_not_only_hot_loop(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
) -> None:
    records = paired_records(baseline, config)
    result = bench.summarize(records, config)
    assert result["retained"] and all(
        item["bitwise_equal"] for item in result["paired_trials"]
    )
    assert (
        result["reductions"]["small_prefix_application_seconds"]
        < config.retention_fraction
    )
    for record in records:
        if record["provider"] == "rust":
            record["measurements"]["full_process_wall_seconds"] = 40.0
    assert not bench.summarize(records, config)["retained"]


def test_retention_rejects_inaccurate_failed_unbalanced_and_different_input_trials(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
) -> None:
    records = paired_records(baseline, config)
    records[0]["accuracy_passed"] = False
    assert not bench.summarize(records, config)["retained"]
    records = paired_records(baseline, config)
    records[1]["arrays"]["observations"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="inputs differ"):
        bench.summarize(records, config)
    records = paired_records(baseline, config)
    records[2]["trial"] = 0
    with pytest.raises(ValueError, match="duplicate or missing"):
        bench.summarize(records, config)


def test_frozen_snapshot_hashes_native_sources_toolchain_and_binary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native_files = [
        "native/Cargo.toml",
        "native/Cargo.lock",
        "native/src/lib.rs",
        "rust-toolchain.toml",
    ]
    for name in native_files:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(name)
    binary = tmp_path / "_native.abi3.so"
    binary.write_bytes(b"native one")
    monkeypatch.setattr(
        bench,
        "source_snapshot",
        lambda *args: {"git_commit": "f" * 40, "source_sha256": {}},
    )

    def command(root: Path, args: list[str]) -> str:
        if args[0] == "git":
            return "" if "--others" in args else "\n".join(native_files)
        return "frozen toolchain"

    monkeypatch.setattr(bench, "_command", command)
    first = bench.frozen_snapshot(tmp_path, tmp_path / "config", binary)
    assert set(first["source_sha256"]) == set(native_files)
    assert first["toolchain"]["rustc"] == "frozen toolchain"
    binary.write_bytes(b"native two")
    assert (
        bench.frozen_snapshot(tmp_path, tmp_path / "config", binary)["native_binary"]
        != first["native_binary"]
    )
    monkeypatch.setattr(bench, "_command", lambda *args: "untracked.rs")
    with pytest.raises(ValueError, match="untracked"):
        bench.frozen_snapshot(tmp_path, tmp_path / "config", binary)


def test_fresh_process_timeouts_are_retained_and_order_alternates(
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(bench, "load_config", lambda *args: config)
    monkeypatch.setattr(bench, "_binary_path", lambda *args: None)
    monkeypatch.setattr(bench, "frozen_snapshot", lambda *args: source)
    monkeypatch.setattr(bench.platform, "platform", lambda: "test platform")

    def timeout(*args: Any, **kwargs: Any) -> Any:
        assert kwargs["timeout"] == config.timeout_seconds
        assert all(kwargs["env"][name] == "1" for name in bench.THREAD_VARIABLES)
        raise subprocess.TimeoutExpired(args[0], config.timeout_seconds)

    monkeypatch.setattr(bench.subprocess, "run", timeout)
    result = bench.run_benchmark(tmp_path / "unused", tmp_path)
    assert result["trial_order"] == [
        ["python", 0],
        ["rust", 0],
        ["rust", 1],
        ["python", 1],
        ["python", 2],
        ["rust", 2],
    ]
    assert all(record["status"] == "timeout" for record in result["records"])
    assert not result["summary"]["retained"]


def test_native_binary_and_compiled_source_revision_are_not_optional_metadata(
    baseline: dict[str, Any],
    config: bench.RustBenchmarkConfig,
) -> None:
    source = {
        "git_commit": "f" * 40,
        "native_binary": {
            "filename": "_native.abi3.so",
            "sha256": "a" * 64,
            "size_bytes": 1,
        },
    }
    record = copy.deepcopy(baseline)
    record.update(
        provider="rust",
        source=source,
        loaded_native_binary=source["native_binary"],
        native_build_info={
            "working_precision_bits": 128,
            "protocol_version": 1,
            "source_revision": "0" * 40,
        },
    )
    record["basis_diagnostics"].update(
        requested_reference_backend="rust", time_evaluation_backend="rust_mpfr_128"
    )
    with pytest.raises(ValueError, match="build identity"):
        bench.validate_worker(record, config, "rust", 0, source)
    record["native_build_info"]["source_revision"] = source["git_commit"]
    record["loaded_native_binary"] = {**source["native_binary"], "sha256": "b" * 64}
    with pytest.raises(ValueError, match="different native binary"):
        bench.validate_worker(record, config, "rust", 0, source)


def test_source_drift_after_workers_prevents_a_canonical_result(
    config: bench.RustBenchmarkConfig,
    source: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(bench, "load_config", lambda *args: config)
    monkeypatch.setattr(bench, "_binary_path", lambda *args: None)
    calls = [0]

    def snapshot(*args: Any) -> dict[str, Any]:
        calls[0] += 1
        return source if calls[0] == 1 else {**source, "git_commit": "0" * 40}

    monkeypatch.setattr(bench, "frozen_snapshot", snapshot)
    monkeypatch.setattr(
        bench.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 1, "", "worker failed"
        ),
    )
    with pytest.raises(ValueError, match="source/build/config changed"):
        bench.run_benchmark(tmp_path / "unused", tmp_path)
