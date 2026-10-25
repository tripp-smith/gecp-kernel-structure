"""Frozen fresh-process comparison of optional Rust reference time evaluation.

Use three alternating fresh trials per backend. The primary application timer
includes construction, observations/fit, native load/initialization and eighty
distinct public time evaluations with ordinary cache capacity eight. An actual
eight-query prefix and repeated cache hits are controls, not universal speedups.
Full process wall additionally includes imports, validation and serialization.
CLI exit zero means the experiment was recorded, not that Rust was retained.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.metadata
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

from experiments.dlr_application_benchmark import (
    array_payload,
    fingerprint,
)
from experiments.dlr_benchmark import THREAD_VARIABLES, _peak_rss, source_snapshot

PROTOCOL = "public_reference128_streaming_rust_v1"
BACKENDS = ("python", "rust")
FAILURES = {"native_unavailable", "runtime_failed", "worker_failed", "timeout"}
ORACLE_FRACTIONS = np.array([0.1234567, 0.4567891, 0.8765432])
FREQUENCY_MODES = np.array([-17, -2, 0, 4, 31], dtype=np.int64)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _positive(value: Any, name: str, *, zero: bool = False) -> float:
    _require(
        type(value) in (int, float)
        and math.isfinite(value)
        and (value >= 0 if zero else value > 0),
        name + " must be finite and " + ("nonnegative" if zero else "positive"),
    )
    return float(value)


@dataclass(frozen=True)
class RustBenchmarkConfig:
    beta: float = 20.0
    wmax: float = 50.0
    tolerance: float = 1e-10
    construction_precision: str = "auto"
    precision_bits: int = 128
    grid_order: int = 8
    max_rank: int = 128
    cache_size: int = 8
    trials: int = 3
    query_count: int = 80
    rows_per_query: int = 64
    small_query_count: int = 8
    warm_hits: int = 100
    timeout_seconds: float = 300.0
    retention_fraction: float = 0.2
    parity_absolute_target: float = 1e-12
    physical_absolute_target: float = 1e-8
    oracle_precision_bits: int = 512

    def __post_init__(self) -> None:
        for name in (
            "beta",
            "wmax",
            "tolerance",
            "timeout_seconds",
            "retention_fraction",
            "parity_absolute_target",
            "physical_absolute_target",
        ):
            _positive(getattr(self, name), name)
        _require(
            1 <= self.beta * self.wmax <= 16384, "finite cutoff in [1,16384] required"
        )
        _require(
            self.tolerance < 1 and self.retention_fraction < 1,
            "tolerance and retention fraction must be <1",
        )
        _require(
            self.construction_precision in {"auto", "reference"},
            "unsupported construction precision",
        )
        _require(
            type(self.precision_bits) is int and self.precision_bits == 128,
            "benchmark requires exactly 128 working bits",
        )
        for name, minimum, maximum in (
            ("grid_order", 3, 128),
            ("max_rank", 1, 128),
            ("cache_size", 1, 128),
            ("trials", 1, 16),
            ("query_count", 1, 188),
            ("rows_per_query", 2, 256),
            ("small_query_count", 1, self.query_count),
            ("warm_hits", 1, 10000),
            ("oracle_precision_bits", 256, 4096),
        ):
            value = getattr(self, name)
            _require(
                type(value) is int and minimum <= value <= maximum,
                name + " is outside its integer budget",
            )


def load_config(path: Path) -> RustBenchmarkConfig:
    raw = tomllib.loads(path.read_text())
    version = raw.pop("schema_version", None)
    _require(type(version) is int and version == 1, "unsupported config schema")
    return RustBenchmarkConfig(**raw)


def _decode_array(payload: Any) -> np.ndarray:
    """Decode complete complex bytes, including signed-zero components."""
    _require(isinstance(payload, dict), "array payload must be object")
    dtype = payload.get("dtype")
    _require(dtype in {"<c16", "<f8", "<i8"}, "unsupported packed array dtype")

    def numeric(value: Any) -> np.ndarray:
        def boolean(item: Any) -> bool:
            return isinstance(item, bool) or (
                isinstance(item, list) and any(boolean(part) for part in item)
            )

        _require(not boolean(value), "boolean numeric array component")
        array = np.asarray(value)
        _require(array.dtype.kind in "iuf", "nonnumeric array component")
        return array

    try:
        if dtype == "<c16":
            real, imag = numeric(payload["real"]), numeric(payload["imag"])
            _require(real.shape == imag.shape, "complex array component shape mismatch")
            array = np.empty(real.shape, dtype=dtype)
            array.real, array.imag = real, imag
        else:
            array = np.asarray(numeric(payload["values"]), dtype=dtype)
        _require(
            array_payload(array) == payload, "array shape/dtype/hash/content mismatch"
        )
        return array
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError("malformed array payload") from exc


def query_grids(config: RustBenchmarkConfig) -> np.ndarray:
    return np.asarray(
        [
            config.beta
            * (np.arange(config.rows_per_query) + 0.25 + 0.004 * index)
            / config.rows_per_query
            for index in range(config.query_count)
        ],
        dtype=np.float64,
    )


def physical_truth(tau: np.ndarray, beta: float) -> np.ndarray:
    """Independent analytic two-pole input; independent MP checks are separate."""
    times = np.asarray(tau, dtype=np.float64)
    value = np.zeros(times.shape)
    with np.errstate(all="ignore"):
        for weight, energy in ((0.4, -0.4), (0.6, 0.31)):
            term = (
                np.exp(-times * energy) / (1 + np.exp(-beta * energy))
                if energy >= 0
                else np.exp((beta - times) * energy) / (1 + np.exp(beta * energy))
            )
            value -= weight * term
    return value[..., None, None] * np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _command(root: Path, arguments: list[str]) -> str:
    result = subprocess.run(
        arguments, cwd=root, capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def _binary_path(root: Path, supplied: Path | None) -> Path | None:
    if supplied is not None:
        _require(supplied.is_file(), "supplied native binary does not exist")
        return supplied.resolve()
    spec = importlib.machinery.PathFinder.find_spec(
        "kernelgecp._native", [str(root / "kernelgecp")]
    )
    return Path(spec.origin).resolve() if spec is not None and spec.origin else None


def frozen_snapshot(
    root: Path, config_path: Path, binary: Path | None
) -> dict[str, Any]:
    result = source_snapshot(root, config_path)
    _require(
        not _command(
            root,
            [
                "git",
                "ls-files",
                "--others",
                "--exclude-standard",
                "--",
                "native",
                "rust-toolchain.toml",
                ".cargo",
            ],
        ),
        "untracked native/build code prevents a frozen run",
    )
    native = _command(
        root, ["git", "ls-files", "--", "native", "rust-toolchain.toml", ".cargo"]
    ).splitlines()
    _require(
        "native/Cargo.toml" in native
        and "native/Cargo.lock" in native
        and "rust-toolchain.toml" in native,
        "locked native sources/toolchain required",
    )
    result["source_sha256"].update({name: _sha(root / name) for name in native})
    result["native_binary"] = (
        {
            "filename": binary.name,
            "sha256": _sha(binary),
            "size_bytes": binary.stat().st_size,
        }
        if binary is not None
        else None
    )
    result["toolchain"] = {
        "rustc": _command(root, ["rustc", "--version", "--verbose"]),
        "cargo": _command(root, ["cargo", "--version", "--verbose"]),
    }
    return result


def _versions() -> dict[str, str]:
    return {
        name: importlib.metadata.version(name)
        for name in ("kernelgecp", "numpy", "scipy", "mpmath")
    }


def _identity(
    config: RustBenchmarkConfig, provider: str, trial: int, source: dict[str, Any]
) -> dict[str, Any]:
    return {
        "protocol": PROTOCOL,
        "provider": provider,
        "trial": trial,
        "config": asdict(config),
        "source": source,
    }


def _oracle(
    beta: float,
    tau_nodes: np.ndarray,
    omega_nodes: np.ndarray,
    coefficients: np.ndarray,
    times: np.ndarray,
    time_output: np.ndarray,
    frequency_output: np.ndarray,
    bits: int,
) -> dict[str, Any]:
    """QR of the exact-physical centered kernel, applied to frozen coefficients."""
    ctx = mp.mp.clone()
    ctx.prec = bits
    digits = math.ceil(bits * math.log10(2)) + 20
    exact_beta = ctx.mpf(float(beta))

    def scalar(value: Any) -> Any:
        return ctx.mpc(float(np.real(value)), float(np.imag(value)))

    def kernel(time: float, omega: float) -> Any:
        t, w = ctx.mpf(float(time)) / exact_beta, ctx.mpf(float(omega))
        return ctx.exp((ctx.mpf(".5") - t) * w) / (2 * ctx.cosh(w / 2))

    core = ctx.matrix([[kernel(t, w) for w in omega_nodes] for t in tau_nodes])
    q, r = ctx.qr(core.T)

    def solve(rhs: Any) -> Any:
        projected = q.T * rhs
        result = ctx.matrix(r.cols, 1)
        for index in range(r.cols - 1, -1, -1):
            result[index] = (
                projected[index]
                - ctx.fsum(r[index, j] * result[j] for j in range(index + 1, r.cols))
            ) / r[index, index]
        return result

    amplitude = np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])
    errors: dict[str, Any] = {}
    for name, points, actual in (
        ("time", times, time_output),
        ("frequency", FREQUENCY_MODES, frequency_output),
    ):
        maximum, physical = ctx.zero, ctx.zero
        predictions, truths = [], []
        for index, point in enumerate(points):
            if name == "time":
                rhs = ctx.matrix([kernel(float(point), w) for w in omega_nodes])
                truth = -ctx.fsum(
                    ctx.mpf(weight)
                    * ctx.exp(-ctx.mpf(float(point)) * ctx.mpf(energy))
                    / (1 + ctx.exp(-exact_beta * ctx.mpf(energy)))
                    for weight, energy in ((0.4, -0.4), (0.6, 0.31))
                )
            else:
                rhs = ctx.matrix(
                    [
                        -exact_beta
                        / (ctx.j * ctx.pi * (2 * int(point) + 1) - ctx.mpf(float(w)))
                        for w in omega_nodes
                    ]
                )
                truth = ctx.fsum(
                    ctx.mpf(weight)
                    / (
                        ctx.j * ctx.pi * (2 * int(point) + 1) / exact_beta
                        - ctx.mpf(energy)
                    )
                    for weight, energy in ((0.4, -0.4), (0.6, 0.31))
                )
            design = solve(rhs)
            prediction_row, truth_row = [], []
            for a, b in np.ndindex(2, 2):
                ideal = ctx.fsum(
                    design[k] * scalar(coefficients[k, a, b])
                    for k in range(len(tau_nodes))
                )
                value = scalar(actual[index, a, b])
                physical_truth_value = truth * scalar(amplitude[a, b])
                maximum = max(maximum, abs(value - ideal))
                physical = max(physical, abs(value - physical_truth_value))
                prediction_row.append(
                    [ctx.nstr(ideal.real, digits), ctx.nstr(ideal.imag, digits)]
                )
                truth_row.append(
                    [
                        ctx.nstr(physical_truth_value.real, digits),
                        ctx.nstr(physical_truth_value.imag, digits),
                    ]
                )
            predictions.append(prediction_row)
            truths.append(truth_row)
        errors[name] = {
            "max_error": float(maximum),
            "canonical_max_error": ctx.nstr(maximum, digits),
            "physical_max_error": float(physical),
            "canonical_physical_max_error": ctx.nstr(physical, digits),
            "ideal_predictions": predictions,
            "physical_truth": truths,
        }
    return {
        "precision_bits": bits,
        "formulation": "centered_exact_physical_core_qr_on_frozen_coefficients",
        "coordinate_scope": "exact_binary64_tau_divided_by_exact_binary64_beta",
        "coefficient_scope": "actual_binary64_fit_coefficients_not_an_ideal_LS_refit",
        "certified": False,
        **errors,
    }


def run_worker(request: dict[str, Any]) -> dict[str, Any]:
    config = RustBenchmarkConfig(**request["config"])
    provider, trial, expected = request["provider"], request["trial"], request["source"]
    _require(
        provider in BACKENDS and type(trial) is int and 0 <= trial < config.trials,
        "invalid worker identity",
    )
    root, config_path = Path(request["root"]), Path(request["config_path"])
    binary = (
        Path(request["native_binary_path"])
        if request.get("native_binary_path")
        else None
    )
    result = {
        **_identity(config, provider, trial, expected),
        "status": "runtime_failed",
        "accuracy_passed": False,
        "protocol_passed": False,
    }
    try:
        _require(
            frozen_snapshot(root, config_path, binary) == expected,
            "worker source/build/config changed before execution",
        )
        if provider == "rust" and expected["native_binary"] is None:
            return {
                **result,
                "status": "native_unavailable",
                "error": "no importable optional native module was found",
            }
        from kernelgecp import FermionicBasis
        from kernelgecp import _reference_backend as adapter

        native_times = {
            "load_seconds": 0.0,
            "setup_seconds": 0.0,
            "load_calls": 0,
            "setup_calls": 0,
        }
        original_load, original_init = (
            adapter._load_native,
            adapter.RustReferenceBackend.__init__,
        )

        def timed_load() -> Any:
            start = time.perf_counter()
            try:
                return original_load()
            finally:
                native_times["load_seconds"] += time.perf_counter() - start
                native_times["load_calls"] += 1

        def timed_init(self: Any, *args: Any, **kwargs: Any) -> None:
            start = time.perf_counter()
            try:
                original_init(self, *args, **kwargs)
            finally:
                native_times["setup_seconds"] += time.perf_counter() - start
                native_times["setup_calls"] += 1

        adapter._load_native, adapter.RustReferenceBackend.__init__ = (
            timed_load,
            timed_init,
        )
        try:
            started = time.perf_counter()
            basis = FermionicBasis(
                config.beta,
                config.wmax,
                tol=config.tolerance,
                precision=config.construction_precision,
                precision_bits=config.precision_bits,
                grid_order=config.grid_order,
                max_rank=config.max_rank,
                cache_size=config.cache_size,
                reference_backend=provider,
            )
            constructed = time.perf_counter()
            calls = {name: 0 for name in ("fit", "stream", "warm", "validation")}
            phase = "fit"
            underlying = basis._numerics.design_tau

            def counted(tau: Any) -> Any:
                calls[phase] += 1
                return underlying(tau)

            basis._numerics.design_tau = counted
            observations = physical_truth(basis.tau_nodes, config.beta)
            fit = basis.fit_tau(observations)
            fitted = time.perf_counter()
            grids = query_grids(config)
            phase = "stream"
            evaluated, prefix_seconds = [], None
            stream_start = time.perf_counter()
            for index, tau in enumerate(grids):
                evaluated.append(fit.evaluate_tau(tau))
                if index + 1 == config.small_query_count:
                    prefix_seconds = time.perf_counter() - started
            finished = time.perf_counter()
            outputs = np.asarray(evaluated)
            workload_rss = _peak_rss()
            phase = "warm"
            warm_start = time.perf_counter()
            for _ in range(config.warm_hits):
                last_warm = fit.evaluate_tau(grids[-1])
            warm_seconds = time.perf_counter() - warm_start
            _require(
                np.array_equal(last_warm, outputs[-1]), "warm cached output changed"
            )
            phase = "validation"
            oracle_tau = config.beta * ORACLE_FRACTIONS
            _require(
                not np.intersect1d(oracle_tau, grids.ravel()).size
                and not np.intersect1d(oracle_tau, basis.tau_nodes).size,
                "independent oracle queries overlap training/workload",
            )
            oracle_time = fit.evaluate_tau(oracle_tau)
            frequency_output = fit.evaluate_matsubara(FREQUENCY_MODES)
            reference = _oracle(
                config.beta,
                basis.tau_nodes,
                basis.omega_nodes,
                fit.coefficients,
                oracle_tau,
                oracle_time,
                frequency_output,
                config.oracle_precision_bits,
            )
            maximum_physical = float(
                np.max(np.abs(outputs - physical_truth(grids, config.beta)))
            )
            diagnostics = asdict(basis.diagnostics)
            native_info = dict(basis.diagnostics.native_backend_info)
            loaded_binary = None
            if provider == "rust":
                module = sys.modules.get("kernelgecp._native")
                _require(
                    module is not None
                    and isinstance(getattr(module, "__file__", None), str),
                    "native module identity missing",
                )
                path = Path(module.__file__)
                loaded_binary = {
                    "filename": path.name,
                    "sha256": _sha(path),
                    "size_bytes": path.stat().st_size,
                }
                _require(
                    loaded_binary == expected["native_binary"],
                    "loaded native binary differs from frozen build",
                )
            protocol_passed = (
                basis.diagnostics.selection_precision_bits == config.precision_bits
                and basis.diagnostics.evaluation_precision_bits == config.precision_bits
                and calls["stream"] == config.query_count
                and calls["warm"] == 0
            )
            accuracy_passed = (
                maximum_physical <= config.physical_absolute_target
                and all(
                    reference[name]["physical_max_error"]
                    <= config.physical_absolute_target
                    and reference[name]["max_error"] <= config.parity_absolute_target
                    for name in ("time", "frequency")
                )
            )
            validation_rss = _peak_rss()
            result.update(
                status="completed"
                if protocol_passed and accuracy_passed
                else ("accuracy_failed" if protocol_passed else "protocol_failed"),
                accuracy_passed=accuracy_passed,
                protocol_passed=protocol_passed,
                rank=basis.rank,
                basis_diagnostics=diagnostics,
                fit_diagnostics=asdict(fit.diagnostics),
                native_build_info=native_info,
                loaded_native_binary=loaded_binary,
                calls=calls,
                thread_controls={
                    name: os.environ.get(name) for name in THREAD_VARIABLES
                },
                thread_control_scope="environment_settings_not_observed_os_thread_counts",
                versions=_versions(),
                oracle=reference,
                physical_stream_max_error=maximum_physical,
                measurements={
                    "application_seconds": finished - started,
                    "construction_seconds": constructed - started,
                    "observations_and_fit_seconds": fitted - constructed,
                    "stream_seconds": finished - stream_start,
                    "small_prefix_application_seconds": prefix_seconds,
                    "warm_cached_seconds": warm_seconds,
                    "native_lazy_load_and_checks_seconds": native_times["load_seconds"],
                    "native_setup_including_load_seconds": native_times[
                        "setup_seconds"
                    ],
                    "native_packing_and_initialization_seconds": max(
                        0.0,
                        native_times["setup_seconds"] - native_times["load_seconds"],
                    ),
                    "rss_after_workload_bytes": workload_rss,
                    "rss_before_payload_serialization_bytes": validation_rss,
                },
                measurement_scopes={
                    "application": (
                        "construction+observations+fit+native_import_setup+"
                        "public_evaluations;common_imports_and_validation_excluded"
                    ),
                    "small_prefix": (
                        "actual_application_timer_after_first_queries;"
                        "not_standalone_process_wall"
                    ),
                    "warm": "public_last_grid_calls;no_underlying_design_expected",
                    "native_setup": "nested_in_construction;load_is_nested",
                    "rss": (
                        "process_HWM_includes_imports_harness_native_and_numerics;"
                        "before_array_payload_and_JSON"
                    ),
                },
                arrays={
                    name: array_payload(value)
                    for name, value in {
                        "tau_nodes": basis.tau_nodes,
                        "omega_nodes": basis.omega_nodes,
                        "observations": observations,
                        "coefficients": fit.coefficients,
                        "query_grids": grids,
                        "outputs": outputs,
                        "oracle_tau": oracle_tau,
                        "oracle_time_output": oracle_time,
                        "frequency_indices": FREQUENCY_MODES,
                        "frequency_output": frequency_output,
                        "warm_output": last_warm,
                    }.items()
                },
            )
        finally:
            adapter._load_native, adapter.RustReferenceBackend.__init__ = (
                original_load,
                original_init,
            )
        _require(
            frozen_snapshot(root, config_path, binary) == expected,
            "worker source/build/config changed during execution",
        )
    except Exception as exc:
        result.update(
            status="runtime_failed",
            accuracy_passed=False,
            protocol_passed=False,
            error=f"{type(exc).__name__}: {exc}",
        )
    return result


def _validate_oracle_payload(
    record: dict[str, Any], arrays: dict[str, np.ndarray], config: RustBenchmarkConfig
) -> None:
    """Recompute errors from complete serialized MP points and frozen outputs."""
    ctx = mp.mp.clone()
    ctx.prec = config.oracle_precision_bits
    amplitude = np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])

    def scalar(value: Any) -> Any:
        return ctx.mpc(float(np.real(value)), float(np.imag(value)))

    for name, points, actual in (
        ("time", arrays["oracle_tau"], arrays["oracle_time_output"]),
        ("frequency", arrays["frequency_indices"], arrays["frequency_output"]),
    ):
        oracle = record["oracle"][name]
        predictions, truths = (
            oracle.get("ideal_predictions"),
            oracle.get("physical_truth"),
        )
        _require(
            isinstance(predictions, list)
            and isinstance(truths, list)
            and len(predictions) == len(truths) == len(points),
            "oracle serialized point count mismatch",
        )
        maximum, physical = ctx.zero, ctx.zero
        for row, point in enumerate(points):
            _require(
                isinstance(predictions[row], list)
                and isinstance(truths[row], list)
                and len(predictions[row]) == len(truths[row]) == 4,
                "oracle complex batch shape mismatch",
            )
            if name == "time":
                truth = -ctx.fsum(
                    ctx.mpf(weight)
                    * ctx.exp(-ctx.mpf(float(point)) * ctx.mpf(energy))
                    / (1 + ctx.exp(-ctx.mpf(config.beta) * ctx.mpf(energy)))
                    for weight, energy in ((0.4, -0.4), (0.6, 0.31))
                )
            else:
                truth = ctx.fsum(
                    ctx.mpf(weight)
                    / (
                        ctx.j * ctx.pi * (2 * int(point) + 1) / ctx.mpf(config.beta)
                        - ctx.mpf(energy)
                    )
                    for weight, energy in ((0.4, -0.4), (0.6, 0.31))
                )
            for a, b in np.ndindex(2, 2):
                index = 2 * a + b
                values = (predictions[row][index], truths[row][index])
                _require(
                    all(
                        isinstance(value, list)
                        and len(value) == 2
                        and all(isinstance(item, str) for item in value)
                        for value in values
                    ),
                    "oracle complex point malformed",
                )
                ideal, stored_truth = (ctx.mpc(*value) for value in values)
                _require(
                    ctx.isfinite(ideal) and ctx.isfinite(stored_truth),
                    "oracle point nonfinite",
                )
                _require(
                    stored_truth == truth * scalar(amplitude[a, b]),
                    "oracle physical truth mismatch",
                )
                value = scalar(actual[row, a, b])
                maximum = max(maximum, abs(value - ideal))
                physical = max(physical, abs(value - stored_truth))
        for field, error in (("max_error", maximum), ("physical_max_error", physical)):
            _require(
                oracle.get(field) == float(error)
                and ctx.mpf(oracle["canonical_" + field]) == error,
                "oracle error disagrees with serialized points",
            )


def validate_worker(
    record: Any,
    config: RustBenchmarkConfig,
    provider: str,
    trial: int,
    source: dict[str, Any],
) -> dict[str, Any]:
    _require(isinstance(record, dict), "worker payload must be object")
    _require(type(record.get("trial")) is int, "worker trial must be an integer")
    _require(isinstance(record.get("config"), dict), "worker config must be object")
    RustBenchmarkConfig(**record["config"])
    for name, expected in _identity(config, provider, trial, source).items():
        _require(record.get(name) == expected, "worker identity mismatch: " + name)
    status = record.get("status")
    _require(
        isinstance(status, str)
        and status in FAILURES | {"completed", "accuracy_failed", "protocol_failed"},
        "unknown worker status",
    )
    for name in ("accuracy_passed", "protocol_passed"):
        _require(type(record.get(name)) is bool, "worker success flags must be bool")
    if status in FAILURES:
        _require(
            not record["accuracy_passed"] and not record["protocol_passed"],
            "failed worker claims success",
        )
        return record
    rank = record.get("rank")
    _require(type(rank) is int and 1 <= rank <= config.max_rank, "invalid worker rank")
    arrays = {
        name: _decode_array(value) for name, value in record.get("arrays", {}).items()
    }
    shapes = {
        "tau_nodes": (rank,),
        "omega_nodes": (rank,),
        "observations": (rank, 2, 2),
        "coefficients": (rank, 2, 2),
        "query_grids": (config.query_count, config.rows_per_query),
        "outputs": (config.query_count, config.rows_per_query, 2, 2),
        "oracle_tau": (3,),
        "oracle_time_output": (3, 2, 2),
        "frequency_indices": (5,),
        "frequency_output": (5, 2, 2),
        "warm_output": (config.rows_per_query, 2, 2),
    }
    _require(set(arrays) == set(shapes), "worker array manifest mismatch")
    complex_names = {
        "observations",
        "coefficients",
        "outputs",
        "oracle_time_output",
        "frequency_output",
        "warm_output",
    }
    for name, shape in shapes.items():
        _require(arrays[name].shape == shape, "worker array shape mismatch: " + name)
        _require(
            arrays[name].dtype.str
            == (
                "<c16"
                if name in complex_names
                else ("<i8" if name == "frequency_indices" else "<f8")
            ),
            "worker array dtype mismatch: " + name,
        )
    _require(
        np.array_equal(arrays["query_grids"], query_grids(config)),
        "worker query configuration mismatch",
    )
    _require(
        np.array_equal(arrays["oracle_tau"], config.beta * ORACLE_FRACTIONS)
        and np.array_equal(arrays["frequency_indices"], FREQUENCY_MODES),
        "worker validation query mismatch",
    )
    _require(
        not np.intersect1d(arrays["oracle_tau"], arrays["query_grids"].ravel()).size
        and not np.intersect1d(arrays["oracle_tau"], arrays["tau_nodes"]).size,
        "worker validation overlaps observations",
    )
    _require(
        np.array_equal(arrays["warm_output"], arrays["outputs"][-1]),
        "warm output differs",
    )
    physical_error = float(
        np.max(
            np.abs(
                arrays["outputs"] - physical_truth(arrays["query_grids"], config.beta)
            )
        )
    )
    _require(
        record.get("physical_stream_max_error") == physical_error,
        "worker physical error mismatch",
    )
    _require(
        record.get("thread_controls") == {name: "1" for name in THREAD_VARIABLES},
        "worker thread control mismatch",
    )
    _require(
        record.get("thread_control_scope")
        == "environment_settings_not_observed_os_thread_counts",
        "invalid thread control scope",
    )
    measurements = record.get("measurements")
    _require(isinstance(measurements, dict), "missing worker measurements")
    for name in (
        "application_seconds",
        "construction_seconds",
        "observations_and_fit_seconds",
        "stream_seconds",
        "small_prefix_application_seconds",
        "warm_cached_seconds",
        "native_lazy_load_and_checks_seconds",
        "native_setup_including_load_seconds",
        "native_packing_and_initialization_seconds",
    ):
        _positive(measurements.get(name), name, zero=True)
    _require(
        measurements["application_seconds"]
        >= measurements["small_prefix_application_seconds"]
        and measurements["application_seconds"] >= measurements["construction_seconds"],
        "inconsistent application timing",
    )
    allowance = (
        64 * np.finfo(np.float64).eps * max(1.0, measurements["application_seconds"])
    )
    _require(
        measurements["construction_seconds"]
        + measurements["observations_and_fit_seconds"]
        + measurements["stream_seconds"]
        <= measurements["application_seconds"] + allowance,
        "component timing exceeds application scope",
    )
    _require(
        measurements["native_setup_including_load_seconds"]
        >= measurements["native_lazy_load_and_checks_seconds"],
        "native timing scopes inconsistent",
    )
    for name in ("rss_after_workload_bytes", "rss_before_payload_serialization_bytes"):
        _require(
            measurements.get(name) is None
            or (type(measurements[name]) is int and measurements[name] > 0),
            "invalid process RSS",
        )
    calls = record.get("calls")
    _require(
        isinstance(calls, dict)
        and set(calls) == {"fit", "stream", "warm", "validation"}
        and all(type(v) is int and v >= 0 for v in calls.values()),
        "invalid backend call counts",
    )
    diagnostic = record.get("basis_diagnostics")
    _require(isinstance(diagnostic, dict), "missing backend diagnostics")
    _require(
        diagnostic.get("requested_reference_backend") == provider,
        "requested backend mismatch",
    )
    active_reference = diagnostic.get("evaluation_precision_bits") == 128
    _require(
        diagnostic.get("time_evaluation_backend")
        == (
            ("rust_mpfr_128" if provider == "rust" else "python_mpmath_reference")
            if active_reference
            else "scipy_float64"
        ),
        "time backend mismatch",
    )
    _require(
        diagnostic.get("frequency_evaluation_backend")
        == ("python_mpmath_reference" if active_reference else "scipy_float64")
        and diagnostic.get("fit_backend") == "numpy_svd_float64_or_complex128",
        "unaccelerated backend metadata mismatch",
    )
    oracle = record.get("oracle")
    _require(
        isinstance(oracle, dict)
        and oracle.get("precision_bits") == config.oracle_precision_bits
        and oracle.get("certified") is False,
        "invalid independent oracle metadata",
    )
    _validate_oracle_payload(record, arrays, config)
    protocol_passed = (
        diagnostic.get("selection_precision_bits") == config.precision_bits
        and diagnostic.get("evaluation_precision_bits") == config.precision_bits
        and calls["stream"] == config.query_count
        and calls["warm"] == 0
    )
    accuracy_passed = physical_error <= config.physical_absolute_target and all(
        oracle[name]["physical_max_error"] <= config.physical_absolute_target
        and oracle[name]["max_error"] <= config.parity_absolute_target
        for name in ("time", "frequency")
    )
    _require(
        record["protocol_passed"] == protocol_passed
        and record["accuracy_passed"] == accuracy_passed,
        "worker false pass/fail flag",
    )
    expected_status = (
        "completed"
        if protocol_passed and accuracy_passed
        else ("accuracy_failed" if protocol_passed else "protocol_failed")
    )
    _require(status == expected_status, "worker status disagrees with recomputed flags")
    if provider == "rust" and active_reference:
        _require(
            record.get("loaded_native_binary") == source["native_binary"],
            "worker loaded a different native binary",
        )
        info = record.get("native_build_info")
        _require(
            isinstance(info, dict)
            and type(info.get("protocol_version")) is int
            and info.get("working_precision_bits") == 128
            and info.get("protocol_version") == 1
            and info.get("source_revision") == source["git_commit"],
            "native build identity/precision mismatch",
        )
    else:
        _require(
            record.get("loaded_native_binary") is None
            and record.get("native_build_info") == {},
            "Python baseline secretly reports a native backend",
        )
    json.dumps(record, allow_nan=False)
    return record


def summarize(
    records: list[dict[str, Any]], config: RustBenchmarkConfig
) -> dict[str, Any]:
    usable = len(records) == 2 * config.trials and all(
        record.get("status") == "completed"
        and record.get("accuracy_passed") is True
        and record.get("protocol_passed") is True
        for record in records
    )
    summary: dict[str, Any] = {
        "retained": False,
        "all_trials_completed": usable,
        "criterion": (
            "declared_median_application_and_fullprocess_reduction;"
            "parity_and_physical_checks_required"
        ),
        "retention_fraction": config.retention_fraction,
    }
    if not usable:
        return {**summary, "reason": "missing_failed_or_inaccurate_trial"}
    by_backend = {
        name: [record for record in records if record["provider"] == name]
        for name in BACKENDS
    }
    _require(
        all(len(value) == config.trials for value in by_backend.values()),
        "unbalanced trial providers",
    )
    _require(
        all(
            {record["trial"] for record in rows} == set(range(config.trials))
            for rows in by_backend.values()
        ),
        "duplicate or missing trial identifiers",
    )
    comparisons = []
    for trial in range(config.trials):
        pair = {
            name: next(
                record for record in by_backend[name] if record["trial"] == trial
            )
            for name in BACKENDS
        }
        for name in (
            "tau_nodes",
            "omega_nodes",
            "observations",
            "query_grids",
            "oracle_tau",
            "frequency_indices",
        ):
            _require(
                pair["python"]["arrays"][name] == pair["rust"]["arrays"][name],
                "paired trial inputs differ: " + name,
            )
        differences = {}
        bitwise = True
        for name in (
            "coefficients",
            "outputs",
            "oracle_time_output",
            "frequency_output",
            "warm_output",
        ):
            left, right = (
                _decode_array(pair[backend]["arrays"][name]) for backend in BACKENDS
            )
            differences[name] = float(np.max(np.abs(left - right)))
            bitwise = bitwise and fingerprint(left) == fingerprint(right)
        oracle_parity = all(
            pair["rust"]["oracle"][name]["max_error"]
            <= max(
                config.parity_absolute_target,
                2 * pair["python"]["oracle"][name]["max_error"],
            )
            for name in ("time", "frequency")
        )
        comparisons.append(
            {
                "trial": trial,
                "bitwise_equal": bitwise,
                "absolute_differences": differences,
                "parity_passed": max(differences.values())
                <= config.parity_absolute_target
                and oracle_parity,
            }
        )
    medians = {
        backend: {
            name: statistics.median(
                record["measurements"][name] for record in by_backend[backend]
            )
            for name in (
                "application_seconds",
                "full_process_wall_seconds",
                "small_prefix_application_seconds",
                "warm_cached_seconds",
            )
        }
        for backend in BACKENDS
    }
    reductions = {
        name: 1 - medians["rust"][name] / medians["python"][name]
        for name in medians["python"]
    }
    retained = (
        all(item["parity_passed"] for item in comparisons)
        and reductions["application_seconds"] >= config.retention_fraction
        and reductions["full_process_wall_seconds"] >= config.retention_fraction
    )
    return {
        **summary,
        "retained": retained,
        "reason": "retained_for_streaming_workload_only"
        if retained
        else "threshold_or_parity_not_met",
        "paired_trials": comparisons,
        "medians": medians,
        "reductions": reductions,
        "small_prefix_scope": (
            "observed_construction_fit_setup_and_first_queries;"
            "not_standalone_process_wall"
        ),
        "warm_scope": "cached_public_calls;no_native_design_acceleration_claim",
    }


def run_benchmark(
    config_path: Path, root: Path, native_binary: Path | None = None
) -> dict[str, Any]:
    config, root, config_path = (
        load_config(config_path),
        root.resolve(),
        config_path.resolve(),
    )
    binary = _binary_path(root, native_binary)
    source = frozen_snapshot(root, config_path, binary)
    records = []
    with tempfile.TemporaryDirectory(prefix="kernelgecp-rust-benchmark-") as temporary:
        for trial in range(config.trials):
            for provider in BACKENDS if trial % 2 == 0 else BACKENDS[::-1]:
                request = {
                    **_identity(config, provider, trial, source),
                    "root": str(root),
                    "config_path": str(config_path),
                    "native_binary_path": str(binary) if binary is not None else None,
                }
                request_path = Path(temporary) / "request.json"
                request_path.write_text(json.dumps(request))
                started = time.perf_counter()
                try:
                    process = subprocess.run(
                        [
                            sys.executable,
                            "-m",
                            "experiments.rust_reference_benchmark",
                            "--worker-request",
                            str(request_path),
                        ],
                        cwd=root,
                        env={**os.environ, **{name: "1" for name in THREAD_VARIABLES}},
                        capture_output=True,
                        text=True,
                        timeout=config.timeout_seconds,
                        check=False,
                    )
                    wall = time.perf_counter() - started
                    if process.returncode:
                        record = {
                            **_identity(config, provider, trial, source),
                            "status": "worker_failed",
                            "accuracy_passed": False,
                            "protocol_passed": False,
                            "exit_code": process.returncode,
                            "error": process.stderr[-4000:],
                        }
                    else:
                        record = validate_worker(
                            json.loads(process.stdout), config, provider, trial, source
                        )
                    record.setdefault("measurements", {})[
                        "full_process_wall_seconds"
                    ] = wall
                except subprocess.TimeoutExpired:
                    record = {
                        **_identity(config, provider, trial, source),
                        "status": "timeout",
                        "accuracy_passed": False,
                        "protocol_passed": False,
                        "measurements": {
                            "full_process_wall_seconds": time.perf_counter() - started
                        },
                    }
                records.append(record)
                print(
                    f"{provider} trial {trial + 1}: {record['status']}",
                    file=sys.stderr,
                    flush=True,
                )
    _require(
        frozen_snapshot(root, config_path, binary) == source,
        "source/build/config changed during benchmark",
    )
    return {
        "schema_version": 1,
        "protocol": PROTOCOL,
        "source": source,
        "config": asdict(config),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "versions": _versions(),
        "trial_order": [[record["provider"], record["trial"]] for record in records],
        "records": records,
        "summary": summarize(records, config),
        "limitations": [
            "Timing and process HWM describe one machine, not a universal speedup.",
            "Application timing includes native loading/setup; full-process wall "
            "also includes imports, provenance, independent checks and JSON.",
            "The eight-query prefix is an observed application-timer control; "
            "warm hits already reuse the public cache.",
            "Only reference time kernels/triangular substitution use Rust. "
            "Frequency evaluation and NumPy SVD/products are unchanged.",
            "Independent physical-coordinate QR is numerical consistency "
            "evidence, not a continuous-domain or Lean certificate.",
            "Thread variables are settings; actual OS counts are not observed.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native-binary", type=Path)
    parser.add_argument("--worker-request", type=Path)
    args = parser.parse_args()
    if args.worker_request is not None:
        print(
            json.dumps(
                run_worker(json.loads(args.worker_request.read_text())), allow_nan=False
            )
        )
        return
    if args.config is None or args.output is None:
        parser.error("--config and --output are required")
    result = run_benchmark(
        args.config, Path(__file__).resolve().parents[1], args.native_binary
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        stream.write(
            json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )
    print(
        f"Wrote Rust experiment to {args.output}; "
        f"retained={result['summary']['retained']}"
    )


if __name__ == "__main__":
    main()
