"""Matched physical time/frequency applications with isolated calibration workers.

Native eps is selected using calibration samples only. Selected construction and
reuse run in another fresh process, so calibration allocations are absent from
its peak RSS. Measurements are machine-dependent observations, not speedups.
"""

from __future__ import annotations

import argparse
import hashlib
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
from dataclasses import asdict, dataclass, is_dataclass
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np

from experiments.dlr_adapters import AdapterUnavailableError, build_adapter
from experiments.dlr_benchmark import (
    METHODS,
    THREAD_VARIABLES,
    _peak_rss,
    _preload_provider,
    _timed,
    _versions,
    provenance_digest,
    source_snapshot,
)
from experiments.dlr_physics import (
    BetheSolverPlan,
    frequency_index_meshes,
    matrix_amplitudes,
    matsubara_green_function,
    physical_tau_meshes,
    semicircle_matsubara,
    semicircle_tau,
    solver_index_mesh,
    solver_validation_index_mesh,
)
from experiments.dlr_workloads import green_function, make_pole_workloads
from kernelgecp import FermionicBasis, LeastSquaresPlan

PROTOCOL = "matched_physical_time_frequency_svd_v1"
CANDIDATE_STATUSES = {
    "accuracy_passed",
    "accuracy_failed",
    "import_unavailable",
    "runtime_failed",
}


def _positive(value: Any, name: str, *, below_one: bool = False) -> None:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0
        or (below_one and value >= 1)
    ):
        raise ValueError(
            f"{name} must be finite and positive" + (" and < 1" if below_one else "")
        )


@dataclass(frozen=True)
class ApplicationCase:
    beta: float
    wmax: float
    target: float = 1e-8
    solver: bool = False

    def __post_init__(self) -> None:
        for name in ("beta", "wmax", "target"):
            _positive(getattr(self, name), name, below_one=name == "target")
        if not math.isfinite(self.beta * self.wmax) or self.beta * self.wmax < 1:
            raise ValueError("finite beta*wmax >= 1 is required")
        if not isinstance(self.solver, bool):
            raise ValueError("solver must be bool")


@dataclass(frozen=True)
class ApplicationConfig:
    cases: tuple[ApplicationCase, ...]
    methods: tuple[str, ...] = METHODS
    native_eps_candidates: tuple[float, ...] = (1e-6, 1e-8, 1e-10, 1e-12)
    grid_order: int = 8
    max_rank: int = 128
    precision_bits: int = 128
    repeats: int = 3
    warmups: int = 1
    fit_rcond: float = 1e-14
    timeout_seconds: float = 600.0
    solver_mixing: float = 0.5
    solver_tolerance: float = 1e-9
    solver_max_iterations: int = 200

    def __post_init__(self) -> None:
        if not self.cases or not all(
            isinstance(c, ApplicationCase) for c in self.cases
        ):
            raise ValueError("cases must be nonempty ApplicationCase records")
        if (
            not self.methods
            or len(set(self.methods)) != len(self.methods)
            or any(method not in METHODS for method in self.methods)
        ):
            raise ValueError("methods must be distinct supported names")
        if not self.native_eps_candidates:
            raise ValueError("native eps candidates must be nonempty")
        for eps in self.native_eps_candidates:
            _positive(eps, "native eps", below_one=True)
        if any(
            a <= b
            for a, b in zip(
                self.native_eps_candidates, self.native_eps_candidates[1:], strict=False
            )
        ):
            raise ValueError("native eps candidates must be strictly descending")
        for name, minimum in (
            ("grid_order", 3),
            ("max_rank", 1),
            ("precision_bits", 64),
            ("repeats", 1),
            ("warmups", 0),
            ("solver_max_iterations", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be integer >= {minimum}")
        if self.precision_bits > 4096:
            raise ValueError("precision_bits must be <= 4096")
        for name in ("fit_rcond", "solver_mixing", "solver_tolerance"):
            _positive(getattr(self, name), name, below_one=True)
        _positive(self.timeout_seconds, "timeout_seconds")


def load_config(path: Path) -> ApplicationConfig:
    raw = tomllib.loads(path.read_text())
    schema_version = raw.pop("schema_version", None)
    if type(schema_version) is not int or schema_version != 1:
        raise ValueError("unsupported application benchmark config schema")
    try:
        cases = tuple(ApplicationCase(**case) for case in raw.pop("cases"))
        for name in ("methods", "native_eps_candidates"):
            if name in raw:
                raw[name] = tuple(raw[name])
        return ApplicationConfig(cases=cases, **raw)
    except (KeyError, TypeError) as exc:
        raise ValueError("invalid application benchmark configuration") from exc


def _array(value: Any) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype.kind not in "ifuc" or not np.all(np.isfinite(array)):
        raise ValueError("array must contain finite numeric values")
    return np.asarray(array, dtype=array.dtype.newbyteorder("<"), order="C")


def fingerprint(value: Any) -> str:
    """Hash shape, dtype and complete complex bytes, independent of byte order."""
    array = _array(value)
    header = json.dumps(
        {"shape": array.shape, "dtype": array.dtype.str}, sort_keys=True
    ).encode()
    return hashlib.sha256(header + array.tobytes()).hexdigest()


def array_payload(value: Any) -> dict[str, Any]:
    array = _array(value)
    result = {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "sha256": fingerprint(array),
    }
    if np.iscomplexobj(array):
        result.update(real=array.real.tolist(), imag=array.imag.tolist())
    else:
        result["values"] = array.tolist()
    return result


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return array_payload(value)
    if is_dataclass(value):
        return _json_value(asdict(value))
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _meshes(case: ApplicationCase, config: ApplicationConfig) -> tuple[np.ndarray, ...]:
    training, calibration, final = physical_tau_meshes(
        beta=case.beta, wmax=case.wmax, order=config.grid_order
    )
    calibration_n, final_n = frequency_index_meshes()
    return training, calibration, final, calibration_n, final_n


def _observations(
    case: ApplicationCase, training: np.ndarray
) -> tuple[tuple[Any, ...], list[np.ndarray]]:
    workloads = make_pole_workloads(beta=case.beta, wmax=case.wmax)
    amplitude = matrix_amplitudes()
    values = [
        green_function(w, training, beta=case.beta, wmax=case.wmax)[:, None, None]
        * amplitude
        for w in workloads
    ]
    return workloads, values


def _inputs(case: ApplicationCase, config: ApplicationConfig) -> dict[str, Any]:
    training, calibration, final, calibration_n, final_n = _meshes(case, config)
    _, observations = _observations(case, training)
    inputs = {
        "training_tau": array_payload(training),
        "calibration_tau": array_payload(calibration),
        "final_tau": array_payload(final),
        "calibration_indices": array_payload(calibration_n),
        "final_indices": array_payload(final_n),
        "matrix_amplitude": array_payload(matrix_amplitudes()),
        "observations": [array_payload(values) for values in observations],
    }
    if case.solver:
        inputs.update(
            solver_projection_indices=array_payload(solver_index_mesh()),
            solver_validation_indices=array_payload(solver_validation_index_mesh()),
        )
    return inputs


def _identity(
    method: str, case: ApplicationCase, config: ApplicationConfig, stage: str
) -> dict[str, Any]:
    if method not in METHODS:
        raise ValueError("unknown method")
    return {
        "method": method,
        "case": asdict(case),
        "config": _json_value(asdict(config)),
        "protocol": PROTOCOL,
        "stage": stage,
        "inputs": _inputs(case, config),
        "certification_status": "not_certified_sampled_physical_validation",
    }


def _build_basis(
    method: str, case: ApplicationCase, config: ApplicationConfig, eps: float
) -> Any:
    if method == "kernelgecp":
        return FermionicBasis(
            case.beta,
            case.wmax,
            eps,
            precision="auto",
            precision_bits=config.precision_bits,
            grid_order=config.grid_order,
            max_rank=config.max_rank,
        )
    return build_adapter(method, case.beta, case.wmax, eps)


def _basis_metadata(basis: Any) -> dict[str, Any]:
    if isinstance(basis, FermionicBasis):
        return {
            "package": "kernelgecp",
            "api": "public_FermionicBasis",
            "native_tolerance_meaning": "sampled_kernel_absolute_max",
            "diagnostics": _json_value(basis.diagnostics),
            "certification_status": "not_certified",
        }
    return _json_value(basis.metadata)


def _apply(design: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    result = np.tensordot(design, coefficients, axes=(1, 0))
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite evaluated values")
    return result


def _error(predicted: np.ndarray, truth: np.ndarray, target: float) -> dict[str, Any]:
    if predicted.shape != truth.shape or not np.all(np.isfinite(truth)):
        raise ValueError("prediction/truth shape or finiteness mismatch")
    absolute = float(np.max(np.abs(predicted - truth)))
    scale = float(np.max(np.abs(truth)))
    if not math.isfinite(absolute) or not math.isfinite(scale):
        raise ValueError("physical error arithmetic is nonfinite")
    return {
        "absolute_error": absolute,
        "relative_error": absolute / scale if scale else None,
        "truth_scale": scale,
        "accuracy_passed": absolute <= target,
        "predicted": array_payload(predicted),
        "truth": array_payload(truth),
    }


def _evaluate_workloads(
    case: ApplicationCase,
    workloads: tuple[Any, ...],
    coefficients: list[np.ndarray],
    tau: np.ndarray,
    indices: np.ndarray,
    tau_design: np.ndarray,
    frequency_design: np.ndarray,
) -> list[dict[str, Any]]:
    amplitude = matrix_amplitudes()
    return [
        {
            "name": w.name,
            "physical_poles": list(w.poles),
            "spectral_weights": list(w.weights),
            "time": _error(
                _apply(tau_design, coeff),
                green_function(w, tau, beta=case.beta, wmax=case.wmax)[:, None, None]
                * amplitude,
                case.target,
            ),
            "matsubara": _error(
                _apply(frequency_design, coeff),
                matsubara_green_function(w, indices, beta=case.beta, wmax=case.wmax)[
                    :, None, None
                ]
                * amplitude,
                case.target,
            ),
        }
        for w, coeff in zip(workloads, coefficients, strict=True)
    ]


def _passed(workloads: list[dict[str, Any]]) -> bool:
    return bool(workloads) and all(
        w[domain]["accuracy_passed"]
        for w in workloads
        for domain in ("time", "matsubara")
    )


def run_calibration(
    method: str,
    case: ApplicationCase,
    config: ApplicationConfig,
    *,
    isolated: bool = False,
) -> dict[str, Any]:
    """Attempt declared eps in order; this function never evaluates final truth."""
    start = time.perf_counter()
    result = _identity(method, case, config, "calibration")
    training, calibration, _, calibration_n, _ = _meshes(case, config)
    workloads, observations = _observations(case, training)
    candidates, candidate_measurements = [], []
    measurements: dict[str, Any] = {
        "rss_scope": "calibration_worker_process"
        if isolated
        else "in_process_diagnostic_not_isolated",
        "rss_observation": "high_water_mark_before_output_json_encoding",
    }
    try:
        _, measurements["provider_import_seconds"] = _timed(
            lambda: _preload_provider(method)
        )
        for eps in config.native_eps_candidates:
            candidate: dict[str, Any] = {"native_eps": eps}
            timing: dict[str, Any] = {"native_eps": eps}
            candidate_start = time.perf_counter()
            try:
                basis, timing["construction_seconds"] = _timed(
                    lambda eps=eps: _build_basis(method, case, config, eps)
                )
                design, timing["training_design_seconds"] = _timed(
                    lambda basis=basis: basis.design_tau(training)
                )
                plan, timing["svd_preparation_seconds"] = _timed(
                    lambda design=design: LeastSquaresPlan(
                        design, rcond=config.fit_rcond
                    )
                )
                coefficients, timing["initial_fit_seconds"] = _timed(
                    lambda plan=plan: [plan.fit(values) for values in observations]
                )
                tau_design, timing["calibration_time_design_seconds"] = _timed(
                    lambda basis=basis: basis.design_tau(calibration)
                )
                frequency_design, timing["calibration_frequency_design_seconds"] = (
                    _timed(lambda basis=basis: basis.design_matsubara(calibration_n))
                )
                records = _evaluate_workloads(
                    case,
                    workloads,
                    coefficients,
                    calibration,
                    calibration_n,
                    tau_design,
                    frequency_design,
                )
                candidate.update(
                    status="accuracy_passed" if _passed(records) else "accuracy_failed",
                    rank=int(basis.rank),
                    basis_metadata=_basis_metadata(basis),
                    fit_diagnostics=_json_value(plan.diagnostics),
                    workloads=records,
                )
            except AdapterUnavailableError as exc:
                candidate.update(status="import_unavailable", error=str(exc))
            except Exception as exc:
                candidate.update(
                    status="runtime_failed", error=f"{type(exc).__name__}: {exc}"
                )
            timing["total_seconds"] = time.perf_counter() - candidate_start
            candidates.append(candidate)
            candidate_measurements.append(timing)
            if candidate["status"] in ("accuracy_passed", "import_unavailable"):
                break
        selected = next(
            (c["native_eps"] for c in candidates if c["status"] == "accuracy_passed"),
            None,
        )
        result.update(
            status="calibration_passed"
            if selected is not None
            else "calibration_failed",
            selected_native_eps=selected,
        )
    except AdapterUnavailableError as exc:
        result.update(
            status="import_unavailable", selected_native_eps=None, error=str(exc)
        )
    except Exception as exc:
        result.update(
            status="runtime_failed",
            selected_native_eps=None,
            error=f"{type(exc).__name__}: {exc}",
        )
    measurements.update(
        candidates=candidate_measurements,
        total_seconds=time.perf_counter() - start,
        peak_process_rss_bytes=_peak_rss(),
    )
    result.update(candidates=candidates, measurements=measurements)
    return result


def selection_digest(calibration: dict[str, Any]) -> str:
    numerical = {
        key: value for key, value in calibration.items() if key != "measurements"
    }
    return hashlib.sha256(
        json.dumps(numerical, sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def _selection_receipt(calibration: dict[str, Any]) -> dict[str, Any]:
    """Compact proof summary; selected workers never load calibration arrays."""
    receipt = {
        key: calibration[key]
        for key in ("method", "case", "config", "protocol", "worker_provenance")
    }
    receipt.update(
        selected_native_eps=calibration["selected_native_eps"],
        selection_digest=selection_digest(calibration),
        candidates=[
            {
                "native_eps": candidate["native_eps"],
                "status": candidate["status"],
                **(
                    {
                        "max_absolute_errors": {
                            domain: max(
                                workload[domain]["absolute_error"]
                                for workload in candidate["workloads"]
                            )
                            for domain in ("time", "matsubara")
                        }
                    }
                    if candidate["status"] in {"accuracy_passed", "accuracy_failed"}
                    else {}
                ),
            }
            for candidate in calibration["candidates"]
        ],
    )
    return receipt


def _validate_receipt(
    receipt: Any, method: str, case: ApplicationCase, config: ApplicationConfig
) -> float:
    if not isinstance(receipt, dict) or any(
        receipt.get(key) != value
        for key, value in {
            "method": method,
            "case": asdict(case),
            "config": _json_value(asdict(config)),
            "protocol": PROTOCOL,
        }.items()
    ):
        raise ValueError("selection receipt identity mismatch")
    digest = receipt.get("selection_digest")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(char not in "0123456789abcdef" for char in digest)
    ):
        raise ValueError("selection receipt digest is malformed")
    candidates = receipt.get("candidates")
    if (
        not isinstance(candidates, list)
        or not candidates
        or len(candidates) > len(config.native_eps_candidates)
    ):
        raise ValueError("selection receipt candidates missing")
    selected = None
    for index, candidate in enumerate(candidates):
        if (
            not isinstance(candidate, dict)
            or candidate.get("native_eps") != config.native_eps_candidates[index]
            or candidate.get("status") not in CANDIDATE_STATUSES
            or selected is not None
        ):
            raise ValueError("selection receipt candidate order/status mismatch")
        status = candidate["status"]
        if status in {"accuracy_passed", "accuracy_failed"}:
            errors = candidate.get("max_absolute_errors", {})
            if set(errors) != {"time", "matsubara"} or any(
                isinstance(error, bool)
                or not isinstance(error, Real)
                or not math.isfinite(error)
                or error < 0
                for error in errors.values()
            ):
                raise ValueError("selection receipt errors invalid")
            if (status == "accuracy_passed") != all(
                error <= case.target for error in errors.values()
            ):
                raise ValueError("selection receipt accuracy mismatch")
        if status == "accuracy_passed":
            selected = candidate["native_eps"]
        if status == "import_unavailable":
            raise ValueError("unavailable calibration cannot select a basis")
    if selected is None or receipt.get("selected_native_eps") != selected:
        raise ValueError("selection receipt does not prove selected eps")
    json.dumps(receipt, allow_nan=False)
    return selected


def _decode_array(payload: Any) -> np.ndarray:
    if not isinstance(payload, dict):
        raise ValueError("array payload must be object")
    try:
        dtype = np.dtype(payload["dtype"])
        if dtype.kind == "c":
            array = np.asarray(payload["real"]) + 1j * np.asarray(payload["imag"])
        else:
            array = np.asarray(payload["values"])
        array = _array(np.asarray(array, dtype=dtype))
        if array_payload(array) != payload:
            raise ValueError("array shape/dtype/hash/content mismatch")
        return array
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError("malformed array payload") from exc


def _validate_error(error: Any, truth: np.ndarray, target: float) -> None:
    if not isinstance(error, dict):
        raise ValueError("missing physical error record")
    predicted = _decode_array(error.get("predicted"))
    if _decode_array(error.get("truth")).shape != truth.shape:
        raise ValueError("physical truth shape mismatch")
    if error != _error(predicted, truth, target):
        raise ValueError("physical error/truth/decision mismatch")


def _validate_workloads(
    records: Any, case: ApplicationCase, config: ApplicationConfig, stage: str
) -> None:
    if not isinstance(records, list) or len(records) != 3:
        raise ValueError("missing physical workload records")
    _, calibration, final, calibration_n, final_n = _meshes(case, config)
    tau, indices = (
        (calibration, calibration_n) if stage == "calibration" else (final, final_n)
    )
    for record, workload in zip(
        records, make_pole_workloads(beta=case.beta, wmax=case.wmax), strict=True
    ):
        if (
            record.get("name") != workload.name
            or record.get("physical_poles") != list(workload.poles)
            or record.get("spectral_weights") != list(workload.weights)
        ):
            raise ValueError("physical workload identity mismatch")
        amplitude = matrix_amplitudes()
        _validate_error(
            record.get("time"),
            green_function(workload, tau, beta=case.beta, wmax=case.wmax)[:, None, None]
            * amplitude,
            case.target,
        )
        _validate_error(
            record.get("matsubara"),
            matsubara_green_function(workload, indices, beta=case.beta, wmax=case.wmax)[
                :, None, None
            ]
            * amplitude,
            case.target,
        )


def _selection(calibration: dict[str, Any], config: ApplicationConfig) -> float | None:
    """Fail closed if chosen eps is not the first fully passing declared trial."""
    candidates = calibration.get("candidates")
    if not isinstance(candidates, list) or len(candidates) > len(
        config.native_eps_candidates
    ):
        raise ValueError("invalid calibration candidate sequence")
    selected = None
    for index, candidate in enumerate(candidates):
        if (
            not isinstance(candidate, dict)
            or candidate.get("native_eps") != config.native_eps_candidates[index]
            or candidate.get("status") not in CANDIDATE_STATUSES
        ):
            raise ValueError("invalid calibration candidate identity/status")
        if selected is not None:
            raise ValueError("calibration continued after a passing candidate")
        if candidate["status"] in {"accuracy_passed", "accuracy_failed"}:
            records = candidate.get("workloads")
            _validate_workloads(
                records, ApplicationCase(**calibration["case"]), config, "calibration"
            )
            if not isinstance(records, list) or len(records) != 3:
                raise ValueError("missing calibration workload errors")
            for record in records:
                for domain in ("time", "matsubara"):
                    error = record.get(domain, {})
                    absolute = error.get("absolute_error")
                    if (
                        isinstance(absolute, bool)
                        or not isinstance(absolute, Real)
                        or not math.isfinite(absolute)
                        or absolute < 0
                        or error.get("accuracy_passed")
                        is not (absolute <= calibration["case"]["target"])
                    ):
                        raise ValueError("invalid calibration accuracy decision")
            if (candidate["status"] == "accuracy_passed") != _passed(records):
                raise ValueError("calibration candidate status disagrees with errors")
        if candidate["status"] == "accuracy_passed":
            selected = candidate["native_eps"]
        if candidate["status"] == "import_unavailable" and index != len(candidates) - 1:
            raise ValueError("calibration continued after unavailable import")
    if calibration.get("selected_native_eps") != selected:
        raise ValueError("selected eps is not proven by calibration")
    if calibration.get("status") == "calibration_passed":
        if selected is None:
            raise ValueError("passing calibration requires selected eps")
    elif selected is not None:
        raise ValueError("failed calibration cannot select eps")
    elif calibration.get("status") == "calibration_failed" and (
        not candidates
        or (
            len(candidates) != len(config.native_eps_candidates)
            and candidates[-1]["status"] != "import_unavailable"
        )
    ):
        raise ValueError("failed calibration omitted declared candidates")
    return selected


def _repeated(call: Any, config: ApplicationConfig) -> dict[str, Any]:
    for _ in range(config.warmups):
        call()
    samples = [_timed(call)[1] for _ in range(config.repeats)]
    return {"samples_seconds": samples, "median_seconds": statistics.median(samples)}


def _uncached_design(basis: Any, kind: str, coordinates: np.ndarray) -> np.ndarray:
    if hasattr(basis, "clear_cache"):
        basis.clear_cache()
    return (
        basis.design_tau(coordinates)
        if kind == "time"
        else basis.design_matsubara(coordinates)
    )


def run_case(
    method: str,
    case: ApplicationCase,
    config: ApplicationConfig,
    *,
    calibration: dict[str, Any] | None = None,
    selection: dict[str, Any] | None = None,
    isolated: bool = False,
) -> dict[str, Any]:
    """Run selected measurements; canonical calibration uses another worker."""
    if selection is not None:
        eps = _validate_receipt(selection, method, case, config)
        digest = selection["selection_digest"]
    else:
        if calibration is None:
            calibration = run_calibration(method, case, config)
        _validate_identity(calibration, method, case, config, "calibration")
        eps = _selection(calibration, config)
        digest = selection_digest(calibration)
    result = _identity(method, case, config, "measurement")
    result.update(
        selected_native_eps=eps,
        selection_digest=digest,
        forward_status="not_run",
        overall_success=False,
    )
    measurements: dict[str, Any] = {
        "scope": "nondeterministic_selected_basis_measurements",
        "rss_scope": "fresh_selected_worker_including_imports_native_and_harness"
        if isolated
        else "in_process_diagnostic_not_isolated",
        "rss_observation": "high_water_mark_before_output_json_encoding",
        "construction_scope": (
            "cold_selected_constructor_after_separately_timed_provider_import"
        ),
        "provider_import_scope": "incremental_after_common_harness_imports",
        "uncached_design_scope": (
            "clear_public_instance_caches_then_design_including_cache_clear_cost"
        ),
        "warmups": config.warmups,
        "repeats": config.repeats,
    }
    start = time.perf_counter()
    if eps is None:
        result.update(status="selection_failed", measurements=measurements)
        return result
    training, _, final, _, final_n = _meshes(case, config)
    workloads, observations = _observations(case, training)
    try:
        _, measurements["provider_import_seconds"] = _timed(
            lambda: _preload_provider(method)
        )
        basis, measurements["construction_seconds"] = _timed(
            lambda eps=eps: _build_basis(method, case, config, eps)
        )
        design, measurements["training_design_seconds"] = _timed(
            lambda basis=basis: basis.design_tau(training)
        )
        plan, measurements["svd_preparation_seconds"] = _timed(
            lambda design=design: LeastSquaresPlan(design, rcond=config.fit_rcond)
        )
        coefficients, measurements["initial_fit_seconds"] = _timed(
            lambda plan=plan: [plan.fit(values) for values in observations]
        )
        # Final forward truth and designs first occur after eps selection is frozen.
        tau_design, measurements["final_time_design_seconds"] = _timed(
            lambda: _uncached_design(basis, "time", final)
        )
        frequency_design, measurements["final_frequency_design_seconds"] = _timed(
            lambda: _uncached_design(basis, "matsubara", final_n)
        )
        records = _evaluate_workloads(
            case, workloads, coefficients, final, final_n, tau_design, frequency_design
        )
        for record, values, coefficient in zip(
            records, observations, coefficients, strict=True
        ):
            record["training_absolute_error"] = float(
                np.max(np.abs(_apply(design, coefficient) - values))
            )
            record["coefficients"] = array_payload(coefficient)
        measurements["workloads"] = [
            {
                "name": workload.name,
                "fit_cached_svd": _repeated(
                    lambda values=values: plan.fit(values), config
                ),
                "time_evaluate_cached": _repeated(
                    lambda coefficient=coefficient: _apply(tau_design, coefficient),
                    config,
                ),
                "frequency_evaluate_cached": _repeated(
                    lambda coefficient=coefficient: _apply(
                        frequency_design, coefficient
                    ),
                    config,
                ),
            }
            for workload, values, coefficient in zip(
                workloads, observations, coefficients, strict=True
            )
        ]
        # Design cost is independent of RHS count. Evaluate every fixture in one
        # common batch for the expensive uncached measurement, rather than
        # multiplying repeated reference-arithmetic design work by three.
        coefficient_batch = np.stack(coefficients, axis=1)
        measurements["batched_evaluation_shape"] = list(coefficient_batch.shape)
        measurements["time_design_uncached"] = _repeated(
            lambda: _uncached_design(basis, "time", final), config
        )
        measurements["frequency_design_uncached"] = _repeated(
            lambda: _uncached_design(basis, "matsubara", final_n), config
        )
        measurements["time_evaluate_uncached_batch"] = _repeated(
            lambda: _apply(_uncached_design(basis, "time", final), coefficient_batch),
            config,
        )
        measurements["frequency_evaluate_uncached_batch"] = _repeated(
            lambda: _apply(
                _uncached_design(basis, "matsubara", final_n), coefficient_batch
            ),
            config,
        )
        result.update(
            status="accuracy_passed" if _passed(records) else "accuracy_failed",
            forward_status="accuracy_passed" if _passed(records) else "accuracy_failed",
            overall_success=_passed(records),
            rank=int(basis.rank),
            basis_metadata=_basis_metadata(basis),
            native_tau_nodes=array_payload(basis.tau_nodes),
            fit_diagnostics=_json_value(plan.diagnostics),
            workloads=records,
        )
        setup_keys = (
            "provider_import_seconds",
            "construction_seconds",
            "training_design_seconds",
            "svd_preparation_seconds",
            "initial_fit_seconds",
            "final_time_design_seconds",
            "final_frequency_design_seconds",
        )
        measurements["selected_pipeline_setup_seconds"] = sum(
            measurements[key] for key in setup_keys
        )
        measurements["pipeline_timing_scope"] = (
            "sum_of_provider_import_constructor_training_design_svd_initial_fits_"
            "final_designs_plus_solver_setup_and_loop_when_present"
        )
        if case.solver:
            solver_pipeline_start = time.perf_counter()
            indices = solver_index_mesh()
            solver_validation_indices = solver_validation_index_mesh()
            solver_design, measurements["solver_frequency_design_seconds"] = _timed(
                lambda: _uncached_design(basis, "matsubara", indices)
            )
            solver_fit, measurements["solver_svd_preparation_seconds"] = _timed(
                lambda: LeastSquaresPlan(solver_design, rcond=config.fit_rcond)
            )
            solver, measurements["solver_plan_preparation_seconds"] = _timed(
                lambda: BetheSolverPlan(
                    basis,
                    indices=indices,
                    hopping=case.wmax / 2,
                    rcond=config.fit_rcond,
                    design=solver_design,
                    fit_plan=solver_fit,
                )
            )

            def solve() -> Any:
                return solver.solve(
                    mixing=config.solver_mixing,
                    tolerance=config.solver_tolerance,
                    max_iterations=config.solver_max_iterations,
                )

            solution, measurements["solver_loop_seconds"] = _timed(solve)
            measurements["solver_setup_and_loop_wall_seconds"] = (
                time.perf_counter() - solver_pipeline_start
            )
            (
                solver_validation_design,
                measurements["solver_validation_frequency_design_seconds"],
            ) = _timed(
                lambda: _uncached_design(basis, "matsubara", solver_validation_indices)
            )
            solver_time = _error(
                _apply(tau_design, solution.coefficients),
                semicircle_tau(final, beta=case.beta, hopping=case.wmax / 2),
                case.target,
            )
            solver_frequency = _error(
                _apply(solver_validation_design, solution.coefficients),
                semicircle_matsubara(
                    solver_validation_indices, beta=case.beta, hopping=case.wmax / 2
                ),
                case.target,
            )
            observation_truth = _error(
                solution.frequency_values,
                semicircle_matsubara(indices, beta=case.beta, hopping=case.wmax / 2),
                case.target,
            )
            result["solver"] = {
                "converged": solution.converged,
                "stop_reason": solution.stop_reason,
                "history": _json_value(solution.history),
                "metadata": _json_value(solution.metadata),
                "fit_diagnostics": _json_value(solution.fit_diagnostics),
                "indices": array_payload(solution.indices),
                "coefficients": array_payload(solution.coefficients),
                "frequency_values": array_payload(solution.frequency_values),
                "time": solver_time,
                "matsubara": solver_frequency,
                "validation_indices": array_payload(solver_validation_indices),
                "independent_truth_on_solver_observations": {
                    "heldout": False,
                    "used_for_selection_or_stopping": False,
                    "error": observation_truth,
                },
                "validation_passed": solver_time["accuracy_passed"]
                and solver_frequency["accuracy_passed"]
                and observation_truth["accuracy_passed"],
            }
            solver_passed = solution.converged and result["solver"]["validation_passed"]
            result["solver"]["status"] = (
                "converged_validation_passed"
                if solver_passed
                else "converged_validation_failed"
                if solution.converged
                else "not_converged"
            )
            result["overall_success"] = _passed(records) and solver_passed
            result["status"] = (
                "accuracy_passed" if result["overall_success"] else "accuracy_failed"
            )
            measurements["solver_reuse_cached_setup"] = _repeated(solve, config)
            measurements["solver_setup_seconds"] = (
                measurements["solver_frequency_design_seconds"]
                + measurements["solver_svd_preparation_seconds"]
                + measurements["solver_plan_preparation_seconds"]
                + measurements["solver_validation_frequency_design_seconds"]
            )
            measurements["selected_pipeline_setup_plus_solver_seconds"] = (
                measurements["selected_pipeline_setup_seconds"]
                + measurements["solver_setup_seconds"]
                + measurements["solver_loop_seconds"]
            )
    except AdapterUnavailableError as exc:
        result.update(
            status="import_unavailable", overall_success=False, error=str(exc)
        )
    except Exception as exc:
        result.update(
            status="runtime_failed",
            overall_success=False,
            error=f"{type(exc).__name__}: {exc}",
        )
    measurements.update(
        total_seconds=time.perf_counter() - start, peak_process_rss_bytes=_peak_rss()
    )
    result["measurements"] = measurements
    return result


def _validate_identity(
    record: Any,
    method: str,
    case: ApplicationCase,
    config: ApplicationConfig,
    stage: str,
) -> None:
    if not isinstance(record, dict):
        raise ValueError("worker payload must be object")
    expected = _identity(method, case, config, stage)
    if any(record.get(key) != value for key, value in expected.items()):
        raise ValueError("worker identity/protocol/inputs mismatch")
    statuses = (
        {
            "calibration_passed",
            "calibration_failed",
            "import_unavailable",
            "runtime_failed",
        }
        if stage == "calibration"
        else {
            "accuracy_passed",
            "accuracy_failed",
            "selection_failed",
            "import_unavailable",
            "runtime_failed",
        }
    )
    if record.get("status") not in statuses:
        raise ValueError("worker status mismatch")
    measurements = record.get("measurements")
    if not isinstance(measurements, dict):
        raise ValueError("worker measurements must be object")
    rss = measurements.get("peak_process_rss_bytes")
    if rss is not None and (
        isinstance(rss, bool) or not isinstance(rss, int) or rss <= 0
    ):
        raise ValueError("worker RSS must be positive integer bytes or unavailable")
    if stage == "measurement" and (
        not isinstance(record.get("overall_success"), bool)
        or record.get("forward_status")
        not in {"not_run", "accuracy_passed", "accuracy_failed"}
        or (
            record["status"]
            in {"selection_failed", "import_unavailable", "runtime_failed"}
            and record["overall_success"]
        )
    ):
        raise ValueError("worker overall/forward status mismatch")
    json.dumps(record, allow_nan=False)


def _validate_solver_convergence(
    solver: dict[str, Any],
    case: ApplicationCase,
    config: ApplicationConfig,
    rank: int,
) -> None:
    """Check saved iteration evidence against the actual returned fixed point.

    The last saved defect must agree with a new binary64 equation evaluation
    within 64 eps times the larger of one, the value scale and target scale.
    This comparison allowance never relaxes the convergence threshold: the
    newly computed defect must itself be <= the configured tolerance.
    """
    history = solver.get("history")
    if not isinstance(history, list) or len(history) > config.solver_max_iterations:
        raise ValueError("solver iteration history has invalid length")
    norm_fields = {"step_norm", "defect_norm", "projection_residual_norm"}
    for index, iteration in enumerate(history, 1):
        if (
            not isinstance(iteration, dict)
            or set(iteration) != {"iteration", *norm_fields}
            or isinstance(iteration.get("iteration"), bool)
            or not isinstance(iteration.get("iteration"), int)
            or iteration["iteration"] != index
        ):
            raise ValueError("solver iteration history is not consecutive 1..N")
        if any(
            isinstance(iteration[name], bool)
            or not isinstance(iteration[name], Real)
            or not math.isfinite(iteration[name])
            or iteration[name] < 0
            for name in norm_fields
        ):
            raise ValueError("solver iteration norms must be finite and nonnegative")
    metadata = solver.get("metadata")
    if not isinstance(metadata, dict):
        raise ValueError("solver convergence metadata missing")
    initial_residual = metadata.get("initial_projection_residual_norm")
    if (
        isinstance(initial_residual, bool)
        or not isinstance(initial_residual, Real)
        or not math.isfinite(initial_residual)
        or initial_residual < 0
    ):
        raise ValueError("solver initial projection norm is invalid")
    for name, expected in {
        "beta": case.beta,
        "hopping": case.wmax / 2,
        "mixing": config.solver_mixing,
        "tolerance": config.solver_tolerance,
    }.items():
        value = metadata.get(name)
        if isinstance(value, bool) or not isinstance(value, Real) or value != expected:
            raise ValueError(f"solver metadata {name} does not match configuration")
    for name, expected in {
        "max_iterations": config.solver_max_iterations,
        "iterations": len(history),
    }.items():
        value = metadata.get(name)
        if isinstance(value, bool) or not isinstance(value, int) or value != expected:
            raise ValueError(
                f"solver metadata {name} does not match iteration evidence"
            )
    expected_metadata = {
        "converged": solver.get("converged"),
        "stop_reason": solver.get("stop_reason"),
        "history": history,
        "fit_diagnostics": solver.get("fit_diagnostics"),
        "arithmetic": "numpy_complex128",
        "precision_bits": 53,
        "index_convention": "ordinary_integer_n;wn=(2*n+1)*pi/beta",
        "initial_state": "frequency_projection_of_1_over_z",
        "convergence_scope": "self_consistency_defect_on_supplied_frequency_grid",
        "certification_status": "not_certified",
        "certification_scope": "no_continuous_input_or_rounding_certificate",
    }
    if any(metadata.get(name) != value for name, value in expected_metadata.items()):
        raise ValueError("solver metadata disagrees with convergence evidence")
    diagnostics = solver.get("fit_diagnostics")
    indices = _decode_array(solver.get("indices"))
    values = _decode_array(solver.get("frequency_values"))
    if values.shape != indices.shape or values.dtype != np.dtype("complex128"):
        raise ValueError("solver frequency values must be a complex128 index vector")
    if (
        not isinstance(diagnostics, dict)
        or diagnostics.get("design_shape") != [indices.size, rank]
        or diagnostics.get("rcond") != config.fit_rcond
        or diagnostics.get("factorization") != "cached_thin_svd"
    ):
        raise ValueError("solver SVD diagnostics do not match configuration")
    # Re-evaluate the physical equation without using the solver's saved defect
    # or a truth resolvent. This uses only its final values and supplied modes.
    with np.errstate(all="ignore"):
        z = 1j * (2.0 * indices + 1.0) * (math.pi / case.beta)
        target = 1.0 / (z - (case.wmax / 2) ** 2 * values)
        actual_defect = float(np.max(np.abs(values - target)))
    if not math.isfinite(actual_defect):
        raise ValueError("solver returned a nonfinite actual fixed-point defect")
    converged = solver.get("converged")
    if not isinstance(converged, bool):
        raise ValueError("solver convergence flag must be bool")
    if metadata.get("converged") is not converged:
        raise ValueError("solver metadata convergence flag must match as bool")
    if history:
        comparison_allowance = (
            64
            * np.finfo(float).eps
            * max(1.0, float(np.max(np.abs(values))), float(np.max(np.abs(target))))
        )
        if abs(history[-1]["defect_norm"] - actual_defect) > comparison_allowance:
            raise ValueError("solver saved final defect disagrees with returned values")
        if any(item["defect_norm"] <= config.solver_tolerance for item in history[:-1]):
            raise ValueError("solver history continued after satisfying tolerance")
        if converged != (history[-1]["defect_norm"] <= config.solver_tolerance):
            raise ValueError("solver convergence flag disagrees with final history")
        if converged != (actual_defect <= config.solver_tolerance):
            raise ValueError("solver convergence flag disagrees with actual defect")
    elif converged:
        raise ValueError("solver cannot converge without an accepted iteration")
    stop_reason = solver.get("stop_reason")
    if converged:
        if stop_reason != "self_consistency_tolerance_reached":
            raise ValueError("solver converged with an inconsistent stop reason")
    elif stop_reason == "iteration_budget_exhausted":
        if len(history) != config.solver_max_iterations:
            raise ValueError("solver budget exhaustion has incomplete history")
    elif stop_reason not in {"nonfinite_iteration", "nonfinite_projection"}:
        raise ValueError("solver nonconvergence stop reason is invalid")
    elif len(history) >= config.solver_max_iterations:
        raise ValueError("solver arithmetic failure followed an exhausted budget")


def _validated_worker(
    record: Any,
    method: str,
    case: ApplicationCase,
    config: ApplicationConfig,
    expected: dict[str, Any],
    stage: str,
    calibration: dict[str, Any] | None = None,
) -> dict[str, Any]:
    _validate_identity(record, method, case, config, stage)
    if record.get("worker_provenance") != expected:
        raise RuntimeError("worker source/configuration does not match frozen parent")
    if stage == "calibration":
        _selection(record, config)
    else:
        if (
            calibration is None
            or record.get("selected_native_eps") != _selection(calibration, config)
            or record.get("selection_digest") != selection_digest(calibration)
        ):
            raise ValueError("measurement selection does not match calibration")
        if record["status"] in {"accuracy_passed", "accuracy_failed"}:
            _validate_workloads(record.get("workloads"), case, config, stage)
            if (record.get("forward_status") == "accuracy_passed") != _passed(
                record.get("workloads", [])
            ):
                raise ValueError("measurement accuracy status mismatch")
            overall_success = _passed(record["workloads"])
            if case.solver:
                solver = record.get("solver", {})
                if _decode_array(solver.get("indices")).tolist() != (
                    solver_index_mesh().tolist()
                ):
                    raise ValueError("solver projection indices mismatch")
                if _decode_array(solver.get("validation_indices")).tolist() != (
                    solver_validation_index_mesh().tolist()
                ):
                    raise ValueError("solver validation indices mismatch")
                _validate_solver_convergence(solver, case, config, record["rank"])
                _, _, final, _, final_n = _meshes(case, config)
                _validate_error(
                    solver.get("time"),
                    semicircle_tau(final, beta=case.beta, hopping=case.wmax / 2),
                    case.target,
                )
                _validate_error(
                    solver.get("matsubara"),
                    semicircle_matsubara(
                        solver_validation_index_mesh(),
                        beta=case.beta,
                        hopping=case.wmax / 2,
                    ),
                    case.target,
                )
                observed = solver.get("independent_truth_on_solver_observations", {})
                if (
                    observed.get("heldout") is not False
                    or observed.get("used_for_selection_or_stopping") is not False
                ):
                    raise ValueError("solver observation truth scope mismatch")
                _validate_error(
                    observed.get("error"),
                    semicircle_matsubara(
                        solver_index_mesh(), beta=case.beta, hopping=case.wmax / 2
                    ),
                    case.target,
                )
                if observed["error"]["predicted"] != solver.get("frequency_values"):
                    raise ValueError("solver observation truth uses different values")
                if not isinstance(solver.get("converged"), bool) or solver.get(
                    "validation_passed"
                ) is not (
                    solver["time"]["accuracy_passed"]
                    and solver["matsubara"]["accuracy_passed"]
                    and observed["error"]["accuracy_passed"]
                ):
                    raise ValueError("solver convergence/validation status mismatch")
                overall_success = (
                    overall_success
                    and solver["converged"]
                    and solver["validation_passed"]
                )
                expected_solver_status = (
                    "converged_validation_passed"
                    if solver["converged"] and solver["validation_passed"]
                    else "converged_validation_failed"
                    if solver["converged"]
                    else "not_converged"
                )
                if solver.get("status") != expected_solver_status:
                    raise ValueError("solver status disagrees with acceptance flags")
            if (
                record.get("overall_success") is not overall_success
                or (record["status"] == "accuracy_passed") != overall_success
            ):
                raise ValueError("overall application status mismatch")
    return record


def run_benchmark(config_path: Path, root: Path) -> dict[str, Any]:
    before = source_snapshot(root, config_path)
    config = load_config(config_path)
    if before != source_snapshot(root, config_path):
        raise RuntimeError("source/configuration changed while loading config")
    environment = os.environ.copy()
    environment.update({name: "1" for name in THREAD_VARIABLES})
    records = []
    with tempfile.TemporaryDirectory(prefix="gecp-dv-workers-") as temporary:
        for index, case in enumerate(config.cases):
            for method in config.methods:
                stages: dict[str, Any] = {}
                for stage in ("calibration", "measurement"):
                    output = Path(temporary) / f"{index}-{method}-{stage}.json"
                    command = [
                        sys.executable,
                        "-m",
                        "experiments.dlr_application_benchmark",
                        "--config",
                        str(config_path.resolve()),
                        "--output",
                        str(output),
                        "--worker",
                        method,
                        "--stage",
                        stage,
                        "--case-index",
                        str(index),
                        "--expected-provenance",
                        provenance_digest(before),
                    ]
                    if stage == "measurement":
                        calibration = stages["calibration"]
                        if calibration.get("status") != "calibration_passed":
                            break
                        receipt_path = (
                            Path(temporary) / f"{index}-{method}-selection.json"
                        )
                        receipt_path.write_text(
                            json.dumps(_selection_receipt(calibration), allow_nan=False)
                        )
                        command += [
                            "--selection",
                            str(receipt_path),
                        ]
                    try:
                        completed = subprocess.run(
                            command,
                            cwd=root,
                            env=environment,
                            timeout=config.timeout_seconds,
                            capture_output=True,
                            text=True,
                        )
                        if completed.returncode or not output.is_file():
                            stages[stage] = {
                                "status": "worker_failed",
                                "returncode": completed.returncode,
                                "stderr": completed.stderr[-4000:],
                            }
                        else:
                            try:
                                stages[stage] = _validated_worker(
                                    json.loads(output.read_text()),
                                    method,
                                    case,
                                    config,
                                    before,
                                    stage,
                                    stages.get("calibration"),
                                )
                            except (
                                ValueError,
                                TypeError,
                                KeyError,
                                RuntimeError,
                            ) as exc:
                                stages[stage] = {
                                    "status": "worker_failed",
                                    "error": f"invalid worker payload: {exc}",
                                }
                    except subprocess.TimeoutExpired:
                        stages[stage] = {
                            "status": "timeout",
                            "timeout_seconds": config.timeout_seconds,
                        }
                measurement = stages.get("measurement")
                status = (
                    measurement["status"]
                    if measurement
                    else stages["calibration"]["status"]
                )
                records.append(
                    {"method": method, "case": asdict(case), "status": status, **stages}
                )
    if before != source_snapshot(root, config_path):
        raise RuntimeError("source/configuration changed during benchmark")
    return {
        "schema_version": 1,
        "protocol": PROTOCOL,
        "claim_level": (
            "matched_sampled_physical_application_evidence_not_speedup_or_certification"
        ),
        "provenance": before,
        "config": _json_value(asdict(config)),
        "sampling_roles": {
            "calibration_indices": array_payload(frequency_index_meshes()[0]),
            "forward_final_indices": array_payload(frequency_index_meshes()[1]),
            "solver_projection_indices": array_payload(solver_index_mesh()),
            "solver_final_indices": array_payload(solver_validation_index_mesh()),
            "final_tau_interval_fractions": ["2/7", "5/7"],
            "solver_observation_truth": "independent_post_loop_not_heldout",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "versions": _versions(),
            "worker_threads": {name: "1" for name in THREAD_VARIABLES},
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("experiments/dlr_application_benchmark.toml"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", choices=METHODS, help=argparse.SUPPRESS)
    parser.add_argument(
        "--stage", choices=("calibration", "measurement"), help=argparse.SUPPRESS
    )
    parser.add_argument("--case-index", type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument("--expected-provenance", help=argparse.SUPPRESS)
    parser.add_argument("--selection", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.worker:
        before = source_snapshot(root, args.config)
        if provenance_digest(before) != args.expected_provenance:
            raise RuntimeError("worker did not start from expected frozen source")
        config = load_config(args.config)
        if args.stage == "calibration":
            result = run_calibration(
                args.worker, config.cases[args.case_index], config, isolated=True
            )
        elif args.stage == "measurement" and args.selection:
            selection = json.loads(args.selection.read_text())
            _validate_receipt(
                selection, args.worker, config.cases[args.case_index], config
            )
            if selection.get("worker_provenance") != before:
                raise RuntimeError("selection receipt provenance mismatch")
            result = run_case(
                args.worker,
                config.cases[args.case_index],
                config,
                selection=selection,
                isolated=True,
            )
        else:
            raise ValueError(
                "worker requires stage and measurement requires selection receipt"
            )
        if before != source_snapshot(root, args.config):
            raise RuntimeError("worker source/configuration changed during execution")
        result["worker_provenance"] = before
        result["measurements"]["worker_process_id"] = os.getpid()
        result["measurements"]["peak_process_rss_bytes"] = _peak_rss()
    else:
        result = run_benchmark(args.config.resolve(), root)
    encoded = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(encoded)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
