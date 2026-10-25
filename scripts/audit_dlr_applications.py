"""Independently audit saved physical DLR applications at 95 decimal digits.

No benchmark or library helper is imported. Pole formulas use exact lifts of
the actual binary64 inputs; the solver oracle uses a causal quadratic root and
selected independent energy integrals. Honest scientific failures are retained.
An inconsistent artifact, hash, source identity, shape, or flag fails the audit.
This is sampled numerical evidence, not an interval or rounding certificate.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

DECIMAL_PRECISION = 95
PROTOCOL = "matched_physical_time_frequency_svd_v1"
THREAD_VARIABLES = {
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "RAYON_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
}
FAILURE_STATUSES = {
    "calibration_failed",
    "import_unavailable",
    "runtime_failed",
    "worker_failed",
    "timeout",
    "selection_failed",
}


class AuditError(ValueError):
    """The saved evidence is inconsistent or cannot be independently checked."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def lift_binary64(value: Any) -> Any:
    """Lift the actual finite binary64 value, not its short decimal spelling."""
    require(not isinstance(value, (bool, np.bool_)), "boolean is not a real input")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise AuditError("invalid binary64 real input") from exc
    require(math.isfinite(number), "nonfinite binary64 input")
    return mp.mpf(number)


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON object key")
        result[key] = value
    return result


def _nonfinite_constant(value: str) -> None:
    raise AuditError(f"nonfinite JSON constant: {value}")


def read_artifact(path: Path, expected_sha256: str) -> tuple[dict[str, Any], str, str]:
    """Verify container identity, decompress if needed, and hash JSON payload."""
    require(
        isinstance(expected_sha256, str)
        and re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is not None,
        "expected SHA256 must be 64 lowercase hexadecimal characters",
    )
    container = path.read_bytes()
    container_hash = hashlib.sha256(container).hexdigest()
    require(container_hash == expected_sha256, "artifact container hash mismatch")
    try:
        payload = (
            gzip.decompress(container) if path.name.endswith(".json.gz") else container
        )
        artifact = json.loads(
            payload, object_pairs_hook=_object, parse_constant=_nonfinite_constant
        )
    except (OSError, EOFError, UnicodeError, json.JSONDecodeError) as exc:
        raise AuditError("malformed JSON/gzip artifact") from exc
    require(isinstance(artifact, dict), "artifact must be a JSON object")
    return artifact, container_hash, hashlib.sha256(payload).hexdigest()


def decode_array(payload: Any) -> np.ndarray:
    """Check the complete dtype/shape/complex-content fingerprint independently."""
    require(isinstance(payload, dict), "array payload must be an object")
    shape = payload.get("shape")
    require(
        isinstance(shape, list)
        and all(type(size) is int and size >= 0 for size in shape),
        "invalid array shape",
    )
    try:
        dtype = np.dtype(payload["dtype"])
        require(dtype.kind in "ifuc", "invalid numeric array dtype")
        require(dtype.str == payload["dtype"], "noncanonical array dtype")
        if dtype.kind == "c":
            real = np.asarray(payload["real"])
            imag = np.asarray(payload["imag"])
            require(
                real.dtype.kind in "ifu" and imag.dtype.kind in "ifu",
                "complex components must contain real numbers",
            )
            require(
                real.shape == imag.shape == tuple(shape),
                "complex component shape mismatch",
            )
            array = real.astype(float) + 1j * imag.astype(float)
        else:
            array = np.asarray(payload["values"])
            require(array.dtype.kind in "ifu", "array must contain real numbers")
            if dtype.kind in "iu":
                require(
                    array.dtype.kind in "iu", "integer payload contains nonintegers"
                )
        array = np.asarray(array, dtype=dtype.newbyteorder("<"), order="C")
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise AuditError("malformed array payload") from exc
    require(list(array.shape) == shape, "array shape mismatch")
    require(np.all(np.isfinite(array)), "nonfinite array payload")
    header = json.dumps(
        {"shape": array.shape, "dtype": array.dtype.str}, sort_keys=True
    ).encode()
    digest = hashlib.sha256(header + array.tobytes()).hexdigest()
    require(digest == payload.get("sha256"), "array content hash mismatch")
    return array


def _all_payloads(value: Any) -> Any:
    if isinstance(value, dict):
        if {"dtype", "shape", "sha256"}.issubset(value):
            yield value
        else:
            for item in value.values():
                yield from _all_payloads(item)
    elif isinstance(value, list):
        for item in value:
            yield from _all_payloads(item)


def _values(payload: Any, shape: tuple[int, ...] | None = None) -> list[Any]:
    array = decode_array(payload)
    require(shape is None or array.shape == shape, "physical array shape mismatch")
    if np.iscomplexobj(array):
        return [
            mp.mpc(lift_binary64(value.real), lift_binary64(value.imag))
            for value in array.flat
        ]
    return [lift_binary64(value) for value in array.flat]


def _indices(payload: Any) -> list[int]:
    array = decode_array(payload)
    require(
        array.ndim == 1 and array.size > 0 and array.dtype.kind in "iu",
        "frequency indices must be a nonempty integer vector",
    )
    result = [int(value) for value in array]
    require(len(set(result)) == len(result), "duplicate frequency indices")
    require(
        all(abs(value) < 2**51 for value in result), "frequency index out of bounds"
    )
    return result


def _times(payload: Any, beta: Any) -> list[Any]:
    array = decode_array(payload)
    require(
        array.ndim == 1 and array.size > 0 and array.dtype.kind == "f",
        "time coordinates must be a nonempty floating vector",
    )
    times = [lift_binary64(value) for value in array]
    require(all(0 <= value <= beta for value in times), "time outside [0,beta]")
    require(
        all(a < b for a, b in zip(times, times[1:], strict=False)),
        "time order/duplicates",
    )
    return times


def _number(value: Any) -> str:
    return str(mp.nstr(value, 32))


def verify_error(
    record: Any,
    truth: list[Any],
    target: Any,
    label: str,
    report: list[dict[str, Any]],
    *,
    expected_shape: tuple[int, ...],
) -> bool:
    """Recompute physical errors; an accurately recorded failure returns False."""
    require(isinstance(record, dict), label + ": missing error record")
    predicted = _values(record.get("predicted"), expected_shape)
    stored_truth = _values(record.get("truth"), expected_shape)
    require(
        len(predicted) == len(truth) == len(stored_truth) > 0,
        label + ": physical point count mismatch",
    )
    absolute = max(abs(a - b) for a, b in zip(predicted, truth, strict=True))
    scale = max(abs(value) for value in truth)
    truth_difference = max(abs(a - b) for a, b in zip(stored_truth, truth, strict=True))
    error_difference = abs(absolute - lift_binary64(record.get("absolute_error")))
    allowance = mp.mpf("2e-14") * max(mp.mpf(1), scale)
    require(truth_difference < allowance, label + ": independent truth disagrees")
    require(error_difference < allowance, label + ": independent error disagrees")
    require(
        abs(scale - lift_binary64(record.get("truth_scale"))) < allowance,
        label + ": truth scale mismatch",
    )
    passed = bool(absolute <= target)
    require(
        type(record.get("accuracy_passed")) is bool
        and record["accuracy_passed"] == passed,
        label + ": accuracy flag mismatch",
    )
    relative = absolute / scale if scale else None
    if relative is None:
        require(record.get("relative_error") is None, label + ": zero truth scale")
    else:
        require(
            abs(relative - lift_binary64(record.get("relative_error")))
            < mp.mpf("2e-13") * max(mp.mpf(1), relative),
            label + ": relative error mismatch",
        )
    report.append(
        {
            "label": label,
            "absolute_error_mp": _number(absolute),
            "stored_error": record["absolute_error"],
            "relative_error_mp": None if relative is None else _number(relative),
            "truth_scale_mp": _number(scale),
            "maximum_truth_difference": _number(truth_difference),
            "absolute_error_difference": _number(error_difference),
            "accuracy_passed": passed,
            "point_count": len(truth),
        }
    )
    return passed


def expected_workloads(case: dict[str, Any]) -> list[dict[str, Any]]:
    """Derive fixed fixture identity from the case, independently of its payload."""
    wmax = float(case["wmax"])
    require(math.isfinite(wmax) and wmax > 0, "invalid physical cutoff")
    return [
        {"name": "central_pole", "physical_poles": [0.0], "spectral_weights": [1.0]},
        {
            "name": "asymmetric_satellites",
            "physical_poles": [wmax * value for value in (-0.75, 0.08, 0.55)],
            "spectral_weights": [0.2, 0.5, 0.3],
        },
        {
            "name": "gapped",
            "physical_poles": [wmax * value for value in (-0.8, -0.3, 0.35, 0.9)],
            "spectral_weights": [0.15, 0.35, 0.3, 0.2],
        },
    ]


def _fixture_identity(record: Any, fixture: dict[str, Any]) -> None:
    require(
        isinstance(record, dict)
        and all(record.get(key) == value for key, value in fixture.items()),
        "physical workload identity mismatch",
    )


def pole_time(tau: Any, poles: list[Any], weights: list[Any], beta: Any) -> Any:
    return -sum(
        weight * mp.exp(-tau * pole) / (1 + mp.exp(-beta * pole))
        for pole, weight in zip(poles, weights, strict=True)
    )


def pole_frequency(index: int, poles: list[Any], weights: list[Any], beta: Any) -> Any:
    z = mp.j * (2 * index + 1) * mp.pi / beta
    return sum(weight / (z - pole) for pole, weight in zip(poles, weights, strict=True))


def semicircle_resolvent(index: int, beta: Any, hopping: Any) -> Any:
    """Causal quadratic root, separate from the benchmark's real hypot formula."""
    z = mp.j * (2 * index + 1) * mp.pi / beta
    root = mp.sqrt(z * z - 4 * hopping * hopping)
    if mp.im(root) * mp.im(z) < 0:
        root = -root
    result = 2 / (z + root)
    tolerance = mp.power(10, -(mp.mp.dps - 10)) * max(1, abs(result))
    require(
        abs(result - 1 / (z - hopping**2 * result)) < tolerance,
        "semicircle root defect",
    )
    require(mp.im(result) * mp.im(z) < 0, "noncausal semicircle root")
    return result


def semicircle_energy_integral(tau: Any, beta: Any, hopping: Any) -> Any:
    """Direct MP integration against the spectral density in physical energy."""
    edge = 2 * hopping

    def integrand(energy: Any) -> Any:
        density = mp.sqrt(max(mp.mpf(0), edge * edge - energy * energy)) / (
            2 * mp.pi * hopping * hopping
        )
        return -density * mp.exp(-tau * energy) / (1 + mp.exp(-beta * energy))

    return mp.quad(integrand, [-edge, -hopping, 0, hopping, edge])


def _git(root: Path, *arguments: str) -> bytes:
    try:
        return subprocess.check_output(
            ["git", *arguments], cwd=root, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as exc:
        raise AuditError("cannot read expected Git source/blob") from exc


def provenance_checks(
    artifact: dict[str, Any], root: Path, expected_source: str
) -> dict[str, Any]:
    """Check immutable blobs, source-era config, and current audit environment."""
    require(
        re.fullmatch(r"[0-9a-f]{40}", expected_source) is not None,
        "source must be the complete expected Git commit SHA",
    )
    provenance = artifact["provenance"]
    commit = provenance["git_commit"]
    require(commit == expected_source, "unexpected source revision")
    require(
        _git(root, "rev-parse", commit + "^{tree}").decode().strip()
        == provenance["git_tree"],
        "Git tree mismatch",
    )
    tree_paths = (
        _git(root, "ls-tree", "-r", "--name-only", commit).decode().splitlines()
    )
    expected_paths = {
        path
        for path in tree_paths
        if path in {"pyproject.toml", "uv.lock"}
        or (
            path.count("/") == 1
            and path.endswith(".py")
            and path.split("/")[0] in {"kernelgecp", "experiments"}
        )
    }
    require(
        set(provenance["source_sha256"]) == expected_paths,
        "source fingerprint omitted or added frozen files",
    )
    for path, expected_hash in provenance["source_sha256"].items():
        require(
            hashlib.sha256(_git(root, "show", commit + ":" + path)).hexdigest()
            == expected_hash,
            "Git blob hash mismatch: " + path,
        )
    config_path = "experiments/dlr_application_benchmark.toml"
    config_bytes = _git(root, "show", commit + ":" + config_path)
    require(
        hashlib.sha256(config_bytes).hexdigest() == provenance["config_sha256"],
        "configuration blob hash mismatch",
    )
    config = tomllib.loads(config_bytes.decode())
    version = config.pop("schema_version")
    require(type(version) is int and version == 1, "invalid config schema")
    for case in config["cases"]:
        case.setdefault("solver", False)
    # Compare against that revision's TOML. In particular, do not inject the
    # later calibration_fraction into cohort A, which predates that field.
    require(config == artifact["config"], "source-era TOML/config identity mismatch")
    recorded_versions = artifact["environment"]["versions"]
    current_versions = {
        name: importlib.metadata.version(name) for name in recorded_versions
    }
    require(current_versions == recorded_versions, "installed versions differ")
    environment = artifact["environment"]
    for key, actual in (
        ("python", platform.python_version()),
        ("machine", platform.machine()),
        ("platform", platform.platform()),
    ):
        require(environment[key] == actual, key + " differs from recorded environment")
    controls = environment["worker_threads"]
    require(
        set(controls) == THREAD_VARIABLES and set(controls.values()) == {"1"},
        "thread controls mismatch",
    )
    benchmark_source = _git(
        root, "show", commit + ":experiments/dlr_application_benchmark.py"
    ).decode()
    require(
        'environment.update({name: "1" for name in THREAD_VARIABLES})'
        in benchmark_source,
        "missing explicit worker environment controls",
    )
    return {
        "source_revision": commit,
        "source_blob_hashes_verified": len(expected_paths),
        "config_sha256": provenance["config_sha256"],
        "installed_versions": current_versions,
        "worker_threads_recorded_and_source_propagation_verified": controls,
        "historical_os_thread_count": "not_recorded_not_independently_recoverable",
        "python_executable": sys.executable,
    }


def _stage_identity(
    stage: dict[str, Any], artifact: dict[str, Any], method: str, case: dict[str, Any]
) -> None:
    if stage["status"] in {"worker_failed", "timeout"}:
        return
    require(
        stage.get("worker_provenance") == artifact["provenance"], "worker provenance"
    )
    require(stage.get("config") == artifact["config"], "worker configuration")
    require(stage.get("method") == method and stage.get("case") == case, "worker case")
    require(stage.get("protocol") == PROTOCOL, "worker protocol mismatch")


def _selection(calibration: dict[str, Any], config: dict[str, Any]) -> Any:
    if calibration["status"] != "calibration_passed":
        require(calibration["status"] in FAILURE_STATUSES, "unknown calibration status")
        return None
    candidates = calibration.get("candidates")
    require(
        isinstance(candidates, list) and candidates, "missing calibration candidates"
    )
    declared = config["native_eps_candidates"]
    require(len(candidates) <= len(declared), "too many calibration candidates")
    fraction = config.get("calibration_fraction", 1.0)
    threshold = calibration["case"]["target"] * fraction
    if "calibration_absolute_target" in calibration:
        require(
            calibration["calibration_absolute_target"] == threshold,
            "calibration target",
        )
    selected = None
    for index, candidate in enumerate(candidates):
        require(
            candidate.get("native_eps") == declared[index], "candidate order mismatch"
        )
        require(selected is None, "calibration continued after passing eps")
        status = candidate.get("status")
        require(
            status in {"accuracy_passed", "accuracy_failed", "runtime_failed"},
            "invalid successful-calibration trial status",
        )
        if status in {"accuracy_passed", "accuracy_failed"}:
            workloads = candidate.get("workloads")
            fixtures = expected_workloads(calibration["case"])
            require(
                isinstance(workloads, list) and len(workloads) == 3,
                "calibration fixtures",
            )
            flags = []
            for workload, fixture in zip(workloads, fixtures, strict=True):
                _fixture_identity(workload, fixture)
                for domain in ("time", "matsubara"):
                    error = workload[domain]
                    absolute = float(error["absolute_error"])
                    require(
                        math.isfinite(absolute) and absolute >= 0, "calibration error"
                    )
                    passed = absolute <= threshold
                    require(
                        error.get("accuracy_passed") is passed, "calibration error flag"
                    )
                    flags.append(passed)
            require((status == "accuracy_passed") == all(flags), "calibration status")
        if status == "accuracy_passed":
            selected = candidate["native_eps"]
    require(
        selected is not None and calibration.get("selected_native_eps") == selected,
        "selection is not the first passing declared eps",
    )
    return selected


def _solver_checks(
    solver: dict[str, Any],
    case: dict[str, Any],
    times: list[Any],
    indices: list[int],
    label: str,
    frequency_report: list[dict[str, Any]],
    time_report: list[dict[str, Any]],
    energy_cache: dict[Any, Any],
) -> dict[str, Any]:
    beta, target = lift_binary64(case["beta"]), lift_binary64(case["target"])
    hopping = lift_binary64(case["wmax"]) / 2
    observed = solver["independent_truth_on_solver_observations"]
    require(
        observed.get("heldout") is False
        and observed.get("used_for_selection_or_stopping") is False,
        "solver observation-truth scope mismatch",
    )
    nodes = _indices(solver["indices"])
    holdout = _indices(solver["validation_indices"])
    require(
        not set(nodes).intersection(holdout), "solver holdout overlaps observations"
    )
    require(set(holdout).issubset(indices), "solver holdout differs from final mesh")
    frequency_flags = []
    for name, grid, error in (
        ("heldout_frequency", holdout, solver["matsubara"]),
        ("observed_frequency", nodes, observed["error"]),
    ):
        truth = [semicircle_resolvent(n, beta, hopping) for n in grid]
        frequency_flags.append(
            verify_error(
                error,
                truth,
                target,
                label + "/" + name,
                frequency_report,
                expected_shape=(len(grid),),
            )
        )
    require(
        observed["error"]["predicted"] == solver["frequency_values"],
        "solver returned values",
    )
    returned = _values(solver["frequency_values"], (len(nodes),))
    defect = max(
        abs(value - 1 / (mp.j * (2 * n + 1) * mp.pi / beta - hopping**2 * value))
        for n, value in zip(nodes, returned, strict=True)
    )
    converged = solver.get("converged")
    require(type(converged) is bool, "solver convergence must be bool")
    history = solver.get("history")
    require(isinstance(history, list), "solver history missing")
    metadata = solver["metadata"]
    tolerance = lift_binary64(metadata["tolerance"])
    if history:
        require(
            abs(defect - lift_binary64(history[-1]["defect_norm"])) < mp.mpf("2e-14"),
            "saved solver defect differs from independent equation",
        )
    require(not converged or bool(history), "convergence without an iteration")
    require(
        converged == bool(history and defect <= tolerance), "solver defect/convergence"
    )
    predicted_time = _values(solver["time"]["predicted"], (len(times),))
    stored_time = _values(solver["time"]["truth"], (len(times),))
    selected = sorted(
        {0, 1, len(times) // 5, len(times) // 2, len(times) - 2, len(times) - 1}
    )
    selected = [index for index in selected if 0 <= index < len(times)]
    for index in selected:
        tau = times[index]
        key = (tau, beta, hopping)
        if key not in energy_cache:
            energy_cache[key] = semicircle_energy_integral(tau, beta, hopping)
        truth = energy_cache[key]
        difference = abs(truth - stored_time[index])
        require(difference < mp.mpf("2e-14"), "semicircle energy integral disagreement")
        time_report.append(
            {
                "label": label,
                "index": index,
                "tau_binary64_exact": _number(tau),
                "energy_integral_95dps": _number(truth),
                "stored_truth_difference": _number(difference),
                "prediction_error_mp": _number(abs(predicted_time[index] - truth)),
                "accuracy_passed": bool(abs(predicted_time[index] - truth) <= target),
            }
        )
    # All saved time samples can be checked for internal error/flag consistency;
    # only the selected energy integrals above are independent time-truth checks.
    saved_time = max(
        abs(a - b) for a, b in zip(predicted_time, stored_time, strict=True)
    )
    saved_time_passed = bool(saved_time <= target)
    require(solver["time"]["accuracy_passed"] is saved_time_passed, "solver time flag")
    require(
        abs(saved_time - lift_binary64(solver["time"]["absolute_error"]))
        < mp.mpf("2e-14"),
        "solver time stored error",
    )
    validated = saved_time_passed and all(frequency_flags)
    require(solver.get("validation_passed") is validated, "solver validation flag")
    status = (
        "converged_validation_passed"
        if converged and validated
        else "converged_validation_failed"
        if converged
        else "not_converged"
    )
    require(solver.get("status") == status, "solver status mismatch")
    return {
        "converged": converged,
        "validation_passed": validated,
        "actual_defect_mp": _number(defect),
        "scope": "equation_defect_on_saved_observed_indices",
    }


def audit_artifact(
    path: Path, root: Path, source: str, expected_sha256: str
) -> dict[str, Any]:
    artifact, container_hash, payload_hash = read_artifact(path, expected_sha256)
    require(
        type(artifact.get("schema_version")) is int and artifact["schema_version"] == 1,
        "unsupported artifact schema",
    )
    require(artifact.get("protocol") == PROTOCOL, "unsupported artifact protocol")
    provenance = provenance_checks(artifact, root, source)
    payloads = list(_all_payloads(artifact))
    for payload in payloads:
        decode_array(payload)
    config = artifact["config"]
    expected = [
        (method, case) for case in config["cases"] for method in config["methods"]
    ]
    records = artifact.get("records")
    require(
        isinstance(records, list) and len(records) == len(expected),
        "missing/extra cases",
    )
    forward: list[dict[str, Any]] = []
    frequencies: list[dict[str, Any]] = []
    selected_time: list[dict[str, Any]] = []
    flags = []
    cache: dict[Any, Any] = {}
    energy_cache: dict[Any, Any] = {}
    shared_inputs = {}
    with mp.workdps(DECIMAL_PRECISION):
        for record, (method, case) in zip(records, expected, strict=True):
            require(
                record.get("method") == method and record.get("case") == case,
                "case order/identity",
            )
            calibration = record.get("calibration")
            require(isinstance(calibration, dict), "missing calibration result")
            _stage_identity(calibration, artifact, method, case)
            selected_eps = _selection(calibration, config)
            measurement = record.get("measurement")
            if measurement is not None:
                require(isinstance(measurement, dict), "invalid measurement result")
                _stage_identity(measurement, artifact, method, case)
            if measurement is None or measurement.get("status") in FAILURE_STATUSES:
                failed_stage = measurement if measurement is not None else calibration
                require(
                    selected_eps is None or measurement is not None,
                    "passing calibration omitted measurement",
                )
                require(
                    failed_stage["status"] in FAILURE_STATUSES,
                    "unexplained missing measurement",
                )
                require(
                    not failed_stage.get("overall_success", False),
                    "failed worker marked successful",
                )
                require(record["status"] == failed_stage["status"], "failed top status")
                flags.append(
                    {
                        "method": method,
                        "case": case,
                        "overall_success": False,
                        "audit_status": "retained_failure_without_physical_predictions",
                        "scientific_status": record["status"],
                    }
                )
                continue
            require(selected_eps is not None, "measurement without passing calibration")
            require(
                measurement.get("selected_native_eps") == selected_eps,
                "measurement eps mismatch",
            )
            require(
                measurement["status"] in {"accuracy_passed", "accuracy_failed"},
                "measurement status",
            )
            inputs = measurement["inputs"]
            case_key = json.dumps(case, sort_keys=True)
            if case_key in shared_inputs:
                require(
                    inputs == shared_inputs[case_key],
                    "methods received different observations/queries",
                )
            shared_inputs[case_key] = inputs
            beta, target = lift_binary64(case["beta"]), lift_binary64(case["target"])
            require(beta > 0 and target > 0, "invalid physical case")
            times = _times(inputs["final_tau"], beta)
            training = _times(inputs["training_tau"], beta)
            calibration_times = _times(inputs["calibration_tau"], beta)
            require(
                not set(times).intersection(training + calibration_times),
                "final time overlap",
            )
            require(
                not set(training).intersection(calibration_times),
                "calibration time overlap",
            )
            indices = _indices(inputs["final_indices"])
            require(
                not set(indices).intersection(_indices(inputs["calibration_indices"])),
                "final frequency overlap",
            )
            amplitude = _values(inputs["matrix_amplitude"], (2, 2))
            require(
                amplitude
                == [
                    mp.mpc(1),
                    mp.mpc("0.25", "0.5"),
                    mp.mpc("0.25", "-0.5"),
                    mp.mpc("0.75"),
                ],
                "fixed matrix amplitude identity mismatch",
            )
            workloads = measurement.get("workloads")
            require(
                isinstance(workloads, list) and len(workloads) == 3,
                "missing physical fixtures",
            )
            passed = True
            for workload, fixture in zip(
                workloads, expected_workloads(case), strict=True
            ):
                _fixture_identity(workload, fixture)
                poles = [lift_binary64(value) for value in fixture["physical_poles"]]
                weights = [
                    lift_binary64(value) for value in fixture["spectral_weights"]
                ]
                label = (
                    f"{method}/beta{case['beta']}/wmax{case['wmax']}/{fixture['name']}"
                )
                key = (beta, tuple(poles), tuple(weights), tuple(times), tuple(indices))
                if key not in cache:
                    cache[key] = (
                        [pole_time(t, poles, weights, beta) for t in times],
                        [pole_frequency(n, poles, weights, beta) for n in indices],
                    )
                for domain, scalar_truth in zip(
                    ("time", "matsubara"), cache[key], strict=True
                ):
                    truth = [
                        scalar * amp for scalar in scalar_truth for amp in amplitude
                    ]
                    domain_passed = verify_error(
                        workload[domain],
                        truth,
                        target,
                        label + "/" + domain,
                        forward,
                        expected_shape=(len(scalar_truth), 2, 2),
                    )
                    passed = passed and domain_passed
            require(
                (measurement.get("forward_status") == "accuracy_passed") == passed,
                "forward status",
            )
            overall = passed
            solver_summary = None
            if case["solver"]:
                solver_summary = _solver_checks(
                    measurement["solver"],
                    case,
                    times,
                    indices,
                    method + "/solver",
                    frequencies,
                    selected_time,
                    energy_cache,
                )
                overall = (
                    overall
                    and solver_summary["converged"]
                    and solver_summary["validation_passed"]
                )
            require(
                measurement.get("overall_success") is overall, "overall flag mismatch"
            )
            require(
                (record.get("status") == "accuracy_passed") == overall,
                "top record status",
            )
            require(record["status"] == measurement["status"], "top measurement status")
            flags.append(
                {
                    "method": method,
                    "case": case,
                    "overall_success": overall,
                    "audit_status": "physical_predictions_independently_audited",
                    "scientific_status": record["status"],
                    "solver": solver_summary,
                }
            )
        truth_max = max(
            (mp.mpf(row["maximum_truth_difference"]) for row in forward),
            default=mp.mpf(0),
        )
        error_max = max(
            (mp.mpf(row["absolute_error_difference"]) for row in forward),
            default=mp.mpf(0),
        )
        return {
            "audit_scope": (
                "independent_95_digit_physical_truth_with_exact_binary64_inputs"
            ),
            "certification_status": "not_interval_or_rounding_certified",
            "artifact_container_sha256": container_hash,
            "artifact_uncompressed_sha256": payload_hash,
            "audit_script_sha256": hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
            "decimal_precision": DECIMAL_PRECISION,
            "provenance": provenance,
            "array_payloads_hash_verified": len(payloads),
            "records_attempted": len(records),
            "forward_workloads_audited": len(forward) // 2,
            "forward_domain_checks": forward,
            "solver_frequency_checks": frequencies,
            "selected_semicircle_energy_integrals": selected_time,
            "overall_flags": flags,
            "overall_pass_count": sum(row["overall_success"] for row in flags),
            "retained_worker_failures": [
                row for row in flags if row["audit_status"].startswith("retained")
            ],
            "failed_forward_domain_checks": [
                row for row in forward if not row["accuracy_passed"]
            ],
            "maximum_forward_stored_truth_difference": _number(truth_max),
            "maximum_forward_reported_error_difference": _number(error_max),
            "limitations": [
                "JSON numbers are lifted as actual binary64, not short decimal reals.",
                "Predictions are checked against physical truth, not rebuilt designs.",
                "All solver frequencies and six selected time integrals are checked.",
                "Thread settings are checked; actual OS thread counts are unavailable.",
                "Retained worker failures have no physical predictions to audit.",
                "The parallel proof route has not delivered continuous certification.",
            ],
        }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--source", required=True, help="full expected source commit SHA"
    )
    parser.add_argument(
        "--expected-sha256", required=True, help="input container SHA256"
    )
    args = parser.parse_args(argv)
    try:
        result = audit_artifact(
            args.input, args.root, args.source, args.expected_sha256
        )
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    except (AuditError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f"DLR application audit failed: {exc}", file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "records_attempted",
                    "forward_workloads_audited",
                    "overall_pass_count",
                    "failed_forward_domain_checks",
                    "retained_worker_failures",
                    "artifact_container_sha256",
                    "artifact_uncompressed_sha256",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
