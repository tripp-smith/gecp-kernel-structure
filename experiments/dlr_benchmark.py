"""Matched physical Green-function basis benchmarks, not performance guarantees.

Run each case/method in a fresh process. Numerical fields are separate from
machine-dependent measurements; all bases use the same cached SVD fitter.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
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
from collections.abc import Callable
from dataclasses import asdict, dataclass
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from experiments.dlr_adapters import (
    AdapterUnavailableError,
    build_adapter,
)
from experiments.dlr_workloads import (
    common_tau_meshes,
    green_function,
    make_pole_workloads,
)

METHODS = ("kernelgecp", "pydlr", "sparse_ir")
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
)


@dataclass(frozen=True)
class BenchmarkCase:
    beta: float
    wmax: float
    tolerance: float = 1e-8

    def __post_init__(self) -> None:
        for value in (self.beta, self.wmax, self.tolerance):
            if (
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, Real)
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError("case parameters must be finite and positive")
        if not math.isfinite(self.beta * self.wmax) or self.beta * self.wmax < 1:
            raise ValueError("benchmark requires finite beta*wmax >= 1")
        if self.tolerance >= 1:
            raise ValueError("tolerance must be less than one")


@dataclass(frozen=True)
class BenchmarkConfig:
    cases: tuple[BenchmarkCase, ...]
    methods: tuple[str, ...] = METHODS
    grid_order: int = 12
    max_rank: int = 128
    repeats: int = 5
    warmups: int = 1
    fit_rcond: float = 1e-14
    timeout_seconds: float = 120.0

    def __post_init__(self) -> None:
        if not self.cases or not all(isinstance(c, BenchmarkCase) for c in self.cases):
            raise ValueError("cases must be nonempty BenchmarkCase records")
        if (
            not self.methods
            or not all(isinstance(method, str) for method in self.methods)
            or len(set(self.methods)) != len(self.methods)
            or any(method not in METHODS for method in self.methods)
        ):
            raise ValueError("methods must be distinct supported names")
        for name, minimum in (
            ("grid_order", 3),
            ("max_rank", 1),
            ("repeats", 1),
            ("warmups", 0),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if not math.isfinite(self.fit_rcond) or not 0 < self.fit_rcond < 1:
            raise ValueError("fit_rcond must lie in (0, 1)")
        if (
            isinstance(self.timeout_seconds, (bool, np.bool_))
            or not isinstance(self.timeout_seconds, Real)
            or not math.isfinite(self.timeout_seconds)
            or self.timeout_seconds <= 0
        ):
            raise ValueError("timeout_seconds must be finite and positive")


def load_config(path: Path) -> BenchmarkConfig:
    raw = tomllib.loads(path.read_text())
    if raw.pop("schema_version", None) != 1:
        raise ValueError("unsupported benchmark config schema")
    try:
        cases = tuple(BenchmarkCase(**case) for case in raw.pop("cases"))
        if "methods" in raw:
            raw["methods"] = tuple(raw["methods"])
        return BenchmarkConfig(cases=cases, **raw)
    except (KeyError, TypeError) as exc:
        raise ValueError("invalid benchmark configuration") from exc


def fingerprint(array: np.ndarray) -> str:
    values = np.asarray(array, dtype="<f8", order="C")
    shape = json.dumps(values.shape).encode()
    return hashlib.sha256(shape + values.tobytes()).hexdigest()


class CachedFitPlan:
    """One SVD for every RHS; no explicit inverse or repeated factorization."""

    def __init__(self, matrix: np.ndarray, rcond: float) -> None:
        if np.iscomplexobj(matrix):
            raise ValueError("design must be real; complex RHS values are supported")
        self.matrix = np.array(matrix, dtype=float, copy=True)
        if (
            self.matrix.ndim != 2
            or min(self.matrix.shape, default=0) == 0
            or not np.all(np.isfinite(self.matrix))
        ):
            raise ValueError("design must be a finite nonempty matrix")
        if not math.isfinite(rcond) or not 0 < rcond < 1:
            raise ValueError("rcond must lie in (0, 1)")
        u, singular, vh = np.linalg.svd(self.matrix, full_matrices=False)
        retained = singular > rcond * singular[0]
        self.u = u[:, retained]
        self.singular = singular[retained]
        self.vh = vh[retained, :]
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            condition = float(singular[0] / singular[-1])
        self.diagnostics = {
            "shape": list(self.matrix.shape),
            "singular_values": singular.tolist(),
            "effective_rank": int(np.sum(retained)),
            "rcond": rcond,
            "precision_bits": 53,
            "coefficient_coordinates": "method_native_not_basis_invariant",
            "condition_number": condition if math.isfinite(condition) else None,
            "condition_status": (
                "estimated"
                if math.isfinite(condition)
                else "singular_or_outside_binary64_range"
            ),
        }

    def fit(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if (
            values.ndim not in (1, 2)
            or values.shape[0] != self.matrix.shape[0]
            or not np.all(np.isfinite(values))
        ):
            raise ValueError("RHS must be finite with matching observation axis")
        projected = self.u.T @ values
        divisor = self.singular if values.ndim == 1 else self.singular[:, None]
        coefficients = np.asarray(self.vh.T @ (projected / divisor))
        if not np.all(np.isfinite(coefficients)):
            raise ValueError("nonfinite fitted coefficients")
        return coefficients


class GECPBasis:
    """Benchmark-local cardinal basis; existing public GECP API is unchanged."""

    def __init__(self, case: BenchmarkCase, config: BenchmarkConfig) -> None:
        from kernelgecp import FermionicKernel, GECPConfig, gecp

        self.beta = case.beta
        self.kernel = FermionicKernel(case.beta * case.wmax)
        result = gecp(
            self.kernel,
            config=GECPConfig(
                pivot="grid",
                precision_bits=53,
                tol=case.tolerance,
                max_rank=config.max_rank,
                grid_order=config.grid_order,
            ),
        )
        if not result.rank:
            raise ValueError("GECP returned an empty basis")
        self.rank = result.rank
        self.tau_nodes = result.t_nodes * case.beta
        self.omega_nodes = result.omega_nodes.copy()
        core = self.kernel(result.t_nodes[:, None], self.omega_nodes[None, :])
        self._core_transpose_lu = lu_factor(core.T)
        self.metadata = {
            "package": "kernelgecp",
            "basis": "selected_cross_cardinal_coordinates",
            "native_tolerance_meaning": "sampled_kernel_absolute_max",
            "native_converged": result.converged,
            "native_stop_reason": result.stop_reason,
            "precision_bits": 53,
            "certification_status": "not_certified",
        }

    def design_matrix(self, tau: np.ndarray) -> np.ndarray:
        columns = self.kernel(tau[:, None] / self.beta, self.omega_nodes[None, :])
        return np.asarray(lu_solve(self._core_transpose_lu, columns.T).T)


def _build_basis(method: str, case: BenchmarkCase, config: BenchmarkConfig) -> Any:
    if method == "kernelgecp":
        return GECPBasis(case, config)
    return build_adapter(method, case.beta, case.wmax, case.tolerance)


def _preload_provider(method: str) -> None:
    try:
        importlib.import_module(method)
    except ImportError as exc:
        raise AdapterUnavailableError(f"{method} import unavailable: {exc}") from exc


def _timed(call: Callable[[], Any]) -> tuple[Any, float]:
    start = time.perf_counter()
    result = call()
    return result, time.perf_counter() - start


def _repeated(call: Callable[[], Any], config: BenchmarkConfig) -> dict[str, Any]:
    for _ in range(config.warmups):
        call()
    samples = [_timed(call)[1] for _ in range(config.repeats)]
    return {"samples_seconds": samples, "median_seconds": statistics.median(samples)}


def _versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _peak_rss() -> int | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(peak if sys.platform == "darwin" else peak * 1024)
    except ImportError:
        return None


def run_case(
    method: str,
    case: BenchmarkCase,
    config: BenchmarkConfig,
    *,
    isolated: bool = False,
) -> dict[str, Any]:
    """Evaluate fixed holdouts only after all construction/fit decisions."""
    if method not in METHODS:
        raise ValueError("unknown method")
    training, validation = common_tau_meshes(
        beta=case.beta, wmax=case.wmax, order=config.grid_order
    )
    workloads = make_pole_workloads(beta=case.beta, wmax=case.wmax)
    observations = [
        green_function(w, training, beta=case.beta, wmax=case.wmax) for w in workloads
    ]
    result: dict[str, Any] = {
        "method": method,
        "grid_order": config.grid_order,
        "case": asdict(case),
        "lambda": case.beta * case.wmax,
        "protocol": "shared_observations_cached_svd_v1",
        "certification_status": "not_certified_sampled_validation",
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "versions": _versions(),
            "threads": {name: os.environ.get(name) for name in THREAD_VARIABLES},
            "numpy_build": np.__config__.show(mode="dicts"),
        },
        "inputs": {
            "training_tau": training.tolist(),
            "validation_tau": validation.tolist(),
            "training_hash": fingerprint(training),
            "validation_hash": fingerprint(validation),
            "observation_hashes": [fingerprint(values) for values in observations],
        },
    }
    measurements: dict[str, Any] = {
        "scope": "nondeterministic_pilot_not_canonical",
        "construction_scope": "actual_adapter_constructor_including_package_defaults",
        "provider_import_scope": "incremental_after_common_harness_not_cold_startup",
        "rss_scope": (
            "fresh_worker_process_including_imports_and_native_libraries"
            if isolated
            else "in_process_diagnostic_not_isolated"
        ),
        "warmups": config.warmups,
        "repeats": config.repeats,
    }
    try:
        _, measurements["provider_import_seconds"] = _timed(
            lambda: _preload_provider(method)
        )
        basis, measurements["construction_seconds"] = _timed(
            lambda: _build_basis(method, case, config)
        )
        design, measurements["training_design_seconds"] = _timed(
            lambda: basis.design_matrix(training)
        )
        plan, measurements["fit_plan_seconds"] = _timed(
            lambda: CachedFitPlan(design, config.fit_rcond)
        )
        fitted = [plan.fit(values) for values in observations]
        # No holdout values are evaluated until the basis and all fits are frozen.
        validation_design, measurements["validation_design_seconds"] = _timed(
            lambda: basis.design_matrix(validation)
        )
        if not np.all(np.isfinite(validation_design)):
            raise ValueError("nonfinite validation design")
        records = []
        measurements["workloads"] = []
        for workload, values, coeffs in zip(
            workloads, observations, fitted, strict=True
        ):
            truth = green_function(workload, validation, beta=case.beta, wmax=case.wmax)
            predicted = validation_design @ coeffs
            if not np.all(np.isfinite(predicted)):
                raise ValueError("nonfinite predicted values")
            error = float(np.max(np.abs(predicted - truth)))
            scale = float(np.max(np.abs(truth)))
            records.append(
                {
                    "name": workload.name,
                    "training_absolute_error": float(
                        np.max(np.abs(design @ coeffs - values))
                    ),
                    "absolute_error": error,
                    "relative_error": error / scale,
                    "accuracy_passed": error <= case.tolerance,
                    "truth_hash": fingerprint(truth),
                    "predicted_hash": fingerprint(predicted),
                    "predicted_values": predicted.tolist(),
                    "truth_values": truth.tolist(),
                    "physical_poles": list(workload.poles),
                    "spectral_weights": list(workload.weights),
                }
            )
            measurements["workloads"].append(
                {
                    "name": workload.name,
                    "fit": _repeated(lambda values=values: plan.fit(values), config),
                    "evaluate_uncached": _repeated(
                        lambda coeffs=coeffs: basis.design_matrix(validation) @ coeffs,
                        config,
                    ),
                    "evaluate_cached": _repeated(
                        lambda coeffs=coeffs: validation_design @ coeffs, config
                    ),
                }
            )
        result.update(
            status=(
                "accuracy_passed"
                if all(record["accuracy_passed"] for record in records)
                else "accuracy_failed"
            ),
            rank=int(basis.rank),
            native_tau_nodes=np.asarray(basis.tau_nodes).tolist(),
            basis_metadata=basis.metadata,
            fit_diagnostics=plan.diagnostics,
            workloads=records,
        )
    except AdapterUnavailableError as exc:
        result.update(status="import_unavailable", error=str(exc))
    except Exception as exc:
        result.update(status="runtime_failed", error=f"{type(exc).__name__}: {exc}")
    measurements["peak_process_rss_bytes"] = _peak_rss()
    result["measurements"] = measurements
    return result


def source_snapshot(root: Path, config_path: Path) -> dict[str, Any]:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=root, check=True, capture_output=True, text=True
        ).stdout.strip()

    if git("status", "--porcelain", "--untracked-files=no"):
        raise ValueError("benchmark requires a clean tracked source revision")
    if git(
        "ls-files", "--others", "--exclude-standard", "--", "experiments", "kernelgecp"
    ):
        raise ValueError("untracked experiment/library code prevents a frozen run")
    paths = [*root.glob("kernelgecp/*.py"), *root.glob("experiments/*.py")]
    paths += [root / "pyproject.toml", root / "uv.lock"]
    return {
        "git_commit": git("rev-parse", "HEAD"),
        "git_tree": git("rev-parse", "HEAD^{tree}"),
        "source_sha256": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)
        },
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
    }


def provenance_digest(snapshot: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(snapshot, sort_keys=True).encode()).hexdigest()


def _validated_worker(
    record: Any,
    method: str,
    case: BenchmarkCase,
    expected: dict[str, Any],
    grid_order: int,
) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise ValueError("worker payload must be an object")
    if record.get("worker_provenance") != expected:
        raise RuntimeError("worker source/configuration does not match frozen parent")
    if (
        record.get("method") != method
        or record.get("case") != asdict(case)
        or record.get("grid_order") != grid_order
        or record.get("protocol") != "shared_observations_cached_svd_v1"
        or record.get("status")
        not in {
            "accuracy_passed",
            "accuracy_failed",
            "import_unavailable",
            "runtime_failed",
        }
    ):
        raise ValueError("worker identity/protocol/status mismatch")
    json.dumps(record, allow_nan=False)
    training, validation = common_tau_meshes(
        beta=case.beta, wmax=case.wmax, order=grid_order
    )
    inputs = record.get("inputs", {})
    if not isinstance(inputs, dict):
        raise ValueError("worker inputs must be an object")
    if (
        inputs.get("training_tau") != training.tolist()
        or inputs.get("validation_tau") != validation.tolist()
        or inputs.get("training_hash") != fingerprint(training)
        or inputs.get("validation_hash") != fingerprint(validation)
        or inputs.get("observation_hashes")
        != [
            fingerprint(green_function(w, training, beta=case.beta, wmax=case.wmax))
            for w in make_pole_workloads(beta=case.beta, wmax=case.wmax)
        ]
    ):
        raise ValueError("worker sampling coordinates do not match the request")
    return record


def run_benchmark(config_path: Path, root: Path) -> dict[str, Any]:
    before = source_snapshot(root, config_path)
    config = load_config(config_path)
    if before != source_snapshot(root, config_path):
        raise RuntimeError("source/configuration changed while loading configuration")
    records: list[dict[str, Any]] = []
    environment = os.environ.copy()
    environment.update({name: "1" for name in THREAD_VARIABLES})
    with tempfile.TemporaryDirectory(prefix="gecp-dlr-workers-") as temporary:
        for index, case in enumerate(config.cases):
            for method in config.methods:
                output = Path(temporary) / f"{index}-{method}.json"
                command = [
                    sys.executable,
                    "-m",
                    "experiments.dlr_benchmark",
                    "--config",
                    str(config_path.resolve()),
                    "--output",
                    str(output),
                    "--worker",
                    method,
                    "--case-index",
                    str(index),
                    "--expected-provenance",
                    provenance_digest(before),
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
                        records.append(
                            {
                                "method": method,
                                "case": asdict(case),
                                "status": "worker_failed",
                                "returncode": completed.returncode,
                                "stderr": completed.stderr[-4000:],
                            }
                        )
                    else:
                        try:
                            records.append(
                                _validated_worker(
                                    json.loads(output.read_text()),
                                    method,
                                    case,
                                    before,
                                    config.grid_order,
                                )
                            )
                        except (ValueError, TypeError) as exc:
                            records.append(
                                {
                                    "method": method,
                                    "case": asdict(case),
                                    "status": "worker_failed",
                                    "error": f"invalid worker payload: {exc}",
                                }
                            )
                except subprocess.TimeoutExpired:
                    records.append(
                        {"method": method, "case": asdict(case), "status": "timeout"}
                    )
    if before != source_snapshot(root, config_path):
        raise RuntimeError("source/configuration changed during benchmark")
    return {
        "schema_version": 1,
        "claim_level": "sampled_forward_pilot_not_speedup_or_certification",
        "provenance": before,
        "config": asdict(config),
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("experiments/dlr_benchmark.toml")
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", choices=METHODS, help=argparse.SUPPRESS)
    parser.add_argument("--case-index", type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument("--expected-provenance", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.worker:
        before = source_snapshot(root, args.config)
        if provenance_digest(before) != args.expected_provenance:
            raise RuntimeError("worker did not start from expected frozen source")
        config = load_config(args.config)
        result = run_case(
            args.worker, config.cases[args.case_index], config, isolated=True
        )
        if before != source_snapshot(root, args.config):
            raise RuntimeError("worker source/configuration changed during execution")
        result["worker_provenance"] = before
    else:
        result = run_benchmark(args.config.resolve(), root)
    encoded = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(encoded)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
