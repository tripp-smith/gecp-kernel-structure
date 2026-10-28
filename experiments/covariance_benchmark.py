"""Frozen, attributed UCI covariance campaign; calibration never reads test rows.

Run from the repository root as ``python -m experiments.covariance_benchmark``.
Each validation case must run in a fresh process. The real-data gates supplement,
and never replace, independent small-matrix tests and the full root verifier.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import resource
import subprocess
import sys
import time
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator, cg
from scipy.spatial.distance import cdist

from kernelgecp import KernelCholeskyFactor

ROOT = Path(__file__).resolve().parents[1]
DATA_URL = (
    "https://archive.ics.uci.edu/static/public/501/"
    "beijing+multi+site+air+quality+data.zip"
)
OUTER_SHA256 = "b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4"
CSV_SHA256 = "a94fbcfc71708b6ffdc033360163d91efb7f675e496f6b5860929aa96b96351b"
FEATURES = ("TEMP", "PRES", "DEWP", "RAIN", "WSPM")
SIZES = (512, 1024, 2048, 4096, 8192, 16384)
RANKS = (32, 64, 128, 256, 512)
TRAIN_END, CAL_END, TEST_END = 16384, 16896, 17408
SOURCE_FILES = (
    "kernelgecp/kernel_cholesky.py",
    "experiments/covariance_benchmark.py",
    "scripts/audit_covariance.py",
    "environments/gp/uv.lock",
    "environments/gp/pyproject.toml",
)
LOCKED_VERSIONS = {"numpy": "2.2.6", "scipy": "1.18.1"}
STD_DEFINITION = {
    "kind": "posterior_predictive_observation_standard_deviation",
    "normalization": "training_target_standard_deviation",
    "includes_observation_noise": True,
    "noise_parameter_is_variance": True,
}
TIMING_SCOPES = {
    "application_workflow": (
        "normalization_factorization_prediction_excludes_imports_and_data_loading"
    ),
    "full_workflow": (
        "compute_dense_CG_validation_NPZ_writing_excludes_imports_loading_final_receipt"
    ),
    "process_wall_seconds": (
        "parent_measured_fresh_child_launch_through_complete_exit_in_manifest"
    ),
}
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def provenance() -> dict[str, Any]:
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    source_hashes = {name: sha256((ROOT / name).read_bytes()) for name in SOURCE_FILES}
    tracked_dirty = bool(
        subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=ROOT,
            text=True,
        ).strip()
    )
    tracked = (
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", *SOURCE_FILES],
            cwd=ROOT,
            capture_output=True,
        ).returncode
        == 0
    )
    return {
        "source_revision": revision,
        "source_hashes": source_hashes,
        "source_clean_and_tracked": tracked and not tracked_dirty,
        "source_binding_scope": "executed_scientific_sources_and_isolated_gp_runtime",
        "development_lock_sha256_context_only": sha256((ROOT / "uv.lock").read_bytes()),
        "versions": {
            name: importlib.metadata.version(name) for name in ("numpy", "scipy")
        },
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "scientific_runtime_lock": "environments/gp/uv.lock",
        "platform": platform.platform(),
        "threads": {key: os.environ.get(key) for key in THREADS},
    }


def load_data(path: Path) -> dict[str, Any]:
    """Verify both archive and station identities before using any observation."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(DATA_URL, timeout=60) as response:
            downloaded = response.read(20 * 1024 * 1024 + 1)
        if len(downloaded) > 20 * 1024 * 1024:
            raise ValueError("unexpectedly large UCI archive")
        if sha256(downloaded) != OUTER_SHA256:
            raise ValueError("downloaded UCI archive hash mismatch")
        path.write_bytes(downloaded)
    raw = path.read_bytes()
    if sha256(raw) != OUTER_SHA256:
        raise ValueError("UCI archive hash mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        nested = archive.read("PRSA2017_Data_20130301-20170228.zip")
    with zipfile.ZipFile(io.BytesIO(nested)) as archive:
        names = [
            name
            for name in archive.namelist()
            if name.endswith("/PRSA_Data_Aotizhongxin_20130301-20170228.csv")
        ]
        if len(names) != 1:
            raise ValueError("expected exactly one Aotizhongxin station CSV")
        station = archive.read(names[0])
    if sha256(station) != CSV_SHA256:
        raise ValueError("UCI station CSV hash mismatch")
    return parse_station(station)


def parse_station(station: bytes) -> dict[str, Any]:
    """Chronological complete cases in just the declared five features and target."""
    records, timestamps, ordinals = [], [], []
    previous: datetime | None = None
    missing = 0
    reader = csv.DictReader(io.StringIO(station.decode("utf-8-sig")))
    for row in reader:
        if row["station"] != "Aotizhongxin":
            raise ValueError("wrong station")
        stamp = datetime(*(int(row[key]) for key in ("year", "month", "day", "hour")))
        if previous is not None and stamp <= previous:
            raise ValueError("observations must be strictly chronological")
        previous = stamp
        entries = [row[key] for key in (*FEATURES, "PM2.5")]
        if any(value == "NA" or value == "" for value in entries):
            missing += 1
            continue
        values = np.array([float(value) for value in entries])
        if not np.all(np.isfinite(values)):
            raise ValueError("nonfinite observation")
        records.append(values)
        timestamps.append(stamp.isoformat())
        ordinals.append(int(row["No"]))
    return {
        "values": np.asarray(records),
        "timestamps": timestamps,
        "ordinals": np.array(ordinals, dtype=np.int64),
        "excluded_missing": missing,
        "outer_sha256": OUTER_SHA256,
        "csv_sha256": CSV_SHA256,
    }


def split_data(data: dict[str, Any], size: int, partition: str) -> dict[str, Any]:
    if (
        type(size) is not int
        or size not in SIZES
        or partition not in ("calibration", "test")
    ):
        raise ValueError("unsupported size or partition")
    if len(data["values"]) < TEST_END:
        raise ValueError("insufficient complete cases for the frozen splits")
    # Each smaller cohort uses only its own training prefix for every scale.
    train = data["values"][:size]
    mean, scale = train[:, :5].mean(axis=0), train[:, :5].std(axis=0)
    y_mean, y_scale = float(train[:, 5].mean()), float(train[:, 5].std())
    if np.any(scale <= 0) or y_scale <= 0:
        raise ValueError("constant training feature/target cannot be standardized")
    start, end = (
        (TRAIN_END, CAL_END) if partition == "calibration" else (CAL_END, TEST_END)
    )
    query = data["values"][start:end]
    return {
        "x": (train[:, :5] - mean) / scale,
        "y": (train[:, 5] - y_mean) / y_scale,
        "query_x": (query[:, :5] - mean) / scale,
        "query_y": (query[:, 5] - y_mean) / y_scale,
        "normalization": {
            "features_mean": mean.tolist(),
            "features_std": scale.tolist(),
            "target_mean": y_mean,
            "target_std": y_scale,
            "fit_rows": size,
        },
        "train_ordinal_sha256": sha256(data["ordinals"][:size].tobytes()),
        "query_ordinal_sha256": sha256(data["ordinals"][start:end].tobytes()),
        "train_first": data["timestamps"][0],
        "train_last": data["timestamps"][size - 1],
        "query_first": data["timestamps"][start],
        "query_last": data["timestamps"][end - 1],
    }


def rbf(left: np.ndarray, right: np.ndarray, lengthscale: float) -> np.ndarray:
    """2D blocked callers; no (n,m,d) broadcast or hidden full-kernel allocation."""
    if not np.isfinite(lengthscale) or lengthscale <= 0:
        raise ValueError("positive finite RBF lengthscale required")
    return np.exp(-0.5 * cdist(left / lengthscale, right / lengthscale, "sqeuclidean"))


def original_operator(
    x: np.ndarray, lengthscale: float, noise: float
) -> LinearOperator:
    def multiply(rhs):
        array = np.asarray(rhs)
        result = np.empty_like(array, dtype=np.float64)
        for start in range(0, len(x), 128):
            end = min(start + 128, len(x))
            result[start:end] = rbf(x[start:end], x, lengthscale) @ array
            result[start:end] += noise * array[start:end]
        return result

    return LinearOperator(
        (len(x), len(x)), matvec=multiply, matmat=multiply, dtype=np.float64
    )


def dense_reference(
    split: dict[str, Any], lengthscale: float, noise: float
) -> tuple[np.ndarray, np.ndarray]:
    x, query, y = split["x"], split["query_x"], split["y"]
    if len(x) > 2048:
        raise ValueError("dense reference is bounded to 2048 training observations")
    matrix = rbf(x, x, lengthscale)
    matrix[np.diag_indices(len(x))] += noise
    chol = cho_factor(matrix, lower=True)
    cross = rbf(query, x, lengthscale)
    mean = cross @ cho_solve(chol, y)
    variance = 1 - np.sum(cross * cho_solve(chol, cross.T).T, axis=1)
    if np.min(variance) < -1e-10:
        raise RuntimeError("negative dense latent variance")
    return mean, np.sqrt(np.maximum(variance, 0) + noise)


def predictions(
    factor: KernelCholeskyFactor,
    split: dict[str, Any],
    lengthscale: float,
    noise: float,
):
    cross = rbf(split["query_x"], split["x"][list(factor.indices)], lengthscale)
    extension = factor.extend(cross, diagonal=np.ones(len(cross)))
    shifted = factor.shifted(noise)
    solved = shifted.solve_with_diagnostics(split["y"])
    mean = extension.features @ (factor.factor.T @ solved.solution)
    # Evaluate the small posterior covariance of features; never form the
    # dense query-by-training cross or a query-by-query covariance.
    small = np.eye(factor.rank) + (factor.factor.T @ factor.factor) / noise
    if factor.rank:
        chol = cho_factor(small, lower=True)
        variance = np.sum(
            extension.features * cho_solve(chol, extension.features.T).T, axis=1
        )
    else:
        variance = np.zeros(len(cross))
    # Residual prior uncertainty is retained, not silently discarded.
    variance += extension.residual_diagonal
    if not np.all(np.isfinite(variance)) or np.min(variance) < -1e-10:
        raise RuntimeError("invalid approximate predictive variance")
    std = np.sqrt(np.maximum(variance, 0) + noise)
    return mean, std, shifted, solved


def accuracy(mean, std, expected_mean, expected_std) -> dict[str, float]:
    return {
        "mean_rms": float(np.sqrt(np.mean((mean - expected_mean) ** 2))),
        "mean_max": float(np.max(np.abs(mean - expected_mean))),
        "std_rms": float(np.sqrt(np.mean((std - expected_std) ** 2))),
        "std_max": float(np.max(np.abs(std - expected_std))),
    }


def calibration(data: dict[str, Any], sizes: tuple[int, ...]) -> dict[str, Any]:
    """Choose calibration likelihood subject to a declared representation budget."""
    if provenance()["versions"] != LOCKED_VERSIONS:
        raise ValueError("calibration requires the same locked GP NumPy/SciPy runtime")
    split = split_data(data, 2048, "calibration")
    candidates = []
    for lengthscale in (0.7, 1.5, 3.0):
        for noise in (0.05, 0.2, 0.8):
            mean, std = dense_reference(split, lengthscale, noise)
            nll = float(
                np.mean(np.log(std) + 0.5 * ((split["query_y"] - mean) / std) ** 2)
            )
            factor = KernelCholeskyFactor.from_columns(
                np.ones(2048),
                lambda j, x=split["x"], ell=lengthscale: rbf(x, x[[j]], ell)[:, 0],
                max_rank=512,
                atol=1e-12,
            )
            approximate_mean, approximate_std, _, _ = predictions(
                factor, split, lengthscale, noise
            )
            metrics = accuracy(approximate_mean, approximate_std, mean, std)
            candidates.append(
                {
                    "lengthscale": lengthscale,
                    "noise": noise,
                    "nll": nll,
                    "rank512_accuracy": metrics,
                    "rank512_admissible": (
                        max(metrics["mean_rms"], metrics["std_rms"]) <= 1e-4
                        and max(metrics["mean_max"], metrics["std_max"]) <= 1e-3
                    ),
                }
            )
            print(
                "calibration_candidate " + json.dumps(candidates[-1], sort_keys=True),
                flush=True,
            )
    admissible = [
        candidate for candidate in candidates if candidate["rank512_admissible"]
    ]
    common = {
        "schema": "kernelgecp.covariance-calibration.v1",
        "provenance": provenance(),
        "dataset": dataset_identity(data),
        "calibration_only": True,
        "test_evaluated": False,
        "sizes": list(sizes),
        "ranks": list(RANKS),
        "standard_deviation_definition": STD_DEFINITION,
        "timing_scopes": TIMING_SCOPES,
        "stochastic_operations": False,
        "random_seed": None,
        "hyperparameter_candidates": candidates,
        "hyperparameter_policy": (
            "best_calibration_nll_subject_to_tenfold_dense_margin_at_n2048_rank512"
        ),
        "rank_policy": (
            "smallest_tenfold_margin_dense_calibration_rank_at_n<=2048_else_512"
        ),
    }
    if not admissible:
        return {
            **common,
            "status": "inconclusive",
            "failure": (
                "all kernel candidates exceeded the declared rank512 "
                "representation budget"
            ),
        }
    selected = min(
        admissible, key=lambda item: (item["nll"], item["lengthscale"], item["noise"])
    )
    lengthscale, noise = selected["lengthscale"], selected["noise"]
    rows, chosen, exhausted = [], {}, []
    for size in sizes:
        split = split_data(data, size, "calibration")
        expected = dense_reference(split, lengthscale, noise) if size <= 2048 else None
        passing = []
        for rank in RANKS:
            start = time.perf_counter()
            factor = KernelCholeskyFactor.from_columns(
                np.ones(size),
                lambda j, x=split["x"]: rbf(x, x[[j]], lengthscale)[:, 0],
                max_rank=rank,
                atol=1e-12,
            )
            mean, std, _, solve = predictions(factor, split, lengthscale, noise)
            metrics = accuracy(mean, std, *expected) if expected is not None else None
            # Tenfold safety margin is fixed before untouched validation.
            eligible = (
                metrics is not None
                and max(metrics["mean_rms"], metrics["std_rms"]) <= 1e-4
                and max(metrics["mean_max"], metrics["std_max"]) <= 1e-3
            )
            if eligible:
                passing.append(rank)
            rows.append(
                {
                    "size": size,
                    "rank_budget": rank,
                    "actual_rank": factor.rank,
                    "accuracy": metrics,
                    "eligible_small_dense_margin": eligible,
                    "predictive_rmse": float(
                        np.sqrt(np.mean((mean - split["query_y"]) ** 2))
                    ),
                    "solve_relative_residual": solve.relative_residual,
                    "seconds_pilot_only": time.perf_counter() - start,
                }
            )
            print(
                "calibration_rank " + json.dumps(rows[-1], sort_keys=True), flush=True
            )
        if size <= 2048 and not passing:
            exhausted.append(size)
        # Beyond the bounded dense oracle, use the largest declared budget.
        # This is explicit fixed-resource policy, not a large-n accuracy claim.
        chosen[str(size)] = min(passing) if passing else (512 if size > 2048 else None)
    return {
        **common,
        "status": "calibrated" if not exhausted else "inconclusive",
        "exhausted_sizes": exhausted,
        "lengthscale": lengthscale,
        "noise": noise,
        "selected_ranks": chosen,
        "rank_calibration": rows,
    }


def dataset_identity(data: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_url": DATA_URL,
        "outer_sha256": OUTER_SHA256,
        "csv_sha256": CSV_SHA256,
        "doi": "10.24432/C5RK5G",
        "attribution": (
            "Chen, S. (2017). Beijing Multi-Site Air Quality. "
            "UCI Machine Learning Repository."
        ),
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "station": "Aotizhongxin",
        "features": list(FEATURES),
        "target": "PM2.5",
        "feature_units": ["degree_Celsius", "hPa", "degree_Celsius", "mm", "m/s"],
        "target_units": "microgram/m^3",
        "complete_cases": len(data["values"]),
        "excluded_missing": data["excluded_missing"],
        "split_complete_case_offsets": [0, TRAIN_END, CAL_END, TEST_END],
        "small_size_policy": "prefix_of_training_only_normalization_per_prefix",
    }


def validate_frozen(frozen: dict[str, Any], *, canonical: bool) -> None:
    if (
        frozen.get("schema") != "kernelgecp.covariance-calibration.v1"
        or frozen.get("status") != "calibrated"
    ):
        raise ValueError("invalid frozen calibration")
    if (
        frozen.get("test_evaluated") is not False
        or frozen.get("calibration_only") is not True
    ):
        raise ValueError("calibration must exclude test observations")
    if (
        frozen["dataset"]["outer_sha256"] != OUTER_SHA256
        or frozen["dataset"]["csv_sha256"] != CSV_SHA256
    ):
        raise ValueError("wrong frozen dataset")
    current = provenance()
    if current["source_hashes"] != frozen["provenance"]["source_hashes"]:
        raise ValueError("benchmark source differs from frozen calibration")
    if current["versions"] != frozen["provenance"]["versions"]:
        raise ValueError("benchmark dependencies differ from frozen calibration")
    if current["versions"] != LOCKED_VERSIONS:
        raise ValueError("benchmark requires the locked GP NumPy/SciPy runtime")
    if frozen.get("standard_deviation_definition") != STD_DEFINITION:
        raise ValueError("ambiguous predictive standard-deviation convention")
    if frozen.get("timing_scopes") != TIMING_SCOPES:
        raise ValueError("ambiguous scientific timing conventions")
    if canonical:
        if not current["source_clean_and_tracked"]:
            raise ValueError("canonical cases require a clean committed implementation")
        if tuple(frozen["sizes"]) != SIZES or tuple(frozen["ranks"]) != RANKS:
            raise ValueError(
                "canonical calibration must cover every accepted size and rank"
            )
        if any(os.environ.get(key) != "1" for key in THREADS):
            raise ValueError(
                "canonical campaign requires fixed single-thread environment"
            )


def run_case(data, frozen, size, repeat, output: Path, *, canonical: bool):
    validate_frozen(frozen, canonical=canonical)
    start_provenance = provenance()
    if (
        type(size) is not int
        or size not in frozen["sizes"]
        or type(repeat) is not int
        or repeat not in (0, 1, 2)
    ):
        raise ValueError("unexpected size or repeat")
    begin = time.perf_counter()
    lengthscale, noise = float(frozen["lengthscale"]), float(frozen["noise"])
    split = split_data(data, size, "test")
    rank = frozen["selected_ranks"][str(size)]
    start = time.perf_counter()
    factor = KernelCholeskyFactor.from_columns(
        np.ones(size),
        lambda j: rbf(split["x"], split["x"][[j]], lengthscale)[:, 0],
        max_rank=rank,
        atol=1e-12,
    )
    construction_seconds = time.perf_counter() - start
    start = time.perf_counter()
    mean, std, shifted, solved = predictions(factor, split, lengthscale, noise)
    prediction_seconds = time.perf_counter() - start
    application_seconds = time.perf_counter() - begin
    expected_mean = expected_std = np.array([], dtype=np.float64)
    metrics = None
    dense_seconds = None
    if size <= 2048:
        start = time.perf_counter()
        expected_mean, expected_std = dense_reference(split, lengthscale, noise)
        dense_seconds = time.perf_counter() - start
        metrics = accuracy(mean, std, expected_mean, expected_std)
    start = time.perf_counter()
    operation = original_operator(split["x"], lengthscale, noise)
    iterations = 0

    def iteration(_):
        nonlocal iterations
        iterations += 1

    solution, info = cg(
        operation,
        split["y"],
        M=shifted.as_preconditioner(),
        rtol=1e-9,
        atol=0,
        maxiter=1000,
        callback=iteration,
    )
    residual = operation @ solution - split["y"]
    cg_error = float(np.linalg.norm(residual) / np.linalg.norm(split["y"]))
    cg_seconds = time.perf_counter() - start
    arrays_path = output.with_suffix(".npz")
    arrays_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        arrays_path,
        root=factor.factor,
        mean=mean,
        std=std,
        expected_mean=expected_mean,
        expected_std=expected_std,
        query_y=split["query_y"],
        cg_solution=solution,
        cg_residual=residual,
        cg_rhs=split["y"],
        pivots=np.array(factor.indices, dtype=np.int64),
        approximate_solution=solved.solution,
    )
    retained_components = {
        "shifted_factor": shifted.storage_nbytes,
        "source_values": data["values"].nbytes,
        "source_ordinals": data["ordinals"].nbytes,
        "normalized_split": sum(
            split[key].nbytes for key in ("x", "y", "query_x", "query_y")
        ),
        "prediction_arrays": mean.nbytes
        + std.nbytes
        + expected_mean.nbytes
        + expected_std.nbytes,
        "solve_arrays": solution.nbytes + residual.nbytes + solved.solution.nbytes,
    }
    retained_bytes = sum(retained_components.values())
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss_bytes = int(rss if sys.platform == "darwin" else rss * 1024)
    accuracy_passed = metrics is None or (
        max(metrics["mean_rms"], metrics["std_rms"]) <= 1e-3
        and max(metrics["mean_max"], metrics["std_max"]) <= 1e-2
    )
    checks = {
        "dense_accuracy": accuracy_passed,
        "approximate_solve": solved.relative_residual <= 1e-10,
        "original_kernel_cg": info == 0 and cg_error <= 1e-8,
        "retained_state": retained_bytes <= 128 * 1024**2,
        "rss": rss_bytes <= 4 * 1024**3,
        "source_clean": (
            provenance()["source_clean_and_tracked"]
            and start_provenance == provenance()
        )
        if canonical
        else False,
    }
    passed = all(checks.values())
    result = {
        "schema": "kernelgecp.covariance-case.v1",
        "status": "passed" if passed else ("failed" if canonical else "pilot"),
        "canonical": canonical,
        "size": size,
        "repeat": repeat,
        "rank_budget": rank,
        "actual_rank": factor.rank,
        "checks": checks,
        "provenance": provenance(),
        "dataset": dataset_identity(data),
        "frozen_calibration_sha256": sha256(
            json.dumps(frozen, sort_keys=True).encode()
        ),
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
        "accuracy": metrics,
        "standard_deviation_definition": STD_DEFINITION,
        "solve_relative_residual": solved.relative_residual,
        "cg_relative_residual": cg_error,
        "cg_info": int(info),
        "cg_iterations": iterations,
        "retained_numerical_bytes": retained_bytes,
        "retained_components_bytes": retained_components,
        "peak_rss_bytes": rss_bytes,
        "factor_metadata": factor.metadata.as_dict(),
        "timing_scopes": TIMING_SCOPES,
        "timing_seconds": {
            "application_workflow": application_seconds,
            "construction": construction_seconds,
            "prediction": prediction_seconds,
            "dense_reference": dense_seconds,
            "original_kernel_cg": cg_seconds,
            "full_workflow": time.perf_counter() - begin,
        },
        "predictive_rmse_original_units": float(
            np.sqrt(np.mean((mean - split["query_y"]) ** 2))
            * split["normalization"]["target_std"]
        ),
        "artifact": {
            "path": arrays_path.name,
            "sha256": sha256(arrays_path.read_bytes()),
        },
        "scope": (
            "finite_real_data_representation_and_linear_algebra_"
            "not_competitive_forecast_or_optimal_placement"
        ),
    }
    write_json(output, result)
    return result


def run_provider(data, frozen, output: Path) -> dict[str, Any]:
    """Actual locked CPU provider on the same real-data frozen n512 factor."""
    import gpytorch
    import torch

    validate_frozen(frozen, canonical=True)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "torch", "gpytorch", "linear_operator")
    }
    expected_versions = {
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "torch": "2.14.1",
        "gpytorch": "1.15.2",
        "linear_operator": "0.6.1",
    }
    if any(
        versions[name].split("+")[0] != expected
        for name, expected in expected_versions.items()
    ):
        raise ValueError("actual GP provider differs from the supported lock")
    split = split_data(data, 512, "test")
    lengthscale, noise = frozen["lengthscale"], frozen["noise"]
    factor = KernelCholeskyFactor.from_columns(
        np.ones(512),
        lambda j: rbf(split["x"], split["x"][[j]], lengthscale)[:, 0],
        max_rank=frozen["selected_ranks"]["512"],
        atol=1e-12,
    )
    mean, std, shifted, _ = predictions(factor, split, lengthscale, noise)
    converted = factor.to_gpytorch(noise)
    rhs = torch.tensor(split["y"], dtype=torch.float64)
    solved = converted.solve(rhs).numpy()
    extension = factor.extend(
        rbf(split["query_x"], split["x"][list(factor.indices)], lengthscale)
    )
    cross = extension.features @ factor.factor.T
    provider_mean = cross @ solved
    provider_variance = (
        1
        + noise
        - np.sum(cross * converted.solve(torch.tensor(cross.T)).numpy().T, axis=1)
    )
    if np.any(provider_variance < 0):
        raise RuntimeError("actual provider produced negative prediction variance")
    provider_std = np.sqrt(provider_variance)
    metrics = accuracy(provider_mean, provider_std, mean, std)
    metrics["solve_max"] = float(np.max(np.abs(solved - shifted.solve(split["y"]))))
    matrix = factor.factor @ factor.factor.T + noise * np.eye(512)
    chol = cho_factor(matrix, lower=True)
    logdet = float(2 * np.log(np.diag(chol[0])).sum())
    logprob = -0.5 * (
        split["y"] @ cho_solve(chol, split["y"]) + logdet + 512 * np.log(2 * np.pi)
    )
    distribution = gpytorch.distributions.MultivariateNormal(
        torch.zeros(512, dtype=torch.float64), converted
    )
    metrics["log_probability_error"] = abs(float(distribution.log_prob(rhs)) - logprob)
    checks = {
        "cpu_float64": converted.device.type == "cpu"
        and converted.dtype == torch.float64,
        "frozen_factor": not converted.requires_grad,
        "numerical": all(
            np.isfinite(value) and value <= 1e-9 for value in metrics.values()
        ),
    }
    arrays_path = output.with_suffix(".npz")
    arrays_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        arrays_path,
        root=factor.factor,
        mean=provider_mean,
        std=provider_std,
        expected_mean=mean,
        expected_std=std,
        solution=solved,
        rhs=split["y"],
    )
    result = {
        "schema": "kernelgecp.covariance-provider.v1",
        "status": "passed" if all(checks.values()) else "failed",
        "size": 512,
        "provider": "gpytorch",
        "versions": versions,
        "provenance": provenance(),
        "dataset": dataset_identity(data),
        "metrics": metrics,
        "standard_deviation_definition": STD_DEFINITION,
        "checks": checks,
        "frozen_calibration_sha256": sha256(
            json.dumps(frozen, sort_keys=True).encode()
        ),
        "artifact": {
            "path": arrays_path.name,
            "sha256": sha256(arrays_path.read_bytes()),
        },
    }
    write_json(output, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("calibrate", "case", "campaign", "provider")
    )
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frozen", type=Path)
    parser.add_argument("--size", type=int)
    parser.add_argument("--repeat", type=int, default=0)
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "calibrate":
            sizes = (512, 1024, 2048) if args.pilot else SIZES
            result = calibration(load_data(args.data), sizes)
            write_json(args.output, result)
            if result["status"] != "calibrated":
                raise RuntimeError(
                    "calibration budget exhausted; failure evidence retained"
                )
        else:
            if args.frozen is None:
                parser.error("--frozen is required")
            frozen = json.loads(args.frozen.read_text())
            if args.command == "provider":
                result = run_provider(load_data(args.data), frozen, args.output)
                if result["status"] != "passed":
                    raise RuntimeError("actual provider gate failed; evidence retained")
            elif args.command == "case":
                result = run_case(
                    load_data(args.data),
                    frozen,
                    args.size,
                    args.repeat,
                    args.output,
                    canonical=not args.pilot,
                )
                if result["status"] == "failed":
                    raise RuntimeError("a covariance gate failed; evidence retained")
            else:
                validate_frozen(frozen, canonical=True)
                args.output.mkdir(parents=True, exist_ok=True)
                case_paths = []
                for size in SIZES:
                    for repeat in range(3):
                        path = args.output / f"n{size}-repeat{repeat}.json"
                        command = [
                            sys.executable,
                            "-m",
                            "experiments.covariance_benchmark",
                            "case",
                            "--data",
                            str(args.data.resolve()),
                            "--frozen",
                            str(args.frozen.resolve()),
                            "--output",
                            str(path.resolve()),
                            "--size",
                            str(size),
                            "--repeat",
                            str(repeat),
                        ]
                        print(f"running n={size} repeat={repeat}", flush=True)
                        process_start = time.perf_counter()
                        with path.with_suffix(".log").open("w") as log:
                            completed = subprocess.run(
                                command,
                                cwd=ROOT,
                                check=False,
                                stdout=log,
                                stderr=subprocess.STDOUT,
                            )
                        case_paths.append(
                            {
                                "path": path.name,
                                "exit_code": completed.returncode,
                                "process_wall_seconds": time.perf_counter()
                                - process_start,
                                "sha256": sha256(path.read_bytes())
                                if path.exists()
                                else None,
                            }
                        )
                write_json(
                    args.output / "campaign.json",
                    {
                        "schema": "kernelgecp.covariance-campaign.v1",
                        "provenance": provenance(),
                        "frozen_path": str(args.frozen.resolve()),
                        "cases": case_paths,
                    },
                )
                if any(case["exit_code"] != 0 for case in case_paths):
                    raise RuntimeError(
                        "campaign has failed processes; evidence retained"
                    )
    except Exception as error:
        if (
            args.command in ("case", "provider", "calibrate")
            and not args.output.exists()
        ):
            write_json(
                args.output,
                {
                    "schema": f"kernelgecp.covariance-{args.command}.v1",
                    "status": "failed",
                    "error": f"{type(error).__name__}: {error}",
                    "provenance": provenance(),
                },
            )
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
