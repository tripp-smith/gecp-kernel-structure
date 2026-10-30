"""Actual fixed-bath CT-HYB spectroscopy with independent chain validation."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import threading
import time
from contextlib import suppress
from pathlib import Path

import numpy as np

from experiments.cthyb_dmft import load_measurement
from experiments.spectroscopy import METHODS, SOURCE_FILES, THREADS, basis, write_json
from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    continue_tau_maxent,
    covariance_modes,
    prepare_spectral_projection,
)

ROOT = Path(__file__).resolve().parents[1]
CALIBRATION_REVISION = "9d28053060329cd20647dc03e038a383c351adf1"
FINAL_REVISION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
CALIBRATION_SEED_BASE = 29100000
FINAL_SEED_BASE = 30100000
OBS = np.arange(0, 3073, 48)
QUERY = 24 + 192 * np.arange(16)
RTOL = 1e-4
IMAG_TOL = 1e-10


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_cov(samples):
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 2 or samples.shape[0] < 2 or not np.isfinite(samples).all():
        raise ValueError(
            "finite independent sample matrix with at least two rows required"
        )
    mean = samples.mean(axis=0)
    residual = samples - mean
    return mean, residual.T @ residual / (len(samples) * (len(samples) - 1))


def load_cohort(folder, receipt_path, binding_path, binding_sha256, stage):
    """Validate audit, raw receipts, hashes, seeds, shape and one common bath."""
    if stage not in ("calibration", "final"):
        raise ValueError("unsupported stage")
    revision = CALIBRATION_REVISION if stage == "calibration" else FINAL_REVISION
    seed_base = CALIBRATION_SEED_BASE if stage == "calibration" else FINAL_SEED_BASE
    if len(binding_sha256) != 64 or sha(binding_path) != binding_sha256:
        raise ValueError("audit binding digest mismatch")
    binding = json.loads(binding_path.read_text())
    if stage == "final" and binding.get("generation_revision") == CALIBRATION_REVISION:
        raise ValueError("consumed final cohort cannot be reused")
    if (
        binding.get("schema") != "kernelgecp.cthyb-audit-binding.v1"
        or binding.get("status") != "passed"
        or binding.get("generation_revision") != revision
        or binding.get("audit_receipt_sha256") != sha(receipt_path)
        or binding.get("result_sha256") != sha(folder / "result.json")
        or binding.get("manifest_sha256") != sha(folder / "manifest.json")
        or binding.get("freeze_sha256") != sha(folder / "freeze.json")
    ):
        raise ValueError("audit binding does not match source")
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("schema") != "kernelgecp.cthyb-audit.v1"
        or receipt.get("status") != "passed"
        or receipt.get("generation_revision") != revision
        or not receipt.get("gates")
        or not all(v is True for v in receipt["gates"].values())
    ):
        raise ValueError("independent EE audit did not pass for expected revision")
    audit_blob = subprocess.check_output(
        ["git", "show", receipt["audit_revision"] + ":experiments/audit_cthyb_dmft.py"],
        cwd=ROOT,
    )
    if hashlib.sha256(audit_blob).hexdigest() != receipt.get("audit_sha256"):
        raise ValueError("audit source hash mismatch")
    identity = json.loads((folder / "manifest.json").read_text())
    result = json.loads((folder / "result.json").read_text())
    config = identity["config"]
    if (
        identity.get("source_revision") != revision
        or result.get("identity") != identity
        or result.get("status") != "passed"
        or result.get("freeze_sha256") != sha(folder / "freeze.json")
        or config.get("seed_base") != seed_base
        or config.get("seed_mapping") != "odd_stride_2_v1"
        or config.get("final_replicas") != 64
        or config.get("n_tau") != 3073
        or config.get("beta") != 10
        or config.get("interaction") != 2
    ):
        raise ValueError("campaign identity, result, or configuration mismatch")
    role, offset = (
        ("kernelgecp-final", 40000)
        if stage == "calibration"
        else ("dense-final", 30000)
    )
    paths = sorted((folder / role).glob("replica-*.json"))
    if len(paths) != 64 or paths != [
        folder / role / f"replica-{i:02d}.json" for i in range(64)
    ]:
        raise ValueError("cohort must have exactly 64 ordered replicas")
    samples, seeds, hashes, bath = [], set(), [], None
    for i, path in enumerate(paths):
        digest = sha(path)
        if (
            result["measurement_archives"].get(str(path.relative_to(folder))) != digest
            or binding.get("measurement_archives", {}).get(
                str(path.relative_to(folder))
            )
            != digest
        ):
            raise ValueError("measurement receipt hash mismatch")
        measurement = load_measurement(path.with_suffix(""), identity)
        if (
            measurement.seed != config["seed_base"] + offset + 1 + 2 * i
            or measurement.seed in seeds
        ):
            raise ValueError("duplicate or incorrect odd-stride seed")
        seeds.add(measurement.seed)
        if (
            measurement.config.beta != 10
            or measurement.config.interaction != 2
            or measurement.config.n_tau != 3073
            or measurement.cycles != 1600000
        ):
            raise ValueError("measurement physics or cycle count mismatch")
        if measurement.solve_status != 0 or measurement.average_sign <= 0.99:
            raise ValueError("native chain solve or sign failed")
        if bath is None:
            bath = measurement.bath_iw.copy()
        elif not np.array_equal(bath, measurement.bath_iw):
            raise ValueError("cohort contains different baths")
        values = measurement.g_tau
        if values.shape != (3073, 2) or not np.isfinite(values).all():
            raise ValueError("invalid G(tau) array")
        if np.max(np.abs(values.imag)) > IMAG_TOL:
            raise ValueError("imaginary G(tau) exceeds declared tolerance")
        samples.append(values.real.mean(axis=1))
        hashes.append(digest)
    return np.array(samples), {
        "stage": stage,
        "role": role,
        "generation_revision": revision,
        "seed_base": seed_base,
        "manifest_sha256": sha(folder / "manifest.json"),
        "freeze_sha256": sha(folder / "freeze.json"),
        "seeds": sorted(seeds),
        "receipt_sha256": sha(receipt_path),
        "binding_sha256": binding_sha256,
        "result_sha256": sha(folder / "result.json"),
        "measurement_receipt_sha256": hashes,
        "spin_operation": "paramagnetic_arithmetic_mean_up_down",
        "imaginary_absolute_tolerance": IMAG_TOL,
    }


def dataset_and_validation(samples, metadata):
    if samples.shape != (64, 3073) or not np.isfinite(samples).all():
        raise ValueError("cohort must contain 64 finite 3073-point chains")
    train, held = samples[:32], samples[32:]
    obs_mean, obs_cov = mean_cov(train[:, OBS])
    query_mean, held_cov = mean_cov(held[:, QUERY])
    tau = np.linspace(0, 10, 3073)
    data = SpectralDataset(
        10, tau[OBS], obs_mean, obs_cov, covariance_rtol=RTOL, provenance=metadata
    )
    return data, tau[QUERY], query_mean, held_cov


def validation_score(covariance, residual):
    _, modes = covariance_modes(covariance, rtol=RTOL)
    rms = modes.whitened_rms(residual)
    return rms, modes, bool(rms <= 2)


def weights(omega, density):
    cdf = np.r_[0.0, np.cumsum(np.diff(omega) * (density[:-1] + density[1:]) / 2)]
    return np.diff(np.interp([-4, -1, 1, 4], omega, cdf))


def mode_record(modes):
    return {
        "eigenvalues": modes.eigenvalues.tolist(),
        "threshold": modes.threshold,
        "relative_cutoff": modes.relative_cutoff,
        "absolute_cutoff": modes.absolute_cutoff,
        "effective_rank": modes.effective_rank,
        "discarded_eigenvalues": modes.eigenvalues[
            modes.eigenvalues <= modes.threshold
        ].tolist(),
    }


# The protocol is fixed before calibration and final validation. Repeats are
# fresh solver invocations on the same cohort, not independent data realizations.
PROTOCOL = {
    "methods": list(METHODS),
    "default_models": ["flat", "broad_gaussian"],
    # One common native-data candidate for every method and model; the ladder
    # is a declared history, not an automatic or method-dependent fallback.
    "basis_control_ladder": [
        {
            "name": "baseline",
            "eps": 1e-8,
            "rcond": 1e-10,
            "outcome": "failed_partial_calibration_pydlr_flat",
            "source_revision": "d5070509b675706e3e18164de96660530ac0489a",
            "failed_whitened_rms": 2.0585741570745255,
            "failed_projection_rank": 20,
            "failed_projected_covariance_rank": 19,
            "progress_sha256": (
                "c84df6e1dc8d2b2afce0f5ea7954ebdc8ee444cfac569d2fb803e582e7de1747"
            ),
            "failed_result_sha256": (
                "03c506bc6c864b8acf7c7f88ba9b08b16c3a17c13c67de97cace5ff28bd5739b"
            ),
            "failed_array_sha256": (
                "ebe689065afd27d657047af795710adb373e96bdb3319d05d55273feb339c523"
            ),
        },
        {"name": "candidate1", "eps": 1e-10, "rcond": 1e-12},
    ],
    "selected_basis_control": "candidate1",
    "eps": 1e-10,
    "rcond": 1e-12,
    "n_omega": 401,
    "n_alpha": 33,
    "covariance_rtol": RTOL,
    "bootstrap_count": 8,
    "bootstrap_repeat": 0,
    "bootstrap_seed": 193,
    "calibration_repeats": 1,
    "final_repeats": 3,
    "validation": "independent_heldout_chain_mean_conditional_on_training_fit",
    "alpha_selection": "normalized_training_only_discrepancy",
    "consumed_final_failure_array_sha256": (
        "407c85405e9f56135442db29bb2a4c75df8ef77b37e3f831a58c754e562817ba"
    ),
    "calibration_generation_revision": CALIBRATION_REVISION,
    "final_generation_revision": FINAL_REVISION,
    "final_seed_base": FINAL_SEED_BASE,
}


def validate_basis_controls():
    """Fail closed if the declared candidate and executed controls diverge."""
    selected = [
        item
        for item in PROTOCOL["basis_control_ladder"]
        if item["name"] == PROTOCOL["selected_basis_control"]
    ]
    if len(selected) != 1 or any(
        selected[0][key] != PROTOCOL[key] for key in ("eps", "rcond")
    ):
        raise ValueError("selected native basis controls do not match protocol")


def snapshot():
    """Require all scientific source blobs to match the committed revision."""
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    tracked = subprocess.check_output(
        ["git", "ls-files", "kernelgecp", "experiments"], cwd=ROOT, text=True
    ).splitlines()
    files = sorted(
        set(SOURCE_FILES)
        | {
            "experiments/cthyb_spectroscopy.py",
            "experiments/cthyb_dmft.py",
            "experiments/audit_cthyb_spectroscopy.py",
        }
        | {p for p in tracked if p.startswith("kernelgecp/") and p.endswith(".py")}
    )
    hashes = {}
    for name in files:
        blob = subprocess.check_output(["git", "show", revision + ":" + name], cwd=ROOT)
        digest = sha(ROOT / name)
        if digest != hashlib.sha256(blob).hexdigest():
            raise ValueError("scientific source differs from committed blob: " + name)
        hashes[name] = digest
    threads = {name: os.environ.get(name) for name in (*THREADS, "BLIS_NUM_THREADS")}
    if any(value != "1" for value in threads.values()):
        raise ValueError("campaign requires declared single-thread environment")
    modules = [
        "triqs_maxent",
        "triqs_maxent.kernels",
        "triqs_maxent.version",
        "triqs",
        "triqs_cthyb",
    ]
    imported = {}
    for name in modules:
        module = importlib.import_module(name)
        path = Path(module.__file__).resolve()
        imported[name] = {"path": str(path), "sha256": sha(path)}
    # Hash installed code and native extensions, excluding generated caches.
    runtime_trees = {}
    for package in (
        "triqs_maxent",
        "triqs",
        "triqs_cthyb",
        "numpy",
        "scipy",
        "mpmath",
        "pydlr",
        "sparse_ir",
        "pylibsparseir",
    ):
        module = importlib.import_module(package)
        root = Path(module.__file__).resolve().parent
        digest = hashlib.sha256()
        count = 0
        for file in sorted(root.rglob("*")):
            if (
                file.is_file()
                and file.suffix in (".py", ".so", ".dylib")
                and "__pycache__" not in file.parts
            ):
                digest.update(str(file.relative_to(root)).encode() + b"\0")
                digest.update(bytes.fromhex(sha(file)))
                count += 1
        if not count:
            raise ValueError(
                "runtime package contains no source or native modules: " + package
            )
        runtime_trees[package] = {
            "root": str(root),
            "files": count,
            "sha256": digest.hexdigest(),
        }
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir")
    }
    version = importlib.import_module("triqs_maxent.version").version
    if str(version) != "4.0.1":
        raise ValueError("requires TRIQS MaxEnt 4.0.1")
    versions["triqs_maxent"] = str(version)
    versions["triqs"] = str(importlib.import_module("triqs.version").version)
    versions["triqs_cthyb"] = str(
        importlib.import_module("triqs_cthyb.version").version
    )
    if versions["triqs"] != "4.0.2" or versions["triqs_cthyb"] != "4.0.1":
        raise ValueError("unsupported TRIQS/CT-HYB provider versions")
    return {
        "revision": revision,
        "sources": hashes,
        "threads": threads,
        "python": sys.version,
        "executable": str(Path(sys.executable).resolve()),
        "executable_sha256": sha(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "packages": versions,
        "provider_modules": imported,
        "runtime_trees": runtime_trees,
    }


def verify_calibration(path, digest, snapshot_value):
    validate_basis_controls()
    if not path or not digest or len(digest) != 64 or sha(path) != digest:
        raise ValueError("final requires pinned passed calibration receipt")
    receipt = json.loads(path.read_text())
    if (
        receipt.get("schema") != "kernelgecp.cthyb-spectroscopy-calibration.v1"
        or receipt.get("status") != "passed"
        or receipt.get("protocol") != PROTOCOL
        or receipt.get("snapshot") != snapshot_value
    ):
        raise ValueError(
            "calibration freeze does not match protocol, source, runtime or EE binding"
        )
    cohort = receipt.get("cohort", {})
    if (
        cohort.get("stage") != "calibration"
        or cohort.get("role") != "kernelgecp-final"
        or cohort.get("generation_revision") != CALIBRATION_REVISION
        or cohort.get("seed_base") != CALIBRATION_SEED_BASE
        or cohort.get("seeds")
        != [CALIBRATION_SEED_BASE + 40001 + 2 * i for i in range(64)]
        or cohort.get("binding_sha256") != receipt.get("binding_sha256")
        or not isinstance(receipt.get("binding_sha256"), str)
        or len(receipt["binding_sha256"]) != 64
    ):
        raise ValueError("calibration cohort binding or identity mismatch")
    expected = {
        f"{method}-{model}-r0.json"
        for method in METHODS
        for model in PROTOCOL["default_models"]
    }
    entries = receipt.get("cases", {})
    if set(entries) != expected:
        raise ValueError("calibration case inventory mismatch")
    for name, hashes in entries.items():
        case = path.parent / name
        archive = case.with_suffix(".npz")
        if (
            sha(case) != hashes["result_sha256"]
            or sha(archive) != hashes["array_sha256"]
        ):
            raise ValueError("calibration case hash mismatch")
        value = json.loads(case.read_text())
        if (
            value.get("status") != "passed"
            or value.get("array_sha256") != sha(archive)
            or value.get("eps") != PROTOCOL["eps"]
            or value.get("rcond") != PROTOCOL["rcond"]
        ):
            raise ValueError("calibration case failed")
    return receipt


def validate_cohort_pair(calibration, final):
    if (
        calibration["generation_revision"] == final["generation_revision"]
        or calibration["binding_sha256"] == final["binding_sha256"]
        or set(calibration["seeds"]) & set(final["seeds"])
    ):
        raise ValueError("calibration and final cohort bindings or seeds overlap")


def verify_calibration_audit(path, digest, calibration_digest, binding_digest):
    if not path or not digest or len(digest) != 64 or sha(path) != digest:
        raise ValueError("final requires a pinned independent calibration audit")
    audit = json.loads(path.read_text())
    if (
        audit.get("schema") != "kernelgecp.cthyb-spectroscopy-independent-audit.v1"
        or audit.get("status") != "passed"
        or audit.get("kind") != "calibration"
        or audit.get("calibration_sha256") != calibration_digest
        or audit.get("binding_sha256") != binding_digest
        or not audit.get("gates")
        or not all(value is True for value in audit["gates"].values())
    ):
        raise ValueError(
            "independent calibration audit did not pass or binding mismatch"
        )
    blob = subprocess.check_output(
        [
            "git",
            "show",
            audit["auditor_revision"] + ":experiments/audit_cthyb_spectroscopy.py",
        ],
        cwd=ROOT,
    )
    if hashlib.sha256(blob).hexdigest() != audit.get("auditor_sha256"):
        raise ValueError("independent calibration auditor hash mismatch")
    return audit


def run(
    folder,
    receipt,
    binding,
    digest,
    stage,
    method,
    model,
    output,
    *,
    bootstrap=0,
    seed=193,
    eps=1e-10,
    rcond=1e-12,
    n_omega=401,
    n_alpha=33,
):
    if method not in METHODS or model not in ("flat", "broad_gaussian"):
        raise ValueError("unsupported method or default model")
    if bootstrap not in (0, *range(2, 65)):
        raise ValueError("bootstrap must be zero or between 2 and 64")
    workflow_start = time.perf_counter()
    samples, metadata = load_cohort(folder, receipt, binding, digest, stage)
    data, query, reference, validation_cov = dataset_and_validation(samples, metadata)
    construction_start = time.perf_counter()
    projection = None
    design = np.eye(len(data.tau))
    if method != "raw":
        b = basis(method, data.beta, eps)
        design = b.design_tau(data.tau)
        projection = prepare_spectral_projection(
            b, data.tau, data.tau, rcond=rcond, metadata={"method": method, "eps": eps}
        )
    construction_seconds = time.perf_counter() - construction_start
    projection_start = time.perf_counter()
    fitted = data if projection is None else projection.apply(data)
    projection_seconds = time.perf_counter() - projection_start
    config = MaxEntConfig(n_omega=n_omega, n_alpha=n_alpha, default_model=model)
    continuation_start = time.perf_counter()
    spectrum = continue_tau_maxent(fitted, config)
    continuation_seconds = time.perf_counter() - continuation_start
    prediction_start = time.perf_counter()
    predicted = spectrum.evaluate_tau(data.beta, query)
    prediction_seconds = time.perf_counter() - prediction_start
    rms, modes, heldout_pass = validation_score(validation_cov, predicted - reference)
    p = np.eye(len(data.tau)) if projection is None else projection.matrix
    from triqs_maxent import LinearOmegaMesh
    from triqs_maxent.kernels import TauKernel

    native_kernel = TauKernel(
        data.tau.copy(),
        LinearOmegaMesh(config.omega_min, config.omega_max, config.n_omega),
        beta=data.beta,
    ).K
    boot_weights, boot_spectra, boot_indices, boot_records = [], [], [], []
    rng = np.random.default_rng(seed)
    bootstrap_start = time.perf_counter()
    for _ in range(bootstrap):
        indices = rng.integers(0, 32, size=32)
        boot_indices.append(indices)
        sample_mean, sample_cov = mean_cov(samples[indices][:, OBS])
        boot_data = SpectralDataset(
            10,
            data.tau,
            sample_mean,
            sample_cov,
            covariance_rtol=RTOL,
            provenance={"bootstrap": "independent_chain_resampling"},
        )
        boot_data = boot_data if projection is None else projection.apply(boot_data)
        trial = continue_tau_maxent(boot_data, config)
        boot_weights.append(weights(trial.omega, trial.density))
        boot_spectra.append(trial.density)
        boot_records.append(
            {
                "selected_index": trial.selected_index,
                "converged": trial.converged,
                "iterations": trial.iterations,
                "effective_rank": trial.effective_rank,
            }
        )
    bootstrap_seconds = time.perf_counter() - bootstrap_start
    gates = {
        "positive": bool(np.min(spectrum.density) >= -1e-6),
        "normalized": bool(abs(spectrum.mass - 1) <= 1e-6),
        "heldout_whitened_rms": heldout_pass,
        "provider_selected_converged": bool(
            spectrum.converged[spectrum.selected_index]
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    archive = output.with_suffix(".npz")
    if output.exists() or archive.exists():
        raise FileExistsError("refusing to overwrite output")
    np.savez_compressed(
        archive,
        omega=spectrum.omega,
        density=spectrum.density,
        observation_tau=data.tau,
        observation_values=data.values,
        covariance=data.covariance,
        projected_values=fitted.values,
        projected_covariance=fitted.covariance,
        projection=p,
        observation_design=design,
        native_kernel=native_kernel,
        projected_kernel=p @ native_kernel,
        query=query,
        query_mean=reference,
        validation_covariance=validation_cov,
        prediction=predicted,
        alphas=spectrum.alphas,
        chi2=spectrum.chi2,
        trial_spectra=spectrum.trial_spectra,
        bootstrap_weights=np.array(boot_weights),
        bootstrap_spectra=np.array(boot_spectra),
        bootstrap_indices=np.array(boot_indices),
    )
    result = {
        "schema": "kernelgecp.cthyb-spectroscopy.v1",
        "status": "passed" if all(gates.values()) else "failed",
        "gates": gates,
        "source": metadata,
        "method": method,
        "config": dataclasses.asdict(config),
        "eps": eps,
        "rcond": rcond,
        "observation_indices": OBS.tolist(),
        "query_indices": QUERY.tolist(),
        "covariance_estimator": "sum_residual_outer_over_N_times_N_minus_1",
        "covariance_rtol": RTOL,
        "validation_covariance": (
            "heldout_chain_mean_covariance_conditional_on_training_fit; "
            "excludes_nonlinear_prediction_uncertainty"
        ),
        "training_covariance_modes": mode_record(data.modes),
        "projected_covariance_modes": mode_record(fitted.modes),
        "validation_covariance_modes": mode_record(modes),
        "projection_svd": None
        if projection is None
        else dataclasses.asdict(projection.diagnostics),
        "provider_diagnostics": spectrum.diagnostics,
        "normalized_training_whitened_rms": spectrum.normalized_training_whitened_rms,
        "units": {"tau": "inverse_energy", "omega": "energy", "G_tau": "dimensionless"},
        "training_rank": data.modes.effective_rank,
        "projected_rank": fitted.modes.effective_rank,
        "validation_rank": modes.effective_rank,
        "whitened_rms": rms,
        "mass": spectrum.mass,
        "raw_mass": spectrum.raw_mass,
        "selected_index": spectrum.selected_index,
        "converged": spectrum.converged,
        "iterations": spectrum.iterations,
        "window_weights": weights(spectrum.omega, spectrum.density).tolist(),
        "bootstrap_count": bootstrap,
        "bootstrap_seed": seed,
        "bootstrap_records": boot_records,
        "bootstrap_window_standard_deviation": np.std(
            boot_weights, axis=0, ddof=1
        ).tolist()
        if bootstrap
        else None,
        "claim_scope": (
            "finite_grid_fixed_bath_empirical_diagnostics_no_native_truth_"
            "no_unique_spectrum_or_coverage_claim"
        ),
        "array_sha256": sha(archive),
        "timing_seconds": {
            "construction": construction_seconds,
            "projection": projection_seconds,
            "continuation": continuation_seconds,
            "bootstrap": bootstrap_seconds,
            "physical_prediction": prediction_seconds,
            "workflow_before_receipt": time.perf_counter() - workflow_start,
        },
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        * (1 if sys.platform == "darwin" else 1024),
        "archive_bytes": archive.stat().st_size,
    }
    write_json(output, result)
    return result


GLOBAL_SECONDS = 330 * 60
CASE_SECONDS = 60 * 60


def launch_worker(command, log_path, timeout, on_started=None):
    """Own a fresh process group, reap it on every exit, and retain its log."""
    started = time.monotonic()
    with log_path.open("xb") as log:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            if on_started is not None:
                on_started(process.pid)
            try:
                code = process.wait(timeout=timeout)
                status = "passed" if code == 0 else "failed"
            except subprocess.TimeoutExpired:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                code, status = process.returncode, "timeout"
        except BaseException:
            if process.poll() is None:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    return {
        "pid": process.pid,
        "command": command,
        "exit_code": code,
        "elapsed_seconds": time.monotonic() - started,
        "status": status,
        "log_file": log_path.name,
        "log_sha256": sha(log_path),
    }


def worker(job_path, job_digest):
    validate_basis_controls()
    if len(job_digest) != 64 or sha(job_path) != job_digest:
        raise ValueError("worker job hash mismatch")
    job = json.loads(job_path.read_text())
    if (
        job.get("schema") != "kernelgecp.cthyb-spectroscopy-worker-job.v1"
        or job.get("protocol") != PROTOCOL
    ):
        raise ValueError("invalid worker job protocol")
    method, model, repeat, stage = (
        job["method"],
        job["model"],
        job["repeat"],
        job["stage"],
    )
    if (
        stage not in ("calibration", "final")
        or method not in METHODS
        or model not in PROTOCOL["default_models"]
        or type(repeat) is not int
        or not 0 <= repeat < PROTOCOL[stage + "_repeats"]
    ):
        raise ValueError("worker case outside frozen inventory")
    if job.get("parent_pid") != os.getppid():
        raise ValueError("worker parent process mismatch")
    output = Path(job["output"])
    if output.name != f"{method}-{model}-r{repeat}.json":
        raise ValueError("worker output does not match case")
    frozen = snapshot()
    if frozen != job["snapshot"]:
        raise ValueError("worker source/runtime snapshot mismatch")
    if os.getpgrp() != os.getpid():
        raise ValueError("worker requires its own process group")
    parent_pid = job["parent_pid"]

    def watch_parent():
        while True:
            time.sleep(1)
            if os.getppid() != parent_pid:
                os.killpg(os.getpgrp(), signal.SIGKILL)

    threading.Thread(target=watch_parent, daemon=True).start()
    if stage == "final":
        calibration_value = verify_calibration(
            Path(job["calibration"]), job["calibration_sha256"], frozen
        )
        verify_calibration_audit(
            Path(job["calibration_audit"]),
            job["calibration_audit_sha256"],
            job["calibration_sha256"],
            calibration_value["binding_sha256"],
        )
    result = run(
        Path(job["folder"]),
        Path(job["audit_receipt"]),
        Path(job["audit_binding"]),
        job["binding_sha256"],
        stage,
        method,
        model,
        output,
        bootstrap=PROTOCOL["bootstrap_count"] if repeat == 0 else 0,
        seed=PROTOCOL["bootstrap_seed"],
        eps=PROTOCOL["eps"],
        rcond=PROTOCOL["rcond"],
        n_omega=PROTOCOL["n_omega"],
        n_alpha=PROTOCOL["n_alpha"],
    )
    # Revalidate every bound raw receipt/archive after the calculation.
    load_cohort(
        Path(job["folder"]),
        Path(job["audit_receipt"]),
        Path(job["audit_binding"]),
        job["binding_sha256"],
        stage,
    )
    if snapshot() != frozen or sha(job_path) != job_digest:
        raise RuntimeError("worker source/runtime or job changed during execution")
    if stage == "final":
        calibration_value = verify_calibration(
            Path(job["calibration"]), job["calibration_sha256"], frozen
        )
        verify_calibration_audit(
            Path(job["calibration_audit"]),
            job["calibration_audit_sha256"],
            job["calibration_sha256"],
            calibration_value["binding_sha256"],
        )
    if result["status"] != "passed":
        raise RuntimeError("worker case gates failed")
    write_json(
        output.with_name(output.stem + "-worker.json"),
        {
            "schema": "kernelgecp.cthyb-spectroscopy-worker.v1",
            "status": "passed",
            "pid": os.getpid(),
            "job_sha256": job_digest,
            "result_sha256": sha(output),
            "array_sha256": sha(output.with_suffix(".npz")),
        },
    )


def write_progress(path, result):
    temporary = path.with_suffix(".tmp")
    write_json(temporary, result)
    os.replace(temporary, path)


def campaign(
    folder,
    receipt,
    binding,
    binding_digest,
    stage,
    output,
    *,
    calibration=None,
    calibration_digest=None,
    calibration_audit=None,
    calibration_audit_digest=None,
):
    started = time.monotonic()
    validate_basis_controls()
    if stage not in ("calibration", "final"):
        raise ValueError("invalid campaign stage")
    frozen = snapshot()
    if stage == "final":
        calibration_value = verify_calibration(calibration, calibration_digest, frozen)
        verify_calibration_audit(
            calibration_audit,
            calibration_audit_digest,
            calibration_digest,
            calibration_value["binding_sha256"],
        )
    _, cohort = load_cohort(folder, receipt, binding, binding_digest, stage)
    if stage == "final":
        validate_cohort_pair(calibration_value["cohort"], cohort)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError("campaign output directory must be empty")
    output = output.resolve()
    target = output / ("calibration.json" if stage == "calibration" else "final.json")
    progress_path = output / "progress.json"
    cases, records, workers = {}, {}, {}
    progress = {
        "schema": "kernelgecp.cthyb-spectroscopy-progress.v1",
        "status": "running",
        "stage": stage,
        "parent_pid": os.getpid(),
        "protocol": PROTOCOL,
        "snapshot": frozen,
        "binding_sha256": binding_digest,
        "cohort": cohort,
        "calibration_cohort": calibration_value["cohort"] if stage == "final" else None,
        "cases": cases,
        "workers": workers,
        "global_timeout_seconds": GLOBAL_SECONDS,
        "case_timeout_seconds": CASE_SECONDS,
    }
    write_progress(progress_path, progress)
    repeats = PROTOCOL[stage + "_repeats"]
    try:
        for repeat in range(repeats):
            for method in METHODS:
                for model in PROTOCOL["default_models"]:
                    name = f"{method}-{model}-r{repeat}.json"
                    case = output / name
                    job_path = output / (case.stem + "-job.json")
                    log_path = output / (case.stem + "-worker.log")
                    job = {
                        "schema": "kernelgecp.cthyb-spectroscopy-worker-job.v1",
                        "stage": stage,
                        "parent_pid": os.getpid(),
                        "method": method,
                        "model": model,
                        "repeat": repeat,
                        "protocol": PROTOCOL,
                        "snapshot": frozen,
                        "folder": str(folder.resolve()),
                        "audit_receipt": str(receipt.resolve()),
                        "audit_binding": str(binding.resolve()),
                        "binding_sha256": binding_digest,
                        "output": str(case),
                        "calibration": str(calibration.resolve())
                        if calibration
                        else None,
                        "calibration_sha256": calibration_digest,
                        "calibration_audit": str(calibration_audit.resolve())
                        if calibration_audit
                        else None,
                        "calibration_audit_sha256": calibration_audit_digest,
                    }
                    # Empty output and exclusive job creation forbid retries.
                    with job_path.open("x") as stream:
                        json.dump(
                            job, stream, sort_keys=True, indent=2, allow_nan=False
                        )
                        stream.write("\n")
                    job_digest = sha(job_path)
                    command = [
                        sys.executable,
                        "-m",
                        "experiments.cthyb_spectroscopy",
                        "--worker",
                        str(job_path),
                        "--job-sha256",
                        job_digest,
                    ]
                    remaining = GLOBAL_SECONDS - (time.monotonic() - started)
                    if remaining <= 0:
                        raise TimeoutError("campaign global time budget exhausted")
                    workers[name] = {
                        "status": "starting",
                        "command": command,
                        "job_file": job_path.name,
                        "job_sha256": job_digest,
                        "log_file": log_path.name,
                    }
                    write_progress(progress_path, progress)

                    def mark_started(pid, case_name=name):
                        workers[case_name].update(pid=pid, status="running")
                        write_progress(progress_path, progress)

                    child = launch_worker(
                        command, log_path, min(CASE_SECONDS, remaining), mark_started
                    )
                    workers[name].update(child)
                    write_progress(progress_path, progress)
                    if child["status"] != "passed":
                        raise RuntimeError(
                            f"worker {name} {child['status']}, "
                            f"exit {child['exit_code']}"
                        )
                    if any(
                        other.get("pid") == child["pid"]
                        for key, other in workers.items()
                        if key != name
                    ):
                        raise RuntimeError("worker process ID unexpectedly reused")
                    worker_receipt = case.with_name(case.stem + "-worker.json")
                    finished = json.loads(worker_receipt.read_text())
                    if (
                        finished.get("status") != "passed"
                        or finished.get("pid") != child["pid"]
                        or finished.get("job_sha256") != job_digest
                        or finished.get("result_sha256") != sha(case)
                        or finished.get("array_sha256") != sha(case.with_suffix(".npz"))
                    ):
                        raise ValueError("worker receipt hash or process mismatch")
                    workers[name]["worker_receipt_file"] = worker_receipt.name
                    workers[name]["worker_receipt_sha256"] = sha(worker_receipt)
                    records[name] = json.loads(case.read_text())
                    cases[name] = {
                        "result_sha256": sha(case),
                        "array_sha256": sha(case.with_suffix(".npz")),
                    }
                    write_progress(progress_path, progress)
    except BaseException as error:
        progress.update(status="failed", error=f"{type(error).__name__}: {error}")
        write_progress(progress_path, progress)
        raise
    try:
        if snapshot() != frozen:
            raise RuntimeError("scientific source or runtime changed during campaign")
        if time.monotonic() - started > GLOBAL_SECONDS:
            raise TimeoutError("campaign global time budget exhausted")
    except BaseException as error:
        progress.update(status="failed", error=f"{type(error).__name__}: {error}")
        write_progress(progress_path, progress)
        raise
    sensitivities = {}
    for repeat in range(repeats):
        for method in METHODS:
            a = records[f"{method}-flat-r{repeat}.json"]
            b = records[f"{method}-broad_gaussian-r{repeat}.json"]
            sensitivities[f"{method}-r{repeat}"] = {
                "broad_minus_flat_window_weights": (
                    np.array(b["window_weights"]) - a["window_weights"]
                ).tolist(),
                "flat_bootstrap_window_std": a["bootstrap_window_standard_deviation"],
                "broad_bootstrap_window_std": b["bootstrap_window_standard_deviation"],
            }
    result = {
        "schema": "kernelgecp.cthyb-spectroscopy-"
        + ("calibration.v1" if stage == "calibration" else "final.v1"),
        "status": "passed"
        if all(x["status"] == "passed" for x in records.values())
        else "failed",
        "protocol": PROTOCOL,
        "snapshot": frozen,
        "binding_sha256": binding_digest,
        "cohort": cohort,
        "calibration_cohort": calibration_value["cohort"] if stage == "final" else None,
        "calibration_sha256": calibration_digest if stage == "final" else None,
        "calibration_audit_sha256": calibration_audit_digest
        if stage == "final"
        else None,
        "cases": cases,
        "workers": workers,
        "parent_pid": os.getpid(),
        "elapsed_seconds": time.monotonic() - started,
        "global_timeout_seconds": GLOBAL_SECONDS,
        "case_timeout_seconds": CASE_SECONDS,
        "sensitivity": sensitivities,
        "repeat_scope": (
            "fresh_solver_invocations_same_fixed_bath_cohort_not_independent_data"
        ),
        "claim_scope": (
            "conditional_heldout_measurement_noise_not_total_prediction_"
            "uncertainty_or_spectral_truth"
        ),
    }
    write_json(target, result)
    progress["status"] = result["status"]
    progress["campaign_sha256"] = sha(target)
    write_progress(progress_path, progress)
    return result


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        parser = argparse.ArgumentParser(
            description="Execute one frozen campaign worker"
        )
        parser.add_argument("--worker", type=Path, required=True)
        parser.add_argument("--job-sha256", required=True)
        args = parser.parse_args()
        worker(args.worker, args.job_sha256)
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--audit-receipt", type=Path, required=True)
    p.add_argument("--audit-binding", type=Path, required=True)
    p.add_argument("--binding-sha256", required=True)
    p.add_argument("--stage", choices=["calibration", "final"], required=True)
    p.add_argument("--calibration", type=Path)
    p.add_argument("--calibration-sha256")
    p.add_argument("--calibration-audit", type=Path)
    p.add_argument("--calibration-audit-sha256")
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    result = campaign(
        args.folder,
        args.audit_receipt,
        args.audit_binding,
        args.binding_sha256,
        args.stage,
        args.output_dir,
        calibration=args.calibration,
        calibration_digest=args.calibration_sha256,
        calibration_audit=args.calibration_audit,
        calibration_audit_digest=args.calibration_audit_sha256,
    )
    if result["status"] != "passed":
        raise SystemExit("campaign gates failed; artifacts retained")


if __name__ == "__main__":
    main()
