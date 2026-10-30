"""Actual TauMaxEnt pilot/campaign; no recovery surrogate or hidden covariance."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import expit, ndtr

from experiments.dlr_adapters import build_adapter
from kernelgecp import FermionicBasis
from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    continue_tau_maxent,
    covariance_modes,
    prepare_spectral_projection,
    save_spectral_dataset,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    "experiments/spectroscopy_campaign.py",
    "experiments/audit_spectroscopy.py",
    "experiments/spectroscopy.py",
    "experiments/dlr_adapters.py",
    "kernelgecp/integrations/maxent.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "environments/locks/triqs-linux-64-conda-explicit.txt",
    "environments/locks/triqs-maxent-pip.txt",
    "environments/locks/dlr-physics-pip.txt",
)
METHODS = ("raw", "kernelgecp", "pydlr", "sparse_ir")
CASES = ("offgrid_atoms_independent", "offgrid_atoms_correlated", "smooth_correlated")
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def provenance():
    return {
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_hashes": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in SOURCE_FILES
        },
        "tracked_clean": not subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=ROOT,
            text=True,
        ).strip(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "pydlr", "sparse-ir")
        },
        "threads": {name: os.environ.get(name) for name in THREADS},
    }


def kernel(beta, tau, omega):
    exponent = np.where(
        omega >= 0, -tau[:, None] * omega, (beta - tau[:, None]) * omega
    )
    return -np.exp(exponent) * expit(beta * np.abs(omega))


def truth(case, beta, tau):
    if case.startswith("offgrid_atoms"):
        poles = np.array([-1.371, 0.413, 2.117])
        masses = np.array([0.2, 0.55, 0.25])
        return kernel(beta, tau, poles) @ masses
    nodes, weights = np.polynomial.legendre.leggauss(512)
    omega = nodes * 4
    density = smooth_density(omega)
    return kernel(beta, tau, omega) @ (weights * 4 * density)


def smooth_density(omega):
    omega = np.asarray(omega, dtype=np.float64)
    means, scales, weights = [-1.1, 0.55, 1.9], [0.32, 0.45, 0.25], [0.2, 0.55, 0.25]
    density = np.zeros_like(omega)
    for mean, sigma, mass in zip(means, scales, weights, strict=True):
        norm = ndtr((4 - mean) / sigma) - ndtr((-4 - mean) / sigma)
        density += (
            mass
            * np.exp(-0.5 * ((omega - mean) / sigma) ** 2)
            / (np.sqrt(2 * np.pi) * sigma * norm)
        )
    return density


def true_cdf(case, omega):
    omega = np.asarray(omega, dtype=np.float64)
    if case.startswith("offgrid_atoms"):
        return sum(
            mass * (omega >= pole)
            for pole, mass in zip(
                [-1.371, 0.413, 2.117], [0.2, 0.55, 0.25], strict=True
            )
        )
    result = np.zeros_like(omega)
    for mean, sigma, mass in zip(
        [-1.1, 0.55, 1.9], [0.32, 0.45, 0.25], [0.2, 0.55, 0.25], strict=True
    ):
        lower = ndtr((-4 - mean) / sigma)
        norm = ndtr((4 - mean) / sigma) - lower
        result += mass * (ndtr((omega - mean) / sigma) - lower) / norm
    return result


def observations(case, seed, *, query_order=24):
    if type(query_order) is not int or query_order not in (24, 36, 48):
        raise ValueError(
            "query_order must be 24, historical 36, or predeclared fresh 48"
        )
    beta = 10.0
    tau = np.linspace(0, beta, 65)
    # The final query nodes are off the observation mesh and do not select alpha.
    query = (np.polynomial.legendre.leggauss(query_order)[0] + 1) * beta / 2
    sigma = 1e-4
    correlation = np.eye(tau.size)
    if case.endswith("correlated"):
        correlation = 0.25 * correlation + 0.75 * np.exp(
            -np.abs(tau[:, None] - tau) / 0.6
        )
    covariance = sigma**2 * correlation
    generator = np.random.default_rng(seed)
    values = truth(case, beta, tau) + np.linalg.cholesky(
        covariance
    ) @ generator.standard_normal(tau.size)
    query_covariance = sigma**2 * np.eye(query.size)
    if case.endswith("correlated"):
        query_covariance = sigma**2 * (
            0.25 * np.eye(query.size)
            + 0.75 * np.exp(-np.abs(query[:, None] - query) / 0.6)
        )
    return (
        SpectralDataset(
            beta,
            tau,
            values,
            covariance,
            provenance={
                "provider": "analytic_synthetic",
                "case": case,
                "seed": seed,
                "noise_sigma": sigma,
            },
        ),
        query,
        query_covariance,
    )


def basis(method, beta, eps):
    if method == "kernelgecp":
        return FermionicBasis(
            beta, 4, eps, precision="float64", grid_order=24, max_rank=96
        )
    return build_adapter(method, beta, 4, eps)


def spectrum_metrics(case, result):
    omega, density = result.omega, result.density
    cdf = np.r_[0.0, np.cumsum(np.diff(omega) * (density[:-1] + density[1:]) / 2)]
    actual = true_cdf(case, omega)
    windows = np.interp([-4, -1, 1, 4], omega, cdf)
    oracle = true_cdf(case, np.array([-4, -1, 1, 4]))
    return {
        "window_weights": np.diff(windows).tolist(),
        "window_weight_errors": (np.diff(windows) - np.diff(oracle)).tolist(),
        "cdf_max_error": float(np.max(np.abs(cdf - actual))),
        "cdf_integrated_absolute_error": float(
            np.trapezoid(np.abs(cdf - actual), omega)
        ),
    }


def run_case(
    case, method, seed, eps, config, output, *, pilot, bootstrap=0, query_order=24
):
    before = provenance()
    if not pilot and (
        not before["tracked_clean"] or any(os.environ.get(v) != "1" for v in THREADS)
    ):
        raise ValueError(
            "canonical spectroscopy requires clean source and fixed single threads"
        )
    data, query, query_covariance = observations(case, seed, query_order=query_order)
    start = time.perf_counter()
    projection = None
    coefficient_fit = None
    observation_design = np.eye(data.tau.size)
    fit_projection_difference = None
    fitted = data
    if method != "raw":
        constructed = basis(method, data.beta, eps)
        observation_design = constructed.design_tau(data.tau)
        projection = prepare_spectral_projection(
            constructed,
            data.tau,
            data.tau,
            rcond=1e-10,
            metadata={"method": method, "eps": eps},
        )
        fitted = projection.apply(data)
        if method == "kernelgecp":
            coefficient_fit = constructed.fit_tau(data.values, data.tau, rcond=1e-10)
            fit_projection_difference = float(
                np.max(np.abs(coefficient_fit.evaluate_tau(data.tau) - fitted.values))
            )
    construction_seconds = time.perf_counter() - start
    config = MaxEntConfig(**config)
    output.parent.mkdir(parents=True, exist_ok=True)
    save_spectral_dataset(
        fitted,
        output.with_name(output.stem + "-dataset.npz"),
        units={"tau": "inverse_energy", "omega": "energy", "G_tau": "dimensionless"},
        provider={"name": method, "continuation": "TRIQS/MaxEnt4.0.1"},
        calculation_settings=dataclasses.asdict(config),
    )
    fit_start = time.perf_counter()
    spectrum = continue_tau_maxent(fitted, config)
    continuation_seconds = time.perf_counter() - fit_start
    prediction_start = time.perf_counter()
    predicted = spectrum.evaluate_tau(data.beta, query)
    prediction_seconds = time.perf_counter() - prediction_start
    reference = truth(case, data.beta, query)
    _, modes = covariance_modes(query_covariance)
    rms = modes.whitened_rms(predicted - reference)
    metrics = spectrum_metrics(case, spectrum)
    bootstrap_weights, bootstrap_cdf, bootstrap_spectra = [], [], []
    bootstrap_values, bootstrap_records = [], []
    generator = np.random.default_rng(seed + 100000)
    for _ in range(bootstrap):
        sampled = SpectralDataset(
            data.beta,
            data.tau,
            data.values
            + np.linalg.cholesky(data.covariance)
            @ generator.standard_normal(data.tau.size),
            data.covariance,
            provenance=data.provenance,
        )
        bootstrap_values.append(sampled.values)
        if projection is not None:
            sampled = projection.apply(sampled)
        resampled = continue_tau_maxent(sampled, config)
        diagnostic = spectrum_metrics(case, resampled)
        bootstrap_weights.append(diagnostic["window_weights"])
        bootstrap_cdf.append(diagnostic["cdf_max_error"])
        bootstrap_spectra.append(resampled.density)
        bootstrap_records.append(
            {
                "selected_index": resampled.selected_index,
                "raw_mass": resampled.raw_mass,
                "mass": resampled.mass,
                "converged": resampled.converged,
                "iterations": resampled.iterations,
                "effective_rank": resampled.effective_rank,
            }
        )
    if bootstrap:
        metrics["bootstrap_window_standard_deviation"] = np.std(
            bootstrap_weights, axis=0, ddof=1
        ).tolist()
        metrics["bootstrap_cdf_error_range"] = [
            float(min(bootstrap_cdf)),
            float(max(bootstrap_cdf)),
        ]
    after = provenance()
    if before != after:
        raise RuntimeError(
            "spectroscopy source or environment changed during generation"
        )
    gates = {
        "positive": bool(np.all(spectrum.density >= 0)),
        "normalized": abs(spectrum.mass - 1) <= 1e-6,
        "heldout_whitened_rms": rms <= 2,
        "provider_selected_converged": spectrum.converged[spectrum.selected_index],
        "coefficient_projection_agreement": fit_projection_difference is None
        or fit_projection_difference <= 1e-10,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output.with_suffix(".npz"),
        omega=spectrum.omega,
        density=spectrum.density,
        tau=data.tau,
        values=data.values,
        covariance=data.covariance,
        projected_values=fitted.values,
        projected_covariance=fitted.covariance,
        projection=np.eye(data.tau.size) if projection is None else projection.matrix,
        observation_design=observation_design,
        query=query,
        query_covariance=query_covariance,
        prediction=predicted,
        reference=reference,
        alphas=spectrum.alphas,
        chi2=spectrum.chi2,
        trial_spectra=spectrum.trial_spectra,
        bootstrap_weights=np.array(bootstrap_weights),
        bootstrap_spectra=np.array(bootstrap_spectra),
        bootstrap_values=np.array(bootstrap_values),
    )
    coefficient_archive = None
    if coefficient_fit is not None:
        archive_path = output.with_name(output.stem + "-gecp-fit.npz")
        coefficient_fit.save(archive_path)
        coefficient_archive = {
            "path": archive_path.name,
            "sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest(),
            "format": "kernelgecp.fit.v2",
            "observation_domain": "tau",
        }
    result = {
        "schema": "kernelgecp.spectroscopy-case.v1",
        "pilot": pilot,
        "status": "passed" if all(gates.values()) else "failed",
        "gates": gates,
        "case": case,
        "method": method,
        "seed": seed,
        "query_order": query_order,
        "eps": eps,
        "config": dataclasses.asdict(config),
        "provenance": before,
        "effective_rank": fitted.modes.effective_rank,
        "whitened_rms": rms,
        "raw_mass": spectrum.raw_mass,
        "mass": spectrum.mass,
        "normalization_training_whitened_rms": (
            spectrum.normalized_training_whitened_rms
        ),
        "selected_index": spectrum.selected_index,
        "provider_diagnostics": spectrum.diagnostics,
        "converged": spectrum.converged,
        "iterations": spectrum.iterations,
        "metrics": metrics,
        "construction_seconds": construction_seconds,
        "continuation_seconds": continuation_seconds,
        "prediction_seconds": prediction_seconds,
        "workflow_before_final_receipt_seconds": time.perf_counter() - start,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        * (1 if sys.platform == "darwin" else 1024),
        "array_storage_bytes": output.with_suffix(".npz").stat().st_size,
        "projection_svd": None
        if projection is None
        else dataclasses.asdict(projection.diagnostics),
        "coefficient_fit_archive": coefficient_archive,
        "fit_projection_max_difference": fit_projection_difference,
        "bootstrap_count": bootstrap,
        "bootstrap_seed": seed + 100000,
        "bootstrap_records": bootstrap_records,
        "uncertainty_claim": "empirical_bootstrap_variation_not_complete_coverage",
        "array_sha256": hashlib.sha256(
            output.with_suffix(".npz").read_bytes()
        ).hexdigest(),
    }
    write_json(output, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--eps", type=float, default=1e-8)
    parser.add_argument("--n-omega", type=int, default=401)
    parser.add_argument("--n-alpha", type=int, default=33)
    parser.add_argument(
        "--default-model", choices=["flat", "broad_gaussian"], default="flat"
    )
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--query-order", type=int, choices=[24, 36, 48], default=24)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.bootstrap <= 64 or args.bootstrap == 1:
        parser.error("bootstrap count must be0 or in[2,64]")
    try:
        config = dataclasses.asdict(
            MaxEntConfig(
                n_omega=args.n_omega,
                n_alpha=args.n_alpha,
                default_model=args.default_model,
            )
        )
        result = run_case(
            args.case,
            args.method,
            args.seed,
            args.eps,
            config,
            args.output,
            pilot=args.pilot,
            bootstrap=args.bootstrap,
            query_order=args.query_order,
        )
        if result["status"] != "passed":
            raise RuntimeError("spectroscopy gates failed; retain failed artifacts")
    except Exception as error:
        if not args.output.exists():
            write_json(
                args.output,
                {"status": "failed", "error": f"{type(error).__name__}: {error}"},
            )
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
