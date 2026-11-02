#!/usr/bin/env python3
"""Diagnostic constrained MaxEnt scout on the consumed training observations only."""

from __future__ import annotations

import argparse
import contextlib
import ctypes
import hashlib
import importlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

FIXTURE = "tests/data/maxent/v2-98b-raw-failed/raw-flat-r0"
ARRAY_SHA256 = "1ce196afc1c642736ccec72487566eef2f075e73ae3120a6719a7cb2370ed02f"
JSON_SHA256 = "f8d13234fda36f4ebc009ad2680db0079839e4ad6be8f8a1c3fc734a2c6fa2fa"
FIXED = {
    "omega_min": -4.0,
    "omega_max": 4.0,
    "n_omega": 401,
    "alpha_min": 1e-4,
    "alpha_max": 1e4,
    "n_alpha": 33,
    "maxiter": 1000,
    "minimizer_mu0": 1e-8,
    "normalize_mass": 1.0,
    "mass_mode": "constrained",
    "covariance_rtol": 1e-4,
    "covariance_atol": 1e-18,
    "beta": 10.0,
    "training_points": 65,
    "expected_rank": 64,
    "models": ["flat", "broad_gaussian"],
    "quadrature_order": 12,
    "threads": 1,
    "provider_version": "4.0.1",
}
FILES = {
    "maxent": "kernelgecp/integrations/maxent.py",
    "mass": "kernelgecp/integrations/_maxent_mass.py",
    "native_mass": "kernelgecp/integrations/_maxent_native_mass.py",
}
TRACKED = [
    *FILES.values(),
    FIXTURE + ".json",
    FIXTURE + ".npz",
    "scripts/providers/installed_wheel.py",
    "scripts/providers/triqs_smoke.py",
    "scripts/providers/maxent_mass_smoke.py",
    ".github/workflows/maxent-science.yml",
    "environments/triqs.yml",
    "environments/locks/triqs-linux-64-conda-explicit.txt",
    "environments/locks/triqs-maxent-pip.txt",
    "scripts/providers/maxent_mass_training_scout.py",
]
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def load_helper(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load source helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bind_source(root, revision):
    hashes = {}
    for name in TRACKED:
        blob = subprocess.check_output(
            ["git", "-C", str(root), "show", revision + ":" + name]
        )
        if (root / name).read_bytes() != blob:
            raise RuntimeError("working source differs from revision: " + name)
        hashes[name] = hashlib.sha256(blob).hexdigest()
    if hashes[FIXTURE + ".npz"] != ARRAY_SHA256:
        raise RuntimeError("fixture archive hash mismatch")
    if hashes[FIXTURE + ".json"] != JSON_SHA256:
        raise RuntimeError("fixture JSON hash mismatch")
    if Path(__file__).read_bytes() != (root / TRACKED[-1]).read_bytes():
        raise RuntimeError("executed scout differs from bound source")
    return hashes


def load_training(np, archive):
    # Do not enumerate, deserialize, or access query/validation members.
    with np.load(archive, allow_pickle=False) as data:
        return tuple(
            data[name].copy()
            for name in ("observation_tau", "observation_values", "covariance")
        )


def forward12(np, beta, tau, omega, spectra):
    nodes, weights = np.polynomial.legendre.leggauss(12)
    widths = np.diff(omega)
    energies = (
        (omega[:-1, None] + omega[1:, None]) / 2 + widths[:, None] * nodes / 2
    ).ravel()
    frac = (nodes + 1) / 2
    amps = spectra[:, :-1, None] * (1 - frac) + spectra[:, 1:, None] * frac
    quadrature = (amps * widths[None, :, None] * weights / 2).reshape(len(spectra), -1)
    exponent = np.where(
        energies >= 0, -tau[:, None] * energies, (beta - tau[:, None]) * energies
    )
    kernel = -np.exp(exponent) / (1 + np.exp(-beta * np.abs(energies)))
    return quadrature @ kernel.T


@contextlib.contextmanager
def capture_fds(log):
    """Capture Python and native writes to both process file descriptors."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    try:
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        ctypes.CDLL(None).fflush(None)
        for target, fd in enumerate(saved, 1):
            os.dup2(fd, target)
            os.close(fd)


def output_hashes(directory):
    names = ["worker.log", "flat-trials.npz", "broad_gaussian-trials.npz"]
    return {
        name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
        for name in names
        if (directory / name).is_file()
    }


def validate_training(np, tau, values, covariance):
    if (
        tau.shape != (65,)
        or values.shape != (65,)
        or covariance.shape != (65, 65)
        or not np.isfinite(tau).all()
        or not np.isfinite(values).all()
        or not np.isfinite(covariance).all()
        or tau[0] != 0.0
        or tau[-1] != 10.0
        or np.any(np.diff(tau) <= 0)
    ):
        raise RuntimeError("invalid finite 65-point beta=10 training observations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    receipt = {
        "schema": "kernelgecp.maxent-training-scout.v1",
        "status": "failed",
        "scope": (
            "consumed development data; training-only diagnostic; no heldout scoring "
            "or scientific acceptance"
        ),
        "fixed_controls": FIXED,
        "cases": {},
    }
    stage = "source"
    with (
        (args.output_dir / "worker.log").open("x") as log,
        capture_fds(log),
    ):
        try:
            root = args.source_root.resolve()
            for name in THREADS:
                os.environ[name] = "1"
            triqs = load_helper(
                root / "scripts/providers/triqs_smoke.py", "_scout_triqs"
            )
            wheel = load_helper(
                root / "scripts/providers/installed_wheel.py", "_scout_wheel"
            )
            smoke = load_helper(
                root / "scripts/providers/maxent_mass_smoke.py", "_scout_smoke"
            )
            identity = triqs.source_identity(root)
            if identity["source_revision"] != args.source_revision:
                raise RuntimeError("source revision mismatch")
            hashes = bind_source(root, args.source_revision)
            metadata = json.loads((root / (FIXTURE + ".json")).read_text())
            if metadata.get("array_sha256") != ARRAY_SHA256:
                raise RuntimeError("fixture metadata hash mismatch")
            receipt.update(
                source_revision=args.source_revision,
                source_sha256=hashes,
                identity=identity,
                environment=triqs.environment_identity(Path(sys.prefix)),
            )
            stage = "imports"
            np = importlib.import_module("numpy")
            pools = importlib.import_module("threadpoolctl")
            importlib.import_module("triqs_maxent")
            version = str(importlib.import_module("triqs_maxent.version").version)
            if version != "4.0.1":
                raise RuntimeError("actual MaxEnt 4.0.1 required")
            modules = {
                k: importlib.import_module(v[:-3].replace("/", "."))
                for k, v in FILES.items()
            }
            paths = {k: Path(v.__file__).resolve() for k, v in modules.items()}
            for key, path in paths.items():
                if (
                    path.is_relative_to(root)
                    or hashlib.sha256(path.read_bytes()).hexdigest()
                    != hashes[FILES[key]]
                ):
                    raise RuntimeError("installed source binding failed: " + key)
            receipt["installation"] = wheel.installed_wheel_identity(paths, FILES)
            smoke.verify_record_hashes(paths)
            receipt["provider_version"] = version
            tau, values, covariance = load_training(np, root / (FIXTURE + ".npz"))
            validate_training(np, tau, values, covariance)
            api = modules["maxent"]
            data = api.SpectralDataset(
                float(tau[-1]), tau, values, covariance, covariance_rtol=1e-4
            )
            if tau.shape != (65,) or data.modes.effective_rank != 64:
                raise RuntimeError("training size or retained covariance rank mismatch")
            receipt["training"] = {
                "points": 65,
                "effective_rank": 64,
                "covariance_threshold": data.modes.threshold,
                "covariance_eigenvalues": data.modes.eigenvalues.tolist(),
            }
            stage = "continuation"
            with pools.threadpool_limits(limits=1):
                if not pools.threadpool_info() or any(
                    x["num_threads"] != 1 for x in pools.threadpool_info()
                ):
                    raise RuntimeError("single-thread limit unavailable")
                for model in ("flat", "broad_gaussian"):
                    try:
                        config = api.MaxEntConfig(
                            n_omega=401,
                            n_alpha=33,
                            maxiter=1000,
                            minimizer_mu0=1e-8,
                            default_model=model,
                            normalize_mass=1.0,
                            mass_mode="constrained",
                        )
                        result = api.continue_tau_maxent(data, config)
                        predictions = forward12(
                            np, data.beta, tau, result.omega, result.trial_spectra
                        )
                        whitened = [
                            float(data.modes.whitened_rms(row - values))
                            for row in predictions
                        ]
                        np.savez_compressed(
                            args.output_dir / (model + "-trials.npz"),
                            omega=result.omega,
                            trial_spectra=result.trial_spectra,
                            alphas=result.alphas,
                            chi2=result.chi2,
                            training_prediction=predictions,
                        )
                        receipt["cases"][model] = {
                            "status": "completed",
                            "selected_index": result.selected_index,
                            "selected_alpha": float(
                                result.alphas[result.selected_index]
                            ),
                            "converged": list(result.converged),
                            "iterations": list(result.iterations),
                            "masses": np.trapezoid(
                                result.trial_spectra, result.omega, axis=1
                            ).tolist(),
                            "training_whitened_rms_12pt": whitened,
                            "diagnostics": result.diagnostics,
                        }
                    except Exception as exc:
                        receipt["cases"][model] = {
                            "status": "failed",
                            "type": type(exc).__name__,
                            "message": str(exc),
                            "partial_trial_limit": (
                                "wrapper reports completed history in error; "
                                "partial spectra unavailable"
                            ),
                        }
                        print(f"{model} failed: {type(exc).__name__}: {exc}")
            receipt["status"] = (
                "diagnostic_completed"
                if all(v["status"] == "completed" for v in receipt["cases"].values())
                else "failed"
            )
        except Exception as exc:
            receipt["error"] = {
                "stage": stage,
                "type": type(exc).__name__,
                "message": str(exc),
            }
            print(f"scout failed at {stage}: {type(exc).__name__}: {exc}")
    receipt["output_sha256"] = output_hashes(args.output_dir)
    (args.output_dir / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return 0 if receipt["status"] == "diagnostic_completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
