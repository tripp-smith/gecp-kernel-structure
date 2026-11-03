#!/usr/bin/env python3
"""Independent, training-only analysis of an authenticated scout receipt.

Does not run a provider, optimize, or access any held-out arrays.  Quadrature
error estimates are not certified bounds; installed identity is self-reported
by the scout and cannot be independently reconstructed on another machine.
"""

import argparse
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
from pathlib import Path

for _thread in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread] = "1"

import numpy as np  # noqa: E402
from scipy.integrate import quad_vec  # noqa: E402
from scipy.linalg import eigh  # noqa: E402
from scipy.special import expit  # noqa: E402

INPUT_SHA = "1ce196afc1c642736ccec72487566eef2f075e73ae3120a6719a7cb2370ed02f"
FIXTURE = "tests/data/maxent/v2-98b-raw-failed/raw-flat-r0.npz"
MODELS = ("flat", "broad_gaussian")
MODULES = (
    "kernelgecp/integrations/maxent.py",
    "kernelgecp/integrations/_maxent_mass.py",
    "kernelgecp/integrations/_maxent_native_mass.py",
)
SOURCE_INVENTORY = frozenset(
    {
        *MODULES,
        FIXTURE,
        "tests/data/maxent/v2-98b-raw-failed/raw-flat-r0.json",
        ".github/workflows/maxent-science.yml",
        "environments/locks/triqs-linux-64-conda-explicit.txt",
        "environments/locks/triqs-maxent-pip.txt",
        "environments/triqs.yml",
        "scripts/providers/installed_wheel.py",
        "scripts/providers/maxent_mass_smoke.py",
        "scripts/providers/maxent_mass_training_scout.py",
        "scripts/providers/maxent_solve_capture.py",
        "scripts/providers/triqs_smoke.py",
    }
)
CONTROLS = {
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
    "mass_curvature": "projected",
    "damping": {"identity": False, "marquardt": True},
    "mass_interior_fraction": 1e-12,
    "covariance_rtol": 1e-4,
    "covariance_atol": 1e-18,
    "beta": 10.0,
    "training_points": 65,
    "expected_rank": 64,
    "models": list(MODELS),
    "quadrature_order": 12,
    "threads": 1,
    "provider_version": "4.0.1",
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def authenticate_receipt(path, expected_sha):
    require(
        not path.is_symlink() and digest(path) == expected_sha,
        "receipt authentication failed",
    )
    return json.loads(path.read_text())


def validate_receipt_header(receipt, revision):
    require(receipt.get("schema") == "kernelgecp.maxent-training-scout.v4", "schema")
    require(
        receipt.get("source_revision") == revision
        and receipt.get("identity", {}).get("source_revision") == revision,
        "source revision",
    )
    require(receipt.get("fixed_controls") == CONTROLS, "fixed controls")


def validate_outputs(directory, hashes):
    expected = {
        "worker.log",
        *(
            d + "-" + m + "-trials.npz"
            for d in ("identity", "marquardt")
            for m in MODELS
        ),
    }
    require(isinstance(hashes, dict) and set(hashes) == expected, "output inventory")
    for name, sha in hashes.items():
        path = directory / name
        require(not path.is_symlink() and digest(path) == sha, "output hash: " + name)


def validate_minimizer_flags(config, expected, native, damping):
    require(
        config == expected
        and native.get("marquardt") is (damping == "marquardt")
        and native.get("J_squared") is False
        and native.get("mu0") == 1e-8,
        "row configuration or native minimizer flags",
    )


def validate_selection(converged, selected, independent_selected):
    require(
        isinstance(converged, list)
        and len(converged) == 33
        and all(type(x) is bool for x in converged),
        "convergence flags",
    )
    require(
        type(selected) is int and 0 <= selected < 33 and converged[selected],
        "selected convergence",
    )
    require(
        selected == independent_selected, "independent training selector disagreement"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", required=True, type=Path)
    parser.add_argument(
        "--receipt-sha256",
        required=True,
        help="authenticated receipt hash, obtained independently",
    )
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(args.artifact_dir.resolve()):
        parser.error("output must be outside the immutable artifact directory")
    output = {
        "status": "inconclusive",
        "cases": {},
        "scope": "consumed training only",
        "quadrature_error_status": "estimated_not_certified",
        "installed_identity_scope": (
            "receipt claims; remote installation not independently observable"
        ),
    }
    try:
        require(
            re.fullmatch(r"[0-9a-f]{40}", args.source_revision) is not None,
            "source revision must be a full lowercase commit SHA",
        )
        root, directory = args.source_root.resolve(), args.artifact_dir.resolve()
        require(not args.output.exists(), "output must not already exist")
        receipt_path = directory / "receipt.json"
        receipt = authenticate_receipt(receipt_path, args.receipt_sha256)
        require(
            receipt.get("schema") == "kernelgecp.maxent-training-scout.v4", "schema"
        )
        output["scout_status"] = receipt.get("status")
        output["receipt_sha256"] = args.receipt_sha256
        if receipt.get("status") != "diagnostic_completed":
            output["reason"] = "scout failed or incomplete; no scientific pass inferred"
            output["scout_error"] = receipt.get("error")
            output["scout_cases"] = receipt.get("cases", {})
            return output
        validate_receipt_header(receipt, args.source_revision)
        source = receipt["source_sha256"]
        require(
            isinstance(source, dict) and set(source) == SOURCE_INVENTORY,
            "source inventory",
        )
        for name, sha in source.items():
            require(
                isinstance(name, str)
                and not name.startswith("/")
                and ".." not in Path(name).parts
                and isinstance(sha, str)
                and re.fullmatch(r"[0-9a-f]{64}", sha),
                "invalid source inventory entry",
            )
            blob = subprocess.check_output(
                ["git", "-C", str(root), "show", args.source_revision + ":" + name]
            )
            require(
                hashlib.sha256(blob).hexdigest() == sha, "source blob mismatch: " + name
            )
        require(all(name in source for name in (*MODULES, FIXTURE)), "source inventory")
        require(source[FIXTURE] == INPUT_SHA, "fixture source hash")
        require(receipt.get("provider_version") == "4.0.1", "provider version")
        install = receipt.get("installation", {})
        require(
            install.get("distribution") == "kernelgecp"
            and install.get("editable") is False
            and install.get("module_count") == 3
            and len(install.get("record_sha256", "")) == 64,
            "installed wheel identity receipt",
        )
        env = receipt.get("environment", {})
        require(
            env.get("conda_package_count", 0) > 0
            and len(env.get("conda_package_set_sha256", "")) == 64
            and all(x == "1" for x in env.get("thread_limits", {}).values())
            and len(env.get("thread_limits", {})) >= 5,
            "environment identity receipt",
        )
        validate_outputs(directory, receipt.get("output_sha256"))
        fixture = root / FIXTURE
        require(digest(fixture) == INPUT_SHA, "fixture hash")
        # Only these three fixture members are accessed; never enumerate the archive.
        with np.load(fixture, allow_pickle=False) as z:
            tau, y, covariance = (
                z[k].copy()
                for k in ("observation_tau", "observation_values", "covariance")
            )
        require(
            tau.shape == y.shape == (65,) and covariance.shape == (65, 65), "shapes"
        )
        require(
            np.isfinite(tau).all()
            and np.isfinite(y).all()
            and np.isfinite(covariance).all(),
            "nonfinite training data",
        )
        eig, vec = eigh((covariance + covariance.T) / 2, driver="evd")
        threshold = max(1e-18, 1e-4 * eig[-1])
        keep = eig > threshold
        require(keep.sum() == 64, "training covariance rank")
        whitening = (vec[:, keep] / np.sqrt(eig[keep])).T
        output["covariance"] = {"rank": int(keep.sum()), "threshold": float(threshold)}
        expected_cases = {
            d + "/" + m for d in ("identity", "marquardt") for m in MODELS
        }
        require(set(receipt["cases"]) == expected_cases, "case inventory")

        def deadline(signum, frame):
            raise TimeoutError("training quadrature time budget exhausted")

        signal.signal(signal.SIGALRM, deadline)
        signal.alarm(120)
        for case in sorted(expected_cases):
            damping, model = case.split("/")
            record = receipt["cases"][case]
            require(record["status"] == "completed", "case incomplete")
            diagnostic = record["diagnostics"]
            parameterization = diagnostic.get("mass_parameterization", {})
            config = diagnostic.get("config", {})
            native = diagnostic.get("native_minimizer", {})
            expected_config = {
                k: CONTROLS[k]
                for k in (
                    "omega_min",
                    "omega_max",
                    "n_omega",
                    "alpha_min",
                    "alpha_max",
                    "n_alpha",
                    "maxiter",
                    "minimizer_mu0",
                    "normalize_mass",
                    "mass_mode",
                    "mass_curvature",
                    "mass_interior_fraction",
                )
            }
            expected_config.update(
                default_model=model, minimizer_marquardt=damping == "marquardt"
            )
            validate_minimizer_flags(config, expected_config, native, damping)
            require(
                diagnostic.get("mass_operation") == "native_mass_constraint"
                and parameterization.get("curvature") == "projected"
                and parameterization.get("d_dv") is False
                and parameterization.get("dA_projection") == 2
                and parameterization.get("chi2_factor") == 1.0
                and parameterization.get("interior_fraction") == 1e-12
                and parameterization.get("floor_policy")
                == "reject_native_entropy_floor_or_overflow"
                and diagnostic.get("selection_rule")
                == "closest_abs_log_chi2_per_retained_rank_to_one",
                "projected curvature, objective, or mass operation mismatch",
            )
            observation = parameterization.get("last_evaluation", {})
            for key in ("min_observed_H", "min_observed_mixture_floor"):
                value = parameterization.get(key)
                require(
                    isinstance(value, (int, float))
                    and np.isfinite(value)
                    and value > 1e-100,
                    "native floor diagnostic: " + key,
                )
            for key in ("min_H", "actual_mixture_floor"):
                value = observation.get(key)
                require(
                    isinstance(value, (int, float))
                    and np.isfinite(value)
                    and value > 1e-100,
                    "last native floor diagnostic: " + key,
                )
            maximum_zeros = parameterization.get("max_latent_zero_count")
            require(
                type(maximum_zeros) is int
                and 0 <= maximum_zeros <= 401
                and type(parameterization.get("evaluation_count")) is int
                and parameterization["evaluation_count"] > 0,
                "latent zero/evaluation diagnostic",
            )
            with np.load(
                directory / (damping + "-" + model + "-trials.npz"), allow_pickle=False
            ) as z:
                omega, spectra, alphas, producer_prediction, native_chi = (
                    z[k].copy()
                    for k in (
                        "omega",
                        "trial_spectra",
                        "alphas",
                        "training_prediction",
                        "chi2",
                    )
                )
            require(
                omega.shape == (401,)
                and spectra.shape == (33, 401)
                and alphas.shape == native_chi.shape == (33,)
                and producer_prediction.shape == (33, 65),
                "trial shapes",
            )
            require(
                np.allclose(omega, np.linspace(-4, 4, 401), rtol=0, atol=1e-13)
                and np.isfinite(spectra).all()
                and (spectra > 0).all(),
                "trial grid or positivity",
            )
            masses = np.trapezoid(spectra, omega, axis=1)
            require(np.allclose(masses, 1, rtol=0, atol=1e-10), "trial mass")
            widths = np.empty(401)
            widths[0], widths[-1] = (
                (omega[1] - omega[0]) / 2,
                (omega[-1] - omega[-2]) / 2,
            )
            widths[1:-1] = (omega[2:] - omega[:-2]) / 2
            exponent = np.where(
                omega[None, :] >= 0,
                -tau[:, None] * omega[None, :],
                (10 - tau[:, None]) * omega[None, :],
            )
            kernel = -np.exp(exponent) * expit(10 * np.abs(omega))[None, :]
            # NormalChi2 receives K.K @ H, where H is density times mesh widths.
            native_residual = ((spectra * widths) @ kernel.T - y) @ whitening.T
            native_recomputed = np.sum(native_residual * native_residual, axis=1)
            require(
                np.isfinite(native_chi).all()
                and np.allclose(native_chi, native_recomputed, rtol=3e-5, atol=2e-5),
                "native K.K on H chi-square disagreement",
            )
            require(
                isinstance(record["converged"], list)
                and len(record["converged"]) == 33
                and all(type(x) is bool for x in record["converged"]),
                "convergence flags",
            )
            converged = np.asarray(record["converged"], dtype=bool)
            selected = record["selected_index"]
            require(
                converged.shape == (33,)
                and type(selected) is int
                and 0 <= selected < 33
                and converged[selected],
                "selected convergence",
            )
            require(
                np.isfinite(alphas).all()
                and np.allclose(
                    np.sort(alphas), np.geomspace(1e-4, 1e4, 33), rtol=1e-12, atol=0
                )
                and np.isclose(record["selected_alpha"], alphas[selected], rtol=1e-12),
                "alpha values or selection",
            )

            def integrand(x, omega=omega, spectra=spectra):
                j = min(max(np.searchsorted(omega, x) - 1, 0), 399)
                fraction = (x - omega[j]) / (omega[j + 1] - omega[j])
                density = spectra[:, j] * (1 - fraction) + spectra[:, j + 1] * fraction
                exponent = np.where(x >= 0, -tau * x, (10 - tau) * x)
                return (
                    density[:, None] * (-np.exp(exponent) * expit(10 * abs(x)))[None, :]
                )

            independent, error = quad_vec(
                integrand, -4, 4, points=omega[1:-1], epsabs=1e-10, epsrel=1e-10
            )
            require(
                np.isfinite(independent).all()
                and np.isfinite(error)
                and np.isfinite(producer_prediction).all(),
                "quadrature failure",
            )
            residual = (independent - y) @ whitening.T
            chi = np.sum(residual * residual, axis=1)
            scores = np.full(33, np.inf)
            scores[converged] = np.abs(
                np.log(np.maximum(chi[converged], np.finfo(float).tiny) / keep.sum())
            )
            recomputed = int(np.argmin(scores))
            pred_error = float(
                np.max(np.abs(independent[selected] - producer_prediction[selected]))
            )
            all_pred_error = float(np.max(np.abs(independent - producer_prediction)))
            reported = record["diagnostics"]["selection_chi2"]
            require(
                isinstance(reported, list)
                and len(reported) == 33
                and all(
                    (
                        isinstance(reported[i], (int, float))
                        and not isinstance(reported[i], bool)
                        and np.isfinite(reported[i])
                    )
                    if converged[i]
                    else reported[i] is None
                    for i in range(33)
                ),
                "selection chi-square inventory",
            )
            reported_error = max(
                abs(float(reported[i]) - chi[i])
                for i in range(33)
                if converged[i] and reported[i] is not None
            )
            output["cases"][case] = {
                "selected_index": selected,
                "independent_selected_index": recomputed,
                "selected_training_chi2": float(chi[selected]),
                "selected_whitened_rms": float(np.sqrt(chi[selected] / keep.sum())),
                "max_mass_error": float(np.max(abs(masses - 1))),
                "native_chi2_max_abs_error": float(
                    np.max(abs(native_chi - native_recomputed))
                ),
                "selected_prediction_max_abs_error": pred_error,
                "all_trial_prediction_max_abs_error": all_pred_error,
                "eligible_selection_chi2_max_abs_error": float(reported_error),
                "quad_vec_estimated_global_error": float(error),
                "converged_count": int(converged.sum()),
                "max_latent_zero_count": maximum_zeros,
                "min_observed_H": parameterization["min_observed_H"],
                "min_observed_mixture_floor": parameterization[
                    "min_observed_mixture_floor"
                ],
            }
            validate_selection(record["converged"], selected, recomputed)
            require(
                all_pred_error <= 1e-8 and reported_error <= 1e-3,
                "independent training forward discrepancy",
            )
        signal.alarm(0)
        validate_outputs(directory, receipt.get("output_sha256"))
        require(
            digest(receipt_path) == args.receipt_sha256
            and digest(fixture) == INPUT_SHA,
            "receipt or fixture changed during audit",
        )
        output["status"] = "training_checks_passed_uncertified"
    except Exception as exc:
        output["status"] = "inconclusive" if isinstance(exc, TimeoutError) else "failed"
        output["reason"] = f"{type(exc).__name__}: {exc}"
    finally:
        signal.alarm(0)
        with args.output.open("x") as file:
            file.write(json.dumps(output, indent=2, allow_nan=False) + "\n")
    return output


if __name__ == "__main__":
    result = main()
    print(json.dumps(result, indent=2, allow_nan=False))
    sys.exit(0 if result["status"] == "training_checks_passed_uncertified" else 1)
