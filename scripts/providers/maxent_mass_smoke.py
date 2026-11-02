#!/usr/bin/env python3
"""Probe the installed native MaxEnt mass interface, not scientific acceptance."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import importlib
import importlib.metadata
import importlib.util
import io
import json
import os
import subprocess
import sys
from pathlib import Path

FILES = {
    "maxent": "kernelgecp/integrations/maxent.py",
    "mass": "kernelgecp/integrations/_maxent_mass.py",
    "native_mass": "kernelgecp/integrations/_maxent_native_mass.py",
}
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


def verify_record_hashes(paths):
    distribution = importlib.metadata.distribution("kernelgecp")
    record = distribution.read_text("RECORD")
    if not record:
        raise RuntimeError("missing wheel RECORD")
    rows = list(csv.reader(io.StringIO(record), strict=True))
    for key, relative in FILES.items():
        matches = [row for row in rows if len(row) == 3 and row[0] == relative]
        if len(matches) != 1:
            raise RuntimeError("missing or duplicate wheel RECORD entry: " + relative)
        digest = (
            base64.urlsafe_b64encode(hashlib.sha256(paths[key].read_bytes()).digest())
            .rstrip(b"=")
            .decode("ascii")
        )
        if matches[0][1] != "sha256=" + digest or matches[0][2] != str(
            paths[key].stat().st_size
        ):
            raise RuntimeError("wheel RECORD hash or size mismatch: " + relative)


def fixture(np):
    """Independent high-order quadrature of a smooth normalized density."""
    beta = 8.0
    tau = np.linspace(0.0, beta, 17)
    nodes, weights = np.polynomial.legendre.leggauss(256)
    omega = 4.0 * nodes
    weights = 4.0 * weights
    density = np.exp(-0.5 * ((omega - 0.4) / 1.0) ** 2)
    density /= weights @ density
    kernel = -np.exp(-tau[:, None] * omega) / (1.0 + np.exp(-beta * omega))
    return beta, tau, kernel @ (weights * density)


def native_cost_checks(np, native):
    """Exercise provider cost chain rules and state refresh on a tiny kernel."""
    provider = importlib.import_module("triqs_maxent")
    kernels = importlib.import_module("triqs_maxent.kernels")
    omega = provider.LinearOmegaMesh(-2.0, 2.0, 9)
    grid = np.linspace(-1.0, 1.0, 9)
    basis, _ = np.linalg.qr(np.column_stack((grid, grid**2, grid**3)))
    matrix = np.diag([3.0, 1.0, 0.1]) @ basis.T
    kernel = kernels.DataKernel(np.arange(3, dtype=float), omega, matrix)
    cost, hidden = native.make_native_mass_cost(1.0)
    cost.K = kernel
    default = provider.DataDefaultModel(np.ones(9), omega)
    cost.D = default
    cost.G = np.array([0.2, -0.1, 0.3])
    cost.err = np.array([0.2, 0.3, 0.4])
    cost.set_alpha(0.7)
    before = hidden.input_size[0]
    kernel.reduce_singular_space(0.5)
    initial = default.D / np.sum(default.D)
    initial_copy = initial.copy()
    initial.setflags(write=False)
    v = hidden(initial).inv()
    after = hidden.input_size[0]
    if not (after < before and after > 0 and v.shape == (after,)):
        raise RuntimeError("singular-space map refresh failed")
    if not np.array_equal(initial, initial_copy):
        raise RuntimeError("inverse mutated its read-only input")
    v = np.linspace(0.1, 0.2, after)
    v_copy = v.copy()
    v.setflags(write=False)
    h = hidden(v).f()
    if abs(float(np.sum(h)) - 1.0) > 1e-12 or np.any(h <= 0):
        raise RuntimeError("native H mass or positivity failed")
    analytic_gradient = cost(v).d()
    analytic_hessian = cost(v).dd()
    step = 1e-5
    identity = np.eye(after)
    gradient = np.array(
        [
            (cost(v + step * e).f() - cost(v - step * e).f()) / (2 * step)
            for e in identity
        ]
    )
    hessian = np.column_stack(
        [
            (cost(v + step * e).d() - cost(v - step * e).d()) / (2 * step)
            for e in identity
        ]
    )
    grad_error = float(
        np.max(
            np.abs(gradient - analytic_gradient)
            / np.maximum(1.0, np.abs(analytic_gradient))
        )
    )
    hess_error = float(
        np.max(
            np.abs(hessian - analytic_hessian)
            / np.maximum(1.0, np.abs(analytic_hessian))
        )
    )
    if (
        not np.isfinite(grad_error)
        or not np.isfinite(hess_error)
        or grad_error > 2e-5
        or hess_error > 5e-4
        or not np.array_equal(v, v_copy)
    ):
        raise RuntimeError("native cost derivative or input-isolation check failed")
    changed = provider.DataDefaultModel(np.exp(grid), omega)
    cost.D = changed
    new_h = hidden(np.zeros(hidden.input_size)).f()
    expected = changed.D / np.sum(changed.D)
    if not np.allclose(new_h, expected, rtol=1e-12, atol=1e-14):
        raise RuntimeError("default-model map refresh failed")
    return {
        "rank_before_reduction": before,
        "rank_after_reduction": after,
        "gradient_max_scaled_error": grad_error,
        "hessian_max_scaled_error": hess_error,
        "gradient_limit": 2e-5,
        "hessian_limit": 5e-4,
        "finite_difference_step": step,
        "default_refresh": True,
        "readonly_input_isolation": True,
    }


def run_probe(root, expected_revision, require_installed):
    result = {
        "schema": "kernelgecp.native-mass-interface-probe.v1",
        "status": "failed",
        "scope": (
            "native interface probe only; not production calibration or acceptance"
        ),
        "cases": {},
    }
    stage = "source"
    try:
        root = root.resolve()
        for name in THREADS:
            os.environ[name] = "1"
        helpers = root / "scripts/providers"
        triqs = load_helper(helpers / "triqs_smoke.py", "_mass_probe_triqs")
        wheel = load_helper(helpers / "installed_wheel.py", "_mass_probe_wheel")
        identity = triqs.source_identity(root)
        if identity["source_revision"] != expected_revision:
            raise RuntimeError("source revision mismatch")
        tracked = [
            *FILES.values(),
            "scripts/providers/triqs_smoke.py",
            "scripts/providers/installed_wheel.py",
            "scripts/providers/maxent_mass_smoke.py",
        ]
        hashes = {}
        for path in tracked:
            blob = subprocess.check_output(
                ["git", "-C", str(root), "show", expected_revision + ":" + path]
            )
            if (root / path).read_bytes() != blob:
                raise RuntimeError("working source differs from revision: " + path)
            hashes[path] = hashlib.sha256(blob).hexdigest()
        if (
            Path(__file__).read_bytes()
            != (helpers / "maxent_mass_smoke.py").read_bytes()
        ):
            raise RuntimeError("executed probe differs from bound source")
        result.update(
            source_revision=expected_revision,
            identities=identity,
            source_sha256=hashes,
            environment=triqs.environment_identity(Path(sys.prefix)),
        )
        stage = "imports"
        with triqs.provider_output_to_stderr():
            np = importlib.import_module("numpy")
            pools = importlib.import_module("threadpoolctl")
            importlib.import_module("triqs_maxent")
            version = importlib.import_module("triqs_maxent.version").version
            if str(version) != "4.0.1":
                raise RuntimeError("actual MaxEnt 4.0.1 required")
            modules = {
                key: importlib.import_module(path[:-3].replace("/", "."))
                for key, path in FILES.items()
            }
            paths = {
                key: Path(module.__file__).resolve() for key, module in modules.items()
            }
            for key, path in paths.items():
                if require_installed and path.is_relative_to(root):
                    raise RuntimeError(
                        "installed wheel required; checkout module imported"
                    )
                if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[FILES[key]]:
                    raise RuntimeError("imported module differs from source: " + key)
            if require_installed:
                result["installation"] = wheel.installed_wheel_identity(paths, FILES)
                verify_record_hashes(paths)
            result["provider_version"] = str(version)
            stage = "native_continuation"
            with pools.threadpool_limits(limits=1):
                info = pools.threadpool_info()
                if not info or any(item["num_threads"] != 1 for item in info):
                    raise RuntimeError("single-thread limit unavailable")
                result["native_cost_checks"] = native_cost_checks(
                    np, modules["native_mass"]
                )
                beta, tau, values = fixture(np)
                sigma = 0.003
                api = modules["maxent"]
                data = api.SpectralDataset(
                    beta, tau, values, np.eye(len(tau)) * sigma**2
                )
                for model in ("flat", "broad_gaussian"):
                    config = api.MaxEntConfig(
                        n_omega=101,
                        n_alpha=5,
                        maxiter=1000,
                        minimizer_mu0=1e-8,
                        default_model=model,
                        normalize_mass=1.0,
                        mass_mode="constrained",
                    )
                    spectrum = api.continue_tau_maxent(data, config)
                    trials = spectrum.trial_spectra
                    masses = np.trapezoid(trials, spectrum.omega, axis=1)
                    diagnostics = spectrum.diagnostics
                    param = diagnostics["mass_parameterization"]
                    if (
                        trials.shape != (5, 101)
                        or not np.isfinite(trials).all()
                        or np.any(trials < 0)
                        or np.max(np.abs(masses - 1)) > 1e-10
                        or not spectrum.converged[spectrum.selected_index]
                        or diagnostics["mass_operation"] != "native_mass_constraint"
                        or diagnostics["config"]["maxiter"] != 1000
                        or diagnostics["config"]["minimizer_mu0"] != 1e-8
                        or not np.isfinite(spectrum.normalized_training_whitened_rms)
                        or not param
                        or not 0 < param["effective_rank"] <= param["original_rank"]
                        or not np.isfinite(param["retained_condition"])
                        or not np.isfinite(spectrum.chi2).all()
                    ):
                        raise RuntimeError(
                            "native constrained trial validation failed: " + model
                        )
                    result["cases"][model] = {
                        "trial_masses": masses.tolist(),
                        "converged": list(spectrum.converged),
                        "iterations": list(spectrum.iterations),
                        "selected_index": spectrum.selected_index,
                        "training_whitened_rms": (
                            spectrum.normalized_training_whitened_rms
                        ),
                        "chi2": spectrum.chi2.tolist(),
                        "mass_parameterization": param,
                        "config": diagnostics["config"],
                        "mass_operation": diagnostics["mass_operation"],
                    }
        result["status"] = "native_mass_interface_probe_passed"
    except Exception as exc:
        result["error"] = {
            "stage": stage,
            "type": type(exc).__name__,
            "message": str(exc),
        }
        print(
            f"mass probe failed at {stage}: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--require-installed", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run_probe(args.source_root, args.source_revision, args.require_installed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return 0 if result["status"] == "native_mass_interface_probe_passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
