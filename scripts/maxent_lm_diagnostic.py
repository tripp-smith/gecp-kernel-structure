#!/usr/bin/env python3
"""Three fixed native LM probes; unit fixture or preexisting calibration archive."""

import argparse
import contextlib
import dataclasses
import hashlib
import importlib
import inspect
import io
import json
import math
import os
import platform
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import scipy
import triqs_maxent
import triqs_maxent.minimizers as minimizers

from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    continue_tau_maxent,
    load_spectral_dataset,
    prepare_spectral_projection,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        return clean(x.item())
    if isinstance(x, float) and not math.isfinite(x):
        return "NaN" if math.isnan(x) else ("+Infinity" if x > 0 else "-Infinity")
    return x


def write_report(path, obj):
    path.write_text(json.dumps(clean(obj), indent=2, allow_nan=False))


parser = argparse.ArgumentParser()
parser.add_argument(
    "--dataset", type=Path, help="existing calibration spectral archive; no generation"
)
parser.add_argument(
    "--model", choices=["broad_gaussian"], help="required for calibration archive"
)
parser.add_argument(
    "--output", type=Path, default=Path("/tmp/gecp-maxent-lm-linux.json")
)
args = parser.parse_args()
if bool(args.dataset) != bool(args.model):
    parser.error("--dataset and --model must be provided together")
if args.dataset:
    projected = load_spectral_dataset(args.dataset)
    config = MaxEntConfig(n_omega=401, n_alpha=33, default_model=args.model)
    manifest = json.loads(args.dataset.with_suffix(".json").read_text())
    if manifest.get("calculation_settings") != dataclasses.asdict(config):
        raise ValueError(
            "archive settings differ from declared calibration configuration"
        )
else:

    class PolynomialDesign:
        beta = 3.0

        def design_tau(self, tau):
            return np.polynomial.legendre.legvander(
                2 * np.asarray(tau) / self.beta - 1, 3
            )

    tau = np.linspace(0, 3, 9)
    values = -np.exp(-0.7 * tau) / (1 + np.exp(-3 * 0.7))
    data = SpectralDataset(3, tau, values, 1e-6 * np.eye(len(tau)))
    projected = prepare_spectral_projection(
        PolynomialDesign(), tau, np.linspace(0, 3, 17)
    ).apply(data)
    config = MaxEntConfig(omega_min=-2, omega_max=2, n_omega=81, n_alpha=9)
presets = [("baseline", {}), ("mu0_1e-12", {"mu0": 1e-12}), ("mu0_1e-8", {"mu0": 1e-8})]
original_class, original_solve = minimizers.LevenbergMinimizer, np.linalg.solve
source = inspect.getsourcefile(original_class)
package_dir = Path(triqs_maxent.__file__).resolve().parent
provider_sources = {
    str(p.relative_to(package_dir)): sha(p) for p in sorted(package_dir.rglob("*.py"))
}
blas = io.StringIO()
with contextlib.redirect_stdout(blas):
    np.show_config()
threads = [
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
]
report = {
    "repo_head": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip(),
    "script_sha256": sha(__file__),
    "integration_sha256": sha("kernelgecp/integrations/maxent.py"),
    "provider_source_hashes": provider_sources,
    "minimizer_source": source,
    "python": sys.version,
    "platform": platform.platform(),
    "processor": platform.processor(),
    "cpu_count": os.cpu_count(),
    "numpy": np.__version__,
    "scipy": scipy.__version__,
    "provider": importlib.import_module("triqs_maxent.version").version,
    "numpy_config": blas.getvalue(),
    "threads": {x: os.environ.get(x) for x in threads},
    "dataset_path": str(args.dataset) if args.dataset else None,
    "dataset_sha256": sha(args.dataset) if args.dataset else None,
    "dataset_manifest_sha256": sha(args.dataset.with_suffix(".json"))
    if args.dataset
    else None,
    "fixture_fingerprint": projected.fingerprint,
    "config": dataclasses.asdict(config),
    "schema": "kernelgecp.maxent-lm-diagnostic.v1",
    "scope": "Consumed calibration fixtures only; not production acceptance.",
    "trials": [],
}
npzpath = args.output.with_suffix(".failure.npz")


def probe_type(opts, hist):
    class ProbeMinimizer(original_class):
        def __init__(self, *a, **kw):
            super().__init__(*a, **(kw | opts))

        def minimize(self, function, v0):
            try:
                value = super().minimize(function, v0)
            except Exception:
                hist.append({"status": "exception"})
                raise
            hist.append(
                {"converged": bool(self.converged), "iterations": int(self.n_iter_last)}
            )
            return value

    return ProbeMinimizer


def solve_capture(failure, npzpath):
    def capture_solve(a, b, *pos, **kw):
        try:
            return original_solve(a, b, *pos, **kw)
        except np.linalg.LinAlgError:
            try:
                frame = inspect.currentframe().f_back
                loc = frame.f_locals
                failure.update(
                    {
                        "caller_file": frame.f_code.co_filename,
                        "caller_line": frame.f_lineno,
                        "caller_function": frame.f_code.co_name,
                        "i": loc.get("i"),
                        "mu": loc.get("mu"),
                    }
                )
                # Copy directly stored alpha; do not evaluate a native function
                # or invoke a property to observe it.
                alpha = loc.get("alpha")
                if isinstance(alpha, (int, float, np.integer, np.floating)):
                    failure["scalar_alpha"] = float(alpha)
                else:
                    function = loc.get("function")
                    stored = getattr(function, "__dict__", {})
                    alpha = stored.get("_alpha") if isinstance(stored, dict) else None
                    failure["scalar_alpha"] = (
                        float(alpha)
                        if isinstance(alpha, (int, float, np.integer, np.floating))
                        else None
                    )
                arrays = {"matrix": np.asarray(a), "rhs": np.asarray(b)}
                for key in ("mu", "J", "f", "Id"):
                    val = loc.get(key)
                    if val is not None:
                        arrays[key] = np.asarray(val)
                if all(
                    x.ndim <= 2
                    and all(d <= 2048 for d in x.shape)
                    and x.size <= 2048**2
                    and x.dtype.kind in "biufc"
                    for x in arrays.values()
                ):
                    np.savez_compressed(npzpath, **arrays)
                    failure["npz"] = str(npzpath)
                    failure["npz_sha256"] = sha(npzpath)
                    failure["matrix_shape"] = list(np.asarray(a).shape)
                    failure["matrix_all_finite"] = bool(np.isfinite(a).all())
                    failure["rhs_all_finite"] = bool(np.isfinite(b).all())
                else:
                    failure["npz"] = None
            except Exception as capture_error:
                failure["capture_error"] = repr(capture_error)
            raise

    return capture_solve


for name, opts in presets:
    hist, failure = [], {}
    row = {"name": name, "overrides": opts, "effective_mu0": opts.get("mu0", 1e-18)}
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    try:
        minimizers.LevenbergMinimizer = probe_type(opts, hist)
        if name == "baseline":
            np.linalg.solve = solve_capture(failure, npzpath)
        result = continue_tau_maxent(projected, config)
        training_rms = result.normalized_training_whitened_rms
        row.update(
            {
                "status": "completed",
                "selected_index": result.selected_index,
                "alphas": result.alphas.tolist(),
                "chi2": result.chi2.tolist(),
                "normalized_physical_chi2": result.diagnostics["selection_chi2"],
                "normalized_training_whitened_rms": training_rms,
                "raw_mass": result.raw_mass,
                "mass": result.mass,
                "converged": result.converged,
                "iterations": result.iterations,
                "spectra": result.trial_spectra.tolist(),
            }
        )
    except Exception as e:
        row.update(
            {
                "status": "exception",
                "exception_chain": "".join(traceback.format_exception(e)),
            }
        )
    finally:
        np.linalg.solve, minimizers.LevenbergMinimizer = original_solve, original_class
        row.update(
            {
                "wall_seconds": time.perf_counter() - start_wall,
                "cpu_seconds": time.process_time() - start_cpu,
                "history": hist,
                "failure": failure,
            }
        )
        report["trials"].append(row)
        write_report(args.output, report)
base = report["trials"][0]
if base["status"] == "completed":
    for row in report["trials"][1:]:
        if row["status"] == "completed":
            row["vs_baseline"] = {
                "max_abs_spectrum": float(
                    np.max(
                        np.abs(np.asarray(row["spectra"]) - np.asarray(base["spectra"]))
                    )
                ),
                "max_abs_native_chi2": float(
                    np.max(np.abs(np.asarray(row["chi2"]) - np.asarray(base["chi2"])))
                ),
                "max_abs_physical_chi2_common_eligible": max(
                    (
                        abs(a - b)
                        for a, b in zip(
                            row["normalized_physical_chi2"],
                            base["normalized_physical_chi2"],
                            strict=True,
                        )
                        if a is not None and b is not None
                    ),
                    default=None,
                ),
            }
runtime_info = io.StringIO()
with contextlib.redirect_stdout(runtime_info):
    try:
        np.show_runtime()
    except Exception as error:
        print(type(error).__name__, str(error))
report["numpy_runtime"] = runtime_info.getvalue()
write_report(args.output, report)
for row in report["trials"]:
    print(
        row["name"],
        row["status"],
        len(row["history"]),
        row.get("normalized_training_whitened_rms"),
    )
