#!/usr/bin/env python3
"""Validate projection matrices without continuation or sampling."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
PACKAGE_FILES = {
    "basis": "kernelgecp/basis.py",
    "fitting": "kernelgecp/fitting.py",
    "maxent": "kernelgecp/integrations/maxent.py",
}


def load_helper(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load source helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_import_path(path, root, require_installed):
    path = Path(path).resolve()
    if require_installed and path.is_relative_to(root.resolve()):
        raise RuntimeError("installed wheel required; checkout module imported")
    return path


def sha_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_sha(a):
    import numpy as np

    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def independent_projector(a):
    import mpmath
    import numpy as np

    ctx = mpmath.mp.clone()
    ctx.dps = 80
    m = ctx.matrix([[ctx.mpf(float(x)) for x in row] for row in a])
    u, singular, _ = ctx.svd(m, full_matrices=False)
    rank = sum(x > ctx.mpf("1e-12") * singular[0] for x in singular)
    retained = u[:, :rank]
    p = retained * retained.T
    return np.array(
        [[float(p[i, j]) for j in range(a.shape[0])] for i in range(a.shape[0])]
    ), int(rank)


def matrix_checks(plan, source, dataset_class, tau):
    import numpy as np

    a = plan.source_design
    oracle, rank = independent_projector(a)
    p = plan.matrix
    delta = np.abs(p - oracle)
    ratio = float(np.max(delta / (2e-8 + 1e-8 * np.abs(oracle))))
    tight = float(delta.max())
    if rank != plan.diagnostics.effective_rank or ratio > 1 or tight > 1e-13:
        raise RuntimeError("independent projector accuracy or rank check failed")
    if (
        np.ascontiguousarray(a, dtype="<f8").tobytes()
        != np.ascontiguousarray(source, dtype="<f8").tobytes()
    ):
        raise RuntimeError("repeated provider design differs from decomposed source")
    rng = np.random.default_rng(49211)
    x = rng.standard_normal((len(tau), len(tau))) / len(tau)
    covariance = x @ x.T + np.eye(len(tau))
    values = rng.standard_normal(len(tau))
    data = dataset_class(10.0, tau, values, covariance)
    projected = plan.apply(data)
    kernel = rng.standard_normal((len(tau), 7)) + 1j * rng.standard_normal(
        (len(tau), 7)
    )
    errors = {
        "projected_values": float(np.max(np.abs(projected.values - oracle @ values))),
        "projected_covariance": float(
            np.max(np.abs(projected.covariance - oracle @ covariance @ oracle.T))
        ),
        "complex_kernel": float(np.max(np.abs(p @ kernel - oracle @ kernel))),
        "complex_linearity": float(
            np.max(np.abs(p @ (1j * kernel) - 1j * (p @ kernel)))
        ),
    }
    if max(errors.values()) > 1e-13:
        raise RuntimeError("linear projection propagation check failed")
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    keep = s > 1e-12 * s[0]
    default_error = float(np.max(np.abs(u[:, keep] @ u[:, keep].T - oracle)))
    return {
        "shape": list(a.shape),
        "source_design_sha256_le_f8_c": canonical_sha(a),
        "repeated_design_sha256_le_f8_c": canonical_sha(source),
        "bitwise_equal": True,
        "hash_layout": "little-endian float64, C order",
        "source_original_dtype": a.dtype.str,
        "source_original_strides": list(a.strides),
        "original_rank": plan.diagnostics.effective_rank,
        "original_condition": plan.diagnostics.condition_number,
        "transformed_rank": plan.transformed_diagnostics.effective_rank,
        "transformed_condition": plan.transformed_diagnostics.condition_number,
        "projection_provenance": json.loads(plan.projection_provenance_json),
        "independent_svd_decimal_digits": 80,
        "independent_rank": rank,
        "projector_max_abs": tight,
        "original_budget_max_ratio": ratio,
        "tight_absolute_limit": 1e-13,
        "propagation_max_abs": errors,
        "binary64_svd_projector_max_abs_vs_oracle": default_error,
    }


def run_probe(root, expected_revision, require_installed):
    root = root.resolve()
    receipt = {
        "schema": "kernelgecp.provider-projection-precision.v1",
        "status": "failed",
        "scope": (
            "provider projection matrix validation only; no MaxEnt continuation, "
            "calibration, or scientific acceptance"
        ),
        "checks": {},
        "cases": {},
    }
    stage = "source"
    try:
        for name in THREADS:
            os.environ[name] = "1"
        helpers = root / "scripts/providers"
        triqs = load_helper(helpers / "triqs_smoke.py", "_projection_triqs_helpers")
        wheel = load_helper(helpers / "installed_wheel.py", "_projection_wheel_helper")
        identity = triqs.source_identity(root)
        if identity["source_revision"] != expected_revision:
            raise RuntimeError("source revision mismatch")
        tracked = [
            *PACKAGE_FILES.values(),
            "experiments/dlr_adapters.py",
            "scripts/providers/triqs_smoke.py",
            "scripts/providers/installed_wheel.py",
            "scripts/providers/projection_precision_smoke.py",
        ]
        source_hashes = {}
        for path in tracked:
            blob = subprocess.check_output(
                ["git", "-C", str(root), "show", expected_revision + ":" + path]
            )
            if (root / path).read_bytes() != blob:
                raise RuntimeError("working source differs from revision: " + path)
            source_hashes[path] = hashlib.sha256(blob).hexdigest()
        if (
            Path(__file__).read_bytes()
            != (root / "scripts/providers/projection_precision_smoke.py").read_bytes()
        ):
            raise RuntimeError("executed probe differs from bound source")
        receipt.update(
            {
                "source_revision": expected_revision,
                "identities": identity,
                "source_sha256": source_hashes,
                "environment": triqs.environment_identity(Path(sys.prefix)),
                "versions": {"python": platform.python_version()},
            }
        )
        stage = "imports"
        with triqs.provider_output_to_stderr():
            np = importlib.import_module("numpy")
            scipy = importlib.import_module("scipy")
            mp = importlib.import_module("mpmath")
            pools = importlib.import_module("threadpoolctl")
            modules = {
                key: importlib.import_module(name)
                for key, name in {
                    "basis": "kernelgecp.basis",
                    "fitting": "kernelgecp.fitting",
                    "maxent": "kernelgecp.integrations.maxent",
                }.items()
            }
            module_paths = {
                k: validate_import_path(v.__file__, root, require_installed)
                for k, v in modules.items()
            }
            for name, path in module_paths.items():
                if sha_file(path) != source_hashes[PACKAGE_FILES[name]]:
                    raise RuntimeError("imported module differs from source: " + name)
            if require_installed:
                receipt["installation"] = wheel.installed_wheel_identity(
                    module_paths, PACKAGE_FILES
                )
            # Only experiments is added: the installed kernelgecp remains authoritative.
            sys.path.insert(0, str(root / "experiments"))
            adapters = importlib.import_module("dlr_adapters")
            if (
                Path(adapters.__file__).resolve()
                != (root / "experiments/dlr_adapters.py").resolve()
            ):
                raise RuntimeError("unexpected comparison adapter")
            for name in ("pydlr", "sparse_ir"):
                importlib.import_module(name)
            versions = importlib.metadata
            receipt["versions"].update(
                {
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "mpmath": mp.__version__,
                    "pydlr": versions.version("pydlr"),
                    "sparse_ir": versions.version("sparse-ir"),
                }
            )
            stage = "projection"
            with pools.threadpool_limits(limits=1):
                info = pools.threadpool_info()
                if not info or any(x["num_threads"] != 1 for x in info):
                    raise RuntimeError("thread pool limit unavailable")
                tau = np.linspace(0, 10, 65)
                bases = {
                    "kernelgecp": modules["basis"].FermionicBasis(
                        10, 4, 1e-10, precision="float64", grid_order=24, max_rank=96
                    ),
                    "pydlr": adapters.build_adapter("pydlr", 10, 4, 1e-10),
                    "sparse_ir": adapters.build_adapter("sparse_ir", 10, 4, 1e-10),
                }
                for name, basis in bases.items():
                    plan = modules["maxent"].prepare_spectral_projection(
                        basis, tau, tau, rcond=1e-12, precision_bits=256
                    )
                    repeat = np.asarray(basis.design_tau(tau), dtype=np.float64)
                    receipt["cases"][name] = matrix_checks(
                        plan, repeat, modules["maxent"].SpectralDataset, tau
                    )
        receipt["checks"] = {
            "source_revision_and_hashes": True,
            "installed_wheel": bool(require_installed),
            "all_three_provider_projection_matrices": True,
            "single_thread": True,
        }
        receipt["status"] = "provider_projection_matrix_validation_passed"
    except Exception as exc:
        receipt["error"] = {
            "stage": stage,
            "type": type(exc).__name__,
            "message": str(exc),
        }
        print(
            f"projection smoke failed at {stage}: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--source-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    p.add_argument("--source-revision", required=True)
    p.add_argument("--require-installed", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run_probe(args.source_root, args.source_revision, args.require_installed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return (
        0 if result["status"] == "provider_projection_matrix_validation_passed" else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
