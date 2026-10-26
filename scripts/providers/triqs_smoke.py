#!/usr/bin/env python3
"""Execute the pinned TRIQS/CT-HYB/MaxEnt runtime, failing closed.

Run inside environments/triqs.yml. Provider output goes to stderr; stdout is
one JSON receipt. This is a runtime calculation probe, not a DMFT or spectral
accuracy benchmark. No optional provider is replaced by a mock or fallback.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib
import json
import os
import platform
import re
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

SCHEMA = "kernelgecp.provider-smoke.v1"
EXPECTED_VERSIONS = {"triqs": "4.0.2", "cthyb": "4.0.1", "maxent": "4.0.1"}
MAXENT_WHEEL_SHA256 = "74c0fb1a4ec4abbc0ee8b46d11f9f01781687d4bfbe74982c0100c0887b4fcac"
THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def require(checks: dict[str, bool], name: str, condition: object) -> None:
    """Record the actual check before rejecting it; never infer success."""
    checks[name] = bool(condition)
    if not checks[name]:
        raise RuntimeError(f"provider check failed: {name}")


def validate_versions(versions: dict[str, str]) -> None:
    for name, expected in EXPECTED_VERSIONS.items():
        if versions.get(name) != expected:
            raise RuntimeError(f"unexpected {name} version; expected {expected}")
    if not re.fullmatch(r"3\.12\.\d+", versions.get("python", "")):
        raise RuntimeError("the supported provider runtime requires Python 3.12")


def module_version(module_name: str) -> str:
    module = importlib.import_module(module_name + ".version")
    value = getattr(module, "version", None)
    if not isinstance(value, str):
        raise RuntimeError("provider does not expose its release version")
    return value


def environment_identity(prefix: Path) -> dict[str, Any]:
    """Hash package identity without exporting prefixes, URLs or user paths."""
    packages = []
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        record = json.loads(path.read_text())
        packages.append(
            {key: record.get(key, "") for key in ("name", "version", "build", "subdir")}
        )
    if not packages:
        raise RuntimeError("a resolved conda package identity is required")
    payload = json.dumps(packages, sort_keys=True, separators=(",", ":")).encode()
    return {
        "system": platform.system(),
        "machine": platform.machine(),
        "python_implementation": platform.python_implementation(),
        "conda_package_count": len(packages),
        "conda_package_set_sha256": hashlib.sha256(payload).hexdigest(),
        "thread_limits": {name: os.environ.get(name) for name in THREAD_VARIABLES},
    }


def source_identity(root: Path) -> dict[str, str]:
    revision = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise RuntimeError("a full source revision is required")
    return {
        "source_revision": revision,
        "probe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "environment_spec_sha256": hashlib.sha256(
            (root / "environments" / "triqs.yml").read_bytes()
        ).hexdigest(),
        "maxent_wheel_sha256": MAXENT_WHEEL_SHA256,
    }


@contextlib.contextmanager
def provider_output_to_stderr() -> Iterator[None]:
    """Also redirect native C++ writes, which bypass redirect_stdout alone."""
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        os.dup2(2, 1)
        with contextlib.redirect_stdout(sys.stderr):
            yield
    finally:
        sys.stderr.flush()
        os.dup2(saved, 1)
        os.close(saved)


def run_probe(expected_revision: str | None = None) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    checks = {"imports": False, "numerical": False}
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "provider": "triqs",
        "status": "failed",
        "checks": checks,
        "versions": {"python": platform.python_version()},
        "metrics": {},
    }
    stage = "environment"
    try:
        for name in THREAD_VARIABLES:
            os.environ[name] = "1"
        receipt["identities"] = source_identity(root)
        receipt["source_revision"] = receipt["identities"]["source_revision"]
        stage = "source_revision"
        require(
            checks,
            "source_revision",
            expected_revision is None
            or (
                re.fullmatch(r"[0-9a-f]{40}", expected_revision) is not None
                and receipt["identities"]["source_revision"] == expected_revision
            ),
        )
        stage = "environment"
        receipt["environment"] = environment_identity(Path(sys.prefix))
        with provider_output_to_stderr():
            stage = "imports"
            np = importlib.import_module("numpy")
            scipy = importlib.import_module("scipy")
            pools = importlib.import_module("threadpoolctl")
            gfs = importlib.import_module("triqs.gfs")
            mesh_types = importlib.import_module("triqs.mesh")
            operators = importlib.import_module("triqs.operators")
            mpi = importlib.import_module("triqs.utility.mpi")
            cthyb = importlib.import_module("triqs_cthyb")
            maxent = importlib.import_module("triqs_maxent")
            checks["imports"] = True
            receipt["versions"].update(
                {
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "triqs": module_version("triqs"),
                    "cthyb": module_version("triqs_cthyb"),
                    "maxent": module_version("triqs_maxent"),
                }
            )
            stage = "versions"
            validate_versions(receipt["versions"])
            checks["pinned_versions"] = True
            require(checks, "one_mpi_rank", mpi.size == 1)
            with pools.threadpool_limits(limits=1):
                stage = "triqs_resolvent"
                beta, energy = 4.0, 0.3
                green = gfs.Gf(
                    mesh=mesh_types.MeshImFreq(beta=beta, statistic="Fermion", n_iw=32),
                    target_shape=[1, 1],
                )
                green << gfs.inverse(gfs.iOmega_n - energy)
                iw = np.array([complex(point.value) for point in green.mesh])
                expected_iw = 1j * (2 * np.arange(-32, 32) + 1) * np.pi / beta
                require(
                    checks,
                    "fermionic_frequency_labels",
                    iw.shape == (64,)
                    and float(np.max(np.abs(iw - expected_iw))) < 1e-12,
                )
                expected = 1.0 / (expected_iw - energy)
                error = float(np.max(np.abs(green.data[:, 0, 0] - expected)))
                require(checks, "triqs_resolvent", np.isfinite(error) and error < 1e-12)
                receipt["metrics"]["triqs_resolvent_max_error"] = error

                stage = "cthyb_sampling"
                solver = cthyb.Solver(
                    beta=beta,
                    gf_struct=[("up", 1), ("down", 1)],
                    n_iw=32,
                    n_tau=257,
                )
                for _, bath in solver.G0_iw:
                    bath << gfs.inverse(gfs.iOmega_n + 0.5 - 0.3**2 * gfs.Wilson(1))
                solver.solve(
                    h_int=operators.n("up", 0) * operators.n("down", 0),
                    n_cycles=200,
                    n_warmup_cycles=50,
                    length_cycle=16,
                    random_seed=20261026,
                    verbosity=0,
                )
                tau_values = [np.asarray(block.data) for _, block in solver.G_tau]
                iw_values = [np.asarray(block.data) for _, block in solver.G_iw]
                sign = complex(solver.average_sign)
                require(
                    checks,
                    "cthyb_sampling",
                    len(tau_values) == 2
                    and len(iw_values) == 2
                    and all(
                        np.all(np.isfinite(array)) for array in tau_values + iw_values
                    )
                    and any(np.max(np.abs(array)) > 0 for array in tau_values)
                    and np.isfinite(sign)
                    and 0 < abs(sign) <= 1 + 1e-12,
                )
                receipt["metrics"].update(
                    {"cthyb_cycles": 200, "cthyb_average_sign_real": float(sign.real)}
                )

                stage = "maxent_minimization"
                tau = np.linspace(0, beta, 33)
                exact = -0.5 * (
                    np.exp(-tau * 0.5) / (1 + np.exp(-beta * 0.5))
                    + np.exp(tau * 0.5) / (1 + np.exp(beta * 0.5))
                )
                rng = np.random.default_rng(20261026)
                observed = exact + rng.normal(0, 1e-4, tau.size)
                continuation = maxent.TauMaxEnt(cost_function="bryan")
                continuation.omega = maxent.LinearOmegaMesh(-2, 2, 81)
                continuation.alpha_mesh = maxent.LogAlphaMesh(1, 10, 3)
                continuation.analyzers = []
                continuation.interactive = False
                continuation.set_G_tau_data(tau, observed)
                continuation.set_error(1e-3)
                result = continuation.run()
                spectra = np.asarray(result.A)
                chi2 = np.asarray(result.chi2)
                require(
                    checks,
                    "maxent_minimization",
                    spectra.shape == (3, 81)
                    and chi2.shape == (3,)
                    and np.all(np.isfinite(spectra))
                    and np.all(spectra >= 0)
                    and np.any(spectra > 0)
                    and np.all(np.isfinite(chi2))
                    and np.all(chi2 >= 0),
                )
                receipt["metrics"].update(
                    {"maxent_alpha_count": 3, "maxent_min_chi2": float(np.min(chi2))}
                )
                stage = "thread_limits"
                threadpools = pools.threadpool_info()
                require(
                    checks,
                    "single_thread",
                    bool(threadpools)
                    and all(pool["num_threads"] == 1 for pool in threadpools),
                )
                receipt["environment"]["observed_threadpools"] = [
                    {key: pool.get(key) for key in ("internal_api", "num_threads")}
                    for pool in threadpools
                ]
        checks["numerical"] = all(checks[key] for key in checks if key != "numerical")
        require(checks, "numerical", checks["numerical"])
        receipt["status"] = "passed"
    except Exception as error:
        receipt["error"] = {"stage": stage, "type": type(error).__name__}
        print(
            f"TRIQS provider probe failed at {stage}: {type(error).__name__}",
            file=sys.stderr,
        )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source-revision")
    args = parser.parse_args()
    receipt = run_probe(args.source_revision)
    rendered = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    sys.stdout.write(rendered)
    return 0 if receipt["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
