#!/usr/bin/env python3
"""Fail-closed actual CT-HYB measured-feedback wiring, not canonical MC accuracy."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import io
import json
import os
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

from triqs_adapter_smoke import installed_wheel_identity
from triqs_smoke import (
    THREAD_VARIABLES,
    environment_identity,
    provider_output_to_stderr,
    source_identity,
)

SCHEMA = "kernelgecp.cthyb-smoke.v1"
FILES = {
    "basis": "kernelgecp/basis.py",
    "adapter": "kernelgecp/integrations/triqs.py",
    "cthyb": "kernelgecp/integrations/cthyb_dmft.py",
    "fitting": "kernelgecp/fitting.py",
    "core": "kernelgecp/dmft.py",
}


def require(condition: object, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _installed(paths: dict[str, Path]) -> dict[str, Any]:
    identity = installed_wheel_identity(
        {name: paths[name] for name in ("basis", "adapter")}
    )
    distribution = metadata.distribution("kernelgecp")
    record = distribution.read_text("RECORD")
    require(
        isinstance(record, str) and distribution.files is not None,
        "incomplete wheel metadata",
    )
    rows = list(csv.reader(io.StringIO(record or ""), strict=True))
    listed = {str(path) for path in distribution.files or ()}
    for name, relative in FILES.items():
        require(
            relative in listed and sum(row[0] == relative for row in rows) == 1,
            "measurement module absent from installed wheel metadata",
        )
        require(
            Path(distribution.locate_file(relative)).resolve() == paths[name],
            "measurement module import is not owned by the installed wheel",
        )
    return {**identity, "module_count": len(FILES)}


def run(revision: str, *, require_installed: bool) -> dict[str, Any]:
    for name in THREAD_VARIABLES:
        require(os.environ.get(name) == "1", f"{name}=1 is required")
    root = Path(__file__).resolve().parents[2]
    source = source_identity(root)
    require(
        source["source_revision"] == revision,
        "requested source revision differs from actual HEAD",
    )
    source["probe_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    modules = {
        name: importlib.import_module(relative[:-3].replace("/", "."))
        for name, relative in FILES.items()
    }
    paths = {name: Path(module.__file__).resolve() for name, module in modules.items()}
    hashes = {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in paths.items()
    }
    for name, relative in FILES.items():
        require(
            hashes[name] == hashlib.sha256((root / relative).read_bytes()).hexdigest(),
            "imported numerical code differs from the requested source",
        )
    installed = _installed(paths) if require_installed else {"required": False}
    if require_installed:
        require(
            all(not path.is_relative_to(root) for path in paths.values()),
            "installed-wheel acceptance requires imports outside the checkout",
        )
    np = importlib.import_module("numpy")
    cthyb = modules["cthyb"]
    basis = modules["basis"].FermionicBasis(10, 8, 1e-7, grid_order=24, max_rank=40)
    config = cthyb.CTHYBConfig(n_iw=32, n_tau=385, warmup_cycles=500)
    observed = np.flatnonzero(np.arange(config.n_tau) % 16 != 7)
    dense = cthyb.CTHYBPlan(config)
    plans = {"kernelgecp": cthyb.CTHYBPlan(config, basis, fit_tau_indices=observed)}
    # Add only the adapter source directory; adding the checkout root would
    # defeat the outside-checkout installed-wheel test by shadowing its package.
    sys.path.insert(0, str(root / "experiments"))
    adapters = importlib.import_module("dlr_adapters")
    for method in ("pydlr", "sparse_ir"):
        plans[method] = cthyb.CTHYBPlan(
            config,
            adapters.build_adapter(method, 10, 8, 1e-7),
            fit_tau_indices=observed,
        )
    bath = dense.initial_bath()
    measurements = [
        dense.sample(bath, cycles=3000, seed=27190001 + i) for i in range(3)
    ]
    m = measurements[0]
    reps = {name: plan.represent(m) for name, plan in plans.items()}
    selected = np.array([0, 1, 3, 7]) + config.n_iw
    discrepancies = {
        name: float(np.max(np.abs(rep.g_iw[selected] - m.g_iw[selected])))
        for name, rep in reps.items()
    }
    native_bath, native_defect = dense.bath_update(m)
    compressed_bath, compressed_defect = plans["kernelgecp"].bath_update(
        m, reps["kernelgecp"]
    )
    dense_second = dense.sample(native_bath, cycles=3000, seed=27191001)
    compressed_second = plans["kernelgecp"].sample(
        compressed_bath, cycles=3000, seed=27192001
    )
    estimate = cthyb.estimate_replicas(
        np.array([dense.observables(item, [0, 1, 3, 7]) for item in measurements])
    )
    expected = m.bath_iw + config.mixing * (
        reps["kernelgecp"].g_iw.mean(axis=1) - m.bath_iw
    )
    feedback_error = float(np.max(np.abs(compressed_bath - expected)))
    moment_error = float(
        np.max(np.abs(m.sigma_moments[0] - config.interaction * m.occupancies[::-1]))
    )
    checks = {
        "imports": True,
        "numerical": True,
        "source_bound": True,
        "actual_fresh_native_chains": all(
            item.solve_status == 0
            for item in (*measurements, dense_second, compressed_second)
        ),
        "density_matrix_moments": moment_error < 1e-12,
        "measured_time_feedback": feedback_error < 1e-12 and compressed_defect > 0,
        "native_feedback_executed": native_defect > 0
        and not np.array_equal(native_bath, bath),
        "all_three_identical_measurement_replays": all(
            rep.rank > 0 and discrepancies[name] < 0.1 for name, rep in reps.items()
        ),
        "independent_replica_covariance": bool(
            np.all(np.isfinite(estimate.covariance))
            and np.max(estimate.standard_error) > 0
        ),
        "owns_outputs": all(
            not arr.flags.writeable
            for arr in (m.g_tau, m.g_iw, compressed_bath, reps["kernelgecp"].g_iw)
        ),
    }
    if require_installed:
        checks["installed_wheel"] = True
    require(all(checks.values()), f"CT-HYB wiring checks failed: {checks}")
    versions = {
        "python": sys.version.split()[0],
        "triqs": importlib.import_module("triqs.version").version,
        "cthyb": importlib.import_module("triqs_cthyb.version").version,
        **{
            name: metadata.version(name)
            for name in ("numpy", "scipy", "pydlr", "sparse-ir", "pylibsparseir")
        },
    }
    require(
        versions["triqs"] == "4.0.2" and versions["cthyb"] == "4.0.1",
        "unsupported actual providers",
    )
    return {
        "schema": SCHEMA,
        "provider": "cthyb",
        "status": "passed",
        **source,
        "versions": versions,
        "checks": checks,
        "environment": environment_identity(Path(sys.prefix)),
        "package_sources": hashes,
        "installed_identity": installed,
        "metrics": {
            "feedback_algebra_error": feedback_error,
            "native_hartree_moment_error": moment_error,
            "selected_replay_discrepancies": discrepancies,
            "density": m.density,
            "double_occupancy": m.double_occupancy,
            "replica_standard_errors": estimate.standard_error.tolist(),
        },
        "sampling": {
            "cycles": 3000,
            "warmup_cycles": 500,
            "length_cycle": 50,
            "replica_seeds": [27190001, 27190002, 27190003],
            "dense_second_seed": 27191001,
            "gecp_second_seed": 27192001,
        },
        "acceptance_scope": (
            "fast_actual_wiring_not_canonical_uncertainty_or_physical_DMFT_acceptance"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--require-installed", action="store_true")
    args = parser.parse_args()
    try:
        with provider_output_to_stderr():
            result = run(args.source_revision, require_installed=args.require_installed)
    except Exception as error:
        result = {
            "schema": SCHEMA,
            "provider": "cthyb",
            "status": "failed",
            "error_type": type(error).__name__,
            "error": str(error),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )
        print(json.dumps(result, sort_keys=True), file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
