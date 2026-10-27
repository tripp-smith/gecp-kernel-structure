#!/usr/bin/env python3
"""Fail-closed actual dense-TRIQS IPT and independent finite-model checks.

This validates the native reference, not compressed-provider benchmark results.
IPT remains an approximate physical model. A missing runtime is never skipped.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import io
import json
import os
import platform
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

from triqs_adapter_smoke import installed_wheel_identity
from triqs_smoke import (
    THREAD_VARIABLES,
    environment_identity,
    provider_output_to_stderr,
    require,
    source_identity,
)

SCHEMA = "kernelgecp.dmft-reference-smoke.v1"
PACKAGE_FILES = {
    "basis": "kernelgecp/basis.py",
    "adapter": "kernelgecp/integrations/triqs.py",
    "core": "kernelgecp/dmft.py",
    "reference": "kernelgecp/integrations/triqs_dmft.py",
    "lehmann": "kernelgecp/lehmann.py",
}


def _wheel_identity(module_paths: dict[str, Path]) -> dict[str, Any]:
    identity = installed_wheel_identity(
        {name: module_paths[name] for name in ("basis", "adapter")}
    )
    distribution = metadata.distribution("kernelgecp")
    files, record = distribution.files, distribution.read_text("RECORD")
    if files is None or not isinstance(record, str):
        raise RuntimeError("installed wheel metadata is incomplete")
    rows = list(csv.reader(io.StringIO(record), strict=True))
    listed = {str(path) for path in files}
    for name, relative in PACKAGE_FILES.items():
        if relative not in listed or sum(row[0] == relative for row in rows) != 1:
            raise RuntimeError("native-reference module is absent from wheel metadata")
        if Path(distribution.locate_file(relative)).resolve() != module_paths[name]:
            raise RuntimeError("native-reference import is not owned by the wheel")
    return {**identity, "module_count": len(PACKAGE_FILES)}


def _pole_oracle(
    np: Any, gfs: Any, meshes: Any, model: Any, interaction: float
) -> dict[str, float]:
    lehmann = importlib.import_module("kernelgecp.lehmann")
    oracle = lehmann.ipt_self_energy_poles(model, interaction)
    labels = np.arange(-512, 512)
    frequency = meshes.MeshImFreq(beta=model.beta, statistic="Fermion", n_iw=512)
    time_mesh = meshes.MeshImTime(beta=model.beta, statistic="Fermion", n_tau=6145)
    giw = gfs.Gf(mesh=frequency.copy(), target_shape=[])
    siw = gfs.Gf(mesh=frequency.copy(), target_shape=[])
    gtau = gfs.Gf(mesh=time_mesh.copy(), target_shape=[])
    stau = gfs.Gf(mesh=time_mesh.copy(), target_shape=[])
    giw.data[:] = model.evaluate_matsubara(labels)
    gtau.set_from_fourier(giw, known_moments=np.r_[0, model.moments(2)].astype(complex))
    tau = np.array([point.value for point in time_mesh])
    stau.data[:] = interaction**2 * gtau.data**2 * gtau.data[::-1]
    # The native tail constraints use only the input Weiss spectral data and
    # endpoint derivatives. The oracle's triple-sum moments are not fed back.
    special = importlib.import_module("scipy.special")
    occupied = special.expit(-model.beta * model.energies)
    empty = special.expit(model.beta * model.energies)
    density = model.density()
    left, right = -(1 - density), -density
    derivative0 = np.sum(model.weights * model.energies * empty)
    derivativeb = np.sum(model.weights * model.energies * occupied)
    first = interaction**2 * density * (1 - density)
    second = interaction**2 * (
        2 * left * derivative0 * right
        - left**2 * derivativeb
        + 2 * right * derivativeb * left
        - right**2 * derivative0
    )
    siw.set_from_fourier(stau, known_moments=np.array([0, first, second], complex))
    iw = 1j * (2 * labels + 1) * np.pi / model.beta
    atom = lehmann.hubbard_atom(model.beta, interaction)
    return {
        "tau_error": float(np.max(np.abs(stau.data - oracle.evaluate_tau(tau)))),
        "iw_error": float(np.max(np.abs(siw.data - oracle.evaluate_matsubara(labels)))),
        "wrong_unreflected_cube_error": float(
            np.max(
                np.abs(
                    interaction**2 * model.evaluate_tau(tau) ** 3
                    - oracle.evaluate_tau(tau)
                )
            )
        ),
        "atomic_dyson_error": float(
            np.max(np.abs(1 / (iw - siw.data) - atom.evaluate_matsubara(labels)))
        ),
    }


def run_probe(
    expected_revision: str, *, require_installed: bool = False
) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    checks = {"imports": False, "numerical": False}
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "provider": "triqs-dense-ipt",
        "status": "failed",
        "versions": {"python": platform.python_version()},
        "checks": checks,
        "metrics": {},
        "runs": [],
        "scope": "native_reference_and_finite_model_checks_not_exact_Hubbard_physics",
    }
    stage = "source_revision"
    try:
        for name in THREAD_VARIABLES:
            os.environ[name] = "1"
        identities = source_identity(root)
        identities["probe_sha256"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        receipt["identities"] = identities
        receipt["source_revision"] = identities["source_revision"]
        require(
            checks, "source_revision", expected_revision == receipt["source_revision"]
        )
        stage = "environment"
        receipt["environment"] = environment_identity(Path(sys.prefix))
        with provider_output_to_stderr():
            stage = "imports"
            np = importlib.import_module("numpy")
            scipy = importlib.import_module("scipy")
            gfs = importlib.import_module("triqs.gfs")
            meshes = importlib.import_module("triqs.mesh")
            mpi = importlib.import_module("triqs.utility.mpi")
            version = importlib.import_module("triqs.version")
            pools = importlib.import_module("threadpoolctl")
            modules = {
                "basis": importlib.import_module("kernelgecp.basis"),
                "adapter": importlib.import_module("kernelgecp.integrations.triqs"),
                "core": importlib.import_module("kernelgecp.dmft"),
                "reference": importlib.import_module(
                    "kernelgecp.integrations.triqs_dmft"
                ),
                "lehmann": importlib.import_module("kernelgecp.lehmann"),
            }
            checks["imports"] = True
            receipt["versions"].update(
                {
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "triqs": version.version,
                }
            )
            require(
                checks,
                "pinned_runtime",
                version.version == "4.0.2"
                and platform.python_version().startswith("3.12."),
            )
            paths = {}
            for name, relative in PACKAGE_FILES.items():
                path = Path(modules[name].__file__).resolve()
                paths[name] = path
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                identities[name + "_sha256"] = digest
                require(
                    checks,
                    name + "_matches_source",
                    digest
                    == hashlib.sha256((root / relative).read_bytes()).hexdigest(),
                )
                if require_installed:
                    require(
                        checks,
                        name + "_outside_checkout",
                        not path.is_relative_to(root),
                    )
            if require_installed:
                receipt["installed_wheel"] = _wheel_identity(paths)
                checks["installed_wheel"] = True
            require(checks, "one_rank", mpi.size == 1)
            with pools.threadpool_limits(limits=1):
                reference, core, lehmann = (
                    modules["reference"],
                    modules["core"],
                    modules["lehmann"],
                )
                stage = "canonical_loops"
                coarse_u2 = fine_u2 = None
                max_g_iw = max_s_iw = max_g_tau = max_s_tau = 0.0
                for interaction in (0, 1, 2):
                    config = core.IPTConfig(interaction=interaction)
                    coarse = reference.solve_dense_ipt(config, record_trajectory=True)
                    fine = reference.solve_dense_ipt(
                        core.IPTConfig(interaction=interaction, n_iw=1024, n_tau=12289)
                    )
                    require(
                        checks,
                        f"U{interaction}_converged",
                        coarse.converged and fine.converged,
                    )
                    require(
                        checks,
                        f"U{interaction}_two_unmixed_defects",
                        len(coarse.history) >= 2
                        and all(
                            row.unmixed_defect < 1e-8 for row in coarse.history[-2:]
                        ),
                    )
                    require(
                        checks,
                        f"U{interaction}_density_causality",
                        abs(coarse.density - 0.5) < 1e-8
                        and coarse.native_tail.positive_frequency_imag_g_max <= 1e-12
                        and coarse.native_tail.positive_frequency_imag_sigma_max
                        <= 1e-12,
                    )
                    if interaction:
                        moment_error = coarse.native_tail.sigma_first_moment_error
                        require(
                            checks,
                            f"U{interaction}_tail_moment",
                            moment_error is not None and moment_error < 1e-6,
                        )
                    max_g_iw = max(
                        max_g_iw,
                        float(np.max(np.abs(coarse.g_iw - fine.g_iw[512:1536]))),
                    )
                    max_s_iw = max(
                        max_s_iw,
                        float(
                            np.max(np.abs(coarse.sigma_iw - fine.sigma_iw[512:1536]))
                        ),
                    )
                    max_g_tau = max(
                        max_g_tau, float(np.max(np.abs(coarse.g_tau - fine.g_tau[::2])))
                    )
                    max_s_tau = max(
                        max_s_tau,
                        float(np.max(np.abs(coarse.sigma_tau - fine.sigma_tau[::2]))),
                    )
                    receipt["runs"].append(
                        {
                            "interaction": interaction,
                            "iterations": len(coarse.history),
                            "converged": coarse.converged,
                            "stop_reason": coarse.stop_reason,
                            "last_two_unmixed_defects": [
                                row.unmixed_defect for row in coarse.history[-2:]
                            ],
                            "density_per_spin": coarse.density,
                            "trajectory_count": len(coarse.trajectory),
                            "sigma_first_moment_error": (
                                coarse.native_tail.sigma_first_moment_error
                            ),
                        }
                    )
                    if interaction == 2:
                        coarse_u2, fine_u2 = coarse, fine
                metrics = receipt["metrics"]
                metrics.update(
                    {
                        "refinement_green_iw": max_g_iw,
                        "refinement_sigma_iw": max_s_iw,
                        "refinement_green_tau": max_g_tau,
                        "refinement_sigma_tau": max_s_tau,
                    }
                )
                require(
                    checks,
                    "native_refinement",
                    max(max_g_iw, max_s_iw, max_g_tau, max_s_tau) < 1e-7,
                )
                stage = "heldout_queries"
                heldout = np.arange(3, fine_u2.tau.size - 1, 4)
                heldout_error = max(
                    float(
                        np.max(
                            np.abs(
                                coarse_u2.evaluate_tau(
                                    fine_u2.tau[heldout], field=field
                                )
                                - values[heldout]
                            )
                        )
                    )
                    for field, values in (
                        ("green", fine_u2.g_tau),
                        ("sigma", fine_u2.sigma_tau),
                    )
                )
                labels = np.array([-701, -601, -513, 512, 600, 700])
                unseen_error = max(
                    float(
                        np.max(
                            np.abs(
                                coarse_u2.evaluate_matsubara(labels, field=field)
                                - fine_u2.evaluate_matsubara(labels, field=field)
                            )
                        )
                    )
                    for field in ("green", "sigma")
                )
                metrics.update(
                    {
                        "heldout_tau_error": heldout_error,
                        "unseen_frequency_error": unseen_error,
                    }
                )
                require(
                    checks,
                    "heldout_native_queries",
                    max(heldout_error, unseen_error) < 1e-7,
                )
                nodes, weights = np.polynomial.legendre.leggauss(256)
                tau_query = coarse_u2.config.beta * (nodes + 1) / 2
                quad_labels = np.array([-8, -2, 0, 6])
                iw = 1j * (2 * quad_labels + 1) * np.pi / coarse_u2.config.beta
                fourier_error = 0.0
                for field in ("green", "sigma"):
                    integral = (
                        np.exp(iw[:, None] * tau_query)
                        * (coarse_u2.config.beta * weights / 2)
                    ) @ coarse_u2.evaluate_tau(tau_query, field=field)
                    fourier_error = max(
                        fourier_error,
                        float(
                            np.max(
                                np.abs(
                                    integral
                                    - coarse_u2.evaluate_matsubara(
                                        quad_labels, field=field
                                    )
                                )
                            )
                        ),
                    )
                metrics["independent_fourier_error"] = fourier_error
                require(checks, "independent_fourier", fourier_error < 1e-7)
                stage = "independent_oracles"
                atomic = _pole_oracle(
                    np, gfs, meshes, lehmann.FermionicPoleModel(10, [0], [1]), 2
                )
                asymmetric = _pole_oracle(
                    np,
                    gfs,
                    meshes,
                    lehmann.FermionicPoleModel(10, [-1.3, 0.2, 1.1], [0.2, 0.5, 0.3]),
                    2,
                )
                metrics.update(
                    {
                        "atomic_dyson_error": atomic["atomic_dyson_error"],
                        "asymmetric_oracle_iw_error": asymmetric["iw_error"],
                        "asymmetric_oracle_tau_error": asymmetric["tau_error"],
                        "wrong_unreflected_cube_error": asymmetric[
                            "wrong_unreflected_cube_error"
                        ],
                    }
                )
                require(checks, "atomic_lehmann", atomic["atomic_dyson_error"] < 1e-8)
                require(
                    checks,
                    "asymmetric_triple_sum",
                    max(asymmetric["iw_error"], asymmetric["tau_error"]) < 1e-8
                    and asymmetric["wrong_unreflected_cube_error"] > 0.1,
                )
                finite_bath = lehmann.impurity_one_bath(10, 2, 0.4)
                native_iw = gfs.Gf(
                    mesh=meshes.MeshImFreq(beta=10, statistic="Fermion", n_iw=512),
                    target_shape=[],
                )
                native_tau = gfs.Gf(
                    mesh=meshes.MeshImTime(beta=10, statistic="Fermion", n_tau=6145),
                    target_shape=[],
                )
                native_iw.data[:] = finite_bath.evaluate_matsubara(np.arange(-512, 512))
                native_tau.set_from_fourier(
                    native_iw,
                    known_moments=np.r_[0, finite_bath.moments(2)].astype(complex),
                )
                finite_bath_error = float(
                    np.max(
                        np.abs(
                            native_tau.data - finite_bath.evaluate_tau(coarse_u2.tau)
                        )
                    )
                )
                metrics["finite_bath_tau_error"] = finite_bath_error
                require(checks, "finite_bath_lehmann", finite_bath_error < 1e-8)
                require(
                    checks,
                    "single_thread",
                    bool(pools.threadpool_info())
                    and all(
                        pool["num_threads"] == 1 for pool in pools.threadpool_info()
                    ),
                )
        require(
            checks,
            "numerical",
            all(value for name, value in checks.items() if name != "numerical"),
        )
        receipt["status"] = "passed"
    except Exception as error:
        receipt["error"] = {"stage": stage, "type": type(error).__name__}
        print(
            f"dense IPT probe failed at {stage}: {type(error).__name__}",
            file=sys.stderr,
        )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--require-installed", action="store_true")
    args = parser.parse_args()
    receipt = run_probe(args.source_revision, require_installed=args.require_installed)
    rendered = json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    sys.stdout.write(rendered)
    return 0 if receipt["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
