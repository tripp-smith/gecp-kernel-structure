#!/usr/bin/env python3
"""Fail-closed actual-TRIQS acceptance for the optional public adapter.

Run using the pinned provider environment, with this checkout installed or on
PYTHONPATH. A missing provider is a failed execution, never a skipped test.
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
from importlib import metadata as package_metadata
from pathlib import Path
from typing import Any

from triqs_smoke import (
    THREAD_VARIABLES,
    environment_identity,
    provider_output_to_stderr,
    require,
    source_identity,
)

SCHEMA = "kernelgecp.triqs-adapter-smoke.v1"
PACKAGE_FILES = {
    "basis": "kernelgecp/basis.py",
    "adapter": "kernelgecp/integrations/triqs.py",
}


def installed_wheel_identity(module_paths: dict[str, Path]) -> dict[str, Any]:
    """Bind imports to actual wheel metadata without exporting private paths."""
    try:
        distribution = package_metadata.distribution("kernelgecp")
    except package_metadata.PackageNotFoundError as error:
        raise RuntimeError(
            "installed kernelgecp distribution metadata is required"
        ) from error
    if set(module_paths) != set(PACKAGE_FILES):
        raise RuntimeError("installed module identity is incomplete")
    files = distribution.files
    record = distribution.read_text("RECORD")
    if files is None or not isinstance(record, str) or not record.strip():
        raise RuntimeError("installed kernelgecp wheel RECORD metadata is required")
    direct_url = distribution.read_text("direct_url.json")
    if direct_url is not None:
        try:
            direct = json.loads(direct_url)
        except (TypeError, ValueError) as error:
            raise RuntimeError("invalid installation metadata") from error
        if not isinstance(direct, dict):
            raise RuntimeError("invalid installation metadata")
        directory = direct.get("dir_info", {})
        if not isinstance(directory, dict) or (
            "editable" in directory and type(directory["editable"]) is not bool
        ):
            raise RuntimeError("invalid installation directory metadata")
        if directory.get("editable") is True:
            raise RuntimeError(
                "editable installations do not establish wheel acceptance"
            )
    try:
        rows = list(csv.reader(io.StringIO(record), strict=True))
    except csv.Error as error:
        raise RuntimeError("invalid wheel RECORD metadata") from error
    if any(len(row) != 3 for row in rows):
        raise RuntimeError("invalid wheel RECORD metadata")
    listed = {str(path) for path in files}
    for name, relative in PACKAGE_FILES.items():
        if relative not in listed or sum(row[0] == relative for row in rows) != 1:
            raise RuntimeError("imported module is absent from wheel file metadata")
        declared_path = Path(distribution.locate_file(relative)).resolve()
        if module_paths[name].resolve() != declared_path:
            raise RuntimeError("imported module does not belong to the installed wheel")
    return {
        "distribution": "kernelgecp",
        "version": distribution.version,
        "record_sha256": hashlib.sha256(record.encode()).hexdigest(),
        "module_count": len(PACKAGE_FILES),
        "editable": False,
    }


def run_probe(
    expected_revision: str, *, require_installed: bool = False
) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    checks = {"imports": False, "numerical": False}
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "provider": "triqs",
        "status": "failed",
        "versions": {"python": platform.python_version()},
        "checks": checks,
        "metrics": {},
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
            basis_module = importlib.import_module("kernelgecp.basis")
            adapter = importlib.import_module("kernelgecp.integrations.triqs")
            checks["imports"] = True
            receipt["versions"].update(
                {
                    "triqs": version.version,
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                }
            )
            require(checks, "pinned_triqs", version.version == "4.0.2")
            require(
                checks,
                "supported_python",
                platform.python_version().startswith("3.12."),
            )
            imported_modules: dict[str, Path] = {}
            for name, module in (("basis", basis_module), ("adapter", adapter)):
                imported_path = Path(module.__file__).resolve()
                imported_modules[name] = imported_path
                actual = hashlib.sha256(imported_path.read_bytes()).hexdigest()
                expected = hashlib.sha256(
                    (root / PACKAGE_FILES[name]).read_bytes()
                ).hexdigest()
                identities[name + "_sha256"] = actual
                require(checks, name + "_matches_source", actual == expected)
                outside_checkout = not imported_path.is_relative_to(root)
                receipt.setdefault("import_origins", {})[name] = (
                    "outside_checkout" if outside_checkout else "checkout"
                )
                if require_installed:
                    require(checks, name + "_installed", outside_checkout)
            if require_installed:
                stage = "installed_distribution"
                checks["installed_distribution"] = False
                receipt["installation"] = installed_wheel_identity(imported_modules)
                require(checks, "installed_distribution", True)
                receipt["import_origins"] = {
                    name: "installed" for name in imported_modules
                }
            require(checks, "one_rank", mpi.size == 1)
            with pools.threadpool_limits(limits=1):
                stage = "analytic_fixtures"
                beta = 4.0
                basis = basis_module.FermionicBasis(
                    beta=beta,
                    wmax=2,
                    tol=1e-12,
                    precision="reference",
                    grid_order=8,
                    max_rank=64,
                )
                energies = np.array([-0.8, 0.6])
                first = np.array([[0.7, 0.15 + 0.2j], [0.15 - 0.2j, 0.3]])
                residues = {
                    "scalar": np.array([0.4, 0.6]),
                    "matrix": np.array([first, np.eye(2) - first]),
                }
                train_mesh = meshes.MeshImTime(beta=beta, statistic="Fermion", n_tau=97)
                query_mesh = meshes.MeshImTime(beta=beta, statistic="Fermion", n_tau=98)
                tau = np.array([point.value for point in train_mesh])
                query = np.array([point.value for point in query_mesh])

                def time_kernel(rows: Any) -> Any:
                    return -np.exp(-rows[:, None] * energies) / (
                        1 + np.exp(-beta * energies)
                    )

                fits, sources = {}, {}
                signed_mesh = meshes.MeshImFreq(beta=beta, statistic="Fermion", n_iw=9)
                positive_mesh = meshes.MeshImFreq(
                    beta=beta,
                    statistic="Fermion",
                    n_iw=9,
                    opt="positive_frequencies_only",
                )
                signed_error = positive_error = 0.0
                for name, weights in residues.items():
                    source = gfs.Gf(
                        mesh=train_mesh, target_shape=list(weights.shape[1:])
                    )
                    source.data[...] = np.tensordot(
                        time_kernel(tau), weights, axes=(1, 0)
                    )
                    fitted = adapter.fit_triqs_tau(basis, source)
                    fits[name], sources[name] = fitted, source
                    output = adapter.to_triqs_tau(fitted, query_mesh)
                    truth = np.tensordot(time_kernel(query), weights, axes=(1, 0))
                    error = float(np.max(np.abs(output.data - truth)))
                    receipt["metrics"][name + "_tau_error"] = error
                    require(checks, name + "_tau", error <= 1e-8)
                    require(checks, name + "_shape", output.data.shape == truth.shape)
                    for mesh, positive in ((signed_mesh, False), (positive_mesh, True)):
                        indices = np.array([point.index for point in mesh])
                        iw = 1j * (2 * indices + 1) * np.pi / beta
                        target = np.tensordot(
                            1 / (iw[:, None] - energies), weights, axes=(1, 0)
                        )
                        predicted = adapter.to_triqs_iw(fitted, mesh)
                        error = float(np.max(np.abs(predicted.data - target)))
                        require(
                            checks,
                            name + ("_positive_iw" if positive else "_signed_iw"),
                            error <= 1e-8,
                        )
                        if positive:
                            positive_error = max(positive_error, error)
                        else:
                            signed_error = max(signed_error, error)
                receipt["metrics"].update(
                    {
                        "signed_iw_error": signed_error,
                        "positive_iw_error": positive_error,
                    }
                )

                stage = "complex_linearity"
                complex_residues = np.array(
                    [
                        [[1 + 2j, 3 + 4j], [-0.6j, 2 - 0.5j]],
                        np.zeros((2, 2)),
                    ]
                )
                arbitrary = gfs.Gf(mesh=train_mesh, target_shape=[2, 2])
                arbitrary.data[...] = np.tensordot(
                    time_kernel(tau), complex_residues, axes=(1, 0)
                )
                arbitrary_fit = adapter.fit_triqs_tau(basis, arbitrary)
                ordinary = np.array([point.index for point in signed_mesh])
                arbitrary_iw = 1j * (2 * ordinary + 1) * np.pi / beta
                arbitrary_truth = np.tensordot(
                    1 / (arbitrary_iw[:, None] - energies),
                    complex_residues,
                    axes=(1, 0),
                )
                complex_error = float(
                    np.max(
                        np.abs(
                            adapter.to_triqs_iw(arbitrary_fit, signed_mesh).data
                            - arbitrary_truth
                        )
                    )
                )
                receipt["metrics"]["complex_linearity_error"] = complex_error
                require(checks, "complex_linearity", complex_error <= 1e-8)

                stage = "endpoints_and_fourier"
                scalar = adapter.to_triqs_tau(fits["scalar"], query_mesh).data
                require(
                    checks,
                    "endpoint_distinction",
                    abs(complex(scalar[0] - scalar[-1])) > 1e-3
                    and abs(complex(scalar[0] + scalar[-1] + 1)) < 1e-8,
                )
                nodes, weights = np.polynomial.legendre.leggauss(96)
                quadrature_tau = beta * (nodes + 1) / 2
                indices = np.arange(-4, 4)
                iw = 1j * (2 * indices + 1) * np.pi / beta
                integral = np.tensordot(
                    np.exp(iw[:, None] * quadrature_tau) * (beta * weights / 2),
                    fits["matrix"].evaluate_tau(quadrature_tau),
                    axes=(1, 0),
                )
                frequencies = adapter.to_triqs_iw(fits["matrix"], signed_mesh).data[
                    5:13
                ]
                fourier_error = float(np.max(np.abs(integral - frequencies)))
                receipt["metrics"]["fourier_error"] = fourier_error
                require(checks, "independent_fourier", fourier_error <= 1e-8)
                require(
                    checks,
                    "matrix_adjoint_not_entrywise_conjugation",
                    np.max(np.abs(frequencies[3] - frequencies[4].conj().T)) < 1e-8
                    and np.max(np.abs(frequencies[3] - frequencies[4].conj())) > 1e-3,
                )

                stage = "block_recipe"
                one = gfs.Gf(mesh=train_mesh, target_shape=[1, 1])
                one.data[:, 0, 0] = sources["scalar"].data
                original = gfs.BlockGf(
                    name_list=["z_one", "a_two"],
                    block_list=[one, sources["matrix"]],
                    make_copies=True,
                )
                blocks = [
                    (name, adapter.fit_triqs_tau(basis, block))
                    for name, block in original
                ]
                result = gfs.BlockGf(
                    name_list=[name for name, _ in blocks],
                    block_list=[
                        adapter.to_triqs_iw(fit, signed_mesh) for _, fit in blocks
                    ],
                    make_copies=True,
                )
                require(
                    checks,
                    "ordered_unequal_blocks",
                    [name for name, _ in result] == ["z_one", "a_two"]
                    and result["z_one"].data.shape == (18, 1, 1)
                    and result["a_two"].data.shape == (18, 2, 2),
                )

                stage = "ownership_and_axes"
                expected = adapter.to_triqs_tau(fits["matrix"], query_mesh).data.copy()
                alternate = basis.fit_tau(
                    np.moveaxis(sources["matrix"].data, 0, 1), tau, axis=1
                )
                require(
                    checks,
                    "observation_axis",
                    np.max(
                        np.abs(
                            adapter.to_triqs_tau(alternate, query_mesh).data - expected
                        )
                    )
                    < 1e-12,
                )
                output = adapter.to_triqs_tau(fits["matrix"], query_mesh)
                sources["matrix"].data[...] = 0
                output.data[...] = 7j
                require(
                    checks,
                    "ownership",
                    np.max(
                        np.abs(
                            adapter.to_triqs_tau(fits["matrix"], query_mesh).data
                            - expected
                        )
                    )
                    == 0,
                )
                stage = "mesh_ownership"
                for domain in ("tau", "iw"):
                    if domain == "tau":
                        caller_mesh = meshes.MeshImTime(
                            beta=beta, statistic="Fermion", n_tau=13
                        )
                        replacement = meshes.MeshImTime(
                            beta=5, statistic="Fermion", n_tau=15
                        )
                        output_replacement = meshes.MeshImTime(
                            beta=6, statistic="Fermion", n_tau=17
                        )
                        exported = adapter.to_triqs_tau(fits["matrix"], caller_mesh)
                    else:
                        caller_mesh = meshes.MeshImFreq(
                            beta=beta, statistic="Fermion", n_iw=5
                        )
                        replacement = meshes.MeshImFreq(
                            beta=5, statistic="Fermion", n_iw=7
                        )
                        output_replacement = meshes.MeshImFreq(
                            beta=6, statistic="Fermion", n_iw=9
                        )
                        exported = adapter.to_triqs_iw(fits["matrix"], caller_mesh)
                    original_count = len(caller_mesh)
                    original_coordinates = np.array(
                        [complex(point.value) for point in caller_mesh]
                    )
                    original_data = exported.data.copy()
                    caller_mesh.copy_from(replacement)
                    detached_from_caller = (
                        exported.mesh.beta == beta
                        and len(exported.mesh) == original_count
                        and np.array_equal(
                            np.array([complex(point.value) for point in exported.mesh]),
                            original_coordinates,
                        )
                        and np.array_equal(exported.data, original_data)
                    )
                    caller_count = len(caller_mesh)
                    caller_coordinates = np.array(
                        [complex(point.value) for point in caller_mesh]
                    )
                    exported.mesh.copy_from(output_replacement)
                    detached_from_output = (
                        caller_mesh.beta == 5
                        and len(caller_mesh) == caller_count
                        and np.array_equal(
                            np.array([complex(point.value) for point in caller_mesh]),
                            caller_coordinates,
                        )
                        and np.array_equal(exported.data, original_data)
                    )
                    require(
                        checks,
                        domain + "_mesh_ownership",
                        detached_from_caller and detached_from_output,
                    )
                checks["mesh_ownership"] = True
                require(
                    checks,
                    "single_thread",
                    bool(pools.threadpool_info())
                    and all(
                        info["num_threads"] == 1 for info in pools.threadpool_info()
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
            f"TRIQS adapter probe failed at {stage}: {type(error).__name__}",
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
