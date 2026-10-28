"""Bind retained native H2 provider artifacts to the supported reader manifest.

This command enriches the explicit settings already saved by green_smoke.py;
it does not run Green, infer geometry from serialized Cell, or certify GF2.
The source provider receipt must bind every original native artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bind_manifest(
    native_directory: Path, receipt_path: Path, output: Path
) -> dict[str, Any]:
    import h5py

    if output.resolve().parent != native_directory.resolve():
        raise ValueError(
            "write the bound manifest beside the retained native artifacts"
        )
    settings = json.loads((native_directory / "run-manifest.json").read_text())
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("schema") != "kernelgecp.provider-smoke.v1"
        or receipt.get("provider") != "green"
        or receipt.get("status") != "passed"
        or settings.get("molecule") != "H2"
        or settings.get("method") != "GF2"
        or settings.get("storage_basis") != "symmetric_lowdin"
    ):
        raise ValueError("a passed actual molecular H2 GF2 provider run is required")
    artifact_hashes = receipt.get("artifacts", {})
    required = {
        "input": "input.h5",
        "results": "sim.h5",
        "native_log": "gf2.log",
        "build_provenance": "build-provenance.json",
        "run_settings": "run-manifest.json",
    }
    files = {}
    for name, filename in required.items():
        actual = sha256(native_directory / filename)
        if not any(
            Path(key).name == filename and digest == actual
            for key, digest in artifact_hashes.items()
        ):
            raise ValueError("source provider receipt does not bind " + filename)
        files[name] = {"path": filename, "sha256": actual}
    destination = native_directory / "provider-receipt.json"
    if destination.exists() and sha256(destination) != sha256(receipt_path):
        raise ValueError("a different provider receipt already exists beside this run")
    if receipt_path.resolve() != destination.resolve():
        shutil.copyfile(receipt_path, destination)
    files["provider_receipt"] = {
        "path": destination.name,
        "sha256": sha256(destination),
    }
    build = json.loads((native_directory / "build-provenance.json").read_text())
    if build.get("sources") != settings.get("source_lock"):
        raise ValueError("run and build source identities differ")
    provider = dict(receipt["versions"])
    with h5py.File(native_directory / "sim.h5", "r") as source:
        provider["green-grids"] = str(source.attrs["__grids_version__"])
    provider["native_executable_sha256"] = build["native_executable_sha256"]
    provider["source_revisions"] = {
        row["name"]: row["revision"] for row in build["sources"]["sources"]
    }
    bond = settings["bond_length_angstrom"]
    manifest = {
        "schema": "kernelgecp.green-molecular-run.v1",
        "molecule": settings["molecule"],
        "orbital_basis": settings["orbital_basis"],
        "storage_basis": settings["storage_basis"],
        "method": settings["method"],
        "restricted": True,
        "beta_inverse_hartree": settings["beta_inverse_hartree"],
        "geometry": [
            {"element": "H", "xyz": [0.0, 0.0, 0.0]},
            {"element": "H", "xyz": [0.0, 0.0, bond]},
        ],
        "units": {
            "energy": "hartree",
            "imaginary_time": "hartree^-1",
            "geometry": "angstrom",
            "green_tau": "dimensionless",
            "dynamic_selfenergy_tau": "hartree^2",
        },
        "provider": provider,
        "files": files,
        "convergence": {
            "iteration": settings["measurements"]["iteration"],
            "energy_change_threshold_hartree": settings["energy_change_threshold"],
            "green_change_threshold": settings["green_change_threshold"],
        },
        "native_grid_lambda": settings["native_grid_lambda"],
        "native_command": settings["native_command"],
        "input_command": settings["input_command"],
        "source_revision": receipt["source_revision"],
        "scope": "native_run_provenance_without_GF2_model_accuracy_certificate",
    }
    output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_directory", type=Path)
    parser.add_argument("provider_receipt", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.native_directory / "bound-manifest.json"
    manifest = bind_manifest(args.native_directory, args.provider_receipt, output)
    print(
        json.dumps(
            {
                "manifest": str(output),
                "sha256": sha256(output),
                "source_revision": manifest["source_revision"],
            }
        )
    )


if __name__ == "__main__":
    main()
