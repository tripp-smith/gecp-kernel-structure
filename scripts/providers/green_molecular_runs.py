"""Run the four bounded native EF H2 GF2 calculations in the pinned image.

These native calculations generate scientific data and actual convergence
receipts. Compression, native-grid accuracy and FCI comparisons are subsequent
independent gates; successful native execution alone never closes EF.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import re
import shutil
import sys
from pathlib import Path
from typing import Any

from scripts.providers.bind_green_manifest import bind_manifest, sha256
from scripts.providers.green_smoke import (
    ENERGY_THRESHOLD,
    GREEN_CHANGE_THRESHOLD,
    _inspect_result,
    _run,
)

BONDS = (0.74, 1.5)
GRIDS = ("1e4", "1e5")
ITERATION_BUDGET = 300


def native_case(root: Path, bond: float, label: str, revision: str) -> dict[str, Any]:
    import h5py
    import numpy as np
    import pyscf

    if (
        pyscf.__version__ != "2.14.0"
        or importlib.metadata.version("green-mbtools") != "1.0.0"
    ):
        raise RuntimeError("the pinned actual Green/PySCF provider is required")
    work = root / f"h2-{bond}-{label}"
    work.mkdir()
    (work / "atom.dat").write_text(f"H 0.0 0.0 0.0\nH 0.0 0.0 {bond}\n")
    build_file = Path("/opt/green-lock/build-provenance.json")
    lock = json.loads(Path("/opt/green-lock/sources.json").read_text())
    if lock.get("release") != "1.0.0":
        raise RuntimeError(
            "a different native provider source cannot satisfy this campaign"
        )
    input_command = [
        sys.executable,
        "/opt/green-src/green-mbpt/python/init_data_mol_df.py",
        "--atom",
        "atom.dat",
        "--basis",
        "sto-3g",
        "--restricted",
        "true",
        "--orth",
        "symmetric_lowdin",
        "--memory",
        "512",
    ]
    _run(input_command, work, "input.log", timeout=900)
    grid = Path("/opt/green-src/green-grids/data/ir") / (label + ".h5")
    native_command = [
        "/opt/green/bin/mbpt.exe",
        "--input_file",
        "input.h5",
        "--results_file",
        "sim.h5",
        "--scf_type",
        "GF2",
        "--dfintegral_hf_file",
        "df_hf_int",
        "--dfintegral_file",
        "df_hf_int",
        "--BETA",
        "20",
        "--grid_file",
        str(grid),
        "--itermax",
        str(ITERATION_BUDGET),
        "--E_thr",
        str(ENERGY_THRESHOLD),
        "--mixing_type",
        "SIGMA_MIXING",
        "--mixing_weight",
        "0.5",
        "--jobs",
        "SC",
    ]
    log = _run(native_command, work, "gf2.log", timeout=1800)
    measurements = _inspect_result(
        work / "sim.h5", log, max_iterations=ITERATION_BUDGET
    )
    shutil.copyfile(build_file, work / "build-provenance.json")
    shutil.copyfile(grid, work / "native-grid.h5")
    settings = {
        "molecule": "H2",
        "bond_length_angstrom": bond,
        "orbital_basis": "sto-3g",
        "storage_basis": "symmetric_lowdin",
        "beta_inverse_hartree": 20.0,
        "method": "GF2",
        "native_grid_lambda": float(label),
        "input_command": input_command,
        "native_command": native_command,
        "energy_change_threshold": ENERGY_THRESHOLD,
        "green_change_threshold": GREEN_CHANGE_THRESHOLD,
        "iteration_budget": ITERATION_BUDGET,
        "source_lock": lock,
        "measurements": measurements,
        "scope": "actual_native_GF2_without_compression_or_FCI_accuracy_claim",
    }
    (work / "run-manifest.json").write_text(
        json.dumps(settings, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    names = (
        "input.log",
        "gf2.log",
        "input.h5",
        "sim.h5",
        "native-grid.h5",
        "run-manifest.json",
        "build-provenance.json",
        "df_hf_int/meta.h5",
        "df_hf_int/VQ_0.h5",
    )
    receipt = {
        "schema": "kernelgecp.provider-smoke.v1",
        "provider": "green",
        "status": "passed",
        "source_revision": revision,
        "versions": {
            "python": platform.python_version(),
            "green": lock["release"],
            "green-mbtools": importlib.metadata.version("green-mbtools"),
            "pyscf": pyscf.__version__,
            "numpy": np.__version__,
            "h5py": h5py.__version__,
        },
        "checks": {
            "imports": True,
            "numerical": True,
            "native_molecular_input_generation": True,
            "native_gf2_self_consistency": True,
            "stored_final_energy_and_green_deltas": True,
            "nonzero_interacting_selfenergy": True,
            "endpoint_normalization_and_particle_count": True,
        },
        "measurements": measurements,
        "artifacts": {"green/" + name: sha256(work / name) for name in names},
    }
    receipt_path = work / "native-receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    manifest_path = work / "bound-manifest.json"
    manifest = bind_manifest(work, receipt_path, manifest_path)
    manifest["files"]["native_grid"] = {
        "path": "native-grid.h5",
        "sha256": sha256(work / "native-grid.h5"),
    }
    for key, filename in (
        ("df_metadata", "df_hf_int/meta.h5"),
        ("df_chunk", "df_hf_int/VQ_0.h5"),
    ):
        manifest["files"][key] = {"path": filename, "sha256": sha256(work / filename)}
    manifest["native_grid_source_revision"] = next(
        row["revision"] for row in lock["sources"] if row["name"] == "green-grids"
    )
    manifest_path.write_text(
        json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    return {
        "bond_length_angstrom": bond,
        "native_grid_lambda": float(label),
        "directory": work.name,
        "manifest_sha256": sha256(manifest_path),
        "provider_receipt_sha256": sha256(work / "provider-receipt.json"),
        "measurements": measurements,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    if not re.fullmatch("[0-9a-f]{40}", args.source_revision):
        parser.error("a full lowercase source revision is required")
    if any(
        os.environ.get(name) != "1"
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    ):
        parser.error("native numerical thread counts must explicitly equal one")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    cases = [
        native_case(args.output_dir, bond, label, args.source_revision)
        for bond in BONDS
        for label in GRIDS
    ]
    report = {
        "schema": "kernelgecp.green-molecular-native.v1",
        "status": "passed",
        "source_revision": args.source_revision,
        "cases": cases,
        "scope": "four_actual_converged_native_runs_without_later_accuracy_claims",
    }
    (args.output_dir / "native-campaign.json").write_text(
        json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"status": "passed", "native_cases": len(cases)}))


if __name__ == "__main__":
    main()
