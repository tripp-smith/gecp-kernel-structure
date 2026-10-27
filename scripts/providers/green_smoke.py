"""Run an actual interacting molecular GF2 calculation in the Green image.

This provider-readiness gate does not benchmark kernelgecp or certify GF2's
physical approximation. A passed receipt is written only after native input
generation, self-consistency and independently inspected final deltas pass.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from numbers import Real
from pathlib import Path
from typing import Any

ENERGY_THRESHOLD = 1e-10
GREEN_CHANGE_THRESHOLD = 1e-6
MAX_ITERATIONS = 100


def validate_convergence(
    log: str,
    *,
    iteration: int,
    energy_change_sum: float,
    green_change: float,
) -> None:
    """Reject native exit-zero/max-iteration and inconsistent convergence claims."""
    if (
        "Simulation Converged" not in log
        or "Reached Maximum number of iterations" in log
    ):
        raise RuntimeError("Green did not report native self-consistency convergence")
    if (
        isinstance(iteration, bool)
        or not isinstance(iteration, int)
        or not 2 <= iteration <= MAX_ITERATIONS
    ):
        raise RuntimeError("a converged pair of native iterations is required")
    for value in (energy_change_sum, green_change):
        if isinstance(value, bool) or not isinstance(value, Real) or value < 0:
            raise RuntimeError("convergence deltas must be nonnegative real numbers")
        try:
            finite = math.isfinite(value)
        except OverflowError:
            finite = False
        if not finite:
            raise RuntimeError("convergence deltas must be finite")
    if energy_change_sum >= ENERGY_THRESHOLD:
        raise RuntimeError("stored energy-change sum exceeds the requested threshold")
    if green_change >= GREEN_CHANGE_THRESHOLD:
        raise RuntimeError("stored final Green-function change is too large")


def _run(command: list[str], work: Path, name: str, *, timeout: int) -> str:
    log_path = work / name
    with log_path.open("w") as output:
        subprocess.run(
            command,
            cwd=work,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
            env={
                **os.environ,
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMBA_NUM_THREADS": "1",
            },
        )
    return log_path.read_text()


def _inspect_result(path: Path, log: str) -> dict[str, Any]:
    import h5py
    import numpy as np

    required = (
        "G_tau/data",
        "G_tau/mesh",
        "Selfenergy/data",
        "Selfenergy/mesh",
        "Sigma1",
        "Energy_1b",
        "Energy_HF",
        "Energy_2b",
        "mu",
    )
    with h5py.File(path, "r") as results:
        iteration = int(results["iter"][()])
        current = results[f"iter{iteration}"]
        previous = results[f"iter{iteration - 1}"]
        if any(key not in current or key not in previous for key in required):
            raise RuntimeError("selected final iteration pair is incomplete")
        green = np.asarray(current["G_tau/data"])
        previous_green = np.asarray(previous["G_tau/data"])
        selfenergy = np.asarray(current["Selfenergy/data"])
        tau = np.asarray(current["G_tau/mesh"])
        if (
            green.ndim != 5
            or green.shape[1:] != (1, 1, 2, 2)
            or green.shape != previous_green.shape
            or green.shape != selfenergy.shape
            or tau.shape != (green.shape[0],)
        ):
            raise RuntimeError("unexpected native restricted H2 molecular data layout")
        if not all(
            np.all(np.isfinite(x)) for x in (green, previous_green, selfenergy, tau)
        ):
            raise RuntimeError("nonfinite native molecular data")
        energy_fields = ("Energy_1b", "Energy_HF", "Energy_2b")
        energy_change = sum(
            abs(float(current[field][()]) - float(previous[field][()]))
            for field in energy_fields
        )
        green_change = float(np.max(np.abs(green - previous_green)))
        validate_convergence(
            log,
            iteration=iteration,
            energy_change_sum=energy_change,
            green_change=green_change,
        )
        correlation = float(current["Energy_2b"][()])
        dynamic_max = float(np.max(np.abs(selfenergy)))
        if (
            not math.isfinite(correlation)
            or abs(correlation) <= 1e-8
            or dynamic_max <= 1e-8
        ):
            raise RuntimeError("native calculation has no interacting GF2 result")
        normalization_error = float(
            np.max(np.abs(green[0, 0, 0] + green[-1, 0, 0] + np.eye(2)))
        )
        electron_count = float(-2 * np.trace(green[-1, 0, 0]).real)
        if normalization_error > 1e-5 or abs(electron_count - 2) > 1e-5:
            raise RuntimeError("H2 endpoint normalization or particle count failed")
        values = {
            "iteration": iteration,
            "energy_change_sum": energy_change,
            "green_change_max": green_change,
            "correlation_energy_hartree": correlation,
            "total_energy_hartree": float(current["Energy_HF"][()]) + correlation,
            "chemical_potential_hartree": float(current["mu"][()]),
            "dynamic_selfenergy_max": dynamic_max,
            "endpoint_normalization_error": normalization_error,
            "electron_count": electron_count,
            "tau_sample_count": len(tau),
            "data_shape": list(green.shape),
        }
        if any(
            not math.isfinite(value)
            for value in values.values()
            if isinstance(value, float)
        ):
            raise RuntimeError("nonfinite native molecular measurements")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.source_revision):
        parser.error("--source-revision must be a full lowercase Git revision")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    # A failed rerun must not leave an older successful receipt at the requested path.
    output.unlink(missing_ok=True)
    work = output.parent / "green"
    work.mkdir(exist_ok=False)

    import green_mbtools.mint  # noqa: F401
    import h5py
    import numpy
    import pyscf

    if (
        pyscf.__version__ != "2.14.0"
        or importlib.metadata.version("green-mbtools") != "1.0.0"
    ):
        raise RuntimeError("the pinned Green/PySCF provider environment is required")
    lock_path = Path("/opt/green-lock/sources.json")
    lock = json.loads(lock_path.read_text())
    executable = Path("/opt/green/bin/mbpt.exe")
    if lock["release"] != "1.0.0" or not executable.is_file():
        raise RuntimeError("the pinned native Green executable is required")
    (work / "atom.dat").write_text("H 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
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
    grid = Path("/opt/green-src/green-grids/data/ir/1e4.h5")
    native_command = [
        str(executable),
        "--input_file",
        "input.h5",
        "--results_file",
        "sim.h5",
        "--scf_type",
        "GF2",
        # The molecular generator writes the bare Coulomb vertices once.
        # Both HF and GF2 consume that same df_hf_int directory.
        "--dfintegral_hf_file",
        "df_hf_int",
        "--dfintegral_file",
        "df_hf_int",
        "--BETA",
        "20",
        "--grid_file",
        str(grid),
        "--itermax",
        str(MAX_ITERATIONS),
        "--E_thr",
        str(ENERGY_THRESHOLD),
        "--mixing_type",
        "SIGMA_MIXING",
        "--mixing_weight",
        "0.5",
        "--jobs",
        "SC",
    ]
    native_log = _run(native_command, work, "gf2.log", timeout=1800)
    measurements = _inspect_result(work / "sim.h5", native_log)
    shutil.copyfile(
        "/opt/green-lock/build-provenance.json", work / "build-provenance.json"
    )
    manifest = {
        "molecule": "H2",
        "bond_length_angstrom": 0.74,
        "orbital_basis": "sto-3g",
        "storage_basis": "symmetric_lowdin",
        "beta_inverse_hartree": 20.0,
        "method": "GF2",
        "native_grid_lambda": 1e4,
        "input_command": input_command,
        "native_command": native_command,
        "energy_change_threshold": ENERGY_THRESHOLD,
        "green_change_threshold": GREEN_CHANGE_THRESHOLD,
        "scope": "native_provider_readiness_without_physical_accuracy_certificate",
        "source_lock": lock,
        "measurements": measurements,
    }
    (work / "run-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    artifacts = {}
    for name in (
        "input.log",
        "gf2.log",
        "input.h5",
        "sim.h5",
        "run-manifest.json",
        "build-provenance.json",
    ):
        file = work / name
        artifacts[f"green/{name}"] = hashlib.sha256(file.read_bytes()).hexdigest()
    receipt = {
        "schema": "kernelgecp.provider-smoke.v1",
        "provider": "green",
        "status": "passed",
        "source_revision": args.source_revision,
        "versions": {
            "python": platform.python_version(),
            "green": lock["release"],
            "green-mbtools": importlib.metadata.version("green-mbtools"),
            "pyscf": pyscf.__version__,
            "numpy": numpy.__version__,
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
        "artifacts": artifacts,
    }
    output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps(receipt, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
