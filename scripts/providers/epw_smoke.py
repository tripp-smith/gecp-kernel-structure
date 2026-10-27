"""Run the pinned, native QE/EPW coarse MgB2 calculation; emit no fake passes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _copy_phonon_data(work: Path) -> int:
    """Collect PH output without upstream pp.py's unchecked shell copies."""
    prefix = "MgB2"
    phsave = work / "_ph0" / f"{prefix}.phsave"
    control = ET.parse(phsave / "control_ph.xml").getroot()
    counts = [
        int(element.text or "0")
        for element in control.iter()
        if element.tag.split("}")[-1] == "NUMBER_OF_Q_POINTS"
    ]
    _require(len(counts) == 1 and counts[0] >= 1, "missing phonon q-point count")
    saved = work / "save"
    saved.mkdir()
    shutil.copytree(phsave, saved / phsave.name)
    for qpoint in range(1, counts[0] + 1):
        dynamical = work / f"{prefix}.dyn{qpoint}"
        _require(dynamical.is_file(), "missing native dynamical matrix")
        shutil.copy2(dynamical, saved / f"{prefix}.dyn_q{qpoint}")
        directory = work / "_ph0"
        if qpoint > 1:
            directory /= f"{prefix}.q_{qpoint}"
        candidates = [
            candidate
            for name in (f"{prefix}.dvscf1", f"{prefix}.dvscf")
            if (candidate := directory / name).is_file()
        ]
        _require(len(candidates) == 1, "missing or ambiguous PH potential file")
        shutil.copy2(candidates[0], saved / f"{prefix}.dvscf_q{qpoint}")
    return counts[0]


def _run_native(
    executable: Path, input_file: str, work: Path, timeout_seconds: int
) -> tuple[str, dict[str, Any]]:
    output = work / input_file.replace(".in", ".out")
    environment = dict(os.environ)
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    # Direct execution of an MPI-enabled executable runs exactly one MPI rank.
    print(f"Running native {executable.name} -in {input_file}", flush=True)
    with output.open("w") as stream:
        process = subprocess.run(
            [str(executable), "-in", input_file],
            cwd=work,
            env=environment,
            stdout=stream,
            stderr=subprocess.STDOUT,
            timeout=timeout_seconds,
            check=False,
        )
    text = output.read_text(errors="replace")
    _require(process.returncode == 0, f"{input_file} failed; inspect {output}")
    _require("Error in routine" not in text, f"native error in {output}")
    return text, {
        "input": input_file,
        "input_sha256": _sha256(work / input_file),
        "log": output.name,
        "log_sha256": _sha256(output),
        "returncode": process.returncode,
        "executable_sha256": _sha256(executable),
        "mpi_ranks": 1,
        "omp_threads": 1,
    }


def _gap_rows(path: Path) -> list[list[float]]:
    lines = path.read_text().splitlines()
    _require(bool(lines) and "znorm" in lines[0], "invalid native gap header")
    rows = [
        [float(value.replace("D", "E")) for value in line.split()]
        for line in lines[1:]
        if line.strip()
    ]
    _require(len(rows) >= 8, "native calculation did not produce a frequency mesh")
    _require(
        all(
            len(row) == 3 and all(math.isfinite(value) for value in row) for row in rows
        ),
        "native imaginary-axis gap output is not finite three-column data",
    )
    _require(all(row[0] > 0 and row[1] > 0 for row in rows), "invalid frequency or Z")
    _require(
        all(left[0] < right[0] for left, right in zip(rows, rows[1:], strict=False)),
        "native frequencies are not strictly increasing",
    )
    _require(rows[0][2] > 0, "MgB2 15 K native gap was not positive")
    return rows


def run(output: Path, source_revision: str, qe_root: Path) -> dict[str, Any]:
    _require(bool(re.fullmatch(r"[0-9a-f]{40}", source_revision)), "invalid source SHA")
    repository = Path(__file__).resolve().parents[2]
    lock_path = repository / "environments/epw/sources.lock.json"
    lock = json.loads(lock_path.read_text())
    source_manifest_path = qe_root / "kernelgecp-source-manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text())
    _require(
        source_manifest["lock_sha256"] == _sha256(lock_path),
        "native image source lock differs from checked-out source lock; rebuild it",
    )
    _require(
        source_manifest["sources"] == lock["sources"]
        and source_manifest["pseudopotentials"] == lock["pseudopotentials"],
        "native source/pseudopotential identities disagree",
    )
    executables = {
        "pw": qe_root / "bin/pw.x",
        "ph": qe_root / "bin/ph.x",
        "epw": qe_root / "EPW/bin/epw.x",
    }
    _require(
        all(
            path.is_file() and os.access(path, os.X_OK) for path in executables.values()
        ),
        "required native executables are unavailable",
    )
    artifact_root = output.parent / "epw-native"
    # Keep the official input's ../../pseudo/ path and its content unchanged.
    work = artifact_root / "test-suite/epw_super"
    _require(not artifact_root.exists(), "native work directory already exists")
    work.mkdir(parents=True)
    pseudo_directory = artifact_root / "pseudo"
    pseudo_directory.mkdir()
    for item in lock["pseudopotentials"]:
        original = qe_root / "pseudo" / item["name"]
        _require(_sha256(original) == item["sha256"], "pseudopotential hash mismatch")
        shutil.copy2(original, pseudo_directory / item["name"])
    for name in lock["calculation"]["inputs"]:
        original = qe_root / lock["calculation"]["upstream_directory"] / name
        _require(
            _sha256(original) == source_manifest["native_input_sha256"][name],
            "upstream native input was modified",
        )
        shutil.copy2(original, work / name)
    shutil.copy2(source_manifest_path, artifact_root / source_manifest_path.name)
    build_provenance: dict[str, str] = {}
    for name in (
        "epw-system-packages.txt",
        "epw-configure.log",
        "epw-build.log",
    ):
        original = Path("/opt") / name
        if original.is_file():
            shutil.copy2(original, artifact_root / name)
            build_provenance[name] = _sha256(original)
    logs: list[dict[str, Any]] = []
    scf, identity = _run_native(executables["pw"], "scf.in", work, 1800)
    logs.append(identity)
    _require("convergence has been achieved" in scf, "SCF did not converge")
    _require("JOB DONE" in scf, "SCF did not finish normally")
    _require(bool(re.search(r"PWSCF\s+v\.7\.5\b", scf)), "unexpected native QE version")
    phonon, identity = _run_native(executables["ph"], "ph.in", work, 7200)
    logs.append(identity)
    _require("JOB DONE" in phonon, "native phonon calculation did not finish")
    phonon_qpoints = _copy_phonon_data(work)
    for name in ("scf_epw.in", "nscf_epw.in"):
        text, identity = _run_native(executables["pw"], name, work, 1800)
        logs.append(identity)
        _require("JOB DONE" in text, f"{name} did not finish normally")
    epw, identity = _run_native(executables["epw"], "epw1.in", work, 5400)
    logs.append(identity)
    _require(bool(re.search(r"EPW\s+v\.6\.0\b", epw)), "unexpected native EPW version")
    _require("Convergence was reached in nsiter" in epw, "EPW solve did not converge")
    matches = re.findall(r"Electron-phonon coupling strength\s*=\s*([0-9.Ee+-]+)", epw)
    _require(
        len(matches) == 1, "native electron-phonon coupling is missing or ambiguous"
    )
    coupling = float(matches[0])
    _require(math.isfinite(coupling) and 0.1 < coupling < 3, "invalid native coupling")
    gap_file = work / "MgB2.imag_iso_015.00"
    rows = _gap_rows(gap_file)
    checks = {
        "imports": True,
        "numerical": True,
        "source_lock_matches_image": True,
        "native_qe_7_5": True,
        "native_epw_6_0": True,
        "scf_converged": True,
        "phonon_calculation_completed": True,
        "official_inputs_unmodified": True,
        "native_eliashberg_converged": True,
        "finite_native_gap_mesh": True,
        "positive_mgb2_gap_at_15_kelvin": True,
    }
    return {
        "schema": "kernelgecp.provider-smoke.v1",
        "provider": "epw",
        "status": "passed",
        "source_revision": source_revision,
        "versions": {"python": platform.python_version(), "qe": "7.5", "epw": "6.0"},
        "checks": checks,
        "source_lock_sha256": _sha256(lock_path),
        "native_source_manifest_sha256": _sha256(source_manifest_path),
        "native_source": source_manifest,
        "build_provenance_sha256": build_provenance,
        "native_executions": logs,
        "calculation": {
            **lock["calculation"],
            "phonon_irreducible_qpoints": phonon_qpoints,
            "temperature_kelvin": 15,
            "electron_phonon_coupling": coupling,
            "native_frequency_count": len(rows),
            "lowest_frequency_gap_ev": rows[0][2],
            "gap_file": str(gap_file.relative_to(output.parent)),
            "gap_sha256": _sha256(gap_file),
        },
        "artifact_directory": str(artifact_root.relative_to(output.parent)),
        "certification_status": (
            "native_runtime_smoke_not_material_accuracy_certificate"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument(
        "--qe-root", type=Path, default=Path(os.environ.get("QE_ROOT", "/opt/qe"))
    )
    arguments = parser.parse_args()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        receipt = run(arguments.output, arguments.source_revision, arguments.qe_root)
    except Exception as error:
        receipt = {
            "schema": "kernelgecp.provider-smoke.v1",
            "provider": "epw",
            "status": "failed",
            "source_revision": arguments.source_revision,
            "versions": {"python": platform.python_version()},
            "checks": {"imports": False, "numerical": False},
            "error": str(error),
        }
        arguments.output.write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n"
        )
        print(f"EPW native smoke failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    arguments.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
