"""Generate authentic official EPW anisotropic fixtures in the pinned image.

This first EG gate runs the solver; it does not close Nambu compression or
native reconstructed-restart acceptance. No local Docker daemon is required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.providers import epw_smoke as native

UPSTREAM_INPUT_SHA256 = (
    "a67e51f747674319d3c5659e402529bcdf3d571b54f12d2802f2531d74ba8a82"
)
UPSTREAM_IO_SHA256 = "416ee21c65fd49b2a04894c2d5e18ffb1a5e8f7c223cab950a28f760cf56ec36"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect_native_file(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    if not lines or "shift(w)" not in lines[0]:
        raise ValueError("expected the full-bandwidth native anisotropic header")
    blocks: list[list[list[float]]] = []
    for line in lines[1:]:
        if len(line) != 100:
            raise ValueError("expected five fixed-width ES20.10 native fields")
        values = [
            float(line[start : start + 20].replace("D", "E"))
            for start in range(0, 100, 20)
        ]
        if (
            not all(math.isfinite(x) for x in values)
            or values[0] == 0
            or values[2] <= 0
        ):
            raise ValueError("invalid frequency or nonfinite/nonpositive native Z")
        if not blocks or values[0] != blocks[-1][0][0]:
            blocks.append([])
        blocks[-1].append(values)
    if len(blocks) < 2 or not blocks[0]:
        raise ValueError("native anisotropic file is empty or frequency-truncated")
    frequencies = [block[0][0] for block in blocks]
    states = [row[1] for row in blocks[0]]
    if len(set(frequencies)) != len(frequencies) or any(
        left >= right for left, right in zip(frequencies, frequencies[1:], strict=False)
    ):
        raise ValueError("native frequency blocks must be unique and increasing")
    if any([row[1] for row in block] != states for block in blocks):
        raise ValueError("state ordinal/order changed between frequency blocks")
    return {
        "filename": path.name,
        "sha256": _sha(path),
        "bytes": path.stat().st_size,
        "frequency_count": len(blocks),
        "state_count": len(states),
        "frequency_ev": frequencies,
        "state_energy_minus_ef0_ev": states,
        "state_ordering": "native_frequency_major_then_ik_then_ibnd_ordinal",
        "columns": ["frequency_ev", "energy_minus_ef0_ev", "Z", "gap_ev", "shift_ev"],
        "positive_only": all(x > 0 for x in frequencies),
    }


def run(output: Path, revision: str, qe_root: Path) -> dict[str, Any]:
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        raise ValueError("a full source revision is required")
    output.parent.mkdir(parents=True, exist_ok=True)
    baseline = native.run(output.parent / "epw-baseline.json", revision, qe_root)
    (output.parent / "epw-baseline.json").write_text(
        json.dumps(baseline, indent=2, sort_keys=True) + "\n"
    )
    upstream = qe_root / "test-suite/epw_super/epw6.in"
    io_source = qe_root / "EPW/src/io/io_supercond.f90"
    if _sha(upstream) != UPSTREAM_INPUT_SHA256 or _sha(io_source) != UPSTREAM_IO_SHA256:
        raise ValueError("anisotropic source/input differs from pinned qe-7.5")
    work = output.parent / "epw-native/test-suite/epw_super"
    shutil.copy2(upstream, work / "epw6.in")
    text, execution = native._run_native(
        qe_root / "EPW/bin/epw.x", "epw6.in", work, 7200
    )
    if len(re.findall("Convergence was reached in nsiter", text)) < 3:
        raise RuntimeError("all three official anisotropic temperatures must converge")
    files = []
    for temperature in (17, 20, 30):
        path = work / f"MgB2.imag_aniso_{temperature:06.2f}"
        files.append({"temperature_kelvin": temperature, **inspect_native_file(path)})
    return {
        "schema": "kernelgecp.epw-native-anisotropic-fixture.v1",
        "status": "passed",
        "source_revision": revision,
        "provider": {"qe": "7.5", "epw": "6.0"},
        "script_sha256": _sha(Path(__file__)),
        "baseline": baseline,
        "native_anisotropic_execution": execution,
        "files": files,
        "upstream_input_sha256": UPSTREAM_INPUT_SHA256,
        "upstream_io_sha256": UPSTREAM_IO_SHA256,
        "chemical_potential_offset_status": "explicit_offset_required_before_Nambu",
        "scope": "actual_coarse_MgB2_fixture_generation_not_full_EG_acceptance",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--qe-root", type=Path, default=Path("/opt/qe"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.unlink(missing_ok=True)
    try:
        result = run(args.output, args.source_revision, args.qe_root)
    except Exception as error:
        args.output.write_text(
            json.dumps(
                {"status": "failed", "error": f"{type(error).__name__}: {error}"}
            )
            + "\n"
        )
        raise SystemExit(1) from error
    args.output.write_text(
        json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
