"""Compare unmodified EPW 6 imag_read/Padé restarts of a frozen 17 K or 20 K input.

This is coarse native interoperability evidence, not an EPW acceleration claim
or an experimental materials prediction. No self-consistency iteration is run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.providers.epw_anisotropic_fixture import (
    SCOUT_INPUT_SHA256,
    UPSTREAM_INPUT_SHA256,
    UPSTREAM_IO_SHA256,
    inspect_native_file,
)

RAW_NAME = "MgB2.imag_aniso_017.00"
OUTPUT_NAME = "MgB2.pade_aniso_017.00"
INPUT_NAME = "epw6.in"
NQSTEP = 500
ATOL = 1e-5
RTOL = 1e-3
CHANNELS = (
    "Z_real",
    "Z_imag",
    "gap_real_ev",
    "gap_imag_ev",
    "shift_real_ev",
    "shift_imag_ev",
)
MISMATCH = "temperature not the same with the input"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def native_names(temperature_kelvin: int) -> tuple[str, str, str]:
    require(
        type(temperature_kelvin) is int and temperature_kelvin in (17, 20),
        "temperature must be 17 or 20 K",
    )
    suffix = f"{temperature_kelvin:03d}.00"
    return (
        f"MgB2.imag_aniso_{suffix}",
        f"MgB2.pade_aniso_{suffix}",
        f"MgB2.pade_aniso_gap0_{suffix}",
    )


def restart_input(text: str, temperature_kelvin: int = 17) -> str:
    native_names(temperature_kelvin)
    changes = (
        (
            r"(?m)^(\s*temps\s*=\s*)17\s+20\s+30(?=\s*!|\s*$)",
            rf"\g<1>{temperature_kelvin}",
        ),
        (r"(?m)^(\s*iverbosity\s*=\s*)0(?=\s|$)", r"\g<1>2"),
    )
    require(not re.search(r"(?im)^\s*imag_read\s*=", text), "imag_read already present")
    for pattern, replacement in changes:
        text, count = re.subn(pattern, replacement, text)
        require(count == 1, "expected exactly one official input setting")
    text, count = re.subn(
        r"(?m)^([ \t]*limag[ \t]*=[ \t]*\.true\.[ \t]*)$",
        r"\1\n  imag_read = .true.",
        text,
        flags=re.IGNORECASE,
    )
    require(count == 1, "expected exactly one enabled limag")
    require(
        bool(re.search(r"(?m)^\s*wscut\s*=\s*1\.0\s*!", text)), "expected 1 eV wscut"
    )
    return text


def wrong_frequency(source: Path, destination: Path, states: int) -> dict[str, Any]:
    lines = source.read_text().splitlines(keepends=True)
    require(len(lines) > states and states > 0, "missing first frequency block")
    index = states
    line = lines[index]
    require(len(line.rstrip("\r\n")) == 100, "malformed negative-control row")
    frequency = float(line[:20].replace("D", "E"))
    replacement = f"{frequency * 1.1:20.10E}"
    require(len(replacement) == 20 and frequency != 0, "invalid frequency corruption")
    lines[index] = replacement + line[20:]
    destination.write_text("".join(lines))
    return {
        "line_1_based": index + 1,
        "state_ordinal_1_based": states,
        "original_frequency_ev": frequency,
        "replacement_frequency_ev": float(replacement),
        "change": "last state of first frequency block, frequency times 1.1",
    }


def expected_mesh(a2f: Path) -> tuple[int, float]:
    lines = a2f.read_text().splitlines()
    require(len(lines) > NQSTEP and "w[meV]" in lines[0], "missing native a2f mesh")
    values = [float(line.split()[0]) for line in lines[1 : NQSTEP + 1]]
    require(all(math.isfinite(v) and v > 0 for v in values), "invalid a2f frequencies")
    require(
        all(b > a for a, b in zip(values, values[1:], strict=False)),
        "unordered a2f frequencies",
    )
    # Pinned supercond.f90 eliashberg_grid: nsw=nqstep*NINT(wscut/wsphmax).
    # The last a2f coordinate prints wsphmax in meV; avoid ambiguous rounding.
    ratio = 1000.0 / values[-1]
    require(
        abs(ratio - (math.floor(ratio) + 0.5)) > 1e-5, "ambiguous native nsw rounding"
    )
    nsw = NQSTEP * math.floor(ratio + 0.5)
    require(nsw >= NQSTEP, "invalid real frequency count")
    return nsw, values[-1] / (1000.0 * NQSTEP)


def execute(executable: Path, work: Path, timeout: int) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.update(
        {
            k: "1"
            for k in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        }
    )
    log = work / "restart.out"
    start = time.monotonic()
    status: dict[str, Any] = {
        "log": str(log),
        "command": [str(executable), "-in", INPUT_NAME],
        "mpi_ranks": 1,
        "threads": 1,
    }
    try:
        with log.open("xb") as stream:
            proc = subprocess.run(
                status["command"],
                cwd=work,
                env=environment,
                stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
        status.update(returncode=proc.returncode, timed_out=False)
    except subprocess.TimeoutExpired:
        status.update(returncode=None, timed_out=True)
    status["wall_seconds"] = time.monotonic() - start
    status["log_sha256"] = sha(log)
    status["log_bytes"] = log.stat().st_size
    return status


def check_execution(execution: dict[str, Any], *, negative: bool) -> None:
    text = Path(execution["log"]).read_text(errors="replace")
    require(
        execution["log_bytes"] > 0 and not execution["timed_out"],
        "native log empty or timeout",
    )
    if negative:
        require(
            execution["returncode"] not in (None, 0), "wrong-frequency run did not fail"
        )
        require(
            "eliashberg_read_aniso_iaxis" in text.lower() and MISMATCH in text.lower(),
            "wrong-frequency error not observed",
        )
    else:
        require(execution["returncode"] == 0, "native process failed")
        require("error in routine" not in text.lower(), "native error in log")
        require(
            re.search(r"EPW\s+v\.6\.0\b", text) is not None, "missing EPW 6 version"
        )
        require(
            "Pade approximant of full-bandwidth anisotropic" in text
            and "Convergence was reached for N =" in text
            and "raxis_pade" in text
            and "Total program execution" in text,
            "native Padé completion missing",
        )


def compare_tables(
    raw: Path, fitted: Path, states: list[float], nfreq: int, step: float
) -> dict[str, Any]:
    require(states and nfreq >= 2 and step > 0, "invalid expected axes")
    metrics = {
        name: {
            "max_absolute_error": 0.0,
            "max_budget_ratio": 0.0,
            "failures": 0,
            "unit": "dimensionless" if name.startswith("Z_") else "eV",
        }
        for name in CHANNELS
    }
    with raw.open() as a, fitted.open() as b:
        headers = [a.readline(), b.readline()]
        required = (
            "Re[znorm(w)]",
            "Im[znorm(w)]",
            "Re[delta(w)]",
            "Im[delta(w)]",
            "Re[shift(w)]",
            "Im[shift(w)]",
        )
        require(
            all(
                h.lstrip().startswith("#") and all(word in h for word in required)
                for h in headers
            ),
            "missing real-axis header",
        )
        axis = []
        for state, energy in enumerate(states):
            for iw in range(nfreq):
                lines = [a.readline().rstrip("\r\n"), b.readline().rstrip("\r\n")]
                require(
                    all(len(line) == 160 for line in lines),
                    f"truncated/malformed real-axis row {state}:{iw}",
                )
                fields = [
                    [line[i : i + 20] for i in range(0, 160, 20)] for line in lines
                ]
                require(fields[0][:2] == fields[1][:2], f"axis mismatch {state}:{iw}")
                rows = [[float(x.replace("D", "E")) for x in row] for row in fields]
                require(
                    all(math.isfinite(v) for row in rows for v in row),
                    "nonfinite real-axis row",
                )
                require(rows[0][1] == energy, f"state order mismatch {state}")
                if state == 0:
                    require(
                        abs(rows[0][0] - (iw + 1) * step) <= 6e-10,
                        f"real frequency grid mismatch {iw}",
                    )
                    axis.append(fields[0][0])
                else:
                    require(
                        fields[0][0] == axis[iw],
                        "real frequency grid changed across states",
                    )
                for j, name in enumerate(CHANNELS, start=2):
                    diff = abs(rows[1][j] - rows[0][j])
                    budget = ATOL + RTOL * abs(rows[0][j])
                    item = metrics[name]
                    item["max_absolute_error"] = max(item["max_absolute_error"], diff)
                    item["max_budget_ratio"] = max(
                        item["max_budget_ratio"], diff / budget
                    )
                    item["failures"] += int(diff > budget)
        require(not a.readline() and not b.readline(), "extra real-axis rows")
    return {
        "state_count": len(states),
        "frequency_count": nfreq,
        "row_count": len(states) * nfreq,
        "atol": ATOL,
        "rtol": RTOL,
        "budget_reference": "absolute raw value",
        "channels": metrics,
        "passed": all(v["failures"] == 0 for v in metrics.values()),
    }


def run(
    workspace: Path,
    raw_file: Path,
    reconstructed_file: Path,
    revision: str,
    qe_root: Path,
    output: Path,
    timeout: int,
    *,
    temperature_kelvin: int = 17,
) -> dict[str, Any]:
    receipt: dict[str, Any] = {
        "schema": "kernelgecp.epw-restart-comparison.v1",
        "status": "failed",
        "source_revision": revision,
        "temperature_kelvin": temperature_kelvin,
        "provider": {"qe": "7.5", "epw": "6.0"},
        "scope": (
            "coarse native imag_read/Pade interoperability, "
            "not material accuracy or acceleration"
        ),
        "self_consistency_iteration": False,
        "native_timeout_seconds_per_case": timeout,
        "tolerances": {
            "atol": ATOL,
            "rtol": RTOL,
            "gap_shift_units": "eV",
            "Z_units": "dimensionless",
        },
        "cases": {},
        "errors": [],
    }
    try:
        raw_name, output_name, gap_name = native_names(temperature_kelvin)
        require(
            re.fullmatch(r"[0-9a-f]{40}", revision) is not None,
            "full source revision required",
        )
        require(1 <= timeout <= 1200, "native timeout must be in [1,1200] seconds")
        workspace = workspace.resolve()
        output = output.resolve()
        qe_root = qe_root.resolve()
        require(
            workspace.is_dir() and not output.is_relative_to(workspace),
            "output must be outside native workspace",
        )
        executable = qe_root / "EPW/bin/epw.x"
        official = qe_root / "test-suite/epw_super/epw6.in"
        io = qe_root / "EPW/src/io/io_supercond.f90"
        require(
            sha(official) == UPSTREAM_INPUT_SHA256 and sha(io) == UPSTREAM_IO_SHA256,
            "pinned QE source mismatch",
        )
        require(
            executable.is_file() and os.access(executable, os.X_OK),
            "missing native executable",
        )
        original = workspace / INPUT_NAME
        require(
            sha(original) == SCOUT_INPUT_SHA256, "expected official wscut=1 eV input"
        )
        source_text = original.read_text()
        modified = restart_input(source_text, temperature_kelvin)
        receipt.update(
            source={
                "binary_sha256": sha(executable),
                "io_source_sha256": sha(io),
                "official_input_sha256": sha(official),
                "script_sha256": sha(Path(__file__)),
            },
            input={
                "original_sha256": sha(original),
                "modified_sha256": hashlib.sha256(modified.encode()).hexdigest(),
                "changes": [
                    f"temps: 17 20 30 -> {temperature_kelvin} K",
                    "imag_read: default false -> true",
                    "iverbosity: 0 -> 2 (output diagnostics)",
                ],
                "wscut_ev": 1.0,
            },
        )
        raw_meta = inspect_native_file(raw_file)
        fitted_meta = inspect_native_file(reconstructed_file)
        for key in (
            "frequency_count",
            "state_count",
            "frequency_ev",
            "state_energy_minus_ef0_ev",
            "columns",
        ):
            require(
                raw_meta[key] == fitted_meta[key],
                f"reconstructed input axis mismatch: {key}",
            )
        require(raw_meta["positive_only"], "expected positive Matsubara frequencies")
        receipt["inputs"] = {"raw": raw_meta, "reconstructed": fitted_meta}
        nfreq, step = expected_mesh(workspace / "MgB2.a2f")
        receipt["expected_real_grid"] = {
            "frequency_count": nfreq,
            "step_ev_derived_from_rounded_a2f": step,
            "a2f_sha256": sha(workspace / "MgB2.a2f"),
            "rule": "nsw = 500*NINT(1 eV/wsphmax), ws(i)=i*wsphmax/500",
        }
        case_root = output.parent / (output.stem + "-restart-cases")
        case_root.mkdir(parents=True, exist_ok=False)
        receipt["case_root"] = str(case_root)
        for name in ("raw", "reconstructed", "wrong-frequency"):
            work = case_root / name
            shutil.copytree(workspace, work)
            # Discard copied outputs so only this execution can supply evidence.
            for stale in (output_name, gap_name, "restart.out"):
                (work / stale).unlink(missing_ok=True)
            (work / INPUT_NAME).write_text(modified)
            target = work / raw_name
            if name == "wrong-frequency":
                corrupt = wrong_frequency(raw_file, target, raw_meta["state_count"])
            else:
                shutil.copyfile(
                    raw_file if name == "raw" else reconstructed_file, target
                )
                corrupt = None
            case = {
                "workdir": str(work),
                "input_sha256": sha(work / INPUT_NAME),
                "restart_file_sha256": sha(target),
                "negative_control": corrupt,
            }
            receipt["cases"][name] = case
            try:
                require(
                    case["input_sha256"] == receipt["input"]["modified_sha256"],
                    "restart input hash mismatch",
                )
                case["execution"] = execute(executable, work, timeout)
                check_execution(case["execution"], negative=name == "wrong-frequency")
                if name != "wrong-frequency":
                    result_path = work / output_name
                    require(
                        result_path.is_file() and result_path.stat().st_size > 0,
                        "missing native Padé table",
                    )
                    case["output"] = {
                        "path": str(result_path),
                        "sha256": sha(result_path),
                        "bytes": result_path.stat().st_size,
                    }
                case["status"] = "passed"
            except Exception as error:
                case.update(status="failed", error=f"{type(error).__name__}: {error}")
                receipt["errors"].append(f"{name}: {type(error).__name__}: {error}")
        if all(
            receipt["cases"][x]["status"] == "passed"
            for x in ("raw", "reconstructed", "wrong-frequency")
        ):
            try:
                comparison = compare_tables(
                    case_root / "raw" / output_name,
                    case_root / "reconstructed" / output_name,
                    raw_meta["state_energy_minus_ef0_ev"],
                    nfreq,
                    step,
                )
                receipt["comparison"] = comparison
                if not comparison["passed"]:
                    receipt["errors"].append(
                        "one or more native output channels exceeded fixed tolerance"
                    )
            except Exception as error:
                receipt["errors"].append(
                    f"table comparison: {type(error).__name__}: {error}"
                )
        if not receipt["errors"]:
            receipt["status"] = "passed"
    except Exception as error:
        receipt["errors"].append(f"{type(error).__name__}: {error}")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--raw-file", type=Path, required=True)
    parser.add_argument("--reconstructed-file", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--qe-root", type=Path, default=Path("/opt/qe"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native-timeout", type=int, default=1200)
    parser.add_argument("--temperature-kelvin", type=int, choices=(17, 20), default=17)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = run(
        args.workspace,
        args.raw_file,
        args.reconstructed_file,
        args.source_revision,
        args.qe_root,
        args.output,
        args.native_timeout,
        temperature_kelvin=args.temperature_kelvin,
    )
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
