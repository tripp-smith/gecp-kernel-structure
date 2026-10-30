"""Development-only native 17 K Padé-order calibration; never final acceptance."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

from scripts.providers import epw_restart as epw

PERCENTAGES = (90, 74, 58, 44, 30, 22, 15)
ORDERS = (98, 80, 64, 48, 32, 24, 16)
BINARY = "63409927dc9e36969d525b96ce4f8751d3294179fbaf1cb27b90248c4042d863"
RAW = "3bb9f9c0574b6001ff27c998a135c94e76c6ef56f2565b5b60ef34c3ea94b872"
FITTED = "6a7ad3ce942c01c1c9a1a7e9b2a4b54d1639c503928cbc9d4fc3bbda4e468258"
_deadline = math.inf


def remaining_timeout():
    remaining = _deadline - time.monotonic()
    if remaining < 1:
        raise TimeoutError("300-minute calibration budget exhausted")
    return min(1200, int(remaining)) if math.isfinite(remaining) else 1200


def checkpoint(out, report):
    path = out / "calibration-progress.json"
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    temporary.replace(path)


def order(percent, count=109):
    n = percent * count // 100
    n += n % 2
    return n - 2 if n > count else n


def extract(text):
    pieces = []
    for name in ("pade_coeff", "pade_eval"):
        found = re.findall(
            r"(?ims)^\s*SUBROUTINE\s+"
            + name
            + r"\s*\(.*?^\s*END SUBROUTINE\s+"
            + name
            + r"\s*$",
            text,
        )
        if len(found) != 1:
            raise ValueError("missing or ambiguous routine " + name)
        pieces.append(found[0])
    return "\n".join(pieces)


WRAPPER = """
module local_stdout
 integer, parameter :: stdout=6
end module
module extracted
contains
{routines}
end module
program check_pade
 use kinds
 use extracted
 use ep_constants, only: kelvin2ev, ryd2ev
 use constants, only: pi
 implicit none
 integer :: n,i,j,k,ios,st
 real(DP) :: v(109,140,5),kt,ws(9)
 complex(DP) :: z(109),u(109),a(109),y
 character(len=4096) :: path,arg,line,mode,gridpath
 call get_command_argument(1,arg)
 read(arg,*) n
 call get_command_argument(2,path)
 call get_command_argument(3,mode)
 if (trim(mode) == 'real') then
  call get_command_argument(4,gridpath)
  open(12,file=trim(gridpath),status='old',iostat=ios)
  if (ios /= 0) stop 4
  read(12,*,iostat=ios) ws
  if (ios /= 0) stop 5
  close(12)
 endif
 open(11,file=trim(path),status='old',iostat=ios)
 if (ios /= 0) stop 1
 read(11,'(a)',iostat=ios) line
 if (ios /= 0) stop 2
 do i=1,109
  do j=1,140
   read(11,'(5ES20.10)',iostat=ios) v(i,j,:)
   if (ios /= 0) stop 3
  enddo
 enddo
 close(11)
 call get_command_argument(5,arg)
 read(arg,*) kt
 kt=kt * kelvin2ev / ryd2ev
 kt=kt * ryd2ev
 do i=1,109
  z(i)=cmplx(0.d0,dble(2*i-1)*pi*kt,kind=DP)
 enddo
 do st=1,140
  do k=3,5
   u(:)=cmplx(v(:,st,k),0.d0,kind=DP)
   call pade_coeff(n,z(1:n),u(1:n),a(1:n))
   if (trim(mode) == 'real') then
    do i=1,9
     call pade_eval(n,z(1:n),a(1:n),cmplx(ws(i),0.d0,kind=DP),y)
     write(*,'(3i6,3ES26.17)') st,k,i,real(y,DP),aimag(y),0.d0
    enddo
   else
    do i=n+1,109
     call pade_eval(n,z(1:n),a(1:n),z(i),y)
     write(*,'(3i6,3ES26.17)') st,k,i,real(y,DP),aimag(y),v(i,st,k)
    enddo
   endif
  enddo
 enddo
end program
"""


def compile_helper(source, out, qe):
    original = extract(source.read_text())
    routines, count = re.subn(
        r"(?im)^(\s*)USE\s+io_global,\s*ONLY\s*:\s*stdout",
        r"\1USE local_stdout, ONLY : stdout",
        original,
    )
    if count != 1:
        raise ValueError("unexpected io_global import in extracted routines")
    src = out / "pade-extracted.f90"
    src.write_text(WRAPPER.format(routines=routines))
    makeinc = qe / "make.inc"
    query = out / "compiler-flags.mk"
    query.write_text(
        "include "
        + str(makeinc)
        + "\n.PHONY: diagnostic_print\ndiagnostic_print:\n"
        + '\t@printf "%s\\n" "$(MPIF90)" "$(F90FLAGS)"\n'
    )
    lines = subprocess.check_output(
        ["make", "-s", "-f", str(query), "diagnostic_print"],
        cwd=qe,
        text=True,
        timeout=remaining_timeout(),
    ).splitlines()
    if len(lines) != 2:
        raise ValueError("cannot resolve QE MPIF90/F90FLAGS")
    compiler = shlex.split(lines[0])
    flags = shlex.split(lines[1])
    if not compiler or not flags:
        raise ValueError("empty compiler or flags")
    modpaths = []
    for name in ("kinds.mod", "ep_constants.mod", "constants.mod"):
        found = list(qe.rglob(name))
        if not found:
            raise ValueError("missing actual QE module " + name)
        # Prefer compiled modules in expected source directories; reject ambiguity.
        expected = [
            x
            for x in found
            if x.parent in (qe / "EPW/src", qe / "Modules", qe / "UtilXlib")
        ]
        choices = expected or found
        if len({epw.sha(x) for x in choices}) != 1:
            raise ValueError("ambiguous QE module " + name)
        modpaths.append(sorted(choices)[0])
    binary = out / "pade-heldout"
    cmd = (
        compiler
        + flags
        + ["-I" + str(x.parent) for x in modpaths]
        + ["-o", str(binary), str(src)]
    )
    subprocess.run(
        cmd,
        cwd=out,
        check=True,
        capture_output=True,
        text=True,
        timeout=remaining_timeout(),
    )
    return {
        "upstream_source_sha256": epw.sha(source),
        "source_bodies_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "extracted_source_sha256": epw.sha(src),
        "binary_sha256": epw.sha(binary),
        "make_inc_sha256": epw.sha(makeinc),
        "module_sha256": {x.name: epw.sha(x) for x in modpaths},
        "command": cmd,
        "compiler": subprocess.check_output(
            compiler + ["--version"], text=True, timeout=remaining_timeout()
        ).splitlines()[0],
        "stdout_rebinding": (
            "local stdout unit 6; coefficients and evaluation bodies "
            "otherwise unchanged"
        ),
        "temperature_sequence": "temps=17*kelvin2eV/ryd2ev; gtemp=temps*ryd2ev",
    }


def heldout(binary, n, source, out, role):
    result = subprocess.run(
        [str(binary), str(n), str(source), "heldout", "unused", "17"],
        cwd=out,
        text=True,
        capture_output=True,
        timeout=remaining_timeout(),
        check=True,
    )
    file = out / f"heldout-{role}-{n}.txt"
    file.write_text(result.stdout)
    errors = [0.0, 0.0, 0.0]
    imaginary = [0.0, 0.0, 0.0]
    count = 0
    seen = set()
    for line in result.stdout.splitlines():
        values = line.split()
        if len(values) != 6:
            raise ValueError("invalid diagnostic line")
        st, ch, i = map(int, values[:3])
        real, imag, truth = map(float, values[3:])
        key = (st, ch, i)
        if key in seen or not (
            1 <= st <= 140
            and 3 <= ch <= 5
            and n < i <= 109
            and all(map(math.isfinite, (real, imag, truth)))
        ):
            raise ValueError("invalid diagnostic data")
        seen.add(key)
        errors[ch - 3] = max(errors[ch - 3], abs(real - truth))
        imaginary[ch - 3] = max(imaginary[ch - 3], abs(imag))
        count += 1
    if count != 140 * 3 * (109 - n):
        raise ValueError("incomplete heldout diagnostic")
    budgets = (1e-5, 1e-6, 1e-6)
    return {
        "max_absolute_error": errors,
        "max_imaginary_defect": imaginary,
        "rows": count,
        "passed": all(
            x <= b and y <= b
            for x, y, b in zip(errors, imaginary, budgets, strict=True)
        ),
        "transcript": str(file),
        "sha256": epw.sha(file),
    }


def compare(a, b, states, nfreq, step):
    return epw.compare_tables(a, b, states, nfreq, step)


def helper_native_check(binary, raw, native, out):
    indices = (0, 1, 9, 49, 99, 499, 999, 1999, 4499)
    expected = {}
    with native.open() as stream:
        next(stream)
        for st in range(1, 141):
            for i in range(4500):
                line = next(stream)
                if i in indices:
                    fields = [float(line[k : k + 20]) for k in range(0, 160, 20)]
                    for ch in range(3):
                        expected[(st, ch + 3, indices.index(i) + 1)] = (
                            fields[0],
                            fields[2 + ch * 2],
                            fields[3 + ch * 2],
                        )
    grid = out / "real-check-grid.txt"
    grid.write_text(" ".join(repr(expected[(1, 3, j)][0]) for j in range(1, 10)))
    proc = subprocess.run(
        [str(binary), "98", str(raw), "real", str(grid), "17"],
        cwd=out,
        text=True,
        capture_output=True,
        timeout=remaining_timeout(),
        check=True,
    )
    transcript = out / "helper-native-check.txt"
    transcript.write_text(proc.stdout)
    seen = set()
    worst = 0.0
    fail = 0
    for line in proc.stdout.splitlines():
        items = line.split()
        if len(items) != 6:
            raise ValueError("invalid helper native check")
        key = tuple(map(int, items[:3]))
        v = list(map(float, items[3:]))
        if key in seen or key not in expected or not all(map(math.isfinite, v)):
            raise ValueError("invalid helper native check data")
        seen.add(key)
        _, real, imag = expected[key]
        ratio = max(
            abs(v[0] - real) / (epw.ATOL + epw.RTOL * abs(real)),
            abs(v[1] - imag) / (epw.ATOL + epw.RTOL * abs(imag)),
        )
        worst = max(worst, ratio)
        fail += int(ratio > 1)
    if len(seen) != 140 * 3 * 9:
        raise ValueError("incomplete helper native check")
    return {
        "passed": fail == 0,
        "failures": fail,
        "max_budget_ratio": worst,
        "rows": len(seen),
        "transcript_sha256": epw.sha(transcript),
    }


def run(artifact, qe, out, source_revision):
    global _deadline
    started = time.monotonic()
    _deadline = started + 300 * 60
    out.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "kernelgecp.epw-pade-calibration.v1",
        "status": "inconclusive",
        "source_revision": source_revision,
        "driver_sha256": epw.sha(Path(__file__)),
        "restart_comparator_sha256": epw.sha(Path(epw.__file__)),
        "scope": "consumed 17 K development only; no final claim",
        "native_printed_convergence": "fixed N message, not convergence test",
        "percentages": list(PERCENTAGES),
        "orders": list(ORDERS),
        "candidates": [],
        "errors": [],
        "budget_seconds": 300 * 60,
    }
    try:
        if re.fullmatch(r"[0-9a-f]{40}", source_revision) is None:
            raise ValueError("full source revision required")
        binary = qe / "EPW/bin/epw.x"
        saved = json.loads((artifact / "restart-comparison.json").read_text())
        raw = artifact / "campaign/case-0/original-replay.dat"
        fit = artifact / "campaign/case-0/fitted-native.dat"
        if epw.sha(binary) != BINARY or saved["source"]["binary_sha256"] != BINARY:
            raise ValueError("native binary mismatch")
        if epw.sha(raw) != RAW or epw.sha(fit) != FITTED:
            raise ValueError("17 K input mismatch")
        if (
            saved["inputs"]["raw"]["sha256"] != RAW
            or saved["inputs"]["reconstructed"]["sha256"] != FITTED
        ):
            raise ValueError("saved receipt input mismatch")
        if tuple(map(order, PERCENTAGES)) != ORDERS:
            raise ValueError("native order mapping mismatch")
        report["helper"] = compile_helper(
            qe / "EPW/src/utilities/utilities.f90", out, qe
        )
        report["binary_sha256"] = BINARY
        report["inputs"] = {"raw": RAW, "reconstructed": FITTED}
        workspace = artifact / "epw-native/test-suite/epw_super"
        if epw.sha(workspace / epw.INPUT_NAME) != epw.SCOUT_INPUT_SHA256:
            raise ValueError("baseline input mismatch")
        modified = epw.restart_input((workspace / epw.INPUT_NAME).read_text())
        if re.search(r"(?im)^\s*npade\s*=", modified):
            raise ValueError("npade already present")
        old = artifact / "restart-comparison-restart-cases/raw" / epw.OUTPUT_NAME
        if epw.sha(old) != saved["cases"]["raw"]["output"]["sha256"]:
            raise ValueError("saved N98 output mismatch")
        states = saved["inputs"]["raw"]["state_energy_minus_ef0_ev"]
        grid = saved["expected_real_grid"]

        def compare_grid(a, b):
            return compare(
                a,
                b,
                states,
                grid["frequency_count"],
                grid["step_ev_derived_from_rounded_a2f"],
            )

        for percent, n in zip(PERCENTAGES, ORDERS, strict=True):
            remaining_timeout()
            row = {"percent": percent, "order": n, "cases": {}, "status": "failed"}
            report["candidates"].append(row)
            checkpoint(out, report)
            for role, source in (("raw", raw), ("reconstructed", fit)):
                work = out / f"p{percent}-{role}"

                def ignore(_path, names):
                    return [
                        x
                        for x in names
                        if "020.00" in x
                        or "030.00" in x
                        or x
                        in (
                            epw.OUTPUT_NAME,
                            "MgB2.pade_aniso_gap0_017.00",
                            "restart.out",
                        )
                    ]

                shutil.copytree(workspace, work, ignore=ignore)
                (work / epw.INPUT_NAME).write_text(
                    modified.replace(
                        "  imag_read = .true.",
                        f"  imag_read = .true.\n  npade = {percent}",
                    )
                )
                shutil.copyfile(source, work / epw.RAW_NAME)
                case = {"input_sha256": epw.sha(work / epw.INPUT_NAME)}
                row["cases"][role] = case
                try:
                    execution = epw.execute(binary, work, remaining_timeout())
                    case["execution"] = execution
                    epw.check_execution(execution, negative=False)
                    if not re.search(
                        r"Convergence was reached for N =\s*" + str(n) + r"\b",
                        (work / "restart.out").read_text(errors="replace"),
                    ):
                        raise ValueError("native N mismatch")
                    case.update(
                        status="passed", output_sha256=epw.sha(work / epw.OUTPUT_NAME)
                    )
                except Exception as error:
                    case.update(
                        status="failed", error=f"{type(error).__name__}: {error}"
                    )
                checkpoint(out, report)
            if any(x["status"] != "passed" for x in row["cases"].values()):
                if percent == 90:
                    raise ValueError("N98 execution failed")
                continue
            try:
                a = out / f"p{percent}-raw" / epw.OUTPUT_NAME
                b = out / f"p{percent}-reconstructed" / epw.OUTPUT_NAME
                if percent == 90:
                    report["baseline_comparison"] = compare_grid(old, a)
                    if not report["baseline_comparison"]["passed"]:
                        raise ValueError("N98 native baseline mismatch")
                    report["helper_native_check"] = helper_native_check(
                        out / "pade-heldout", raw, old, out
                    )
                    if not report["helper_native_check"]["passed"]:
                        raise ValueError("N98 extracted helper mismatch")
                row["comparison"] = compare_grid(a, b)
                row["heldout_raw"] = heldout(out / "pade-heldout", n, raw, out, "raw")
                row["heldout_reconstructed"] = heldout(
                    out / "pade-heldout", n, fit, out, "reconstructed"
                )
                if (
                    row["comparison"]["passed"]
                    and row["heldout_raw"]["passed"]
                    and row["heldout_reconstructed"]["passed"]
                ):
                    row["status"] = "passed"
            except Exception as error:
                row["error"] = f"{type(error).__name__}: {error}"
                if percent == 90:
                    raise
            checkpoint(out, report)
        report["passing_percentages"] = [
            r["percent"] for r in report["candidates"] if r["status"] == "passed"
        ]
        report["status"] = (
            "passed_calibration" if report["passing_percentages"] else "inconclusive"
        )
        if report["passing_percentages"]:
            report["selected_percent"] = report["passing_percentages"][0]
            report["selection_policy"] = (
                "first passing candidate in predeclared descending order; no final data"
            )
            freeze = {
                "schema": "kernelgecp.epw-pade-calibration-freeze.v1",
                "status": "passed_calibration",
                "percent": report["selected_percent"],
                "order": order(report["selected_percent"]),
                "binary_sha256": BINARY,
                "raw_sha256": RAW,
                "fitted_sha256": FITTED,
                "helper": report["helper"],
                "driver_sha256": epw.sha(Path(__file__)),
                "selection_policy": report["selection_policy"],
            }
            target = out / "calibration-freeze.json"
            target.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
            report["freeze_sha256"] = epw.sha(target)
    except Exception as error:
        report["errors"].append(f"{type(error).__name__}: {error}")
    report["elapsed_seconds"] = time.monotonic() - started
    checkpoint(out, report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--artifact", type=Path, required=True)
    p.add_argument("--qe-root", type=Path, default=Path("/opt/qe"))
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--source-revision", required=True)
    args = p.parse_args()
    report = run(
        args.artifact.resolve(),
        args.qe_root.resolve(),
        args.output.resolve(),
        args.source_revision,
    )
    (args.output / "calibration.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return 0 if report["status"] == "passed_calibration" else 1


if __name__ == "__main__":
    sys.exit(main())
