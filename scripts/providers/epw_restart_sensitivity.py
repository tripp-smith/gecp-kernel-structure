"""Truth-assisted development-only sensitivity for the consumed 17 K restart."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import time
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import math

from scripts.providers import epw_pade_calibration as calibration
from scripts.providers import epw_restart as epw

SCALES = (1.0, 0.1, 0.01, 0.001, 0.0)
BUDGET = 300 * 60


def remaining(deadline):
    seconds = deadline - time.monotonic()
    if seconds < 1:
        raise TimeoutError("300-minute sensitivity budget exhausted")
    return min(1200, int(seconds))


def checkpoint(out, report):
    temp = out / "sensitivity-progress.tmp"
    temp.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temp.replace(out / "sensitivity-progress.json")


FIELD = re.compile(r" *[+-]?\d\.\d{10}[EeDd][+-]\d{2}")


def prepare(raw, fitted, context, target, scale):
    """Format only the three fields and measure the resulting rounded errors."""
    if scale not in SCALES:
        raise ValueError("scale is not predeclared")
    raw_meta = epw.inspect_native_file(raw)
    fit_meta = epw.inspect_native_file(fitted)
    for key in (
        "frequency_count",
        "state_count",
        "frequency_ev",
        "state_energy_minus_ef0_ev",
        "columns",
    ):
        epw.require(raw_meta[key] == fit_meta[key], "input axes mismatch")
    epw.require(
        context.get("provider") == {"qe": "7.5", "epw": "6.0"}
        and context.get("temperature_kelvin") == 17
        and context.get("state_count") == raw_meta["state_count"]
        and context.get("matsubara_indices") == list(range(raw_meta["frequency_count"]))
        and context.get("positive_only") is True
        and context.get("full_bandwidth") is True
        and context.get("chi_convention") == "native_fbw_column"
        and context.get("state_order") == "frequency_major_ik_ibnd_ordinal"
        and context.get("units") == "eV"
        and raw_meta["positive_only"],
        "context mismatch",
    )
    source = raw.read_bytes()
    raw_lines = source.decode("ascii").splitlines(keepends=True)
    fit_lines = fitted.read_bytes().decode("ascii").splitlines(keepends=True)
    epw.require(
        len(raw_lines) == len(fit_lines) and raw_lines[0] == fit_lines[0],
        "input header or row count mismatch",
    )
    output = [raw_lines[0]]
    metrics = {
        name: {
            "max_absolute_error": 0.0,
            "equal_value_count": 0,
            "equal_text_count": 0,
            "field_count": 0,
        }
        for name in ("Z", "gap_ev", "shift_ev")
    }
    for raw_line, fit_line in zip(raw_lines[1:], fit_lines[1:], strict=True):
        left, right = raw_line.rstrip("\r\n"), fit_line.rstrip("\r\n")
        epw.require(len(left) == len(right) == 100, "malformed native row")
        epw.require(left[:40] == right[:40], "input axes text mismatch")
        raw_fields = [left[i : i + 20] for i in range(0, 100, 20)]
        fit_fields = [right[i : i + 20] for i in range(0, 100, 20)]
        epw.require(
            all(FIELD.fullmatch(f) for f in raw_fields + fit_fields),
            "malformed ES20.10 field",
        )
        row = left[:40]
        for col, name in enumerate(("Z", "gap_ev", "shift_ev"), start=2):
            original = float(raw_fields[col].replace("D", "E").replace("d", "e"))
            fitted_value = float(fit_fields[col].replace("D", "E").replace("d", "e"))
            # Explicit endpoints avoid cancellation changing a rounded fitted value.
            value = (
                original
                if scale == 0
                else fitted_value
                if scale == 1
                else original + scale * (fitted_value - original)
            )
            epw.require(math.isfinite(value), "nonfinite perturbed field")
            formatted = raw_fields[col] if scale == 0 else f"{value:20.10E}"
            epw.require(
                len(formatted) == 20 and FIELD.fullmatch(formatted),
                "perturbed value is not representable in ES20.10",
            )
            actual = float(formatted.replace("D", "E").replace("d", "e"))
            item = metrics[name]
            item["max_absolute_error"] = max(
                item["max_absolute_error"], abs(actual - original)
            )
            item["equal_value_count"] += int(actual == original)
            item["equal_text_count"] += int(formatted == raw_fields[col])
            item["field_count"] += 1
            row += formatted
        output.append(row + raw_line[len(left) :])
    payload = source if scale == 0 else "".join(output).encode("ascii")
    target.write_bytes(payload)
    epw.inspect_native_file(target)
    return {
        "sha256": epw.sha(target),
        "frequency_and_state_text_unchanged": True,
        "fields": metrics,
    }


def run(artifact, qe, out, revision, *, budget_seconds=BUDGET):
    started = time.monotonic()
    deadline = started + budget_seconds
    out.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "kernelgecp.epw-restart-sensitivity.v1",
        "status": "inconclusive",
        "source_revision": revision,
        "scope": (
            "Truth-assisted consumed 17 K development feasibility only; "
            "no fitted-restart acceptance or new temperature validation."
        ),
        "provider": {"qe": "7.5", "epw": "6.0"},
        "npade_percent": 90,
        "native_order": 98,
        "scales": list(SCALES),
        "budget_seconds": budget_seconds,
        "cases": {
            "raw": {"status": "not_run"},
            **{f"scale-{x:g}": {"status": "not_run", "scale": x} for x in SCALES},
        },
        "errors": [],
        "driver_sha256": epw.sha(Path(__file__)),
        "comparator_sha256": epw.sha(Path(epw.__file__)),
    }
    checkpoint(out, report)
    try:
        epw.require(
            0 < budget_seconds <= BUDGET, "budget exceeds 300 minutes or is invalid"
        )
        epw.require(
            re.fullmatch(r"[0-9a-f]{40}", revision) is not None,
            "full source revision required",
        )
        saved = json.loads((artifact / "restart-comparison.json").read_text())
        epw.require(saved["provider"] == report["provider"], "saved provider mismatch")
        context = json.loads((artifact / "campaign/case-0/case.json").read_text())[
            "fitted_context"
        ]
        raw, fit = (
            artifact / "campaign/case-0" / x
            for x in ("original-replay.dat", "fitted-native.dat")
        )
        binary, official, io = (
            qe / "EPW/bin/epw.x",
            qe / "test-suite/epw_super/epw6.in",
            qe / "EPW/src/io/io_supercond.f90",
        )
        workspace = artifact / "epw-native/test-suite/epw_super"
        epw.require(
            epw.sha(binary) == calibration.BINARY == saved["source"]["binary_sha256"],
            "binary mismatch",
        )
        epw.require(
            epw.sha(raw) == calibration.RAW == saved["inputs"]["raw"]["sha256"]
            and epw.sha(fit)
            == calibration.FITTED
            == saved["inputs"]["reconstructed"]["sha256"],
            "input hash mismatch",
        )
        epw.require(
            epw.sha(official) == epw.UPSTREAM_INPUT_SHA256
            and epw.sha(io) == epw.UPSTREAM_IO_SHA256,
            "pinned source mismatch",
        )
        original = workspace / epw.INPUT_NAME
        epw.require(
            epw.sha(original) == epw.SCOUT_INPUT_SHA256, "workspace input mismatch"
        )
        epw.require(
            context["temperature_kelvin"] == 17
            and context["provider"] == report["provider"],
            "context mismatch",
        )
        modified = epw.restart_input(original.read_text())
        epw.require(
            not re.search(r"(?im)^\s*npade\s*=", modified), "npade already present"
        )
        old = artifact / "restart-comparison-restart-cases/raw" / epw.OUTPUT_NAME
        epw.require(
            epw.sha(old) == saved["cases"]["raw"]["output"]["sha256"],
            "saved baseline mismatch",
        )
        meta, fitted_meta = epw.inspect_native_file(raw), epw.inspect_native_file(fit)
        for key in (
            "frequency_count",
            "state_count",
            "frequency_ev",
            "state_energy_minus_ef0_ev",
            "columns",
        ):
            epw.require(meta[key] == fitted_meta[key], "input axes mismatch")
        nfreq, step = epw.expected_mesh(workspace / "MgB2.a2f")
        epw.require(
            saved["expected_real_grid"]["a2f_sha256"]
            == epw.sha(workspace / "MgB2.a2f"),
            "mesh mismatch",
        )
        report["source"] = {
            "binary_sha256": epw.sha(binary),
            "io_sha256": epw.sha(io),
            "official_input_sha256": epw.sha(official),
            "workspace_input_sha256": epw.sha(original),
            "raw_sha256": epw.sha(raw),
            "fitted_sha256": epw.sha(fit),
            "saved_raw_output_sha256": epw.sha(old),
            "a2f_sha256": epw.sha(workspace / "MgB2.a2f"),
            "context_sha256": hashlib.sha256(
                json.dumps(context, sort_keys=True).encode()
            ).hexdigest(),
            "saved_receipt_sha256": epw.sha(artifact / "restart-comparison.json"),
            "case_metadata_sha256": epw.sha(artifact / "campaign/case-0/case.json"),
        }
        report["tolerances"] = {
            "atol": epw.ATOL,
            "rtol": epw.RTOL,
            "reference": "absolute raw output",
        }
        checkpoint(out, report)

        def ignore(_path, names):
            return [
                x
                for x in names
                if "020.00" in x
                or "030.00" in x
                or x in (epw.OUTPUT_NAME, "MgB2.pade_aniso_gap0_017.00", "restart.out")
            ]

        def execute_case(name, scale=None):
            remaining(deadline)
            case = {"status": "failed", "scale": scale}
            report["cases"][name] = case
            checkpoint(out, report)
            try:
                work = out / name
                shutil.copytree(workspace, work, ignore=ignore)
                (work / epw.INPUT_NAME).write_text(modified)
                target = work / epw.RAW_NAME
                if scale is None:
                    shutil.copyfile(raw, target)
                    epw.require(epw.sha(target) == calibration.RAW, "raw copy mismatch")
                else:
                    case["perturbation"] = prepare(raw, fit, context, target, scale)
                case["input_sha256"] = epw.sha(work / epw.INPUT_NAME)
                epw.require(
                    case["input_sha256"]
                    == hashlib.sha256(modified.encode()).hexdigest(),
                    "input copy mismatch",
                )
                case["execution"] = epw.execute(binary, work, remaining(deadline))
                epw.check_execution(case["execution"], negative=False)
                epw.require(
                    re.search(
                        r"Convergence was reached for N =\s*98\b",
                        (work / "restart.out").read_text(errors="replace"),
                    )
                    is not None,
                    "native N mismatch",
                )
                result = work / epw.OUTPUT_NAME
                epw.require(
                    result.is_file() and result.stat().st_size > 0,
                    "missing native output",
                )
                case.update(output_sha256=epw.sha(result), status="executed")
            except Exception as error:
                case["error"] = f"{type(error).__name__}: {error}"
                report["errors"].append(f"{name}: {case['error']}")
            checkpoint(out, report)

        execute_case("raw")
        epw.require(
            report["cases"]["raw"]["status"] == "executed", "raw baseline failed"
        )
        baseline = out / "raw" / epw.OUTPUT_NAME
        report["baseline_comparison"] = epw.compare_tables(
            old, baseline, meta["state_energy_minus_ef0_ev"], nfreq, step
        )
        epw.require(
            report["baseline_comparison"]["passed"], "saved baseline comparison failed"
        )
        for scale in SCALES:
            remaining(deadline)
            name = f"scale-{scale:g}"
            execute_case(name, scale)
            case = report["cases"][name]
            if case["status"] == "executed":
                try:
                    case["comparison"] = epw.compare_tables(
                        baseline,
                        out / name / epw.OUTPUT_NAME,
                        meta["state_energy_minus_ef0_ev"],
                        nfreq,
                        step,
                    )
                    case["status"] = (
                        "passed"
                        if case["comparison"]["passed"]
                        else "comparison_failed"
                    )
                except Exception as error:
                    case.update(
                        status="failed", error=f"{type(error).__name__}: {error}"
                    )
                    report["errors"].append(f"{name}: {case['error']}")
            checkpoint(out, report)
        report["status"] = (
            "completed"
            if not report["errors"] and report["cases"]["scale-0"]["status"] == "passed"
            else "inconclusive"
        )
    except Exception as error:
        report["errors"].append(f"{type(error).__name__}: {error}")
    report["elapsed_seconds"] = time.monotonic() - started
    checkpoint(out, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--qe-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    return (
        0
        if run(args.artifact, args.qe_root, args.output, args.source_revision)["status"]
        == "completed"
        else 1
    )


if __name__ == "__main__":
    sys.exit(main())
