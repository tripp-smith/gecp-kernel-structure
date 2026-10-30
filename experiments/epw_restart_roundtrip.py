"""Frozen 20 K model-plus-XOR experiment; no native solver acceptance is implied."""

from __future__ import annotations

import argparse
import gzip
import importlib.metadata
import os
import platform
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from experiments.epw_frequency import (
    BUDGET,
    ROOT,
    SOURCE_FILES,
    THREADS,
    fit,
    load,
    measure,
    measured_fit_timings,
    passed,
    sha,
    write,
)
from kernelgecp.integrations.epw import read_epw_anisotropic, write_epw_anisotropic
from kernelgecp.integrations.epw_restart_archive import EPWRestartArchive
from scripts.providers.epw_anisotropic_fixture import (
    SCOUT_INPUT_SHA256,
    UPSTREAM_INPUT_SHA256,
    UPSTREAM_IO_SHA256,
)

HOLDOUT = (45, 49, 53)
LIMIT = 330 * 60
CASE_LIMIT = 1200
SOURCES = tuple(
    dict.fromkeys(
        (
            *SOURCE_FILES,
            "experiments/epw_restart_roundtrip.py",
            "kernelgecp/integrations/epw_restart_archive.py",
            "scripts/providers/epw_restart.py",
            "scripts/providers/epw_smoke.py",
            "tests/test_epw_restart_roundtrip.py",
            "experiments/audit_epw_restart_roundtrip.py",
            "tests/test_epw_restart_roundtrip_audit.py",
        )
    )
)
CONTROLS = {
    "temperature_kelvin": 20,
    "rank": 60,
    "dictionary_size": 1201,
    "pole_cutoff_ev": 13,
    "precision_bits": 256,
    "rcond": 1e-12,
    "max_evaluation_condition": 1e4,
    "real_nambu_parity": True,
    "holdout_rows": list(HOLDOUT),
    "budget": BUDGET,
    "policy": "lossless_native_xor_v1",
}


def source_environment():
    if subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True
    ).strip():
        raise ValueError("tracked source must be clean")
    for path in SOURCES:
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", path],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    if any(os.environ.get(k) != "1" for k in THREADS):
        raise ValueError("single-thread environment required")
    versions = {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "mpmath")}
    if versions != {
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "mpmath": "1.3.0",
    } or sys.version_info[:2] != (3, 12):
        raise ValueError("scientific runtime versions mismatch")
    if (
        Path(sys.prefix).resolve()
        != (ROOT / "environments/epw-science/.venv").resolve()
    ):
        raise ValueError("isolated EPW runtime required")
    return {
        "revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_hashes": {p: sha(ROOT / p) for p in SOURCES},
        "threads": {k: os.environ[k] for k in THREADS},
        "versions": versions,
        "python": platform.python_version(),
        "prefix": str(Path(sys.prefix).resolve()),
        "executable": str(Path(sys.executable).resolve()),
        "platform": platform.platform(),
    }


def native_context(native_root, offset, temperature):
    native_root = Path(native_root).resolve()
    receipt_path = native_root / "epw-anisotropic.json"
    rec = load(receipt_path)
    execution = rec.get("native_anisotropic_execution", {})
    if rec.get("script_sha256") != sha(
        ROOT / "scripts/providers/epw_anisotropic_fixture.py"
    ):
        raise ValueError("producer source hash mismatch")
    if not all(rec.get("baseline", {}).get("checks", {}).values()) or not rec.get(
        "baseline", {}
    ).get("checks"):
        raise ValueError("native baseline checks missing or failed")
    if (
        rec.get("schema") != "kernelgecp.epw-native-anisotropic-fixture.v1"
        or rec.get("status") != "passed"
        or rec.get("provider") != {"qe": "7.5", "epw": "6.0"}
        or rec.get("wscut_ev") != 1.0
        or execution.get("returncode") != 0
        or rec.get("upstream_input_sha256") != UPSTREAM_INPUT_SHA256
        or rec.get("upstream_io_sha256") != UPSTREAM_IO_SHA256
        or rec.get("executed_input_sha256") != SCOUT_INPUT_SHA256
    ):
        raise ValueError("native receipt unsupported or failed")
    baseline_executables = [
        x.get("executable_sha256")
        for x in rec["baseline"].get("native_executions", [])
        if x.get("input") == "epw1.in"
    ]
    if (
        len(baseline_executables) != 1
        or execution.get("executable_sha256") != baseline_executables[0]
    ):
        raise ValueError("native executable identity mismatch")
    work = native_root / "epw-native/test-suite/epw_super"
    inp, log = work / "epw6.in", work / "epw6.out"
    if (
        sha(inp) != rec.get("executed_input_sha256")
        or sha(inp) != execution.get("input_sha256")
        or sha(log) != execution.get("log_sha256")
    ):
        raise ValueError("native input/log hash mismatch")
    text = log.read_text()
    number = r"([+-]?\d+\.\d+[Ee][+-]\d+)"
    fermi = re.findall(r"Fermi level \(eV\)\s*=\s*" + number, text)
    mu = re.findall(r"Chemical potential \(itemp\s*=\s*2\)\s*=\s*" + number, text)
    if (
        len(fermi) != 1
        or len(mu) != 1
        or len(re.findall("Convergence was reached in nsiter", text)) < 3
    ):
        raise ValueError(
            "missing/ambiguous native Fermi, chemical potential, or convergence"
        )
    if not np.isfinite(offset) or abs(float(mu[0]) - float(fermi[0]) - offset) > 1e-10:
        raise ValueError(
            "explicit chemical offset does not match printed native values"
        )
    files = rec.get("files", [])
    selected = [f for f in files if f.get("temperature_kelvin") == temperature]
    if (
        len(files) != 3
        or {f.get("temperature_kelvin") for f in files} != {17, 20, 30}
        or len(selected) != 1
    ):
        raise ValueError("native temperature receipt mismatch")
    for f in files:
        if (
            Path(f["filename"]).name != f["filename"]
            or f["filename"] != f"MgB2.imag_aniso_{f['temperature_kelvin']:06.2f}"
        ):
            raise ValueError("native filename mismatch")
        p = work / f["filename"]
        if not p.is_file() or sha(p) != f["sha256"] or p.stat().st_size != f["bytes"]:
            raise ValueError("native file hash mismatch")
    f = selected[0]
    path = work / f["filename"]
    if (
        f.get("columns")
        != [
            "frequency_ev",
            "energy_minus_ef0_ev",
            "Z",
            "gap_ev",
            "shift_ev",
        ]
        or not f.get("positive_only")
        or f.get("state_ordering") != "native_frequency_major_then_ik_then_ibnd_ordinal"
    ):
        raise ValueError("expected positive full-bandwidth table")
    context = {
        "schema": "kernelgecp.epw-anisotropic-input.v1",
        "provider": rec["provider"],
        "temperature_kelvin": temperature,
        "units": "eV",
        "state_order": "frequency_major_ik_ibnd_ordinal",
        "state_count": f["state_count"],
        "matsubara_indices": list(range(f["frequency_count"])),
        "positive_only": True,
        "full_bandwidth": True,
        "chi_convention": "native_fbw_column",
        "reference_fermi_energy_ev": float(fermi[0]),
        "chemical_potential_offset_ev": offset,
        "native_sha256": sha(path),
        "native_filename": path.name,
        "source_provenance": {
            "receipt_sha256": sha(receipt_path),
            "input_sha256": sha(inp),
            "log_sha256": sha(log),
            "executable_sha256": execution["executable_sha256"],
            "source_revision": rec["source_revision"],
            "source_input_sha256": rec["upstream_input_sha256"],
            "source_io_sha256": rec["upstream_io_sha256"],
            "printed_fermi_ev": fermi[0],
            "printed_mu_ev": mu[0],
            "corpus_sha256": {f["filename"]: f["sha256"] for f in files},
        },
    }
    return context, path, f


def read_checked(path, context, f):
    data = read_epw_anisotropic(path, context)
    if not np.array_equal(
        data.frequencies_ev, np.asarray(f["frequency_ev"])
    ) or not np.array_equal(
        data.state_energies_ev, np.asarray(f["state_energy_minus_ef0_ev"])
    ):
        raise ValueError("receipt frequency/state values mismatch")
    return data


def same_native(a, b):
    for name in ("frequencies_ev", "state_energies_ev", "z", "delta_ev", "chi_ev"):
        x, y = getattr(a, name), getattr(b, name)
        if x.shape != y.shape or not np.array_equal(
            np.asarray(x, dtype="<f8").view("<u8"),
            np.asarray(y, dtype="<f8").view("<u8"),
        ):
            raise ValueError("native bits mismatch: " + name)
    if not np.array_equal(a.indices, b.indices):
        raise ValueError("native frequency ordering mismatch")


def _case(native_root, offset, output, repeat, frozen):
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    record = {
        "schema": "kernelgecp.epw-restart-roundtrip-case.v1",
        "status": "running",
        "repeat": repeat,
        "freeze_sha256": sha(frozen),
        "native_acceptance": "not_tested",
        "corrected_values_are_fit_predictions": False,
    }
    try:
        freeze = load(frozen)
        if source_environment() != freeze["environment"]:
            raise ValueError("frozen source/runtime changed")
        context, path, table = native_context(native_root, offset, 20)
        if context != freeze["context"] or freeze["controls"] != CONTROLS:
            raise ValueError("frozen controls or native provenance changed")
        data = read_checked(path, context, table)
        n = len(data.indices)
        if n <= max(HOLDOUT):
            raise ValueError("insufficient frequency rows")
        training = [i for i in range(n) if i not in HOLDOUT]
        t = time.monotonic()
        model = fit(data, training, 60)
        fit_time = time.monotonic() - t
        t = time.monotonic()
        metrics, arrays = {}, {}
        for name, rows in (
            ("heldout", HOLDOUT),
            ("training", training),
            ("full", list(range(n))),
        ):
            metrics[name], vals = measure(data, model, rows)
            arrays.update({name + "_" + k: v for k, v in vals.items()})
        evaluation = time.monotonic() - t
        # Record the frozen pure-model result before any truth-assisted corrections.
        record.update(
            pure_model_metrics=metrics,
            pure_model_heldout_passed=passed(metrics["heldout"], 1),
            diagnostics=model.diagnostics,
            fit_rows=training,
        )
        write(output / "case.json", record)
        if not record["pure_model_heldout_passed"]:
            raise ValueError("pure-model heldout exceeds precommitted budgets")
        model.save(output / "model.json")
        t = time.monotonic()
        archive = EPWRestartArchive.from_data(data, model, policy=CONTROLS["policy"])
        archive.save(output / "bundle.json")
        encode = time.monotonic() - t
        t = time.monotonic()
        restored = EPWRestartArchive.load(output / "bundle.json")
        reconstructed = restored.reconstruct_native()
        same_native(data, reconstructed)
        decode = time.monotonic() - t
        t = time.monotonic()
        replay = output / "raw-replay.dat"
        write_epw_anisotropic(data, replay)
        if replay.read_bytes() != path.read_bytes():
            raise ValueError("source replay differs")
        restored_context = restored.write_native_restart(output / "restored-native.dat")
        same_native(
            data, read_epw_anisotropic(output / "restored-native.dat", restored_context)
        )
        export = time.monotonic() - t
        arrays.update(
            coefficients=model.coefficients,
            frequencies_ev=data.frequencies_ev,
            state_energies_ev=data.state_energies_ev,
            restored_z=reconstructed.z,
            restored_gap=reconstructed.delta_ev,
            restored_chi=reconstructed.chi_ev,
        )
        np.savez_compressed(output / "scientific-arrays.npz", **arrays)
        if (
            source_environment() != freeze["environment"]
            or native_context(native_root, offset, 20)[0] != context
        ):
            raise ValueError("source or native provenance changed during case")
        raw_binary = sum(
            getattr(data, k).nbytes
            for k in ("frequencies_ev", "state_energies_ev", "z", "delta_ev", "chi_ev")
        )
        record.update(
            status="passed",
            lossless_native_bits_equal=True,
            printed_native_fields_equal=True,
            bytes={
                "raw_ascii": path.stat().st_size,
                "gzip_ascii": len(gzip.compress(path.read_bytes(), mtime=0)),
                "raw_binary_arrays": raw_binary,
                "model_only": sum(
                    (output / x).stat().st_size for x in ("model.json", "model.npz")
                ),
                "model_plus_xor_bundle": sum(
                    (output / x).stat().st_size for x in ("bundle.json", "bundle.npz")
                ),
            },
            timing_seconds={
                "construction_and_fit": fit_time,
                **measured_fit_timings(model),
                "evaluation": evaluation,
                "encode": encode,
                "decode": decode,
                "export": export,
                "whole": time.monotonic() - start,
            },
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            * (1 if sys.platform == "darwin" else 1024),
            files_sha256={
                p.name: sha(p)
                for p in output.iterdir()
                if p.is_file() and p.name != "case.json"
            },
            scope=(
                "Finite 20 K pure-model diagnostics and truth-assisted all-field "
                "XOR reconstruction; no solver acceptance or continuous theorem."
            ),
        )
    except Exception as e:
        record.update(
            status="failed",
            error=f"{type(e).__name__}: {e}",
            elapsed_seconds=time.monotonic() - start,
        )
        write(output / "case.json", record)
        raise
    write(output / "case.json", record)
    return record


def campaign(native_root, temperature, offset, output):
    if temperature != 20:
        raise ValueError("only precommitted 20 K supported")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    record = {
        "schema": "kernelgecp.epw-restart-roundtrip-campaign.v1",
        "status": "running",
        "native_acceptance": "not_tested",
        "children": [],
    }
    try:
        environment = source_environment()
        # Hash and validate native metadata before reading or evaluating fields.
        context, _, _ = native_context(native_root, offset, 20)
        frozen = output / "freeze.json"
        write(
            frozen,
            {
                "schema": "kernelgecp.epw-restart-roundtrip-freeze.v1",
                "environment": environment,
                "controls": CONTROLS,
                "context": context,
                "native_root": str(Path(native_root).resolve()),
                "offset": offset,
            },
        )
        record["freeze_sha256"] = sha(frozen)
        write(output / "campaign.json", record)
        for repeat in range(3):
            remaining = LIMIT - (time.monotonic() - start)
            if remaining <= 0:
                raise TimeoutError("campaign deadline")
            case = output / f"case-{repeat}"
            log = output / f"case-{repeat}.log"
            cmd = [
                sys.executable,
                "-m",
                "experiments.epw_restart_roundtrip",
                "case",
                "--native-root",
                str(Path(native_root).resolve()),
                "--temperature-kelvin",
                "20",
                "--chemical-potential-offset-ev",
                repr(offset),
                "--output",
                str(case),
                "--freeze",
                str(frozen),
                "--repeat",
                str(repeat),
            ]
            child = {"repeat": repeat, "command": cmd}
            record["children"].append(child)
            write(output / "campaign.json", record)
            with log.open("xb") as stream:
                result = subprocess.run(
                    cmd,
                    cwd=ROOT,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=min(remaining, CASE_LIMIT),
                )
            child.update(exit_code=result.returncode, log_sha256=sha(log))
            receipt = case / "case.json"
            if receipt.exists():
                child["receipt_sha256"] = sha(receipt)
            if result.returncode or not receipt.exists():
                raise ValueError("child failed or omitted receipt")
            rec = load(receipt)
            if (
                rec["status"] != "passed"
                or rec["repeat"] != repeat
                or rec["freeze_sha256"] != sha(frozen)
            ):
                raise ValueError("child receipt mismatch")
            for name, digest in rec["files_sha256"].items():
                if Path(name).name != name or sha(case / name) != digest:
                    raise ValueError("child artifact mismatch")
        with np.load(
            output / "case-0/scientific-arrays.npz", allow_pickle=False
        ) as first:
            for i in (1, 2):
                with np.load(
                    output / f"case-{i}/scientific-arrays.npz", allow_pickle=False
                ) as other:
                    if set(first.files) != set(other.files) or any(
                        not np.array_equal(first[k], other[k]) for k in first.files
                    ):
                        raise ValueError("repeated numeric arrays differ")
        if (
            source_environment() != environment
            or native_context(native_root, offset, 20)[0] != context
        ):
            raise ValueError("source or native provenance changed during campaign")
        record.update(status="passed", repeat_arrays_identical=True)
    except Exception as e:
        record.update(
            status="inconclusive"
            if isinstance(e, (TimeoutError, subprocess.TimeoutExpired))
            else "failed",
            error=f"{type(e).__name__}: {e}",
        )
    record["process_wall_seconds"] = time.monotonic() - start
    write(output / "campaign.json", record)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("case", "campaign"))
    p.add_argument("--native-root", type=Path, required=True)
    p.add_argument("--temperature-kelvin", type=int, required=True)
    p.add_argument("--chemical-potential-offset-ev", type=float, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--freeze", type=Path)
    p.add_argument("--repeat", type=int, choices=(0, 1, 2))
    a = p.parse_args()
    try:
        if a.temperature_kelvin != 20:
            raise ValueError("only 20 K supported")
        if a.command == "case":
            if a.freeze is None or a.repeat is None:
                raise ValueError("case requires freeze and repeat")
            f = load(a.freeze)
            if (
                f["native_root"] != str(a.native_root.resolve())
                or f["offset"] != a.chemical_potential_offset_ev
            ):
                raise ValueError("freeze invocation mismatch")
            result = _case(
                a.native_root,
                a.chemical_potential_offset_ev,
                a.output,
                a.repeat,
                a.freeze,
            )
        else:
            if a.freeze is not None or a.repeat is not None:
                raise ValueError("campaign disallows freeze and repeat")
            result = campaign(
                a.native_root,
                a.temperature_kelvin,
                a.chemical_potential_offset_ev,
                a.output,
            )
        if result["status"] != "passed":
            raise ValueError(result.get("error", "campaign failed"))
    except Exception as e:
        print(f"{type(e).__name__}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
