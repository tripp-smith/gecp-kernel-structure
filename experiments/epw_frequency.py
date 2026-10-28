"""Frozen EPW real-parity frequency-cross experiment; no EPW solver execution."""

from __future__ import annotations

import argparse
import dataclasses
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from kernelgecp.integrations.epw import (
    nambu_green,
    read_epw_anisotropic,
    write_epw_anisotropic,
)
from kernelgecp.integrations.epw_frequency import (
    EPWFrequencyCrossFit,
    fit_epw_frequency_cross,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = (
    "experiments/epw_frequency.py",
    "experiments/audit_epw_frequency.py",
    "kernelgecp/integrations/epw.py",
    "kernelgecp/integrations/epw_frequency.py",
    "kernelgecp/integrations/_frequency_cross.py",
    "kernelgecp/_basis_io.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "kernelgecp/types.py",
    "kernelgecp/__init__.py",
    "kernelgecp/integrations/__init__.py",
    "kernelgecp/_basis_numerics.py",
    "kernelgecp/_reference_backend.py",
    "kernelgecp/gecp.py",
    "kernelgecp/kernels.py",
    "scripts/providers/epw_anisotropic_fixture.py",
    "environments/epw-science/pyproject.toml",
    "environments/epw-science/uv.lock",
)
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
CAL = (36, 40, 44)
FINAL = (45, 49, 53)
RANKS = (48, 60, 72)
BUDGET = {"green_scaled": 1e-6, "z": 1e-5, "gap_ev": 1e-6, "chi_ev": 1e-6}
LIMIT = 330 * 60


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    os.replace(tmp, path)


def load(path):
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict) or not data:
        raise ValueError("invalid metadata")
    return data


def source_environment(*, pilot=False):
    files = {p: sha(ROOT / p) for p in SOURCE_FILES}
    clean = not subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True
    ).strip()
    # Untracked source cannot be bound to a clean revision.
    for p in SOURCE_FILES:
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", p],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
    if not pilot and not clean:
        raise ValueError("tracked source must be clean")
    if any(os.environ.get(k) != "1" for k in THREADS):
        raise ValueError("single-thread environment required")
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "pytest")
    }
    required = {
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "mpmath": "1.3.0",
        "pytest": "9.0.3",
    }
    if versions != required or sys.version_info[:2] != (3, 12):
        raise ValueError("scientific runtime versions mismatch")
    prefix = (ROOT / "environments/epw-science/.venv").resolve()
    if not pilot and Path(sys.prefix).resolve() != prefix:
        raise ValueError("isolated EPW runtime required")
    return {
        "revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_hashes": files,
        "tracked_clean": clean,
        "python": platform.python_version(),
        "python_executable": str(Path(sys.executable).resolve()),
        "python_prefix": str(Path(sys.prefix).resolve()),
        "versions": versions,
        "platform": platform.platform(),
        "threads": {k: os.environ[k] for k in THREADS},
        "development_lock_sha256_context_only": sha(ROOT / "uv.lock"),
    }


def stable(env):
    return {
        k: v
        for k, v in env.items()
        if k
        not in ("revision", "tracked_clean", "development_lock_sha256_context_only")
    }


def native_context(native_root, offset):
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
    mu = re.findall(r"Chemical potential \(itemp\s*=\s*1\)\s*=\s*" + number, text)
    if (
        len(fermi) != 1
        or len(mu) != 1
        or len(re.findall("Convergence was reached in nsiter", text)) < 3
    ):
        raise ValueError(
            "missing/ambiguous native Fermi, chemical potential, or convergence"
        )
    if (
        not np.isfinite(offset)
        or offset != 0
        or abs(float(mu[0]) - float(fermi[0]) - offset) > 1e-10
    ):
        raise ValueError(
            "explicit zero chemical offset does not match printed native values"
        )
    files = rec.get("files", [])
    selected = [f for f in files if f.get("temperature_kelvin") == 17]
    if (
        len(files) != 3
        or {f.get("temperature_kelvin") for f in files} != {17, 20, 30}
        or len(selected) != 1
    ):
        raise ValueError("native temperature receipt mismatch")
    for f in files:
        p = work / f["filename"]
        if not p.is_file() or sha(p) != f["sha256"] or p.stat().st_size != f["bytes"]:
            raise ValueError("native file hash mismatch")
    f = selected[0]
    path = work / f["filename"]
    if f.get("columns") != [
        "frequency_ev",
        "energy_minus_ef0_ev",
        "Z",
        "gap_ev",
        "shift_ev",
    ] or not f.get("positive_only"):
        raise ValueError("expected positive full-bandwidth table")
    context = {
        "schema": "kernelgecp.epw-anisotropic-input.v1",
        "provider": rec["provider"],
        "temperature_kelvin": 17,
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
    data = read_epw_anisotropic(path, context)
    if not np.array_equal(
        data.frequencies_ev, np.asarray(f["frequency_ev"])
    ) or not np.array_equal(
        data.state_energies_ev, np.asarray(f["state_energy_minus_ef0_ev"])
    ):
        raise ValueError("receipt frequency/state values mismatch")
    return data, context, path


def partition(n):
    if n <= max(FINAL):
        raise ValueError("not enough native frequencies")
    return {
        "train": [i for i in range(n) if i not in CAL and i not in FINAL],
        "calibration": list(CAL),
        "final": list(FINAL),
    }


def fit(data, rows, rank):
    return fit_epw_frequency_cross(
        data,
        real_nambu_parity=True,
        fit_rows=rows,
        rank=rank,
        pole_cutoff_ev=13,
        dictionary_size=1201,
        precision_bits=256,
        rcond=1e-12,
        max_evaluation_condition=1e4,
    )


def measure(data, model, rows):
    rows = np.asarray(rows, dtype=int)
    truth = nambu_green(data, rows=rows)
    pred = model.evaluate_matsubara(data.indices[rows])
    scale = np.max(np.abs(truth), axis=(-2, -1))
    err = np.abs(pred - truth)
    if np.any(scale == 0):
        raise ValueError("zero matrix scale")
    # Independent direct solve recovers fields using the printed native frequency.
    inv = np.linalg.solve(pred, np.broadcast_to(np.eye(2), pred.shape))
    z = (inv[..., 0, 0] + inv[..., 1, 1]) / (2j * data.frequencies_ev[rows, None])
    chi = (
        (inv[..., 1, 1] - inv[..., 0, 0]) / 2
        - data.state_energies_ev[None, :]
        + data.context["chemical_potential_offset_ev"]
    )
    gap = -(inv[..., 0, 1] + inv[..., 1, 0]) / (2 * z)
    imaginary = max(float(np.max(np.abs(x.imag))) for x in (z, chi, gap))
    if not np.isfinite(imaginary) or imaginary > 1e-8:
        raise ValueError("non-real recovered field")
    metrics = {
        "green_scaled": float(np.max(err / scale[..., None, None])),
        "green_global_scaled": float(np.max(err) / np.max(np.abs(truth))),
        "z": float(np.max(np.abs(z.real - data.z[rows]))),
        "gap_ev": float(np.max(np.abs(gap.real - data.delta_ev[rows]))),
        "chi_ev": float(np.max(np.abs(chi.real - data.chi_ev[rows]))),
        "imaginary_defect": imaginary,
    }
    return metrics, {
        "truth_green": truth,
        "predicted_green": pred,
        "truth_z": data.z[rows],
        "truth_gap": data.delta_ev[rows],
        "truth_chi": data.chi_ev[rows],
        "predicted_z": z.real,
        "predicted_gap": gap.real,
        "predicted_chi": chi.real,
    }


def passed(metrics, factor):
    return (
        all(
            np.isfinite(metrics[k]) and metrics[k] <= factor * v
            for k, v in BUDGET.items()
        )
        and metrics["imaginary_defect"] <= 1e-8
    )


def measured_fit_timings(model):
    times = model.diagnostics["timing_seconds"]
    keys = (
        "basis_construction",
        "nambu_observation_construction",
        "least_squares_preparation",
        "coefficient_fit",
        "training_replay",
    )
    if set(times) != set(keys) or any(
        not np.isfinite(times[k]) or times[k] < 0 for k in keys
    ):
        raise ValueError("missing or invalid measured fit timings")
    return {
        **times,
        "fit_total": times["least_squares_preparation"] + times["coefficient_fit"],
    }


def calibrate(native_root, output, freeze, offset, pilot=False):
    output = Path(output)
    freeze = Path(freeze)
    if output.exists() or freeze.exists():
        raise FileExistsError("calibration output exists")
    output.mkdir(parents=True)
    started = time.monotonic()
    env = source_environment(pilot=pilot)
    data, context, path = native_context(native_root, offset)
    parts = partition(len(data.indices))
    attempts = []
    selected = None
    report = {
        "schema": "kernelgecp.epw-frequency-calibration.v1",
        "status": "running",
        "pilot_noncanonical": pilot,
        "provenance": env,
        "context": context,
        "partitions": parts,
        "attempts": attempts,
    }
    try:
        for rank in RANKS:
            if time.monotonic() - started >= LIMIT:
                raise TimeoutError("resource deadline")
            attempt = {"rank": rank}
            attempts.append(attempt)
            try:
                t = time.monotonic()
                model = fit(data, parts["train"], rank)
                attempt["construction_and_fit_seconds"] = time.monotonic() - t
                attempt["timing_seconds"] = measured_fit_timings(model)
                attempt["diagnostics"] = model.diagnostics
                attempt["training"], training = measure(data, model, parts["train"])
                attempt["calibration"], calibration = measure(
                    data, model, parts["calibration"]
                )
                arrays = {
                    **{"training_" + k: v for k, v in training.items()},
                    **{"calibration_" + k: v for k, v in calibration.items()},
                }
                np.savez_compressed(output / f"rank-{rank}.npz", **arrays)
                attempt["arrays_sha256"] = sha(output / f"rank-{rank}.npz")
                attempt["passed"] = passed(attempt["training"], 0.1) and passed(
                    attempt["calibration"], 0.1
                )
                if attempt["passed"]:
                    selected = rank
                    break
            except Exception as error:
                attempt.update(passed=False, error=f"{type(error).__name__}: {error}")
            finally:
                write(output / "calibration.json", report)
        if selected is None:
            raise ValueError("no rank passed predeclared calibration margins")
        if source_environment(pilot=pilot) != env:
            raise ValueError("source or runtime changed")
        report.update(
            status="passed",
            selected_rank=selected,
            process_wall_seconds=time.monotonic() - started,
        )
        write(output / "calibration.json", report)
        frozen = {
            "schema": "kernelgecp.epw-frequency-freeze.v1",
            "status": "passed",
            "pilot_noncanonical": pilot,
            "provenance": env,
            "context": context,
            "partitions": parts,
            "native_root": str(Path(native_root).resolve()),
            "native_sha256": sha(path),
            "calibration_path": str((output / "calibration.json").resolve()),
            "calibration_sha256": sha(output / "calibration.json"),
            "selected_rank": selected,
            "rank_ladder": list(RANKS),
            "budget": BUDGET,
            "calibration_margin": 0.1,
            "controls": {
                "pole_cutoff_ev": 13,
                "dictionary_size": 1201,
                "precision_bits": 256,
                "rcond": 1e-12,
                "max_evaluation_condition": 1e4,
            },
            "scope": (
                "final rows neither constructed nor scored during calibration; "
                "no experimental prediction or speedup claim"
            ),
        }
        write(freeze, frozen)
    except Exception as error:
        report.update(
            status="inconclusive" if isinstance(error, TimeoutError) else "failed",
            error=f"{type(error).__name__}: {error}",
        )
        write(output / "calibration.json", report)
        raise
    return report


def checked_freeze(freeze_path, native_root, offset):
    frozen = load(freeze_path)
    env = source_environment()
    if (
        frozen.get("schema") != "kernelgecp.epw-frequency-freeze.v1"
        or frozen.get("status") != "passed"
        or frozen.get("pilot_noncanonical") is not False
        or stable(frozen.get("provenance", {})) != stable(env)
    ):
        raise ValueError("freeze source/runtime mismatch")
    if sha(frozen["calibration_path"]) != frozen.get("calibration_sha256"):
        raise ValueError("calibration receipt hash mismatch")
    calibration = load(frozen["calibration_path"])
    if (
        calibration.get("status") != "passed"
        or calibration.get("selected_rank") != frozen.get("selected_rank")
        or calibration.get("context") != frozen.get("context")
        or calibration.get("partitions") != frozen.get("partitions")
        or stable(calibration.get("provenance", {})) != stable(env)
    ):
        raise ValueError("calibration receipt disagrees with freeze")
    data, context, path = native_context(native_root, offset)
    if (
        frozen.get("context") != context
        or frozen.get("native_sha256") != sha(path)
        or frozen.get("partitions") != partition(len(data.indices))
        or frozen.get("rank_ladder") != list(RANKS)
        or frozen.get("budget") != BUDGET
    ):
        raise ValueError("freeze corpus or controls mismatch")
    if frozen.get("selected_rank") not in RANKS or frozen.get("controls") != {
        "pole_cutoff_ev": 13,
        "dictionary_size": 1201,
        "precision_bits": 256,
        "rcond": 1e-12,
        "max_evaluation_condition": 1e4,
    }:
        raise ValueError("unsupported frozen selection")
    return frozen, env, data, path


def run_case(native_root, freeze_path, offset, repeat, output):
    if repeat not in (0, 1, 2):
        raise ValueError("repeat must be 0, 1 or 2")
    output = Path(output)
    if output.exists():
        raise FileExistsError("case output exists")
    output.mkdir(parents=True)
    started = time.monotonic()
    record = {
        "schema": "kernelgecp.epw-frequency-case.v1",
        "status": "running",
        "repeat": repeat,
    }
    try:
        record["freeze_sha256"] = sha(freeze_path)
        frozen, env, data, path = checked_freeze(freeze_path, native_root, offset)
        parts = frozen["partitions"]
        record.update(provenance=env, partitions=parts, rank=frozen["selected_rank"])
        t = time.monotonic()
        model = fit(data, parts["train"], frozen["selected_rank"])
        fit_seconds = time.monotonic() - t
        t = time.monotonic()
        arrays = {}
        metrics = {}
        for name, rows in (
            ("training", parts["train"]),
            ("calibration", parts["calibration"]),
            ("final", parts["final"]),
            ("full", list(range(len(data.indices)))),
        ):
            m, a = measure(data, model, rows)
            metrics[name] = m
            arrays.update({name + "_" + k: v for k, v in a.items()})
            if not passed(m, 1):
                raise ValueError(f"{name} exceeds fixed acceptance budgets")
        evaluation_seconds = time.monotonic() - t
        t = time.monotonic()
        archive = output / "fit.json"
        model.save(archive)
        restored = EPWFrequencyCrossFit.load(archive)
        for name in ("calibration", "final", "full"):
            rows = parts[name] if name != "full" else list(range(len(data.indices)))
            if not np.array_equal(
                model.evaluate_matsubara(data.indices[rows]),
                restored.evaluate_matsubara(data.indices[rows]),
            ):
                raise ValueError("archive reload is not exactly deterministic")
        recovery = restored.reconstruct_native()
        replay = output / "original-replay.dat"
        write_epw_anisotropic(data, replay)
        if replay.read_bytes() != path.read_bytes():
            raise ValueError("native lossless replay mismatch")
        fitted = output / "fitted-native.dat"
        fitted_context = write_epw_anisotropic(
            recovery.data, fitted, preserve_source=False
        )
        fitted_data = read_epw_anisotropic(fitted, fitted_context)
        if fitted_data.z.shape != data.z.shape:
            raise ValueError("fitted native export mismatch")
        arrays.update(
            coefficients=model.coefficients,
            frequencies_ev=data.frequencies_ev,
            state_energies_ev=data.state_energies_ev,
        )
        np.savez_compressed(output / "scientific-arrays.npz", **arrays)
        archive_seconds = time.monotonic() - t
        if source_environment() != env:
            raise ValueError("source changed during case")
        raw_binary = sum(
            a.nbytes
            for a in (
                data.frequencies_ev,
                data.state_energies_ev,
                data.z,
                data.delta_ev,
                data.chi_ev,
            )
        )
        record.update(
            status="passed",
            metrics=metrics,
            diagnostics=model.diagnostics,
            recovery=dataclasses.asdict(recovery.diagnostics),
            fitted_context=fitted_context,
            timing_seconds={
                "construction_and_fit": fit_seconds,
                **measured_fit_timings(model),
                "evaluation": evaluation_seconds,
                "archive_and_export": archive_seconds,
                "total": time.monotonic() - started,
            },
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            * (1 if sys.platform == "darwin" else 1024),
            bytes={
                "raw_ascii": path.stat().st_size,
                "raw_binary_arrays": raw_binary,
                "gzip_ascii": len(gzip.compress(path.read_bytes(), mtime=0)),
                "fit_npz_plus_json": (output / "fit.npz").stat().st_size
                + archive.stat().st_size,
            },
            files_sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()},
            scope=(
                "finite native representation only; no experimental prediction "
                "or speedup claim"
            ),
        )
    except Exception as error:
        record.update(
            status="failed",
            error=f"{type(error).__name__}: {error}",
            timing_seconds={"total": time.monotonic() - started},
        )
        write(output / "case.json", record)
        raise
    write(output / "case.json", record)
    return record


def campaign(native_root, freeze, offset, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("campaign output exists; no overwrite or resume")
    output.mkdir(parents=True)
    start = time.monotonic()
    manifest = {
        "schema": "kernelgecp.epw-frequency-campaign.v1",
        "status": "running",
        "children": [],
    }
    try:
        manifest["freeze_sha256"] = sha(freeze)
        frozen, env, _, _ = checked_freeze(freeze, native_root, offset)
        manifest["provenance"] = env
        for repeat in range(3):
            remaining = LIMIT - (time.monotonic() - start)
            if remaining <= 0:
                raise TimeoutError("resource deadline")
            case = output / f"case-{repeat}"
            log = output / f"case-{repeat}.log"
            before = time.monotonic()
            cmd = [
                sys.executable,
                "-m",
                "experiments.epw_frequency",
                "case",
                "--native-root",
                str(native_root),
                "--freeze",
                str(freeze),
                "--chemical-potential-offset-ev",
                str(offset),
                "--repeat",
                str(repeat),
                "--output",
                str(case),
            ]
            child = {"repeat": repeat, "command": cmd}
            manifest["children"].append(child)
            write(output / "campaign.json", manifest)
            try:
                with log.open("xb") as stream:
                    result = subprocess.run(
                        cmd,
                        stdout=stream,
                        stderr=subprocess.STDOUT,
                        cwd=ROOT,
                        timeout=remaining,
                    )
                child.update(
                    exit_code=result.returncode,
                    process_wall_seconds=time.monotonic() - before,
                    log_sha256=sha(log),
                )
                receipt = case / "case.json"
                if receipt.exists():
                    child["receipt_sha256"] = sha(receipt)
                if result.returncode != 0 or not receipt.exists():
                    raise ValueError("case failed or omitted receipt")
                rec = load(receipt)
                if (
                    rec.get("status") != "passed"
                    or rec.get("freeze_sha256") != sha(freeze)
                    or rec.get("repeat") != repeat
                    or stable(rec.get("provenance", {})) != stable(env)
                ):
                    raise ValueError("invalid child receipt")
                for name, digest in rec["files_sha256"].items():
                    if Path(name).name != name or sha(case / name) != digest:
                        raise ValueError("child artifact hash mismatch")
                child["status"] = "passed"
            except subprocess.TimeoutExpired as error:
                child.update(
                    status="inconclusive",
                    process_wall_seconds=time.monotonic() - before,
                    log_sha256=sha(log),
                )
                raise TimeoutError("case deadline") from error
            except Exception as error:
                child.update(status="failed", error=str(error))
                raise
            finally:
                write(output / "campaign.json", manifest)
        with np.load(
            output / "case-0/scientific-arrays.npz", allow_pickle=False
        ) as first:
            for repeat in (1, 2):
                with np.load(
                    output / f"case-{repeat}/scientific-arrays.npz", allow_pickle=False
                ) as other:
                    if set(first.files) != set(other.files) or any(
                        not np.array_equal(first[k], other[k]) for k in first.files
                    ):
                        raise ValueError("repeated scientific arrays differ")
        if source_environment() != env:
            raise ValueError("source changed during campaign")
        manifest.update(status="passed", repeat_arrays_identical=True)
    except Exception as error:
        manifest.update(
            status="inconclusive" if isinstance(error, TimeoutError) else "failed",
            error=f"{type(error).__name__}: {error}",
        )
    manifest["process_wall_seconds"] = time.monotonic() - start
    write(output / "campaign.json", manifest)
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("calibration", "case", "campaign"))
    p.add_argument("--native-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--freeze", type=Path, required=True)
    p.add_argument("--chemical-potential-offset-ev", type=float, required=True)
    p.add_argument("--repeat", type=int, choices=(0, 1, 2))
    p.add_argument("--pilot", action="store_true")
    a = p.parse_args()
    try:
        if a.command == "calibration":
            result = calibrate(
                a.native_root,
                a.output,
                a.freeze,
                a.chemical_potential_offset_ev,
                a.pilot,
            )
        elif a.command == "case":
            if a.pilot or a.repeat is None:
                raise ValueError("case requires repeat and disallows pilot")
            result = run_case(
                a.native_root,
                a.freeze,
                a.chemical_potential_offset_ev,
                a.repeat,
                a.output,
            )
        else:
            if a.pilot or a.repeat is not None:
                raise ValueError("campaign disallows pilot/repeat")
            result = campaign(
                a.native_root, a.freeze, a.chemical_potential_offset_ev, a.output
            )
        if result["status"] != "passed":
            raise ValueError(result.get("error", "inconclusive"))
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
