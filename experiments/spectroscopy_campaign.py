"""Frozen, fresh-process MaxEnt calibration and reproducibility campaign.

This orchestrator never searches controls or retries failed cases. The provider
and independent scientific audit are separate obligations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from experiments import spectroscopy as science

LIMIT_SECONDS = 330 * 60
DEFAULTS = ("flat", "broad_gaussian")
REPEATS = 3


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, path)


def read(path):
    with open(path) as stream:
        value = json.load(stream)
    if not isinstance(value, dict) or not value:
        raise ValueError(f"empty or invalid metadata: {path}")
    return value


def environment(python):
    # A different interpreter is not silently substituted, even if its version matches.
    if Path(python).resolve() != Path(sys.executable).resolve():
        raise ValueError("--python must be the interpreter running this driver")
    p = science.provenance()
    if not p.get("tracked_clean") or any(
        p.get("threads", {}).get(k) != "1" for k in science.THREADS
    ):
        raise ValueError("clean tracked source and single-thread environment required")
    p["python_executable"] = str(Path(sys.executable).resolve())
    return p


def stable_environment(p):
    return {k: v for k, v in p.items() if k not in ("source_revision", "tracked_clean")}


def controls():
    return {
        "cases": list(science.CASES),
        "methods": list(science.METHODS),
        "defaults": list(DEFAULTS),
        "eps": 1e-8,
        "n_omega": 401,
        "n_alpha": 33,
        "minimizer_mu0": 1e-8,
        "calibration_seeds": [101, 102, 103],
        "calibration_query_order": 24,
        "final_seeds": [3001, 3002, 3003],
        "final_query_order": 64,
        "bootstrap_first_repeat": 8,
        "repeats": REPEATS,
    }


def jobs(calibration):
    for ci, case in enumerate(science.CASES):
        for method in science.METHODS:
            for model in DEFAULTS:
                for repeat in range(1 if calibration else REPEATS):
                    yield {
                        "case": case,
                        "method": method,
                        "default_model": model,
                        "seed": (101 if calibration else 3001) + ci,
                        "eps": 1e-8,
                        "n_omega": 401,
                        "n_alpha": 33,
                        "minimizer_mu0": 1e-8,
                        "query_order": 24 if calibration else 64,
                        "bootstrap": 0 if calibration or repeat else 8,
                        "pilot": calibration,
                        "repeat": repeat,
                    }


def stem(job):
    return f"{job['case']}-{job['method']}-{job['default_model']}-r{job['repeat']}"


def validate_child(job, result, env, directory):
    if (
        result.get("status") != "passed"
        or not isinstance(result.get("gates"), dict)
        or not result["gates"]
        or not all(v is True for v in result["gates"].values())
    ):
        raise ValueError("failed or missing child gates/status")
    for key in ("case", "method", "seed", "eps", "pilot"):
        if result.get(key) != job[key]:
            raise ValueError(f"child {key} mismatch")
    if (
        result.get("bootstrap_count") != job["bootstrap"]
        or result.get("query_order") != job["query_order"]
    ):
        raise ValueError("child query/bootstrap mismatch")
    config = result.get("config", {})
    for k in ("default_model", "n_omega", "n_alpha", "minimizer_mu0"):
        if config.get(k) != job[k]:
            raise ValueError(f"child config {k} mismatch")
    if result.get("provenance") != {
        k: v for k, v in env.items() if k != "python_executable"
    }:
        raise ValueError("child provenance mismatch")
    npz = directory / (stem(job) + ".npz")
    if (
        not npz.is_file()
        or not npz.stat().st_size
        or sha(npz) != result.get("array_sha256")
    ):
        raise ValueError("missing or corrupt child arrays")
    with np.load(npz, allow_pickle=False) as arrays:
        if not arrays.files:
            raise ValueError("empty child arrays")
        for k in arrays.files:
            arrays[k]  # force decompression and checksum validation
    dataset = directory / (stem(job) + "-dataset.npz")
    if not dataset.is_file() or not dataset.stat().st_size:
        raise ValueError("missing dataset archive")
    dataset_json = dataset.with_suffix(".json")
    if not dataset_json.is_file() or not dataset_json.stat().st_size:
        raise ValueError("missing dataset application manifest")
    fit = result.get("coefficient_fit_archive")
    if job["method"] == "kernelgecp":
        if not isinstance(fit, dict) or fit.get("path") != stem(job) + "-gecp-fit.npz":
            raise ValueError("missing GECP fit archive metadata")
        if not (directory / fit["path"]).is_file() or sha(
            directory / fit["path"]
        ) != fit.get("sha256"):
            raise ValueError("corrupt GECP fit archive")
    return {
        "json_sha256": sha(directory / (stem(job) + ".json")),
        "npz_sha256": sha(npz),
        "dataset_sha256": sha(dataset),
        "dataset_json_sha256": sha(dataset_json),
        "fit_sha256": fit["sha256"] if fit else None,
    }


def child_command(python, job, output):
    cmd = [
        python,
        "-m",
        "experiments.spectroscopy",
        "--case",
        job["case"],
        "--method",
        job["method"],
        "--seed",
        str(job["seed"]),
        "--eps",
        str(job["eps"]),
        "--n-omega",
        str(job["n_omega"]),
        "--n-alpha",
        str(job["n_alpha"]),
        "--minimizer-mu0",
        str(job["minimizer_mu0"]),
        "--default-model",
        job["default_model"],
        "--query-order",
        str(job["query_order"]),
        "--bootstrap",
        str(job["bootstrap"]),
        "--output",
        str(output),
    ]
    return cmd + (["--pilot"] if job["pilot"] else [])


def run_jobs(directory, manifest, env, python, deadline):
    directory.mkdir(parents=True, exist_ok=True)
    records = manifest["children"]
    for job in jobs(manifest["kind"] == "calibration"):
        name = stem(job)
        paths = list(directory.glob(name + ".*")) + list(
            directory.glob(name + "-*.npz")
        )
        existing = records.get(name)
        if existing is not None:
            try:
                if (
                    existing.get("job") != job
                    or existing.get("exit_code") != 0
                    or existing.get("status") != "passed"
                ):
                    raise ValueError("failed or mismatched previous child; no retry")
                hashes = validate_child(
                    job, read(directory / (name + ".json")), env, directory
                )
                if hashes != existing.get("hashes") or sha(
                    directory / (name + ".log")
                ) != existing.get("log_sha256"):
                    raise ValueError("previous child hash mismatch")
            except Exception as error:
                manifest.update(status="failed", error=f"resume {name}: {error}")
                return
            continue
        if paths:
            manifest.update(status="failed", error=f"orphan child artifacts: {name}")
            return
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            manifest.update(
                status="inconclusive", error="330 minute resource deadline reached"
            )
            return
        out = directory / (name + ".json")
        log = directory / (name + ".log")
        record = {"job": job, "status": "running", "exit_code": None}
        records[name] = record
        write(directory / "campaign.json", manifest)
        started = time.monotonic()
        try:
            with log.open("xb") as stream:
                proc = subprocess.run(
                    child_command(python, job, out),
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    timeout=remaining,
                    cwd=science.ROOT,
                )
            record["exit_code"] = proc.returncode
            record["process_wall_seconds"] = time.monotonic() - started
            record["log_sha256"] = sha(log)
            if proc.returncode != 0:
                raise ValueError(f"child exit {proc.returncode}")
            record["hashes"] = validate_child(job, read(out), env, directory)
            if environment(python) != env:
                raise ValueError("source or runtime changed during child")
            record["status"] = "passed"
        except subprocess.TimeoutExpired:
            record.update(
                status="inconclusive",
                process_wall_seconds=time.monotonic() - started,
                log_sha256=sha(log),
            )
            manifest.update(status="inconclusive", error=f"deadline during {name}")
            return
        except Exception as error:
            record.update(
                status="failed",
                error=str(error),
                process_wall_seconds=time.monotonic() - started,
            )
            if log.exists():
                record["log_sha256"] = sha(log)
            record["available_file_sha256"] = {
                f.name: sha(f)
                for f in [
                    out,
                    out.with_suffix(".npz"),
                    out.with_name(out.stem + "-dataset.npz"),
                    out.with_name(out.stem + "-gecp-fit.npz"),
                ]
                if f.is_file()
            }
            manifest.update(status="failed", error=f"{name}: {error}")
            return
        finally:
            write(directory / "campaign.json", manifest)
    manifest["status"] = "passed"


def compare_repeats(directory, manifest):
    skipped = []
    for job in jobs(False):
        if job["repeat"] == 0:
            continue
        reference = dict(job, repeat=0)
        with (
            np.load(directory / (stem(reference) + ".npz"), allow_pickle=False) as a,
            np.load(directory / (stem(job) + ".npz"), allow_pickle=False) as b,
        ):
            if set(a.files) != set(b.files):
                raise ValueError("repeat array keys differ")
            for name in a.files:
                if name.startswith("bootstrap_"):
                    skipped.append(
                        {
                            "case": stem(job),
                            "array": name,
                            "reason": "bootstrap only in repeat 0",
                        }
                    )
                elif not np.array_equal(a[name], b[name]):
                    raise ValueError(f"repeat array differs: {stem(job)} {name}")
    manifest["repeat_comparison"] = {"status": "passed", "skipped": skipped}


def execute(kind, directory, freeze_path, python=sys.executable):
    directory = Path(directory).resolve()
    freeze_path = Path(freeze_path).resolve()
    env = environment(python)
    cfg = controls()
    started = time.monotonic()
    mp = directory / "campaign.json"
    if kind == "calibration":
        if mp.exists() or freeze_path.exists():
            raise ValueError("calibration outputs already exist")
        manifest = {
            "schema": "kernelgecp.spectroscopy-campaign.v1",
            "kind": kind,
            "status": "running",
            "controls": cfg,
            "provenance": env,
            "children": {},
        }
    else:
        frozen = read(freeze_path)
        if (
            frozen.get("status") != "passed"
            or frozen.get("controls") != cfg
            or stable_environment(frozen.get("provenance", {}))
            != stable_environment(env)
        ):
            raise ValueError("freeze controls, source or runtime mismatch")
        if mp.exists():
            manifest = read(mp)
            if (
                manifest.get("kind") != kind
                or manifest.get("controls") != cfg
                or manifest.get("freeze_sha256") != sha(freeze_path)
                or stable_environment(manifest.get("provenance", {}))
                != stable_environment(env)
            ):
                raise ValueError("resume configuration/source mismatch")
        else:
            manifest = {
                "schema": "kernelgecp.spectroscopy-campaign.v1",
                "kind": kind,
                "status": "running",
                "controls": cfg,
                "provenance": env,
                "freeze_sha256": sha(freeze_path),
                "children": {},
            }
    manifest["status"] = "running"
    try:
        run_jobs(directory, manifest, env, python, started + LIMIT_SECONDS)
        if manifest["status"] == "passed" and kind == "campaign":
            compare_repeats(directory, manifest)
        if manifest["status"] == "passed" and kind == "calibration":
            freeze = {
                "schema": "kernelgecp.spectroscopy-freeze.v1",
                "status": "passed",
                "controls": cfg,
                "provenance": env,
                "calibration_manifest_sha256": None,
                "scope": (
                    "fixed controls; final seeds and query order differ "
                    "from calibration; "
                    "no data-driven control search"
                ),
            }
            manifest["process_wall_seconds"] = time.monotonic() - started
            write(mp, manifest)
            freeze["calibration_manifest_sha256"] = sha(mp)
            write(freeze_path, freeze)
    except Exception as error:
        manifest.update(status="failed", error=str(error))
    manifest["process_wall_seconds"] = time.monotonic() - started
    # Preserve exact calibration manifest bytes referenced by freeze.
    if not (kind == "calibration" and manifest["status"] == "passed"):
        write(mp, manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("calibration", "campaign"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()
    try:
        result = execute(args.kind, args.output, args.freeze, args.python)
        if result["status"] != "passed":
            raise RuntimeError(result.get("error", result["status"]))
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
