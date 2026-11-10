"""Frozen, fresh-process MaxEnt calibration and reproducibility campaign.

This orchestrator never searches controls or retries failed cases. The provider
and independent scientific audit are separate obligations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from experiments import spectroscopy as science

LIMIT_SECONDS = 330 * 60
SHARD_LIMIT_SECONDS = 300 * 60
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


def environment(python, profile="v2"):
    # A different interpreter is not silently substituted, even if its version matches.
    if Path(python).resolve() != Path(sys.executable).resolve():
        raise ValueError("--python must be the interpreter running this driver")
    p = science.provenance() if profile == "v2" else science.provenance(profile)
    if not p.get("tracked_clean") or any(
        p.get("threads", {}).get(k) != "1" for k in science.THREADS
    ):
        raise ValueError("clean tracked source and single-thread environment required")
    p["python_executable"] = str(Path(sys.executable).resolve())
    return p


def stable_environment(p):
    return {k: v for k, v in p.items() if k not in ("source_revision", "tracked_clean")}


def controls(profile="v2"):
    if profile not in ("v2", "v3", "v4"):
        raise ValueError("unsupported spectroscopy profile")
    result = {
        "profile": "corrective_projection_v2_synthetic_binary64",
        "rcond": 1e-10,
        "precision_bits": 53,
        "cases": list(science.CASES),
        "methods": list(science.METHODS),
        "defaults": list(DEFAULTS),
        "eps": 1e-8,
        "n_omega": 401,
        "n_alpha": 33,
        "minimizer_mu0": 1e-8,
        "calibration_seeds": [101, 102, 103],
        "calibration_query_order": 24,
        "final_seeds": [4001, 4002, 4003],
        "final_query_order": 64,
        "bootstrap_first_repeat": 8,
        "repeats": REPEATS,
    }
    if profile in ("v3", "v4"):
        result.update(
            profile="constrained_mass_v3_synthetic_binary64",
            final_seeds=[5001, 5002, 5003],
            final_query_order=80,
            mass_mode="constrained",
            mass_curvature="projected",
            mass_interior_fraction=1e-12,
            minimizer_marquardt=True,
            maxiter=1000,
            calibration_scope="consumed_development",
        )
    if profile == "v4":
        result.update(
            profile="constrained_mass_v4_synthetic_binary64",
            mass_curvature_floor_relative=1e-10,
        )
    return result


def jobs(calibration, profile="v2"):
    controls(profile)
    for ci, case in enumerate(science.CASES):
        for method in science.METHODS:
            for model in DEFAULTS:
                for repeat in range(1 if calibration else REPEATS):
                    job = {
                        "case": case,
                        "method": method,
                        "default_model": model,
                        "seed": (
                            101
                            if calibration
                            else (5001 if profile in ("v3", "v4") else 4001)
                        )
                        + ci,
                        "eps": 1e-8,
                        "rcond": 1e-10,
                        "precision_bits": 53,
                        "n_omega": 401,
                        "n_alpha": 33,
                        "minimizer_mu0": 1e-8,
                        "query_order": 24
                        if calibration
                        else (80 if profile in ("v3", "v4") else 64),
                        "bootstrap": 0 if calibration or repeat else 8,
                        "pilot": calibration,
                        "repeat": repeat,
                    }
                    if profile in ("v3", "v4"):
                        job["profile"] = profile
                    yield job


def stem(job):
    return f"{job['case']}-{job['method']}-{job['default_model']}-r{job['repeat']}"


def validate_child(job, result, env, directory):
    if (
        result.get("schema")
        != "kernelgecp.spectroscopy-case." + job.get("profile", "v2")
        or result.get("status") != "passed"
        or not isinstance(result.get("gates"), dict)
        or not result["gates"]
        or not all(v is True for v in result["gates"].values())
    ):
        raise ValueError("failed or missing child gates/status")
    for key in ("case", "method", "seed", "eps", "pilot", "rcond", "precision_bits"):
        if result.get(key) != job[key]:
            raise ValueError(f"child {key} mismatch")
    if type(result.get("precision_bits")) is not int:
        raise ValueError("child precision_bits must be an integer")
    if (
        result.get("bootstrap_count") != job["bootstrap"]
        or result.get("query_order") != job["query_order"]
    ):
        raise ValueError("child query/bootstrap mismatch")
    config = result.get("config", {})
    for k in ("default_model", "n_omega", "n_alpha", "minimizer_mu0"):
        if config.get(k) != job[k]:
            raise ValueError(f"child config {k} mismatch")
    if job.get("profile") in ("v3", "v4") and (
        config.get("mass_mode") != "constrained"
        or config.get("mass_curvature") != "projected"
        or config.get("mass_interior_fraction") != 1e-12
        or config.get("minimizer_marquardt") is not True
        or config.get("maxiter") != 1000
    ):
        raise ValueError("child v3 mass controls mismatch")
    if (
        job.get("profile") == "v4"
        and config.get("mass_curvature_floor_relative") != 1e-10
    ) or (
        job.get("profile", "v2") != "v4" and "mass_curvature_floor_relative" in config
    ):
        raise ValueError("child curvature floor mismatch")
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
        "--rcond",
        str(job["rcond"]),
        "--precision-bits",
        str(job["precision_bits"]),
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
    if job.get("profile") in ("v3", "v4"):
        cmd += ["--profile", job["profile"]]
    return cmd + (["--pilot"] if job["pilot"] else [])


def run_jobs(directory, manifest, env, python, deadline, profile="v2", selected=None):
    directory.mkdir(parents=True, exist_ok=True)
    records = manifest["children"]
    for job in (
        selected
        if selected is not None
        else jobs(manifest["kind"] == "calibration", profile)
    ):
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
            if (
                environment(python) if profile == "v2" else environment(python, profile)
            ) != env:
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
                    out.with_name(out.stem + "-failed-solve.npz"),
                ]
                if f.is_file()
            }
            manifest.update(status="failed", error=f"{name}: {error}")
            return
        finally:
            write(directory / "campaign.json", manifest)
    manifest["status"] = "passed"


def compare_repeats(directory, manifest, profile="v2", selected=None):
    skipped = []
    for job in selected if selected is not None else jobs(False, profile):
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


def execute(kind, directory, freeze_path, python=sys.executable, *, profile="v2"):
    directory = Path(directory).resolve()
    freeze_path = Path(freeze_path).resolve()
    env = environment(python) if profile == "v2" else environment(python, profile)
    cfg = controls(profile)
    started = time.monotonic()
    mp = directory / "campaign.json"
    if kind == "calibration":
        if mp.exists() or freeze_path.exists():
            raise ValueError("calibration outputs already exist")
        manifest = {
            "schema": "kernelgecp.spectroscopy-campaign." + profile,
            "kind": kind,
            "status": "running",
            "controls": cfg,
            "provenance": env,
            "children": {},
        }
    else:
        frozen = read(freeze_path)
        if (
            frozen.get("schema") != "kernelgecp.spectroscopy-freeze." + profile
            or frozen.get("status") != "passed"
            or frozen.get("controls") != cfg
            or stable_environment(frozen.get("provenance", {}))
            != stable_environment(env)
        ):
            raise ValueError("freeze controls, source or runtime mismatch")
        if mp.exists():
            manifest = read(mp)
            if (
                manifest.get("schema") != "kernelgecp.spectroscopy-campaign." + profile
                or manifest.get("kind") != kind
                or manifest.get("controls") != cfg
                or manifest.get("freeze_sha256") != sha(freeze_path)
                or stable_environment(manifest.get("provenance", {}))
                != stable_environment(env)
            ):
                raise ValueError("resume configuration/source mismatch")
        else:
            manifest = {
                "schema": "kernelgecp.spectroscopy-campaign." + profile,
                "kind": kind,
                "status": "running",
                "controls": cfg,
                "provenance": env,
                "freeze_sha256": sha(freeze_path),
                "children": {},
            }
    manifest["status"] = "running"
    try:
        run_jobs(directory, manifest, env, python, started + LIMIT_SECONDS, profile)
        if manifest["status"] == "passed" and kind == "campaign":
            compare_repeats(directory, manifest, profile)
        if manifest["status"] == "passed" and kind == "calibration":
            freeze = {
                "schema": "kernelgecp.spectroscopy-freeze." + profile,
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


def shard_jobs(ordinal):
    if type(ordinal) is not int or not 0 <= ordinal < 6:
        raise ValueError("shard ordinal must be 0 through 5")
    all_jobs = list(jobs(False, "v4"))
    if len(all_jobs) != 72:
        raise ValueError("unexpected v4 job partition")
    return all_jobs[ordinal * 12 : (ordinal + 1) * 12]


def checked_freeze(freeze_path, env):
    frozen = read(freeze_path)
    if (
        frozen.get("schema") != "kernelgecp.spectroscopy-freeze.v4"
        or frozen.get("status") != "passed"
        or frozen.get("controls") != controls("v4")
        or frozen.get("provenance") != env
    ):
        raise ValueError("freeze controls, source or runtime mismatch")
    return sha(freeze_path)


def expected_files(selected):
    names = set()
    for job in selected:
        name = stem(job)
        names.update(
            name + suffix
            for suffix in (".json", ".npz", ".log", "-dataset.json", "-dataset.npz")
        )
        if job["method"] == "kernelgecp":
            names.add(name + "-gecp-fit.npz")
    return names


def inventory(directory, selected):
    names = expected_files(selected)
    actual = {p.name for p in directory.iterdir() if p.name != "campaign.json"}
    if actual != names or any(
        not (directory / n).is_file() or (directory / n).is_symlink() for n in names
    ):
        raise ValueError("shard file inventory mismatch")
    return {name: sha(directory / name) for name in sorted(names)}


def execute_shard(directory, freeze_path, ordinal, python=sys.executable):
    directory, freeze_path = Path(directory).resolve(), Path(freeze_path).resolve()
    selected = shard_jobs(ordinal)
    env = environment(python, "v4")
    digest = checked_freeze(freeze_path, env)
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    manifest = {
        "schema": "kernelgecp.spectroscopy-campaign-shard.v4",
        "kind": "campaign",
        "ordinal": ordinal,
        "status": "running",
        "controls": controls("v4"),
        "provenance": env,
        "freeze_sha256": digest,
        "children": {},
    }
    write(directory / "campaign.json", manifest)
    try:
        run_jobs(
            directory,
            manifest,
            env,
            python,
            started + SHARD_LIMIT_SECONDS,
            "v4",
            selected,
        )
        if manifest["status"] == "passed":
            compare_repeats(directory, manifest, "v4", selected)
            manifest["inventory"] = inventory(directory, selected)
    except Exception as error:
        manifest.update(status="failed", error=str(error))
    manifest["process_wall_seconds"] = time.monotonic() - started
    write(directory / "campaign.json", manifest)
    return manifest


def assemble_shards(directory, freeze_path, shard_root, python=sys.executable):
    directory, freeze_path, shard_root = map(
        lambda p: Path(p).resolve(), (directory, freeze_path, shard_root)
    )
    env = environment(python, "v4")
    digest = checked_freeze(freeze_path, env)
    started = time.monotonic()
    if {p.name for p in shard_root.iterdir()} != {f"shard-{i}" for i in range(6)}:
        raise ValueError("missing or extra shard directories")
    receipts = []
    for ordinal in range(6):
        folder = shard_root / f"shard-{ordinal}"
        if not folder.is_dir() or folder.is_symlink():
            raise ValueError("invalid shard directory")
        receipt_path = folder / "campaign.json"
        if not receipt_path.is_file() or receipt_path.is_symlink():
            raise ValueError("invalid shard receipt file")
        receipt_digest = sha(receipt_path)
        receipt = read(receipt_path)
        selected = shard_jobs(ordinal)
        expected = {stem(j): j for j in selected}
        if (
            receipt.get("schema") != "kernelgecp.spectroscopy-campaign-shard.v4"
            or receipt.get("kind") != "campaign"
            or type(receipt.get("ordinal")) is not int
            or receipt["ordinal"] != ordinal
            or receipt.get("status") != "passed"
            or receipt.get("controls") != controls("v4")
            or receipt.get("provenance") != env
            or receipt.get("freeze_sha256") != digest
            or set(receipt.get("children", {})) != set(expected)
            or receipt.get("repeat_comparison", {}).get("status") != "passed"
            or receipt.get("inventory") != inventory(folder, selected)
        ):
            raise ValueError(f"invalid shard receipt {ordinal}")
        for name, job in expected.items():
            record = receipt["children"][name]
            if (
                record.get("job") != job
                or record.get("status") != "passed"
                or type(record.get("exit_code")) is not int
                or record.get("exit_code") != 0
                or record.get("hashes")
                != validate_child(job, read(folder / (name + ".json")), env, folder)
                or record.get("log_sha256") != sha(folder / (name + ".log"))
            ):
                raise ValueError(f"invalid shard child {name}")
        receipts.append((folder, receipt, receipt_digest))
    directory.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema": "kernelgecp.spectroscopy-campaign.v4",
        "kind": "campaign",
        "status": "running",
        "controls": controls("v4"),
        "provenance": env,
        "freeze_sha256": digest,
        "children": {},
        "shard_receipts": {},
    }
    write(directory / "campaign.json", manifest)
    try:
        (directory / "shard-receipts").mkdir()
        for ordinal, (folder, receipt, receipt_digest) in enumerate(receipts):
            for name, digest_file in receipt["inventory"].items():
                shutil.copyfile(folder / name, directory / name)
                if sha(directory / name) != digest_file:
                    raise ValueError("copied child hash mismatch")
            target = directory / "shard-receipts" / f"shard-{ordinal}.json"
            shutil.copyfile(folder / "campaign.json", target)
            if sha(target) != receipt_digest:
                raise ValueError("copied shard receipt hash mismatch")
            manifest["shard_receipts"][str(ordinal)] = {
                "path": str(target.relative_to(directory)),
                "sha256": receipt_digest,
            }
            manifest["children"].update(receipt["children"])
        compare_repeats(directory, manifest, "v4")
        if checked_freeze(freeze_path, environment(python, "v4")) != digest:
            raise ValueError("freeze changed during assembly")
        manifest["status"] = "passed"
    except Exception as error:
        manifest.update(status="failed", error=str(error))
    manifest["assembly_wall_seconds"] = time.monotonic() - started
    write(directory / "campaign.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("calibration", "campaign"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--profile", choices=["v2", "v3", "v4"], default="v2")
    parser.add_argument("--shard", type=int)
    parser.add_argument("--assemble-shards", type=Path)
    args = parser.parse_args()
    try:
        if args.shard is not None or args.assemble_shards is not None:
            if (
                args.kind != "campaign"
                or args.profile != "v4"
                or (args.shard is not None and args.assemble_shards is not None)
            ):
                raise ValueError(
                    "sharding requires v4 campaign and exactly one shard mode"
                )
            if args.shard is not None:
                result = execute_shard(
                    args.output, args.freeze, args.shard, args.python
                )
            else:
                result = assemble_shards(
                    args.output, args.freeze, args.assemble_shards, args.python
                )
        else:
            result = execute(
                args.kind, args.output, args.freeze, args.python, profile=args.profile
            )
        if result["status"] != "passed":
            raise RuntimeError(result.get("error", result["status"]))
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
