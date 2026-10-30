"""Frozen, bounded fixed-bath CT-HYB acquisition; not an EE DMFT campaign.

Freeze once, then run bounded acquire segments. Each chain is attempted at most
once. An interrupted attempt or inconsistent archive requires investigation,
not deletion or automatic replay. The independent audit is a separate gate.
"""

from __future__ import annotations

import argparse
import dataclasses
import fcntl
import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from experiments import cthyb_spectroscopy as spectroscopy
from experiments.cthyb_dmft import (
    ROOT,
    THREADS,
    _read,
    _sha,
    _write,
    load_measurement,
    save_measurement,
)
from kernelgecp.integrations.cthyb_dmft import CTHYBPlan

BINDING_SHA = "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
REVISION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
SEEDS = [32100001 + 2 * i for i in range(256)]
SCHEMA = "kernelgecp.cthyb-fixed-bath.v1"
RESULT_SCHEMA = "kernelgecp.cthyb-fixed-bath-result.v1"
TIMING_SCHEMA = "kernelgecp.cthyb-fixed-bath-timing.v1"


def role(i: int) -> str:
    return "train" if i < 128 else "heldout"


def _bath_hash(bath: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(bath).tobytes()).hexdigest()


def parent_data(parent: Path, binding_path: Path):
    if _sha(binding_path) != BINDING_SHA:
        raise ValueError("accepted parent binding hash mismatch")
    binding = _read(binding_path)
    if (
        binding.get("schema") != "kernelgecp.cthyb-audit-binding.v1"
        or binding.get("status") != "passed"
        or binding.get("generation_revision") != REVISION
    ):
        raise ValueError("parent binding is not the accepted passed audit")
    files = {
        "audit": "independent-audit.json",
        "manifest": "manifest.json",
        "freeze": "freeze.json",
        "result": "result.json",
    }
    for key, name in files.items():
        if (
            _sha(parent / name)
            != binding["audit_receipt_sha256" if key == "audit" else key + "_sha256"]
        ):
            raise ValueError("parent " + key + " hash mismatch")
    if _read(parent / files["audit"]).get("status") != "passed":
        raise ValueError("parent independent audit did not pass")
    archives = binding["measurement_archives"]
    prior_seeds = set()
    for name, digest in archives.items():
        path = parent / name
        if path.is_symlink() or _sha(path) != digest:
            raise ValueError("parent receipt hash mismatch")
        receipt = _read(path)
        arr = path.with_suffix(".npz")
        if arr.is_symlink() or _sha(arr) != receipt["arrays_sha256"]:
            raise ValueError("parent raw archive hash mismatch")
        prior_seeds.add(receipt["seed"])
    if prior_seeds.intersection(SEEDS):
        raise ValueError("new seed allocation overlaps parent")
    receipt_path = parent / "dense-final/replica-00.json"
    receipt = _read(receipt_path)
    parent_manifest = _read(parent / "manifest.json")
    if receipt["identity"] != parent_manifest:
        raise ValueError("parent bath receipt identity differs from manifest")
    measurement = load_measurement(receipt_path.with_suffix(""), parent_manifest)
    config = measurement.config
    if config.beta != 10 or config.interaction != 2 or measurement.cycles != 1600000:
        raise ValueError("parent native configuration is not the required bath")
    if any(
        dataclasses.asdict(config)[k] != parent_manifest["config"][k]
        for k in dataclasses.asdict(config)
    ):
        raise ValueError("parent manifest and raw native configuration mismatch")
    parent_info = {
        "binding_sha256": _sha(binding_path),
        **{key + "_sha256": _sha(parent / name) for key, name in files.items()},
        "bath_receipt_sha256": _sha(receipt_path),
        "bath_archive_sha256": _sha(receipt_path.with_suffix(".npz")),
        "bath_bytes_sha256": _bath_hash(measurement.bath_iw),
    }
    return parent_info, measurement.bath_iw, config


def provenance():
    if any(os.environ.get(name) != "1" for name in THREADS):
        raise ValueError("all thread environment variables must equal 1")
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()

    # Freeze every repository Python/TOML/lock source, including analysis controls.
    def relevant(path):
        return path.endswith((".py", ".toml", ".lock", ".txt"))

    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "HEAD"], cwd=ROOT, text=True
    ).splitlines()
    paths = sorted(p for p in paths if relevant(p))
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard"], cwd=ROOT, text=True
    ).splitlines()
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", "HEAD"], cwd=ROOT, text=True
    ).splitlines()
    if any(relevant(p) for p in untracked + changed) or any(
        not (ROOT / p).is_file() for p in paths
    ):
        raise ValueError("frozen sources must exactly match committed inventory")
    sources = {p: _sha(ROOT / p) for p in paths}
    # Bind installed provider Python sources and native extensions.
    modules = (
        "numpy",
        "scipy",
        "triqs",
        "triqs_cthyb",
        "pydlr",
        "sparse_ir",
        "pylibsparseir",
    )
    provider_files = {}
    versions = {}
    for name in modules:
        mod = importlib.import_module(name)
        path = Path(mod.__file__).resolve()
        base = path.parent
        for item in sorted(base.rglob("*")):
            if item.is_file() and (
                item.suffix in (".py", ".so", ".dylib", ".pyd") or ".so." in item.name
            ):
                provider_files[str(item)] = _sha(item)
        try:
            versions[name] = importlib.metadata.version(name.replace("_", "-"))
        except importlib.metadata.PackageNotFoundError:
            versions[name] = str(getattr(mod, "__version__", "unknown"))
            if versions[name] == "unknown" and name in ("triqs", "triqs_cthyb"):
                versions[name] = str(importlib.import_module(name + ".version").version)
    # Linux maps contain the actual dynamically loaded provider libraries,
    # including dependencies installed outside the package tree.
    maps = Path("/proc/self/maps")
    if not maps.exists():
        raise ValueError(
            "canonical native binary capture requires Linux /proc/self/maps"
        )
    for line in maps.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 6 and parts[-1].startswith("/"):
            path = Path(parts[-1])
            if path.is_file() and (".so" in path.name or path.suffix == ".dylib"):
                provider_files[str(path.resolve())] = _sha(path)
    return {
        "revision": revision,
        "sources": sources,
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "versions": versions,
            "provider_files": provider_files,
        },
        "threads": {n: os.environ[n] for n in THREADS},
    }


def protocol(config):
    return {
        "seeds": SEEDS,
        "train_count": 128,
        "cycles": 1600000,
        "config": dataclasses.asdict(config),
        "segment_seconds": 14400,
        "memory_limit_bytes": 4294967296,
    }


CALIBRATION_GATES = frozenset(
    {
        "ee_source_and_raw_chains",
        "frozen_protocol_and_hashes",
        "covariance_and_projection",
        "native_trials_and_quadrature",
        "heldout_rms",
        "bootstrap_and_defaults",
        "repeat_arrays",
        "fresh_process_receipts",
    }
)


def calibration_binding(
    receipt: Path, receipt_sha256: str, audit_path: Path, audit_sha256: str
):
    """Pin calibration to the exact current analysis code and provider runtime."""
    snapshot = spectroscopy.snapshot()
    cal = spectroscopy.verify_calibration(receipt, receipt_sha256, snapshot)
    audit = spectroscopy.verify_calibration_audit(
        audit_path, audit_sha256, receipt_sha256, cal["binding_sha256"]
    )
    gates = audit.get("gates")
    if (
        not isinstance(gates, dict)
        or set(gates) != CALIBRATION_GATES
        or not all(v is True for v in gates.values())
        or audit.get("cases") != 8
        or audit.get("campaign_sha256") != receipt_sha256
        or audit.get("auditor_revision") != snapshot["revision"]
        or audit.get("auditor_sha256")
        != snapshot["sources"].get("experiments/audit_cthyb_spectroscopy.py")
    ):
        raise ValueError("calibration independent audit gates or source mismatch")
    return {
        "receipt_sha256": receipt_sha256,
        "audit_sha256": audit_sha256,
        "binding_sha256": cal["binding_sha256"],
        "auditor_revision": audit["auditor_revision"],
        "auditor_sha256": audit["auditor_sha256"],
    }


def freeze(
    parent: Path,
    binding: Path,
    output: Path,
    calibration: Path,
    calibration_sha256: str,
    calibration_audit: Path,
    calibration_audit_sha256: str,
):
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("freeze requires an empty output directory")
    info, _, config = parent_data(parent, binding)
    CTHYBPlan(
        config
    )  # Load the exact native provider before recording shared libraries.
    cal = calibration_binding(
        calibration, calibration_sha256, calibration_audit, calibration_audit_sha256
    )
    manifest = {
        "schema": SCHEMA,
        "protocol": protocol(config),
        "parent": info,
        "calibration": cal,
        "provenance": provenance(),
    }
    output.mkdir(parents=True, exist_ok=True)
    _write(output / "manifest.json", manifest)
    _write(
        output / "result.json",
        {
            "schema": RESULT_SCHEMA,
            "manifest_sha256": _sha(output / "manifest.json"),
            "chains": [],
            "status": "in_progress",
        },
    )


def validate(
    output: Path,
    parent: Path,
    binding: Path,
    calibration: Path,
    calibration_sha256: str,
    calibration_audit: Path,
    calibration_audit_sha256: str,
):
    manifest = _read(output / "manifest.json")
    info, bath, config = parent_data(parent, binding)
    CTHYBPlan(config)
    cal = calibration_binding(
        calibration, calibration_sha256, calibration_audit, calibration_audit_sha256
    )
    if manifest != {
        "schema": SCHEMA,
        "protocol": protocol(config),
        "parent": info,
        "calibration": cal,
        "provenance": provenance(),
    }:
        raise ValueError("frozen source, runtime, protocol or parent changed")
    digest = _sha(output / "manifest.json")
    result = _read(output / "result.json")
    entries = result.get("chains")
    if (
        result.get("schema") != RESULT_SCHEMA
        or result.get("manifest_sha256") != digest
        or not isinstance(entries, list)
        or len(entries) > 256
        or result.get("status")
        != ("complete" if len(entries) == 256 else "in_progress")
    ):
        raise ValueError("invalid fixed-bath result")
    folder = output / "chains"
    expected = set()
    for i, entry in enumerate(entries):
        stem = folder / f"chain-{i:03d}"
        identity = {"manifest_sha256": digest, "index": i, "role": role(i)}
        m = load_measurement(stem, identity)
        timing = _read(stem.with_suffix(".timing.json"))
        expected_entry = {
            "index": i,
            "seed": SEEDS[i],
            "role": role(i),
            "receipt_sha256": _sha(stem.with_suffix(".json")),
            "arrays_sha256": _sha(stem.with_suffix(".npz")),
            "timing_sha256": _sha(stem.with_suffix(".timing.json")),
        }
        if (
            entry != expected_entry
            or m.seed != SEEDS[i]
            or m.cycles != 1600000
            or m.config != config
            or not np.array_equal(m.bath_iw, bath)
            or not 0.99 < m.average_sign <= 1
            or m.solve_status != 0
        ):
            raise ValueError("stored chain does not match frozen protocol")
        if (
            timing.get("schema") != TIMING_SCHEMA
            or timing.get("seed") != m.seed
            or timing.get("solve_status") != 0
            or timing.get("average_sign") != m.average_sign
            or timing.get("receipt_sha256") != entry["receipt_sha256"]
            or timing.get("arrays_sha256") != entry["arrays_sha256"]
            or not 0 <= timing.get("elapsed_seconds", -1) <= 14400
            or not 0 < timing.get("peak_rss_bytes", 0) <= 4294967296
        ):
            raise ValueError("invalid timing or native diagnostics")
        expected.update(
            {stem.name + suffix for suffix in (".json", ".npz", ".timing.json")}
        )
    actual = {p.name for p in folder.iterdir()} if folder.exists() else set()
    if actual != expected:
        raise ValueError("orphan, interrupted attempt or unexpected chain file")
    return manifest, result, bath, config


def worker(
    output,
    parent,
    binding,
    calibration,
    calibration_sha256,
    calibration_audit,
    calibration_audit_sha256,
    index,
    lock_fd,
):
    if (
        lock_fd < 0
        or os.fstat(lock_fd).st_ino != (output / ".acquire.lock").stat().st_ino
    ):
        raise ValueError("worker requires the active acquisition lock")
    manifest, result, bath, config = validate(
        output,
        parent,
        binding,
        calibration,
        calibration_sha256,
        calibration_audit,
        calibration_audit_sha256,
    )
    if index != len(result["chains"]):
        raise ValueError("chain must be the next predeclared index")
    folder = output / "chains"
    folder.mkdir(exist_ok=True)
    stem = folder / f"chain-{index:03d}"
    marker = stem.with_suffix(".attempt")
    with marker.open("x") as f:
        f.write(json.dumps({"index": index, "seed": SEEDS[index]}) + "\n")
        f.flush()
        os.fsync(f.fileno())
    start = time.monotonic()
    try:
        m = CTHYBPlan(config).sample(bath, cycles=1600000, seed=SEEDS[index])
    except BaseException as error:
        _write(
            stem.with_suffix(".failure.json"),
            {
                "schema": "kernelgecp.cthyb-fixed-bath-failure.v1",
                "index": index,
                "seed": SEEDS[index],
                "elapsed_seconds": time.monotonic() - start,
                "exception_type": type(error).__name__,
                "message": str(error),
            },
        )
        raise
    elapsed = time.monotonic() - start
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss = int(rss if platform.system() == "Darwin" else rss * 1024)
    # Save the returned raw data before validating diagnostics, so a bad sign,
    # resource overrun or changed source cannot discard an attempted chain.
    identity = {
        "manifest_sha256": result["manifest_sha256"],
        "index": index,
        "role": role(index),
    }
    save_measurement(stem, m, identity)
    if (
        m.seed != SEEDS[index]
        or m.config != config
        or m.cycles != 1600000
        or not np.array_equal(m.bath_iw, bath)
        or m.solve_status != 0
        or not 0.99 < m.average_sign <= 1
        or rss > 4294967296
        or elapsed > 14400
    ):
        raise ValueError("native solve or resource validation failed")
    # An interrupted or invalid chain remains accompanied by its attempt marker.
    if (
        calibration_binding(
            calibration, calibration_sha256, calibration_audit, calibration_audit_sha256
        )
        != manifest["calibration"]
        or manifest["provenance"] != provenance()
        or parent_data(parent, binding)[0] != manifest["parent"]
    ):
        raise ValueError("frozen execution identity changed during chain")
    entry = {
        "index": index,
        "seed": m.seed,
        "role": role(index),
        "receipt_sha256": _sha(stem.with_suffix(".json")),
        "arrays_sha256": _sha(stem.with_suffix(".npz")),
    }
    timing = {
        "schema": TIMING_SCHEMA,
        "seed": m.seed,
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": rss,
        "solve_status": m.solve_status,
        "average_sign": m.average_sign,
        "receipt_sha256": entry["receipt_sha256"],
        "arrays_sha256": entry["arrays_sha256"],
    }
    _write(stem.with_suffix(".timing.json"), timing)
    entry["timing_sha256"] = _sha(stem.with_suffix(".timing.json"))
    result["chains"].append(entry)
    result["status"] = "complete" if len(result["chains"]) == 256 else "in_progress"
    _write(output / "result.json", result)
    marker.unlink()


def acquire(
    parent,
    binding,
    output,
    calibration,
    calibration_sha256,
    calibration_audit,
    calibration_audit_sha256,
    max_new,
):
    if type(max_new) is not int or not 1 <= max_new <= 256:
        raise ValueError("max-new-chains must be between 1 and 256")
    # The lock is outside the chain inventory and held across the entire segment.
    with (output / ".acquire.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        deadline = time.monotonic() + 14400
        for _ in range(max_new):
            _, result, _, _ = validate(
                output,
                parent,
                binding,
                calibration,
                calibration_sha256,
                calibration_audit,
                calibration_audit_sha256,
            )
            if len(result["chains"]) == 256:
                break
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            cmd = [
                sys.executable,
                "-m",
                "experiments.cthyb_fixed_bath",
                "--parent",
                str(parent),
                "--binding",
                str(binding),
                "--output",
                str(output),
                "--calibration",
                str(calibration),
                "--calibration-sha256",
                calibration_sha256,
                "--calibration-audit",
                str(calibration_audit),
                "--calibration-audit-sha256",
                calibration_audit_sha256,
                "worker",
                "--index",
                str(len(result["chains"])),
                "--lock-fd",
                str(lock.fileno()),
            ]
            subprocess.run(
                cmd, cwd=ROOT, check=True, timeout=remaining, pass_fds=(lock.fileno(),)
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--binding", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--calibration-sha256", required=True)
    parser.add_argument("--calibration-audit", type=Path, required=True)
    parser.add_argument("--calibration-audit-sha256", required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("freeze")
    a = sub.add_parser("acquire")
    a.add_argument("--max-new-chains", type=int, required=True)
    w = sub.add_parser("worker")
    w.add_argument("--index", type=int, required=True)
    w.add_argument("--lock-fd", type=int, required=True)
    args = parser.parse_args()
    calargs = (
        args.calibration,
        args.calibration_sha256,
        args.calibration_audit,
        args.calibration_audit_sha256,
    )
    if args.command == "freeze":
        freeze(args.parent, args.binding, args.output, *calargs)
    elif args.command == "acquire":
        acquire(args.parent, args.binding, args.output, *calargs, args.max_new_chains)
    else:
        worker(
            args.output, args.parent, args.binding, *calargs, args.index, args.lock_fd
        )


if __name__ == "__main__":
    main()
