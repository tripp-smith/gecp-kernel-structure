"""Untouched 64-chain EE validation with immutable inconclusive parent data."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import signal
import subprocess
import tarfile
import tempfile
import tomllib
from pathlib import Path, PurePosixPath
from typing import Any

from experiments.cthyb_dmft import (
    ROOT,
    _bytes,
    _identity,
    _read,
    _sha,
    _write,
    run_final,
)

PARENT_REVISION = "45d42c73ee526c5824344eed99785ce329c96b52"
PARENT_HASHES = {
    "manifest.json": "4cf53b0d639232612e77a14e58658c3559fef0e334253a128bb0dab72d03996c",
    "freeze.json": "836ece257b8569381c648e913192f0585282440ad28b0f9a4aaa5bc43a5b7957",
    "result.json": "8051f86e5038f109c66b02c6bf1df2728ee95cb1c3f50d70864c4bcd8eb7e462",
}
NEW_SOURCES = (
    "experiments/revalidate_cthyb_dmft.py",
    "experiments/cthyb_revalidation.toml",
)


def inventory(folder: Path) -> dict[str, str]:
    result = {}
    for path in sorted(folder.rglob("*")):
        if path.is_symlink():
            raise ValueError("parent evidence cannot contain symlinks")
        if path.is_file() and path.suffix in (".json", ".npz"):
            result[str(path.relative_to(folder))] = _sha(path)
    if not result:
        raise ValueError("empty parent evidence")
    return result


def validate_parent(parent: Path, audit_path: Path) -> dict[str, Any]:
    for relative, expected in PARENT_HASHES.items():
        if _sha(parent / relative) != expected:
            raise ValueError("the retained inconclusive parent was changed")
    manifest = _read(parent / "manifest.json")
    if manifest.get("source_revision") != PARENT_REVISION:
        raise ValueError("wrong parent sampling revision")
    frozen = _read(parent / "freeze.json")
    result = _read(parent / "result.json")
    for files in (frozen["calibration_archives"], result["measurement_archives"]):
        for relative, digest in files.items():
            name = PurePosixPath(relative)
            if (
                name.is_absolute()
                or ".." in name.parts
                or _sha(parent.joinpath(*name.parts)) != digest
            ):
                raise ValueError("retained parent raw-file hash changed")
    for receipt_path in parent.glob("*/*.json"):
        receipt = _read(receipt_path)
        if receipt.get("schema") == "kernelgecp.cthyb-measurement.v1" and (
            receipt.get("identity") != manifest
            or _sha(receipt_path.with_suffix(".npz")) != receipt.get("arrays_sha256")
        ):
            raise ValueError(
                "retained parent chain data or generation identity changed"
            )
    audit = _read(audit_path)
    expected_gates = {
        "loop_controllers",
        "standard_errors",
        "physical_observable_agreement",
        "native_and_represented_self_consistency",
        "all_provider_replay_and_cutoff",
        "independent_native_refinement",
        "independent_low_frequency_quadrature",
    }
    gates = audit.get("gates", {})
    if (
        audit.get("schema") != "kernelgecp.cthyb-audit.v1"
        or audit.get("status") != "inconclusive"
        or audit.get("generation_revision") != PARENT_REVISION
        or set(gates) != expected_gates
        or gates.get("standard_errors") is not False
        or not all(
            value is True for key, value in gates.items() if key != "standard_errors"
        )
        or not isinstance(audit.get("max_standard_error"), float)
        or not 3e-4 < audit["max_standard_error"] < 6e-4
    ):
        raise ValueError("parent requires its independently verified uncertainty miss")
    return {
        "generation_revision": PARENT_REVISION,
        "receipt_sha256": dict(PARENT_HASHES),
        "independent_audit_sha256": _sha(audit_path),
        "inventory": inventory(parent),
    }


def pack_parent(parent: Path, audit_path: Path, destination: Path) -> None:
    """Deterministically retain raw parent chains and their independent audit."""
    proof = validate_parent(parent, audit_path)
    if destination.exists():
        raise FileExistsError("refusing to replace retained parent evidence")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with (
        destination.open("xb") as raw,
        gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed,
        tarfile.open(fileobj=compressed, mode="w|") as archive,
    ):
        paths = {
            "parent/" + relative: parent / relative for relative in proof["inventory"]
        }
        paths["parent-audit.json"] = audit_path
        for name, path in sorted(paths.items()):
            contents = path.read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(contents)
            info.mode = 0o644
            info.mtime = 0
            archive.addfile(info, io.BytesIO(contents))


def extract_parent(
    bundle: Path, destination: Path, expected_hash: str
) -> tuple[Path, Path]:
    if _sha(bundle) != expected_hash:
        raise ValueError("parent bundle hash mismatch")
    if destination.exists():
        return destination / "parent", destination / "parent-audit.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        staging = Path(temporary) / "bundle"
        staging.mkdir()
        seen = set()
        total = 0
        with tarfile.open(bundle, "r:gz") as archive:
            for entry in archive:
                name = PurePosixPath(entry.name)
                total += entry.size
                if (
                    not entry.isfile()
                    or name.is_absolute()
                    or ".." in name.parts
                    or "\\" in entry.name
                    or len(entry.name) > 256
                    or entry.name in seen
                    or entry.size > 64 * 1024**2
                    or total > 512 * 1024**2
                    or not (
                        entry.name.startswith("parent/")
                        or entry.name == "parent-audit.json"
                    )
                ):
                    raise ValueError("unsafe or oversized parent bundle member")
                seen.add(entry.name)
                target = staging.joinpath(*name.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                stream = archive.extractfile(entry)
                if stream is None:
                    raise ValueError("missing bundle member payload")
                with target.open("xb") as output:
                    shutil.copyfileobj(stream, output)
        validate_parent(staging / "parent", staging / "parent-audit.json")
        os.replace(staging, destination)
    return destination / "parent", destination / "parent-audit.json"


def load_protocol(path: Path) -> dict[str, Any]:
    protocol = tomllib.loads(path.read_text())
    fixed = {
        "schema_version": 1,
        "cycles": 400000,
        "cutoff": 8.0,
        "final_replicas": 64,
        "refinement_replicas": 16,
        "seed_base": 28100000,
        "max_standard_error": 0.0003,
        "agreement_absolute": 0.001,
        "agreement_sigma": 3.0,
        "job_minutes": 300,
    }
    if any(
        type(protocol.get(key)) is not type(value) or protocol[key] != value
        for key, value in fixed.items()
    ):
        raise ValueError(
            "revalidation protocol differs from the predeclared scientific contract"
        )
    digest = protocol.get("parent_bundle_sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(x not in "0123456789abcdef" for x in digest)
    ):
        raise ValueError("a SHA256-bound parent bundle is required")
    return protocol


def prepare(
    parent: Path, audit_path: Path, output: Path, protocol: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    proof = validate_parent(parent, audit_path)
    parent_manifest = _read(parent / "manifest.json")
    if parent_manifest["config"].get("seed_mapping") != "odd_stride_2_v1":
        raise ValueError(
            "parent calibration used colliding native RNG streams; fresh calibration required"
        )
    data = dict(parent_manifest["config"])
    for name in ("final_replicas", "refinement_replicas", "seed_base"):
        data[name] = protocol[name]
    if data["seed_base"] == parent_manifest["config"]["seed_base"]:
        raise ValueError("new validation cannot reuse old seed families")
    identity = _identity(data)
    revision = identity["source_revision"]
    for name in NEW_SOURCES:
        digest = _sha(ROOT / name)
        blob = subprocess.check_output(["git", "show", revision + ":" + name], cwd=ROOT)
        if hashlib.sha256(blob).hexdigest() != digest:
            raise ValueError("revalidation implementation/protocol must be committed")
        identity["sources"][name] = digest
    identity["revalidation"] = {
        "schema": "kernelgecp.cthyb-revalidation.v1",
        "parent_generation_revision": PARENT_REVISION,
        "parent_receipt_sha256": proof["receipt_sha256"],
        "parent_independent_audit_sha256": proof["independent_audit_sha256"],
        "parent_inventory_sha256": hashlib.sha256(
            _bytes(proof["inventory"])
        ).hexdigest(),
        "parent_status": "inconclusive_uncertainty_not_relabelled",
        "cycles_and_cutoff": "inherited_parent_calibration",
        "new_validation": "fresh_loops_and_64_new_chains_per_path_no_pooling",
    }
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "manifest.json"
    if manifest.exists():
        prior = _read(manifest)
        if {k: v for k, v in prior.items() if k != "source_revision"} != {
            k: v for k, v in identity.items() if k != "source_revision"
        }:
            raise ValueError("revalidation source/config/runtime or parent changed")
        identity = prior
    else:
        if any(
            (output / name).exists()
            for name in (
                "dense-loop",
                "kernelgecp-loop",
                "dense-final",
                "kernelgecp-final",
            )
        ):
            raise ValueError("new freeze cannot follow preexisting final data")
        _write(manifest, identity)
    retained = output / "parent"
    if not retained.exists():
        shutil.copytree(parent, retained)
    if inventory(retained) != proof["inventory"]:
        raise ValueError("retained raw parent differs from original provenance")
    retained_audit = output / "parent-audit.json"
    if not retained_audit.exists():
        shutil.copyfile(audit_path, retained_audit)
    if _sha(retained_audit) != proof["independent_audit_sha256"]:
        raise ValueError("retained parent audit changed")
    frozen = {
        "identity": identity,
        "status": "frozen_calibration_not_final_acceptance",
        "cycles": protocol["cycles"],
        "cutoff": protocol["cutoff"],
        "selection_scope": "inherited_parent_calibration_new_replica_design",
        "parent_proof": proof,
        "protocol": protocol,
        "final_data_untouched_at_freeze": True,
        "calibration_archives": {
            "parent/" + name: digest for name, digest in proof["inventory"].items()
        }
        | {"parent-audit.json": proof["independent_audit_sha256"]},
    }
    path = output / "freeze.json"
    if path.exists():
        if _read(path) != json.loads(_bytes(frozen)):
            raise ValueError("revalidation freeze changed")
    else:
        _write(path, frozen)
    return data, identity, frozen


def main() -> int:
    import fcntl

    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    pack = commands.add_parser("pack-parent")
    pack.add_argument("--parent", type=Path, required=True)
    pack.add_argument("--parent-audit", type=Path, required=True)
    pack.add_argument("--output", type=Path, required=True)
    run = commands.add_parser("run")
    run.add_argument("--parent", type=Path)
    run.add_argument("--parent-audit", type=Path)
    run.add_argument("--parent-bundle", type=Path)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument(
        "--protocol", type=Path, default=ROOT / "experiments/cthyb_revalidation.toml"
    )
    args = parser.parse_args()
    if args.command == "pack-parent":
        pack_parent(args.parent, args.parent_audit, args.output)
        print(_sha(args.output))
        return 0
    protocol = load_protocol(args.protocol)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / ".campaign.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.parent_bundle is not None:
            parent, audit_path = extract_parent(
                args.parent_bundle,
                args.output / "parent-bundle",
                protocol["parent_bundle_sha256"],
            )
        elif args.parent is not None and args.parent_audit is not None:
            parent, audit_path = args.parent, args.parent_audit
        else:
            raise ValueError("provide a parent bundle or both original parent paths")
        data, identity, frozen = prepare(parent, audit_path, args.output, protocol)

        def expired(signum: int, frame: Any) -> None:
            raise TimeoutError("checkpointed scientific job budget exhausted")

        signal.signal(signal.SIGALRM, expired)
        signal.alarm(60 * protocol["job_minutes"])
        try:
            result = run_final(data, args.output, identity, frozen)
        except TimeoutError:
            _write(
                args.output / "budget-exhausted.json",
                {"status": "inconclusive_job_budget", "identity": identity},
            )
            return 2
        finally:
            signal.alarm(0)
    print(json.dumps({"status": result["status"], "checks": result["checks"]}))
    return 0 if result["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
