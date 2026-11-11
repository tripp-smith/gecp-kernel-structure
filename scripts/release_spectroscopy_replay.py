#!/usr/bin/env python3
"""Authenticate and replay the canonical saved spectroscopy evidence."""

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from scripts import release_science_payloads as payloads
from scripts import release_science_source_equivalence as science_equivalence
from scripts import release_source_equivalence as version_equivalence
from scripts.providers.verify_receipt import read_receipt
from scripts.release_science_fetch import fetch_spectroscopy_origin
from scripts.release_science_runs import REPOSITORY

PIN = Path("data/release-science-origin-1.2.0.json")
SCHEMA = "kernelgecp.release-science-origin-pin.v1"
COUNTS = {
    "synthetic_calibration": 24,
    "synthetic_final": 72,
    "measured_calibration": 8,
    "measured_final": 24,
    "raw_chains": 256,
}


def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(pin):
    if (
        not isinstance(pin, dict)
        or set(pin) != {"schema", "revision", "run_id"}
        or pin["schema"] != SCHEMA
        or not isinstance(pin["revision"], str)
        or not re.fullmatch(r"[0-9a-f]{40}", pin["revision"])
        or type(pin["run_id"]) is not int
        or pin["run_id"] <= 0
    ):
        raise ValueError("invalid origin pin")
    return pin["revision"], pin["run_id"]


def verify_checkout(repo, candidate):
    if not isinstance(candidate, str) or not re.fullmatch(r"[0-9a-f]{40}", candidate):
        raise ValueError("candidate must be a full lowercase commit ID")
    head = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=all"],
        text=True,
    )
    if head != candidate or status:
        raise ValueError("candidate checkout must be exact and clean")


def compare_candidate(repo, source, candidate):
    version = subprocess.check_output(
        ["git", "-C", str(repo), "show", candidate + ":pyproject.toml"], text=True
    )
    import tomllib

    value = tomllib.loads(version)["project"]["version"]
    if value == "1.1.0":
        return science_equivalence.compare(repo, source, candidate)
    if value == "1.2.0":
        parents = subprocess.check_output(
            ["git", "-C", str(repo), "rev-list", "--parents", "-n", "1", candidate],
            text=True,
        ).split()
        if len(parents) != 2:
            raise ValueError("release candidate must have one parent")
        producer = parents[1]
        result = science_equivalence.compare(repo, source, producer)
        guard = version_equivalence.compare(str(repo), candidate, producer)
        if guard["version_transition"] != "1.1.0 -> 1.2.0":
            raise ValueError("wrong release version transition")
        return {"science_to_producer": result, "producer_to_candidate": guard}
    raise ValueError("unsupported candidate version")


def replay(repo, candidate, source, run_id, output, python):
    output = Path(output)
    if output.exists() or output.is_symlink() or not output.parent.is_dir():
        raise ValueError("output must be new")
    if Path(repo).resolve() in output.resolve().parents:
        raise ValueError("output must be outside checkout")
    verify_checkout(repo, candidate)
    equivalence = compare_candidate(repo, source, candidate)
    origin = fetch_spectroscopy_origin(REPOSITORY, source, run_id)
    artifact = origin["artifact"]
    output.mkdir()
    (output / "origin.json").write_text(
        json.dumps(origin, indent=2, sort_keys=True) + "\n"
    )
    (output / "equivalence.json").write_text(
        json.dumps(equivalence, indent=2, sort_keys=True) + "\n"
    )
    archive = output / (artifact["name"] + ".zip")
    payloads.download(artifact["id"], archive, output / "download.log")
    if (
        archive.stat().st_size != artifact["size_in_bytes"]
        or sha256_file(archive) != artifact["digest"][7:]
    ):
        raise ValueError("download digest or size mismatch")
    env = dict(os.environ)
    env.pop("GH_TOKEN", None)
    env.pop("GITHUB_TOKEN", None)
    command = [
        str(python),
        "-m",
        "scripts.release_spectroscopy_v4",
        "--archive",
        str(archive),
        "--archive-sha256",
        artifact["digest"][7:],
        "--archive-size",
        str(artifact["size_in_bytes"]),
        "--run-id",
        str(run_id),
        "--source-revision",
        source,
        "--triqs-python",
        str(python),
        "--output",
        str(output / "audit"),
    ]
    subprocess.run(command, cwd=repo, env=env, check=True, timeout=18000)
    receipt = read_receipt(output / "audit/receipt.json")
    if (
        not isinstance(receipt, dict)
        or receipt.get("schema") != "kernelgecp.spectroscopy-v4-saved-replay.v1"
        or receipt.get("status") != "passed"
        or receipt.get("source_revision") != source
        or receipt.get("run_id") != str(run_id)
        or receipt.get("archive_sha256") != artifact["digest"][7:]
        or receipt.get("archive_size") != artifact["size_in_bytes"]
        or receipt.get("counts") != COUNTS
    ):
        raise ValueError("invalid wrapper receipt")
    hashes = receipt.get("receipt_sha256")
    names = {
        "synthetic-cal",
        "synthetic-final",
        "measured-cal",
        "raw",
        "measured-final",
    }
    if not isinstance(hashes, dict) or set(hashes) != names:
        raise ValueError("invalid wrapper receipt inventory")
    for name, digest in hashes.items():
        if (
            not isinstance(digest, str)
            or not re.fullmatch(r"[0-9a-f]{64}", digest)
            or sha256_file(output / "audit" / (name + ".json")) != digest
        ):
            raise ValueError("wrapper receipt digest mismatch")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--science-revision")
    parser.add_argument("--run-id", type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    args = parser.parse_args()
    try:
        if (args.science_revision is None) != (args.run_id is None):
            raise ValueError("both preflight inputs are required")
        verify_checkout(Path.cwd(), args.candidate)
        pin = (
            {"schema": SCHEMA, "revision": args.science_revision, "run_id": args.run_id}
            if args.science_revision is not None
            else read_receipt(PIN)
        )
        if args.science_revision is None:
            committed = subprocess.check_output(
                ["git", "show", f"{args.candidate}:{PIN}"]
            )
            if PIN.read_bytes() != committed:
                raise ValueError("origin pin must match committed bytes")
        source, run_id = identity(pin)
        replay(Path.cwd(), args.candidate, source, run_id, args.output, args.python)
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"spectroscopy replay failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
