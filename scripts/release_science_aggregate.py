"""Bind authenticated producer runs and six fresh replays to a release revision."""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

from scripts.providers.verify_receipt import read_receipt
from scripts.release_science_fetch import fetch_science_runs, fetch_spectroscopy_origin
from scripts.release_science_gate import origin_check, validate_wrapper
from scripts.release_science_payloads import validate_selection
from scripts.release_science_runs import PRODUCTION_PROFILE, REPOSITORY, profile_roles
from scripts.release_science_source_equivalence import compare as compare_science
from scripts.release_source_equivalence import compare as compare_release

PLAN_SCHEMA = "kernelgecp.release-science-plan.v1"
RECEIPT_SCHEMA = "kernelgecp.release-science-acceptance.v1"
ORIGIN_PIN = Path("data/release-science-origin-1.2.0.json")
ROLES = frozenset(profile_roles(PRODUCTION_PROFILE)) | {"spectroscopy"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write(path, record):
    with path.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def git(*args):
    return subprocess.check_output(["git", *args], text=True, timeout=60).strip()


def prepare(output, version, revision):
    # Workflows pass the verified collect outputs for a tag. A dispatch instead
    # audits main before the version-only commit exists.
    event, ref, head = (
        os.environ[key] for key in ("GITHUB_EVENT_NAME", "GITHUB_REF", "GITHUB_SHA")
    )
    require(os.environ["GITHUB_REPOSITORY"] == REPOSITORY, "unexpected repository")
    require(
        git("rev-parse", "HEAD") == head, "checkout does not match workflow revision"
    )
    require(
        not git("status", "--porcelain", "--untracked-files=all"), "checkout is dirty"
    )
    if event == "workflow_dispatch":
        require(
            ref == "refs/heads/main" and not (version or revision),
            "baseline dispatch requires main",
        )
        producer, candidate, mode = head, head, "baseline"
    else:
        require(
            event == "push" and version in ("1.1.0", "1.2.0"),
            "unsupported publication event or version",
        )
        require(
            ref == "refs/tags/v" + version and head == revision, "tag identity mismatch"
        )
        if version == "1.1.0":
            return "legacy", head
        parents = git("show", "-s", "--format=%P", head).split()
        require(len(parents) == 1, "version-only release must have one parent")
        producer, candidate, mode = parents[0], head, "tag"

    pin = read_receipt(ORIGIN_PIN)
    require(
        isinstance(pin, dict)
        and set(pin) == {"schema", "revision", "run_id"}
        and pin["schema"] == "kernelgecp.release-science-origin-pin.v1"
        and isinstance(pin["revision"], str)
        and type(pin["run_id"]) is int
        and pin["run_id"] > 0,
        "invalid canonical spectroscopy origin pin",
    )
    science = pin["revision"]
    equivalence = compare_release(".", candidate, baseline=producer)
    require(
        equivalence["version_transition"]
        == ("1.1.0 -> 1.2.0" if mode == "tag" else "1.1.0 -> 1.1.0"),
        "unexpected version transition",
    )
    scientific_equivalence = compare_science(".", science, producer)
    selection = fetch_science_runs(REPOSITORY, producer, profile=PRODUCTION_PROFILE)
    origin = fetch_spectroscopy_origin(REPOSITORY, science, pin["run_id"])

    output.mkdir()
    write(output / "selected-runs.json", selection)
    write(output / "spectroscopy-origin.json", origin)
    write(
        output / "plan.json",
        {
            "schema": PLAN_SCHEMA,
            "repository": REPOSITORY,
            "run_id": os.environ["GITHUB_RUN_ID"],
            "run_attempt": os.environ["GITHUB_RUN_ATTEMPT"],
            "producer_revision": producer,
            "candidate_revision": candidate,
            "science_revision": science,
            "mode": mode,
            "selection_sha256": digest(output / "selected-runs.json"),
            "origin_sha256": digest(output / "spectroscopy-origin.json"),
            "release_equivalence": equivalence,
            "scientific_equivalence": scientific_equivalence,
        },
    )
    return mode, producer


def aggregate(plan_dir, roles_dir, output, *, run_id, run_attempt, candidate):
    plan = read_receipt(plan_dir / "plan.json")
    selected_path = plan_dir / "selected-runs.json"
    origin_path = plan_dir / "spectroscopy-origin.json"
    selection, origin = (
        read_receipt(selected_path),
        origin_check(read_receipt(origin_path)),
    )
    validate_selection(selection)
    require(
        plan["schema"] == PLAN_SCHEMA
        and plan["repository"] == REPOSITORY
        and plan["run_id"] == run_id
        and plan["run_attempt"] == run_attempt
        and plan["candidate_revision"] == candidate
        and plan["mode"] in ("baseline", "tag")
        and selection.get("profile") == PRODUCTION_PROFILE
        and plan["producer_revision"] == selection["revision"]
        and plan["science_revision"] == origin["revision"]
        and origin["repository"] == REPOSITORY
        and origin["attempt"] == 1
        and plan["selection_sha256"] == digest(selected_path)
        and plan["origin_sha256"] == digest(origin_path),
        "wrong invocation, source, or authenticated plan digests",
    )
    require(
        not roles_dir.is_symlink()
        and {entry.name for entry in roles_dir.iterdir()} == ROLES,
        "missing or unknown replay roles",
    )
    receipts = {}
    for role in sorted(ROLES):
        folder = roles_dir / role
        require(
            not folder.is_symlink() and folder.is_dir(),
            f"{role}: invalid replay directory",
        )
        receipt_path = folder / "role-receipt.json"
        wrapper_path = folder / "wrapper-receipt.json"
        receipt, wrapper = read_receipt(receipt_path), read_receipt(wrapper_path)
        source = origin["revision"] if role == "spectroscopy" else selection["revision"]
        require(
            receipt["schema"] == "kernelgecp.science-role-replay.v1"
            and receipt["status"] == "passed"
            and receipt["role"] == role
            and receipt["producer_revision"] == plan["producer_revision"]
            and receipt["source_revision"] == source
            and receipt["selection_sha256"] == plan["selection_sha256"]
            and receipt["origin_sha256"] == plan["origin_sha256"]
            and receipt["wrapper_receipt_sha256"] == digest(wrapper_path)
            and type(receipt["command_exit"]) is int
            and receipt["command_exit"] == 0,
            f"{role}: invalid replay receipt binding",
        )
        validate_wrapper(role, wrapper, source, origin)
        expected = (
            [origin["artifact"]]
            if role == "spectroscopy"
            else selection["runs"][role]["artifacts"]
        )
        actual = [
            {key: record[key] for key in ("id", "name", "digest", "size_in_bytes")}
            for record in receipt["artifacts"]
        ]
        require(actual == expected, f"{role}: replay artifact identity mismatch")
        receipts[role] = {
            "role_receipt_sha256": digest(receipt_path),
            "wrapper_receipt_sha256": digest(wrapper_path),
            "inputs": origin if role == "spectroscopy" else selection["runs"][role],
        }
    write(
        output,
        {
            **plan,
            "schema": RECEIPT_SCHEMA,
            "status": "passed",
            "plan_sha256": digest(plan_dir / "plan.json"),
            "replays": receipts,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    pre = commands.add_parser("prepare")
    pre.add_argument("--output", required=True, type=Path)
    pre.add_argument("--version", default="")
    pre.add_argument("--revision", default="")
    post = commands.add_parser("aggregate")
    post.add_argument("--plan", required=True, type=Path)
    post.add_argument("--roles", required=True, type=Path)
    post.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            mode, producer = prepare(args.output, args.version, args.revision)
            with Path(os.environ["GITHUB_OUTPUT"]).open("a") as stream:
                stream.write(f"mode={mode}\nproducer={producer}\n")
        else:
            aggregate(
                args.plan,
                args.roles,
                args.output,
                run_id=os.environ["GITHUB_RUN_ID"],
                run_attempt=os.environ["GITHUB_RUN_ATTEMPT"],
                candidate=os.environ["GITHUB_SHA"],
            )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.SubprocessError,
    ) as error:
        parser.exit(1, f"release science gate rejected: {error}\n")


if __name__ == "__main__":
    main()
