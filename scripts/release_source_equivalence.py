#!/usr/bin/env python3
"""Compare committed source only; this does not establish scientific acceptance.

An upstream eligibility gate must derive the baseline from a fresh authenticated
scientific run or receipt and independently audit its payloads. Supplying a
baseline here does not establish its trustworthiness or approve a release.
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tomllib

BASELINE = "feaac669eb49d700f0059f5be755a9e4403a1719"
OLD, NEW = b"1.1.0", b"1.2.0"
VERSION_FILES = {
    "kernelgecp/__init__.py",
    "pyproject.toml",
    "uv.lock",
    "native/Cargo.toml",
    "native/Cargo.lock",
}
# Compare all tracked files in these trees, including workflows, fixtures,
# scientific data and verification scripts, regardless of filename extension.
# Root prose and release documents are outside this scientific-source check;
# their release binding and installed-package acceptance are separate gates.
ROOT_INPUTS = {
    ".dockerignore",
    ".python-version",
    "pyproject.toml",
    "uv.lock",
    "rust-toolchain.toml",
    ".gitattributes",
}


def git(repo, *args):
    return subprocess.check_output(["git", "-C", repo, *args], stderr=subprocess.PIPE)


def commit(repo, sha):
    if not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise ValueError("identity must be a full lowercase SHA-1 commit ID")
    if git(repo, "rev-parse", "--verify", sha + "^{commit}").decode().strip() != sha:
        raise ValueError("identity must name the commit itself")
    return sha


def critical(path):
    return (
        path in ROOT_INPUTS
        or path.endswith((".lock", ".lock.json"))
        or path.startswith(
            (
                ".github/",
                "data/",
                "environments/",
                "examples/",
                "experiments/",
                "kernelgecp/",
                "native/",
                "scripts/",
                "tests/",
            )
        )
    )


def tree(repo, sha):
    entries = {}
    for item in git(repo, "ls-tree", "-rz", "--full-tree", sha).split(b"\0"):
        if not item:
            continue
        meta, raw = item.split(b"\t", 1)
        mode, typ, oid = meta.decode().split()
        path = raw.decode("utf-8", "strict")
        if critical(path):
            if typ != "blob" or mode not in ("100644", "100755"):
                raise ValueError(f"non-regular tracked source: {path}")
            entries[path] = (mode, oid)
    return entries


def blob(repo, oid):
    return git(repo, "cat-file", "blob", oid)


def normalized(path, data, expected):
    if path == "kernelgecp/__init__.py":
        pattern = rb'(?m)^__version__ = "([^"\r\n]+)"$'
        matches = list(re.finditer(pattern, data))
        if len(matches) != 1 or matches[0].group(1) != expected:
            raise ValueError("invalid unique __version__ assignment")
        a, b = matches[0].span(1)
        return data[:a] + b"<VERSION>" + data[b:]
    doc = tomllib.loads(data.decode("utf-8"))
    if path == "pyproject.toml":
        section, name, table = b"[project]", "kernelgecp", doc["project"]
    elif path == "native/Cargo.toml":
        section, name, table = b"[package]", "kernelgecp-native", doc["package"]
    else:
        section = b"[[package]]"
        name = "kernelgecp" if path == "uv.lock" else "kernelgecp-native"
        matches = [p for p in doc["package"] if p.get("name") == name]
        if len(matches) != 1:
            raise ValueError(f"non-unique package: {path}")
        table = matches[0]
    if table.get("name") != name or table.get("version") != expected.decode():
        raise ValueError(f"incorrect package identity/version: {path}")
    # Section boundaries include all TOML tables so unrelated nested fields cannot
    # be mistaken for the package's own direct version assignment.
    headers = list(re.finditer(rb"(?m)^\[\[?[^\r\n]+\]\]?\s*$", data))
    candidates = []
    for i, h in enumerate(headers):
        if h.group().strip() != section:
            continue
        end = headers[i + 1].start() if i + 1 < len(headers) else len(data)
        block = data[h.end() : end]
        names = list(re.finditer(rb'(?m)^name = "([^"\r\n]+)"$', block))
        if len(names) == 1 and names[0].group(1).decode() == name:
            versions = list(re.finditer(rb'(?m)^version = "([^"\r\n]+)"$', block))
            if len(versions) != 1 or versions[0].group(1) != expected:
                raise ValueError(f"invalid version stanza: {path}")
            a, b = versions[0].span(1)
            candidates.append((h.end() + a, h.end() + b))
    if len(candidates) != 1:
        raise ValueError(f"ambiguous stanza: {path}")
    a, b = candidates[0]
    return data[:a] + b"<VERSION>" + data[b:]


def compare(repo, candidate, baseline=None):
    baseline = BASELINE if baseline is None else baseline
    commit(repo, baseline)
    commit(repo, candidate)
    subprocess.run(
        ["git", "-C", repo, "merge-base", "--is-ancestor", baseline, candidate],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    left, right = tree(repo, baseline), tree(repo, candidate)
    additions = right.keys() - left.keys()
    if additions:
        raise ValueError(f"unexpected critical additions: {sorted(additions)}")
    if not left.keys() >= VERSION_FILES or not right.keys() >= VERSION_FILES:
        raise ValueError("missing version file")
    if left.keys() != right.keys():
        raise ValueError(
            "critical tracked path mismatch: "
            f"removed={sorted(left.keys() - right.keys())}, "
            f"added={sorted(right.keys() - left.keys())}"
        )
    version = None
    for path in sorted(left):
        mode_a, oid_a = left[path]
        mode_b, oid_b = right[path]
        if mode_a != mode_b:
            raise ValueError(f"mode changed: {path}")
        if path in VERSION_FILES:
            source = blob(repo, oid_a)
            target = blob(repo, oid_b)
            canonical = normalized(path, source, OLD)
            if target == source:
                current = OLD
            elif canonical == normalized(path, target, NEW):
                current = NEW
            else:
                raise ValueError(f"non-version change: {path}")
            if version is not None and current != version:
                raise ValueError("version fields do not share a single transition")
            version = current
        elif oid_a != oid_b:
            raise ValueError(f"critical blob changed: {path}")

    def inventory_digest(entries):
        inventory = "\n".join(
            f"{p}\t{entries[p][0]}\t{entries[p][1]}" for p in sorted(entries)
        ).encode()
        return hashlib.sha256(inventory).hexdigest()

    return {
        "baseline_commit": baseline,
        "candidate_commit": candidate,
        "ancestor_verified": True,
        "protected_paths": len(left),
        "baseline_inventory_sha256": inventory_digest(left),
        "candidate_inventory_sha256": inventory_digest(right),
        "candidate_protected_entries": [
            {"path": p, "mode": right[p][0], "blob_id": right[p][1]}
            for p in sorted(right)
        ],
        "version_transition": "1.1.0 -> " + version.decode(),
        "allowed_additions": sorted(additions),
        "scope": (
            "committed source equivalence only; "
            "not scientific or installed-package acceptance or baseline trust; "
            "the upstream eligibility gate must derive the baseline from a fresh "
            "authenticated scientific run or receipt and independently audit payloads"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument(
        "--candidate", required=True, help="caller-supplied full expected commit SHA"
    )
    parser.add_argument(
        "--baseline",
        default=None,
        help=(
            "full baseline commit SHA; caller selection does not confer trust "
            "or approval"
        ),
    )
    args = parser.parse_args()
    try:
        print(
            json.dumps(
                compare(args.repo, args.candidate, baseline=args.baseline),
                sort_keys=True,
            )
        )
    except (
        ValueError,
        KeyError,
        UnicodeError,
        tomllib.TOMLDecodeError,
        subprocess.CalledProcessError,
    ) as e:
        print(f"equivalence rejected: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
