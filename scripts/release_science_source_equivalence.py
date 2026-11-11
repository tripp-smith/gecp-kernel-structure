#!/usr/bin/env python3
"""Compare frozen science source S with committed producer source P.

This is an input identity check, not artifact authentication, replay, scientific
acceptance, or the separate P-to-K release version comparison. The caller must
authenticate both identities and must not infer acceptance from this check.
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys

# Exact paths only. These are release plumbing,
# not science providers or production code. Their contents are not approved by
# this comparison; other release gates must validate and execute them.
RELEASE_ONLY = frozenset(
    {
        ".github/workflows/publish.yml",
        ".github/workflows/release-science.yml",
        "scripts/release_science_fetch.py",
        "scripts/release_science_runs.py",
        "scripts/release_science_payloads.py",
        "scripts/release_science_gate.py",
        "scripts/release_science_aggregate.py",
        "data/release-science-origin-1.2.0.json",
        "scripts/release_spectroscopy_v4.py",
        "scripts/release_science_source_equivalence.py",
        "tests/test_release_science_fetch.py",
        "tests/test_release_science_runs.py",
        "tests/test_release_science_payloads.py",
        "tests/test_release_science_gate.py",
        "tests/test_release_science_aggregate.py",
        "tests/test_release_spectroscopy_origin.py",
        "tests/test_release_spectroscopy_v4.py",
        "tests/test_release_science_source_equivalence.py",
    }
)

# Root prose can evolve independently of computational inputs. No directory or
# filename-pattern exclusions: settings, source, workflows, evidence, fixtures,
# lockfiles and unknown tracked paths are all protected, including their mode.
ROOT_PROSE = frozenset(
    {
        "APPLICATION.md",
        "FINDINGS.md",
        "PRODUCTIZATION.md",
        "README.md",
        "RESEARCH.md",
        "SPEC.md",
        "SPECTROSCOPY.md",
        "FINAL_HANDOFF.md",
        "MATHLIB.md",
        "CHANGELOG.md",
        "CONTRIBUTING.md",
    }
)


def git(repo, *args):
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], stderr=subprocess.PIPE
    )


def commit(repo, identity):
    if not re.fullmatch(r"[0-9a-f]{40}", identity):
        raise ValueError("identity must be a full lowercase SHA-1 commit ID")
    if (
        git(repo, "rev-parse", "--verify", identity + "^{commit}").decode().strip()
        != identity
    ):
        raise ValueError("identity must name the commit itself")
    return identity


def tree(repo, identity):
    entries = {}
    for item in git(repo, "ls-tree", "-rz", "--full-tree", identity).split(b"\0"):
        if not item:
            continue
        meta, raw = item.split(b"\t", 1)
        mode, typ, oid = meta.decode("ascii").split()
        path = raw.decode("utf-8", "strict")
        if typ != "blob" or mode not in ("100644", "100755"):
            raise ValueError(f"non-regular tracked input: {path}")
        entries[path] = (mode, oid)
    return entries


def inventory_digest(entries):
    # NUL-delimited fields avoid ambiguous records for unusual tracked names.
    encoded = b"".join(
        p.encode("utf-8") + b"\0" + mode.encode() + b"\0" + oid.encode() + b"\0"
        for p, (mode, oid) in sorted(entries.items())
    )
    return hashlib.sha256(encoded).hexdigest()


def compare(repo, science, producer):
    commit(repo, science)
    commit(repo, producer)
    subprocess.run(
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", science, producer],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    source, target = tree(repo, science), tree(repo, producer)
    removed = source.keys() - target.keys() - ROOT_PROSE
    added = target.keys() - source.keys()
    forbidden = added - RELEASE_ONLY - ROOT_PROSE
    if removed or forbidden:
        raise ValueError(
            f"tracked input mismatch: removed={sorted(removed)}, "
            f"added={sorted(forbidden)}"
        )
    for path in sorted((source.keys() & target.keys()) - ROOT_PROSE - RELEASE_ONLY):
        if source[path] != target[path]:
            raise ValueError(f"tracked input mode or blob changed: {path}")
    return {
        "science_commit": science,
        "producer_commit": producer,
        "protected_paths": len(source),
        "science_inventory_sha256": inventory_digest(source),
        "producer_inventory_sha256": inventory_digest(target),
        "ancestor_verified": True,
        "allowed_additions": sorted(added & RELEASE_ONLY),
        "changed_release_only": sorted(
            p
            for p in source.keys() & target.keys() & RELEASE_ONLY
            if source[p] != target[p]
        ),
        "scope": (
            "committed input identity only; no artifact authentication or "
            "scientific acceptance"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument(
        "--science", required=True, help="caller-authenticated full science commit"
    )
    parser.add_argument(
        "--producer", required=True, help="authenticated full producer commit"
    )
    args = parser.parse_args()
    try:
        print(
            json.dumps(compare(args.repo, args.science, args.producer), sort_keys=True)
        )
    except (ValueError, UnicodeError, subprocess.CalledProcessError) as exc:
        print(f"science source equivalence rejected: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
