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

# Exact paths only. These are release plumbing or release metadata,
# not science providers or production code. Their contents are not approved by
# this comparison; other release gates must validate and execute them.
RELEASE_ONLY = frozenset(
    {
        ".github/workflows/publish.yml",
        ".github/workflows/release-science.yml",
        ".github/workflows/spectroscopy-evidence.yml",
        "CITATION.cff",
        "scripts/release_spectroscopy_replay.py",
        "tests/test_release_spectroscopy_replay.py",
        "scripts/release_science_fetch.py",
        "scripts/release_science_runs.py",
        "scripts/release_science_payloads.py",
        "scripts/release_science_gate.py",
        "scripts/release_science_aggregate.py",
        "data/release-science-origin-1.2.0.json",
        "scripts/release_spectroscopy_v4.py",
        "scripts/release_cthyb_campaign_scope.py",
        "tests/test_release_cthyb_campaign_scope.py",
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
        "PACKAGE_README.md",
        "METHODS.md",
        "INTERACTING_DMFT.md",
        "CTHYB_DMFT.md",
        "MOLECULAR_GREEN.md",
        "EPW_INTEROPERABILITY.md",
        "COVARIANCE.md",
        "INSTALLATION_COMPATIBILITY.md",
        "MIGRATION.md",
        "PROVIDER_NOTICES.md",
        "RELEASE_NOTES.md",
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


CTHYB_WORKFLOW = ".github/workflows/cthyb-science.yml"


def workflow_outside_scope(repo, oid):
    """Read a committed blob and isolate the uniquely delimited scope job."""
    content = git(repo, "cat-file", "blob", oid)
    roots = list(re.finditer(rb"(?m)^jobs:\r?$", content))
    if len(roots) != 1:
        raise ValueError("malformed cthyb jobs delimiter")
    job_region = content[roots[0].end() :]
    # Reject alternate or malformed two-space job keys instead of trying to
    # interpret YAML aliases, quoted keys, or duplicate mappings.
    headers = list(re.finditer(rb"(?m)^  [^ \t#\r\n][^\r\n]*", job_region))
    if any(
        not re.fullmatch(rb"  [A-Za-z_][A-Za-z0-9_-]*:", m.group()) for m in headers
    ):
        raise ValueError("malformed cthyb job header")
    scopes = [i for i, match in enumerate(headers) if match.group() == b"  scope:"]
    calibrations = [
        i for i, match in enumerate(headers) if match.group() == b"  calibration:"
    ]
    if len(scopes) != 1 or len(calibrations) != 1 or calibrations[0] != scopes[0] + 1:
        raise ValueError("malformed cthyb scope job delimiters")
    offset = roots[0].end()
    scope_body = job_region[headers[scopes[0]].end() : headers[calibrations[0]].start()]
    if re.search(rb"(?m)^ {0,3}[^ \t#\r\n]", scope_body):
        raise ValueError("cthyb scope escapes its job")
    return content[: offset + headers[scopes[0]].start()], content[
        offset + headers[calibrations[0]].start() :
    ]


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
        if source[path] == target[path]:
            continue
        if (
            path != CTHYB_WORKFLOW
            or source[path][0] != target[path][0]
            or workflow_outside_scope(repo, source[path][1])
            != workflow_outside_scope(repo, target[path][1])
        ):
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
