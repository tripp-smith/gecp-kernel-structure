"""Select the original CT-HYB preflight for a narrowly scoped release bump."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from experiments.replay_cthyb_acceptance import campaign_needed
from scripts import release_source_equivalence as equivalence


def eligible(repo: Path, event: str, candidate: str, base: str) -> bool:
    """Fail closed to ordinary candidate preflight for any unverified identity."""
    if event != "pull_request":
        return False
    try:
        equivalence.commit(str(repo), base)
        equivalence.commit(str(repo), candidate)
        if (
            equivalence.git(str(repo), "rev-parse", "HEAD").decode().strip()
            != candidate
        ):
            return False
        if equivalence.git(str(repo), "status", "--porcelain", "--untracked-files=all"):
            return False
        parents = (
            equivalence.git(str(repo), "rev-list", "--parents", "-n", "1", candidate)
            .decode()
            .split()
        )
        if parents != [candidate, base]:
            return False
        paths = (
            equivalence.git(
                str(repo),
                "diff",
                "--no-renames",
                "--name-only",
                "-z",
                base,
                candidate,
            )
            .decode()
            .split("\0")
        )
        if set(paths[:-1]) != equivalence.VERSION_FILES or paths[-1] != "":
            return False
        result = equivalence.compare(str(repo), candidate, baseline=base)
        return result["version_transition"] == "1.1.0 -> 1.2.0"
    except (ValueError, KeyError, UnicodeError, OSError, subprocess.CalledProcessError):
        return False


def decide(repo: Path, event: str, candidate: str, base: str) -> bool:
    if event == "workflow_dispatch":
        return campaign_needed(repo, force=True)
    if not eligible(repo, event, candidate, base):
        return campaign_needed(repo)
    # A detached clone provides both the exact baseline files and git history that
    # the original preflight checks, without changing the candidate checkout.
    with tempfile.TemporaryDirectory(prefix="cthyb-scope-base-") as temp:
        subprocess.run(
            ["git", "clone", "--quiet", "--no-local", "--no-checkout", str(repo), temp],
            check=True,
        )
        subprocess.run(
            ["git", "-C", temp, "checkout", "--quiet", "--detach", base],
            check=True,
        )
        return campaign_needed(Path(temp))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--base", default="")
    args = parser.parse_args()
    print(
        json.dumps(
            {
                "requires_campaign": decide(
                    Path.cwd(), args.event, args.candidate, args.base
                )
            }
        )
    )


if __name__ == "__main__":
    main()
