"""Audit provider receipts and independently replay their saved DMFT campaigns.

This component does not decide overall scientific release eligibility.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath

from scripts.providers.verify_receipt import read_receipt, validate_receipt

ROOT = Path(__file__).resolve().parents[1]
LANES = {
    "provider-triqs-ubuntu-24.04": "triqs",
    "provider-triqs-macos-15": "triqs",
    "provider-gp-ubuntu-24.04": "gp",
    "provider-gp-macos-15": "gp",
    "provider-green-linux": "green",
    "provider-epw-linux": "epw",
}
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe(root, relative):
    require(isinstance(relative, str) and relative, "missing artifact path")
    p = PurePosixPath(relative)
    require(
        not p.is_absolute()
        and all(x not in ("", ".", "..") for x in relative.split("/"))
        and "\\" not in relative,
        "unsafe artifact path",
    )
    path = root
    for part in p.parts:
        path = path / part
        require(not path.is_symlink(), "symlink in artifact path")
    require(
        path.is_file() and path.stat().st_nlink == 1,
        f"missing or hardlinked artifact: {relative}",
    )
    return path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check_bindings(root, receipt, provider):
    bindings = []
    if provider == "green":
        bindings.extend(receipt["artifacts"].items())
    if provider == "epw":
        base = receipt["artifact_directory"]
        require(base == "epw-native", "unexpected EPW artifact directory")
        work = base + "/test-suite/epw_super/"
        for execution in receipt["native_executions"]:
            bindings.extend(
                (work + execution[k], execution[k + "_sha256"])
                for k in ("input", "log")
            )
        bindings.append(
            (receipt["calculation"]["gap_file"], receipt["calculation"]["gap_sha256"])
        )
        bindings.append(
            (
                base + "/kernelgecp-source-manifest.json",
                receipt["native_source_manifest_sha256"],
            )
        )
        bindings.extend(
            (base + "/" + name, sha)
            for name, sha in receipt["build_provenance_sha256"].items()
        )
    seen = set()
    for name, expected in bindings:
        require(name not in seen, "duplicate bound artifact path")
        seen.add(name)
        require(
            isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
            "invalid artifact hash",
        )
        require(digest(safe(root, name)) == expected, f"artifact hash mismatch: {name}")
    return {name: expected for name, expected in bindings}


@contextmanager
def source_snapshot(repository, revision):
    """A private, detached checkout: no untracked current-worktree imports."""
    with tempfile.TemporaryDirectory(prefix="gecp-provider-source-") as tmp:
        snapshot = Path(tmp) / "source"
        subprocess.run(
            [
                "git",
                "clone",
                "--quiet",
                "--no-checkout",
                "--shared",
                str(repository),
                str(snapshot),
            ],
            check=True,
            capture_output=True,
            timeout=60,
        )
        subprocess.run(
            ["git", "checkout", "--quiet", "--detach", revision],
            cwd=snapshot,
            check=True,
            capture_output=True,
            timeout=60,
        )
        head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=snapshot, text=True, timeout=60
        ).strip()
        require(head == revision, "snapshot revision mismatch")
        require(
            not subprocess.check_output(
                ["git", "status", "--porcelain"], cwd=snapshot, text=True, timeout=60
            ).strip(),
            "snapshot is dirty",
        )
        yield snapshot


def audit(
    artifacts_root,
    revision,
    triqs_python,
    output,
    *,
    runner=subprocess.run,
    repository=ROOT,
):
    result = {
        "schema": "kernelgecp.provider-payload-audit.v1",
        "status": "failed",
        "scope": (
            "provider runtime receipts and independent saved ED campaign audits "
            "only; not release eligibility"
        ),
        "source_revision": revision,
        "lanes": {},
        "limitations": [
            "EPW executable hashes refer to binaries in the provider image, "
            "not files in the uploaded artifact",
            "Authenticated artifact manifest binding is required upstream",
            "Artifact extraction must not create hardlinks",
        ],
    }
    output = Path(output).absolute()
    logs = output.parent / (output.name + ".logs")
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        logs.mkdir(exist_ok=True)
        require(re.fullmatch(r"[0-9a-f]{40}", revision) is not None, "invalid revision")
        head = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True
        ).strip()
        require(head == revision, "source HEAD differs from requested revision")
        dirty = subprocess.check_output(
            [
                "git",
                "status",
                "--porcelain",
                "--untracked-files=no",
                "--",
                "experiments",
                "kernelgecp",
                "scripts/providers",
                "environments/locks",
            ],
            cwd=repository,
            text=True,
        )
        require(not dirty.strip(), "scientific input sources are dirty")
        root = Path(artifacts_root).absolute()
        require(
            root.is_dir() and not any(p.is_symlink() for p in (root, *root.parents)),
            "invalid artifacts root",
        )
        for lane, provider in LANES.items():
            folder = root / lane
            require(
                not folder.is_symlink() and folder.is_dir(),
                f"missing or unsafe lane: {lane}",
            )
            receipt_path = safe(folder, provider + ".json")
            receipt = validate_receipt(read_receipt(receipt_path), provider, revision)
            bound = check_bindings(folder, receipt, provider)
            result["lanes"][lane] = {
                "receipt_valid": True,
                "receipt_sha256": digest(receipt_path),
                "bound_payload_sha256": bound,
            }
        for lane in (x for x, p in LANES.items() if p == "triqs"):
            folder = root / lane / "interacting-dmft"
            require(
                not folder.is_symlink() and folder.is_dir(), "missing DMFT campaign"
            )
            report = read_receipt(safe(folder, "report.json"))
            require(isinstance(report, dict), "DMFT report must be an object")
            require(
                report.get("schema") == "kernelgecp.interacting-dmft.v1"
                and report.get("pilot") is False
                and report.get("all_mandatory_passed") is True
                and report.get("source_revision") == revision,
                "noncanonical DMFT campaign",
            )
            require(
                len(report.get("references", {})) == 3
                and len(report.get("validation", [])) == 27
                and len(report.get("dense_validation", [])) == 9,
                "wrong DMFT cohort",
            )
            # Disallow symlinks anywhere in the auditor input tree.
            require(
                not any(p.is_symlink() for p in folder.rglob("*")),
                "symlink in DMFT campaign",
            )
            target = logs / (lane + "-dmft.json")
            logfile = logs / (lane + ".log")
            target.unlink(missing_ok=True)
            env = dict(os.environ)
            env.update({key: "1" for key in THREADS})
            env.pop("PYTHONPATH", None)
            try:
                with source_snapshot(repository, revision) as snapshot:
                    result["execution_source"] = {
                        "revision": revision,
                        "mode": "private detached Git checkout",
                        "isolated_python": True,
                    }
                    process = runner(
                        [
                            str(triqs_python),
                            "-I",
                            "-c",
                            (
                                "import sys,runpy; sys.path.insert(0,sys.argv.pop(1)); "
                                "runpy.run_module('experiments.audit_interacting_dmft',"
                                "run_name='__main__',alter_sys=True)"
                            ),
                            str(snapshot),
                            str(folder),
                            "--output",
                            str(target),
                            "--repository",
                            str(snapshot),
                        ],
                        cwd=snapshot,
                        env=env,
                        text=True,
                        capture_output=True,
                        check=False,
                    )
            except OSError as exc:
                logfile.write_text(str(exc) + "\n")
                raise ValueError(f"DMFT runtime failed: {exc}") from exc
            logfile.write_text((process.stdout or "") + "\n" + (process.stderr or ""))
            require(process.returncode == 0, f"DMFT auditor failed: {lane}")
            audited = read_receipt(target)
            require(isinstance(audited, dict), "DMFT audit must be an object")
            require(
                isinstance(audited.get("checks"), dict)
                and audited.get("schema") == "kernelgecp.interacting-dmft-audit.v1"
                and audited.get("passed") is True
                and audited.get("source_revision") == revision
                and audited.get("campaign_report_sha256")
                == digest(folder / "report.json")
                and audited.get("audited_records") == {"compressed": 27, "dense": 9}
                and all(
                    audited.get("checks", {}).get(x) is True
                    for x in (
                        "source_bindings",
                        "selection_and_environment",
                        "reference",
                        "records",
                        "repetitions",
                        "finite_pole_oracle",
                    )
                ),
                "invalid independent DMFT audit",
            )
            result["lanes"][lane]["independent_dmft"] = str(target)
        result["status"] = "passed"
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        subprocess.SubprocessError,
    ) as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    try:
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    except OSError as exc:
        print(f"cannot write provider audit receipt {output}: {exc}", file=sys.stderr)
        result["status"] = "failed"
        result["error"] = f"receipt write failed: {exc}"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["audit"])
    parser.add_argument("--artifacts-root", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--triqs-python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(
        args.artifacts_root, args.source_revision, args.triqs_python, args.output
    )
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
