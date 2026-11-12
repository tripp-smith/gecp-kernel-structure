"""Replay caller-authenticated v4 evidence without sampling or origin checks."""

import argparse
import json
import os
import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path

from scripts.release_epw_audits import as_text, bounded_runner
from scripts.release_epw_native import checked_path, inventory, read_json, safe
from scripts.release_provider_audits import digest, require, source_snapshot
from scripts.release_science_payloads import (
    CHUNK,
    preflight,
)
from scripts.release_science_payloads import (
    inventory as zip_inventory,
)

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = Path("/home/runner/work/gecp-kernel-structure/gecp-kernel-structure")
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
GATES = {
    "ee_source_and_raw_chains",
    "frozen_protocol_and_hashes",
    "covariance_and_projection",
    "native_trials_and_quadrature",
    "heldout_rms",
    "bootstrap_and_defaults",
    "repeat_arrays",
    "fresh_process_receipts",
}
TIMEOUT = 5400
BOOTSTRAP = """
import importlib, importlib.metadata as md, json, os, platform, runpy, sys

assert sys.platform == "linux" and platform.machine() in ("x86_64", "AMD64")
names = (
    "triqs",
    "triqs_cthyb",
    "triqs_maxent",
    "pydlr",
    "sparse_ir",
    "numpy",
    "scipy",
    "pylibsparseir",
    "mpmath",
)
mods = {n: importlib.import_module(n) for n in names}
versions = {
    n: md.version(n)
    for n in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir")
}
versions.update(
    {
        n: str(importlib.import_module(n + ".version").version)
        for n in ("triqs", "triqs_cthyb", "triqs_maxent")
    }
)
expected = {
    "triqs": "4.0.2",
    "triqs_cthyb": "4.0.1",
    "triqs_maxent": "4.0.1",
    "numpy": "2.2.6",
    "scipy": "1.18.1",
    "mpmath": "1.3.0",
    "pydlr": "1.0.1",
    "sparse-ir": "2.1.6",
    "pylibsparseir": "0.12.0",
}
threads = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
assert (
    platform.python_version() == "3.12.16"
    and versions == expected
    and all(os.environ.get(n) == "1" for n in threads)
)
prefix = os.path.realpath(sys.prefix)
paths = {n: os.path.realpath(m.__file__) for n, m in mods.items()}
assert all(os.path.commonpath((prefix, p)) == prefix for p in paths.values())
with open(sys.argv.pop(1), "x") as f:
    json.dump(
        {
            "python": sys.version,
            "executable": sys.executable,
            "prefix": prefix,
            "versions": versions,
            "modules": paths,
            "threads": {n: os.environ[n] for n in threads},
        },
        f,
        sort_keys=True,
    )
source = sys.argv.pop(1)
module = sys.argv.pop(1)
sys.path.insert(0, source)
runpy.run_module(module, run_name="__main__", alter_sys=True)
"""


def extract(archive, destination):
    count = preflight(archive)
    with zipfile.ZipFile(archive) as z:
        entries, total = zip_inventory(z)
        require(len(entries) == count, "ZIP entry count mismatch")
        size = 0
        for entry in entries:
            target = destination.joinpath(*entry.filename.rstrip("/").split("/"))
            if entry.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            n = 0
            with z.open(entry) as src, target.open("xb") as dst:
                while chunk := src.read(CHUNK):
                    n += len(chunk)
                    size += len(chunk)
                    require(n <= entry.file_size and size <= total, "ZIP size mismatch")
                    dst.write(chunk)
            require(n == entry.file_size, "ZIP member size mismatch")
        require(size == total, "ZIP total mismatch")


def shard_jobs(ordinal):
    """The 6caf85f v4 campaign's declared order and controls."""
    jobs = []
    for ci, case in enumerate(
        ("offgrid_atoms_independent", "offgrid_atoms_correlated", "smooth_correlated")
    ):
        for method in ("raw", "kernelgecp", "pydlr", "sparse_ir"):
            for model in ("flat", "broad_gaussian"):
                for repeat in range(3):
                    job = {
                        "case": case,
                        "method": method,
                        "default_model": model,
                        "seed": 5001 + ci,
                        "eps": 1e-8,
                        "rcond": 1e-10,
                        "precision_bits": 53,
                        "n_omega": 401,
                        "n_alpha": 33,
                        "minimizer_mu0": 1e-8,
                        "query_order": 80,
                        "bootstrap": 8 if repeat == 0 else 0,
                        "pilot": False,
                        "repeat": repeat,
                        "profile": "v4",
                    }
                    name = f"{case}-{method}-{model}-r{repeat}"
                    jobs.append((name, job))
    return dict(jobs[ordinal * 12 : (ordinal + 1) * 12])


def check_shards(directory, manifest, freeze_hash):
    refs = manifest.get("shard_receipts")
    children = manifest.get("children")
    require(
        isinstance(refs, dict)
        and set(refs) == {str(i) for i in range(6)}
        and isinstance(children, dict),
        "shard receipt inventory",
    )
    seen = set()
    for ordinal in range(6):
        ref = refs[str(ordinal)]
        path = f"shard-receipts/shard-{ordinal}.json"
        require(
            isinstance(ref, dict)
            and set(ref) == {"path", "sha256"}
            and ref["path"] == path,
            "shard receipt path",
        )
        receipt_path = safe(directory, path)
        require(digest(receipt_path) == ref["sha256"], "shard receipt hash")
        r = read_json(receipt_path)
        expected = shard_jobs(ordinal)
        part = r.get("children")
        require(
            r.get("schema") == "kernelgecp.spectroscopy-campaign-shard.v4"
            and r.get("kind") == "campaign"
            and type(r.get("ordinal")) is int
            and r["ordinal"] == ordinal
            and r.get("status") == "passed"
            and r.get("controls") == manifest.get("controls")
            and r.get("provenance") == manifest.get("provenance")
            and r.get("freeze_sha256") == freeze_hash
            and r.get("repeat_comparison", {}).get("status") == "passed"
            and isinstance(part, dict)
            and set(part) == set(expected),
            "invalid shard receipt or partition",
        )
        require(not (seen & set(part)), "overlapping shard children")
        seen.update(part)
        filenames = set()
        for name, job in expected.items():
            record = part[name]
            require(
                isinstance(record, dict)
                and record.get("job") == job
                and record.get("status") == "passed"
                and type(record.get("exit_code")) is int
                and record.get("exit_code") == 0
                and record == children.get(name),
                "shard child mismatch",
            )
            filenames.update(
                name + suffix
                for suffix in (".json", ".npz", ".log", "-dataset.json", "-dataset.npz")
            )
            if job["method"] == "kernelgecp":
                filenames.add(name + "-gecp-fit.npz")
        inv = r.get("inventory")
        require(isinstance(inv, dict) and set(inv) == filenames, "shard file inventory")
        for name, expected_hash in inv.items():
            require(digest(safe(directory, name)) == expected_hash, "shard file hash")
    require(seen == set(children), "assembled shard union")


def evidence_checks(e, source, revision):
    def f(name):
        return safe(e, name)

    def j(name):
        return read_json(f(name))

    def sha(name):
        return digest(f(name))

    freeze = j("freeze.json")
    require(
        freeze.get("schema") == "kernelgecp.spectroscopy-freeze.v4"
        and freeze.get("status") == "passed"
        and freeze.get("provenance", {}).get("source_revision") == revision,
        "freeze profile or source",
    )
    for folder, audit, kind, count in [
        ("calibration", "calibration-audit.json", "calibration", 24),
        ("campaign", "campaign-audit.json", "campaign", 72),
    ]:
        m = j(folder + "/campaign.json")
        a = j(audit)
        require(
            m.get("schema") == "kernelgecp.spectroscopy-campaign.v4"
            and m.get("status") == "passed"
            and m.get("kind") == kind
            and m.get("provenance", {}).get("source_revision") == revision
            and isinstance(m.get("children"), dict)
            and len(m["children"]) == count,
            "synthetic inventory",
        )
        synthetic(a, kind, count, sha(folder + "/campaign.json"), sha("freeze.json"))
        if kind == "campaign":
            check_shards(e / folder, m, sha("freeze.json"))
    cal = "native-calibration/calibration.json"
    final = "native-final/final.json"
    for name, kind in [(cal, "calibration"), (final, "final")]:
        m = j(name)
        require(
            m.get("schema") == "kernelgecp.cthyb-spectroscopy-" + kind + ".v4"
            and m.get("snapshot", {}).get("revision") == revision,
            "measured identity",
        )
    ee_binding = digest(safe(source, "data/cthyb/audit-binding-e83d195.json"))
    for a, kind, count, campaign, binding in [
        ("native-calibration-audit.json", "calibration", 8, cal, ee_binding),
        ("native-final-audit.json", "final", 24, final, sha("fixed-bath-binding.json")),
    ]:
        measured(
            j(a),
            kind,
            count,
            sha(campaign),
            sha(cal),
            binding,
            revision,
            digest(safe(source, "experiments/audit_cthyb_spectroscopy.py")),
        )
    r, b = j("fixed-bath-audit.json"), j("fixed-bath-binding.json")
    raw(r, b, e, source, revision, sha(cal), sha("native-calibration-audit.json"))
    manifest = j("fixed-bath/manifest.json")
    require(
        manifest.get("provenance", {}).get("revision") == revision
        and manifest.get("provenance", {}).get("sources") == sources(source),
        "raw source revision or inventory",
    )
    return {
        "synthetic_calibration": 24,
        "synthetic_final": 72,
        "measured_calibration": 8,
        "measured_final": 24,
        "raw_chains": 256,
    }


def synthetic(a, kind, count, manifest, freeze):
    require(
        a.get("status") == "passed"
        and a.get("kind") == kind
        and type(a.get("cases")) is int
        and a["cases"] == count
        and a.get("independently_audited_cases") == 24
        and a.get("manifest_sha256") == manifest
        and a.get("freeze_sha256") == freeze,
        "synthetic audit",
    )


def measured(a, kind, count, campaign, calibration, binding, revision, auditor):
    require(
        a.get("schema") == "kernelgecp.cthyb-spectroscopy-independent-audit.v4"
        and a.get("status") == "passed"
        and a.get("kind") == kind
        and type(a.get("cases")) is int
        and a["cases"] == count
        and a.get("campaign_sha256") == campaign
        and a.get("calibration_sha256") == calibration
        and a.get("binding_sha256") == binding
        and a.get("auditor_revision") == revision
        and a.get("auditor_sha256") == auditor
        and isinstance(a.get("gates"), dict)
        and set(a["gates"]) == GATES
        and all(v is True for v in a["gates"].values()),
        "measured audit",
    )


def raw(r, b, e, source, revision, cal_hash, cal_audit_hash):
    require(
        r.get("schema") == "kernelgecp.cthyb-fixed-bath-audit.v4"
        and b.get("schema") == "kernelgecp.cthyb-fixed-bath-audit-binding.v4"
        and r.get("status") == b.get("status") == "passed"
        and (r.get("chain_count"), r.get("train_count"), r.get("heldout_count"))
        == (256, 128, 128),
        "raw audit",
    )
    require(
        r.get("audit_revision") == revision
        and r.get("audit_source_sha256")
        == digest(safe(source, "experiments/audit_cthyb_fixed_bath.py")),
        "raw auditor identity",
    )
    require(
        b.get("receipt_sha256") == digest(safe(e, "fixed-bath-audit.json")),
        "raw receipt hash",
    )
    for k in (
        "audit_revision",
        "audit_source_sha256",
        "calibration",
        "manifest_sha256",
        "parent_binding_sha256",
        "result_sha256",
    ):
        require(b.get(k) == r.get(k), "raw binding " + k)
    require(
        r.get("manifest_sha256") == digest(safe(e, "fixed-bath/manifest.json"))
        and r.get("result_sha256") == digest(safe(e, "fixed-bath/result.json"))
        and r.get("parent_binding_sha256")
        == digest(safe(source, "data/cthyb/linux-audit-binding-7200ec0.json")),
        "raw inputs",
    )
    require(
        r.get("calibration", {}).get("receipt_sha256") == cal_hash
        and r["calibration"].get("audit_sha256") == cal_audit_hash,
        "raw calibration",
    )
    expected = {
        f"chains/chain-{i:03d}{ext}"
        for i in range(256)
        for ext in (".json", ".npz", ".timing.json")
    }
    inv = b.get("inventory")
    require(isinstance(inv, dict) and set(inv) == expected, "raw chain inventory")
    for name, h in inv.items():
        require(digest(safe(e, "fixed-bath/" + name)) == h, "raw chain hash")


def sources(source):
    names = subprocess.check_output(
        ["git", "ls-files"], cwd=source, text=True, timeout=60
    ).splitlines()
    return {
        p: digest(safe(source, p))
        for p in names
        if p.endswith((".py", ".toml", ".lock", ".txt"))
    }


def audit(
    archive,
    archive_sha256,
    archive_size,
    run_id,
    revision,
    python,
    output,
    *,
    repository=ROOT,
    workspace=WORKSPACE,
    runner=bounded_runner,
):
    archive = checked_path(archive)
    output = checked_path(output)
    workspace = checked_path(workspace)
    require(
        re.fullmatch("[0-9a-f]{64}", archive_sha256) is not None
        and re.fullmatch("[0-9a-f]{40}", revision) is not None
        and re.fullmatch("[1-9][0-9]*", str(run_id)) is not None
        and type(archive_size) is int
        and archive_size > 0,
        "invalid authenticated identity",
    )
    require(
        archive.is_file()
        and archive.stat().st_nlink == 1
        and archive.stat().st_size == archive_size
        and digest(archive) == archive_sha256,
        "archive authentication",
    )
    require(
        archive.name == f"maxent-fixed-bath-final-{run_id}-1.zip", "artifact identity"
    )
    require(
        workspace == WORKSPACE
        and workspace.is_dir()
        and not output.exists()
        and output.parent.is_dir()
        and workspace not in output.parents
        and archive not in output.parents,
        "unsafe workspace/output",
    )
    evidence = workspace / "artifacts"
    require(
        not evidence.exists() and not evidence.is_symlink(),
        "evidence path must be absent",
    )
    output.mkdir()
    receipt = {
        "schema": "kernelgecp.spectroscopy-v4-saved-replay.v1",
        "status": "failed",
        "source_revision": revision,
        "archive_sha256": archive_sha256,
        "archive_size": archive_size,
        "run_id": str(run_id),
        "scope": (
            "caller-authenticated saved evidence only; origin not verified; "
            "not release eligibility"
        ),
    }
    before = None
    try:
        evidence.mkdir()
        extract(archive, evidence)
        before = inventory(evidence)
        with source_snapshot(repository, revision) as source:
            source_before = sources(source)
            counts = evidence_checks(evidence, source, revision)
            env = dict(os.environ)
            env.update(dict.fromkeys(THREADS, "1"))
            env.pop("PYTHONPATH", None)
            ee = workspace / "data/cthyb"
            for binding in (
                "independent-audit-e83d195.json",
                "audit-binding-e83d195.json",
                "linux-audit-binding-7200ec0.json",
            ):
                require(
                    digest(safe(ee, binding))
                    == digest(safe(source, "data/cthyb/" + binding)),
                    "workspace data differs from frozen source",
                )
            cal = evidence / "native-calibration/calibration.json"
            ca = evidence / "native-calibration-audit.json"
            fb = evidence / "fixed-bath"
            freeze = evidence / "freeze.json"
            jobs = [
                (
                    "synthetic-cal",
                    "experiments.audit_spectroscopy",
                    [
                        "--profile",
                        "v4",
                        "--manifest",
                        evidence / "calibration/campaign.json",
                        "--freeze",
                        freeze,
                        "--output",
                        output / "synthetic-cal.json",
                    ],
                ),
                (
                    "synthetic-final",
                    "experiments.audit_spectroscopy",
                    [
                        "--profile",
                        "v4",
                        "--manifest",
                        evidence / "campaign/campaign.json",
                        "--freeze",
                        freeze,
                        "--output",
                        output / "synthetic-final.json",
                    ],
                ),
                (
                    "measured-cal",
                    "experiments.audit_cthyb_spectroscopy",
                    [
                        "--profile",
                        "v4",
                        cal,
                        "--ee-root",
                        evidence / "cthyb-input/gecp-ee-odd-9d28053",
                        "--ee-audit",
                        ee / "independent-audit-e83d195.json",
                        "--ee-binding",
                        ee / "audit-binding-e83d195.json",
                        "--output",
                        output / "measured-cal.json",
                    ],
                ),
                (
                    "raw",
                    "experiments.audit_cthyb_fixed_bath",
                    [
                        "--profile",
                        "v4",
                        "--folder",
                        fb,
                        "--parent-dir",
                        evidence / "parent",
                        "--parent-binding",
                        ee / "linux-audit-binding-7200ec0.json",
                        "--receipt",
                        output / "raw.json",
                        "--binding",
                        output / "raw-binding.json",
                        "--calibration",
                        cal,
                        "--calibration-audit",
                        ca,
                    ],
                ),
                (
                    "measured-final",
                    "experiments.audit_cthyb_spectroscopy",
                    [
                        "--profile",
                        "v4",
                        evidence / "native-final/final.json",
                        "--ee-root",
                        fb,
                        "--ee-audit",
                        evidence / "fixed-bath-audit.json",
                        "--ee-binding",
                        evidence / "fixed-bath-binding.json",
                        "--calibration",
                        cal,
                        "--calibration-audit",
                        ca,
                        "--output",
                        output / "measured-final.json",
                    ],
                ),
            ]
            deadline = time.monotonic() + TIMEOUT
            for name, module, args in jobs:
                cmd = [
                    str(Path(python).absolute()),
                    "-I",
                    "-c",
                    BOOTSTRAP,
                    str(output / (name + "-runtime.json")),
                    str(source),
                    module,
                    *map(str, args),
                ]
                try:
                    remaining = deadline - time.monotonic()
                    require(remaining > 0, "audit deadline")
                    proc = runner(
                        cmd,
                        cwd=source,
                        env=env,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=remaining,
                    )
                    (output / (name + ".log")).write_text(
                        as_text(proc.stdout) + "\n" + as_text(proc.stderr)
                    )
                    require(proc.returncode == 0, "auditor failed: " + name)
                except subprocess.TimeoutExpired as error:
                    (output / (name + ".log")).write_text(
                        as_text(error.stdout)
                        + "\n"
                        + as_text(error.stderr)
                        + "\nTimeoutExpired\n"
                    )
                    raise ValueError("audit timed out") from error
                read_json(output / (name + "-runtime.json"))
            for name, kind, count, manifest in [
                ("synthetic-cal", "calibration", 24, "calibration/campaign.json"),
                ("synthetic-final", "campaign", 72, "campaign/campaign.json"),
            ]:
                synthetic(
                    read_json(output / (name + ".json")),
                    kind,
                    count,
                    digest(safe(evidence, manifest)),
                    digest(freeze),
                )
            auditor = digest(safe(source, "experiments/audit_cthyb_spectroscopy.py"))
            measured(
                read_json(output / "measured-cal.json"),
                "calibration",
                8,
                digest(cal),
                digest(cal),
                digest(ee / "audit-binding-e83d195.json"),
                revision,
                auditor,
            )
            measured(
                read_json(output / "measured-final.json"),
                "final",
                24,
                digest(evidence / "native-final/final.json"),
                digest(cal),
                digest(evidence / "fixed-bath-binding.json"),
                revision,
                auditor,
            )
            original = read_json(evidence / "fixed-bath-audit.json")
            recomputed = read_json(output / "raw.json")
            rebinding = read_json(output / "raw-binding.json")
            require(
                original == recomputed
                and rebinding == read_json(evidence / "fixed-bath-binding.json"),
                "recomputed raw audit mismatch",
            )
            require(sources(source) == source_before, "source changed")
            require(
                inventory(evidence) == before
                and digest(archive) == archive_sha256
                and archive.stat().st_size == archive_size,
                "evidence changed",
            )
            receipt.update(
                status="passed",
                counts=counts,
                source_sha256=source_before,
                receipt_sha256={
                    name: digest(output / (name + ".json")) for name, _, _ in jobs
                },
            )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.SubprocessError,
        zipfile.BadZipFile,
    ) as exc:
        receipt["error"] = f"{type(exc).__name__}: {exc}"
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    if before is not None:
        try:
            require(
                inventory(evidence) == before and digest(archive) == archive_sha256,
                "final input inventory changed",
            )
        except (OSError, ValueError) as exc:
            receipt.update(status="failed", error=str(exc))
            (output / "receipt.json").write_text(
                json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
            )
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", required=True, type=Path)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--archive-size", required=True, type=int)
    p.add_argument("--run-id", required=True)
    p.add_argument("--source-revision", required=True)
    p.add_argument("--triqs-python", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    try:
        result = audit(
            a.archive,
            a.archive_sha256,
            a.archive_size,
            a.run_id,
            a.source_revision,
            a.triqs_python,
            a.output,
        )
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
