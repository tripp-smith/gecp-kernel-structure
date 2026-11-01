#!/usr/bin/env python3
"""EH release orchestration. Run in the clean frozen EH checkout and pinned runtime.

Evidence root layout follows the maxent-science workflow artifacts; parent and
cthyb-input are extracted beforehand using the project's safe extraction procedure.
The candidate is a descendant commit, not an alternate auditor checkout.
Receipt preflight alone is not scientific acceptance; the CLI always reruns audits.
"""

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = None
LOGS = None
FROZEN = "feaac669eb49d700f0059f5be755a9e4403a1719"


def fail(message):
    raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, stderr=subprocess.PIPE)


def load(path):
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 4 * 2**20:
        fail("invalid or oversized JSON receipt")

    def pairs(xs):
        d = {}
        for k, v in xs:
            if k in d:
                fail("duplicate JSON key: " + k)
            d[k] = v
        return d

    def constant(x):
        fail("nonfinite JSON constant: " + x)

    def finite_float(x):
        value = float(x)
        if not math.isfinite(value):
            fail("nonfinite JSON number: " + x)
        return value

    return json.loads(
        path.read_text(),
        object_pairs_hook=pairs,
        parse_constant=constant,
        parse_float=finite_float,
    )


def checked_file(p):
    if p.is_symlink() or not p.is_file():
        fail("missing or unsafe file: " + str(p))
    return p


def source_check(candidate):
    if git("rev-parse", "HEAD").decode().strip() != FROZEN:
        fail("auditors must execute at frozen EH HEAD")
    if git("status", "--porcelain", "--untracked-files=no").strip():
        fail("dirty auditor checkout")
    if not re.fullmatch("[0-9a-f]{40}", candidate):
        fail("full candidate commit required")
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", FROZEN, candidate],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    # Inventory is constructed from the frozen Git tree, never from evidence.
    names = git("ls-tree", "-r", "--name-only", FROZEN).decode().splitlines()
    selected = {
        p
        for p in names
        if (p.startswith("kernelgecp/") and p.endswith(".py"))
        or p.startswith("environments/locks/")
        or p.startswith("native/src/")
        or p in ("native/Cargo.toml", "native/Cargo.lock")
        or (p.startswith("experiments/") and p.endswith(".toml"))
    }
    selected |= {
        "pyproject.toml",
        "uv.lock",
        "experiments/spectroscopy.py",
        "experiments/spectroscopy_campaign.py",
        "experiments/audit_spectroscopy.py",
        "experiments/dlr_adapters.py",
        "experiments/dlr_physics.py",
        "experiments/dlr_workloads.py",
        "experiments/cthyb_dmft.py",
        "experiments/cthyb_loop_state.py",
        "experiments/cthyb_fixed_bath.py",
        "experiments/cthyb_spectroscopy.py",
        "experiments/audit_cthyb_dmft.py",
        "experiments/audit_cthyb_fixed_bath.py",
        "experiments/audit_cthyb_spectroscopy.py",
        "experiments/replay_cthyb_acceptance.py",
    }
    if not selected <= set(names):
        fail("frozen source inventory incomplete")
    # Adding Python modules under the scientific package is a code change too.
    candidate_names = set(
        git("ls-tree", "-r", "--name-only", candidate).decode().splitlines()
    )

    def guarded_inventory(paths):
        return {
            x
            for x in paths
            if (x.startswith("kernelgecp/") and x.endswith(".py"))
            or x.startswith("native/src/")
            or x.startswith("environments/locks/")
        }

    if guarded_inventory(candidate_names) != guarded_inventory(names):
        fail("scientific package or lock inventory changed")
    for path in sorted(selected):
        old = git("show", FROZEN + ":" + path)
        new = git("show", candidate + ":" + path)
        if old != new:
            fail("scientific source or lock changed: " + path)


def execute(cmd, name):
    with (LOGS / (name + ".log")).open("wb") as log:
        subprocess.run(
            cmd,
            cwd=ROOT,
            env={**os.environ, "PYTHONPATH": str(ROOT)},
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def run(module, args, output):
    cmd = [sys.executable, "-m", module, *map(str, args), "--output", str(output)]
    execute(cmd, output.stem)
    result = load(checked_file(output))
    shutil.copyfile(output, LOGS / output.name)
    return result


def sha(path):
    return hashlib.sha256(checked_file(path).read_bytes()).hexdigest()


def audit_summary(path, *, kind, cases, schema=None):
    value = load(checked_file(path))
    if (
        type(value) is not dict
        or value.get("status") != "passed"
        or value.get("kind") != kind
        or type(value.get("cases")) is not int
        or value["cases"] != cases
    ):
        fail("invalid independent audit summary: " + str(path))
    if schema is not None and value.get("schema") != schema:
        fail("wrong audit schema")
    return value


def synthetic_summary(path, manifest, freeze, kind, cases):
    value = audit_summary(path, kind=kind, cases=cases)
    if (
        value.get("independently_audited_cases") != 24
        or type(value.get("independently_audited_cases")) is not int
        or value.get("manifest_sha256") != sha(manifest)
        or value.get("freeze_sha256") != sha(freeze)
    ):
        fail("invalid synthetic audit inventory or binding")
    children = load(manifest).get("children")
    if type(children) is not dict or len(children) != cases:
        fail("synthetic children incomplete")
    return value


def receipt_preflight(e):
    """Structural/hash checks only. This function never establishes acceptance."""
    if e.is_symlink() or not e.is_dir():
        fail("invalid evidence root")
    for p in e.rglob("*"):
        if p.is_symlink():
            fail("symlink in evidence")
        if p.is_file() and p.suffix == ".json":
            load(p)

    def f(p):
        return checked_file(e / p)

    freeze, sc, sf = (
        f("freeze.json"),
        f("calibration/campaign.json"),
        f("campaign/campaign.json"),
    )
    synthetic_summary(f("calibration-audit.json"), sc, freeze, "calibration", 24)
    synthetic_summary(f("campaign-audit.json"), sf, freeze, "campaign", 72)
    cal, ca, final = (
        f("native-calibration/calibration.json"),
        f("native-calibration-audit.json"),
        f("native-final/final.json"),
    )
    schema = "kernelgecp.cthyb-spectroscopy-independent-audit.v1"
    ac = audit_summary(ca, kind="calibration", cases=8, schema=schema)
    af = audit_summary(
        f("native-final-audit.json"), kind="final", cases=24, schema=schema
    )
    raw, binding = f("fixed-bath-audit.json"), f("fixed-bath-binding.json")
    r, b = load(raw), load(binding)
    expected_gates = {
        "bootstrap_and_defaults",
        "covariance_and_projection",
        "ee_source_and_raw_chains",
        "fresh_process_receipts",
        "frozen_protocol_and_hashes",
        "heldout_rms",
        "native_trials_and_quadrature",
        "repeat_arrays",
    }
    for a in (ac, af):
        if (
            a.get("auditor_revision") != FROZEN
            or a.get("auditor_sha256")
            != hashlib.sha256(
                git("show", FROZEN + ":experiments/audit_cthyb_spectroscopy.py")
            ).hexdigest()
            or type(a.get("gates")) is not dict
            or set(a["gates"]) != expected_gates
            or any(x is not True for x in a["gates"].values())
        ):
            fail("measured auditor identity or gates mismatch")
    if (
        ac.get("campaign_sha256") != sha(cal)
        or ac.get("calibration_sha256") != sha(cal)
        or ac.get("binding_sha256")
        != sha(ROOT / "data/cthyb/audit-binding-e83d195.json")
    ):
        fail("calibration binding mismatch")
    if (
        af.get("campaign_sha256") != sha(final)
        or af.get("calibration_sha256") != sha(cal)
        or af.get("binding_sha256") != sha(binding)
    ):
        fail("final binding mismatch")
    if (
        r.get("schema") != "kernelgecp.cthyb-fixed-bath-audit.v1"
        or b.get("schema") != "kernelgecp.cthyb-fixed-bath-audit-binding.v1"
        or r.get("status") != "passed"
        or b.get("status") != "passed"
    ):
        fail("raw audit schema/status")
    if (
        r.get("audit_revision") != FROZEN
        or r.get("audit_source_sha256")
        != hashlib.sha256(
            git("show", FROZEN + ":experiments/audit_cthyb_fixed_bath.py")
        ).hexdigest()
    ):
        fail("raw auditor identity")
    if (r.get("chain_count"), r.get("train_count"), r.get("heldout_count")) != (
        256,
        128,
        128,
    ):
        fail("raw count mismatch")
    if b.get("receipt_sha256") != sha(raw):
        fail("raw receipt binding mismatch")
    for key in (
        "audit_revision",
        "audit_source_sha256",
        "calibration",
        "manifest_sha256",
        "parent_binding_sha256",
        "result_sha256",
    ):
        if b.get(key) != r.get(key):
            fail("raw binding mismatch: " + key)
    if r.get("manifest_sha256") != sha(f("fixed-bath/manifest.json")) or r.get(
        "result_sha256"
    ) != sha(f("fixed-bath/result.json")):
        fail("raw artifact mismatch")
    if r.get("calibration", {}).get("receipt_sha256") != sha(cal) or r.get(
        "calibration", {}
    ).get("audit_sha256") != sha(ca):
        fail("raw calibration mismatch")
    expected = {
        f"chains/chain-{i:03d}{ext}"
        for i in range(256)
        for ext in (".json", ".npz", ".timing.json")
    }
    inv = b.get("inventory")
    if type(inv) is not dict or set(inv) != expected:
        fail("raw inventory mismatch")
    for name, digest in inv.items():
        if digest != sha(e / "fixed-bath" / name):
            fail("raw inventory hash mismatch: " + name)
    return {
        "status": "structurally_consistent",
        "synthetic_calibration": 24,
        "synthetic_final": 72,
        "measured_calibration": 8,
        "measured_final": 24,
        "raw_chains": 256,
    }


def check(e, candidate):
    source_check(candidate)
    receipt_preflight(e)

    def f(p):
        return checked_file(e / p)

    freeze = f("freeze.json")
    sc = f("calibration/campaign.json")
    sf = f("campaign/campaign.json")
    cal = f("native-calibration/calibration.json")
    ca = f("native-calibration-audit.json")
    final = f("native-final/final.json")
    eeaudit = checked_file(ROOT / "data/cthyb/independent-audit-e83d195.json")
    eebind = checked_file(ROOT / "data/cthyb/audit-binding-e83d195.json")
    parentbinding = checked_file(ROOT / "data/cthyb/linux-audit-binding-7200ec0.json")
    with tempfile.TemporaryDirectory(prefix="eh-audit-") as t:
        t = Path(t)
        for manifest, audit, kind, count in [
            (sc, "calibration-audit.json", "calibration", 24),
            (sf, "campaign-audit.json", "campaign", 72),
        ]:
            run(
                "experiments.audit_spectroscopy",
                ["--manifest", manifest, "--freeze", freeze],
                t / (kind + ".json"),
            )
            synthetic_summary(t / (kind + ".json"), manifest, freeze, kind, count)
            synthetic_summary(e / audit, manifest, freeze, kind, count)
        run(
            "experiments.audit_cthyb_spectroscopy",
            [
                cal,
                "--ee-root",
                e / "cthyb-input/gecp-ee-odd-9d28053",
                "--ee-audit",
                eeaudit,
                "--ee-binding",
                eebind,
            ],
            t / "measured-cal.json",
        )
        schema = "kernelgecp.cthyb-spectroscopy-independent-audit.v1"
        audit_summary(
            t / "measured-cal.json", kind="calibration", cases=8, schema=schema
        )
        oldcal = audit_summary(ca, kind="calibration", cases=8, schema=schema)
        if (
            oldcal.get("campaign_sha256") != sha(cal)
            or oldcal.get("calibration_sha256") != sha(cal)
            or oldcal.get("binding_sha256") != sha(eebind)
        ):
            fail("measured calibration binding mismatch")
        raw_a, raw_b = t / "raw-audit.json", t / "raw-binding.json"
        execute(
            [
                sys.executable,
                "-m",
                "experiments.audit_cthyb_fixed_bath",
                "--folder",
                str(e / "fixed-bath"),
                "--parent-dir",
                str(e / "parent"),
                "--parent-binding",
                str(parentbinding),
                "--receipt",
                str(raw_a),
                "--binding",
                str(raw_b),
                "--calibration",
                str(cal),
                "--calibration-audit",
                str(ca),
            ],
            "fixed-bath",
        )
        for raw_output in (raw_a, raw_b):
            load(checked_file(raw_output))
            shutil.copyfile(raw_output, LOGS / raw_output.name)
        if load(checked_file(raw_a)).get("chain_count") != 256:
            fail("raw cohort incomplete")
        # The final auditor consumes original raw receipt/binding bytes: final
        # campaign embeds their digests. Recomputed output is a separate check.
        original_raw, original_binding = (
            f("fixed-bath-audit.json"),
            f("fixed-bath-binding.json"),
        )
        oldraw, oldbinding = load(original_raw), load(original_binding)
        freshraw = load(raw_a)
        for key in (
            "schema",
            "status",
            "claim_scope",
            "chain_count",
            "train_count",
            "heldout_count",
            "manifest_sha256",
            "result_sha256",
            "parent_binding_sha256",
            "calibration",
            "audit_revision",
            "audit_source_sha256",
        ):
            if oldraw.get(key) != freshraw.get(key):
                fail("original raw audit mismatch: " + key)
        if oldbinding.get("receipt_sha256") != sha(original_raw):
            fail("original raw receipt binding mismatch")
        run(
            "experiments.audit_cthyb_spectroscopy",
            [
                final,
                "--ee-root",
                e / "fixed-bath",
                "--ee-audit",
                original_raw,
                "--ee-binding",
                original_binding,
                "--calibration",
                cal,
                "--calibration-audit",
                ca,
            ],
            t / "final.json",
        )
        audit_summary(t / "final.json", kind="final", cases=24, schema=schema)
        oldfinal = audit_summary(
            e / "native-final-audit.json", kind="final", cases=24, schema=schema
        )
        if (
            oldfinal.get("campaign_sha256") != sha(final)
            or oldfinal.get("calibration_sha256") != sha(cal)
            or oldfinal.get("binding_sha256") != sha(original_binding)
        ):
            fail("measured final binding mismatch")
    return {
        "schema": "kernelgecp.spectroscopy-evidence-replay.v1",
        "status": "passed",
        "scope": "EH finite scientific campaign; frozen-data replay, no new sampling",
        "audit_receipts": {
            path.name: sha(path) for path in sorted(LOGS.glob("*.json"))
        },
        "candidate": candidate,
        "auditor_revision": FROZEN,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("--candidate", required=True)
    p.add_argument("--frozen-checkout", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    global ROOT, LOGS
    try:
        ROOT = a.frozen_checkout.resolve(strict=True)
        if not sys.platform.startswith("linux"):
            fail("scientific audit replay requires Linux")
        a.output.parent.mkdir(parents=True, exist_ok=True)
        LOGS = a.output.parent / (a.output.name + ".logs")
        LOGS.mkdir(exist_ok=False)
        result = check(a.evidence, a.candidate)
        code = 0
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.SubprocessError,
    ) as exc:
        result = {"status": "failed", "error": str(exc)}
        code = 1
    try:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(
            json.dumps(result, allow_nan=False, sort_keys=True, indent=2) + "\n"
        )
    except OSError as exc:
        print("unable to write result: " + str(exc), file=sys.stderr)
        return 1
    return code


if __name__ == "__main__":
    raise SystemExit(main())
