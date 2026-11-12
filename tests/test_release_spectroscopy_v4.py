import zipfile
from pathlib import Path

import pytest

from scripts import release_spectroscopy_v4 as replay


def test_calibration_only_archive_fails_closed(tmp_path):
    archive = Path("/tmp/gecp-eh-0e21-calibration-passed.zip")
    if not archive.is_file():
        pytest.skip("optional authenticated development artifact unavailable")
    replay.extract(archive, tmp_path)
    with pytest.raises((ValueError, OSError)):
        replay.evidence_checks(tmp_path, tmp_path, "0" * 40)


def test_unsafe_zip_rejected(tmp_path):
    z = tmp_path / "bad.zip"
    with zipfile.ZipFile(z, "w") as archive:
        archive.writestr("../escape", "bad")
    with pytest.raises(ValueError):
        replay.extract(z, tmp_path / "out")
    assert not (tmp_path.parent / "escape").exists()


def test_synthetic_audit_requires_exact_bindings():
    good = dict(
        status="passed",
        kind="campaign",
        cases=72,
        independently_audited_cases=24,
        manifest_sha256="a",
        freeze_sha256="b",
    )
    replay.synthetic(good, "campaign", 72, "a", "b")
    for key, value in [("cases", 24), ("freeze_sha256", "bad"), ("status", "failed")]:
        with pytest.raises(ValueError):
            replay.synthetic(dict(good, **{key: value}), "campaign", 72, "a", "b")


def test_measured_audit_rejects_failed_or_missing_gate():
    good = dict(
        schema="kernelgecp.cthyb-spectroscopy-independent-audit.v4",
        status="passed",
        kind="final",
        cases=24,
        campaign_sha256="a",
        calibration_sha256="b",
        binding_sha256="c",
        auditor_revision="d",
        auditor_sha256="e",
        gates=dict.fromkeys(replay.GATES, True),
    )
    replay.measured(good, "final", 24, "a", "b", "c", "d", "e")
    for gates in ({}, dict.fromkeys(replay.GATES, False)):
        with pytest.raises(ValueError):
            replay.measured(
                dict(good, gates=gates), "final", 24, "a", "b", "c", "d", "e"
            )


def test_authentication_rejects_wrong_digest(tmp_path, monkeypatch):
    monkeypatch.setattr(replay, "checked_path", lambda p: Path(p).absolute())
    archive = tmp_path / "maxent-fixed-bath-final-1-1.zip"
    archive.write_bytes(b"not an archive")
    with pytest.raises(ValueError, match="authentication"):
        replay.audit(
            archive,
            "0" * 64,
            archive.stat().st_size,
            "1",
            "0" * 40,
            Path("/usr/bin/python3"),
            tmp_path / "output",
        )
    assert not (tmp_path / "output").exists()


def shard_fixture(root):
    import json

    manifest = {
        "controls": {"profile": "v4"},
        "provenance": {"source_revision": "a"},
        "children": {},
        "shard_receipts": {},
    }
    (root / "shard-receipts").mkdir()
    jobs = []
    for ci, case in enumerate(
        ("offgrid_atoms_independent", "offgrid_atoms_correlated", "smooth_correlated")
    ):
        for method in ("raw", "kernelgecp", "pydlr", "sparse_ir"):
            for model in ("flat", "broad_gaussian"):
                for repeat in range(3):
                    name = f"{case}-{method}-{model}-r{repeat}"
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
                    jobs.append((name, job))
    for ordinal in range(6):
        part = {}
        inv = {}
        for name, job in jobs[ordinal * 12 : (ordinal + 1) * 12]:
            part[name] = {"job": job, "status": "passed", "exit_code": 0}
            suffixes = [".json", ".npz", ".log", "-dataset.json", "-dataset.npz"]
            if job["method"] == "kernelgecp":
                suffixes.append("-gecp-fit.npz")
            for suffix in suffixes:
                file = root / (name + suffix)
                file.write_text(name + suffix)
                inv[file.name] = replay.digest(file)
        manifest["children"].update(part)
        receipt = {
            "schema": "kernelgecp.spectroscopy-campaign-shard.v4",
            "kind": "campaign",
            "ordinal": ordinal,
            "status": "passed",
            "controls": manifest["controls"],
            "provenance": manifest["provenance"],
            "freeze_sha256": "frozen",
            "repeat_comparison": {"status": "passed"},
            "children": part,
            "inventory": inv,
        }
        path = root / "shard-receipts" / f"shard-{ordinal}.json"
        path.write_text(json.dumps(receipt))
        manifest["shard_receipts"][str(ordinal)] = {
            "path": f"shard-receipts/shard-{ordinal}.json",
            "sha256": replay.digest(path),
        }
    return manifest


def test_shard_structure_accepts_complete_partition(tmp_path):
    manifest = shard_fixture(tmp_path)
    replay.check_shards(tmp_path, manifest, "frozen")


@pytest.mark.parametrize(
    "attack", ["missing", "digest", "ordinal", "cohort", "inventory"]
)
def test_shard_structure_rejects_tampering(tmp_path, attack):
    import json

    manifest = shard_fixture(tmp_path)
    ref = manifest["shard_receipts"]["0"]
    path = tmp_path / ref["path"]
    if attack == "missing":
        del manifest["shard_receipts"]["5"]
    elif attack == "digest":
        ref["sha256"] = "0" * 64
    elif attack == "inventory":
        receipt = json.loads(path.read_text())
        victim = next(iter(receipt["inventory"]))
        (tmp_path / victim).write_text("tampered")
    else:
        receipt = json.loads(path.read_text())
        if attack == "ordinal":
            receipt["ordinal"] = 1
        else:
            wrong = next(iter(manifest["children"].keys() - receipt["children"].keys()))
            receipt["children"][wrong] = manifest["children"][wrong]
        path.write_text(json.dumps(receipt))
        ref["sha256"] = replay.digest(path)
    with pytest.raises(ValueError):
        replay.check_shards(tmp_path, manifest, "frozen")
