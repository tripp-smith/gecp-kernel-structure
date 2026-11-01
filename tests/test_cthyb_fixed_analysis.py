import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from experiments import audit_cthyb_spectroscopy as audit
from experiments import cthyb_spectroscopy as m


def test_fixed_split_statistics_and_overlap():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(256, 3073))
    data, _, held, cov = m.dataset_and_validation(x, {})
    np.testing.assert_allclose(data.values, x[:128, m.OBS].mean(axis=0))
    np.testing.assert_allclose(
        data.covariance, np.cov(x[:128, m.OBS], rowvar=False) / 128
    )
    np.testing.assert_allclose(held, x[128:, m.QUERY].mean(axis=0))
    np.testing.assert_allclose(cov, np.cov(x[128:, m.QUERY], rowvar=False) / 128)
    with pytest.raises(ValueError, match="overlap"):
        m.validate_cohort_pair(
            {
                "generation_revision": "a",
                "binding_sha256": "a",
                "seeds": [m.FIXED_SEEDS[0]],
            },
            {"generation_revision": "b", "binding_sha256": "b", "seeds": m.FIXED_SEEDS},
        )


def fixture(tmp_path, monkeypatch):
    chains = tmp_path / "chains"
    chains.mkdir()
    bath = np.array([1 + 0j])
    cfg = {"beta": 10}
    calibration = {
        "receipt_sha256": "a" * 64,
        "audit_sha256": "b" * 64,
        "binding_sha256": "c" * 64,
        "auditor_revision": "1" * 40,
        "auditor_sha256": "d" * 64,
    }
    manifest = {
        "schema": "kernelgecp.cthyb-fixed-bath.v1",
        "protocol": {
            "seeds": m.FIXED_SEEDS,
            "train_count": 128,
            "cycles": 1600000,
            "config": cfg,
        },
        "parent": {
            "binding_sha256": m.PARENT_BINDING,
            "bath_bytes_sha256": hashlib.sha256(bath.tobytes()).hexdigest(),
        },
        "calibration": calibration,
        "provenance": {
            "revision": "1" * 40,
            "sources": {"experiments/audit_cthyb_spectroscopy.py": "d" * 64},
        },
    }

    def write(path, obj):
        path.write_text(json.dumps(obj))

    write(tmp_path / "manifest.json", manifest)
    mh = m.sha(tmp_path / "manifest.json")
    entries, inventory = [], {}
    for i, seed in enumerate(m.FIXED_SEEDS):
        role = "train" if i < 128 else "heldout"
        stem = chains / f"chain-{i:03d}"
        np.savez(
            stem.with_suffix(".npz"), g_tau=np.full((3073, 2), i + 0j), bath_iw=bath
        )
        write(
            stem.with_suffix(".json"),
            {
                "schema": "kernelgecp.cthyb-measurement.v1",
                "identity": {"manifest_sha256": mh, "index": i, "role": role},
                "seed": seed,
                "config": cfg,
                "cycles": 1600000,
                "solve_status": 0,
                "average_sign": 1,
                "arrays_sha256": m.sha(stem.with_suffix(".npz")),
            },
        )
        write(stem.with_suffix(".timing.json"), {})
        paths = [stem.with_suffix(x) for x in (".json", ".npz", ".timing.json")]
        hs = [m.sha(p) for p in paths]
        inventory.update(
            {str(p.relative_to(tmp_path)): h for p, h in zip(paths, hs, strict=True)}
        )
        entries.append(
            dict(
                index=i,
                seed=seed,
                role=role,
                receipt_sha256=hs[0],
                arrays_sha256=hs[1],
                timing_sha256=hs[2],
            )
        )
    write(
        tmp_path / "result.json",
        {
            "schema": "kernelgecp.cthyb-fixed-bath-result.v1",
            "status": "complete",
            "manifest_sha256": mh,
            "chains": entries,
        },
    )
    source_hash = hashlib.sha256(b"audit").hexdigest()
    receipt = tmp_path / "audit.json"
    write(
        receipt,
        {
            "schema": "kernelgecp.cthyb-fixed-bath-audit.v1",
            "status": "passed",
            "parent_binding_sha256": m.PARENT_BINDING,
            "calibration": calibration,
            "manifest_sha256": mh,
            "result_sha256": m.sha(tmp_path / "result.json"),
            "chain_count": 256,
            "train_count": 128,
            "heldout_count": 128,
            "audit_revision": "2" * 40,
            "audit_source_sha256": source_hash,
        },
    )
    binding = tmp_path / "binding.json"
    write(
        binding,
        {
            "schema": "kernelgecp.cthyb-fixed-bath-audit-binding.v1",
            "status": "passed",
            "calibration": calibration,
            "receipt_sha256": m.sha(receipt),
            "manifest_sha256": mh,
            "result_sha256": m.sha(tmp_path / "result.json"),
            "parent_binding_sha256": m.PARENT_BINDING,
            "inventory": inventory,
            "audit_revision": "2" * 40,
            "audit_source_sha256": source_hash,
        },
    )
    monkeypatch.setattr(m.subprocess, "check_output", lambda *a, **k: b"audit")
    monkeypatch.setattr(
        m,
        "load_measurement",
        lambda stem, identity: SimpleNamespace(
            seed=m.FIXED_SEEDS[identity["index"]],
            cycles=1600000,
            solve_status=0,
            average_sign=1,
            config=SimpleNamespace(),
            bath_iw=bath,
            g_tau=np.full((3073, 2), identity["index"] + 0j),
        ),
    )
    monkeypatch.setattr(m.dataclasses, "asdict", lambda cfg: {"beta": 10})
    return receipt, binding


def test_fixed_bath_reader_and_independent_audit(tmp_path, monkeypatch):
    receipt, binding = fixture(tmp_path, monkeypatch)
    samples, cohort = m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    stats, independent = audit.fixed_raw_cohort(tmp_path, receipt, binding)
    assert samples.shape == (256, 3073)
    assert cohort == independent
    np.testing.assert_array_equal(samples, stats[0])
    (tmp_path / "chains/chain-015.npz").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="inventory hash"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    with pytest.raises(ValueError, match="inventory hash"):
        audit.fixed_raw_cohort(tmp_path, receipt, binding)


def test_wrong_split_and_missing_rejected(tmp_path, monkeypatch):
    receipt, binding = fixture(tmp_path, monkeypatch)
    (tmp_path / "chains/chain-255.timing.json").unlink()
    with pytest.raises((ValueError, FileNotFoundError)):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    with pytest.raises((ValueError, FileNotFoundError)):
        audit.fixed_raw_cohort(tmp_path, receipt, binding)


@pytest.mark.parametrize("change", ["split", "seeds", "cohort", "parent"])
def test_rebound_invalid_fixed_manifest_rejected(tmp_path, monkeypatch, change):
    receipt, binding = fixture(tmp_path, monkeypatch)
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if change == "split":
        manifest["protocol"]["train_count"] = 127
    elif change == "seeds":
        manifest["protocol"]["seeds"][0] = 29140001
    elif change == "cohort":
        manifest["schema"] = "kernelgecp.cthyb-fixed-bath-other.v1"
    else:
        manifest["parent"]["binding_sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest))
    manifest_hash = m.sha(manifest_path)
    result_path = tmp_path / "result.json"
    result = json.loads(result_path.read_text())
    result["manifest_sha256"] = manifest_hash
    result_path.write_text(json.dumps(result))
    a = json.loads(receipt.read_text())
    a["manifest_sha256"] = manifest_hash
    a["result_sha256"] = m.sha(result_path)
    receipt.write_text(json.dumps(a))
    b = json.loads(binding.read_text())
    b["manifest_sha256"] = manifest_hash
    b["result_sha256"] = m.sha(result_path)
    b["receipt_sha256"] = m.sha(receipt)
    binding.write_text(json.dumps(b))
    with pytest.raises(ValueError, match="manifest, split or result"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    with pytest.raises(ValueError, match="manifest, split or result"):
        audit.fixed_raw_cohort(tmp_path, receipt, binding)


@pytest.mark.parametrize("part", ["missing", "receipt", "binding", "source"])
def test_frozen_calibration_raw_binding_rejected(tmp_path, monkeypatch, part):
    receipt, binding = fixture(tmp_path, monkeypatch)
    b = json.loads(binding.read_text())
    if part == "missing":
        b.pop("calibration")
    elif part == "binding":
        b["calibration"]["receipt_sha256"] = "e" * 64
    elif part == "receipt":
        a = json.loads(receipt.read_text())
        a["calibration"]["audit_sha256"] = "e" * 64
        receipt.write_text(json.dumps(a))
        b["receipt_sha256"] = m.sha(receipt)
    else:
        b["calibration"]["auditor_sha256"] = "e" * 64
    binding.write_text(json.dumps(b))
    with pytest.raises(ValueError, match="frozen calibration binding"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    with pytest.raises(ValueError, match="frozen calibration binding"):
        audit.fixed_raw_cohort(tmp_path, receipt, binding)


def test_final_arguments_must_match_frozen_calibration(tmp_path):
    cal_path, audit_path = tmp_path / "cal.json", tmp_path / "audit.json"
    cal_path.write_text("calibration")
    audit_path.write_text("independent audit")
    cal = {"binding_sha256": "c" * 64}
    independent = {"auditor_revision": "1" * 40, "auditor_sha256": "d" * 64}
    ref = {
        "receipt_sha256": m.sha(cal_path),
        "audit_sha256": m.sha(audit_path),
        "binding_sha256": cal["binding_sha256"],
        **independent,
    }
    m.check_fixed_calibration(
        {"calibration": ref}, cal_path, audit_path, cal, independent
    )
    audit_path.write_text("changed audit")
    with pytest.raises(ValueError, match="pre-acquisition freeze"):
        m.check_fixed_calibration(
            {"calibration": ref}, cal_path, audit_path, cal, independent
        )
