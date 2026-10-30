import json
from pathlib import Path

import numpy as np
import pytest

from experiments import audit_epw_restart_roundtrip as audit


def test_archived_training_coordinates_and_nearby_condition_boundary():
    fixture = Path(__file__).parent / "data/epw/restart20-a9c2c39/bundle.npz"
    with np.load(fixture, allow_pickle=False) as archive:
        rows = archive["fit_rows"]
        archived = archive["observation_x"]
        nodes = archive["nodes"]
        poles = archive["poles"]
        weights = archive["weights"]

    checked = audit.training_coordinates(archived, rows)
    design, condition = audit.direct_design(checked, nodes, poles, weights)
    assert condition < 1e4
    assert rows[0] == 0 and archived[0] in nodes
    assert design[0, np.flatnonzero(nodes == archived[0])[0]] == 1
    nearby = archived[0]
    for _ in range(3):
        nearby = np.nextafter(nearby, 0.0)
    assert nearby == ((np.pi * audit.KB * 20) ** 2)
    with pytest.raises(ValueError, match="binary64 mismatch"):
        audit.training_coordinates(np.array([nearby]), [0])
    with pytest.raises(ValueError, match="cardinal evaluation condition"):
        audit.direct_design(np.array([nearby]), nodes, poles, weights)


def test_direct_product_interpolates_nodes_and_rejects_duplicate_nodes():
    nodes = np.array([1.0, 2.0])
    poles = np.array([0.0, -1.0])
    weights = np.array([-2.0, 6.0])
    design, _ = audit.direct_design(nodes, nodes, poles, weights)
    np.testing.assert_array_equal(design, np.eye(2))
    with pytest.raises(ValueError):
        audit.direct_design(np.array([1.0]), np.array([1.0, 1.0]), poles, weights)


def test_native_nambu_algebra_recovers_synthetic_fields():
    data = {
        "frequency": np.array([0.01, 0.03]),
        "energy": np.array([-0.2, 0.1]),
        "fields": {
            0: [[1.2, 0.003, 0.01], [1.3, 0.004, -0.02]],
            1: [[1.1, 0.002, 0.03], [1.4, 0.005, 0.01]],
        },
        "context": {"chemical_potential_offset_ev": 0.002},
    }
    g, z, d, c = audit.truth(data, [0, 1])
    rz, rd, rc, _, _ = audit.recover(g, data["frequency"], data["energy"], 0.002)
    np.testing.assert_allclose(rz, z)
    np.testing.assert_allclose(rd, d)
    np.testing.assert_allclose(rc, c, atol=1e-14)


def test_hash_and_xor_little_endian_bitwise():
    target = np.array([1.0, -0.0, 1e-200])
    baseline = np.array([2.0, 0.0, 1e-100])
    correction = audit.bits(target) ^ audit.bits(baseline)
    assert correction.dtype == np.dtype("<u8")
    np.testing.assert_array_equal(audit.bits(baseline) ^ correction, audit.bits(target))
    assert audit.hash_bits(audit.bits(target)) != audit.hash_bits(audit.bits(baseline))


def test_duplicate_json_and_missing_campaign_fail_closed(tmp_path):
    p = tmp_path / "x.json"
    p.write_text('{"a":1,"a":2}')
    with pytest.raises(ValueError, match="duplicate"):
        audit.read(p)
    p.write_text(json.dumps({"schema": "wrong"}))
    with pytest.raises((KeyError, ValueError, FileNotFoundError)):
        audit.audit_campaign(p)


def test_heldout_budget_is_fit_only():
    m = {
        "green_scaled": 1e-7,
        "z": 1e-6,
        "gap_ev": 1e-7,
        "chi_ev": 1e-7,
        "imaginary_defect": 0.0,
    }
    assert audit.passed(m, 1)
    m["gap_ev"] = 1.1e-6
    assert not audit.passed(m, 1)


def test_synthetic_bundle_rejects_inventory_and_hash_tampering(tmp_path):
    import hashlib

    payload = tmp_path / "bundle.npz"
    np.savez(
        payload,
        __metadata__=np.frombuffer(b"{}", dtype=np.uint8),
        corrections=np.zeros((2, 1, 3), dtype="<u8"),
    )
    record = {
        "schema": "kernelgecp.epw-restart-xor-application.v1",
        "policy": "lossless_native_xor_v1",
        "units": "eV",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "state_order": "frequency_major_ik_ibnd_ordinal",
        "extra_encoded_information": (
            "uint64 XOR corrections for all native Z, gap, and chi fields"
        ),
        "archive": {
            "file": "bundle.npz",
            "sha256": hashlib.sha256(payload.read_bytes()).hexdigest(),
        },
    }
    record["record_sha256"] = audit.digest_json(record)
    path = tmp_path / "bundle.json"
    path.write_text(json.dumps(record))
    data = {
        "context": {
            "provider": record["provider"],
            "state_order": record["state_order"],
        }
    }
    with pytest.raises(ValueError, match="inventory"):
        audit.bundle(path, data)
    payload.write_bytes(payload.read_bytes() + b"tamper")
    with pytest.raises(ValueError, match="hash"):
        audit.bundle(path, data)
