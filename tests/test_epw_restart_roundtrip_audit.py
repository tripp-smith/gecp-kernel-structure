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


def portable_fixture(tmp_path):
    """Minimal structural fixture, not a completed scientific campaign."""
    freeze = {
        "schema": "kernelgecp.epw-restart-roundtrip-freeze.v1",
        "controls": audit.CONTROLS,
        "offset": 0.0,
        "context": {"chemical_potential_offset_ev": 0.0},
        "environment": {
            "executable": "/old/venv/bin/python",
            "revision": "a" * 40,
            "source_hashes": {"experiments/epw_restart_roundtrip.py": "legacy"},
        },
        "native_root": "/old/artifacts",
    }
    old_root = "/old/artifacts/restart20"
    children = []
    for i in range(3):
        folder = tmp_path / f"case-{i}"
        folder.mkdir()
        (folder / "case.json").write_text("{}")
        log = tmp_path / f"case-{i}.log"
        log.write_text("worker")
        children.append(
            {
                "repeat": i,
                "exit_code": 0,
                "log_sha256": audit.sha(log),
                "receipt_sha256": audit.sha(folder / "case.json"),
                "command": [
                    "/old/venv/bin/python",
                    "-m",
                    "experiments.epw_restart_roundtrip",
                    "case",
                    "--native-root",
                    "/old/artifacts",
                    "--temperature-kelvin",
                    "20",
                    "--chemical-potential-offset-ev",
                    "0.0",
                    "--output",
                    f"{old_root}/case-{i}",
                    "--freeze",
                    f"{old_root}/freeze.json",
                    "--repeat",
                    str(i),
                ],
            }
        )
    freeze_path = tmp_path / "freeze.json"
    freeze_path.write_text(json.dumps(freeze))
    campaign = {
        "schema": "kernelgecp.epw-restart-roundtrip-campaign.v1",
        "status": "passed",
        "native_acceptance": "not_tested",
        "repeat_arrays_identical": True,
        "freeze_sha256": audit.sha(freeze_path),
        "children": children,
    }
    manifest = tmp_path / "campaign.json"
    manifest.write_text(json.dumps(campaign))
    return manifest, freeze, campaign


def test_relocation_uses_override_without_rewriting_records(tmp_path, monkeypatch):
    manifest, freeze, _ = portable_fixture(tmp_path)
    before = {
        str(p.relative_to(tmp_path)): audit.sha(p)
        for p in tmp_path.rglob("*")
        if p.is_file()
    }
    native_root = tmp_path / "relocated-native"
    native_root.mkdir()
    calls = []
    monkeypatch.setattr(audit, "provenance", lambda env: calls.append(env["revision"]))

    def fake_native(root, context):
        assert root == native_root and context == freeze["context"]
        calls.append("relocated")
        return {"frequency": [0.1]}

    monkeypatch.setattr(audit, "native", fake_native)
    # The stubbed fixture intentionally stops at the first incomplete case receipt.
    with pytest.raises(KeyError, match="schema"):
        audit.audit_campaign(manifest, native_root_override=native_root)
    assert calls == ["a" * 40, "relocated"]
    after = {
        str(p.relative_to(tmp_path)): audit.sha(p)
        for p in tmp_path.rglob("*")
        if p.is_file()
    }
    assert after == before


@pytest.mark.parametrize(
    "change", ["argv", "traversal", "mixed", "freeze", "native_root", "executable"]
)
def test_relocation_rejects_bad_original_bindings(tmp_path, change):
    _, freeze, campaign = portable_fixture(tmp_path)
    children = campaign["children"]
    if change == "argv":
        children[1]["command"][7] = "17"
    elif change == "traversal":
        children[1]["command"][11] = "/old/artifacts/restart20/../restart20/case-1"
    elif change == "mixed":
        children[2]["command"][13] = "/different/restart20/freeze.json"
    elif change == "freeze":
        children[0]["command"][13] = "/old/artifacts/restart20/not-freeze.json"
    elif change == "native_root":
        freeze["native_root"] = "/old/artifacts/../artifacts"
    else:
        freeze["environment"]["executable"] = "/old/./venv/bin/python"
    with pytest.raises(ValueError):
        audit.validate_relocated_commands(children, freeze)


def test_relocation_rejects_altered_frozen_hash_before_source_access(
    tmp_path, monkeypatch
):
    manifest, _, _ = portable_fixture(tmp_path)
    (tmp_path / "freeze.json").write_text("{}")
    monkeypatch.setattr(
        audit, "native", lambda *_: pytest.fail("read native prematurely")
    )
    with pytest.raises(ValueError, match="campaign"):
        audit.audit_campaign(manifest, native_root_override=tmp_path)


def test_relocation_preserves_provenance_and_native_context_checks(
    tmp_path, monkeypatch
):
    manifest, freeze, _ = portable_fixture(tmp_path)

    def reject_provenance(env):
        assert env["revision"] == freeze["environment"]["revision"]
        raise ValueError("producer source hash mismatch")

    monkeypatch.setattr(audit, "provenance", reject_provenance)
    with pytest.raises(ValueError, match="producer source hash mismatch"):
        audit.audit_campaign(manifest, native_root_override=tmp_path)
    monkeypatch.setattr(audit, "provenance", lambda _: None)

    def reject_context(root, context):
        assert root == tmp_path and context == freeze["context"]
        raise ValueError("native context or file hash mismatch")

    monkeypatch.setattr(audit, "native", reject_context)
    with pytest.raises(ValueError, match="native context or file hash mismatch"):
        audit.audit_campaign(manifest, native_root_override=tmp_path)


def test_native_override_rechecks_producer_blob_context_and_file_hashes(
    tmp_path, monkeypatch
):
    import hashlib

    producer = b"synthetic producer source"
    monkeypatch.setattr(audit, "blob", lambda revision, path: producer)
    work = tmp_path / "epw-native/test-suite/epw_super"
    work.mkdir(parents=True)
    (work / "epw6.in").write_bytes(b"native input")
    (work / "epw6.out").write_text(
        "Fermi level (eV) = 1.0E+00\n"
        "Chemical potential (itemp = 2) = 1.0E+00\n"
        + "Convergence was reached in nsiter\n"
        * 3
    )
    files = []
    for temp in (17, 20, 30):
        name = f"MgB2.imag_aniso_{temp:06.2f}"
        (work / name).write_bytes(b"synthetic native data")
        files.append(
            {
                "temperature_kelvin": temp,
                "filename": name,
                "sha256": audit.sha(work / name),
                "bytes": (work / name).stat().st_size,
                "frequency_count": 54,
                "state_count": 1,
            }
        )
    binary = "b" * 64
    receipt = {
        "schema": "kernelgecp.epw-native-anisotropic-fixture.v1",
        "status": "passed",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "wscut_ev": 1.0,
        "source_revision": "a" * 40,
        "script_sha256": hashlib.sha256(producer).hexdigest(),
        "baseline": {
            "checks": {"native": True},
            "native_executions": [{"input": "epw1.in", "executable_sha256": binary}],
        },
        "native_anisotropic_execution": {
            "returncode": 0,
            "executable_sha256": binary,
            "input_sha256": audit.sha(work / "epw6.in"),
            "log_sha256": audit.sha(work / "epw6.out"),
        },
        "upstream_input_sha256": (
            "a67e51f747674319d3c5659e402529bcdf3d571b54f12d2802f2531d74ba8a82"
        ),
        "upstream_io_sha256": (
            "416ee21c65fd49b2a04894c2d5e18ffb1a5e8f7c223cab950a28f760cf56ec36"
        ),
        "executed_input_sha256": (
            "7a74f36be9a1604ec8b4184e1bb985f14d237835ec43b19c6ac7e9f4dc2937e1"
        ),
        "files": files,
    }
    # The real native input has a pinned hash; simulate only that constant via sha
    # while leaving the other file checks genuine.
    real_sha = audit.sha

    def fixture_sha(path):
        if path == work / "epw6.in":
            return receipt["executed_input_sha256"]
        return real_sha(path)

    receipt["native_anisotropic_execution"]["input_sha256"] = receipt[
        "executed_input_sha256"
    ]
    monkeypatch.setattr(audit, "sha", fixture_sha)
    (tmp_path / "epw-anisotropic.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="context"):
        audit.native(tmp_path, {"chemical_potential_offset_ev": 0.0})
    (work / files[0]["filename"]).write_bytes(b"altered")
    with pytest.raises(ValueError, match="corpus hash"):
        audit.native(tmp_path, {"chemical_potential_offset_ev": 0.0})
    monkeypatch.setattr(audit, "blob", lambda revision, path: b"changed")
    with pytest.raises(ValueError, match="native producer blob"):
        audit.native(tmp_path, {"chemical_potential_offset_ev": 0.0})


def interpreter_record():
    return {
        "schema": "kernelgecp.epw-interpreter-identity.v1",
        "invocation": "/old/venv/bin/python-link",
        "resolved_executable": "/old/venv/bin/python",
        "binary_sha256": "b" * 64,
        "invocation_is_symlink": True,
        "invocation_link_target": "python",
    }


def test_relocated_new_identity_binds_raw_and_resolved(tmp_path):
    _, freeze, campaign = portable_fixture(tmp_path)
    identity = interpreter_record()
    freeze["environment"]["interpreter_identity"] = identity
    for child in campaign["children"]:
        child["command"][0] = identity["invocation"]
    assert audit.validate_relocated_commands(campaign["children"], freeze) == identity


@pytest.mark.parametrize(
    "field",
    [
        "schema",
        "binary_sha256",
        "resolved_executable",
        "invocation",
        "invocation_link_target",
        "invocation_is_symlink",
    ],
)
def test_relocated_new_identity_rejects_malformed_or_mismatched(tmp_path, field):
    _, freeze, campaign = portable_fixture(tmp_path)
    identity = interpreter_record()
    freeze["environment"]["interpreter_identity"] = identity
    for child in campaign["children"]:
        child["command"][0] = identity["invocation"]
    bad = {
        "schema": "wrong",
        "binary_sha256": "bad",
        "resolved_executable": "/wrong",
        "invocation": "/wrong",
        "invocation_link_target": None,
        "invocation_is_symlink": 1,
    }
    identity[field] = bad[field]
    with pytest.raises(ValueError):
        audit.validate_relocated_commands(campaign["children"], freeze)


def test_legacy_alias_is_explicitly_inconclusive(tmp_path):
    _, freeze, campaign = portable_fixture(tmp_path)
    for child in campaign["children"]:
        child["command"][0] = "/old/venv/bin/python-link"
    with pytest.raises(ValueError, match="inconclusive_legacy_executable_mapping"):
        audit.validate_relocated_commands(campaign["children"], freeze)


def test_relocated_rejects_child_identity_mismatch(tmp_path, monkeypatch):
    manifest, freeze, campaign = portable_fixture(tmp_path)
    identity = interpreter_record()
    freeze["environment"]["interpreter_identity"] = identity
    freeze["environment"]["source_hashes"] = {
        "experiments/epw_restart_roundtrip.py": next(iter(audit.IDENTITY_PRODUCERS))
    }
    for child in campaign["children"]:
        child["command"][0] = identity["invocation"]
    case = tmp_path / "case-0/case.json"
    case.write_text(
        json.dumps({"interpreter_identity": {**identity, "binary_sha256": "c" * 64}})
    )
    campaign["children"][0]["receipt_sha256"] = audit.sha(case)
    freeze_path = tmp_path / "freeze.json"
    freeze_path.write_text(json.dumps(freeze))
    campaign["freeze_sha256"] = audit.sha(freeze_path)
    manifest.write_text(json.dumps(campaign))
    monkeypatch.setattr(audit, "provenance", lambda _: None)
    monkeypatch.setattr(audit, "native", lambda *_: {"frequency": [0.1]})
    with pytest.raises(ValueError, match="child interpreter identity mismatch"):
        audit.audit_campaign(manifest, native_root_override=tmp_path)


def test_present_null_identity_is_not_legacy(tmp_path):
    _, freeze, campaign = portable_fixture(tmp_path)
    freeze["environment"]["interpreter_identity"] = None
    with pytest.raises(ValueError, match="interpreter identity"):
        audit.validate_relocated_commands(campaign["children"], freeze)


def test_legacy_producer_cannot_claim_invented_interpreter_identity(
    tmp_path, monkeypatch
):
    manifest, freeze, campaign = portable_fixture(tmp_path)
    identity = interpreter_record()
    freeze["environment"]["interpreter_identity"] = identity
    freeze["environment"]["source_hashes"] = {
        "experiments/epw_restart_roundtrip.py": (
            "465e1cd646fbb4461985fa21a8aba9e75cddaff1ed6a1d26a027f6236e0f3046"
        )
    }
    for child in campaign["children"]:
        child["command"][0] = identity["invocation"]
    freeze_path = tmp_path / "freeze.json"
    freeze_path.write_text(json.dumps(freeze))
    campaign["freeze_sha256"] = audit.sha(freeze_path)
    manifest.write_text(json.dumps(campaign))
    monkeypatch.setattr(audit, "provenance", lambda _: None)
    monkeypatch.setattr(audit, "native", lambda *_: pytest.fail("native accessed"))
    for override in (None, tmp_path):
        with pytest.raises(ValueError, match="producer does not support"):
            audit.audit_campaign(manifest, native_root_override=override)


@pytest.mark.parametrize("relocated", [False, True])
def test_identity_capable_source_cannot_downgrade_to_legacy(
    tmp_path, monkeypatch, relocated
):
    manifest, freeze, campaign = portable_fixture(tmp_path)
    freeze["environment"]["source_hashes"]["experiments/epw_restart_roundtrip.py"] = (
        next(iter(audit.IDENTITY_PRODUCERS))
    )
    freeze_path = tmp_path / "freeze.json"
    freeze_path.write_text(json.dumps(freeze))
    campaign["freeze_sha256"] = audit.sha(freeze_path)
    manifest.write_text(json.dumps(campaign))
    monkeypatch.setattr(audit, "provenance", lambda _: None)
    monkeypatch.setattr(audit, "native", lambda *_: pytest.fail("native accessed"))
    with pytest.raises(ValueError, match="identity-capable producer omitted"):
        audit.audit_campaign(
            manifest, native_root_override=tmp_path if relocated else None
        )
