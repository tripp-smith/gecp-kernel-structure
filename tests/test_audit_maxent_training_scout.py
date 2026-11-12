import hashlib

import pytest

from scripts import audit_maxent_training_scout as audit


def test_header_rejects_bad_schema_source_and_controls():
    revision = "a" * 40
    valid = {
        "schema": "kernelgecp.maxent-training-scout.v4",
        "source_revision": revision,
        "identity": {"source_revision": revision},
        "fixed_controls": audit.CONTROLS,
    }
    audit.validate_receipt_header(valid, revision)
    for key, value in [
        ("schema", "v3"),
        ("source_revision", "b" * 40),
        ("fixed_controls", {}),
    ]:
        with pytest.raises(ValueError):
            audit.validate_receipt_header({**valid, key: value}, revision)


def test_authenticated_inventory_rejects_missing_or_changed_files(tmp_path):
    hashes = {}
    for name in (
        "worker.log",
        *(
            f"{d}-{m}-trials.npz"
            for d in ("identity", "marquardt")
            for m in audit.MODELS
        ),
    ):
        (tmp_path / name).write_bytes(name.encode())
        hashes[name] = hashlib.sha256(name.encode()).hexdigest()
    audit.validate_outputs(tmp_path, hashes)
    with pytest.raises(ValueError, match="inventory"):
        audit.validate_outputs(tmp_path, {})
    (tmp_path / "worker.log").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="hash"):
        audit.validate_outputs(tmp_path, hashes)


def test_rejects_wrong_minimizer_flags_and_selector():
    native = {"marquardt": False, "J_squared": False, "mu0": 1e-8}
    audit.validate_minimizer_flags({"mode": 1}, {"mode": 1}, native, "identity")
    with pytest.raises(ValueError, match="flags"):
        audit.validate_minimizer_flags({"mode": 1}, {"mode": 1}, native, "marquardt")
    flags = [True] * 33
    audit.validate_selection(flags, 16, 16)
    with pytest.raises(ValueError, match="selector"):
        audit.validate_selection(flags, 16, 17)
    with pytest.raises(ValueError, match="flags"):
        audit.validate_selection([1] * 33, 16, 16)
    flags[16] = False
    with pytest.raises(ValueError, match="convergence"):
        audit.validate_selection(flags, 16, 16)


def test_bad_receipt_hash_is_rejected(tmp_path):
    path = tmp_path / "receipt.json"
    path.write_text('{"schema":"v4"}')
    with pytest.raises(ValueError, match="authentication"):
        audit.authenticate_receipt(path, "0" * 64)
