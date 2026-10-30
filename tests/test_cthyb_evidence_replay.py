"""Fail-closed checks for the frozen CT-HYB evidence replay."""

import importlib.util
import json
import zipfile
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "replay_cthyb_acceptance",
    Path(__file__).parents[1] / "experiments/replay_cthyb_acceptance.py",
)
assert spec and spec.loader
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def valid_audit():
    return {
        "status": "passed",
        "generation_revision": replay.GENERATION,
        "audit_sha256": replay.AUDITOR_SHA,
        "gates": dict.fromkeys(replay.GATES, True),
    }


def test_digest_mismatch_fails_before_providers(tmp_path):
    archive = tmp_path / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    archive.parent.mkdir(parents=True)
    archive.write_bytes(b"altered")
    with pytest.raises(ValueError, match="archive digest mismatch"):
        replay.replay(tmp_path, tmp_path / "output.json")


def test_unsafe_archive_fails(tmp_path):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("../escape", b"x")
        for i in range(466):
            z.writestr(f"{i}.txt", b"x")
    with pytest.raises(ValueError, match="unsafe archive"):
        replay.extract(archive, tmp_path / "out")
    assert not (tmp_path / "escape").exists()


@pytest.mark.parametrize("mode", ["missing", "failed"])
def test_gate_must_be_present_and_true(mode):
    audit = valid_audit()
    key = next(iter(replay.GATES))
    if mode == "missing":
        del audit["gates"][key]
    else:
        audit["gates"][key] = False
    with pytest.raises(ValueError, match="gate"):
        replay.check_gates(audit)


def test_unavailable_provider_fails(monkeypatch):
    def missing(_name):
        raise ModuleNotFoundError("provider missing")

    monkeypatch.setattr(replay.importlib, "import_module", missing)
    with pytest.raises(ModuleNotFoundError, match="provider missing"):
        replay.require_providers()


def test_changed_scientific_source_fails(tmp_path, monkeypatch):
    folder = tmp_path / "data"
    folder.mkdir()
    source = tmp_path / "source.py"
    source.write_text("changed")
    monkeypatch.setattr(replay, "SOURCES", {"source.py"})
    manifest = {
        "source_revision": replay.GENERATION,
        "sources": {"source.py": "0" * 64},
    }
    (folder / "manifest.json").write_text(json.dumps(manifest))
    binding = {"generation_revision": replay.GENERATION}
    for key, name in (
        ("manifest_sha256", "manifest.json"),
        ("freeze_sha256", "freeze.json"),
        ("result_sha256", "result.json"),
        ("audit_receipt_sha256", "independent-audit.json"),
    ):
        path = folder / name
        if not path.exists():
            path.write_text("{}")
        binding[key] = replay.digest(path)
    with pytest.raises(ValueError, match="scientific source changed"):
        replay.verify_tree(folder, binding, tmp_path)


def test_scope_accepted_docs_head_skips_campaign(tmp_path, monkeypatch):
    _scope_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(replay, "sources_unchanged", lambda _root, _sources: True)
    assert replay.campaign_needed(tmp_path) is False


def test_scope_changed_science_requires_campaign(tmp_path, monkeypatch):
    _scope_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(replay, "sources_unchanged", lambda _root, _sources: False)
    assert replay.campaign_needed(tmp_path) is True


def test_scope_forced_dispatch_requires_campaign(tmp_path):
    assert replay.campaign_needed(tmp_path, force=True) is True


def test_scope_corrupt_archive_fails_closed(tmp_path):
    archive = tmp_path / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    binding = tmp_path / "data/cthyb/linux-audit-binding-7200ec0.json"
    archive.parent.mkdir(parents=True)
    binding.parent.mkdir(parents=True)
    archive.write_bytes(b"corrupt")
    binding.write_text("{}")
    with pytest.raises(ValueError, match="archive digest mismatch"):
        replay.campaign_needed(tmp_path)


def test_scope_missing_baseline_requires_campaign(tmp_path):
    assert replay.campaign_needed(tmp_path) is True


def _scope_fixture(tmp_path, monkeypatch):
    archive = tmp_path / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    binding = tmp_path / "data/cthyb/linux-audit-binding-7200ec0.json"
    archive.parent.mkdir(parents=True)
    binding.parent.mkdir(parents=True)
    archive.write_bytes(b"fixture")
    binding.write_text("{}")
    monkeypatch.setattr(replay, "ZIP_SHA", replay.digest(archive))
    monkeypatch.setattr(replay, "BINDING_SHA", replay.digest(binding))

    def extract(_archive, dest):
        (dest / "manifest.json").write_text('{"sources": {}}')

    monkeypatch.setattr(replay, "extract", extract)
    monkeypatch.setattr(replay, "verify_tree", lambda *_args, **_kwargs: None)
