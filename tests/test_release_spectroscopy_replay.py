import hashlib
import json
from types import SimpleNamespace

import pytest

from scripts import release_spectroscopy_replay as replay


def test_pin_requires_exact_schema_and_types():
    valid = {"schema": replay.SCHEMA, "revision": "a" * 40, "run_id": 123}
    assert replay.identity(valid) == ("a" * 40, 123)
    for bad in (
        {**valid, "run_id": True},
        {**valid, "extra": 1},
        {**valid, "revision": "a" * 39},
        {**valid, "schema": "old"},
    ):
        with pytest.raises(ValueError):
            replay.identity(bad)


def test_replay_checks_archive_and_receipts_and_strips_token(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    output = tmp_path / "result"
    content = b"saved archive"
    digest = hashlib.sha256(content).hexdigest()
    origin = {
        "artifact": {
            "name": "maxent-fixed-bath-final-123-1",
            "id": 42,
            "digest": "sha256:" + digest,
            "size_in_bytes": len(content),
        }
    }
    monkeypatch.setattr(replay, "verify_checkout", lambda *a: None)
    monkeypatch.setattr(replay, "compare_candidate", lambda *a: {})
    monkeypatch.setattr(replay, "fetch_spectroscopy_origin", lambda *a: origin)
    monkeypatch.setenv("GH_TOKEN", "secret")
    monkeypatch.setenv("GITHUB_TOKEN", "secret")

    def download(ident, target, log):
        assert ident == 42
        target.write_bytes(content)
        log.write_text("")

    monkeypatch.setattr(replay.payloads, "download", download)

    def run(command, **kwargs):
        assert "GH_TOKEN" not in kwargs["env"] and "GITHUB_TOKEN" not in kwargs["env"]
        audit = output / "audit"
        audit.mkdir()
        hashes = {}
        for name in (
            "synthetic-cal",
            "synthetic-final",
            "measured-cal",
            "raw",
            "measured-final",
        ):
            path = audit / (name + ".json")
            path.write_text(name)
            hashes[name] = replay.sha256_file(path)
        receipt = {
            "schema": "kernelgecp.spectroscopy-v4-saved-replay.v1",
            "status": "passed",
            "source_revision": "a" * 40,
            "run_id": "123",
            "archive_sha256": digest,
            "archive_size": len(content),
            "counts": replay.COUNTS,
            "receipt_sha256": hashes,
        }
        (audit / "receipt.json").write_text(json.dumps(receipt))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(replay.subprocess, "run", run)
    assert (
        replay.replay(repo, "b" * 40, "a" * 40, 123, output, "python")["status"]
        == "passed"
    )
    origin["artifact"]["digest"] = "sha256:" + "0" * 64
    with pytest.raises(ValueError, match="digest or size"):
        replay.replay(repo, "b" * 40, "a" * 40, 123, tmp_path / "bad", "python")


def test_candidate_release_requires_one_parent_and_exact_transition(
    monkeypatch, tmp_path
):
    candidate, source, parent = "c" * 40, "a" * 40, "b" * 40
    monkeypatch.setattr(
        replay.science_equivalence,
        "compare",
        lambda repo, s, p: {"source": s, "producer": p},
    )
    monkeypatch.setattr(
        replay.version_equivalence,
        "compare",
        lambda repo, c, p: {"version_transition": "1.1.0 -> 1.2.0"},
    )

    def git(command, **kwargs):
        if "show" in command:
            return '[project]\nversion = "1.2.0"\n'
        return f"{candidate} {parent}\n"

    monkeypatch.setattr(replay.subprocess, "check_output", git)
    result = replay.compare_candidate(tmp_path, source, candidate)
    assert result["science_to_producer"] == {"source": source, "producer": parent}
    monkeypatch.setattr(
        replay.version_equivalence,
        "compare",
        lambda *a: {"version_transition": "1.1.0 -> 1.1.0"},
    )
    with pytest.raises(ValueError, match="transition"):
        replay.compare_candidate(tmp_path, source, candidate)
    monkeypatch.setattr(
        replay.subprocess,
        "check_output",
        lambda command, **kw: (
            '[project]\nversion = "1.2.0"\n'
            if "show" in command
            else f"{candidate} {parent} {'d' * 40}\n"
        ),
    )
    with pytest.raises(ValueError, match="one parent"):
        replay.compare_candidate(tmp_path, source, candidate)


def test_checkout_identity_and_cleanliness(monkeypatch, tmp_path):
    candidate = "a" * 40

    def git(command, **kwargs):
        return candidate + "\n" if "rev-parse" in command else ""

    monkeypatch.setattr(replay.subprocess, "check_output", git)
    replay.verify_checkout(tmp_path, candidate)
    with pytest.raises(ValueError, match="lowercase"):
        replay.verify_checkout(tmp_path, candidate.upper())
    monkeypatch.setattr(
        replay.subprocess,
        "check_output",
        lambda command, **kw: (
            candidate + "\n" if "rev-parse" in command else "?? unknown\n"
        ),
    )
    with pytest.raises(ValueError, match="clean"):
        replay.verify_checkout(tmp_path, candidate)


def test_published_descendant_tree_gate(tmp_path, monkeypatch):
    import subprocess

    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(tmp_path), *args], text=True
        ).strip()

    git("init", "-q")
    git("config", "user.email", "test@example.com")
    git("config", "user.name", "Test")
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "1.1.0"\n')
    (tmp_path / "science.py").write_text("original\n")
    git("add", ".")
    git("commit", "-qm", "producer")
    producer = git("rev-parse", "HEAD")
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "1.2.0"\n')
    git("commit", "-qam", "release")
    anchor = git("rev-parse", "HEAD")
    source = "a" * 40
    monkeypatch.setattr(replay, "PUBLISHED_RELEASE", anchor)
    monkeypatch.setattr(replay, "FROZEN_SCIENCE", source)

    def science(repo, s, p):
        if p != producer or s != source:
            raise ValueError("unexpected producer")
        return {"source": s, "producer": p}

    monkeypatch.setattr(replay.science_equivalence, "compare", science)
    monkeypatch.setattr(
        replay.version_equivalence,
        "compare",
        lambda repo, c, p: {"version_transition": "1.1.0 -> 1.2.0"},
    )
    (tmp_path / "PROJECT_COST.md").write_text("cost\n")
    git("add", ".")
    git("commit", "-qm", "documentation")
    candidate = git("rev-parse", "HEAD")
    result = replay.compare_candidate(tmp_path, source, candidate)
    provenance = result["published_descendant"]
    assert provenance["anchor_commit"] == anchor
    assert provenance["candidate_commit"] == candidate
    assert provenance["changed_paths"] == ["PROJECT_COST.md"]
    assert provenance["compared_path_count"] == 3
    assert (
        result["published_release_guard"]["science_to_producer"]["producer"] == producer
    )

    with pytest.raises(ValueError, match="frozen science source"):
        replay.compare_candidate(tmp_path, "b" * 40, candidate)

    (tmp_path / "science.py").write_text("modified\n")
    git("commit", "-qam", "science mutation")
    with pytest.raises(ValueError, match="science.py"):
        replay.compare_candidate(tmp_path, source, git("rev-parse", "HEAD"))

    git("checkout", "-q", producer)
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "1.2.0"\n')
    git("commit", "-qam", "unrelated release")
    unrelated = git("rev-parse", "HEAD")
    # Even an unrelated 1.2.0 commit cannot use the descendant exception.
    monkeypatch.setattr(
        replay.version_equivalence,
        "compare",
        lambda *a: {"version_transition": "wrong"},
    )
    with pytest.raises(ValueError, match="transition"):
        replay.compare_candidate(tmp_path, source, unrelated)
