import subprocess
from pathlib import Path

import pytest

from scripts import release_cthyb_campaign_scope as scope
from scripts import release_source_equivalence as eq

ROOT = Path(__file__).resolve().parents[1]


def git(repo, *args):
    return (
        subprocess.check_output(["git", "-C", str(repo), *args], stderr=subprocess.PIPE)
        .decode()
        .strip()
    )


@pytest.fixture
def release(tmp_path):
    git(tmp_path, "init", "-q")
    git(tmp_path, "config", "user.name", "Fixture")
    git(tmp_path, "config", "user.email", "fixture@example.invalid")
    for path in eq.VERSION_FILES:
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(eq.git(str(ROOT), "show", f"{eq.BASELINE}:{path}"))
    git(tmp_path, "add", ".")
    git(tmp_path, "commit", "-qm", "base")
    base = git(tmp_path, "rev-parse", "HEAD")
    for path in eq.VERSION_FILES:
        target = tmp_path / path
        target.write_bytes(
            eq.normalized(path, target.read_bytes(), eq.OLD).replace(
                b"<VERSION>", eq.NEW
            )
        )
    git(tmp_path, "add", ".")
    git(tmp_path, "commit", "-qm", "bump")
    return tmp_path, base, git(tmp_path, "rev-parse", "HEAD")


def test_version_candidate_runs_preflight_at_detached_base(release, monkeypatch):
    repo, base, candidate = release

    def check(root, *, force=False):
        assert root != repo and not force
        assert git(root, "rev-parse", "HEAD") == base
        assert git(root, "status", "--porcelain") == ""
        assert (root / "pyproject.toml").read_bytes() == eq.git(
            str(repo), "show", f"{base}:pyproject.toml"
        )
        return False

    monkeypatch.setattr(scope, "campaign_needed", check)
    assert scope.decide(repo, "pull_request", candidate, base) is False


def test_mixed_and_dirty_changes_use_candidate(release, monkeypatch):
    repo, base, candidate = release
    monkeypatch.setattr(scope, "campaign_needed", lambda root, **kw: root == repo)
    (repo / "extra.txt").write_text("change")
    assert scope.decide(repo, "pull_request", candidate, base) is True
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "mixed")
    assert (
        scope.decide(repo, "pull_request", git(repo, "rev-parse", "HEAD"), base) is True
    )
    assert scope.decide(repo, "pull_request", candidate, "0" * 40) is True


def test_dispatch_forces_original_preflight(release, monkeypatch):
    repo, base, candidate = release
    monkeypatch.setattr(
        scope, "campaign_needed", lambda root, *, force=False: root == repo and force
    )
    assert scope.decide(repo, "workflow_dispatch", candidate, base) is True


def test_missing_baseline_evidence_requires_campaign(release):
    repo, base, candidate = release
    assert scope.decide(repo, "pull_request", candidate, base) is True


def test_changed_baseline_evidence_fails_closed(release):
    repo, base, _ = release
    git(repo, "reset", "--hard", "-q", base)
    archive = repo / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    binding = repo / "data/cthyb/linux-audit-binding-7200ec0.json"
    archive.parent.mkdir(parents=True)
    binding.parent.mkdir(parents=True)
    archive.write_bytes(b"changed")
    binding.write_text("{}")
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "evidence")
    base = git(repo, "rev-parse", "HEAD")
    for path in eq.VERSION_FILES:
        target = repo / path
        target.write_bytes(
            eq.normalized(path, target.read_bytes(), eq.OLD).replace(
                b"<VERSION>", eq.NEW
            )
        )
    git(repo, "add", ".")
    git(repo, "commit", "-qm", "bump")
    with pytest.raises(ValueError, match="archive digest mismatch"):
        scope.decide(repo, "pull_request", git(repo, "rev-parse", "HEAD"), base)
