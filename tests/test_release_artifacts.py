"""Publication gates reject stale producers, altered payloads and public conflicts."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import subprocess
import urllib.error
import zipfile
from pathlib import Path
from typing import Any

import pytest

from scripts import release_artifacts as release


def _git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def _run(role: str, revision: str, run_id: int) -> dict[str, Any]:
    return {
        "id": run_id,
        "head_sha": revision,
        "head_branch": "main",
        "event": "push",
        "status": "completed",
        "conclusion": "success",
        "run_attempt": 1,
        "path": f".github/workflows/{'ci' if role == 'ci' else 'native'}.yml",
        "repository": {"full_name": release.REPOSITORY},
        "head_repository": {"full_name": release.REPOSITORY},
        "html_url": f"https://github.com/{release.REPOSITORY}/actions/runs/{run_id}",
    }


def _producers(revision: str = "a" * 40) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    runs = {"ci": [_run("ci", revision, 10)], "native": [_run("native", revision, 20)]}
    artifacts = [
        {
            "id": index + 30,
            "name": name,
            "expired": False,
            "workflow_run": {"id": 20, "head_sha": revision, "head_branch": "main"},
        }
        for index, name in enumerate(sorted(release.ARTIFACTS))
    ]
    return runs, artifacts


def test_selection_accepts_only_complete_exact_main_receipts() -> None:
    runs, artifacts = _producers()
    selected = release.select_runs(release.REPOSITORY, "a" * 40, runs, artifacts)
    assert selected["ci_run"]["id"] == 10 and selected["native_run"]["id"] == 20
    assert set(selected["artifact_ids"]) == release.ARTIFACTS


@pytest.mark.parametrize("role", ["ci", "native"])
@pytest.mark.parametrize(
    "field,value",
    [
        ("head_sha", "b" * 40),
        ("head_branch", "feature"),
        ("event", "pull_request"),
        ("status", "in_progress"),
        ("conclusion", "failure"),
        ("path", ".github/workflows/other.yml"),
        ("id", True),
        ("run_attempt", 0),
        ("repository", {"full_name": "other/repository"}),
        ("head_repository", {"full_name": "fork/repository"}),
    ],
)
def test_selection_rejects_wrong_or_unfinished_producer(
    role: str, field: str, value: Any
) -> None:
    runs, artifacts = _producers()
    runs[role][0][field] = value
    with pytest.raises(release.ReleaseError):
        release.select_runs(release.REPOSITORY, "a" * 40, runs, artifacts)


def test_later_failed_main_run_is_not_hidden_by_an_earlier_success() -> None:
    runs, artifacts = _producers()
    failed = _run("ci", "a" * 40, 11)
    failed["conclusion"] = "failure"
    runs["ci"].append(failed)
    with pytest.raises(release.ReleaseError, match="latest"):
        release.select_runs(release.REPOSITORY, "a" * 40, runs, artifacts)


@pytest.mark.parametrize(
    "mutation",
    [
        "expired",
        "wrong_run",
        "wrong_sha",
        "wrong_branch",
        "bool_id",
        "missing",
        "extra",
        "duplicate",
    ],
)
def test_selection_rejects_stale_or_ambiguous_artifacts(mutation: str) -> None:
    runs, artifacts = _producers()
    if mutation == "expired":
        artifacts[0]["expired"] = True
    elif mutation == "wrong_run":
        artifacts[0]["workflow_run"]["id"] = 21
    elif mutation == "wrong_sha":
        artifacts[0]["workflow_run"]["head_sha"] = "b" * 40
    elif mutation == "wrong_branch":
        artifacts[0]["workflow_run"]["head_branch"] = "feature"
    elif mutation == "bool_id":
        artifacts[0]["id"] = True
    elif mutation == "missing":
        artifacts.pop()
    elif mutation == "extra":
        artifacts.append(copy.deepcopy(artifacts[0]))
    else:
        artifacts[0]["name"] = artifacts[1]["name"]
    with pytest.raises(release.ReleaseError):
        release.select_runs(release.REPOSITORY, "a" * 40, runs, artifacts)


@pytest.fixture
def project(tmp_path: Path) -> Path:
    root = tmp_path / "project"
    (root / "native").mkdir(parents=True)
    (root / "kernelgecp").mkdir()
    (root / "pyproject.toml").write_text(
        '[project]\nname="kernelgecp"\nversion="1.1.0"\n'
    )
    (root / "native/Cargo.toml").write_text(
        '[package]\nname="kernelgecp-native"\nversion="1.1.0"\n'
    )
    (root / "kernelgecp/__init__.py").write_text('__version__ = "1.1.0"\n')
    for name in ("METHODS.md", "CITATION.cff", "RELEASE_NOTES.md"):
        (root / name).write_text(f"Fixture {name}\n")
    _git(root, "init", "-q", "-b", "main")
    _git(root, "add", ".")
    _git(
        root,
        "-c",
        "user.name=Release fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "-c",
        "commit.gpgsign=false",
        "commit",
        "-qm",
        "fixture",
    )
    _git(root, "update-ref", "refs/remotes/origin/main", "HEAD")
    _git(root, "tag", "v1.1.0")
    return root


def test_tag_selection_checks_real_git_identity_before_github(
    project: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    revision = _git(project, "rev-parse", "HEAD")
    runs, artifacts = _producers(revision)

    def api(path: str, *query: str) -> dict[str, Any]:
        if path.endswith("/artifacts"):
            return {"artifacts": artifacts}
        return {"workflow_runs": runs["ci" if "ci.yml" in path else "native"]}

    monkeypatch.setattr(release, "_gh", api)
    selected = release.select_release(
        project, "refs/tags/v1.1.0", revision, release.REPOSITORY
    )
    assert selected["source_revision"] == revision and selected["version"] == "1.1.0"


@pytest.mark.parametrize(
    "ref", ["refs/heads/main", "refs/tags/v01.1.0", "refs/tags/v1.1.0rc1"]
)
def test_nonstable_refs_never_query_github(
    project: Path, ref: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        release, "_gh", lambda *args: pytest.fail("must reject before network")
    )
    with pytest.raises(release.ReleaseError):
        release.select_release(
            project, ref, _git(project, "rev-parse", "HEAD"), release.REPOSITORY
        )


@pytest.mark.parametrize("bad", ["event", "cargo", "runtime", "outside_main"])
def test_tag_version_commit_and_ancestry_mismatch_stop_before_network(
    project: Path, bad: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    revision = _git(project, "rev-parse", "HEAD")
    if bad in ("cargo", "runtime", "outside_main"):
        file = project / (
            "native/Cargo.toml" if bad == "cargo" else "kernelgecp/__init__.py"
        )
        file.write_text(file.read_text().replace("1.1.0", "1.2.0"))
        _git(project, "add", ".")
        _git(
            project,
            "-c",
            "user.name=Fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "-c",
            "commit.gpgsign=false",
            "commit",
            "-qm",
            "changed",
        )
        _git(project, "tag", "-f", "v1.1.0")
        revision = _git(project, "rev-parse", "HEAD")
        if bad != "outside_main":
            _git(project, "update-ref", "refs/remotes/origin/main", revision)
    monkeypatch.setattr(
        release, "_gh", lambda *args: pytest.fail("must reject before network")
    )
    with pytest.raises(release.ReleaseError):
        release.select_release(
            project,
            "refs/tags/v1.1.0",
            "b" * 40 if bad == "event" else revision,
            release.REPOSITORY,
        )


@pytest.fixture
def collected(
    tmp_path: Path, project: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, dict[str, Any], Path]:
    """Exercise composition; complete archive inspectors have separate fixtures."""
    revision = _git(project, "rev-parse", "HEAD")
    runs, items = _producers(revision)
    selected = {
        "schema": "kernelgecp.selected-release.v1",
        "repository": release.REPOSITORY,
        "version": "1.1.0",
        "tag": "v1.1.0",
        "source_revision": revision,
        "source_tree": _git(project, "rev-parse", "HEAD^{tree}"),
        **release.select_runs(release.REPOSITORY, revision, runs, items),
    }
    artifacts = tmp_path / "artifacts"
    base = artifacts / "kernelgecp-release-source"
    for folder in ("source", "dist", "validation/sdist-wheel"):
        (base / folder).mkdir(parents=True)
    report = {
        "schema": "kernelgecp.corresponding-source.receipt.v1",
        "version": "1.1.0",
        "source_revision": revision,
        "source_tree": selected["source_tree"],
        "combined_native_license_expression": "Apache-2.0 AND LGPL-3.0-or-later",
    }
    for key, filename, payload in (
        (
            "archive",
            "kernelgecp-1.1.0-native-source.tar.gz",
            b"controlled archive fixture",
        ),
        ("license_bundle", "kernelgecp-1.1.0-native-licenses.txt", b"fixture licenses"),
        ("sbom", "kernelgecp-1.1.0-native-dependencies.cyclonedx.json", b"{}"),
    ):
        (base / "source" / filename).write_bytes(payload)
        report[key] = filename
        report["sha256" if key == "archive" else key + "_sha256"] = hashlib.sha256(
            payload
        ).hexdigest()
    (base / "source/source-report.json").write_text(json.dumps(report))
    for path in (
        base / "dist/kernelgecp-1.1.0-py3-none-any.whl",
        base / "dist/kernelgecp-1.1.0.tar.gz",
        base / "validation/sdist-wheel/kernelgecp-1.1.0-py3-none-any.whl",
    ):
        path.write_bytes(b"controlled portable fixture")
    inventory = {
        "source_revision": revision,
        "source_tree": selected["source_tree"],
        "source_asset": {"filename": report["archive"], "sha256": report["sha256"]},
        "license_bundle_sha256": report["license_bundle_sha256"],
        "sbom_sha256": report["sbom_sha256"],
    }
    for role, tag in (
        ("linux", "manylinux_2_28_x86_64"),
        ("macos", "macosx_11_0_arm64"),
    ):
        folder = artifacts / f"kernelgecp-release-native-{role}"
        folder.mkdir()
        with zipfile.ZipFile(
            folder / f"kernelgecp-1.1.0-cp312-abi3-{tag}.whl", "w"
        ) as archive:
            archive.writestr("kernelgecp/_native_inventory.json", json.dumps(inventory))
    monkeypatch.setattr(release, "check_archives", lambda **kwargs: {})
    monkeypatch.setattr(
        release, "inspect_corresponding_source", lambda *args, **kwargs: report
    )
    return artifacts, selected, tmp_path / "release"


def _collect(
    project: Path, collected: tuple[Path, dict[str, Any], Path]
) -> tuple[Path, dict[str, Any]]:
    artifacts, selected, output = collected
    manifest = release.collect_artifacts(
        project_root=project, selected=selected, artifacts=artifacts, output=output
    )
    return output, manifest


def test_collector_separates_four_uploads_from_rebuild_and_source_assets(
    project: Path, collected: tuple[Path, dict[str, Any], Path]
) -> None:
    output, manifest = _collect(project, collected)
    assert len(manifest["dist"]) == 4 and len(list((output / "dist").iterdir())) == 4
    assert not (output / "dist/kernelgecp-1.1.0-native-source.tar.gz").exists()
    assert release.verify_payload(output) == manifest


@pytest.mark.parametrize(
    "bad", ["extra", "wrong_target", "source_revision", "source_hash", "symlink"]
)
def test_collector_rejects_ambiguous_stale_or_changed_downloads(
    project: Path, collected: tuple[Path, dict[str, Any], Path], bad: str
) -> None:
    artifacts, selected, output = collected
    if bad == "extra":
        (artifacts / "unexpected").mkdir()
    elif bad == "wrong_target":
        next((artifacts / "kernelgecp-release-native-linux").iterdir()).rename(
            artifacts / "kernelgecp-release-native-linux/wrong.whl"
        )
    elif bad == "source_revision":
        selected["source_revision"] = "b" * 40
    elif bad == "source_hash":
        (
            artifacts
            / "kernelgecp-release-source/source/kernelgecp-1.1.0-native-source.tar.gz"
        ).write_bytes(b"tampered")
    else:
        target = artifacts / "kernelgecp-release-source/dist/extra.whl"
        try:
            target.symlink_to(project / "pyproject.toml")
        except OSError:
            pytest.skip("symbolic links unavailable")
    with pytest.raises(ValueError):
        release.collect_artifacts(
            project_root=project, selected=selected, artifacts=artifacts, output=output
        )
    assert not output.exists()


@pytest.mark.parametrize(
    "bad", ["hash", "extra", "missing", "path", "duplicate", "roles", "identity"]
)
def test_frozen_payload_rejects_changes(
    project: Path, collected: tuple[Path, dict[str, Any], Path], bad: str
) -> None:
    output, manifest = _collect(project, collected)
    candidate = output / "dist" / manifest["dist"][0]["filename"]
    if bad == "hash":
        candidate.write_bytes(b"tampered")
    elif bad == "extra":
        (output / "extra").write_bytes(b"unexpected")
    elif bad == "missing":
        candidate.unlink()
    elif bad == "path":
        manifest["dist"][0]["filename"] = "../outside"
    elif bad == "duplicate":
        manifest["dist"].append(copy.deepcopy(manifest["dist"][0]))
    elif bad == "roles":
        manifest["dist"][0]["role"] = "sdist"
    else:
        manifest["tag"] = "v2.0.0"
    (output / "release-manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        release.verify_payload(output)


@pytest.mark.parametrize("bad", [None, "hash", "filename", "missing"])
def test_public_pypi_requires_exact_filename_hash_set(
    project: Path,
    collected: tuple[Path, dict[str, Any], Path],
    monkeypatch: pytest.MonkeyPatch,
    bad: str | None,
) -> None:
    output, manifest = _collect(project, collected)
    urls = [
        {"filename": x["filename"], "digests": {"sha256": x["sha256"]}}
        for x in manifest["dist"]
    ]
    if bad == "hash":
        urls[0]["digests"]["sha256"] = "b" * 64
    elif bad == "filename":
        urls[0]["filename"] = "other.whl"
    elif bad == "missing":
        urls.pop()
    payload = {"info": {"name": "kernelgecp", "version": "1.1.0"}, "urls": urls}
    monkeypatch.setattr(
        release.urllib.request,
        "urlopen",
        lambda *args, **kwargs: io.BytesIO(json.dumps(payload).encode()),
    )
    if bad is None:
        assert release.verify_pypi(output)["status"] == "verified"
    else:
        with pytest.raises(release.ReleaseError):
            release.verify_pypi(output)


def test_public_absence_is_bounded_and_never_success(
    project: Path,
    collected: tuple[Path, dict[str, Any], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output, _ = _collect(project, collected)
    calls = []

    def missing(url: str, **kwargs: Any) -> Any:
        calls.append(url)
        raise urllib.error.HTTPError(url, 404, "not found", {}, None)

    monkeypatch.setattr(release.urllib.request, "urlopen", missing)
    delays = []
    monkeypatch.setattr(release.time, "sleep", delays.append)
    with pytest.raises(release.ReleaseError):
        release.verify_pypi(output)
    assert len(calls) == 6
    assert len(delays) == 5 and max(delays) <= 60
    assert 100 < sum(delays) <= 135


def test_public_metadata_visibility_recovers_without_relaxing_hash_checks(
    project: Path,
    collected: tuple[Path, dict[str, Any], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output, manifest = _collect(project, collected)
    payload = {
        "info": {"name": "kernelgecp", "version": "1.1.0"},
        "urls": [
            {"filename": item["filename"], "digests": {"sha256": item["sha256"]}}
            for item in manifest["dist"]
        ],
    }
    calls = []
    delays = []

    def eventually_visible(url: str, **kwargs: Any) -> Any:
        calls.append(url)
        if len(calls) <= 2:
            raise urllib.error.HTTPError(url, 404, "not yet visible", {}, None)
        return io.BytesIO(json.dumps(payload).encode())

    monkeypatch.setattr(release.urllib.request, "urlopen", eventually_visible)
    monkeypatch.setattr(release.time, "sleep", delays.append)
    result = release.verify_pypi(output)
    assert result["status"] == "verified"
    assert result["files"] == {
        item["filename"]: item["sha256"] for item in manifest["dist"]
    }
    assert len(calls) == 3 and len(delays) == 2
    assert 0 < delays[0] < delays[1] <= 60
