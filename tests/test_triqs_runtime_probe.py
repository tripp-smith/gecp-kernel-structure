"""Failure/identity checks for the optional actual-provider runtime probe."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from scripts.providers import triqs_smoke as probe


def test_failed_check_is_recorded_before_rejection() -> None:
    checks: dict[str, bool] = {}
    with pytest.raises(RuntimeError, match="cthyb_sampling"):
        probe.require(checks, "cthyb_sampling", False)
    assert checks == {"cthyb_sampling": False}


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("python", "3.13.2"),
        ("triqs", "3.3.3"),
        ("cthyb", "4.0.2"),
        ("maxent", "4.0.1.dev1"),
        ("maxent", ""),
    ],
)
def test_runtime_requires_the_declared_release_tuple(name: str, value: str) -> None:
    versions = {"python": "3.12.16", **probe.EXPECTED_VERSIONS}
    versions[name] = value
    with pytest.raises(RuntimeError):
        probe.validate_versions(versions)


def test_runtime_accepts_the_declared_release_tuple() -> None:
    probe.validate_versions({"python": "3.12.16", **probe.EXPECTED_VERSIONS})


def test_environment_identity_excludes_paths_and_urls(tmp_path: Path) -> None:
    metadata = tmp_path / "conda-meta"
    metadata.mkdir()
    record = {
        "name": "triqs",
        "version": "4.0.2",
        "build": "py312_example_0",
        "subdir": "osx-arm64",
        "url": "https://private.invalid/token=secret",
        "prefix": str(tmp_path),
    }
    path = metadata / "triqs.json"
    path.write_text(json.dumps(record))
    first = probe.environment_identity(tmp_path)
    record["url"] = "https://another.invalid/other-token"
    record["prefix"] = "/another/private/path"
    path.write_text(json.dumps(record))
    second = probe.environment_identity(tmp_path)
    assert first == second
    assert first["conda_package_count"] == 1
    rendered = json.dumps(first)
    assert str(tmp_path) not in rendered
    assert "secret" not in rendered
    assert "https://" not in rendered
    record["build"] = "changed_build"
    path.write_text(json.dumps(record))
    assert (
        probe.environment_identity(tmp_path)["conda_package_set_sha256"]
        != first["conda_package_set_sha256"]
    )


def test_an_unidentified_runtime_cannot_pass(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="conda package identity"):
        probe.environment_identity(tmp_path)


def test_missing_provider_does_not_become_an_import_only_success(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        probe, "source_identity", lambda root: {"source_revision": "a" * 40}
    )
    monkeypatch.setattr(probe, "environment_identity", lambda prefix: {})
    for name in probe.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")

    def missing(name: str) -> Any:
        raise ModuleNotFoundError("deliberately unavailable test provider")

    monkeypatch.setattr(probe.importlib, "import_module", missing)
    receipt = probe.run_probe()
    assert receipt["schema"] == "kernelgecp.provider-smoke.v1"
    assert receipt["provider"] == "triqs"
    assert receipt["status"] == "failed"
    assert receipt["checks"] == {
        "imports": False,
        "numerical": False,
        "source_revision": True,
    }
    assert receipt["error"] == {"stage": "imports", "type": "ModuleNotFoundError"}


def test_cli_failure_has_nonzero_exit_and_a_failed_receipt(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    receipt = {
        "schema": probe.SCHEMA,
        "provider": "triqs",
        "status": "failed",
        "checks": {"imports": False, "numerical": False},
        "versions": {},
    }
    destination = tmp_path / "nested" / "triqs.json"
    monkeypatch.setattr(probe, "run_probe", lambda revision: receipt)
    monkeypatch.setattr(
        probe.sys, "argv", ["triqs_smoke", "--output", str(destination)]
    )
    assert probe.main() == 1
    assert json.loads(destination.read_text()) == receipt
    assert json.loads(capsys.readouterr().out) == receipt


@pytest.mark.parametrize("expected", ["bad-revision", "a" * 40])
def test_declared_source_must_match_actual_head(
    expected: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        probe, "source_identity", lambda root: {"source_revision": "b" * 40}
    )
    for name in probe.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    receipt = probe.run_probe(expected)
    assert receipt["status"] == "failed"
    assert receipt["checks"]["source_revision"] is False
    assert receipt["checks"]["imports"] is False
    assert receipt["error"]["stage"] == "source_revision"
