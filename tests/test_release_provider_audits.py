"""Mocked rejection tests; these do not establish actual provider acceptance."""

import json
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_provider_receipts import _receipt

from scripts import release_provider_audits as audit

SHA = "a" * 40


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def setup(tmp_path, monkeypatch):
    root = tmp_path / "artifacts"
    for lane, provider in audit.LANES.items():
        receipt = _receipt(provider)
        folder = root / lane
        if provider == "green":
            for name in receipt["artifacts"]:
                path = folder / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(name.encode())
                receipt["artifacts"][name] = audit.digest(path)
        if provider == "epw":
            receipt["artifact_directory"] = "epw-native"
            receipt["build_provenance_sha256"] = {}
            base = folder / "epw-native/test-suite/epw_super"
            base.mkdir(parents=True)
            for i, item in enumerate(receipt["native_executions"]):
                for kind in ("input", "log"):
                    name = f"{kind}{i}"
                    (base / name).write_text(name)
                    item[kind] = name
                    item[kind + "_sha256"] = audit.digest(base / name)
            gap = base / "gap"
            gap.write_text("gap")
            receipt["calculation"].update(
                gap_file=str(gap.relative_to(folder)), gap_sha256=audit.digest(gap)
            )
            manifest = folder / "epw-native/kernelgecp-source-manifest.json"
            manifest.write_text("manifest")
            receipt["native_source_manifest_sha256"] = audit.digest(manifest)
        write(folder / (provider + ".json"), receipt)
        if provider == "triqs":
            write(
                root / lane / "interacting-dmft/report.json",
                {
                    "schema": "kernelgecp.interacting-dmft.v1",
                    "pilot": False,
                    "all_mandatory_passed": True,
                    "source_revision": SHA,
                    "references": {str(i): {} for i in range(3)},
                    "validation": [{}] * 27,
                    "dense_validation": [{}] * 9,
                },
            )

    @contextmanager
    def snapshot(repository, revision):
        yield tmp_path

    monkeypatch.setattr(audit, "source_snapshot", snapshot)
    monkeypatch.setattr(
        audit.subprocess,
        "check_output",
        lambda cmd, **kw: SHA if "rev-parse" in cmd else "",
    )

    def run(cmd, **kwargs):
        report = Path(cmd[5]) / "report.json"
        assert cmd[1] == "-I"
        assert "PYTHONPATH" not in kwargs["env"]
        write(
            Path(cmd[cmd.index("--output") + 1]),
            {
                "schema": "kernelgecp.interacting-dmft-audit.v1",
                "passed": True,
                "source_revision": SHA,
                "campaign_report_sha256": audit.digest(report),
                "audited_records": {"compressed": 27, "dense": 9},
                "checks": dict.fromkeys(
                    (
                        "source_bindings",
                        "selection_and_environment",
                        "reference",
                        "records",
                        "repetitions",
                        "finite_pole_oracle",
                    ),
                    True,
                ),
            },
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    def invoke(runner=run, revision=SHA):
        return audit.audit(
            root,
            revision,
            Path("python"),
            tmp_path / "out.json",
            runner=runner,
            repository=tmp_path,
        )

    return root, invoke, run


def test_mock_success(setup):
    _, invoke, _ = setup
    result = invoke()
    assert result["status"] == "passed"
    for lane, provider in audit.LANES.items():
        assert result["lanes"][lane]["receipt_sha256"] == audit.digest(
            setup[0] / lane / (provider + ".json")
        )
    green = result["lanes"]["provider-green-linux"]["bound_payload_sha256"]
    assert len(green) == 6
    assert green["green/input.log"] == audit.digest(
        setup[0] / "provider-green-linux/green/input.log"
    )


@pytest.mark.parametrize(
    "case",
    [
        "missing",
        "symlink",
        "revision",
        "invalid",
        "noncanonical",
        "missing_report",
        "failed",
        "runtime",
        "array_report",
        "array_output",
        "tampered",
        "hardlink",
    ],
)
def test_rejections(setup, monkeypatch, tmp_path, case):
    root, invoke, run = setup
    lane = root / "provider-triqs-ubuntu-24.04"
    if case == "missing":
        (root / "provider-gp-macos-15/gp.json").unlink()
    elif case == "symlink":
        (root / "provider-gp-macos-15/gp.json").unlink()
        (root / "provider-gp-macos-15/gp.json").symlink_to(lane / "triqs.json")
    elif case == "invalid":

        def reject(*args):
            raise ValueError("invalid receipt")

        monkeypatch.setattr(audit, "validate_receipt", reject)
    elif case == "array_report":
        write(lane / "interacting-dmft/report.json", [])
    elif case == "array_output":

        def run(cmd, **kwargs):
            write(Path(cmd[cmd.index("--output") + 1]), [])
            return SimpleNamespace(returncode=0, stdout="", stderr="")
    elif case == "tampered":
        (root / "provider-green-linux/green/input.log").write_text("tampered")
    elif case == "hardlink":
        target = root / "provider-gp-macos-15/gp.json"
        (tmp_path / "hardlink").hardlink_to(target)
    elif case == "noncanonical":
        write(lane / "interacting-dmft/report.json", {"pilot": True})
    elif case == "missing_report":
        (lane / "interacting-dmft/report.json").unlink()
    elif case == "failed":

        def run(*a, **k):
            return SimpleNamespace(returncode=1, stdout="", stderr="failure")
    elif case == "runtime":

        def run(*args, **kwargs):
            raise FileNotFoundError("runtime missing")

    result = invoke(run, "b" * 40 if case == "revision" else SHA)
    assert result["status"] == "failed"
    assert json.loads((tmp_path / "out.json").read_text())["error"]


def test_hash_and_path_rejections(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (root / "file").write_text("correct")
    receipt = {"artifacts": {"file": "0" * 64}}
    with pytest.raises(ValueError, match="hash mismatch"):
        audit.check_bindings(root, receipt, "green")
    for path in ("../escape", "/absolute", "a/../file"):
        with pytest.raises(ValueError, match="unsafe"):
            audit.safe(root, path)
    (root / "link").symlink_to(root / "file")
    with pytest.raises(ValueError, match="symlink"):
        audit.safe(root, "link")


def test_snapshot_excludes_shadow_modules(tmp_path):
    repository = tmp_path / "repository"
    repository.mkdir()

    def git(*args):
        return subprocess.check_output(
            ["git", *args], cwd=repository, text=True
        ).strip()

    git("init", "-q")
    git("config", "user.email", "test@example.com")
    git("config", "user.name", "Test")
    for name in ("numpy", "mpmath"):
        (repository / name).mkdir()
        (repository / name / "__init__.py").write_text("TRUSTED = True\n")
    git("add", ".")
    git("commit", "-qm", "fixture")
    revision = git("rev-parse", "HEAD")
    for name in ("numpy", "mpmath"):
        (repository / (name + ".py")).write_text("raise RuntimeError('shadow')\n")
    with audit.source_snapshot(repository, revision) as snapshot:
        code = (
            "import sys; sys.path.insert(0,sys.argv[1]); import numpy,mpmath; "
            "assert numpy.TRUSTED and mpmath.TRUSTED"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", code, str(snapshot)],
            cwd=repository,
            check=True,
        )
