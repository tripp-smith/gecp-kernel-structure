"""Structural fixtures use synthetic bytes; no scientific acceptance is claimed."""

import hashlib
import json
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import release_molecular_audits as module

SHA = "a" * 40


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data if isinstance(data, bytes) else json.dumps(data).encode())
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(root):
    cases = []
    references = {}
    for bond in (0.74, 1.5):
        for grid, label in ((1e4, "1e4"), (1e5, "1e5")):
            name = f"h2-{bond}-{label}"
            folder = root / "native" / name
            data_hash = put(folder / "input.h5", b"synthetic")
            provider_hash = put(
                folder / "provider-receipt.json",
                {
                    "schema": "kernelgecp.provider-smoke.v1",
                    "provider": "green",
                    "status": "passed",
                    "source_revision": SHA,
                },
            )
            manifest_hash = put(
                folder / "bound-manifest.json",
                {
                    "source_revision": SHA,
                    "files": {"input": {"path": "input.h5", "sha256": data_hash}},
                },
            )
            cases.append(
                {
                    "directory": name,
                    "bond_length_angstrom": bond,
                    "native_grid_lambda": grid,
                    "manifest_sha256": manifest_hash,
                    "provider_receipt_sha256": provider_hash,
                }
            )
            references[name] = {
                key: put(root / "comparison" / name / key, b"synthetic")
                for key in ("calibration.npz", "validation.npz")
            }
    native_hash = put(
        root / "native" / "native-campaign.json",
        {
            "schema": "kernelgecp.green-molecular-native.v1",
            "status": "passed",
            "source_revision": SHA,
            "cases": cases,
        },
    )
    selection_hash = put(root / "comparison" / "selection.json", [])
    put(
        root / "comparison" / "report.json",
        {
            "schema": "kernelgecp.molecular-green.v2",
            "pilot": False,
            "source_revision": SHA,
            "native_source_revision": SHA,
            "native_campaign_sha256": native_hash,
            "selection_sha256": selection_hash,
            "reference_artifacts": references,
            "calibration": [],
            "cutoff_refinement_attempts": [],
            "validation": [],
        },
    )


@contextmanager
def snapshot(repo, revision):
    assert revision == SHA
    yield Path(repo)


def runner(command, **kwargs):
    assert command[1] == "-I" and "PYTHONPATH" not in kwargs["env"]
    result = {
        "schema": "kernelgecp.molecular-green-audit.v1",
        "status": "passed",
        "generation_revision": SHA,
        "native_source_revision": SHA,
        "validated_identity_count": 36,
        "gates": dict.fromkeys(module.GATES, True),
    }
    for name, n in (
        ("final_trials", 36),
        ("native", 4),
        ("native_grid_refinement", 2),
        ("cutoff_controls", 12),
        ("first_passing_calibration", 12),
        ("model_differences", 4),
    ):
        result[name] = [{}] * n
    put(Path(command[-1]), result)
    return SimpleNamespace(returncode=0, stdout="test", stderr="")


def execute(tmp_path, run=runner):
    return module.audit(
        tmp_path / "payload",
        SHA,
        "/python",
        tmp_path / "out",
        repository=tmp_path,
        snapshot=snapshot,
        runner=run,
    )


def test_success_binds_synthetic_payload_and_isolated_auditor(tmp_path):
    fixture(tmp_path / "payload")
    result = execute(tmp_path)
    assert result["status"] == "passed"
    assert (tmp_path / "out" / "auditor.log").exists()
    with pytest.raises(ValueError, match="output must be new"):
        execute(tmp_path)


@pytest.mark.parametrize(
    "damage",
    [
        "orphan",
        "link",
        "hardlink",
        "tamper",
        "wrong_source",
        "bad_audit",
        "failed_process",
    ],
)
def test_damage_fails_closed_and_preserves_attempt(tmp_path, damage):
    root = tmp_path / "payload"
    fixture(root)
    if damage == "orphan":
        put(root / "comparison" / "orphan", b"x")
    elif damage == "link":
        (root / "comparison" / "link").symlink_to(
            root / "comparison" / "selection.json"
        )
    elif damage == "hardlink":
        (root / "comparison" / "hard").hardlink_to(
            root / "comparison" / "selection.json"
        )
    elif damage == "tamper":
        put(root / "native" / "h2-0.74-1e4" / "input.h5", b"altered")
    elif damage == "wrong_source":
        path = root / "comparison" / "report.json"
        report = json.loads(path.read_text())
        report["source_revision"] = "b" * 40
        put(path, report)

    def bad_runner(command, **kwargs):
        result = runner(command, **kwargs)
        if damage == "failed_process":
            return SimpleNamespace(returncode=1, stdout="failed", stderr="detail")
        path = Path(command[-1])
        data = json.loads(path.read_text())
        data["gates"]["representation"] = False
        put(path, data)
        return result

    result = execute(
        tmp_path, bad_runner if damage in ("bad_audit", "failed_process") else runner
    )
    assert result["status"] == "failed"
    assert json.loads((tmp_path / "out" / "receipt.json").read_text())["error"]


@pytest.mark.parametrize("bypass", [False, True])
def test_output_inside_payload_does_not_mutate_input(tmp_path, bypass):
    root = tmp_path / "payload"
    fixture(root)
    before = module.inventory(root)
    output = root / "attempt"
    if bypass:
        output = root / "comparison" / ".." / "attempt"
    with pytest.raises(ValueError, match="inside payload"):
        module.audit(root, SHA, "/python", output)
    assert module.inventory(root) == before
    assert not (root / "attempt").exists()


def test_cli_rejects_symlink_traversal_without_mutation(tmp_path):
    root = tmp_path / "payload"
    fixture(root)
    before = module.inventory(root)
    alias = tmp_path / "alias"
    alias.symlink_to(root, target_is_directory=True)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.release_molecular_audits",
            "--payload",
            str(root),
            "--source-revision",
            SHA,
            "--molecular-python",
            "/python",
            "--output",
            str(alias / "attempt"),
        ],
        cwd=module.ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 1 and "symlink" in result.stderr
    assert "Traceback" not in result.stderr
    assert module.inventory(root) == before


def test_timeout_preserves_output_and_failed_receipt(tmp_path):
    fixture(tmp_path / "payload")

    def timeout(command, **kwargs):
        assert kwargs["timeout"] == 300
        raise subprocess.TimeoutExpired(
            command, 300, output=b"partial out", stderr=b"partial err"
        )

    result = execute(tmp_path, timeout)
    assert result["status"] == "failed" and "TimeoutExpired" in result["error"]
    log = (tmp_path / "out" / "auditor.log").read_text()
    assert "partial out" in log and "partial err" in log
    assert (
        json.loads((tmp_path / "out" / "receipt.json").read_text())["status"]
        == "failed"
    )
