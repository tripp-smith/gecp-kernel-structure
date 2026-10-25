"""Native wheel staging rejects forged asset identities and metadata."""

import hashlib
import json
import tomllib
from pathlib import Path

import pytest

from scripts.build_native_wheel import (
    NativeWheelError,
    overlay_project_metadata,
    read_source_report,
)


def _report(tmp_path: Path) -> Path:
    report = {
        "schema": "kernelgecp.corresponding-source.receipt.v1",
        "version": "1.0.0",
        "source_revision": "a" * 40,
        "source_tree": "b" * 40,
        "archive": "kernelgecp-1.0.0-native-source.tar.gz",
        "license_bundle": "kernelgecp-1.0.0-native-licenses.txt",
        "sbom": "kernelgecp-1.0.0-native-dependencies.cyclonedx.json",
        "combined_native_license_expression": "Apache-2.0 AND LGPL-3.0-or-later",
    }
    for name, hash_name in (
        ("archive", "sha256"),
        ("license_bundle", "license_bundle_sha256"),
        ("sbom", "sbom_sha256"),
    ):
        payload = (name + " deterministic test data").encode()
        (tmp_path / report[name]).write_bytes(payload)
        report[hash_name] = hashlib.sha256(payload).hexdigest()
    path = tmp_path / "source-report.json"
    path.write_text(json.dumps(report))
    return path


def test_report_checks_all_external_artifact_hashes(tmp_path: Path) -> None:
    path = _report(tmp_path)
    value = read_source_report(path)
    assert value["source_revision"] == "a" * 40
    (tmp_path / value["license_bundle"]).write_bytes(b"changed notices")
    with pytest.raises(NativeWheelError, match="hash mismatch"):
        read_source_report(path)


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("schema", "unknown", "schema"),
        ("source_revision", "HEAD", "source_revision"),
        ("source_tree", "a" * 39, "source_tree"),
        ("archive", "../outside.tar.gz", "basename"),
        ("license_bundle", "/tmp/licenses.txt", "basename"),
        ("sbom", "folder\\inventory.json", "basename"),
        ("sha256", "not a digest", "sha256"),
        ("combined_native_license_expression", "Apache-2.0", "license expression"),
        (
            "combined_native_license_expression",
            "Apache-2.0\n[project]\n",
            "license expression",
        ),
    ],
)
def test_report_rejects_forgery(
    tmp_path: Path, field: str, value: str, message: str
) -> None:
    path = _report(tmp_path)
    data = json.loads(path.read_text())
    data[field] = value
    path.write_text(json.dumps(data))
    with pytest.raises(NativeWheelError, match=message):
        read_source_report(path)


def test_report_rejects_missing_or_linked_external_assets(tmp_path: Path) -> None:
    path = _report(tmp_path)
    data = json.loads(path.read_text())
    archive = tmp_path / data["archive"]
    payload = archive.read_bytes()
    archive.unlink()
    with pytest.raises(NativeWheelError, match="missing archive"):
        read_source_report(path)
    target = tmp_path / "other.tar.gz"
    target.write_bytes(payload)
    try:
        archive.symlink_to(target)
    except OSError:
        pytest.skip("platform does not permit local symbolic links")
    with pytest.raises(NativeWheelError, match="missing archive"):
        read_source_report(path)


def test_metadata_overlay_changes_only_staged_native_license_and_inventory() -> None:
    original = (
        '[project]\nname = "kernelgecp"\nversion = "1.0.0"\n'
        'license = { text = "Apache-2.0" }\n'
        "[tool.maturin.sbom]\nrust = false\nauditwheel = false\n"
    )
    expression = "Apache-2.0 AND LGPL-3.0-or-later AND MIT"
    data = tomllib.loads(overlay_project_metadata(original, expression))
    assert data["project"]["license"] == expression
    assert data["project"]["license-files"] == ["LICENSE", "native/_wheel_licenses.txt"]
    assert data["tool"]["maturin"]["sbom"] == {
        "rust": False,
        "auditwheel": False,
        "include": ["native/_wheel_dependencies.cyclonedx.json"],
    }
    assert tomllib.loads(original)["project"]["license"] == {"text": "Apache-2.0"}


def test_metadata_overlay_rejects_automatic_host_sbom() -> None:
    text = (
        '[project]\nlicense = { text = "Apache-2.0" }\n'
        "[tool.maturin.sbom]\nrust = true\nauditwheel = false\n"
    )
    with pytest.raises(NativeWheelError):
        overlay_project_metadata(text, "Apache-2.0 AND LGPL-3.0-or-later")
