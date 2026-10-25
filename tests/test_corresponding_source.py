"""Corresponding-source archives bind Git sources and complete vendor payloads."""

from __future__ import annotations

import copy
import gzip
import hashlib
import io
import json
import subprocess
import tarfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from scripts.prepare_corresponding_source import (
    SourceAssetError,
    inspect_corresponding_source,
    prepare_corresponding_source,
)

_CRATES = {"pyo3": "0.27.2", "rug": "1.30.0", "gmp-mpfr-sys": "1.7.1"}
_REQUIRED = {
    "pyproject.toml": '[project]\nname = "kernelgecp"\nversion = "1.0.0"\n',
    "uv.lock": "version = 1\n",
    "LICENSE": "Apache License, Version 2.0\n",
    "README.md": "Fixture project\n",
    "APPLICATION.md": "Fixture application guide\n",
    "rust-toolchain.toml": '[toolchain]\nchannel = "1.92.0"\n',
    "scripts/prepare_corresponding_source.py": '"""Fixture source tool."""\n',
    "scripts/offline_native_build.sh": "#!/bin/sh\nset -eu\n",
    "native/Cargo.toml": (
        '[package]\nname = "kernelgecp-native"\nversion = "1.0.0"\n'
        'edition = "2021"\nlicense = "Apache-2.0"\n'
    ),
    "native/build.rs": "fn main() {}\n",
    "native/src/lib.rs": "pub fn fixture() {}\n",
    "native/README.md": "Native fixture\n",
    "native/SOURCE_BUILD.md": "Build instructions\n",
    "native/THIRD_PARTY_NOTICES.md": "Vendored dependencies and notices\n",
    "native/licenses/LGPL-3.0.txt": "Fixture LGPL notice\n",
    "native/licenses/GPL-3.0.txt": "Fixture GPL notice\n",
    "kernelgecp/__init__.py": '__version__ = "1.0.0"\n',
    "kernelgecp/_native.pyi": "NATIVE_PROTOCOL_VERSION: int\n",
    "kernelgecp/py.typed": "",
}


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _git(root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *arguments], text=True, stderr=subprocess.STDOUT
    ).strip()


def _commit(root: Path) -> str:
    _git(root, "add", "--all")
    _git(root, "update-index", "--chmod=+x", "scripts/offline_native_build.sh")
    _git(
        root,
        "-c",
        "user.name=Source fixture",
        "-c",
        "user.email=source-fixture@example.invalid",
        "-c",
        "commit.gpgsign=false",
        "-c",
        f"core.hooksPath={root.parent / 'nonexistent-hooks'}",
        "commit",
        "-qm",
        "Freeze fixture",
    )
    return _git(root, "rev-parse", "HEAD")


def _refresh_vendor(crate: Path, package_checksum: str | None = None) -> None:
    checksum_path = crate / ".cargo-checksum.json"
    if package_checksum is None:
        package_checksum = json.loads(checksum_path.read_text())["package"]
    files = {
        path.relative_to(crate).as_posix(): _sha(path.read_bytes())
        for path in sorted(crate.rglob("*"))
        if path.is_file() and path != checksum_path
    }
    checksum_path.write_text(
        json.dumps({"package": package_checksum, "files": files}, sort_keys=True)
    )


def _vendor_crate(vendor: Path, name: str, version: str) -> Path:
    crate = vendor / f"{name}-{version}"
    (crate / "src").mkdir(parents=True)
    license_expression = "MIT" if name == "pyo3" else "LGPL-3.0+"
    (crate / "Cargo.toml").write_text(
        f'[package]\nname = "{name}"\nversion = "{version}"\n'
        f'license = "{license_expression}"\n'
    )
    (crate / "LICENSE").write_text(f"Nonempty fixture license for {name}\n")
    (crate / "src/lib.rs").write_text("pub fn fixture_dependency() {}\n")
    _refresh_vendor(crate, _sha(f"registry package {name}-{version}".encode()))
    return crate


@dataclass(frozen=True)
class _Fixture:
    root: Path
    vendor: Path
    revision: str
    output: Path

    def prepare(self, *, output: Path | None = None) -> dict[str, Any]:
        return prepare_corresponding_source(
            project_root=self.root,
            revision=self.revision,
            vendor=self.vendor,
            output=output or self.output,
        )

    def archive(self, report: dict[str, Any], *, output: Path | None = None) -> Path:
        name = Path(report["archive"])
        assert str(name) == name.name and not name.is_absolute()
        return (output or self.output) / name


@pytest.fixture
def source_fixture(tmp_path: Path) -> _Fixture:
    root, vendor = tmp_path / "project", tmp_path / "vendor"
    root.mkdir()
    vendor.mkdir()
    _git(root, "init", "--object-format=sha1", "-q", "-b", "main")
    _git(root, "config", "core.autocrlf", "false")
    for relative, content in _REQUIRED.items():
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(content)
    lock = [
        "version = 4\n",
        '[[package]]\nname = "kernelgecp-native"\nversion = "1.0.0"\n',
    ]
    for name, version in _CRATES.items():
        crate = _vendor_crate(vendor, name, version)
        checksum = json.loads((crate / ".cargo-checksum.json").read_text())["package"]
        lock.append(
            f'[[package]]\nname = "{name}"\nversion = "{version}"\n'
            'source = "registry+https://github.com/rust-lang/crates.io-index"\n'
            f'checksum = "{checksum}"\n'
        )
    (root / "native/Cargo.lock").write_text("\n".join(lock))
    return _Fixture(root, vendor, _commit(root), tmp_path / "output")


def _members(archive: Path) -> list[tuple[tarfile.TarInfo, bytes | None]]:
    records = []
    with tarfile.open(archive, "r:gz") as stream:
        for member in stream.getmembers():
            source = stream.extractfile(member) if member.isfile() else None
            records.append((copy.copy(member), source.read() if source else None))
    return records


def _rewrite(
    archive: Path,
    directory: Path,
    mutation: Callable[
        [list[tuple[tarfile.TarInfo, bytes | None]]],
        list[tuple[tarfile.TarInfo, bytes | None]],
    ],
) -> Path:
    records = mutation(_members(archive))
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w", format=tarfile.PAX_FORMAT) as stream:
        for member, payload in records:
            if payload is not None:
                member.size = len(payload)
            stream.addfile(member, io.BytesIO(payload) if payload is not None else None)
    directory.mkdir()
    target = directory / archive.name
    target.write_bytes(gzip.compress(raw.getvalue(), mtime=0))
    return target


def _rewrite_manifest(
    archive: Path, directory: Path, mutation: Callable[[dict[str, Any]], None]
) -> Path:
    def mutate(
        records: list[tuple[tarfile.TarInfo, bytes | None]],
    ) -> list[tuple[tarfile.TarInfo, bytes | None]]:
        result = []
        found = 0
        for member, payload in records:
            if member.name.endswith("/SOURCE_MANIFEST.json"):
                assert payload is not None
                manifest = json.loads(payload)
                mutation(manifest)
                payload = json.dumps(manifest, sort_keys=True).encode()
                found += 1
            result.append((member, payload))
        assert found == 1
        return result

    return _rewrite(archive, directory, mutate)


def test_source_asset_success_is_deterministic_and_git_bound(
    source_fixture: _Fixture, tmp_path: Path
) -> None:
    first = source_fixture.prepare()
    second = source_fixture.prepare(output=tmp_path / "repeat")
    archive = source_fixture.archive(first)
    assert archive.name == "kernelgecp-1.0.0-native-source.tar.gz"
    assert (
        archive.read_bytes()
        == source_fixture.archive(second, output=tmp_path / "repeat").read_bytes()
    )
    assert first == second
    assert first["schema"] == "kernelgecp.corresponding-source.receipt.v1"
    assert archive.read_bytes()[4:8] == bytes(4)
    assert first["sha256"] == _sha(archive.read_bytes())
    assert first["source_revision"] == source_fixture.revision
    assert first["source_tree"] == _git(source_fixture.root, "rev-parse", "HEAD^{tree}")
    assert first["version"] == "1.0.0"
    assert first["vendored_crates"] == 3
    assert first["archive_bytes"] == archive.stat().st_size
    assert isinstance(first["combined_native_license_expression"], str)
    assert "Apache-2.0" in first["combined_native_license_expression"]
    assert "LGPL-3.0-or-later" in first["combined_native_license_expression"]
    assert "LGPL-3.0+" not in first["combined_native_license_expression"]
    for name_key, hash_key in (
        ("license_bundle", "license_bundle_sha256"),
        ("sbom", "sbom_sha256"),
    ):
        name = Path(first[name_key])
        assert str(name) == name.name and not name.is_absolute()
        assert _sha((source_fixture.output / name).read_bytes()) == first[hash_key]
    inspected = inspect_corresponding_source(
        archive,
        expected_sha256=first["sha256"],
        project_root=source_fixture.root,
        revision=source_fixture.revision,
    )
    assert inspected["source_revision"] == source_fixture.revision
    assert inspected["payload_files"] == first["payload_files"]
    portable = inspect_corresponding_source(archive, expected_sha256=first["sha256"])
    assert portable["source_revision"] == source_fixture.revision
    records = _members(archive)
    manifests = [
        payload
        for member, payload in records
        if member.name.endswith("/SOURCE_MANIFEST.json")
    ]
    assert len(manifests) == 1 and manifests[0] is not None
    manifest = json.loads(manifests[0])
    assert manifest["schema"] == "kernelgecp.corresponding-source.v1"
    assert manifest["lock_sha256"] == _sha(
        (source_fixture.root / "native/Cargo.lock").read_bytes()
    )
    paths = {item["path"] for item in manifest["files"]}
    assert "SOURCE_MANIFEST.json" not in paths
    assert {".cargo/config.toml", "NATIVE_DEPENDENCIES.cyclonedx.json"} <= paths
    assert {item["name"] for item in manifest["packages"]} == set(_CRATES)
    assert all(
        item["license"] and item["license_files"] for item in manifest["packages"]
    )
    assert str(tmp_path) not in json.dumps(manifest)
    script = next(
        member
        for member, _ in records
        if member.name.endswith("/scripts/offline_native_build.sh")
    )
    assert script.mode & 0o777 == 0o755


def test_source_asset_refuses_overwrite(source_fixture: _Fixture) -> None:
    first = source_fixture.prepare()
    original = source_fixture.archive(first).read_bytes()
    with pytest.raises(SourceAssetError):
        source_fixture.prepare()
    assert source_fixture.archive(first).read_bytes() == original


def test_source_asset_rejects_a_file_as_output_directory(
    source_fixture: _Fixture,
) -> None:
    source_fixture.output.write_text("existing caller file")
    with pytest.raises(SourceAssetError):
        source_fixture.prepare()
    assert source_fixture.output.read_text() == "existing caller file"


@pytest.mark.parametrize("revision", ["HEAD", "a" * 39, "not-a-revision", "0" * 40])
def test_prepare_rejects_nonfull_or_unknown_revision(
    source_fixture: _Fixture, revision: str
) -> None:
    with pytest.raises(SourceAssetError):
        prepare_corresponding_source(
            project_root=source_fixture.root,
            revision=revision,
            vendor=source_fixture.vendor,
            output=source_fixture.output,
        )


def test_prepare_rejects_dirty_selected_source(source_fixture: _Fixture) -> None:
    (source_fixture.root / "README.md").write_text("Uncommitted source\n")
    with pytest.raises(SourceAssetError):
        source_fixture.prepare()


def test_inspection_uses_git_blobs_with_a_dirty_checkout(
    source_fixture: _Fixture,
) -> None:
    report = source_fixture.prepare()
    (source_fixture.root / "README.md").write_text("Later local changes\n")
    result = inspect_corresponding_source(
        source_fixture.archive(report),
        expected_sha256=report["sha256"],
        project_root=source_fixture.root,
        revision=source_fixture.revision,
    )
    assert result["source_revision"] == source_fixture.revision


def test_prepare_rejects_project_native_version_mismatch(
    source_fixture: _Fixture,
) -> None:
    path = source_fixture.root / "native/Cargo.toml"
    path.write_text(path.read_text().replace('version = "1.0.0"', 'version = "9.9.9"'))
    revision = _commit(source_fixture.root)
    with pytest.raises(SourceAssetError):
        prepare_corresponding_source(
            project_root=source_fixture.root,
            revision=revision,
            vendor=source_fixture.vendor,
            output=source_fixture.output,
        )


def test_prepare_rejects_missing_required_committed_source(
    source_fixture: _Fixture,
) -> None:
    (source_fixture.root / "native/SOURCE_BUILD.md").unlink()
    revision = _commit(source_fixture.root)
    with pytest.raises(SourceAssetError):
        prepare_corresponding_source(
            project_root=source_fixture.root,
            revision=revision,
            vendor=source_fixture.vendor,
            output=source_fixture.output,
        )


@pytest.mark.parametrize(
    "failure",
    [
        "file_checksum",
        "package_checksum",
        "extra_file",
        "missing_file",
        "missing_checksum",
        "malformed_checksum",
        "checksum_json_list",
        "unsafe_checksum_path",
        "missing_license",
        "empty_license_declaration",
        "extra_crate",
    ],
)
def test_prepare_rejects_incomplete_or_unbound_vendor(
    source_fixture: _Fixture, failure: str
) -> None:
    crate = source_fixture.vendor / "rug-1.30.0"
    if failure == "file_checksum":
        (crate / "src/lib.rs").write_text("tampered\n")
    elif failure == "package_checksum":
        _refresh_vendor(crate, "0" * 64)
    elif failure == "extra_file":
        (crate / "unlisted.txt").write_text("unbound\n")
    elif failure == "missing_file":
        (crate / "src/lib.rs").unlink()
    elif failure == "missing_checksum":
        (crate / ".cargo-checksum.json").unlink()
    elif failure == "malformed_checksum":
        (crate / ".cargo-checksum.json").write_text("not-json")
    elif failure == "checksum_json_list":
        (crate / ".cargo-checksum.json").write_text("[]")
    elif failure == "unsafe_checksum_path":
        path = crate / ".cargo-checksum.json"
        checksums = json.loads(path.read_text())
        checksums["files"]["../outside"] = "0" * 64
        path.write_text(json.dumps(checksums))
    elif failure == "missing_license":
        (crate / "LICENSE").unlink()
        _refresh_vendor(crate)
    elif failure == "empty_license_declaration":
        path = crate / "Cargo.toml"
        path.write_text(
            path.read_text().replace('license = "LGPL-3.0+"', 'license = ""')
        )
        _refresh_vendor(crate)
    else:
        _vendor_crate(source_fixture.vendor, "unlocked", "0.1.0")
    with pytest.raises(SourceAssetError):
        source_fixture.prepare()


def test_inspection_rejects_wrong_outer_sha(source_fixture: _Fixture) -> None:
    report = source_fixture.prepare()
    with pytest.raises(SourceAssetError):
        inspect_corresponding_source(
            source_fixture.archive(report), expected_sha256="0" * 64
        )


def test_prepare_rejects_checksum_bound_vendor_symlink(
    source_fixture: _Fixture,
) -> None:
    crate = source_fixture.vendor / "rug-1.30.0"
    try:
        (crate / "src/linked.rs").symlink_to(crate / "src/lib.rs")
    except OSError as error:
        pytest.skip(f"test environment cannot create a symlink: {error}")
    _refresh_vendor(crate)
    with pytest.raises(SourceAssetError):
        source_fixture.prepare()


@pytest.mark.parametrize(
    "failure",
    [
        "schema",
        "extra_manifest_field",
        "inventory_sha",
        "inventory_missing",
        "duplicate_inventory",
        "inventory_size",
        "inventory_mode",
        "extra_file_field",
        "source_blob",
        "source_blob_sha",
        "extra_blob_field",
        "source_revision",
        "package_checksum",
        "license_inventory",
    ],
)
def test_inspection_rejects_manifest_tampering_after_rehash(
    source_fixture: _Fixture, tmp_path: Path, failure: str
) -> None:
    report = source_fixture.prepare()

    def mutate(manifest: dict[str, Any]) -> None:
        if failure == "schema":
            manifest["schema"] = "unrecognized"
        elif failure == "extra_manifest_field":
            manifest["unrecognized"] = "unbound metadata"
        elif failure == "inventory_sha":
            manifest["files"][0]["sha256"] = "0" * 64
        elif failure == "inventory_missing":
            manifest["files"].pop()
        elif failure == "duplicate_inventory":
            manifest["files"].append(manifest["files"][0].copy())
        elif failure == "inventory_size":
            manifest["files"][0]["size"] = -1
        elif failure == "inventory_mode":
            manifest["files"][0]["mode"] = "invalid"
        elif failure == "extra_file_field":
            manifest["files"][0]["unrecognized"] = "unbound metadata"
        elif failure == "source_blob":
            manifest["source_blobs"][0]["git_blob"] = "0" * 40
        elif failure == "source_blob_sha":
            manifest["source_blobs"][0]["sha256"] = "0" * 64
        elif failure == "extra_blob_field":
            manifest["source_blobs"][0]["unrecognized"] = "unbound metadata"
        elif failure == "source_revision":
            manifest["source_revision"] = "bad-revision"
        elif failure == "package_checksum":
            manifest["packages"][0]["registry_checksum"] = "0" * 64
        else:
            manifest["packages"][0]["license_files"] = []

    tampered = _rewrite_manifest(
        source_fixture.archive(report), tmp_path / "tampered", mutate
    )
    with pytest.raises(SourceAssetError):
        inspect_corresponding_source(
            tampered,
            expected_sha256=_sha(tampered.read_bytes()),
            project_root=source_fixture.root,
            revision=source_fixture.revision,
        )


@pytest.mark.parametrize(
    "failure",
    [
        "unsafe_path",
        "symlink",
        "duplicate_member",
        "missing_manifest",
        "payload_tamper",
    ],
)
def test_inspection_rejects_unsafe_or_inconsistent_tar_after_rehash(
    source_fixture: _Fixture, tmp_path: Path, failure: str
) -> None:
    report = source_fixture.prepare()

    def mutate(
        records: list[tuple[tarfile.TarInfo, bytes | None]],
    ) -> list[tuple[tarfile.TarInfo, bytes | None]]:
        first = copy.copy(records[0][0])
        root = first.name.split("/", 1)[0]
        if failure == "unsafe_path":
            first.name = root + "/../outside.txt"
            return records + [(first, b"unsafe")]
        if failure == "symlink":
            first.name = root + "/linked"
            first.type = tarfile.SYMTYPE
            first.linkname = "/outside"
            first.size = 0
            return records + [(first, None)]
        if failure == "duplicate_member":
            return records + [(copy.copy(records[0][0]), records[0][1])]
        if failure == "missing_manifest":
            return [
                (member, payload)
                for member, payload in records
                if not member.name.endswith("/SOURCE_MANIFEST.json")
            ]
        changed = False
        result = []
        for member, payload in records:
            if member.name.endswith("/README.md") and not changed:
                payload = b"changed payload"
                changed = True
            result.append((member, payload))
        assert changed
        return result

    tampered = _rewrite(source_fixture.archive(report), tmp_path / "tampered", mutate)
    with pytest.raises(SourceAssetError):
        inspect_corresponding_source(
            tampered, expected_sha256=_sha(tampered.read_bytes())
        )


def test_inspection_rejects_malformed_compressed_input(tmp_path: Path) -> None:
    archive = tmp_path / "kernelgecp-1.0.0-native-source.tar.gz"
    archive.write_bytes(b"not gzip or tar")
    with pytest.raises(SourceAssetError):
        inspect_corresponding_source(
            archive, expected_sha256=_sha(archive.read_bytes())
        )
