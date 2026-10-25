"""Release gates reject archive leakage, omissions and false installed checks."""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

from scripts.check_distribution import (
    DistributionError,
    check_archives,
    check_installed,
)

DIST_INFO = "kernelgecp-1.0.0.dist-info"
METADATA = b"""\
Metadata-Version: 2.4
Name: kernelgecp
Version: 1.0.0
Requires-Python: >=3.12
License-Expression: Apache-2.0

Package description.
"""


def _fixture_archives(
    tmp_path: Path,
    *,
    wheel_changes: dict[str, bytes | None] | None = None,
    sdist_changes: dict[str, bytes | None] | None = None,
) -> tuple[Path, Path, Path]:
    root = tmp_path / "checkout"
    package = root / "kernelgecp"
    package.mkdir(parents=True)
    (root / "pyproject.toml").write_text(
        '[project]\nname = "kernelgecp"\nversion = "1.0.0"\n'
        'requires-python = ">=3.12"\n',
        encoding="utf-8",
    )
    (root / "LICENSE").write_text("Apache License\nVersion 2.0\n", encoding="utf-8")
    (root / "APPLICATION.md").write_text("Application guide.\n", encoding="utf-8")
    (package / "__init__.py").write_text('__version__ = "1.0.0"\n', encoding="utf-8")
    (package / "basis.py").write_text('"""Public basis."""\n', encoding="utf-8")
    (package / "py.typed").touch()
    sources = {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in package.iterdir()
    }
    wheel_files = {
        **sources,
        f"{DIST_INFO}/METADATA": METADATA,
        f"{DIST_INFO}/WHEEL": (
            b"Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n"
        ),
        f"{DIST_INFO}/RECORD": b"",
        f"{DIST_INFO}/licenses/LICENSE": (root / "LICENSE").read_bytes(),
    }
    sdist_files = {
        **sources,
        "pyproject.toml": (root / "pyproject.toml").read_bytes(),
        "APPLICATION.md": (root / "APPLICATION.md").read_bytes(),
        "LICENSE": (root / "LICENSE").read_bytes(),
        "PKG-INFO": METADATA,
    }
    for files, changes in (
        (wheel_files, wheel_changes or {}),
        (sdist_files, sdist_changes or {}),
    ):
        for name, payload in changes.items():
            if payload is None:
                files.pop(name, None)
            else:
                files[name] = payload
    wheel = tmp_path / "kernelgecp-1.0.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, payload in wheel_files.items():
            archive.writestr(name, payload)
    sdist = tmp_path / "kernelgecp-1.0.0.tar.gz"
    with tarfile.open(sdist, "w:gz") as archive:
        for name, payload in sdist_files.items():
            member = tarfile.TarInfo(f"kernelgecp-1.0.0/{name}")
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))
    return root, wheel, sdist


def test_distribution_gate_accepts_complete_minimal_archives(tmp_path: Path) -> None:
    root, wheel, sdist = _fixture_archives(tmp_path)
    report = check_archives(project_root=root, wheel=wheel, sdist=sdist)
    assert report["name"] == "kernelgecp"
    assert report["version"] == "1.0.0"
    assert report["package_source_files"] == 3
    assert all(len(archive["sha256"]) == 64 for archive in report["archives"])
    json.dumps(report, allow_nan=False)


def _native_fixture(
    tmp_path: Path, *, changes: dict[str, bytes | None] | None = None
) -> tuple[Path, Path, Path]:
    notices = b"Complete third-party notices and LGPL/GPL source instructions.\n"
    sbom = b'{"bomFormat":"CycloneDX","specVersion":"1.5","components":[]}\n'
    expression = "Apache-2.0 AND LGPL-3.0-or-later"
    inventory = {
        "schema": "kernelgecp.native-wheel.v1",
        "version": "1.0.0",
        "source_revision": "a" * 40,
        "source_tree": "b" * 40,
        "source_asset": {
            "filename": "kernelgecp-1.0.0-native-source.tar.gz",
            "sha256": "c" * 64,
            "url": "https://github.com/tripp-smith/gecp-kernel-structure/releases/download/v1.0.0/kernelgecp-1.0.0-native-source.tar.gz",
        },
        "license_bundle_sha256": hashlib.sha256(notices).hexdigest(),
        "sbom_sha256": hashlib.sha256(sbom).hexdigest(),
        "license_expression": expression,
    }
    stub = b"def build_info() -> dict[str, str | int | bool]: ...\n"
    payload = {
        "kernelgecp/_native.abi3.so": b"\x7fELF" + bytes(60),
        "kernelgecp/_native.pyi": stub,
        "kernelgecp/_native_inventory.json": json.dumps(inventory).encode(),
        f"{DIST_INFO}/METADATA": METADATA.replace(
            b"Apache-2.0", expression.encode()
        ).replace(
            b"\n\nPackage description.",
            b"\nLicense-File: LICENSE\nLicense-File: native/_wheel_licenses.txt"
            b"\n\nPackage description.",
        ),
        f"{DIST_INFO}/WHEEL": (
            b"Wheel-Version: 1.0\nRoot-Is-Purelib: false\n"
            b"Tag: cp312-abi3-manylinux_2_28_x86_64\n"
        ),
        f"{DIST_INFO}/licenses/native/_wheel_licenses.txt": notices,
        f"{DIST_INFO}/sboms/_wheel_dependencies.cyclonedx.json": sbom,
    }
    payload.update(changes or {})
    root, wheel, sdist = _fixture_archives(
        tmp_path,
        wheel_changes=payload,
        sdist_changes={"kernelgecp/_native.pyi": stub},
    )
    (root / "kernelgecp/_native.pyi").write_bytes(stub)
    native_wheel = wheel.with_name(
        "kernelgecp-1.0.0-cp312-abi3-manylinux_2_28_x86_64.whl"
    )
    wheel.rename(native_wheel)
    return root, native_wheel, sdist


def test_native_gate_accepts_explicit_platform_wheel_and_identical_sources(
    tmp_path: Path,
) -> None:
    root, wheel, sdist = _native_fixture(tmp_path)
    result = check_archives(project_root=root, wheel=wheel, sdist=sdist, kind="native")
    assert result["kind"] == "native" and result["package_source_files"] == 4
    with pytest.raises(DistributionError, match="universal"):
        check_archives(project_root=root, wheel=wheel, sdist=sdist)


@pytest.mark.parametrize(
    "changes,message",
    [
        ({"kernelgecp/_native_inventory.json": None}, "missing inventory"),
        ({"kernelgecp/_native.abi3.so": b"not a binary"}, "malformed"),
        (
            {
                "kernelgecp/_native.abi3.so": b"\x7fELF"
                + bytes(30)
                + b"/home/runner/build"
            },
            "host build path",
        ),
        ({"kernelgecp/arbitrary.so": b"\x7fELF"}, "unexpected wheel"),
        ({f"{DIST_INFO}/licenses/native/_wheel_licenses.txt": None}, "missing notices"),
        (
            {f"{DIST_INFO}/licenses/native/_wheel_licenses.txt": b"changed"},
            "hash mismatch",
        ),
        (
            {f"{DIST_INFO}/sboms/_wheel_dependencies.cyclonedx.json": b"{}"},
            "hash mismatch",
        ),
        ({f"{DIST_INFO}/METADATA": METADATA}, "license metadata"),
        (
            {
                f"{DIST_INFO}/WHEEL": (
                    b"Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n"
                )
            },
            "platform tags",
        ),
        ({"kernelgecp/_native.pyi": None}, "sources differ"),
    ],
)
def test_native_gate_rejects_missing_or_forged_native_payload(
    tmp_path: Path, changes: dict[str, bytes | None], message: str
) -> None:
    root, wheel, sdist = _native_fixture(tmp_path, changes=changes)
    with pytest.raises(DistributionError, match=message):
        check_archives(project_root=root, wheel=wheel, sdist=sdist, kind="native")


@pytest.mark.parametrize(
    "field,value",
    [
        ("version", "9.0.0"),
        ("source_revision", "HEAD"),
        ("license_expression", "Apache-2.0"),
    ],
)
def test_native_inventory_rejects_wrong_identity(
    tmp_path: Path, field: str, value: str
) -> None:
    root, wheel, sdist = _native_fixture(tmp_path)
    with zipfile.ZipFile(wheel) as archive:
        files = {name: archive.read(name) for name in archive.namelist()}
    data = json.loads(files["kernelgecp/_native_inventory.json"])
    data[field] = value
    files["kernelgecp/_native_inventory.json"] = json.dumps(data).encode()
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, payload in files.items():
            archive.writestr(name, payload)
    with pytest.raises(DistributionError):
        check_archives(project_root=root, wheel=wheel, sdist=sdist, kind="native")


def test_native_sdist_must_include_current_rust_sources(tmp_path: Path) -> None:
    root, wheel, sdist = _native_fixture(tmp_path)
    (root / "native").mkdir()
    (root / "native/Cargo.toml").write_text('[package]\nname="kernelgecp-native"\n')
    with pytest.raises(DistributionError, match="native sources differ"):
        check_archives(project_root=root, wheel=wheel, sdist=sdist, kind="native")


@pytest.mark.parametrize(
    ("name", "payload", "message"),
    [
        ("FINDINGS 2.md", b"do not read", "excluded"),
        ("experiments/data/run.json", b"{}", "excluded"),
        ("kernelgecp/data.json", b"{}", "unexpected wheel payload"),
        ("../outside.py", b"", "unsafe"),
        (f"{DIST_INFO}/private-note.txt", b"", "unexpected wheel payload"),
        ("kernelgecp/py.typed", None, "sources differ"),
        ("kernelgecp/basis.py", None, "sources differ"),
        (f"{DIST_INFO}/licenses/LICENSE", b"other license", "license differs"),
        (f"{DIST_INFO}/METADATA", None, "exactly one METADATA"),
        (f"{DIST_INFO}/WHEEL", b"Tag: py3-none-any\n", "pure universal"),
    ],
)
def test_distribution_gate_rejects_bad_wheel_payload(
    tmp_path: Path, name: str, payload: bytes | None, message: str
) -> None:
    root, wheel, sdist = _fixture_archives(tmp_path, wheel_changes={name: payload})
    with pytest.raises(DistributionError, match=message):
        check_archives(project_root=root, wheel=wheel, sdist=sdist)


@pytest.mark.parametrize(
    ("name", "payload", "message"),
    [
        ("MATHLIB 2.md", b"do not read", "excluded"),
        ("RESEARCH.md", b"research document", "unexpected sdist payload"),
        ("experiments/data/result.json", b"{}", "excluded"),
        ("kernelgecp/py.typed", None, "sources differ"),
        ("APPLICATION.md", None, "sdist missing APPLICATION"),
        ("LICENSE", b"wrong", "LICENSE differs"),
        ("PKG-INFO", None, "sdist missing PKG-INFO"),
    ],
)
def test_distribution_gate_rejects_bad_sdist_payload(
    tmp_path: Path, name: str, payload: bytes | None, message: str
) -> None:
    root, wheel, sdist = _fixture_archives(tmp_path, sdist_changes={name: payload})
    with pytest.raises(DistributionError, match=message):
        check_archives(project_root=root, wheel=wheel, sdist=sdist)


@pytest.mark.parametrize(
    "metadata",
    [
        b"",
        METADATA.replace(b"Version: 1.0.0\n", b"Version: 2.0.0\n"),
        METADATA.replace(b"Name: kernelgecp\n", b"Name: another-project\n"),
        METADATA.replace(b"Requires-Python: >=3.12\n", b"Requires-Python: >=3.13\n"),
        METADATA.replace(b"License-Expression: Apache-2.0\n", b""),
        METADATA.replace(b"Name: kernelgecp\n", b"Name: kernelgecp\nName: other\n"),
    ],
)
@pytest.mark.parametrize("archive_kind", ["wheel", "sdist"])
def test_distribution_gate_rejects_malformed_metadata(
    tmp_path: Path, metadata: bytes, archive_kind: str
) -> None:
    kwargs = (
        {"wheel_changes": {f"{DIST_INFO}/METADATA": metadata}}
        if archive_kind == "wheel"
        else {"sdist_changes": {"PKG-INFO": metadata}}
    )
    root, wheel, sdist = _fixture_archives(tmp_path, **kwargs)
    with pytest.raises(DistributionError):
        check_archives(project_root=root, wheel=wheel, sdist=sdist)


def test_distribution_gate_cli_returns_nonzero_on_failure(tmp_path: Path) -> None:
    root, wheel, sdist = _fixture_archives(
        tmp_path, wheel_changes={"kernelgecp/py.typed": None}
    )
    script = Path(__file__).parents[1] / "scripts/check_distribution.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "inspect",
            "--project-root",
            str(root),
            "--wheel",
            str(wheel),
            "--sdist",
            str(sdist),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "distribution check failed" in result.stderr
    assert "sha256" not in result.stdout


def test_installed_gate_rejects_running_inside_checkout(tmp_path: Path) -> None:
    with pytest.raises(DistributionError, match="outside the checkout"):
        check_installed(project_root=Path.cwd())


def test_installed_gate_rejects_checkout_import(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = Path(__file__).parents[1]
    monkeypatch.chdir(tmp_path)
    with pytest.raises(DistributionError, match="checkout source"):
        check_installed(project_root=root)
