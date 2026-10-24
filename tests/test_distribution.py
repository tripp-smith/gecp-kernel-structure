"""Release gates reject archive leakage, omissions and false installed checks."""

from __future__ import annotations

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
