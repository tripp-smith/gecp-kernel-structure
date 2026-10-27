"""Replayable installed conda locks reject private or unidentified artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from scripts.providers import export_conda_lock as exporter


def _record(name: str = "example", *, subdir: str = "linux-64") -> dict[str, Any]:
    return {
        "name": name,
        "version": "1.2.3",
        "build": "py312_example_0",
        "subdir": subdir,
        "sha256": "a" * 64,
        "url": (
            f"https://conda.anaconda.org/conda-forge/{subdir}/"
            f"{name}-1.2.3-py312_example_0.conda"
        ),
        "depends": ["python >=3.12,<3.13", "numpy >=2,<3"],
        "license": "BSD-3-Clause",
        "link": {"source": "/private/build/cache", "prefix": "/private/environment"},
    }


def _write(
    prefix: Path, record: dict[str, Any], filename: str = "example.json"
) -> None:
    directory = prefix / "conda-meta"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / filename).write_text(json.dumps(record))


def test_lock_is_exact_hash_pinned_public_and_deterministic(tmp_path: Path) -> None:
    _write(tmp_path, _record("z_native"), "z.json")
    _write(tmp_path, _record("a_noarch", subdir="noarch"), "a.json")
    explicit, inventory = exporter.build_lock(tmp_path)
    lines = explicit.splitlines()
    assert lines[0] == "@EXPLICIT"
    assert len(lines) == 3
    assert lines[1].startswith(
        "https://conda.anaconda.org/conda-forge/noarch/a_noarch-"
    )
    assert all(line.endswith("#" + "a" * 64) for line in lines[1:])
    assert inventory["schema"] == "kernelgecp.conda-runtime-lock.v1"
    assert inventory["platform"] == "linux-64"
    assert inventory["package_count"] == 2
    assert inventory["explicit_sha256"] == hashlib.sha256(explicit.encode()).hexdigest()
    package = inventory["packages"][0]
    assert package["depends"] == ["numpy >=2,<3", "python >=3.12,<3.13"]
    assert package["license"] == "BSD-3-Clause"
    assert package["build"] == "py312_example_0"
    assert "/private" not in json.dumps(inventory)
    assert str(tmp_path) not in explicit + json.dumps(inventory)
    assert exporter.build_lock(tmp_path) == (explicit, inventory)


@pytest.mark.parametrize(
    "key", ["sha256", "url", "version", "build", "license", "depends"]
)
def test_missing_metadata_never_produces_unverified_lock(
    tmp_path: Path, key: str
) -> None:
    record = _record()
    del record[key]
    _write(tmp_path, record)
    with pytest.raises(exporter.CondaLockError):
        exporter.build_lock(tmp_path)


@pytest.mark.parametrize(
    "url",
    [
        "http://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda",
        "https://user:password@conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda",
        "https://private.invalid/conda-forge/linux-64/example-1.2.3-py312_example_0.conda",
        "file:///private/cache/example-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/private-channel/linux-64/example-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda?token=secret",
        "https://conda.anaconda.org/t/secret/conda-forge/linux-64/example-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda#token=secret",
        "https://conda.anaconda.org/conda-forge/linux-64/other-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/conda-forge/osx-arm64/example-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda\n",
        "https://conda.anaconda.org/conda-forge/\nlinux-64/example-1.2.3-py312_example_0.conda",
        "https://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda#",
        "https://conda.anaconda.org/conda-forge/linux-64/example-1.2.3-py312_example_0.conda?",
    ],
)
def test_private_token_or_mismatched_urls_are_rejected(
    tmp_path: Path, url: str
) -> None:
    record = _record()
    record["url"] = url
    _write(tmp_path, record)
    with pytest.raises(exporter.CondaLockError):
        exporter.build_lock(tmp_path)


@pytest.mark.parametrize("digest", [None, "", "a" * 32, "z" * 64, 123, True])
def test_sha256_is_required_not_an_md5_or_placeholder(
    tmp_path: Path, digest: Any
) -> None:
    record = _record()
    record["sha256"] = digest
    _write(tmp_path, record)
    with pytest.raises(exporter.CondaLockError, match="SHA256"):
        exporter.build_lock(tmp_path)


@pytest.mark.parametrize("number", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_values_even_in_unused_metadata_are_rejected(
    tmp_path: Path, number: float
) -> None:
    record = _record()
    record["unused"] = number
    _write(tmp_path, record)
    with pytest.raises(exporter.CondaLockError, match="nonfinite"):
        exporter.build_lock(tmp_path)


def test_duplicate_keys_names_and_mixed_platforms_are_rejected(tmp_path: Path) -> None:
    _write(tmp_path, _record())
    path = tmp_path / "conda-meta" / "example.json"
    original = path.read_text()
    path.write_text(original[:-1] + ', "sha256": "' + "b" * 64 + '"}')
    with pytest.raises(exporter.CondaLockError, match="duplicate metadata key"):
        exporter.build_lock(tmp_path)
    path.write_text(original)
    _write(tmp_path, _record(), "duplicate.json")
    with pytest.raises(exporter.CondaLockError, match="duplicate installed"):
        exporter.build_lock(tmp_path)
    (path.parent / "duplicate.json").unlink()
    _write(tmp_path, _record("mac", subdir="osx-arm64"), "mac.json")
    with pytest.raises(exporter.CondaLockError, match="one native platform"):
        exporter.build_lock(tmp_path)


def test_json_exponent_overflow_is_not_a_finite_metadata_number(tmp_path: Path) -> None:
    _write(tmp_path, _record())
    path = tmp_path / "conda-meta" / "example.json"
    path.write_text(path.read_text()[:-1] + ', "unused": 1e999}')
    with pytest.raises(exporter.CondaLockError, match="nonfinite"):
        exporter.build_lock(tmp_path)


@pytest.mark.parametrize("field", ["depends", "license"])
def test_inventory_text_cannot_leak_paths_or_tokens(tmp_path: Path, field: str) -> None:
    record = _record()
    record[field] = ["/private/dependency"] if field == "depends" else "token=secret"
    _write(tmp_path, record)
    with pytest.raises(exporter.CondaLockError):
        exporter.build_lock(tmp_path)


def test_empty_or_symlinked_metadata_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(exporter.CondaLockError):
        exporter.build_lock(tmp_path)
    directory = tmp_path / "conda-meta"
    directory.mkdir()
    with pytest.raises(exporter.CondaLockError):
        exporter.build_lock(tmp_path)
    source = tmp_path / "elsewhere.json"
    source.write_text(json.dumps(_record()))
    (directory / "linked.json").symlink_to(source)
    with pytest.raises(exporter.CondaLockError, match="regular JSON"):
        exporter.build_lock(tmp_path)


def test_cli_outputs_no_files_when_validation_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    record = _record()
    del record["sha256"]
    _write(tmp_path, record)
    explicit, inventory = tmp_path / "lock.txt", tmp_path / "inventory.json"
    monkeypatch.setattr(
        exporter.sys,
        "argv",
        [
            "export_conda_lock",
            "--prefix",
            str(tmp_path),
            "--explicit",
            str(explicit),
            "--inventory",
            str(inventory),
        ],
    )
    with pytest.raises(SystemExit) as result:
        exporter.main()
    assert result.value.code == 1
    assert not explicit.exists()
    assert not inventory.exists()
