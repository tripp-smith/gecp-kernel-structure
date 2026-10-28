"""Freeze boundaries and safe immutable parent retention for fresh QMC data."""

from __future__ import annotations

import hashlib
import io
import tarfile
from pathlib import Path

import pytest

from experiments.revalidate_cthyb_dmft import (
    ROOT,
    extract_parent,
    inventory,
    load_protocol,
    validate_parent,
)


def test_protocol_keeps_targets_and_declares_completely_new_64_chains() -> None:
    protocol = load_protocol(ROOT / "experiments/cthyb_revalidation.toml")
    assert protocol["final_replicas"] == 64
    assert protocol["cycles"] == 400000
    assert protocol["seed_base"] == 28100000
    assert protocol["max_standard_error"] == 3e-4
    assert protocol["agreement_absolute"] == 1e-3 and protocol["agreement_sigma"] == 3


@pytest.mark.parametrize(
    "before,after",
    [
        ("final_replicas = 64", "final_replicas = 48"),
        ("seed_base = 28100000", "seed_base = 27100000"),
        ("max_standard_error = 0.0003", "max_standard_error = 0.0004"),
        ("job_minutes = 300", "job_minutes = 600"),
    ],
)
def test_protocol_rejects_posthoc_tuning(
    tmp_path: Path, before: str, after: str
) -> None:
    text = (ROOT / "experiments/cthyb_revalidation.toml").read_text()
    path = tmp_path / "changed.toml"
    path.write_text(text.replace(before, after))
    with pytest.raises(ValueError, match="predeclared"):
        load_protocol(path)


def test_hash_bound_parent_extract_and_original_miss(tmp_path: Path) -> None:
    protocol = load_protocol(ROOT / "experiments/cthyb_revalidation.toml")
    parent, audit = extract_parent(
        ROOT / "experiments/data/cthyb_45d42c7_inconclusive.tar.gz",
        tmp_path / "unpacked",
        protocol["parent_bundle_sha256"],
    )
    proof = validate_parent(parent, audit)
    assert proof["generation_revision"].startswith("45d42c7")
    assert len(proof["inventory"]) > 100
    # Repeated extraction cannot replace the original raw files.
    previous = inventory(parent)
    extract_parent(
        ROOT / "experiments/data/cthyb_45d42c7_inconclusive.tar.gz",
        tmp_path / "unpacked",
        protocol["parent_bundle_sha256"],
    )
    assert inventory(parent) == previous
    (parent / "result.json").write_text("{}")
    with pytest.raises(ValueError, match="changed"):
        validate_parent(parent, audit)


@pytest.mark.parametrize(
    "name,kind",
    [
        ("../outside", "file"),
        ("/absolute", "file"),
        ("parent/link", "symlink"),
        ("other/file", "file"),
    ],
)
def test_bundle_rejects_unsafe_members(tmp_path: Path, name: str, kind: str) -> None:
    bundle = tmp_path / "bad.tar.gz"
    with tarfile.open(bundle, "w:gz") as archive:
        info = tarfile.TarInfo(name)
        if kind == "symlink":
            info.type = tarfile.SYMTYPE
            info.linkname = "/outside"
            archive.addfile(info)
        else:
            info.size = 1
            archive.addfile(info, io.BytesIO(b"x"))
    digest = hashlib.sha256(bundle.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="unsafe"):
        extract_parent(bundle, tmp_path / "unpacked", digest)
    assert not (tmp_path / "outside").exists()
