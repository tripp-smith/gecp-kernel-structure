"""Reject corrupted development evidence; this is not production acceptance."""

import hashlib
import json
import shutil
import tarfile
from pathlib import Path

import pytest

from scripts import audit_maxent_lm_diagnostic as audit

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def reports(tmp_path_factory):
    target = tmp_path_factory.mktemp("lm-evidence")
    metadata = json.loads(
        (ROOT / "data/spectroscopy/lm-linux-eeb5523.json").read_text()
    )
    archive = ROOT / metadata["archive"]["path"]
    assert (
        hashlib.sha256(archive.read_bytes()).hexdigest()
        == metadata["archive"]["sha256"]
    )
    with tarfile.open(archive) as source:
        source.extractall(target, filter="data")
    return target / "reports"


def test_independent_capture_analysis_records_failures_without_claiming_acceptance(
    reports,
):
    result = audit.audit(reports)
    assert "not production acceptance" in result["scope"]
    captures = [
        row["captured_system"] for row in result["runs"] if "captured_system" in row
    ]
    assert len(captures) == 4
    for row in captures:
        assert row["exact_binary_matrix_rank"] == row["shape"][0]
        assert row["exact_binary_augmented_rank"] == row["shape"][0]


@pytest.mark.parametrize("mutation", ["missing", "nonstandard_json", "capture_bytes"])
def test_incomplete_or_corrupted_reports_are_rejected(reports, tmp_path, mutation):
    copied = tmp_path / "reports"
    shutil.copytree(reports, copied)
    if mutation == "missing":
        shutil.rmtree(next(copied.iterdir()))
    elif mutation == "nonstandard_json":
        path = next(copied.glob("*/result.json"))
        path.write_text('{"bad": NaN}')
    else:
        path = next(copied.glob("*/result.failure.npz"))
        path.write_bytes(path.read_bytes()[:-3] + b"bad")
    with pytest.raises(ValueError):
        audit.audit(copied)
