import hashlib
import io
import json
import tarfile

import numpy as np
import pytest

from scripts import diagnose_spectral_projection as diag


def test_wrong_archive_hash_writes_failure(tmp_path, monkeypatch):
    archive = tmp_path / "wrong.tar.gz"
    archive.write_bytes(b"incorrect")
    monkeypatch.setattr(diag, "ARCHIVE_BYTES", archive.stat().st_size)
    output = tmp_path / "result.json"
    assert diag.main(["--archive", str(archive), "--output", str(output)]) == 1
    assert json.loads(output.read_text())["reason"] == "archive SHA-256 mismatch"


def test_rejects_unsafe_archive_members(tmp_path, monkeypatch):
    archive = tmp_path / "unsafe.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        entry = tarfile.TarInfo("../escape")
        entry.size = 1
        tar.addfile(entry, io.BytesIO(b"x"))
    monkeypatch.setattr(diag, "ARCHIVE_BYTES", archive.stat().st_size)
    monkeypatch.setattr(
        diag, "ARCHIVE_SHA256", hashlib.sha256(archive.read_bytes()).hexdigest()
    )
    with pytest.raises(ValueError, match="unsafe archive"):
        diag.load_archive(archive)


def test_rejects_nonfinite_saved_matrix():
    data = io.BytesIO()
    a = np.ones((65, 23))
    a[0, 0] = np.nan
    np.savez(data, observation_design=a, projection=np.eye(65))
    with pytest.raises(ValueError, match="invalid design"):
        diag.read_matrices(data.getvalue(), 23)


def test_high_precision_recovers_ill_conditioned_span():
    # Orthogonal complement is known exactly despite a tiny second singular value.
    a = np.array([[1.0, 1.0], [0.0, 1e-10], [0.0, 0.0]])
    p = np.diag([1.0, 1.0, 0.0])
    result = diag.high_precision(a, p)
    assert result["rank"] == 2
    assert result["within_tolerance"]
    assert result["max_abs"] < 1e-14
