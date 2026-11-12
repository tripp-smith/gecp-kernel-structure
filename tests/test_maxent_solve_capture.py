import hashlib

import numpy as np
import pytest

from scripts.providers.maxent_solve_capture import capture_failed_solve


def _failure():
    J = np.array([[1.0, 1.0], [1.0, 1.0]])
    Id = np.eye(2)
    mu = 0.0
    v = np.array([2.0, 3.0])
    alpha = 0.1  # noqa: F841
    i = 7  # noqa: F841
    return np.linalg.solve(J + mu * Id, v)  # noqa: F841


def test_success_delegates_and_restores(tmp_path):
    previous = np.linalg.solve
    metadata = {}
    path = tmp_path / "capture.npz"
    with capture_failed_solve(path, metadata):
        np.testing.assert_array_equal(
            np.linalg.solve(np.eye(2), [2.0, 3.0]), [2.0, 3.0]
        )
    assert np.linalg.solve is previous and not path.exists() and not metadata


def test_singular_capture(tmp_path):
    metadata = {}
    path = tmp_path / "capture.npz"
    with capture_failed_solve(path, metadata), pytest.raises(np.linalg.LinAlgError):
        _failure()
    assert metadata["status"] == "captured"
    assert metadata["iteration"] == 7 and metadata["scalar_alpha"] == 0.1
    assert metadata["npz_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(metadata["caller_source_sha256"]) == 64
    with np.load(path, allow_pickle=False) as arrays:
        assert set(arrays.files) == {"matrix", "rhs", "J", "Id", "mu", "v"}
        for key in arrays.files:
            assert metadata["arrays"][key]["finite"]
            assert (
                metadata["arrays"][key]["sha256"]
                == hashlib.sha256(
                    np.ascontiguousarray(arrays[key]).tobytes()
                ).hexdigest()
            )


def test_cap_and_dtype_failure_keep_original(monkeypatch, tmp_path):
    original = np.linalg.solve

    def bad(*args, **kwargs):
        raise np.linalg.LinAlgError("original")

    monkeypatch.setattr(np.linalg, "solve", bad)
    for index, matrix in enumerate((np.zeros((129, 129)), np.array([[object()]]))):
        info = {}
        with (
            capture_failed_solve(tmp_path / f"{index}.npz", info),
            pytest.raises(np.linalg.LinAlgError, match="original"),
        ):
            np.linalg.solve(matrix, np.ones(len(matrix)))
        assert (
            info["status"] == "capture_failed"
            and "unsupported" in info["capture_error"]
        )
        assert np.linalg.solve is bad
    monkeypatch.setattr(np.linalg, "solve", original)


def test_nested_context_and_write_failure(tmp_path):
    prior = np.linalg.solve
    target = tmp_path / "existing.npz"
    target.write_bytes(b"existing")
    outer, inner = {}, {}
    with capture_failed_solve(tmp_path / "outer.npz", outer):
        outer_solve = np.linalg.solve
        with capture_failed_solve(target, inner), pytest.raises(np.linalg.LinAlgError):
            _failure()
        assert np.linalg.solve is outer_solve
    assert np.linalg.solve is prior
    assert target.read_bytes() == b"existing"
    assert (
        inner["status"] == "capture_failed"
        and "FileExistsError" in inner["capture_error"]
    )


def test_total_byte_cap_checked_before_write(monkeypatch, tmp_path):
    def bad(*args, **kwargs):
        raise np.linalg.LinAlgError("singular sentinel")

    def call_with_large_locals():
        J = np.zeros((128, 128), dtype=np.complex128)
        Id = J.copy()
        v = J.copy()  # noqa: F841
        mu = J.copy()  # noqa: F841
        return np.linalg.solve(J, Id)  # noqa: F841

    monkeypatch.setattr(np.linalg, "solve", bad)
    info = {}
    path = tmp_path / "too-big.npz"
    with (
        capture_failed_solve(path, info),
        pytest.raises(np.linalg.LinAlgError, match="singular sentinel"),
    ):
        call_with_large_locals()
    assert info["status"] == "capture_failed"
    assert "byte cap" in info["capture_error"]
    assert not path.exists()
