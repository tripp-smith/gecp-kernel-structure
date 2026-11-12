"""Independent bounded saved-factor and structural failure checks."""

import copy
import importlib.util
import io
import json
import struct
import zipfile
from pathlib import Path

import numpy as np
import pytest

from scripts import release_covariance_provider as replay


def fixture():
    k01 = np.exp(-0.5)
    k02 = np.exp(-2.0)
    diagonal = np.sqrt(1 - k01**2)
    root = np.array([[1.0, 0.0], [k01, diagonal], [k02, (k01 - k02 * k01) / diagonal]])
    x = np.array([[0.0], [1.0], [2.0]])
    query = np.array([[0.5], [1.5]])
    y = np.array([0.3, -0.2, 0.4])
    noise = 0.3
    matrix = root @ root.T + noise * np.eye(3)
    kernel = np.exp(-0.5 * (query - x[:2].T) ** 2)
    features = np.linalg.solve(root[:2], kernel.T).T
    cross = features @ root.T
    solution = np.linalg.solve(matrix, y)
    mean = cross @ solution
    std = np.sqrt(
        1 + noise - np.sum(cross * np.linalg.solve(matrix, cross.T).T, axis=1)
    )
    case = {
        "root": root,
        "pivots": np.array([0, 1]),
        "cg_rhs": y,
        "approximate_solution": solution,
        "mean": mean,
        "std": std,
    }
    provider = {
        "root": root.copy(),
        "rhs": y.copy(),
        "solution": solution.copy(),
        "mean": mean.copy(),
        "std": std.copy(),
        "expected_mean": mean.copy(),
        "expected_std": std.copy(),
    }
    return provider, case, {"x": x, "query_x": query, "y": y}, noise


def test_independent_dense_reference_and_corruption():
    provider, case, split, noise = fixture()
    result = replay.calculate_reference(provider, case, split, noise, 1.0)
    matrix = case["root"] @ case["root"].T + noise * np.eye(3)
    independent = -0.5 * (
        split["y"] @ np.linalg.solve(matrix, split["y"])
        + np.linalg.slogdet(matrix)[1]
        + 3 * np.log(2 * np.pi)
    )
    assert result["log_probability"] == pytest.approx(independent, abs=1e-12)
    for key in (
        "solution",
        "mean",
        "std",
        "expected_mean",
        "expected_std",
        "rhs",
        "root",
    ):
        damaged = copy.deepcopy(provider)
        damaged[key].flat[0] += 0.01
        with pytest.raises(ValueError):
            replay.calculate_reference(damaged, case, split, noise, 1.0)


def test_duplicate_pivots_and_nonpositive_noise_rejected():
    provider, case, split, noise = fixture()
    case["pivots"][1] = 0
    with pytest.raises(ValueError, match="pivots"):
        replay.calculate_reference(provider, case, split, noise, 1.0)
    with pytest.raises(ValueError, match="controls"):
        replay.calculate_reference(provider, case, split, 0.0, 1.0)


def test_strict_json_and_npz(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"a":1,"a":2}')
    with pytest.raises(ValueError, match="duplicate"):
        replay.read_json(path)
    path.write_text(json.dumps({"a": float("nan")}))
    with pytest.raises(ValueError, match="nonfinite"):
        replay.read_json(path)
    archive = tmp_path / "arrays.npz"
    np.savez_compressed(archive, value=np.ones(2), extra=np.ones(1))
    with pytest.raises(ValueError, match="keyset"):
        replay.read_arrays(archive, {"value"})
    np.savez_compressed(archive, value=np.array([np.inf]))
    with pytest.raises(ValueError, match="invalid"):
        replay.read_arrays(archive, {"value"})


def test_missing_input_does_not_download(tmp_path, monkeypatch):
    import urllib.request

    monkeypatch.setattr(
        urllib.request, "urlopen", lambda *a, **k: pytest.fail("network used")
    )
    with pytest.raises(ValueError):
        replay.replay(
            tmp_path, tmp_path, tmp_path / "missing.zip", tmp_path / "out.json"
        )


def test_optional_provider_is_lazy():
    path = Path(replay.__file__)
    source = path.read_text()
    assert 'importlib.import_module("torch")' in source
    assert "import torch\n" not in source


@pytest.mark.skipif(
    importlib.util.find_spec("gpytorch") is None,
    reason="optional real GPyTorch runtime",
)
def test_actual_provider_small_fixture():
    provider, case, split, noise = fixture()
    result = replay.actual_provider(
        replay.calculate_reference(provider, case, split, noise, 1.0), noise
    )
    assert result["device"] == "cpu"
    assert max(result["errors"].values()) <= 1e-9


def test_forged_npy_shape_rejected_before_numpy_load(tmp_path, monkeypatch):
    header = io.BytesIO()
    np.lib.format.write_array_header_1_0(
        header, {"descr": "<f8", "fortran_order": False, "shape": (2**40,)}
    )
    archive = tmp_path / "forged.npz"
    with zipfile.ZipFile(archive, "w") as out:
        out.writestr("value.npy", header.getvalue())
    monkeypatch.setattr(np, "load", lambda *a, **k: pytest.fail("allocation attempted"))
    with pytest.raises(ValueError, match="dimensions|allocation"):
        replay.read_arrays(archive, {"value"})


def test_npy_truncation_rejected_before_numpy_load(tmp_path, monkeypatch):
    header = io.BytesIO()
    np.lib.format.write_array_header_1_0(
        header, {"descr": "<f8", "fortran_order": False, "shape": (2,)}
    )
    archive = tmp_path / "truncated.npz"
    with zipfile.ZipFile(archive, "w") as out:
        out.writestr("value.npy", header.getvalue() + bytes(8))
    monkeypatch.setattr(np, "load", lambda *a, **k: pytest.fail("allocation attempted"))
    with pytest.raises(ValueError, match="payload size"):
        replay.read_arrays(archive, {"value"})


def test_npy_huge_header_rejected_before_parser_or_load(tmp_path, monkeypatch):
    archive = tmp_path / "header.npz"
    with zipfile.ZipFile(archive, "w") as out:
        out.writestr("value.npy", b"\x93NUMPY\x02\x00" + struct.pack("<I", 2**32 - 1))
    monkeypatch.setattr(np, "load", lambda *a, **k: pytest.fail("allocation attempted"))
    monkeypatch.setattr(
        np.lib.format,
        "read_array_header_2_0",
        lambda *a, **k: pytest.fail("header parser called"),
    )
    with pytest.raises(ValueError, match="header size"):
        replay.read_arrays(archive, {"value"})
