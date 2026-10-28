"""Dataset chronology, leakage isolation, and blocked covariance controls."""

from __future__ import annotations

import copy

import numpy as np
import pytest
from scipy.linalg import cho_factor, cho_solve

from experiments import covariance_benchmark as benchmark
from kernelgecp import KernelCholeskyFactor


def _synthetic_data():
    rng = np.random.default_rng(219)
    return {
        "values": rng.normal(size=(benchmark.TEST_END + 10, 6)),
        "timestamps": [str(index) for index in range(benchmark.TEST_END + 10)],
        "ordinals": np.arange(benchmark.TEST_END + 10, dtype=np.int64),
    }


def test_prefix_normalization_excludes_later_training_calibration_and_test_data():
    data = _synthetic_data()
    before = benchmark.split_data(data, 512, "calibration")
    changed = copy.deepcopy(data)
    changed["values"][512 : benchmark.TRAIN_END] *= 100
    changed["values"][benchmark.CAL_END :] = 1e100
    after = benchmark.split_data(changed, 512, "calibration")
    assert before["normalization"] == after["normalization"]
    np.testing.assert_array_equal(before["x"], after["x"])
    np.testing.assert_array_equal(before["y"], after["y"])
    np.testing.assert_array_equal(before["query_x"], after["query_x"])
    assert before["query_first"] == str(benchmark.TRAIN_END)
    test = benchmark.split_data(data, 512, "test")
    assert test["query_first"] == str(benchmark.CAL_END)
    np.testing.assert_allclose(before["x"].mean(axis=0), 0, atol=1e-15)
    np.testing.assert_allclose(before["x"].std(axis=0), 1, atol=1e-15)


def test_parser_keeps_complete_cases_only_in_declared_features_and_target():
    header = "No,year,month,day,hour,PM2.5,TEMP,PRES,DEWP,RAIN,WSPM,station,unused\n"
    rows = (
        "1,2013,3,1,0,4,-1,1023,-18,0,4,Aotizhongxin,NA\n"
        "2,2013,3,1,1,NA,-1,1023,-18,0,4,Aotizhongxin,1\n"
        "3,2013,3,1,2,8,-2,1024,-17,0,5,Aotizhongxin,2\n"
    )
    parsed = benchmark.parse_station((header + rows).encode())
    assert parsed["excluded_missing"] == 1
    np.testing.assert_array_equal(parsed["ordinals"], [1, 3])
    np.testing.assert_array_equal(parsed["values"][:, -1], [4, 8])
    with pytest.raises(ValueError, match="chronological"):
        benchmark.parse_station((header + rows + rows.splitlines()[0]).encode())
    with pytest.raises(ValueError, match="station"):
        benchmark.parse_station(
            (header + rows.replace("Aotizhongxin", "other")).encode()
        )


def test_wrong_archive_hash_fails_before_parsing_or_network(tmp_path):
    archive = tmp_path / "wrong.zip"
    archive.write_bytes(b"not the attributed archive")
    with pytest.raises(ValueError, match="hash mismatch"):
        benchmark.load_data(archive)


def test_blocked_original_operator_matches_independent_pairwise_formula(monkeypatch):
    rng = np.random.default_rng(183)
    x = rng.normal(size=(257, 5))
    expected = np.zeros((257, 257))
    for dimension in range(5):
        expected += (x[:, dimension, None] - x[None, :, dimension]) ** 2
    expected = np.exp(-0.5 * expected / 1.7**2) + 0.2 * np.eye(257)
    original = benchmark.rbf
    calls = []

    def checked(left, right, lengthscale):
        assert len(left) <= 128
        assert left.ndim == right.ndim == 2
        calls.append(len(left))
        return original(left, right, lengthscale)

    monkeypatch.setattr(benchmark, "rbf", checked)
    operation = benchmark.original_operator(x, 1.7, 0.2)
    rhs = rng.normal(size=(257, 3))
    np.testing.assert_allclose(operation @ rhs, expected @ rhs, atol=4e-14)
    assert calls == [128, 128, 1]


def test_predictive_standard_deviation_keeps_unrepresented_prior_uncertainty():
    rng = np.random.default_rng(89)
    x = rng.normal(size=(19, 5))
    query = rng.normal(size=(11, 5))
    split = {"x": x, "query_x": query, "y": rng.normal(size=19)}
    factor = KernelCholeskyFactor.from_columns(
        np.ones(19), lambda j: benchmark.rbf(x, x[[j]], 1.5)[:, 0], max_rank=5
    )
    mean, std, _, _ = benchmark.predictions(factor, split, 1.5, 0.2)
    extension = factor.extend(benchmark.rbf(query, x[list(factor.indices)], 1.5))
    cross = extension.features @ factor.factor.T
    dense = factor.factor @ factor.factor.T + 0.2 * np.eye(19)
    chol = cho_factor(dense, lower=True)
    expected_mean = cross @ cho_solve(chol, split["y"])
    expected_variance = 1.2 - np.sum(cross * cho_solve(chol, cross.T).T, axis=1)
    np.testing.assert_allclose(mean, expected_mean, atol=2e-14)
    np.testing.assert_allclose(std**2, expected_variance, atol=2e-14)
    assert np.max(extension.features**2) < np.max(std**2)


def test_frozen_validation_refuses_test_contaminated_or_changed_source(monkeypatch):
    monkeypatch.setattr(
        benchmark,
        "provenance",
        lambda: {
            "source_hashes": {"source": "abc"},
            "versions": {"numpy": "2.2.6"},
        },
    )
    frozen = {
        "schema": "kernelgecp.covariance-calibration.v1",
        "status": "calibrated",
        "test_evaluated": True,
        "calibration_only": True,
    }
    with pytest.raises(ValueError, match="test observations"):
        benchmark.validate_frozen(frozen, canonical=False)
    frozen["test_evaluated"] = False
    frozen["dataset"] = {
        "outer_sha256": benchmark.OUTER_SHA256,
        "csv_sha256": benchmark.CSV_SHA256,
    }
    frozen["provenance"] = {"source_hashes": {"source": "changed"}}
    with pytest.raises(ValueError, match="source differs"):
        benchmark.validate_frozen(frozen, canonical=False)


def test_scientific_freeze_binds_executed_runtime_not_unrelated_development_lock(
    tmp_path, monkeypatch
):
    import subprocess

    monkeypatch.setattr(benchmark, "ROOT", tmp_path)
    for name in (*benchmark.SOURCE_FILES, "uv.lock"):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"frozen {name}")
    monkeypatch.setattr(
        benchmark.subprocess,
        "check_output",
        lambda args, **kwargs: "a" * 40 if args[1] == "rev-parse" else "",
    )
    monkeypatch.setattr(
        benchmark.subprocess,
        "run",
        lambda args, **kwargs: subprocess.CompletedProcess(args, 0),
    )
    monkeypatch.setattr(
        benchmark.importlib.metadata, "version", benchmark.LOCKED_VERSIONS.__getitem__
    )
    original = benchmark.provenance()
    frozen = {
        "schema": "kernelgecp.covariance-calibration.v1",
        "status": "calibrated",
        "test_evaluated": False,
        "calibration_only": True,
        "dataset": {
            "outer_sha256": benchmark.OUTER_SHA256,
            "csv_sha256": benchmark.CSV_SHA256,
        },
        "provenance": original,
        "standard_deviation_definition": benchmark.STD_DEFINITION,
        "timing_scopes": benchmark.TIMING_SCOPES,
    }
    (tmp_path / "uv.lock").write_text("unrelated development-only update")
    current = benchmark.provenance()
    assert current["source_hashes"] == original["source_hashes"]
    assert (
        current["development_lock_sha256_context_only"]
        != original["development_lock_sha256_context_only"]
    )
    benchmark.validate_frozen(frozen, canonical=False)
    (tmp_path / "environments/gp/uv.lock").write_text(
        "changed executed scientific runtime"
    )
    with pytest.raises(ValueError, match="source differs"):
        benchmark.validate_frozen(frozen, canonical=False)
