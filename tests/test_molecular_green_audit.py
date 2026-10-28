"""Independent molecular oracle checks and intentional acceptance rejection."""

from __future__ import annotations

import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import expit

import experiments.audit_molecular_green as auditor


def test_matrix_causality_uses_eigenvalues_not_diagonal_signs() -> None:
    imaginary = np.array([[-1.0, 2.0], [2.0, -1.0]])
    values = (1j * imaginary)[None, None, None]
    assert auditor.causality(values, np.array([0])) == pytest.approx(1.0)
    assert auditor.causality(
        -1j * np.eye(2)[None, None, None], np.array([0])
    ) == pytest.approx(-1.0)


def test_restricted_density_preserves_two_spin_factor() -> None:
    labels = np.array([-1, 0])
    tau = np.array([1.0, 19.0])
    fields = {name: np.zeros((2, 1, 1, 2, 2), complex) for name in auditor.FIELDS}
    endpoints = np.zeros((2, 1, 1, 2, 2), complex)
    endpoints[-1, 0, 0] = np.diag([-0.75, -0.25])
    native = SimpleNamespace(
        density=2 + 0j,
        time=lambda q, f: np.zeros((len(q), 1, 1, 2, 2)),
        frequency=lambda q, f: np.zeros((len(q), 1, 1, 2, 2)),
    )
    metrics = auditor.worker_metrics(
        {**fields, "tau": tau, "indices": labels, "green_endpoints": endpoints}, native
    )
    assert metrics["density"] == 2
    assert metrics["density_error"] == 0


def test_finite_thermal_trace_matches_independent_noninteracting_formula() -> None:
    energies = np.array([0.2, 0.7])
    mu = 0.3
    trace = auditor.FiniteTrace(np.diag(energies), np.zeros((2, 2, 2, 2)), mu, 0.4)
    occupied = expit(-20 * (energies - mu))
    assert trace.density == pytest.approx(2 * occupied.sum(), abs=1e-13)
    assert trace.energy == pytest.approx(
        0.4 + 2 * np.sum(energies * occupied), abs=1e-13
    )
    tau = np.array([0.0, 0.1, 3.0, 19.9, 20.0])
    expected = -np.exp(-tau[:, None] * (energies - mu)) * expit(20 * (energies - mu))
    np.testing.assert_allclose(
        np.diagonal(trace.time(tau)[:, 0, 0], axis1=-2, axis2=-1), expected, atol=2e-14
    )
    indices = np.array([-32, -1, 0, 31])
    expected = 1 / (1j * (2 * indices[:, None] + 1) * np.pi / 20 - (energies - mu))
    np.testing.assert_allclose(
        np.diagonal(trace.frequency(indices)[:, 0, 0], axis1=-2, axis2=-1),
        expected,
        atol=3e-13,
    )
    for first in range(4):
        for second in range(4):
            a, b = trace.operators[first], trace.operators[second]
            np.testing.assert_array_equal(
                a @ b.T + b.T @ a, np.eye(16) if first == second else np.zeros((16, 16))
            )


@pytest.mark.parametrize("payload", ['{"x":NaN}', '{"x":1e999}', '{"x":1,"x":2}', "[]"])
def test_receipt_parser_rejects_nonfinite_duplicate_and_nonobject(
    tmp_path: Path, payload: str
) -> None:
    path = tmp_path / "bad.json"
    path.write_text(payload)
    with pytest.raises(ValueError):
        auditor.read(path)


@pytest.mark.parametrize("kind", ["object", "nonfinite", "missing"])
def test_numeric_archive_rejects_unsafe_data(tmp_path: Path, kind: str) -> None:
    path = tmp_path / "bad.npz"
    if kind == "object":
        values = np.array([object()], dtype=object)
    elif kind == "nonfinite":
        values = np.array([np.nan])
    else:
        values = np.array([1.0])
    np.savez(path, wrong=values) if kind == "missing" else np.savez(
        path, expected=values
    )
    with pytest.raises(ValueError):
        auditor.numeric_archive(path, {"expected"})


def test_audit_rejects_pilot_without_consulting_reported_pass(tmp_path: Path) -> None:
    (tmp_path / "report.json").write_text(
        json.dumps(
            {
                "schema": "kernelgecp.molecular-green.v1",
                "pilot": True,
                "all_mandatory_passed": True,
            }
        )
    )
    with pytest.raises(ValueError, match="pilot"):
        auditor.audit(tmp_path, tmp_path)


def test_auditor_never_imports_producer_truth_or_trusts_its_gate() -> None:
    source = inspect.getsource(auditor)
    for forbidden in (
        "from experiments.molecular_green import",
        "from experiments.molecular_native import",
        "from experiments.molecular_reference import",
        "read_green_tau(",
        'report["all_mandatory_passed"]',
        'report["deterministic_repetitions"]',
    ):
        assert forbidden not in source


@pytest.mark.parametrize("name", ["../outside", "/absolute", "x\\y"])
def test_artifact_paths_cannot_escape_bundle(tmp_path: Path, name: str) -> None:
    with pytest.raises(ValueError, match="unsafe"):
        auditor.relative(tmp_path, name)


def test_fresh_final_time_grid_excludes_retired_queries() -> None:
    native = SimpleNamespace(tau=np.linspace(0, 20, 101))
    fresh = auditor.expected_time(native, "validation")
    old = auditor.expected_time(native, "regression")
    assert fresh.shape == (96,)
    assert not np.intersect1d(fresh, old).size


def test_audit_rejects_reused_final_frequency_before_artifact_access(
    tmp_path: Path,
) -> None:
    config = {
        "schema_version": 2,
        "validation_tau_rule": "gauss_legendre_96",
        "regression_tau_rule": "native_mod4_3_plus_gauss64",
        "retired_validation_source": "6050ddd33af30deafe0d99057ea5d8a277c411a9",
        "calibration_indices": [0],
        "regression_indices": [1],
        "validation_indices": [1],
    }
    (tmp_path / "report.json").write_text(
        json.dumps(
            {
                "schema": "kernelgecp.molecular-green.v2",
                "pilot": False,
                "config": config,
            }
        )
    )
    with pytest.raises(ValueError, match="final frequency query overlap"):
        auditor.audit(tmp_path, tmp_path)
