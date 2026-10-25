"""The executable same-workload example is exercised, not just its protocol."""

import importlib.util

import numpy as np
import pytest

from examples.dlr_interoperability import run_demo
from examples.sparse_ir_physical_recipe import (
    _exact_iw,
    _exact_tau,
    _spectral_fixture,
    run_recipe,
)


def test_public_gecp_physical_interoperability_example():
    record = run_demo(("kernelgecp",))["records"][0]
    assert record["method"] == "kernelgecp"
    assert record["accuracy_passed"]
    assert record["solver_converged"]
    assert record["solver_iterations"] > 1


@pytest.mark.parametrize("method", ["pydlr", "sparse_ir"])
def test_actual_optional_physical_interoperability_example(method):
    if importlib.util.find_spec(method) is None:
        pytest.skip(f"optional {method} is not installed")
    record = run_demo((method,))["records"][0]
    assert record["method"] == method
    assert record["accuracy_passed"]
    assert record["solver_converged"]


@pytest.mark.parametrize("methods", [(), ("unknown",), ("kernelgecp", "kernelgecp")])
def test_example_rejects_invalid_provider_selection(methods):
    with pytest.raises(ValueError):
        run_demo(methods)


def test_standalone_upstream_recipe_uses_actual_sparse_ir():
    if importlib.util.find_spec("sparse_ir") is None:
        pytest.skip("optional sparse_ir is not installed")
    result = run_recipe()
    assert result["accuracy_passed"]
    assert result["rank"] > 0
    assert result["time_error"] <= 1e-8
    assert result["frequency_error"] <= 1e-8


def test_upstream_recipe_has_normalized_positive_matrix_spectrum():
    beta, wmax, energies, residues = _spectral_fixture()
    assert np.all(np.abs(energies) <= wmax)
    np.testing.assert_array_equal(residues, residues.conj().transpose(0, 2, 1))
    assert np.min(np.linalg.eigvalsh(residues)) > 0
    np.testing.assert_array_equal(residues.sum(axis=0), np.eye(2))
    endpoints = _exact_tau(np.array([0.0, beta]), beta, energies, residues)
    np.testing.assert_allclose(endpoints.sum(axis=0), -np.eye(2), rtol=0, atol=1e-15)


def test_upstream_recipe_time_truth_integrates_to_frequency_truth():
    beta, _, energies, residues = _spectral_fixture()
    nodes, weights = np.polynomial.legendre.leggauss(96)
    tau = 0.5 * beta * (nodes + 1)
    reduced_n = np.array([-7, -3, -1, 1, 3, 7])
    phase = np.exp(1j * reduced_n[:, None] * np.pi * tau / beta)
    integrated = np.einsum(
        "t,nt,tab->nab",
        0.5 * beta * weights,
        phase,
        _exact_tau(tau, beta, energies, residues),
    )
    expected = _exact_iw(reduced_n, beta, energies, residues)
    np.testing.assert_allclose(integrated, expected, rtol=0, atol=3e-14)
