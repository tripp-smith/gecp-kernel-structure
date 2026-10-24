"""Public physical facade: independent truth, reuse, and persistence boundaries."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import lstsq

from kernelgecp import FermionicBasis, FermionicFit, GECPConfig, NumericalBreakdownError
from kernelgecp._basis_io import read_archive, write_archive


@pytest.fixture(scope="module")
def basis():
    return FermionicBasis(beta=3.0, wmax=4.0, tol=1e-9, grid_order=6)


def green(tau, beta=3.0):
    eps = np.array([-2.1, 0.4, 3.1])
    weights = np.array([0.2, 0.5, 0.3])
    return (
        -(np.exp(-np.asarray(tau)[:, None] * eps) / (1 + np.exp(-beta * eps))) @ weights
    )


def green_iw(n, beta=3.0):
    eps = np.array([-2.1, 0.4, 3.1])
    weights = np.array([0.2, 0.5, 0.3])
    iw = 1j * np.pi * (2 * np.asarray(n) + 1) / beta
    return (1 / (iw[:, None] - eps)) @ weights


def test_independent_time_and_physical_frequency(basis):
    fit = basis.fit_tau(green(basis.tau_nodes))
    tau = (np.arange(301) + 0.314159) * basis.beta / 302
    n = np.array([-100, -10, -3, -1, 0, 2, 13, 99])
    assert np.max(np.abs(fit.evaluate_tau(tau) - green(tau))) < 2e-8
    assert np.max(np.abs(fit.evaluate_matsubara(n) - green_iw(n))) < 2e-8
    assert fit.diagnostics.training_residual_max < 1e-12
    assert fit.diagnostics.certification_status == "not_certified"
    assert basis.diagnostics.convergence_scope == "finite_construction_grid_only"


@pytest.mark.parametrize("axis", [0, 1, -1])
def test_batched_complex_matrix_layout_and_independent_qr(basis, axis):
    tau = np.linspace(0, basis.beta, 81)
    g = green(tau)
    weights = np.array([[1, 2j], [3 - 2j, -0.4]])
    samples = g[:, None, None] * weights
    samples = np.moveaxis(samples, 0, axis)
    fit = basis.fit_tau(samples, tau, axis=axis)
    design = basis.design_tau(tau)
    independent = lstsq(
        design, (g[:, None, None] * weights).reshape(81, 4), lapack_driver="gelsy"
    )[0]
    assert np.allclose(
        np.moveaxis(fit.coefficients, fit.axis, 0).reshape(basis.rank, 4),
        independent,
        atol=1e-12,
    )
    query = np.linspace(0.003, basis.beta - 0.003, 52)
    expected = np.moveaxis(green(query)[:, None, None] * weights, 0, axis)
    assert fit.evaluate_tau(query).shape == expected.shape
    assert np.max(np.abs(fit.evaluate_tau(query) - expected)) < 1e-7
    n = np.array([-5, -1, 0, 4])
    expected_iw = np.moveaxis(green_iw(n)[:, None, None] * weights, 0, axis)
    assert np.max(np.abs(fit.evaluate_matsubara(n) - expected_iw)) < 1e-7


def test_svd_is_cached_and_cache_owns_coordinates(basis, monkeypatch):
    basis.clear_cache()
    original = np.linalg.svd
    calls = []

    def observed(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "svd", observed)
    tau = np.linspace(0, basis.beta, 55)
    expected = tau.copy()
    first = basis.prepare_tau(tau)
    tau[:] = 0
    second = basis.prepare_tau(expected.copy())
    assert first is second
    basis.fit_tau(green(expected), expected)
    basis.fit_tau(2 * green(expected), expected)
    assert len(calls) == 1
    assert basis.cache_info["fit_plan_entries"] == 1
    with pytest.raises(ValueError):
        basis.prepare_tau(expected, rcond=True)


def test_bounded_cache_and_immutable_arrays():
    b = FermionicBasis(2.0, 2.0, grid_order=4, cache_size=2)
    for count in range(30, 35):
        b.prepare_tau(np.linspace(0, b.beta, count))
    assert b.cache_info == {
        "design_entries": 2,
        "fit_plan_entries": 2,
        "capacity_each": 2,
    }
    fit = b.fit_tau(-np.ones(b.rank))
    for array in [
        b.tau_nodes,
        b.omega_nodes,
        b.design_tau([0.1, 0.2]),
        fit.coefficients,
    ]:
        with pytest.raises(ValueError):
            array.setflags(write=True)
    b.clear_cache()
    assert b.cache_info["design_entries"] == b.cache_info["fit_plan_entries"] == 0


@pytest.mark.parametrize("precision", ["float64", "reference", "auto"])
def test_arithmetic_metadata_and_native_research_defaults(precision):
    defaults = GECPConfig()
    b = FermionicBasis(2, 3, tol=1e-7, precision=precision, grid_order=4)
    bits = 128 if precision == "reference" else 53
    assert b.diagnostics.selection_precision_bits == bits
    assert b.diagnostics.evaluation_precision_bits == bits
    assert b.diagnostics.export_precision_bits == b.diagnostics.fit_precision_bits == 53
    assert b.diagnostics.coordinate_precision_bits == 53
    assert GECPConfig() == defaults
    assert defaults.pivot == "adaptive" and defaults.precision_bits == 128


def test_auto_falls_back_only_for_arithmetic_breakdown(monkeypatch):
    original = FermionicBasis._select
    calls = []

    def observed(self, bits):
        calls.append(bits)
        if bits == 53:
            raise NumericalBreakdownError("controlled range failure")
        return original(self, bits)

    monkeypatch.setattr(FermionicBasis, "_select", observed)
    b = FermionicBasis(2, 2, grid_order=4)
    assert calls == [53, 128]
    assert b.diagnostics.attempts[0].outcome == "arithmetic_failure"
    assert b.diagnostics.evaluation_precision_bits == 128
    with pytest.raises(NumericalBreakdownError):
        FermionicBasis(2, 2, grid_order=4, precision="float64")

    def invalid(self, bits):
        raise ValueError("controlled input error")

    monkeypatch.setattr(FermionicBasis, "_select", invalid)
    with pytest.raises(ValueError, match="controlled input"):
        FermionicBasis(2, 2, grid_order=4)


def test_nonconvergence_is_not_silently_accepted():
    with pytest.raises(NumericalBreakdownError, match="did not converge"):
        FermionicBasis(3, 4, max_rank=1, grid_order=4)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"beta": True},
        {"beta": 0},
        {"beta": float("inf")},
        {"wmax": -1},
        {"beta": 0.1, "wmax": 0.1},
        {"beta": 1e300, "wmax": 1e300},
        {"tol": 0},
        {"tol": 1},
        {"tol": float("nan")},
        {"tol": "1e-8"},
        {"precision": "guess"},
        {"precision_bits": 53},
        {"grid_order": True},
        {"grid_order": 1},
        {"max_rank": 0},
        {"cache_size": 0},
        {"max_rank": 1.5},
    ],
)
def test_invalid_configuration(kwargs):
    params = {"beta": 2, "wmax": 3, **kwargs}
    with pytest.raises(ValueError):
        FermionicBasis(**params)


@pytest.mark.parametrize(
    "tau", [[-0.1], [3.01], [np.nan], [np.inf], [[0.1]], [True], [1j], "abc"]
)
def test_invalid_tau(basis, tau):
    with pytest.raises(ValueError):
        basis.design_tau(tau)


@pytest.mark.parametrize(
    "indices", [[0.5], [True], [[0]], [np.nan], [2**51], [-(2**51)], ["0"]]
)
def test_invalid_matsubara_indices(basis, indices):
    with pytest.raises(ValueError):
        basis.design_matsubara(indices)


def test_evaluation_shape_overflow_and_empty_queries(basis):
    with pytest.raises(ValueError, match="coefficient axis"):
        basis.evaluate_tau(1.0, [0.2])
    with pytest.raises(ValueError, match="basis rank"):
        basis.evaluate_tau(np.ones(basis.rank + 1), [0.2])
    with pytest.raises(ValueError, match="axis"):
        basis.evaluate_tau(np.ones(basis.rank), [0.2], axis=True)
    assert basis.evaluate_tau(np.ones(basis.rank), []).shape == (0,)
    assert basis.evaluate_tau(np.empty((basis.rank, 0)), [0.1]).shape == (1, 0)


@pytest.mark.parametrize("precision", ["float64", "reference"])
def test_safe_roundtrip_no_reselection(tmp_path, monkeypatch, precision):
    b = FermionicBasis(3, 4, grid_order=4, precision=precision)
    samples = green(b.tau_nodes)[:, None] * np.array([1 + 2j, -3j])
    fit = b.fit_tau(samples)
    basis_path = tmp_path / "basis.npz"
    fit_path = tmp_path / "fit.npz"
    b.save(basis_path)
    fit.save(fit_path)

    def forbidden(*args, **kwargs):
        pytest.fail("restoring must not rerun GECP")

    monkeypatch.setattr(FermionicBasis, "_select", forbidden)
    restored_basis = FermionicBasis.load(basis_path)
    restored_fit = FermionicFit.load(fit_path)
    tau = np.linspace(0, b.beta, 25)
    n = np.array([-4, -1, 0, 2])
    assert np.array_equal(restored_basis.tau_nodes, b.tau_nodes)
    assert np.array_equal(restored_fit.coefficients, fit.coefficients)
    assert np.array_equal(restored_fit.evaluate_tau(tau), fit.evaluate_tau(tau))
    assert np.array_equal(restored_fit.evaluate_matsubara(n), fit.evaluate_matsubara(n))
    assert restored_fit.diagnostics is None
    assert restored_basis.diagnostics.history_status == "restored_unverified_history"
    with pytest.raises(ValueError, match="standalone basis"):
        FermionicBasis.load(fit_path)
    with pytest.raises(ValueError, match="fitted function"):
        FermionicFit.load(basis_path)


@pytest.mark.parametrize(
    "change",
    ["version", "conventions", "nodes", "coefficient_shape", "precision", "history"],
)
def test_corrupt_semantic_archives_rejected(tmp_path, basis, change):
    original = tmp_path / "original.npz"
    basis.fit_tau(green(basis.tau_nodes)).save(original)
    metadata, arrays = read_archive(original)
    if change == "version":
        metadata["version"] = 999
    elif change == "conventions":
        metadata["conventions"]["fourier"] = "negative sign"
    elif change == "nodes":
        arrays["tau_nodes"] = np.zeros(basis.rank)
    elif change == "coefficient_shape":
        arrays["coefficients"] = np.ones(basis.rank + 1)
    elif change == "precision":
        metadata["construction"]["evaluation_precision_bits"] = 4097
    else:
        metadata["construction"]["converged"] = False
    corrupted = tmp_path / "corrupted.npz"
    write_archive(corrupted, metadata, arrays)
    with pytest.raises(ValueError):
        FermionicFit.load(corrupted)


def test_restored_history_cannot_claim_certification(tmp_path, basis):
    path = tmp_path / "basis.npz"
    basis.save(path)
    metadata, arrays = read_archive(path)
    metadata["construction"]["certification_status"] = "certified"
    write_archive(path, metadata, arrays)
    restored = FermionicBasis.load(path)
    assert restored.diagnostics.certification_status == "not_certified"
    assert replace(restored.diagnostics).converged


def test_large_finite_rcond_has_portable_diagnostics(tmp_path, basis):
    tau = np.linspace(0, basis.beta, 31)
    fit = basis.fit_tau(green(tau), tau, rcond=1e308)
    assert fit.diagnostics.least_squares.rank_threshold is None
    assert fit.diagnostics.least_squares.effective_rank == 0
    path = tmp_path / "zero_rank.npz"
    fit.save(path)
    restored = FermionicFit.load(path)
    assert np.array_equal(restored.coefficients, fit.coefficients)


def test_selection_isolates_harmless_underflow_without_weakening_other_traps(
    monkeypatch,
):
    import kernelgecp.basis as basis_module

    original = basis_module.gecp
    observed = []

    def checked(*args, **kwargs):
        observed.append(np.geterr())
        return original(*args, **kwargs)

    monkeypatch.setattr(basis_module, "gecp", checked)
    with np.errstate(all="raise"):
        FermionicBasis(2, 2, grid_order=4)
        assert np.geterr() == dict.fromkeys(
            ["divide", "over", "under", "invalid"], "raise"
        )
    assert observed == [
        {"divide": "raise", "over": "raise", "under": "ignore", "invalid": "raise"}
    ]
