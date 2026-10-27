"""Frequency fitting: independent physical/QR checks and archive boundaries."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import lstsq

from kernelgecp import FermionicBasis, FermionicFit
from kernelgecp._basis_io import read_archive, write_archive


@pytest.fixture(scope="module")
def basis():
    return FermionicBasis(2.0, 2.0, tol=1e-9, grid_order=6)


def _frequency_truth(indices, beta):
    energies = np.array([-0.7, 0.2, 1.1])
    masses = np.array([0.2, 0.5, 0.3])
    z = 1j * (2 * np.asarray(indices) + 1) * np.pi / beta
    return (1 / (z[:, None] - energies)) @ masses


def _time_truth(tau, beta):
    energies = np.array([-0.7, 0.2, 1.1])
    masses = np.array([0.2, 0.5, 0.3])
    return (
        -(np.exp(-np.asarray(tau)[:, None] * energies) / (1 + np.exp(-beta * energies)))
        @ masses
    )


@pytest.mark.parametrize("axis", [0, 1, 2, -1])
def test_signed_complex_matrix_physical_truth_and_independent_qr(basis, axis):
    indices = np.arange(-40, 40)
    amplitude = np.array([[1 + 0.5j, 0.2 - 0.3j], [-0.7j, 0.4 + 1j]])
    values0 = _frequency_truth(indices, basis.beta)[:, None, None] * amplitude
    values = np.moveaxis(values0, 0, axis)
    fit = basis.fit_matsubara(values, indices, axis=axis)
    independent = lstsq(
        basis.design_matsubara(indices),
        values0.reshape(indices.size, -1),
        lapack_driver="gelsy",
    )[0]
    np.testing.assert_allclose(
        np.moveaxis(fit.coefficients, fit.axis, 0).reshape(basis.rank, -1),
        independent,
        rtol=1e-12,
        atol=1e-12,
    )
    tau = np.concatenate(
        ([0.0, basis.beta], (np.arange(137) + 0.271828) * basis.beta / 138)
    )
    heldout = np.array([-100, -51, 44, 80])
    expected_tau = np.moveaxis(
        _time_truth(tau, basis.beta)[:, None, None] * amplitude, 0, axis
    )
    expected_iw = np.moveaxis(
        _frequency_truth(heldout, basis.beta)[:, None, None] * amplitude, 0, axis
    )
    assert fit.evaluate_tau(tau).shape == expected_tau.shape
    assert np.max(np.abs(fit.evaluate_tau(tau) - expected_tau)) < 1e-8
    assert np.max(np.abs(fit.evaluate_matsubara(heldout) - expected_iw)) < 1e-8
    assert fit.observation_domain == fit.diagnostics.observation_domain == "matsubara"
    assert fit.diagnostics.training_residual_max < 1e-9
    assert fit.diagnostics.certification_status == "not_certified"
    assert fit.diagnostics.least_squares.effective_rank == basis.rank


def test_positive_only_complex_minimum_norm_does_not_complete_conjugates(basis):
    indices = np.array([0, 2, 5])
    values = np.array([1 + 2j, -0.7 + 1.1j, 3 - 4j])
    fit = basis.fit_matsubara(values, indices)
    design = basis.design_matsubara(indices)
    independent = lstsq(design, values, lapack_driver="gelsy")[0]
    np.testing.assert_allclose(fit.coefficients, independent, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fit.evaluate_matsubara(indices), values, atol=1e-12)
    assert fit.diagnostics.least_squares.design_shape == (3, basis.rank)
    assert fit.diagnostics.least_squares.effective_rank == 3
    # This is an arbitrary complex signal. Symmetry was neither supplied nor
    # imposed, even though indices happen to be positive only.
    negative = fit.evaluate_matsubara(-indices - 1)
    assert np.max(np.abs(negative - values.conj())) > 0.1
    scaled = basis.fit_matsubara((0.3 - 0.8j) * values, indices)
    np.testing.assert_allclose(
        scaled.coefficients, (0.3 - 0.8j) * fit.coefficients, atol=1e-12
    )


def test_positive_only_training_success_is_not_cross_domain_accuracy(basis):
    indices = np.arange(40)
    fit = basis.fit_matsubara(_frequency_truth(indices, basis.beta), indices)
    assert fit.diagnostics.training_residual_max < 1e-8
    assert fit.diagnostics.least_squares.condition_number > 1e10
    # This retained analytic fixture demonstrates the observation/conditioning
    # boundary. It is not a promise that every positive-only fit is inaccurate.
    tau = (np.arange(137) + 0.271828) * basis.beta / 138
    error = np.max(np.abs(fit.evaluate_tau(tau) - _time_truth(tau, basis.beta)))
    assert error > 1e-7


def test_repeated_observations_have_direct_rank_one_minimum_norm_solution(basis):
    indices = np.array([3, 3, 3, 3])
    values = np.array([1 + 2j, 2 - 1j, -3 + 0.5j, 4 - 2j])
    fit = basis.fit_matsubara(values, indices, rcond=1e-12)
    row = basis.design_matsubara([3])[0]
    expected = row.conj() * np.mean(values) / np.sum(np.abs(row) ** 2)
    np.testing.assert_allclose(fit.coefficients, expected, atol=1e-13, rtol=1e-13)
    assert fit.diagnostics.least_squares.effective_rank == 1
    assert fit.diagnostics.least_squares.conditioning_status == "rank_truncated"
    assert fit.diagnostics.training_residual_max == pytest.approx(
        np.max(np.abs(values - np.mean(values)))
    )


def test_ill_conditioned_frequency_mesh_exposes_rank_truncation(basis):
    indices = 2**40 + np.arange(24, dtype=np.int64)
    values = _frequency_truth(indices, basis.beta)
    fit = basis.fit_matsubara(values, indices, rcond=1e-10)
    diagnostic = fit.diagnostics.least_squares
    assert diagnostic.effective_rank < basis.rank
    assert diagnostic.conditioning_status == "rank_truncated"
    assert np.all(np.isfinite(fit.coefficients))
    assert fit.diagnostics.certification_status == "not_certified"


def test_caches_reuse_svd_own_coordinates_and_separate_domains(basis, monkeypatch):
    basis.clear_cache()
    original = np.linalg.svd
    calls = []

    def observed(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "svd", observed)
    indices = np.array([0, 1, 2])
    expected = indices.copy()
    frequency = basis.prepare_matsubara(indices)
    indices[:] = 5
    assert basis.prepare_matsubara(expected.copy()) is frequency
    basis.fit_matsubara(np.ones(3, dtype=complex), expected)
    basis.fit_matsubara(2j * np.ones(3), expected)
    time = basis.prepare_tau(expected)
    assert time is not frequency
    assert np.iscomplexobj(frequency.design)
    assert not np.iscomplexobj(time.design)
    assert len(calls) == 2
    assert basis.cache_info["fit_plan_entries"] == 2
    truncated = basis.prepare_matsubara(expected, rcond=2)
    assert truncated is not frequency
    assert truncated.diagnostics.effective_rank == 0
    assert basis.prepare_matsubara(expected) is frequency
    assert len(calls) == 3
    for array in (frequency.design, time.design):
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_frequency_plan_lru_is_bounded_across_both_domains():
    b = FermionicBasis(2, 2, tol=1e-6, grid_order=4, cache_size=2)
    original = b.prepare_matsubara([0, 1])
    b.prepare_tau([0.1, 0.3])
    b.prepare_matsubara([2, 3])
    replacement = b.prepare_matsubara([0, 1])
    assert replacement is not original
    assert b.cache_info == {
        "design_entries": 2,
        "fit_plan_entries": 2,
        "capacity_each": 2,
    }
    b.clear_cache()
    assert b.cache_info["design_entries"] == b.cache_info["fit_plan_entries"] == 0


def test_frequency_fit_owns_data_preserves_empty_batches_and_order(basis):
    indices = np.array([5, -1, 0, -4, 2, -3])
    data = np.arange(24).reshape(2, 6, 2) * (1 + 2j)
    fit = basis.fit_matsubara(data, indices, axis=1)
    expected = fit.coefficients.copy()
    data[:] = 0
    indices[:] = 17
    np.testing.assert_array_equal(fit.coefficients, expected)
    with pytest.raises(ValueError):
        fit.coefficients.setflags(write=True)
    with pytest.raises(AttributeError):
        fit.observation_domain = "tau"
    empty = basis.fit_matsubara(np.empty((2, 0, 6, 3)), [5, -1, 0, -4, 2, -3], axis=2)
    assert empty.coefficients.shape == (2, 0, basis.rank, 3)
    assert empty.evaluate_tau([0.1, 0.2]).shape == (2, 0, 2, 3)
    assert empty.evaluate_matsubara([]).shape == (2, 0, 0, 3)
    assert empty.diagnostics.training_residual_max == 0


@pytest.mark.parametrize(
    "indices",
    [
        [],
        None,
        [True],
        [0, False],
        [[0]],
        [0.0],
        [1j],
        [np.nan],
        [2**51],
        [-(2**51)],
        [2**80],
        np.array([2**63], dtype=np.uint64),
        ["0"],
    ],
)
def test_invalid_frequency_observations_rejected(basis, indices):
    with pytest.raises(ValueError):
        basis.prepare_matsubara(indices)


@pytest.mark.parametrize(
    "rcond", [True, np.bool_(False), -1, np.nan, np.inf, 1j, "0", [0], 10**1000]
)
def test_invalid_rank_cutoff_rejected_even_for_cache_hits(basis, rcond):
    basis.prepare_matsubara([0, 1], rcond=0)
    with pytest.raises(ValueError):
        basis.prepare_matsubara([0, 1], rcond=rcond)


@pytest.mark.parametrize(
    "values,axis",
    [
        (1, 0),
        ([1], 0),
        ([True, False], 0),
        ([1, np.nan], 0),
        ([1j, np.inf], 0),
        ([1, 2], True),
        ([1, 2], 1),
    ],
)
def test_invalid_values_and_axes_rejected(basis, values, axis):
    with pytest.raises(ValueError):
        basis.fit_matsubara(values, [0, 1], axis=axis)


@pytest.mark.parametrize("domain", [None, "frequency", "", True, 2, [], {}])
def test_manual_fit_domain_is_explicit_and_validated(basis, domain):
    with pytest.raises(ValueError, match="observation_domain"):
        FermionicFit(basis, np.ones(basis.rank), observation_domain=domain)


def test_manual_default_time_domain_and_diagnostic_consistency(basis):
    manual = FermionicFit(basis, np.ones(basis.rank))
    assert manual.observation_domain == "tau"
    time_fit = basis.fit_tau(-np.ones(basis.rank))
    assert time_fit.diagnostics.observation_domain == "tau"
    with pytest.raises(ValueError, match="domains must agree"):
        FermionicFit(
            basis,
            np.ones(basis.rank),
            diagnostics=time_fit.diagnostics,
            observation_domain="matsubara",
        )
    with pytest.raises(ValueError, match="observation_domain"):
        FermionicFit(
            basis,
            np.ones(basis.rank),
            diagnostics=replace(time_fit.diagnostics, observation_domain="invalid"),
        )


@pytest.mark.parametrize("domain", ["tau", "matsubara"])
def test_version_two_fit_and_basis_roundtrip(tmp_path, basis, monkeypatch, domain):
    indices = np.arange(-16, 16)
    amplitude = np.array([[1], [2j]])
    if domain == "matsubara":
        data = _frequency_truth(indices, basis.beta)[None, :] * amplitude
        fit = basis.fit_matsubara(data, indices, axis=-1)
    else:
        times = np.linspace(0, basis.beta, indices.size)
        data = _time_truth(times, basis.beta)[None, :] * amplitude
        fit = basis.fit_tau(data, times, axis=-1)
    path = tmp_path / "frequency.npz"
    basis_path = tmp_path / "basis.npz"
    fit.save(path)
    basis.save(basis_path)
    metadata, _ = read_archive(path)
    assert metadata["version"] == 2
    assert metadata["fit"]["observation_domain"] == domain
    assert metadata["fit"]["diagnostics"]["observation_domain"] == domain
    assert read_archive(basis_path)[0]["version"] == 2

    def forbidden(*args, **kwargs):
        pytest.fail("archive load must not reselect pivots or fit an SVD")

    monkeypatch.setattr(FermionicBasis, "_select", forbidden)
    monkeypatch.setattr(FermionicBasis, "prepare_tau", forbidden)
    monkeypatch.setattr(FermionicBasis, "prepare_matsubara", forbidden)
    restored = FermionicFit.load(path)
    restored_basis = FermionicBasis.load(basis_path)
    assert restored.observation_domain == domain
    assert restored.diagnostics is None
    assert restored.axis == fit.axis
    np.testing.assert_array_equal(restored.coefficients, fit.coefficients)
    np.testing.assert_array_equal(restored_basis.tau_nodes, basis.tau_nodes)
    for query in ([0, 1], [-4, 3]):
        np.testing.assert_array_equal(
            restored.evaluate_matsubara(query), fit.evaluate_matsubara(query)
        )


@pytest.mark.parametrize("domain", [None, "iw", "", True, 2, [], {}, "missing"])
def test_version_two_archive_domain_rejected_before_core_restore(
    tmp_path, basis, monkeypatch, domain
):
    import kernelgecp.basis as basis_module

    path = tmp_path / "bad_domain.npz"
    basis.fit_matsubara([1 + 2j, 3j], [0, 1]).save(path)
    metadata, arrays = read_archive(path)
    if domain == "missing":
        del metadata["fit"]["observation_domain"]
    else:
        metadata["fit"]["observation_domain"] = domain
    write_archive(path, metadata, arrays)

    def forbidden(*args, **kwargs):
        pytest.fail("invalid observation domain reached core allocation")

    monkeypatch.setattr(basis_module, "BasisNumerics", forbidden)
    with pytest.raises(ValueError):
        FermionicFit.load(path)
