"""Adversarial checks of facade ownership and untrusted historical metadata."""

from __future__ import annotations

import copy

import numpy as np
import pytest

from kernelgecp import FermionicBasis, FermionicFit
from kernelgecp._basis_io import read_archive, write_archive


@pytest.fixture(scope="module")
def small_basis():
    return FermionicBasis(2.0, 2.0, tol=1e-6, grid_order=4, precision="float64")


@pytest.fixture
def fit_archive(tmp_path, small_basis):
    path = tmp_path / "fit.npz"
    coefficients = np.arange(small_basis.rank)[:, None] * np.array([1 + 2j, -3j])
    FermionicFit(small_basis, coefficients).save(path)
    metadata, arrays = read_archive(path)
    return path, metadata, arrays


def test_fit_owns_complex_batches_and_cache_eviction_does_not_change_them(small_basis):
    generator = np.random.default_rng(721)
    coefficients = generator.normal(
        size=(2, 3, small_basis.rank)
    ) + 1j * generator.normal(size=(2, 3, small_basis.rank))
    original = coefficients.copy()
    fit = FermionicFit(small_basis, coefficients, axis=-1)
    coefficients[:] = 0
    tau = np.array([0.03, 0.22, 1.73])
    expected = np.moveaxis(
        small_basis.design_tau(tau)
        @ np.moveaxis(original, -1, 0).reshape(small_basis.rank, -1),
        0,
        -1,
    ).reshape(2, 3, tau.size)
    np.testing.assert_allclose(fit.evaluate_tau(tau), expected, rtol=1e-13, atol=1e-13)
    small_basis.clear_cache()
    np.testing.assert_allclose(fit.evaluate_tau(tau), expected, rtol=1e-13, atol=1e-13)
    with pytest.raises(ValueError):
        fit.coefficients.setflags(write=True)
    assert fit.evaluate_matsubara([]).shape == (2, 3, 0)


def test_zero_batch_and_negative_observation_axis_roundtrip(tmp_path, small_basis):
    tau = np.linspace(0, small_basis.beta, 19)
    samples = np.empty((2, 0, tau.size), dtype=np.complex128)
    fit = small_basis.fit_tau(samples, tau, axis=-1)
    assert fit.coefficients.shape == (2, 0, small_basis.rank)
    assert fit.diagnostics.training_residual_max == 0
    path = tmp_path / "empty_fit.npz"
    fit.save(path)
    restored = FermionicFit.load(path)
    assert restored.axis == 2
    assert restored.coefficients.shape == fit.coefficients.shape
    assert restored.coefficients.dtype == np.complex128
    assert restored.evaluate_tau([0.2, 0.3]).shape == (2, 0, 2)


def test_fit_plan_cache_separates_rank_thresholds_and_rejects_bool_hits(small_basis):
    tau = np.linspace(0, small_basis.beta, 15)
    full = small_basis.prepare_tau(tau, rcond=0)
    truncated = small_basis.prepare_tau(tau.copy(), rcond=1)
    assert full is not truncated
    assert full.diagnostics.effective_rank > 0
    assert truncated.diagnostics.effective_rank == 0
    assert small_basis.prepare_tau(tau, rcond=0) is full
    for invalid in (True, np.bool_(False), 10**1000):
        with pytest.raises(ValueError):
            small_basis.prepare_tau(tau, rcond=invalid)


@pytest.mark.parametrize("precision", [np.array("auto"), np.array(["auto"]), True])
def test_precision_policy_requires_a_string(precision):
    with pytest.raises(ValueError):
        FermionicBasis(2, 2, grid_order=4, precision=precision)


@pytest.mark.parametrize(
    "field,value",
    [
        ("rank", 0),
        ("requested_policy", "reference"),
        ("reference_precision_bits", 129),
        ("selection_precision_bits", 128),
        ("coordinate_precision_bits", 128),
        ("export_precision_bits", 128),
        ("fit_precision_bits", 128),
    ],
)
def test_restore_rejects_contradictory_precision_and_rank_history(
    fit_archive, field, value
):
    path, metadata, arrays = fit_archive
    metadata["construction"][field] = value
    write_archive(path, metadata, arrays)
    with pytest.raises(ValueError):
        FermionicFit.load(path)


@pytest.mark.parametrize("corruption", ["precision", "outcome", "count"])
def test_restore_rejects_incoherent_attempt_history(fit_archive, corruption):
    path, metadata, arrays = fit_archive
    attempts = metadata["construction"]["attempts"]
    if corruption == "precision":
        attempts[-1]["precision_bits"] = 4097
    elif corruption == "outcome":
        attempts[-1]["outcome"] = "contradicts_final_stop_reason"
    else:
        attempts.insert(
            0, {"precision_bits": 53, "outcome": "arithmetic_failure", "detail": ""}
        )
    write_archive(path, metadata, arrays)
    with pytest.raises(ValueError):
        FermionicFit.load(path)


@pytest.mark.parametrize(
    "corruption", ["extra_metadata", "extra_fit", "version_bool", "kind", "extra_array"]
)
def test_restore_rejects_unknown_or_mismatched_schema_payloads(fit_archive, corruption):
    path, metadata, arrays = fit_archive
    if corruption == "extra_metadata":
        metadata["unexpected"] = 123
    elif corruption == "extra_fit":
        metadata["fit"]["unexpected"] = 123
    elif corruption == "version_bool":
        metadata["version"] = True
    elif corruption == "kind":
        metadata["kind"] = "basis"
    else:
        arrays["unexpected"] = np.array([123])
    write_archive(path, metadata, arrays)
    with pytest.raises(ValueError):
        FermionicFit.load(path)


def test_restored_history_cannot_expand_current_evidence_scope(fit_archive):
    path, metadata, arrays = fit_archive
    history = metadata["construction"]
    history["certification_status"] = "certified"
    history["certification_scope"] = "continuous_error_proof"
    history["convergence_scope"] = "all_continuous_inputs"
    history["history_status"] = "independently_verified"
    metadata["fit"]["diagnostics"] = {
        "training_residual_max": 0,
        "certification_status": "certified",
    }
    write_archive(path, metadata, arrays)
    restored = FermionicFit.load(path)
    diagnostic = restored.basis.diagnostics
    assert diagnostic.history_status == "restored_unverified_history"
    assert diagnostic.convergence_scope == "finite_construction_grid_only"
    assert diagnostic.certification_status == "not_certified"
    assert (
        diagnostic.certification_scope == "no_continuous_input_or_rounding_certificate"
    )
    assert restored.diagnostics is None


@pytest.mark.parametrize("count,bits", [(257, 53), (65, 4096)])
def test_restore_limits_core_work_before_numerical_allocation(
    fit_archive, monkeypatch, count, bits
):
    import kernelgecp.basis as basis_module

    path, metadata, arrays = fit_archive
    metadata = copy.deepcopy(metadata)
    parameters = metadata["parameters"]
    parameters["max_rank"] = count
    parameters["precision"] = "float64" if bits == 53 else "reference"
    parameters["precision_bits"] = 128 if bits == 53 else bits
    history = metadata["construction"]
    history.update(
        rank=count,
        requested_policy=parameters["precision"],
        reference_precision_bits=parameters["precision_bits"],
        selection_precision_bits=bits,
        evaluation_precision_bits=bits,
        tie_counts=[0] * count,
        attempts=[
            {"precision_bits": bits, "outcome": history["stop_reason"], "detail": ""}
        ],
    )
    arrays["tau_nodes"] = np.linspace(0, parameters["beta"], count)
    cutoff = parameters["beta"] * parameters["wmax"]
    arrays["omega_nodes"] = np.linspace(-cutoff, cutoff, count)
    arrays["coefficients"] = np.zeros(count)
    write_archive(path, metadata, arrays)

    def forbidden_numerics(*args, **kwargs):
        pytest.fail("unsafe restored core reached numerical allocation")

    monkeypatch.setattr(basis_module, "BasisNumerics", forbidden_numerics)
    with pytest.raises(ValueError):
        FermionicFit.load(path)
