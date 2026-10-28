"""Independent resolvent, expm and actual PySCF-sector molecular checks."""

import importlib.util

import numpy as np
import pytest
from scipy.linalg import expm

from experiments.molecular_reference import (
    MODE_ORDER,
    build_h2_reference,
    molecular_hamiltonian,
    reference_from_integrals,
)
from kernelgecp.lehmann import fock_annihilation_operators
from kernelgecp.types import NumericalBreakdownError


def test_noninteracting_matrix_resolvent_and_restricted_axes():
    h = np.array([[-0.7, 0.2], [0.2, 0.3]])
    model = reference_from_integrals(h, np.zeros((2, 2, 2, 2)), beta=7, mu=-0.1)
    labels = np.array([-25, -3, -1, 0, 2, 24])
    z = 1j * (2 * labels + 1) * np.pi / 7
    expected = np.stack(
        [np.linalg.solve((w - 0.1) * np.eye(2) - h, np.eye(2)) for w in z]
    )
    actual = model.evaluate_matsubara(labels)
    assert actual.shape == (len(labels), 1, 1, 2, 2)
    np.testing.assert_allclose(actual[:, 0, 0], expected, atol=3e-15)
    endpoints = model.evaluate_tau([0, 7])
    np.testing.assert_allclose(endpoints.sum(axis=0)[0, 0], -np.eye(2), atol=3e-15)
    assert MODE_ORDER == (
        "orbital0_up",
        "orbital0_down",
        "orbital1_up",
        "orbital1_down",
    )


def test_interacting_hamiltonian_and_thermal_trace_are_independent_oracles():
    u, beta, mu = 1.4, 3.2, 0.13
    h = np.array([[-0.2, 0.17], [0.17, 0.4]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = eri[1, 1, 1, 1] = u
    model = reference_from_integrals(h, eri, beta=beta, mu=mu)
    operators = fock_annihilation_operators(4)
    number = [c.T @ c for c in operators]
    expected_h = sum(
        h[p, q] * operators[2 * p + s].T @ operators[2 * q + s]
        for p, q, s in np.ndindex(2, 2, 2)
    )
    expected_h += u * (number[0] @ number[1] + number[2] @ number[3])
    np.testing.assert_allclose(model.hamiltonian, expected_h, atol=1e-16)
    grand = expected_h - mu * sum(number)
    # Shift only the independent expm calculation for bounded exponentials.
    grand -= np.eye(16) * np.linalg.eigvalsh(grand)[0]
    partition = np.trace(expm(-beta * grand))
    times = np.linspace(0, beta, 7)
    observed = model.evaluate_spin_tau(times)
    for row, tau in enumerate(times):
        for p, q in np.ndindex(4, 4):
            expected = (
                -np.trace(
                    expm(-(beta - tau) * grand)
                    @ operators[p]
                    @ expm(-tau * grand)
                    @ operators[q].T
                )
                / partition
            )
            assert observed[row, p, q] == pytest.approx(expected, abs=2e-14)
    density = -np.trace(observed[-1]).real
    assert density == pytest.approx(model.electron_count, abs=2e-14)
    assert abs(model.electron_count - 2) > 0.01  # supplied mu was not neutralized
    neutral = reference_from_integrals(h, eri, beta=beta)
    assert neutral.electron_count == pytest.approx(2, abs=2e-13)


def test_matrix_fourier_quadrature_has_independent_sign_and_normalization():
    h = np.array([[-0.5, 0.2], [0.2, 0.1]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = 0.8
    model = reference_from_integrals(h, eri, beta=6)
    nodes, weights = np.polynomial.legendre.leggauss(100)
    times = model.beta * (nodes + 1) / 2
    labels = np.array([-5, -1, 0, 4])
    phases = np.exp(1j * (2 * labels[:, None] + 1) * np.pi * times / model.beta)
    integrated = np.einsum(
        "nt,tpq,t->npq",
        phases,
        model.evaluate_spin_tau(times),
        weights * model.beta / 2,
    )
    np.testing.assert_allclose(
        model.evaluate_spin_matsubara(labels), integrated, atol=2e-14
    )
    for array in (
        model.hamiltonian,
        model.lowdin_x,
        model.spectral_weights,
        integrated[:0],
    ):
        if array.size:
            with pytest.raises(ValueError):
                array.setflags(write=True)


@pytest.mark.parametrize("beta", [True, 0, -1, np.nan, np.inf, 10**1000])
def test_bad_beta_fails_before_optional_import(beta):
    with pytest.raises(ValueError):
        build_h2_reference(0.74, beta=beta)


@pytest.mark.parametrize(
    "field,value",
    [
        ("mu", True),
        ("mu", np.nan),
        ("mu", 10**1000),
        ("nuclear_energy", True),
        ("nuclear_energy", np.inf),
    ],
)
def test_bad_thermodynamic_parameters_are_rejected(field, value):
    with pytest.raises(ValueError):
        reference_from_integrals(np.eye(2), np.zeros((2, 2, 2, 2)), **{field: value})


def test_invalid_integrals_and_resource_arithmetic_fail_clearly():
    eri = np.zeros((2, 2, 2, 2))
    for h in ([[0, True], [1, 0]], [[0, 1], [0, 0]], [[0, np.inf], [0, 0]]):
        with pytest.raises(ValueError):
            molecular_hamiltonian(h, eri)
    bad = eri.copy()
    bad[0, 1, 0, 0] = 0.1
    with pytest.raises(ValueError, match="permutation"):
        molecular_hamiltonian(np.eye(2), bad)
    with pytest.raises(NumericalBreakdownError):
        reference_from_integrals(np.eye(2), eri, mu=1e308)


@pytest.mark.skipif(
    importlib.util.find_spec("pyscf") is None,
    reason="actual optional PySCF runtime required",
)
@pytest.mark.parametrize("bond", [0.74, 1.5])
def test_actual_pyscf_all_fci_sectors_and_lowdin_model(bond):
    model = build_h2_reference(bond)
    assert len(model.sector_energy_errors) == 9
    assert max(error for _, _, error in model.sector_energy_errors) < 3e-10
    np.testing.assert_allclose(
        model.lowdin_x.T @ model.overlap_ao @ model.lowdin_x, np.eye(2), atol=3e-15
    )
    assert model.electron_count == pytest.approx(2, abs=2e-12)
    assert model.exact_hf_energy_hartree < -0.9
    assert model.thermal_energy_hartree < model.exact_hf_energy_hartree
    returned = model.evaluate_tau([0, 10, 20])
    np.testing.assert_allclose(
        returned[[0, -1]].sum(axis=0)[0, 0], -np.eye(2), atol=3e-14
    )
    supplied = build_h2_reference(bond, mu=0.0462183952331543)
    assert supplied.mu == 0.0462183952331543
