"""Independent finite-model physics and arithmetic checks for ED references."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import expit

from kernelgecp.lehmann import (
    FermionicPoleModel,
    finite_fock_lehmann,
    fock_annihilation_operators,
    hubbard_atom,
    impurity_one_bath,
    ipt_self_energy_poles,
    thermal_trace_tau,
)
from kernelgecp.types import NumericalBreakdownError


@pytest.mark.parametrize("modes", [1, 2, 3, 4])
def test_all_fermionic_canonical_anticommutation_relations(modes):
    operators = fock_annihilation_operators(modes)
    identity = np.eye(2**modes)
    for i, left in enumerate(operators):
        for j, right in enumerate(operators):
            np.testing.assert_array_equal(left @ right + right @ left, 0 * identity)
            np.testing.assert_array_equal(
                left @ right.T + right.T @ left,
                identity if i == j else 0 * identity,
            )
        with pytest.raises(ValueError):
            left.setflags(write=True)


@pytest.mark.parametrize("beta", [0.2, 10.0, 1000.0])
@pytest.mark.parametrize("interaction", [0.0, 1.0, 2.0])
def test_fock_atom_matches_independent_two_pole_resolvent(beta, interaction):
    atom = hubbard_atom(beta, interaction)
    indices = np.array([-29, -3, -1, 0, 2, 31])
    z = 1j * (2 * indices + 1) * np.pi / beta
    expected = 0.5 / (z - interaction / 2) + 0.5 / (z + interaction / 2)
    np.testing.assert_allclose(atom.evaluate_matsubara(indices), expected, atol=2e-13)
    assert atom.mass == pytest.approx(1, abs=2e-14)
    assert atom.density() == pytest.approx(0.5, abs=2e-14)
    np.testing.assert_allclose(atom.moments(), [1, 0, interaction**2 / 4], atol=2e-14)
    times = np.array([0, beta * 0.13, beta * 0.5, beta])
    exact = FermionicPoleModel(beta, [-interaction / 2, interaction / 2], [0.5, 0.5])
    np.testing.assert_allclose(
        atom.evaluate_tau(times), exact.evaluate_tau(times), atol=2e-14
    )


def _independent_impurity_hamiltonian(interaction, hybridization):
    operators = fock_annihilation_operators(4)
    numbers = [c.T @ c for c in operators]
    h = interaction * (numbers[0] - 0.5 * np.eye(16)) @ (numbers[1] - 0.5 * np.eye(16))
    for impurity, bath in ((0, 2), (1, 3)):
        h = h + hybridization * (
            operators[impurity].T @ operators[bath]
            + operators[bath].T @ operators[impurity]
        )
    return h, operators[0]


@pytest.mark.parametrize("beta", [0.4, 3.0, 50.0])
def test_interacting_finite_bath_against_independent_matrix_exponentials(beta):
    h, c = _independent_impurity_hamiltonian(2.0, 0.7)
    poles = impurity_one_bath(beta, 2.0, 0.7)
    times = np.linspace(0, beta, 11)
    trace = thermal_trace_tau(beta, h, c, times)
    np.testing.assert_allclose(poles.evaluate_tau(times), trace, atol=3e-13, rtol=2e-12)
    assert poles.mass == pytest.approx(1, abs=2e-14)
    assert poles.density() == pytest.approx(0.5, abs=2e-14)
    np.testing.assert_allclose(poles.moments(), [1, 0, 1 + 0.7**2], atol=3e-13)
    np.testing.assert_allclose(
        poles.evaluate_tau(times), poles.evaluate_tau(beta - times), atol=2e-14
    )
    np.testing.assert_allclose(
        poles.evaluate_matsubara([-3, -1]),
        poles.evaluate_matsubara([2, 0]).conj(),
        atol=2e-14,
    )
    positive = poles.evaluate_matsubara(np.arange(50))
    assert np.all(positive.imag < 0)


def test_complex_rotated_fock_basis_preserves_scalar_thermal_reference():
    first, second = fock_annihilation_operators(2)
    h = 0.6 * first.T @ first - 0.9 * second.T @ second
    rng = np.random.default_rng(4837)
    q, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
    rotated_h = q @ h @ q.conj().T
    # Explicitly construct the Hermitian matrix; no approximate projection
    # occurs inside the reference helper's Hamiltonian validation.
    rotated_h = np.triu(rotated_h) + np.triu(rotated_h, 1).conj().T
    rotated_h[np.diag_indices(4)] = rotated_h.diagonal().real
    c = q @ first @ q.conj().T
    model = finite_fock_lehmann(2.3, rotated_h, c)
    expected = FermionicPoleModel(2.3, [0.6], [1.0])
    tau = np.linspace(0, 2.3, 13)
    np.testing.assert_allclose(
        model.evaluate_tau(tau), expected.evaluate_tau(tau), atol=2e-14
    )
    np.testing.assert_allclose(
        thermal_trace_tau(2.3, rotated_h, c, tau),
        expected.evaluate_tau(tau),
        atol=2e-14,
    )


def test_noninteracting_bath_and_decoupled_atom_limits_are_independent_oracles():
    beta, v = 10.0, 0.7
    indices = np.arange(-8, 8)
    z = 1j * (2 * indices + 1) * np.pi / beta
    np.testing.assert_allclose(
        impurity_one_bath(beta, 0, v).evaluate_matsubara(indices),
        z / (z**2 - v**2),
        atol=2e-14,
    )
    np.testing.assert_allclose(
        impurity_one_bath(beta, 2, 0).evaluate_matsubara(indices),
        hubbard_atom(beta, 2).evaluate_matsubara(indices),
        atol=2e-14,
    )
    np.testing.assert_allclose(
        impurity_one_bath(beta, 2, -v).evaluate_matsubara(indices),
        impurity_one_bath(beta, 2, v).evaluate_matsubara(indices),
        atol=2e-14,
    )


def test_ipt_triples_match_reflected_product_and_reject_the_asymmetric_cube():
    weiss = FermionicPoleModel(10, [-1.3, 0.2, 2.0], [0.2, 0.5, 0.3])
    sigma = ipt_self_energy_poles(weiss, 2)
    tau = np.linspace(0, weiss.beta, 31)
    g = weiss.evaluate_tau(tau)
    reflected = 4 * g**2 * weiss.evaluate_tau(weiss.beta - tau)
    np.testing.assert_allclose(sigma.evaluate_tau(tau), reflected, atol=2e-15)
    assert np.max(np.abs(sigma.evaluate_tau(tau) - 4 * g**3)) > 0.1
    n0 = weiss.density()
    assert sigma.mass == pytest.approx(4 * n0 * (1 - n0), rel=2e-15)
    assert np.all(sigma.weights >= 0)
    assert sigma.spectral_cutoff <= 3 * weiss.spectral_cutoff
    assert sigma.mass != pytest.approx(1)
    # Direct quadrature provides a separate frequency check of the reflection.
    nodes, weights = np.polynomial.legendre.leggauss(160)
    points = (nodes + 1) * weiss.beta / 2
    values = (
        4 * weiss.evaluate_tau(points) ** 2 * weiss.evaluate_tau(weiss.beta - points)
    )
    indices = np.array([-4, -1, 0, 3])
    phase = np.exp(1j * (2 * indices[:, None] + 1) * np.pi * points / weiss.beta)
    integrated = phase @ (weights * values * weiss.beta / 2)
    np.testing.assert_allclose(
        sigma.evaluate_matsubara(indices), integrated, atol=2e-14
    )


def test_atomic_ipt_limit_is_exact_but_interacting_finite_bath_is_a_different_model():
    beta, u = 10.0, 2.0
    sigma = ipt_self_energy_poles(FermionicPoleModel(beta, [0], [1]), u)
    indices = np.arange(-12, 12)
    z = 1j * (2 * indices + 1) * np.pi / beta
    np.testing.assert_allclose(sigma.evaluate_matsubara(indices), u**2 / (4 * z))
    np.testing.assert_allclose(
        1 / (z - sigma.evaluate_matsubara(indices)),
        hubbard_atom(beta, u).evaluate_matsubara(indices),
        atol=2e-14,
    )
    weiss = FermionicPoleModel(beta, [-0.7, 0.7], [0.5, 0.5])
    finite_sigma = ipt_self_energy_poles(weiss, u).evaluate_matsubara(indices)
    approximate = 1 / (1 / weiss.evaluate_matsubara(indices) - finite_sigma)
    exact = impurity_one_bath(beta, u, 0.7).evaluate_matsubara(indices)
    assert np.max(np.abs(approximate - exact)) > 1e-4


def test_stable_thermal_extremes_and_nonunit_zero_mass_semantics():
    poles = FermionicPoleModel(100, [-1000, 0, 1000], [0.2, 0.3, 0.5])
    with np.errstate(all="raise"):
        values = poles.evaluate_tau([0, 1e-12, 50, 100 - 1e-12, 100])
    assert np.all(np.isfinite(values))
    np.testing.assert_allclose(values[[0, -1]], [-0.65, -0.35], atol=1e-15)
    assert poles.density() == pytest.approx(0.35)
    assert poles.evaluate_tau([]).shape == (0,)
    assert poles.evaluate_matsubara([]).shape == (0,)
    nonunit = FermionicPoleModel(2, [-1, 1], [2, 3])
    assert nonunit.mass == 5
    assert np.sum(nonunit.evaluate_tau([0, 2])) == pytest.approx(-5)
    zero = ipt_self_energy_poles(FermionicPoleModel(2, [0], [1]), 0)
    np.testing.assert_array_equal(zero.evaluate_tau([0, 1, 2]), [0, 0, 0])
    np.testing.assert_array_equal(zero.moments(), [0, 0, 0])
    np.testing.assert_array_equal(zero.evaluate_matsubara([-1, 0]), [0j, 0j])
    # Tiny holes must not disappear through subtraction from a rounded one.
    tiny = ipt_self_energy_poles(FermionicPoleModel(1, [-40], [1]), 2)
    assert tiny.mass == pytest.approx(4 * expit(-40) * expit(40), rel=2e-15)


def test_model_owns_immutable_sorted_arrays_and_records_explicit_cutoff():
    energies = np.array([1.0, -2.0])
    weights = np.array([0.3, 0.7])
    model = FermionicPoleModel(2, energies, weights, cutoff=3)
    energies[:] = 0
    weights[:] = 0
    np.testing.assert_array_equal(model.energies, [-2, 1])
    np.testing.assert_array_equal(model.weights, [0.7, 0.3])
    assert model.declared_cutoff == 3
    assert model.spectral_cutoff == 2
    for array in (
        model.energies,
        model.weights,
        model.moments(),
        model.evaluate_tau([0]),
    ):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(AttributeError):
        model.beta = 5


@pytest.mark.parametrize("beta", [True, 0, -1, np.nan, np.inf, 10**1000])
def test_invalid_temperature_is_rejected(beta):
    with pytest.raises(ValueError):
        FermionicPoleModel(beta, [0], [1])


@pytest.mark.parametrize(
    "energies,weights",
    [
        ([], []),
        ([0], [0.5, 0.5]),
        ([True], [1]),
        ([1j], [1]),
        ([np.inf], [1]),
        ([0], [-1]),
        ([0], [True]),
        ([0], [np.nan]),
        ([[0]], [1]),
        ([0, 1], [1e308, 1e308]),
    ],
)
def test_invalid_pole_data_is_rejected(energies, weights):
    with pytest.raises(ValueError):
        FermionicPoleModel(2, energies, weights)


@pytest.mark.parametrize("cutoff", [True, -1, np.inf, np.nan, 0.5])
def test_invalid_or_exceeded_explicit_support_bound_is_rejected(cutoff):
    with pytest.raises(ValueError):
        FermionicPoleModel(2, [-1, 1], [0.5, 0.5], cutoff=cutoff)


@pytest.mark.parametrize("indices", [[True], [0.0], [1j], [[0]], [2**51], [2**80]])
def test_invalid_frequency_indices_are_not_coerced(indices):
    with pytest.raises(ValueError):
        FermionicPoleModel(2, [0], [1]).evaluate_matsubara(indices)


@pytest.mark.parametrize("tau", [[True], [-0.1], [2.1], [np.nan], [1j], [[0]]])
def test_invalid_time_coordinates_are_rejected(tau):
    with pytest.raises(ValueError):
        FermionicPoleModel(2, [0], [1]).evaluate_tau(tau)


@pytest.mark.parametrize("order", [True, -1, 2.0, 17])
def test_invalid_moment_order_is_rejected(order):
    with pytest.raises(ValueError):
        FermionicPoleModel(2, [0], [1]).moments(order)


def test_unrepresentable_frequency_and_moment_arithmetic_fail_closed():
    with pytest.raises(NumericalBreakdownError):
        FermionicPoleModel(1e-308, [0], [1]).evaluate_matsubara([0])
    with pytest.raises(NumericalBreakdownError):
        FermionicPoleModel(2, [1e200], [1]).moments(2)
    np.testing.assert_array_equal(
        FermionicPoleModel(2, [1e200], [0]).moments(), [0, 0, 0]
    )


def test_unrepresentable_fock_energy_shift_fails_without_overflow_warning():
    c = fock_annihilation_operators(1)[0]
    h = np.diag([-1e308, 1e308])
    with np.errstate(all="raise"):
        with pytest.raises(NumericalBreakdownError):
            finite_fock_lehmann(1, h, c)
        with pytest.raises(NumericalBreakdownError, match="energy shift"):
            thermal_trace_tau(1, h, c, [0, 0.5, 1])


def test_invalid_fock_and_ipt_inputs_fail_before_physical_claims():
    c = fock_annihilation_operators(1)[0]
    for h, operator in (
        (np.ones((3, 3)), np.eye(3)),
        (np.array([[0, 1], [0, 0]]), c),
        (np.eye(2), 2 * c),
        (np.eye(2), np.eye(2)),
    ):
        with pytest.raises(ValueError):
            finite_fock_lehmann(2, h, operator)
    for count in (True, 0, -1, 5, 2.0):
        with pytest.raises(ValueError):
            fock_annihilation_operators(count)
    for interaction in (True, -1, np.nan, np.inf):
        with pytest.raises(ValueError):
            hubbard_atom(2, interaction)
        with pytest.raises(ValueError):
            ipt_self_energy_poles(FermionicPoleModel(2, [0], [1]), interaction)
    with pytest.raises(ValueError, match="unit-mass"):
        ipt_self_energy_poles(FermionicPoleModel(2, [0], [2]), 1)
