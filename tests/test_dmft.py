"""Independent pole/QR oracles and controller boundaries for projected IPT."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import lstsq

from kernelgecp import FermionicBasis, NumericalBreakdownError
from kernelgecp.dmft import IPTConfig, IPTIteration, IPTPlan
from kernelgecp.lehmann import (
    FermionicPoleModel,
    hubbard_atom,
    ipt_self_energy_poles,
)


class PoleBasis:
    """Physical spectral coordinates, independent of GECP cardinal bases."""

    def __init__(self, energies, *, beta=3.0, wmax=None):
        self.energies = np.asarray(energies, dtype=float)
        self.beta = beta
        self.wmax = (
            max(1.0, float(np.max(np.abs(self.energies)))) if wmax is None else wmax
        )
        self.rank = len(self.energies)
        self.designs = []

    def design_tau(self, tau):
        times = np.asarray(tau)
        # All test poles and beta are small: this direct formula intentionally
        # differs from the stable occupied/empty oracle in lehmann.py.
        result = -np.exp(-times[:, None] * self.energies) / (
            1 + np.exp(-self.beta * self.energies)
        )
        self.designs.append(result)
        return result

    def design_matsubara(self, indices):
        z = 1j * (2 * np.asarray(indices) + 1) * np.pi / self.beta
        result = 1 / (z[:, None] - self.energies)
        self.designs.append(result)
        return result


def config(**kwargs):
    return replace(IPTConfig(), beta=3.0, n_iw=16, n_tau=65, **kwargs)


def pole_plan(cfg=None, indices=None):
    cfg = config() if cfg is None else cfg
    return IPTPlan(
        PoleBasis([-1, 0, 1], beta=cfg.beta, wmax=1),
        PoleBasis(np.arange(-3, 4), beta=cfg.beta, wmax=3),
        config=cfg,
        frequency_fit_indices=indices,
    )


def frequencies(plan):
    return 1j * (2 * plan.indices + 1) * np.pi / plan.config.beta


def test_config_defaults_and_primitive_normalization():
    actual = IPTConfig()
    assert actual == IPTConfig(10, 1, 2, 512, 6145, 0.5, 200, 1e-8, 1e-13)
    normalized = IPTConfig(n_iw=np.int64(3), beta=np.float64(2))
    assert type(normalized.n_iw) is int
    assert type(normalized.beta) is float
    with pytest.raises(FrozenInstanceError):
        actual.beta = 5


@pytest.mark.parametrize(
    "field,value",
    [
        ("beta", 0),
        ("beta", np.inf),
        ("beta", True),
        ("beta", "3"),
        ("hopping", 0),
        ("hopping", -1),
        ("hopping", 1e308),
        ("hopping", 1e-300),
        ("interaction", -1),
        ("interaction", 1e308),
        ("interaction", 1j),
        ("mixing", 0),
        ("mixing", 1.0001),
        ("tolerance", 0),
        ("tolerance", np.nan),
        ("rcond", -1),
        ("rcond", False),
        ("rcond", np.inf),
        ("n_iw", 0),
        ("n_iw", 2.0),
        ("n_iw", True),
        ("n_tau", 1),
        ("max_iterations", 0),
    ],
)
def test_invalid_configuration_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        IPTConfig(**{field: value})


def test_reflected_product_and_sigma_match_independent_asymmetric_poles_and_qr():
    cfg = config(interaction=1.7)
    weiss = FermionicPoleModel(cfg.beta, [-0.7, 1.2], [0.65, 0.35])
    exact_sigma = ipt_self_energy_poles(weiss, cfg.interaction)
    green = PoleBasis(weiss.energies, beta=cfg.beta, wmax=1.3)
    sigma = PoleBasis(np.unique(exact_sigma.energies), beta=cfg.beta, wmax=3 * 1.3)
    plan = IPTPlan(green, sigma, config=cfg)
    z = frequencies(plan)
    truth_g0 = weiss.evaluate_matsubara(plan.indices)
    bath = z - 1 / truth_g0
    update = plan.update(bath)
    np.testing.assert_allclose(update.weiss_coefficients, weiss.weights, atol=2e-14)
    np.testing.assert_allclose(
        update.weiss_tau, weiss.evaluate_tau(plan.tau), atol=1e-14
    )
    np.testing.assert_allclose(
        update.weiss_reflected_tau,
        weiss.evaluate_tau(cfg.beta - plan.tau),
        atol=1e-14,
    )
    np.testing.assert_allclose(
        update.sigma_tau_product, exact_sigma.evaluate_tau(plan.tau), atol=2e-14
    )
    assert (
        np.max(abs(update.sigma_tau_product - cfg.interaction**2 * update.weiss_tau**3))
        > 0.1
    )
    np.testing.assert_allclose(
        update.sigma_iw, exact_sigma.evaluate_matsubara(plan.indices), atol=3e-14
    )
    np.testing.assert_allclose(update.sigma_tau, update.sigma_tau_product, atol=2e-14)
    expected_g = 1 / (1 / truth_g0 - exact_sigma.evaluate_matsubara(plan.indices))
    np.testing.assert_allclose(update.g_iw, expected_g, atol=2e-14)
    # Pivoted QR supplies an independent coefficient oracle, rather than
    # comparing the cached SVD implementation with another call to itself.
    qr_coefficients = lstsq(
        green.design_matsubara(plan.indices), expected_g, lapack_driver="gelsy"
    )[0]
    np.testing.assert_allclose(update.g_coefficients, qr_coefficients, atol=2e-14)
    np.testing.assert_allclose(
        update.g_tau, green.design_tau(plan.tau) @ qr_coefficients, atol=2e-14
    )
    assert update.diagnostics.sigma_time_fit.effective_rank == sigma.rank
    assert update.diagnostics.certification_status == "not_certified"
    assert (
        update.diagnostics.hartree_convention
        == "mu_U_over_2_cancels_static_sigma_U_over_2"
    )


def test_atomic_weiss_dynamic_sigma_and_candidate_green_match_exact_fock_atom():
    plan = pole_plan(config(interaction=2))
    update = plan.update(np.zeros(plan.indices.size, complex))
    z = frequencies(plan)
    atom = hubbard_atom(plan.config.beta, 2)
    np.testing.assert_allclose(update.sigma_iw, 1 / z, atol=3e-14)
    np.testing.assert_allclose(
        update.sigma_tau, -np.ones(plan.tau.size) / 2, atol=2e-14
    )
    np.testing.assert_allclose(
        update.g_iw, atom.evaluate_matsubara(plan.indices), atol=2e-14
    )
    np.testing.assert_allclose(update.g_tau, atom.evaluate_tau(plan.tau), atol=3e-14)
    # This is one algebraic atomic update; t>0 has not magically made it a
    # self-consistent Hubbard/Bethe solution.
    assert update.unmixed_defect > 0.1


def test_noninteracting_semicircle_quadrature_and_two_passes():
    plan = pole_plan(config(interaction=0))
    result = plan.solve(record_trajectory=True)
    assert result.converged
    assert result.stop_reason == "self_consistency_tolerance_reached"
    assert len(result.history) == 2
    assert result.final_update is not None
    assert result.certification_status == "not_certified"
    assert (
        result.convergence_scope
        == "unprojected_bethe_defect_on_frequency_observations_only"
    )
    assert np.max(abs(result.sigma_iw)) == 0

    def rho(energy):
        return np.sqrt(4 - energy**2) / (2 * np.pi)

    for position in (0, 14, 15, 16, 19, 31):
        z = frequencies(plan)[position]
        truth = quad(lambda energy, z=z: (rho(energy) / (z - energy)).real, -2, 2)[0]
        truth += (
            1j * quad(lambda energy, z=z: (rho(energy) / (z - energy)).imag, -2, 2)[0]
        )
        assert abs(result.g_iw[position] - truth) < 2e-13
    np.testing.assert_allclose(result.bath_iw, result.g_iw, atol=1e-14)
    assert result.g_iw is result.final_update.g_iw
    assert result.sigma_tau_product is result.final_update.sigma_tau_product
    assert result.g_coefficients is result.final_update.g_coefficients
    assert result.diagnostics is result.final_update.diagnostics
    assert result.unmixed_defect == result.final_update.unmixed_defect
    assert (
        result.full_grid_unmixed_defect == result.final_update.full_grid_unmixed_defect
    )
    failed_budget = pole_plan(config(interaction=0, max_iterations=1)).solve()
    assert not failed_budget.converged
    assert failed_budget.stop_reason == "iteration_budget_exhausted"
    assert len(failed_budget.history) == 1


def test_tiny_mixed_step_cannot_hide_large_unmixed_defect():
    plan = pole_plan(config(mixing=1e-12, max_iterations=3, tolerance=1e-8))
    result = plan.solve(np.zeros(plan.indices.size, complex))
    assert not result.converged
    assert all(step.unmixed_defect > 0.1 for step in result.history)
    assert all(step.mixed_step < 1e-8 for step in result.history)


def test_stop_requires_strict_threshold_and_two_consecutive_passes(monkeypatch):
    plan = pole_plan(config(interaction=0, tolerance=0.125, max_iterations=6))
    real_update = plan.update
    defects = iter([0.0, 0.125, 0.0, 0.2, 0.0, 0.0])

    def prescribed(bath):
        return replace(real_update(bath), unmixed_defect=next(defects))

    monkeypatch.setattr(plan, "update", prescribed)
    result = plan.solve()
    assert result.converged
    assert len(result.history) == 6
    assert [step.unmixed_defect for step in result.history] == [0, 0.125, 0, 0.2, 0, 0]


def test_subset_holdouts_never_enter_fitting_controller_or_later_feedback():
    observed = np.array([0, -1, 1, -2, 3, -4, 6, -7])
    plan = pole_plan(config(max_iterations=4), observed)
    observed[:] = 12
    np.testing.assert_array_equal(
        plan.frequency_fit_indices, [0, -1, 1, -2, 3, -4, 6, -7]
    )
    positions = plan.frequency_fit_indices + plan.config.n_iw
    omitted = np.ones(plan.indices.size, dtype=bool)
    omitted[positions] = False
    bath = np.zeros(plan.indices.size, complex)
    changed_bath = bath.copy()
    changed_bath[omitted] = 0.3 - 0.2j
    first, changed = plan.update(bath), plan.update(changed_bath)
    for field in (
        "weiss_coefficients",
        "weiss_tau",
        "weiss_reflected_tau",
        "sigma_coefficients",
        "sigma_iw",
        "sigma_tau_product",
        "g_coefficients",
        "g_tau",
    ):
        np.testing.assert_array_equal(getattr(first, field), getattr(changed, field))
    assert first.unmixed_defect == changed.unmixed_defect
    assert np.max(abs(first.g_iw[omitted] - changed.g_iw[omitted])) > 1e-3
    results = [plan.solve(b, record_trajectory=True) for b in (bath, changed_bath)]
    assert results[0].converged == results[1].converged
    assert len(results[0].history) == len(results[1].history)
    for original, perturbed in zip(
        results[0].bath_trajectory, results[1].bath_trajectory, strict=True
    ):
        np.testing.assert_array_equal(original[positions], perturbed[positions])
    assert [h.unmixed_defect for h in results[0].history] == [
        h.unmixed_defect for h in results[1].history
    ]


def test_positive_only_observations_are_preserved_without_conjugate_completion():
    plan = pole_plan(indices=[3, 0, 1])
    update = plan.update(np.zeros(plan.indices.size, complex))
    np.testing.assert_array_equal(update.frequency_fit_indices, [3, 0, 1])
    assert update.diagnostics.weiss_frequency_fit.design_shape == (3, 3)


@pytest.mark.parametrize("attribute", ["config", "green_basis", "selfenergy_basis"])
def test_prepared_plan_identity_cannot_be_reassigned(attribute):
    plan = pole_plan()
    before = plan.update(np.zeros(plan.indices.size, complex))
    replacement = (
        replace(plan.config, beta=4)
        if attribute == "config"
        else PoleBasis([2], beta=4, wmax=6)
    )
    with pytest.raises(AttributeError):
        setattr(plan, attribute, replacement)
    after = plan.update(np.zeros(plan.indices.size, complex))
    assert after.config is before.config
    np.testing.assert_array_equal(after.sigma_iw, before.sigma_iw)


def test_full_grid_defect_is_reported_but_never_used_for_observed_stop():
    plan = pole_plan(config(interaction=0, max_iterations=2), indices=[-1, 0])
    bath = plan.solve().bath_iw.copy()
    omitted = np.ones(plan.indices.size, dtype=bool)
    omitted[plan.frequency_fit_indices + plan.config.n_iw] = False
    bath[omitted] = 0
    result = plan.solve(bath)
    assert result.converged
    assert all(h.unmixed_defect < plan.config.tolerance for h in result.history)
    assert all(h.full_grid_unmixed_defect > 0.01 for h in result.history)
    assert result.full_grid_unmixed_defect > plan.config.tolerance
    assert result.convergence_scope.endswith("frequency_observations_only")


def test_rank_deficient_complex_fit_retains_minimum_norm_and_explicit_rank():
    cfg = config(interaction=0)
    # Two identical columns and one complex observation have a direct
    # minimum-norm solution: divide the coefficient equally between columns.
    green = PoleBasis([0, 0], beta=cfg.beta, wmax=1)
    sigma = PoleBasis([0, 0], beta=cfg.beta, wmax=3)
    plan = IPTPlan(green, sigma, config=cfg, frequency_fit_indices=[0])
    bath = np.zeros(plan.indices.size, complex)
    bath[cfg.n_iw] = 0.2 + 0.1j
    update = plan.update(bath)
    z = frequencies(plan)[cfg.n_iw]
    coefficient = z / (z - bath[cfg.n_iw])
    np.testing.assert_allclose(update.weiss_coefficients, [coefficient / 2] * 2)
    assert np.max(abs(update.weiss_coefficients.imag)) > 0.01
    assert update.diagnostics.weiss_frequency_fit.effective_rank == 1
    assert update.diagnostics.weiss_frequency_fit.design_shape == (1, 2)
    assert update.diagnostics.sigma_time_fit.effective_rank == 1
    assert update.diagnostics.certification_status == "not_certified"
    truncated = IPTPlan(green, sigma, config=replace(cfg, rcond=2))
    candidate = truncated.update(bath)
    assert candidate.diagnostics.weiss_frequency_fit.effective_rank == 0
    assert np.max(abs(candidate.g_coefficients)) == 0
    assert candidate.diagnostics.green_projection_residual > 0.1
    # A large discarded residual remains explicit even if an algebraic raw
    # noninteracting candidate happens to satisfy a controller elsewhere.
    assert candidate.diagnostics.sigma_time_fit.certification_status == "not_certified"


@pytest.mark.parametrize(
    "indices", [[], [0, 0], [16], [-17], [0.0], [True], [[0]], [np.nan]]
)
def test_invalid_observation_subset_rejected(indices):
    with pytest.raises(ValueError):
        pole_plan(indices=indices)


def test_two_svds_are_cached_and_designs_input_and_all_records_are_owned(monkeypatch):
    calls = []
    original = np.linalg.svd

    def observed(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "svd", observed)
    plan = pole_plan(config(max_iterations=3))
    assert len(calls) == 2
    bath = np.zeros(plan.indices.size, complex)
    first = plan.update(bath)
    # Provider-returned design arrays and subsequent caller mutation cannot
    # rewrite cached arithmetic or a previous receipt.
    for basis in (plan.green_basis, plan.selfenergy_basis):
        for design in basis.designs:
            design[:] = 0
    again = plan.update(bath)
    np.testing.assert_array_equal(first.sigma_iw, again.sigma_iw)
    result = plan.solve(bath, record_trajectory=True)
    assert len(calls) == 2
    bath[:] = 99
    assert np.max(abs(first.bath_iw)) == 0
    for record in (first, result):
        for name in record.__dataclass_fields__:
            value = getattr(record, name)
            if isinstance(value, np.ndarray):
                assert not value.flags.writeable
                with pytest.raises(ValueError):
                    value.setflags(write=True)
    for frame in result.bath_trajectory:
        with pytest.raises(ValueError):
            frame.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        first.unmixed_defect = 0
    with pytest.raises(FrozenInstanceError):
        result.converged = True
    assert isinstance(result.history[0], IPTIteration)
    assert pole_plan(config(max_iterations=1)).solve().bath_trajectory == ()


def test_trajectory_records_actual_inputs_and_terminal_proposal_is_not_applied():
    plan = pole_plan(config(max_iterations=4, mixing=0.3))
    result = plan.solve(record_trajectory=True)
    assert not result.converged
    assert len(result.history) == len(result.bath_trajectory) == 4
    for i, bath in enumerate(result.bath_trajectory):
        update = plan.update(bath)
        proposed = (1 - 0.3) * bath + 0.3 * update.g_iw
        assert result.history[i].mixed_step == pytest.approx(
            np.max(abs(proposed - bath))
        )
        if i < 3:
            np.testing.assert_array_equal(result.bath_trajectory[i + 1], proposed)
    np.testing.assert_array_equal(result.bath_iw, result.bath_trajectory[-1])
    np.testing.assert_array_equal(result.bath_iw, result.final_update.bath_iw)
    assert np.max(abs(result.bath_iw - proposed)) > 1e-5


@pytest.mark.parametrize(
    "value",
    [[0], np.zeros((2, 16)), [True] * 32, ["0"] * 32, [np.inf] * 32, [np.nan] * 32],
)
def test_invalid_caller_bath_raises_instead_of_returning_failed_iteration(value):
    plan = pole_plan()
    with pytest.raises(ValueError):
        plan.update(value)
    with pytest.raises(ValueError):
        plan.solve(value)


def test_first_numerical_failure_retains_failure_and_no_substitute_solution():
    plan = pole_plan()
    singular_bath = frequencies(plan)
    with pytest.raises(NumericalBreakdownError, match="Weiss"):
        plan.update(singular_bath)
    result = plan.solve(singular_bath, record_trajectory=True)
    assert not result.converged
    assert result.stop_reason == "numerical_breakdown"
    assert result.final_update is None
    assert result.history == result.bath_trajectory == ()
    assert result.g_iw is result.g_tau is result.sigma_iw is result.sigma_tau is None
    assert result.g_coefficients is result.sigma_coefficients is None
    assert result.weiss_iw is result.weiss_tau is result.weiss_coefficients is None
    assert (
        result.diagnostics
        is result.unmixed_defect
        is result.full_grid_unmixed_defect
        is None
    )
    assert "iteration 1: NumericalBreakdownError" in result.failure_detail
    np.testing.assert_array_equal(result.failed_bath_iw, singular_bath)
    with pytest.raises(ValueError):
        result.failed_bath_iw.setflags(write=True)


def test_later_breakdown_preserves_last_valid_candidate_and_failed_input(monkeypatch):
    plan = pole_plan()
    actual = plan.update
    attempted = []

    def fail_second(bath):
        attempted.append(np.array(bath))
        if len(attempted) == 2:
            raise NumericalBreakdownError("independent injected backend failure")
        return actual(bath)

    monkeypatch.setattr(plan, "update", fail_second)
    result = plan.solve(record_trajectory=True)
    assert not result.converged
    assert result.stop_reason == "numerical_breakdown"
    assert len(result.history) == len(result.bath_trajectory) == 1
    np.testing.assert_array_equal(result.bath_iw, attempted[0])
    np.testing.assert_array_equal(result.failed_bath_iw, attempted[1])
    np.testing.assert_array_equal(result.g_iw, actual(attempted[0]).g_iw)
    assert "iteration 2" in result.failure_detail


def test_basis_validation_budget_and_arithmetic_context_fail_closed():
    green, sigma = PoleBasis([0]), PoleBasis([0], wmax=3)
    with pytest.raises(TypeError, match="config"):
        IPTPlan(green, sigma, config={})
    with pytest.raises(TypeError, match="PhysicalBasis"):
        IPTPlan(object(), sigma, config=config())
    with pytest.raises(ValueError, match="beta"):
        IPTPlan(green, sigma, config=replace(config(), beta=4))
    with pytest.raises(ValueError, match="three times"):
        IPTPlan(green, green, config=config())
    sigma.rank = 2**26
    with pytest.raises(ValueError, match="preparation budget"):
        IPTPlan(green, sigma, config=config())
    assert green.designs == sigma.designs == []
    sigma.rank = 2
    with pytest.raises(ValueError, match="shape"):
        IPTPlan(green, sigma, config=config())
    sigma.rank = 1
    plan = IPTPlan(green, sigma, config=config())
    with pytest.raises(ValueError, match="boolean"):
        plan.solve(record_trajectory=1)
    with np.errstate(all="raise"):
        before = np.geterr()
        plan.solve(np.zeros(plan.indices.size, complex))
        assert np.geterr() == before


@pytest.mark.parametrize("bad_design", [np.nan, np.inf, "bad", True])
def test_invalid_provider_designs_are_rejected(bad_design):
    class InvalidBasis(PoleBasis):
        def design_tau(self, tau):
            return np.full((len(tau), self.rank), bad_design)

    with pytest.raises(ValueError):
        IPTPlan(InvalidBasis([0]), PoleBasis([0], wmax=3), config=config())


def test_real_public_fermionic_bases_can_run_a_bounded_interacting_loop():
    cfg = config(interaction=1, max_iterations=100, tolerance=1e-8)
    green = FermionicBasis(beta=cfg.beta, wmax=4, tol=1e-10, grid_order=8)
    sigma = FermionicBasis(beta=cfg.beta, wmax=12, tol=1e-10, grid_order=8)
    result = IPTPlan(green, sigma, config=cfg).solve()
    assert result.converged
    assert len(result.history) > 2
    assert result.history[-1].unmixed_defect < cfg.tolerance
    assert result.history[-2].unmixed_defect < cfg.tolerance
    assert np.max(result.g_iw[cfg.n_iw :].imag) < 0
    assert np.max(result.sigma_iw[cfg.n_iw :].imag) < 0
    # No conjugacy projection is implemented: this arises from this PH fixture.
    np.testing.assert_allclose(
        result.g_iw[: cfg.n_iw], result.g_iw[cfg.n_iw :][::-1].conj(), atol=2e-11
    )
