import numpy as np

from kernelgecp import (
    FermionicKernel,
    dlr_rank_bound,
    dyadic_taylor_rank_bound,
    exp_family_rank_bound,
    fermionic_dyadic_taylor_approximation,
    fermionic_separated_approximation,
)


def test_explicit_rank_counts() -> None:
    assert exp_family_rank_bound(16, 1 / 256) == 25
    assert dlr_rank_bound(16, 1e-6) > 50
    assert dyadic_taylor_rank_bound(16, 1e-6) == 1600


def test_verified_dyadic_taylor_construction_at_band_boundaries() -> None:
    approximation = fermionic_dyadic_taylor_approximation(64, 1e-6)
    t = np.asarray([0.0, 0.125, 0.5, 1.0])[:, None]
    omega = np.asarray([-64.0, -32.0, -16.0, -1.0, 0.0, 1.0, 16.0, 32.0, 64.0])[None, :]
    approximate = approximation.evaluate(t, omega)
    exact = FermionicKernel(64)(t, omega)
    error = float(np.max(np.abs(approximate - exact)))
    assert approximation.rank == 16 * 20 * 7
    assert error <= approximation.error_bound * 1.01 + 1e-14


def test_verified_dyadic_taylor_extreme_cutoff_and_validation() -> None:
    approximation = fermionic_dyadic_taylor_approximation(1e6, 1e-5)
    t = np.asarray([0.0, 0.5, 1.0])[:, None]
    omega = np.asarray([-1e6, -1.0, 0.0, 1.0, 1e6])[None, :]
    approximate = approximation.evaluate(t, omega)
    exact = FermionicKernel(1e6)(t, omega)
    assert np.max(np.abs(approximate - exact)) <= approximation.error_bound + 1e-14
    with np.testing.assert_raises(ValueError):
        approximation.evaluate(0.5, 1e6 + 1)


def test_composite_interpolant_matches_nodes_and_dense_values() -> None:
    approximation = fermionic_separated_approximation(16, order=12)
    t = np.linspace(0, 1, 21)
    at_nodes = approximation.evaluate(t, approximation.omega_nodes)
    exact_nodes = approximation.kernel(t[:, None], approximation.omega_nodes[None, :])
    assert np.max(np.abs(at_nodes - exact_nodes)) < 2e-12
    omega = np.linspace(-16, 16, 501)
    approximate = approximation.evaluate(t, omega)
    exact = FermionicKernel(16)(t[:, None], omega[None, :])
    assert np.max(np.abs(approximate - exact)) < 2e-7
