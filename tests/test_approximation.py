import numpy as np

from kernelgecp import (
    FermionicKernel,
    dlr_rank_bound,
    exp_family_rank_bound,
    fermionic_separated_approximation,
)


def test_explicit_rank_counts() -> None:
    assert exp_family_rank_bound(16, 1 / 256) == 25
    assert dlr_rank_bound(16, 1e-6) > 50


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
