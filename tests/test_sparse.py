import numpy as np

from kernelgecp import FermionicKernel, sparse_representation


def test_two_delta_fixture_uses_two_atoms() -> None:
    kernel = FermionicKernel(10)
    t = np.linspace(0, 1, 121)
    frequencies = np.array([-2.5, 3.0])
    weights = np.array([0.8, -0.35])
    values = kernel(t[:, None], frequencies[None, :]) @ weights
    result = sparse_representation(
        t,
        values,
        cutoff=10,
        tolerance=1e-8,
        max_atoms=4,
        candidate_frequencies=frequencies,
    )
    assert result.converged
    assert result.frequencies.size == 2
    assert result.validation_error < 1e-8
