"""Run with ``uv run python -m examples.reusable_basis`` from the checkout."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from kernelgecp import FermionicBasis, FermionicFit


def main() -> None:
    beta = 4.0
    energy = 0.7
    basis = FermionicBasis(beta=beta, wmax=3.0, tol=1e-9)
    tau = np.linspace(0, beta, 101)
    atomic_green = -np.exp(-tau * energy) / (1 + np.exp(-beta * energy))
    # One observation axis plus two matrix axes. Complex values are supported.
    amplitudes = np.array([[1.0, 0.2j], [-0.2j, 2.0]])
    values = atomic_green[:, None, None] * amplitudes
    fit = basis.fit_tau(values, tau)
    # Another RHS on the identical grid reuses the prepared SVD.
    basis.fit_tau(2 * values, tau)
    assert basis.cache_info["fit_plan_entries"] == 1

    holdout = (np.arange(80) + 0.5) * beta / 80
    expected = -np.exp(-holdout * energy) / (1 + np.exp(-beta * energy))
    time_error = float(
        np.max(np.abs(fit.evaluate_tau(holdout) - expected[:, None, None] * amplitudes))
    )
    # n is the ordinary integer index, NOT the odd number 2*n+1.
    n = np.array([-10, -1, 0, 2, 9])
    expected_iw = 1 / (1j * (2 * n + 1) * np.pi / beta - energy)
    frequency_error = float(
        np.max(
            np.abs(fit.evaluate_matsubara(n) - expected_iw[:, None, None] * amplitudes)
        )
    )
    assert time_error < 1e-7 and frequency_error < 1e-7

    with TemporaryDirectory(prefix="kernelgecp-tutorial-") as directory:
        path = Path(directory) / "green.npz"
        fit.save(path)
        restored = FermionicFit.load(path)
        assert np.array_equal(restored.evaluate_matsubara(n), fit.evaluate_matsubara(n))

    print(f"rank={basis.rank}, time_error={time_error:.3e}")
    print(f"frequency_error={frequency_error:.3e}, {basis.cache_info}")
    print(basis.diagnostics)
    print("Errors above are sampled checks for this fixture, not certificates.")


if __name__ == "__main__":
    main()
