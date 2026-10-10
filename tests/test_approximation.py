import itertools

import mpmath as mp
import numpy as np

from kernelgecp import (
    DyadicTaylorApproximation,
    FermionicKernel,
    dlr_rank_bound,
    dyadic_taylor_rank_bound,
    exp_family_rank_bound,
    fermionic_dyadic_taylor_approximation,
    fermionic_separated_approximation,
)


def _stable_masked_eighth_order_error(t: mp.mpf, omega: mp.mpf) -> mp.mpf:
    """Evaluate the p=1 masked error from its convergent exponential tail."""

    x = t * omega
    if omega <= 1:
        active = True
    else:
        band = int(mp.ceil(mp.log(omega, 2)))
        active = t * mp.mpf(2) ** (band - 1) <= 1
    if not active:
        return mp.exp(-x)
    return mp.nsum(lambda k: (-x) ** k / mp.factorial(k), [8, mp.inf])


def _eighth_order_beta_moment(x: mp.mpf) -> mp.mpf:
    return mp.quad(lambda u: (1 - u) ** 7 * mp.exp(-u * x), [0, 1])


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


def test_masked_dyadic_remainder_reverses_order_two_sign() -> None:
    with mp.workdps(100):
        approximation = DyadicTaylorApproximation(
            cutoff=4, accuracy_bits=1, scale=2, precision_bits=384
        )
        rows = [mp.mpf("0.5"), mp.mpf("0.6")]
        columns = [mp.mpf("2"), mp.mpf("2.5")]

        stable = [
            [_stable_masked_eighth_order_error(t, omega) for omega in columns]
            for t in rows
        ]
        implemented = [
            [
                mp.exp(-t * omega)
                - approximation._positive_value(t, omega) * (1 + mp.exp(-omega))
                for omega in columns
            ]
            for t in rows
        ]
        for stable_row, implemented_row in zip(stable, implemented, strict=True):
            for expected, actual in zip(stable_row, implemented_row, strict=True):
                assert mp.almosteq(expected, actual)

        determinant = stable[0][0] * stable[1][1] - stable[0][1] * stable[1][0]
        assert mp.mpf("4.96e-6") < determinant < mp.mpf("4.97e-6")

        for lower_time in map(mp.mpf, ["0.45", "0.5"]):
            for upper_time in map(mp.mpf, ["0.55", "0.6", "0.7"]):
                for lower_frequency in map(mp.mpf, ["1.8", "2"]):
                    for upper_frequency in map(mp.mpf, ["2.1", "2.5", "3"]):
                        diagonal = _stable_masked_eighth_order_error(
                            lower_time, lower_frequency
                        ) * _stable_masked_eighth_order_error(
                            upper_time, upper_frequency
                        )
                        off_diagonal = _stable_masked_eighth_order_error(
                            lower_time, upper_frequency
                        ) * _stable_masked_eighth_order_error(
                            upper_time, lower_frequency
                        )
                        assert diagonal - off_diagonal > 0


def test_unmasked_eighth_order_tail_beta_structure_and_order_two_sign() -> None:
    with mp.workdps(100):
        for x in map(mp.mpf, ["0.001", "0.01", "0.1", "0.5", "1", "2"]):
            beta_tail = x**8 / mp.factorial(7) * _eighth_order_beta_moment(x)
            convergent_tail = mp.nsum(
                lambda k, x=x: (-x) ** k / mp.factorial(k), [8, mp.inf]
            )
            direct_tail = mp.exp(-x) - mp.fsum(
                (-x) ** k / mp.factorial(k) for k in range(8)
            )
            assert mp.almosteq(beta_tail, convergent_tail)
            assert mp.almosteq(beta_tail, direct_tail)

        rows = list(map(mp.mpf, ["0.05", "0.1", "0.25", "0.5", "0.75", "1"]))
        columns = list(map(mp.mpf, ["0.05", "0.1", "0.25", "0.5", "1", "1.5", "2"]))
        moments = {
            (time, frequency): _eighth_order_beta_moment(time * frequency)
            for time in rows
            for frequency in columns
        }
        approximation = DyadicTaylorApproximation(
            cutoff=2, accuracy_bits=1, scale=1, precision_bits=384
        )
        smallest_signed_margin = mp.inf
        for time, frequency in moments:
            implemented_error = mp.exp(-time * frequency) - (
                approximation._positive_value(time, frequency)
                * (1 + mp.exp(-frequency))
            )
            beta_error = (
                (time * frequency) ** 8 / mp.factorial(7) * moments[time, frequency]
            )
            assert mp.almosteq(beta_error, implemented_error)

        for lower_time, upper_time in itertools.combinations(rows, 2):
            for lower_frequency, upper_frequency in itertools.combinations(columns, 2):
                determinant = (
                    moments[lower_time, lower_frequency]
                    * moments[upper_time, upper_frequency]
                    - moments[lower_time, upper_frequency]
                    * moments[upper_time, lower_frequency]
                )
                assert determinant < 0
                smallest_signed_margin = min(smallest_signed_margin, -determinant)
        assert smallest_signed_margin > mp.mpf("4.3e-6")


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
