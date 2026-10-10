import itertools
import math
from fractions import Fraction

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
    return _eighth_order_beta_power_moment(0, x)


def _eighth_order_beta_power_moment(power: int, x: mp.mpf) -> mp.mpf:
    return mp.quad(lambda u: u**power * (1 - u) ** 7 * mp.exp(-u * x), [0, 1])


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


def test_tilted_beta_moment_gap_and_variance_bound() -> None:
    unweighted_gap = sum(
        Fraction((-1) ** degree * math.comb(7, degree))
        * (Fraction(1, degree + 2) - 2 * Fraction(1, degree + 3))
        for degree in range(8)
    )
    assert unweighted_gap == Fraction(1, 120)

    with mp.workdps(100):
        points = {mp.mpf(index) / 50 for index in range(101)}
        points.update(map(mp.mpf, ["1e-12", "1e-9", "1e-6", "0.001", "0.01", "0.1"]))
        smallest_moment_gap = mp.inf
        smallest_variance_margin = mp.inf
        for x in sorted(points):
            moment_zero = _eighth_order_beta_power_moment(0, x)
            moment_one = _eighth_order_beta_power_moment(1, x)
            moment_two = _eighth_order_beta_power_moment(2, x)
            moment_gap = moment_one - 2 * moment_two
            variance_margin = moment_one * moment_zero - x * (
                moment_two * moment_zero - moment_one**2
            )
            assert moment_gap > 0
            assert variance_margin > 0
            smallest_moment_gap = min(smallest_moment_gap, moment_gap)
            smallest_variance_margin = min(smallest_variance_margin, variance_margin)
        assert smallest_moment_gap > mp.mpf("0.0062")
        assert smallest_variance_margin > mp.mpf("0.0008")


def test_beta_moment_derivatives_and_elasticity_monotonicity() -> None:
    with mp.workdps(100):
        derivative_points = map(mp.mpf, ["0", "1e-12", "1e-6", "0.1", "1", "2"])
        for x in derivative_points:
            moment_one = _eighth_order_beta_power_moment(1, x)
            moment_two = _eighth_order_beta_power_moment(2, x)
            differentiated_zero = mp.diff(
                lambda y: _eighth_order_beta_power_moment(0, y), x
            )
            differentiated_one = mp.diff(
                lambda y: _eighth_order_beta_power_moment(1, y), x
            )
            assert mp.almosteq(differentiated_zero, -moment_one)
            assert mp.almosteq(differentiated_one, -moment_two)

        points = sorted(
            {mp.mpf(index) / 50 for index in range(101)}
            | {mp.mpf(value) for value in ["1e-12", "1e-9", "1e-6"]}
        )
        elasticities: list[mp.mpf] = []
        smallest_negative_derivative = mp.inf
        for x in points:
            moment_zero = _eighth_order_beta_power_moment(0, x)
            moment_one = _eighth_order_beta_power_moment(1, x)
            moment_two = _eighth_order_beta_power_moment(2, x)
            elasticity = -x * moment_one / moment_zero
            elasticity_derivative = (
                x * (moment_two * moment_zero - moment_one**2)
                - moment_one * moment_zero
            ) / moment_zero**2
            assert elasticity_derivative < 0
            elasticities.append(elasticity)
            smallest_negative_derivative = min(
                smallest_negative_derivative, -elasticity_derivative
            )

        assert all(upper < lower for lower, upper in itertools.pairwise(elasticities))
        assert smallest_negative_derivative > mp.mpf("0.078")


def test_beta_moment_scale_ratio_and_local_minor_sign() -> None:
    with mp.workdps(100):
        smallest_ratio_derivative = mp.inf
        for lower_time, upper_time in [
            (mp.mpf("0.05"), mp.mpf("1")),
            (mp.mpf("0.2"), mp.mpf("0.75")),
            (mp.mpf("0.49"), mp.mpf("0.51")),
        ]:
            frequencies = [
                mp.mpf("0.03"),
                mp.mpf("0.1"),
                mp.mpf("0.5"),
                mp.mpf("1"),
                2 / upper_time,
            ]
            ratios: list[mp.mpf] = []
            for frequency in frequencies:
                lower_moment = _eighth_order_beta_power_moment(
                    0, lower_time * frequency
                )
                upper_moment = _eighth_order_beta_power_moment(
                    0, upper_time * frequency
                )
                lower_first = _eighth_order_beta_power_moment(1, lower_time * frequency)
                upper_first = _eighth_order_beta_power_moment(1, upper_time * frequency)
                analytic_derivative = (
                    -lower_time * lower_first * upper_moment
                    + upper_time * lower_moment * upper_first
                ) / upper_moment**2
                differentiated_ratio = mp.diff(
                    lambda omega, lower_time=lower_time, upper_time=upper_time: (
                        _eighth_order_beta_power_moment(0, lower_time * omega)
                        / _eighth_order_beta_power_moment(0, upper_time * omega)
                    ),
                    frequency,
                )
                assert mp.almosteq(analytic_derivative, differentiated_ratio)
                assert analytic_derivative > 0
                ratios.append(lower_moment / upper_moment)
                smallest_ratio_derivative = min(
                    smallest_ratio_derivative, analytic_derivative
                )
            assert all(upper > lower for lower, upper in itertools.pairwise(ratios))
        assert smallest_ratio_derivative > mp.mpf("0.0015")

        rows = list(map(mp.mpf, ["0.07", "0.21", "0.57", "0.93"]))
        columns = list(map(mp.mpf, ["0.03", "0.13", "0.41", "1.1", "2"]))
        smallest_beta_margin = mp.inf
        smallest_tail_margin = mp.inf
        for lower_time, upper_time in itertools.combinations(rows, 2):
            for lower_frequency, upper_frequency in itertools.combinations(columns, 2):
                beta_determinant = _eighth_order_beta_moment(
                    lower_time * lower_frequency
                ) * _eighth_order_beta_moment(
                    upper_time * upper_frequency
                ) - _eighth_order_beta_moment(
                    lower_time * upper_frequency
                ) * _eighth_order_beta_moment(upper_time * lower_frequency)
                tail_determinant = (
                    (lower_time * upper_time * lower_frequency * upper_frequency) ** 8
                    / mp.factorial(7) ** 2
                    * beta_determinant
                )
                assert beta_determinant < 0
                assert tail_determinant < 0
                smallest_beta_margin = min(smallest_beta_margin, -beta_determinant)
                smallest_tail_margin = min(smallest_tail_margin, -tail_determinant)
        assert smallest_beta_margin > mp.mpf("2.4e-5")
        assert smallest_tail_margin > mp.mpf("1e-46")


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
