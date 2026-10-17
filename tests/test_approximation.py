import itertools
import math
import random
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
from kernelgecp.surrogates import fraction_determinant


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


def test_truncated_power_composition_and_higher_order_minor_sweep() -> None:
    with mp.workdps(100):
        for time, frequency, cutoff in [
            (mp.mpf("0.2"), mp.mpf("0.1"), mp.mpf("2")),
            (mp.mpf("1.3"), mp.mpf("0.7"), mp.mpf("1.5")),
            (mp.mpf("2"), mp.mpf("1"), mp.mpf("1.2")),
            (mp.mpf("0.05"), mp.mpf("2"), mp.mpf("3")),
        ]:
            beta_moment = _eighth_order_beta_moment(time * frequency)
            scaled_beta = frequency**8 * beta_moment
            truncated_integral = mp.quad(
                lambda value, frequency=frequency, time=time: (
                    (frequency - value) ** 7 * mp.exp(-time * value)
                ),
                [0, frequency],
            )
            common_cutoff_integral = mp.quad(
                lambda value, frequency=frequency, time=time: (
                    max(frequency - value, mp.mpf(0)) ** 7 * mp.exp(-time * value)
                ),
                [0, frequency, cutoff],
            )
            tail = (time * frequency) ** 8 / mp.factorial(7) * beta_moment
            composition_tail = time**8 / mp.factorial(7) * common_cutoff_integral
            assert mp.almosteq(scaled_beta, truncated_integral)
            assert mp.almosteq(truncated_integral, common_cutoff_integral)
            assert mp.almosteq(tail, composition_tail)

        grid = list(map(mp.mpf, ["0.02", "0.08", "0.2", "0.5", "1", "1.4"]))
        moments = {
            (time, frequency): mp.hyp1f1(1, 9, -time * frequency) / 8
            for time in grid
            for frequency in grid
        }
        expected_margins = {
            2: mp.mpf("6e-6"),
            3: mp.mpf("4e-12"),
            4: mp.mpf("1.7e-18"),
            5: mp.mpf("2e-24"),
        }
        for order, expected_margin in expected_margins.items():
            signature = (-1) ** (order * (order - 1) // 2)
            smallest_signed_margin = mp.inf
            for rows in itertools.combinations(grid, order):
                for columns in itertools.combinations(grid, order):
                    if rows[-1] * columns[-1] > 2:
                        continue
                    determinant = mp.det(
                        mp.matrix(
                            [
                                [moments[row, column] for column in columns]
                                for row in rows
                            ]
                        )
                    )
                    signed_determinant = signature * determinant
                    assert signed_determinant > 0
                    smallest_signed_margin = min(
                        smallest_signed_margin, signed_determinant
                    )
            assert smallest_signed_margin > expected_margin


def test_truncated_power_order_two_support_geometry_exactly() -> None:
    def kernel(row: Fraction, column: Fraction) -> Fraction:
        return max(column - row, Fraction()) ** 7

    def minor(
        rows: tuple[Fraction, Fraction], columns: tuple[Fraction, Fraction]
    ) -> Fraction:
        return kernel(rows[0], columns[0]) * kernel(rows[1], columns[1]) - kernel(
            rows[0], columns[1]
        ) * kernel(rows[1], columns[0])

    assert minor((Fraction(2), Fraction(3)), (Fraction(0), Fraction(1))) == 0
    assert minor((Fraction(0), Fraction(2)), (Fraction(0), Fraction(3))) == 0
    assert minor((Fraction(0), Fraction(2)), (Fraction(1), Fraction(3))) == 1
    assert minor((Fraction(0), Fraction(1)), (Fraction(1), Fraction(2))) == 1
    assert minor((Fraction(0), Fraction(1)), (Fraction(2), Fraction(3))) == 14197

    points = list(map(Fraction, range(-3, 5)))
    for rows in itertools.combinations(points, 2):
        for columns in itertools.combinations(points, 2):
            determinant = minor(rows, columns)
            paired_support = rows[0] < columns[0] and rows[1] < columns[1]
            assert determinant >= 0
            assert (determinant > 0) is paired_support


def test_normalized_truncated_power_convolution_exactly() -> None:
    def normalized(degree: int, row: Fraction, column: Fraction) -> Fraction:
        if row >= column:
            return Fraction()
        return (column - row) ** degree / math.factorial(degree)

    def expanded_polynomial_integral(
        degree: int, row: Fraction, column: Fraction
    ) -> Fraction:
        # Independently expand (column - z)^degree before integrating termwise.
        return sum(
            (
                Fraction(math.comb(degree, exponent))
                * column ** (degree - exponent)
                * (-1) ** exponent
                * (column ** (exponent + 1) - row ** (exponent + 1))
                / (exponent + 1)
            )
            for exponent in range(degree + 1)
        ) / math.factorial(degree)

    cases = [
        (Fraction(-5), Fraction(-2), Fraction(3), Fraction(9)),
        (Fraction(-7, 3), Fraction(-1, 2), Fraction(5, 4), Fraction(11, 3)),
        (Fraction(0), Fraction(0), Fraction(0), Fraction(4)),
        (Fraction(-3), Fraction(2, 5), Fraction(2, 5), Fraction(8)),
    ]
    for degree in range(8):
        for lower, row, column, upper in cases:
            assert lower <= row <= column <= upper
            integral = expanded_polynomial_integral(degree, row, column)
            assert integral == normalized(degree + 1, row, column)

            # The two factors vanish off (row, column), so enlarging either
            # common-domain endpoint contributes exactly zero.
            assert normalized(0, row, lower) == 0
            assert normalized(degree, upper, column) == 0

    row, column = Fraction(-3, 2), Fraction(7, 3)
    assert math.factorial(7) * normalized(7, row, column) == (column - row) ** 7


def test_all_normalized_truncated_power_minors_are_nonnegative_exactly() -> None:
    def normalized(degree: int, row: Fraction, column: Fraction) -> Fraction:
        if row >= column:
            return Fraction()
        return (column - row) ** degree / math.factorial(degree)

    points = list(map(Fraction, range(-2, 4)))
    checked = 0
    for degree in range(8):
        for order in range(1, 5):
            for rows in itertools.combinations(points, order):
                for columns in itertools.combinations(points, order):
                    matrix = [
                        [normalized(degree, row, column) for column in columns]
                        for row in rows
                    ]
                    assert fraction_determinant(matrix) >= 0
                    checked += 1

    assert checked == 7088


def test_andreief_identity_for_three_monomials_exactly() -> None:
    order = 3
    permutations = list(itertools.permutations(range(order)))

    def sign(permutation: tuple[int, ...]) -> int:
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(order)
            for right in range(left + 1, order)
        )
        return (-1) ** inversions

    # The pair-integral matrix is the 3 x 3 Hilbert matrix:
    # integral_0^1 x^(i+j) dx = 1/(i+j+1).
    hilbert = [
        [Fraction(1, row + column + 1) for column in range(order)]
        for row in range(order)
    ]
    determinant_side = math.factorial(order) * fraction_determinant(hilbert)

    # Independently expand both sampled determinants and integrate each
    # coordinate monomial over the unit cube.
    product_measure_side = sum(
        Fraction(sign(left_permutation) * sign(right_permutation))
        * math.prod(
            Fraction(
                1,
                left_permutation[coordinate] + right_permutation[coordinate] + 1,
            )
            for coordinate in range(order)
        )
        for left_permutation in permutations
        for right_permutation in permutations
    )

    assert determinant_side == Fraction(1, 360)
    assert product_measure_side == determinant_side


def test_common_sample_permutation_signs_cancel_exactly() -> None:
    ordered_sample = (Fraction(1, 4), Fraction(3, 4), Fraction(3, 2))
    columns = (Fraction(1), Fraction(2), Fraction(3))
    products: set[Fraction] = set()

    for sample in itertools.permutations(ordered_sample):
        left_matrix = [[value**degree for value in sample] for degree in range(3)]
        right_matrix = [
            [max(column - value, Fraction()) ** 3 for column in columns]
            for value in sample
        ]
        products.add(
            fraction_determinant(left_matrix) * fraction_determinant(right_matrix)
        )

    assert products == {Fraction(119085, 131072)}

    repeated_sample = (Fraction(1, 4), Fraction(3, 4), Fraction(3, 4))
    repeated_right = [
        [max(column - value, Fraction()) ** 3 for column in columns]
        for value in repeated_sample
    ]
    assert fraction_determinant(repeated_right) == 0


def test_truncated_power_interlacing_witness_is_strict_at_high_orders() -> None:
    for order in range(1, 13):
        columns = tuple(Fraction(index + 1) for index in range(order))
        sample = (columns[0] / 2,) + tuple(
            (columns[index - 1] + columns[index]) / 2 for index in range(1, order)
        )
        matrix = [
            [max(column - value, Fraction()) ** 7 for column in columns]
            for value in sample
        ]
        determinant = fraction_determinant(matrix)

        assert all(sample[index] < columns[index] for index in range(order))
        assert all(
            matrix[row][column] == 0 for row in range(order) for column in range(row)
        )
        assert determinant == Fraction(1, 2 ** (7 * order))
        assert determinant > 0


def test_positive_tail_has_factorial_free_gecp_envelope() -> None:
    mp.mp.dps = 120
    nodes = [mp.mpf(index) / 8 for index in range(1, 9)]
    bound = mp.mpf(1) / math.factorial(8)

    def taylor(x: mp.mpf) -> mp.mpf:
        return mp.fsum((-x) ** degree / math.factorial(degree) for degree in range(8))

    def tail(t: mp.mpf, omega: mp.mpf) -> mp.mpf:
        denominator = 1 + mp.exp(-omega)
        exact = mp.exp(-t * omega) / denominator
        approximation = taylor(t * omega) / denominator
        remainder = (mp.exp(-t * omega) - taylor(t * omega)) / denominator
        assert mp.almosteq(exact, approximation + remainder)
        assert abs(remainder) <= bound
        return remainder

    original = mp.matrix([[tail(t, omega) for omega in nodes] for t in nodes])
    residual = mp.matrix(original)
    selected_rows: list[int] = []
    selected_columns: list[int] = []

    for rank in range(1, len(nodes) + 1):
        pivot_row, pivot_column = max(
            itertools.product(range(len(nodes)), repeat=2),
            key=lambda index: abs(residual[index]),
        )
        pivot = residual[pivot_row, pivot_column]
        assert pivot != 0
        assert abs(pivot) <= bound
        selected_rows.append(pivot_row)
        selected_columns.append(pivot_column)

        core = mp.matrix(
            [
                [original[row, column] for column in selected_columns]
                for row in selected_rows
            ]
        )
        assert abs(mp.det(core)) <= bound**rank * (1 + mp.mpf("1e-90"))

        residual = mp.matrix(
            [
                [
                    residual[row, column]
                    - residual[row, pivot_column] * residual[pivot_row, column] / pivot
                    for column in range(len(nodes))
                ]
                for row in range(len(nodes))
            ]
        )


def test_positive_series_split_aligns_mixed_determinant_coefficients() -> None:
    mp.mp.dps = 140

    def negative_taylor(x: mp.mpf) -> mp.mpf:
        return mp.fsum((-x) ** degree / math.factorial(degree) for degree in range(8))

    def positive_taylor(x: mp.mpf) -> mp.mpf:
        return mp.fsum(x**degree / math.factorial(degree) for degree in range(8))

    def mixed_coefficient(
        rows: list[mp.mpf],
        columns: list[mp.mpf],
        tail_columns: int,
        approximation,
        remainder,
    ) -> mp.mpf:
        order = len(rows)
        return mp.fsum(
            mp.det(
                mp.matrix(
                    [
                        [
                            remainder(row * column)
                            if column_index in selected
                            else approximation(row * column)
                            for column_index, column in enumerate(columns)
                        ]
                        for row in rows
                    ]
                )
            )
            for selected_tuple in itertools.combinations(range(order), tail_columns)
            for selected in [set(selected_tuple)]
        )

    def negative_remainder(x: mp.mpf) -> mp.mpf:
        return mp.exp(-x) - negative_taylor(x)

    obstruction_nodes = [mp.mpf("0.5"), mp.mpf(1)]
    negative_mixed = mixed_coefficient(
        obstruction_nodes,
        obstruction_nodes,
        1,
        negative_taylor,
        negative_remainder,
    )
    # A 2 x 2 negative-exponential determinant should be negative, but its
    # one-tail-column coefficient is positive: the direct split is not
    # coefficientwise sign compatible.
    assert negative_mixed > 0

    def positive_remainder(x: mp.mpf) -> mp.mpf:
        return mp.exp(x) - positive_taylor(x)

    for order in range(1, 7):
        nodes = [mp.mpf(index) / (order + 1) for index in range(1, order + 1)]
        denominator_product = mp.fprod(1 + mp.exp(node) for node in nodes)
        for tail_columns in range(order + 1):
            coefficient = mixed_coefficient(
                nodes,
                nodes,
                tail_columns,
                positive_taylor,
                positive_remainder,
            )
            assert coefficient > 0
            assert coefficient / denominator_product > 0


def test_block_feature_and_mixed_column_coefficients_exactly() -> None:
    rows = (Fraction(1), Fraction(2), Fraction(3))
    columns = (Fraction(1), Fraction(2), Fraction(4))
    order = len(rows)
    low_terms = 3
    terms = 6

    def feature_sum(row: Fraction, column: Fraction, start: int, stop: int) -> Fraction:
        return sum(
            (row * column) ** degree / math.factorial(degree)
            for degree in range(start, stop)
        )

    low = [
        [feature_sum(row, column, 0, low_terms) for column in columns] for row in rows
    ]
    high = [
        [feature_sum(row, column, low_terms, terms) for column in columns]
        for row in rows
    ]

    multilinear_coefficients = []
    for high_count in range(order + 1):
        coefficient = Fraction()
        for selected_tuple in itertools.combinations(range(order), high_count):
            selected = set(selected_tuple)
            coefficient += fraction_determinant(
                [
                    [
                        high[row][column] if column in selected else low[row][column]
                        for column in range(order)
                    ]
                    for row in range(order)
                ]
            )
        multilinear_coefficients.append(coefficient)

    cauchy_binet_coefficients = [Fraction() for _ in range(order + 1)]
    for exponents in itertools.combinations(range(terms), order):
        left = [
            [row**degree / math.factorial(degree) for degree in exponents]
            for row in rows
        ]
        right = [[column**degree for column in columns] for degree in exponents]
        high_count = sum(degree >= low_terms for degree in exponents)
        cauchy_binet_coefficients[high_count] += fraction_determinant(
            left
        ) * fraction_determinant(right)

    assert multilinear_coefficients == cauchy_binet_coefficients
    assert all(coefficient > 0 for coefficient in cauchy_binet_coefficients)


def test_positive_tail_mixed_coefficients_converge() -> None:
    mp.mp.dps = 120
    rows = (mp.mpf("0.2"), mp.mpf("0.55"), mp.mpf("0.9"))
    columns = (mp.mpf("0.15"), mp.mpf("0.6"), mp.mpf("0.95"))
    order = len(rows)

    def taylor_entry(row: mp.mpf, column: mp.mpf, start: int, stop: int) -> mp.mpf:
        product = row * column
        return mp.fsum(
            product**degree / mp.factorial(degree) for degree in range(start, stop)
        )

    low = mp.matrix(
        [[taylor_entry(row, column, 0, 8) for column in columns] for row in rows]
    )
    exact_tail = mp.matrix(
        [
            [mp.exp(row * column) - low[i, j] for j, column in enumerate(columns)]
            for i, row in enumerate(rows)
        ]
    )

    def coefficient(tail: mp.matrix, degree: int) -> mp.mpf:
        return mp.fsum(
            mp.det(
                mp.matrix(
                    [
                        [
                            tail[i, j] if j in selected else low[i, j]
                            for j in range(order)
                        ]
                        for i in range(order)
                    ]
                )
            )
            for selected_tuple in itertools.combinations(range(order), degree)
            for selected in [set(selected_tuple)]
        )

    exact_coefficients = [
        coefficient(exact_tail, degree) for degree in range(order + 1)
    ]
    assert all(value > 0 for value in exact_coefficients)

    errors = []
    for terms in (16, 24, 40, 64):
        finite_tail = mp.matrix(
            [
                [taylor_entry(row, column, 8, terms) for column in columns]
                for row in rows
            ]
        )
        errors.append(
            max(
                abs(coefficient(finite_tail, degree) - exact_coefficients[degree])
                for degree in range(1, order + 1)
            )
        )

    assert errors[-1] < mp.mpf("1e-80")
    assert all(later < earlier for earlier, later in itertools.pairwise(errors))


def test_positive_tail_rank_amplification() -> None:
    mp.mp.dps = 350
    smallest_relative_slack = mp.inf

    for order in range(9, 17):
        nodes = [mp.mpf(index) / (order + 1) for index in range(1, order + 1)]
        low = mp.matrix(
            [
                [
                    mp.fsum(
                        (row * column) ** degree / mp.factorial(degree)
                        for degree in range(8)
                    )
                    for column in nodes
                ]
                for row in nodes
            ]
        )
        full = mp.matrix([[mp.exp(row * column) for column in nodes] for row in nodes])
        tail = full - low
        full_det = mp.det(full)
        assert full_det > 0

        for weight in (mp.mpf(2), mp.mpf(8), mp.mpf(40320)):
            amplified_det = mp.det(low + weight * tail)
            lower_bound = weight ** (order - 8) * full_det
            slack = amplified_det - lower_bound
            assert slack > 0
            smallest_relative_slack = min(
                smallest_relative_slack, slack / amplified_det
            )

    assert smallest_relative_slack > 0


def test_factorial_amplified_positive_kernel_envelope_and_minors() -> None:
    mp.mp.dps = 250
    weight = mp.mpf(mp.factorial(8))
    largest_entry = mp.mpf(0)
    smallest_determinant = mp.inf

    def amplified_entry(row: mp.mpf, column: mp.mpf) -> mp.mpf:
        product = row * column
        low = mp.fsum(product**degree / mp.factorial(degree) for degree in range(8))
        tail = mp.exp(product) - low
        return (low + weight * tail) / (1 + mp.exp(column))

    endpoint_grid = [mp.mpf("1e-20")] + [mp.mpf(index) / 16 for index in range(1, 17)]
    for row in endpoint_grid:
        for column in endpoint_grid:
            value = amplified_entry(row, column)
            assert value > 0
            assert value <= 2
            largest_entry = max(largest_entry, value)

    for order in range(1, 9):
        rows = [mp.mpf(index) / (order + 1) for index in range(1, order + 1)]
        columns = [mp.mpf(index + 1) / (order + 2) for index in range(1, order + 1)]
        determinant = mp.det(
            mp.matrix(
                [[amplified_entry(row, column) for column in columns] for row in rows]
            )
        )
        assert determinant > 0
        assert determinant <= mp.mpf(2) ** order
        smallest_determinant = min(smallest_determinant, determinant)

    assert largest_entry < 2
    assert smallest_determinant > 0


def test_exact_total_positive_hadamard_envelope_through_order_eight() -> None:
    order = 8
    terms = 16
    rows = tuple(Fraction(index, order + 1) for index in range(1, order + 1))
    columns = tuple(Fraction(index, order + 2) for index in range(1, order + 1))

    def truncated_exp(row: Fraction, column: Fraction) -> Fraction:
        product = row * column
        return sum(
            (product**degree / math.factorial(degree) for degree in range(terms)),
            Fraction(),
        )

    matrix = [[truncated_exp(row, column) for column in columns] for row in rows]

    # The finite positive-feature kernel is strictly totally positive on this
    # sample: check every one of its 12,869 nonempty ordered square minors.
    checked_minors = 0
    for minor_order in range(1, order + 1):
        for selected_rows in itertools.combinations(range(order), minor_order):
            for selected_columns in itertools.combinations(range(order), minor_order):
                minor = [
                    [matrix[i][j] for j in selected_columns] for i in selected_rows
                ]
                assert fraction_determinant(minor) > 0
                checked_minors += 1
    assert checked_minors == math.comb(2 * order, order) - 1

    # Exact scalar Schur elimination mirrors the Lean induction: positivity is
    # retained and each Schur diagonal can only decrease.
    schur = matrix
    while schur:
        determinant = fraction_determinant(schur)
        diagonal_product = math.prod(schur[index][index] for index in range(len(schur)))
        assert 0 < determinant <= diagonal_product
        if len(schur) == 1:
            break
        pivot = schur[0][0]
        next_schur = [
            [
                schur[i][j] - schur[i][0] * schur[0][j] / pivot
                for j in range(1, len(schur))
            ]
            for i in range(1, len(schur))
        ]
        for index in range(len(next_schur)):
            assert 0 < next_schur[index][index] <= schur[index + 1][index + 1]
        schur = next_schur


def test_positive_series_ten_pivot_half_contraction() -> None:
    mp.mp.dps = 180
    determinant_bound = mp.mpf(2) ** 10 / mp.factorial(8) ** 2

    for size in (16, 24, 32):
        rows = [mp.mpf(index) / size for index in range(1, size + 1)]
        columns = [mp.mpf(index) / size for index in range(1, size + 1)]
        original = mp.matrix(
            [
                [mp.exp(row * column) / (1 + mp.exp(column)) for column in columns]
                for row in rows
            ]
        )
        residual = mp.matrix(original)
        selected_rows: list[int] = []
        selected_columns: list[int] = []

        for _ in range(10):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[row, column]),
                    -row,
                    -column,
                    row,
                    column,
                )
                for row in range(size)
                for column in range(size)
            )
            pivot = residual[pivot_row, pivot_column]
            assert pivot != 0
            selected_rows.append(pivot_row)
            selected_columns.append(pivot_column)
            pivot_column_values = [residual[row, pivot_column] for row in range(size)]
            pivot_row_values = [residual[pivot_row, column] for column in range(size)]
            residual = mp.matrix(
                [
                    [
                        residual[row, column]
                        - pivot_column_values[row] * pivot_row_values[column] / pivot
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
            )

        selected_core = mp.matrix(
            [
                [original[row, column] for column in selected_columns]
                for row in selected_rows
            ]
        )
        sorted_core = mp.matrix(
            [
                [original[row, column] for column in sorted(selected_columns)]
                for row in sorted(selected_rows)
            ]
        )
        assert mp.almosteq(abs(mp.det(selected_core)), mp.det(sorted_core))
        assert abs(mp.det(selected_core)) <= determinant_bound
        residual_maximum = max(
            abs(residual[row, column]) for row in range(size) for column in range(size)
        )
        assert residual_maximum <= mp.mpf("0.5")


def test_signed_unit_band_twelve_pivot_half_contraction() -> None:
    mp.mp.dps = 180
    determinant_bound = mp.mpf(2) ** 12 / mp.mpf(128) ** 4

    for size in (16, 24, 32):
        rows = [mp.mpf(index) / (size - 1) for index in range(size)]
        columns = [-1 + 2 * mp.mpf(index) / (size - 1) for index in range(size)]
        original = mp.matrix(
            [
                [mp.exp(-row * column) / (1 + mp.exp(-column)) for column in columns]
                for row in rows
            ]
        )
        residual = mp.matrix(original)
        selected_rows: list[int] = []
        selected_columns: list[int] = []

        for _ in range(12):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[row, column]),
                    -row,
                    -column,
                    row,
                    column,
                )
                for row in range(size)
                for column in range(size)
            )
            pivot = residual[pivot_row, pivot_column]
            assert pivot != 0
            selected_rows.append(pivot_row)
            selected_columns.append(pivot_column)
            pivot_column_values = [residual[row, pivot_column] for row in range(size)]
            pivot_row_values = [residual[pivot_row, column] for column in range(size)]
            residual = mp.matrix(
                [
                    [
                        residual[row, column]
                        - pivot_column_values[row] * pivot_row_values[column] / pivot
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
            )

        selected_core = mp.matrix(
            [
                [original[row, column] for column in selected_columns]
                for row in selected_rows
            ]
        )
        sorted_core = mp.matrix(
            [
                [original[row, column] for column in sorted(selected_columns)]
                for row in sorted(selected_rows)
            ]
        )
        assert mp.almosteq(abs(mp.det(selected_core)), mp.det(sorted_core))
        assert abs(mp.det(selected_core)) <= determinant_bound
        residual_maximum = max(
            abs(residual[row, column]) for row in range(size) for column in range(size)
        )
        assert residual_maximum <= mp.mpf("0.5")


def test_translated_unit_bands_have_uniform_twelve_pivot_contraction() -> None:
    mp.mp.dps = 180
    size = 16
    weight = mp.mpf(128)
    determinant_bound = mp.mpf(2) ** 12 / weight**4
    rows = [mp.mpf(index) / (size - 1) for index in range(size)]

    for center in map(mp.mpf, ("-8", "-2", "3", "9")):
        columns = [center - 1 + 2 * mp.mpf(index) / (size - 1) for index in range(size)]

        for row in rows:
            for column in columns:
                shifted = center + 1 - column
                product = row * shifted
                low = mp.fsum(
                    product**degree / mp.factorial(degree) for degree in range(8)
                )
                tail = mp.exp(product) - low
                amplified = (
                    mp.exp(-row * (center + 1))
                    * (low + weight * tail)
                    / (1 + mp.exp(shifted - (center + 1)))
                )
                assert 0 < amplified <= 2

        original = mp.matrix(
            [
                [mp.exp(-row * column) / (1 + mp.exp(-column)) for column in columns]
                for row in rows
            ]
        )
        residual = mp.matrix(original)
        selected_rows: list[int] = []
        selected_columns: list[int] = []

        for _ in range(12):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[row, column]),
                    -row,
                    -column,
                    row,
                    column,
                )
                for row in range(size)
                for column in range(size)
            )
            pivot = residual[pivot_row, pivot_column]
            assert pivot != 0
            selected_rows.append(pivot_row)
            selected_columns.append(pivot_column)
            pivot_column_values = [residual[row, pivot_column] for row in range(size)]
            pivot_row_values = [residual[pivot_row, column] for column in range(size)]
            residual = mp.matrix(
                [
                    [
                        residual[row, column]
                        - pivot_column_values[row] * pivot_row_values[column] / pivot
                        for column in range(size)
                    ]
                    for row in range(size)
                ]
            )

        selected_core = mp.matrix(
            [
                [original[row, column] for column in selected_columns]
                for row in selected_rows
            ]
        )
        sorted_core = mp.matrix(
            [
                [original[row, column] for column in sorted(selected_columns)]
                for row in sorted(selected_rows)
            ]
        )
        assert mp.almosteq(abs(mp.det(selected_core)), mp.det(sorted_core))
        assert abs(mp.det(selected_core)) <= determinant_bound
        residual_maximum = max(
            abs(residual[row, column]) for row in range(size) for column in range(size)
        )
        assert residual_maximum <= mp.mpf("0.5")


def test_residual_amplification_bridge_and_frozen_column_obstruction() -> None:
    rows = (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1))
    columns = (Fraction(1, 4), Fraction(1, 2), Fraction(1), Fraction(2))

    def finite_exp(x: Fraction, terms: int = 6) -> Fraction:
        return sum(
            (x**degree / math.factorial(degree) for degree in range(terms)),
            Fraction(),
        )

    full = [[finite_exp(row * column) for column in columns] for row in rows]
    frozen_tail = [
        [Fraction() if j == 3 else full[i][j] - 1 for j in range(4)] for i in range(4)
    ]
    frozen_low = [[full[i][j] - frozen_tail[i][j] for j in range(4)] for i in range(4)]
    amplified_frozen = [
        [128 * frozen_tail[i][j] + frozen_low[i][j] for j in range(4)] for i in range(4)
    ]
    original_det = fraction_determinant(full)
    amplified_frozen_det = fraction_determinant(amplified_frozen)
    assert original_det == Fraction(91063777, 85899345920)
    assert amplified_frozen_det == -Fraction(32176499, 983040)
    assert amplified_frozen_det < 128 ** (4 - 2) * original_det

    # A one-pivot Schur residual of a finite positive-feature kernel has the
    # compatible split that the formal bridge requires.  The selected point
    # is x=y=0, so the exact Schur residual simply removes the constant feature.
    fresh = (Fraction(1, 4), Fraction(1, 2), Fraction(1))
    residual = [[finite_exp(row * column) - 1 for column in fresh] for row in fresh]
    low = [[row * column for column in fresh] for row in fresh]
    tail = [[residual[i][j] - low[i][j] for j in range(3)] for i in range(3)]

    coefficients: list[Fraction] = []
    for degree in range(4):
        coefficient = Fraction()
        for selected_tuple in itertools.combinations(range(3), degree):
            selected = set(selected_tuple)
            coefficient += fraction_determinant(
                [
                    [tail[i][j] if j in selected else low[i][j] for j in range(3)]
                    for i in range(3)
                ]
            )
        coefficients.append(coefficient)
    assert coefficients[0] == coefficients[1] == 0
    assert all(coefficient >= 0 for coefficient in coefficients)

    weight = 128
    amplified_residual = [
        [weight * tail[i][j] + low[i][j] for j in range(3)] for i in range(3)
    ]
    assert weight ** (3 - 1) * fraction_determinant(residual) <= fraction_determinant(
        amplified_residual
    )

    augmented_nodes = (Fraction(0), *fresh)
    augmented = [
        [finite_exp(row * column) for column in augmented_nodes]
        for row in augmented_nodes
    ]
    assert fraction_determinant(augmented) == fraction_determinant(residual)


def test_schur_feature_elimination_preserves_right_feature_rank_exactly() -> None:
    nodes = tuple(Fraction(index, 4) for index in range(5))
    columns = tuple(Fraction(index, 5) for index in range(1, 6))

    def finite_exp(x: Fraction) -> Fraction:
        return sum(
            (x**degree / math.factorial(degree) for degree in range(6)),
            Fraction(),
        )

    original = [[finite_exp(row * column) for column in columns] for row in nodes]
    left = [[Fraction(1), row] for row in nodes]
    right = [[Fraction(1) for _ in columns], list(columns)]
    tail = [
        [
            original[i][j]
            - sum(left[i][feature] * right[feature][j] for feature in range(2))
            for j in range(5)
        ]
        for i in range(5)
    ]

    def row_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - pivot_matrix[i][pivot_column] * companion[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    residual = original
    conditioned_tail = tail
    conditioned_left = left
    for step in range(3):
        _, _, _, pivot_row, pivot_column = max(
            (
                abs(residual[i][j]),
                -i,
                -j,
                i,
                j,
            )
            for i in range(5)
            for j in range(5)
        )
        next_residual = row_condition(residual, residual, pivot_row, pivot_column)
        conditioned_tail = row_condition(
            residual, conditioned_tail, pivot_row, pivot_column
        )
        conditioned_left = row_condition(
            residual, conditioned_left, pivot_row, pivot_column
        )
        residual = next_residual

        transported = [
            [
                conditioned_tail[i][j]
                + sum(
                    conditioned_left[i][feature] * right[feature][j]
                    for feature in range(2)
                )
                for j in range(5)
            ]
            for i in range(5)
        ]
        assert transported == residual

        if step == 0:
            fresh_rows = [i for i in range(5) if i != pivot_row]
            fresh_columns = [j for j in range(5) if j != pivot_column]
            conditioned_low_now = [
                [
                    sum(
                        conditioned_left[i][feature] * right[feature][j]
                        for feature in range(2)
                    )
                    for j in fresh_columns
                ]
                for i in fresh_rows
            ]
            conditioned_tail_now = [
                [conditioned_tail[i][j] for j in fresh_columns] for i in fresh_rows
            ]
            mixed_coefficients: list[Fraction] = []
            for degree in range(5):
                coefficient = Fraction()
                for tail_columns_tuple in itertools.combinations(range(4), degree):
                    tail_columns = set(tail_columns_tuple)
                    coefficient += fraction_determinant(
                        [
                            [
                                conditioned_tail_now[i][j]
                                if j in tail_columns
                                else conditioned_low_now[i][j]
                                for j in range(4)
                            ]
                            for i in range(4)
                        ]
                    )
                mixed_coefficients.append(coefficient)
            assert mixed_coefficients == [
                Fraction(),
                Fraction(),
                Fraction(777, 3338240000000),
                -Fraction(846159, 4172800000000000),
                Fraction(81, 260800000000000000),
            ]

    # The prefix may be arbitrarily long, but the transported low block still
    # factors through the original two features, so every 3 x 3 minor vanishes.
    conditioned_low = [
        [
            sum(
                conditioned_left[i][feature] * right[feature][j] for feature in range(2)
            )
            for j in range(5)
        ]
        for i in range(5)
    ]
    for selected_rows in itertools.combinations(range(5), 3):
        for selected_columns in itertools.combinations(range(5), 3):
            minor = [
                [conditioned_low[i][j] for j in selected_columns] for i in selected_rows
            ]
            assert fraction_determinant(minor) == 0


def test_finite_weight_norming_replaces_conditioned_sign_and_root_assumptions() -> None:
    nodes = tuple(56 + 9 * index for index in range(9))
    denominator = 387420489
    weights = tuple(
        Fraction(numerator, denominator)
        for numerator in (
            1271136797320,
            -8739065481575,
            26815762573600,
            -47745138240800,
            53778864502000,
            -39151013357456,
            17959180439200,
            -4739832125600,
            550492313800,
        )
    )
    for degree in range(9):
        assert (
            sum(
                (
                    weight * node**degree
                    for weight, node in zip(weights, nodes, strict=True)
                ),
                Fraction(),
            )
            == 1
        )
    assert sum(map(abs, weights), Fraction()) < 518173
    assert Fraction(518173) * 4**18 < Fraction(56) ** 10

    def finite_exp(x: Fraction, terms: int = 9) -> Fraction:
        return sum(
            (x**degree / math.factorial(degree) for degree in range(terms)),
            Fraction(),
        )

    def multiply(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(right))),
                    Fraction(),
                )
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    def row_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - pivot_matrix[i][pivot_column] * companion[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def column_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - companion[i][pivot_column] * pivot_matrix[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def bi_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        return column_condition(
            pivot_matrix,
            row_condition(pivot_matrix, companion, pivot_row, pivot_column),
            pivot_row,
            pivot_column,
        )

    def conditioned_pencil_coefficients(
        row_numerators: tuple[int, ...],
        column_numerators: tuple[int, ...],
        grid_denominator: int,
        low_terms: int,
        prefix_length: int,
        fresh_rows: tuple[int, ...],
        fresh_columns: tuple[int, ...],
    ) -> list[Fraction]:
        rows = tuple(Fraction(value, grid_denominator) for value in row_numerators)
        columns = tuple(
            Fraction(value, grid_denominator) for value in column_numerators
        )
        size = len(rows)
        original = [[finite_exp(x * y) for y in columns] for x in rows]
        left = [
            [x**degree / math.factorial(degree) for degree in range(low_terms)]
            for x in rows
        ]
        right = [[y**degree for y in columns] for degree in range(low_terms)]
        low = multiply(left, right)
        tail = [[original[i][j] - low[i][j] for j in range(size)] for i in range(size)]
        residual = original
        conditioned_tail = tail
        conditioned_left = left
        selected_rows: set[int] = set()
        selected_columns: set[int] = set()
        for _ in range(prefix_length):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[i][j]),
                    -i,
                    -j,
                    i,
                    j,
                )
                for i in range(size)
                for j in range(size)
                if i not in selected_rows and j not in selected_columns
            )
            conditioned_tail = row_condition(
                residual, conditioned_tail, pivot_row, pivot_column
            )
            conditioned_left = row_condition(
                residual, conditioned_left, pivot_row, pivot_column
            )
            residual = row_condition(residual, residual, pivot_row, pivot_column)
            selected_rows.add(pivot_row)
            selected_columns.add(pivot_column)
        conditioned_low = multiply(conditioned_left, right)
        order = len(fresh_rows)
        coefficients: list[Fraction] = []
        for degree in range(order + 1):
            coefficient = Fraction()
            for selected_tuple in itertools.combinations(range(order), degree):
                selected = set(selected_tuple)
                coefficient += fraction_determinant(
                    [
                        [
                            conditioned_tail[fresh_rows[i]][fresh_columns[j]]
                            if j in selected
                            else conditioned_low[fresh_rows[i]][fresh_columns[j]]
                            for j in range(order)
                        ]
                        for i in range(order)
                    ]
                )
            coefficients.append(coefficient)
        return coefficients

    # A genuine conditioned positive-Taylor pencil can violate amplification
    # at the single weight 128.
    coefficients = conditioned_pencil_coefficients(
        (0, 3, 4, 6, 7, 10, 13),
        (1, 2, 3, 4, 6, 8, 13),
        13,
        1,
        1,
        (0, 3, 5),
        (0, 1, 3),
    )

    def evaluate(values: list[Fraction], value: int) -> Fraction:
        return sum(
            (coefficient * value**degree for degree, coefficient in enumerate(values)),
            Fraction(),
        )

    assert abs(evaluate(coefficients, 128)) < 128**2 * abs(evaluate(coefficients, 1))
    norming_values = [abs(evaluate(coefficients, node)) for node in nodes]
    assert 56**2 * abs(evaluate(coefficients, 1)) <= 518173 * max(norming_values)

    # Real-rootedness also fails: after removing the forced quadratic zero,
    # this exact quadratic has negative discriminant.
    nonreal_coefficients = conditioned_pencil_coefficients(
        (0, 2, 3, 4, 6, 7, 11),
        (3, 5, 6, 7, 8, 9, 11),
        11,
        2,
        2,
        (1, 2, 3, 4),
        (1, 2, 4, 5),
    )
    assert nonreal_coefficients[:2] == [0, 0]
    quadratic = nonreal_coefficients[2:]
    assert len(quadratic) == 3
    assert quadratic[1] ** 2 - 4 * quadratic[0] * quadratic[2] < 0


def test_two_sided_projection_and_plucker_factorization_exactly() -> None:
    def finite_exp(x: Fraction, terms: int = 9) -> Fraction:
        return sum(
            (x**degree / math.factorial(degree) for degree in range(terms)),
            Fraction(),
        )

    def multiply(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(right))),
                    Fraction(),
                )
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    def add(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [left[i][j] + right[i][j] for j in range(len(left[0]))]
            for i in range(len(left))
        ]

    def row_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - pivot_matrix[i][pivot_column] * companion[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def column_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - companion[i][pivot_column] * pivot_matrix[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def bi_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        return column_condition(
            pivot_matrix,
            row_condition(pivot_matrix, companion, pivot_row, pivot_column),
            pivot_row,
            pivot_column,
        )

    def mixed_coefficients(
        tail: list[list[Fraction]],
        low: list[list[Fraction]],
        rows: tuple[int, ...],
        columns: tuple[int, ...],
    ) -> list[Fraction]:
        order = len(rows)
        coefficients: list[Fraction] = []
        for degree in range(order + 1):
            coefficient = Fraction()
            for selected_tuple in itertools.combinations(range(order), degree):
                selected = set(selected_tuple)
                coefficient += fraction_determinant(
                    [
                        [
                            tail[rows[i]][columns[j]]
                            if j in selected
                            else low[rows[i]][columns[j]]
                            for j in range(order)
                        ]
                        for i in range(order)
                    ]
                )
            coefficients.append(coefficient)
        return coefficients

    def conditioned_data(
        row_values: tuple[int, ...],
        column_values: tuple[int, ...],
        denominator: int,
        low_terms: int,
        prefix_length: int,
        *,
        two_sided: bool,
    ) -> tuple[
        list[list[Fraction]],
        list[list[Fraction]],
        set[int],
        set[int],
    ]:
        rows = tuple(Fraction(value, denominator) for value in row_values)
        columns = tuple(Fraction(value, denominator) for value in column_values)
        size = len(rows)
        original = [[finite_exp(x * y) for y in columns] for x in rows]
        left = [
            [x**degree / math.factorial(degree) for degree in range(low_terms)]
            for x in rows
        ]
        right = [[y**degree for y in columns] for degree in range(low_terms)]
        low = multiply(left, right)
        tail = [[original[i][j] - low[i][j] for j in range(size)] for i in range(size)]
        direct_low = [line[:] for line in low]
        residual = [line[:] for line in original]
        selected_rows: set[int] = set()
        selected_columns: set[int] = set()
        for _ in range(prefix_length):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[i][j]),
                    -i,
                    -j,
                    i,
                    j,
                )
                for i in range(size)
                for j in range(size)
                if i not in selected_rows and j not in selected_columns
            )
            if two_sided:
                tail = bi_condition(residual, tail, pivot_row, pivot_column)
                direct_low = bi_condition(residual, direct_low, pivot_row, pivot_column)
                left = row_condition(residual, left, pivot_row, pivot_column)
                right = column_condition(residual, right, pivot_row, pivot_column)
                assert direct_low == multiply(left, right)
            else:
                tail = row_condition(residual, tail, pivot_row, pivot_column)
                direct_low = row_condition(
                    residual, direct_low, pivot_row, pivot_column
                )
            residual = row_condition(residual, residual, pivot_row, pivot_column)
            selected_rows.add(pivot_row)
            selected_columns.add(pivot_column)
            assert residual == add(tail, direct_low)
            if two_sided:
                assert all(not tail[pivot_row][j] for j in range(size))
                assert all(not tail[i][pivot_column] for i in range(size))
                assert all(not direct_low[pivot_row][j] for j in range(size))
                assert all(not direct_low[i][pivot_column] for i in range(size))
        return tail, direct_low, selected_rows, selected_columns

    witness = ((1, 2, 6, 8, 9, 10, 11), (2, 3, 4, 10, 15, 19, 20))
    one_tail, one_low, _, _ = conditioned_data(*witness, 11, 2, 2, two_sided=False)
    bi_tail, bi_low, _, _ = conditioned_data(*witness, 11, 2, 2, two_sided=True)
    sample = (1, 2, 3, 4)
    one_signs = [
        0 if value == 0 else 1 if value > 0 else -1
        for value in mixed_coefficients(one_tail, one_low, sample, sample)
    ]
    bi_signs = [
        0 if value == 0 else 1 if value > 0 else -1
        for value in mixed_coefficients(bi_tail, bi_low, sample, sample)
    ]
    assert one_signs == [0, 0, 1, -1, 1]
    assert bi_signs == [0, 0, 1, 1, 1]

    rng = random.Random(0xB17B1)
    checked_samples = 0
    for _ in range(36):
        size = 6
        denominator = rng.choice((7, 9, 11, 13))
        row_values = tuple(sorted(rng.sample(range(denominator + 1), size)))
        column_values = tuple(sorted(rng.sample(range(1, 2 * denominator), size)))
        low_terms = rng.choice((1, 2, 3))
        prefix_length = rng.choice((1, 2))
        tail, low, selected_rows, selected_columns = conditioned_data(
            row_values,
            column_values,
            denominator,
            low_terms,
            prefix_length,
            two_sided=True,
        )
        fresh_rows = tuple(index for index in range(size) if index not in selected_rows)
        fresh_columns = tuple(
            index for index in range(size) if index not in selected_columns
        )
        for rows in itertools.combinations(fresh_rows, 3):
            for columns in itertools.combinations(fresh_columns, 3):
                coefficients = mixed_coefficients(tail, low, rows, columns)
                signs = {value > 0 for value in coefficients if value}
                assert len(signs) <= 1
                checked_samples += 1
    assert checked_samples == 2760

    def check_plucker_factorizations(
        row_values: tuple[int, ...],
        column_values: tuple[int, ...],
        denominator: int,
        low_terms: int,
        prefix_length: int,
    ) -> int:
        rows = tuple(Fraction(value, denominator) for value in row_values)
        columns = tuple(Fraction(value, denominator) for value in column_values)
        size = len(rows)
        original = [[finite_exp(x * y) for y in columns] for x in rows]
        original_left = [
            [x**degree / math.factorial(degree) for degree in range(low_terms)]
            for x in rows
        ]
        original_right = [[y**degree for y in columns] for degree in range(low_terms)]
        conditioned_left = [line[:] for line in original_left]
        conditioned_right = [line[:] for line in original_right]
        residual = [line[:] for line in original]
        selected_rows: list[int] = []
        selected_columns: list[int] = []
        for _ in range(prefix_length):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[i][j]),
                    -i,
                    -j,
                    i,
                    j,
                )
                for i in range(size)
                for j in range(size)
                if i not in selected_rows and j not in selected_columns
            )
            conditioned_left = row_condition(
                residual, conditioned_left, pivot_row, pivot_column
            )
            conditioned_right = column_condition(
                residual, conditioned_right, pivot_row, pivot_column
            )
            residual = row_condition(residual, residual, pivot_row, pivot_column)
            selected_rows.append(pivot_row)
            selected_columns.append(pivot_column)

        core_det = fraction_determinant(
            [[original[i][j] for j in selected_columns] for i in selected_rows]
        )
        assert core_det
        fresh_rows = tuple(index for index in range(size) if index not in selected_rows)
        fresh_columns = tuple(
            index for index in range(size) if index not in selected_columns
        )
        minor_order = min(2, low_terms)
        checked = 0
        for sampled_rows in itertools.combinations(fresh_rows, minor_order):
            augmented_rows = (*selected_rows, *sampled_rows)
            for sampled_columns in itertools.combinations(fresh_columns, minor_order):
                augmented_columns = (*selected_columns, *sampled_columns)
                for features in itertools.combinations(range(low_terms), minor_order):
                    left_core = [
                        [
                            original[i][selected_columns[j]]
                            if j < prefix_length
                            else original_left[i][features[j - prefix_length]]
                            for j in range(prefix_length + minor_order)
                        ]
                        for i in augmented_rows
                    ]
                    right_core = [
                        [
                            original[selected_rows[i]][j]
                            if i < prefix_length
                            else original_right[features[i - prefix_length]][j]
                            for j in augmented_columns
                        ]
                        for i in range(prefix_length + minor_order)
                    ]
                    left_minor = [
                        [conditioned_left[i][degree] for degree in features]
                        for i in sampled_rows
                    ]
                    right_minor = [
                        [conditioned_right[degree][j] for j in sampled_columns]
                        for degree in features
                    ]
                    left_det = fraction_determinant(left_core)
                    right_det = fraction_determinant(right_core)
                    conditioned_left_det = fraction_determinant(left_minor)
                    conditioned_right_det = fraction_determinant(right_minor)
                    assert left_det == core_det * conditioned_left_det
                    assert right_det == core_det * conditioned_right_det
                    assert left_det * right_det == core_det**2 * (
                        conditioned_left_det * conditioned_right_det
                    )
                    checked += 1
        return checked

    checked_plucker = check_plucker_factorizations(*witness, 11, 2, 2)
    plucker_rng = random.Random(0xB1C0)
    for _ in range(18):
        denominator = plucker_rng.choice((7, 9, 11, 13))
        checked_plucker += check_plucker_factorizations(
            tuple(sorted(plucker_rng.sample(range(denominator + 1), 6))),
            tuple(sorted(plucker_rng.sample(range(1, 2 * denominator), 6))),
            denominator,
            plucker_rng.choice((1, 2, 3)),
            plucker_rng.choice((1, 2)),
        )
    assert checked_plucker == 2382


def test_graded_plucker_pairing_survives_termwise_sign_obstruction() -> None:
    def finite_exp(value: Fraction, terms: int = 14) -> Fraction:
        return sum(
            (value**degree / math.factorial(degree) for degree in range(terms)),
            Fraction(),
        )

    def multiply(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(right))),
                    Fraction(),
                )
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    def solve(
        coefficient: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        order = len(coefficient)
        columns = len(right[0])
        upper = [row[:] for row in coefficient]
        transformed = [row[:] for row in right]
        for pivot_index in range(order):
            pivot = upper[pivot_index][pivot_index]
            assert pivot
            for i in range(pivot_index + 1, order):
                factor = upper[i][pivot_index] / pivot
                for j in range(pivot_index, order):
                    upper[i][j] -= factor * upper[pivot_index][j]
                for j in range(columns):
                    transformed[i][j] -= factor * transformed[pivot_index][j]
        solution = [[Fraction() for _ in range(columns)] for _ in range(order)]
        for i in reversed(range(order)):
            for j in range(columns):
                solution[i][j] = (
                    transformed[i][j]
                    - sum(
                        (upper[i][k] * solution[k][j] for k in range(i + 1, order)),
                        Fraction(),
                    )
                ) / upper[i][i]
        return solution

    def row_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - pivot_matrix[i][pivot_column] * companion[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def column_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - companion[i][pivot_column] * pivot_matrix[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def bi_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        return column_condition(
            pivot_matrix,
            row_condition(pivot_matrix, companion, pivot_row, pivot_column),
            pivot_row,
            pivot_column,
        )

    def mixed_coefficient(
        high: list[list[Fraction]],
        low: list[list[Fraction]],
        degree: int,
    ) -> Fraction:
        order = len(high)
        return sum(
            (
                fraction_determinant(
                    [
                        [
                            high[i][j] if j in selected else low[i][j]
                            for j in range(order)
                        ]
                        for i in range(order)
                    ]
                )
                for selected_tuple in itertools.combinations(range(order), degree)
                for selected in [set(selected_tuple)]
            ),
            Fraction(),
        )

    denominator = 17
    row_values = (3, 5, 6, 12, 13, 15, 17)
    column_values = (1, 7, 11, 12, 13, 16, 17)
    rows = tuple(Fraction(value, denominator) for value in row_values)
    columns = tuple(Fraction(value, denominator) for value in column_values)
    size = len(rows)
    feature_count = 8
    prefix_length = 3

    original = [[finite_exp(x * y) for y in columns] for x in rows]
    left = [
        [x**degree / math.factorial(degree) for degree in range(feature_count)]
        for x in rows
    ]
    right = [[y**degree for y in columns] for degree in range(feature_count)]
    residual = [line[:] for line in original]
    selected_rows: list[int] = []
    selected_columns: list[int] = []
    for _ in range(prefix_length):
        _, _, _, pivot_row, pivot_column = max(
            (
                abs(residual[i][j]),
                -i,
                -j,
                i,
                j,
            )
            for i in range(size)
            for j in range(size)
            if i not in selected_rows and j not in selected_columns
        )
        residual = row_condition(residual, residual, pivot_row, pivot_column)
        selected_rows.append(pivot_row)
        selected_columns.append(pivot_column)

    assert selected_rows == [6, 0, 3]
    assert selected_columns == [6, 0, 2]
    fresh_rows = (1, 2, 4)
    fresh_columns = (3, 4, 5)
    augmented_rows = (*selected_rows, *fresh_rows)
    augmented_columns = (*selected_columns, *fresh_columns)

    products: list[Fraction] = []
    for features in itertools.combinations(range(feature_count), len(fresh_rows)):
        left_core = [
            [
                original[i][selected_columns[j]]
                if j < prefix_length
                else left[i][features[j - prefix_length]]
                for j in range(prefix_length + len(features))
            ]
            for i in augmented_rows
        ]
        right_core = [
            [
                original[selected_rows[i]][j]
                if i < prefix_length
                else right[features[i - prefix_length]][j]
                for j in augmented_columns
            ]
            for i in range(prefix_length + len(features))
        ]
        products.append(
            fraction_determinant(left_core) * fraction_determinant(right_core)
        )

    assert len(products) == math.comb(feature_count, len(fresh_rows)) == 56
    assert sum(product < 0 for product in products) == 2
    assert sum(product > 0 for product in products) == 54

    fixed_left = [[original[i][j] for j in selected_columns] for i in augmented_rows]
    fixed_right = [[original[i][j] for j in augmented_columns] for i in selected_rows]
    feature_left = [
        [left[i][degree] for degree in range(feature_count)] for i in augmented_rows
    ]
    feature_right = [
        [right[degree][j] for j in augmented_columns] for degree in range(feature_count)
    ]
    fixed_product = multiply(fixed_left, fixed_right)
    feature_product = multiply(feature_left, feature_right)
    graded_pairing = sum(products, Fraction())

    assert graded_pairing == mixed_coefficient(
        fixed_product, feature_product, prefix_length
    )
    assert graded_pairing > 0

    # Linearize every finite weight as one bordered determinant.  The signed
    # four-cycle also rules out repairing these matrices by a fixed row/column
    # sign gauge into an entrywise-nonnegative (hence totally nonnegative)
    # matrix.
    finite_weights = (1, 56, 65, 74, 83, 92, 101, 110, 119, 128)
    order = len(augmented_rows)
    global_low = multiply(left, right)
    global_tail = [
        [original[i][j] - global_low[i][j] for j in range(size)] for i in range(size)
    ]
    core = [[original[i][j] for j in selected_columns] for i in selected_rows]
    core_det = fraction_determinant(core)
    assert core_det
    for weight in finite_weights:
        rank_update = [
            [feature_product[i][j] + weight * fixed_product[i][j] for j in range(order)]
            for i in range(order)
        ]
        bordered = [
            [
                Fraction(i == j)
                if j < prefix_length
                else fixed_right[i][j - prefix_length]
                for j in range(prefix_length + order)
            ]
            for i in range(prefix_length)
        ] + [
            [
                -weight * fixed_left[i][j]
                if j < prefix_length
                else feature_product[i][j - prefix_length]
                for j in range(prefix_length + order)
            ]
            for i in range(order)
        ]
        assert fraction_determinant(bordered) == fraction_determinant(rank_update)
        sign_cycle = (
            bordered[0][0]
            * bordered[0][prefix_length]
            * bordered[prefix_length][0]
            * bordered[prefix_length][prefix_length]
        )
        assert sign_cycle < 0

        companion = [
            [global_low[i][j] + weight * global_tail[i][j] for j in range(size)]
            for i in range(size)
        ]
        pivot_kernel = [row[:] for row in original]
        conditioned = [row[:] for row in companion]
        for pivot_row, pivot_column in zip(
            selected_rows, selected_columns, strict=True
        ):
            conditioned = bi_condition(
                pivot_kernel, conditioned, pivot_row, pivot_column
            )
            pivot_kernel = row_condition(
                pivot_kernel, pivot_kernel, pivot_row, pivot_column
            )

        kernel_left = [[original[i][j] for j in selected_columns] for i in fresh_rows]
        kernel_right = [[original[i][j] for j in fresh_columns] for i in selected_rows]
        companion_core = [
            [companion[i][j] for j in selected_columns] for i in selected_rows
        ]
        companion_left = [
            [companion[i][j] for j in selected_columns] for i in fresh_rows
        ]
        companion_right = [
            [companion[i][j] for j in fresh_columns] for i in selected_rows
        ]
        companion_fresh = [[companion[i][j] for j in fresh_columns] for i in fresh_rows]
        core_solve_companion_right = solve(core, companion_right)
        core_solve_kernel_right = solve(core, kernel_right)
        core_solve_mixed = solve(
            core, multiply(companion_core, core_solve_kernel_right)
        )
        schur = [
            [
                companion_fresh[i][j]
                - multiply(kernel_left, core_solve_companion_right)[i][j]
                - multiply(companion_left, core_solve_kernel_right)[i][j]
                + multiply(kernel_left, core_solve_mixed)[i][j]
                for j in range(len(fresh_columns))
            ]
            for i in range(len(fresh_rows))
        ]
        conditioned_fresh = [
            [conditioned[i][j] for j in fresh_columns] for i in fresh_rows
        ]
        assert schur == conditioned_fresh

        zero_core = [[Fraction() for _ in selected_rows] for _ in selected_rows]
        double_bordered = (
            [
                [*core[i], *companion_core[i], *companion_right[i]]
                for i in range(prefix_length)
            ]
            + [
                [*zero_core[i], *core[i], *kernel_right[i]]
                for i in range(prefix_length)
            ]
            + [
                [*kernel_left[i], *companion_left[i], *companion_fresh[i]]
                for i in range(len(fresh_rows))
            ]
        )
        assert fraction_determinant(double_bordered) == core_det**2 * (
            fraction_determinant(conditioned_fresh)
        )
        assert any(
            fraction_determinant(
                [[double_bordered[i][j] for j in columns] for i in rows]
            )
            < 0
            for rows in itertools.combinations(range(len(double_bordered)), 2)
            for columns in itertools.combinations(range(len(double_bordered)), 2)
        )

    # Exercise the arbitrary-run closed form independently of the fixed
    # finite-Taylor witness.  Cauchy kernels give exact nonsingular pivot
    # cores, while deterministic random companions, weights, prefix lengths,
    # and fresh blocks probe the full rectangular projection before the
    # square determinant specialization.  Every core action is evaluated by
    # exact Gaussian solves rather than by forming an inverse.
    closed_form_rng = random.Random(0xC105EDF0)
    checked_closed_forms = 0
    for _ in range(12):
        random_rows = tuple(sorted(closed_form_rng.sample(range(1, 45), 6)))
        random_columns = tuple(sorted(closed_form_rng.sample(range(47, 101), 6)))
        random_kernel = [
            [Fraction(1, random_rows[i] + random_columns[j]) for j in range(6)]
            for i in range(6)
        ]
        companion_base = [
            [
                Fraction(
                    closed_form_rng.randint(-11, 11),
                    closed_form_rng.randint(1, 9),
                )
                for _ in range(6)
            ]
            for _ in range(6)
        ]
        companion_direction = [
            [
                Fraction(closed_form_rng.randint(-9, 9), closed_form_rng.randint(1, 7))
                for _ in range(6)
            ]
            for _ in range(6)
        ]

        pivot_kernel = [line[:] for line in random_kernel]
        random_selected_rows: list[int] = []
        random_selected_columns: list[int] = []
        for _ in range(3):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(pivot_kernel[i][j]),
                    -i,
                    -j,
                    i,
                    j,
                )
                for i in range(6)
                for j in range(6)
                if i not in random_selected_rows and j not in random_selected_columns
            )
            assert pivot_kernel[pivot_row][pivot_column]
            random_selected_rows.append(pivot_row)
            random_selected_columns.append(pivot_column)
            pivot_kernel = row_condition(
                pivot_kernel, pivot_kernel, pivot_row, pivot_column
            )

        for prefix_length in (1, 2, 3):
            prefix_rows = random_selected_rows[:prefix_length]
            prefix_columns = random_selected_columns[:prefix_length]
            available_rows = [i for i in range(6) if i not in prefix_rows]
            available_columns = [j for j in range(6) if j not in prefix_columns]
            fresh_order = closed_form_rng.randint(1, min(3, len(available_rows)))
            sampled_rows = tuple(
                sorted(closed_form_rng.sample(available_rows, fresh_order))
            )
            sampled_columns = tuple(
                sorted(closed_form_rng.sample(available_columns, fresh_order))
            )
            core = [[random_kernel[i][j] for j in prefix_columns] for i in prefix_rows]
            core_det = fraction_determinant(core)
            assert core_det

            for weight in (-3, 0, 2, 7):
                random_companion = [
                    [
                        companion_base[i][j] + weight * companion_direction[i][j]
                        for j in range(6)
                    ]
                    for i in range(6)
                ]
                recursive = [line[:] for line in random_companion]
                recursive_kernel = [line[:] for line in random_kernel]
                for pivot_row, pivot_column in zip(
                    prefix_rows, prefix_columns, strict=True
                ):
                    recursive = bi_condition(
                        recursive_kernel, recursive, pivot_row, pivot_column
                    )
                    recursive_kernel = row_condition(
                        recursive_kernel,
                        recursive_kernel,
                        pivot_row,
                        pivot_column,
                    )

                kernel_left = [
                    [random_kernel[i][j] for j in prefix_columns] for i in sampled_rows
                ]
                kernel_right = [
                    [random_kernel[i][j] for j in sampled_columns] for i in prefix_rows
                ]
                companion_core = [
                    [random_companion[i][j] for j in prefix_columns]
                    for i in prefix_rows
                ]
                companion_left = [
                    [random_companion[i][j] for j in prefix_columns]
                    for i in sampled_rows
                ]
                companion_right = [
                    [random_companion[i][j] for j in sampled_columns]
                    for i in prefix_rows
                ]
                companion_fresh = [
                    [random_companion[i][j] for j in sampled_columns]
                    for i in sampled_rows
                ]
                left_right_correction = multiply(
                    kernel_left, solve(core, companion_right)
                )
                right_left_correction = multiply(
                    companion_left, solve(core, kernel_right)
                )
                mixed_correction = multiply(
                    kernel_left,
                    solve(core, multiply(companion_core, solve(core, kernel_right))),
                )
                closed = [
                    [
                        companion_fresh[i][j]
                        - left_right_correction[i][j]
                        - right_left_correction[i][j]
                        + mixed_correction[i][j]
                        for j in range(fresh_order)
                    ]
                    for i in range(fresh_order)
                ]
                recursive_fresh = [
                    [recursive[i][j] for j in sampled_columns] for i in sampled_rows
                ]
                assert recursive_fresh == closed

                zero_core = [
                    [Fraction() for _ in range(prefix_length)]
                    for _ in range(prefix_length)
                ]
                double_bordered = (
                    [
                        [*core[i], *companion_core[i], *companion_right[i]]
                        for i in range(prefix_length)
                    ]
                    + [
                        [*zero_core[i], *core[i], *kernel_right[i]]
                        for i in range(prefix_length)
                    ]
                    + [
                        [*kernel_left[i], *companion_left[i], *companion_fresh[i]]
                        for i in range(fresh_order)
                    ]
                )
                assert fraction_determinant(double_bordered) == core_det**2 * (
                    fraction_determinant(recursive_fresh)
                )
                checked_closed_forms += 1

    assert checked_closed_forms == 144


def test_upper_tangent_lift_and_oriented_coefficients_exactly() -> None:
    def multiply(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(right))),
                    Fraction(),
                )
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    def upper_tangent(
        base: list[list[Fraction]], companion: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        row_count = len(base)
        column_count = len(base[0])
        return [
            [
                base[i][j]
                if row_layer == column_layer
                else companion[i][j]
                if row_layer == 0 and column_layer == 1
                else Fraction()
                for column_layer in range(2)
                for j in range(column_count)
            ]
            for row_layer in range(2)
            for i in range(row_count)
        ]

    def block_diagonal(
        matrix: list[list[Fraction]],
    ) -> list[list[Fraction]]:
        row_count = len(matrix)
        column_count = len(matrix[0])
        return [
            [
                matrix[i][j] if row_layer == column_layer else Fraction()
                for column_layer in range(2)
                for j in range(column_count)
            ]
            for row_layer in range(2)
            for i in range(row_count)
        ]

    def solve_vector(
        coefficient: list[list[Fraction]], right: list[Fraction]
    ) -> list[Fraction]:
        order = len(coefficient)
        upper = [line[:] for line in coefficient]
        transformed = right[:]
        for pivot_index in range(order):
            pivot_row = next(
                i for i in range(pivot_index, order) if upper[i][pivot_index]
            )
            if pivot_row != pivot_index:
                upper[pivot_index], upper[pivot_row] = (
                    upper[pivot_row],
                    upper[pivot_index],
                )
                transformed[pivot_index], transformed[pivot_row] = (
                    transformed[pivot_row],
                    transformed[pivot_index],
                )
            pivot = upper[pivot_index][pivot_index]
            for i in range(pivot_index + 1, order):
                factor = upper[i][pivot_index] / pivot
                for j in range(pivot_index, order):
                    upper[i][j] -= factor * upper[pivot_index][j]
                transformed[i] -= factor * transformed[pivot_index]
        solution = [Fraction() for _ in range(order)]
        for i in reversed(range(order)):
            solution[i] = (
                transformed[i]
                - sum(
                    (upper[i][j] * solution[j] for j in range(i + 1, order)),
                    Fraction(),
                )
            ) / upper[i][i]
        return solution

    def row_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - pivot_matrix[i][pivot_column] * companion[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def column_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_matrix[pivot_row][pivot_column]
        assert pivot
        return [
            [
                companion[i][j]
                - companion[i][pivot_column] * pivot_matrix[pivot_row][j] / pivot
                for j in range(len(companion[0]))
            ]
            for i in range(len(companion))
        ]

    def bi_condition(
        pivot_matrix: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        return column_condition(
            pivot_matrix,
            row_condition(pivot_matrix, companion, pivot_row, pivot_column),
            pivot_row,
            pivot_column,
        )

    tangent_rng = random.Random(0x7A6E17)
    checked_factorizations = 0
    checked_triple_expansions = 0
    for _ in range(24):
        row_count = tangent_rng.randint(2, 4)
        feature_count = tangent_rng.randint(1, 4)
        column_count = tangent_rng.randint(2, 4)
        left = [
            [
                Fraction(tangent_rng.randint(-5, 7), tangent_rng.randint(1, 5))
                for _ in range(feature_count)
            ]
            for _ in range(row_count)
        ]
        weight = [
            [
                Fraction(tangent_rng.randint(-4, 6), tangent_rng.randint(1, 5))
                for _ in range(feature_count)
            ]
            for _ in range(feature_count)
        ]
        right = [
            [
                Fraction(tangent_rng.randint(-6, 8), tangent_rng.randint(1, 5))
                for _ in range(column_count)
            ]
            for _ in range(feature_count)
        ]
        base = multiply(left, right)
        companion = multiply(multiply(left, weight), right)
        identity = [
            [Fraction(i == j) for j in range(feature_count)]
            for i in range(feature_count)
        ]
        middle = upper_tangent(identity, weight)
        factored = multiply(
            multiply(block_diagonal(left), middle), block_diagonal(right)
        )
        lifted = upper_tangent(base, companion)
        assert factored == lifted

        # Exercise the exact ordered three-factor Cauchy--Binet expansion on
        # arbitrary (usually nonmonotone) outer samples.  Only the two feature
        # flags are enumerated in increasing order.
        sample_order = tangent_rng.randint(
            1,
            min(3, 2 * row_count, 2 * feature_count, 2 * column_count),
        )
        sampled_lifted_rows = tuple(
            tangent_rng.sample(range(2 * row_count), sample_order)
        )
        sampled_lifted_columns = tuple(
            tangent_rng.sample(range(2 * column_count), sample_order)
        )
        lifted_sample_det = fraction_determinant(
            [
                [lifted[i][j] for j in sampled_lifted_columns]
                for i in sampled_lifted_rows
            ]
        )
        lifted_left = block_diagonal(left)
        lifted_right = block_diagonal(right)
        triple_sum = sum(
            (
                fraction_determinant(
                    [
                        [lifted_left[i][j] for j in left_flag]
                        for i in sampled_lifted_rows
                    ]
                )
                * fraction_determinant(
                    [[middle[i][j] for j in right_flag] for i in left_flag]
                )
                * fraction_determinant(
                    [
                        [lifted_right[i][j] for j in sampled_lifted_columns]
                        for i in right_flag
                    ]
                )
                for left_flag in itertools.combinations(
                    range(2 * feature_count), sample_order
                )
                for right_flag in itertools.combinations(
                    range(2 * feature_count), sample_order
                )
            ),
            Fraction(),
        )
        assert triple_sum == lifted_sample_det
        checked_triple_expansions += 1

        selected_row = tangent_rng.randrange(row_count)
        selected_column = tangent_rng.randrange(column_count)
        fresh_row = next(i for i in range(row_count) if i != selected_row)
        fresh_column = next(j for j in range(column_count) if j != selected_column)
        lifted_rows = (selected_row, row_count + selected_row, fresh_row)
        lifted_columns = (
            selected_column,
            column_count + selected_column,
            column_count + fresh_column,
        )
        distinguished = [[lifted[i][j] for j in lifted_columns] for i in lifted_rows]
        double_bordered = [
            [
                base[selected_row][selected_column],
                companion[selected_row][selected_column],
                companion[selected_row][fresh_column],
            ],
            [
                Fraction(),
                base[selected_row][selected_column],
                base[selected_row][fresh_column],
            ],
            [
                base[fresh_row][selected_column],
                companion[fresh_row][selected_column],
                companion[fresh_row][fresh_column],
            ],
        ]
        assert distinguished == double_bordered
        checked_factorizations += 1
    assert checked_factorizations == 24

    # The full lift is not totally nonnegative, but the distinguished
    # complete-prefix minors exhibit a narrower coefficientwise orientation.
    # Interpolate their exact tail-weight polynomials and retain this as
    # evidence for a future matrix-valued total-positivity theorem.
    coefficient_rng = random.Random(0x7A6E)
    checked_coefficients = 0
    found_negative_ordinary_minor = False
    for _ in range(8):
        size = 6
        denominator = coefficient_rng.choice((11, 13, 17))
        row_nodes = tuple(
            sorted(
                Fraction(value, denominator)
                for value in coefficient_rng.sample(range(1, denominator + 1), size)
            )
        )
        column_nodes = tuple(
            sorted(
                Fraction(value, denominator)
                for value in coefficient_rng.sample(range(1, 2 * denominator), size)
            )
        )
        terms = 8
        low_terms = 3
        kernel = [
            [
                sum(
                    (
                        x**degree * y**degree / math.factorial(degree)
                        for degree in range(terms)
                    ),
                    Fraction(),
                )
                for y in column_nodes
            ]
            for x in row_nodes
        ]
        residual = [line[:] for line in kernel]
        selected_rows: list[int] = []
        selected_columns: list[int] = []
        for _ in range(3):
            _, _, _, pivot_row, pivot_column = max(
                (
                    abs(residual[i][j]),
                    -i,
                    -j,
                    i,
                    j,
                )
                for i in range(size)
                for j in range(size)
                if i not in selected_rows and j not in selected_columns
            )
            selected_rows.append(pivot_row)
            selected_columns.append(pivot_column)
            residual = row_condition(residual, residual, pivot_row, pivot_column)

        for prefix_length in (1, 2, 3):
            prefix_rows = selected_rows[:prefix_length]
            prefix_columns = selected_columns[:prefix_length]
            available_rows = [i for i in range(size) if i not in prefix_rows]
            available_columns = [j for j in range(size) if j not in prefix_columns]
            for order in (1, 2, 3):
                for _ in range(5):
                    sampled_rows = tuple(
                        sorted(coefficient_rng.sample(available_rows, order))
                    )
                    sampled_columns = tuple(
                        sorted(coefficient_rng.sample(available_columns, order))
                    )
                    values: list[Fraction] = []
                    for tail_weight in range(order + 1):
                        companion = [
                            [
                                sum(
                                    (
                                        (1 if degree < low_terms else tail_weight)
                                        * x**degree
                                        * y**degree
                                        / math.factorial(degree)
                                        for degree in range(terms)
                                    ),
                                    Fraction(),
                                )
                                for y in column_nodes
                            ]
                            for x in row_nodes
                        ]
                        pivot_kernel = [line[:] for line in kernel]
                        for pivot_row, pivot_column in zip(
                            prefix_rows, prefix_columns, strict=True
                        ):
                            companion = bi_condition(
                                pivot_kernel,
                                companion,
                                pivot_row,
                                pivot_column,
                            )
                            pivot_kernel = row_condition(
                                pivot_kernel,
                                pivot_kernel,
                                pivot_row,
                                pivot_column,
                            )
                        values.append(
                            fraction_determinant(
                                [
                                    [companion[i][j] for j in sampled_columns]
                                    for i in sampled_rows
                                ]
                            )
                        )
                    vandermonde = [
                        [Fraction(node) ** degree for degree in range(order + 1)]
                        for node in range(order + 1)
                    ]
                    coefficients = solve_vector(vandermonde, values)
                    cross_inversions = sum(
                        selected > fresh
                        for selected in prefix_rows
                        for fresh in sampled_rows
                    ) + sum(
                        selected > fresh
                        for selected in prefix_columns
                        for fresh in sampled_columns
                    )
                    orientation = -1 if cross_inversions % 2 else 1
                    assert all(
                        orientation * coefficient >= 0 for coefficient in coefficients
                    )
                    assert any(
                        orientation * coefficient > 0 for coefficient in coefficients
                    )
                    assert orientation * values[1] > 0
                    checked_coefficients += 1

        if not found_negative_ordinary_minor:
            weighted = [
                [
                    sum(
                        (
                            (1 if degree < low_terms else 7)
                            * x**degree
                            * y**degree
                            / math.factorial(degree)
                            for degree in range(terms)
                        ),
                        Fraction(),
                    )
                    for y in column_nodes
                ]
                for x in row_nodes
            ]
            lifted = upper_tangent(kernel, weighted)
            found_negative_ordinary_minor = any(
                fraction_determinant([[lifted[i][j] for j in columns] for i in rows])
                < 0
                for rows in itertools.combinations(range(2 * size), 2)
                for columns in itertools.combinations(range(2 * size), 2)
            )

    assert checked_coefficients == 360
    assert checked_triple_expansions == 24
    assert found_negative_ordinary_minor


def test_tangent_factorwise_orientation_obstruction_exactly() -> None:
    def multiply(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(right))),
                    Fraction(),
                )
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    def witness(
        feature_count: int,
    ) -> tuple[
        list[tuple[tuple[int, ...], tuple[int, ...], Fraction]],
        Fraction,
    ]:
        row_nodes = (Fraction(1, 3), Fraction(2, 3))
        column_nodes = (Fraction(1, 4), Fraction(3, 4))
        weights = tuple(Fraction(degree + 2) for degree in range(feature_count))
        left_base = [
            [node**degree / math.factorial(degree) for degree in range(feature_count)]
            for node in row_nodes
        ]
        right_base = [
            [node**degree for node in column_nodes] for degree in range(feature_count)
        ]
        zero_block = [Fraction() for _ in range(feature_count)]
        left = [
            left_base[0] + zero_block,
            zero_block + left_base[0],
            left_base[1] + zero_block,
        ]
        middle = [
            [
                Fraction(i == j)
                + (
                    weights[i]
                    if i < feature_count and j == feature_count + i
                    else Fraction()
                )
                for j in range(2 * feature_count)
            ]
            for i in range(2 * feature_count)
        ]
        right = [
            [right_base[i][0], Fraction(), Fraction()] for i in range(feature_count)
        ] + [
            [Fraction(), right_base[i][0], right_base[i][1]]
            for i in range(feature_count)
        ]

        terms: list[tuple[tuple[int, ...], tuple[int, ...], Fraction]] = []
        for left_flag in itertools.combinations(range(2 * feature_count), 3):
            left_minor = fraction_determinant(
                [[left[i][j] for j in left_flag] for i in range(3)]
            )
            for right_flag in itertools.combinations(range(2 * feature_count), 3):
                term = (
                    left_minor
                    * fraction_determinant(
                        [[middle[i][j] for j in right_flag] for i in left_flag]
                    )
                    * fraction_determinant(
                        [[right[i][j] for j in range(3)] for i in right_flag]
                    )
                )
                if term:
                    terms.append((left_flag, right_flag, term))

        complete_det = fraction_determinant(multiply(multiply(left, middle), right))
        assert sum((term for _, _, term in terms), Fraction()) == complete_det
        base = multiply(left_base, right_base)
        companion = multiply(
            [
                [left_base[i][k] * weights[k] for k in range(feature_count)]
                for i in range(len(left_base))
            ],
            right_base,
        )
        core = base[0][0]
        conditioned = (
            companion[1][1]
            - base[1][0] * companion[0][1] / core
            - companion[1][0] * base[0][1] / core
            + base[1][0] * companion[0][0] * base[0][1] / core**2
        )
        assert complete_det == core**2 * conditioned
        return terms, complete_det

    two_feature_terms, two_feature_det = witness(2)
    assert [term for _, _, term in two_feature_terms] == [
        Fraction(1, 2),
        Fraction(1, 36),
    ]
    assert two_feature_det == Fraction(19, 36)

    three_feature_terms, three_feature_det = witness(3)
    three_feature_values = [term for _, _, term in three_feature_terms]
    assert len(three_feature_values) == 12
    assert sum(term > 0 for term in three_feature_values) == 10
    assert sum(term < 0 for term in three_feature_values) == 2
    assert (
        (0, 1, 5),
        (0, 4, 5),
        Fraction(-1, 192),
    ) in three_feature_terms
    assert three_feature_det == Fraction(24929, 27648)
    assert three_feature_det > 0


def test_conditioned_tangent_nonnegative_weight_boundary_exactly() -> None:
    def bi_condition(
        pivot_kernel: list[list[Fraction]],
        companion: list[list[Fraction]],
        pivot_row: int,
        pivot_column: int,
    ) -> list[list[Fraction]]:
        pivot = pivot_kernel[pivot_row][pivot_column]
        return [
            [
                companion[i][j]
                - pivot_kernel[i][pivot_column] * companion[pivot_row][j] / pivot
                - companion[i][pivot_column] * pivot_kernel[pivot_row][j] / pivot
                + pivot_kernel[i][pivot_column]
                * companion[pivot_row][pivot_column]
                * pivot_kernel[pivot_row][j]
                / pivot**2
                for j in range(len(pivot_kernel[0]))
            ]
            for i in range(len(pivot_kernel))
        ]

    def conditioned_minor(
        degrees: tuple[int, ...],
        weights: tuple[Fraction, ...],
        row_nodes: tuple[Fraction, ...],
        column_nodes: tuple[Fraction, ...],
        selected_rows: tuple[int, ...],
        selected_columns: tuple[int, ...],
        fresh_rows: tuple[int, ...],
        fresh_columns: tuple[int, ...],
    ) -> Fraction:
        left = [
            [node**degree / math.factorial(degree) for degree in degrees]
            for node in row_nodes
        ]
        right = [[node**degree for node in column_nodes] for degree in degrees]
        kernel = [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(len(degrees))),
                    Fraction(),
                )
                for j in range(len(column_nodes))
            ]
            for i in range(len(row_nodes))
        ]
        companion = [
            [
                sum(
                    (
                        left[i][k] * weights[k] * right[k][j]
                        for k in range(len(degrees))
                    ),
                    Fraction(),
                )
                for j in range(len(column_nodes))
            ]
            for i in range(len(row_nodes))
        ]
        for pivot_row, pivot_column in zip(
            selected_rows, selected_columns, strict=True
        ):
            companion = bi_condition(kernel, companion, pivot_row, pivot_column)
            kernel = bi_condition(kernel, kernel, pivot_row, pivot_column)
        return fraction_determinant(
            [[companion[i][j] for j in fresh_columns] for i in fresh_rows]
        )

    rng = random.Random(0xCE11)
    positive_cases = 0
    for _ in range(5000):
        node_count = rng.randrange(3, 7)
        feature_count = rng.randrange(2, 8)
        degrees = tuple(sorted(rng.sample(range(15), feature_count)))
        row_denominator = rng.randrange(node_count + 1, 18)
        column_denominator = rng.randrange(node_count + 1, 18)
        row_nodes = tuple(
            Fraction(value, row_denominator)
            for value in sorted(rng.sample(range(1, row_denominator), node_count))
        )
        column_nodes = tuple(
            Fraction(value, column_denominator)
            for value in sorted(rng.sample(range(1, column_denominator), node_count))
        )
        weights = tuple(
            Fraction(rng.randrange(8), rng.randrange(1, 6))
            for _ in range(feature_count)
        )
        prefix_length = rng.randrange(1, min(3, node_count))
        selected_rows = tuple(rng.sample(range(node_count), prefix_length))
        selected_columns = tuple(rng.sample(range(node_count), prefix_length))
        available_rows = [i for i in range(node_count) if i not in selected_rows]
        available_columns = [j for j in range(node_count) if j not in selected_columns]
        fresh_order = rng.randrange(1, min(3, len(available_rows)) + 1)
        fresh_rows = tuple(sorted(rng.sample(available_rows, fresh_order)))
        fresh_columns = tuple(sorted(rng.sample(available_columns, fresh_order)))
        determinant = conditioned_minor(
            degrees,
            weights,
            row_nodes,
            column_nodes,
            selected_rows,
            selected_columns,
            fresh_rows,
            fresh_columns,
        )
        cross_inversions = sum(
            selected > fresh for selected in selected_rows for fresh in fresh_rows
        ) + sum(
            selected > fresh for selected in selected_columns for fresh in fresh_columns
        )
        orientation = -1 if cross_inversions % 2 else 1
        assert orientation * determinant >= 0
        positive_cases += orientation * determinant > 0

    assert positive_cases == 3964

    signed_weight_failure = conditioned_minor(
        (0, 1),
        (Fraction(-2), Fraction()),
        (Fraction(1, 3), Fraction(2, 3)),
        (Fraction(1, 4), Fraction(3, 4)),
        (0,),
        (0,),
        (1,),
        (1,),
    )
    assert signed_weight_failure == Fraction(-4, 169)


def test_interleaved_tangent_middle_minors_exactly() -> None:
    def inversion_orientation(values: tuple[int, ...]) -> int:
        inversions = sum(
            values[i] > values[j]
            for i in range(len(values))
            for j in range(i + 1, len(values))
        )
        return -1 if inversions % 2 else 1

    checked_feature_major_minors = 0
    checked_layer_major_minors = 0
    found_negative_layer_major_minor = False
    for feature_count in range(1, 6):
        # Mix positive and zero edge weights.  The formal theorem permits both,
        # and zero weights exercise the weak (rather than strict) sign boundary.
        weights = tuple(
            Fraction() if feature % 3 == 1 else Fraction(feature + 1, feature + 2)
            for feature in range(feature_count)
        )
        order = 2 * feature_count
        feature_major = [
            [
                Fraction(1)
                if i == j
                else weights[i // 2]
                if i % 2 == 0 and j == i + 1
                else Fraction()
                for j in range(order)
            ]
            for i in range(order)
        ]

        def interleave(layer_major_index: int, count: int = feature_count) -> int:
            if layer_major_index < count:
                return 2 * layer_major_index
            return 2 * (layer_major_index - count) + 1

        for minor_order in range(order + 1):
            flags = tuple(itertools.combinations(range(order), minor_order))
            for rows in flags:
                for columns in flags:
                    feature_major_det = fraction_determinant(
                        [[feature_major[i][j] for j in columns] for i in rows]
                    )
                    assert feature_major_det >= 0
                    checked_feature_major_minors += 1

                    interleaved_rows = tuple(interleave(i) for i in rows)
                    interleaved_columns = tuple(interleave(j) for j in columns)
                    layer_major_det = fraction_determinant(
                        [
                            [feature_major[i][j] for j in interleaved_columns]
                            for i in interleaved_rows
                        ]
                    )
                    orientation = inversion_orientation(
                        interleaved_rows
                    ) * inversion_orientation(interleaved_columns)
                    assert orientation * layer_major_det >= 0
                    found_negative_layer_major_minor |= layer_major_det < 0
                    checked_layer_major_minors += 1

    expected = sum(
        math.comb(4 * feature_count, 2 * feature_count) for feature_count in range(1, 6)
    )
    assert checked_feature_major_minors == expected
    assert checked_layer_major_minors == expected
    assert found_negative_layer_major_minor


def test_noncrossing_support_kills_nonidentity_matchings_exactly() -> None:
    weights = (Fraction(2, 3), Fraction(), Fraction(4, 5))
    order = 2 * len(weights)
    matrix = [
        [
            Fraction(1)
            if i == j
            else weights[i // 2]
            if i % 2 == 0 and j == i + 1
            else Fraction()
            for j in range(order)
        ]
        for i in range(order)
    ]

    for x1 in range(order):
        for x2 in range(x1 + 1, order):
            for y1 in range(order):
                for y2 in range(y1 + 1, order):
                    assert matrix[x2][y1] == 0 or matrix[x1][y2] == 0

    checked_nonidentity_terms = 0
    for minor_order in range(order + 1):
        flags = tuple(itertools.combinations(range(order), minor_order))
        identity = tuple(range(minor_order))
        for rows in flags:
            for columns in flags:
                for permutation in itertools.permutations(range(minor_order)):
                    if permutation == identity:
                        continue
                    term = math.prod(
                        matrix[rows[permutation[j]]][columns[j]]
                        for j in range(minor_order)
                    )
                    assert term == 0
                    checked_nonidentity_terms += 1

    expected = sum(
        math.comb(order, minor_order) ** 2 * (math.factorial(minor_order) - 1)
        for minor_order in range(order + 1)
    )
    assert checked_nonidentity_terms == expected


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
