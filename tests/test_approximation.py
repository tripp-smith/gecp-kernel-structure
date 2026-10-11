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


def test_block_feature_cauchy_binet_coefficients_exactly() -> None:
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
