import itertools
import math
import random
from fractions import Fraction

import mpmath as mp

from kernelgecp.surrogates import (
    exact_gecp,
    exact_gecp_pivot_sign_coherence,
    first_exact_gecp_sign_coherence_failure,
    fraction_determinant,
    geometric_surrogate,
    inspect_surrogate,
)


def test_fraction_determinant_and_exact_pivots() -> None:
    matrix = geometric_surrogate(4, Fraction(2, 3))
    assert (
        fraction_determinant([[Fraction(2), Fraction(1)], [Fraction(1), Fraction(3)]])
        == 5
    )
    rows, columns, pivots = exact_gecp(matrix)
    assert len(rows) == len(columns) == len(pivots) == 4
    assert all(pivot != 0 for pivot in pivots)


def test_bordered_core_determinant_identity_exactly() -> None:
    matrix = geometric_surrogate(4, Fraction(2, 3))
    pivot_rows, pivot_columns, pivots = exact_gecp(matrix)
    residual = [row.copy() for row in matrix]
    selected_rows: list[int] = []
    selected_columns: list[int] = []
    for step, pivot in enumerate(pivots):
        core = [
            [matrix[row][column] for column in selected_columns]
            for row in selected_rows
        ]
        core_det = fraction_determinant(core)
        for row in range(4):
            for column in range(4):
                bordered_rows = [*selected_rows, row]
                bordered_columns = [*selected_columns, column]
                bordered = [
                    [matrix[i][j] for j in bordered_columns] for i in bordered_rows
                ]
                assert (
                    fraction_determinant(bordered) == core_det * residual[row][column]
                )
        pivot_row = pivot_rows[step]
        pivot_column = pivot_columns[step]
        assert residual[pivot_row][pivot_column] == pivot
        pivot_column_values = [residual[row][pivot_column] for row in range(4)]
        pivot_row_values = residual[pivot_row].copy()
        for row in range(4):
            for column in range(4):
                residual[row][column] -= (
                    pivot_column_values[row] * pivot_row_values[column] / pivot
                )
        selected_rows.append(pivot_row)
        selected_columns.append(pivot_column)


def test_low_rank_perturbation_determinant_bound_exactly() -> None:
    rng = random.Random(0xC0FFEE)

    def product(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum((left[i][a] * right[a][j] for a in range(len(right))), Fraction())
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    for size in range(2, 6):
        for rank in range(1, size):
            for _trial in range(3):
                left = [
                    [Fraction(rng.randint(-2, 2)) for _ in range(rank)]
                    for _ in range(size)
                ]
                right = [
                    [Fraction(rng.randint(-2, 2)) for _ in range(size)]
                    for _ in range(rank)
                ]
                background = product(left, right)
                error = [
                    [Fraction(rng.randint(-2, 2), 7) for _ in range(size)]
                    for _ in range(size)
                ]
                error[0][0] = Fraction(1, 7)
                matrix = [
                    [background[i][j] + error[i][j] for j in range(size)]
                    for i in range(size)
                ]
                epsilon = max(abs(value) for row in error for value in row)
                column_bound = max(
                    Fraction(1),
                    max(abs(value) for row in background for value in row),
                )
                determinant_sum = Fraction()
                mixed_bound_sum = Fraction()
                for mask in range(1 << size):
                    error_columns = mask.bit_count()
                    mixed = [
                        [
                            error[i][j] if mask & (1 << j) else background[i][j]
                            for j in range(size)
                        ]
                        for i in range(size)
                    ]
                    mixed_det = fraction_determinant(mixed)
                    determinant_sum += mixed_det
                    if rank + error_columns < size:
                        assert mixed_det == 0
                    term_bound = (
                        math.factorial(size)
                        * epsilon**error_columns
                        * column_bound ** (size - error_columns)
                    )
                    assert abs(mixed_det) <= term_bound
                    if rank + error_columns >= size:
                        mixed_bound_sum += term_bound

                determinant = fraction_determinant(matrix)
                assert determinant == determinant_sum
                assert abs(determinant) <= mixed_bound_sum
                assert abs(determinant) <= (
                    2**size
                    * math.factorial(size)
                    * epsilon ** (size - rank)
                    * column_bound**rank
                )


def test_exceptional_columns_shift_low_rank_determinant_threshold_exactly() -> None:
    rng = random.Random(0xECCE)

    def product(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                sum((left[i][a] * right[a][j] for a in range(len(right))), Fraction())
                for j in range(len(right[0]))
            ]
            for i in range(len(left))
        ]

    size = 6
    rank = 2
    left = [[Fraction(rng.randint(-3, 3)) for _ in range(rank)] for _ in range(size)]
    right = [[Fraction(rng.randint(-3, 3)) for _ in range(size)] for _ in range(rank)]
    low_rank = product(left, right)

    for exceptional_count in (0, 2):
        exceptional = set(range(exceptional_count))
        exact_background = [row.copy() for row in low_rank]
        for i in range(size):
            for j in exceptional:
                exact_background[i][j] = Fraction(rng.randint(-5, 5))

        error = [[Fraction() for _ in range(size)] for _ in range(size)]
        for i in range(size):
            for j in range(size):
                if j not in exceptional:
                    error[i][j] = Fraction(rng.randint(-2, 2), 11)
        error[0][exceptional_count] = Fraction(1, 11)
        matrix = [
            [exact_background[i][j] + error[i][j] for j in range(size)]
            for i in range(size)
        ]

        effective_rank = rank + exceptional_count
        epsilon = max(abs(value) for row in error for value in row)
        column_bound = max(
            Fraction(1),
            max(abs(value) for row in exact_background for value in row),
        )
        determinant_sum = Fraction()
        for mask in range(1 << size):
            error_columns = mask.bit_count()
            mixed = [
                [
                    error[i][j] if mask & (1 << j) else exact_background[i][j]
                    for j in range(size)
                ]
                for i in range(size)
            ]
            mixed_det = fraction_determinant(mixed)
            determinant_sum += mixed_det
            if effective_rank + error_columns < size:
                assert mixed_det == 0

        determinant = fraction_determinant(matrix)
        assert determinant == determinant_sum
        assert abs(determinant) <= (
            2**size
            * math.factorial(size)
            * epsilon ** (size - effective_rank)
            * column_bound**effective_rank
        )


def test_exceptional_columns_compose_with_exact_gecp_selected_core() -> None:
    size = 6
    rank = 2
    exceptional_columns = {0, 1}
    left = [[Fraction(i + 1), Fraction((i + 1) ** 2)] for i in range(size)]
    right = [
        [Fraction(j + 2) for j in range(size)],
        [Fraction((j + 2) ** 2) for j in range(size)],
    ]
    low_rank = [
        [
            sum((left[i][a] * right[a][j] for a in range(rank)), Fraction())
            for j in range(size)
        ]
        for i in range(size)
    ]
    exact_background = [row.copy() for row in low_rank]
    for i in range(size):
        exact_background[i][0] = Fraction(7 * i + 3)
        exact_background[i][1] = Fraction(5 * i * i - 2 * i + 1)
    error = [[Fraction() for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if j not in exceptional_columns:
                error[i][j] = Fraction(int(i == j), 13)
    matrix = [
        [exact_background[i][j] + error[i][j] for j in range(size)] for i in range(size)
    ]

    pivot_rows, pivot_columns, pivots = exact_gecp(matrix)
    assert len(pivots) == size
    selected_core = [[matrix[i][j] for j in pivot_columns] for i in pivot_rows]
    selected_exceptions = {
        index
        for index, column in enumerate(pivot_columns)
        if column in exceptional_columns
    }
    assert len(selected_exceptions) == 2
    assert fraction_determinant(selected_core) == math.prod(pivots)

    epsilon = max(abs(value) for row in error for value in row)
    column_bound = max(
        Fraction(1),
        max(abs(value) for row in exact_background for value in row),
    )
    effective_rank = rank + len(selected_exceptions)
    assert abs(fraction_determinant(selected_core)) <= (
        2**size
        * math.factorial(size)
        * epsilon ** (size - effective_rank)
        * column_bound**effective_rank
    )


def test_prescribed_two_corner_run_has_exactly_two_exceptional_columns() -> None:
    size = 7
    matrix = geometric_surrogate(size, Fraction(2, 3))
    residual = [row.copy() for row in matrix]
    selected_rows: list[int] = []
    selected_columns: list[int] = []
    pivots: list[Fraction] = []

    def update(pivot_row: int, pivot_column: int) -> None:
        pivot = residual[pivot_row][pivot_column]
        assert pivot != 0
        selected_rows.append(pivot_row)
        selected_columns.append(pivot_column)
        pivots.append(pivot)
        pivot_column_values = [residual[i][pivot_column] for i in range(size)]
        pivot_row_values = residual[pivot_row].copy()
        for i in range(size):
            for j in range(size):
                residual[i][j] -= pivot_column_values[i] * pivot_row_values[j] / pivot

    update(0, size - 1)
    update(size - 1, 0)
    central_columns = range(1, size - 1)
    for _step in range(size - 2):
        pivot_row, pivot_column = max(
            itertools.product(range(size), central_columns),
            key=lambda coordinate: (
                abs(residual[coordinate[0]][coordinate[1]]),
                -coordinate[0],
                -coordinate[1],
            ),
        )
        update(pivot_row, pivot_column)

    exceptional_indices = {0, 1}
    assert selected_columns[:2] == [size - 1, 0]
    assert len(exceptional_indices) == 2
    assert all(
        1 <= column < size - 1
        for index, column in enumerate(selected_columns)
        if index not in exceptional_indices
    )
    selected_core = [[matrix[i][j] for j in selected_columns] for i in selected_rows]
    assert fraction_determinant(selected_core) == math.prod(pivots)


def test_above_quarter_trajectory_prefix_reconstructs_final_residual_exactly() -> None:
    size = 7
    matrix = [
        [64 * value for value in row]
        for row in geometric_surrogate(size, Fraction(2, 3))
    ]
    residual = [row.copy() for row in matrix]

    def update(pivot_row: int, pivot_column: int) -> Fraction:
        pivot = residual[pivot_row][pivot_column]
        assert pivot != 0
        pivot_column_values = [residual[i][pivot_column] for i in range(size)]
        pivot_row_values = residual[pivot_row].copy()
        for i in range(size):
            for j in range(size):
                residual[i][j] -= pivot_column_values[i] * pivot_row_values[j] / pivot
        return pivot

    update(0, size - 1)
    update(size - 1, 0)
    initial_continuation_residual = [row.copy() for row in residual]
    trajectory = [[row.copy() for row in residual]]
    prefix_rows: list[int] = []
    prefix_columns: list[int] = []
    prefix_pivots: list[Fraction] = []

    while True:
        pivot_row, pivot_column = max(
            itertools.product(range(size), range(size)),
            key=lambda coordinate: (
                abs(residual[coordinate[0]][coordinate[1]]),
                -coordinate[0],
                -coordinate[1],
            ),
        )
        pivot = residual[pivot_row][pivot_column]
        if abs(pivot) <= Fraction(1, 4):
            break
        assert 1 <= pivot_column < size - 1
        assert abs(pivot) == max(abs(value) for row in residual for value in row)
        prefix_rows.append(pivot_row)
        prefix_columns.append(pivot_column)
        prefix_pivots.append(update(pivot_row, pivot_column))
        trajectory.append([row.copy() for row in residual])

    assert len(prefix_pivots) == 3
    assert all(abs(pivot) > Fraction(1, 4) for pivot in prefix_pivots)
    assert max(abs(value) for row in residual for value in row) <= Fraction(1, 4)

    reconstructed = [row.copy() for row in initial_continuation_residual]
    for step, (pivot_row, pivot_column, pivot) in enumerate(
        zip(prefix_rows, prefix_columns, prefix_pivots, strict=True)
    ):
        assert reconstructed == trajectory[step]
        assert reconstructed[pivot_row][pivot_column] == pivot
        pivot_column_values = [reconstructed[i][pivot_column] for i in range(size)]
        pivot_row_values = reconstructed[pivot_row].copy()
        for i in range(size):
            for j in range(size):
                reconstructed[i][j] -= (
                    pivot_column_values[i] * pivot_row_values[j] / pivot
                )
    assert reconstructed == trajectory[len(prefix_pivots)] == residual


def test_above_quarter_stopping_budget_exactly() -> None:
    transition_neighbors = {
        neighbor
        for exponent in range(1, 17)
        for neighbor in (2**exponent - 1, 2**exponent, 2**exponent + 1)
    }
    for scale_plus_one in sorted(set(range(1, 257)) | transition_neighbors):
        approximation_order = 24 * scale_plus_one
        continuation_length = 2048 * scale_plus_one**2
        binary_bits = 2 * scale_plus_one + 10
        effective_rank = 16 * scale_plus_one * approximation_order
        assert effective_rank <= continuation_length
        assert continuation_length + 2 <= 2**binary_bits
        assert (binary_bits + 3) * (
            continuation_length + 2
        ) + effective_rank <= approximation_order * (
            continuation_length - effective_rank
        )


def test_quadratic_block_satisfies_exact_half_contraction_condition() -> None:
    for scale in range(33):
        accuracy_order = 16 * (scale + 1)
        assert 128 * (2 * accuracy_order + 1) * (scale + 1) <= 2**accuracy_order


def test_logarithmic_block_satisfies_exact_half_contraction_condition() -> None:
    dyadic_neighbors = {
        neighbor
        for exponent in range(1, 65)
        for neighbor in (2**exponent - 1, 2**exponent, 2**exponent + 1)
    }
    scales_plus_one = set(range(1, 1025)) | dyadic_neighbors
    for scale_plus_one in sorted(scales_plus_one):
        ceil_log_two = (scale_plus_one - 1).bit_length()
        accuracy_order = 16 + 2 * ceil_log_two
        assert 128 * (2 * accuracy_order + 1) * scale_plus_one <= 2**accuracy_order


def test_accuracy_block_satisfies_exact_dyadic_target_condition() -> None:
    dyadic_neighbors = {
        neighbor
        for exponent in range(1, 65)
        for neighbor in (2**exponent - 1, 2**exponent, 2**exponent + 1)
    }
    scales_plus_one = set(range(1, 1025)) | dyadic_neighbors
    for scale_plus_one in sorted(scales_plus_one):
        ceil_log_two = (scale_plus_one - 1).bit_length()
        for target_order in range(65):
            accuracy_order = 16 + 2 * ceil_log_two + 2 * target_order
            assert (
                64 * (2 * accuracy_order + 1) * scale_plus_one * 2**target_order
                <= 2**accuracy_order
            )


def test_sylvester_hadamard_determinant_obstruction_exactly() -> None:
    def kronecker(
        left: list[list[Fraction]], right: list[list[Fraction]]
    ) -> list[list[Fraction]]:
        return [
            [
                left_value * right_value
                for left_value in left_row
                for right_value in right_row
            ]
            for left_row in left
            for right_row in right
        ]

    hadamard_two = [[Fraction(1), Fraction(1)], [Fraction(1), Fraction(-1)]]
    matrix = hadamard_two
    for exponent in range(1, 7):
        order = 2**exponent
        assert len(matrix) == order
        assert all(abs(value) == 1 for row in matrix for value in row)
        determinant = fraction_determinant(matrix)
        assert determinant**2 == order**order
        matrix = kronecker(matrix, hadamard_two)

    for constant_base in range(257):
        exponent = 2 * constant_base + 2
        assert constant_base**2 < 2**exponent


def test_fermionic_corner_cross_ratio_localization_high_precision() -> None:
    with mp.workdps(100):
        tolerance = mp.mpf("1e-90")

        def kernel(time: mp.mpf, frequency: mp.mpf) -> mp.mpf:
            return mp.exp(-time * frequency) / (1 + mp.exp(-frequency))

        def update(
            pivot_time: mp.mpf,
            pivot_frequency: mp.mpf,
            time: mp.mpf,
            frequency: mp.mpf,
        ) -> mp.mpf:
            return kernel(time, frequency) - (
                kernel(time, pivot_frequency)
                * kernel(pivot_time, frequency)
                / kernel(pivot_time, pivot_frequency)
            )

        arbitrary_cases = [
            ("0.17", "-1.3", "0.81", "2.4"),
            ("0.72", "3.1", "0.08", "-0.9"),
            ("0", "8", "0.375", "0.4"),
            ("1", "-8", "0.625", "-0.4"),
        ]
        for raw_case in arbitrary_cases:
            pivot_time, pivot_frequency, time, frequency = map(mp.mpf, raw_case)
            direct = update(pivot_time, pivot_frequency, time, frequency)
            factorized = kernel(time, frequency) * (
                1 - mp.exp((time - pivot_time) * (frequency - pivot_frequency))
            )
            assert abs(direct - factorized) < tolerance

        logarithmic_half_area = mp.log(2)
        grid = [mp.mpf(index) / 16 for index in range(17)]
        for pivot_time in grid:
            for time in grid:
                if time < pivot_time:
                    continue
                for pivot_frequency in map(mp.mpf, ("-2", "0", "3")):
                    for frequency_gap in map(mp.mpf, ("0", "0.125", "0.5", "1", "2")):
                        frequency = pivot_frequency - frequency_gap
                        area = (time - pivot_time) * frequency_gap
                        residual = update(pivot_time, pivot_frequency, time, frequency)
                        assert residual >= -tolerance
                        assert abs(residual) <= area + tolerance
                        if area <= logarithmic_half_area:
                            assert abs(residual) <= mp.mpf("0.5") + tolerance
                        reflected_residual = update(
                            1 - pivot_time,
                            -pivot_frequency,
                            1 - time,
                            -frequency,
                        )
                        assert reflected_residual >= -tolerance
                        assert abs(reflected_residual) <= area + tolerance
                        if area <= logarithmic_half_area:
                            assert abs(reflected_residual) <= mp.mpf("0.5") + tolerance


def test_fermionic_two_corner_power_secant_tent_high_precision() -> None:
    with mp.workdps(100):
        tolerance = mp.mpf("1e-85")

        def kernel(time: mp.mpf, frequency: mp.mpf) -> mp.mpf:
            return mp.exp(-time * frequency) / (1 + mp.exp(-frequency))

        def first_residual(
            upper_frequency: mp.mpf,
            time: mp.mpf,
            frequency: mp.mpf,
        ) -> mp.mpf:
            return kernel(time, frequency) - (
                kernel(time, upper_frequency)
                * kernel(0, frequency)
                / kernel(0, upper_frequency)
            )

        def two_corner_residual(
            lower_frequency: mp.mpf,
            upper_frequency: mp.mpf,
            time: mp.mpf,
            frequency: mp.mpf,
        ) -> mp.mpf:
            pivot = first_residual(upper_frequency, 1, lower_frequency)
            return first_residual(upper_frequency, time, frequency) - (
                first_residual(upper_frequency, time, lower_frequency)
                * first_residual(upper_frequency, 1, frequency)
                / pivot
            )

        bands = [
            (mp.mpf("-3"), mp.mpf("2")),
            (mp.mpf("-0.75"), mp.mpf("4.5")),
            (mp.mpf("-2"), mp.mpf("3")),
            (mp.mpf("-4.5"), mp.mpf("0.75")),
        ]
        grid = [mp.mpf(index) / 16 for index in range(17)]
        for lower_frequency, upper_frequency in bands:
            a = mp.exp(-upper_frequency)
            b = mp.exp(-lower_frequency)
            for time in grid:
                for frequency_weight in grid:
                    frequency = lower_frequency + frequency_weight * (
                        upper_frequency - lower_frequency
                    )
                    z = mp.exp(-frequency)
                    secant = ((b - z) * a**time + (z - a) * b**time) / (b - a)
                    formula = (z**time - secant) / (1 + z)
                    residual = two_corner_residual(
                        lower_frequency, upper_frequency, time, frequency
                    )
                    upper_area = time * (upper_frequency - frequency)
                    lower_area = (1 - time) * (frequency - lower_frequency)

                    assert abs(residual - formula) < tolerance
                    assert residual >= -tolerance
                    assert residual <= min(upper_area, lower_area) + tolerance
                    if time in (0, 1) or frequency_weight in (0, 1):
                        assert abs(residual) < tolerance


def test_two_corner_half_contraction_obstruction_exactly() -> None:
    endpoint_base = Fraction(24)
    interior_base = Fraction(5, 4)
    lower_endpoint = endpoint_base**-3
    upper_endpoint = endpoint_base**3
    coordinate = interior_base**3
    secant = (upper_endpoint - coordinate) / (
        upper_endpoint - lower_endpoint
    ) * endpoint_base**-2 + (coordinate - lower_endpoint) / (
        upper_endpoint - lower_endpoint
    ) * endpoint_base**2
    residual = (interior_base**2 - secant) / (1 + coordinate)
    initial_pivot = upper_endpoint / (1 + upper_endpoint)

    assert secant == Fraction(31569, 379832)
    assert residual == Fraction(4495348, 8973531)
    assert residual - initial_pivot / 2 == Fraction(222676, 224338275)

    for integer_base in range(24, 65):
        base = Fraction(integer_base)
        lower_endpoint = base**-3
        upper_endpoint = base**3
        secant = (upper_endpoint - coordinate) / (
            upper_endpoint - lower_endpoint
        ) * base**-2 + (coordinate - lower_endpoint) / (
            upper_endpoint - lower_endpoint
        ) * base**2
        residual = (interior_base**2 - secant) / (1 + coordinate)
        initial_pivot = upper_endpoint / (1 + upper_endpoint)
        assert residual > initial_pivot / 2

    with mp.workdps(100):
        tolerance = mp.mpf("1e-90")
        cutoff = 3 * mp.log(24)
        time = mp.mpf(2) / 3
        frequency = -3 * mp.log(mp.mpf(5) / 4)

        def kernel(sample_time: mp.mpf, sample_frequency: mp.mpf) -> mp.mpf:
            return mp.exp(-sample_time * sample_frequency) / (
                1 + mp.exp(-sample_frequency)
            )

        def first_residual(sample_time: mp.mpf, sample_frequency: mp.mpf) -> mp.mpf:
            return kernel(sample_time, sample_frequency) - (
                kernel(sample_time, cutoff)
                * kernel(0, sample_frequency)
                / kernel(0, cutoff)
            )

        nested_residual = first_residual(time, frequency) - (
            first_residual(time, -cutoff)
            * first_residual(1, frequency)
            / first_residual(1, -cutoff)
        )
        exact_residual = mp.mpf(4495348) / 8973531
        exact_margin = mp.mpf(222676) / 224338275

        assert -cutoff < frequency < cutoff
        assert abs(nested_residual - exact_residual) < tolerance
        assert abs(nested_residual - kernel(0, cutoff) / 2 - exact_margin) < tolerance


def test_third_pivot_frequency_localization_high_precision() -> None:
    with mp.workdps(100):
        tolerance = mp.mpf("1e-80")

        def kernel(time: mp.mpf, frequency: mp.mpf) -> mp.mpf:
            return mp.exp(-time * frequency) / (1 + mp.exp(-frequency))

        def two_corner_residual(
            cutoff: mp.mpf, time: mp.mpf, frequency: mp.mpf
        ) -> mp.mpf:
            def first_residual(sample_time: mp.mpf, sample_frequency: mp.mpf) -> mp.mpf:
                return kernel(sample_time, sample_frequency) - (
                    kernel(sample_time, cutoff)
                    * kernel(0, sample_frequency)
                    / kernel(0, cutoff)
                )

            return first_residual(time, frequency) - (
                first_residual(time, -cutoff)
                * first_residual(1, frequency)
                / first_residual(1, -cutoff)
            )

        cutoffs = [2 * mp.log(4), mp.mpf(3), mp.mpf(4), mp.mpf(8), mp.mpf(16)]
        times = [mp.mpf(index) / 40 for index in range(41)]
        frequency_weights = [mp.mpf(index) / 80 for index in range(81)]
        for cutoff in cutoffs:
            center = two_corner_residual(cutoff, mp.mpf("0.5"), mp.mpf(0))
            center_formula = mp.mpf("0.5") - 1 / (2 * mp.cosh(cutoff / 2))
            assert abs(center - center_formula) < tolerance
            assert center > mp.mpf("0.25")

            sampled_maximum = -mp.inf
            sampled_frequency = mp.mpf(0)
            for time in times:
                for weight in frequency_weights:
                    frequency = -cutoff + 2 * cutoff * weight
                    residual = two_corner_residual(cutoff, time, frequency)
                    if abs(frequency) >= cutoff / 2:
                        assert residual <= mp.mpf("0.25") + tolerance
                    if residual > sampled_maximum:
                        sampled_maximum = residual
                        sampled_frequency = frequency

            assert sampled_maximum >= center - tolerance
            assert abs(sampled_frequency) < cutoff / 2


def test_outer_half_quarter_bound_survives_later_pivots() -> None:
    quarter = Fraction(1, 4)
    mixed_bound_cases = [
        (Fraction(2), quarter, Fraction(1), quarter),
        (Fraction(2), -quarter, Fraction(1), -quarter),
        (Fraction(-2), quarter, Fraction(1), -quarter),
        (Fraction(-2), -quarter, Fraction(1), quarter),
    ]
    for pivot, value, selected_column, selected_row in mixed_bound_cases:
        assert abs(value) <= quarter
        assert abs(selected_row) <= quarter
        assert abs(selected_column) <= abs(pivot)
        assert (value * pivot) * (selected_column * selected_row) >= 0
        updated = value - selected_column * selected_row / pivot
        assert abs(updated) <= quarter

    with mp.workdps(80):
        tolerance = mp.mpf("1e-60")

        def kernel(time: mp.mpf, frequency: mp.mpf) -> mp.mpf:
            return mp.exp(-time * frequency) / (1 + mp.exp(-frequency))

        for cutoff in map(mp.mpf, ("4", "8", "16")):
            times = [mp.mpf(index) / 30 for index in range(31)]
            frequencies = [
                -cutoff + 2 * cutoff * mp.mpf(index) / 120 for index in range(121)
            ]

            first_pivot = kernel(0, cutoff)
            first = [
                [
                    kernel(time, frequency)
                    - kernel(time, cutoff) * kernel(0, frequency) / first_pivot
                    for frequency in frequencies
                ]
                for time in times
            ]
            second_pivot = first[-1][0]
            residual = [
                [
                    first[i][j] - first[i][0] * first[-1][j] / second_pivot
                    for j in range(len(frequencies))
                ]
                for i in range(len(times))
            ]

            for _step in range(5):
                outer_maximum = max(
                    abs(residual[i][j])
                    for i in range(len(times))
                    for j, frequency in enumerate(frequencies)
                    if abs(frequency) >= cutoff / 2
                )
                assert outer_maximum <= mp.mpf("0.25") + tolerance

                pivot_i, pivot_j = max(
                    (
                        (i, j)
                        for i in range(len(times))
                        for j in range(len(frequencies))
                    ),
                    key=lambda index: abs(residual[index[0]][index[1]]),
                )
                pivot = residual[pivot_i][pivot_j]
                assert abs(pivot) > tolerance
                if abs(pivot) > mp.mpf("0.25") + tolerance:
                    assert abs(frequencies[pivot_j]) < cutoff / 2

                residual = [
                    [
                        residual[i][j]
                        - residual[i][pivot_j] * residual[pivot_i][j] / pivot
                        for j in range(len(frequencies))
                    ]
                    for i in range(len(times))
                ]


def test_geometric_minors_have_expected_sign() -> None:
    for size in range(2, 6):
        record = inspect_surrogate(size, Fraction(3, 4))
        assert not any(record.minor_sign_mismatches.values())
        assert not any(record.zero_minors.values())


def test_geometric_minor_sign_tracks_row_and_column_orientation_exactly() -> None:
    def permutation_sign(permutation: tuple[int, ...]) -> int:
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(len(permutation))
            for j in range(i + 1, len(permutation))
        )
        return -1 if inversions % 2 else 1

    for size in range(2, 5):
        matrix = geometric_surrogate(size, Fraction(2, 3))
        ordered_sign = -1 if (size * (size - 1) // 2) % 2 else 1
        permutations = list(itertools.permutations(range(size)))
        for row_permutation in permutations:
            for column_permutation in permutations:
                permuted = [
                    [matrix[row][column] for column in column_permutation]
                    for row in row_permutation
                ]
                determinant = fraction_determinant(permuted)
                expected_sign = (
                    ordered_sign
                    * permutation_sign(row_permutation)
                    * permutation_sign(column_permutation)
                )
                assert determinant != 0
                assert (1 if determinant > 0 else -1) == expected_sign


def test_exp_taylor_principal_determinant_matches_vandermonde_product() -> None:
    def vandermonde_product(values: list[Fraction]) -> Fraction:
        return math.prod(
            values[j] - values[i]
            for i in range(len(values))
            for j in range(i + 1, len(values))
        )

    for size in range(1, 7):
        rows = [Fraction(2 * i + 1, size + 1) for i in range(size)]
        columns = [Fraction(3 * j + 2, size + 2) for j in range(size)]
        principal = [
            [
                sum(
                    rows[i] ** k * Fraction(1, math.factorial(k)) * columns[j] ** k
                    for k in range(size)
                )
                for j in range(size)
            ]
            for i in range(size)
        ]
        expected = (
            vandermonde_product(rows)
            * math.prod(Fraction(1, math.factorial(k)) for k in range(size))
            * vandermonde_product(columns)
        )
        assert fraction_determinant(principal) == expected
        assert expected > 0


def test_exp_taylor_rectangular_feature_determinants_are_positive_exactly() -> None:
    for size in range(1, 6):
        rows = [Fraction(i, size) for i in range(size)]
        columns = [Fraction(2 * j, size + 1) for j in range(size)]
        for terms in range(size, size + 4):
            truncated = [
                [
                    sum(
                        rows[i] ** k * Fraction(1, math.factorial(k)) * columns[j] ** k
                        for k in range(terms)
                    )
                    for j in range(size)
                ]
                for i in range(size)
            ]
            assert fraction_determinant(truncated) > 0


def test_exp_matrix_determinant_retains_principal_lower_bound() -> None:
    with mp.workdps(100):
        for size in range(1, 6):
            rows = [mp.mpf(i) / size for i in range(size)]
            columns = [mp.mpf(2 * j) / (size + 1) for j in range(size)]
            matrix = mp.matrix(
                [
                    [mp.exp(rows[i] * columns[j]) for j in range(size)]
                    for i in range(size)
                ]
            )
            row_vandermonde = mp.fprod(
                rows[j] - rows[i] for i in range(size) for j in range(i + 1, size)
            )
            column_vandermonde = mp.fprod(
                columns[j] - columns[i] for i in range(size) for j in range(i + 1, size)
            )
            principal = (
                row_vandermonde
                * mp.fprod(mp.mpf(1) / math.factorial(k) for k in range(size))
                * column_vandermonde
            )
            determinant = mp.det(matrix)
            assert determinant > 0
            assert determinant >= principal * (1 - mp.mpf("1e-80"))


def test_exp_kernel_determinants_have_all_ordered_signatures() -> None:
    with mp.workdps(100):
        for size in range(1, 7):
            rows = [mp.mpf(3 * i - size) / (size + 1) for i in range(size)]
            columns = [mp.mpf(2 * j - size) / (size + 2) for j in range(size)]
            matrix = mp.matrix(
                [
                    [mp.exp(-rows[i] * columns[j]) for j in range(size)]
                    for i in range(size)
                ]
            )
            determinant = mp.det(matrix)
            expected_sign = -1 if math.comb(size, 2) % 2 else 1
            assert determinant != 0
            assert mp.sign(determinant) == expected_sign


def test_generalized_vandermonde_determinants_are_positive_exactly() -> None:
    for size in range(1, 6):
        nodes = [Fraction(2 * i + 1, size + 2) for i in range(size)]
        for exponents in itertools.combinations(range(size + 3), size):
            matrix = [[node**exponent for exponent in exponents] for node in nodes]
            assert fraction_determinant(matrix) > 0


def test_generalized_vandermonde_zero_node_boundary_exactly() -> None:
    for size in range(1, 6):
        nodes = [Fraction(i, size) for i in range(size)]
        for exponents in itertools.combinations(range(size + 3), size):
            matrix = [[node**exponent for exponent in exponents] for node in nodes]
            determinant = fraction_determinant(matrix)
            assert determinant >= 0
            assert (determinant > 0) == (exponents[0] == 0)


def test_geometric_surrogate_selected_crosses_are_sign_coherent() -> None:
    for size in range(2, 9):
        for q in (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)):
            checks = exact_gecp_pivot_sign_coherence(geometric_surrogate(size, q))
            assert len(checks) == size
            assert all(checks)
            assert (
                first_exact_gecp_sign_coherence_failure(geometric_surrogate(size, q))
                is None
            )


def test_exact_sign_coherence_failure_retains_minimized_witness() -> None:
    failure = first_exact_gecp_sign_coherence_failure(
        [
            [Fraction(-4), Fraction(-1)],
            [Fraction(-3), Fraction(4)],
        ]
    )
    assert failure is not None
    assert (
        failure.step,
        failure.pivot_row,
        failure.pivot_column,
        failure.row,
        failure.column,
        failure.cross_product,
    ) == (0, 0, 0, 1, 1, Fraction(-48))


def test_minimized_sign_regular_pivot_obstruction() -> None:
    left = [[Fraction(4), Fraction(3)], [Fraction(2), Fraction(1)]]
    right = [[Fraction(1), Fraction(2)], [Fraction(3), Fraction(4)]]
    assert fraction_determinant(left) == fraction_determinant(right) == -2
    left_rows, left_columns, _ = exact_gecp(left)
    right_rows, right_columns, _ = exact_gecp(right)
    assert (left_rows[0], left_columns[0]) == (0, 0)
    assert (right_rows[0], right_columns[0]) == (1, 1)
