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
