import itertools
from fractions import Fraction

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
