"""Independent exact tests of the all-orders, matched-rank Fischer bridge."""

import itertools
import math
from fractions import Fraction


def _det(matrix: list[list[Fraction]]) -> Fraction:
    work = [line.copy() for line in matrix]
    value = Fraction(1)
    for k in range(len(work)):
        p = next((i for i in range(k, len(work)) if work[i][k]), None)
        if p is None:
            return Fraction()
        if p != k:
            work[k], work[p] = work[p], work[k]
            value = -value
        pivot = work[k][k]
        value *= pivot
        for i in range(k + 1, len(work)):
            multiplier = work[i][k] / pivot
            for j in range(k + 1, len(work)):
                work[i][j] -= multiplier * work[k][j]
    return value


def _minor(
    matrix: list[list[Fraction]], rows: tuple[int, ...], columns: tuple[int, ...]
) -> Fraction:
    return _det([[matrix[i][j] for j in columns] for i in rows])


def _kernel(size: int, case: int) -> list[list[Fraction]]:
    # Positive features and increasing positive nodes give strict total
    # positivity through every order <= size by Cauchy--Binet.
    return [
        [
            Fraction(2 ** ((i * 3 + case) % 7), 3 ** ((j * 2 + case) % 5))
            * sum(
                (
                    Fraction((i + 1) * (j + 2), size * (size + 1)) ** k
                    / math.factorial(k)
                    for k in range(size + 2)
                ),
                Fraction(),
            )
            for j in range(size)
        ]
        for i in range(size)
    ]


def test_fischer_and_column_maximal_pivot_exactly() -> None:
    checked = 0
    off_diagonal = 0
    for case in range(4):
        matrix = _kernel(5, case)
        for order in range(1, 6):
            tuples = tuple(itertools.combinations(range(5), order))
            for rows, columns in itertools.product(tuples, repeat=2):
                block = [[matrix[i][j] for j in columns] for i in rows]
                determinant = _det(block)
                assert determinant > 0
                for q in range(order):
                    complement = tuple(i for i in range(order) if i != q)
                    old_minor = _minor(block, complement, complement)
                    assert determinant <= block[q][q] * old_minor
                    p = max(range(order), key=lambda i: (block[i][q], -i))
                    pivot = block[p][q]
                    fresh_rows = tuple(i for i in range(order) if i != p)
                    residual = [
                        [
                            block[i][j] - block[i][q] * block[p][j] / pivot
                            for j in complement
                        ]
                        for i in fresh_rows
                    ]
                    assert abs(_det(residual)) == determinant / pivot
                    assert abs(_det(residual)) <= old_minor
                    checked += 1
                    off_diagonal += p != q
    assert checked == 2520
    assert off_diagonal > 1000


def test_local_minor_envelope_through_exact_prefixes() -> None:
    # This tests, but does NOT formally prove, the next all-prefix gauge lemma.
    # Maximize over ALL old row tuples, while holding each local column tuple
    # fixed; selecting just the pivot-deleted old row tuple is insufficient.
    checked = 0
    for case in range(4):
        matrix = _kernel(6, case)
        rows, columns = tuple(range(6)), tuple(range(6))
        row_sign, column_sign = [1] * 6, [1] * 6
        while len(rows) > 1:
            p, q = max(
                itertools.product(rows, columns),
                key=lambda ij: (abs(matrix[ij[0]][ij[1]]), -ij[0], -ij[1]),
            )
            pivot = matrix[p][q]
            fresh_rows = tuple(i for i in rows if i != p)
            fresh_columns = tuple(j for j in columns if j != q)
            residual = [
                [v - matrix[i][q] * matrix[p][j] / pivot for j, v in enumerate(line)]
                for i, line in enumerate(matrix)
            ]
            row_sign = [
                value * (-1 if i < p else 1) for i, value in enumerate(row_sign)
            ]
            column_sign = [
                value * (-1 if j < q else 1) for j, value in enumerate(column_sign)
            ]
            for order in range(1, min(4, len(fresh_rows)) + 1):
                for local_columns in itertools.combinations(fresh_columns, order):
                    old_max = max(
                        abs(_minor(matrix, r, local_columns))
                        for r in itertools.combinations(rows, order)
                    )
                    for local_rows in itertools.combinations(fresh_rows, order):
                        value = _minor(residual, local_rows, local_columns)
                        gauge = math.prod(row_sign[i] for i in local_rows) * math.prod(
                            column_sign[j] for j in local_columns
                        )
                        assert gauge * value > 0
                        assert abs(value) <= old_max
                        checked += 1
            rows, columns, matrix = fresh_rows, fresh_columns, residual
    assert checked == 1376
