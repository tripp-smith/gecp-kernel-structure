"""Exact rational checks for geometric sign-regular surrogate matrices."""

from __future__ import annotations

import itertools
from dataclasses import asdict, dataclass
from fractions import Fraction

RationalMatrix = list[list[Fraction]]


@dataclass(frozen=True, slots=True)
class SurrogateRecord:
    size: int
    q: str
    pivot_rows: list[int]
    pivot_columns: list[int]
    pivots: list[str]
    minor_counts: dict[str, int]
    minor_sign_mismatches: dict[str, int]
    zero_minors: dict[str, int]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def geometric_surrogate(size: int, q: Fraction) -> RationalMatrix:
    """Return the exact matrix with entry ``q ** (i * j)``."""

    if size < 1 or not 0 < q < 1:
        raise ValueError("size must be positive and q must lie in (0, 1)")
    return [[q ** (i * j) for j in range(size)] for i in range(size)]


def fraction_determinant(matrix: RationalMatrix) -> Fraction:
    """Compute an exact determinant by fraction-preserving elimination."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    if size == 0:
        return Fraction(1)
    work = [row.copy() for row in matrix]
    determinant = Fraction(1)
    for column in range(size):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column] != 0), None
        )
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            determinant = -determinant
        pivot = work[column][column]
        determinant *= pivot
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot
            for entry in range(column + 1, size):
                work[row][entry] -= multiplier * work[column][entry]
    return determinant


def exact_gecp(matrix: RationalMatrix) -> tuple[list[int], list[int], list[Fraction]]:
    """Run lexicographic complete pivoting using exact ``Fraction`` arithmetic."""

    rows = len(matrix)
    columns = len(matrix[0]) if rows else 0
    if rows == 0 or any(len(row) != columns for row in matrix):
        raise ValueError("matrix must be nonempty and rectangular")
    residual = [row.copy() for row in matrix]
    pivot_rows: list[int] = []
    pivot_columns: list[int] = []
    pivots: list[Fraction] = []
    for _ in range(min(rows, columns)):
        row, column = max(
            itertools.product(range(rows), range(columns)),
            key=lambda coordinate: (
                abs(residual[coordinate[0]][coordinate[1]]),
                -coordinate[0],
                -coordinate[1],
            ),
        )
        pivot = residual[row][column]
        if pivot == 0:
            break
        pivot_rows.append(row)
        pivot_columns.append(column)
        pivots.append(pivot)
        pivot_column = [residual[i][column] for i in range(rows)]
        pivot_row_values = residual[row].copy()
        for i in range(rows):
            for j in range(columns):
                residual[i][j] -= pivot_column[i] * pivot_row_values[j] / pivot
    return pivot_rows, pivot_columns, pivots


def inspect_surrogate(size: int, q: Fraction) -> SurrogateRecord:
    """Enumerate every square minor and the exact GECP path."""

    matrix = geometric_surrogate(size, q)
    counts: dict[str, int] = {}
    mismatches: dict[str, int] = {}
    zeros: dict[str, int] = {}
    for order in range(1, size + 1):
        key = str(order)
        counts[key] = 0
        mismatches[key] = 0
        zeros[key] = 0
        expected_sign = -1 if (order * (order - 1) // 2) % 2 else 1
        for row_indices in itertools.combinations(range(size), order):
            for column_indices in itertools.combinations(range(size), order):
                minor = [[matrix[i][j] for j in column_indices] for i in row_indices]
                determinant = fraction_determinant(minor)
                counts[key] += 1
                if determinant == 0:
                    zeros[key] += 1
                elif (1 if determinant > 0 else -1) != expected_sign:
                    mismatches[key] += 1
    rows, columns, pivots = exact_gecp(matrix)
    return SurrogateRecord(
        size=size,
        q=f"{q.numerator}/{q.denominator}",
        pivot_rows=rows,
        pivot_columns=columns,
        pivots=[f"{pivot.numerator}/{pivot.denominator}" for pivot in pivots],
        minor_counts=counts,
        minor_sign_mismatches=mismatches,
        zero_minors=zeros,
    )


def exact_surrogate_census(
    sizes: range = range(2, 9),
    q_values: tuple[Fraction, ...] = (
        Fraction(1, 2),
        Fraction(2, 3),
        Fraction(3, 4),
    ),
) -> list[SurrogateRecord]:
    """Run the fixed Milestone E exact-surrogate protocol."""

    return [inspect_surrogate(size, q) for size in sizes for q in q_values]
