"""Exact coverage obstructions and separate actual-exponential evidence."""

import itertools
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from experiments.transfer_nodes import (
    conditioned_factors,
    cutoff_coefficients,
    determinant_two,
    matched_coordinates,
    transfer_numerator,
)
from scripts.generate_coverage_certificate import generate_certificate

ROWS = tuple(Fraction(v, 29) for v in (7, 11, 12, 18, 19, 21, 22, 23, 27))
COLUMNS = tuple(Fraction(v, 63) for v in (8, 9, 11, 13, 24, 33, 39, 57, 60))


def _mixed(a: list[list[Fraction]], b: list[list[Fraction]]) -> Fraction:
    return a[0][0] * b[1][1] + b[0][0] * a[1][1] - a[0][1] * b[1][0] - b[0][1] * a[1][0]


def test_integer_certificate_is_reproducible() -> None:
    source, report = generate_certificate()
    root = Path(__file__).resolve().parents[1]
    target = root / "GECPKernelStructure/Fermionic/AdaptiveCoverageObstructionData.lean"
    assert target.read_text() == source
    assert report["selected_rows"] == [8, 0, 3, 2, 7, 5, 1]
    assert report["selected_columns"] == [8, 0, 6, 4, 7, 3, 5]
    assert report["feature_scale_decimal_digits"] == 305
    assert report["stop_reason"] == "completed"


def test_coherent_unit_square_neighboring_coverage_fails_exactly() -> None:
    features = 13
    left, right, residual, sr, sc = conditioned_factors(
        ROWS, COLUMNS, features, 7, features
    )
    assert sr == [8, 0, 3, 2, 7, 5, 1]
    assert sc == [8, 0, 6, 4, 7, 3, 5]
    rows, columns = (6, 4), (1, 2)

    def block(degrees: range) -> list[list[Fraction]]:
        return [
            [
                sum((left[i][k] * right[k][j] for k in degrees), Fraction())
                for j in columns
            ]
            for i in rows
        ]

    for i, j in itertools.product(rows, columns):
        assert (
            sum((left[i][k] * right[k][j] for k in range(features)), Fraction())
            == residual[i][j]
        )
    coordinates = matched_coordinates(left, right, rows, columns, 1)
    for cutoff in range(2, features - 1):
        a, b, c = cutoff_coefficients(coordinates, features, cutoff)
        assert min(a, b, c) > 0
        assert b * b > 4 * a * c
    for moved in (7, 8):
        low, singleton, upper = (
            block(range(moved)),
            block(range(moved, moved + 1)),
            block(range(moved + 1, features)),
        )
        high = [[singleton[i][j] + upper[i][j] for j in range(2)] for i in range(2)]
        l_transfer, r_transfer = _mixed(low, singleton), _mixed(singleton, upper)
        numerator = (
            determinant_two(low) * r_transfer**2
            - _mixed(low, high) * l_transfer * r_transfer
            + determinant_two(high) * l_transfer**2
        )
        assert (numerator, l_transfer, r_transfer) == transfer_numerator(
            coordinates, features, moved
        )
        assert numerator > 0 and l_transfer > 0 and r_transfer > 0
        for matrix in (low, high, upper):
            assert determinant_two(matrix) > 0


def test_outer_prefix_loses_positive_transfer_sums_exactly() -> None:
    rows = tuple(Fraction(v, 13) for v in (4, 6, 8, 12))
    columns = (Fraction(1, 5), Fraction(2, 5), Fraction(13, 15), Fraction(22, 3))
    features = 12
    denominators = [
        1 + sum((y**k / math.factorial(k) for k in range(features)), Fraction())
        for y in columns
    ]
    left = [[x**k / math.factorial(k) for k in range(features)] for x in rows]
    right = [
        [y**k / denominators[j] for j, y in enumerate(columns)] for k in range(features)
    ]
    matrix = [
        [
            sum((left[i][k] * right[k][j] for k in range(features)), Fraction())
            for j in range(4)
        ]
        for i in range(4)
    ]
    sr: list[int] = []
    sc: list[int] = []
    for expected in ((3, 2), (2, 3)):
        _, _, _, row, column = max(
            (abs(matrix[i][j]), -i, -j, i, j)
            for i in range(4)
            for j in range(4)
            if i not in sr and j not in sc
        )
        assert (row, column) == expected
        pivot = matrix[row][column]
        assert pivot != 0
        left = [
            [
                value - matrix[i][column] * left[row][k] / pivot
                for k, value in enumerate(line)
            ]
            for i, line in enumerate(left)
        ]
        right = [
            [
                value - right[k][column] * matrix[row][j] / pivot
                for j, value in enumerate(line)
            ]
            for k, line in enumerate(right)
        ]
        matrix = [
            [
                value - matrix[i][column] * matrix[row][j] / pivot
                for j, value in enumerate(line)
            ]
            for i, line in enumerate(matrix)
        ]
        sr.append(row)
        sc.append(column)
    for i, j in itertools.product((0, 1), repeat=2):
        assert (
            sum((left[i][k] * right[k][j] for k in range(features)), Fraction())
            == matrix[i][j]
        )
    coordinates = matched_coordinates(left, right, (0, 1), (0, 1), 1)
    for moved in (5, 6):
        numerator, l_transfer, r_transfer = transfer_numerator(
            coordinates, features, moved
        )
        assert numerator > 0 and l_transfer < 0 and r_transfer > 0
    a, b, c = cutoff_coefficients(coordinates, features, 6)
    assert min(a, b, c) > 0 and b * b > 4 * a * c


def _exponential_record(
    rows: tuple[Fraction, ...],
    columns: tuple[Fraction, ...],
    prefix: tuple[tuple[int, int], ...],
    cutoffs: tuple[int, ...],
    weighted: bool,
) -> tuple[str, ...]:
    xs = [mp.mpf(x.numerator) / x.denominator for x in rows]
    ys = [mp.mpf(y.numerator) / y.denominator for y in columns]
    weights = [1 / (1 + mp.exp(y)) if weighted else mp.mpf(1) for y in ys]
    size, terms = len(rows), max(cutoffs) + 1
    left = mp.matrix([[x**k / mp.factorial(k) for k in range(terms)] for x in xs])
    right = mp.matrix(
        [[y**k * weights[j] for j, y in enumerate(ys)] for k in range(terms)]
    )
    kernel = mp.matrix(
        [[mp.exp(x * y) * weights[j] for j, y in enumerate(ys)] for x in xs]
    )
    sr: list[int] = []
    sc: list[int] = []
    for expected in prefix:
        candidates = sorted(
            (abs(kernel[i, j]), -i, -j, i, j)
            for i in range(size)
            for j in range(size)
            if i not in sr and j not in sc
        )
        peak, _, _, row, column = candidates[-1]
        assert (row, column) == expected
        assert (peak - candidates[-2][0]) / peak > mp.mpf("1e-12")
        pivot = kernel[row, column]
        assert pivot != 0
        left = mp.matrix(
            [
                [
                    left[i, k] - kernel[i, column] * left[row, k] / pivot
                    for k in range(terms)
                ]
                for i in range(size)
            ]
        )
        right = mp.matrix(
            [
                [
                    right[k, j] - right[k, column] * kernel[row, j] / pivot
                    for j in range(size)
                ]
                for k in range(terms)
            ]
        )
        kernel = mp.matrix(
            [
                [
                    kernel[i, j] - kernel[i, column] * kernel[row, j] / pivot
                    for j in range(size)
                ]
                for i in range(size)
            ]
        )
        sr.append(row)
        sc.append(column)
    fresh_rows = [i for i in range(size) if i not in sr]
    fresh_columns = [j for j in range(size) if j not in sc]
    orientation = (-1) ** (
        sum(i > j for i in sr for j in fresh_rows)
        + sum(i > j for i in sc for j in fresh_columns)
    )
    residual = mp.matrix([[kernel[i, j] for j in fresh_columns] for i in fresh_rows])

    def mixed(a: mp.matrix, b: mp.matrix) -> mp.mpf:
        return (
            a[0, 0] * b[1, 1]
            + b[0, 0] * a[1, 1]
            - a[0, 1] * b[1, 0]
            - b[0, 1] * a[1, 0]
        )

    record: list[str] = []
    for moved in cutoffs:
        low = mp.matrix(
            [
                [
                    sum(left[i, k] * right[k, j] for k in range(moved))
                    for j in fresh_columns
                ]
                for i in fresh_rows
            ]
        )
        singleton = mp.matrix(
            [
                [left[i, moved] * right[moved, j] for j in fresh_columns]
                for i in fresh_rows
            ]
        )
        high = (
            residual - low
        )  # The full exponential tail, not a finite feature truncation.
        a, b, c = (
            orientation * v for v in (mp.det(low), mixed(low, high), mp.det(high))
        )
        l_transfer, r_transfer = (
            orientation * v
            for v in (mixed(low, singleton), mixed(singleton, high - singleton))
        )
        assert min(a, b, c, r_transfer) > 0
        assert l_transfer < 0 if weighted else l_transfer > 0
        numerator = a * r_transfer**2 - b * l_transfer * r_transfer + c * l_transfer**2
        scale = (
            abs(a * r_transfer**2)
            + abs(b * l_transfer * r_transfer)
            + abs(c * l_transfer**2)
        )
        assert numerator > 0 and b * b > 4 * a * c
        record.extend(
            (
                mp.nstr(numerator / scale, 45),
                mp.nstr((b * b - 4 * a * c) / (b * b + 4 * a * c), 45),
            )
        )
    return tuple(record)


def test_actual_exponential_coverage_failures_numerically() -> None:
    records: dict[int, tuple[str, ...]] = {}
    for digits in (160, 320):
        with mp.workdps(digits):
            unit = _exponential_record(
                ROWS,
                COLUMNS,
                ((8, 8), (0, 0), (3, 6), (2, 4), (7, 7), (5, 3), (1, 5)),
                (7, 8),
                False,
            )
            outer = _exponential_record(
                tuple(Fraction(v, 13) for v in (4, 6, 8, 12)),
                (Fraction(1, 5), Fraction(2, 5), Fraction(13, 15), Fraction(22, 3)),
                ((3, 2), (2, 3)),
                (5, 6),
                True,
            )
            records[digits] = unit + outer
    assert records[160] == records[320]
    assert records[160][0] == "0.00733642663708844273392090806344167651521625274"
    assert records[160][2] == "0.0154696286973923495855634025862919067083240956"
