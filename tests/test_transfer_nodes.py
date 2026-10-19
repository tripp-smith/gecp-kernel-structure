import itertools
from fractions import Fraction

import mpmath as mp
import pytest

from experiments.transfer_nodes import (
    audit,
    conditioned_factors,
    cutoff_coefficients,
    determinant_two,
    matched_coordinates,
    transfer_numerator,
)


def test_interior_node_failure_and_alternate_certificates_exactly() -> None:
    rows = tuple(Fraction(value, 21) for value in (17, 18, 19, 20))
    columns = tuple(Fraction(value, 20) for value in (15, 16, 17, 19))
    features = 9
    left, right, residual, selected_rows, selected_columns = conditioned_factors(
        rows, columns, features, 2, features
    )
    assert selected_rows == selected_columns == [3, 0]
    coordinates = matched_coordinates(left, right, (1, 2), (1, 2), 1)
    for i, j in itertools.product((1, 2), repeat=2):
        assert (
            sum(
                (left[i][degree] * right[degree][j] for degree in range(features)),
                Fraction(),
            )
            == residual[i][j]
        )

    def block(degrees: range) -> list[list[Fraction]]:
        return [
            [
                sum((left[i][k] * right[k][j] for k in degrees), Fraction())
                for j in (1, 2)
            ]
            for i in (1, 2)
        ]

    for cutoff in range(1, features):
        low = block(range(cutoff))
        high = block(range(cutoff, features))
        constant, mixed, quadratic = cutoff_coefficients(coordinates, features, cutoff)
        assert determinant_two(low) == constant
        assert determinant_two(high) == quadratic
        assert mixed == (
            determinant_two(
                [[low[i][j] + high[i][j] for j in range(2)] for i in range(2)]
            )
            - constant
            - quadratic
        )
        assert constant >= 0 and mixed > 0 and quadratic >= 0
        assert mixed**2 > 4 * constant * quadratic
        if cutoff in (3, 4):
            assert constant > 0 and quadratic > 0

    numerator, left_transfer, right_transfer = transfer_numerator(
        coordinates, features, 3
    )
    assert left_transfer > 0 and right_transfer > 0
    assert numerator > 0
    constant, mixed, quadratic = cutoff_coefficients(coordinates, features, 3)
    scale = (
        constant * right_transfer**2
        + mixed * left_transfer * right_transfer
        + quadratic * left_transfer**2
    )
    assert Fraction(9, 100000) < numerator / scale < Fraction(1, 10000)
    assert transfer_numerator(coordinates, features, 2)[0] < 0
    assert transfer_numerator(coordinates, features, 4)[0] < 0


def test_adaptive_transfer_coverage_exactly() -> None:
    report = audit(cases=100)
    assert report["coverage_holds_on_samples"] is True
    assert report["transfers"] == 54845
    assert report["interior_cutoffs_covered"] == 47496
    assert report["positive_node_numerators"] == 4
    assert report["stop_reason"] == "completed"
    with pytest.raises(ValueError, match="positive"):
        audit(cases=0)


def test_interior_node_failure_persists_for_exponential_numerically() -> None:
    records: dict[int, tuple[str, ...]] = {}
    for digits in (80, 160):
        with mp.workdps(digits):
            rows = [mp.mpf(v) / 21 for v in (17, 18, 19, 20)]
            columns = [mp.mpf(v) / 20 for v in (15, 16, 17, 19)]
            kernel = mp.matrix([[mp.exp(x * y) for y in columns] for x in rows])
            left = mp.matrix(
                [
                    [x**degree / mp.factorial(degree) for degree in range(5)]
                    for x in rows
                ]
            )
            right = mp.matrix([[y**degree for y in columns] for degree in range(5)])
            selected_rows: list[int] = []
            selected_columns: list[int] = []
            for expected in ((3, 3), (0, 0)):
                _, _, _, row, column = max(
                    (abs(kernel[i, j]), -i, -j, i, j)
                    for i in range(4)
                    for j in range(4)
                    if i not in selected_rows and j not in selected_columns
                )
                assert (row, column) == expected
                pivot = kernel[row, column]
                left = mp.matrix(
                    [
                        [
                            left[i, k] - kernel[i, column] * left[row, k] / pivot
                            for k in range(5)
                        ]
                        for i in range(4)
                    ]
                )
                right = mp.matrix(
                    [
                        [
                            right[k, j] - right[k, column] * kernel[row, j] / pivot
                            for j in range(4)
                        ]
                        for k in range(5)
                    ]
                )
                kernel = mp.matrix(
                    [
                        [
                            kernel[i, j] - kernel[i, column] * kernel[row, j] / pivot
                            for j in range(4)
                        ]
                        for i in range(4)
                    ]
                )
                selected_rows.append(row)
                selected_columns.append(column)

            residual = mp.matrix([[kernel[i, j] for j in (1, 2)] for i in (1, 2)])

            def block(
                degrees: range,
                left_factor: mp.matrix = left,
                right_factor: mp.matrix = right,
            ) -> mp.matrix:
                return mp.matrix(
                    [
                        [
                            sum(left_factor[i, k] * right_factor[k, j] for k in degrees)
                            for j in (1, 2)
                        ]
                        for i in (1, 2)
                    ]
                )

            def mixed(first: mp.matrix, second: mp.matrix) -> mp.mpf:
                return (
                    first[0, 0] * second[1, 1]
                    + second[0, 0] * first[1, 1]
                    - first[0, 1] * second[1, 0]
                    - second[0, 1] * first[1, 0]
                )

            values: list[str] = []
            for degree in (2, 3, 4):
                low = block(range(degree))
                moved = block(range(degree, degree + 1))
                upper = residual - low - moved
                constant = mp.det(low)
                linear = mixed(low, moved + upper)
                quadratic = mp.det(moved + upper)
                lsum = mixed(low, moved)
                rsum = mixed(moved, upper)
                numerator = (
                    constant * rsum**2 - linear * lsum * rsum + quadratic * lsum**2
                )
                scale = (
                    abs(constant * rsum**2)
                    + abs(linear * lsum * rsum)
                    + abs(quadratic * lsum**2)
                )
                normalized = numerator / scale
                if degree == 3:
                    assert mp.mpf("0.0001") < normalized < mp.mpf("0.00011")
                else:
                    assert normalized < 0
                values.append(mp.nstr(normalized, 35))
            records[digits] = tuple(values)
    # Precision agreement is numerical evidence, separate from the finite
    # rational obstruction formally proved in Lean.
    assert records[80] == records[160]
