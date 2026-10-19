"""Independent numerical checks of the analytic one-pivot sign theorem.

The decimal comparisons are evidence, not proof or runtime certification.
The depth-ten regression separately rejects a proposed multi-pivot shortcut.
"""

import itertools
from fractions import Fraction

import mpmath as mp

from experiments.transfer_nodes import (
    conditioned_factors,
    cutoff_coefficients,
    matched_coordinates,
    transfer_numerator,
)


def test_actual_exponential_one_pivot_sign_schedule() -> None:
    records: dict[int, tuple[str, ...]] = {}
    triples = ((0, 1, 2, 2), (1, 2, 3, 4), (97, 98, 99, 100))
    for digits in (80, 160):
        values: list[str] = []
        with mp.workdps(digits):
            for rate_numerator, rate_denominator in ((1, 100), (1, 2), (1, 1)):
                rate = mp.mpf(rate_numerator) / rate_denominator
                for nx, ny, nz, denominator in triples:
                    nodes = [mp.mpf(t) / denominator for t in (nx, ny, nz)]
                    for i, j in itertools.combinations(range(9), 2):
                        matrix = mp.matrix(
                            [[mp.exp(rate * t), t**i, t**j] for t in nodes]
                        )
                        raw = mp.det(matrix)
                        normalized = mp.matrix(
                            [
                                [1, t**i * mp.exp(-rate * t), t**j * mp.exp(-rate * t)]
                                for t in nodes
                            ]
                        )
                        assert mp.almosteq(
                            raw,
                            mp.exp(rate * sum(nodes)) * mp.det(normalized),
                            rel_eps=mp.mpf("1e-60"),
                        )
                        sign = -1 if i == 0 and j >= 2 else 1
                        assert sign * raw > 0
                        values.append(mp.nstr(raw, 35))
        records[digits] = tuple(values)
    assert len(records[80]) == 324
    assert records[80] == records[160]


def test_independent_one_pivot_augmented_minor_products() -> None:
    with mp.workdps(100):
        for i, j in itertools.combinations(range(12), 2):
            determinants = []
            for rate, nodes in (
                (mp.mpf(1), (mp.mpf(0), mp.mpf("0.5"), mp.mpf(1))),
                (mp.mpf("0.3"), (mp.mpf("0.11"), mp.mpf("0.12"), mp.mpf("0.99"))),
            ):
                # Direct elimination is independent of the derivative-ratio proof.
                determinants.append(
                    mp.det(mp.matrix([[mp.exp(rate * t), t**i, t**j] for t in nodes]))
                )
            assert determinants[0] * determinants[1] > 0


def test_bad_node_need_not_be_pivot_depth_plus_one_exactly() -> None:
    rows = tuple(Fraction(t, 34) for t in (3, 4, 7, 12, 13, 14, 18, 20, 22, 25, 27, 31))
    columns = tuple(
        Fraction(t, 33) for t in (1, 2, 3, 6, 11, 14, 16, 19, 21, 26, 28, 31)
    )
    depth, features = 10, 23
    left, right, residual, selected_rows, selected_columns = conditioned_factors(
        rows, columns, features, depth, features
    )
    assert selected_rows == [11, 0, 6, 9, 2, 3, 10, 1, 8, 5]
    assert selected_columns == [11, 0, 7, 4, 9, 3, 10, 2, 6, 8]
    fresh_rows, fresh_columns = (4, 7), (1, 5)
    # All retained features reconstruct the same kernel used for pivoting.
    for i, j in itertools.product(fresh_rows, fresh_columns):
        assert (
            sum((left[i][k] * right[k][j] for k in range(features)), Fraction())
            == residual[i][j]
        )
    orientation = (-1) ** (
        sum(i > j for i in selected_rows for j in fresh_rows)
        + sum(i > j for i in selected_columns for j in fresh_columns)
    )
    coordinates = matched_coordinates(
        left, right, fresh_rows, fresh_columns, orientation
    )
    numerator, l_transfer, r_transfer = transfer_numerator(coordinates, features, 10)
    assert depth + 1 != 10
    assert numerator > 0 and l_transfer > 0 and r_transfer > 0
    constant, mixed, quadratic = cutoff_coefficients(coordinates, features, 10)
    scale = (
        constant * r_transfer**2
        + mixed * l_transfer * r_transfer
        + quadratic * l_transfer**2
    )
    assert Fraction(24, 10000) < numerator / scale < Fraction(25, 10000)
    assert transfer_numerator(coordinates, features, 9)[0] < 0
    assert transfer_numerator(coordinates, features, 11)[0] < 0
    for cutoff in range(2, features - 1):
        constant, mixed, quadratic = cutoff_coefficients(coordinates, features, cutoff)
        assert constant > 0 and mixed > 0 and quadratic > 0
        assert mixed**2 > 4 * constant * quadratic
        assert (
            min(
                transfer_numerator(coordinates, features, cutoff - 1)[0],
                transfer_numerator(coordinates, features, cutoff)[0],
            )
            <= 0
        )
