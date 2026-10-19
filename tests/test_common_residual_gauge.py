"""Exact signed-frequency/boundary tests of one common residual sign gauge."""

import itertools
import math
from fractions import Fraction

import mpmath as mp

from kernelgecp.surrogates import fraction_determinant


def _vandermonde(values: tuple[Fraction, ...]) -> Fraction:
    return fraction_determinant([[x**k for k in range(len(values))] for x in values])


def test_vandermonde_border_factorization_without_sign_assumptions() -> None:
    samples = (
        (),
        (Fraction(0),),
        (Fraction(2), Fraction(-1)),
        (Fraction(1, 2), Fraction(3), Fraction(-2)),
        (Fraction(1), Fraction(1)),
    )
    for selected, fresh in itertools.product(samples, repeat=2):
        cross = math.prod(y - x for y in fresh for x in selected)
        assert _vandermonde(selected + fresh) == (
            _vandermonde(selected) * cross * _vandermonde(fresh)
        )


def _actual_rational_sample() -> list[list[Fraction]]:
    # Exact values of exp(t*y)/(1+exp(y)) at t=i/5 and y=5*j*log(2).
    # Time includes 0 and 1; frequency labels include negative, zero, positive.
    return [
        [Fraction(2) ** (i * j) / (1 + Fraction(2) ** (5 * j)) for j in range(-3, 3)]
        for i in range(6)
    ]


def _minor(matrix, rows, columns):
    return fraction_determinant([[matrix[i][j] for j in columns] for i in rows])


def _step(matrix, p, q):
    pivot = matrix[p][q]
    assert pivot
    return [
        [v - matrix[i][q] * matrix[p][j] / pivot for j, v in enumerate(line)]
        for i, line in enumerate(matrix)
    ]


def test_common_gauge_for_complete_and_noncomplete_prefixes() -> None:
    checked = 0
    noncomplete = 0
    for mode in ("complete", "smallest", "interior"):
        matrix = _actual_rational_sample()
        rows, columns = tuple(range(6)), tuple(range(6))
        sr, sc = [], []
        while True:
            row_sign = {i: math.prod(-1 if i < p else 1 for p in sr) for i in rows}
            col_sign = {j: math.prod(-1 if j < q else 1 for q in sc) for j in columns}
            gauged = [
                [
                    row_sign.get(i, 0) * v * col_sign.get(j, 0)
                    for j, v in enumerate(line)
                ]
                for i, line in enumerate(matrix)
            ]
            for order in range(len(rows) + 1):
                for r, c in itertools.product(
                    itertools.combinations(rows, order),
                    itertools.combinations(columns, order),
                ):
                    raw = _minor(matrix, r, c)
                    transformed = _minor(gauged, r, c)
                    assert transformed > 0
                    assert transformed == abs(raw)
                    checked += 1
                    # The absolute identity also covers tuples not in increasing
                    # order, unlike positivity itself.
                    assert abs(_minor(gauged, r[::-1], c)) == abs(raw)
            if not rows:
                break
            candidates = tuple(itertools.product(rows, columns))
            cp = max(
                candidates, key=lambda ij: (abs(matrix[ij[0]][ij[1]]), -ij[0], -ij[1])
            )
            if mode == "complete":
                p, q = cp
            elif mode == "smallest":
                p, q = min(candidates, key=lambda ij: (abs(matrix[ij[0]][ij[1]]), ij))
            else:
                p, q = rows[len(rows) // 2], columns[len(columns) // 2]
            noncomplete += (p, q) != cp
            matrix = _step(matrix, p, q)
            sr.append(p)
            sc.append(q)
            rows = tuple(i for i in rows if i != p)
            columns = tuple(j for j in columns if j != q)
    assert checked == 3825
    assert noncomplete >= 8


def test_gauged_and_actual_schur_steps_have_same_absolute_minors() -> None:
    checked = 0
    matrix = _actual_rational_sample()
    rows, columns = tuple(range(6)), tuple(range(6))
    sr, sc = [], []
    while len(rows) > 1:
        row_sign = {i: math.prod(-1 if i < p else 1 for p in sr) for i in rows}
        col_sign = {j: math.prod(-1 if j < q else 1 for q in sc) for j in columns}
        gauged = [
            [row_sign.get(i, 0) * v * col_sign.get(j, 0) for j, v in enumerate(line)]
            for i, line in enumerate(matrix)
        ]
        p, q = max(
            itertools.product(rows, columns),
            key=lambda ij: (abs(matrix[ij[0]][ij[1]]), -ij[0], -ij[1]),
        )
        actual_next, gauged_next = _step(matrix, p, q), _step(gauged, p, q)
        rr, cc = tuple(i for i in rows if i != p), tuple(j for j in columns if j != q)
        for order in range(len(rr) + 1):
            for r, c in itertools.product(
                itertools.combinations(rr, order), itertools.combinations(cc, order)
            ):
                assert abs(_minor(actual_next, r, c)) == abs(_minor(gauged_next, r, c))
                # In augmented ordered positions, delete the row having the
                # SAME RANK as the pivot column, not necessarily the pivot row.
                augmented_rows = tuple(sorted(r + (p,)))
                augmented_columns = tuple(sorted(c + (q,)))
                matched_row = augmented_rows[augmented_columns.index(q)]
                old_rows = tuple(i for i in augmented_rows if i != matched_row)
                assert abs(_minor(actual_next, r, c)) <= abs(
                    _minor(matrix, old_rows, c)
                )
                checked += 1
        matrix, rows, columns = actual_next, rr, cc
        sr.append(p)
        sc.append(q)
    assert checked == 350


def test_common_gauge_absolute_identity_with_repeated_coordinates() -> None:
    matrix = _actual_rational_sample()
    matrix = _step(matrix, 2, 3)
    rows, columns = (0, 0, 5), (1, 4, 4)
    raw = _minor(matrix, rows, columns)
    assert raw == 0
    gauged = [
        [(-1 if i < 2 else 1) * v * (-1 if j < 3 else 1) for j, v in enumerate(line)]
        for i, line in enumerate(matrix)
    ]
    assert abs(_minor(gauged, rows, columns)) == abs(raw)


def test_actual_exponential_signed_boundary_sample_at_independent_precisions() -> None:
    exact = _actual_rational_sample()
    selected = ((2, 3), (4, 1), (0, 5))  # Deliberately not a complete-pivot rule.
    for p, q in selected:
        exact = _step(exact, p, q)
    records = []
    for digits in (160, 320):
        with mp.workdps(digits):
            matrix = [
                [
                    mp.exp(mp.mpf(i) / 5 * (5 * j * mp.log(2)))
                    / (1 + mp.exp(5 * j * mp.log(2)))
                    for j in range(-3, 3)
                ]
                for i in range(6)
            ]
            for p, q in selected:
                matrix = _step(matrix, p, q)
            rows, columns = (1, 3, 5), (0, 2, 4)
            values = []
            for order in (1, 2, 3):
                for r, c in itertools.product(
                    itertools.combinations(rows, order),
                    itertools.combinations(columns, order),
                ):
                    observed = mp.det(mp.matrix([[matrix[i][j] for j in c] for i in r]))
                    value = _minor(exact, r, c)
                    reference = mp.mpf(value.numerator) / value.denominator
                    assert abs((observed - reference) / reference) < mp.mpf("1e-140")
                    values.append(mp.nstr(observed, 45))
            records.append(values)
    assert records[0] == records[1]


def test_frequency_reflection_transports_the_original_run_exactly() -> None:
    reflected = _actual_rational_sample()
    original = [line[::-1] for line in reflected]
    # Original frequencies are increasing; the reflected list has reversed
    # indices. Transport the CHOSEN original pivots, not a new tie policy.
    rows, columns = tuple(range(6)), tuple(range(6))
    for _ in range(6):
        p, q = max(
            itertools.product(rows, columns),
            key=lambda ij: (abs(original[ij[0]][ij[1]]), -ij[0], -ij[1]),
        )
        assert original[p][q] == reflected[p][5 - q]
        original = _step(original, p, q)
        reflected = _step(reflected, p, 5 - q)
        assert original == [line[::-1] for line in reflected]
        rows = tuple(i for i in rows if i != p)
        columns = tuple(j for j in columns if j != q)
