"""Independent finite witnesses; the continuous rate itself is proved in Lean."""

import itertools
import math
from fractions import Fraction

import mpmath as mp

from kernelgecp.surrogates import fraction_determinant


def _sample():
    # Actual positive presentation at t=i/5, y=5*j*log(2), including endpoints.
    return [
        [Fraction(2) ** (i * j) / (1 + Fraction(2) ** (5 * j)) for j in range(-3, 3)]
        for i in range(6)
    ]


def _minor(matrix, rows, columns):
    return fraction_determinant([[matrix[i][j] for j in columns] for i in rows])


def _trace():
    matrix = _sample()
    states, pivots = [matrix], []
    rows = columns = tuple(range(6))
    while rows:
        p, q = max(
            itertools.product(rows, columns),
            key=lambda ij: (abs(matrix[ij[0]][ij[1]]), -ij[0], -ij[1]),
        )
        value = matrix[p][q]
        assert value
        pivots.append((p, q, value))
        matrix = [
            [v - matrix[i][q] * matrix[p][j] / value for j, v in enumerate(line)]
            for i, line in enumerate(matrix)
        ]
        states.append(matrix)
        rows = tuple(i for i in rows if i != p)
        columns = tuple(j for j in columns if j != q)
    return states, pivots


def _witness(states, pivots, marked):
    rows = columns = ()
    product = Fraction(1)
    inside_unmarked = rank_changes = 0
    # Include some unmarked pivot columns in the band: ignoring an in-band
    # pivot is essential when selecting exactly the budget order of marks.
    band = {pivots[i][1] for i in marked} | {pivots[0][1]}
    for k in reversed(range(len(pivots))):
        p, q, pivot = pivots[k]
        new = abs(_minor(states[k + 1], rows, columns))
        assert new > 0
        if k in marked:
            rows, columns = (p,) + rows, (q,) + columns
            assert abs(_minor(states[k], rows, columns)) == abs(pivot) * new
            product *= abs(pivot)
        else:
            ar, ac = sorted((p,) + rows), sorted((q,) + columns)
            rank = ac.index(q)
            rank_changes += ar[rank] != p
            rows = tuple(i for j, i in enumerate(ar) if j != rank)
            columns = tuple(i for j, i in enumerate(ac) if j != rank)
            assert abs(_minor(states[k], rows, columns)) >= new
            inside_unmarked += q in band
        assert len(rows) == len(columns) == len({i for i in marked if i >= k})
        assert set(columns) == {pivots[i][1] for i in marked if i >= k}
        assert set(columns) <= band
        assert not set(rows).intersection(p for p, _, _ in pivots[:k])
        assert not set(columns).intersection(q for _, q, _ in pivots[:k])
        assert abs(_minor(states[k], rows, columns)) >= product
    return rows, columns, product, inside_unmarked, rank_changes


def test_backwards_marked_witness_for_every_mark_set_exactly() -> None:
    states, pivots = _trace()
    inside = changed = 0
    for mask in range(64):
        marked = {k for k in range(6) if mask & (1 << k)}
        rows, columns, product, count, changes = _witness(states, pivots, marked)
        assert len(rows) == len(columns) == len(marked)
        assert abs(_minor(_sample(), rows, columns)) >= product
        inside += count
        changed += changes
    assert 0 < inside < 192  # Both in-band and outside-band unmarked steps.
    assert changed > 0


def test_panel_count_from_independently_enumerated_initial_minors() -> None:
    states, pivots = _trace()
    checked = 0
    for width in (2, 3, 4):
        for band in itertools.combinations(range(6), width):
            for order in range(1, min(3, width) + 1):
                maximum = max(
                    abs(_minor(states[0], r, c))
                    for r in itertools.combinations(range(6), order)
                    for c in itertools.combinations(band, order)
                )
                # Choose the smallest dyadic threshold whose power bounds
                # every initial local minor. No residual-budget assumption.
                epsilon = Fraction(1, 2**32)
                while epsilon**order < maximum:
                    epsilon *= 2
                large = [
                    k
                    for k, (_, q, v) in enumerate(pivots)
                    if q in band and abs(v) > epsilon
                ]
                assert len(large) < order
                checked += 1
    assert checked == 135


def test_actual_exponential_witnesses_at_two_precisions() -> None:
    states, pivots = _trace()
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            matrix = [
                [
                    mp.exp(mp.mpf(i) / 5 * (5 * j * mp.log(2)))
                    / (1 + mp.exp(5 * j * mp.log(2)))
                    for j in range(-3, 3)
                ]
                for i in range(6)
            ]
            for marked in ({0, 2, 4}, {1, 3}, set(range(6))):
                r, c, product, _, _ = _witness(states, pivots, marked)
                determinant = abs(
                    mp.det(mp.matrix([[matrix[i][j] for j in c] for i in r]))
                )
                exact = abs(_minor(states[0], r, c))
                expected = mp.mpf(exact.numerator) / exact.denominator
                assert mp.isfinite(determinant) and determinant > 0, (
                    "Uncertified determinant"
                )
                assert mp.almosteq(
                    determinant, expected, rel_eps=mp.mpf(10) ** (-digits + 25)
                )
                assert expected >= mp.mpf(product.numerator) / product.denominator
                values.append(mp.nstr(determinant, 45))
        records.append(values)
    assert records[0] == records[1]


def test_finite_panel_partition_and_global_block_arithmetic() -> None:
    for s in range(12):
        labels = tuple(itertools.product((False, True), range(s + 1)))
        nodes = {Fraction(0), Fraction(1, 2), Fraction(-1, 2)}
        for j in range(s + 1):
            for sign in (-1, 1):
                nodes.update((Fraction(sign * 2**j), Fraction(sign * 3 * 2**j, 4)))
        for w in nodes:
            if abs(w) > 2**s:
                continue
            membership = []
            for negative, j in labels:
                u = -w if negative else w
                member = (not negative or w < 0) and (
                    0 <= u <= 1 if j == 0 else 2 ** (j - 1) < u <= 2**j
                )
                if member:
                    membership.append((negative, j))
            assert len(membership) == 1
        assert len(labels) == 2 * (s + 1)
        for q in range(65):
            order = 16 * (33 + 4 * q)
            count = 32 * (s + 1) * (33 + 4 * q)
            assert count == len(labels) * order
            # If every panel has < order large pivots, there are < count.
            assert math.prod((len(labels), order - 1)) < count
