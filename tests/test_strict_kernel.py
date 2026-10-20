"""Exact independent strict-SPD complete-pivot and fresh-Schur checks."""

from fractions import Fraction
from itertools import combinations, product

from kernelgecp.surrogates import fraction_determinant


def _step(matrix, p):
    assert matrix[p][p] > 0
    return [
        [v - matrix[i][p] * matrix[p][j] / matrix[p][p] for j, v in enumerate(row)]
        for i, row in enumerate(matrix)
    ]


def _check_fresh_congruence(matrix, pivot, fresh):
    others = [i for i in fresh if i != pivot]
    indices = [pivot, *others]
    gram = [[matrix[i][j] for j in indices] for i in indices]
    elimination = [[-matrix[pivot][j] / matrix[pivot][pivot] for j in others]]
    elimination += [
        [Fraction(i == j) for j in range(len(others))] for i in range(len(others))
    ]
    # The lower identity block proves injectivity without solving or inversion.
    assert elimination[1:] == [
        [Fraction(i == j) for j in range(len(others))] for i in range(len(others))
    ]
    updated = _step(matrix, pivot)
    for i, j in product(range(len(others)), repeat=2):
        congruence = sum(
            (
                elimination[a][i] * gram[a][b] * elimination[b][j]
                for a, b in product(range(len(indices)), repeat=2)
            ),
            Fraction(),
        )
        assert congruence == updated[others[i]][others[j]]


def _all_complete_branches(matrix, selected=()):
    n = len(matrix)
    fresh = [i for i in range(n) if i not in selected]
    assert matrix == list(map(list, zip(*matrix, strict=True)))
    for i in selected:
        assert all(v == 0 for v in matrix[i])
    for order in range(1, len(fresh) + 1):
        for indices in combinations(fresh, order):
            assert (
                fraction_determinant([[matrix[i][j] for j in indices] for i in indices])
                > 0
            )
    for i, j in product(fresh, repeat=2):
        if i != j:
            assert abs(matrix[i][j]) ** 2 < matrix[i][i] * matrix[j][j]

    maximum = max(abs(v) for row in matrix for v in row)
    if maximum == 0:
        assert not fresh and len(selected) == n
        return [selected]
    complete = {
        (i, j) for i, j in product(range(n), repeat=2) if abs(matrix[i][j]) == maximum
    }
    cholesky = {(i, i) for i in range(n) if matrix[i][i] == maximum}
    assert complete == cholesky  # necessity, not diagonal tie preference
    traces = []
    for p, q in sorted(complete):
        assert p == q and p in fresh and matrix[p][p] > 0
        _check_fresh_congruence(matrix, p, fresh)
        traces.extend(_all_complete_branches(_step(matrix, p), (*selected, p)))
    return traces


def test_strict_feature_kernel_every_complete_branch_and_prefix():
    nodes = [Fraction(i, 2) for i in range(-2, 3)]
    matrix = [
        [sum((x**d * y**d / (d + 1) for d in range(5)), Fraction()) for y in nodes]
        for x in nodes
    ]
    # Injective Vandermonde features have full rank; check all exact principal
    # minors as an independent strict-positivity certificate for this sample.
    traces = _all_complete_branches(matrix)
    assert len(traces) >= 2  # reflection ties do not undermine diagonal necessity
    assert all(set(trace) == set(range(5)) for trace in traces)


def test_all_diagonal_ties_are_kept_not_forced_to_one_trace():
    matrix = [[Fraction(i == j) for j in range(3)] for i in range(3)]
    traces = _all_complete_branches(matrix)
    assert len(traces) == 6
    assert len(set(traces)) == 6


def test_strictness_cannot_be_dropped_from_pivot_necessity():
    matrix = [[Fraction(1), Fraction(-1)], [Fraction(-1), Fraction(1)]]
    assert fraction_determinant(matrix) == 0
    assert all(matrix[i][i] > 0 for i in range(2))
    assert all(abs(v) <= abs(matrix[0][1]) for row in matrix for v in row)
    assert matrix[0][1] != 0  # a successful off-diagonal complete pivot in weak PSD
