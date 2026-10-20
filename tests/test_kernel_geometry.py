"""Exact, independent checks of kernel PSD order and original-kernel geometry."""

from fractions import Fraction
from itertools import combinations

from kernelgecp.surrogates import fraction_determinant


def _kernel(x, y):
    return sum((x**d * y**d / (d + 1) for d in range(5)), Fraction())


def _step(matrix, pivot):
    value = matrix[pivot][pivot]
    assert value > 0
    return [
        [
            entry - matrix[i][pivot] * matrix[pivot][j] / value
            for j, entry in enumerate(row)
        ]
        for i, row in enumerate(matrix)
    ]


def _assert_psd(matrix):
    # For a real symmetric finite matrix, every principal minor is nonnegative
    # iff it is PSD. This includes singular and repeated-sample cases.
    assert matrix == list(map(list, zip(*matrix, strict=True)))
    for order in range(1, len(matrix) + 1):
        for indices in combinations(range(len(matrix)), order):
            minor = [[matrix[i][j] for j in indices] for i in indices]
            assert fraction_determinant(minor) >= 0


def test_every_diagonal_prefix_original_energy_and_fill_bound():
    nodes = [Fraction(i, 16) for i in range(17)]
    original = [[_kernel(x, y) for y in nodes] for x in nodes]
    # d*x^(d-1)*y^d/(d+1) <= d/(d+1) on the entire interval [0,1].
    lipschitz = sum((Fraction(d, d + 1) for d in range(1, 5)), Fraction())
    for i, x in enumerate(nodes):
        for j, s in enumerate(nodes):
            for y in range(len(nodes)):
                assert abs(original[i][y] - original[j][y]) <= lipschitz * abs(x - s)

    records = 0
    for order in ((0, 16, 8, 4, 12), (16, 15, 14, 13, 12)):
        # Neither sequence is required to maximize the diagonal.
        residual, selected = original, []
        for length in range(6):
            removed = [
                [v - residual[i][j] for j, v in enumerate(row)]
                for i, row in enumerate(original)
            ]
            sample = (0, 4, 8, 16, 4)  # repeated nodes must remain valid
            for matrix in (residual, removed):
                _assert_psd([[matrix[i][j] for j in sample] for i in sample])
            diagonal_max = max(residual[i][i] for i in range(len(nodes)))
            assert max(map(abs, (v for row in residual for v in row))) == diagonal_max
            for i, x in enumerate(nodes):
                for s in selected:
                    assert all(v == 0 for v in residual[s])
                    energy = original[i][i] - 2 * original[i][s] + original[s][s]
                    assert residual[i][i] <= energy <= 2 * lipschitz * abs(x - nodes[s])
                    records += 1
            if selected:
                fill = max(min(abs(x - nodes[s]) for s in selected) for x in nodes)
                assert all(
                    abs(v) <= 2 * lipschitz * fill for row in residual for v in row
                )
            if length < 5:
                pivot = order[length]
                residual = _step(residual, pivot)
                selected.append(pivot)
        assert all(v == 0 for row in residual for v in row)
    assert records == 510


def test_weak_psd_off_diagonal_complete_tie_is_not_a_diagonal_pivot():
    matrix = [[Fraction(1) for _ in range(2)] for _ in range(2)]
    _assert_psd(matrix)
    p, q = 0, 1
    assert p != q and matrix[p][q] != 0
    assert all(abs(v) <= abs(matrix[p][q]) for row in matrix for v in row)
    # Strict SPD matters for necessity; weak PSD only guarantees availability.
    assert fraction_determinant(matrix) == 0


def test_original_lipschitz_constant_cannot_replace_twice_constant():
    # The PSD kernel 1-|x-y| on [0,1], sampled at 0 and 1/4, has L=1.
    # After the pivot at 0, R(1/4,1/4)=7/16 exceeds L*h=1/4.
    x = Fraction(1, 4)
    matrix = [[Fraction(1), 1 - x], [1 - x, Fraction(1)]]
    _assert_psd(matrix)
    residual = _step(matrix, 0)
    assert residual[1][1] == Fraction(7, 16)
    assert x < residual[1][1] <= 2 * x
