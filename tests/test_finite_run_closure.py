"""Finite-only exact witnesses, solve-based crosses, and exact-stop checks."""

import itertools
from fractions import Fraction

import mpmath as mp

from kernelgecp.surrogates import fraction_determinant


def _sample(sign=1):
    return [
        [
            Fraction(2) ** (sign * i * j) / (1 + Fraction(2) ** (sign * 5 * j))
            for j in range(-3, 3)
        ]
        for i in range(6)
    ]


def _step(matrix, p, q):
    pivot = matrix[p][q]
    assert pivot  # A zero complete maximum is a stop, never an update.
    return [
        [v - matrix[i][q] * matrix[p][j] / pivot for j, v in enumerate(line)]
        for i, line in enumerate(matrix)
    ]


def _trace(matrix):
    states, pivots = [matrix], []
    for _ in range(min(len(matrix), len(matrix[0]))):
        p, q = max(
            itertools.product(range(len(matrix)), range(len(matrix[0]))),
            key=lambda ij: (abs(matrix[ij[0]][ij[1]]), -ij[0], -ij[1]),
        )
        if matrix[p][q] == 0:
            break
        pivots.append((p, q, matrix[p][q]))
        matrix = _step(matrix, p, q)
        states.append(matrix)
    return states, pivots


def _minor(matrix, rows, columns):
    return fraction_determinant([[matrix[i][j] for j in columns] for i in rows])


def _solve(matrix, rhs):
    """Independent exact elimination/back substitution; no inverse is formed."""
    work = [list(line) + [value] for line, value in zip(matrix, rhs, strict=True)]
    n = len(work)
    for k in range(n):
        p = next(i for i in range(k, n) if work[i][k])
        work[k], work[p] = work[p], work[k]
        for i in range(k + 1, n):
            ratio = work[i][k] / work[k][k]
            for j in range(k, n + 1):
                work[i][j] -= ratio * work[k][j]
    result = [Fraction()] * n
    for k in reversed(range(n)):
        result[k] = (
            work[k][n]
            - sum((work[k][j] * result[j] for j in range(k + 1, n)), Fraction())
        ) / work[k][k]
    return result


def _cross(matrix, pivots):
    rows = [p for p, _, _ in pivots]
    columns = [q for _, q, _ in pivots]
    core = [[matrix[p][q] for q in columns] for p in rows]
    result = [[Fraction() for _ in matrix[0]] for _ in matrix]
    for y in range(len(matrix[0])):
        coefficients = _solve(core, [matrix[p][y] for p in rows])
        for x in range(len(matrix)):
            result[x][y] = sum(
                (matrix[x][q] * c for q, c in zip(columns, coefficients, strict=True)),
                Fraction(),
            )
    return result


def test_every_prefix_and_mark_set_without_future_pivots() -> None:
    matrix = _sample()
    states, pivots = _trace(matrix)
    witnesses = steps = 0
    for length in range(len(pivots) + 1):
        # All subsequent calculations receive only these finite slices.
        prefix, prefix_states = pivots[:length], states[: length + 1]
        for mask in range(2**length):
            rows = columns = ()
            product = Fraction(1)
            for k in reversed(range(length)):
                p, q, pivot = prefix[k]
                before = abs(_minor(prefix_states[k + 1], rows, columns))
                assert before > 0
                if mask & (1 << k):
                    rows, columns = (p,) + rows, (q,) + columns
                    product *= abs(pivot)
                    assert (
                        abs(_minor(prefix_states[k], rows, columns))
                        == abs(pivot) * before
                    )
                else:
                    ar, ac = sorted((p,) + rows), sorted((q,) + columns)
                    rank = ac.index(q)
                    rows = tuple(v for i, v in enumerate(ar) if i != rank)
                    columns = tuple(v for i, v in enumerate(ac) if i != rank)
                    assert abs(_minor(prefix_states[k], rows, columns)) >= before
                assert abs(_minor(prefix_states[k], rows, columns)) >= product
                steps += 1
            assert len(rows) == len(columns) == mask.bit_count()
            assert set(columns) == {
                prefix[i][1] for i in range(length) if mask & (1 << i)
            }
            assert abs(_minor(matrix, rows, columns)) >= product
            witnesses += 1
    assert witnesses == 127 and steps == 642


def test_exact_selected_cross_residual_and_green_sums_at_every_prefix() -> None:
    rho = [Fraction((-1) ** j * (j + 1), 7) for j in range(6)]
    checked = 0
    for sign in (-1, 1):
        matrix = _sample(sign)
        states, pivots = _trace(matrix)
        for length in range(len(pivots) + 1):
            cross = _cross(matrix, pivots[:length])
            epsilon = max(abs(v) for row in states[length] for v in row)
            for x in range(6):
                actual = sum((rho[y] * matrix[x][y] for y in range(6)), Fraction())
                approximation = sum(
                    (rho[y] * cross[x][y] for y in range(6)), Fraction()
                )
                integrated_error = sum(
                    (rho[y] * states[length][x][y] for y in range(6)), Fraction()
                )
                assert actual - approximation == integrated_error
                assert abs(integrated_error) <= epsilon * sum(map(abs, rho))
                for y in range(6):
                    assert matrix[x][y] - cross[x][y] == states[length][x][y]
                    checked += 1
    assert checked == 504


def test_rank_deficient_exact_stop_does_not_attempt_zero_update() -> None:
    for matrix, expected_rank in (
        ([[Fraction(0) for _ in range(5)] for _ in range(4)], 0),
        ([[Fraction((i + 1) * (j + 2), 7) for j in range(5)] for i in range(4)], 1),
    ):
        states, pivots = _trace(matrix)
        assert len(pivots) == expected_rank
        assert len(states) == expected_rank + 1
        assert all(pivot for _, _, pivot in pivots)
        assert all(v == 0 for row in states[-1] for v in row)
        assert _cross(matrix, pivots) == matrix


def test_actual_exponential_finite_crosses_at_two_precisions() -> None:
    exact = _sample(-1)
    _, pivots = _trace(exact)
    records = []
    green_checks = 0
    rho = [Fraction((-1) ** j * (j + 1), 7) for j in range(6)]
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            matrix = [
                [
                    mp.exp(-mp.mpf(i) / 5 * (5 * j * mp.log(2)))
                    / (1 + mp.exp(-5 * j * mp.log(2)))
                    for j in range(-3, 3)
                ]
                for i in range(6)
            ]
            weights = [mp.mpf(v.numerator) / v.denominator for v in rho]
            tolerance = mp.mpf(10) ** (-digits + 25)
            for length in range(7):
                prefix = pivots[:length]
                residual = [line.copy() for line in matrix]
                for p, q, _ in prefix:
                    residual = _step(residual, p, q)
                rows, columns = [p for p, _, _ in prefix], [q for _, q, _ in prefix]
                reference = _cross(exact, prefix)
                cross_matrix = [[mp.mpf(0) for _ in range(6)] for _ in range(6)]
                if length:
                    core = mp.matrix([[matrix[p][q] for q in columns] for p in rows])
                    determinant = abs(mp.det(core))
                    assert mp.isfinite(determinant) and determinant > 0, (
                        "Uncertified core"
                    )
                    values.append(mp.nstr(determinant, 45))
                for y in range(6):
                    solution = (
                        mp.lu_solve(core, mp.matrix([matrix[p][y] for p in rows]))
                        if length
                        else []
                    )
                    for x in range(6):
                        cross = mp.fsum(
                            matrix[x][q] * solution[j] for j, q in enumerate(columns)
                        )
                        target = (
                            mp.mpf(reference[x][y].numerator)
                            / reference[x][y].denominator
                        )
                        assert mp.isfinite(cross)
                        assert mp.almosteq(
                            cross, target, rel_eps=tolerance, abs_eps=tolerance
                        )
                        assert mp.almosteq(
                            matrix[x][y] - cross,
                            residual[x][y],
                            rel_eps=tolerance,
                            abs_eps=tolerance,
                        )
                        values.append(mp.nstr(cross, 45))
                        cross_matrix[x][y] = cross
                for x in range(6):
                    original_green = mp.fsum(
                        weights[y] * matrix[x][y] for y in range(6)
                    )
                    approximate_green = mp.fsum(
                        weights[y] * cross_matrix[x][y] for y in range(6)
                    )
                    residual_green = mp.fsum(
                        weights[y] * residual[x][y] for y in range(6)
                    )
                    exact_error = sum(
                        (rho[y] * (exact[x][y] - reference[x][y]) for y in range(6)),
                        Fraction(),
                    )
                    target_error = (
                        mp.mpf(exact_error.numerator) / exact_error.denominator
                    )
                    observed_error = original_green - approximate_green
                    assert mp.isfinite(observed_error)
                    assert mp.almosteq(
                        observed_error,
                        target_error,
                        rel_eps=tolerance,
                        abs_eps=tolerance,
                    )
                    assert mp.almosteq(
                        observed_error,
                        residual_green,
                        rel_eps=tolerance,
                        abs_eps=tolerance,
                    )
                    green_checks += 1
        records.append(values)
    assert len(records[0]) == 258
    assert records[0] == records[1]
    assert green_checks == 84
