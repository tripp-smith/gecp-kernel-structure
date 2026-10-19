"""Exact panel structure/arithmetic and separate actual-exponential evidence."""

import itertools
import math
from fractions import Fraction

import mpmath as mp


def _contains(negative: bool, panel: int, omega: Fraction) -> bool:
    if negative:
        if omega >= 0:
            return False
        omega = -omega
    if panel == 0:
        return 0 <= omega <= 1
    return 2 ** (panel - 1) < omega <= 2**panel


def _coordinates(p: int, panel: int, negative: bool):
    lower = 2 ** max(0, panel - 1)
    cutoff = min(Fraction(1), Fraction(p, lower))
    times = (
        Fraction(0),
        Fraction(1),
        cutoff,
        max(Fraction(), cutoff - Fraction(1, 4 * lower)),
        min(Fraction(1), cutoff + Fraction(1, 4 * lower)),
    )
    frequencies = (
        (Fraction(0), Fraction(1, 2), Fraction(1))
        if panel == 0
        else (
            Fraction(lower) + Fraction(1, 4),
            Fraction(3 * lower, 2),
            Fraction(2 * lower),
        )
    )
    if negative:
        times = tuple(1 - t for t in times)
        frequencies = tuple(-w for w in frequencies if w != 0)
    return times, frequencies


def _mp(value: Fraction):
    return mp.mpf(value.numerator) / value.denominator


def _kernel(t, omega):
    # Both branches avoid exponentially large intermediate values.
    if omega >= 0:
        return mp.exp(-t * omega) / (1 + mp.exp(-omega))
    return mp.exp((1 - t) * omega) / (1 + mp.exp(omega))


def _approximation(t: Fraction, omega: Fraction, p: int, panel: int, negative: bool):
    if negative:
        t, omega = 1 - t, -omega
    if panel and t * 2 ** (panel - 1) > p:
        return mp.mpf(0)
    x = _mp(t) * _mp(omega)
    return mp.fsum((-x) ** k / mp.factorial(k) for k in range(8 * p)) / (
        1 + mp.exp(-_mp(omega))
    )


def test_masked_panel_feature_factorization_exactly() -> None:
    checked = 0
    for p, panel, negative in itertools.product(
        (1, 2, 4), (0, 1, 5, 12), (False, True)
    ):
        times, frequencies = _coordinates(p, panel, negative)
        for t, omega in itertools.product(times, frequencies):
            assert _contains(negative, panel, omega)
            u, w = (1 - t, -omega) if negative else (t, omega)
            active = not panel or u * 2 ** (panel - 1) <= p
            # A rational column weight isolates exact finite-factor algebra.
            # The actual exponential denominator is checked separately below.
            weight = 1 / (1 + w * w)
            left = [
                (-u) ** k / math.factorial(k) if active else Fraction()
                for k in range(8 * p)
            ]
            right = [w**k * weight for k in range(8 * p)]
            value = sum((a * b for a, b in zip(left, right, strict=True)), Fraction())
            direct = (
                sum(
                    ((-u * w) ** k / math.factorial(k) for k in range(8 * p)),
                    Fraction(),
                )
                * weight
                if active
                else Fraction()
            )
            assert len(left) == len(right) == 8 * p
            assert value == direct
            checked += 1
    assert checked == 345


def test_actual_panel_error_bound_at_independent_precisions() -> None:
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            for p, panel, negative in itertools.product(
                (1, 2, 4), (0, 1, 4, 12, 20), (False, True)
            ):
                times, frequencies = _coordinates(p, panel, negative)
                for t, omega in itertools.product(times, frequencies):
                    assert _contains(negative, panel, omega)
                    error = abs(
                        _kernel(_mp(t), _mp(omega))
                        - _approximation(t, omega, p, panel, negative)
                    )
                    assert mp.isfinite(error)
                    assert error <= mp.mpf(2) ** (-p)
                    values.append(mp.nstr(error, 45))
        records.append(values)
    assert len(records[0]) == 435
    assert records[0] == records[1]


def test_small_order_actual_panel_determinants_at_independent_precisions() -> None:
    # These small samples check implementation behavior, not the enormous
    # target-order determinant. The all-orders budget is proved in Lean.
    records = []
    for digits in (320, 640):
        values = []
        with mp.workdps(digits):
            for p, panel, negative in itertools.product(
                (1, 2), (0, 1, 8), (False, True)
            ):
                order = 16 * p
                lower = 2 ** max(0, panel - 1)
                times = [Fraction(i, (order - 1) * lower) for i in range(order)]
                frequencies = (
                    [Fraction(j, order - 1) for j in range(order)]
                    if panel == 0
                    else [
                        lower + Fraction(lower * (j + 1), order) for j in range(order)
                    ]
                )
                if negative:
                    times = [1 - t for t in times]
                    frequencies = [-w for w in frequencies]
                    if panel == 0:
                        # Zero is allocated to the positive central panel only.
                        frequencies = [-Fraction(j + 1, order) for j in range(order)]
                assert all(_contains(negative, panel, w) for w in frequencies)
                matrix = mp.matrix(
                    [[_kernel(_mp(t), _mp(w)) for w in frequencies] for t in times]
                )
                determinant = abs(mp.det(matrix))
                assert mp.isfinite(determinant) and determinant > 0, (
                    "Uncertified determinant"
                )
                bound = (
                    2**order
                    * math.factorial(order)
                    * mp.mpf(2) ** (-p * (order - 8 * p) + 8 * p)
                )
                assert determinant <= bound
                values.append(mp.nstr(determinant, 45))
        records.append(values)
    assert len(records[0]) == 12
    assert records[0] == records[1]


def test_local_budget_exponent_and_accuracy_arithmetic_exactly() -> None:
    for m in range(65):
        p = 2 * m + 1
        order, rank = 16 * p, 8 * p
        assert p * (order - rank) == order * m + rank
        assert math.factorial(order) <= order**order
    for q in range(1001):
        assert 32 * (33 + 4 * q) * 2**q <= 2 ** (16 + 2 * q)
    for q in range(17):
        p = 33 + 4 * q
        order, rank = 16 * p, 8 * p
        # Clear positive dyadic denominators; no floating-point arithmetic.
        assert math.factorial(order) * 2 ** (order + rank + q * order) <= 2 ** (
            p * (order - rank)
        )


def test_signed_panel_membership_boundaries_do_not_overlap() -> None:
    # Finite exact regression only; the universal partition/count proof is
    # intentionally left to the next research phase.
    for scale in range(9):
        positive = {Fraction(0), Fraction(1, 2), Fraction(1)}
        for j in range(scale):
            positive.update(
                (Fraction(2**j), Fraction(3 * 2**j, 2), Fraction(2 ** (j + 1)))
            )
        samples = positive | {-w for w in positive}
        for omega in samples:
            hits = sum(
                _contains(negative, panel, omega)
                for negative in (False, True)
                for panel in range(scale + 1)
            )
            assert hits == 1
