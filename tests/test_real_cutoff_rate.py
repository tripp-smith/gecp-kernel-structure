"""Exact count arithmetic and independent analytic/first-pivot observations."""

import itertools
from fractions import Fraction

import mpmath as mp


def _mp(x: Fraction):
    return mp.mpf(x.numerator) / x.denominator


def _cover(cutoff: Fraction) -> int:
    # No floating logarithm determines which side of a dyadic boundary wins.
    s = 0
    while 2**s < cutoff:
        s += 1
    return s


def _kernel(t, w):
    if w >= 0:
        return mp.exp(-t * w) / (1 + mp.exp(-w))
    return mp.exp((1 - t) * w) / (1 + mp.exp(w))


def test_every_rank_block_conversion_exactly() -> None:
    checked = 0
    for width in range(1, 65):
        block = 128 * width
        for j in range(25):
            for delta in (-1, 0, 1, block - 1):
                k = max(0, block * j + delta)
                quotient = k // block
                q = max(0, quotient - 9)
                if quotient >= 9:
                    assert 32 * width * (33 + 4 * q) <= k
                else:
                    assert q == 0  # The initial error bound is used here.
                # Equivalent to the exponent comparison, without logs/floats.
                assert Fraction(k, block) < q + 10
                checked += 1
    assert checked == 6400


def test_real_cutoff_logarithmic_covers_at_two_precisions() -> None:
    cutoffs = {Fraction(1, 32), Fraction(3, 4), Fraction(10**6)}
    for s in range(17):
        cutoffs.update(
            (
                Fraction(2**s) - Fraction(1, 4),
                Fraction(2**s),
                Fraction(2**s) + Fraction(1, 4),
            )
        )
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            for cutoff in sorted(cutoffs):
                s = _cover(cutoff)
                assert cutoff <= 2**s
                assert not s or 2 ** (s - 1) < cutoff
                lam = _mp(cutoff)
                assert (s + 1) * mp.log(2) <= 2 * mp.log(2 + lam)
                if cutoff >= 1:
                    assert (s + 1) * mp.log(2) <= 2 * mp.log(1 + lam)
                for k in (0, 1, 1056 * (s + 1), 4096 * (s + 1)):
                    cover = 1024 * mp.exp(-mp.log(2) * k / (128 * (s + 1)))
                    envelope = 1024 * mp.exp(
                        -(mp.log(2) ** 2) * k / (256 * mp.log(2 + lam))
                    )
                    assert 0 < cover <= envelope
                    values.append(mp.nstr(envelope, 45))
        records.append(values)
    assert len(cutoffs) == 53
    assert len(records[0]) == 212
    assert records[0] == records[1]


def test_real_tolerance_ceiling_including_large_tolerances() -> None:
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            for cutoff, epsilon in itertools.product(
                (Fraction(1, 32), Fraction(1), Fraction(3, 2), Fraction(10**6)),
                (
                    Fraction(2048),
                    Fraction(1024),
                    Fraction(2),
                    Fraction(1),
                    Fraction(1, 2),
                    Fraction(1, 10**6),
                    Fraction(1, 10**20),
                ),
            ):
                lam, eps = _mp(cutoff), _mp(epsilon)
                threshold = 256 / mp.log(2) ** 2 * mp.log(2 + lam) * mp.log(1024 / eps)
                rank = max(0, int(mp.ceil(threshold)))
                assert rank >= threshold
                bound = 1024 * mp.exp(
                    -(mp.log(2) ** 2) * rank / (256 * mp.log(2 + lam))
                )
                assert mp.isfinite(bound) and 0 < bound <= eps
                values.append((rank, mp.nstr(bound, 45)))
        records.append(values)
    assert len(records[0]) == 28
    assert records[0] == records[1]


def test_small_cutoff_obstruction_from_actual_first_complete_pivot() -> None:
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            for C, c in itertools.product(
                (Fraction(1, 16), Fraction(1), Fraction(1024)),
                (Fraction(1, 8), Fraction(1), Fraction(4)),
            ):
                cutoff = min(Fraction(1), c * c / (16 * C))
                assert 0 < cutoff <= 1 and 16 * C * cutoff <= c * c
                lam = _mp(cutoff)
                pivot = _kernel(mp.mpf(0), lam)
                for t, w in itertools.product(
                    (mp.mpf(0), mp.mpf("0.5"), mp.mpf(1)), (-lam, mp.mpf(0), lam)
                ):
                    assert 0 < _kernel(t, w) <= pivot
                actual = (
                    _kernel(mp.mpf(1), -lam)
                    - _kernel(mp.mpf(1), lam) * _kernel(mp.mpf(0), -lam) / pivot
                )
                exact_form = -mp.expm1(-lam)
                assert mp.isfinite(actual) and actual > 0, "Uncertified residual"
                assert mp.almosteq(
                    actual, exact_form, rel_eps=mp.mpf(10) ** (-digits + 25)
                )
                proposed = _mp(C) * mp.exp(-_mp(c) / mp.log1p(lam))
                assert mp.isfinite(proposed) and proposed > 0
                assert proposed <= lam / 4 < lam / 2 <= actual
                values.append((mp.nstr(actual, 45), mp.nstr(proposed, 45)))
        records.append(values)
    assert len(records[0]) == 9
    assert records[0] == records[1]
