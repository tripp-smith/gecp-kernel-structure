from fractions import Fraction

import mpmath as mp


def test_real_d1_parameter_choices_and_rank_constant() -> None:
    records = []
    for digits in (160, 320):
        values = []
        with mp.workdps(digits):
            for cutoff in map(Fraction, (1, 2, 3, 16, 10**6)):
                s = 0
                while 2**s < cutoff:
                    s += 1
                for epsilon in (
                    Fraction(1, 2),
                    Fraction(1, 3),
                    Fraction(1, 10**6),
                    Fraction(1, 2**40),
                    1 - Fraction(1, 2**20),
                ):
                    p = 1
                    while Fraction(1, 2**p) > epsilon:
                        p += 1
                    assert p > 0 and Fraction(1, 2**p) <= epsilon
                    rank = 16 * p * (s + 1)
                    eps = mp.mpf(epsilon.numerator) / epsilon.denominator
                    bound = (
                        128
                        * mp.log(1 + cutoff.numerator / cutoff.denominator)
                        * mp.log(128 / eps)
                    )
                    assert rank <= bound
                    values.append((rank, mp.nstr(bound, 45)))
        records.append(values)
    assert len(records[0]) == 25
    assert records[0] == records[1]
