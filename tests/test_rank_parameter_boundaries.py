"""Exact checks of binary64 parameters at dyadic and underflow boundaries."""

import math
from fractions import Fraction

from kernelgecp import (
    DyadicTaylorApproximation,
    dlr_rank_bound,
    dyadic_taylor_rank_bound,
    exp_family_rank_bound,
    fermionic_dyadic_taylor_approximation,
)


def test_dyadic_factory_parameters_cover_actual_cutoff_and_tolerance() -> None:
    cutoff = math.nextafter(16.0, math.inf)
    tolerance = math.nextafter(2.0**-40, 0.0)
    result = fermionic_dyadic_taylor_approximation(cutoff, tolerance)
    assert Fraction(cutoff) <= 2**result.scale
    assert Fraction(1, 2**result.accuracy_bits) <= Fraction(tolerance)
    assert result.error_bound <= tolerance


def test_all_representable_sides_of_dyadic_transitions_have_minimal_parameters() -> (
    None
):
    cutoffs = [
        value
        for exponent in (0, 1, 4, 8, 32, 100, 500, 1023)
        for value in (
            math.nextafter(2.0**exponent, 0.0),
            2.0**exponent,
            math.nextafter(2.0**exponent, math.inf),
        )
        if value >= 1
    ]
    tolerances = [
        value
        for exponent in (1, 2, 20, 40, 100, 500, 1022, 1074)
        for value in (
            math.nextafter(2.0**-exponent, 0.0),
            2.0**-exponent,
            math.nextafter(2.0**-exponent, math.inf),
        )
        if 0 < value < 1
    ]
    for cutoff in cutoffs:
        for tolerance in tolerances:
            result = fermionic_dyadic_taylor_approximation(cutoff, tolerance)
            assert Fraction(cutoff) <= 2**result.scale
            if result.scale:
                assert 2 ** (result.scale - 1) < Fraction(cutoff)
            error = Fraction(1, 2**result.accuracy_bits)
            assert error <= Fraction(tolerance) < 2 * error
            assert error <= Fraction(result.error_bound) <= Fraction(tolerance)
            assert result.rank == dyadic_taylor_rank_bound(cutoff, tolerance)


def test_rank_helpers_handle_subnormal_tolerances_with_exact_counts() -> None:
    for cutoff in (1.0, math.nextafter(16.0, 0.0), 16.0, 1e308):
        for tolerance in (0.9, math.nextafter(2.0**-40, 0.0), 2.0**-1074):
            scales = 1
            while 2**scales <= Fraction(cutoff):
                scales += 1
            orders = 1
            while 4**orders <= 1 / Fraction(tolerance):
                orders += 1
            assert exp_family_rank_bound(cutoff, tolerance) == scales * orders
            ratio = 3 / Fraction(tolerance)
            band_orders = 1
            while 4**band_orders <= ratio:
                band_orders += 1
            central = dlr_rank_bound(cutoff, tolerance) - 2 * scales * band_orders
            # This exact binary cover is stronger than the natural-log count.
            assert 2 ** (central - 1) >= ratio
            assert 2 ** (central - 2) < ratio


def test_manual_subfloat_error_bound_is_exported_outward_not_as_zero() -> None:
    result = DyadicTaylorApproximation(cutoff=1.0, accuracy_bits=2000, scale=0)
    assert result.error_bound == math.nextafter(0.0, math.inf)
    assert Fraction(result.error_bound) >= Fraction(1, 2**2000)
    assert result.metadata.configuration["error_bound"] > 0
    result.metadata.canonical_json()
