"""Identity objectives expose binary64 cover and enclosure rounding errors."""

import math
from fractions import Fraction

import mpmath as mp
import pytest

from kernelgecp import (
    NumericalBreakdownError,
    certified_pivot,
    interval_certified_pivot,
)


def test_lipschitz_upper_encloses_one_ulp_interval_identity() -> None:
    maximum = math.nextafter(1.0, math.inf)
    result = certified_pivot(
        lambda t, omega: t,
        t_bounds=(1.0, maximum),
        omega_bounds=(0.0, 1.0),
        lipschitz=(1.0, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.upper_bound >= maximum
    assert result.lower_bound <= result.t <= maximum
    assert not result.certified
    assert result.termination_reason == "cell_budget"


@pytest.mark.parametrize("name", ["abs_tol", "rel_tol"])
@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_lipschitz_rejects_nonfinite_tolerances(name: str, value: float) -> None:
    with pytest.raises(ValueError, match="tolerance"):
        certified_pivot(
            lambda t, omega: t,
            omega_bounds=(0.0, 1.0),
            lipschitz=1.0,
            **{name: value},
        )


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("sign", [-1.0, 1.0])
@pytest.mark.parametrize("negative", [False, True])
def test_identity_enclosures_cover_both_signs_and_axes(
    axis: int, sign: float, negative: bool
) -> None:
    endpoint = math.nextafter(1.0, math.inf)
    bounds = (-endpoint, -1.0) if negative else (1.0, endpoint)
    result = certified_pivot(
        lambda t, omega: sign * (t if axis == 0 else omega),
        t_bounds=bounds if axis == 0 else (0.0, 1.0),
        omega_bounds=bounds if axis == 1 else (0.0, 1.0),
        lipschitz=(1.0, 0.0) if axis == 0 else (0.0, 1.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.upper_bound >= endpoint
    assert result.lower_bound == abs(result.value) <= endpoint
    assert not result.certified


@pytest.mark.parametrize("negative", [False, True])
def test_subnormal_radius_is_not_rounded_to_zero(negative: bool) -> None:
    endpoint = math.nextafter(0.0, math.inf)
    result = certified_pivot(
        lambda t, omega: t,
        t_bounds=(-endpoint, 0.0) if negative else (0.0, endpoint),
        omega_bounds=(0.0, 1.0),
        lipschitz=(1.0, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.upper_bound >= endpoint
    assert result.lower_bound == 0.0
    assert result.eta == 0.0 and not result.certified


def test_large_finite_width_is_bounded_without_intermediate_overflow() -> None:
    endpoint = float.fromhex("0x1.fffffffffffffp+1023")
    result = certified_pivot(
        lambda t, omega: t,
        t_bounds=(-endpoint, endpoint),
        omega_bounds=(0.0, 1.0),
        lipschitz=(1.0, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.upper_bound == endpoint
    assert result.lower_bound == 0.0 and not result.certified


def test_positive_lipschitz_product_underflow_is_enclosed_outward() -> None:
    constant = math.nextafter(0.0, math.inf)
    result = certified_pivot(
        lambda t, omega: 0.0,
        omega_bounds=(0.0, 1.0),
        lipschitz=(constant, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert Fraction(result.upper_bound) >= Fraction(constant) / 2 > 0
    assert not result.certified


def test_eta_is_a_downward_bound_on_the_reported_gap_ratio() -> None:
    result = certified_pivot(
        lambda t, omega: 1.0,
        omega_bounds=(0.0, 1.0),
        lipschitz=(0.5, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.upper_bound == 1.25
    assert Fraction(result.eta) <= Fraction(result.lower_bound) / Fraction(
        result.upper_bound
    )


def test_gap_comparison_does_not_round_the_allowed_tolerance_up() -> None:
    absolute = math.nextafter(0.25, 0.0)
    relative = 0.75 * (0.25 - absolute)
    assert absolute + relative == 0.25  # Rounded sum would falsely pass.
    result = certified_pivot(
        lambda t, omega: 1.0,
        omega_bounds=(0.0, 1.0),
        lipschitz=(0.5, 0.0),
        abs_tol=absolute,
        rel_tol=relative,
        max_cells=1,
    )
    assert Fraction(result.upper_bound) - Fraction(result.lower_bound) > (
        Fraction(absolute) + Fraction(relative) * Fraction(result.lower_bound)
    )
    assert not result.certified


def test_sampled_ties_are_recorded_without_a_hidden_callback_reevaluation() -> None:
    calls: list[tuple[float, float]] = []

    def constant(t: float, omega: float) -> float:
        calls.append((t, omega))
        return -1.0

    result = certified_pivot(
        constant,
        omega_bounds=(0.0, 1.0),
        lipschitz=(1.0, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=5,
    )
    assert len(calls) == result.cells_evaluated == 5
    assert (result.t, result.omega) == min(calls)
    assert result.candidate_count == result.exact_tie_count == len(set(calls))
    assert result.near_tie_count == result.exact_tie_count
    assert result.value == -1.0 and result.lower_bound == 1.0
    assert result.tie_rtol == 0.0
    assert result.tie_scope == "distinct_sampled_exact_callback_magnitudes"
    assert (
        result.objective_scope == "exact_point_values_with_valid_cell_lipschitz_bounds"
    )


def test_resolution_exhaustion_preserves_the_unresolved_upper_bound() -> None:
    endpoint = math.nextafter(1.0, math.inf)
    result = certified_pivot(
        lambda t, omega: t,
        t_bounds=(1.0, endpoint),
        omega_bounds=(0.0, math.nextafter(0.0, math.inf)),
        lipschitz=(1.0, 0.0),
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=100,
    )
    assert not result.certified and result.termination_reason == "coordinate_resolution"
    assert result.cells_evaluated == 1
    assert result.upper_bound >= endpoint > result.lower_bound


def test_zero_lipschitz_constant_can_still_certify_an_exact_constant() -> None:
    endpoint = float.fromhex("0x1.fffffffffffffp+1023")
    result = certified_pivot(
        lambda t, omega: 3.0,
        t_bounds=(-endpoint, endpoint),
        omega_bounds=(-endpoint, endpoint),
        lipschitz=0.0,
        abs_tol=0.0,
        rel_tol=0.0,
        max_cells=1,
    )
    assert result.certified and result.lower_bound == result.upper_bound == 3.0
    assert result.eta == 1.0 and result.candidate_count == 1


@pytest.mark.parametrize("budget", [0, -1, True, 1.5, math.inf, math.nan])
def test_cell_budget_must_be_a_positive_integer(budget: object) -> None:
    with pytest.raises(ValueError, match="budget"):
        certified_pivot(
            lambda t, omega: t,
            omega_bounds=(0.0, 1.0),
            lipschitz=1.0,
            max_cells=budget,  # type: ignore[arg-type]
        )


def test_unrepresentable_upper_bound_is_an_explicit_arithmetic_failure() -> None:
    with pytest.raises(NumericalBreakdownError, match="upper bound.*range"):
        certified_pivot(
            lambda t, omega: 1e308,
            t_bounds=(0.0, 1e308),
            omega_bounds=(0.0, 1.0),
            lipschitz=(1e308, 0.0),
        )


def test_standalone_lipschitz_metadata_exposes_all_evidence_categories() -> None:
    result = certified_pivot(
        lambda t, omega: 2.0,
        omega_bounds=(0.0, 1.0),
        lipschitz=0.0,
    )
    metadata = result.metadata
    assert metadata.precision_bits == metadata.requested_precision_bits == 53
    assert metadata.arithmetic == "binary64_points_exact_rational_lipschitz_bounds"
    assert metadata.configuration["enclosure_arithmetic"] == "exact_binary_rational"
    assert metadata.configuration["enclosure_precision_bits"] is None
    assert metadata.configuration["bound_precision_bits"] == 53
    assert metadata.configuration["value_precision_bits"] == 53
    assert metadata.configuration["coordinate_precision_bits"] == 53
    assert metadata.conditioning_status == "not_applicable"
    assert metadata.converged is result.certified is True
    assert metadata.stop_reason == result.termination_reason
    assert metadata.tie_counts == [1]
    assert metadata.certification_scope == result.objective_scope
    metadata.canonical_json()


@pytest.mark.parametrize("original_precision", [79, 80, 81, 82])
@pytest.mark.parametrize("fail", [False, True])
def test_interval_metadata_and_exact_context_restoration(
    original_precision: int, fail: bool
) -> None:
    previous = mp.iv.prec
    mp.iv.prec = original_precision

    def objective(t: object, omega: object) -> object:
        if fail:
            raise ValueError("deliberate callback failure")
        return t

    try:
        if fail:
            with pytest.raises(ValueError, match="callback failure"):
                interval_certified_pivot(
                    lambda t, omega: t,
                    objective,
                    omega_bounds=(0.0, 1.0),
                    precision_bits=96,
                )
        else:
            result = interval_certified_pivot(
                lambda t, omega: t,
                objective,
                omega_bounds=(0.0, 1.0),
                precision_bits=96,
                abs_tol=0.01,
                max_cells=1,
            )
            metadata = result.metadata
            assert metadata.arithmetic == "mpmath_interval_with_float64_exports"
            assert metadata.precision_bits == metadata.requested_precision_bits == 96
            assert result.enclosure_precision_bits > 96
            assert metadata.configuration["enclosure_precision_bits"] == (
                result.enclosure_precision_bits
            )
            assert result.bound_precision_bits == result.coordinate_precision_bits == 53
            assert metadata.converged is result.certified
            assert metadata.tie_counts == [result.near_tie_count]
            assert metadata.conditioning_status == "not_applicable"
            metadata.canonical_json()
        assert mp.iv.prec == original_precision
    finally:
        mp.iv.prec = previous
