"""Independent regressions for numerical defects exposed by the final audit."""

import json
from dataclasses import asdict

import mpmath as mp
import numpy as np
import pytest

from kernelgecp import (
    FermionicKernel,
    GECPConfig,
    NumericalBreakdownError,
    gecp,
    gecp_matrix,
    interval_certified_pivot,
)
from kernelgecp.certified import _make_interval_cell
from kernelgecp.gecp import _FermionicResidualEvaluator
from kernelgecp.research import certify_fermionic_block_contraction


def test_cancellation_point_estimate_never_becomes_a_certified_lower_bound() -> None:
    def point(t: object, omega: object) -> mp.mpf:
        return (mp.mpf(1) + mp.mpf("1.3e-16")) - mp.mpf(1)

    def interval(t: object, omega: object) -> object:
        return (mp.iv.mpf(1) + mp.iv.mpf("1.3e-16")) - mp.iv.mpf(1)

    certificate = interval_certified_pivot(
        point,
        interval,
        omega_bounds=(-1.0, 1.0),
        precision_bits=53,
        abs_tol=1e-20,
        rel_tol=0,
        max_cells=2,
    )
    with mp.workdps(80):
        truth = mp.mpf("1.3e-16")
        assert mp.mpf(certificate.lower_bound) <= truth <= certificate.upper_bound
        assert certificate.value > certificate.upper_bound
        assert mp.mpf(certificate.eta) <= truth / certificate.upper_bound
    assert certificate.certified
    json.dumps(asdict(certificate), allow_nan=False)


@pytest.mark.parametrize("value", ["1e-1000", "-1e-1000"])
def test_interval_underflow_cannot_certify_a_zero_objective(value: str) -> None:
    certificate = interval_certified_pivot(
        lambda t, omega: mp.mpf(value),
        lambda t, omega: mp.iv.mpf(value),
        omega_bounds=(1.0, 2.0),
        abs_tol=0,
        rel_tol=0,
        max_cells=1,
    )
    assert certificate.lower_bound == 0
    assert certificate.upper_bound == np.nextafter(0.0, np.inf)
    assert not certificate.certified
    assert certificate.eta == 0


def test_zero_crossing_singleton_enclosure_has_zero_absolute_lower_bound() -> None:
    certificate = interval_certified_pivot(
        lambda t, omega: mp.mpf(2),
        lambda t, omega: mp.iv.mpf([-1, 1]),
        omega_bounds=(1.0, 2.0),
        abs_tol=0,
        rel_tol=0,
        max_cells=1,
    )
    assert certificate.lower_bound == 0
    assert not certificate.certified


def test_mean_value_bound_uses_actual_rounded_midpoint_radius() -> None:
    left, right = 1.0, float(np.nextafter(1.0, np.inf))

    def objective(t: object, omega: object) -> object:
        return t - 1

    old_dps = mp.iv.dps
    mp.iv.dps = 60
    try:
        cell = _make_interval_cell(
            objective,
            objective,
            (left, right, left, right),
            0,
            128,
            lambda t, omega: (mp.iv.mpf(1), mp.iv.mpf(0)),
        )
    finally:
        mp.iv.dps = old_dps
    assert cell.upper >= right - left
    certificate = interval_certified_pivot(
        objective,
        objective,
        t_bounds=(left, right),
        omega_bounds=(left, right),
        abs_tol=0,
        rel_tol=0,
        max_cells=20,
        gradient_interval_function=lambda t, omega: (mp.iv.mpf(1), mp.iv.mpf(0)),
    )
    assert certificate.lower_bound <= right - left <= certificate.upper_bound
    assert not certificate.certified
    assert certificate.termination_reason == "coordinate_resolution"


def test_inconsistent_gradient_enclosures_fail_closed() -> None:
    with pytest.raises(ValueError, match="inconsistent.*bounds"):
        interval_certified_pivot(
            lambda t, omega: t,
            lambda t, omega: t,
            omega_bounds=(1.0, 2.0),
            gradient_interval_function=lambda t, omega: (mp.iv.mpf(0), mp.iv.mpf(0)),
            max_cells=1,
        )


def test_sampled_ties_are_recorded_without_continuous_multiplicity_claim() -> None:
    certificate = interval_certified_pivot(
        lambda t, omega: 1 + mp.mpf("1e-10") * t,
        lambda t, omega: 1 + mp.iv.mpf("1e-10") * t,
        omega_bounds=(-1.0, 1.0),
        abs_tol=1e-8,
        tie_rtol=1e-8,
        max_cells=2,
    )
    assert certificate.certified
    assert certificate.candidate_count == 11
    assert 1 <= certificate.exact_tie_count < certificate.near_tie_count
    assert certificate.near_tie_count == certificate.candidate_count
    assert certificate.tie_scope == "distinct_sampled_singleton_lower_bounds_only"


def test_adaptive_result_records_actual_ties_and_rounded_model_scope() -> None:
    result = gecp(
        FermionicKernel(1),
        config=GECPConfig(
            max_rank=1,
            pivot="adaptive",
            precision_bits=96,
            certificate_abs_tol=1e-2,
            max_cells=2_000,
        ),
    )
    assert result.rank == 1
    certificate = result.pivot_certificates[0]
    assert certificate.near_tie_count >= 2
    assert result.tie_counts == [certificate.near_tie_count]
    assert type(result.tie_counts[0]) is int
    assert certificate.objective_scope == "frozen_rounded_fermionic_residual_model"
    metadata = result.metadata
    assert metadata.tie_status == "recorded_sampled_candidates"
    assert metadata.tie_counts == result.tie_counts
    assert metadata.certification_scope.endswith("frozen_rounded_residual_model_only")
    assert not metadata.configuration["exact_selected_cross_certified"]
    assert not metadata.configuration["float64_evaluate_residual_certified"]
    metadata.canonical_json()


def test_block_diagnostic_uses_the_same_binary_coordinates_as_its_certificates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: list[tuple[mp.mpf, mp.mpf]] = []
    original = _FermionicResidualEvaluator.add_pivot

    def record(
        evaluator: _FermionicResidualEvaluator,
        t: mp.mpf,
        omega: mp.mpf,
        pivot: mp.mpf,
    ) -> None:
        captured.append((t, omega))
        original(evaluator, t, omega, pivot)

    monkeypatch.setattr(_FermionicResidualEvaluator, "add_pivot", record)
    evidence = certify_fermionic_block_contraction(
        1.3, 1, precision_bits=96, abs_tol=1e-2, max_cells=2000
    )
    assert evidence.completed_steps == 1
    assert captured[0] == (mp.mpf(0), mp.mpf(1.3))
    assert all(
        c.objective_scope == "frozen_rounded_fermionic_residual_model"
        for c in evidence.pivot_certificates
    )


def test_representable_large_rank_one_update_avoids_outer_product_overflow() -> None:
    matrix = np.full((2, 2), 1e308)
    with np.errstate(over="raise", invalid="raise"):
        result = gecp_matrix(
            matrix,
            config=GECPConfig(precision_bits=53, pivot="grid", max_rank=2, tol=1e294),
        )
    assert result.rank == 1
    assert result.converged
    assert result.residual_matrix is not None
    assert np.max(np.abs(result.residual_matrix)) == 0
    result.metadata.canonical_json()


def test_unrepresentable_float64_update_is_an_explicit_breakdown() -> None:
    with pytest.raises(NumericalBreakdownError, match="GECP update"):
        gecp_matrix(
            [[1e308, 1e308], [1e308, -1e308]],
            config=GECPConfig(precision_bits=53, pivot="grid", max_rank=2),
        )


@pytest.mark.parametrize("transpose", [False, True])
def test_mixed_scales_do_not_underflow_a_representable_schur_term(
    transpose: bool,
) -> None:
    matrix = np.array([[1e200, 1e200], [1e-200, 0.0]])
    if transpose:
        matrix = matrix.T
    result = gecp_matrix(
        matrix,
        config=GECPConfig(
            precision_bits=53, pivot="grid", max_rank=2, tol=1e-210, tie_rtol=0
        ),
    )
    assert result.converged and result.rank == 2
    assert result.row_indices == [0, 1]
    assert result.column_indices == [0, 1]
    assert abs(result.pivots[1]) == pytest.approx(1e-200, rel=2e-15, abs=0)
    # Direct binary64 products suffice for this independent exact determinant
    # check; LAPACK's unscaled determinant estimate can lose the tiny row.
    determinant = matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]
    assert determinant == pytest.approx(-1.0, rel=2e-15)
    assert np.prod(result.pivots) == pytest.approx(determinant, rel=2e-15)


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_near_tie_tolerance_never_reselects_an_annihilated_coordinate(
    precision_bits: int,
) -> None:
    result = gecp_matrix(
        [[1.0, 0.0], [0.0, 1e-20]],
        config=GECPConfig(
            precision_bits=precision_bits,
            pivot="grid",
            max_rank=2,
            tol=1e-25,
            tie_rtol=1e-8,
        ),
    )
    assert result.converged
    assert result.row_indices == result.column_indices == [0, 1]
    assert result.tie_counts == [1, 1]
    assert np.array_equal(result.residual_matrix, np.zeros((2, 2)))


def test_unrepresentable_tiny_schur_term_is_an_explicit_breakdown() -> None:
    with pytest.raises(NumericalBreakdownError, match="GECP update"):
        gecp_matrix(
            [[1.0, 1e-200], [1e-200, 0.0]],
            config=GECPConfig(precision_bits=53, pivot="grid", max_rank=2, tol=1e-310),
        )


@pytest.mark.parametrize("transpose", [False, True])
def test_roundoff_does_not_promote_an_exact_duplicate_core(transpose: bool) -> None:
    row = [0.7088646663651996, 0.9626470984418726]
    matrix = np.array([row, row])
    if transpose:
        matrix = matrix.T
    with pytest.raises(NumericalBreakdownError, match="exact duplicate"):
        gecp_matrix(
            matrix,
            config=GECPConfig(
                precision_bits=53, pivot="grid", max_rank=2, tol=1e-100, tie_rtol=0
            ),
        )


def test_failed_float64_core_diagnostic_is_an_explicit_breakdown(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed_svd(*args: object, **kwargs: object) -> object:
        raise np.linalg.LinAlgError("independent diagnostic failure")

    monkeypatch.setattr(np.linalg, "svd", failed_svd)
    with pytest.raises(NumericalBreakdownError, match="core diagnostic failed"):
        gecp_matrix([[1.0]], config=GECPConfig(precision_bits=53, pivot="grid"))


@pytest.mark.parametrize("t", [0.0, 0.125, 0.5, 0.875, 1.0])
@pytest.mark.parametrize("omega", [-1000.0, -19.0, -0.25, 0.0, 0.25, 19.0, 1000.0])
def test_centered_kernel_formula_against_independent_high_precision_expression(
    t: float, omega: float
) -> None:
    kernel = FermionicKernel(1000)
    with mp.workdps(100):
        x, w = mp.mpf(t), mp.mpf(omega)
        expected = float(mp.exp(-(x - mp.mpf("0.5")) * w) / (2 * mp.cosh(w / 2)))
    assert kernel(t, omega) == pytest.approx(expected, rel=2e-14, abs=0.0)
    assert kernel.centered(t, omega) == pytest.approx(expected, rel=2e-14, abs=0.0)
