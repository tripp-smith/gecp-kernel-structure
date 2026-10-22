"""Independent checks of numerical evidence and explicit non-certification scopes."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

from kernelgecp import (
    FermionicKernel,
    GECPConfig,
    NumericalBreakdownError,
    gecp,
    gecp_matrix,
    pivoted_cholesky,
    sparse_representation,
)
from kernelgecp.types import FloatArray, NumericalMetadata


def test_cholesky_reports_actual_precision_and_independent_core_conditioning() -> None:
    matrix = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]])
    result = pivoted_cholesky(
        matrix, config=GECPConfig(precision_bits=128, max_rank=3, tol=1e-12)
    )
    metadata = result.metadata
    assert result.converged
    assert metadata.arithmetic == "float64"
    assert metadata.precision_bits == 53
    assert metadata.requested_precision_bits == 128
    assert metadata.precision_status == "float64_only"
    assert metadata.conditioning_status == "estimated"
    assert metadata.certification_status == "not_certified"
    assert metadata.convergence_scope == (
        "original_sampled_matrix_reconstruction_infinity_error"
    )
    assert len(metadata.conditioning_history) == result.rank
    for count, record in enumerate(metadata.conditioning_history, start=1):
        selected = result.indices[:count]
        singular = np.linalg.svd(matrix[np.ix_(selected, selected)], compute_uv=False)
        assert record.stage == "selected_core"
        assert record.status == "estimated"
        assert record.columns == record.effective_rank == count
        assert record.sigma_min == pytest.approx(singular[-1])
        assert record.condition_number == pytest.approx(singular[0] / singular[-1])
    assert result.factor @ result.factor.T + result.residual_matrix == pytest.approx(
        matrix, abs=1e-14
    )
    actual_error = np.max(np.abs(matrix - result.factor @ result.factor.T))
    assert result.metadata.configuration["returned_error"] == actual_error


def test_cholesky_small_symmetry_defect_is_reported_not_hidden_by_projection() -> None:
    matrix = np.array([[2.0, 1.0], [1.0 + 2e-13, 3.0]])
    result = pivoted_cholesky(matrix, config=GECPConfig(max_rank=2, tol=1e-14))
    actual = matrix - result.factor @ result.factor.T
    assert np.array_equal(result.residual_matrix, actual)
    assert result.metadata.configuration["symmetry_error"] == np.max(
        np.abs(matrix - matrix.T)
    )
    assert result.metadata.configuration["returned_error"] == np.max(np.abs(actual))
    assert not result.converged
    assert result.stop_reason == "roundoff_or_symmetry_limit"


def test_cholesky_tiny_tolerance_uses_returned_error_not_working_diagonal() -> None:
    matrix = np.array([[2.0, 1.0], [1.0, 2.0]])
    result = pivoted_cholesky(matrix, config=GECPConfig(max_rank=2, tol=1e-30))
    actual_error = float(np.max(np.abs(matrix - result.factor @ result.factor.T)))
    assert result.converged is (actual_error <= 1e-30)
    assert result.metadata.configuration["returned_error"] == actual_error
    assert result.residual_matrix == pytest.approx(
        matrix - result.factor @ result.factor.T, abs=0
    )


def test_cholesky_large_finite_spd_input_avoids_symmetrization_overflow() -> None:
    matrix = np.eye(2) * 1e308
    tolerance = 1e294
    result = pivoted_cholesky(
        matrix,
        config=GECPConfig(max_rank=2, tol=tolerance, precision_bits=53),
    )
    actual_error = float(np.max(np.abs(matrix - result.factor @ result.factor.T)))
    assert np.all(np.isfinite(result.factor))
    assert np.all(np.isfinite(result.residual_matrix))
    assert result.converged is (actual_error <= tolerance)
    assert result.converged
    assert result.metadata.configuration["returned_error"] == actual_error
    result.metadata.canonical_json()


@pytest.mark.parametrize("tie_rtol,first,count", [(1e-8, 0, 2), (0.0, 1, 1)])
def test_cholesky_near_ties_respect_explicit_zero_tolerance(
    tie_rtol: float, first: int, count: int
) -> None:
    result = pivoted_cholesky(
        np.diag([1.0, 1.0 + 1e-10]),
        config=GECPConfig(max_rank=1, tie_rtol=tie_rtol, precision_bits=53),
    )
    assert result.indices == [first]
    assert result.metadata.tie_counts == [count]
    assert result.metadata.tie_rtol == tie_rtol
    assert result.metadata.precision_status == "as_requested"
    assert not result.converged
    assert result.stop_reason == "max_rank"


@pytest.mark.parametrize(
    "matrix,max_rank,converged,reason",
    [
        (np.empty((0, 0)), 3, True, "tolerance"),
        (np.zeros((3, 3)), 3, True, "tolerance"),
        (np.eye(3), 0, False, "max_rank"),
    ],
)
def test_cholesky_rank_zero_and_empty_results_have_explicit_metadata(
    matrix: FloatArray, max_rank: int, converged: bool, reason: str
) -> None:
    result = pivoted_cholesky(matrix, config=GECPConfig(max_rank=max_rank))
    assert result.rank == 0
    assert result.indices == result.pivots == []
    assert result.factor.shape == (matrix.shape[0], 0)
    assert result.converged is converged
    assert result.stop_reason == reason
    assert result.metadata.conditioning_status == "not_applicable"
    assert result.metadata.conditioning_history == []
    assert result.metadata.tie_counts == []
    assert result.metadata.certification_status == "not_certified"
    result.metadata.canonical_json()


@pytest.mark.parametrize(
    "matrix", [np.array([[1.0, 2.0], [2.0, 1.0]]), np.array([[0.0, 1.0], [1.0, 0.0]])]
)
def test_cholesky_rejects_materially_indefinite_residual(matrix: FloatArray) -> None:
    with pytest.raises(NumericalBreakdownError, match="negative diagonal|PSD"):
        pivoted_cholesky(matrix, config=GECPConfig(max_rank=2))


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_cholesky_rejects_nonfinite_input(value: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        pivoted_cholesky(np.diag([1.0, value]))


def test_conditioning_estimator_failure_is_not_a_factorization_certificate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed_svd(*_args: object, **_kwargs: object) -> None:
        raise np.linalg.LinAlgError("deliberate diagnostic failure")

    monkeypatch.setattr("kernelgecp.diagnostics.np.linalg.svd", failed_svd)
    result = pivoted_cholesky(np.eye(2), config=GECPConfig(max_rank=2))
    assert result.converged
    assert result.metadata.conditioning_status == "failed"
    assert all(
        record.status == "failed" for record in result.metadata.conditioning_history
    )
    assert all(
        record.sigma_min is None for record in result.metadata.conditioning_history
    )
    assert result.metadata.certification_status == "not_certified"
    result.metadata.canonical_json()


def test_sparse_zero_target_uses_no_atoms_or_optimizer() -> None:
    result = sparse_representation(
        np.linspace(0, 1, 9), np.zeros(9), cutoff=2, candidate_frequencies=[-1, 1]
    )
    assert result.converged
    assert result.stop_reason == "tolerance"
    assert result.iterations == result.frequencies.size == result.weights.size == 0
    assert result.validation_error == result.training_error == 0
    assert result.metadata.tie_status == "not_applicable"
    assert result.metadata.tie_counts == []
    assert result.metadata.conditioning_status == "not_applicable"
    assert result.metadata.refinement.status == "not_applicable"
    assert result.metadata.certification_status == "not_certified"
    result.metadata.canonical_json()


def test_sparse_exhausts_library_without_reselecting_an_atom() -> None:
    t = np.linspace(0, 1, 21)
    candidates = np.array([-1.0, 1.0])
    values = np.sin(8 * t)
    result = sparse_representation(
        t,
        values,
        cutoff=2,
        tolerance=1e-14,
        max_atoms=10,
        candidate_frequencies=candidates,
        refine=False,
    )
    assert not result.converged
    assert result.stop_reason == "candidate_exhausted"
    assert result.iterations == 2
    assert np.unique(result.frequencies).size == 2
    assert sorted(result.frequencies) == sorted(candidates)
    assert len(result.metadata.tie_counts) == result.iterations
    design = FermionicKernel(2)(t[:, None], result.frequencies[None, :])
    singular = np.linalg.svd(design, compute_uv=False)
    record = result.metadata.conditioning_history[-1]
    assert record.stage == "returned_design"
    assert record.sigma_min == pytest.approx(singular[-1])
    assert record.condition_number == pytest.approx(singular[0] / singular[-1])
    assert result.training_error == pytest.approx(
        np.max(np.abs(design @ result.weights - values))
    )
    assert result.metadata.refinement.status == "not_requested"


def test_sparse_rank_deficient_design_has_no_infinite_condition_number() -> None:
    result = sparse_representation(
        [0.0, 0.0, 0.0],
        [-1.0, 0.0, 1.0],
        cutoff=2,
        max_atoms=2,
        candidate_frequencies=[-1.0, 1.0],
        refine=False,
    )
    assert not result.converged
    assert result.metadata.conditioning_status == "rank_deficient"
    record = result.metadata.conditioning_history[-1]
    assert record.status == "rank_deficient"
    assert record.columns == 2
    assert record.effective_rank == 1
    assert record.condition_number is None
    result.metadata.canonical_json()


@pytest.mark.parametrize("tie_rtol", [0.0, 1e-8])
def test_sparse_exact_ties_record_deterministic_first_available_atom(
    tie_rtol: float,
) -> None:
    result = sparse_representation(
        [0.5],
        [1.0],
        cutoff=2,
        candidate_frequencies=[-2.0, 2.0],
        refine=False,
        tie_rtol=tie_rtol,
    )
    assert result.frequencies.tolist() == [-2.0]
    assert result.metadata.tie_counts == [2]
    assert result.metadata.tie_rtol == tie_rtol


@pytest.fixture
def off_grid_sparse_problem() -> tuple[FloatArray, FloatArray, FloatArray]:
    t = np.linspace(0, 1, 31)
    values = np.asarray(FermionicKernel(2)(t, 0.2)) * 0.8
    return t, values, np.array([-1.0, 1.0])


def test_sparse_optimizer_budget_and_representation_convergence_are_distinct(
    monkeypatch: pytest.MonkeyPatch,
    off_grid_sparse_problem: tuple[FloatArray, FloatArray, FloatArray],
) -> None:
    t, values, candidates = off_grid_sparse_problem

    def budget_exhausted(
        _objective: object, parameters: FloatArray, **kwargs: object
    ) -> SimpleNamespace:
        assert kwargs["max_nfev"] == 1
        assert parameters.size == 2
        return SimpleNamespace(
            x=np.array([0.2, 0.8]),
            success=False,
            status=0,
            nfev=1,
            message="maximum function evaluations exceeded",
            jac=np.zeros((t.size, 2)),
        )

    monkeypatch.setattr("kernelgecp.sparse.least_squares", budget_exhausted)
    result = sparse_representation(
        t,
        values,
        cutoff=2,
        max_atoms=1,
        tolerance=1e-12,
        candidate_frequencies=candidates,
        refinement_max_nfev=1,
    )
    assert result.converged
    assert result.stop_reason == "tolerance"
    refinement = result.metadata.refinement
    assert refinement.status == "accepted"
    assert refinement.accepted
    assert refinement.optimizer_success is False
    assert refinement.optimizer_status == 0
    assert refinement.evaluations == refinement.max_evaluations == 1
    assert result.training_error <= 1e-12
    assert result.metadata.certification_status == "not_certified"


@pytest.mark.parametrize("mode", ["worse", "nonfinite", "failure"])
def test_sparse_rejects_bad_refinement_and_reports_returned_error(
    monkeypatch: pytest.MonkeyPatch,
    off_grid_sparse_problem: tuple[FloatArray, FloatArray, FloatArray],
    mode: str,
) -> None:
    t, values, candidates = off_grid_sparse_problem
    baseline = sparse_representation(
        t, values, cutoff=2, max_atoms=1, candidate_frequencies=candidates, refine=False
    )

    def bad_refinement(
        _objective: object, parameters: FloatArray, **_kwargs: object
    ) -> SimpleNamespace:
        if mode == "failure":
            raise ValueError("deliberate optimizer failure")
        proposed = parameters.copy()
        proposed[-1] = 0 if mode == "worse" else np.nan
        return SimpleNamespace(
            x=proposed,
            success=True,
            status=1,
            nfev=2,
            message="optimizer status alone is not an error guarantee",
            jac=np.zeros((t.size, parameters.size)),
        )

    monkeypatch.setattr("kernelgecp.sparse.least_squares", bad_refinement)
    result = sparse_representation(
        t, values, cutoff=2, max_atoms=1, candidate_frequencies=candidates
    )
    assert result.frequencies == pytest.approx(baseline.frequencies)
    assert result.weights == pytest.approx(baseline.weights)
    assert result.training_error == pytest.approx(baseline.training_error)
    assert not result.metadata.refinement.accepted
    assert result.metadata.refinement.status == (
        "failed" if mode == "failure" else "rejected_worse_or_nonfinite"
    )
    actual_error = float(
        np.max(
            np.abs(
                FermionicKernel(2)(t[:, None], result.frequencies[None, :])
                @ result.weights
                - values
            )
        )
    )
    assert result.validation_error == result.training_error == actual_error
    assert result.converged is bool(
        actual_error <= 1e-8 * max(1, np.max(np.abs(values)))
    )
    result.metadata.canonical_json()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"t": [np.nan]},
        {"t": [np.inf]},
        {"t": [-0.1]},
        {"values": [np.nan]},
        {"values": [np.inf]},
        {"cutoff": np.nan},
        {"cutoff": np.inf},
        {"tolerance": np.nan},
        {"tolerance": np.inf},
        {"tie_rtol": np.nan},
        {"tie_rtol": np.inf},
        {"tie_rtol": -1.0},
        {"refinement_max_nfev": 0},
        {"candidate_frequencies": [np.nan]},
        {"candidate_frequencies": [np.inf]},
        {"candidate_frequencies": [3.0]},
        {"candidate_frequencies": [1.0, 1.0]},
        {"candidate_frequencies": []},
    ],
)
def test_sparse_rejects_invalid_inputs(kwargs: dict[str, object]) -> None:
    arguments: dict[str, object] = {
        "t": [0.5],
        "values": [1.0],
        "cutoff": 2,
        "candidate_frequencies": [-1.0, 1.0],
    }
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        sparse_representation(**arguments)  # type: ignore[arg-type]


def test_gecp_metadata_view_distinguishes_grid_from_continuous_evidence() -> None:
    grid = gecp_matrix(
        [[2.0, 0.0], [0.0, 1.0]],
        config=GECPConfig(max_rank=2, precision_bits=96),
    )
    metadata = grid.metadata
    assert metadata.algorithm == "gecp_grid"
    assert metadata.arithmetic == "mpmath"
    assert metadata.precision_bits == metadata.requested_precision_bits == 96
    assert metadata.tie_status == "recorded"
    assert metadata.tie_counts == grid.tie_counts
    assert metadata.tie_rtol == 2.0 ** (5 - 96)
    assert metadata.convergence_scope == "sampled_grid"
    assert metadata.certification_status == "not_certified"
    assert "no_continuous_certificate" in metadata.certification_scope
    assert all(
        record.sigma_min_decimal is not None for record in metadata.conditioning_history
    )
    adaptive = gecp(
        FermionicKernel(1),
        config=GECPConfig(
            max_rank=0,
            tol=2,
            pivot="adaptive",
            precision_bits=96,
            certificate_abs_tol=1e-2,
            max_cells=2_000,
        ),
    )
    assert adaptive.converged
    assert adaptive.metadata.algorithm == "gecp_interval_adaptive"
    assert adaptive.metadata.tie_status == "not_estimated"
    assert adaptive.metadata.tie_counts == []
    assert adaptive.metadata.convergence_scope == "continuous_interval_bound"
    assert adaptive.metadata.certification_status == "certified"
    assert (
        adaptive.metadata.certification_scope == "outward_interval_residual_upper_bound"
    )


@pytest.mark.parametrize("precision,first,count", [(53, 0, 2), (96, 1, 1)])
def test_gecp_metadata_records_backend_effective_default_tie_tolerance(
    precision: int, first: int, count: int
) -> None:
    result = gecp_matrix(
        [[1.0, 0.0], [0.0, 1.0 + 5e-15]],
        config=GECPConfig(max_rank=1, precision_bits=precision, pivot="grid"),
    )
    assert result.row_indices == result.column_indices == [first]
    assert result.metadata.tie_counts == [count]
    effective_tolerance = (
        32 * np.finfo(float).eps if precision == 53 else 2.0 ** (5 - precision)
    )
    assert result.metadata.tie_rtol == effective_tolerance
    assert result.metadata.configuration["tie_rtol"] is None


def test_uncertified_adaptive_metadata_does_not_promote_budget_exhaustion() -> None:
    result = gecp(
        FermionicKernel(1),
        config=GECPConfig(
            max_rank=1,
            pivot="adaptive",
            precision_bits=96,
            certificate_abs_tol=0,
            certificate_rel_tol=0,
            max_cells=1,
        ),
    )
    assert not result.converged
    assert result.stop_reason == "uncertified_pivot"
    assert result.metadata.certification_status == "not_certified"
    assert result.metadata.tie_status == "not_estimated"
    assert result.metadata.conditioning_status == "not_applicable"
    result.metadata.canonical_json()


def test_metadata_serialization_is_finite_deterministic_and_input_sensitive() -> None:
    t = np.linspace(0, 1, 11)
    options = {"cutoff": 2, "candidate_frequencies": [-1, 1], "refine": False}
    first = sparse_representation(t, np.sin(t), **options)
    repeated = sparse_representation(t, np.sin(t), **options)
    changed = sparse_representation(t, np.cos(t), **options)
    encoded = first.metadata.canonical_json()
    assert encoded == repeated.metadata.canonical_json()
    assert encoded == json.dumps(
        first.metadata.as_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    assert json.loads(encoded)["precision_bits"] == 53
    assert (
        first.metadata.configuration["target_fingerprint"]
        != changed.metadata.configuration["target_fingerprint"]
    )
    assert (
        first.metadata.configuration["candidate_fingerprint"]
        == (changed.metadata.configuration["candidate_fingerprint"])
    )
    with pytest.raises(ValueError):
        NumericalMetadata(configuration={"invalid": np.nan}).canonical_json()


def test_manual_metadata_defaults_do_not_claim_measured_execution() -> None:
    metadata = NumericalMetadata()
    assert metadata.algorithm == metadata.arithmetic == "not_recorded"
    assert metadata.precision_bits is None
    assert metadata.conditioning_status == metadata.tie_status == "not_recorded"
    assert metadata.certification_status == "not_certified"
    assert metadata.refinement.status == "not_applicable"
