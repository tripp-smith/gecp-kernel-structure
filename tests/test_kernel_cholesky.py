"""Independent dense and arithmetic checks for the column-oracle API."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import cg

from kernelgecp import GECPConfig, pivoted_cholesky
from kernelgecp.kernel_cholesky import (
    KernelCholeskyExtension,
    KernelCholeskyFactor,
    KernelSolveResult,
)
from kernelgecp.types import NumericalBreakdownError


def _from_dense(matrix, **kwargs):
    return KernelCholeskyFactor.from_columns(
        np.diag(matrix), lambda j: matrix[:, j], **kwargs
    )


def test_column_factor_matches_independent_dense_pivoting_and_error():
    # Deliberately asymmetric geometry/amplitude gives resolved pivot margins.
    # Symmetric near-tie trajectories can differ between binary64 algorithms.
    points = np.sort(np.random.default_rng(713).uniform(-2, 2, 31))
    amplitude = np.linspace(0.8, 1.2, points.size)
    matrix = (
        amplitude[:, None]
        * np.exp(-((points[:, None] - points[None, :]) ** 2))
        * amplitude
    )
    calls = []

    def column(index):
        calls.append(index)
        return matrix[:, index]

    factor = KernelCholeskyFactor.from_columns(
        np.diag(matrix), column, max_rank=12, atol=1e-12
    )
    dense = pivoted_cholesky(
        matrix, config=GECPConfig(max_rank=12, tol=1e-12, precision_bits=53)
    )
    assert factor.indices == tuple(dense.indices)
    np.testing.assert_allclose(factor.factor, dense.factor, atol=1e-12)
    np.testing.assert_allclose(
        factor.residual_diagonal,
        np.diag(matrix - factor.factor @ factor.factor.T),
        atol=5e-16,
    )
    assert calls == list(factor.indices)
    assert factor.diagnostics.column_evaluations == factor.rank
    assert factor.diagnostics.stop_reason == "max_rank"
    assert not factor.diagnostics.converged
    assert factor.metadata.certification_status == "not_certified"
    assert factor.metadata.precision_bits == 53


def test_symmetric_roundoff_ties_are_admissible_against_dense_schur_solves():
    # Linux NumPy and the independent dense implementation select opposite
    # members of this symmetric pair after nine pivots. Check each choice by
    # direct principal-block solves rather than requiring bit-identical paths.
    points = np.linspace(-2, 2, 31)
    matrix = np.exp(-((points[:, None] - points[None, :]) ** 2))
    factor = _from_dense(matrix, max_rank=12, atol=1e-12)
    selected = []
    for pivot in factor.indices:
        diagonal = np.diag(matrix).copy()
        if selected:
            columns = matrix[:, selected]
            solved = np.linalg.solve(matrix[np.ix_(selected, selected)], columns.T)
            diagonal -= np.sum(columns * solved.T, axis=1)
        diagonal[selected] = -np.inf
        # Independent Schur subtraction accumulates absolute roundoff relative
        # to the original matrix; it need not reproduce a tiny residual's ulps.
        assert np.max(diagonal) - diagonal[pivot] <= 512 * np.finfo(float).eps
        selected.append(pivot)
    residual = matrix - factor.factor @ factor.factor.T
    np.testing.assert_allclose(np.diag(residual), factor.residual_diagonal, atol=5e-16)
    assert np.linalg.eigvalsh(residual)[0] >= -2e-14
    assert np.max(np.abs(residual)) <= np.max(factor.residual_diagonal) + 2e-14


def test_rank_deficiency_duplicates_and_zero_rank_have_explicit_stopping():
    features = np.array([[1, 0], [1, 0], [0, 1], [1, 1], [0, 0.0]])
    matrix = features @ features.T
    factor = _from_dense(matrix, max_rank=5, atol=1e-13)
    assert factor.rank == 2
    assert factor.diagnostics.converged
    np.testing.assert_allclose(factor.factor @ factor.factor.T, matrix, atol=1e-14)
    assert factor.indices == (3, 0)
    assert factor.diagnostics.tie_counts == (1, 3)
    for diagonal in ([], [0, 0]):
        empty = KernelCholeskyFactor.from_columns(
            diagonal, lambda _: pytest.fail("no column should be queried"), max_rank=0
        )
        assert empty.rank == 0
        assert empty.diagnostics.converged
        assert empty.factor.shape == (len(diagonal), 0)
    capped = _from_dense(np.eye(3), max_rank=0)
    assert capped.rank == 0
    assert not capped.diagnostics.converged
    assert capped.diagnostics.stop_reason == "max_rank"


def test_near_ties_choose_first_index_but_small_scale_is_not_an_absolute_tie():
    diagonal = np.array([1 - 2e-15, 1, 0.2])
    factor = _from_dense(np.diag(diagonal), max_rank=1)
    assert factor.indices == (0,)
    assert factor.diagnostics.tie_counts == (2,)
    small = _from_dense(np.diag([1e-20, 2e-20]), atol=0, max_rank=1)
    assert small.indices == (1,)


def test_factor_owns_inputs_even_when_callback_mutates_caller_diagonal():
    diagonal = np.ones(3)
    returned_column = np.ones(3)

    def column(_):
        diagonal[:] = 0
        return returned_column

    factor = KernelCholeskyFactor.from_columns(diagonal, column, max_rank=3)
    returned_column[:] = 0
    np.testing.assert_array_equal(factor.diagonal, 1)
    np.testing.assert_array_equal(factor.factor, np.ones((3, 1)))
    for array in (factor.diagonal, factor.residual_diagonal, factor.factor):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(AttributeError):
        factor.indices = ()
    metadata = factor.metadata
    metadata.tie_counts.clear()
    assert factor.metadata.tie_counts == [3]


def test_public_extension_and_solve_record_constructors_take_immutable_ownership():
    features = np.ones((3, 2))
    diagonal = np.ones(3)
    extension = KernelCholeskyExtension(features, diagonal)
    solution = KernelSolveResult(diagonal, 1e-14, 0, "independent")
    features[:] = 0
    diagonal[:] = 0
    np.testing.assert_array_equal(extension.features, 1)
    np.testing.assert_array_equal(extension.residual_diagonal, 1)
    np.testing.assert_array_equal(solution.solution, 1)
    for array in (extension.features, extension.residual_diagonal, solution.solution):
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_clipping_negative_roundoff_does_not_establish_tighter_convergence():
    factor = _from_dense(np.diag([2.0]), max_rank=1, atol=0)
    assert factor.residual_diagonal[0] == 0
    assert factor.diagnostics.returned_diagonal_error > 0
    assert not factor.converged
    assert factor.stop_reason == "roundoff_limit"


def test_out_of_sample_features_and_gp_predictions_match_independent_dense_scipy():
    rng = np.random.default_rng(941)
    train = rng.normal(size=(17, 3))
    test = rng.normal(size=(7, 3))
    matrix, cross = train @ train.T, test @ train.T
    factor = _from_dense(matrix, max_rank=6, atol=1e-12)
    extension = factor.extend(
        cross[:, factor.indices], diagonal=np.sum(test**2, axis=1)
    )
    np.testing.assert_allclose(extension.features @ factor.factor.T, cross, atol=2e-14)
    np.testing.assert_allclose(extension.residual_diagonal, 0, atol=2e-14)
    y = rng.normal(size=17)
    noise = 0.3
    chol = cho_factor(matrix + noise * np.eye(17), lower=True)
    expected_mean = cross @ cho_solve(chol, y)
    expected_variance = np.sum(test**2, axis=1) - np.sum(
        cross * cho_solve(chol, cross.T).T, axis=1
    )
    shifted = factor.shifted(noise)
    new_cross = extension.features @ factor.factor.T
    actual_mean = new_cross @ shifted.solve(y)
    actual_variance = np.sum(test**2, axis=1) - np.sum(
        new_cross * shifted.solve(new_cross.T).T, axis=1
    )
    np.testing.assert_allclose(actual_mean, expected_mean, atol=1e-13)
    np.testing.assert_allclose(actual_variance, expected_variance, atol=2e-13)
    with pytest.raises(ValueError):
        extension.features.setflags(write=True)
    with pytest.raises(NumericalBreakdownError, match="out-of-sample"):
        factor.extend(100 * cross[:, factor.indices], diagonal=np.ones(7))


@pytest.mark.parametrize("vector", [False, True])
def test_shifted_solve_operator_and_preconditioner_match_dense_scipy(vector):
    rng = np.random.default_rng(2143)
    root = rng.normal(size=(23, 5))
    factor = _from_dense(root @ root.T, max_rank=7, atol=1e-12)
    noise = np.linspace(0.1, 0.7, 23)
    rhs = rng.normal(size=(23, 3))
    if vector:
        rhs = rhs[:, 0]
    matrix = root @ root.T + np.diag(noise)
    expected = cho_solve(cho_factor(matrix, lower=True), rhs)
    shifted = factor.shifted(noise)
    result = shifted.solve_with_diagnostics(rhs)
    np.testing.assert_allclose(result.solution, expected, atol=3e-12)
    assert result.relative_residual < 1e-10
    operator = factor.as_linear_operator(noise)
    np.testing.assert_allclose(operator @ rhs, matrix @ rhs, atol=2e-13)
    np.testing.assert_allclose(operator.T @ rhs, matrix @ rhs, atol=2e-13)
    np.testing.assert_allclose(shifted.as_preconditioner() @ rhs, expected, atol=3e-12)
    np.testing.assert_array_equal(shifted.solve(np.zeros_like(rhs)), 0)


def test_original_kernel_preconditioned_cg_uses_the_column_factor():
    points = np.linspace(-2, 2, 73)
    kernel = np.exp(-((points[:, None] - points[None, :]) ** 2))
    factor = _from_dense(kernel, max_rank=7, atol=0)
    matrix = kernel + 0.01 * np.eye(73)
    rhs = np.sin(points)
    solution, info = cg(
        matrix, rhs, M=factor.shifted(0.01).as_preconditioner(), rtol=1e-10, atol=0
    )
    assert info == 0
    assert np.linalg.norm(matrix @ solution - rhs) / np.linalg.norm(rhs) < 1e-8


def test_tiny_noise_fallback_avoids_woodbury_cancellation_without_dense_state():
    factor = _from_dense(np.diag([1.0, 2.0, 3.0]), max_rank=3, atol=0)
    result = factor.shifted(1e-24, max_refinements=0).solve_with_diagnostics([1, 4, 9])
    assert result.method == "implicit_householder_svd"
    np.testing.assert_allclose(result.solution, [1, 2, 3], atol=2e-15)
    assert result.relative_residual < 1e-14


def test_unrepresentable_final_subnormal_solution_does_not_pass_normalized_gate():
    factor = _from_dense(np.zeros((1, 1)), max_rank=0)
    with pytest.raises(NumericalBreakdownError, match="rescaling"):
        factor.shifted(2).solve([np.nextafter(0, 1)])


def test_constructor_does_not_request_all_columns_or_allocate_dense_kernel(monkeypatch):
    size, rank = 1009, 4
    points = np.linspace(-2, 2, size)
    calls = []
    zeros = np.zeros

    def guarded_zeros(shape, *args, **kwargs):
        assert shape != (size, size), "attempted full covariance/residual allocation"
        return zeros(shape, *args, **kwargs)

    def column(j):
        calls.append(j)
        return np.exp(-((points - points[j]) ** 2))

    monkeypatch.setattr(np, "zeros", guarded_zeros)
    factor = KernelCholeskyFactor.from_columns(np.ones(size), column, max_rank=rank)
    assert len(calls) == factor.rank == rank
    assert factor.factor.shape == (size, rank)
    assert factor.storage_nbytes < 8 * size * (rank + 3)
    assert factor.shifted(0.1).storage_nbytes < 8 * size * (rank + 4)


@pytest.mark.parametrize("diagonal", [[True], [np.nan], [np.inf], [-1], [[1]], [1j]])
def test_invalid_diagonal_is_rejected(diagonal):
    with pytest.raises(ValueError):
        KernelCholeskyFactor.from_columns(diagonal, lambda _: [1])


@pytest.mark.parametrize(
    "column", [[1], [1, np.nan], [True, False], [1j, 0], [[1], [0]]]
)
def test_invalid_callback_result_is_rejected(column):
    with pytest.raises(ValueError):
        KernelCholeskyFactor.from_columns([1, 1], lambda _: column, max_rank=1)


def test_inconsistent_callbacks_and_non_psd_columns_fail_explicitly():
    with pytest.raises(ValueError, match="diagonal"):
        KernelCholeskyFactor.from_columns([1, 1], lambda _: [0, 0])
    with pytest.raises(NumericalBreakdownError, match="PSD entry"):
        KernelCholeskyFactor.from_columns([1, 1], lambda _: [1, 2])
    nonsymmetric = np.array([[1, 0.1], [0.5, 1]])
    with pytest.raises(ValueError, match="symmetric"):
        _from_dense(nonsymmetric, max_rank=2, atol=0)
    indefinite = np.array([[1, 0.9, 0.9], [0.9, 1, -0.9], [0.9, -0.9, 1]])
    with pytest.raises(NumericalBreakdownError, match="negative residual"):
        _from_dense(indefinite, max_rank=3, atol=0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_rank": True},
        {"max_rank": -1},
        {"max_rank": 1.5},
        {"atol": True},
        {"atol": -1},
        {"rtol": np.inf},
        {"tie_rtol": 1},
        {"psd_rtol": 1},
        {"psd_rtol": -1},
    ],
)
def test_invalid_factor_controls_are_not_coerced(kwargs):
    with pytest.raises(ValueError):
        _from_dense(np.eye(2), **kwargs)


@pytest.mark.parametrize("shift", [True, 0, -1, np.inf, [1], [1, 0], [1, np.nan]])
def test_invalid_shifted_solve_inputs_are_rejected(shift):
    factor = _from_dense(np.eye(2), max_rank=2)
    with pytest.raises(ValueError):
        factor.shifted(shift)


def test_invalid_extension_rhs_and_arithmetic_fail_clearly():
    factor = _from_dense(np.eye(2), max_rank=1)
    for cross in ([1], [[1, 2]], [[True]], [[np.nan]]):
        with pytest.raises(ValueError):
            factor.extend(cross)
    for rhs in ([True, False], [1], np.ones((2, 1, 1)), [1j, 0]):
        with pytest.raises(ValueError):
            factor.shifted(1).solve(rhs)
    with pytest.raises(NumericalBreakdownError):
        factor.shifted(np.nextafter(0, 1))


def test_default_package_import_does_not_load_optional_gp_provider():
    code = (
        "import kernelgecp; import kernelgecp.kernel_cholesky; import sys; "
        "assert 'torch' not in sys.modules; assert 'gpytorch' not in sys.modules; "
        "assert 'linear_operator' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def test_actual_cpu_gpytorch_uses_identical_factor_and_owns_conversion():
    torch = pytest.importorskip("torch")
    gpytorch = pytest.importorskip("gpytorch")
    pytest.importorskip("linear_operator")
    root = np.array([[1, 0.2], [0.4, -0.3], [-0.2, 0.8], [0.7, 0.1]])
    factor = _from_dense(root @ root.T, max_rank=2, atol=1e-12)
    shift = np.array([0.1, 0.2, 0.15, 0.3])
    rhs = np.array([0.3, -0.8, 1.1, 0.7])
    converted = factor.to_gpytorch(shift)
    assert converted.device.type == "cpu"
    assert converted.dtype == torch.float64
    assert not converted.requires_grad
    np.testing.assert_allclose(
        converted.to_dense().numpy(), factor.as_linear_operator(shift) @ np.eye(4)
    )
    np.testing.assert_allclose(
        converted.solve(torch.tensor(rhs)).numpy(),
        factor.shifted(shift).solve(rhs),
        atol=1e-11,
    )
    distribution = gpytorch.distributions.MultivariateNormal(
        torch.zeros(4, dtype=torch.float64), converted
    )
    expected_matrix = root @ root.T + np.diag(shift)
    expected_logdet = np.linalg.slogdet(expected_matrix)[1]
    expected_logprob = -0.5 * (
        rhs @ np.linalg.solve(expected_matrix, rhs)
        + expected_logdet
        + 4 * np.log(2 * np.pi)
    )
    assert float(distribution.log_prob(torch.tensor(rhs))) == pytest.approx(
        expected_logprob, abs=1e-11
    )
    before = factor.factor.copy()
    converted.representation()[0][0, 0] = 99
    np.testing.assert_array_equal(factor.factor, before)
