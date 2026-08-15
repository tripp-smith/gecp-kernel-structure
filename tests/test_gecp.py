import numpy as np
import pytest

from kernelgecp import (
    FermionicKernel,
    GECPConfig,
    cross_approximation,
    evaluate_residual,
    gecp,
    gecp_matrix,
)


def test_rank_two_interpolation_and_determinant_product() -> None:
    matrix = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]])
    result = gecp_matrix(matrix, config=GECPConfig(tol=1e-14, max_rank=3))
    assert result.converged
    assert result.rank == 3
    assert result.residual_matrix is not None
    assert np.max(np.abs(result.residual_matrix[result.row_indices, :])) < 1e-12
    assert np.max(np.abs(result.residual_matrix[:, result.column_indices])) < 1e-12
    products = np.cumprod(result.pivots)
    assert result.core_det_history == pytest.approx(products)


def test_iterative_and_direct_cross_residual_agree() -> None:
    kernel = FermionicKernel(20)
    t = np.linspace(0, 1, 17)
    omega = np.linspace(-20, 20, 31)
    result = gecp(
        kernel,
        config=GECPConfig(tol=1e-9, max_rank=8, pivot="grid"),
        t_grid=t,
        omega_grid=omega,
    )
    direct = evaluate_residual(kernel, result, t, omega)
    assert result.residual_matrix is not None
    assert direct == pytest.approx(result.residual_matrix, abs=2e-12)
    approximation = cross_approximation(
        kernel, result.t_nodes, result.omega_nodes, t, omega
    )
    assert approximation.shape == (t.size, omega.size)


def test_rank_budget_is_reported() -> None:
    result = gecp_matrix(
        np.eye(4), config=GECPConfig(tol=1e-12, max_rank=2, pivot="grid")
    )
    assert not result.converged
    assert result.stop_reason == "max_rank"
    assert result.rank == 2


def test_ties_are_lexicographic_and_recorded() -> None:
    result = gecp_matrix(np.eye(3), config=GECPConfig(max_rank=1, pivot="grid"))
    assert result.row_indices == [0]
    assert result.column_indices == [0]
    assert result.tie_counts == [3]
