import numpy as np
import pytest

from kernelgecp import GECPConfig, gecp_matrix, pivoted_cholesky


def test_psd_gecp_matches_canonical_cholesky() -> None:
    points = np.array([-1.7, -0.4, 0.2, 1.1])
    gram = np.exp(-((points[:, None] - points[None, :]) ** 2))
    config = GECPConfig(tol=1e-12, max_rank=4, pivot="grid")
    gecp_result = gecp_matrix(gram, config=config)
    cholesky = pivoted_cholesky(gram, config=config)
    assert gecp_result.row_indices == cholesky.indices
    assert gecp_result.column_indices == cholesky.indices
    assert gecp_result.pivots == pytest.approx(cholesky.pivots)
    assert cholesky.factor @ cholesky.factor.T == pytest.approx(gram, abs=1e-12)


def test_cholesky_rejects_nonsymmetric_input() -> None:
    with pytest.raises(ValueError):
        pivoted_cholesky(np.array([[1.0, 2.0], [0.0, 1.0]]))


def test_psd_ties_use_the_same_diagonal_preferring_rule() -> None:
    gram = np.ones((3, 3))
    config = GECPConfig(tol=1e-14, max_rank=3, pivot="grid")
    gecp_result = gecp_matrix(gram, config=config)
    cholesky = pivoted_cholesky(gram, config=config)
    assert gecp_result.row_indices == [0]
    assert gecp_result.column_indices == [0]
    assert cholesky.indices == [0]
    assert gecp_result.pivots == pytest.approx(cholesky.pivots)
    assert np.linalg.eigvalsh(cholesky.residual_matrix).min() >= -1e-14
