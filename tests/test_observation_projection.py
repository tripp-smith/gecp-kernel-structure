"""Independent observation-space oracle for the cached retained SVD."""

import numpy as np
import pytest

from kernelgecp.integrations.maxent import _SpectralLeastSquaresPlan


def test_complex_projection_matches_independent_qr_and_is_owned():
    rng = np.random.default_rng(541)
    design = rng.normal(size=(13, 4)) + 1j * rng.normal(size=(13, 4))
    plan = _SpectralLeastSquaresPlan(design)
    projection = plan.observation_projection()
    q, _ = np.linalg.qr(design)
    np.testing.assert_allclose(projection, q @ q.conj().T, atol=1e-15)
    np.testing.assert_allclose(projection.conj().T, projection, atol=1e-16)
    np.testing.assert_allclose(projection @ projection, projection, atol=1e-15)
    values = rng.normal(size=(13, 2, 3)) + 1j * rng.normal(size=(13, 2, 3))
    projected = np.einsum("ij,jab->iab", projection, values)
    predicted = np.einsum("ij,jab->iab", design, plan.fit(values))
    np.testing.assert_allclose(projected, predicted, atol=3e-15)
    with pytest.raises(ValueError):
        projection.setflags(write=True)
    design[:] = 0
    np.testing.assert_allclose(plan.observation_projection(), projection)


def test_zero_rank_projection_and_full_matrix_resource_limit():
    plan = _SpectralLeastSquaresPlan(np.zeros((5, 2), complex))
    np.testing.assert_array_equal(plan.observation_projection(), np.zeros((5, 5)))
    bounded = _SpectralLeastSquaresPlan(np.ones((2049, 1)))
    with pytest.raises(ValueError, match="2048"):
        bounded.observation_projection()
