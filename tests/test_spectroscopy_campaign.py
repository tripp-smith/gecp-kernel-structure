"""Independent conventions and negative acceptance controls for the campaign."""

import numpy as np
import pytest
from scipy.integrate import quad

from experiments.spectroscopy import observations, smooth_density, true_cdf, truth


def test_smooth_density_and_cdf_accept_integer_window_edges():
    edges = np.array([-4, -1, 1, 4])
    cdf = true_cdf("smooth_correlated", edges)
    assert cdf[0] == pytest.approx(0, abs=1e-15)
    assert cdf[-1] == pytest.approx(1, abs=1e-15)
    for lower, upper, difference in zip(
        edges[:-1], edges[1:], np.diff(cdf), strict=True
    ):
        mass = quad(lambda w: float(smooth_density(w)), lower, upper, epsabs=1e-13)[0]
        assert difference == pytest.approx(mass, abs=1e-13)


def test_synthetic_forward_sign_and_fermionic_endpoint_rule():
    for case in ("offgrid_atoms_independent", "smooth_correlated"):
        values = truth(case, 10, np.array([0, 10]))
        assert values[0] + values[1] == pytest.approx(-1, abs=2e-13)
        assert np.max(values) < 0


def test_seeds_and_observation_query_isolation_are_explicit():
    data, query, covariance = observations("offgrid_atoms_correlated", 111)
    repeated, repeated_query, _ = observations("offgrid_atoms_correlated", 111)
    other, _, _ = observations("offgrid_atoms_correlated", 112)
    np.testing.assert_array_equal(data.values, repeated.values)
    np.testing.assert_array_equal(query, repeated_query)
    assert not np.array_equal(data.values, other.values)
    assert not np.any(np.isin(query, data.tau))
    assert np.linalg.norm(covariance - np.diag(np.diag(covariance))) > 0
    assert data.modes.effective_rank == len(data.tau)
