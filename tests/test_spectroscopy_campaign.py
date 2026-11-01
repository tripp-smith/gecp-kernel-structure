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
    final, final_query, _ = observations(
        "offgrid_atoms_correlated", 1001, query_order=36
    )
    assert len(final_query) == 36
    assert not np.any(np.isin(query, final_query))
    assert not np.array_equal(final.values, data.values)


def test_new_final_noise_and_query_nodes_are_disjoint_from_consumed_validation():
    from experiments.spectroscopy_campaign import controls, jobs

    frozen = controls()
    assert frozen["final_seeds"] == [3001, 3002, 3003]
    assert frozen["final_query_order"] == 64
    assert frozen["minimizer_mu0"] == 1e-8
    for offset, case in enumerate(frozen["cases"]):
        old, old_query, _ = observations(case, 2001 + offset, query_order=48)
        new, new_query, _ = observations(case, 3001 + offset, query_order=64)
        assert not np.array_equal(old.values, new.values)
        assert not np.intersect1d(old_query, new_query).size
        assert not np.intersect1d(new.tau, new_query).size
    assert all(job["minimizer_mu0"] == 1e-8 for job in jobs(True))
    assert all(job["minimizer_mu0"] == 1e-8 for job in jobs(False))
