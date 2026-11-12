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
    assert frozen["final_seeds"] == [4001, 4002, 4003]
    assert frozen["final_query_order"] == 64
    assert frozen["minimizer_mu0"] == 1e-8
    assert frozen["rcond"] == 1e-10 and frozen["precision_bits"] == 53
    assert frozen["profile"] == "corrective_projection_v2_synthetic_binary64"
    # Check reserved identities without generating the untouched final noise.
    assert not set(frozen["final_seeds"]) & {
        101,
        102,
        103,
        1001,
        1002,
        1003,
        2001,
        2002,
        2003,
        3001,
        3002,
        3003,
    }
    fixture, query, _ = observations("offgrid_atoms_correlated", 111, query_order=64)
    assert not np.intersect1d(fixture.tau, query).size
    assert all(job["minimizer_mu0"] == 1e-8 for job in jobs(True))
    assert all(job["minimizer_mu0"] == 1e-8 for job in jobs(False))


def test_case_serializes_single_source_design_with_binary64_controls(
    tmp_path, monkeypatch
):
    from experiments import spectroscopy as science

    calls = []

    class Basis:
        beta = 10

        def design_tau(self, tau):
            calls.append(1)
            return np.polynomial.legendre.legvander(2 * tau / 10 - 1, 3)

    class Spectrum:
        omega = np.linspace(-4, 4, 9)
        density = np.ones(9) / 8
        mass = raw_mass = 1.0
        alphas = np.array([1.0])
        chi2 = np.array([1.0])
        trial_spectra = density[None, :]
        selected_index = 0
        converged = (True,)
        iterations = (1,)
        effective_rank = 4
        normalized_training_whitened_rms = 0.1
        diagnostics = {}

        def evaluate_tau(self, beta, tau):
            return science.truth("offgrid_atoms_correlated", beta, tau)

    monkeypatch.setattr(science, "basis", lambda *args: Basis())
    monkeypatch.setattr(science, "provenance", lambda: {})
    monkeypatch.setattr(science, "continue_tau_maxent", lambda *args: Spectrum())
    output = tmp_path / "case.json"
    result = science.run_case(
        "offgrid_atoms_correlated",
        "pydlr",
        101,
        1e-8,
        {},
        output,
        pilot=True,
        rcond=1e-10,
        precision_bits=53,
    )
    assert calls == [1]
    assert result["schema"] == "kernelgecp.spectroscopy-case.v2"
    assert result["rcond"] == result["projection_svd"]["rcond"] == 1e-10
    assert result["precision_bits"] == 53
    with np.load(output.with_suffix(".npz")) as archive:
        expected = np.polynomial.legendre.legvander(2 * archive["tau"] / 10 - 1, 3)
        np.testing.assert_array_equal(archive["observation_design"], expected)
    with pytest.raises(ValueError, match="precision_bits=53"):
        science.run_case(
            "offgrid_atoms_correlated",
            "pydlr",
            101,
            1e-8,
            {},
            tmp_path / "bad.json",
            pilot=True,
            precision_bits=256,
        )

    class DisagreeingFit:
        def evaluate_tau(self, tau):
            return np.full(len(tau), 100.0)

        def save(self, path):
            path.write_bytes(b"fixture")

    class GecpBasis(Basis):
        def fit_tau(self, values, tau, *, rcond):
            assert rcond == 1e-10
            return DisagreeingFit()

    monkeypatch.setattr(science, "basis", lambda *args: GecpBasis())
    failed = science.run_case(
        "offgrid_atoms_correlated",
        "kernelgecp",
        101,
        1e-8,
        {},
        tmp_path / "gecp.json",
        pilot=True,
    )
    assert failed["status"] == "failed"
    assert failed["gates"]["coefficient_projection_agreement"] is False
