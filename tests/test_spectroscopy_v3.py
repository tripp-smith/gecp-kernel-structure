import numpy as np
import pytest

from experiments import audit_spectroscopy as audit
from experiments import spectroscopy as science
from experiments import spectroscopy_campaign as campaign


def test_v3_profiles_are_separate_and_independently_match_jobs():
    old = campaign.controls()
    new = campaign.controls("v3")
    assert old["final_seeds"] == [4001, 4002, 4003]
    assert old["final_query_order"] == 64
    assert new["final_seeds"] == [5001, 5002, 5003]
    assert new["final_query_order"] == 80
    assert new["calibration_scope"] == "consumed_development"
    assert new["minimizer_marquardt"] is True
    for calibration in (True, False):
        generated = list(campaign.jobs(calibration, "v3"))
        checked = dict(audit.jobs(calibration, "v3"))
        assert {campaign.stem(job): job for job in generated} == checked
        assert all(job["profile"] == "v3" for job in generated)
        assert all(
            job["query_order"] == (24 if calibration else 80) for job in generated
        )
    assert "profile" not in next(campaign.jobs(False))
    assert campaign.controls() == old


def test_v3_query_nodes_consumed_development_seed_only():
    _, query, covariance = science.observations(science.CASES[0], 111, query_order=80)
    expected = (np.polynomial.legendre.leggauss(80)[0] + 1) * 5
    np.testing.assert_array_equal(query, expected)
    assert covariance.shape == (80, 80)
    assert not np.any(np.isclose(query[:, None], np.linspace(0, 10, 65), atol=1e-12))


def test_v3_command_and_native_audit_reject_wrong_flags():
    job = next(campaign.jobs(False, "v3"))
    cmd = campaign.child_command("python", job, "result.json")
    assert cmd[cmd.index("--profile") + 1] == "v3"
    d = {
        "provider": "TRIQS/MaxEnt",
        "version": "4.0.1",
        "mass_operation": "native_mass_constraint",
        "mass_parameterization": {
            "curvature": "projected",
            "d_dv": False,
            "dA_projection": 2,
            "interior_fraction": 1e-12,
            "chi2_factor": 1.0,
            "floor_policy": "reject_native_entropy_floor_or_overflow",
            "min_observed_H": 1e-15,
            "min_observed_mixture_floor": 1e-15,
        },
        "native_minimizer": {
            "class": "LevenbergMinimizer",
            "maxiter": 1000,
            "marquardt": True,
            "J_squared": False,
            "mu0": 1e-8,
        },
    }
    audit.audit_v3_diagnostics(d)
    d["native_minimizer"]["marquardt"] = False
    with pytest.raises(ValueError, match="native controls"):
        audit.audit_v3_diagnostics(d)
    d["native_minimizer"]["marquardt"] = True
    d["mass_operation"] = "explicit_rescale"
    with pytest.raises(ValueError, match="native controls"):
        audit.audit_v3_diagnostics(d)


def test_v3_rejects_a_v2_freeze_before_reading_numerical_payload(tmp_path):
    import json

    manifest, freeze = tmp_path / "campaign.json", tmp_path / "freeze.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.spectroscopy-campaign.v3",
                "kind": "campaign",
                "status": "passed",
            }
        )
    )
    freeze.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.spectroscopy-freeze.v2",
                "status": "passed",
            }
        )
    )
    with pytest.raises(ValueError, match="freeze status/schema"):
        audit.audit_campaign(manifest, freeze, profile="v3")
    with pytest.raises(ValueError, match="campaign schema/status"):
        audit.audit_campaign(manifest, freeze)


def test_v3_independent_controls_reject_post_fit_rescale_configuration():
    config = audit.solver_config("flat", "v3")
    assert config["mass_mode"] == "constrained"
    assert config["minimizer_marquardt"] is True
    legacy = audit.solver_config("flat")
    assert "mass_mode" not in legacy
    assert config != legacy
    # This identity check is made before accepting a repeat as an identical run.
    row = next(campaign.jobs(False, "v3"))
    result = {
        "schema": "kernelgecp.spectroscopy-case.v2",
        "status": "passed",
        "gates": {"positive": True},
    }
    with pytest.raises(ValueError, match="child gates/status"):
        campaign.validate_child(row, result, {}, None)
