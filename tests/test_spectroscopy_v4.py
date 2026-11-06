import copy

import pytest

from experiments import audit_cthyb_spectroscopy as native_audit
from experiments import audit_spectroscopy as synthetic_audit
from experiments import cthyb_fixed_bath as fixed
from experiments import cthyb_spectroscopy as native
from experiments import spectroscopy_campaign as campaign


def test_v4_controls_and_reservations_are_versioned():
    p = native.protocol_for("v4")
    assert p["profile"] == "native_shifted_curvature_v4"
    assert p["mass_curvature_floor_relative"] == 1e-10
    assert fixed.seeds_for("v4") == [34100001 + 2 * i for i in range(256)]
    assert native.schema_for("-worker", "v4").endswith(".v4")
    assert "mass_curvature_floor_relative" not in native.protocol_for("v3")
    c = campaign.controls("v4")
    assert c["profile"] == "constrained_mass_v4_synthetic_binary64"
    assert c["mass_curvature_floor_relative"] == 1e-10
    assert c["final_seeds"] == [5001, 5002, 5003]
    assert c["final_query_order"] == 80
    assert "mass_curvature_floor_relative" not in campaign.controls("v3")
    assert (
        synthetic_audit.solver_config("flat", "v4")["mass_curvature_floor_relative"]
        == 1e-10
    )
    for calibration in (False, True):
        jobs = list(campaign.jobs(calibration, "v4"))
        assert {campaign.stem(j): j for j in jobs} == dict(
            synthetic_audit.jobs(calibration, "v4")
        )
        assert all(j["profile"] == "v4" for j in jobs)
        command = campaign.child_command("python", jobs[0], "out.json")
        assert command[command.index("--profile") + 1] == "v4"


def metric():
    return {
        "evaluation_count": 1,
        "certified": False,
        "last_evaluation": {
            "floor_relative": 1e-10,
            "dimension": 2,
            "scale": 2.0,
            "normalized_min_eigenvalue": 0.0,
            "normalized_max_eigenvalue": 1.0,
            "normalized_shift": 2e-10,
            "shift": 4e-10,
            "shifted_normalized_min_eigenvalue": 2e-10,
            "shifted_normalized_max_eigenvalue": 1.0000000002,
            "condition_estimate": 1.0000000002 / 2e-10,
            "certified": False,
        },
    }


@pytest.mark.parametrize(
    "audit",
    [native_audit.audit_curvature_metric, synthetic_audit.audit_curvature_metric],
)
def test_metric_audits_reject_missing_wrong_and_legacy(audit):
    config = {"mass_curvature_floor_relative": 1e-10}
    parameter = {"curvature_metric": metric()}
    audit(parameter, config)
    for mutate in (
        lambda p: p.clear(),
        lambda p: p["curvature_metric"]["last_evaluation"].update(floor_relative=1e-9),
        lambda p: p["curvature_metric"].update(evaluation_count=0),
        lambda p: p["curvature_metric"]["last_evaluation"].update(shift=0),
    ):
        bad = copy.deepcopy(parameter)
        mutate(bad)
        with pytest.raises((ValueError, AssertionError)):
            audit(bad, config)
    with pytest.raises((ValueError, AssertionError)):
        audit(parameter, {})
