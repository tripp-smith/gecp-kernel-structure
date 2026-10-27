"""Independent covariance/acceptance arithmetic and fail-closed archive gates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from experiments.audit_cthyb_dmft import _agreement, _statistics, audit


def test_independent_scalar_covariance_closed_form() -> None:
    rows = np.array([[1.0, 3.0], [3.0, 1.0], [5.0, 5.0]])
    result = _statistics(rows)
    np.testing.assert_array_equal(result["mean"], [3.0, 3.0])
    np.testing.assert_array_equal(result["covariance"], [[4.0, 2.0], [2.0, 4.0]])
    np.testing.assert_allclose(
        result["covariance_of_mean"], np.array([[4.0, 2.0], [2.0, 4.0]]) / 3
    )
    np.testing.assert_allclose(result["standard_error"], np.sqrt([4 / 3, 4 / 3]))
    assert result["replicas"] == 3


@pytest.mark.parametrize("rows", [np.ones((1, 3)), np.array([[1.0], [float("nan")]])])
def test_independent_statistics_refuse_invalid_chain_populations(
    rows: np.ndarray,
) -> None:
    with pytest.raises(ValueError):
        _statistics(rows)


def test_agreement_uses_three_combined_standard_errors_per_observable() -> None:
    first = {"mean": np.array([0.0, 0.0]), "standard_error": np.array([0.0003, 0.0001])}
    second = {
        "mean": np.array([0.002, 0.002]),
        "standard_error": np.array([0.0003, 0.0001]),
    }
    data = {"agreement_absolute": 0.001, "agreement_sigma": 3}
    result = _agreement(first, second, data)
    np.testing.assert_allclose(
        result["limits"], 0.001 + 3 * np.sqrt([2 * 0.0003**2, 2 * 0.0001**2])
    )
    assert result["passed"] is False
    second["mean"][1] = 0.0012
    assert _agreement(first, second, data)["passed"] is True


def test_audit_refuses_tuning_pilot_as_canonical_evidence(tmp_path: Path) -> None:
    identity = {
        "schema": "kernelgecp.cthyb-campaign.v1",
        "config": {
            "canonical": False,
            "cycle_budgets": [1000, 4000, 16000],
        },
    }
    (tmp_path / "manifest.json").write_text(json.dumps(identity))
    with pytest.raises(ValueError, match="tuning pilot"):
        audit(tmp_path)


def test_audit_does_not_depend_on_campaign_estimator_or_acceptance_booleans() -> None:
    # Static import/use guard supplements the independent numerical examples.
    # This keeps accidental refactoring from turning the auditor into an alias.
    import inspect

    import experiments.audit_cthyb_dmft as auditor

    source = inspect.getsource(auditor)
    assert "from kernelgecp.integrations.cthyb_dmft import" not in source
    assert "estimate_replicas(" not in source
    assert 'saved["checks"]' not in source
    assert 'saved["status"]' not in source
