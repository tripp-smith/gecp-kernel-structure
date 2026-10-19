"""Precision-independent numerical evidence is explicitly not certification."""

import mpmath as mp
import pytest

from experiments import companion_envelope
from experiments.companion_envelope import WEIGHTS, audit_case, experiment


def test_original_kernel_companion_stress_evidence() -> None:
    report = audit_case(2, 16)
    assert report["precision_agreement_45_significant_digits"]
    assert report["stop_reason"] == "precision_agreement"
    for run in report["runs"]:
        assert len(run["selected_rows"]) == len(run["selected_columns"]) == 2
        assert len(set(run["sample_columns"])) == 18
        assert (
            run["selected_columns"][0] >= 20
        )  # The first pivot is outside the local band.
        assert run["weights"] == list(WEIGHTS)
        assert run["sampled_envelope_holds"]
        assert run["certification"] == "Not interval-certified or a theorem"
        with mp.workdps(80):
            assert all(0 < mp.mpf(v) < 1 for v in run["determinant_ratios"])
            assert 0 < mp.mpf(run["relative_pivot_gap"]) < 1


def test_companion_stress_rejects_invalid_configuration() -> None:
    for args in ((-1, 16, 320), (2, 0, 320), (2, 16, 53)):
        with pytest.raises(ValueError):
            experiment(*args)


def test_companion_stress_precision_disagreement_is_uncertified(monkeypatch) -> None:
    def fake_experiment(depth, cutoff, digits):
        return {
            "selected_rows": [digits],
            "selected_columns": [],
            "sample_columns": [],
            "determinant_ratios": ["1e-99"],
            "sampled_envelope_holds": True,
        }

    monkeypatch.setattr(companion_envelope, "experiment", fake_experiment)
    report = audit_case(2, 16)
    assert not report["precision_agreement_45_significant_digits"]
    assert report["stop_reason"] == "uncertified_precision_disagreement"


def test_companion_stress_zero_determinant_is_rejected(monkeypatch) -> None:
    monkeypatch.setattr(mp, "det", lambda matrix: mp.mpf(0))
    with pytest.raises(ArithmeticError, match="Uncertified zero determinant"):
        experiment(0, 1, 100)


def test_companion_stress_does_not_round_away_envelope_violation(monkeypatch) -> None:
    monkeypatch.setattr(mp, "det", lambda matrix: 2**18 * (1 + mp.mpf("1e-60")))
    report = audit_case(0, 1, digits=100)
    assert report["precision_agreement_45_significant_digits"]
    assert report["stop_reason"] == "sampled_envelope_violation"
    assert not report["runs"][0]["sampled_envelope_holds"]
