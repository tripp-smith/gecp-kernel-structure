"""Campaign contract regressions; no optional provider substitutes are used."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from experiments import interacting_dmft as campaign

CONFIG = Path(__file__).resolve().parents[1] / "experiments/interacting_dmft.toml"


def test_canonical_config_reserves_unobserved_queries():
    config = campaign.load_config(CONFIG)
    observed = campaign._fit_indices(config)
    assert config["repeats"] == 3
    assert not np.any(np.isin(config["validation_indices"], observed))
    assert not np.any(np.isin(config["calibration_indices"], observed))
    assert np.all(np.isin(np.arange(-18, 18), observed))
    assert config["validation_tau_rule"] == "gauss_legendre_64"


@pytest.mark.parametrize(
    "old,new",
    [
        ("repeats = 3", "repeats = 1"),
        ("n_iw = 512", "n_iw = true"),
        ("mixing = 0.5", "mixing = 2.0"),
        ("target = 1e-6", "target = nan"),
        ("reference_n_tau = 12289", "reference_n_tau = 12290"),
        ("interactions = [0.0, 1.0, 2.0]", "interactions = [0.0, 0.0]"),
        ("gauss_legendre_64", "reuse_training_rows"),
        ("validation_indices = [-947", "validation_indices = [-401"),
    ],
)
def test_invalid_campaign_contract_fails(tmp_path, old, new):
    source = CONFIG.read_text()
    assert old in source
    changed = tmp_path / "invalid.toml"
    changed.write_text(source.replace(old, new))
    with pytest.raises(ValueError):
        campaign.load_config(changed)


def test_complex_scalar_serialization_and_nonfinite_rejection(tmp_path):
    output = tmp_path / "record.json"
    campaign._write(output, {"value": np.complex128(1 - 2j)})
    assert json.loads(output.read_text()) == {"value": {"real": 1, "imag": -2}}
    with pytest.raises(ValueError):
        campaign._write(output, {"invalid": np.nan})
    assert json.loads(output.read_text()) == {"value": {"real": 1, "imag": -2}}


def test_canonical_cli_cannot_pass_a_failed_determinism_gate(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["campaign", "--output-dir", str(tmp_path)])
    monkeypatch.setattr(
        campaign,
        "run",
        lambda *args, **kwargs: {
            "all_mandatory_passed": False,
            "campaign_checks_passed": True,
        },
    )
    with pytest.raises(SystemExit) as error:
        campaign.main()
    assert error.value.code == 1
