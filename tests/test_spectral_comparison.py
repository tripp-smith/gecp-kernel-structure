"""Fairness, returned-model metrics and deterministic small comparison checks."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pytest

from kernelgecp.comparison import (
    METHODS,
    ComparisonConfig,
    _fixtures,
    assess_proposal,
    run_spectral_comparison,
)
from kernelgecp.kernels import FermionicKernel
from kernelgecp.spectral_common import AtomProposal
from kernelgecp.types import NumericalBreakdownError, NumericalMetadata


@pytest.mark.parametrize(
    "field,value",
    [
        ("cutoff", float("nan")),
        ("noise_sigma", float("inf")),
        ("native_tolerance", 0),
        ("training_count", True),
        ("matsubara_modes", 0),
        ("nls_starts", 0),
        ("atom_budgets", (2, 1)),
        ("atom_budgets", (1, True)),
        ("tolerances", (1e-8, 1e-8)),
        ("tolerances", (float("nan"),)),
    ],
)
def test_comparison_config_validation(field: str, value: object) -> None:
    with pytest.raises(ValueError):
        replace(ComparisonConfig(), **{field: value})


def test_validation_values_cannot_change_returned_fit() -> None:
    t = np.linspace(0, 1, 31)
    held_t = np.linspace(0.01, 0.99, 18)
    omega = np.array([-2.2, 0.8])
    values = FermionicKernel(8)(t[:, None], omega) @ np.array([0.4, 0.6])
    proposal = AtomProposal(omega, NumericalMetadata(), {})
    first = assess_proposal(
        proposal,
        method="test",
        budget=2,
        training_t=t,
        observations=values,
        validation_t=held_t,
        validation_values=np.zeros(18),
        cutoff=8,
        tolerances=(1e-8,),
    )
    second = assess_proposal(
        proposal,
        method="test",
        budget=2,
        training_t=t,
        observations=values,
        validation_t=held_t,
        validation_values=np.ones(18),
        cutoff=8,
        tolerances=(1e-8,),
    )
    assert first["weights"] == second["weights"]
    assert first["frequencies"] == second["frequencies"]
    assert first["training_error"] == second["training_error"]
    assert first["held_out_error"] != second["held_out_error"]
    assert not first["certified"]


def test_small_comparison_same_inputs_metrics_and_repeat_bytes() -> None:
    config = ComparisonConfig(
        tolerances=(1e-4,),
        atom_budgets=(1, 2),
        training_count=33,
        validation_count=32,
        candidate_count=101,
        density_count=101,
        matsubara_modes=3,
        grid_order=4,
        precision_bits=53,
        nls_max_evaluations=8,
        nls_starts=1,
        vf_max_iterations=3,
    )
    first = run_spectral_comparison(config, git_commit="small-test")
    second = run_spectral_comparison(config, git_commit="small-test")
    assert json.dumps(first, sort_keys=True, allow_nan=False) == json.dumps(
        second, sort_keys=True, allow_nan=False
    )
    assert not np.intersect1d(first["training_times"], first["validation_times"]).size
    assert first["methods"] == list(METHODS)
    for fixture in first["fixtures"]:
        assert len(fixture["results"]) == len(METHODS) * len(config.atom_budgets)
        for result in fixture["results"]:
            assert result["training_input_hash"] == fixture["training_input_hash"]
            assert result["atom_count"] <= result["requested_atom_budget"]
            assert all(abs(w) <= config.cutoff for w in result["frequencies"])
            assert len(result["weights"]) == result["atom_count"]
            assert result["passes"][0]["passed"] == (result["held_out_error"] <= 1e-4)
            assert not result["certified"]
        assert fixture["transform"]["oracle_scope"].startswith("evaluation_only")
    assert first["known_library_control"]["candidate_strategy"].endswith("not_blind")


def test_nonphysical_proposal_rejected() -> None:
    with pytest.raises(ValueError, match="physical"):
        assess_proposal(
            AtomProposal(np.array([9.0]), NumericalMetadata(), {}),
            method="test",
            budget=1,
            training_t=np.array([0.0, 1.0]),
            observations=np.zeros(2),
            validation_t=np.array([0.5]),
            validation_values=np.zeros(1),
            cutoff=8,
            tolerances=(1e-8,),
        )


def test_blind_primary_truth_is_not_inserted_in_fixed_dictionary() -> None:
    config = ComparisonConfig()
    dictionary = np.linspace(-config.cutoff, config.cutoff, config.candidate_count)
    for name, frequencies, _, _ in _fixtures(config)[:2]:
        assert name in ("clean_off_grid", "noisy_near_fermi")
        distances = np.min(np.abs(frequencies[:, None] - dictionary), axis=1)
        assert np.all(distances > 1e-4)


def test_algorithm_breakdown_is_retained_instead_of_dropped(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed(*args: object, **kwargs: object) -> None:
        raise NumericalBreakdownError("intentional Prony breakdown")

    monkeypatch.setattr("kernelgecp.comparison.prony_atoms", failed)
    result = run_spectral_comparison(
        ComparisonConfig(
            tolerances=(1e-4,),
            atom_budgets=(1,),
            training_count=33,
            validation_count=32,
            candidate_count=101,
            density_count=101,
            matsubara_modes=3,
            grid_order=4,
            precision_bits=53,
            nls_max_evaluations=3,
            nls_starts=1,
            vf_max_iterations=1,
        )
    )
    for fixture in result["fixtures"]:
        prony = next(r for r in fixture["results"] if r["method"] == "prony")
        assert prony["atom_count"] == 0
        assert not prony["native_metadata"]["converged"]
        assert prony["native_metadata"]["stop_reason"] == "numerical_failure"
        assert "intentional" in prony["native_diagnostics"]["failure"]
        assert not prony["passes"][0]["passed"]
    json.dumps(result, allow_nan=False)
