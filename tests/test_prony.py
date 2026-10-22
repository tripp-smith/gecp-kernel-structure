"""Prony recovery independent of a supplied generating order or dictionary."""

import json

import numpy as np
import pytest

from kernelgecp.prony import prony_atoms


def test_prony_data_rank_recovers_both_signs_below_budget() -> None:
    t = np.linspace(0, 1, 257)
    frequencies = np.array([-5.4, -2.2, 0.8, 4.5])
    # Construct exponentials independently of the fermionic kernel evaluator.
    values = np.exp(-t[:, None] * frequencies) @ np.array([0.01, 0.1, 0.4, 0.2])
    result = prony_atoms(t, values, cutoff=8, max_atoms=8, tolerance=1e-10)
    assert result.metadata.converged
    assert result.diagnostics["detected_rank"] == 4
    assert np.allclose(result.frequencies, frequencies, atol=1e-6, rtol=0)
    assert result.diagnostics["hankel_stride"] > 1
    assert result.diagnostics["hankel_shape"][1] == 9
    assert result.diagnostics["hankel_shape"][0] < len(t)
    assert not result.diagnostics["native_fit_weights_are_spectral_weights"]
    assert result.metadata.refinement.status == "not_requested"
    assert result.metadata.certification_status == "not_certified"
    json.dumps(result.diagnostics, allow_nan=False)
    result.metadata.canonical_json()


def test_prony_zero_short_and_budget_exhaustion() -> None:
    t = np.linspace(0, 1, 33)
    zero = prony_atoms(t, np.zeros(33), cutoff=4, max_atoms=8)
    assert zero.metadata.converged and not zero.frequencies.size
    limited = prony_atoms(t, np.exp(-t), cutoff=4, max_atoms=0)
    assert not limited.metadata.converged
    assert limited.metadata.stop_reason == "atom_budget_exhausted"
    short = prony_atoms([0, 1], [1, 0.5], cutoff=4, max_atoms=4)
    assert not short.metadata.converged
    assert short.metadata.stop_reason == "insufficient_samples"


@pytest.mark.parametrize("case", ["outside", "complex"])
def test_nonphysical_prony_nodes_rejected_without_clipping(case: str) -> None:
    t = np.linspace(0, 1, 65)
    values = np.exp(-10 * t) if case == "outside" else np.cos(5 * t)
    result = prony_atoms(t, values, cutoff=4, max_atoms=4, tolerance=1e-12)
    assert not result.frequencies.size
    assert result.diagnostics["rejected_root_count"] > 0
    assert not result.metadata.converged
    assert all(not item["accepted"] for item in result.diagnostics["raw_roots"])


def test_noisy_clustered_data_retains_diagnostics_and_repeats() -> None:
    t = np.linspace(0, 1, 129)
    frequencies = np.array([-3.8, -0.45, 0.35, 3.2])
    clean = np.exp(-t[:, None] * frequencies) @ np.array([0.02, 0.3, 0.4, 0.1])
    observed = clean + np.random.default_rng(17).normal(0, 1e-8, len(t))
    first = prony_atoms(t, observed, cutoff=8, max_atoms=8, tolerance=1e-10)
    second = prony_atoms(t, observed, cutoff=8, max_atoms=8, tolerance=1e-10)
    assert first.diagnostics["rank_threshold_normalized"] > 0
    assert first.frequencies.size <= 8
    assert first.metadata.converged == (first.diagnostics["native_fit_error"] <= 1e-10)
    assert first.metadata.canonical_json() == second.metadata.canonical_json()
    assert np.array_equal(first.frequencies, second.frequencies)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_atoms": True},
        {"max_atoms": -1},
        {"cutoff": np.nan},
        {"tolerance": 0},
        {"rank_rtol": -1},
        {"pole_imag_tolerance": np.inf},
    ],
)
def test_prony_rejects_bad_controls(kwargs: dict) -> None:
    controls = {"max_atoms": 2, "cutoff": 4, **kwargs}
    with pytest.raises(ValueError):
        prony_atoms(np.linspace(0, 1, 9), np.ones(9), **controls)


@pytest.mark.parametrize(
    "t,values",
    [
        ([], []),
        ([0, 0.4, 1], [1, 2, 3]),
        ([0, 0], [1, 2]),
        ([0, 1], [1]),
        ([0, 1], [1, np.inf]),
        ([-1, 0, 1], [1, 2, 3]),
    ],
)
def test_prony_rejects_invalid_time_data(t: object, values: object) -> None:
    with pytest.raises(ValueError):
        prony_atoms(t, values, cutoff=4, max_atoms=2)
