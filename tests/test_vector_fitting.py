"""Independent rational fixtures for genuine, support-blind pole relocation."""

import json

import numpy as np
import pytest

from kernelgecp.vector_fitting import _real_basis, vector_fitting_atoms


def _points(count: int = 40) -> np.ndarray:
    return 1j * np.pi * (2 * np.arange(count) + 1)


def _rational(
    points: np.ndarray, poles: np.ndarray, residues: np.ndarray
) -> np.ndarray:
    return (1 / (points[:, None] - poles[None, :])) @ residues


def test_vector_fitting_recovers_positive_and_negative_offgrid_poles() -> None:
    points = _points()
    poles = np.array([-5.4, -2.2, 0.8, 4.5])
    residues = np.array([0.3, 0.5, 0.6, -0.2])
    values = _rational(points, poles, residues)
    result = vector_fitting_atoms(
        points, values, cutoff=8, max_atoms=4, tolerance=1e-12
    )
    assert result.metadata.converged
    assert result.metadata.stop_reason == "native_tolerance"
    assert np.allclose(result.frequencies, poles, atol=2e-9, rtol=0)
    assert np.allclose(np.array(result.diagnostics["raw_residues"])[:, 0], residues)
    # A fixed-pole least-squares fit cannot pass this fixture.  The history
    # proves actual eig(A-bc) relocation from a truth-independent pole set.
    assert result.diagnostics["history"][0]["native_relative_error"] > 1e-3
    assert result.diagnostics["relocations"] >= 1
    assert result.diagnostics["history"][0]["maximum_sorted_pole_movement"] > 1
    assert result.diagnostics["native_relative_error"] < 1e-12
    assert result.metadata.configuration["polynomial_terms"].startswith("none")
    assert result.metadata.configuration["stability_flipping"] is False
    assert result.metadata.refinement.optimizer_success is None
    assert result.metadata.certification_status == "not_certified"
    assert result.metadata.conditioning_history


def test_zero_iterations_is_a_visible_fixed_pole_budget_failure() -> None:
    points = _points()
    values = 1 / (points - 2.3)
    result = vector_fitting_atoms(
        points, values, cutoff=4, max_atoms=1, max_iterations=0
    )
    assert not result.metadata.converged
    assert result.metadata.stop_reason == "iteration_budget"
    assert result.diagnostics["relocations"] == 0
    assert result.frequencies.tolist() == [0]
    assert result.diagnostics["native_relative_error"] > 1e-3


def test_one_iteration_relocates_single_offgrid_pole() -> None:
    points = _points()
    result = vector_fitting_atoms(
        points,
        -0.3 / (points + 1.75),
        cutoff=4,
        max_atoms=1,
        max_iterations=1,
        tolerance=1e-13,
    )
    assert result.metadata.converged
    assert np.allclose(result.frequencies, [-1.75], atol=1e-12, rtol=0)
    assert result.diagnostics["relocations"] == 1


def test_native_convergence_does_not_mean_nonphysical_poles_are_atoms() -> None:
    points = _points()
    poles = np.array([1.2 - 0.7j, 1.2 + 0.7j])
    residues = np.array([0.4 + 0.2j, 0.4 - 0.2j])
    result = vector_fitting_atoms(
        points,
        _rational(points, poles, residues),
        cutoff=4,
        max_atoms=2,
        tolerance=1e-12,
    )
    assert result.metadata.converged
    assert result.frequencies.size == 0
    assert not result.diagnostics["physical_native_converged"]
    assert result.diagnostics["rejected_pole_count"] == 2
    assert all(
        item["reason"] == "nonreal" for item in result.diagnostics["rejected_poles"]
    )
    assert "native_rational" in result.metadata.convergence_scope


def test_out_of_band_pole_is_rejected_without_clipping() -> None:
    points = _points()
    result = vector_fitting_atoms(
        points, 1 / (points - 10), cutoff=4, max_atoms=1, tolerance=1e-12
    )
    assert result.metadata.converged
    assert result.frequencies.size == 0
    assert result.diagnostics["rejected_pole_count"] == 1
    assert result.diagnostics["rejected_poles"][0]["reason"] == "outside_cutoff"
    assert result.diagnostics["raw_poles"][0][0] == pytest.approx(10)


def test_real_pair_basis_and_state_relocation_match_complex_formula() -> None:
    """A hand-constructed complex-pair check, independent of the fit loop."""
    points = _points(5)
    poles = np.array([-2 + 0j, 0.5 - 1.2j, 0.5 + 1.2j])
    basis, state, control, residue_map = _real_basis(points, poles)
    coefficients = np.array([0.2, -0.3, 0.7])
    complex_residues = residue_map @ coefficients
    assert np.allclose(
        basis @ coefficients, _rational(points, poles, complex_residues), atol=1e-15
    )
    expected = np.linalg.eigvals(
        np.diag(poles) - np.outer(np.ones(3), complex_residues)
    )
    actual = np.linalg.eigvals(state - np.outer(control, coefficients))
    assert np.allclose(np.sort_complex(actual), np.sort_complex(expected), atol=1e-13)


@pytest.mark.parametrize("max_atoms", [0, 4])
def test_zero_data_returns_zero_model_without_spurious_atoms(max_atoms: int) -> None:
    points = _points(6)
    result = vector_fitting_atoms(
        points, np.zeros(points.size), cutoff=4, max_atoms=max_atoms
    )
    assert result.metadata.converged
    assert result.metadata.stop_reason == "zero_model_tolerance"
    assert result.frequencies.size == 0
    assert result.diagnostics["relocations"] == 0
    assert result.diagnostics["native_relative_error"] == 0
    assert result.metadata.conditioning_status == "not_applicable"


def test_nonzero_data_with_no_atoms_is_not_converged() -> None:
    points = _points(4)
    result = vector_fitting_atoms(points, 1 / (points - 2), cutoff=4, max_atoms=0)
    assert result.frequencies.size == 0
    assert not result.metadata.converged
    assert result.metadata.stop_reason == "atom_budget"


def test_one_sample_reports_sample_limited_requested_budget() -> None:
    points = _points(1)
    result = vector_fitting_atoms(
        points, 0.7 / (points + 1.4), cutoff=4, max_atoms=9, tolerance=1e-12
    )
    assert result.metadata.converged
    assert np.allclose(result.frequencies, [-1.4], atol=1e-12, rtol=0)
    assert result.metadata.configuration["effective_atom_budget"] == 1
    assert result.diagnostics["sample_limited"]


def test_numerical_failure_is_visible_and_keeps_no_fake_convergence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def failed(*args: object, **kwargs: object) -> None:
        raise np.linalg.LinAlgError("intentional least-squares failure")

    monkeypatch.setattr(np.linalg, "lstsq", failed)
    points = _points(4)
    result = vector_fitting_atoms(points, 1 / (points - 2), cutoff=4, max_atoms=1)
    assert not result.metadata.converged
    assert result.metadata.stop_reason == "numerical_breakdown"
    assert "intentional" in result.diagnostics["failure_message"]
    assert result.frequencies.size == 0
    json.dumps(result.diagnostics, allow_nan=False)
    result.metadata.canonical_json()


def test_determinism_and_json_safe_metadata() -> None:
    points = _points(12)
    values = 1 / (points + 1.7) - 0.3 / (points - 2.1)
    first = vector_fitting_atoms(points, values, cutoff=4, max_atoms=2)
    second = vector_fitting_atoms(points, values, cutoff=4, max_atoms=2)
    assert np.array_equal(first.frequencies, second.frequencies)
    assert first.metadata.canonical_json() == second.metadata.canonical_json()
    assert json.dumps(first.diagnostics, sort_keys=True, allow_nan=False) == json.dumps(
        second.diagnostics, sort_keys=True, allow_nan=False
    )


@pytest.mark.parametrize(
    "arguments",
    [
        {"max_atoms": -1},
        {"max_atoms": 2.5},
        {"max_atoms": True},
        {"max_iterations": -1},
        {"max_iterations": 2.5},
        {"max_iterations": False},
        {"cutoff": 0},
        {"cutoff": np.inf},
        {"tolerance": 0},
        {"tolerance": np.nan},
        {"pole_imag_tolerance": -1},
        {"pole_imag_tolerance": np.inf},
    ],
)
def test_invalid_configuration(arguments: dict) -> None:
    config = {"cutoff": 4, "max_atoms": 2, **arguments}
    with pytest.raises(ValueError):
        vector_fitting_atoms(_points(4), np.ones(4), **config)


@pytest.mark.parametrize(
    ("points", "values"),
    [
        ([], []),
        ([[1j, 2j]], [[1, 2]]),
        ([1j, 2j], [1]),
        ([1j, 1j], [1, 2]),
        ([1.0, 2j], [1, 2]),
        ([1j, np.nan], [1, 2]),
        ([1j, 2j], [1, np.inf]),
    ],
)
def test_invalid_samples(points: object, values: object) -> None:
    with pytest.raises(ValueError):
        vector_fitting_atoms(points, values, cutoff=4, max_atoms=2)
