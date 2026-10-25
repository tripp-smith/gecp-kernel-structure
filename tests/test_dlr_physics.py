"""Independent physical oracles, validation isolation and cached solver checks."""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError

import mpmath as mp
import numpy as np
import pytest

from experiments import dlr_physics
from experiments.dlr_physics import (
    BetheSolverPlan,
    frequency_index_meshes,
    matrix_amplitudes,
    matsubara_frequencies,
    matsubara_green_function,
    physical_tau_meshes,
    semicircle_matsubara,
    semicircle_tau,
    solver_index_mesh,
    solver_validation_index_mesh,
)
from experiments.dlr_workloads import common_tau_meshes, make_pole_workloads
from kernelgecp import FermionicBasis, LeastSquaresPlan


@pytest.mark.parametrize("beta,wmax", [(0.25, 2.0), (10.0, 2.0), (100.0, 10.24)])
def test_time_meshes_reuse_training_and_are_actually_pairwise_disjoint(
    beta: float, wmax: float
) -> None:
    training, calibration, final = physical_tau_meshes(beta=beta, wmax=wmax, order=8)
    expected_training, expected_midpoints = common_tau_meshes(
        beta=beta, wmax=wmax, order=8
    )
    np.testing.assert_array_equal(training, expected_training)
    np.testing.assert_array_equal(calibration, expected_midpoints)
    left = training[:-1]
    widths = training[1:] - left
    np.testing.assert_array_equal(
        final,
        np.column_stack((left + widths / 5.0, left + (4.0 / 5.0) * widths)).ravel(),
    )
    assert len(final) == 2 * (len(training) - 1)
    for mesh in (training, calibration, final):
        assert np.all(np.diff(mesh) > 0)
        assert mesh.min() >= 0 and mesh.max() <= beta
        with pytest.raises(ValueError):
            mesh.setflags(write=True)
    for left, right in (
        (training, calibration),
        (training, final),
        (calibration, final),
    ):
        assert not np.intersect1d(left, right).size
    assert np.all(training[:-1] < final[::2])
    assert np.all(final[::2] < calibration)
    assert np.all(calibration < final[1::2])
    assert np.all(final[1::2] < training[1:])


def test_meshes_reject_actual_floating_holdouts_that_collapse(monkeypatch) -> None:
    left = 0.5
    middle = np.nextafter(left, np.inf)
    right = np.nextafter(middle, np.inf)
    monkeypatch.setattr(
        dlr_physics,
        "common_tau_meshes",
        lambda **kwargs: (np.array([left, right]), np.array([middle])),
    )
    with pytest.raises(ValueError, match="disjoint|increasing"):
        physical_tau_meshes(beta=1, wmax=1)


def test_fixed_frequency_grids_isolate_forward_and_solver_holdouts() -> None:
    calibration, final = frequency_index_meshes()
    solver_final = solver_validation_index_mesh()
    solver = solver_index_mesh()
    assert calibration.dtype == final.dtype == solver.dtype == np.int64
    np.testing.assert_array_equal(
        final,
        [-801, -401, -171, -71, -38, -22, -9, -7, 6, 8, 21, 37, 70, 170, 400, 800],
    )
    np.testing.assert_array_equal(
        calibration,
        np.setdiff1d(
            np.unique(np.r_[np.arange(-32, 32), [-256, -128, -64, 63, 127, 255]]),
            final,
        ),
    )
    assert 4 in calibration and -5 in calibration
    assert not np.intersect1d(calibration, final).size
    assert not np.intersect1d(solver, solver_final).size
    np.testing.assert_array_equal(solver_final, final[np.abs(2 * final + 1) >= 37])
    assert set(range(-18, 18)).issubset(solver)
    np.testing.assert_array_equal(
        solver, np.setdiff1d(np.arange(-128, 128), solver_final)
    )
    for mesh in (calibration, final, solver, solver_final):
        np.testing.assert_array_equal(mesh, -mesh[::-1] - 1)
        assert not mesh.flags.writeable


def test_matrix_amplitude_is_complex_hermitian_positive_and_unit_maximum() -> None:
    amplitude = matrix_amplitudes()
    assert amplitude.shape == (2, 2)
    assert amplitude.dtype == np.complex128
    np.testing.assert_array_equal(amplitude, amplitude.conj().T)
    assert np.min(np.linalg.eigvalsh(amplitude)) > 0
    assert np.max(np.abs(amplitude)) == 1.0
    assert amplitude[0, 1].imag != 0
    assert not amplitude.flags.writeable


def test_finite_pole_frequency_truth_matches_independent_expression() -> None:
    beta, wmax = 7.5, 2.25
    indices = np.array([-200, -6, -2, -1, 0, 1, 5, 199])
    for workload in make_pole_workloads(beta=beta, wmax=wmax):
        with mp.workdps(90):
            expected = np.array(
                [
                    complex(
                        mp.fsum(
                            mp.mpf(weight)
                            / (
                                1j * (2 * int(index) + 1) * mp.pi / mp.mpf(beta)
                                - mp.mpf(pole)
                            )
                            for pole, weight in zip(
                                workload.poles, workload.weights, strict=True
                            )
                        )
                    )
                    for index in indices
                ]
            )
        actual = matsubara_green_function(workload, indices, beta=beta, wmax=wmax)
        np.testing.assert_allclose(actual, expected, rtol=5e-16, atol=2e-17)
        np.testing.assert_allclose(actual[:4], actual[-1:-5:-1].conj(), atol=2e-17)


@pytest.mark.parametrize("indices", [[], 0, [1.0], [True], [False, 0], [0, 0], [2**51]])
def test_frequency_truth_rejects_invalid_ordinary_indices(indices: object) -> None:
    with pytest.raises(ValueError):
        matsubara_frequencies(indices, beta=10)


def test_semicircle_root_is_causal_normalized_and_solves_quadratic() -> None:
    beta, hopping = 10.0, 1.0
    indices = np.array([-(10**8), -128, -2, -1, 0, 1, 127, 10**8 - 1])
    omega = matsubara_frequencies(indices, beta=beta)
    values = semicircle_matsubara(indices, beta=beta, hopping=hopping)
    np.testing.assert_allclose(values, values[::-1].conj(), rtol=0, atol=0)
    assert np.all(values.imag * omega < 0)
    np.testing.assert_allclose(
        values, 1 / (1j * omega - hopping**2 * values), rtol=4e-16
    )
    assert (1j * omega[-1] * values[-1]).real == pytest.approx(1.0, abs=1e-15)


def test_semicircle_root_retains_nonzero_subnormal_when_quadratic_scale_overflows() -> (
    None
):
    beta, hopping = 2e-308, 8e307
    actual = semicircle_matsubara([-1, 0], beta=beta, hopping=hopping)
    with mp.workdps(90):
        frequency = mp.pi / mp.mpf(beta)
        expected = float(
            2 / (frequency + mp.sqrt(frequency**2 + 4 * mp.mpf(hopping) ** 2))
        )
    assert expected > 0 and actual[1].imag < 0
    np.testing.assert_allclose(
        actual, [1j * expected, -1j * expected], rtol=1e-14, atol=np.nextafter(0.0, 1.0)
    )


def test_semicircle_quadrature_matches_independent_high_precision_energy_integral() -> (
    None
):
    beta, hopping = 10.0, 1.0
    tau = np.array([0.0, 0.13, 1.7, 5.0, 9.87, 10.0])
    with mp.workdps(70):
        expected = np.array(
            [
                float(
                    mp.quad(
                        lambda energy, sample_time=time: (
                            -mp.sqrt(4 - energy**2)
                            / (2 * mp.pi)
                            * mp.exp(-mp.mpf(sample_time) * energy)
                            / (1 + mp.exp(-mp.mpf(beta) * energy))
                        ),
                        [-2, 0, 2],
                    )
                )
                for time in tau
            ]
        )
    actual = semicircle_tau(tau, beta=beta, hopping=hopping)
    np.testing.assert_allclose(actual, expected, rtol=8e-15, atol=3e-16)
    assert actual[0] + actual[-1] == pytest.approx(-1.0, abs=6e-16)
    np.testing.assert_allclose(
        actual, semicircle_tau(beta - tau, beta=beta, hopping=hopping)
    )


class _TinyBasis:
    beta = 2.0
    rank = 4

    def __init__(self) -> None:
        self.calls = 0
        self.matrix = np.array(
            [[1, 0.1j, 0, 0], [0, 1, 0.2, 0], [0, 0, 1, -0.2j], [0.1, 0, 0, 1]],
            dtype=np.complex128,
        )

    def design_matsubara(self, indices: object) -> np.ndarray:
        self.calls += 1
        return self.matrix

    def design_tau(self, tau: object) -> np.ndarray:
        raise AssertionError("the scalar frequency solver does not evaluate time truth")


def test_solver_caches_one_design_and_svd_and_matches_direct_iteration(
    monkeypatch,
) -> None:
    basis = _TinyBasis()
    indices = np.array([-2, -1, 0, 1])
    svds = 0
    original_svd = np.linalg.svd

    def counted_svd(*args, **kwargs):
        nonlocal svds
        svds += 1
        return original_svd(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "svd", counted_svd)
    plan = BetheSolverPlan(basis, indices=indices, hopping=0.7)
    first = plan.solve(tolerance=1e-12)
    second = plan.solve(tolerance=1e-12)
    assert basis.calls == 1 and svds == 1
    assert first.converged and first.stop_reason == "self_consistency_tolerance_reached"
    assert first.iterations == second.iterations
    np.testing.assert_array_equal(first.coefficients, second.coefficients)
    # Independently evolve point values without projection or any SVD.
    z = 1j * (2 * indices + 1) * np.pi / basis.beta
    direct = 1 / z
    direct_iterations = 0
    for _step in range(1, 201):
        direct_iterations += 1
        direct = 0.5 * direct + 0.5 / (z - 0.7**2 * direct)
        direct_defect = float(np.max(np.abs(direct - 1 / (z - 0.7**2 * direct))))
        if direct_defect <= 1e-12:
            break
    assert first.iterations == direct_iterations
    np.testing.assert_allclose(first.frequency_values, direct, rtol=2e-14, atol=1e-15)
    # QR supplies an independent coefficient solve, different from cached SVD.
    q, r = np.linalg.qr(basis.matrix)
    expected_coefficients = np.linalg.solve(r, q.conj().T @ direct)
    np.testing.assert_allclose(first.coefficients, expected_coefficients, rtol=2e-14)
    assert first.history[-1].defect_norm <= 1e-12
    assert max(item.projection_residual_norm for item in first.history) < 2e-15
    assert first.certification_status == "not_certified"
    json.dumps(first.metadata, allow_nan=False)
    for array in (first.coefficients, first.frequency_values, first.indices):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        first.converged = False  # type: ignore[misc]


def test_solver_accepts_precomputed_design_and_plan_without_refactor(
    monkeypatch,
) -> None:
    basis = _TinyBasis()
    fitting = LeastSquaresPlan(basis.matrix, rcond=2e-14)
    monkeypatch.setattr(
        np.linalg,
        "svd",
        lambda *args, **kwargs: pytest.fail("unexpected factorization"),
    )
    plan = BetheSolverPlan(
        basis, indices=[-2, -1, 0, 1], design=basis.matrix, fit_plan=fitting
    )
    assert basis.calls == 0
    assert plan.fit_plan is fitting
    assert plan.solve().fit_diagnostics.rcond == 2e-14


def test_solver_keeps_budget_exhaustion_and_rank_loss_visible() -> None:
    basis = _TinyBasis()
    plan = BetheSolverPlan(basis, indices=[-2, -1, 0, 1])
    result = plan.solve(max_iterations=1, tolerance=1e-15)
    assert not result.converged
    assert result.stop_reason == "iteration_budget_exhausted"
    assert result.iterations == 1 and result.history[0].defect_norm > result.tolerance
    lost = BetheSolverPlan(basis, indices=[-2, -1, 0, 1], rcond=2).solve(
        max_iterations=3
    )
    assert not lost.converged
    assert lost.fit_diagnostics.effective_rank == 0
    assert lost.fit_diagnostics.conditioning_status == "zero_rank"
    assert lost.history[-1].step_norm == 0
    assert lost.history[-1].defect_norm > lost.tolerance
    assert lost.history[-1].projection_residual_norm > 0


@pytest.mark.parametrize(
    "options",
    [
        {"mixing": 0},
        {"mixing": 1.1},
        {"mixing": True},
        {"mixing": float("nan")},
        {"tolerance": 0},
        {"tolerance": True},
        {"max_iterations": 0},
        {"max_iterations": 1.5},
        {"max_iterations": True},
    ],
)
def test_solver_rejects_invalid_controls(options: dict) -> None:
    plan = BetheSolverPlan(_TinyBasis(), indices=[-2, -1, 0, 1])
    with pytest.raises(ValueError):
        plan.solve(**options)


@pytest.mark.parametrize("hopping", [0, -1, True, float("nan"), 1e300, 1e-300])
def test_solver_rejects_invalid_or_unrepresentable_hopping(hopping: object) -> None:
    with pytest.raises(ValueError):
        BetheSolverPlan(_TinyBasis(), indices=[-2, -1, 0, 1], hopping=hopping)


def test_solver_rejects_inconsistent_prepared_design_or_fit() -> None:
    basis = _TinyBasis()
    for design in (np.eye(3), np.full((4, 4), np.nan), np.full((4, 4), "bad")):
        with pytest.raises(ValueError):
            BetheSolverPlan(basis, indices=[-2, -1, 0, 1], design=design)
    with pytest.raises(ValueError, match="fit_plan"):
        BetheSolverPlan(
            basis,
            indices=[-2, -1, 0, 1],
            fit_plan=LeastSquaresPlan(np.eye(4)),
        )


def test_public_basis_solver_matches_independent_truth_on_observations() -> None:
    basis = FermionicBasis(10, 2, tol=1e-10, precision="float64", grid_order=12)
    plan = BetheSolverPlan(basis)
    result = plan.solve()
    assert result.converged
    assert result.history[-1].defect_norm <= 1e-9
    assert result.fit_diagnostics.effective_rank == basis.rank
    # Truth is consulted only after the loop. These are observed samples,
    # explicitly not an untouched holdout or the fresh canonical final split.
    expected = semicircle_matsubara(result.indices, beta=10, hopping=1)
    np.testing.assert_allclose(result.frequency_values, expected, rtol=0, atol=1e-8)
    # Reuse already explored pilot time coordinates; reserve fresh fifths for
    # frozen canonical generation, so this regression does not inspect them.
    pilot_tau = np.array([0.0, 0.13, 1.7, 5.0, 9.87, 10.0])
    time_values = basis.design_tau(pilot_tau) @ result.coefficients
    np.testing.assert_allclose(
        time_values, semicircle_tau(pilot_tau, beta=10, hopping=1), rtol=0, atol=1e-8
    )


def test_low_frequency_hole_pilot_preserves_failed_observability_witness() -> None:
    basis = FermionicBasis(10, 2, tol=1e-10, precision="float64", grid_order=12)
    old_holes = [-12, -6, -3, -2, 1, 2, 5, 11]
    old_observations = np.setdiff1d(np.arange(-128, 128), old_holes)
    result = BetheSolverPlan(basis, indices=old_observations).solve()
    assert result.converged
    assert result.history[-1].defect_norm <= 1e-9
    assert result.fit_diagnostics.condition_number > 1e8
    # This failed pilot motivated keeping essential low modes observed. It is
    # retained as failed evidence, independent of the untouched fresh holdouts.
    low_holes = np.array([-2, 1])
    prediction = basis.design_matsubara(low_holes) @ result.coefficients
    truth = semicircle_matsubara(low_holes, beta=10, hopping=1)
    assert np.max(np.abs(prediction - truth)) > 1e-6
