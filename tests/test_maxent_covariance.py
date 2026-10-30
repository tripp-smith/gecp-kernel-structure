"""Independent covariance/forward checks; actual optional provider is separate."""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from scipy.integrate import quad

from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    continue_tau_maxent,
    covariance_modes,
    diagonal_spectral_datasets,
    load_spectral_dataset,
    prepare_spectral_projection,
    save_spectral_dataset,
    spectral_forward,
)


class PolynomialDesign:
    beta = 3.0

    def design_tau(self, tau):
        return np.polynomial.legendre.legvander(2 * np.asarray(tau) / self.beta - 1, 3)


def dataset():
    tau = np.linspace(0, 3, 9)
    correlation = np.exp(-np.abs(tau[:, None] - tau) / 0.8)
    return SpectralDataset(3, tau, -0.5 + 0.1 * np.cos(tau), 1e-6 * correlation)


def test_complete_projection_matches_independent_qr_covariance_and_noise_modes():
    data = dataset()
    target = np.linspace(0, 3, 17)
    basis = PolynomialDesign()
    plan = prepare_spectral_projection(basis, data.tau, target)
    q, r = np.linalg.qr(basis.design_tau(data.tau))
    independent = np.linalg.solve(r.T, basis.design_tau(target).T).T @ q.T
    projected = plan.apply(data)
    np.testing.assert_allclose(plan.matrix, independent, atol=2e-15)
    np.testing.assert_allclose(projected.values, independent @ data.values, atol=1e-15)
    expected = independent @ data.covariance @ independent.T
    np.testing.assert_allclose(projected.covariance, expected, atol=3e-21)
    assert projected.modes.effective_rank == 4
    assert projected.source_fingerprint == data.fingerprint
    # Treating interpolated points as independent changes the noise model.
    assert np.linalg.norm(expected - np.diag(np.diag(expected))) > 1e-6
    residual = np.linspace(-1, 1, target.size) * 1e-4
    eigen, vectors = np.linalg.eigh(expected)
    keep = eigen > projected.modes.threshold
    oracle = np.sqrt(np.mean((vectors[:, keep].T @ residual) ** 2 / eigen[keep]))
    assert projected.modes.whitened_rms(residual) == pytest.approx(oracle, abs=2e-15)


def test_covariance_explicit_rank_rejection_and_roundoff_symmetrization():
    matrix, modes = covariance_modes(np.diag([1.0, 1e-9, 0.0]), rtol=1e-8)
    assert modes.effective_rank == 1
    assert modes.whitened_rms([2, 300, -500]) == 2
    np.testing.assert_array_equal(matrix, np.diag([1, 1e-9, 0]))
    for invalid in (np.diag([1, -1e-4]), [[1, 0.1], [0, 1]], np.zeros((2, 2))):
        with pytest.raises(ValueError):
            covariance_modes(invalid)
    for invalid in (0, -1, 1, True, np.nan):
        with pytest.raises(ValueError):
            covariance_modes(np.eye(3), rtol=invalid)
    close = np.array([[1, 0.2 + 1e-15], [0.2, 1]])
    _, modes = covariance_modes(close)
    assert modes.symmetry_correction > 0


def test_same_mesh_projection_is_stable_for_an_ill_scaled_coefficient_basis():
    class BadlyScaledDesign:
        beta = 3.0

        def design_tau(self, tau):
            x = np.asarray(tau) / self.beta
            return np.column_stack((np.ones_like(x), x, 1e-9 * x**2))

    data = dataset()
    plan = prepare_spectral_projection(
        BadlyScaledDesign(), data.tau, data.tau, rcond=1e-12
    )
    # Independent unscaled QR spans the same polynomial observations.
    q, _ = np.linalg.qr(np.column_stack((np.ones(9), data.tau, data.tau**2)))
    np.testing.assert_allclose(plan.matrix, q @ q.T, atol=2e-15)
    np.testing.assert_allclose(plan.matrix, plan.matrix.T, atol=1e-16)
    np.testing.assert_allclose(plan.matrix @ plan.matrix, plan.matrix, atol=1e-15)


def test_chained_projection_preserves_complete_original_observation_operator():
    data = dataset()
    first = prepare_spectral_projection(
        PolynomialDesign(), data.tau, np.linspace(0, 3, 17)
    )
    intermediate = first.apply(data)
    second = prepare_spectral_projection(PolynomialDesign(), intermediate.tau, data.tau)
    result = second.apply(intermediate)
    expected = second.matrix @ first.matrix
    np.testing.assert_allclose(result.projection, expected, atol=1e-15)
    np.testing.assert_array_equal(result.projection_source_tau, data.tau)
    np.testing.assert_allclose(result.values, expected @ data.values, atol=1e-15)
    np.testing.assert_allclose(
        result.covariance, expected @ data.covariance @ expected.T, atol=1e-21
    )
    assert result.source_fingerprint == data.fingerprint


def test_owned_data_projection_and_no_implicit_complex_projection():
    data = dataset()
    for array in (data.values, data.tau, data.covariance, data.modes.retained_vectors):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        data.beta = 7
    with pytest.raises(ValueError, match="real-part"):
        SpectralDataset(3, data.tau, data.values.astype(complex), data.covariance)
    for bad in (data.tau[1:], data.tau[::-1], np.zeros(9)):
        with pytest.raises(ValueError):
            SpectralDataset(3, bad, data.values, data.covariance)
    plan = prepare_spectral_projection(PolynomialDesign(), data.tau, data.tau)
    with pytest.raises(ValueError, match="exact prepared"):
        plan.apply(
            SpectralDataset(4, np.linspace(0, 4, 9), data.values, data.covariance)
        )


def test_diagonal_channels_extract_declared_time_major_marginals():
    data = dataset()
    values = np.column_stack((data.values, data.values * 2))
    joint = np.kron(data.covariance, np.array([[1, 0.25], [0.25, 2]]))
    first, second = diagonal_spectral_datasets(
        3, data.tau, values, joint, ("up", "down"), provenance={"run": "native"}
    )
    np.testing.assert_array_equal(first.values, data.values)
    np.testing.assert_allclose(second.covariance, 2 * data.covariance)
    assert first.provenance["channel"] == "up"
    assert first.provenance["source"] == {"run": "native"}
    assert second.provenance["ordering"] == ["up", "down"]
    with pytest.raises(ValueError):
        diagonal_spectral_datasets(3, data.tau, values, joint, ("up", "up"))


def test_spectral_forward_independent_adaptive_quadrature_and_endpoint_sum_rule():
    omega = np.linspace(-2, 2, 41)
    density = np.maximum(1 - np.abs(omega) / 2, 0) / 2
    beta = 10.0
    tau = np.array([0, 0.2, 1.7, 5, 8.4, 10])
    actual = spectral_forward(beta, tau, omega, density)
    expected = [
        quad(
            lambda w, t=t: (
                -np.exp(-t * w) / (1 + np.exp(-beta * w)) * np.interp(w, omega, density)
            ),
            -2,
            2,
            points=omega,
            epsabs=1e-13,
        )[0]
        for t in tau
    ]
    np.testing.assert_allclose(actual, expected, atol=1e-14)
    assert actual[0] + actual[-1] == pytest.approx(-np.trapezoid(density, omega))
    assert np.all(np.isfinite(spectral_forward(1e4, [0, 1e4], omega, density)))


def test_archive_safe_roundtrip_and_tamper_rejection(tmp_path):
    data = prepare_spectral_projection(
        PolynomialDesign(), dataset().tau, np.linspace(0, 3, 17)
    ).apply(dataset())
    target = tmp_path / "observations.npz"
    manifest = save_spectral_dataset(
        data,
        target,
        units={"tau": "Hartree^-1"},
        provider={"name": "test", "version": "1"},
    )
    result = load_spectral_dataset(target)
    assert result.fingerprint == data.fingerprint
    np.testing.assert_array_equal(result.projection, data.projection)
    with pytest.raises(FileExistsError):
        save_spectral_dataset(data, target, units={"tau": "1"}, provider={"name": "x"})
    document = json.loads(manifest.read_text())
    document["effective_rank"] += 1
    manifest.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="rank mismatch"):
        load_spectral_dataset(target)
    manifest.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError, match="duplicate"):
        load_spectral_dataset(target)
    document["effective_rank"] -= 1
    manifest.write_text(json.dumps(document))
    target.write_bytes(target.read_bytes() + b"tamper")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_spectral_dataset(target)


@pytest.mark.parametrize(
    "settings",
    [
        dict(n_alpha=True),
        dict(n_omega=2),
        dict(alpha_min=2, alpha_max=1),
        dict(default_model="silent_fallback"),
        dict(normalize_mass=0),
        *[
            dict(minimizer_mu0=value)
            for value in (True, np.bool_(False), 0, -1, np.nan, np.inf, -np.inf, "1e-8")
        ],
    ],
)
def test_resource_and_configuration_contract(settings):
    with pytest.raises(ValueError):
        MaxEntConfig(**settings)


def test_no_optional_provider_needed_and_missing_provider_fails(monkeypatch):
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import kernelgecp; "
            "import kernelgecp.integrations.maxent; import sys; "
            "assert 'triqs_maxent' not in sys.modules",
        ],
        check=True,
    )
    original = importlib.import_module

    def missing(name, *args, **kwargs):
        if name.startswith("triqs_maxent"):
            raise ImportError("optional provider absent")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", missing)
    with pytest.raises(ImportError, match="optional"):
        continue_tau_maxent(dataset())


def test_actual_tau_maxent_provider_positive_and_rank_aware():
    pytest.importorskip("triqs_maxent")
    beta = 10
    tau = np.linspace(0, beta, 33)
    pole = 0.7
    values = -(np.exp(-pole * tau) + np.exp(-pole * (beta - tau))) / (
        2 * (1 + np.exp(-beta * pole))
    )
    data = SpectralDataset(beta, tau, values, 1e-8 * np.eye(33))
    settings = {
        "omega_min": -2,
        "omega_max": 2,
        "n_omega": 81,
        "n_alpha": 9,
        "minimizer_mu0": 1e-8,
    }
    result = continue_tau_maxent(data, MaxEntConfig(**settings))
    assert result.diagnostics["config"]["minimizer_mu0"] == 1e-8
    assert result.diagnostics["version"] == "4.0.1"
    assert result.effective_rank == 33
    assert result.converged[result.selected_index]
    diagnostics = result.diagnostics
    assert diagnostics["selection_chi2"][
        result.selected_index
    ] / result.effective_rank == pytest.approx(
        result.normalized_training_whitened_rms**2
    )
    diagnostics["selection_chi2"][result.selected_index] = -1
    assert result.diagnostics["selection_chi2"][result.selected_index] >= 0
    assert np.min(result.density) >= 0
    assert abs(result.mass - 1) < 1e-12
    assert abs(result.raw_mass - 1) < 0.002
    assert result.normalized_training_whitened_rms < 3
    # This small integration case does not close the independent science gate.
    assert np.max(np.abs(result.evaluate_tau(beta, tau) - values)) < 1e-3


def test_actual_provider_projects_forward_kernel_with_complete_observation_map(
    monkeypatch,
):
    provider = pytest.importorskip("triqs_maxent")
    created = []
    actual = provider.TauMaxEnt

    def capture(*args, **kwargs):
        instance = actual(*args, **kwargs)
        created.append(instance)
        return instance

    monkeypatch.setattr(provider, "TauMaxEnt", capture)
    tau = dataset().tau
    values = -np.exp(-0.7 * tau) / (1 + np.exp(-3 * 0.7))
    data = SpectralDataset(3, tau, values, 1e-6 * np.eye(len(tau)))
    plan = prepare_spectral_projection(PolynomialDesign(), tau, np.linspace(0, 3, 17))
    projected = plan.apply(data)
    result = continue_tau_maxent(
        projected,
        MaxEntConfig(
            omega_min=-2, omega_max=2, n_omega=81, n_alpha=9, minimizer_mu0=1e-8
        ),
    )
    solver = created[0]
    independent_kernel = -np.exp(-tau[:, None] * result.omega) / (
        1 + np.exp(-3 * result.omega)
    )
    expected = plan.matrix @ independent_kernel
    target_kernel = -np.exp(-projected.tau[:, None] * result.omega) / (
        1 + np.exp(-3 * result.omega)
    )
    assert np.max(np.abs(expected - target_kernel)) > 1e-3
    np.testing.assert_allclose(
        solver.K.K_delta, expected * solver.omega.delta, atol=1e-16
    )
    assert result.effective_rank == 4
    assert result.diagnostics["forward_operator"] == (
        "complete_projection_times_native_source_tau_kernel"
    )


def test_selection_scores_delivered_spectrum_and_projection_not_provider_chi2():
    from kernelgecp.integrations.maxent import _selection_metrics
    from kernelgecp.types import NumericalBreakdownError

    omega = np.linspace(-2, 2, 81)
    tau = np.linspace(0, 3, 9)
    good = np.exp(-(((omega - 0.6) / 0.5) ** 2))
    good /= np.trapezoid(good, omega)
    bad = np.exp(-(((omega + 1.0) / 0.35) ** 2))
    bad /= np.trapezoid(bad, omega)
    values = spectral_forward(3, tau, omega, good) + 0.001
    data = SpectralDataset(3, tau, values, 1e-6 * np.eye(9))
    trials = np.array([1.02 * bad, good, np.zeros_like(good), -good, good])
    native_chi = np.array([9.0, 100.0, 9.0, 9.0, 9.0])
    converged = np.array([True, True, True, True, False])
    for training in (
        data,
        prepare_spectral_projection(
            PolynomialDesign(), tau, np.linspace(0, 3, 17)
        ).apply(data),
    ):
        selected, masses, metrics = _selection_metrics(
            training, omega, trials, native_chi, converged, 1.0
        )
        assert selected == 1
        assert masses[0] == pytest.approx(1.02)
        assert metrics[2:] == [None, None, None]
        prediction = spectral_forward(3, tau, omega, good)
        if training.projection is not None:
            prediction = training.projection @ prediction
        expected = training.modes.whitened_rms(prediction - training.values) ** 2
        assert metrics[1] / training.modes.effective_rank == pytest.approx(expected)
        assert (
            _selection_metrics(training, omega, trials, native_chi, converged, None)[0]
            == 0
        )
        with pytest.raises(NumericalBreakdownError, match="finite-mass"):
            _selection_metrics(
                training, omega, trials[2:4], native_chi[2:4], converged[2:4], 1.0
            )


@pytest.mark.parametrize("mu0", [None, 2e-7])
def test_provider_failure_retains_cause_and_completed_trial_history(monkeypatch, mu0):
    from types import SimpleNamespace

    from kernelgecp.types import NumericalBreakdownError

    constructor_calls = []

    class Minimizer:
        def __init__(self, **kwargs):
            constructor_calls.append(kwargs)

        def minimize(self):
            self.converged, self.n_iter_last = True, 7

    class Solver:
        def __init__(self, **kwargs):
            self.minimizer = kwargs["minimizer"]

        def set_G_tau_data(self, tau, values):
            pass

        def set_cov(self, covariance):
            self.err = np.ones(len(covariance))

        def run(self):
            self.minimizer.minimize()
            np.seterr(all="ignore")
            raise np.linalg.LinAlgError("singular training system")

    provider = SimpleNamespace(
        TauMaxEnt=Solver,
        LinearOmegaMesh=lambda lo, hi, n: np.linspace(lo, hi, n),
        LogAlphaMesh=lambda *args: None,
        DataDefaultModel=lambda *args: None,
    )
    modules = {
        "triqs_maxent": provider,
        "triqs_maxent.minimizers": SimpleNamespace(LevenbergMinimizer=Minimizer),
        "triqs_maxent.version": SimpleNamespace(version="4.0.1"),
    }
    original = importlib.import_module
    monkeypatch.setattr(
        importlib,
        "import_module",
        lambda name: modules[name] if name in modules else original(name),
    )
    previous_policy = np.geterr()
    with pytest.raises(NumericalBreakdownError) as caught:
        settings = {"maxiter": 123}
        if mu0 is not None:
            settings["minimizer_mu0"] = mu0
        continue_tau_maxent(dataset(), MaxEntConfig(**settings))
    assert constructor_calls == [{"maxiter": 123, "mu0": 1e-18 if mu0 is None else mu0}]
    assert np.geterr() == previous_policy
    assert "LinAlgError: singular training system" in str(caught.value)
    assert "completed trials: 1" in str(caught.value)
    assert "history (converged, iterations): [(True, 7)]" in str(caught.value)
    assert isinstance(caught.value.__cause__, np.linalg.LinAlgError)


def test_existing_positional_configuration_keeps_its_meaning():
    config = MaxEntConfig(-4, 4, 401, 1e-4, 1e4, 33, 1000, "broad_gaussian", None)
    assert config.default_model == "broad_gaussian"
    assert config.normalize_mass is None
    assert config.minimizer_mu0 == 1e-18


def test_historical_spectral_dataset_load_preserves_original_manifest():
    import json
    from pathlib import Path

    path = (
        Path(__file__).parent
        / "data/maxent/linux-singular-e8ac40e"
        / "offgrid_atoms_correlated-pydlr-broad_gaussian-r0-dataset.npz"
    )
    original = path.with_suffix(".json").read_bytes()
    manifest = json.loads(original)
    assert "minimizer_mu0" not in manifest["calculation_settings"]
    restored = load_spectral_dataset(path)
    assert restored.fingerprint == manifest["fingerprint"]
    assert path.with_suffix(".json").read_bytes() == original
