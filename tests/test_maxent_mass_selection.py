"""Mass-mode validation and training-only selection; no provider execution."""

import numpy as np
import pytest

from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    _selection_metrics,
    spectral_forward,
)
from kernelgecp.types import NumericalBreakdownError


def test_legacy_mass_mode_and_explicit_constraint():
    assert MaxEntConfig().mass_mode == "rescale"
    assert MaxEntConfig(normalize_mass=None).mass_mode == "rescale"
    assert MaxEntConfig(mass_mode="constrained").normalize_mass == 1.0
    assert "mass_mode" not in MaxEntConfig().as_record()
    assert "mass_curvature" not in MaxEntConfig().as_record()
    assert "mass_interior_fraction" not in MaxEntConfig().as_record()
    assert "minimizer_marquardt" not in MaxEntConfig().as_record()
    assert MaxEntConfig(mass_mode="constrained").as_record()["mass_mode"] == (
        "constrained"
    )
    with pytest.raises(ValueError, match="explicit normalize_mass"):
        MaxEntConfig(mass_mode="constrained", normalize_mass=None)
    for invalid in ("unknown", [], True):
        with pytest.raises(ValueError, match="mass_mode"):
            MaxEntConfig(mass_mode=invalid)


def test_native_curvature_is_explicit_and_never_ignored():
    assert MaxEntConfig(mass_mode="constrained").mass_curvature == "exact"
    projected = MaxEntConfig(mass_mode="constrained", mass_curvature="projected")
    assert projected.as_record()["mass_curvature"] == "projected"
    for invalid in ("unknown", [], True):
        with pytest.raises(ValueError, match="mass_curvature"):
            MaxEntConfig(mass_mode="constrained", mass_curvature=invalid)
    with pytest.raises(ValueError, match="requires constrained"):
        MaxEntConfig(mass_curvature="projected")


def test_interior_fraction_is_explicit_and_never_ignored():
    assert MaxEntConfig().mass_interior_fraction == 0.0
    mixture = MaxEntConfig(mass_mode="constrained", mass_interior_fraction=1e-12)
    assert mixture.as_record()["mass_interior_fraction"] == 1e-12
    for invalid in (-1, 1, True, [], np.nan, 1e-30):
        with pytest.raises(ValueError, match="mass_interior_fraction"):
            MaxEntConfig(mass_mode="constrained", mass_interior_fraction=invalid)
    with pytest.raises(ValueError, match="requires constrained"):
        MaxEntConfig(mass_interior_fraction=1e-12)


def test_native_damping_flag_is_explicit_and_boolean():
    assert MaxEntConfig().minimizer_marquardt is False
    assert (
        MaxEntConfig(minimizer_marquardt=True).as_record()["minimizer_marquardt"]
        is True
    )
    for invalid in (0, 1, "true", None):
        with pytest.raises(ValueError, match="minimizer_marquardt"):
            MaxEntConfig(minimizer_marquardt=invalid)


def test_constraint_rejects_mass_violation_instead_of_rescaling():
    omega = np.array([-1.0, 0.0, 1.0])
    spectra = np.full((1, 3), 0.55)
    data = SpectralDataset(2.0, [0.0, 1.0, 2.0], np.zeros(3), np.eye(3))
    args = (data, omega, spectra, np.array([1.0]), np.array([True]), 1.0)
    with pytest.raises(NumericalBreakdownError, match="no positive finite-mass"):
        _selection_metrics(*args, mass_mode="constrained")
    selected, masses, metrics = _selection_metrics(*args)
    assert selected == 0 and masses[0] == pytest.approx(1.1)
    assert metrics[0] is not None


def test_constraint_selection_does_not_modify_rounding_scale():
    omega = np.array([-1.0, 0.0, 1.0])
    spectra = np.full((1, 3), 0.5 * (1 + 5e-12))
    tau = np.array([0.0, 1.0, 2.0])
    observed = spectral_forward(2.0, tau, omega, spectra[0])
    data = SpectralDataset(2.0, tau, observed, np.eye(3) * 1e-16)
    args = (data, omega, spectra, np.array([1.0]), np.array([True]), 1.0)
    _, _, constrained = _selection_metrics(*args, mass_mode="constrained")
    _, _, rescaled = _selection_metrics(*args)
    assert constrained[0] == 0.0
    assert rescaled[0] is not None and rescaled[0] > 0
