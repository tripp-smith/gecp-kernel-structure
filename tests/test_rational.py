"""Independent transform identities, resolution checks and physical AAA poles."""

import json

import mpmath as mp
import numpy as np
import pytest

from kernelgecp.kernels import FermionicKernel
from kernelgecp.rational import aaa_atoms, matsubara_transform


def test_fermionic_matsubara_identity_independent_quadrature() -> None:
    with mp.workdps(60):
        for omega in (-8, -0.35, 0, 1.7, 8):
            for mode in (-3, 0, 2):
                z = mp.j * mp.pi * (2 * mode + 1)
                integral = -mp.quad(
                    lambda t, z=z, omega=omega: (
                        mp.exp(z * t) * mp.exp(-t * omega) / (1 + mp.exp(-omega))
                    ),
                    [0, 1],
                )
                assert abs(integral - 1 / (z - omega)) < mp.mpf("1e-55")


def test_shared_observation_transform_ladder_not_oracle() -> None:
    poles = np.array([-5.4, -2.2, 0.8, 4.5])
    weights = np.array([0.18, 0.31, 0.29, 0.22])
    errors = []
    for count in (33, 65, 129, 257):
        t = np.linspace(0, 1, count)
        observed = FermionicKernel(8)(t[:, None], poles) @ weights
        z, transformed, diagnostics = matsubara_transform(t, observed, modes=8)
        exact = (1 / (z[:, None] - poles)) @ weights
        errors.append(float(np.max(np.abs(transformed - exact))))
        assert diagnostics["sample_count"] == count
        assert diagnostics["certification_status"] == "not_certified"
        assert diagnostics["source_ladder_max_abs_difference"] > 0
        assert diagnostics["conjugate_symmetry_max_abs_error"] < 1e-14
        assert diagnostics["coarse_fine_max_abs_difference"] < 1e-10
    assert all(new < old / 5 for old, new in zip(errors, errors[1:], strict=False))
    assert errors[-1] < 1e-10


def test_aaa_recovers_actual_spectral_poles_from_rational_data() -> None:
    z = 1j * np.pi * (2 * np.arange(-16, 16) + 1)
    poles = np.array([-5.4, -2.2, 0.8, 4.5])
    values = (1 / (z[:, None] - poles)) @ np.array([0.18, 0.31, 0.29, 0.22])
    result = aaa_atoms(z, values, cutoff=8, max_atoms=4, tolerance=1e-12)
    assert result.metadata.converged
    assert np.allclose(result.frequencies, poles, atol=1e-7, rtol=0)
    assert result.diagnostics["barycentric_support_count"] <= 5
    assert result.diagnostics["rejected_pole_count"] == 0
    assert "native_complex" in result.metadata.convergence_scope
    assert result.metadata.certification_status == "not_certified"
    result.metadata.canonical_json()
    json.dumps(result.diagnostics, allow_nan=False)


@pytest.mark.parametrize("pole,reason", [(10, "outside_cutoff"), (1 + 0.5j, "nonreal")])
def test_aaa_nonphysical_poles_are_not_clipped(pole: complex, reason: str) -> None:
    z = 1j * np.arange(1, 17)
    result = aaa_atoms(z, 1 / (z - pole), cutoff=4, max_atoms=1, tolerance=1e-12)
    assert result.metadata.converged
    assert not result.frequencies.size
    assert result.diagnostics["rejected_poles"][0]["reason"] == reason


def test_aaa_zero_budget_zero_data_and_repeat() -> None:
    z = 1j * np.arange(1, 9)
    empty = aaa_atoms(z, np.zeros(8), cutoff=4, max_atoms=4)
    assert empty.metadata.converged and not empty.frequencies.size
    limited = aaa_atoms(z, 1 / (z - 1.3), cutoff=4, max_atoms=0)
    assert not limited.metadata.converged and not limited.frequencies.size
    first = aaa_atoms(z, 1 / (z - 1.3), cutoff=4, max_atoms=1)
    second = aaa_atoms(z, 1 / (z - 1.3), cutoff=4, max_atoms=1)
    assert first.metadata.canonical_json() == second.metadata.canonical_json()
    assert np.array_equal(first.frequencies, second.frequencies)


@pytest.mark.parametrize(
    "t,values",
    [
        ([0, 0.2, 0.8, 1], [1, 2, 3, 4]),
        ([0, 0.5, 1], [1, 2, 3]),
        ([0, 1 / 3, 2 / 3, 1], [1, np.nan, 3, 4]),
        ([0, 1 / 3, 2 / 3, 1], [1, 2, 3]),
    ],
)
def test_transform_rejects_invalid_shared_samples(t: object, values: object) -> None:
    with pytest.raises(ValueError):
        matsubara_transform(t, values)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_atoms": True},
        {"max_atoms": -1},
        {"cutoff": np.inf},
        {"tolerance": np.nan},
        {"pole_imag_tolerance": -1},
    ],
)
def test_aaa_rejects_bad_controls(kwargs: dict) -> None:
    controls = {"max_atoms": 2, "cutoff": 4, **kwargs}
    with pytest.raises(ValueError):
        aaa_atoms([1j, 2j], [1, 2], **controls)


def test_aaa_rejects_nonfinite_and_duplicate_samples() -> None:
    for z, values in (([1j, 1j], [1, 2]), ([1j, 2j], [np.inf, 2])):
        with pytest.raises(ValueError):
            aaa_atoms(z, values, cutoff=4, max_atoms=1)
