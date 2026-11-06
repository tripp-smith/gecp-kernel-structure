"""Independent arithmetic checks for the opt-in projected metric."""

import numpy as np
import pytest

from kernelgecp.integrations._maxent_curvature import shifted_projected_metric
from kernelgecp.integrations._maxent_native_mass import make_native_mass_cost
from kernelgecp.integrations.maxent import MaxEntConfig
from kernelgecp.types import NumericalBreakdownError


def test_shift_independent_two_by_two_formula_and_ownership():
    matrix = np.array([[2.0, 1.2], [0.8, 0.5]])
    old = matrix.copy()
    floor = 0.01
    metric, info = shifted_projected_metric(matrix, floor)
    # Eigenvalues of [[1, .5], [.5, .25]], divided by scale 2.
    shift = floor * 2
    expected = np.array([[2.0, 1.0], [1.0, 0.5]]) + 2 * shift * np.eye(2)
    np.testing.assert_allclose(metric, expected, rtol=1e-14, atol=0)
    assert np.array_equal(matrix, old) and not np.shares_memory(matrix, metric)
    assert info["normalized_shift"] == pytest.approx(shift)
    assert info["condition_estimate"] == pytest.approx(1.27 / 0.02)
    assert info["certified"] is False


@pytest.mark.parametrize(
    "matrix",
    [
        np.zeros((2, 2)),
        np.array([[np.inf]]),
        np.array([[np.nan]]),
        np.array([[np.nextafter(0.0, 1.0)]]),
        np.array(
            [
                [np.finfo(float).max, np.finfo(float).max],
                [np.finfo(float).max, -np.finfo(float).max],
            ]
        ),
    ],
)
def test_unrepresentable_inputs_rejected(matrix):
    with pytest.raises(NumericalBreakdownError):
        shifted_projected_metric(matrix, 1e-2)


@pytest.mark.parametrize(
    "matrix",
    [np.ones((129, 129)), np.ones((2, 3)), np.array([[1j]]), np.array([[True]])],
)
def test_invalid_shape_or_type_rejected(matrix):
    with pytest.raises(ValueError):
        shifted_projected_metric(matrix, 1e-10)


@pytest.mark.parametrize(
    "floor", [True, np.bool_(True), -1, 0.1, np.nan, np.inf, "1e-10"]
)
def test_invalid_floor_rejected_before_provider_import(floor):
    with pytest.raises(ValueError, match="mass_curvature_floor_relative"):
        MaxEntConfig(mass_curvature_floor_relative=floor)
    with pytest.raises(ValueError, match="mass_curvature_floor_relative"):
        make_native_mass_cost(
            1.0, curvature="projected", curvature_floor_relative=floor
        )


def test_modes_and_legacy_records():
    for kwargs in (
        {},
        {"mass_mode": "constrained"},
        {"mass_mode": "constrained", "mass_curvature": "projected"},
    ):
        assert (
            MaxEntConfig(**kwargs).as_record()
            == MaxEntConfig(**kwargs, mass_curvature_floor_relative=0.0).as_record()
        )
        assert "mass_curvature_floor_relative" not in MaxEntConfig(**kwargs).as_record()
    for kwargs in ({}, {"mass_mode": "constrained"}):
        with pytest.raises(ValueError, match="constrained projected"):
            MaxEntConfig(**kwargs, mass_curvature_floor_relative=1e-10)
    config = MaxEntConfig(
        mass_mode="constrained",
        mass_curvature="projected",
        mass_curvature_floor_relative=1e-10,
    )
    assert config.as_record()["mass_curvature_floor_relative"] == 1e-10


def test_mixed_scale_entries_and_equal_subnormal_are_preserved():
    tiny = np.nextafter(0.0, 1.0)
    for offdiagonal in (1e-20, 3 * tiny):
        matrix = np.array([[1e308, offdiagonal], [offdiagonal, 1e308]])
        metric, info = shifted_projected_metric(matrix, 1e-10)
        assert info["normalized_shift"] == 0
        assert metric[0, 1] == metric[1, 0] == offdiagonal


def test_unequal_unrepresentable_half_sum_fails_closed():
    tiny = np.nextafter(0.0, 1.0)
    matrix = np.array([[1.0, 3 * tiny], [tiny, 1.0]])
    with pytest.raises(NumericalBreakdownError, match="half sum"):
        shifted_projected_metric(matrix, 1e-10)
