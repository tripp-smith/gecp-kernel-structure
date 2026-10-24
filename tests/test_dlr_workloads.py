"""Independent physical-pole truth and isolated common-mesh checks."""

from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from experiments.dlr_workloads import (
    PoleWorkload,
    common_tau_meshes,
    green_function,
    make_pole_workloads,
    validate_tau_meshes,
)
from kernelgecp.grids import time_grid


def _mp_direct_truth(workload: PoleWorkload, tau: np.ndarray, beta: float) -> list:
    """Use the literal physical expression, including large positive exponents."""

    with mp.workdps(100):
        exact_beta = mp.mpf(beta)
        return [
            -mp.fsum(
                mp.mpf(weight)
                * mp.exp(-mp.mpf(time) * mp.mpf(pole))
                / (1 + mp.exp(-exact_beta * mp.mpf(pole)))
                for pole, weight in zip(workload.poles, workload.weights, strict=True)
            )
            for time in tau
        ]


@pytest.mark.parametrize("beta,wmax", [(0.25, 2.0), (8.0, 3.0), (1000.0, 100.0)])
def test_finite_pole_values_match_independent_high_precision_expression(
    beta: float, wmax: float
) -> None:
    tau = beta * np.array([0.0, 1e-8, 0.013, 0.29, 0.5, 0.91, 1.0])
    for workload in make_pole_workloads(beta=beta, wmax=wmax):
        expected = np.array(
            [float(value) for value in _mp_direct_truth(workload, tau, beta)]
        )
        # The independent expression never calls a repository kernel or fit.
        actual = green_function(workload, tau, beta=beta, wmax=wmax)
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=1e-16)
        assert actual.dtype == np.float64
        assert actual.shape == tau.shape
        assert np.all(actual <= 0)


def test_exact_central_pole_and_normalized_endpoint_identity() -> None:
    for beta in (1e-3, 1.0, 1e6):
        fixtures = make_pole_workloads(beta=beta, wmax=2.0)
        central = green_function(fixtures[0], [0.0, beta / 3, beta], beta=beta, wmax=2)
        np.testing.assert_array_equal(central, [-0.5, -0.5, -0.5])
        for workload in fixtures:
            endpoints = green_function(workload, [0.0, beta], beta=beta, wmax=2)
            assert float(endpoints.sum()) == pytest.approx(-1.0, rel=0, abs=5e-16)


def test_extreme_scales_and_numpy_underflow_traps_are_safe() -> None:
    workload = PoleWorkload("extreme", (-1e-300, 1e-300), (0.4, 0.6))
    times = np.array([0.0, 5e299, 1e300])
    actual = green_function(workload, times, beta=1e300, wmax=1e-300)
    expected = np.array(
        [float(value) for value in _mp_direct_truth(workload, times, 1e300)]
    )
    np.testing.assert_allclose(actual, expected, rtol=3e-16)
    with np.errstate(all="raise"):
        trapped = green_function(
            PoleWorkload("underflow", (-100.0, 100.0), (0.5, 0.5)),
            [0.0, 5.0, 10.0],
            beta=10,
            wmax=100,
        )
    assert np.all(np.isfinite(trapped))


def test_fixture_supports_weights_and_scaling_are_fixed() -> None:
    first = make_pole_workloads(beta=1, wmax=2)
    second = make_pole_workloads(beta=100, wmax=4)
    assert [workload.name for workload in first] == [
        "central_pole",
        "asymmetric_satellites",
        "gapped",
    ]
    for left, right in zip(first, second, strict=True):
        assert right.weights == left.weights
        np.testing.assert_array_equal(right.poles, 2 * np.asarray(left.poles))
        assert all(abs(pole) <= 2 for pole in left.poles)
        assert all(weight > 0 for weight in left.weights)
        assert sum(left.weights) == pytest.approx(1.0, abs=2e-16)
    assert all(abs(pole) >= 0.3 * 2 for pole in first[2].poles)


@pytest.mark.parametrize(
    "poles,weights",
    [
        ((), ()),
        (0.0, (1.0,)),
        (((0.0,),), (1.0,)),
        ((0.0,), 1.0),
        ((0.0,), (0.4, 0.6)),
        ((0.0, 0.0), (0.4, 0.6)),
        ((float("nan"),), (1.0,)),
        ((float("inf"),), (1.0,)),
        ((1j,), (1.0,)),
        ((0.0,), (0.0,)),
        ((0.0, 1.0), (-0.2, 1.2)),
        ((0.0,), (float("nan"),)),
        ((0.0,), (0.9,)),
        ((0.0, 1.0), (1e308, 1e308)),
    ],
)
def test_workload_rejects_malformed_shapes_and_nonphysical_measure(
    poles: object, weights: object
) -> None:
    with pytest.raises(ValueError):
        PoleWorkload("bad", poles, weights)  # type: ignore[arg-type]


@pytest.mark.parametrize("name", ["", " ", 1])
def test_workload_requires_nonempty_name(name: object) -> None:
    with pytest.raises(ValueError, match="name"):
        PoleWorkload(name, (0.0,), (1.0,))  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "value", [0, -1, float("nan"), float("inf"), True, [1.0], np.array(1.0)]
)
@pytest.mark.parametrize("parameter", ["beta", "wmax"])
def test_physical_scalar_validation(value: object, parameter: str) -> None:
    parameters = {"beta": 1.0, "wmax": 1.0, parameter: value}
    for function in (make_pole_workloads, common_tau_meshes):
        with pytest.raises(ValueError):
            function(**parameters)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        green_function(
            PoleWorkload("central", (0.0,), (1.0,)),
            [0.0],
            **parameters,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("beta,wmax", [(1e308, 1e308), (1e-300, 1e-300)])
def test_unrepresentable_dimensionless_cutoff_is_rejected(
    beta: float, wmax: float
) -> None:
    with pytest.raises(ValueError, match="beta \\* wmax"):
        make_pole_workloads(beta=beta, wmax=wmax)
    with pytest.raises(ValueError, match="beta \\* wmax"):
        common_tau_meshes(beta=beta, wmax=wmax)


@pytest.mark.parametrize(
    "tau",
    [0.0, np.array(0.0), [[0.0]], [], [float("nan")], [-1e-8], [1.00000001], [1j]],
)
def test_green_function_rejects_invalid_time_shapes_and_domain(tau: object) -> None:
    with pytest.raises(ValueError):
        green_function(PoleWorkload("central", (0.0,), (1.0,)), tau, beta=1, wmax=1)


def test_green_function_rejects_poles_outside_declared_physical_cutoff() -> None:
    with pytest.raises(ValueError, match="poles"):
        green_function(PoleWorkload("outside", (1.1,), (1.0,)), [0.5], beta=1, wmax=1)


@pytest.mark.parametrize("beta,wmax", [(0.2, 1.0), (1.0, 16.0), (100.0, 10.24)])
def test_common_meshes_are_canonical_endpoint_aware_disjoint_and_repeatable(
    beta: float, wmax: float
) -> None:
    training, validation = common_tau_meshes(beta=beta, wmax=wmax, order=8)
    repeated = common_tau_meshes(beta=beta, wmax=wmax, order=8)
    assert training.tobytes() == repeated[0].tobytes()
    assert validation.tobytes() == repeated[1].tobytes()
    assert training[0] == 0
    assert training[-1] == beta
    expected = time_grid(8, max(1, beta * wmax)) * beta
    expected[0], expected[-1] = 0, beta
    np.testing.assert_array_equal(training, expected)
    assert np.all(training[:-1] < validation)
    assert np.all(validation < training[1:])
    assert not np.intersect1d(training, validation).size
    assert training[1] < min(beta / 2, 1 / wmax)
    validate_tau_meshes(training, validation, beta=beta)


@pytest.mark.parametrize("order", [True, 1, 0, 2.5, np.array(8)])
def test_common_meshes_reject_noninteger_or_too_small_order(order: object) -> None:
    with pytest.raises(ValueError, match="order"):
        common_tau_meshes(beta=1, wmax=1, order=order)  # type: ignore[arg-type]


def test_nonbinary_panel_boundaries_do_not_duplicate_or_overlap_holdout() -> None:
    training, validation = common_tau_meshes(beta=1, wmax=10, order=8)
    # Four left panels [0,.1], [.1,.2], [.2,.4], [.4,.5] and reflection.
    assert len(training) == 2 * 4 * (8 - 1) + 1
    assert np.count_nonzero(training == 0.1) == 1
    assert np.count_nonzero(training == 0.9) == 1
    assert np.all(training[:-1] < validation)
    assert np.all(validation < training[1:])
    assert not np.intersect1d(training, validation).size


@pytest.mark.parametrize(
    "training,validation",
    [
        ([0.0, 0.0, 1.0], [0.5]),
        ([0.0, 1.0], [0.5, 0.5]),
        ([0.0, 0.5, 1.0], [0.5]),
        ([1.0, 0.0], [0.5]),
        ([0.0, 1.0], [1.01]),
        ([0.0, float("nan")], [0.5]),
        ([0.0, 1.0], 0.5),
        (0.0, [0.5]),
        ([[0.0, 1.0]], [0.5]),
        ([0.0, 1.0], []),
    ],
)
def test_mesh_validation_rejects_duplicate_overlap_and_invalid_shape(
    training: object, validation: object
) -> None:
    with pytest.raises(ValueError):
        validate_tau_meshes(training, validation, beta=1)


def test_common_meshes_reject_float64_coordinate_resolution_exhaustion() -> None:
    with pytest.raises(ValueError, match="duplicate|disjoint|increasing"):
        common_tau_meshes(beta=1, wmax=1e20, order=16)
