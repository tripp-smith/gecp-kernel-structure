"""Independent physical identities and precision boundaries for frozen bases."""

from __future__ import annotations

from typing import Any

import mpmath as mp
import numpy as np
import pytest

import kernelgecp._basis_numerics as numerics
from kernelgecp._basis_numerics import BasisNumerics
from kernelgecp.types import NumericalBreakdownError


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_physical_pole_identity_nonunit_beta_and_both_frequency_signs(
    precision_bits: int,
) -> None:
    beta = 3.25
    poles = np.array([-0.8, 0.4])
    weights = np.array([0.35, 0.65])
    nodes = np.array([0.0, beta])
    basis = BasisNumerics(beta, 4.0, nodes, beta * poles, precision_bits=precision_bits)

    def physical(times: np.ndarray) -> np.ndarray:
        return (
            -(np.exp(-times[:, None] * poles) / (1 + np.exp(-beta * poles))) @ weights
        )

    coefficients = physical(nodes)
    target = np.linspace(0, beta, 19)
    np.testing.assert_allclose(
        basis.design_tau(target) @ coefficients,
        physical(target),
        rtol=2e-14,
        atol=2e-15,
    )
    modes = np.array([-4, -1, 0, 3], dtype=np.int64)
    physical_frequencies = 1j * np.pi * (2 * modes + 1) / beta
    exact = (1 / (physical_frequencies[:, None] - poles)) @ weights
    np.testing.assert_allclose(
        basis.design_matsubara(modes) @ coefficients, exact, rtol=2e-14, atol=2e-15
    )


def test_reference_nodal_stability_on_ill_conditioned_cutoff1000_fixture() -> None:
    beta = 20.0
    basis = BasisNumerics(
        beta,
        1000.0,
        beta * np.arange(6) * 0.001,
        np.arange(900.0, 1001.0, 20.0),
        precision_bits=128,
    )
    assert basis.exported_core_condition is not None
    assert basis.exported_core_condition > 1e10
    assert basis.core_condition == basis.exported_core_condition
    assert 0 < basis.nodal_error < 1e-24
    observed = np.max(np.abs(basis.design_tau(basis.tau_nodes) - np.eye(6)))
    assert basis.nodal_error == observed


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_factorization_is_cached_for_repeated_real_and_complex_evaluations(
    monkeypatch: pytest.MonkeyPatch, precision_bits: int
) -> None:
    calls = []
    if precision_bits == 53:
        factor = numerics.lu_factor

        def float_factor(*args: Any, **kwargs: Any) -> Any:
            calls.append(1)
            return factor(*args, **kwargs)

        monkeypatch.setattr(numerics, "lu_factor", float_factor)
    else:
        context_type = type(mp.mp)
        mp_factor = context_type.lu

        def reference_factor(context: Any, *args: Any, **kwargs: Any) -> Any:
            calls.append(1)
            return mp_factor(context, *args, **kwargs)

        monkeypatch.setattr(context_type, "lu", reference_factor)
    basis = BasisNumerics(
        2.0, 3.0, [0.0, 2.0], [-2.0, 2.0], precision_bits=precision_bits
    )
    for _ in range(2):
        basis.design_tau(np.linspace(0, 2, 7))
        basis.design_matsubara(np.arange(-3, 3))
    assert len(calls) == 1


def test_owned_reference_context_preserves_precision_and_results() -> None:
    with mp.workprec(71):
        basis = BasisNumerics(2.0, 3.0, [0.0, 2.0], [-2.0, 2.0], precision_bits=128)
        first = basis.design_matsubara(np.arange(-3, 3))
        assert mp.mp.prec == 71
        with mp.workprec(53):
            second = basis.design_matsubara(np.arange(-3, 3))
            assert mp.mp.prec == 53
        assert mp.mp.prec == 71
    np.testing.assert_array_equal(first, second)


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_nodes_are_copied_immutable_and_empty_queries_have_correct_shape(
    precision_bits: int,
) -> None:
    tau, omega = np.array([0.0, 2.0]), np.array([-2.0, 2.0])
    basis = BasisNumerics(2.0, 3.0, tau, omega, precision_bits=precision_bits)
    before = basis.design_tau([0.25, 1.5])
    tau[:] = 1.0
    omega[:] = 0.0
    np.testing.assert_array_equal(before, basis.design_tau([0.25, 1.5]))
    for nodes in (basis.tau_nodes, basis.omega_nodes):
        assert not nodes.flags.writeable
        with pytest.raises(ValueError):
            nodes.setflags(write=True)
    assert basis.design_tau([]).shape == (0, 2)
    assert basis.design_matsubara(np.array([], dtype=int)).shape == (0, 2)


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_kernel_tail_underflow_does_not_depend_on_numpy_caller_traps(
    precision_bits: int,
) -> None:
    with np.errstate(all="raise"):
        basis = BasisNumerics(
            20.0, 1000.0, [0.0], [1000.0], precision_bits=precision_bits
        )
        observed = basis.design_tau([0.0, 20.0])
        transformed = basis.design_matsubara(np.array([-1, 0]))
    assert observed[0, 0] == 1.0
    assert observed[1, 0] == 0.0
    assert np.all(np.isfinite(transformed))


@pytest.mark.parametrize(
    "changes",
    [
        {"tau_nodes": []},
        {"tau_nodes": [0.0]},
        {"tau_nodes": [[0.0, 2.0]]},
        {"tau_nodes": [0.0, float("nan")]},
        {"tau_nodes": [-1.0, 2.0]},
        {"tau_nodes": [0.0, 3.0]},
        {"omega_nodes": [-4.0, 2.0]},
        {"omega_nodes": [-2.0, float("inf")]},
        {"omega_nodes": [-2.0, 2j]},
        {"precision_bits": True},
        {"precision_bits": 52},
        {"precision_bits": 128.0},
        {"beta": 0.0},
        {"cutoff": 0.5},
    ],
)
def test_invalid_construction_is_rejected(changes: dict[str, Any]) -> None:
    controls = dict(
        beta=2.0,
        cutoff=3.0,
        tau_nodes=[0.0, 2.0],
        omega_nodes=[-2.0, 2.0],
        precision_bits=128,
    )
    with pytest.raises(ValueError):
        BasisNumerics(**(controls | changes))


@pytest.mark.parametrize("precision_bits", [53, 128])
def test_singular_selected_core_is_explicit_failure(precision_bits: int) -> None:
    with pytest.raises(NumericalBreakdownError, match="core LU failed"):
        BasisNumerics(2.0, 3.0, [0.0, 0.0], [-2.0, 2.0], precision_bits=precision_bits)


@pytest.mark.parametrize(
    "values", [0.5, [[0.5]], [float("inf")], [1j], [True], [-0.1], [3]]
)
def test_invalid_time_queries_are_rejected(values: Any) -> None:
    basis = BasisNumerics(2.0, 3.0, [0.0], [0.0], precision_bits=128)
    with pytest.raises(ValueError):
        basis.design_tau(values)


@pytest.mark.parametrize("values", [0, [[0]], [0.0], [float("nan")], [1j], [True]])
def test_invalid_matsubara_queries_are_rejected(values: Any) -> None:
    basis = BasisNumerics(2.0, 3.0, [0.0], [0.0], precision_bits=128)
    with pytest.raises(ValueError):
        basis.design_matsubara(values)
