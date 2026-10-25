"""Independent finite-query precision enclosures and failure boundary checks."""

from __future__ import annotations

import json
import math
from dataclasses import FrozenInstanceError
from fractions import Fraction
from typing import Any

import mpmath as mp
import numpy as np
import pytest
from mpmath.ctx_iv import MPIntervalContext

from kernelgecp import FermionicBasis, FermionicFit
from kernelgecp.precision import (
    IntervalArray,
    _certificate_hash,
    _directed,
    _endpoint,
    _exact,
    _exact_dot_error,
    _IntervalLU,
    _Uncertified,
    _upper_abs,
    certify_fit_precision,
)


@pytest.fixture(scope="module")
def basis() -> FermionicBasis:
    return FermionicBasis(3.25, 2.0, tol=1e-5, grid_order=6, precision="auto")


def _kernel(time: Any, frequency: Any, beta: Any) -> Any:
    return mp.exp(-time * frequency / beta) / (1 + mp.exp(-frequency))


def _oracle(
    basis: FermionicBasis,
    training: np.ndarray,
    data: np.ndarray,
    times: np.ndarray,
    modes: np.ndarray,
) -> tuple[Any, Any, Any, Any, Any]:
    beta = mp.mpf(basis.beta)
    frequencies = [mp.mpf(float(value)) for value in basis.omega_nodes]
    nodes = [mp.mpf(float(value)) for value in basis.tau_nodes]
    core = mp.matrix([[_kernel(time, w, beta) for w in frequencies] for time in nodes])

    def tau_rows(values: np.ndarray) -> Any:
        return mp.matrix(
            [
                list(
                    mp.lu_solve(
                        core.T,
                        mp.matrix(
                            [_kernel(mp.mpf(float(t)), w, beta) for w in frequencies]
                        ),
                    )
                )
                for t in values
            ]
        )

    design = tau_rows(training)
    right = mp.matrix([[mp.mpc(complex(value)) for value in row] for row in data])
    columns = [
        mp.qr_solve(design, right[:, column])[0] for column in range(data.shape[1])
    ]
    coefficients = mp.matrix(
        [[column[row] for column in columns] for row in range(basis.rank)]
    )
    time_design = tau_rows(times)
    frequency_design = mp.matrix(
        [
            list(
                mp.lu_solve(
                    core.T,
                    mp.matrix(
                        [
                            -beta / (mp.j * mp.pi * (2 * int(n) + 1) - w)
                            for w in frequencies
                        ]
                    ),
                )
            )
            for n in modes
        ]
    )
    return (
        design,
        coefficients,
        time_design * coefficients,
        frequency_design * coefficients,
        core,
    )


def _contains(
    enclosure: IntervalArray, exact: Any, shape: tuple[int, ...], axis: int
) -> None:
    for name, component in (("real", mp.re), ("imag", mp.im)):
        lower = np.moveaxis(getattr(enclosure, f"{name}_lower"), axis, 0).reshape(shape)
        upper = np.moveaxis(getattr(enclosure, f"{name}_upper"), axis, 0).reshape(shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                assert mp.mpf(float(lower[i, j])) <= component(exact[i, j])
                assert component(exact[i, j]) <= mp.mpf(float(upper[i, j]))


@pytest.mark.parametrize("axis", [0, 1, -1])
def test_public_complex_matrix_fit_encloses_independent_high_precision_qr(
    basis: FermionicBasis, axis: int
) -> None:
    training = np.linspace(0, basis.beta, 2 * basis.rank + 3)
    matrix = np.array([[1, 0.2 + 0.35j], [0.2 - 0.35j, 0.7]])
    values0 = (np.exp(-0.3 * training) - 0.2 * np.exp(-1.3 * training))[
        :, None, None
    ] * matrix
    values = np.moveaxis(values0, 0, axis)
    fit = basis.fit_tau(values, training, axis=axis, rcond=1e-14)
    times = np.array([0.1, 0.7, 2.3])
    modes = np.array([-7, -1, 0, 4])
    certificate = certify_fit_precision(
        fit, values, training, query_tau=times, matsubara_indices=modes
    )
    assert certificate.certified, certificate.termination_reason
    assert certificate.coefficients is not None
    assert certificate.training_design is not None
    assert certificate.time is not None and certificate.matsubara is not None
    assert certificate.core is not None
    with mp.workdps(90):
        design, coefficients, ideal_time, ideal_frequency, core = _oracle(
            basis, training, values0.reshape(training.size, -1), times, modes
        )
        _contains(
            certificate.training_design.ideal, design, (training.size, basis.rank), 0
        )
        _contains(
            certificate.coefficients.ideal, coefficients, (basis.rank, 4), fit.axis
        )
        _contains(
            certificate.time.ideal_prediction, ideal_time, (times.size, 4), fit.axis
        )
        _contains(
            certificate.matsubara.ideal_prediction,
            ideal_frequency,
            (modes.size, 4),
            fit.axis,
        )
        _contains(certificate.core.ideal, core, (basis.rank, basis.rank), 0)
        actual_coeff = np.moveaxis(fit.coefficients, fit.axis, 0).reshape(
            basis.rank, -1
        )
        coefficient_bounds = np.moveaxis(
            certificate.coefficients.actual_error_bound, fit.axis, 0
        ).reshape(basis.rank, -1)
        for i in range(basis.rank):
            for j in range(4):
                assert abs(
                    mp.mpc(complex(actual_coeff[i, j])) - coefficients[i, j]
                ) <= mp.mpf(float(coefficient_bounds[i, j]))
        for result, actual, oracle in (
            (certificate.time, fit.evaluate_tau(times), ideal_time),
            (certificate.matsubara, fit.evaluate_matsubara(modes), ideal_frequency),
        ):
            actual = np.moveaxis(actual, fit.axis, 0).reshape(oracle.rows, 4)
            direct = np.moveaxis(result.direct_error_bound, fit.axis, 0).reshape(
                oracle.rows, 4
            )
            total = np.moveaxis(result.total_error_bound, fit.axis, 0).reshape(
                oracle.rows, 4
            )
            for i in range(oracle.rows):
                for j in range(4):
                    error = abs(mp.mpc(complex(actual[i, j])) - oracle[i, j])
                    assert error <= mp.mpf(float(direct[i, j]))
                    assert error <= mp.mpf(float(total[i, j]))
    assert certificate.coefficients.max_error < 1e-8
    assert certificate.time.max_error < 1e-8
    assert certificate.matsubara.max_error < 1e-8
    assert certificate.core.normalization_error_bound.max() > 0
    json.dumps(certificate.to_dict(), allow_nan=False)


def test_interval_lu_cramer_and_underresolved_pivot() -> None:
    ctx = MPIntervalContext()
    ctx.prec = 256
    a, b, c, d = Fraction(2), Fraction(1, 3), Fraction(-1, 5), Fraction(4, 7)
    rhs = [Fraction(1, 11), Fraction(-2, 9)]

    def convert(q: Fraction) -> Any:
        return ctx.mpf(q.numerator) / ctx.mpf(q.denominator)

    lu = _IntervalLU([[convert(a), convert(b)], [convert(c), convert(d)]], ctx, "test")
    solution = lu.solve([convert(value) for value in rhs])
    determinant = a * d - b * c
    expected = [
        (rhs[0] * d - b * rhs[1]) / determinant,
        (a * rhs[1] - rhs[0] * c) / determinant,
    ]
    for enclosed, exact in zip(solution, expected, strict=True):
        assert enclosed.a <= convert(exact) <= enclosed.b
    upper = float(1 + 2**-27)
    lower = float(1 - 2**-27)
    for bits in (53, 256):
        current = MPIntervalContext()
        current.prec = bits
        matrix = [
            [_exact(current, x) for x in row] for row in [[1.0, upper], [lower, 1.0]]
        ]
        if bits == 53:
            with pytest.raises(_Uncertified, match="pivot_contains_zero"):
                _IntervalLU(matrix, current, "small_precision")
        else:
            assert _IntervalLU(matrix, current, "resolved").solve(
                [current.mpf(1), current.mpf(2)]
            )


def test_outward_exports_cover_exact_adjacent_endpoints_subnormals_and_maximum() -> (
    None
):
    ctx = MPIntervalContext()
    ctx.prec = 256
    for value in [
        0.0,
        1.0,
        np.nextafter(1.0, 2.0),
        np.nextafter(0.0, 1.0),
        np.finfo(float).max,
    ]:
        assert _endpoint(_exact(ctx, value), True) == value
        assert _endpoint(_exact(ctx, value), False) == value
    smallest = math.nextafter(0.0, math.inf)
    tiny = ctx.mpf(2) ** -2000
    assert _upper_abs(tiny) == smallest
    assert _upper_abs(-tiny) == smallest
    assert _endpoint(-tiny, False) == -smallest
    assert _endpoint(-tiny, True) == 0.0
    assert _upper_abs(ctx.mpf(2) ** -(10**20)) == smallest
    for value in [ctx.inf, ctx.nan, ctx.mpf(2) ** 2048]:
        with pytest.raises(_Uncertified):
            _upper_abs(value)
    fraction = Fraction(1) + Fraction(1, 2**60)
    assert _directed(fraction, True) == np.nextafter(1.0, 2.0)
    assert _directed(fraction, False) == 1.0


def test_exact_binary_dot_error_keeps_cancellation_and_underflow() -> None:
    ctx = MPIntervalContext()
    ctx.prec = 256
    design = np.array([[1e100 + 2j, 1e100 - 2j, 1 + 3j]])
    coefficients = np.array([[1 + 1j], [-1 - 1j], [0.1 - 0.2j]])
    actual = design @ coefficients
    bound = _exact_dot_error(design, coefficients, actual, ctx)
    with mp.workdps(170):
        true_dot = sum(
            mp.mpc(complex(design[0, i])) * mp.mpc(complex(coefficients[i, 0]))
            for i in range(3)
        )
        assert (
            0
            < abs(true_dot - mp.mpc(complex(actual[0, 0])))
            <= mp.mpf(float(bound[0, 0]))
        )
    tiny = np.array([[math.nextafter(0.0, math.inf)]])
    with np.errstate(under="ignore"):
        result = tiny @ np.array([[0.5]])
    assert result[0, 0] == 0
    assert _exact_dot_error(tiny, np.array([[0.5]]), result, ctx)[0, 0] == tiny[0, 0]


def test_certificates_are_owned_immutable_content_bound_and_repeatable(
    basis: FermionicBasis,
) -> None:
    training = np.linspace(0, basis.beta, basis.rank + 5)
    values = np.exp(-training)
    fit = basis.fit_tau(values, training)
    times = np.array([0.2, 0.8])
    first = certify_fit_precision(fit, values, training, query_tau=times)
    second = certify_fit_precision(fit, values, training, query_tau=times)
    assert first.certified and first.to_dict() == second.to_dict()
    assert first.time is not None and first.coefficients is not None
    before = first.to_dict()
    values[:] = 0
    training[:] = 0
    times[:] = 0
    assert first.to_dict() == before
    for array in (
        first.time.direct_error_bound,
        first.time.ideal_prediction.real_lower,
        first.coefficients.actual_error_bound,
        first.coefficients.ideal.imag_upper,
    ):
        assert not array.flags.writeable
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        first.certified = False
    changed = certify_fit_precision(
        fit, values, np.linspace(0, basis.beta, values.size), query_tau=[0.2, 0.8]
    )
    assert changed.content_sha256 != first.content_sha256
    assert _certificate_hash(np.array([1 + 2j])) != _certificate_hash(
        np.array([1 - 2j])
    )
    assert _certificate_hash(np.array([1 + 2j])) == _certificate_hash(
        np.array([1 + 2j], dtype=">c16")
    )


def test_unrelated_and_truncated_fits_can_have_valid_large_discrepancy(
    basis: FermionicBasis,
) -> None:
    times = np.linspace(0, basis.beta, basis.rank + 5)
    values = np.exp(-times)
    fits = [
        FermionicFit(basis, np.zeros(basis.rank)),
        basis.fit_tau(values, times, rcond=2),
    ]
    for fit in fits:
        certificate = certify_fit_precision(fit, values, times, query_tau=[0.1, 0.5])
        assert certificate.certified
        assert (
            certificate.coefficients is not None
            and certificate.coefficients.max_error > 0.5
        )
        assert certificate.time is not None and certificate.time.max_error > 0.5
        assert certificate.physical_truth_certificate is False
        assert certificate.continuous_certificate is False


def test_rank_work_and_singular_normal_system_are_explicit_inconclusive(
    basis: FermionicBasis,
) -> None:
    fit = FermionicFit(basis, np.ones(basis.rank))
    training = np.linspace(0, basis.beta, basis.rank + 2)
    values = np.ones(training.size)
    for kwargs, reason in (
        ({"max_rank": 1}, "rank_budget"),
        ({"max_work": 1}, "work_budget"),
    ):
        result = certify_fit_precision(fit, values, training, query_tau=[0.5], **kwargs)
        assert not result.certified and result.termination_reason == reason
        assert result.coefficients is None and result.time is None
    under = certify_fit_precision(fit, [1.0], [0.0], query_tau=[0.5])
    assert (
        not under.certified
        and under.termination_reason == "underdetermined_training_design"
    )
    singular = certify_fit_precision(
        fit, values, np.zeros(training.size), query_tau=[0.5]
    )
    assert (
        not singular.certified
        and singular.termination_reason == "normal_pivot_contains_zero"
    )
    assert singular.coefficients is None


def test_owned_interval_context_and_numpy_traps_are_preserved(
    basis: FermionicBasis,
) -> None:
    times = np.linspace(0, basis.beta, basis.rank + 3)
    values = np.exp(-times) + 0.3j
    fit = basis.fit_tau(values, times)
    old_mp, old_iv = mp.mp.prec, mp.iv.prec
    old_numpy = np.geterr()
    try:
        mp.mp.prec, mp.iv.prec = 71, 83
        with np.errstate(all="raise"):
            result = certify_fit_precision(
                fit,
                values,
                times,
                query_tau=[0.0, 0.1, basis.beta],
                matsubara_indices=[-1, 0],
            )
            assert result.certified, result.termination_reason
            assert mp.mp.prec == 71 and mp.iv.prec == 83
            assert all(mode == "raise" for mode in np.geterr().values())
        assert np.geterr() == old_numpy
    finally:
        mp.mp.prec, mp.iv.prec = old_mp, old_iv


@pytest.mark.parametrize(
    "changes",
    [
        {"precision_bits": True},
        {"precision_bits": 52},
        {"precision_bits": 4097},
        {"max_rank": 0},
        {"max_work": 0},
        {"max_work": 2**31},
        {"query_tau": [-0.1]},
        {"query_tau": [1j]},
        {"query_tau": [True]},
        {"matsubara_indices": [True]},
        {"matsubara_indices": [1.0]},
        {"matsubara_indices": [2**51]},
        {"query_tau": None},
    ],
)
def test_invalid_controls_and_queries_are_rejected(
    basis: FermionicBasis, changes: dict[str, Any]
) -> None:
    fit = FermionicFit(basis, np.ones(basis.rank))
    kwargs = {"query_tau": [0.5], **changes}
    with pytest.raises(ValueError):
        certify_fit_precision(fit, np.ones(basis.rank), **kwargs)


@pytest.mark.parametrize(
    "values", [True, np.ones((2, 3)), np.full(4, np.nan), [], ["1"]]
)
def test_invalid_observations_are_rejected(basis: FermionicBasis, values: Any) -> None:
    fit = FermionicFit(basis, np.ones(basis.rank))
    with pytest.raises(ValueError):
        certify_fit_precision(fit, values, query_tau=[0.5])


def test_frequency_only_and_default_training_grid(basis: FermionicBasis) -> None:
    values = np.sin(basis.tau_nodes) + 0.2j
    fit = basis.fit_tau(values)
    result = certify_fit_precision(fit, values, matsubara_indices=[-2, -1, 0, 1])
    assert result.certified and result.time is None and result.matsubara is not None
    assert (
        result.training_design is not None
        and result.training_design.max_solve_residual >= 0
    )
    assert result.observation_axis == 0


def test_uncertain_three_by_three_late_row_swap_encloses_all_exact_vertices() -> None:
    from itertools import product

    ctx = MPIntervalContext()
    ctx.prec = 256
    centers = [Fraction(value) for value in [4, 2, 1, 1, 0.5, 3, 2, 5, 4]]
    radius = Fraction(1, 1024)

    def interval(value: Fraction) -> Any:
        lower, upper = value - radius, value + radius
        return ctx.mpf(
            [
                ctx.mpf(lower.numerator) / ctx.mpf(lower.denominator),
                ctx.mpf(upper.numerator) / ctx.mpf(upper.denominator),
            ]
        )

    def determinant(matrix: list[list[Fraction]]) -> Fraction:
        a, b, c = matrix
        return (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )

    matrix = [
        [interval(value) for value in centers[row * 3 : row * 3 + 3]]
        for row in range(3)
    ]
    factor = _IntervalLU(matrix, ctx, "uncertain_three_by_three")
    assert factor.permutation == [0, 2, 1]
    rhs = [Fraction(1), Fraction(-2), Fraction(3)]
    enclosed = factor.solve([ctx.mpf(int(value)) for value in rhs])
    for signs in product((-1, 1), repeat=9):
        entries = [
            center + sign * radius for center, sign in zip(centers, signs, strict=True)
        ]
        vertex = [entries[row * 3 : row * 3 + 3] for row in range(3)]
        divisor = determinant(vertex)
        assert divisor != 0
        for column in range(3):
            replaced = [row.copy() for row in vertex]
            for row in range(3):
                replaced[row][column] = rhs[row]
            exact = determinant(replaced) / divisor
            expected = ctx.mpf(exact.numerator) / ctx.mpf(exact.denominator)
            assert enclosed[column].a <= expected <= enclosed[column].b
