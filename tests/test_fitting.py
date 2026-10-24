from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from scipy.linalg import lstsq

from kernelgecp.fitting import LeastSquaresPlan
from kernelgecp.types import NumericalBreakdownError


@pytest.mark.parametrize("axis", [0, 1, 2, -1])
@pytest.mark.parametrize("complex_design", [False, True])
def test_batched_fits_match_independent_qr(axis: int, complex_design: bool) -> None:
    generator = np.random.default_rng(81)
    design = generator.normal(size=(13, 5))
    if complex_design:
        design = design + 1j * generator.normal(size=design.shape)
    values = generator.normal(size=(13, 2, 3, 4)) + 1j * generator.normal(
        size=(13, 2, 3, 4)
    )
    supplied = np.moveaxis(values, 0, axis)
    original = supplied.copy()
    plan = LeastSquaresPlan(design)
    coefficients = plan.fit(supplied, axis=axis)
    expected, _, rank, _ = lstsq(
        design,
        values.reshape(13, -1),
        cond=plan.diagnostics.rcond,
        lapack_driver="gelsy",
    )
    assert rank == 5
    assert coefficients.shape == np.moveaxis(np.empty((5, 2, 3, 4)), 0, axis).shape
    np.testing.assert_allclose(
        np.moveaxis(coefficients, axis, 0).reshape(5, -1),
        expected,
        rtol=2e-13,
        atol=2e-13,
    )
    expected_residual = (design @ expected - values.reshape(13, -1)).reshape(
        values.shape
    )
    np.testing.assert_allclose(
        plan.residual(supplied, coefficients, axis=axis),
        np.moveaxis(expected_residual, 0, axis),
        rtol=2e-13,
        atol=2e-13,
    )
    np.testing.assert_array_equal(supplied, original)


def test_vector_and_underdetermined_minimum_norm_fit() -> None:
    design = np.array([[1.0, 2.0, -1.0], [2.0, -1.0, 3.0]])
    values = np.array([2.0, -3.0])
    plan = LeastSquaresPlan(design)
    expected, _, _, _ = lstsq(design, values, lapack_driver="gelsy")
    coefficients = plan.fit(values)
    assert coefficients.shape == (3,)
    assert coefficients.dtype == np.float64
    np.testing.assert_allclose(coefficients, expected, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(plan.residual(values, coefficients), 0, atol=2e-14)


def test_default_cutoff_and_binary64_promotion() -> None:
    plan = LeastSquaresPlan(np.eye(3, dtype=np.float32))
    assert plan.diagnostics.rcond == np.finfo(np.float64).eps * 3
    assert plan.diagnostics.design_dtype == "float64"
    assert plan.design.dtype == np.float64
    assert plan.fit(np.array([1, 2, 3], dtype=np.float32)).dtype == np.float64
    complex_values = np.array([1, 2j, 3], dtype=np.complex64)
    assert plan.fit(complex_values).dtype == np.complex128


def test_rank_truncation_is_explicit_and_strict() -> None:
    plan = LeastSquaresPlan(np.diag([1.0, 0.25, 0.0]), rcond=0.25)
    np.testing.assert_array_equal(plan.fit([2, 3, 4]), [2, 0, 0])
    diagnostic = plan.diagnostics
    assert diagnostic.effective_rank == 1
    assert diagnostic.rcond == 0.25
    assert diagnostic.rank_threshold == 0.25
    assert diagnostic.singular_values == (1.0, 0.25, 0.0)
    assert diagnostic.conditioning_status == "rank_truncated"
    assert diagnostic.condition_number is None
    assert diagnostic.retained_condition_number == 1.0
    assert diagnostic.precision_bits == 53
    assert diagnostic.certification_status == "not_certified"
    assert "supplied_training_values" in diagnostic.residual_scope
    with pytest.raises(FrozenInstanceError):
        diagnostic.effective_rank = 3


def test_zero_rank_and_empty_batches_are_supported() -> None:
    plan = LeastSquaresPlan(np.zeros((4, 3)))
    values = np.arange(8).reshape(2, 4) + 1j
    fitted = plan.fit(values, axis=1)
    assert fitted.shape == (2, 3)
    assert fitted.dtype == np.complex128
    np.testing.assert_array_equal(fitted, 0)
    np.testing.assert_array_equal(plan.residual(values, fitted, axis=1), -values)
    assert plan.diagnostics.effective_rank == 0
    assert plan.diagnostics.conditioning_status == "zero_rank"
    empty = np.empty((2, 4, 0))
    assert plan.fit(empty, axis=1).shape == (2, 3, 0)
    assert plan.residual(empty, np.empty((2, 3, 0)), axis=1).shape == empty.shape


def test_design_is_copied_and_immutable_and_fits_reuse_one_svd(monkeypatch) -> None:
    design = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
    original = design.copy()
    calls = 0
    original_svd = np.linalg.svd

    def counted_svd(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original_svd(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "svd", counted_svd)
    plan = LeastSquaresPlan(design)
    design[:] = 100
    np.testing.assert_array_equal(plan.design, original)
    assert not plan.design.flags.writeable
    with pytest.raises(ValueError):
        plan.design[0, 0] = 99
    with pytest.raises(ValueError):
        plan.design.setflags(write=True)
    for scale in [1, 3, 0.2]:
        values = original @ np.array([2.0, -1.0]) * scale
        fitted = plan.fit(values)
        np.testing.assert_allclose(fitted, np.array([2, -1]) * scale, atol=1e-14)
        np.testing.assert_allclose(plan.residual(values, fitted), 0, atol=1e-14)
    assert calls == 1


@pytest.mark.parametrize(
    "design",
    [
        [1, 2],
        np.ones((1, 1, 1)),
        np.empty((0, 2)),
        np.empty((2, 0)),
        [[np.inf]],
        [[np.nan + 1j]],
        [[True]],
        [[1.0, True]],
        [np.array([1.0]), np.array([True])],
        [["1"]],
        np.array([[1]], dtype=object),
    ],
)
def test_design_validation(design) -> None:
    with pytest.raises(ValueError):
        LeastSquaresPlan(design)


@pytest.mark.parametrize("rcond", [-1, np.inf, np.nan, True, 1j, "0.1", 10**1_000])
def test_rcond_validation(rcond) -> None:
    with pytest.raises(ValueError):
        LeastSquaresPlan(np.eye(2), rcond=rcond)


@pytest.mark.parametrize("axis", [True, np.bool_(False), 0.5, 1, -2])
def test_fit_axis_validation(axis) -> None:
    with pytest.raises(ValueError):
        LeastSquaresPlan(np.eye(2)).fit([1, 2], axis=axis)


@pytest.mark.parametrize(
    "values", [1.0, [1.0], [np.inf, 1], [np.nan, 2j], [True, False], [1, True]]
)
def test_fit_value_validation(values) -> None:
    with pytest.raises(ValueError):
        LeastSquaresPlan(np.eye(2)).fit(values)


@pytest.mark.parametrize("coefficients", [[1], [[1, 2]], [np.nan, 2], [True, False]])
def test_residual_coefficient_validation(coefficients) -> None:
    with pytest.raises(ValueError):
        LeastSquaresPlan(np.eye(2)).residual([1, 2], coefficients)


def test_nonfinite_computation_is_not_a_successful_fit() -> None:
    plan = LeastSquaresPlan([[1e-300]])
    with pytest.raises(NumericalBreakdownError, match="coefficients are nonfinite"):
        plan.fit([1e300])
    with pytest.raises(NumericalBreakdownError, match="residual is nonfinite"):
        LeastSquaresPlan([[1e300]]).residual([0], [1e300])


def test_zero_rcond_retains_positive_subnormal_singular_values() -> None:
    plan = LeastSquaresPlan(np.diag([1.0, 1e-310]), rcond=0.0)
    assert plan.diagnostics.effective_rank == 2
    assert plan.diagnostics.condition_number is None
    with np.errstate(all="raise"):
        np.testing.assert_allclose(plan.fit([1, 1e-310]), [1, 1], atol=1e-14)


def test_local_arithmetic_context_preserves_caller_traps() -> None:
    plan = LeastSquaresPlan([[1e-100]])
    values = [1e-300]
    with np.errstate(all="raise"):
        previous = np.geterr()
        coefficients = plan.fit(values)
        plan.residual(values, coefficients)
        assert np.geterr() == previous
