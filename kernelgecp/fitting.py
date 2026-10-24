"""Reusable binary64 least-squares fits with explicit rank truncation.

The factorization concerns only the supplied, rounded design matrix. Neither
the fitted coefficients nor residuals certify a continuous-domain error bound.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from numbers import Real
from typing import cast

import numpy as np
import numpy.typing as npt

from .types import NumericalBreakdownError

type NumericArray = npt.NDArray[np.float64] | npt.NDArray[np.complex128]


@dataclass(frozen=True, slots=True)
class LeastSquaresDiagnostics:
    """Immutable diagnostics of one cached thin SVD.

    ``condition_number`` describes the untruncated design and is unavailable
    when a singular value is zero or its ratio overflows binary64.
    ``retained_condition_number`` describes only the retained singular vectors.
    The rank rule is strictly ``singular_value > rank_threshold``. A null
    threshold means the finite rcond times the largest singular value exceeds
    binary64 range; all singular values are then discarded.
    """

    design_shape: tuple[int, int]
    design_dtype: str
    singular_values: tuple[float, ...]
    effective_rank: int
    rcond: float
    rank_threshold: float | None
    condition_number: float | None
    retained_condition_number: float | None
    conditioning_status: str
    factorization: str = "cached_thin_svd"
    arithmetic: str = "numpy_float64_or_complex128"
    precision_bits: int = 53
    converged: bool = True
    stop_reason: str = "svd_factorization_completed"
    convergence_scope: str = "numerical_factorization_of_supplied_design"
    solve_scope: str = "truncated_svd_minimum_norm_fit_on_supplied_samples"
    residual_scope: str = "design_times_coefficients_minus_supplied_training_values"
    certification_status: str = "not_certified"
    certification_scope: str = "no_analytic_or_rounding_certificate"


def _contains_boolean_sequence(values: object) -> bool:
    if isinstance(values, (bool, np.bool_)):
        return True
    if isinstance(values, np.ndarray):
        return bool(values.dtype.kind == "b")
    if isinstance(values, (list, tuple)):
        return any(_contains_boolean_sequence(item) for item in values)
    return False


def _numeric_array(values: npt.ArrayLike, name: str) -> NumericArray:
    """Validate before conversion so booleans/strings are not silently coerced."""

    if _contains_boolean_sequence(values):
        raise ValueError(f"{name} must contain real or complex numbers, not booleans")
    try:
        array = np.asarray(values)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a regular numeric array") from error
    if array.dtype.kind not in "iufc":
        raise ValueError(f"{name} must contain real or complex numbers, not booleans")
    dtype = np.complex128 if np.iscomplexobj(array) else np.float64
    with np.errstate(all="ignore"):
        result = np.asarray(array, dtype=dtype)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite binary64 values")
    return result


def _immutable_copy(array: NumericArray) -> NumericArray:
    """Copy into an immutable buffer, including protection against setflags."""

    return cast(
        NumericArray,
        np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape),
    )


def _normalize_axis(axis: int, ndim: int) -> int:
    if isinstance(axis, (bool, np.bool_)):
        raise ValueError("axis must be an integer, not a boolean")
    try:
        integer_axis = operator.index(axis)
    except TypeError as error:
        raise ValueError("axis must be an integer") from error
    if integer_axis < -ndim or integer_axis >= ndim:
        raise ValueError(f"axis is out of range for a {ndim}-dimensional array")
    return integer_axis % ndim


def _condition_ratio(largest: float, smallest: float) -> float | None:
    if smallest == 0:
        return None
    ratio = largest / smallest
    return ratio if np.isfinite(ratio) else None


class LeastSquaresPlan:
    """Factor a fixed design once and fit any number of independent batches.

    ``design`` has shape ``(observations, coefficients)``. ``fit(values, axis=a)``
    replaces the observation axis of ``values`` by the coefficient axis and
    preserves every other dimension. Real and complex inputs are supported;
    arithmetic uses float64 or complex128 regardless of input precision.

    A thin SVD is cached at construction. Singular values at or below
    ``rcond * largest_singular_value`` are discarded. The default ``rcond`` is
    machine epsilon times the larger design dimension. Rank-deficient and
    underdetermined systems return the truncated minimum-norm solution.
    """

    __slots__ = ("_design", "_diagnostics", "_singular_values", "_u", "_vh")

    def __init__(self, design: npt.ArrayLike, *, rcond: float | None = None) -> None:
        matrix = _numeric_array(design, "design")
        if matrix.ndim != 2 or min(matrix.shape) == 0:
            raise ValueError("design must be a nonempty two-dimensional matrix")
        if rcond is None:
            relative_cutoff = float(np.finfo(np.float64).eps * max(matrix.shape))
        else:
            if isinstance(rcond, (bool, np.bool_)) or not isinstance(rcond, Real):
                raise ValueError("rcond must be a finite nonnegative real number")
            try:
                relative_cutoff = float(rcond)
            except OverflowError as error:
                raise ValueError("rcond must be finite in binary64") from error
            if not np.isfinite(relative_cutoff) or relative_cutoff < 0:
                raise ValueError("rcond must be a finite nonnegative real number")
        self._design = _immutable_copy(matrix)
        u, singular_values, vh = np.linalg.svd(self._design, full_matrices=False)
        if not all(np.all(np.isfinite(part)) for part in (u, singular_values, vh)):
            raise NumericalBreakdownError("SVD factors are nonfinite in binary64")
        self._u = _immutable_copy(u)
        self._singular_values = _immutable_copy(
            np.asarray(singular_values, dtype=np.float64)
        )
        self._vh = _immutable_copy(vh)
        largest = float(singular_values[0])
        threshold = relative_cutoff * largest
        rank = int(np.count_nonzero(singular_values > threshold))
        smallest = float(singular_values[-1])
        self._diagnostics = LeastSquaresDiagnostics(
            design_shape=(matrix.shape[0], matrix.shape[1]),
            design_dtype=str(matrix.dtype),
            singular_values=tuple(float(value) for value in singular_values),
            effective_rank=rank,
            rcond=relative_cutoff,
            rank_threshold=threshold if np.isfinite(threshold) else None,
            condition_number=_condition_ratio(largest, smallest),
            retained_condition_number=(
                _condition_ratio(largest, float(singular_values[rank - 1]))
                if rank
                else None
            ),
            conditioning_status=(
                "zero_rank"
                if rank == 0
                else "rank_truncated"
                if rank < singular_values.size
                else "full_rank"
            ),
        )

    @property
    def design(self) -> NumericArray:
        """Return the plan-owned, immutable binary64 design matrix."""

        return self._design

    @property
    def diagnostics(self) -> LeastSquaresDiagnostics:
        """Return immutable factorization, conditioning and arithmetic metadata."""

        return self._diagnostics

    def _observations(
        self, values: npt.ArrayLike, axis: int
    ) -> tuple[NumericArray, int]:
        array = _numeric_array(values, "values")
        if array.ndim == 0:
            raise ValueError("values must have an observation axis")
        normalized_axis = _normalize_axis(axis, array.ndim)
        if array.shape[normalized_axis] != self._design.shape[0]:
            raise ValueError("values observation axis must match the design rows")
        return array, normalized_axis

    def fit(self, values: npt.ArrayLike, *, axis: int = 0) -> NumericArray:
        """Fit independently along ``axis`` without refactoring the design."""

        array, normalized_axis = self._observations(values, axis)
        moved = np.moveaxis(array, normalized_axis, 0)
        right_hand_side = moved.reshape(self._design.shape[0], -1)
        rank = self._diagnostics.effective_rank
        with np.errstate(all="ignore"):
            if rank:
                projected = self._u[:, :rank].conj().T @ right_hand_side
                projected = projected / self._singular_values[:rank, None]
                coefficients = self._vh[:rank, :].conj().T @ projected
            else:
                coefficients = np.zeros(
                    (self._design.shape[1], right_hand_side.shape[1]),
                    dtype=np.result_type(self._design.dtype, array.dtype),
                )
        if not np.all(np.isfinite(coefficients)):
            raise NumericalBreakdownError(
                "least-squares coefficients are nonfinite in binary64"
            )
        shaped = coefficients.reshape((self._design.shape[1], *moved.shape[1:]))
        return np.moveaxis(shaped, 0, normalized_axis)

    def residual(
        self,
        values: npt.ArrayLike,
        coefficients: npt.ArrayLike,
        *,
        axis: int = 0,
    ) -> NumericArray:
        """Return ``design @ coefficients - values`` on the supplied samples.

        The returned array has the original ``values`` shape. This residual
        does not estimate error between observations or at held-out points.
        """

        array, normalized_axis = self._observations(values, axis)
        coefficient_array = _numeric_array(coefficients, "coefficients")
        expected_shape = list(array.shape)
        expected_shape[normalized_axis] = self._design.shape[1]
        if coefficient_array.shape != tuple(expected_shape):
            raise ValueError("coefficients shape must match the fitted values shape")
        moved_values = np.moveaxis(array, normalized_axis, 0)
        moved_coefficients = np.moveaxis(coefficient_array, normalized_axis, 0)
        with np.errstate(all="ignore"):
            residual = self._design @ moved_coefficients.reshape(
                self._design.shape[1], -1
            ) - moved_values.reshape(self._design.shape[0], -1)
        if not np.all(np.isfinite(residual)):
            raise NumericalBreakdownError("training residual is nonfinite in binary64")
        return np.moveaxis(residual.reshape(moved_values.shape), 0, normalized_axis)
