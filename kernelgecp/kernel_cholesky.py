"""Column-only finite PSD-kernel approximation and shifted linear algebra.

The caller supplies a deterministic symmetric PSD kernel through its diagonal
and selected columns. Checks on requested columns are necessary consistency
checks, not a proof of the unobserved matrix's symmetry or positivity.
"""

from __future__ import annotations

import importlib
import math
import operator
from collections.abc import Callable
from dataclasses import dataclass
from numbers import Real
from typing import Any, cast

import numpy as np
import numpy.typing as npt
from scipy.linalg import cho_factor, cho_solve, qr, solve_triangular, svdvals
from scipy.linalg.lapack import dormqr
from scipy.sparse.linalg import LinearOperator

from .fitting import _immutable_copy, _numeric_array
from .types import (
    ConditioningDiagnostic,
    FloatArray,
    NumericalBreakdownError,
    NumericalMetadata,
)

type ColumnOracle = Callable[[int], npt.ArrayLike]
_EPS = np.finfo(np.float64).eps


def _real(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be finite in float64") from error
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError(
            f"{name} must be finite and {'positive' if positive else '>=0'}"
        )
    return result


def _count(value: object, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a nonnegative integer")
    try:
        result = operator.index(value)  # type: ignore[arg-type]
    except TypeError as error:
        raise ValueError(f"{name} must be a nonnegative integer") from error
    if result < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return result


def _array(values: npt.ArrayLike, name: str, *, ndim: int | None = None) -> FloatArray:
    result = _numeric_array(values, name)
    if np.iscomplexobj(result) or (ndim is not None and result.ndim != ndim):
        raise ValueError(f"{name} must be a real array of the required shape")
    return cast(FloatArray, result)


def _owned(values: npt.ArrayLike) -> FloatArray:
    return cast(FloatArray, _immutable_copy(np.asarray(values, dtype=np.float64)))


def _shift(values: npt.ArrayLike, size: int, *, positive: bool) -> FloatArray:
    result = _array(values, "shift")
    if result.ndim == 0:
        scalar = _real(float(result), "shift", positive=positive)
        result = np.full(size, scalar, dtype=np.float64)
    if result.shape != (size,) or np.any(result < 0):
        raise ValueError("shift must be a nonnegative scalar or length-n diagonal")
    if positive and np.any(result <= 0):
        raise ValueError("shifted solve requires a strictly positive diagonal")
    return _owned(result)


def _finite(values: FloatArray, scope: str) -> None:
    if not np.all(np.isfinite(values)):
        raise NumericalBreakdownError(f"nonfinite {scope}")


@dataclass(frozen=True, slots=True)
class KernelCholeskyDiagnostics:
    """Numerical finite-diagonal evidence; never a PSD or rounding certificate."""

    converged: bool
    stop_reason: str
    threshold: float
    roundoff_slack: float
    tie_rtol: float
    tie_counts: tuple[int, ...]
    residual_history: tuple[float, ...]
    clipped_negative_diagonals: int
    returned_diagonal_error: float
    root_sigma_min: float | None
    root_condition_number: float | None
    column_evaluations: int


@dataclass(frozen=True, slots=True)
class KernelCholeskyExtension:
    """Owned out-of-sample features and optional checked residual diagonal."""

    features: FloatArray
    residual_diagonal: FloatArray | None

    def __post_init__(self) -> None:
        features = _array(self.features, "features", ndim=2)
        object.__setattr__(self, "features", _owned(features))
        if self.residual_diagonal is not None:
            diagonal = _array(self.residual_diagonal, "residual_diagonal", ndim=1)
            if diagonal.shape != (features.shape[0],) or np.any(diagonal < 0):
                raise ValueError(
                    "residual_diagonal must be a nonnegative length-m vector"
                )
            object.__setattr__(self, "residual_diagonal", _owned(diagonal))


@dataclass(frozen=True, slots=True, init=False, eq=False)
class KernelCholeskyFactor:
    """Immutable float64 factor constructed with one column per selected pivot.

    ``column(j)`` returns the full column K[:, j]. Only selected columns are
    queried. The callback is not retained. Ties choose the first unselected
    positive index within ``tie_rtol * maximum_residual_diagonal``. Stopping
    measures a clipped finite residual diagonal, assuming the kernel is PSD.
    The existing dense ``pivoted_cholesky`` API is independent and unchanged.
    """

    factor: FloatArray
    diagonal: FloatArray
    residual_diagonal: FloatArray
    indices: tuple[int, ...]
    pivots: tuple[float, ...]
    diagnostics: KernelCholeskyDiagnostics

    def __init__(
        self,
        diagonal: npt.ArrayLike,
        column: ColumnOracle,
        *,
        max_rank: int = 256,
        atol: float = 1e-10,
        rtol: float = 0.0,
        tie_rtol: float = 32 * _EPS,
        psd_rtol: float = 128 * _EPS,
    ) -> None:
        diag = _owned(_array(diagonal, "diagonal", ndim=1))
        if np.any(diag < 0) or not callable(column):
            raise ValueError("nonnegative diagonal and callable column oracle required")
        limit = min(_count(max_rank, "max_rank"), diag.size)
        absolute, relative = _real(atol, "atol"), _real(rtol, "rtol")
        tie, psd = _real(tie_rtol, "tie_rtol"), _real(psd_rtol, "psd_rtol")
        if tie >= 1 or psd >= 1:
            raise ValueError("relative tie and PSD allowances must be smaller than one")
        scale = float(np.max(diag, initial=0))
        threshold, slack = absolute + relative * scale, psd * scale
        if not math.isfinite(threshold):
            raise ValueError("stopping threshold must be finite in float64")
        root = np.zeros((diag.size, limit), dtype=np.float64)
        residual = diag.copy()
        selected: list[int] = []
        pivots: list[float] = []
        ties: list[int] = []
        history: list[float] = []
        clipped = 0
        converged, stop = False, "max_rank"
        for rank in range(limit + 1):
            _finite(residual, "residual diagonal")
            if np.any(residual < -slack):
                raise NumericalBreakdownError("materially negative residual diagonal")
            clipped += int(np.count_nonzero(residual < 0))
            np.maximum(residual, 0, out=residual)
            maximum = float(np.max(residual, initial=0))
            history.append(maximum)
            if maximum <= threshold:
                converged, stop = True, "tolerance"
                break
            if rank == limit:
                break
            if maximum <= slack:
                stop = "roundoff_limit"
                break
            eligible = (residual > 0) & (maximum - residual <= tie * maximum)
            eligible[selected] = False
            candidates = np.flatnonzero(eligible)
            if not candidates.size:
                stop = "roundoff_limit"
                break
            pivot = int(candidates[0])
            values = _array(column(pivot), "oracle column", ndim=1)
            if values.shape != diag.shape:
                raise ValueError("oracle column must have shape (n,)")
            if abs(float(values[pivot]) - float(diag[pivot])) > slack:
                raise ValueError("oracle column disagrees with supplied diagonal")
            with np.errstate(all="ignore"):
                bounds = np.sqrt(diag) * math.sqrt(float(diag[pivot]))
                excess = np.abs(values) - bounds
            if np.any(excess > slack):
                raise NumericalBreakdownError("oracle column violates PSD entry bound")
            previous = root[:, :rank]
            with np.errstate(all="ignore"):
                reconstructed = previous @ previous[pivot]
                remainder = values - reconstructed
            _finite(remainder, "column elimination")
            if selected and np.any(np.abs(remainder[selected]) > slack):
                raise ValueError(
                    "oracle columns disagree on selected symmetric entries"
                )
            magnitude = float(residual[pivot])
            with np.errstate(all="ignore"):
                new_column = remainder / math.sqrt(magnitude)
                new_column[selected] = 0
                new_column[pivot] = math.sqrt(magnitude)
                root[:, rank] = new_column
                residual -= new_column * new_column
            _finite(new_column, "Cholesky column")
            residual[pivot] = 0
            selected.append(pivot)
            pivots.append(magnitude)
            ties.append(int(candidates.size))
        returned = root[:, : len(selected)]
        with np.errstate(all="ignore"):
            residual = diag - np.einsum("ij,ij->i", returned, returned)
        _finite(residual, "returned residual diagonal")
        if np.any(residual < -slack):
            raise NumericalBreakdownError(
                "returned factor has negative residual diagonal"
            )
        diagonal_error = float(np.max(np.abs(residual), initial=0))
        clipped += int(np.count_nonzero(residual < 0))
        residual = np.maximum(residual, 0)
        actual_maximum = float(np.max(residual, initial=0))
        history[-1] = actual_maximum
        if diagonal_error > threshold and converged:
            converged, stop = False, "roundoff_limit"
        elif diagonal_error <= threshold:
            converged, stop = True, "tolerance"
        smallest, condition = None, None
        if selected:
            singular = svdvals(returned[selected])
            smallest = float(singular[-1])
            if smallest > 0:
                with np.errstate(all="ignore"):
                    ratio = float(singular[0] / singular[-1])
                condition = ratio if math.isfinite(ratio) else None
        object.__setattr__(self, "factor", _owned(returned))
        object.__setattr__(self, "diagonal", _owned(diag))
        object.__setattr__(self, "residual_diagonal", _owned(residual))
        object.__setattr__(self, "indices", tuple(selected))
        object.__setattr__(self, "pivots", tuple(pivots))
        object.__setattr__(
            self,
            "diagnostics",
            KernelCholeskyDiagnostics(
                converged,
                stop,
                threshold,
                slack,
                tie,
                tuple(ties),
                tuple(history),
                clipped,
                diagonal_error,
                smallest,
                condition,
                len(selected),
            ),
        )

    @classmethod
    def from_columns(
        cls,
        diagonal: npt.ArrayLike,
        column: ColumnOracle,
        *,
        max_rank: int = 256,
        atol: float = 1e-10,
        rtol: float = 0.0,
        tie_rtol: float = 32 * _EPS,
        psd_rtol: float = 128 * _EPS,
    ) -> KernelCholeskyFactor:
        """Construct without requesting or retaining an n-by-n kernel matrix."""
        return cls(
            diagonal,
            column,
            max_rank=max_rank,
            atol=atol,
            rtol=rtol,
            tie_rtol=tie_rtol,
            psd_rtol=psd_rtol,
        )

    @property
    def rank(self) -> int:
        return int(self.factor.shape[1])

    @property
    def converged(self) -> bool:
        return self.diagnostics.converged

    @property
    def stop_reason(self) -> str:
        return self.diagnostics.stop_reason

    @property
    def shape(self) -> tuple[int, int]:
        return (self.diagonal.size, self.diagonal.size)

    @property
    def storage_nbytes(self) -> int:
        """Owned numerical arrays and scalar histories; excludes Python overhead."""
        return (
            self.factor.nbytes
            + self.diagonal.nbytes
            + self.residual_diagonal.nbytes
            + 8
            * (
                len(self.indices)
                + len(self.pivots)
                + len(self.diagnostics.residual_history)
                + len(self.diagnostics.tie_counts)
            )
        )

    @property
    def metadata(self) -> NumericalMetadata:
        """Return a fresh metadata snapshot; mutating it cannot change the factor."""
        d = self.diagnostics
        return NumericalMetadata(
            algorithm="column_oracle_pivoted_cholesky",
            arithmetic="float64",
            precision_bits=53,
            requested_precision_bits=53,
            precision_status="as_requested",
            conditioning_status="recorded" if self.rank else "not_applicable",
            conditioning_method="selected_root_svd_once_at_termination",
            conditioning_history=[
                ConditioningDiagnostic(
                    stage="selected_root",
                    status="measured",
                    columns=self.rank,
                    sigma_min=d.root_sigma_min,
                    condition_number=d.root_condition_number,
                )
            ]
            if self.rank
            else [],
            tie_status="recorded" if self.rank else "not_applicable",
            tie_policy="first_unselected_positive_index_within_relative_tolerance",
            tie_rtol=d.tie_rtol,
            tie_tolerance_rule="tie_rtol*maximum_residual_diagonal",
            tie_counts=list(d.tie_counts),
            converged=d.converged,
            stop_reason=d.stop_reason,
            convergence_scope="returned_finite_max_absolute_diagonal_reconstruction_error",
            certification_scope="callback_psd_assumed_not_verified_no_rounding_certificate",
            configuration={
                "threshold": d.threshold,
                "roundoff_slack": d.roundoff_slack,
                "returned_diagonal_error": d.returned_diagonal_error,
                "column_evaluations": d.column_evaluations,
                "retained_numerical_bytes": self.storage_nbytes,
            },
        )

    def extend(
        self, cross_columns: npt.ArrayLike, *, diagonal: npt.ArrayLike | None = None
    ) -> KernelCholeskyExtension:
        """Extend using K[new, selected pivots], in the recorded pivot order.

        Optional diagonal data checks the new residual's nonnegativity. This
        neither queries nor reconstructs a full training/test covariance.
        """
        cross = _array(cross_columns, "cross_columns", ndim=2)
        if cross.shape[1] != self.rank:
            raise ValueError(
                "cross_columns must have exactly rank columns in pivot order"
            )
        if self.rank:
            with np.errstate(all="ignore"):
                features = solve_triangular(
                    self.factor[list(self.indices)],
                    cross.T,
                    lower=True,
                    check_finite=False,
                ).T
            _finite(features, "out-of-sample extension")
        else:
            features = np.zeros_like(cross)
        remaining = None
        if diagonal is not None:
            diag = _array(diagonal, "new diagonal", ndim=1)
            if diag.shape != (cross.shape[0],) or np.any(diag < 0):
                raise ValueError("new diagonal must be a nonnegative length-m vector")
            with np.errstate(all="ignore"):
                remaining = diag - np.einsum("ij,ij->i", features, features)
            _finite(remaining, "out-of-sample residual diagonal")
            allowance = 128 * _EPS * float(np.max(diag, initial=0))
            if np.any(remaining < -allowance):
                raise NumericalBreakdownError(
                    "negative out-of-sample residual diagonal"
                )
            remaining = _owned(np.maximum(remaining, 0))
        return KernelCholeskyExtension(_owned(features), remaining)

    def as_linear_operator(self, shift: npt.ArrayLike = 0.0) -> Any:
        """Return SciPy's operator for L L.T + diag(shift), without dense state."""
        diagonal = _shift(shift, self.shape[0], positive=False)

        def multiply(value: npt.ArrayLike) -> FloatArray:
            rhs = _array(value, "operator input")
            if rhs.ndim not in (1, 2) or rhs.shape[0] != self.shape[0]:
                raise ValueError("operator input must have leading dimension n")
            with np.errstate(all="ignore"):
                result = self.factor @ (self.factor.T @ rhs)
                result += diagonal * rhs if rhs.ndim == 1 else diagonal[:, None] * rhs
            _finite(result, "operator multiplication")
            return result

        return LinearOperator(
            self.shape,
            matvec=multiply,
            matmat=multiply,
            rmatvec=multiply,
            dtype=np.float64,
        )

    def shifted(
        self,
        shift: npt.ArrayLike,
        *,
        residual_tolerance: float = 1e-10,
        max_refinements: int = 3,
    ) -> ShiftedKernelCholesky:
        """Prepare solve/preconditioner state for a strictly positive diagonal."""
        return ShiftedKernelCholesky(
            self,
            shift,
            residual_tolerance=residual_tolerance,
            max_refinements=max_refinements,
        )

    def to_gpytorch(self, shift: npt.ArrayLike = 0.0) -> Any:
        """Copy into a CPU float64 linear_operator, with no selection autograd.

        The provider object owns its tensor; mutating it cannot modify this
        factor. Optional imports occur only on this explicit call.
        """
        try:
            torch = importlib.import_module("torch")
            operators = importlib.import_module("linear_operator.operators")
        except ImportError as error:
            raise ImportError(
                "CPU GP conversion requires torch and linear_operator"
            ) from error
        diagonal = _shift(shift, self.shape[0], positive=False)
        tensor = torch.tensor(self.factor.copy(), dtype=torch.float64, device="cpu")
        if self.rank:
            result = operators.LowRankRootLinearOperator(tensor)
        else:
            result = operators.DiagLinearOperator(
                torch.zeros(self.shape[0], dtype=torch.float64, device="cpu")
            )
        return result.add_diagonal(
            torch.tensor(diagonal.copy(), dtype=torch.float64, device="cpu")
        )


@dataclass(frozen=True, slots=True)
class KernelSolveResult:
    """Measured residual of the approximate shifted system, not the true kernel."""

    solution: FloatArray
    relative_residual: float
    refinements: int
    method: str

    def __post_init__(self) -> None:
        solution = _array(self.solution, "solution")
        if solution.ndim not in (1, 2):
            raise ValueError("solution must be a vector or matrix")
        object.__setattr__(self, "solution", _owned(solution))
        object.__setattr__(
            self,
            "relative_residual",
            _real(self.relative_residual, "relative_residual"),
        )
        object.__setattr__(self, "refinements", _count(self.refinements, "refinements"))
        if not isinstance(self.method, str) or not self.method:
            raise ValueError("method must be a nonempty string")


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ShiftedKernelCholesky:
    """Small solve-based shifted state, usable as a SciPy CG preconditioner.

    Retains a rank-by-rank Cholesky system, not another n-by-r basis. If
    subtraction in Woodbury loses accuracy, a temporary Householder solve
    avoids it. Every returned solution meets the requested measured residual.
    """

    factor: KernelCholeskyFactor
    shift: FloatArray
    residual_tolerance: float
    max_refinements: int
    _small_cholesky: FloatArray

    def __init__(
        self,
        factor: KernelCholeskyFactor,
        shift: npt.ArrayLike,
        *,
        residual_tolerance: float = 1e-10,
        max_refinements: int = 3,
    ) -> None:
        if not isinstance(factor, KernelCholeskyFactor):
            raise TypeError("factor must be a KernelCholeskyFactor")
        diagonal = _shift(shift, factor.shape[0], positive=True)
        tolerance = _real(residual_tolerance, "residual_tolerance", positive=True)
        refinements = _count(max_refinements, "max_refinements")
        if tolerance >= 1:
            raise ValueError("residual_tolerance must be smaller than one")
        with np.errstate(all="ignore"):
            scaled = factor.factor / np.sqrt(diagonal)[:, None]
            small = np.eye(factor.rank) + scaled.T @ scaled
        _finite(small, "shifted small system")
        try:
            chol = (
                cho_factor(small, lower=True, check_finite=False)[0]
                if factor.rank
                else small
            )
        except np.linalg.LinAlgError as error:
            raise NumericalBreakdownError(
                "shifted small Cholesky system is unresolved in float64"
            ) from error
        _finite(chol, "shifted Cholesky system")
        object.__setattr__(self, "factor", factor)
        object.__setattr__(self, "shift", diagonal)
        object.__setattr__(self, "residual_tolerance", tolerance)
        object.__setattr__(self, "max_refinements", refinements)
        object.__setattr__(self, "_small_cholesky", _owned(chol))

    @property
    def storage_nbytes(self) -> int:
        """Includes referenced factor once, shift, and retained small solve state."""
        return (
            self.factor.storage_nbytes + self.shift.nbytes + self._small_cholesky.nbytes
        )

    def _woodbury(self, rhs: FloatArray) -> FloatArray:
        with np.errstate(all="ignore"):
            scaled = rhs / self.shift[:, None]
            if not self.factor.rank:
                return scaled
            correction = cho_solve(
                (self._small_cholesky, True),
                self.factor.factor.T @ scaled,
                check_finite=False,
            )
            return cast(
                FloatArray,
                scaled - (self.factor.factor @ correction) / self.shift[:, None],
            )

    def _householder(self, rhs: FloatArray) -> FloatArray:
        # Q is stored as n*r reflectors. No complete n*n Q is materialized.
        with np.errstate(all="ignore"):
            root = self.factor.factor / np.sqrt(self.shift)[:, None]
            scaled_rhs = rhs / np.sqrt(self.shift)[:, None]
        _finite(root, "scaled root")
        _finite(scaled_rhs, "scaled right-hand side")
        (reflectors, tau), triangular = qr(root, mode="raw", check_finite=False)

        def apply(values: FloatArray, trans: bytes) -> FloatArray:
            _, workspace, info = dormqr(b"L", trans, reflectors, tau, values, -1)
            if info != 0:
                raise NumericalBreakdownError("Householder workspace query failed")
            transformed, _, info = dormqr(
                b"L", trans, reflectors, tau, values, int(workspace[0])
            )
            if info != 0:
                raise NumericalBreakdownError("Householder application failed")
            return cast(FloatArray, transformed)

        transformed = apply(scaled_rhs, b"T")
        # A small SVD avoids a second potentially ill-conditioned Gram solve.
        u, singular, _ = np.linalg.svd(triangular, full_matrices=False)
        with np.errstate(all="ignore"):
            denominators = 1 + singular * singular
            transformed[: self.factor.rank] = u @ (
                (u.T @ transformed[: self.factor.rank]) / denominators[:, None]
            )
            result = apply(transformed, b"N") / np.sqrt(self.shift)[:, None]
        _finite(result, "Householder solution")
        return result

    def solve_with_diagnostics(self, rhs: npt.ArrayLike) -> KernelSolveResult:
        """Solve one or multiple real RHS; reject inaccurate/unrepresentable output."""
        values = _array(rhs, "rhs")
        if values.ndim not in (1, 2) or values.shape[0] != self.factor.shape[0]:
            raise ValueError("rhs must be a vector or matrix with leading dimension n")
        matrix = values[:, None] if values.ndim == 1 else values
        if matrix.size == 0 or not np.any(matrix):
            return KernelSolveResult(_owned(values * 0), 0.0, 0, "zero_rhs")
        operation = self.factor.as_linear_operator(self.shift)
        scales = np.max(np.abs(matrix), axis=0)
        # Per-column scaling prevents vector-norm overflow and avoids hiding
        # a small RHS behind a much larger neighboring RHS.
        active = scales > 0
        normalized = np.zeros_like(matrix)
        normalized[:, active] = matrix[:, active] / scales[active]
        norms = np.linalg.norm(normalized, axis=0)

        def residual(solution: FloatArray) -> tuple[FloatArray, float]:
            difference = normalized - operation.matmat(solution)
            _finite(difference, "shifted solve residual")
            errors = np.linalg.norm(difference[:, active], axis=0) / norms[active]
            return difference, float(np.max(errors, initial=0))

        solution = self._woodbury(normalized)
        method, refinements = "woodbury_cholesky", 0
        for attempt in range(self.max_refinements + 1):
            if np.all(np.isfinite(solution)):
                difference, error = residual(solution)
                if error <= self.residual_tolerance:
                    break
            else:
                error = math.inf
            if attempt == self.max_refinements:
                break
            if not math.isfinite(error):
                break
            solution += self._woodbury(difference)
            refinements += 1
        if error > self.residual_tolerance and self.factor.rank:
            solution = self._householder(normalized)
            method = "implicit_householder_svd"
            _, error = residual(solution)
        if not math.isfinite(error) or error > self.residual_tolerance:
            raise NumericalBreakdownError(
                f"shifted solve residual {error:.3g} exceeds "
                f"{self.residual_tolerance:.3g}"
            )
        with np.errstate(all="ignore"):
            result = solution * scales
        _finite(result, "shifted solution")
        # Rescaling can itself underflow. Gate the returned binary64 values,
        # not only their pre-rescaling intermediate representation.
        normalized_output = np.zeros_like(result)
        normalized_output[:, active] = result[:, active] / scales[active]
        _, error = residual(normalized_output)
        if error > self.residual_tolerance:
            raise NumericalBreakdownError(
                "returned shifted solve loses accuracy on rescaling"
            )
        output = result[:, 0] if values.ndim == 1 else result
        return KernelSolveResult(_owned(output), error, refinements, method)

    def solve(self, rhs: npt.ArrayLike) -> FloatArray:
        """Return the owning solution after the measured residual gate passes."""
        return self.solve_with_diagnostics(rhs).solution

    def as_preconditioner(self) -> Any:
        """Return M=(L L.T + D)^-1 for scipy.sparse.linalg.cg's M argument."""
        return LinearOperator(
            self.factor.shape,
            matvec=self.solve,
            matmat=self.solve,
            rmatvec=self.solve,
            dtype=np.float64,
        )
