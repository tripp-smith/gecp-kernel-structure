"""Finite-query interval certificates for actual fermionic fit arithmetic.

The ideal target is full-column-rank least squares on the supplied exact binary
observations, using exact physical tau/beta. It is not unknown physical truth.
Interval arithmetic is a trusted numerical component, not a verified program.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, fields, is_dataclass
from fractions import Fraction
from typing import Any, cast

import numpy as np
import numpy.typing as npt
from mpmath.ctx_iv import MPIntervalContext  # type: ignore[import-untyped]

from .basis import FermionicBasis, FermionicFit, _indices_array, _tau_array
from .fitting import NumericArray, _immutable_copy, _numeric_array
from .kernels import FermionicKernel
from .types import FloatArray, NumericalBreakdownError


class _Uncertified(ArithmeticError):
    """An enclosure could not be completed within the declared arithmetic."""


def _bounded_integer(value: int, name: str, lower: int, upper: int) -> int:
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or not lower <= value <= upper
    ):
        raise ValueError(f"{name} must be an integer in [{lower}, {upper}]")
    return int(value)


def _certificate_hash(array: npt.ArrayLike) -> str:
    values = np.asarray(array)
    if values.dtype.kind not in "iufc" or not np.all(np.isfinite(values)):
        raise ValueError("certificate hashes require finite numeric arrays")
    values = np.asarray(values, dtype=values.dtype.newbyteorder("<"), order="C")
    header = json.dumps(
        {"dtype": values.dtype.str, "shape": values.shape}, sort_keys=True
    ).encode()
    return hashlib.sha256(header + values.tobytes()).hexdigest()


def _exact(ctx: Any, value: Any) -> Any:
    """Import a numeric scalar as its exact binary rational, not repr(float)."""
    if hasattr(value, "_mpi_") or hasattr(value, "_mpci_"):
        return value
    if isinstance(value, (complex, np.complexfloating)):
        number = complex(value)
        return ctx.mpc(_exact(ctx, number.real), _exact(ctx, number.imag))
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return ctx.mpf(int(value))
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("interval input must be finite")
    numerator, denominator = number.as_integer_ratio()
    return ctx.mpf(numerator) / ctx.mpf(denominator)


def _directed(value: Fraction, upward: bool) -> float:
    try:
        result = float(value)
    except OverflowError as exc:
        raise _Uncertified("bound_export_overflow") from exc
    if not math.isfinite(result):
        raise _Uncertified("bound_export_overflow")
    if (Fraction(result) < value) if upward else (Fraction(result) > value):
        result = math.nextafter(result, math.inf if upward else -math.inf)
    if not math.isfinite(result):
        raise _Uncertified("bound_export_overflow")
    return result


def _endpoint(value: Any, upward: bool) -> float:
    sign, mantissa, exponent, bitcount = value._mpi_[1 if upward else 0]
    if bitcount < 0:
        raise _Uncertified("bound_export_overflow")
    if not mantissa:
        return 0.0
    if exponent + bitcount > 1024:
        raise _Uncertified("bound_export_overflow")
    if exponent + bitcount < -1074:
        tiny = math.nextafter(0.0, math.inf)
        return (-tiny if not upward else 0.0) if sign else (tiny if upward else 0.0)
    numerator = -mantissa if sign else mantissa
    exact = (
        Fraction(numerator << exponent)
        if exponent >= 0
        else Fraction(numerator, 1 << -exponent)
    )
    return _directed(exact, upward)


def _upper_abs(value: Any) -> float:
    return max(0.0, _endpoint(abs(value), True))


def _lower_abs(value: Any) -> float:
    return max(0.0, _endpoint(abs(value), False))


def _real_imag(value: Any, ctx: Any) -> tuple[Any, Any]:
    if hasattr(value, "_mpci_"):
        return value.real, value.imag
    return value, ctx.mpf(0)


class _IntervalLU:
    """Small solve-based interval PLU with an explicit zero-exclusion gate.

    Pivot order uses interval lower magnitudes. Every division is made only
    after the selected real pivot interval excludes zero. The factorization
    encloses the ordinary elimination of every matrix represented by inputs.
    """

    def __init__(self, matrix: list[list[Any]], ctx: Any, stage: str) -> None:
        self.ctx = ctx
        self.stage = stage
        self.rank = len(matrix)
        if not self.rank or any(len(row) != self.rank for row in matrix):
            raise ValueError("interval LU requires a nonempty square matrix")
        self.factors = [row.copy() for row in matrix]
        self.permutation = list(range(self.rank))
        self.pivot_lower_bounds: list[float] = []
        for column in range(self.rank):
            chosen = max(
                range(column, self.rank),
                key=lambda row: abs(self.factors[row][column]).a,
            )
            pivot = self.factors[chosen][column]
            if not (pivot.a > 0 or pivot.b < 0):
                raise _Uncertified(f"{stage}_pivot_contains_zero")
            self.factors[column], self.factors[chosen] = (
                self.factors[chosen],
                self.factors[column],
            )
            self.permutation[column], self.permutation[chosen] = (
                self.permutation[chosen],
                self.permutation[column],
            )
            pivot = self.factors[column][column]
            self.pivot_lower_bounds.append(_lower_abs(pivot))
            for row in range(column + 1, self.rank):
                multiplier = self.factors[row][column] / pivot
                self.factors[row][column] = multiplier
                for index in range(column + 1, self.rank):
                    self.factors[row][index] -= multiplier * self.factors[column][index]

    def solve(self, rhs: list[Any]) -> list[Any]:
        if len(rhs) != self.rank:
            raise ValueError("interval RHS must match core rank")
        forward: list[Any] = []
        for row in range(self.rank):
            forward.append(
                rhs[self.permutation[row]]
                - sum(
                    (self.factors[row][index] * forward[index] for index in range(row)),
                    self.ctx.mpf(0),
                )
            )
        result = [self.ctx.mpf(0) for _ in range(self.rank)]
        for row in range(self.rank - 1, -1, -1):
            result[row] = (
                forward[row]
                - sum(
                    (
                        self.factors[row][index] * result[index]
                        for index in range(row + 1, self.rank)
                    ),
                    self.ctx.mpf(0),
                )
            ) / self.factors[row][row]
        return result


class _IntervalBasis:
    """Exact-physical core and cached interval transpose solves.

    Helpers also accept same-context interval query times/frequencies for
    separate bounded consumers. This class itself makes no domain-wide claim.
    """

    def __init__(self, basis: FermionicBasis, precision_bits: int = 256) -> None:
        self.ctx = MPIntervalContext()
        self.ctx.prec = _bounded_integer(precision_bits, "precision_bits", 53, 4096)
        self.c = self.ctx
        self.beta = _exact(self.ctx, basis.beta)
        self.rank = basis.rank
        self.omega = [_exact(self.ctx, value) for value in basis.omega_nodes]
        self.core = [
            [self.kernel(float(time), omega) for omega in self.omega]
            for time in basis.tau_nodes
        ]
        self.transpose_lu = _IntervalLU(
            [list(column) for column in zip(*self.core, strict=True)], self.ctx, "core"
        )

    def kernel(self, tauphysical: Any, omega: Any) -> Any:
        time = _exact(self.ctx, tauphysical) / self.beta
        frequency = _exact(self.ctx, omega)
        if frequency.a >= 0:
            return self.ctx.exp(-time * frequency) / (1 + self.ctx.exp(-frequency))
        if frequency.b <= 0:
            return self.ctx.exp((1 - time) * frequency) / (1 + self.ctx.exp(frequency))
        # Stable sign-split hull for a frequency cell crossing zero.
        negative = self.ctx.mpf([frequency.a, 0])
        positive = self.ctx.mpf([0, frequency.b])
        left = self.ctx.exp((1 - time) * negative) / (1 + self.ctx.exp(negative))
        right = self.ctx.exp(-time * positive) / (1 + self.ctx.exp(-positive))
        return self.ctx.mpf([min(left.a, right.a), max(left.b, right.b)])

    def tau_rows(self, tau: Any) -> list[list[Any]]:
        return [
            self.transpose_lu.solve([self.kernel(time, omega) for omega in self.omega])
            for time in tau
        ]

    def frequency_columns(self, indices: Any) -> list[list[Any]]:
        return [
            [
                -self.beta / (self.ctx.j * self.ctx.pi * (2 * int(index) + 1) - omega)
                for omega in self.omega
            ]
            for index in indices
        ]

    def frequency_rows(self, indices: Any) -> list[list[Any]]:
        return [self.transpose_lu.solve(row) for row in self.frequency_columns(indices)]


def _readonly(values: npt.ArrayLike) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    return cast(
        FloatArray,
        np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape),
    )


def _maximum(array: FloatArray) -> float:
    return float(np.max(array, initial=0.0))


def _json_data(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if is_dataclass(value) and not isinstance(value, type):
        return {
            item.name: _json_data(getattr(value, item.name)) for item in fields(value)
        }
    if isinstance(value, (list, tuple)):
        return [_json_data(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


@dataclass(frozen=True, slots=True)
class IntervalArray:
    """Immutable outward binary64 exports of an interval array."""

    real_lower: FloatArray
    real_upper: FloatArray
    imag_lower: FloatArray
    imag_upper: FloatArray
    bound_precision_bits: int = 53

    @property
    def shape(self) -> tuple[int, ...]:
        return self.real_lower.shape


@dataclass(frozen=True, slots=True)
class CorePrecisionBounds:
    """Physical core versus rounded normalization and a binary64 kernel proxy.

    ``binary64_kernel_error_bound`` already includes normalization error; these
    two informational bounds must not be treated as disjoint additive stages.
    Reference factor entries are not reconstructed by the binary64 proxy.
    """

    ideal: IntervalArray
    normalization_error_bound: FloatArray
    binary64_kernel_error_bound: FloatArray
    transpose_lu_pivot_lower_bounds: tuple[float, ...]
    comparison_scope: str = (
        "reconstructed_binary64_kernel_not_internal_reference_factors"
    )


@dataclass(frozen=True, slots=True)
class DesignPrecisionBounds:
    """Ideal rows versus actual exports, and their true-core backward residual.

    The direct design discrepancy includes normalization, kernel, factorization,
    solve and export effects together. The residual is evaluated independently
    against the enclosed ideal core and RHS, without assuming LU stability.
    """

    ideal: IntervalArray
    actual_error_bound: FloatArray
    core_solve_residual_bound: FloatArray

    @property
    def max_error(self) -> float:
        return _maximum(self.actual_error_bound)

    @property
    def max_solve_residual(self) -> float:
        return _maximum(self.core_solve_residual_bound)


@dataclass(frozen=True, slots=True)
class CoefficientPrecisionBounds:
    """Full-rank ideal LS coefficients versus any supplied frozen fit.

    A successful interval Gram solve proves full column rank of the ideal time
    design. It does not assert numerical rank equivalence with the actual SVD.
    """

    ideal: IntervalArray
    actual_error_bound: FloatArray
    normal_solve_residual_bound: FloatArray
    gram_lu_pivot_lower_bounds: tuple[float, ...]

    @property
    def max_error(self) -> float:
        return _maximum(self.actual_error_bound)


@dataclass(frozen=True, slots=True)
class QueryPrecisionBounds:
    """Direct prediction enclosure and a separately justified additive budget."""

    design: DesignPrecisionBounds
    ideal_prediction: IntervalArray
    direct_error_bound: FloatArray
    coefficient_contribution: FloatArray
    design_contribution: FloatArray
    dot_rounding_bound: FloatArray
    total_error_bound: FloatArray

    @property
    def max_error(self) -> float:
        return _maximum(self.direct_error_bound)


@dataclass(frozen=True, slots=True)
class FitPrecisionCertificate:
    """Finite-query arithmetic evidence, including valid large-error reports.

    ``certified`` means enclosures completed. It is not an accuracy tolerance
    pass, a statement about physical truth, or a continuous error certificate.
    Missing bounds in an inconclusive report are None, never invented zeros.
    """

    certified: bool
    termination_reason: str
    content_sha256: str
    input_hashes: tuple[tuple[str, str], ...]
    rank: int
    precision_bits: int
    estimated_work: int
    observation_axis: int
    core: CorePrecisionBounds | None = None
    training_design: DesignPrecisionBounds | None = None
    coefficients: CoefficientPrecisionBounds | None = None
    time: QueryPrecisionBounds | None = None
    matsubara: QueryPrecisionBounds | None = None
    objective_scope: str = "finite_query_exact_physical_full_column_rank_least_squares"
    arithmetic: str = "owned_mpmath_interval_lu_gram_and_exact_binary_rational_dot"
    coordinate_precision_bits: int = 53
    observation_precision_bits: int = 53
    coefficient_precision_bits: int = 53
    bound_precision_bits: int = 53
    continuous_certificate: bool = False
    physical_truth_certificate: bool = False
    conditioning_status: str = "interval_core_and_gram_pivot_gates"
    conditioning_scope: str = (
        "regularity_when_certified_not_a_norm_condition_number_bound"
    )

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], _json_data(self))


def _reshape(values: FloatArray, shape: tuple[int, ...], axis: int) -> FloatArray:
    return _readonly(np.moveaxis(values.reshape(shape), 0, axis))


def _enclosure(
    matrix: list[list[Any]], ctx: Any, shape: tuple[int, ...], axis: int = 0
) -> IntervalArray:
    components: list[list[float]] = [[], [], [], []]
    for row in matrix:
        for value in row:
            real, imag = _real_imag(value, ctx)
            for target, interval, upward in (
                (components[0], real, False),
                (components[1], real, True),
                (components[2], imag, False),
                (components[3], imag, True),
            ):
                target.append(_endpoint(interval, upward))
    arrays = [_reshape(np.asarray(values), shape, axis) for values in components]
    return IntervalArray(arrays[0], arrays[1], arrays[2], arrays[3])


def _error_matrix(ideal: list[list[Any]], actual: NumericArray, ctx: Any) -> FloatArray:
    return _readonly(
        [
            [
                _upper_abs(value - _exact(ctx, actual[i, j]))
                for j, value in enumerate(row)
            ]
            for i, row in enumerate(ideal)
        ]
    )


def _design_bounds(
    exact: _IntervalBasis,
    ideal: list[list[Any]],
    actual: NumericArray,
    rhs: list[list[Any]],
) -> DesignPrecisionBounds:
    residual = [
        [
            sum(
                (
                    exact.core[k][j] * _exact(exact.ctx, actual[i, k])
                    for k in range(exact.rank)
                ),
                exact.ctx.mpf(0),
            )
            - rhs[i][j]
            for j in range(exact.rank)
        ]
        for i in range(len(ideal))
    ]
    return DesignPrecisionBounds(
        _enclosure(ideal, exact.ctx, actual.shape),
        _error_matrix(ideal, actual, exact.ctx),
        _readonly([[_upper_abs(value) for value in row] for row in residual]),
    )


def _fraction_complex(value: Any) -> tuple[Fraction, Fraction]:
    return Fraction(float(np.real(value))), Fraction(float(np.imag(value)))


def _exact_dot_error(
    design: NumericArray, coefficients: NumericArray, actual: NumericArray, ctx: Any
) -> FloatArray:
    result = np.empty(actual.shape, dtype=np.float64)
    for row in range(design.shape[0]):
        for column in range(coefficients.shape[1]):
            real, imag = Fraction(0), Fraction(0)
            for index in range(design.shape[1]):
                ar, ai = _fraction_complex(design[row, index])
                cr, ci = _fraction_complex(coefficients[index, column])
                real += ar * cr - ai * ci
                imag += ar * ci + ai * cr
            vr, vi = _fraction_complex(actual[row, column])
            dr, di = real - vr, imag - vi
            if not dr and not di:
                result[row, column] = 0.0
            else:
                squared = dr * dr + di * di
                squared_interval = ctx.mpf(squared.numerator) / ctx.mpf(
                    squared.denominator
                )
                result[row, column] = _upper_abs(ctx.sqrt(squared_interval))
    return _readonly(result)


def _query_bounds(
    exact: _IntervalBasis,
    ideal_design: list[list[Any]],
    actual_design: NumericArray,
    rhs: list[list[Any]],
    ideal_coefficients: list[list[Any]],
    actual_coefficients: NumericArray,
    actual_output: NumericArray,
    batch_shape: tuple[int, ...],
    axis: int,
) -> QueryPrecisionBounds:
    ctx, rank = exact.ctx, exact.rank
    rows, batch = actual_output.shape
    ideal_output = [
        [
            sum(
                (ideal_design[i][k] * ideal_coefficients[k][j] for k in range(rank)),
                ctx.mpf(0),
            )
            for j in range(batch)
        ]
        for i in range(rows)
    ]
    coefficient = np.empty((rows, batch))
    design = np.empty((rows, batch))
    total = np.empty((rows, batch))
    dot = _exact_dot_error(actual_design, actual_coefficients, actual_output, ctx)
    for i in range(rows):
        for j in range(batch):
            coefficient_term = sum(
                (
                    abs(ideal_design[i][k])
                    * abs(
                        ideal_coefficients[k][j]
                        - _exact(ctx, actual_coefficients[k, j])
                    )
                    for k in range(rank)
                ),
                ctx.mpf(0),
            )
            design_term = sum(
                (
                    abs(ideal_design[i][k] - _exact(ctx, actual_design[i, k]))
                    * abs(_exact(ctx, actual_coefficients[k, j]))
                    for k in range(rank)
                ),
                ctx.mpf(0),
            )
            coefficient[i, j] = _upper_abs(coefficient_term)
            design[i, j] = _upper_abs(design_term)
            total[i, j] = _directed(
                Fraction(coefficient[i, j])
                + Fraction(design[i, j])
                + Fraction(dot[i, j]),
                True,
            )
    shape = (rows, *batch_shape)
    return QueryPrecisionBounds(
        _design_bounds(exact, ideal_design, actual_design, rhs),
        _enclosure(ideal_output, ctx, shape, axis),
        _reshape(_error_matrix(ideal_output, actual_output, ctx), shape, axis),
        _reshape(coefficient, shape, axis),
        _reshape(design, shape, axis),
        _reshape(dot, shape, axis),
        _reshape(total, shape, axis),
    )


def certify_fit_precision(
    fit: FermionicFit,
    values: npt.ArrayLike,
    tau: npt.ArrayLike | None = None,
    *,
    query_tau: npt.ArrayLike | None = None,
    matsubara_indices: npt.ArrayLike | None = None,
    precision_bits: int = 256,
    max_rank: int = 64,
    max_work: int = 2**22,
) -> FitPrecisionCertificate:
    """Enclose actual fit/evaluation discrepancy from an ideal physical LS fit.

    All supplied numeric inputs denote their exact binary64 values after the
    ordinary public API conversion. Data use ``fit.axis`` as their observation
    axis. The supplied fit can be truncated or manually constructed; a large
    valid discrepancy does not become an accuracy pass. At least one nonempty
    query domain is required. Budgets are checked before interval allocation.
    """
    if not isinstance(fit, FermionicFit) or not isinstance(fit.basis, FermionicBasis):
        raise TypeError("fit must be a FermionicFit with a FermionicBasis")
    bits = _bounded_integer(precision_bits, "precision_bits", 53, 4096)
    rank_budget = _bounded_integer(max_rank, "max_rank", 1, 256)
    work_budget = _bounded_integer(max_work, "max_work", 1, 2**30)
    basis = fit.basis
    training = _readonly(
        _tau_array(basis.tau_nodes if tau is None else tau, basis.beta)
    )
    observations = _immutable_copy(_numeric_array(values, "values"))
    coefficients = _immutable_copy(_numeric_array(fit.coefficients, "coefficients"))
    axis = fit.axis
    if observations.ndim != coefficients.ndim or observations.ndim <= axis:
        raise ValueError(
            "values must have the fit's batch dimensions and observation axis"
        )
    expected_shape = list(coefficients.shape)
    expected_shape[axis] = training.size
    if observations.shape != tuple(expected_shape) or not observations.size:
        raise ValueError(
            "values shape must match training tau and fit batch dimensions"
        )
    times = _readonly(_tau_array([] if query_tau is None else query_tau, basis.beta))
    modes = np.array(
        _indices_array([] if matsubara_indices is None else matsubara_indices),
        copy=True,
    )
    if not times.size and not modes.size:
        raise ValueError("at least one finite query domain must be nonempty")
    moved = np.moveaxis(observations, axis, 0)
    batch_shape = moved.shape[1:]
    batch = math.prod(batch_shape)
    data_matrix = moved.reshape(training.size, batch)
    coefficient_matrix = np.moveaxis(coefficients, axis, 0).reshape(basis.rank, batch)
    rank, m, q = basis.rank, training.size, times.size + modes.size
    work = int(
        (
            2 * rank**3
            + (2 * m + q) * rank**2
            + (m * rank + rank**2 + 2 * q * rank) * batch
        )
        * ((bits + 63) // 64)
    )
    hashes = [
        ("basis_parameters", _certificate_hash([basis.beta, basis.wmax, basis.cutoff])),
        ("basis_tau_nodes", _certificate_hash(basis.tau_nodes)),
        ("basis_omega_nodes", _certificate_hash(basis.omega_nodes)),
        ("training_tau", _certificate_hash(training)),
        ("observations", _certificate_hash(observations)),
        ("coefficients", _certificate_hash(coefficients)),
        ("query_tau", _certificate_hash(times)),
        ("matsubara_indices", _certificate_hash(modes)),
    ]

    def report(reason: str, **bounds: Any) -> FitPrecisionCertificate:
        metadata = {
            "input_hashes": hashes,
            "rank": rank,
            "precision_bits": bits,
            "observation_axis": axis,
            "estimated_work": work,
        }
        digest = hashlib.sha256(
            json.dumps(metadata, sort_keys=True).encode()
        ).hexdigest()
        return FitPrecisionCertificate(
            reason == "enclosed",
            reason,
            digest,
            tuple(hashes),
            rank,
            bits,
            work,
            axis,
            **bounds,
        )

    if rank > rank_budget:
        return report("rank_budget")
    if work > work_budget:
        return report("work_budget")
    if m < rank:
        return report("underdetermined_training_design")
    try:
        with np.errstate(all="ignore"):
            actual_training = _immutable_copy(basis.design_tau(training))
            actual_time = (
                _immutable_copy(basis.design_tau(times)) if times.size else None
            )
            actual_frequency = (
                _immutable_copy(basis.design_matsubara(modes)) if modes.size else None
            )
            output_time = (
                _immutable_copy(fit.evaluate_tau(times)) if times.size else None
            )
            output_frequency = (
                _immutable_copy(fit.evaluate_matsubara(modes)) if modes.size else None
            )
        for name, actual in (
            ("actual_training_design", actual_training),
            ("actual_time_design", actual_time),
            ("actual_frequency_design", actual_frequency),
            ("actual_time_output", output_time),
            ("actual_frequency_output", output_frequency),
        ):
            if actual is not None:
                hashes.append((name, _certificate_hash(actual)))
        exact = _IntervalBasis(basis, bits)
        ctx = exact.ctx
        ideal_training = exact.tau_rows(training)
        training_rhs = [
            [exact.kernel(float(time), omega) for omega in exact.omega]
            for time in training
        ]
        training_bounds = _design_bounds(
            exact, ideal_training, actual_training, training_rhs
        )
        gram = [
            [
                sum(
                    (ideal_training[k][i] * ideal_training[k][j] for k in range(m)),
                    ctx.mpf(0),
                )
                for j in range(rank)
            ]
            for i in range(rank)
        ]
        rhs = [
            [
                sum(
                    (
                        ideal_training[k][i] * _exact(ctx, data_matrix[k, j])
                        for k in range(m)
                    ),
                    ctx.mpf(0),
                )
                for j in range(batch)
            ]
            for i in range(rank)
        ]
        normal = _IntervalLU(gram, ctx, "normal")
        columns = [normal.solve([rhs[i][j] for i in range(rank)]) for j in range(batch)]
        ideal_coefficients = [list(row) for row in zip(*columns, strict=True)]
        coefficient_residual = [
            [
                sum(
                    (
                        gram[i][k] * _exact(ctx, coefficient_matrix[k, j])
                        for k in range(rank)
                    ),
                    ctx.mpf(0),
                )
                - rhs[i][j]
                for j in range(batch)
            ]
            for i in range(rank)
        ]
        coefficient_shape = (rank, *batch_shape)
        coefficient_bounds = CoefficientPrecisionBounds(
            _enclosure(ideal_coefficients, ctx, coefficient_shape, axis),
            _reshape(
                _error_matrix(ideal_coefficients, coefficient_matrix, ctx),
                coefficient_shape,
                axis,
            ),
            _reshape(
                _readonly(
                    [
                        [_upper_abs(value) for value in row]
                        for row in coefficient_residual
                    ]
                ),
                coefficient_shape,
                axis,
            ),
            tuple(normal.pivot_lower_bounds),
        )
        with np.errstate(all="ignore"):
            normalized = np.asarray(basis.tau_nodes / basis.beta)
        normalized_core = [
            [
                ctx.exp(-_exact(ctx, time) * omega) / (1 + ctx.exp(-omega))
                if omega.a >= 0
                else ctx.exp((1 - _exact(ctx, time)) * omega) / (1 + ctx.exp(omega))
                for omega in exact.omega
            ]
            for time in normalized
        ]
        with np.errstate(all="ignore"):
            rounded_core = FermionicKernel(basis.cutoff)(
                normalized[:, None], basis.omega_nodes[None, :]
            )
        hashes.append(("reconstructed_binary64_core", _certificate_hash(rounded_core)))
        core_bounds = CorePrecisionBounds(
            _enclosure(exact.core, ctx, (rank, rank)),
            _readonly(
                [
                    [
                        _upper_abs(exact.core[i][j] - normalized_core[i][j])
                        for j in range(rank)
                    ]
                    for i in range(rank)
                ]
            ),
            _error_matrix(exact.core, np.asarray(rounded_core), ctx),
            tuple(exact.transpose_lu.pivot_lower_bounds),
        )
        time_bounds, frequency_bounds = None, None
        if actual_time is not None and output_time is not None:
            ideal_time = exact.tau_rows(times)
            time_rhs = [
                [exact.kernel(float(time), omega) for omega in exact.omega]
                for time in times
            ]
            time_bounds = _query_bounds(
                exact,
                ideal_time,
                actual_time,
                time_rhs,
                ideal_coefficients,
                coefficient_matrix,
                np.moveaxis(output_time, axis, 0).reshape(times.size, batch),
                batch_shape,
                axis,
            )
        if actual_frequency is not None and output_frequency is not None:
            frequency_bounds = _query_bounds(
                exact,
                exact.frequency_rows(modes),
                actual_frequency,
                exact.frequency_columns(modes),
                ideal_coefficients,
                coefficient_matrix,
                np.moveaxis(output_frequency, axis, 0).reshape(modes.size, batch),
                batch_shape,
                axis,
            )
        return report(
            "enclosed",
            core=core_bounds,
            training_design=training_bounds,
            coefficients=coefficient_bounds,
            time=time_bounds,
            matsubara=frequency_bounds,
        )
    except _Uncertified as exc:
        return report(str(exc))
    except (ArithmeticError, ValueError, NumericalBreakdownError) as exc:
        return report(f"arithmetic_inconclusive:{type(exc).__name__}")
