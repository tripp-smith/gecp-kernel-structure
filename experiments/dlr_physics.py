"""Shared physical truth and a small projected Bethe self-consistency workload.

All frequency indices are ordinary fermionic integers, with
``z = 1j * (2*n+1)*pi/beta``. Truth evaluations use finite spectral measures or
the independent semicircle solution, never a fitted basis. Numerical convergence,
projection quality and independent physical accuracy are separate observations;
none is a continuous-domain or rounding certificate.
"""

from __future__ import annotations

import math
import operator
from dataclasses import asdict, dataclass
from numbers import Real
from typing import Any, Protocol, cast

import numpy as np
import numpy.typing as npt

from experiments.dlr_workloads import (
    PoleWorkload,
    _physical_parameters,
    _positive_scalar,
    _real_vector,
    common_tau_meshes,
    validate_tau_meshes,
)
from kernelgecp.fitting import LeastSquaresDiagnostics, LeastSquaresPlan
from kernelgecp.types import NumericalBreakdownError

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]
IndexArray = npt.NDArray[np.int64]


class PhysicalBasis(Protocol):
    """The common public surface needed by the physical solver."""

    @property
    def rank(self) -> int: ...

    @property
    def beta(self) -> float: ...

    def design_tau(self, tau: npt.ArrayLike) -> npt.NDArray[Any]: ...

    def design_matsubara(self, indices: npt.ArrayLike) -> npt.NDArray[Any]: ...


def _immutable(array: npt.NDArray[Any]) -> npt.NDArray[Any]:
    # A bytes-backed array also prevents callers from re-enabling write access.
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


def _integer(name: str, value: Any, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be an integer >= {minimum}") from error
    if result < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return result


def _contains_boolean(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return True
    if isinstance(value, np.ndarray):
        return value.dtype.kind == "b"
    return isinstance(value, (list, tuple)) and any(
        _contains_boolean(item) for item in value
    )


def _indices(indices: Any) -> IndexArray:
    if _contains_boolean(indices):
        raise ValueError("Matsubara indices must not contain booleans")
    raw = np.asarray(indices)
    if raw.ndim != 1 or not raw.size or raw.dtype.kind not in "iu":
        raise ValueError("Matsubara indices must be a nonempty integer vector")
    limit = 2**51 - 1
    if np.any(raw > limit) or (raw.dtype.kind == "i" and np.any(raw < -limit)):
        raise ValueError("Matsubara index magnitude must be smaller than 2**51")
    if np.unique(raw).size != raw.size:
        raise ValueError("Matsubara indices must not contain duplicates")
    return np.asarray(raw, dtype=np.int64)


def matsubara_frequencies(indices: Any, *, beta: float) -> FloatArray:
    """Return finite, nonzero physical frequencies for ordinary integer indices."""

    beta = _positive_scalar("beta", beta)
    nodes = _indices(indices)
    with np.errstate(all="ignore"):
        frequencies = (2.0 * nodes + 1.0) * (math.pi / beta)
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies == 0):
        raise ValueError("Matsubara frequencies must be finite and nonzero in float64")
    return np.asarray(frequencies, dtype=np.float64)


def physical_tau_meshes(
    *, beta: float, wmax: float, order: int = 16
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Training nodes, midpoint calibration and untouched interval fifths.

    All three actual binary64 meshes are checked for ordering, domain validity,
    duplicates and pairwise overlap after scaling. Resolution exhaustion fails
    explicitly instead of silently merging nominally different coordinates.
    """

    training, calibration = common_tau_meshes(beta=beta, wmax=wmax, order=order)
    left = training[:-1]
    widths = training[1:] - left
    final = np.column_stack((left + widths / 5.0, left + (4.0 / 5.0) * widths))
    final = np.asarray(final.ravel(), dtype=np.float64)
    validate_tau_meshes(training, final, beta=beta)
    validate_tau_meshes(calibration, final, beta=beta)
    return cast(
        tuple[FloatArray, FloatArray, FloatArray],
        tuple(_immutable(array) for array in (training, calibration, final)),
    )


def frequency_index_meshes() -> tuple[IndexArray, IndexArray]:
    """Fixed disjoint calibration and final-validation ordinary indices.

    Forward final frequencies are independent of the time-fit calibration.
    Only the middle and higher subset is withheld from the solver, whose
    observations must retain every low mode through ordinary n=17.
    Conjugate pairs use ``n -> -n-1``. No truth evaluation selects these lists.
    """

    final = np.array(
        [-801, -401, -171, -71, -38, -22, -9, -7, 6, 8, 21, 37, 70, 170, 400, 800],
        dtype=np.int64,
    )
    calibration = np.setdiff1d(
        np.unique(
            np.concatenate(
                (np.arange(-32, 32, dtype=np.int64), [-256, -128, -64, 63, 127, 255])
            )
        ),
        final,
    )
    return cast(IndexArray, _immutable(calibration)), cast(
        IndexArray, _immutable(final)
    )


def solver_validation_index_mesh() -> IndexArray:
    """The disjoint solver holdout subset, keeping essential low modes observed.

    Withholding very low Matsubara modes made the physical fit unobservable in
    a failed pilot. The fixed solver design retains all ``n`` in [-18, 17];
    final solver validation uses only ``abs(2*n+1) >= 37`` from the fresh forward
    holdout. Truth on the complete observed grid is checked separately after
    iteration and is explicitly independent truth on observations, not a holdout.
    """

    _, final = frequency_index_meshes()
    return cast(IndexArray, _immutable(final[np.abs(2 * final + 1) >= 37]))


def solver_index_mesh() -> IndexArray:
    """Ordinary indices [-128, 127], excluding middle and higher solver holdouts."""

    final = solver_validation_index_mesh()
    nodes = np.setdiff1d(np.arange(-128, 128, dtype=np.int64), final)
    return cast(IndexArray, _immutable(nodes))


def matrix_amplitudes() -> ComplexArray:
    """A fixed Hermitian positive matrix with maximum absolute entry one."""

    return cast(
        ComplexArray,
        _immutable(np.array([[1.0, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])),
    )


def matsubara_green_function(
    workload: PoleWorkload, indices: Any, *, beta: float, wmax: float
) -> ComplexArray:
    """Evaluate the independent finite-pole resolvent ``sum w/(z-p)``."""

    beta, wmax, _ = _physical_parameters(beta, wmax)
    frequencies = matsubara_frequencies(indices, beta=beta)
    poles = np.asarray(workload.poles)
    if np.any(np.abs(poles) > wmax):
        raise ValueError("poles must lie in the physical interval [-wmax, wmax]")
    with np.errstate(all="ignore"):
        values = np.sum(
            np.asarray(workload.weights)[None, :]
            / (1j * frequencies[:, None] - poles[None, :]),
            axis=1,
        )
    if not np.all(np.isfinite(values)):
        raise NumericalBreakdownError("finite-pole frequency truth is nonfinite")
    return np.asarray(values, dtype=np.complex128)


def semicircle_matsubara(indices: Any, *, beta: float, hopping: float) -> ComplexArray:
    """Causal semicircle resolvent, using a cancellation-free quadratic root.

    For imaginary ``z=iw`` the physical root is
    ``-2j*sign(w)/(abs(w)+hypot(abs(w),2*hopping))``. This is algebraically
    independent of the self-consistency iteration and has ``G(z) ~ 1/z``.
    """

    hopping = _positive_scalar("hopping", hopping)
    frequencies = matsubara_frequencies(indices, beta=beta)
    if not math.isfinite(2.0 * hopping):
        raise ValueError("twice hopping must be finite in float64")
    absolute = np.abs(frequencies)
    # Normalize before the square root and sum: even a true hypot exceeding
    # float64 range can have a representable, nonzero reciprocal resolvent.
    scale = np.maximum(absolute, 2.0 * hopping)
    with np.errstate(all="ignore"):
        normalized_frequency = absolute / scale
        normalized_root = np.hypot(normalized_frequency, 2.0 * hopping / scale)
        magnitude = (2.0 / (normalized_root + normalized_frequency)) / scale
        result = -1j * np.sign(frequencies) * magnitude
    if not np.all(np.isfinite(result)):
        raise NumericalBreakdownError("semicircle frequency truth is nonfinite")
    return np.asarray(result, dtype=np.complex128)


def semicircle_tau(
    tau: Any, *, beta: float, hopping: float, quadrature_order: int = 512
) -> FloatArray:
    """Semicircle time truth from Gauss-Chebyshev quadrature of the second kind.

    The normalized measure is ``sqrt(4*hopping**2-w**2)/(2*pi*hopping**2)``
    on ``[-2*hopping, 2*hopping]``. This fixed spectral quadrature neither calls
    the solver nor uses any basis coefficients. Its error is not certified.
    """

    hopping = _positive_scalar("hopping", hopping)
    beta, _, _ = _physical_parameters(beta, 2.0 * hopping)
    order = _integer("quadrature_order", quadrature_order, 1)
    times = _real_vector("tau", tau)
    if np.any(times < 0) or np.any(times > beta):
        raise ValueError("tau must lie in the physical interval [0, beta]")
    theta = np.arange(1, order + 1, dtype=np.float64) * (math.pi / (order + 1))
    poles = 2.0 * hopping * np.cos(theta)
    weights = (2.0 / (order + 1)) * np.sin(theta) ** 2
    positive = poles >= 0
    terms = np.empty((times.size, order), dtype=np.float64)
    with np.errstate(under="ignore"):
        terms[:, positive] = np.exp(-times[:, None] * poles[positive]) / (
            1.0 + np.exp(-beta * poles[positive])
        )
        terms[:, ~positive] = np.exp((beta - times[:, None]) * poles[~positive]) / (
            1.0 + np.exp(beta * poles[~positive])
        )
    values = -(terms @ weights)
    if not np.all(np.isfinite(values)):
        raise NumericalBreakdownError("semicircle time quadrature is nonfinite")
    return np.asarray(values, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class BetheIteration:
    iteration: int
    step_norm: float
    defect_norm: float
    projection_residual_norm: float


@dataclass(frozen=True, slots=True)
class BetheSolverResult:
    """Immutable numerical solver evidence; physical truth is evaluated elsewhere."""

    coefficients: ComplexArray
    frequency_values: ComplexArray
    indices: IndexArray
    history: tuple[BetheIteration, ...]
    converged: bool
    stop_reason: str
    fit_diagnostics: LeastSquaresDiagnostics
    beta: float
    hopping: float
    mixing: float
    tolerance: float
    max_iterations: int
    initial_projection_residual_norm: float
    certification_status: str = "not_certified"
    convergence_scope: str = "self_consistency_defect_on_supplied_frequency_grid"
    certification_scope: str = "no_continuous_input_or_rounding_certificate"

    @property
    def iterations(self) -> int:
        return len(self.history)

    @property
    def metadata(self) -> dict[str, Any]:
        """A fresh JSON-friendly diagnostics record, excluding array payloads."""

        return {
            "beta": self.beta,
            "hopping": self.hopping,
            "mixing": self.mixing,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "iterations": self.iterations,
            "converged": self.converged,
            "stop_reason": self.stop_reason,
            "initial_projection_residual_norm": self.initial_projection_residual_norm,
            "history": [asdict(item) for item in self.history],
            "fit_diagnostics": asdict(self.fit_diagnostics),
            "arithmetic": "numpy_complex128",
            "precision_bits": 53,
            "index_convention": "ordinary_integer_n;wn=(2*n+1)*pi/beta",
            "initial_state": "frequency_projection_of_1_over_z",
            "convergence_scope": self.convergence_scope,
            "certification_status": self.certification_status,
            "certification_scope": self.certification_scope,
        }


class BetheSolverPlan:
    """Cache one frequency design and SVD for repeated scalar Bethe solves.

    ``design`` and ``fit_plan`` can be supplied after separately timing their
    preparation. They must describe this basis on these ordinary indices.
    The iteration starts from the projected noninteracting atomic ``1/z``.
    Each step mixes the bath equation with the current projected values and
    projects the mixed values through the same cached SVD. Convergence requires
    the actual equation defect, rather than merely a small projected step.
    """

    def __init__(
        self,
        basis: PhysicalBasis,
        *,
        indices: Any = None,
        hopping: float = 1.0,
        rcond: float = 1e-14,
        design: npt.ArrayLike | None = None,
        fit_plan: LeastSquaresPlan | None = None,
    ) -> None:
        self.beta = _positive_scalar("basis.beta", basis.beta)
        self.rank = _integer("basis.rank", basis.rank, 1)
        self.hopping = _positive_scalar("hopping", hopping)
        self._hopping_squared = self.hopping * self.hopping
        if not math.isfinite(self._hopping_squared) or self._hopping_squared == 0:
            raise ValueError("hopping squared must be finite and positive in float64")
        self.indices = cast(
            IndexArray,
            _immutable(_indices(solver_index_mesh() if indices is None else indices)),
        )
        self._z = 1j * matsubara_frequencies(self.indices, beta=self.beta)
        raw_design = basis.design_matsubara(self.indices) if design is None else design
        raw = np.asarray(raw_design)
        if raw.shape != (self.indices.size, self.rank) or raw.dtype.kind not in "iufc":
            raise ValueError("frequency design must have shape (indices, basis.rank)")
        self.design = cast(
            ComplexArray, _immutable(np.asarray(raw, dtype=np.complex128))
        )
        if not np.all(np.isfinite(self.design)):
            raise ValueError("frequency design must contain only finite values")
        if fit_plan is None:
            self.fit_plan = LeastSquaresPlan(self.design, rcond=rcond)
        else:
            if not isinstance(fit_plan, LeastSquaresPlan):
                raise ValueError("fit_plan must be a LeastSquaresPlan")
            if not np.array_equal(fit_plan.design, self.design):
                raise ValueError("fit_plan must factor the supplied frequency design")
            # A precomputed plan's own cutoff is authoritative and is reported.
            self.fit_plan = fit_plan

    def solve(
        self,
        *,
        mixing: float = 0.5,
        tolerance: float = 1e-9,
        max_iterations: int = 200,
    ) -> BetheSolverResult:
        if isinstance(mixing, (bool, np.bool_)) or not isinstance(mixing, Real):
            raise ValueError("mixing must be a finite real number in (0, 1]")
        mixing = float(mixing)
        if not math.isfinite(mixing) or mixing <= 0 or mixing > 1:
            raise ValueError("mixing must be a finite real number in (0, 1]")
        tolerance = _positive_scalar("tolerance", tolerance)
        limit = _integer("max_iterations", max_iterations, 1)
        initial = 1.0 / self._z
        coefficients = np.asarray(self.fit_plan.fit(initial), dtype=np.complex128)
        values = self.design @ coefficients
        initial_residual = float(np.max(np.abs(values - initial)))
        history: list[BetheIteration] = []
        converged = False
        stop_reason = "iteration_budget_exhausted"
        for iteration in range(1, limit + 1):
            with np.errstate(all="ignore"):
                target = 1.0 / (self._z - self._hopping_squared * values)
                mixed = (1.0 - mixing) * values + mixing * target
            if not np.all(np.isfinite(mixed)):
                stop_reason = "nonfinite_iteration"
                break
            try:
                next_coefficients = np.asarray(
                    self.fit_plan.fit(mixed), dtype=np.complex128
                )
            except NumericalBreakdownError:
                stop_reason = "nonfinite_projection"
                break
            with np.errstate(all="ignore"):
                projected = self.design @ next_coefficients
                defect = projected - 1.0 / (self._z - self._hopping_squared * projected)
            if not np.all(np.isfinite(projected)) or not np.all(np.isfinite(defect)):
                stop_reason = "nonfinite_iteration"
                break
            history.append(
                BetheIteration(
                    iteration=iteration,
                    step_norm=float(np.max(np.abs(projected - values))),
                    defect_norm=float(np.max(np.abs(defect))),
                    projection_residual_norm=float(np.max(np.abs(projected - mixed))),
                )
            )
            coefficients, values = next_coefficients, projected
            if history[-1].defect_norm <= tolerance:
                converged = True
                stop_reason = "self_consistency_tolerance_reached"
                break
        return BetheSolverResult(
            coefficients=cast(ComplexArray, _immutable(coefficients)),
            frequency_values=cast(ComplexArray, _immutable(values)),
            indices=self.indices,
            history=tuple(history),
            converged=converged,
            stop_reason=stop_reason,
            fit_diagnostics=self.fit_plan.diagnostics,
            beta=self.beta,
            hopping=self.hopping,
            mixing=mixing,
            tolerance=tolerance,
            max_iterations=limit,
            initial_projection_residual_norm=initial_residual,
        )
