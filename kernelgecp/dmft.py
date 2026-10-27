"""Projected half-filled single-band Bethe IPT with explicit numerical scopes.

The Hartree constant U/2 cancels the half-filled chemical potential analytically.
The supplied bases represent only decaying Green functions and dynamic self-
energy. Convergence of the discrete projected update is not physical accuracy
or an exact solution of the interacting Hubbard impurity problem.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from numbers import Real
from typing import Protocol, cast

import numpy as np
import numpy.typing as npt

from .basis import _indices_array
from .fitting import (
    LeastSquaresDiagnostics,
    LeastSquaresPlan,
    NumericArray,
    _immutable_copy,
    _numeric_array,
)
from .types import FloatArray, NumericalBreakdownError

type ComplexArray = npt.NDArray[np.complex128]
type IndexArray = npt.NDArray[np.int64]

_MAX_PLAN_ELEMENTS = 2**24


def _real(value: object, name: str, *, zero_allowed: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be finite in binary64") from error
    if not np.isfinite(result) or (result < 0 if zero_allowed else result <= 0):
        qualifier = "nonnegative" if zero_allowed else "positive"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return result


def _integer(value: object, name: str, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    try:
        result = operator.index(value)  # type: ignore[arg-type]
    except TypeError as error:
        raise ValueError(f"{name} must be an integer >= {minimum}") from error
    if result < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return result


def _complex_vector(value: npt.ArrayLike, name: str, size: int) -> ComplexArray:
    array = _numeric_array(value, name)
    if array.shape != (size,):
        raise ValueError(f"{name} must be a one-dimensional array of length {size}")
    return np.asarray(array, dtype=np.complex128)


def _freeze_complex(array: npt.ArrayLike) -> ComplexArray:
    return cast(ComplexArray, _immutable_copy(np.asarray(array, dtype=np.complex128)))


def _freeze_indices(array: IndexArray) -> IndexArray:
    return np.frombuffer(array.tobytes(), dtype=np.int64).reshape(array.shape)


def _checked(array: NumericArray, name: str) -> NumericArray:
    if not np.all(np.isfinite(array)):
        raise NumericalBreakdownError(f"nonfinite {name} in the IPT update")
    return array


def _maximum(array: NumericArray, name: str) -> float:
    with np.errstate(all="ignore"):
        value = float(np.max(np.abs(array), initial=0))
    if not np.isfinite(value):
        raise NumericalBreakdownError(f"nonfinite {name} norm in the IPT update")
    return value


@dataclass(frozen=True, slots=True)
class IPTConfig:
    """A bounded numerical half-filled paramagnetic Bethe IPT calculation."""

    beta: float = 10.0
    hopping: float = 1.0
    interaction: float = 2.0
    n_iw: int = 512
    n_tau: int = 6145
    mixing: float = 0.5
    max_iterations: int = 200
    tolerance: float = 1e-8
    rcond: float = 1e-13

    def __post_init__(self) -> None:
        for name in ("beta", "hopping", "mixing", "tolerance"):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        for name in ("interaction", "rcond"):
            object.__setattr__(
                self, name, _real(getattr(self, name), name, zero_allowed=True)
            )
        if self.mixing > 1:
            raise ValueError("mixing must lie in (0, 1]")
        for name, minimum in (("n_iw", 1), ("n_tau", 2), ("max_iterations", 1)):
            object.__setattr__(self, name, _integer(getattr(self, name), name, minimum))
        if not np.isfinite(self.hopping * self.hopping) or self.hopping**2 == 0:
            raise ValueError("hopping squared must be finite and positive")
        if not np.isfinite(self.interaction * self.interaction):
            raise ValueError("interaction squared must be finite")


_DEFAULT_CONFIG = IPTConfig()


class PhysicalBasis(Protocol):
    """Physical time/frequency designs sharing one coefficient coordinate."""

    @property
    def rank(self) -> int: ...

    @property
    def beta(self) -> float: ...

    @property
    def wmax(self) -> float: ...

    def design_tau(self, tau: npt.ArrayLike) -> NumericArray: ...

    def design_matsubara(self, indices: npt.ArrayLike) -> NumericArray: ...


@dataclass(frozen=True, slots=True)
class IPTDiagnostics:
    """Actual fitting diagnostics; none certifies unsampled physical error."""

    weiss_frequency_fit: LeastSquaresDiagnostics
    sigma_time_fit: LeastSquaresDiagnostics
    green_frequency_fit: LeastSquaresDiagnostics
    weiss_projection_residual: float
    sigma_projection_residual: float
    green_projection_residual: float
    arithmetic: str = "binary64_complex_svd_and_scalar_dyson"
    residual_scope: str = "supplied_frequency_observations_and_time_training_grid"
    hartree_convention: str = "mu_U_over_2_cancels_static_sigma_U_over_2"
    certification_status: str = "not_certified"


@dataclass(frozen=True, slots=True)
class IPTUpdate:
    """One candidate, computed at bath_iw before any bath mixing.

    g_iw is the unprojected scalar Dyson candidate. g_tau and sigma_tau are
    reconstructions from their fitted coefficients. sigma_tau_product retains
    the unprojected reflected product used as self-energy training data.
    """

    config: IPTConfig
    indices: IndexArray
    tau: FloatArray
    frequency_fit_indices: IndexArray
    bath_iw: ComplexArray
    weiss_iw: ComplexArray
    weiss_tau: ComplexArray
    weiss_reflected_tau: ComplexArray
    sigma_tau_product: ComplexArray
    sigma_tau: ComplexArray
    sigma_iw: ComplexArray
    g_iw: ComplexArray
    g_tau: ComplexArray
    weiss_coefficients: ComplexArray
    sigma_coefficients: ComplexArray
    g_coefficients: ComplexArray
    diagnostics: IPTDiagnostics
    unmixed_defect: float
    full_grid_unmixed_defect: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", _freeze_indices(self.indices))
        object.__setattr__(
            self, "frequency_fit_indices", _freeze_indices(self.frequency_fit_indices)
        )
        object.__setattr__(self, "tau", _immutable_copy(self.tau))
        for name in (
            "bath_iw",
            "weiss_iw",
            "weiss_tau",
            "weiss_reflected_tau",
            "sigma_tau_product",
            "sigma_tau",
            "sigma_iw",
            "g_iw",
            "g_tau",
            "weiss_coefficients",
            "sigma_coefficients",
            "g_coefficients",
        ):
            object.__setattr__(self, name, _freeze_complex(getattr(self, name)))


@dataclass(frozen=True, slots=True)
class IPTIteration:
    """Observed defect and proposed mixed-step norm, never a physics verdict.

    mixed_step is the norm of the proposed mixed bath minus its input bath.
    A terminal proposal is not applied after convergence or budget exhaustion.
    """

    iteration: int
    unmixed_defect: float
    mixed_step: float
    weiss_projection_residual: float | None = None
    sigma_projection_residual: float | None = None
    green_projection_residual: float | None = None
    full_grid_unmixed_defect: float | None = None


@dataclass(frozen=True, slots=True)
class IPTResult:
    """Last valid candidate, actual control history, and explicit stop status.

    bath_iw is the input that generated final_update, not a post-stop proposal.
    A failure before the first valid update leaves numerical output properties
    as None; it never substitutes the noninteracting solution or zero arrays.
    """

    config: IPTConfig
    indices: IndexArray
    tau: FloatArray
    frequency_fit_indices: IndexArray
    bath_iw: ComplexArray
    final_update: IPTUpdate | None
    history: tuple[IPTIteration, ...]
    converged: bool
    stop_reason: str
    bath_trajectory: tuple[ComplexArray, ...] = ()
    failure_detail: str | None = None
    failed_bath_iw: ComplexArray | None = None
    convergence_scope: str = "unprojected_bethe_defect_on_frequency_observations_only"
    certification_status: str = "not_certified"

    def __post_init__(self) -> None:
        object.__setattr__(self, "indices", _freeze_indices(self.indices))
        object.__setattr__(self, "tau", _immutable_copy(self.tau))
        object.__setattr__(
            self, "frequency_fit_indices", _freeze_indices(self.frequency_fit_indices)
        )
        object.__setattr__(self, "bath_iw", _freeze_complex(self.bath_iw))
        object.__setattr__(self, "history", tuple(self.history))
        object.__setattr__(
            self,
            "bath_trajectory",
            tuple(_freeze_complex(b) for b in self.bath_trajectory),
        )
        if self.failed_bath_iw is not None:
            object.__setattr__(
                self, "failed_bath_iw", _freeze_complex(self.failed_bath_iw)
            )

    @property
    def g_iw(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.g_iw

    @property
    def sigma_iw(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.sigma_iw

    @property
    def weiss_iw(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.weiss_iw

    @property
    def weiss_tau(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.weiss_tau

    @property
    def g_tau(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.g_tau

    @property
    def sigma_tau(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.sigma_tau

    @property
    def sigma_tau_product(self) -> ComplexArray | None:
        return (
            None if self.final_update is None else self.final_update.sigma_tau_product
        )

    @property
    def weiss_coefficients(self) -> ComplexArray | None:
        return (
            None if self.final_update is None else self.final_update.weiss_coefficients
        )

    @property
    def sigma_coefficients(self) -> ComplexArray | None:
        return (
            None if self.final_update is None else self.final_update.sigma_coefficients
        )

    @property
    def g_coefficients(self) -> ComplexArray | None:
        return None if self.final_update is None else self.final_update.g_coefficients

    @property
    def diagnostics(self) -> IPTDiagnostics | None:
        return None if self.final_update is None else self.final_update.diagnostics

    @property
    def unmixed_defect(self) -> float | None:
        return None if self.final_update is None else self.final_update.unmixed_defect

    @property
    def full_grid_unmixed_defect(self) -> float | None:
        return (
            None
            if self.final_update is None
            else self.final_update.full_grid_unmixed_defect
        )


class IPTPlan:
    """Cache two SVDs for scalar dynamic-self-energy Bethe IPT.

    The Green-frequency plan is shared by Weiss and candidate-G fits. The
    self-energy is fitted independently on the common uniform time grid.
    Its cutoff must be at least three times the Green cutoff for the reflected
    three-factor product; this is not a compact-support claim about the exact
    interacting model. No symmetry, tail or moment repair is performed.

    frequency_fit_indices selects unique ordinary indices from the full grid.
    Omitted modes participate in neither fitting nor stopping. Full-grid Dyson
    values and defects remain available as observations, not controller inputs.
    """

    def __init__(
        self,
        green_basis: PhysicalBasis,
        selfenergy_basis: PhysicalBasis,
        *,
        config: IPTConfig = _DEFAULT_CONFIG,
        frequency_fit_indices: npt.ArrayLike | None = None,
    ) -> None:
        if not isinstance(config, IPTConfig):
            raise TypeError("config must be IPTConfig")
        self._config = config
        self._green_basis = green_basis
        self._selfenergy_basis = selfenergy_basis
        ranks = []
        for name, basis in (("Green", green_basis), ("self-energy", selfenergy_basis)):
            try:
                rank = _integer(basis.rank, f"{name} rank", 1)
                beta = _real(basis.beta, f"{name} beta")
                _real(basis.wmax, f"{name} wmax")
            except AttributeError as error:
                raise TypeError(
                    "bases must implement the PhysicalBasis interface"
                ) from error
            if beta != config.beta:
                raise ValueError("basis beta must equal config.beta")
            ranks.append(rank)
        if selfenergy_basis.wmax < 3 * green_basis.wmax:
            raise ValueError("self-energy wmax must be at least three times Green wmax")
        supplied = (
            None
            if frequency_fit_indices is None
            else _indices_array(frequency_fit_indices)
        )
        if supplied is not None and (
            not supplied.size
            or np.unique(supplied).size != supplied.size
            or np.any(supplied < -config.n_iw)
            or np.any(supplied >= config.n_iw)
        ):
            raise ValueError(
                "frequency_fit_indices must be a nonempty unique grid subset"
            )
        observed_count = 2 * config.n_iw if supplied is None else supplied.size
        elements = ranks[0] * (observed_count + 2 * config.n_tau) + ranks[1] * (
            2 * config.n_iw + config.n_tau
        )
        if elements > _MAX_PLAN_ELEMENTS:
            raise ValueError("IPT designs exceed the 2**24-element preparation budget")
        self._indices = _freeze_indices(
            np.arange(-config.n_iw, config.n_iw, dtype=np.int64)
        )
        self._frequency_fit_indices = (
            self._indices if supplied is None else _freeze_indices(supplied)
        )
        self._positions = self._frequency_fit_indices + config.n_iw
        self._tau = cast(
            FloatArray, _immutable_copy(np.linspace(0, config.beta, config.n_tau))
        )
        with np.errstate(all="ignore"):
            frequency = (2.0 * self._indices + 1) * (np.pi / config.beta)
        if not np.all(np.isfinite(frequency)) or np.any(frequency == 0):
            raise ValueError(
                "physical Matsubara frequencies must be finite and nonzero"
            )
        self._z = 1j * frequency
        self._hopping_squared = config.hopping * config.hopping
        self._interaction_squared = config.interaction * config.interaction
        green_iw = self._design(
            green_basis.design_matsubara(self._frequency_fit_indices),
            observed_count,
            ranks[0],
            "Green frequency",
        )
        self._green_tau = self._design(
            green_basis.design_tau(self._tau), config.n_tau, ranks[0], "Green time"
        )
        self._green_reflected_tau = self._design(
            green_basis.design_tau(config.beta - self._tau),
            config.n_tau,
            ranks[0],
            "reflected Green time",
        )
        self._sigma_tau = self._design(
            selfenergy_basis.design_tau(self._tau),
            config.n_tau,
            ranks[1],
            "self-energy time",
        )
        self._sigma_iw = self._design(
            selfenergy_basis.design_matsubara(self._indices),
            2 * config.n_iw,
            ranks[1],
            "self-energy frequency",
        )
        self._green_fit = LeastSquaresPlan(green_iw, rcond=config.rcond)
        self._sigma_fit = LeastSquaresPlan(self._sigma_tau, rcond=config.rcond)

    @staticmethod
    def _design(value: npt.ArrayLike, rows: int, rank: int, name: str) -> NumericArray:
        array = _numeric_array(value, name + " design")
        if array.shape != (rows, rank):
            raise ValueError(f"{name} design must have shape ({rows}, {rank})")
        return _immutable_copy(array)

    @property
    def config(self) -> IPTConfig:
        return self._config

    @property
    def green_basis(self) -> PhysicalBasis:
        return self._green_basis

    @property
    def selfenergy_basis(self) -> PhysicalBasis:
        return self._selfenergy_basis

    @property
    def indices(self) -> IndexArray:
        return self._indices

    @property
    def tau(self) -> FloatArray:
        return self._tau

    @property
    def frequency_fit_indices(self) -> IndexArray:
        return self._frequency_fit_indices

    def update(self, bath_iw: npt.ArrayLike) -> IPTUpdate:
        """Return the actual candidate and all projection diagnostics.

        The algebraic update accepts asymmetric scalar test data too, enabling
        independent checks of the beta-tau reflection. It does not thereby
        implement a doped IPT model or change the half-filled Hartree convention.
        """
        bath = _complex_vector(bath_iw, "bath_iw", self._indices.size)
        with np.errstate(all="ignore"):
            weiss_iw = cast(
                ComplexArray, _checked(1 / (self._z - bath), "Weiss values")
            )
            weiss_coefficients = self._green_fit.fit(weiss_iw[self._positions])
            weiss_tau = _checked(
                self._green_tau @ weiss_coefficients, "Weiss time values"
            )
            reflected = _checked(
                self._green_reflected_tau @ weiss_coefficients, "reflected Weiss values"
            )
            product = (
                np.zeros_like(weiss_tau)
                if self._interaction_squared == 0
                else self._interaction_squared * weiss_tau**2 * reflected
            )
            _checked(product, "self-energy product")
            sigma_coefficients = self._sigma_fit.fit(product)
            sigma_tau = _checked(
                self._sigma_tau @ sigma_coefficients, "self-energy time"
            )
            sigma_iw = _checked(
                self._sigma_iw @ sigma_coefficients, "self-energy frequency"
            )
            g_iw = _checked(1 / (self._z - bath - sigma_iw), "Dyson candidate")
            g_coefficients = self._green_fit.fit(g_iw[self._positions])
            g_tau = _checked(self._green_tau @ g_coefficients, "candidate Green time")
            difference = _checked(self._hopping_squared * g_iw - bath, "bath defect")
        diagnostics = IPTDiagnostics(
            self._green_fit.diagnostics,
            self._sigma_fit.diagnostics,
            self._green_fit.diagnostics,
            _maximum(
                self._green_fit.residual(weiss_iw[self._positions], weiss_coefficients),
                "Weiss projection residual",
            ),
            _maximum(
                self._sigma_fit.residual(product, sigma_coefficients), "Sigma residual"
            ),
            _maximum(
                self._green_fit.residual(g_iw[self._positions], g_coefficients),
                "Green projection residual",
            ),
        )
        return IPTUpdate(
            self.config,
            self.indices,
            self.tau,
            self.frequency_fit_indices,
            bath,
            weiss_iw,
            cast(ComplexArray, weiss_tau),
            cast(ComplexArray, reflected),
            cast(ComplexArray, product),
            cast(ComplexArray, sigma_tau),
            cast(ComplexArray, sigma_iw),
            cast(ComplexArray, g_iw),
            cast(ComplexArray, g_tau),
            cast(ComplexArray, weiss_coefficients),
            cast(ComplexArray, sigma_coefficients),
            cast(ComplexArray, g_coefficients),
            diagnostics,
            _maximum(difference[self._positions], "observed bath defect"),
            _maximum(difference, "full-grid bath defect"),
        )

    def _initial_bath(self) -> ComplexArray:
        """Cancellation-free independent noninteracting semicircle solution."""
        frequency = self._z.imag
        scale = np.maximum(np.abs(frequency), 2 * self.config.hopping)
        with np.errstate(all="ignore"):
            normalized = np.abs(frequency) / scale
            root = np.hypot(normalized, 2 * self.config.hopping / scale)
            magnitude = (2 / (normalized + root)) / scale
            bath = -1j * np.sign(frequency) * self._hopping_squared * magnitude
        return cast(ComplexArray, _checked(bath, "initial semicircle bath"))

    def solve(
        self,
        initial_bath: npt.ArrayLike | None = None,
        *,
        record_trajectory: bool = False,
    ) -> IPTResult:
        """Iterate until two consecutive unprojected observed defects pass.

        Invalid or nonfinite caller input raises ValueError. Iterative numerical
        breakdown returns a failed result retaining the last valid candidate,
        history and failed bath. Only nonterminal proposals advance the bath.
        """
        if type(record_trajectory) is not bool:
            raise ValueError("record_trajectory must be a boolean")
        bath = (
            self._initial_bath()
            if initial_bath is None
            else _complex_vector(
                initial_bath, "initial_bath", self._indices.size
            ).copy()
        )
        history: list[IPTIteration] = []
        trajectory: list[ComplexArray] = []
        final: IPTUpdate | None = None
        converged = False
        stop_reason = "iteration_budget_exhausted"
        failure_detail = None
        failed_bath = None
        consecutive = 0
        for iteration in range(1, self.config.max_iterations + 1):
            try:
                update = self.update(bath)
                with np.errstate(all="ignore"):
                    proposed = cast(
                        ComplexArray,
                        _checked(
                            (1 - self.config.mixing) * bath
                            + self.config.mixing * self._hopping_squared * update.g_iw,
                            "mixed bath proposal",
                        ),
                    )
                    mixed_step = _maximum(proposed - bath, "mixed bath step")
                diagnostic = update.diagnostics
                history.append(
                    IPTIteration(
                        iteration,
                        update.unmixed_defect,
                        mixed_step,
                        diagnostic.weiss_projection_residual,
                        diagnostic.sigma_projection_residual,
                        diagnostic.green_projection_residual,
                        update.full_grid_unmixed_defect,
                    )
                )
                if record_trajectory:
                    trajectory.append(update.bath_iw)
                final = update
                consecutive = (
                    consecutive + 1
                    if update.unmixed_defect < self.config.tolerance
                    else 0
                )
                if consecutive >= 2:
                    converged = True
                    stop_reason = "self_consistency_tolerance_reached"
                    break
                if iteration < self.config.max_iterations:
                    bath = proposed
            except (
                NumericalBreakdownError,
                np.linalg.LinAlgError,
                FloatingPointError,
            ) as error:
                stop_reason = "numerical_breakdown"
                failure_detail = (
                    f"iteration {iteration}: {type(error).__name__}: {error}"
                )
                failed_bath = bath
                break
        return IPTResult(
            config=self.config,
            indices=self.indices,
            tau=self.tau,
            frequency_fit_indices=self.frequency_fit_indices,
            bath_iw=final.bath_iw if final is not None else bath,
            final_update=final,
            history=tuple(history),
            converged=converged,
            stop_reason=stop_reason,
            bath_trajectory=tuple(trajectory),
            failure_detail=failure_detail,
            failed_bath_iw=failed_bath,
        )
