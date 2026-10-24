"""Reusable fermionic cardinal bases in explicit physical coordinates.

Construction stops on a finite grid; ``tol`` is not a promise about arbitrary
input Green functions or unsampled points. Reference arithmetic applies to
kernel evaluation and the cardinal solve, but exported designs and subsequent
least-squares fits use binary64. No continuous or rounding certificate is made.
"""

from __future__ import annotations

import json
import operator
from collections import OrderedDict
from dataclasses import asdict, dataclass
from numbers import Real
from os import PathLike
from typing import Literal, cast

import numpy as np
import numpy.typing as npt

from ._basis_io import read_archive, write_archive
from ._basis_numerics import BasisNumerics
from .fitting import (
    LeastSquaresDiagnostics,
    LeastSquaresPlan,
    NumericArray,
    _contains_boolean_sequence,
    _immutable_copy,
    _normalize_axis,
    _numeric_array,
)
from .gecp import gecp
from .kernels import FermionicKernel
from .types import FloatArray, GECPConfig, GECPResult, NumericalBreakdownError

type PrecisionPolicy = Literal["auto", "float64", "reference"]
type ArchivePath = str | PathLike[str]

_CONVENTIONS: dict[str, object] = {
    "statistics": "fermionic",
    "tau_domain": "closed_[0,beta]",
    "dimensionless_coordinates": "t=tau/beta;Omega=beta*epsilon",
    "kernel": "positive_exp(-t*Omega)/(1+exp(-Omega))",
    "green_function": "negative_kernel_spectral_integral",
    "fourier": "integral_0^beta_exp(+i*wn*tau)_G(tau)_dtau",
    "matsubara_indices": "ordinary_integer_n;wn=(2*n+1)*pi/beta",
    "coefficients": "cardinal_values_at_selected_physical_tau_nodes",
}
_MAX_RESTORE_RANK = 256
_MAX_RESTORE_WORK = 2**24  # rank squared times working-precision bits


@dataclass(frozen=True, slots=True)
class ConstructionAttempt:
    """One explicit construction attempt, including a failed fast path."""

    precision_bits: int
    outcome: str
    detail: str


@dataclass(frozen=True, slots=True)
class BasisDiagnostics:
    """Construction evidence, arithmetic stages and limits of the result."""

    rank: int
    requested_policy: PrecisionPolicy
    reference_precision_bits: int
    selection_precision_bits: int
    evaluation_precision_bits: int
    coordinate_precision_bits: int
    export_precision_bits: int
    fit_precision_bits: int
    sampled_residual: float
    converged: bool
    stop_reason: str
    tie_counts: tuple[int, ...]
    attempts: tuple[ConstructionAttempt, ...]
    exported_core_condition: float | None
    nodal_error: float
    history_status: str = "constructed_in_this_process"
    convergence_scope: str = "finite_construction_grid_only"
    conditioning_scope: str = "binary64_exported_core_2norm_estimate"
    certification_status: str = "not_certified"
    certification_scope: str = "no_continuous_input_or_rounding_certificate"


@dataclass(frozen=True, slots=True)
class FitDiagnostics:
    """Training evidence only, not an independent accuracy estimate."""

    observation_axis: int
    training_residual_max: float
    least_squares: LeastSquaresDiagnostics | None
    evidence_status: str = "computed_on_supplied_training_samples"
    certification_status: str = "not_certified"


def _positive_real(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be finite in binary64") from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive real number")
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


def _tau_array(tau: npt.ArrayLike, beta: float) -> FloatArray:
    values = _numeric_array(tau, "tau")
    if np.iscomplexobj(values) or values.ndim != 1:
        raise ValueError("tau must be a one-dimensional real array")
    if np.any(values < 0) or np.any(values > beta):
        raise ValueError("tau must lie in the closed interval [0, beta]")
    return cast(FloatArray, values)


def _indices_array(indices: npt.ArrayLike) -> npt.NDArray[np.int64]:
    if _contains_boolean_sequence(indices):
        raise ValueError("Matsubara indices must not contain booleans")
    values = np.asarray(indices)
    if values.ndim == 1 and values.size == 0:
        return np.array([], dtype=np.int64)
    if values.ndim != 1 or values.dtype.kind not in "iu":
        raise ValueError("Matsubara indices must be a one-dimensional integer array")
    # At this bound odd integer 2*n+1 is still exactly representable in binary64.
    limit = 2**51 - 1
    if np.any(values > limit) or (values.dtype.kind == "i" and np.any(values < -limit)):
        raise ValueError("Matsubara index magnitude must be smaller than 2**51")
    return np.asarray(values, dtype=np.int64)


def _cache_key(
    kind: str, values: npt.NDArray[np.generic], rcond: float | None = None
) -> tuple[str, str, bytes, float | None]:
    return kind, values.dtype.str, values.tobytes(), rcond


def _mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or not all(isinstance(k, str) for k in value):
        raise ValueError(f"invalid archive {name}")
    return cast(dict[str, object], value)


def _json_record(value: object) -> object:
    """Convert dataclass tuple fields to the archive's strict JSON tree."""
    return cast(object, json.loads(json.dumps(value, allow_nan=False)))


class FermionicBasis:
    """Select and cache a reusable physical imaginary-time basis.

    ``beta`` is inverse temperature, ``wmax`` is the spectral energy cutoff and
    ``beta*wmax >= 1``. The new facade uses finite-grid GECP explicitly; it does
    not alter the adaptive/reference defaults of the research ``gecp`` API.

    ``precision='float64'`` is a strict fast path. ``'reference'`` constructs,
    evaluates kernels and solves the selected core at ``precision_bits``.
    ``'auto'`` first tries float64 construction, then retries reference only on
    ``NumericalBreakdownError``. The chosen application arithmetic matches the
    successful construction. All public arrays and least-squares fits are
    float64/complex128. Inspect ``diagnostics`` and validate your own workload.

    Caches belong to this instance, are bounded by ``cache_size`` entries each,
    and own their coordinates. Repeated fitting on one grid reuses one SVD.
    Instances are not intended for concurrent cache mutation across threads.
    """

    def __init__(
        self,
        beta: float,
        wmax: float,
        tol: float = 1e-8,
        *,
        precision: PrecisionPolicy = "auto",
        precision_bits: int = 128,
        grid_order: int = 12,
        max_rank: int = 256,
        cache_size: int = 8,
    ) -> None:
        self._set_parameters(
            beta, wmax, tol, precision, precision_bits, grid_order, max_rank, cache_size
        )
        attempts: list[ConstructionAttempt] = []
        chosen_bits = 53 if precision != "reference" else self._reference_bits
        try:
            result = self._select(chosen_bits)
        except NumericalBreakdownError as error:
            attempts.append(
                ConstructionAttempt(chosen_bits, "arithmetic_failure", str(error))
            )
            if precision != "auto":
                raise
            chosen_bits = self._reference_bits
            result = self._select(chosen_bits)
        attempts.append(ConstructionAttempt(chosen_bits, result.stop_reason, ""))
        if not result.converged or result.rank == 0:
            raise NumericalBreakdownError(
                f"basis construction did not converge: {result.stop_reason}, "
                f"rank={result.rank}"
            )
        with np.errstate(all="ignore"):
            tau_nodes = np.asarray(result.t_nodes * self.beta, dtype=np.float64)
        self._set_nodes(tau_nodes, result.omega_nodes, chosen_bits)
        self._diagnostics = BasisDiagnostics(
            rank=self.rank,
            requested_policy=precision,
            reference_precision_bits=self._reference_bits,
            selection_precision_bits=chosen_bits,
            evaluation_precision_bits=chosen_bits,
            coordinate_precision_bits=53,
            export_precision_bits=53,
            fit_precision_bits=53,
            sampled_residual=float(result.residual_history[-1]),
            converged=result.converged,
            stop_reason=result.stop_reason,
            tie_counts=tuple(result.tie_counts),
            attempts=tuple(attempts),
            exported_core_condition=self._numerics.core_condition,
            nodal_error=self._numerics.nodal_error,
        )

    def _set_parameters(
        self,
        beta: object,
        wmax: object,
        tol: object,
        precision: object,
        precision_bits: object,
        grid_order: object,
        max_rank: object,
        cache_size: object,
    ) -> None:
        self._beta = _positive_real(beta, "beta")
        self._wmax = _positive_real(wmax, "wmax")
        self._tol = _positive_real(tol, "tol")
        if self._tol >= 1:
            raise ValueError("tol must be smaller than 1")
        with np.errstate(all="ignore"):
            self._cutoff = self._beta * self._wmax
        if not np.isfinite(self._cutoff) or self._cutoff < 1:
            raise ValueError("beta*wmax must be finite and >= 1")
        if not isinstance(precision, str) or precision not in (
            "auto",
            "float64",
            "reference",
        ):
            raise ValueError("precision must be 'auto', 'float64' or 'reference'")
        self._precision = cast(PrecisionPolicy, precision)
        self._reference_bits = _integer(precision_bits, "precision_bits", 64)
        self._grid_order = _integer(grid_order, "grid_order", 2)
        self._max_rank = _integer(max_rank, "max_rank", 1)
        self._cache_size = _integer(cache_size, "cache_size", 1)
        if self._reference_bits > 4096:
            raise ValueError("precision_bits must be <= 4096")
        self._design_cache: OrderedDict[
            tuple[str, str, bytes, float | None], NumericArray
        ] = OrderedDict()
        self._plan_cache: OrderedDict[
            tuple[str, str, bytes, float | None], LeastSquaresPlan
        ] = OrderedDict()

    def _select(self, bits: int) -> GECPResult:
        # Tiny positive kernel tails may round to zero in a binary64 sample.
        # The research update's own stricter traps still fail closed on unsafe
        # arithmetic. Do not inherit or change the caller's global trap state.
        with np.errstate(under="ignore"):
            return gecp(
                FermionicKernel(self.cutoff),
                config=GECPConfig(
                    tol=self._tol,
                    max_rank=self._max_rank,
                    pivot="grid",
                    precision_bits=bits,
                    grid_order=self._grid_order,
                ),
            )

    def _set_nodes(
        self, tau_nodes: npt.ArrayLike, omega_nodes: npt.ArrayLike, bits: int
    ) -> None:
        tau = _tau_array(tau_nodes, self.beta)
        omega = _numeric_array(omega_nodes, "omega_nodes")
        if (
            tau.size == 0
            or tau.size > self._max_rank
            or omega.ndim != 1
            or np.iscomplexobj(omega)
            or omega.size != tau.size
            or np.unique(tau).size != tau.size
            or np.unique(omega).size != omega.size
            or np.any(np.abs(omega) > self.cutoff)
        ):
            raise ValueError("invalid, repeated, or out-of-domain selected nodes")
        self._tau_nodes = cast(FloatArray, _immutable_copy(tau))
        self._omega_nodes = cast(FloatArray, _immutable_copy(omega))
        self._numerics = BasisNumerics(
            self.beta,
            self.cutoff,
            self._tau_nodes,
            self._omega_nodes,
            precision_bits=bits,
        )

    @property
    def beta(self) -> float:
        """Inverse temperature in the caller's physical units."""
        return self._beta

    @property
    def wmax(self) -> float:
        """Spectral energy cutoff in the reciprocal units of beta."""
        return self._wmax

    @property
    def cutoff(self) -> float:
        """Dimensionless cutoff beta*wmax."""
        return self._cutoff

    @property
    def rank(self) -> int:
        """Number of selected cardinal coefficients."""
        return int(self._tau_nodes.size)

    @property
    def tau_nodes(self) -> FloatArray:
        """Immutable physical imaginary-time nodes, in selection order."""
        return self._tau_nodes

    @property
    def omega_nodes(self) -> FloatArray:
        """Immutable dimensionless spectral nodes Omega=beta*epsilon."""
        return self._omega_nodes

    @property
    def diagnostics(self) -> BasisDiagnostics:
        """Immutable construction, precision and conditioning evidence."""
        return self._diagnostics

    @property
    def cache_info(self) -> dict[str, int]:
        """Current entry counts; caches are bounded separately, not by bytes."""
        return {
            "design_entries": len(self._design_cache),
            "fit_plan_entries": len(self._plan_cache),
            "capacity_each": self._cache_size,
        }

    def clear_cache(self) -> None:
        """Release cached designs and fitting plans, retaining the core factor."""
        self._design_cache.clear()
        self._plan_cache.clear()

    def _design(self, kind: str, values: npt.NDArray[np.generic]) -> NumericArray:
        key = _cache_key(kind, values)
        if key in self._design_cache:
            self._design_cache.move_to_end(key)
            return self._design_cache[key]
        result = _immutable_copy(
            self._numerics.design_tau(values)
            if kind == "tau"
            else self._numerics.design_matsubara(values)
        )
        self._design_cache[key] = result
        while len(self._design_cache) > self._cache_size:
            self._design_cache.popitem(last=False)
        return result

    def design_tau(self, tau: npt.ArrayLike) -> NumericArray:
        """Return immutable cardinal rows at physical tau in [0,beta]."""
        return self._design("tau", _tau_array(tau, self.beta))

    def design_matsubara(self, indices: npt.ArrayLike) -> NumericArray:
        """Return physical Fourier rows for ordinary integer fermionic n.

        The transform is integral exp(+i*wn*tau)*G(tau) d tau with
        wn=(2*n+1)*pi/beta. Indices are not the reduced odd integers 2*n+1.
        """
        return self._design("matsubara", _indices_array(indices))

    def prepare_tau(
        self, tau: npt.ArrayLike | None = None, *, rcond: float | None = None
    ) -> LeastSquaresPlan:
        """Prepare or reuse one SVD for a physical observation grid."""
        nodes = self.tau_nodes if tau is None else _tau_array(tau, self.beta)
        # Validate rcond even on a cache hit; bool hashes like zero/one.
        if rcond is not None:
            if isinstance(rcond, (bool, np.bool_)) or not isinstance(rcond, Real):
                raise ValueError("rcond must be a finite nonnegative real number")
            try:
                rcond = float(rcond)
            except OverflowError as error:
                raise ValueError("rcond must be finite in binary64") from error
            if not np.isfinite(rcond) or rcond < 0:
                raise ValueError("rcond must be a finite nonnegative real number")
        key = _cache_key("tau", nodes, rcond)
        if key in self._plan_cache:
            self._plan_cache.move_to_end(key)
            return self._plan_cache[key]
        plan = LeastSquaresPlan(self.design_tau(nodes), rcond=rcond)
        self._plan_cache[key] = plan
        while len(self._plan_cache) > self._cache_size:
            self._plan_cache.popitem(last=False)
        return plan

    def fit_tau(
        self,
        values: npt.ArrayLike,
        tau: npt.ArrayLike | None = None,
        *,
        axis: int = 0,
        rcond: float | None = None,
    ) -> FermionicFit:
        """Fit real/complex scalar or tensor observations with a cached SVD.

        The observation axis is replaced by the coefficient axis. All other
        axes are preserved. Supplying tau=None means selected tau_nodes.
        """
        plan = self.prepare_tau(tau, rcond=rcond)
        coefficients = plan.fit(values, axis=axis)
        residual = plan.residual(values, coefficients, axis=axis)
        normalized_axis = _normalize_axis(axis, coefficients.ndim)
        return FermionicFit(
            self,
            coefficients,
            axis=normalized_axis,
            diagnostics=FitDiagnostics(
                observation_axis=normalized_axis,
                training_residual_max=float(np.max(np.abs(residual), initial=0)),
                least_squares=plan.diagnostics,
            ),
        )

    def _evaluate(
        self, design: NumericArray, coefficients: npt.ArrayLike, axis: int
    ) -> NumericArray:
        array = _numeric_array(coefficients, "coefficients")
        if array.ndim == 0:
            raise ValueError("coefficients must have a coefficient axis")
        normalized_axis = _normalize_axis(axis, array.ndim)
        if array.shape[normalized_axis] != self.rank:
            raise ValueError("coefficient axis must equal basis rank")
        moved = np.moveaxis(array, normalized_axis, 0)
        with np.errstate(all="ignore"):
            result = design @ moved.reshape(self.rank, -1)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("basis evaluation is nonfinite in binary64")
        return cast(
            NumericArray,
            np.moveaxis(
                result.reshape((design.shape[0], *moved.shape[1:])), 0, normalized_axis
            ),
        )

    def evaluate_tau(
        self, coefficients: npt.ArrayLike, tau: npt.ArrayLike, *, axis: int = 0
    ) -> NumericArray:
        """Evaluate coefficient batches, replacing their rank axis by tau."""
        return self._evaluate(self.design_tau(tau), coefficients, axis)

    def evaluate_matsubara(
        self, coefficients: npt.ArrayLike, indices: npt.ArrayLike, *, axis: int = 0
    ) -> NumericArray:
        """Evaluate coefficient batches at ordinary integer fermionic n."""
        return self._evaluate(self.design_matsubara(indices), coefficients, axis)

    def _payload(
        self, kind: str
    ) -> tuple[dict[str, object], dict[str, npt.NDArray[np.generic]]]:
        metadata: dict[str, object] = {
            "format": "kernelgecp.fermionic_basis",
            "version": 1,
            "kind": kind,
            "conventions": dict(_CONVENTIONS),
            "parameters": {
                "beta": self.beta,
                "wmax": self.wmax,
                "tol": self._tol,
                "precision": self._precision,
                "precision_bits": self._reference_bits,
                "grid_order": self._grid_order,
                "max_rank": self._max_rank,
                "cache_size": self._cache_size,
            },
            "construction": _json_record(asdict(self.diagnostics)),
        }
        return metadata, {"tau_nodes": self.tau_nodes, "omega_nodes": self.omega_nodes}

    def save(self, path: ArchivePath) -> None:
        """Atomically write a versioned NPZ/JSON archive, never pickle.

        Caches are not persisted. Loading validates nodes and rebuilds factors
        without selecting new pivots. Historical stopping evidence is retained
        as unverified history, not promoted to a mathematical certificate.
        """
        write_archive(path, *self._payload("basis"))

    @classmethod
    def _restore(
        cls, metadata: dict[str, object], arrays: dict[str, npt.NDArray[np.generic]]
    ) -> FermionicBasis:
        expected_envelope = {
            "format",
            "version",
            "kind",
            "conventions",
            "parameters",
            "construction",
        }
        if metadata.get("kind") == "fit":
            expected_envelope.add("fit")
        if (
            set(metadata) != expected_envelope
            or metadata.get("kind") not in ("basis", "fit")
            or metadata.get("format") != "kernelgecp.fermionic_basis"
            or type(metadata.get("version")) is not int
            or metadata["version"] != 1
            or metadata.get("conventions") != _CONVENTIONS
        ):
            raise ValueError("unsupported or incompatible fermionic basis archive")
        parameters = _mapping(metadata.get("parameters"), "parameters")
        expected = {
            "beta",
            "wmax",
            "tol",
            "precision",
            "precision_bits",
            "grid_order",
            "max_rank",
            "cache_size",
        }
        if set(parameters) != expected:
            raise ValueError("invalid archive parameters")
        result = cls.__new__(cls)
        result._set_parameters(**parameters)
        history = _mapping(metadata.get("construction"), "construction")
        bits = _integer(
            history.get("evaluation_precision_bits"), "evaluation_precision_bits", 53
        )
        if bits not in (53, result._reference_bits):
            raise ValueError("archive arithmetic disagrees with precision policy")
        if result._precision == "float64" and bits != 53:
            raise ValueError(
                "archive fast-path arithmetic disagrees with precision policy"
            )
        if result._precision == "reference" and bits != result._reference_bits:
            raise ValueError(
                "archive reference arithmetic disagrees with precision policy"
            )
        if set(arrays) not in (
            {"tau_nodes", "omega_nodes"},
            {"tau_nodes", "omega_nodes", "coefficients"},
        ):
            raise ValueError("invalid archive array payload")
        rank = arrays["tau_nodes"].size
        if rank > _MAX_RESTORE_RANK or rank * rank * bits > _MAX_RESTORE_WORK:
            raise ValueError(
                "archive selected core exceeds the restoration work budget"
            )
        expected_history: dict[str, object] = {
            "rank": rank,
            "requested_policy": result._precision,
            "reference_precision_bits": result._reference_bits,
            "selection_precision_bits": bits,
            "coordinate_precision_bits": 53,
            "export_precision_bits": 53,
            "fit_precision_bits": 53,
        }
        for key, expected_value in expected_history.items():
            if (
                type(history.get(key)) is not type(expected_value)
                or history[key] != expected_value
            ):
                raise ValueError(f"inconsistent archive construction field: {key}")
        result._set_nodes(arrays["tau_nodes"], arrays["omega_nodes"], bits)
        # Loaded histories are not trusted as proof, including any hand-edited
        # certification field. Fixed scopes always come from this code.
        residual = history.get("sampled_residual")
        if (
            isinstance(residual, bool)
            or not isinstance(residual, Real)
            or not np.isfinite(float(residual))
            or float(residual) < 0
        ):
            raise ValueError("invalid archive sampled residual")
        ties = history.get("tie_counts")
        if not isinstance(ties, list) or len(ties) != result.rank:
            raise ValueError("invalid archive tie history")
        tie_counts = tuple(_integer(item, "tie_count", 0) for item in ties)
        attempts_value = history.get("attempts")
        if not isinstance(attempts_value, list) or not 1 <= len(attempts_value) <= 2:
            raise ValueError("invalid archive construction attempts")
        attempts: list[ConstructionAttempt] = []
        for value in attempts_value:
            attempt = _mapping(value, "attempt")
            if (
                set(attempt) != {"precision_bits", "outcome", "detail"}
                or not isinstance(attempt["outcome"], str)
                or not isinstance(attempt["detail"], str)
            ):
                raise ValueError("invalid archive construction attempt")
            attempts.append(
                ConstructionAttempt(
                    _integer(attempt["precision_bits"], "attempt precision", 53),
                    attempt["outcome"],
                    attempt["detail"],
                )
            )
        stop_reason = history.get("stop_reason")
        if (
            not isinstance(stop_reason, str)
            or history.get("converged") is not True
            or float(residual) > result._tol
        ):
            raise ValueError("archive does not report converged sampled construction")
        if attempts[-1] != ConstructionAttempt(bits, stop_reason, ""):
            raise ValueError("archive final construction attempt is inconsistent")
        if len(attempts) == 2 and (
            result._precision != "auto"
            or bits != result._reference_bits
            or attempts[0].precision_bits != 53
            or attempts[0].outcome != "arithmetic_failure"
        ):
            raise ValueError("archive fallback construction history is inconsistent")
        result._diagnostics = BasisDiagnostics(
            rank=result.rank,
            requested_policy=result._precision,
            reference_precision_bits=result._reference_bits,
            selection_precision_bits=bits,
            evaluation_precision_bits=bits,
            coordinate_precision_bits=53,
            export_precision_bits=53,
            fit_precision_bits=53,
            sampled_residual=float(residual),
            converged=True,
            stop_reason=stop_reason,
            tie_counts=tie_counts,
            attempts=tuple(attempts),
            exported_core_condition=result._numerics.core_condition,
            nodal_error=result._numerics.nodal_error,
            history_status="restored_unverified_history",
        )
        return result

    @classmethod
    def load(cls, path: ArchivePath) -> FermionicBasis:
        """Validate and restore a basis without recomputing pivot selection."""
        metadata, arrays = read_archive(path)
        if metadata.get("kind") != "basis" or set(arrays) != {
            "tau_nodes",
            "omega_nodes",
        }:
            raise ValueError("archive does not contain a standalone basis")
        return cls._restore(metadata, arrays)


class FermionicFit:
    """Immutable coefficient batches associated with a particular basis."""

    def __init__(
        self,
        basis: FermionicBasis,
        coefficients: npt.ArrayLike,
        *,
        axis: int = 0,
        diagnostics: FitDiagnostics | None = None,
    ) -> None:
        array = _numeric_array(coefficients, "coefficients")
        if array.ndim == 0:
            raise ValueError("coefficients must have a coefficient axis")
        self._axis = _normalize_axis(axis, array.ndim)
        if array.shape[self._axis] != basis.rank:
            raise ValueError("coefficient axis must equal basis rank")
        self._basis = basis
        self._coefficients = _immutable_copy(array)
        self._diagnostics = diagnostics

    @property
    def basis(self) -> FermionicBasis:
        """Basis owning the physical convention and core factor."""
        return self._basis

    @property
    def coefficients(self) -> NumericArray:
        """Owned immutable coefficient batches in cardinal coordinates."""
        return self._coefficients

    @property
    def axis(self) -> int:
        """Coefficient axis, preserved as the query axis in evaluations."""
        return self._axis

    @property
    def diagnostics(self) -> FitDiagnostics | None:
        """Training diagnostics, absent for manually supplied coefficients."""
        return self._diagnostics

    def evaluate_tau(self, tau: npt.ArrayLike) -> NumericArray:
        """Evaluate physical imaginary time with the stored batch layout."""
        return self.basis.evaluate_tau(self.coefficients, tau, axis=self.axis)

    def evaluate_matsubara(self, indices: npt.ArrayLike) -> NumericArray:
        """Evaluate ordinary integer fermionic Matsubara indices."""
        return self.basis.evaluate_matsubara(self.coefficients, indices, axis=self.axis)

    def save(self, path: ArchivePath) -> None:
        """Persist the basis and complex coefficients without executable data."""
        metadata, arrays = self.basis._payload("fit")
        metadata["fit"] = {
            "axis": self.axis,
            "diagnostics": (
                _json_record(asdict(self.diagnostics)) if self.diagnostics else None
            ),
        }
        arrays["coefficients"] = self.coefficients
        write_archive(path, metadata, arrays)

    @classmethod
    def load(cls, path: ArchivePath) -> FermionicFit:
        """Restore coefficients and factors, without rerunning GECP or SVD.

        Historical training diagnostics are not independently revalidated and
        are intentionally not restored as current evidence.
        """
        metadata, arrays = read_archive(path)
        if metadata.get("kind") != "fit" or "coefficients" not in arrays:
            raise ValueError("archive does not contain a fitted function")
        fit = _mapping(metadata.get("fit"), "fit")
        if set(fit) != {"axis", "diagnostics"}:
            raise ValueError("invalid fitted-function archive envelope")
        axis = _integer(fit.get("axis"), "axis", 0)
        basis = FermionicBasis._restore(metadata, arrays)
        return cls(basis, arrays["coefficients"], axis=axis)
