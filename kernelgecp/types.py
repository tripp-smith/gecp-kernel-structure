"""Typed configuration and result records shared by numerical algorithms."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


class NumericalBreakdownError(RuntimeError):
    """Raised when an algorithm cannot safely continue above its tolerance."""


@dataclass(frozen=True, slots=True)
class GECPConfig:
    """Configuration for deterministic finite-grid or adaptive GECP."""

    tol: float = 1e-10
    max_rank: int = 256
    pivot: Literal["grid", "adaptive"] = "adaptive"
    precision_bits: int = 128
    grid_order: int = 24
    certificate_abs_tol: float = 1e-12
    certificate_rel_tol: float = 1e-8
    max_cells: int = 100_000
    tie_rtol: float | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.tol) or self.tol <= 0:
            raise ValueError("tol must be finite and positive")
        if self.max_rank < 0:
            raise ValueError("max_rank must be nonnegative")
        if self.precision_bits < 53:
            raise ValueError("precision_bits must be at least 53")
        if self.grid_order < 2:
            raise ValueError("grid_order must be at least 2")
        if any(
            not np.isfinite(value) or value < 0
            for value in (self.certificate_abs_tol, self.certificate_rel_tol)
        ):
            raise ValueError("certificate tolerances must be finite and nonnegative")
        if self.max_cells < 1:
            raise ValueError("max_cells must be positive")
        if self.tie_rtol is not None and (
            not np.isfinite(self.tie_rtol) or self.tie_rtol < 0
        ):
            raise ValueError("tie_rtol must be finite and nonnegative")

    def normalized(self) -> dict[str, Any]:
        """Return a deterministic JSON-compatible configuration mapping."""

        return asdict(self)


DEFAULT_GECP_CONFIG = GECPConfig()


@dataclass(frozen=True, slots=True)
class PivotCertificate:
    """Bounds for the caller's interval objective, with an approximate point value.

    ``value`` is an ordinary numerical evaluation; only ``lower_bound`` and
    ``upper_bound`` are enclosure-derived bounds in interval mode. Candidate
    tie counts concern the distinct sampled coordinates, not all continuous
    maximizers. ``objective_scope`` identifies the function being enclosed.
    """

    t: float
    omega: float
    value: float
    lower_bound: float
    upper_bound: float
    eta: float
    cells_evaluated: int
    certified: bool
    termination_reason: str
    candidate_count: int = 0
    exact_tie_count: int = 0
    near_tie_count: int = 0
    tie_rtol: float | None = None
    tie_scope: str = "not_recorded"
    objective_scope: str = "caller_supplied_objective"
    arithmetic: str = "not_recorded"
    precision_bits: int | None = None
    requested_precision_bits: int | None = None
    enclosure_precision_bits: int | None = None
    bound_precision_bits: int | None = None
    coordinate_precision_bits: int | None = None

    @property
    def metadata(self) -> NumericalMetadata:
        """Expose scalar-search evidence without inventing matrix conditioning."""
        interval = self.arithmetic == "mpmath_interval_with_float64_exports"
        return NumericalMetadata(
            algorithm="interval_certified_pivot"
            if interval
            else "lipschitz_pivot"
            if self.precision_bits is not None
            else "not_recorded",
            arithmetic=self.arithmetic,
            precision_bits=self.precision_bits,
            requested_precision_bits=self.requested_precision_bits,
            precision_status="mixed_precision_explicit_in_configuration"
            if self.precision_bits is not None
            else "not_recorded",
            conditioning_status="not_applicable",
            conditioning_method="scalar_objective_no_matrix_factorization",
            tie_status="recorded_sampled_candidates"
            if self.candidate_count
            else "not_recorded",
            tie_policy="lexicographic_sampled_singleton_lower_bound"
            if interval
            else "lexicographic_sampled_absolute_point_value",
            tie_rtol=self.tie_rtol,
            tie_counts=[self.near_tie_count] if self.candidate_count else [],
            tie_tolerance_rule=self.tie_scope,
            certification_status="certified" if self.certified else "not_certified",
            certification_scope=self.objective_scope,
            converged=self.certified,
            stop_reason=self.termination_reason,
            convergence_scope="declared_objective_pivot_gap_under_callback_obligations",
            configuration={
                "enclosure_precision_bits": self.enclosure_precision_bits,
                "enclosure_arithmetic": "mpmath_interval"
                if interval
                else "exact_binary_rational"
                if self.precision_bits is not None
                else "not_recorded",
                "bound_precision_bits": self.bound_precision_bits,
                "value_precision_bits": self.bound_precision_bits,
                "coordinate_precision_bits": self.coordinate_precision_bits,
                "point_precision_bits": self.precision_bits,
                "candidate_count": self.candidate_count,
                "exact_tie_count": self.exact_tie_count,
                "near_tie_count": self.near_tie_count,
                "cells_evaluated": self.cells_evaluated,
                "lower_bound": self.lower_bound,
                "upper_bound": self.upper_bound,
                "eta": self.eta,
            },
        )


@dataclass(frozen=True, slots=True)
class HighPrecisionData:
    """Canonical decimal encodings from an arbitrary-precision GECP run."""

    t_nodes: list[str]
    omega_nodes: list[str]
    pivots: list[str]
    residual_history: list[str]
    core_det_history: list[str]
    core_sigma_min_history: list[str]


@dataclass(frozen=True, slots=True)
class ConditioningDiagnostic:
    """Numerical conditioning estimate; unavailable values are null, not infinity."""

    stage: str
    status: str
    columns: int
    effective_rank: int | None = None
    sigma_min: float | None = None
    condition_number: float | None = None
    sigma_min_decimal: str | None = None


@dataclass(frozen=True, slots=True)
class RefinementDiagnostic:
    """Optimizer termination is distinct from representation convergence."""

    status: str = "not_applicable"
    optimizer_success: bool | None = None
    optimizer_status: int | None = None
    evaluations: int = 0
    max_evaluations: int | None = None
    accepted: bool = False
    message: str = ""


@dataclass(frozen=True, slots=True)
class NumericalMetadata:
    """Explicit evidence scope shared by the numerical result records.

    Defaults describe an unrecorded/manual result, never a measured execution.
    Algorithm entry points populate all applicable fields.
    """

    schema_version: int = 1
    algorithm: str = "not_recorded"
    arithmetic: str = "not_recorded"
    precision_bits: int | None = None
    requested_precision_bits: int | None = None
    precision_status: str = "not_recorded"
    conditioning_status: str = "not_recorded"
    conditioning_method: str = "not_recorded"
    conditioning_history: list[ConditioningDiagnostic] = field(default_factory=list)
    tie_status: str = "not_recorded"
    tie_policy: str = "not_recorded"
    tie_rtol: float | None = None
    tie_tolerance_rule: str = "not_recorded"
    tie_counts: list[int] = field(default_factory=list)
    certification_status: str = "not_certified"
    certification_scope: str = "no_analytic_or_rounding_certificate"
    converged: bool | None = None
    stop_reason: str = "not_recorded"
    convergence_scope: str = "not_recorded"
    configuration: dict[str, Any] = field(default_factory=dict)
    refinement: RefinementDiagnostic = field(default_factory=RefinementDiagnostic)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def canonical_json(self) -> str:
        return json.dumps(
            self.as_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
        )


@dataclass(slots=True)
class GECPResult:
    """A complete, auditable record of a finite-grid GECP run."""

    t_nodes: FloatArray
    omega_nodes: FloatArray
    row_indices: list[int]
    column_indices: list[int]
    pivots: list[float]
    residual_history: list[float]
    core_det_history: list[float]
    core_sigma_min_history: list[float]
    tie_counts: list[int]
    pivot_certificates: list[PivotCertificate] = field(default_factory=list)
    rank: int = 0
    converged: bool = False
    stop_reason: str = "not_started"
    config: dict[str, Any] = field(default_factory=dict)
    precision_bits: int = 53
    residual_matrix: FloatArray | None = None
    high_precision: HighPrecisionData | None = None

    @property
    def metadata(self) -> NumericalMetadata:
        """A scoped view over stored GECP evidence, without recomputing estimates."""
        adaptive = (
            self.residual_matrix is None and self.config.get("pivot") == "adaptive"
        )
        default_tie = (
            2.0 ** (5 - self.precision_bits)
            if self.precision_bits > 53
            else 32 * np.finfo(float).eps
        )
        effective_tie = self.config.get("tie_rtol")
        if effective_tie is None:
            effective_tie = default_tie if default_tie else None
        conditioning = [
            ConditioningDiagnostic(
                stage="selected_core",
                status="estimated" if np.isfinite(sigma) else "not_representable",
                columns=i + 1,
                sigma_min=sigma if np.isfinite(sigma) else None,
                sigma_min_decimal=(
                    self.high_precision.core_sigma_min_history[i]
                    if self.high_precision is not None
                    else None
                ),
            )
            for i, sigma in enumerate(self.core_sigma_min_history)
        ]
        return NumericalMetadata(
            algorithm="gecp_interval_adaptive" if adaptive else "gecp_grid",
            arithmetic="mpmath" if adaptive or self.precision_bits > 53 else "float64",
            precision_bits=self.precision_bits,
            requested_precision_bits=int(
                self.config.get("precision_bits", self.precision_bits)
            ),
            precision_status="as_requested",
            conditioning_status=(
                "not_representable"
                if any(record.status == "not_representable" for record in conditioning)
                else "estimated"
                if conditioning
                else "not_applicable"
            ),
            conditioning_method="selected_core_smallest_singular_value_only",
            conditioning_history=conditioning,
            tie_status=(
                "recorded_sampled_candidates"
                if self.pivot_certificates
                else "not_applicable"
            )
            if adaptive
            else "recorded",
            tie_policy="lexicographic_best_singleton_interval_lower_bound"
            if adaptive
            else "lexicographic_near_maximum",
            tie_rtol=(
                self.pivot_certificates[0].tie_rtol if self.pivot_certificates else None
            )
            if adaptive
            else effective_tie,
            tie_tolerance_rule=(
                "sampled_lower_bound_gap<=tie_rtol*max(1,best_lower_bound); "
                "continuous_tie_multiplicity_not_estimated"
                if adaptive
                else "nonzero_unselected_coordinates; gap<=tie_rtol*max(1,max); "
                "default 32*float64_eps or 2^(5-mp_bits)"
            ),
            tie_counts=list(self.tie_counts),
            certification_status=(
                "certified"
                if self.converged and self.stop_reason == "certified_tolerance"
                else "not_certified"
            ),
            certification_scope=(
                "outward_interval_frozen_rounded_residual_model_only"
                if adaptive
                else "sampled_grid_only_no_continuous_certificate"
            ),
            convergence_scope="continuous_interval_bound"
            if adaptive
            else "sampled_grid",
            converged=self.converged,
            stop_reason=self.stop_reason,
            configuration={
                **self.config,
                **(
                    {
                        "interval_objective_scope": (
                            "kernel_minus_approximation_defined_by_frozen_rounded_"
                            "pivot_and_row_column_coefficients"
                        ),
                        "exact_selected_cross_certified": False,
                        "float64_evaluate_residual_certified": False,
                        "tie_scope": "distinct_sampled_coordinates_only",
                    }
                    if adaptive
                    else {}
                ),
            },
        )


@dataclass(slots=True)
class CholeskyResult:
    """Canonical diagonal-pivoted Cholesky output."""

    indices: list[int]
    pivots: list[float]
    residual_history: list[float]
    rank: int
    converged: bool
    stop_reason: str
    factor: FloatArray
    residual_matrix: FloatArray
    metadata: NumericalMetadata = field(default_factory=NumericalMetadata)


@dataclass(slots=True)
class SparseResult:
    """Problem-specific sparse spectral representation."""

    frequencies: FloatArray
    weights: FloatArray
    validation_error: float
    converged: bool
    stop_reason: str
    iterations: int
    metadata: NumericalMetadata = field(default_factory=NumericalMetadata)

    @property
    def training_error(self) -> float:
        """The legacy validation_error is computed on the supplied training samples."""
        return self.validation_error


@dataclass(frozen=True, slots=True)
class PivotGeometry:
    """Physical distances and nearest selected reflected node, not a symmetry proof."""

    t_boundary_distances: tuple[str, str]
    omega_boundary_distances: tuple[str, str]
    dyadic_scale: int
    dyadic_panel: str
    reflection_t: str
    reflection_omega: str
    reflection_nearest_index: int
    reflection_distance: str
    reflection_matched: bool


@dataclass(frozen=True, slots=True)
class CensusRecord:
    """Schema-stable canonical record for one census run."""

    schema_version: int
    git_commit: str
    config_hash: str
    cutoff: str
    tolerance: str
    precision_bits: int
    algorithm: str
    rank: int
    residual: str
    pivots: list[str]
    t_nodes: list[str]
    omega_nodes: list[str]
    core_sigma_min: list[str]
    tie_counts: list[int]
    converged: bool
    stop_reason: str
    residual_history: list[str] = field(default_factory=list)
    core_det_history: list[str] = field(default_factory=list)
    pivot_magnitudes: list[str] = field(default_factory=list)
    pivot_geometry: list[PivotGeometry] = field(default_factory=list)
    grid_shape: tuple[int, int] = (0, 0)
    grid_coordinate_precision_bits: int | None = 53
    grid_generation: str = "float64_composite_chebyshev_lobatto_dyadic_endpoint_panels"
    reflection_match_tolerance: str = "1e-12"
    reflection_distance_metric: str = "max(abs(dt),abs(domega)/cutoff)"
    dyadic_boundary_convention: str = (
        "central_abs_omega_le_1; scale_s:2^(s-1)<abs_omega<=2^s"
    )
    git_dirty: bool = False
    environment: dict[str, str] = field(default_factory=dict)
    metadata: NumericalMetadata = field(default_factory=NumericalMetadata)
