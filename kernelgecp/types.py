"""Typed configuration and result records shared by numerical algorithms."""

from __future__ import annotations

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
        if self.certificate_abs_tol < 0 or self.certificate_rel_tol < 0:
            raise ValueError("certificate tolerances must be nonnegative")
        if self.max_cells < 1:
            raise ValueError("max_cells must be positive")
        if self.tie_rtol is not None and self.tie_rtol < 0:
            raise ValueError("tie_rtol must be nonnegative")

    def normalized(self) -> dict[str, Any]:
        """Return a deterministic JSON-compatible configuration mapping."""

        return asdict(self)


DEFAULT_GECP_CONFIG = GECPConfig()


@dataclass(frozen=True, slots=True)
class PivotCertificate:
    """Bounds witnessing the quality of one approximate continuous pivot."""

    t: float
    omega: float
    value: float
    lower_bound: float
    upper_bound: float
    eta: float
    cells_evaluated: int
    certified: bool
    termination_reason: str


@dataclass(frozen=True, slots=True)
class HighPrecisionData:
    """Canonical decimal encodings from an arbitrary-precision GECP run."""

    t_nodes: list[str]
    omega_nodes: list[str]
    pivots: list[str]
    residual_history: list[str]
    core_det_history: list[str]
    core_sigma_min_history: list[str]


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


@dataclass(slots=True)
class SparseResult:
    """Problem-specific sparse spectral representation."""

    frequencies: FloatArray
    weights: FloatArray
    validation_error: float
    converged: bool
    stop_reason: str
    iterations: int


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
