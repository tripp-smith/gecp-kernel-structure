"""Realistic synthetic applications of the repository's numerical results."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt

from .approximation import dyadic_taylor_rank_bound
from .cholesky import pivoted_cholesky
from .gecp import cross_approximation, gecp, gecp_matrix
from .kernels import FermionicKernel
from .sparse import sparse_representation
from .types import FloatArray, GECPConfig

SpectrumFixture = Literal["hubbard_three_peak", "gapped_two_band"]


@dataclass(frozen=True, slots=True)
class SyntheticApplicationConfig:
    """Deterministic controls for the post-v1 synthetic application suite."""

    cutoff: float = 8.0
    gecp_tolerance: float = 1e-8
    precision_bits: int = 128
    grid_order: int = 12
    max_rank: int = 64
    quadrature_order: int = 2_001
    validation_time_order: int = 161
    sparse_time_order: int = 121
    sparse_candidate_order: int = 1_601
    sparse_max_atoms: int = 8
    sparse_tolerance: float = 1e-8
    noise_sigma: float = 1e-8
    psd_point_count: int = 72
    psd_max_rank: int = 40
    psd_tolerance: float = 1e-6
    seed: int = 20_260_815

    def __post_init__(self) -> None:
        if not np.isfinite(self.cutoff) or self.cutoff < 6:
            raise ValueError("cutoff must be finite and at least six")
        if not np.isfinite(self.gecp_tolerance) or not 0 < self.gecp_tolerance < 1:
            raise ValueError("gecp_tolerance must lie strictly between zero and one")
        if self.precision_bits < 53:
            raise ValueError("precision_bits must be at least 53")
        if self.grid_order < 4 or self.max_rank < 1:
            raise ValueError("grid_order and max_rank are too small")
        if self.quadrature_order < 101 or self.validation_time_order < 21:
            raise ValueError("validation grids are too small")
        if self.sparse_time_order < 21 or self.sparse_candidate_order < 101:
            raise ValueError("sparse-recovery grids are too small")
        if self.sparse_max_atoms < 4:
            raise ValueError("sparse_max_atoms must permit the four-atom fixtures")
        if not np.isfinite(self.sparse_tolerance) or self.sparse_tolerance <= 0:
            raise ValueError("sparse_tolerance must be finite and positive")
        if not np.isfinite(self.noise_sigma) or self.noise_sigma < 0:
            raise ValueError("noise_sigma must be finite and nonnegative")
        if self.psd_point_count < 16 or self.psd_max_rank < 1:
            raise ValueError("PSD application sizes are too small")
        if not np.isfinite(self.psd_tolerance) or self.psd_tolerance <= 0:
            raise ValueError("psd_tolerance must be finite and positive")

    def normalized(self) -> dict[str, float | int]:
        """Return a stable JSON-compatible configuration mapping."""

        return asdict(self)


@dataclass(frozen=True, slots=True)
class GreenCompressionResult:
    """Universal-basis compression result for one continuous spectrum."""

    fixture: str
    description: str
    cutoff: float
    density_l1: float
    gecp_rank: int
    gecp_converged: bool
    gecp_stop_reason: str
    selected_core_sigma_min: float
    kernel_validation_error: float
    green_validation_error: float
    transfer_bound: float
    transfer_ratio: float
    transfer_verified: bool
    dyadic_rank_bound: int
    frequency_nodes: int
    basis_compression_factor: float


@dataclass(frozen=True, slots=True)
class SparseRecoveryResult:
    """Sparse spectral recovery result with held-out, noiseless validation."""

    fixture: str
    description: str
    candidate_strategy: str
    candidate_count: int
    noise_sigma: float
    true_frequencies: list[float]
    true_weights: list[float]
    recovered_frequencies: list[float]
    recovered_weights: list[float]
    true_atom_count: int
    recovered_atom_count: int
    training_error: float
    held_out_error: float
    max_frequency_error: float | None
    converged: bool
    stop_reason: str
    iterations: int


@dataclass(frozen=True, slots=True)
class PSDLandmarkResult:
    """GECP/Cholesky landmark-selection result for a spatial covariance."""

    point_count: int
    selected_rank: int
    pivot_paths_agree: bool
    gecp_converged: bool
    cholesky_converged: bool
    stop_reason: str
    cross_max_error: float
    cholesky_max_error: float
    relative_frobenius_error: float
    residual_min_eigenvalue: float
    selected_core_sigma_min: float
    compression_factor: float
    selected_indices: list[int]
    selected_points: list[list[float]]


@dataclass(frozen=True, slots=True)
class SyntheticApplicationSuite:
    """Schema-stable collection of all realistic synthetic demonstrations."""

    schema_version: int
    git_commit: str
    config_hash: str
    config: dict[str, float | int]
    green_compression: list[GreenCompressionResult]
    sparse_recovery: list[SparseRecoveryResult]
    psd_landmarks: PSDLandmarkResult

    def as_dict(self) -> dict[str, object]:
        """Return a recursively JSON-compatible dictionary."""

        return asdict(self)

    def canonical_json(self) -> str:
        """Serialize without timestamps or platform-dependent whitespace."""

        return json.dumps(
            self.as_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
        )


def quadrature_weights(nodes: npt.ArrayLike) -> FloatArray:
    """Return composite-trapezoid weights for an increasing nonuniform grid."""

    values = np.asarray(nodes, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ValueError("nodes must be a finite one-dimensional grid")
    differences = np.diff(values)
    if np.any(differences <= 0):
        raise ValueError("nodes must be strictly increasing")
    weights = np.empty_like(values)
    weights[0] = differences[0] / 2.0
    weights[-1] = differences[-1] / 2.0
    weights[1:-1] = (differences[:-1] + differences[1:]) / 2.0
    return weights


def synthetic_spectral_density(
    fixture: SpectrumFixture, omega: npt.ArrayLike
) -> FloatArray:
    """Evaluate a normalized-shape, nonnegative synthetic spectral density.

    Normalization over a finite cutoff is performed by the application runner,
    so this function returns the underlying shape rather than a grid-dependent
    normalized vector.
    """

    frequencies = np.asarray(omega, dtype=np.float64)
    if not np.all(np.isfinite(frequencies)):
        raise ValueError("omega must contain only finite values")

    def gaussian(center: float, width: float) -> FloatArray:
        scaled = (frequencies - center) / width
        return np.asarray(
            np.exp(-0.5 * scaled * scaled) / (width * math.sqrt(2.0 * math.pi)),
            dtype=np.float64,
        )

    if fixture == "hubbard_three_peak":
        density = (
            0.39 * gaussian(-4.0, 0.7)
            + 0.22 * gaussian(0.0, 0.24)
            + 0.39 * gaussian(4.0, 0.7)
        )
    elif fixture == "gapped_two_band":
        half_width = 1.35

        def semicircle(center: float) -> FloatArray:
            coordinate = (frequencies - center) / half_width
            return np.asarray(
                2.0
                * np.sqrt(np.maximum(0.0, 1.0 - coordinate * coordinate))
                / (math.pi * half_width),
                dtype=np.float64,
            )

        density = 0.5 * semicircle(-3.7) + 0.5 * semicircle(3.7)
    else:
        raise ValueError(f"unknown spectral fixture: {fixture!r}")
    return np.asarray(density, dtype=np.float64)


def fermionic_green_from_density(
    t: npt.ArrayLike,
    omega: npt.ArrayLike,
    density: npt.ArrayLike,
    *,
    cutoff: float,
) -> FloatArray:
    """Integrate a sampled spectral density against the fermionic kernel."""

    time_values = np.atleast_1d(np.asarray(t, dtype=np.float64))
    frequencies = np.asarray(omega, dtype=np.float64)
    density_values = np.asarray(density, dtype=np.float64)
    if frequencies.ndim != 1 or density_values.shape != frequencies.shape:
        raise ValueError("omega and density must be matching one-dimensional arrays")
    if not np.all(np.isfinite(density_values)):
        raise ValueError("density must contain only finite values")
    weights = quadrature_weights(frequencies) * density_values
    kernel = FermionicKernel(cutoff)
    return np.asarray(
        kernel(time_values[:, None], frequencies[None, :]) @ weights,
        dtype=np.float64,
    )


def synthetic_sensor_points(count: int, seed: int) -> FloatArray:
    """Generate clustered two-dimensional sensor locations reproducibly."""

    if count < 16:
        raise ValueError("count must be at least sixteen")
    generator = np.random.default_rng(seed)
    centers = np.asarray(
        [[0.18, 0.22], [0.76, 0.23], [0.28, 0.78], [0.78, 0.76]],
        dtype=np.float64,
    )
    labels = np.arange(count) % centers.shape[0]
    points = centers[labels] + generator.normal(0.0, 0.075, size=(count, 2))
    return np.asarray(np.clip(points, 0.0, 1.0), dtype=np.float64)


def _config_hash(config: SyntheticApplicationConfig) -> str:
    canonical = json.dumps(config.normalized(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _normalized_density_weights(
    fixture: SpectrumFixture, omega: FloatArray
) -> tuple[FloatArray, FloatArray]:
    density = synthetic_spectral_density(fixture, omega)
    integration_weights = quadrature_weights(omega)
    mass = float(integration_weights @ density)
    if not np.isfinite(mass) or mass <= 0:
        raise ValueError("synthetic density has invalid mass")
    normalized = density / mass
    return normalized, integration_weights * normalized


def _green_compression_results(
    config: SyntheticApplicationConfig,
) -> list[GreenCompressionResult]:
    kernel = FermionicKernel(config.cutoff)
    result = gecp(
        kernel,
        config=GECPConfig(
            tol=config.gecp_tolerance,
            max_rank=config.max_rank,
            pivot="grid",
            precision_bits=config.precision_bits,
            grid_order=config.grid_order,
        ),
    )
    omega = np.linspace(
        -config.cutoff, config.cutoff, config.quadrature_order, dtype=np.float64
    )
    t = np.linspace(0.0, 1.0, config.validation_time_order, dtype=np.float64)
    exact_kernel = np.asarray(kernel(t[:, None], omega[None, :]), dtype=np.float64)
    selected_cross = cross_approximation(
        kernel, result.t_nodes, result.omega_nodes, t, omega
    )
    kernel_error = float(np.max(np.abs(exact_kernel - selected_cross)))
    sigma_min = (
        result.core_sigma_min_history[-1] if result.core_sigma_min_history else 1.0
    )
    descriptions = {
        "hubbard_three_peak": (
            "narrow quasiparticle peak with two broad Hubbard-like sidebands"
        ),
        "gapped_two_band": "two compact semicircular bands separated by a hard gap",
    }
    outputs: list[GreenCompressionResult] = []
    for fixture in ("hubbard_three_peak", "gapped_two_band"):
        normalized_density, spectral_weights = _normalized_density_weights(
            fixture, omega
        )
        exact_green = exact_kernel @ spectral_weights
        approximate_green = selected_cross @ spectral_weights
        green_error = float(np.max(np.abs(exact_green - approximate_green)))
        density_l1 = float(quadrature_weights(omega) @ np.abs(normalized_density))
        transfer_bound = kernel_error * density_l1
        transfer_ratio = 0.0 if transfer_bound == 0 else green_error / transfer_bound
        outputs.append(
            GreenCompressionResult(
                fixture=fixture,
                description=descriptions[fixture],
                cutoff=config.cutoff,
                density_l1=density_l1,
                gecp_rank=result.rank,
                gecp_converged=result.converged,
                gecp_stop_reason=result.stop_reason,
                selected_core_sigma_min=float(sigma_min),
                kernel_validation_error=kernel_error,
                green_validation_error=green_error,
                transfer_bound=transfer_bound,
                transfer_ratio=transfer_ratio,
                transfer_verified=bool(
                    green_error <= transfer_bound + 64 * np.finfo(float).eps
                ),
                dyadic_rank_bound=dyadic_taylor_rank_bound(
                    config.cutoff, config.gecp_tolerance
                ),
                frequency_nodes=config.quadrature_order,
                basis_compression_factor=config.quadrature_order / max(1, result.rank),
            )
        )
    return outputs


def _sparse_recovery_result(
    *,
    fixture: str,
    description: str,
    true_frequencies: FloatArray,
    true_weights: FloatArray,
    noise_sigma: float,
    restrict_candidates: bool,
    config: SyntheticApplicationConfig,
    generator: np.random.Generator,
) -> SparseRecoveryResult:
    kernel = FermionicKernel(config.cutoff)
    t_train = np.linspace(0.0, 1.0, config.sparse_time_order, dtype=np.float64)
    t_held_out = (
        np.arange(config.sparse_time_order, dtype=np.float64) + 0.5
    ) / config.sparse_time_order
    exact_train = np.asarray(
        kernel(t_train[:, None], true_frequencies[None, :]) @ true_weights,
        dtype=np.float64,
    )
    observed = exact_train + generator.normal(0.0, noise_sigma, size=exact_train.shape)
    candidates = (
        true_frequencies.copy()
        if restrict_candidates
        else np.unique(
            np.concatenate(
                (
                    np.linspace(
                        -config.cutoff,
                        config.cutoff,
                        config.sparse_candidate_order,
                        dtype=np.float64,
                    ),
                    true_frequencies,
                )
            )
        )
    )
    tolerance = max(config.sparse_tolerance, 5.0 * noise_sigma)
    recovered = sparse_representation(
        t_train,
        observed,
        cutoff=config.cutoff,
        tolerance=tolerance,
        max_atoms=config.sparse_max_atoms,
        candidate_frequencies=candidates,
        refine=True,
    )
    exact_held_out = np.asarray(
        kernel(t_held_out[:, None], true_frequencies[None, :]) @ true_weights,
        dtype=np.float64,
    )
    recovered_held_out = np.asarray(
        kernel(t_held_out[:, None], recovered.frequencies[None, :]) @ recovered.weights,
        dtype=np.float64,
    )
    held_out_error = float(np.max(np.abs(exact_held_out - recovered_held_out)))
    frequency_error: float | None = None
    if recovered.frequencies.size == true_frequencies.size:
        frequency_error = float(
            np.max(np.abs(np.sort(recovered.frequencies) - np.sort(true_frequencies)))
        )
    return SparseRecoveryResult(
        fixture=fixture,
        description=description,
        candidate_strategy=(
            "known_transition_library" if restrict_candidates else "dense_blind_scan"
        ),
        candidate_count=int(candidates.size),
        noise_sigma=noise_sigma,
        true_frequencies=true_frequencies.tolist(),
        true_weights=true_weights.tolist(),
        recovered_frequencies=recovered.frequencies.tolist(),
        recovered_weights=recovered.weights.tolist(),
        true_atom_count=int(true_frequencies.size),
        recovered_atom_count=int(recovered.frequencies.size),
        training_error=recovered.validation_error,
        held_out_error=held_out_error,
        max_frequency_error=frequency_error,
        converged=recovered.converged,
        stop_reason=recovered.stop_reason,
        iterations=recovered.iterations,
    )


def _sparse_recovery_results(
    config: SyntheticApplicationConfig,
) -> list[SparseRecoveryResult]:
    generator = np.random.default_rng(config.seed)
    return [
        _sparse_recovery_result(
            fixture="quasiparticle_satellites",
            description=(
                "four clean quasiparticle and satellite lines with unequal weights"
            ),
            true_frequencies=np.asarray([-5.4, -2.2, 0.8, 4.5]),
            true_weights=np.asarray([0.12, 0.23, 0.41, 0.24]),
            noise_sigma=0.0,
            restrict_candidates=True,
            config=config,
            generator=generator,
        ),
        _sparse_recovery_result(
            fixture="noisy_near_fermi_pair",
            description=(
                "two near-Fermi lines plus incoherent side peaks with "
                "deterministic noise"
            ),
            true_frequencies=np.asarray([-3.8, -0.45, 0.35, 3.2]),
            true_weights=np.asarray([0.18, 0.32, 0.28, 0.22]),
            noise_sigma=config.noise_sigma,
            restrict_candidates=False,
            config=config,
            generator=generator,
        ),
    ]


def _psd_landmark_result(config: SyntheticApplicationConfig) -> PSDLandmarkResult:
    points = synthetic_sensor_points(config.psd_point_count, config.seed + 1)
    differences = points[:, None, :] - points[None, :, :]
    squared_distances = np.sum(differences * differences, axis=2)
    length_scale = 0.4
    covariance = np.exp(-0.5 * squared_distances / (length_scale * length_scale))
    covariance += 1e-12 * np.eye(config.psd_point_count)
    algorithm_config = GECPConfig(
        tol=config.psd_tolerance,
        max_rank=config.psd_max_rank,
        pivot="grid",
        precision_bits=53,
    )
    complete = gecp_matrix(covariance, config=algorithm_config)
    cholesky = pivoted_cholesky(covariance, config=algorithm_config)
    indices = complete.row_indices
    if indices:
        core = covariance[np.ix_(indices, indices)]
        selected_columns = covariance[:, indices]
        cross = selected_columns @ np.linalg.solve(core, covariance[indices, :])
        sigma_min = float(np.linalg.svd(core, compute_uv=False)[-1])
    else:
        cross = np.zeros_like(covariance)
        sigma_min = 1.0
    cholesky_approximation = cholesky.factor @ cholesky.factor.T
    cross_error = covariance - cross
    cholesky_error = covariance - cholesky_approximation
    denominator = float(np.linalg.norm(covariance))
    return PSDLandmarkResult(
        point_count=config.psd_point_count,
        selected_rank=complete.rank,
        pivot_paths_agree=(
            complete.row_indices == complete.column_indices == cholesky.indices
        ),
        gecp_converged=complete.converged,
        cholesky_converged=cholesky.converged,
        stop_reason=complete.stop_reason,
        cross_max_error=float(np.max(np.abs(cross_error))),
        cholesky_max_error=float(np.max(np.abs(cholesky_error))),
        relative_frobenius_error=(
            0.0
            if denominator == 0
            else float(np.linalg.norm(cross_error)) / denominator
        ),
        residual_min_eigenvalue=float(
            np.min(np.linalg.eigvalsh(cholesky.residual_matrix))
        ),
        selected_core_sigma_min=sigma_min,
        compression_factor=config.psd_point_count / max(1, complete.rank),
        selected_indices=indices,
        selected_points=points[indices].tolist(),
    )


def run_synthetic_application_suite(
    config: SyntheticApplicationConfig | None = None,
    *,
    git_commit: str = "unknown",
) -> SyntheticApplicationSuite:
    """Run all deterministic realistic synthetic demonstrations."""

    effective_config = SyntheticApplicationConfig() if config is None else config
    return SyntheticApplicationSuite(
        schema_version=1,
        git_commit=git_commit,
        config_hash=_config_hash(effective_config),
        config=effective_config.normalized(),
        green_compression=_green_compression_results(effective_config),
        sparse_recovery=_sparse_recovery_results(effective_config),
        psd_landmarks=_psd_landmark_result(effective_config),
    )
