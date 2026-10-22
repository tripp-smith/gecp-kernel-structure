"""Shared-observation, non-oracle comparisons of spectral representations.

All errors are finite-sample evidence, not analytic continuation, support
identification, global sparsity optimality, or continuous certificates.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import asdict, dataclass, replace
from functools import partial
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from .applications import quadrature_weights, synthetic_spectral_density
from .diagnostics import (
    array_fingerprint,
    conditioning_diagnostic,
    conditioning_summary,
)
from .gecp import gecp
from .kernels import FermionicKernel
from .prony import prony_atoms
from .rational import aaa_atoms, matsubara_transform
from .sparse import sparse_representation
from .spectral_common import AtomProposal
from .types import (
    FloatArray,
    GECPConfig,
    GECPResult,
    NumericalBreakdownError,
    NumericalMetadata,
)
from .vector_fitting import vector_fitting_atoms

METHODS = (
    "universal_gecp",
    "aaa",
    "vector_fitting",
    "prony",
    "nonlinear_least_squares",
    "greedy_pursuit",
)


@dataclass(frozen=True, slots=True)
class ComparisonConfig:
    """Frozen deterministic budgets, never the generating support order."""

    cutoff: float = 8.0
    tolerances: tuple[float, ...] = (1e-6, 1e-8)
    atom_budgets: tuple[int, ...] = (1, 2, 3, 4, 6, 8, 10, 12)
    training_count: int = 513
    validation_count: int = 1024
    candidate_count: int = 1601
    density_count: int = 2001
    matsubara_modes: int = 12
    noise_sigma: float = 1e-8
    seed: int = 20261022
    grid_order: int = 12
    precision_bits: int = 128
    universal_max_rank: int = 64
    native_tolerance: float = 1e-11
    vf_max_iterations: int = 30
    nls_max_evaluations: int = 300
    nls_starts: int = 3

    def __post_init__(self) -> None:
        if (
            isinstance(self.cutoff, bool)
            or not np.isfinite(self.cutoff)
            or self.cutoff < 6
        ):
            raise ValueError("cutoff must be finite and at least six")
        for name, minimum in (
            ("training_count", 17),
            ("validation_count", 8),
            ("candidate_count", 17),
            ("density_count", 101),
            ("matsubara_modes", 2),
            ("seed", 0),
            ("grid_order", 4),
            ("precision_bits", 53),
            ("universal_max_rank", 1),
            ("vf_max_iterations", 1),
            ("nls_max_evaluations", 1),
            ("nls_starts", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer at least {minimum}")
        if (
            not self.tolerances
            or any(
                isinstance(v, bool) or not np.isfinite(v) or not 0 < v < 1
                for v in self.tolerances
            )
            or len(set(self.tolerances)) != len(self.tolerances)
        ):
            raise ValueError("tolerances must be distinct finite values in (0,1)")
        if (
            not self.atom_budgets
            or any(
                isinstance(v, bool) or not isinstance(v, int) or v < 1
                for v in self.atom_budgets
            )
            or tuple(sorted(set(self.atom_budgets))) != self.atom_budgets
        ):
            raise ValueError(
                "atom_budgets must be strictly increasing positive integers"
            )
        if max(self.atom_budgets) * 2 + 1 > self.training_count:
            raise ValueError("training samples must support every tested atom budget")
        if self.matsubara_modes * 4 + 1 > self.training_count:
            raise ValueError("training grid does not resolve requested Matsubara modes")
        if (
            isinstance(self.noise_sigma, bool)
            or not np.isfinite(self.noise_sigma)
            or self.noise_sigma < 0
        ):
            raise ValueError("noise_sigma must be finite and nonnegative")
        if (
            isinstance(self.native_tolerance, bool)
            or not np.isfinite(self.native_tolerance)
            or not 0 < self.native_tolerance < 1
        ):
            raise ValueError("native_tolerance must be finite and in (0,1)")

    def normalized(self) -> dict[str, Any]:
        return asdict(self)


def _design(t: FloatArray, omega: FloatArray, cutoff: float) -> FloatArray:
    return np.asarray(FermionicKernel(cutoff)(t[:, None], omega[None, :]), dtype=float)


def _fit_weights(
    t: FloatArray, values: FloatArray, omega: FloatArray, cutoff: float
) -> FloatArray:
    if not omega.size:
        return np.empty(0, dtype=float)
    return np.asarray(
        np.linalg.lstsq(_design(t, omega, cutoff), values, rcond=None)[0], dtype=float
    )


def _nls_proposal(
    t: FloatArray, values: FloatArray, budget: int, config: ComparisonConfig
) -> AtomProposal:
    """Bounded variable-projection nonlinear LS, independent deterministic starts."""
    history: list[dict[str, Any]] = []
    best = np.empty(0, dtype=float)
    best_error = float("inf")
    for start in range(config.nls_starts):
        # Shifted uniform initializations depend only on budget, not true atoms.
        phase = (start + 0.5) / config.nls_starts
        initial = (
            -config.cutoff + 2 * config.cutoff * (np.arange(budget) + phase) / budget
        )

        def objective(omega: FloatArray) -> FloatArray:
            return (
                _design(t, omega, config.cutoff)
                @ _fit_weights(t, values, omega, config.cutoff)
                - values
            )

        try:
            result = least_squares(
                objective,
                initial,
                bounds=(-config.cutoff, config.cutoff),
                max_nfev=config.nls_max_evaluations,
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
            )
            error = float(np.max(np.abs(objective(result.x))))
            valid = bool(np.isfinite(error) and np.all(np.isfinite(result.x)))
            history.append(
                {
                    "start": start,
                    "success": bool(result.success),
                    "status": int(result.status),
                    "evaluations": int(result.nfev),
                    "training_error": error if valid else None,
                    "message": str(result.message),
                }
            )
            # Selection uses training samples only, never held-out errors.
            if valid and error < best_error:
                best_error = error
                best = np.sort(np.asarray(result.x, dtype=float))
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            history.append({"start": start, "success": False, "message": str(exc)})
    converged = bool(
        best_error <= config.native_tolerance * max(1.0, float(np.max(np.abs(values))))
    )
    return AtomProposal(
        best,
        NumericalMetadata(
            algorithm="bounded_variable_projection_nonlinear_least_squares",
            arithmetic="float64",
            precision_bits=53,
            precision_status="float64_only",
            conditioning_status="reported_by_common_returned_design",
            conditioning_method="common_returned_design_svd",
            tie_status="deterministic",
            tie_policy="first_start_with_strictly_best_training_error",
            certification_scope="sampled_training_only_no_continuous_or_support_certificate",
            convergence_scope="training_infinity_error_not_optimizer_termination",
            converged=converged,
            stop_reason="tolerance"
            if converged
            else "optimizer_budget_or_stationarity",
            configuration={
                "max_nfev_per_start": config.nls_max_evaluations,
                "starts": config.nls_starts,
                "requested_atoms": budget,
                "linear_weights": "least_squares_each_objective_evaluation",
            },
        ),
        {
            "starts": history,
            "native_error": best_error if np.isfinite(best_error) else None,
        },
    )


def assess_proposal(
    proposal: AtomProposal,
    *,
    method: str,
    budget: int,
    training_t: FloatArray,
    observations: FloatArray,
    validation_t: FloatArray,
    validation_values: FloatArray,
    cutoff: float,
    tolerances: tuple[float, ...],
) -> dict[str, Any]:
    """Fit on observations; held-out values affect errors, never fit selection."""
    for coordinates, values in (
        (training_t, observations),
        (validation_t, validation_values),
    ):
        if (
            coordinates.ndim != 1
            or values.shape != coordinates.shape
            or coordinates.size == 0
            or not np.all(np.isfinite(coordinates))
            or not np.all(np.isfinite(values))
            or np.any((coordinates < 0) | (coordinates > 1))
        ):
            raise ValueError("expected nonempty paired finite physical time samples")
    if (
        not tolerances
        or any(not np.isfinite(eps) or eps <= 0 for eps in tolerances)
        or not np.isfinite(cutoff)
        or cutoff <= 0
    ):
        raise ValueError("expected finite positive cutoff and tolerances")
    if np.iscomplexobj(proposal.frequencies):
        raise ValueError("proposal must contain real physical frequencies")
    omega = np.asarray(proposal.frequencies, dtype=float)
    if (
        omega.ndim != 1
        or omega.size > budget
        or not np.all(np.isfinite(omega))
        or np.any(np.abs(omega) > cutoff)
    ):
        raise ValueError("proposal must contain finite physical frequencies")
    weights = _fit_weights(training_t, observations, omega, cutoff)
    training_error = float(
        np.max(np.abs(_design(training_t, omega, cutoff) @ weights - observations))
    )
    held_out_error = float(
        np.max(
            np.abs(_design(validation_t, omega, cutoff) @ weights - validation_values)
        )
    )
    validation_scale = float(np.max(np.abs(validation_values)))
    training_scale = max(1.0, float(np.max(np.abs(observations))))
    converged = bool(training_error <= min(tolerances) * training_scale)
    return {
        "method": method,
        "requested_atom_budget": budget,
        "atom_count": int(omega.size),
        "frequencies": omega.tolist(),
        "weights": weights.tolist(),
        "negative_weight_count": int(np.count_nonzero(weights < 0)),
        "training_error": training_error,
        "held_out_error": held_out_error,
        "converged": converged,
        "stop_reason": "training_tolerance"
        if converged
        else "tested_atom_or_algorithm_budget",
        "convergence_scope": "returned_training_error_at_strictest_requested_tolerance",
        "held_out_relative_error": (
            held_out_error / validation_scale if validation_scale > 0 else None
        ),
        "passes": [
            {"tolerance": epsilon, "passed": bool(held_out_error <= epsilon)}
            for epsilon in tolerances
        ],
        "conditioning": asdict(
            conditioning_diagnostic(
                _design(training_t, omega, cutoff), "returned_time_design"
            )
        ),
        "native_metadata": proposal.metadata.as_dict(),
        "native_diagnostics": proposal.diagnostics,
        "fit_scope": "weights_least_squares_on_shared_observations_only",
        "error_scope": (
            "noiseless_held_out_time_samples_not_continuous_or_support_certificate"
        ),
        "training_input_hash": array_fingerprint(
            np.column_stack((training_t, observations))
        ),
        "precision_bits": 53,
        "certified": False,
    }


def _safe_proposal(call: Callable[[], AtomProposal]) -> AtomProposal:
    """Retain a numerical failure in the cohort, never rename it convergence."""
    try:
        return call()
    except (NumericalBreakdownError, np.linalg.LinAlgError, FloatingPointError) as exc:
        return AtomProposal(
            np.empty(0),
            NumericalMetadata(
                algorithm="failed_proposal_see_method",
                arithmetic="float64",
                precision_bits=53,
                precision_status="float64_only",
                converged=False,
                stop_reason="numerical_failure",
                convergence_scope="native_algorithm_failed",
                conditioning_status="failed",
                conditioning_method="failed",
                tie_status="not_applicable",
                tie_policy="not_applicable",
                certification_scope="no_continuous_or_support_certificate",
            ),
            {"failure": f"{type(exc).__name__}: {exc}"},
        )


def _universal_prefix(
    result: GECPResult, budget: int, tolerance: float
) -> AtomProposal:
    """Metadata describes this prefix, not the successful full generated run."""
    rank = min(budget, result.rank)
    base = result.metadata
    history = base.conditioning_history[:rank]
    converged = result.residual_history[rank] <= tolerance
    return AtomProposal(
        result.omega_nodes[:rank].copy(),
        replace(
            base,
            converged=converged,
            stop_reason="prefix_kernel_tolerance"
            if converged
            else "tested_prefix_budget",
            convergence_scope="selected_prefix_sampled_kernel_grid",
            conditioning_history=history,
            conditioning_status=conditioning_summary(history),
            tie_counts=base.tie_counts[:rank],
            configuration={
                **base.configuration,
                "prefix_rank": rank,
                "full_generated_rank": result.rank,
            },
        ),
        {
            "kernel_rank": result.rank,
            "prefix_rank": rank,
            "prefix_kernel_error": result.residual_history[rank],
            "prefix_kernel_error_decimal": (
                result.high_precision.residual_history[rank]
                if result.high_precision is not None
                else str(result.residual_history[rank])
            ),
            "weights_use_shared_observations": True,
            "node_selection": "kernel_only_complete_pivoting",
            "additional_refinement": False,
        },
    )


def _fixtures(
    config: ComparisonConfig,
) -> list[tuple[str, FloatArray, FloatArray, float]]:
    result = [
        (
            "clean_off_grid",
            np.array([-5.417, -2.231, 0.817, 4.537]),
            np.array([0.18, 0.31, 0.29, 0.22]),
            0.0,
        ),
        (
            "noisy_near_fermi",
            np.array([-3.813, -0.457, 0.367, 3.217]),
            np.array([0.22, 0.28, 0.31, 0.19]),
            config.noise_sigma,
        ),
    ]
    grid = np.linspace(-config.cutoff, config.cutoff, config.density_count)
    for name in ("hubbard_three_peak", "gapped_two_band"):
        density = synthetic_spectral_density(name, grid)
        masses = quadrature_weights(grid) * density
        result.append((name, grid, masses / np.sum(masses), 0.0))
    return result


def run_spectral_comparison(
    config: ComparisonConfig, *, git_commit: str = "not_recorded"
) -> dict[str, Any]:
    """Run all six methods on common observations and preserve every tested fit."""
    t = np.linspace(0.0, 1.0, config.training_count)
    validation_t = (
        np.arange(config.validation_count, dtype=float) + 0.5
    ) / config.validation_count
    if np.intersect1d(t, validation_t).size:
        # Arbitrary allowed counts may share rational coordinates. Keep only
        # truly held-out samples and record their actual count below.
        validation_t = validation_t[~np.isin(validation_t, t)]
    if not validation_t.size:
        raise ValueError("validation coordinates must contain held-out samples")
    candidates = np.linspace(-config.cutoff, config.cutoff, config.candidate_count)
    universal = gecp(
        FermionicKernel(config.cutoff),
        config=GECPConfig(
            tol=min(config.tolerances),
            max_rank=config.universal_max_rank,
            pivot="grid",
            precision_bits=config.precision_bits,
            grid_order=config.grid_order,
        ),
    )
    universal_ranks = []
    for epsilon in config.tolerances:
        indices = [
            k for k, error in enumerate(universal.residual_history) if error <= epsilon
        ]
        rank = indices[0] if indices else universal.rank
        universal_ranks.append(
            {
                "tolerance": epsilon,
                "rank": rank,
                "sampled_kernel_error": universal.residual_history[rank],
                "sampled_kernel_error_decimal": (
                    universal.high_precision.residual_history[rank]
                    if universal.high_precision is not None
                    else str(universal.residual_history[rank])
                ),
                "converged": bool(indices),
            }
        )
    fixtures: list[dict[str, Any]] = []
    for fixture_index, (name, true_omega, true_weights, noise_sigma) in enumerate(
        _fixtures(config)
    ):
        clean = _design(t, true_omega, config.cutoff) @ true_weights
        noise = np.random.default_rng(config.seed + fixture_index).normal(
            0.0, noise_sigma, size=t.size
        )
        observed = clean + noise
        held_out = _design(validation_t, true_omega, config.cutoff) @ true_weights
        z, transformed, transform_diagnostics = matsubara_transform(
            t, observed, modes=config.matsubara_modes
        )
        # Oracle is evaluation-only, never supplied to the rational algorithms.
        exact_transform = np.sum(
            true_weights[None, :] / (z[:, None] - true_omega[None, :]), axis=1
        )
        transform_diagnostics["oracle_validation_error"] = float(
            np.max(np.abs(transformed - exact_transform))
        )
        transform_diagnostics["oracle_scope"] = (
            "evaluation_only_discrete_spectral_truth_not_supplied_to_fit"
        )
        results = []
        for budget in config.atom_budgets:
            common: dict[str, Any] = dict(
                training_t=t,
                observations=observed,
                validation_t=validation_t,
                validation_values=held_out,
                cutoff=config.cutoff,
                tolerances=config.tolerances,
                budget=budget,
            )
            proposals = {
                "universal_gecp": _universal_prefix(
                    universal, budget, min(config.tolerances)
                ),
                "aaa": _safe_proposal(
                    partial(
                        aaa_atoms,
                        z,
                        transformed,
                        cutoff=config.cutoff,
                        max_atoms=budget,
                        tolerance=config.native_tolerance,
                    )
                ),
                "vector_fitting": _safe_proposal(
                    partial(
                        vector_fitting_atoms,
                        z,
                        transformed,
                        cutoff=config.cutoff,
                        max_atoms=budget,
                        tolerance=config.native_tolerance,
                        max_iterations=config.vf_max_iterations,
                    )
                ),
                "prony": _safe_proposal(
                    partial(
                        prony_atoms,
                        t,
                        observed,
                        cutoff=config.cutoff,
                        max_atoms=budget,
                        tolerance=config.native_tolerance,
                    )
                ),
                "nonlinear_least_squares": _nls_proposal(t, observed, budget, config),
            }
            omp = sparse_representation(
                t,
                observed,
                cutoff=config.cutoff,
                max_atoms=budget,
                candidate_frequencies=candidates,
                refine=False,
                tolerance=config.native_tolerance,
            )
            proposals["greedy_pursuit"] = AtomProposal(
                omp.frequencies,
                omp.metadata,
                {
                    "candidate_strategy": "independent_fixed_uniform_dictionary",
                    "candidate_count": candidates.size,
                    "candidate_hash": array_fingerprint(candidates),
                    "native_error": omp.validation_error,
                    "additional_refinement": False,
                },
            )
            for method in METHODS:
                results.append(
                    assess_proposal(proposals[method], method=method, **common)
                )
        smallest = []
        for method in METHODS:
            for epsilon in config.tolerances:
                passing = [
                    r
                    for r in results
                    if r["method"] == method and r["held_out_error"] <= epsilon
                ]
                count = min((r["atom_count"] for r in passing), default=None)
                universal_rank = next(
                    r["rank"] for r in universal_ranks if r["tolerance"] == epsilon
                )
                smallest.append(
                    {
                        "method": method,
                        "tolerance": epsilon,
                        "smallest_tested_passing_atoms": count,
                        "universal_kernel_rank": universal_rank,
                        "support_to_universal_rank": count / universal_rank
                        if count is not None and universal_rank
                        else None,
                    }
                )
        fixtures.append(
            {
                "name": name,
                "noise_sigma": noise_sigma,
                "noise_max": float(np.max(np.abs(noise))),
                "input_cohort": "same_uniform_time_observations",
                "observations": observed.tolist(),
                "training_input_hash": array_fingerprint(
                    np.column_stack((t, observed))
                ),
                "validation_values_hash": array_fingerprint(held_out),
                "transform": transform_diagnostics,
                "results": results,
                "smallest_tested": smallest,
                "true_atom_count": int(true_omega.size)
                if name in ("clean_off_grid", "noisy_near_fermi")
                else None,
                "true_frequencies": true_omega.tolist()
                if name in ("clean_off_grid", "noisy_near_fermi")
                else None,
                "true_weights": true_weights.tolist()
                if name in ("clean_off_grid", "noisy_near_fermi")
                else None,
                "truth_scope": "generating_four_atoms"
                if name in ("clean_off_grid", "noisy_near_fermi")
                else "continuous_density_discretized_by_fixed_trapezoid",
            }
        )
    # Explicitly separate the prior-informed control from support-blind comparisons.
    known_omega = np.array([-5.4, -2.2, 0.8, 4.5])
    known_values = _design(t, known_omega, config.cutoff) @ np.array(
        [0.18, 0.31, 0.29, 0.22]
    )
    known = sparse_representation(
        t,
        known_values,
        cutoff=config.cutoff,
        candidate_frequencies=known_omega,
        refine=False,
        tolerance=config.native_tolerance,
        max_atoms=4,
    )
    config_json = json.dumps(config.normalized(), sort_keys=True, separators=(",", ":"))
    return {
        "schema_version": 1,
        "git_commit": git_commit,
        "config_hash": hashlib.sha256(config_json.encode()).hexdigest(),
        "config": config.normalized(),
        "methods": list(METHODS),
        "training_times": t.tolist(),
        "validation_times": validation_t.tolist(),
        "validation_count_actual": validation_t.size,
        "universal_ranks": universal_ranks,
        "universal_metadata": universal.metadata.as_dict(),
        "universal_high_precision": (
            asdict(universal.high_precision)
            if universal.high_precision is not None
            else None
        ),
        "fixtures": fixtures,
        "known_library_control": {
            "candidate_strategy": "explicit_generating_transition_library_not_blind",
            "method": "greedy_pursuit",
            "atom_count": known.frequencies.size,
            "training_error": known.validation_error,
            "metadata": known.metadata.as_dict(),
        },
        "claim_scope": (
            "finite_sample_comparison_smallest_tested_support_"
            "not_global_optimality_or_unique_recovery"
        ),
    }
