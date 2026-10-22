"""Problem-specific sparse spectral recovery for fixed Green's functions."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.optimize import least_squares

from .diagnostics import (
    array_fingerprint,
    conditioning_diagnostic,
    conditioning_summary,
)
from .kernels import FermionicKernel
from .types import (
    FloatArray,
    NumericalBreakdownError,
    NumericalMetadata,
    RefinementDiagnostic,
    SparseResult,
)


def _design(kernel: FermionicKernel, t: FloatArray, omega: FloatArray) -> FloatArray:
    return np.asarray(kernel(t[:, None], omega[None, :]), dtype=np.float64)


def sparse_representation(
    t: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    cutoff: float,
    tolerance: float = 1e-8,
    max_atoms: int = 16,
    candidate_frequencies: npt.ArrayLike | None = None,
    refine: bool = True,
    tie_rtol: float | None = None,
    refinement_max_nfev: int = 2_000,
) -> SparseResult:
    """OMP plus optional bounded refinement; convergence is training-only.

    Refinement is accepted only when finite and nonworsening in training infinity
    error. Optimizer termination and approximation convergence are separate facts.
    """
    t_values = np.asarray(t, dtype=np.float64)
    targets = np.asarray(values, dtype=np.float64)
    if t_values.ndim != 1 or targets.shape != t_values.shape or t_values.size == 0:
        raise ValueError("t and values must be nonempty one-dimensional arrays")
    if (
        not np.all(np.isfinite(t_values))
        or np.any((t_values < 0) | (t_values > 1))
        or not np.all(np.isfinite(targets))
    ):
        raise ValueError("invalid sample coordinates or values")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if not isinstance(max_atoms, (int, np.integer)) or max_atoms < 1:
        raise ValueError("max_atoms must be a positive integer")
    max_atoms = int(max_atoms)
    if tie_rtol is not None and (not np.isfinite(tie_rtol) or tie_rtol < 0):
        raise ValueError("tie_rtol must be finite and nonnegative")
    if (
        not isinstance(refinement_max_nfev, (int, np.integer))
        or refinement_max_nfev < 1
    ):
        raise ValueError("refinement_max_nfev must be a positive integer")
    refinement_max_nfev = int(refinement_max_nfev)
    kernel = FermionicKernel(cutoff)
    candidates = (
        np.linspace(-cutoff, cutoff, 4_001, dtype=np.float64)
        if candidate_frequencies is None
        else np.asarray(candidate_frequencies, dtype=np.float64)
    )
    if (
        candidates.ndim != 1
        or candidates.size == 0
        or not np.all(np.isfinite(candidates))
        or np.any(np.abs(candidates) > cutoff)
        or np.unique(candidates).size != candidates.size
    ):
        raise ValueError(
            "candidate frequencies must be finite, unique, nonempty and in domain"
        )
    design = _design(kernel, t_values, candidates)
    norms = np.linalg.norm(design, axis=0)
    if np.any(norms == 0) or not np.all(np.isfinite(norms)):
        raise ValueError("candidate design contains an unusable column")

    selected: list[int] = []
    weights = np.empty(0, dtype=np.float64)
    residual = targets.copy()
    scale = max(1.0, float(np.max(np.abs(targets))))
    effective_tie_rtol = tie_rtol if tie_rtol is not None else 32 * np.finfo(float).eps
    conditioning = []
    tie_counts: list[int] = []
    error = float(np.max(np.abs(residual)))
    rank_limit = min(max_atoms, int(candidates.size))
    for _ in range(rank_limit):
        if error <= tolerance * scale:
            break
        with np.errstate(over="raise", invalid="raise"):
            try:
                correlations = np.abs(design.T @ residual) / norms
            except FloatingPointError as exc:
                raise NumericalBreakdownError("nonfinite sparse correlations") from exc
        available = np.ones(candidates.size, dtype=bool)
        available[selected] = False
        maximum = float(np.max(correlations[available]))
        near = np.flatnonzero(
            available
            & (np.abs(correlations - maximum) <= effective_tie_rtol * max(1.0, maximum))
        )
        index = int(near[0])
        tie_counts.append(int(near.size))
        selected.append(index)
        active = design[:, selected]
        try:
            weights = np.asarray(np.linalg.lstsq(active, targets, rcond=None)[0])
        except np.linalg.LinAlgError as exc:
            raise NumericalBreakdownError("sparse least-squares solve failed") from exc
        residual = targets - active @ weights
        if not np.all(np.isfinite(residual)) or not np.all(np.isfinite(weights)):
            raise NumericalBreakdownError("nonfinite sparse least-squares result")
        error = float(np.max(np.abs(residual)))
        conditioning.append(conditioning_diagnostic(active, "omp_design"))

    frequencies = candidates[selected].copy()
    refinement = RefinementDiagnostic(
        status="not_requested" if not refine else "not_applicable"
    )
    if refine and selected:
        count = len(selected)

        def objective(parameters: FloatArray) -> FloatArray:
            return (
                _design(kernel, t_values, parameters[:count]) @ parameters[count:]
                - targets
            )

        initial = np.concatenate((frequencies, weights))
        lower = np.concatenate((np.full(count, -cutoff), np.full(count, -np.inf)))
        upper = np.concatenate((np.full(count, cutoff), np.full(count, np.inf)))
        try:
            refined = least_squares(
                objective,
                initial,
                bounds=(lower, upper),
                xtol=1e-13,
                ftol=1e-13,
                gtol=1e-13,
                max_nfev=refinement_max_nfev,
            )
            parameters = np.asarray(refined.x, dtype=np.float64)
            valid = bool(
                np.all(np.isfinite(parameters))
                and np.all(np.abs(parameters[:count]) <= cutoff)
            )
            new_error = (
                float(np.max(np.abs(objective(parameters)))) if valid else float("inf")
            )
            accepted = bool(np.isfinite(new_error) and new_error <= error)
            if accepted:
                frequencies = parameters[:count].copy()
                weights = parameters[count:].copy()
                error = new_error
            refinement = RefinementDiagnostic(
                status="accepted" if accepted else "rejected_worse_or_nonfinite",
                optimizer_success=bool(refined.success),
                optimizer_status=int(refined.status),
                evaluations=int(refined.nfev),
                max_evaluations=refinement_max_nfev,
                accepted=accepted,
                message=str(refined.message),
            )
            conditioning.append(
                conditioning_diagnostic(
                    np.asarray(refined.jac, dtype=np.float64),
                    "refinement_jacobian_candidate",
                )
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            refinement = RefinementDiagnostic(
                status="failed", max_evaluations=refinement_max_nfev, message=str(exc)
            )

    returned_design = _design(kernel, t_values, frequencies)
    error = float(np.max(np.abs(returned_design @ weights - targets)))
    if not np.isfinite(error):
        raise NumericalBreakdownError("nonfinite returned sparse error")
    conditioning.append(conditioning_diagnostic(returned_design, "returned_design"))
    converged = bool(error <= tolerance * scale)
    stop_reason = (
        "tolerance"
        if converged
        else (
            "candidate_exhausted" if len(selected) == candidates.size else "max_atoms"
        )
    )
    return SparseResult(
        frequencies=frequencies,
        weights=weights,
        validation_error=error,
        converged=converged,
        stop_reason=stop_reason,
        iterations=len(selected),
        metadata=NumericalMetadata(
            algorithm="omp_bounded_nonlinear_refinement" if refine else "omp",
            arithmetic="float64",
            precision_bits=53,
            precision_status="float64_only",
            conditioning_status=conditioning_summary(conditioning),
            conditioning_method="float64_design_and_candidate_jacobian_svd",
            conditioning_history=conditioning,
            tie_status="recorded" if selected else "not_applicable",
            tie_policy="first_available_index_within_correlation_tolerance",
            tie_rtol=effective_tie_rtol,
            tie_tolerance_rule="tie_rtol*max(1,maximum_normalized_correlation)",
            tie_counts=tie_counts,
            certification_scope="training_samples_only_no_continuous_or_support_certificate",
            convergence_scope="training_samples_relative_infinity_error",
            converged=converged,
            stop_reason=stop_reason,
            configuration={
                "cutoff": cutoff,
                "tolerance": tolerance,
                "scale": scale,
                "max_atoms": max_atoms,
                "rank_limit": rank_limit,
                "candidate_count": int(candidates.size),
                "candidate_fingerprint": array_fingerprint(candidates),
                "time_fingerprint": array_fingerprint(t_values),
                "target_fingerprint": array_fingerprint(targets),
                "refine": refine,
                "refinement_max_nfev": refinement_max_nfev,
            },
            refinement=refinement,
        ),
    )
