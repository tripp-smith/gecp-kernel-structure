"""Problem-specific sparse spectral recovery for fixed Green's functions."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.optimize import least_squares

from .kernels import FermionicKernel
from .types import FloatArray, SparseResult


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
) -> SparseResult:
    """Recover a sparse signed measure by OMP and bounded frequency refinement."""

    t_values = np.asarray(t, dtype=np.float64)
    targets = np.asarray(values, dtype=np.float64)
    if t_values.ndim != 1 or targets.shape != t_values.shape or t_values.size == 0:
        raise ValueError("t and values must be nonempty one-dimensional arrays")
    if np.any((t_values < 0) | (t_values > 1)) or not np.all(np.isfinite(targets)):
        raise ValueError("invalid sample coordinates or values")
    if tolerance <= 0 or max_atoms < 1:
        raise ValueError("tolerance and max_atoms must be positive")
    kernel = FermionicKernel(cutoff)
    candidates = (
        np.linspace(-cutoff, cutoff, 4_001, dtype=np.float64)
        if candidate_frequencies is None
        else np.asarray(candidate_frequencies, dtype=np.float64)
    )
    if candidates.ndim != 1 or np.any(np.abs(candidates) > cutoff):
        raise ValueError("candidate frequencies must lie in the kernel domain")
    design = _design(kernel, t_values, candidates)
    norms = np.linalg.norm(design, axis=0)
    if np.any(norms == 0):
        raise ValueError("candidate design contains a zero column")

    selected: list[int] = []
    weights = np.empty(0, dtype=np.float64)
    residual = targets.copy()
    stop_reason = "max_atoms"
    converged = False
    scale = max(1.0, float(np.max(np.abs(targets))))

    for _ in range(max_atoms):
        correlations = np.abs(design.T @ residual) / norms
        correlations[selected] = -np.inf
        index = int(np.argmax(correlations))
        selected.append(index)
        active = design[:, selected]
        weights = np.asarray(np.linalg.lstsq(active, targets, rcond=None)[0])
        residual = targets - active @ weights
        if float(np.max(np.abs(residual))) <= tolerance * scale:
            converged = True
            stop_reason = "tolerance"
            break

    frequencies = candidates[selected].copy()
    if refine and selected:
        count = len(selected)

        def objective(parameters: FloatArray) -> FloatArray:
            current_frequencies = parameters[:count]
            current_weights = parameters[count:]
            prediction = _design(kernel, t_values, current_frequencies)
            return prediction @ current_weights - targets

        initial = np.concatenate((frequencies, weights))
        lower = np.concatenate((np.full(count, -cutoff), np.full(count, -np.inf)))
        upper = np.concatenate((np.full(count, cutoff), np.full(count, np.inf)))
        refined = least_squares(
            objective,
            initial,
            bounds=(lower, upper),
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
            max_nfev=2_000,
        )
        frequencies = np.asarray(refined.x[:count], dtype=np.float64)
        weights = np.asarray(refined.x[count:], dtype=np.float64)

    error = float(
        np.max(np.abs(_design(kernel, t_values, frequencies) @ weights - targets))
    )
    if error <= tolerance * scale:
        converged = True
        stop_reason = "tolerance"
    return SparseResult(
        frequencies=frequencies,
        weights=weights,
        validation_error=error,
        converged=converged,
        stop_reason=stop_reason,
        iterations=len(selected),
    )
