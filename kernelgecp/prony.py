"""Uniform-time, data-ranked Prony-family spectral proposals.

The rectangular Hankel ESPRIT construction is adapted from Potts and Tasche,
``Parameter estimation for nonincreasing exponential sums by Prony-like
methods`` (Linear Algebra and its Applications, 2013). Here exponential nodes
larger than one are allowed: negative physical frequencies are not unstable
poles to be reflected or discarded. This is a numerical proposal, not an
identifiability, minimal-support, or continuous-error certificate.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from .diagnostics import (
    array_fingerprint,
    conditioning_diagnostic,
    conditioning_summary,
)
from .spectral_common import AtomProposal
from .types import (
    ConditioningDiagnostic,
    FloatArray,
    NumericalBreakdownError,
    NumericalMetadata,
    RefinementDiagnostic,
)


def prony_atoms(
    t: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    cutoff: float,
    max_atoms: int,
    tolerance: float = 1e-10,
    rank_rtol: float | None = None,
    pole_imag_tolerance: float = 1e-5,
) -> AtomProposal:
    """Propose physical frequencies using a strided Hankel SVD and ESPRIT.

    ``max_atoms`` is an upper budget, never the supplied true order. Numerical
    rank is selected from the data by an explicit singular-value threshold.
    A strided column construction keeps storage linear in sample count times
    budget, while spanning approximately half the time interval even on dense
    grids. Overlapping rows collectively use every input sample. The eigenvalue
    relation is ``node = exp(-h * omega)`` for the recorded *Hankel* step ``h``.

    Only finite positive-root frequencies with imaginary part at most the
    absolute ``pole_imag_tolerance`` and real part in ``[-cutoff, cutoff]`` are
    retained. Raw complex roots, small-imaginary-part projection, and every
    rejection remain in diagnostics. No out-of-band clipping, pole reflection,
    nonlinear refinement, inverse, or unstable amplitude-to-weight conversion
    is performed. Callers must fit original-kernel weights separately.

    Native convergence refers to a stable exponential-basis least-squares fit
    on these same time samples, with absolute infinity tolerance ``tolerance``.
    It does not assert recovery of the generating spectrum or held-out success.
    """
    if np.iscomplexobj(t) or np.iscomplexobj(values):
        raise ValueError("t and values must be real")
    times = np.asarray(t, dtype=np.float64)
    targets = np.asarray(values, dtype=np.float64)
    if times.ndim != 1 or targets.shape != times.shape or times.size < 2:
        raise ValueError(
            "t and values must be matching 1D arrays with at least 2 samples"
        )
    if (
        not np.all(np.isfinite(times))
        or np.any((times < 0) | (times > 1))
        or not np.all(np.isfinite(targets))
    ):
        raise ValueError("sample coordinates and values must be finite and t in [0, 1]")
    differences = np.diff(times)
    if np.any(differences <= 0):
        raise ValueError("time samples must be strictly increasing")
    sample_step = float((times[-1] - times[0]) / (times.size - 1))
    uniformity_error = float(np.max(np.abs(differences - sample_step)))
    uniformity_tolerance = (
        64 * np.finfo(float).eps * max(float(np.max(np.abs(times))), sample_step)
    )
    if uniformity_error > uniformity_tolerance:
        raise ValueError("Prony requires uniformly spaced time samples")
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if (
        isinstance(max_atoms, (bool, np.bool_))
        or not isinstance(max_atoms, (int, np.integer))
        or max_atoms < 0
    ):
        raise ValueError("max_atoms must be a nonnegative integer")
    max_atoms = int(max_atoms)
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if rank_rtol is not None and (not np.isfinite(rank_rtol) or rank_rtol < 0):
        raise ValueError("rank_rtol must be finite and nonnegative")
    if not np.isfinite(pole_imag_tolerance) or pole_imag_tolerance < 0:
        raise ValueError("pole_imag_tolerance must be finite and nonnegative")

    count = int(times.size)
    effective_budget = min(max_atoms, (count - 1) // 2)
    scale = float(np.max(np.abs(targets)))
    configuration: dict[str, Any] = {
        "cutoff": float(cutoff),
        "max_atoms": max_atoms,
        "effective_atom_budget": effective_budget,
        "tolerance": float(tolerance),
        "rank_rtol": rank_rtol,
        "pole_imag_tolerance": float(pole_imag_tolerance),
        "sample_count": count,
        "sample_step": sample_step,
        "t_fingerprint": array_fingerprint(times),
        "values_fingerprint": array_fingerprint(targets),
    }
    diagnostics: dict[str, Any] = {
        "variant": "rectangular_strided_hankel_esprit",
        "input_domain": "uniform_imaginary_time_samples",
        "sample_count": count,
        "sample_step": sample_step,
        "uniformity_error": uniformity_error,
        "uniformity_tolerance": uniformity_tolerance,
        "effective_atom_budget": effective_budget,
        "sample_limited_budget": effective_budget < max_atoms,
        "normalization_scale": scale,
        "hankel_shape": [0, 0],
        "hankel_stride": None,
        "hankel_step": None,
        "singular_values_normalized": [],
        "rank_threshold_normalized": None,
        "effective_rank_rtol": None,
        "rank_threshold_rule": "not_applied",
        "detected_rank": 0,
        "selected_rank": 0,
        "rank_exceeds_budget": False,
        "raw_roots": [],
        "rejected_root_count": 0,
        "accepted_root_count": 0,
        "small_imaginary_projection_count": 0,
        "hankel_rank_approximation_max_error_normalized": None,
        "subspace_shift_max_error": None,
        "native_fit_error": scale,
        "native_fit_relative_error": 0.0 if scale == 0 else 1.0,
        "native_fit_basis": "endpoint_scaled_real_exponentials",
        "native_fit_weights_are_spectral_weights": False,
        "nonlinear_refinement": "not_requested",
        "frequency_order": "ascending",
    }
    conditioning: list[ConditioningDiagnostic] = []

    def result(frequencies: FloatArray, stop_reason: str) -> AtomProposal:
        error = float(diagnostics["native_fit_error"])
        return AtomProposal(
            frequencies,
            NumericalMetadata(
                algorithm="prony_hankel_esprit_atoms",
                arithmetic="numpy_float64_svd_lstsq_eigvals",
                precision_bits=53,
                requested_precision_bits=53,
                precision_status="fixed_float64",
                conditioning_status=conditioning_summary(conditioning),
                conditioning_method="float64_retained_hankel_shift_and_fit_svd",
                conditioning_history=conditioning,
                tie_status="not_applicable",
                tie_policy="ascending_real_frequency_no_pivot_selection",
                certification_status="not_certified",
                certification_scope="no_spectrum_identifiability_or_continuous_certificate",
                converged=error <= tolerance,
                stop_reason=stop_reason,
                convergence_scope="sampled_time_exponential_fit_absolute_infinity_error",
                configuration=configuration,
                refinement=RefinementDiagnostic(status="not_requested"),
            ),
            diagnostics,
        )

    empty = np.empty(0, dtype=np.float64)
    if scale <= tolerance:
        return result(empty, "zero_model_within_tolerance")
    if effective_budget == 0:
        return result(
            empty, "atom_budget_exhausted" if max_atoms == 0 else "insufficient_samples"
        )

    # Each column is shifted by stride observations. The overlap between the
    # long rows covers the full input sequence without a quadratic-size Hankel.
    stride = max(1, (count - 1) // (2 * effective_budget))
    hankel_step = stride * sample_step
    rows = count - stride * effective_budget
    normalized = targets / scale
    hankel = np.column_stack(
        [
            normalized[stride * j : stride * j + rows]
            for j in range(effective_budget + 1)
        ]
    )
    try:
        left, singular, right_h = np.linalg.svd(hankel, full_matrices=False)
    except np.linalg.LinAlgError as exc:
        raise NumericalBreakdownError("Prony Hankel SVD failed") from exc
    if not np.all(np.isfinite(singular)):
        raise NumericalBreakdownError("Prony Hankel singular values are nonfinite")
    machine_threshold = max(hankel.shape) * np.finfo(float).eps * float(singular[0])
    if rank_rtol is None:
        threshold = max(machine_threshold, tolerance / scale * np.sqrt(hankel.size))
        threshold_rule = (
            "max(matrix_dimension*eps*sigma_max,tolerance/input_scale*sqrt(entries))"
        )
    else:
        threshold = max(machine_threshold, rank_rtol * float(singular[0]))
        threshold_rule = "max(matrix_dimension*eps*sigma_max,rank_rtol*sigma_max)"
    detected_rank = int(np.count_nonzero(singular > threshold))
    rank = min(detected_rank, effective_budget)
    diagnostics.update(
        {
            "hankel_shape": list(hankel.shape),
            "hankel_stride": stride,
            "hankel_step": hankel_step,
            "singular_values_normalized": singular.tolist(),
            "rank_threshold_normalized": float(threshold),
            "effective_rank_rtol": float(threshold / singular[0]),
            "rank_threshold_rule": threshold_rule,
            "detected_rank": detected_rank,
            "selected_rank": rank,
            "rank_exceeds_budget": detected_rank > effective_budget,
        }
    )
    if rank == 0:
        return result(empty, "rank_threshold_excludes_all_modes")
    conditioning.append(
        ConditioningDiagnostic(
            "prony_retained_hankel_subspace",
            "estimated",
            rank,
            rank,
            float(singular[rank - 1]),
            float(singular[0] / singular[rank - 1]),
        )
    )
    basis = right_h[:rank, :].T
    first, shifted = basis[:-1, :], basis[1:, :]
    try:
        shift_matrix = np.linalg.lstsq(first, shifted, rcond=None)[0]
        roots = np.linalg.eigvals(shift_matrix)
    except np.linalg.LinAlgError as exc:
        raise NumericalBreakdownError(
            "Prony subspace shift solve or eigenvalues failed"
        ) from exc
    if not np.all(np.isfinite(shift_matrix)) or not np.all(np.isfinite(roots)):
        raise NumericalBreakdownError("Prony subspace shift or roots are nonfinite")
    conditioning.append(conditioning_diagnostic(first, "prony_subspace_shift_design"))
    diagnostics["subspace_shift_max_error"] = float(
        np.max(np.abs(first @ shift_matrix - shifted))
    )
    diagnostics["hankel_rank_approximation_max_error_normalized"] = float(
        np.max(np.abs(hankel - (left[:, :rank] * singular[:rank]) @ right_h[:rank, :]))
    )

    raw: list[dict[str, Any]] = []
    accepted: list[float] = []
    for root in sorted(roots, key=lambda z: (float(np.real(z)), float(np.imag(z)))):
        node = complex(root)
        record: dict[str, Any] = {
            "node_real": node.real,
            "node_imag": node.imag,
            "frequency_real": None,
            "frequency_imag": None,
            "accepted": False,
            "projected_to_real": False,
            "rejection_reason": None,
        }
        if node == 0:
            record["rejection_reason"] = "zero_node"
        else:
            frequency = -np.log(node) / hankel_step
            real, imag = float(frequency.real), float(frequency.imag)
            if not np.isfinite(real) or not np.isfinite(imag):
                record["rejection_reason"] = "nonfinite_frequency"
            else:
                record["frequency_real"], record["frequency_imag"] = real, imag
                if node.real <= 0:
                    record["rejection_reason"] = "nonpositive_real_node"
                elif abs(imag) > pole_imag_tolerance:
                    record["rejection_reason"] = "nonreal_frequency"
                elif abs(real) > cutoff:
                    record["rejection_reason"] = "frequency_outside_cutoff"
                elif real in accepted:
                    record["rejection_reason"] = "duplicate_frequency"
                else:
                    accepted.append(real)
                    record["accepted"] = True
                    record["projected_to_real"] = imag != 0.0
        raw.append(record)
    frequencies = np.asarray(sorted(accepted), dtype=np.float64)
    diagnostics.update(
        {
            "raw_roots": raw,
            "accepted_root_count": int(frequencies.size),
            "rejected_root_count": len(raw) - int(frequencies.size),
            "small_imaginary_projection_count": sum(
                record["projected_to_real"] for record in raw
            ),
        }
    )
    if frequencies.size:
        # Both signs stay bounded by one by scaling each exponential at its
        # appropriate endpoint. These amplitudes are not physical weights.
        anchors = np.where(frequencies >= 0, times[0], times[-1])
        design = np.exp(-(times[:, None] - anchors[None, :]) * frequencies[None, :])
        try:
            amplitudes = np.linalg.lstsq(design, normalized, rcond=None)[0]
        except np.linalg.LinAlgError as exc:
            raise NumericalBreakdownError(
                "Prony native exponential fit failed"
            ) from exc
        relative_error = float(np.max(np.abs(design @ amplitudes - normalized)))
        native_error = relative_error * scale
        if not np.isfinite(native_error) or not np.all(np.isfinite(amplitudes)):
            raise NumericalBreakdownError("Prony native exponential fit is nonfinite")
        diagnostics["native_fit_error"] = native_error
        diagnostics["native_fit_relative_error"] = relative_error
        conditioning.append(
            conditioning_diagnostic(design, "prony_native_exponential_fit")
        )
    if diagnostics["native_fit_error"] <= tolerance:
        return result(frequencies, "native_time_tolerance_met")
    if diagnostics["rejected_root_count"]:
        return result(frequencies, "physical_roots_rejected")
    if detected_rank > effective_budget:
        return result(frequencies, "atom_budget_exhausted")
    if effective_budget < max_atoms:
        return result(frequencies, "sample_limited_atom_budget")
    return result(frequencies, "rank_selected_native_tolerance_not_met")
