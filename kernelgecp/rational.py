"""Shared-observation Matsubara data and genuine AAA spectral proposals.

The transform is a numerical preprocessing operation, not extra observations.
For a finite spectral measure and fermionic ``z = i*pi*(2*n+1)``,
``-integral(exp(z*t)*G(t), t=0..1) = sum(g/(z-omega))``.  Consequently
real poles of the transformed rational approximation are candidate spectral
frequencies. Poles of a rational fit directly to time data have no such meaning.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.interpolate import AAA, CubicSpline
from scipy.special import roots_legendre

from .diagnostics import array_fingerprint, conditioning_summary
from .spectral_common import AtomProposal
from .types import ConditioningDiagnostic, NumericalMetadata

ComplexArray = npt.NDArray[np.complex128]


def _positive_integer(value: int, name: str, *, allow_zero: bool = False) -> int:
    minimum = 0 if allow_zero else 1
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or value < minimum
    ):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def matsubara_transform(
    t: npt.ArrayLike, values: npt.ArrayLike, *, modes: int = 16
) -> tuple[ComplexArray, ComplexArray, dict[str, Any]]:
    """Transform uniform real observations on [0,1], without density access.

    A not-a-knot cubic spline through the supplied observations is integrated
    with composite Gauss-Legendre quadrature (orders four and eight) over each
    source interval, subdivided to have at least eight intervals per positive
    mode. Spline knots are integration boundaries. Returned conjugate modes
    number ``2*modes``. The discrepancies below are empirical diagnostics,
    not bounds on interpolation, aliasing, measurement noise, or continuous
    transform error. Increasing the integration order adds no observations.
    """
    mode_count = _positive_integer(modes, "modes")
    if np.iscomplexobj(t) or np.iscomplexobj(values):
        raise ValueError("time coordinates and observations must be real")
    times = np.asarray(t, dtype=np.float64)
    targets = np.asarray(values, dtype=np.float64)
    if times.ndim != 1 or times.size < 4 or targets.shape != times.shape:
        raise ValueError("at least four paired one-dimensional observations required")
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(targets)):
        raise ValueError("coordinates and observations must be finite")
    expected = np.linspace(0.0, 1.0, times.size)
    coordinate_tolerance = 16 * np.finfo(float).eps
    if not np.allclose(times, expected, rtol=0.0, atol=coordinate_tolerance):
        raise ValueError("observations must use a uniform grid including 0 and 1")
    # Canonical times remove insignificant coordinate-rounding differences only.
    spline = CubicSpline(expected, targets, bc_type="not-a-knot")
    indices = np.arange(-mode_count, mode_count, dtype=np.float64)
    points = np.asarray(1j * np.pi * (2 * indices + 1), dtype=np.complex128)
    subdivisions = max(1, int(np.ceil(8 * mode_count / (times.size - 1))))
    intervals = (times.size - 1) * subdivisions

    def integrate(interpolant: Any, order: int) -> ComplexArray:
        abscissae, weights = roots_legendre(order)
        midpoints = (np.arange(intervals) + 0.5) / intervals
        grid = (midpoints[:, None] + abscissae[None, :] / (2 * intervals)).ravel()
        combined_weights = np.tile(weights / (2 * intervals), intervals)
        integrand = np.exp(points[:, None] * grid) * interpolant(grid)[None, :]
        return np.asarray(-(integrand @ combined_weights), dtype=np.complex128)

    coarse = integrate(spline, 4)
    transformed = integrate(spline, 8)
    difference = float(np.max(np.abs(transformed - coarse)))
    # Every-other-source comparison measures a different interpolant, unlike
    # oversampling the same spline. Include the final endpoint for even N.
    coarser_indices = np.unique(np.append(np.arange(0, times.size, 2), times.size - 1))
    source_difference: float | None = None
    if coarser_indices.size >= 4:
        coarser = CubicSpline(
            expected[coarser_indices], targets[coarser_indices], bc_type="not-a-knot"
        )
        source_difference = float(np.max(np.abs(transformed - integrate(coarser, 8))))
    if not np.all(np.isfinite(transformed)):
        raise ValueError("nonfinite numerical Matsubara transform")
    diagnostics: dict[str, Any] = {
        "algorithm": "shared_observations_cubic_spline_gauss_matsubara",
        "identity": "-integral_0^1_exp(z*t)*G(t)_dt=sum_g/(z-omega)",
        "sample_count": int(times.size),
        "positive_mode_count": mode_count,
        "mode_count": int(points.size),
        "source_time_fingerprint": array_fingerprint(times),
        "source_target_fingerprint": array_fingerprint(targets),
        "interpolation": "not_a_knot_cubic_spline_of_same_observations",
        "quadrature": "composite_gauss_legendre_on_spline_intervals",
        "quadrature_intervals": intervals,
        "coarse_order": 4,
        "fine_order": 8,
        "coarse_quadrature_nodes": 4 * intervals,
        "fine_quadrature_nodes": 8 * intervals,
        "coarse_fine_max_abs_difference": difference,
        "coarser_source_sample_count": int(coarser_indices.size),
        "source_ladder_max_abs_difference": source_difference,
        "conjugate_symmetry_max_abs_error": float(
            np.max(np.abs(transformed - transformed[::-1].conj()))
        ),
        "arithmetic": "complex128",
        "precision_bits": 53,
        "certification_status": "not_certified",
        "error_scope": (
            "empirical_quadrature_and_interpolation_ladders_not_error_bounds"
        ),
        "limitation": "interpolation_and_observation_noise_are_not_removed",
    }
    return points, transformed, diagnostics


def _complex_conditioning(matrix: ComplexArray) -> ConditioningDiagnostic:
    columns = int(matrix.shape[1])
    if min(matrix.shape) == 0:
        return ConditioningDiagnostic(
            "aaa_loewner_weight_matrix", "not_applicable", columns
        )
    try:
        singular = np.linalg.svd(matrix, compute_uv=False)
    except np.linalg.LinAlgError:
        return ConditioningDiagnostic("aaa_loewner_weight_matrix", "failed", columns)
    if not np.all(np.isfinite(singular)):
        return ConditioningDiagnostic("aaa_loewner_weight_matrix", "failed", columns)
    threshold = max(matrix.shape) * np.finfo(float).eps * float(singular[0])
    rank = int(np.count_nonzero(singular > threshold))
    minimum = float(singular[-1]) if matrix.shape[0] >= columns else 0.0
    condition = float(singular[0] / singular[-1]) if minimum > 0 else None
    return ConditioningDiagnostic(
        "aaa_loewner_weight_matrix",
        "rank_deficient" if rank < columns else "estimated",
        columns,
        rank,
        minimum,
        condition if condition is not None and np.isfinite(condition) else None,
    )


def aaa_atoms(
    z: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    cutoff: float,
    max_atoms: int,
    tolerance: float = 1e-10,
    pole_imag_tolerance: float = 1e-5,
) -> AtomProposal:
    """Return finite, nearly-real, in-cutoff poles from SciPy's AAA algorithm.

    ``tolerance`` is relative to the native data infinity norm. ``max_atoms``
    permits at most ``max_atoms+1`` barycentric support terms, hence at most
    ``max_atoms`` poles. An accepted pole satisfies
    ``abs(imag) <= pole_imag_tolerance*max(1,abs(real))``; its tiny imaginary
    component is discarded explicitly in diagnostics. Out-of-band and complex
    poles are rejected, not clipped or repaired. Native AAA convergence does
    not imply a converged original-time-domain spectral representation.
    """
    budget = _positive_integer(max_atoms, "max_atoms", allow_zero=True)
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if not np.isfinite(pole_imag_tolerance) or pole_imag_tolerance < 0:
        raise ValueError("pole_imag_tolerance must be finite and nonnegative")
    points = np.asarray(z, dtype=np.complex128)
    targets = np.asarray(values, dtype=np.complex128)
    if points.ndim != 1 or points.size < 2 or targets.shape != points.shape:
        raise ValueError("at least two paired one-dimensional samples required")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(targets)):
        raise ValueError("rational samples must be finite")
    if np.unique(points).size != points.size:
        raise ValueError("rational sample points must be distinct")
    scale = float(np.max(np.abs(targets)))
    if not np.isfinite(scale):
        raise ValueError("target magnitude exceeds float64 range")
    max_terms = min(budget + 1, int(points.size))
    conditioning: list[ConditioningDiagnostic] = []
    caught: list[str] = []
    rejected: list[dict[str, Any]] = []
    accepted: list[dict[str, Any]] = []
    poles = np.empty(0, dtype=np.complex128)
    support_count = 0
    error_history: list[float] = []
    native_error = scale
    failure: str | None = None
    if scale > 0 and budget > 0:
        try:
            with warnings.catch_warnings(record=True) as records:
                warnings.simplefilter("always", RuntimeWarning)
                approximation = AAA(
                    points,
                    targets,
                    rtol=tolerance,
                    max_terms=max_terms,
                    clean_up=True,
                )
                poles = np.asarray(approximation.poles(), dtype=np.complex128)
                predictions = np.asarray(approximation(points), dtype=np.complex128)
            caught = [str(record.message) for record in records]
            support = np.asarray(approximation.support_points, dtype=np.complex128)
            support_values = np.asarray(
                approximation.support_values, dtype=np.complex128
            )
            support_count = int(support.size)
            error_history = [float(error) for error in approximation.errors]
            native_error = float(np.max(np.abs(predictions - targets)))
            if not np.isfinite(native_error):
                failure = "nonfinite_native_prediction"
            mask = ~np.isin(points, support)
            if np.any(mask) and support_count:
                loewner = (targets[mask, None] - support_values[None, :]) / (
                    points[mask, None] - support[None, :]
                )
                conditioning.append(_complex_conditioning(loewner))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            failure = f"{type(exc).__name__}: {exc}"
    if scale == 0:
        stop_reason = "zero_data"
    elif budget == 0:
        stop_reason = "max_atoms"
    elif failure is not None:
        stop_reason = "numerical_failure"
    elif native_error <= tolerance * scale:
        stop_reason = "native_tolerance"
    else:
        stop_reason = "max_terms_or_native_stagnation"
    for pole in sorted(poles, key=lambda value: (float(value.real), float(value.imag))):
        raw = [
            float(part) if np.isfinite(part) else None
            for part in (pole.real, pole.imag)
        ]
        if not np.isfinite(pole):
            rejected.append({"pole": raw, "reason": "nonfinite"})
        elif abs(pole.imag) > pole_imag_tolerance * max(1.0, abs(pole.real)):
            rejected.append({"pole": raw, "reason": "nonreal"})
        elif abs(pole.real) > cutoff:
            rejected.append({"pole": raw, "reason": "outside_cutoff"})
        else:
            accepted.append(
                {
                    "pole": raw,
                    "frequency": float(pole.real),
                    "imaginary_component_discarded": float(pole.imag),
                }
            )
    frequencies = np.asarray([item["frequency"] for item in accepted], dtype=np.float64)
    converged = failure is None and native_error <= tolerance * scale
    native_relative_error = native_error / scale if scale else 0.0
    diagnostics = {
        "native_max_abs_error": native_error if np.isfinite(native_error) else None,
        "native_relative_error": (
            native_relative_error if np.isfinite(native_relative_error) else None
        ),
        "native_error_history": [
            error if np.isfinite(error) else None for error in error_history
        ],
        "native_error_history_scope": "before_scipy_froissart_cleanup",
        "sample_count": int(points.size),
        "barycentric_support_count": support_count,
        "barycentric_max_terms": max_terms,
        "raw_pole_count": int(poles.size),
        "accepted_pole_count": int(frequencies.size),
        "rejected_pole_count": len(rejected),
        "accepted_poles": accepted,
        "rejected_poles": rejected,
        "pole_acceptance_rule": "abs(imag)<=tol*max(1,abs(real));abs(real)<=cutoff",
        "scipy_warnings": caught,
        "failure": failure,
        "froissart_cleanup": True,
        "nonlinear_refinement": False,
        "native_domain": "supplied_complex_frequency_samples",
        "spectral_interpretation": (
            "real_poles_only;weights_must_be_fit_in_original_kernel"
        ),
    }
    return AtomProposal(
        frequencies=frequencies,
        metadata=NumericalMetadata(
            algorithm="scipy_aaa_matsubara_spectral_proposal",
            arithmetic="complex128",
            precision_bits=53,
            precision_status="float64_only",
            conditioning_status=conditioning_summary(conditioning),
            conditioning_method="complex_loewner_weight_matrix_svd_null_direction_expected",
            conditioning_history=conditioning,
            tie_status="not_estimated",
            tie_policy="scipy_first_maximum_after_unique_sorted_sample_points",
            certification_status="not_certified",
            certification_scope="no_continuous_transform_or_spectral_support_certificate",
            converged=converged,
            stop_reason=stop_reason,
            convergence_scope="native_complex_sample_relative_infinity_error_only",
            configuration={
                "cutoff": cutoff,
                "max_atoms": budget,
                "tolerance": tolerance,
                "pole_imag_tolerance": pole_imag_tolerance,
                "native_scale": scale,
                "max_terms": max_terms,
            },
        ),
        diagnostics=diagnostics,
    )
