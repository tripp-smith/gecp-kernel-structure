"""Scalar, real-rational vector fitting for spectral-frequency proposals.

This is a direct implementation of the pole-relocation equations, not a copy
of the SINTEF software.  Each iteration solves
``[Phi, -f*Phi] [r,c] = f`` and relocates poles to ``eig(A-b*c)``.
There are no constant/polynomial terms, electrical stability flips, nonlinear
optimizers, or hidden frequency-domain oracle observations.
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
from .types import ConditioningDiagnostic, FloatArray, NumericalMetadata

ComplexArray = npt.NDArray[np.complex128]


def _encoded(values: ComplexArray) -> list[list[float]]:
    return [[float(value.real), float(value.imag)] for value in values]


def _real_basis(
    z: ComplexArray, poles: ComplexArray
) -> tuple[ComplexArray, FloatArray, FloatArray, ComplexArray]:
    """Real-parameter basis, state realization, and residue conversion map.

    Conjugate poles need conjugate residues, not individually real residues.
    For ``a=u+iv``, the two basis columns are ``phi_a+phi_conj(a)`` and
    ``i*(phi_a-phi_conj(a))``.  Their state block is
    ``A=[[u,v],[-v,u]], b=[2,0]``.  Thus ``c @ (z*I-A)^-1 @ b`` has the
    same denominator residues as the fitted basis, without forming an inverse.
    """
    size = poles.size
    basis = np.empty((z.size, size), dtype=np.complex128)
    state = np.zeros((size, size), dtype=np.float64)
    control = np.zeros(size, dtype=np.float64)
    residue_map = np.zeros((size, size), dtype=np.complex128)
    column = 0
    for index, pole in enumerate(poles):
        if pole.imag < 0:
            continue
        if pole.imag == 0:
            basis[:, column] = 1 / (z - pole)
            state[column, column] = pole.real
            control[column] = 1
            residue_map[index, column] = 1
            column += 1
            continue
        partner = int(np.argmin(np.abs(poles - pole.conjugate())))
        if abs(poles[partner] - pole.conjugate()) > 1e-12 * max(1, abs(pole)):
            raise FloatingPointError("pole set lost conjugate symmetry")
        first = 1 / (z - pole)
        second = 1 / (z - pole.conjugate())
        basis[:, column] = first + second
        basis[:, column + 1] = 1j * (first - second)
        state[column : column + 2, column : column + 2] = [
            [pole.real, pole.imag],
            [-pole.imag, pole.real],
        ]
        control[column] = 2
        residue_map[index, column : column + 2] = [1, 1j]
        residue_map[partner, column : column + 2] = [1, -1j]
        column += 2
    if column != size or not np.all(np.isfinite(basis)):
        raise FloatingPointError("invalid or singular pole basis")
    return basis, state, control, residue_map


def _scaled_solve(
    matrix: ComplexArray,
    targets: ComplexArray,
    stage: str,
    history: list[ConditioningDiagnostic],
) -> FloatArray:
    """Real/imaginary stacked, column-normalized real least-squares solve."""
    stacked = np.vstack((matrix.real, matrix.imag))
    right = np.concatenate((targets.real, targets.imag))
    if not np.all(np.isfinite(stacked)):
        raise FloatingPointError("nonfinite least-squares design")
    # Avoid overflow when computing the Euclidean column norm.
    maxima = np.max(np.abs(stacked), axis=0)
    if np.any(maxima == 0):
        raise FloatingPointError("zero least-squares column")
    norms = maxima * np.linalg.norm(stacked / maxima, axis=0)
    if not np.all(np.isfinite(norms)):
        raise FloatingPointError("unrepresentable least-squares column norm")
    scaled = stacked / norms
    history.append(conditioning_diagnostic(stacked, f"{stage}_unscaled"))
    history.append(conditioning_diagnostic(scaled, f"{stage}_scaled"))
    solution = np.asarray(
        np.linalg.lstsq(scaled, right, rcond=None)[0] / norms,
        dtype=np.float64,
    )
    if not np.all(np.isfinite(solution)):
        raise FloatingPointError("nonfinite least-squares solution")
    return solution


def vector_fitting_atoms(
    z: npt.ArrayLike,
    values: npt.ArrayLike,
    *,
    cutoff: float,
    max_atoms: int,
    tolerance: float = 1e-10,
    max_iterations: int = 30,
    pole_imag_tolerance: float = 1e-5,
) -> AtomProposal:
    """Propose real spectral atoms by genuine iterative vector fitting.

    ``z`` contains distinct off-real sample coordinates (normally fermionic
    Matsubara points); ``values`` are the rational responses at those points.
    All rational coefficients are real in a real/conjugate-pair basis.  Initial
    poles are equally spaced in ``[-cutoff, cutoff]``; a one-pole fit starts at
    zero.  The pole count is capped at the number of complex observations.

    Native convergence is a sampled rational error, scaled by
    ``max(1,max(abs(values)))``.  It does *not* certify a physical spectral fit.
    Final poles with ``abs(imag) <= pole_imag_tolerance*max(1,abs(real))``
    may use their real part; more complex or out-of-band poles are rejected.
    The returned proposal must subsequently have its weights fitted and its
    time-domain error checked independently.  No positivity is imposed.
    """
    points = np.asarray(z, dtype=np.complex128)
    targets = np.asarray(values, dtype=np.complex128)
    if points.ndim != 1 or points.size == 0 or targets.shape != points.shape:
        raise ValueError("z and values must be nonempty one-dimensional arrays")
    if (
        not np.all(np.isfinite(points))
        or not np.all(np.isfinite(targets))
        or np.any(points.imag == 0)
        or np.unique(points).size != points.size
    ):
        raise ValueError("samples must be finite and z distinct and off the real axis")
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    if not np.isfinite(pole_imag_tolerance) or pole_imag_tolerance < 0:
        raise ValueError("pole_imag_tolerance must be finite and nonnegative")
    for label, number in (("max_atoms", max_atoms), ("max_iterations", max_iterations)):
        if (
            isinstance(number, (bool, np.bool_))
            or not isinstance(number, (int, np.integer))
            or number < 0
        ):
            raise ValueError(f"{label} must be a nonnegative integer")
    max_atoms, max_iterations = int(max_atoms), int(max_iterations)
    target_maximum = float(np.max(np.abs(targets)))
    if not np.isfinite(target_maximum):
        raise ValueError("target magnitude exceeds float64 range")
    scale = max(1.0, target_maximum)
    normalized = targets / scale
    rank_limit = min(max_atoms, int(points.size))
    initial = (
        np.zeros(1, dtype=np.complex128)
        if rank_limit == 1
        else np.linspace(-cutoff, cutoff, rank_limit).astype(np.complex128)
    )
    poles = initial.copy()
    selected = np.empty(0, dtype=np.complex128)
    selected_residues = np.empty(0, dtype=np.complex128)
    conditioning: list[ConditioningDiagnostic] = []
    iterations: list[dict[str, Any]] = []
    best_error = float(np.max(np.abs(normalized)))
    stop_reason = "atom_budget" if rank_limit == 0 else "iteration_budget"
    failure_message: str | None = None
    relocations = 0
    if best_error <= tolerance:
        stop_reason = "zero_model_tolerance"
    elif rank_limit:
        try:
            for iteration in range(max_iterations + 1):
                basis, state, control, residue_map = _real_basis(points, poles)
                coefficients = _scaled_solve(
                    basis, normalized, f"iteration_{iteration}_residues", conditioning
                )
                predicted = basis @ coefficients
                error = float(np.max(np.abs(predicted - normalized)))
                if not np.isfinite(error):
                    raise FloatingPointError("nonfinite native error")
                residues = residue_map @ coefficients * scale
                if not np.all(np.isfinite(residues)):
                    raise FloatingPointError("unrepresentable native residues")
                iterations.append(
                    {
                        "iteration": iteration,
                        "poles": _encoded(poles),
                        "native_relative_error": error,
                    }
                )
                if error < best_error:
                    best_error = error
                    selected, selected_residues = poles.copy(), residues.copy()
                if error <= tolerance:
                    stop_reason = "native_tolerance"
                    break
                if iteration == max_iterations:
                    break
                relocation_design = np.hstack((basis, -normalized[:, None] * basis))
                relocation = _scaled_solve(
                    relocation_design,
                    normalized,
                    f"iteration_{iteration}_relocation",
                    conditioning,
                )
                denominator = relocation[rank_limit:]
                # This real-block realization is similar to diag(a)-ones*c.T,
                # but preserves conjugate symmetry without projecting roots.
                relocated = np.sort_complex(
                    np.asarray(
                        np.linalg.eigvals(state - np.outer(control, denominator)),
                        dtype=np.complex128,
                    )
                )
                if not np.all(np.isfinite(relocated)):
                    raise FloatingPointError("nonfinite relocated poles")
                movement = float(np.max(np.abs(relocated - np.sort_complex(poles))))
                iterations[-1]["relocated_poles"] = _encoded(relocated)
                iterations[-1]["maximum_sorted_pole_movement"] = movement
                poles = relocated
                relocations += 1
        except (np.linalg.LinAlgError, FloatingPointError) as exc:
            stop_reason = "numerical_breakdown"
            failure_message = str(exc)

    frequencies: list[float] = []
    rejected: list[dict[str, Any]] = []
    projected = 0
    duplicate_tolerance = 32 * np.finfo(float).eps * max(1.0, cutoff)
    for pole in np.sort_complex(selected):
        reason: str | None = None
        if abs(pole.imag) > pole_imag_tolerance * max(1.0, abs(pole.real)):
            reason = "nonreal"
        elif abs(pole.real) > cutoff:
            reason = "outside_cutoff"
        elif frequencies and abs(pole.real - frequencies[-1]) <= duplicate_tolerance:
            reason = "duplicate_within_float64_tolerance"
        if reason:
            rejected.append({"pole": _encoded(np.array([pole]))[0], "reason": reason})
        else:
            frequencies.append(float(pole.real))
            projected += int(pole.imag != 0)
    returned = np.asarray(frequencies, dtype=np.float64)
    physical_error: float | None = None
    physical_residues: list[float] = []
    try:
        if returned.size:
            physical_basis = 1 / (points[:, None] - returned[None, :])
            physical_coefficients = _scaled_solve(
                physical_basis, normalized, "physical_residues", conditioning
            )
            physical_error = float(
                np.max(np.abs(physical_basis @ physical_coefficients - normalized))
            )
            physical_residues = (physical_coefficients * scale).tolist()
        else:
            physical_error = float(np.max(np.abs(normalized)))
    except (np.linalg.LinAlgError, FloatingPointError) as exc:
        failure_message = (
            f"{failure_message or ''}; physical residue fit: {exc}"
        ).lstrip("; ")
    native_absolute_error = best_error * scale
    return AtomProposal(
        frequencies=returned,
        metadata=NumericalMetadata(
            algorithm="scalar_real_rational_vector_fitting",
            arithmetic="float64_complex128",
            precision_bits=53,
            precision_status="float64_only",
            conditioning_status=conditioning_summary(conditioning),
            conditioning_method="real_imag_stacked_design_svd_before_and_after_scaling",
            conditioning_history=conditioning,
            tie_status="not_applicable",
            tie_policy="equally_spaced_initial_poles_sorted_eigenvalues_best_native_error",
            converged=best_error <= tolerance,
            stop_reason=stop_reason,
            convergence_scope="native_rational_samples_scaled_infinity_error_only",
            certification_scope="no_continuous_time_or_support_certificate",
            configuration={
                "cutoff": cutoff,
                "max_atoms": max_atoms,
                "effective_atom_budget": rank_limit,
                "tolerance": tolerance,
                "scale": scale,
                "max_iterations": max_iterations,
                "pole_imag_tolerance": pole_imag_tolerance,
                "duplicate_tolerance": duplicate_tolerance,
                "polynomial_terms": "none_strictly_proper_rational_model",
                "stability_flipping": False,
                "coefficient_constraint": "real_or_conjugate_pair_residues",
                "sample_fingerprint": array_fingerprint(
                    np.column_stack((points.real, points.imag))
                ),
                "target_fingerprint": array_fingerprint(
                    np.column_stack((targets.real, targets.imag))
                ),
            },
        ),
        diagnostics={
            "initial_poles": _encoded(initial),
            "raw_poles": _encoded(selected),
            "raw_residues": _encoded(selected_residues),
            "native_relative_error": best_error,
            "native_absolute_error": native_absolute_error
            if np.isfinite(native_absolute_error)
            else None,
            "physical_native_relative_error": physical_error,
            "physical_native_converged": physical_error is not None
            and physical_error <= tolerance,
            "physical_native_residues": physical_residues,
            "relocations": relocations,
            "history": iterations,
            "rejected_poles": rejected,
            "rejected_pole_count": len(rejected),
            "projected_small_imaginary_count": projected,
            "raw_pole_count": int(selected.size),
            "physical_pole_count": int(returned.size),
            "failure_message": failure_message,
            "iteration_selection": "lowest_native_error_including_zero_model",
            "sample_limited": max_atoms > rank_limit,
            "projection_rule": "abs(imag)<=pole_imag_tolerance*max(1,abs(real))",
            "nonclaims": "native convergence is not time-fit or support certification",
        },
    )
