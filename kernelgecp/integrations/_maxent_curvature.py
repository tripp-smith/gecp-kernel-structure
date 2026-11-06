"""Bounded, explicitly shifted projected metric; estimates are uncertified."""

from __future__ import annotations

import math
from numbers import Real
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

from kernelgecp.types import NumericalBreakdownError


def validate_floor_relative(value: object) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(
            "mass_curvature_floor_relative must be finite and in [0, 0.01]"
        )
    numeric = cast(float, value)
    if not math.isfinite(numeric) or not 0 <= numeric <= 1e-2:
        raise ValueError(
            "mass_curvature_floor_relative must be finite and in [0, 0.01]"
        )
    converted = float(value)
    if numeric > 0 and converted == 0:
        raise ValueError("mass_curvature_floor_relative is not representable")
    return converted


def shifted_projected_metric(
    matrix: NDArray[np.float64], floor_relative: float
) -> tuple[NDArray[np.float64], dict[str, Any]]:
    """Return the declared symmetric shifted metric and floating estimates."""
    rho = validate_floor_relative(floor_relative)
    if rho == 0:
        raise ValueError("shifted metric requires a positive relative floor")
    raw = np.asarray(matrix)
    if (
        raw.ndim != 2
        or raw.shape[0] != raw.shape[1]
        or not 1 <= raw.shape[0] <= 128
        or raw.dtype.kind not in "iuf"
    ):
        raise ValueError("curvature must be a real square matrix of dimension 1..128")
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        a = np.asarray(raw, dtype=np.float64)
        half = a * 0.5
        symmetric = np.where(a == a.T, a, half + half.T)
        scale = float(np.max(np.abs(symmetric)))
    unequal = a != a.T
    if np.any(unequal & ((half * 2 != a) | (half.T * 2 != a.T))):
        raise NumericalBreakdownError("curvature half sum is not representable")
    if (
        not np.isfinite(a).all()
        or not math.isfinite(scale)
        or scale < np.finfo(np.float64).tiny
    ):
        raise NumericalBreakdownError("curvature has zero or nonfinite scale")
    try:
        normalized = symmetric / scale
        eigenvalues = np.linalg.eigvalsh(normalized)
        minimum, maximum = float(eigenvalues[0]), float(eigenvalues[-1])
        shift = max(0.0, rho * len(a) - minimum)
        shifted = normalized.copy()
        diagonal = np.diag_indices_from(shifted)
        shifted[diagonal] += shift
        if shift > 0 and np.any(shifted[diagonal] == normalized[diagonal]):
            raise NumericalBreakdownError("curvature shift is not representable")
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            result = symmetric.copy()
            result[diagonal] = shifted[diagonal] * scale
            physical_shift = shift * scale
        if (
            not np.isfinite(result).all()
            or not math.isfinite(physical_shift)
            or (
                shift > 0
                and (
                    physical_shift == 0
                    or np.any(result[diagonal] == symmetric[diagonal])
                )
            )
        ):
            raise NumericalBreakdownError("curvature rescaling is not representable")
        # Check the actual returned matrix, including rounding and underflow.
        returned_eigenvalues = np.linalg.eigvalsh(result / scale)
        actual_minimum = float(returned_eigenvalues[0])
        actual_maximum = float(returned_eigenvalues[-1])
        if not np.isfinite(returned_eigenvalues).all() or actual_minimum <= 0:
            raise NumericalBreakdownError(
                "shifted curvature is not numerically positive"
            )
        condition = actual_maximum / actual_minimum
        if not math.isfinite(condition):
            raise NumericalBreakdownError("curvature condition estimate is nonfinite")
    except np.linalg.LinAlgError as error:
        raise NumericalBreakdownError("curvature eigensolve failed") from error
    return result, {
        "floor_relative": rho,
        "dimension": len(a),
        "scale": scale,
        "normalized_min_eigenvalue": minimum,
        "normalized_max_eigenvalue": maximum,
        "normalized_shift": shift,
        "shift": physical_shift,
        "shifted_normalized_min_eigenvalue": actual_minimum,
        "shifted_normalized_max_eigenvalue": actual_maximum,
        "condition_estimate": condition,
        "certified": False,
    }
