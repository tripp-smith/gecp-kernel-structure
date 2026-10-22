"""Shared float64 conditioning estimates and deterministic input fingerprints."""

from __future__ import annotations

import hashlib

import numpy as np

from .types import ConditioningDiagnostic, FloatArray


def conditioning_summary(history: list[ConditioningDiagnostic]) -> str:
    if not history or all(record.status == "not_applicable" for record in history):
        return "not_applicable"
    for status in ("failed", "rank_deficient", "not_representable"):
        if any(record.status == status for record in history):
            return status
    return "estimated"


def array_fingerprint(array: FloatArray) -> str:
    """Hash shape and canonical little-endian float64 bytes, without storing data."""
    canonical = np.asarray(array, dtype="<f8", order="C")
    return hashlib.sha256(
        str(canonical.shape).encode() + canonical.tobytes()
    ).hexdigest()


def conditioning_diagnostic(matrix: FloatArray, stage: str) -> ConditioningDiagnostic:
    """Estimate two-norm conditioning; this is not an outward-rounded certificate."""
    rows, columns = matrix.shape
    if min(rows, columns) == 0:
        return ConditioningDiagnostic(stage, "not_applicable", columns)
    try:
        singular = np.linalg.svd(matrix, compute_uv=False)
    except np.linalg.LinAlgError:
        return ConditioningDiagnostic(stage, "failed", columns)
    if not np.all(np.isfinite(singular)):
        return ConditioningDiagnostic(stage, "failed", columns)
    threshold = max(rows, columns) * np.finfo(np.float64).eps * float(singular[0])
    rank = int(np.count_nonzero(singular > threshold))
    minimum = float(singular[-1]) if rows >= columns else 0.0
    if rank < columns:
        return ConditioningDiagnostic(stage, "rank_deficient", columns, rank, minimum)
    condition = float(singular[0] / singular[-1])
    return ConditioningDiagnostic(
        stage,
        "estimated",
        columns,
        rank,
        minimum,
        condition if np.isfinite(condition) else None,
    )
