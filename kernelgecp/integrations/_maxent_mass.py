"""Positive exponential-family coordinates with finite-precision fixed mass."""

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any

import numpy as np
from scipy.special import logsumexp

from ..types import FloatArray, NumericalBreakdownError

_MAX_ROWS = 2048
_MAX_COLUMNS = 128
_MAX_HESSIAN_BYTES = 128 * 1024 * 1024


def _immutable(value: Any) -> FloatArray:
    array = np.ascontiguousarray(value, dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)


def _shape(value: Any, dimensions: int) -> tuple[int, ...]:
    """Inspect supported inputs without allocating a coerced float array."""
    if isinstance(value, np.ndarray):
        if value.dtype.kind not in "iuf":
            raise ValueError("arrays must have a real numeric dtype")
        return tuple(value.shape)
    if not isinstance(value, (list, tuple)):
        raise ValueError("expected an ndarray, list, or tuple")
    if dimensions == 1:
        return (len(value),)
    if not value:
        return (0, 0)
    if not isinstance(value[0], (list, tuple, np.ndarray)):
        raise ValueError("features must be a rectangular matrix")
    columns = len(value[0])
    if any(
        not isinstance(row, (list, tuple, np.ndarray)) or len(row) != columns
        for row in value
    ):
        raise ValueError("features must be a rectangular matrix")
    return (len(value), columns)


def _real_array(value: Any) -> FloatArray:
    array = np.asarray(value)
    if array.dtype.kind not in "iuf":
        raise ValueError("arrays must have a real numeric dtype")
    return np.asarray(array, dtype=float)


def _hessian_peak(n: int, rank: int) -> int:
    # Immutable returned bytes duplicate the tensor; include explicit centered,
    # weighted-centered, covariance, h, and p arrays (not library/RSS costs).
    return 16 * n * rank * rank + 16 * n * rank + 8 * rank * rank + 16 * n


class NormalizedMassMap:
    """Map gauge-free coordinates to positive mass contributions.

    Features have one row per contribution. Constant feature directions are
    removed before an SVD; retained left singular vectors form the coordinates.
    The inverse accepts only positive, representable members of this family.
    """

    def __init__(
        self, default_weights: Any, features: Any, mass: float, *, rcond: float = 1e-12
    ) -> None:
        d_shape = _shape(default_weights, 1)
        feature_shape = _shape(features, 2)
        if (
            len(d_shape) != 1
            or len(feature_shape) != 2
            or feature_shape[0] != d_shape[0]
        ):
            raise ValueError("default weights and feature dimensions do not agree")
        n, m = feature_shape
        if not n or not m or n > _MAX_ROWS or m > _MAX_COLUMNS:
            raise ValueError("feature dimensions exceed resource caps or are empty")
        if _hessian_peak(n, min(n, m)) > _MAX_HESSIAN_BYTES:
            raise ValueError("Hessian estimated explicit-array peak cap exceeded")
        d = _real_array(default_weights)
        features = _real_array(features)
        if not np.isfinite(d).all() or np.any(d <= 0):
            raise ValueError("default weights must be finite and positive")
        if not np.isfinite(features).all():
            raise ValueError("features must be finite")
        if isinstance(mass, (bool, np.bool_)) or not isinstance(
            mass, (int, float, np.integer, np.floating)
        ):
            raise ValueError("mass must be a real scalar")
        if not np.isfinite(mass) or mass < np.finfo(float).tiny:
            raise ValueError("mass must be finite and positive")
        if not np.isfinite(rcond) or not 0 < rcond < 1:
            raise ValueError("rcond must be finite and between zero and one")
        centered = features - features[0, :]
        centered -= centered.mean(axis=0)
        if not np.isfinite(centered).all():
            raise NumericalBreakdownError("feature centering is nonfinite")
        try:
            u, singular, _ = np.linalg.svd(centered, full_matrices=False)
            original_singular = np.linalg.svd(features, compute_uv=False)
        except np.linalg.LinAlgError as error:
            raise NumericalBreakdownError("feature SVD failed") from error
        # Relative to centered variation only: constant offsets are gauge.
        rank = int(np.count_nonzero(singular > rcond * singular[0]))
        if not rank:
            raise ValueError("centered features have zero effective rank")
        self._basis = _immutable(u[:, :rank])
        self._default = _immutable(d)
        self._singular_values = _immutable(singular)
        self._mass = float(mass)
        self._log_default = _immutable(np.log(d))
        self._diagnostics = MappingProxyType(
            {
                "original_rows": n,
                "original_columns": m,
                "original_rank": int(
                    np.count_nonzero(original_singular > rcond * original_singular[0])
                ),
                "effective_rank": rank,
                "rcond": float(rcond),
                "singular_values": tuple(float(s) for s in singular),
                "retained_condition": float(singular[0] / singular[rank - 1]),
                "hessian_bytes": n * rank * rank * 8,
                "estimated_explicit_array_peak_bytes": _hessian_peak(n, rank),
                "rank_truncation": "heuristic_relative_centered_svd",
                "certification": "not_certified",
            }
        )
        # Check that even the default distribution is representable.
        self.f(np.zeros(rank))

    @property
    def rank(self) -> int:
        return int(self._basis.shape[1])

    @property
    def input_size(self) -> int:
        return self.rank

    @property
    def basis(self) -> FloatArray:
        return self._basis

    @property
    def default_weights(self) -> FloatArray:
        return self._default

    @property
    def singular_values(self) -> FloatArray:
        return self._singular_values

    @property
    def mass(self) -> float:
        return self._mass

    @property
    def diagnostics(self) -> Mapping[str, Any]:
        return self._diagnostics

    def _vector(self, value: Any, length: int, name: str) -> FloatArray:
        result = _real_array(value)
        if result.shape != (length,) or not np.isfinite(result).all():
            raise ValueError(f"{name} must be a finite vector of length {length}")
        return result

    def f(self, v: Any) -> FloatArray:
        v = self._vector(v, self.rank, "coordinates")
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            logits = self._log_default + self._basis @ v
            logits -= np.max(logits)
            result = self._mass * np.exp(logits - logsumexp(logits))
        if (
            not np.isfinite(result).all()
            or np.any(result <= 0)
            or not np.isfinite(result.sum())
            or abs(result.sum() / self._mass - 1) > 1e-12
        ):
            raise NumericalBreakdownError("normalized mass map underflow or overflow")
        return _immutable(result)

    def d(self, v: Any) -> FloatArray:
        h = self.f(v)
        p = h / self._mass
        centered = self._basis - p @ self._basis
        result = h[:, None] * centered
        if not np.isfinite(result).all():
            raise NumericalBreakdownError("normalized mass Jacobian is nonfinite")
        return _immutable(result)

    def dd(self, v: Any) -> FloatArray:
        h = self.f(v)
        p = h / self._mass
        centered = self._basis - p @ self._basis
        covariance = centered.T @ (p[:, None] * centered)
        result = np.empty((h.size, self.rank, self.rank))
        np.multiply(centered[:, :, None], centered[:, None, :], out=result)
        result -= covariance
        result *= h[:, None, None]
        if not np.isfinite(result).all():
            raise NumericalBreakdownError("normalized mass Hessian is nonfinite")
        return _immutable(result)

    def inv(self, h: Any) -> FloatArray:
        h = self._vector(h, self._default.size, "mass contributions")
        if np.any(h <= 0) or not np.isclose(h.sum(), self._mass, rtol=1e-10, atol=0):
            raise ValueError("mass contributions must be positive and sum to mass")
        log_ratio = np.log(h) - self._log_default
        log_ratio -= log_ratio.mean()
        try:
            v, _, _, _ = np.linalg.lstsq(self._basis, log_ratio, rcond=None)
        except np.linalg.LinAlgError as error:
            raise NumericalBreakdownError(
                "inverse least-squares solve failed"
            ) from error
        if not np.isfinite(v).all() or not np.allclose(
            self.f(v), h, rtol=1e-9, atol=self._mass * 1e-13 / h.size
        ):
            raise ValueError("mass contributions are outside the represented family")
        return _immutable(v)
