"""Finite-grid Cauchy cross selection and binary64 cardinal evaluation.

All diagnostics describe the rounded finite computation, not an error certificate.
"""

from __future__ import annotations

import math
import operator
from numbers import Real
from typing import Any, cast

import mpmath as mp
import numpy as np
import numpy.typing as npt

from ..fitting import LeastSquaresPlan, NumericArray, _immutable_copy, _numeric_array
from ..types import NumericalBreakdownError

FloatArray = npt.NDArray[np.float64]


def _integer(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        result = operator.index(value)  # type: ignore[arg-type]
    except TypeError as error:
        raise ValueError(f"{name} must be an integer") from error
    if not low <= result <= high:
        raise ValueError(f"{name} must be between {low} and {high}")
    return result


def _positive(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_, str, bytes, complex)):
        raise ValueError(f"{name} must be a positive finite real")
    try:
        number = float(value)  # type: ignore[arg-type]
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{name} must be a positive finite real") from error
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be a positive finite real")
    return number


def _points(value: npt.ArrayLike, name: str, *, unique: bool) -> FloatArray:
    raw = _numeric_array(value, name)
    if raw.ndim != 1 or raw.size == 0 or np.iscomplexobj(raw):
        raise ValueError(f"{name} must be a nonempty real vector")
    x = cast(FloatArray, raw)
    if np.any(x <= 0) or (unique and np.unique(x).size != x.size):
        raise ValueError(f"{name} must be positive" + (" and unique" if unique else ""))
    return x


class CauchyFrequencyPlan:
    """Select a real Cauchy dictionary and fit values at observation frequencies.

    Inputs are positive squared frequencies, poles are nonpositive squared energies.
    Selection uses log Schur products. Evaluation uses stored binary64 weights;
    its condition limit is a pointwise numerical guard, not a certificate.
    """

    __slots__ = (
        "rank",
        "precision_bits",
        "pole_cutoff_ev",
        "max_evaluation_condition",
        "dictionary_size",
        "observation_x",
        "pivot_rows",
        "nodes",
        "poles",
        "log_residual_history",
        "tie_counts",
        "near_tie_counts",
        "weights",
        "_plans",
        "_sealed",
    )

    def __setattr__(self, name: str, value: object) -> None:
        if getattr(self, "_sealed", False):
            raise AttributeError("CauchyFrequencyPlan is immutable")
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        if getattr(self, "_sealed", False):
            raise AttributeError("CauchyFrequencyPlan is immutable")
        object.__delattr__(self, name)

    def __init__(
        self,
        observation_x: npt.ArrayLike,
        *,
        pole_cutoff_ev: float = 13,
        rank: int = 60,
        dictionary_size: int = 1201,
        precision_bits: int = 256,
        max_evaluation_condition: float = 1e4,
    ) -> None:
        x = _points(observation_x, "observation_x", unique=True)
        n = x.size
        size = _integer(dictionary_size, "dictionary_size", 2, 10000)
        self.rank = _integer(rank, "rank", 1, min(n, size, 256))
        self.precision_bits = _integer(precision_bits, "precision_bits", 53, 4096)
        if n > 10000 or n * size > 5_000_000:
            raise ValueError("observation/dictionary resource limit exceeded")
        self.pole_cutoff_ev = _positive(pole_cutoff_ev, "pole_cutoff_ev")
        self.max_evaluation_condition = _positive(
            max_evaluation_condition, "max_evaluation_condition"
        )
        if self.max_evaluation_condition < 1:
            raise ValueError("max_evaluation_condition must be at least 1")
        cutoff2 = self.pole_cutoff_ev * self.pole_cutoff_ev
        if not math.isfinite(cutoff2) or cutoff2 <= 1e-10:
            raise ValueError("squared pole cutoff must be finite and exceed 1e-10")
        self.dictionary_size = size
        poles = -np.r_[0.0, np.geomspace(1e-10, cutoff2, size - 1)]
        if np.unique(poles).size != size:
            raise ValueError("pole dictionary contains duplicate binary64 coordinates")
        with np.errstate(over="ignore", invalid="ignore"):
            distances = x[:, None] - poles
        if not np.all(np.isfinite(distances)):
            raise NumericalBreakdownError("Cauchy distance overflows binary64")
        base = -np.log(distances)
        row_log = np.zeros(n)
        column_log = np.zeros(size)
        rows: list[int] = []
        columns: list[int] = []
        log_residual: list[float] = []
        tie_counts: list[int] = []
        near_tie_counts: list[int] = []
        for _ in range(self.rank):
            score = base + row_log[:, None] + column_log
            flat = int(np.argmax(score))
            i, j = np.unravel_index(flat, score.shape)
            best = float(score[i, j])
            if not math.isfinite(best):
                raise NumericalBreakdownError("no finite Cauchy pivot remains")
            rows.append(int(i))
            columns.append(int(j))
            log_residual.append(best)
            tie_counts.append(int(np.count_nonzero(score == best)))
            near_tie_counts.append(int(np.count_nonzero(score >= best - 1e-12)))
            with np.errstate(divide="ignore", invalid="ignore"):
                row_log += np.log(np.abs(x - x[i])) - np.log(x - poles[j])
                column_log += np.log(np.abs(poles - poles[j])) - np.log(x[i] - poles)
            row_log[i] = -np.inf
            column_log[j] = -np.inf
        self.observation_x = cast(FloatArray, _immutable_copy(x))
        self.pivot_rows = tuple(rows)
        self.nodes = cast(FloatArray, _immutable_copy(x[rows]))
        self.poles = cast(FloatArray, _immutable_copy(poles[columns]))
        self.log_residual_history = tuple(log_residual)
        self.tie_counts = tuple(tie_counts)
        self.near_tie_counts = tuple(near_tie_counts)
        self.weights = self._weights()
        self._plans: dict[float | None, LeastSquaresPlan] = {}
        self._sealed = True

    def _weights(self) -> FloatArray:
        with mp.workprec(self.precision_bits):
            nodes = [mp.mpf(float(v)) for v in self.nodes]
            poles = [mp.mpf(float(v)) for v in self.poles]
            weights = [
                mp.fprod(z - p for p in poles)
                / mp.fprod(z - other for k, other in enumerate(nodes) if k != j)
                for j, z in enumerate(nodes)
            ]
            try:
                rounded = np.asarray([float(w) for w in weights], dtype=np.float64)
            except (OverflowError, ValueError) as error:
                raise NumericalBreakdownError(
                    "cardinal weights unrepresentable"
                ) from error
        if not np.all(np.isfinite(rounded)) or np.any(rounded == 0):
            raise NumericalBreakdownError("cardinal weights unrepresentable")
        return cast(FloatArray, _immutable_copy(rounded))

    def evaluation_diagnostics(
        self, x: npt.ArrayLike
    ) -> tuple[FloatArray, dict[str, Any]]:
        query = _points(x, "x", unique=False)
        result = np.empty((query.size, self.rank), dtype=np.float64)
        conditions: list[float] = []
        max_weight = float(np.max(np.abs(self.weights)))
        normalized = self.weights / max_weight
        if np.any(normalized == 0):
            raise NumericalBreakdownError("scaled cardinal weights underflow")
        for i, z in enumerate(query):
            delta = z - self.nodes
            exact = np.flatnonzero(delta == 0)
            if exact.size:
                result[i] = 0
                result[i, exact[0]] = 1
                conditions.append(1.0)
                continue
            nearest = float(np.min(np.abs(delta)))
            with np.errstate(all="ignore"):
                terms = normalized * (nearest / delta)
                constant = nearest / max_weight
            denominator = math.fsum([constant, *terms])
            numerator = math.fsum([abs(constant), *(abs(t) for t in terms)])
            if not math.isfinite(denominator) or denominator == 0 or numerator == 0:
                raise NumericalBreakdownError(
                    "cardinal denominator is zero or nonfinite"
                )
            condition = numerator / abs(denominator)
            if constant == 0:
                # The omitted positive constant is less than one subnormal ULP.
                # Reject only if its conservative effect could cross the cap.
                tiny = float(np.nextafter(0.0, 1.0))
                lower = float(np.nextafter(abs(denominator) - tiny, 0.0))
                upper = float(np.nextafter(numerator + tiny, np.inf))
                if lower <= 0 or upper / lower > self.max_evaluation_condition:
                    raise NumericalBreakdownError(
                        "underflow prevents conditioning assessment"
                    )
            if (
                not math.isfinite(condition)
                or condition > self.max_evaluation_condition
            ):
                raise NumericalBreakdownError(
                    "cardinal evaluation condition limit exceeded"
                )
            result[i] = terms / denominator
            if not np.all(np.isfinite(result[i])):
                raise NumericalBreakdownError("cardinal evaluation is nonfinite")
            conditions.append(condition)
        return cast(FloatArray, _immutable_copy(result)), {
            "condition_max": max(conditions),
            "conditions": tuple(conditions),
            "precision_bits": 53,
            "certification_status": "not_certified",
        }

    def design(self, x: npt.ArrayLike) -> FloatArray:
        return self.evaluation_diagnostics(x)[0]

    def prepare(self, *, rcond: float | None = None) -> LeastSquaresPlan:
        if rcond is not None:
            if isinstance(rcond, (bool, np.bool_)) or not isinstance(rcond, Real):
                raise ValueError("rcond must be a finite nonnegative real number")
            try:
                rcond = float(rcond)
            except OverflowError as error:
                raise ValueError("rcond must be finite in binary64") from error
            if not math.isfinite(rcond) or rcond < 0:
                raise ValueError("rcond must be a finite nonnegative real number")
        if rcond not in self._plans:
            plan = LeastSquaresPlan(self.design(self.observation_x), rcond=rcond)
            if len(self._plans) >= 8:
                del self._plans[next(iter(self._plans))]
            self._plans[rcond] = plan
        return self._plans[rcond]

    def fit(
        self, values: npt.ArrayLike, *, axis: int = 0, rcond: float | None = None
    ) -> NumericArray:
        return _immutable_copy(self.prepare(rcond=rcond).fit(values, axis=axis))

    def to_payload(self) -> tuple[dict[str, Any], dict[str, npt.NDArray[Any]]]:
        metadata = {
            "version": 1,
            "pole_cutoff_ev": self.pole_cutoff_ev,
            "rank": self.rank,
            "dictionary_size": self.dictionary_size,
            "precision_bits": self.precision_bits,
            "max_evaluation_condition": self.max_evaluation_condition,
            "pivot_rows": list(self.pivot_rows),
            "log_residual_history": list(self.log_residual_history),
            "tie_counts": list(self.tie_counts),
            "near_tie_counts": list(self.near_tie_counts),
            "certification_status": "not_certified",
        }
        arrays: dict[str, npt.NDArray[Any]] = {
            "observation_x": self.observation_x.copy(),
            "nodes": self.nodes.copy(),
            "poles": self.poles.copy(),
            "weights": self.weights.copy(),
        }
        return metadata, arrays

    @classmethod
    def from_payload(
        cls, metadata: dict[str, Any], arrays: dict[str, npt.NDArray[Any]]
    ) -> CauchyFrequencyPlan:
        expected_meta = {
            "version",
            "pole_cutoff_ev",
            "rank",
            "dictionary_size",
            "precision_bits",
            "max_evaluation_condition",
            "pivot_rows",
            "log_residual_history",
            "tie_counts",
            "near_tie_counts",
            "certification_status",
        }
        if (
            set(metadata) != expected_meta
            or set(arrays) != {"observation_x", "nodes", "poles", "weights"}
            or type(metadata["version"]) is not int
            or metadata["version"] != 1
        ):
            raise ValueError("invalid Cauchy payload keys or version")
        if metadata["certification_status"] != "not_certified":
            raise ValueError("invalid certification status")
        for key, value in arrays.items():
            if not isinstance(value, np.ndarray) or value.dtype != np.float64:
                raise ValueError(f"Cauchy payload {key} must be a binary64 array")
        x = _points(arrays["observation_x"], "observation_x", unique=True)
        n = x.size
        size = _integer(metadata["dictionary_size"], "dictionary_size", 2, 10000)
        rank = _integer(metadata["rank"], "rank", 1, min(n, size, 256))
        bits = _integer(metadata["precision_bits"], "precision_bits", 53, 4096)
        if n > 10000 or n * size > 5_000_000:
            raise ValueError("observation/dictionary resource limit exceeded")
        cutoff = _positive(metadata["pole_cutoff_ev"], "pole_cutoff_ev")
        limit = _positive(
            metadata["max_evaluation_condition"], "max_evaluation_condition"
        )
        cutoff2 = cutoff * cutoff
        if not math.isfinite(cutoff2) or cutoff2 <= 1e-10 or limit < 1:
            raise ValueError("invalid cutoff or condition limit")
        for key in ("nodes", "poles", "weights"):
            if arrays[key].shape != (rank,) or not np.all(np.isfinite(arrays[key])):
                raise ValueError(f"invalid Cauchy payload {key}")
        nodes = _points(arrays["nodes"], "nodes", unique=True)
        poles, weights = arrays["poles"], arrays["weights"]
        if (
            np.any(poles > 0)
            or np.any(poles < -cutoff2)
            or np.unique(poles).size != rank
            or np.any(weights == 0)
        ):
            raise ValueError("invalid Cauchy payload poles or weights")
        rows_raw = metadata["pivot_rows"]
        logs_raw = metadata["log_residual_history"]
        ties_raw = metadata["tie_counts"]
        near_raw = metadata["near_tie_counts"]
        for field in (rows_raw, logs_raw, ties_raw, near_raw):
            if not isinstance(field, list) or len(field) != rank:
                raise ValueError("invalid Cauchy payload history")
        rows = tuple(_integer(v, "pivot row", 0, n - 1) for v in rows_raw)
        if len(set(rows)) != rank or not np.array_equal(nodes, x[list(rows)]):
            raise ValueError("Cauchy payload nodes do not match pivot rows")
        logs: list[float] = []
        for value in logs_raw:
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
                raise ValueError("invalid Cauchy log residual")
            converted = float(value)
            if not math.isfinite(converted):
                raise ValueError("invalid Cauchy log residual")
            logs.append(converted)
        ties = tuple(_integer(v, "tie count", 1, n * size) for v in ties_raw)
        near = tuple(_integer(v, "near-tie count", 1, n * size) for v in near_raw)
        if any(a > b for a, b in zip(ties, near, strict=True)):
            raise ValueError("invalid Cauchy tie counts")
        obj = cls.__new__(cls)
        obj.rank = rank
        obj.precision_bits = bits
        obj.pole_cutoff_ev = cutoff
        obj.max_evaluation_condition = limit
        obj.dictionary_size = size
        obj.observation_x = cast(FloatArray, _immutable_copy(x))
        obj.pivot_rows = rows
        obj.nodes = cast(FloatArray, _immutable_copy(nodes))
        obj.poles = cast(FloatArray, _immutable_copy(poles))
        obj.log_residual_history = tuple(logs)
        obj.tie_counts = ties
        obj.near_tie_counts = near
        expected = obj._weights()
        if not np.all(
            np.isclose(weights, expected, rtol=8 * np.finfo(float).eps, atol=0)
        ):
            raise ValueError(
                "Cauchy payload weights inconsistent with stored nodes/poles"
            )
        obj.weights = cast(FloatArray, _immutable_copy(weights))
        obj._plans = {}
        obj._sealed = True
        return obj
