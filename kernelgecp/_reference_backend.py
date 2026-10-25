"""Explicit, bounded bridge to the optional 128-bit reference accelerator.

Importing this module does not import a native extension.  Portable archives
store coordinates and coefficients; native handles are always reconstructed.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable
from typing import Any, Literal, Protocol, cast

import numpy as np
import numpy.typing as npt

from .types import FloatArray, NumericalBreakdownError

type ReferenceBackend = Literal["python", "rust"]
type NativeInfo = tuple[tuple[str, str | int | bool], ...]

_PROTOCOL_VERSION = 1
_WORKING_BITS = 128
_MAX_RANK = 128
_BATCH_ROWS = 256
_MAX_CELLS = 2**20
_MAX_FACTOR_BYTES = 2 * _MAX_RANK**2 * 96


class _NativeHandle(Protocol):
    @property
    def rank(self) -> int: ...

    def design_tau(self, normalized_times: list[float]) -> bytes: ...


def _native_info(value: object) -> NativeInfo:
    if not isinstance(value, dict) or not value or len(value) > 32:
        raise NumericalBreakdownError(
            "Rust reference backend returned invalid build info"
        )
    copied: dict[str, str | int | bool] = {}
    for key, item in value.items():
        if (
            not isinstance(key, str)
            or not key
            or len(key) > 64
            or type(item) not in (str, int, bool)
            or (isinstance(item, str) and len(item) > 256)
            or (type(item) is int and not -(2**63) <= item < 2**63)
        ):
            raise NumericalBreakdownError(
                "Rust reference backend returned invalid build info"
            )
        copied[key] = item
    for key, expected in (
        ("protocol_version", _PROTOCOL_VERSION),
        ("working_precision_bits", _WORKING_BITS),
    ):
        if type(copied.get(key)) is not int or copied[key] != expected:
            raise NumericalBreakdownError(
                f"Rust reference backend build info disagrees on {key}"
            )
    return tuple(sorted(copied.items()))


def _load_native() -> Any:
    try:
        module = importlib.import_module("kernelgecp._native")
    except Exception as error:
        raise NumericalBreakdownError(
            "Rust reference backend is unavailable; install a compatible kernelgecp "
            "platform wheel or use reference_backend='python'"
        ) from error
    for key, expected in (
        ("NATIVE_PROTOCOL_VERSION", _PROTOCOL_VERSION),
        ("WORKING_PRECISION_BITS", _WORKING_BITS),
    ):
        actual = getattr(module, key, None)
        if type(actual) is not int or actual != expected:
            raise NumericalBreakdownError(
                f"Rust reference backend has an incompatible {key}"
            )
    return module


class RustReferenceBackend:
    """Own a checked native handle for one existing mpmath transpose LU.

    Times supplied here are already normalized by the existing NumPy path.
    The native implementation must not repeat physical coordinate division.
    """

    def __init__(
        self,
        omega_nodes: FloatArray,
        ordering: list[int],
        lower: Any,
        upper: Any,
        context: Any,
    ) -> None:
        rank = int(omega_nodes.size)
        if context.prec != _WORKING_BITS:
            raise ValueError("Rust reference backend requires exactly 128 working bits")
        if not 1 <= rank <= _MAX_RANK:
            raise ValueError("Rust reference backend supports ranks 1 through 128")
        if (
            omega_nodes.ndim != 1
            or not np.all(np.isfinite(omega_nodes))
            or np.any(np.abs(omega_nodes) > 16384)
            or np.unique(omega_nodes).size != rank
            or len(ordering) != rank
            or any(type(item) is not int for item in ordering)
            or sorted(ordering) != list(range(rank))
            or lower.rows != rank
            or lower.cols != rank
            or upper.rows != rank
            or upper.cols != rank
        ):
            raise ValueError("invalid Rust reference factor dimensions or coordinates")
        total_bytes = 0

        def pack(matrix: Any) -> list[list[str]]:
            nonlocal total_bytes
            rows: list[list[str]] = []
            for i in range(rank):
                row: list[str] = []
                for j in range(rank):
                    value = matrix[i, j]
                    if not context.isfinite(value):
                        raise ValueError("Rust reference factors must be finite")
                    token = str(context.nstr(value, 60))
                    total_bytes += len(token.encode("utf-8"))
                    if len(token) > 96 or total_bytes > _MAX_FACTOR_BYTES:
                        raise ValueError(
                            "Rust reference factor transfer exceeds its budget"
                        )
                    if context.mpf(token) != value:
                        raise NumericalBreakdownError(
                            "Rust reference factor transfer changed a 128-bit value"
                        )
                    row.append(token)
                rows.append(row)
            return rows

        packed_lower, packed_upper = pack(lower), pack(upper)
        module = _load_native()
        try:
            self.info = _native_info(module.build_info())
            factory = cast(
                Callable[
                    [int, list[float], list[int], list[list[str]], list[list[str]]],
                    _NativeHandle,
                ],
                module.ReferenceCardinal128,
            )
            handle = factory(
                _WORKING_BITS,
                [float(value) for value in omega_nodes],
                ordering.copy(),
                packed_lower,
                packed_upper,
            )
            if type(handle.rank) is not int or handle.rank != rank:
                raise NumericalBreakdownError(
                    "Rust reference backend returned a wrong rank"
                )
        except NumericalBreakdownError:
            raise
        except Exception as error:
            raise NumericalBreakdownError(
                "Rust reference factor initialization failed; "
                "use reference_backend='python' to retain the Python path"
            ) from error
        self.rank = rank
        self._handle = handle

    def design_tau(self, normalized_times: FloatArray) -> FloatArray:
        """Return checked binary64 rows, with bounded native-call allocation."""
        if (
            normalized_times.ndim != 1
            or normalized_times.dtype.kind != "f"
            or not np.all(np.isfinite(normalized_times))
            or np.any(normalized_times < 0)
            or np.any(normalized_times > 1)
        ):
            raise ValueError("normalized times must be a finite float vector in [0,1]")
        if normalized_times.size == 0:
            return np.empty((0, self.rank), dtype=np.float64)
        batch_size = min(_BATCH_ROWS, _MAX_CELLS // self.rank)
        result = np.empty((normalized_times.size, self.rank), dtype=np.float64)
        for start in range(0, normalized_times.size, batch_size):
            stop = min(start + batch_size, normalized_times.size)
            try:
                payload = self._handle.design_tau(
                    [float(value) for value in normalized_times[start:stop]]
                )
            except Exception as error:
                raise NumericalBreakdownError(
                    "Rust reference time evaluation failed; "
                    "no Python fallback was substituted"
                ) from error
            expected = (stop - start) * self.rank * 8
            if type(payload) is not bytes or len(payload) != expected:
                raise NumericalBreakdownError(
                    "Rust reference backend returned an invalid array payload"
                )
            batch = np.frombuffer(payload, dtype="<f8").reshape(stop - start, self.rank)
            if not np.all(np.isfinite(batch)):
                raise NumericalBreakdownError(
                    "Rust reference backend returned nonfinite binary64 values"
                )
            result[start:stop] = batch
        return cast(npt.NDArray[np.float64], result)
