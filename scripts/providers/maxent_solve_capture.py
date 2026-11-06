"""Capture bounded native solve failures without changing solver behavior."""

from __future__ import annotations

import hashlib
import inspect
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np

_MAX_BYTES = 1024 * 1024
_MAX_DIMENSION = 128


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def _scalar(value: Any) -> int | float | None:
    if isinstance(value, (bool, np.bool_)):
        return None
    if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
        return int(value) if isinstance(value, (int, np.integer)) else float(value)
    return None


def _capture(path: Path, metadata: dict[str, Any], a: Any, b: Any, frame: Any) -> None:
    loc = frame.f_locals
    filename = Path(frame.f_code.co_filename)
    metadata.update(
        caller_file=str(filename),
        caller_source_sha256=_sha(filename),
        caller_line=frame.f_lineno,
        caller_function=frame.f_code.co_name,
        iteration=_scalar(loc.get("i", loc.get("iteration"))),
        scalar_alpha=_scalar(loc.get("alpha")),
    )
    if metadata["scalar_alpha"] is None:
        function = loc.get("function")
        try:
            stored = object.__getattribute__(function, "__dict__")
        except (AttributeError, TypeError):
            stored = {}
        if isinstance(stored, dict):
            metadata["scalar_alpha"] = _scalar(stored.get("_alpha"))
    values = {"matrix": a, "rhs": b}
    for name in ("J", "Id", "mu", "v"):
        if name in loc and loc[name] is not None:
            values[name] = loc[name]
    arrays: dict[str, np.ndarray[Any, Any]] = {}
    total = 0
    for name, value in values.items():
        if not isinstance(value, (np.ndarray, np.number, int, float, complex)):
            raise ValueError(f"unsupported capture array type: {name}")
        array = np.asarray(value)
        if (
            array.dtype.kind not in "iufc"
            or array.ndim > 2
            or any(size > _MAX_DIMENSION for size in array.shape)
        ):
            raise ValueError(f"unsupported capture array shape or dtype: {name}")
        total += array.nbytes
        if total > _MAX_BYTES:
            raise ValueError("capture total array byte cap exceeded")
        arrays[name] = array
    metadata["arrays"] = {
        name: {
            "shape": list(array.shape),
            "dtype": str(array.dtype),
            "bytes": array.nbytes,
            "finite": bool(np.isfinite(array).all()),
            "sha256": hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest(),
        }
        for name, array in arrays.items()
    }
    metadata["array_bytes"] = total
    # Exclusive creation prevents replacing earlier evidence.
    with path.open("xb") as file:
        np.savez_compressed(file, **arrays)
    metadata.update(npz=str(path), npz_sha256=_sha(path), status="captured")


@contextmanager
def capture_failed_solve(path: str | Path, metadata: dict[str, Any]):
    """Observe a failed solve once, restoring the previous callable on exit.

    The captured data are diagnostic only. The original solve result and
    exceptions are passed through, including when evidence capture itself fails.
    This process-global patch is intended for a single-threaded native worker.
    """
    previous = np.linalg.solve

    def wrapped(a: Any, b: Any, *args: Any, **kwargs: Any) -> Any:
        try:
            return previous(a, b, *args, **kwargs)
        except np.linalg.LinAlgError:
            try:
                frame = inspect.currentframe()
                if frame is None or frame.f_back is None:
                    raise RuntimeError("native caller frame unavailable")
                metadata["status"] = "capture_failed"
                _capture(Path(path), metadata, a, b, frame.f_back)
            except BaseException as error:
                try:
                    metadata["capture_error"] = f"{type(error).__name__}: {error}"
                    metadata["status"] = "capture_failed"
                except BaseException:
                    pass
            raise

    np.linalg.solve = wrapped
    try:
        yield
    finally:
        np.linalg.solve = previous
