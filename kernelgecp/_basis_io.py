"""Bounded, pickle-free NPZ containers for the public basis facade.

The facade owns metadata versions and semantic schemas. This module enforces
the transport boundary: primitive JSON, numeric arrays, bounded resources,
and transactional replacement of an existing archive.
"""

from __future__ import annotations

import io
import json
import math
import os
import re
import struct
import tempfile
import zipfile
import zlib
from collections.abc import Mapping
from pathlib import Path
from typing import Any, BinaryIO, cast

import numpy as np
import numpy.typing as npt

_METADATA_KEY = "__metadata__"
_MAX_MEMBERS = 64
_MAX_UNCOMPRESSED_BYTES = 64 * 1024 * 1024
_MAX_METADATA_BYTES = 1024 * 1024
_MAX_HEADER_BYTES = 16 * 1024
_MAX_COMPRESSION_RATIO = 1000
_MAX_JSON_DEPTH = 64
_MAX_KEY_LENGTH = 128
_MAX_CENTRAL_DIRECTORY_BYTES = 64 * 1024
_MAX_CONTAINER_BYTES = _MAX_UNCOMPRESSED_BYTES + 1024 * 1024
_MEMBER_KEY = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z", flags=re.ASCII)


def _validate_key(key: object, *, metadata: bool = False) -> str:
    if (
        not isinstance(key, str)
        or len(key) > _MAX_KEY_LENGTH
        or not _MEMBER_KEY.fullmatch(key)
        or (key == _METADATA_KEY and not metadata)
    ):
        raise ValueError(
            "archive array names must be identifiers excluding __metadata__"
        )
    return key


def _validate_json(value: object, *, depth: int = 0) -> None:
    if depth > _MAX_JSON_DEPTH:
        raise ValueError("archive metadata is nested too deeply")
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float:
        if not math.isfinite(value):
            raise ValueError("archive metadata numbers must be finite")
        return
    if type(value) is list:
        for item in value:
            _validate_json(item, depth=depth + 1)
        return
    if type(value) is dict:
        for key, item in value.items():
            if type(key) is not str:
                raise ValueError("archive metadata object keys must be strings")
            _validate_json(item, depth=depth + 1)
        return
    raise ValueError("archive metadata must contain only primitive JSON values")


def _metadata_bytes(metadata: dict[str, object]) -> bytes:
    if type(metadata) is not dict:
        raise ValueError("archive metadata root must be a JSON object")
    _validate_json(metadata)
    try:
        encoded = json.dumps(
            metadata,
            allow_nan=False,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    except (ValueError, TypeError, RecursionError) as error:
        raise ValueError(
            "archive metadata cannot be encoded as primitive JSON"
        ) from error
    if len(encoded) > _MAX_METADATA_BYTES:
        raise ValueError("archive metadata exceeds the size limit")
    return encoded


def _validate_dtype(dtype: np.dtype[Any]) -> None:
    if (
        dtype.fields is not None
        or dtype.subdtype is not None
        or dtype.kind not in "biufc"
    ):
        raise ValueError("archive arrays must have primitive numeric or boolean dtypes")


def _npy_bytes(array: npt.NDArray[Any]) -> bytes:
    if not isinstance(array, np.ndarray):
        raise ValueError("archive payloads must be NumPy arrays")
    _validate_dtype(array.dtype)
    if array.ndim > 64 or array.nbytes > _MAX_UNCOMPRESSED_BYTES:
        raise ValueError("archive array exceeds the resource limits")
    stream = io.BytesIO()
    np.lib.format.write_array(stream, array, allow_pickle=False)
    return stream.getvalue()


def write_archive(
    path: str | os.PathLike[str],
    metadata: dict[str, object],
    arrays: Mapping[str, npt.NDArray[Any]],
) -> None:
    """Atomically write deterministic NPZ bytes without modifying caller data.

    The parent directory must exist. Any validation or serialization failure
    leaves an existing destination intact. Metadata schema/version validation
    belongs to the calling facade.
    """

    if not isinstance(arrays, Mapping):
        raise ValueError("archive arrays must be a mapping")
    if len(arrays) + 1 > _MAX_MEMBERS:
        raise ValueError("archive contains too many arrays")
    encoded_metadata = _metadata_bytes(metadata)
    members = {
        _METADATA_KEY: _npy_bytes(np.frombuffer(encoded_metadata, dtype=np.uint8))
    }
    total_bytes = len(members[_METADATA_KEY])
    for key in sorted(_validate_key(key) for key in arrays):
        member = _npy_bytes(arrays[key])
        total_bytes += len(member)
        if total_bytes > _MAX_UNCOMPRESSED_BYTES:
            raise ValueError("archive exceeds the total uncompressed size limit")
        members[key] = member
    target = Path(path)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b",
            prefix=f".{target.name}.",
            suffix=".tmp",
            dir=target.parent,
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            with zipfile.ZipFile(
                stream, mode="w", compression=zipfile.ZIP_STORED
            ) as archive:
                for key, member in members.items():
                    info = zipfile.ZipInfo(
                        f"{key}.npy", date_time=(1980, 1, 1, 0, 0, 0)
                    )
                    info.compress_type = zipfile.ZIP_STORED
                    info.create_system = 3
                    info.external_attr = 0o600 << 16
                    archive.writestr(info, member)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("archive metadata contains duplicate JSON keys")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> object:
    raise ValueError(f"archive metadata contains a nonfinite JSON number: {value}")


def _read_metadata(array: npt.NDArray[Any]) -> dict[str, object]:
    if array.dtype != np.dtype(np.uint8) or array.ndim != 1:
        raise ValueError("archive metadata must be one uint8 UTF-8 JSON array")
    if array.size > _MAX_METADATA_BYTES:
        raise ValueError("archive metadata exceeds the size limit")
    try:
        metadata = json.loads(
            array.tobytes().decode("utf-8"),
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_json_constant,
        )
    except (ValueError, RecursionError) as error:
        raise ValueError("archive metadata is invalid UTF-8 JSON") from error
    if type(metadata) is not dict:
        raise ValueError("archive metadata root must be a JSON object")
    _validate_json(metadata)
    return cast(dict[str, object], metadata)


def _read_npy(member: bytes) -> npt.NDArray[Any]:
    """Validate the NPY header and byte count before NumPy can allocate data."""

    stream = io.BytesIO(member)
    try:
        version = np.lib.format.read_magic(stream)
        if version == (1, 0):
            shape, _, dtype = np.lib.format.read_array_header_1_0(
                stream, max_header_size=_MAX_HEADER_BYTES
            )
        elif version == (2, 0):
            shape, _, dtype = np.lib.format.read_array_header_2_0(
                stream, max_header_size=_MAX_HEADER_BYTES
            )
        else:
            raise ValueError("archive arrays require NPY format version 1 or 2")
        _validate_dtype(dtype)
        if len(shape) > 64 or any(
            type(dimension) is not int
            or dimension < 0
            or dimension > np.iinfo(np.intp).max
            for dimension in shape
        ):
            raise ValueError("archive array has an invalid shape")
        expected_bytes = math.prod(shape) * dtype.itemsize
        if (
            expected_bytes > _MAX_UNCOMPRESSED_BYTES
            or expected_bytes != len(member) - stream.tell()
        ):
            raise ValueError("archive array byte count does not match its shape")
        array = np.load(
            io.BytesIO(member), allow_pickle=False, max_header_size=_MAX_HEADER_BYTES
        )
        if not isinstance(array, np.ndarray):
            raise ValueError("archive payload is not a single NumPy array")
        return array
    except (ValueError, EOFError, OverflowError) as error:
        raise ValueError("archive contains an invalid NumPy array") from error


def _preflight_zip(stream: BinaryIO) -> None:
    """Bound the central directory before ZipFile allocates its member table."""

    stream.seek(0, os.SEEK_END)
    size = stream.tell()
    if size > _MAX_CONTAINER_BYTES:
        raise ValueError("archive container exceeds the size limit")
    stream.seek(max(0, size - 65_557))
    tail = stream.read()
    offset = tail.rfind(b"PK\x05\x06")
    if offset < 0 or len(tail) - offset < 22:
        raise ValueError("archive is missing its ZIP end record")
    record = struct.unpack("<4s4H2LH", tail[offset : offset + 22])
    (
        _,
        disk,
        directory_disk,
        disk_count,
        count,
        directory_size,
        directory_offset,
        comment,
    ) = record
    if (
        disk != 0
        or directory_disk != 0
        or disk_count != count
        or count < 1
        or count > _MAX_MEMBERS
        or directory_size > _MAX_CENTRAL_DIRECTORY_BYTES
        or directory_offset + directory_size > size
        or offset + 22 + comment != len(tail)
        or tail[max(0, offset - 20) : offset - 16] == b"PK\x06\x07"
    ):
        raise ValueError("archive ZIP directory exceeds the supported resource limits")
    stream.seek(0)


def read_archive(
    path: str | os.PathLike[str],
) -> tuple[dict[str, object], dict[str, npt.NDArray[Any]]]:
    """Read bounded numeric arrays and primitive metadata with pickle disabled.

    Duplicate/path-like ZIP members, unsupported compression, excessive sizes,
    and compression ratios above 1000 are rejected before loading arrays.
    Returned metadata is untrusted until the facade validates its schema.
    """

    try:
        with Path(path).open("rb") as source:
            _preflight_zip(source)
            return _read_zip(source)
    except (zipfile.BadZipFile, EOFError, RuntimeError, zlib.error) as error:
        raise ValueError("invalid or corrupt basis archive") from error


def _read_zip(
    source: BinaryIO,
) -> tuple[dict[str, object], dict[str, npt.NDArray[Any]]]:
    with zipfile.ZipFile(source, mode="r") as archive:
        entries = archive.infolist()
        if not entries or len(entries) > _MAX_MEMBERS:
            raise ValueError("archive has an invalid number of members")
        names = [entry.filename for entry in entries]
        if len(set(names)) != len(names):
            raise ValueError("archive contains duplicate ZIP members")
        total_bytes = 0
        keyed_entries: dict[str, zipfile.ZipInfo] = {}
        for entry in entries:
            if not entry.filename.endswith(".npy"):
                raise ValueError("archive members must be named NumPy arrays")
            key = _validate_key(entry.filename[:-4], metadata=True)
            if (
                entry.is_dir()
                or entry.flag_bits & 1
                or entry.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED)
                or entry.file_size < 0
                or entry.compress_size < 0
            ):
                raise ValueError("archive contains an unsupported ZIP member")
            if entry.file_size > _MAX_COMPRESSION_RATIO * max(entry.compress_size, 1):
                raise ValueError("archive member exceeds the compression ratio limit")
            total_bytes += entry.file_size
            if total_bytes > _MAX_UNCOMPRESSED_BYTES:
                raise ValueError("archive exceeds the total uncompressed size limit")
            keyed_entries[key] = entry
        if _METADATA_KEY not in keyed_entries:
            raise ValueError("archive is missing its metadata member")
        arrays: dict[str, npt.NDArray[Any]] = {}
        metadata: dict[str, object] | None = None
        for key, entry in keyed_entries.items():
            with archive.open(entry, mode="r") as stream:
                member = stream.read(_MAX_UNCOMPRESSED_BYTES + 1)
            if len(member) != entry.file_size or len(member) > _MAX_UNCOMPRESSED_BYTES:
                raise ValueError("archive ZIP member size is inconsistent")
            array = _read_npy(member)
            if key == _METADATA_KEY:
                metadata = _read_metadata(array)
            else:
                arrays[key] = array
        if metadata is None:
            raise ValueError("archive is missing its metadata member")
        return metadata, arrays
