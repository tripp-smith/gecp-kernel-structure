import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from kernelgecp import _basis_io
from kernelgecp._basis_io import read_archive, write_archive


def _npy(array: np.ndarray, *, allow_pickle: bool = False) -> bytes:
    stream = io.BytesIO()
    np.save(stream, array, allow_pickle=allow_pickle)
    return stream.getvalue()


def _metadata(value: object = None) -> bytes:
    encoded = json.dumps({"schema_version": 1} if value is None else value).encode()
    return _npy(np.frombuffer(encoded, dtype=np.uint8))


def _zip(path: Path, members: list[tuple[str, bytes]], *, compression=0) -> None:
    with zipfile.ZipFile(path, "w", compression=compression) as archive:
        for name, content in members:
            archive.writestr(name, content)


def test_archive_round_trip_preserves_complex_dtypes_and_is_deterministic(
    tmp_path,
) -> None:
    metadata = {"schema_version": 1, "nested": [True, None, "λ", {"a": 1.25}]}
    arrays = {
        "complex": np.array([[1 + 2j, 3 - 4j]], dtype=">c16"),
        "single": np.array([2 + 3j], dtype=np.complex64),
        "fortran": np.asfortranarray(np.arange(12).reshape(3, 4)),
        "empty": np.empty((2, 0, 3)),
        "flag": np.array([True, False]),
        "scalar": np.array(2.5),
    }
    path = tmp_path / "basis.npz"
    second_path = tmp_path / "second.npz"
    write_archive(path, metadata, arrays)
    write_archive(second_path, metadata, dict(reversed(list(arrays.items()))))
    assert path.read_bytes() == second_path.read_bytes()
    restored_metadata, restored = read_archive(path)
    assert restored_metadata == metadata
    assert restored.keys() == arrays.keys()
    for key, array in arrays.items():
        assert restored[key].dtype == array.dtype
        np.testing.assert_array_equal(restored[key], array)
    assert restored["fortran"].flags.f_contiguous


@pytest.mark.parametrize(
    "metadata",
    [
        {"bad": np.float32(1)},
        {"bad": float("nan")},
        {"bad": float("inf")},
        {"bad": (1, 2)},
        {1: "bad"},
        [1, 2],
        {"bad": object()},
    ],
)
def test_writer_rejects_nonprimitive_metadata_without_overwrite(
    tmp_path, metadata
) -> None:
    path = tmp_path / "basis.npz"
    path.write_bytes(b"original")
    with pytest.raises(ValueError):
        write_archive(path, metadata, {"x": np.arange(3)})
    assert path.read_bytes() == b"original"
    assert sorted(item.name for item in tmp_path.iterdir()) == ["basis.npz"]


@pytest.mark.parametrize(
    "name", ["__metadata__", "../x", "x/y", "x.npy", "", "λ", "a" * 129, 3]
)
def test_writer_rejects_invalid_array_names(tmp_path, name) -> None:
    with pytest.raises(ValueError):
        write_archive(tmp_path / "bad.npz", {}, {name: np.arange(2)})


@pytest.mark.parametrize(
    "array",
    [
        np.array([object()], dtype=object),
        np.array(["a"]),
        np.array([b"a"]),
        np.array([(1, 2)], dtype=[("a", "i4"), ("b", "i4")]),
        np.array(["2026-01-01"], dtype="datetime64[D]"),
        [1, 2],
    ],
)
def test_writer_rejects_non_numeric_arrays(tmp_path, array) -> None:
    with pytest.raises(ValueError):
        write_archive(tmp_path / "bad.npz", {}, {"x": array})


def test_failed_replace_preserves_destination_and_cleans_temporary(
    tmp_path, monkeypatch
):
    path = tmp_path / "basis.npz"
    path.write_bytes(b"original")

    def fail_replace(*args):
        raise OSError("simulated replacement failure")

    monkeypatch.setattr(_basis_io.os, "replace", fail_replace)
    with pytest.raises(OSError, match="simulated"):
        write_archive(path, {}, {"x": np.arange(2)})
    assert path.read_bytes() == b"original"
    assert sorted(item.name for item in tmp_path.iterdir()) == ["basis.npz"]


def test_reader_rejects_duplicate_zip_members(tmp_path) -> None:
    path = tmp_path / "bad.npz"
    with pytest.warns(UserWarning, match="Duplicate name"):
        _zip(path, [("__metadata__.npy", _metadata())] * 2)
    with pytest.raises(ValueError, match="duplicate ZIP"):
        read_archive(path)


@pytest.mark.parametrize("name", ["../x.npy", "x/y.npy", "x.txt", "x.npy.npy"])
def test_reader_rejects_non_array_and_path_members(tmp_path, name) -> None:
    path = tmp_path / "bad.npz"
    _zip(path, [("__metadata__.npy", _metadata()), (name, _npy(np.arange(2)))])
    with pytest.raises(ValueError):
        read_archive(path)


@pytest.mark.parametrize(
    "contents",
    [
        b"not JSON",
        b"\xff",
        b"[]",
        b'{"schema_version":1,"schema_version":2}',
        b'{"bad":NaN}',
        b'{"bad":Infinity}',
        b'{"bad":1e999}',
    ],
)
def test_reader_rejects_invalid_metadata(tmp_path, contents) -> None:
    path = tmp_path / "bad.npz"
    _zip(path, [("__metadata__.npy", _npy(np.frombuffer(contents, dtype=np.uint8)))])
    with pytest.raises(ValueError):
        read_archive(path)


@pytest.mark.parametrize(
    "array", [np.array([123]), np.array([[123]], dtype=np.uint8), np.array(["{}"])]
)
def test_reader_rejects_metadata_with_wrong_array_type(tmp_path, array) -> None:
    path = tmp_path / "bad.npz"
    _zip(path, [("__metadata__.npy", _npy(array))])
    with pytest.raises(ValueError):
        read_archive(path)


def test_reader_rejects_object_array_before_pickle_load(tmp_path, monkeypatch) -> None:
    path = tmp_path / "bad.npz"
    _zip(
        path,
        [
            (
                "payload.npy",
                _npy(np.array([{"x": 1}], dtype=object), allow_pickle=True),
            ),
            ("__metadata__.npy", _metadata()),
        ],
    )

    def forbidden_load(*args, **kwargs):
        pytest.fail("object payload reached the NumPy loader")

    monkeypatch.setattr(_basis_io.np, "load", forbidden_load)
    with pytest.raises(ValueError, match="invalid NumPy"):
        read_archive(path)


def test_reader_checks_npy_shape_before_allocation(tmp_path, monkeypatch) -> None:
    path = tmp_path / "bad.npz"
    header = io.BytesIO()
    np.lib.format.write_array_header_1_0(
        header, {"descr": "<f8", "fortran_order": False, "shape": (2**60,)}
    )
    _zip(path, [("payload.npy", header.getvalue()), ("__metadata__.npy", _metadata())])

    def forbidden_load(*args, **kwargs):
        pytest.fail("oversized shape reached the NumPy loader")

    monkeypatch.setattr(_basis_io.np, "load", forbidden_load)
    with pytest.raises(ValueError, match="invalid NumPy"):
        read_archive(path)


def test_reader_disables_pickle_for_every_array(tmp_path, monkeypatch) -> None:
    path = tmp_path / "basis.npz"
    write_archive(path, {}, {"x": np.array([1 + 2j])})
    actual_load = np.load
    calls = []

    def checked_load(*args, **kwargs):
        calls.append(kwargs["allow_pickle"])
        return actual_load(*args, **kwargs)

    monkeypatch.setattr(_basis_io.np, "load", checked_load)
    read_archive(path)
    assert calls == [False, False]


def test_size_entry_and_compression_limits(tmp_path, monkeypatch) -> None:
    path = tmp_path / "basis.npz"
    write_archive(path, {}, {"x": np.arange(10)})
    monkeypatch.setattr(_basis_io, "_MAX_UNCOMPRESSED_BYTES", 100)
    with pytest.raises(ValueError, match="uncompressed"):
        read_archive(path)
    with pytest.raises(ValueError, match="resource|uncompressed"):
        write_archive(tmp_path / "oversized.npz", {}, {"x": np.arange(100)})
    monkeypatch.undo()
    _zip(path, [(f"x{index}.npy", b"") for index in range(65)])
    with pytest.raises(ValueError, match="resource"):
        read_archive(path)
    _zip(
        path,
        [("__metadata__.npy", _metadata()), ("bomb.npy", b"0" * 2_000_000)],
        compression=zipfile.ZIP_DEFLATED,
    )
    with pytest.raises(ValueError, match="compression ratio"):
        read_archive(path)


@pytest.mark.parametrize("corruption", [b"", b"not a ZIP file", b"PK\x03\x04"])
def test_reader_rejects_corrupt_containers(tmp_path, corruption) -> None:
    path = tmp_path / "bad.npz"
    path.write_bytes(corruption)
    with pytest.raises(ValueError):
        read_archive(path)


def test_reader_rejects_missing_metadata_and_trailing_npy_bytes(tmp_path) -> None:
    path = tmp_path / "bad.npz"
    _zip(path, [("x.npy", _npy(np.arange(2)))])
    with pytest.raises(ValueError, match="missing"):
        read_archive(path)
    _zip(path, [("__metadata__.npy", _metadata() + b"unexplained")])
    with pytest.raises(ValueError, match="invalid NumPy"):
        read_archive(path)
