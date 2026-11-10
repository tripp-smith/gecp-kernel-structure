import hashlib
import io
import stat
import subprocess
import zipfile

import pytest

from scripts import release_science_payloads as payloads
from scripts.release_science_runs import REPOSITORY, ROLES, SCOPE


def zip_bytes(entries):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as archive:
        for name, content in entries:
            archive.writestr(name, content)
    return buf.getvalue()


def selection(blob):
    runs = {}
    for index, (role, (workflow, _, templates)) in enumerate(ROLES.items(), 1):
        runs[role] = {
            "workflow": workflow,
            "id": index,
            "attempt": 1,
            "artifacts": [],
            "auxiliary_artifacts": [],
        }
        for template in templates:
            name = template.format(id=index, attempt=1)
            runs[role]["artifacts"].append(
                {
                    "id": sum(len(r["artifacts"]) for r in runs.values()) + 1,
                    "name": name,
                    "size_in_bytes": len(blob),
                    "digest": "sha256:" + hashlib.sha256(blob).hexdigest(),
                }
            )
    return {
        "repository": REPOSITORY,
        "revision": "a" * 40,
        "scope": SCOPE,
        "runs": runs,
    }


def fake(blob):
    def get(_ident, path, log):
        path.write_bytes(blob)
        log.write_text("")

    return get


def test_materialize(tmp_path):
    blob = zip_bytes([("folder/file", "hello")])
    selected = selection(blob)
    result = payloads.materialize(
        selected, tmp_path / "out", tmp_path / "receipt.json", fake(blob)
    )
    assert result["status"] == "materialized"
    assert len(result["artifacts"]) == sum(len(x[2]) for x in ROLES.values())
    assert (
        tmp_path / "out/providers/provider-green-linux/folder/file"
    ).read_text() == "hello"
    assert result["provenance"].startswith("caller must")


@pytest.mark.parametrize("mutation", ["digest", "size", "missing", "duplicate"])
def test_metadata_and_integrity_fail(tmp_path, mutation):
    blob = zip_bytes([("file", "x")])
    data = selection(blob)
    records = data["runs"]["providers"]["artifacts"]
    if mutation == "digest":
        records[0]["digest"] = "sha256:" + "0" * 64
    elif mutation == "size":
        records[0]["size_in_bytes"] += 1
    elif mutation == "missing":
        records.pop()
    else:
        records[1]["id"] = records[0]["id"]
    with pytest.raises(ValueError):
        payloads.materialize(data, tmp_path / "out", tmp_path / "receipt", fake(blob))
    if mutation in ("digest", "size"):
        assert (tmp_path / "out/failure.log").exists()
        with pytest.raises(ValueError, match="must be new"):
            payloads.materialize(
                data, tmp_path / "out", tmp_path / "receipt", fake(blob)
            )


@pytest.mark.parametrize(
    "entries",
    [
        [("../x", "x")],
        [("/x", "x")],
        [("a\\b", "x")],
        [("a", "x"), ("a/b", "x")],
        [("a", "x"), ("a", "y")],
    ],
)
def test_unsafe_inventory(tmp_path, entries):
    blob = zip_bytes(entries)
    with pytest.raises(ValueError):
        payloads.materialize(
            selection(blob), tmp_path / "out", tmp_path / "receipt", fake(blob)
        )
    assert (tmp_path / "out/failure.log").exists()
    assert not (tmp_path / "out/providers/provider-triqs-ubuntu-24.04").exists()


def test_symlink_and_limits(monkeypatch):
    entry = zipfile.ZipInfo("link")
    entry.create_system = 3
    entry.external_attr = (stat.S_IFLNK | 0o777) << 16
    blob = zip_bytes([(entry, "target")])
    with (
        zipfile.ZipFile(io.BytesIO(blob)) as archive,
        pytest.raises(ValueError, match="nonregular"),
    ):
        payloads.inventory(archive)
    blob = zip_bytes([("file", "hello")])
    monkeypatch.setattr(payloads, "MAX_MEMBER", 4)
    with (
        zipfile.ZipFile(io.BytesIO(blob)) as archive,
        pytest.raises(ValueError, match="budget"),
    ):
        payloads.inventory(archive)


def test_download_failure(tmp_path, monkeypatch):
    class Process:
        stdout = io.BytesIO(b"partial")

        def wait(self):
            return 1

        def kill(self):
            pass

    monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: Process())
    with pytest.raises(ValueError, match="download failed"):
        payloads.download(1, tmp_path / "x.zip", tmp_path / "log")
    assert (tmp_path / "x.zip").read_bytes() == b"partial"


def test_bad_name_cli(tmp_path):
    import json
    import sys

    data = selection(zip_bytes([("file", "x")]))
    data["runs"]["providers"]["artifacts"][0]["name"] = []
    selected = tmp_path / "selected.json"
    selected.write_text(json.dumps(data))
    run = subprocess.run(
        [
            sys.executable,
            "scripts/release_science_payloads.py",
            "materialize",
            "--selected-runs",
            str(selected),
            "--output-root",
            str(tmp_path / "out"),
            "--receipt",
            str(tmp_path / "receipt"),
        ],
        capture_output=True,
        text=True,
    )
    assert run.returncode != 0
    assert "invalid artifact name" in run.stderr
    assert "Traceback" not in run.stderr


def test_unsupported_method(tmp_path):
    blob = bytearray(zip_bytes([("file", "x")]))
    central = blob.index(b"PK\x01\x02")
    blob[8:10] = (99).to_bytes(2, "little")
    blob[central + 10 : central + 12] = (99).to_bytes(2, "little")
    with pytest.raises(ValueError, match="unsupported ZIP compression"):
        payloads.materialize(
            selection(blob), tmp_path / "out", tmp_path / "receipt", fake(blob)
        )
    assert "unsupported ZIP compression" in (tmp_path / "out/failure.log").read_text()


def test_preflight_bounds_before_zipfile(tmp_path, monkeypatch):
    import struct

    blob = bytearray(zip_bytes([("file", "x")]))
    end = blob.index(b"PK\x05\x06")
    monkeypatch.setattr(
        payloads.zipfile, "ZipFile", lambda *_: pytest.fail("ZipFile called")
    )
    for field, value in ((8, 20001), (12, payloads.MAX_CENTRAL + 1)):
        modified = bytearray(blob)
        if field == 8:
            struct.pack_into("<HH", modified, end + field, value, value)
        else:
            struct.pack_into("<I", modified, end + field, value)
        path = tmp_path / str(field)
        path.write_bytes(modified)
        with pytest.raises(ValueError, match="too many|central directory"):
            payloads.preflight(path)


def test_valid_zip64(tmp_path):
    import struct

    blob = zip_bytes([("file", "x")])
    end = blob.index(b"PK\x05\x06")
    fields = struct.unpack_from("<4s4H2IH", blob, end)
    _, _, _, _, count, cd_size, cd_offset, _ = fields
    record = struct.pack(
        "<4sQ2H2I4Q", b"PK\x06\x06", 44, 45, 45, 0, 0, count, count, cd_size, cd_offset
    )
    locator = struct.pack("<4sIQI", b"PK\x06\x07", 0, end, 1)
    footer = struct.pack(
        "<4s4H2IH", b"PK\x05\x06", 0, 0, 0xFFFF, 0xFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0
    )
    archive = blob[:end] + record + locator + footer
    path = tmp_path / "zip64.zip"
    path.write_bytes(archive)
    assert payloads.preflight(path) == 1
    with zipfile.ZipFile(path) as z:
        assert payloads.inventory(z)[1] == 1


def test_hidden_zip64_locator_rejected_before_zipfile(tmp_path, monkeypatch):
    import struct

    blob = bytearray(zip_bytes([("file", "x")]))
    central = blob.index(b"PK\x01\x02")
    end = blob.index(b"PK\x05\x06")
    _, _, _, _, count, cd_size, cd_offset, _ = struct.unpack_from("<4s4H2IH", blob, end)
    record = struct.pack(
        "<4sQ2H2I4Q",
        b"PK\x06\x06",
        44,
        45,
        45,
        0,
        0,
        count,
        count,
        cd_size,
        cd_offset,
    )
    locator = struct.pack("<4sIQI", b"PK\x06\x07", 0, end, 1)
    inserted = record + locator
    struct.pack_into("<H", blob, central + 32, len(inserted))
    struct.pack_into("<I", blob, end + 12, cd_size + len(inserted))
    modified = blob[:end] + inserted + blob[end:]
    path = tmp_path / "hidden.zip"
    path.write_bytes(modified)
    monkeypatch.setattr(
        payloads.zipfile, "ZipFile", lambda *_: pytest.fail("ZipFile called")
    )
    with pytest.raises(ValueError, match="unexpected ZIP64 locator"):
        payloads.preflight(path)


def test_corrupt_deflate_is_clear_failure(tmp_path):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("file", b"actual compressed payload")
    blob = bytearray(buffer.getvalue())
    # Reserved DEFLATE block type: metadata and selected archive hash remain
    # internally consistent, but decoding the actual member must fail.
    blob[30 + len("file")] = 0x07
    with pytest.raises(ValueError, match="invalid ZIP archive"):
        payloads.materialize(
            selection(blob), tmp_path / "out", tmp_path / "receipt", fake(blob)
        )
    assert "invalid ZIP archive" in (tmp_path / "out/failure.log").read_text()
