"""Materialize selected science artifacts; input provenance is the caller's duty.

The caller must freshly fetch and authenticate selection metadata. This module
checks consistency and ZIP integrity; it does not establish scientific eligibility.
"""

import argparse
import hashlib
import json
import re
import stat
import struct
import subprocess
import zipfile
import zlib
from pathlib import Path

try:
    from scripts.providers.verify_receipt import read_receipt
    from scripts.release_science_runs import REPOSITORY, SCOPE, profile_roles
except ModuleNotFoundError:
    from providers.verify_receipt import read_receipt
    from release_science_runs import REPOSITORY, SCOPE, profile_roles

MAX_ZIP = 2 * 2**30
MAX_MEMBER = 2 * 2**30
MAX_TOTAL = 4 * 2**30
MAX_ENTRIES = 20_000
MAX_CENTRAL = 16 * 2**20
MAX_BUNDLE = 8 * 2**30
CHUNK = 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def positive(value):
    return type(value) is int and value > 0


def validate_selection(data):
    require(
        isinstance(data, dict) and data.get("repository") == REPOSITORY,
        "unexpected repository",
    )
    require(data.get("scope") == SCOPE, "unexpected selection scope")
    require(
        isinstance(data.get("revision"), str)
        and re.fullmatch(r"[0-9a-f]{40}", data["revision"]),
        "invalid revision",
    )
    roles = profile_roles(data.get("profile"))
    if "profile" in data:
        require(data["profile"] is not None, "invalid profile")
    runs = data.get("runs")
    require(isinstance(runs, dict) and set(runs) == set(roles), "wrong roles")
    ids, names = set(), set()
    for role, (workflow, _, templates) in roles.items():
        run = runs[role]
        require(
            isinstance(run, dict)
            and run.get("workflow") == workflow
            and positive(run.get("id"))
            and positive(run.get("attempt")),
            "invalid run",
        )
        expected = {s.format(id=run["id"], attempt=run["attempt"]) for s in templates}
        artifacts, auxiliary = run.get("artifacts"), run.get("auxiliary_artifacts")
        require(
            isinstance(artifacts, list) and isinstance(auxiliary, list),
            "invalid artifacts",
        )
        require(
            all(isinstance(a, dict) for a in artifacts + auxiliary),
            "invalid artifact record",
        )
        require(
            all(isinstance(a.get("name"), str) for a in artifacts + auxiliary),
            "invalid artifact name",
        )
        require(
            len(artifacts) == len(expected)
            and {a.get("name") for a in artifacts} == expected,
            "wrong artifact names",
        )
        limit = {"providers": 2, "molecular": 1, "epw": 1}.get(role, 0)
        require(len(auxiliary) <= limit, "too many auxiliary artifacts")
        for item in artifacts + auxiliary:
            name, ident, size, digest = (
                item.get(k) for k in ("name", "id", "size_in_bytes", "digest")
            )
            require(
                isinstance(name, str)
                and positive(ident)
                and positive(size)
                and isinstance(digest, str)
                and re.fullmatch(r"sha256:[0-9a-f]{64}", digest),
                "invalid artifact metadata",
            )
            require(
                ident not in ids and name not in names, "duplicate artifact ID or name"
            )
            ids.add(ident)
            names.add(name)
            if item in auxiliary:
                require(
                    re.fullmatch(
                        r"tripp-smith~gecp-kernel-structure~[A-Za-z0-9_-]+\.dockerbuild",
                        name,
                    ),
                    "invalid auxiliary artifact",
                )
            else:
                require(size <= MAX_ZIP, "ZIP exceeds byte budget")
    return data


def download(ident, target, log):
    """Stream gh output to disk with a finite cap, preserving diagnostic stderr."""
    with target.open("xb") as out, log.open("xb") as err:
        proc = subprocess.Popen(
            ["gh", "api", f"repos/{REPOSITORY}/actions/artifacts/{ident}/zip"],
            stdout=subprocess.PIPE,
            stderr=err,
        )
        try:
            size = 0
            while chunk := proc.stdout.read(CHUNK):
                size += len(chunk)
                if size > MAX_ZIP:
                    raise ValueError("ZIP exceeds byte budget")
                out.write(chunk)
            require(proc.wait() == 0, "artifact download failed; see download log")
        except BaseException:
            proc.kill()
            proc.wait()
            raise


def preflight(path):
    """Bound the central directory before ZipFile allocates its inventory."""
    size = path.stat().st_size
    require(size <= MAX_ZIP and size >= 22, "invalid ZIP size")
    with path.open("rb") as stream:
        tail_size = min(size, 22 + 65535)
        stream.seek(size - tail_size)
        tail = stream.read(tail_size)
        marker = tail.rfind(b"PK\x05\x06")
        require(marker >= 0 and marker + 22 <= len(tail), "missing ZIP end record")
        end = struct.unpack_from("<4s4H2IH", tail, marker)
        _, disk, cd_disk, disk_count, count, cd_size, cd_offset, comment = end
        require(marker + 22 + comment == len(tail), "invalid ZIP end record")
        end_offset = size - tail_size + marker
        require(
            disk == 0 and cd_disk == 0 and disk_count == count,
            "multi-disk ZIP unsupported",
        )
        is64 = count == 0xFFFF or cd_size == 0xFFFFFFFF or cd_offset == 0xFFFFFFFF
        locator = b""
        if end_offset >= 20:
            stream.seek(end_offset - 20)
            locator = stream.read(20)
        has_locator = locator[:4] == b"PK\x06\x07"
        require(is64 or not has_locator, "unexpected ZIP64 locator")
        if is64:
            require(has_locator, "missing ZIP64 locator")
            sig, disk64, offset64, disks = struct.unpack("<4sIQI", locator)
            require(
                sig == b"PK\x06\x07" and disk64 == 0 and disks == 1,
                "invalid ZIP64 locator",
            )
            require(offset64 + 56 <= end_offset - 20, "invalid ZIP64 end offset")
            stream.seek(offset64)
            record = stream.read(56)
            fields = struct.unpack("<4sQ2H2I4Q", record)
            (
                sig,
                length,
                _,
                _,
                disk64,
                cd_disk64,
                disk_count64,
                count64,
                cd_size64,
                cd_offset64,
            ) = fields
            require(
                sig == b"PK\x06\x06"
                and length == 44
                and offset64 + 56 == end_offset - 20,
                "unsupported ZIP64 end record",
            )
            require(
                disk64 == 0 and cd_disk64 == 0 and disk_count64 == count64,
                "multi-disk ZIP64 unsupported",
            )
            require(
                count in (0xFFFF, count64)
                and cd_size in (0xFFFFFFFF, cd_size64)
                and cd_offset in (0xFFFFFFFF, cd_offset64),
                "inconsistent ZIP64 footer",
            )
            count, cd_size, cd_offset, end_offset = (
                count64,
                cd_size64,
                cd_offset64,
                offset64,
            )
        require(count <= MAX_ENTRIES, "too many ZIP entries")
        require(cd_size <= MAX_CENTRAL, "ZIP central directory exceeds byte budget")
        require(
            cd_offset <= end_offset and cd_size == end_offset - cd_offset,
            "invalid ZIP central directory bounds",
        )
        stream.seek(cd_offset)
        actual = 0
        while stream.tell() < end_offset:
            header = stream.read(46)
            require(
                len(header) == 46 and header[:4] == b"PK\x01\x02",
                "invalid ZIP central directory entry",
            )
            name_len, extra_len, comment_len = struct.unpack_from("<3H", header, 28)
            next_offset = stream.tell() + name_len + extra_len + comment_len
            require(
                next_offset <= end_offset, "invalid ZIP central directory entry bounds"
            )
            actual += 1
            require(actual <= MAX_ENTRIES, "too many ZIP entries")
            stream.seek(next_offset)
        require(actual == count, "ZIP entry count mismatch")
    return count


def inventory(archive):
    entries = archive.infolist()
    require(len(entries) <= MAX_ENTRIES, "too many ZIP entries")
    paths, files, directories = set(), set(), set()
    total = 0
    for entry in entries:
        name = entry.filename
        require("\x00" not in entry.orig_filename, "unsafe ZIP path")
        require(
            entry.compress_type in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED),
            "unsupported ZIP compression method",
        )
        require(
            name
            and "\\" not in name
            and not name.startswith("/")
            and "\x00" not in name,
            "unsafe ZIP path",
        )
        parts = name.rstrip("/").split("/")
        require(
            all(p not in ("", ".", "..") for p in parts)
            and not re.match(r"^[A-Za-z]:", parts[0]),
            "unsafe ZIP path",
        )
        key = "/".join(parts)
        require(key not in paths, "duplicate ZIP path")
        paths.add(key)
        require(not entry.flag_bits & 1, "encrypted ZIP member")
        mode = entry.external_attr >> 16
        kind = stat.S_IFMT(mode)
        is_dir = entry.is_dir()
        require(not is_dir or entry.file_size == 0, "directory contains data")
        require(
            kind in (0, stat.S_IFDIR if is_dir else stat.S_IFREG),
            "nonregular ZIP member",
        )
        require(entry.file_size <= MAX_MEMBER, "ZIP member exceeds byte budget")
        total += entry.file_size
        require(total <= MAX_TOTAL, "ZIP uncompressed byte budget exceeded")
        parents = {"/".join(parts[:i]) for i in range(1, len(parts))}
        require(not parents.intersection(files), "ZIP directory/file collision")
        if is_dir:
            require(key not in files, "ZIP directory/file collision")
            directories.add(key)
        else:
            require(key not in directories, "ZIP directory/file collision")
            files.add(key)
        directories.update(parents)
    return entries, total


def materialize(selected, output_root, receipt, downloader=download, *, role=None):
    validate_selection(selected)
    roles = profile_roles(selected.get("profile"))
    if role is not None:
        require(role in roles, "unknown selected role")
        roles = {role: roles[role]}
    output_root, receipt = Path(output_root), Path(receipt)
    require(
        not output_root.exists() and not output_root.is_symlink(),
        "output root must be new",
    )
    require(not receipt.exists() and not receipt.is_symlink(), "receipt must be new")
    output_root.mkdir(parents=False)
    records, total_bundle = [], 0
    try:
        for role in roles:
            for item in selected["runs"][role]["artifacts"]:
                role_dir = output_root / role
                role_dir.mkdir(exist_ok=True)
                zip_path = role_dir / (item["name"] + ".zip")
                log = role_dir / (item["name"] + ".download.log")
                downloader(item["id"], zip_path, log)
                h = hashlib.sha256()
                size = 0
                with zip_path.open("rb") as stream:
                    while chunk := stream.read(CHUNK):
                        size += len(chunk)
                        require(size <= MAX_ZIP, "ZIP exceeds byte budget")
                        h.update(chunk)
                require(
                    size == item["size_in_bytes"]
                    and "sha256:" + h.hexdigest() == item["digest"],
                    "ZIP size or digest mismatch",
                )
                dest = role_dir / item["name"]
                count = preflight(zip_path)
                with zipfile.ZipFile(zip_path) as archive:
                    entries, uncompressed = inventory(archive)
                    require(len(entries) == count, "ZIP entry count mismatch")
                    prior_bundle = total_bundle
                    total_bundle += uncompressed
                    require(
                        total_bundle <= MAX_BUNDLE,
                        "bundle uncompressed byte budget exceeded",
                    )
                    dest.mkdir()
                    extracted = 0
                    for entry in entries:
                        path = dest.joinpath(*entry.filename.rstrip("/").split("/"))
                        if entry.is_dir():
                            path.mkdir(parents=True, exist_ok=True)
                        else:
                            path.parent.mkdir(parents=True, exist_ok=True)
                            with (
                                archive.open(entry) as source,
                                path.open("xb") as target,
                            ):
                                member_bytes = 0
                                while chunk := source.read(CHUNK):
                                    member_bytes += len(chunk)
                                    extracted += len(chunk)
                                    require(
                                        member_bytes <= entry.file_size
                                        and member_bytes <= MAX_MEMBER,
                                        "ZIP member byte count exceeded",
                                    )
                                    require(
                                        extracted <= uncompressed
                                        and extracted <= MAX_TOTAL
                                        and prior_bundle + extracted <= MAX_BUNDLE,
                                        "ZIP extracted byte budget exceeded",
                                    )
                                    target.write(chunk)
                                require(
                                    member_bytes == entry.file_size,
                                    "ZIP member byte count mismatch",
                                )
                    require(
                        extracted == uncompressed, "ZIP extracted byte count mismatch"
                    )
                records.append(
                    {
                        "role": role,
                        "name": item["name"],
                        "id": item["id"],
                        "archive": str(zip_path),
                        "directory": str(dest),
                        "digest": item["digest"],
                        "size_in_bytes": size,
                        "entries": len(entries),
                        "uncompressed_bytes": uncompressed,
                    }
                )
        result = {
            "status": "materialized",
            "scope": (
                "archive integrity and safe extraction only; "
                "scientific eligibility not established"
            ),
            "provenance": (
                "caller must freshly fetch and authenticate selected metadata"
            ),
            "repository": REPOSITORY,
            "revision": selected["revision"],
            "artifacts": records,
            "auxiliary_artifacts": {
                role: selected["runs"][role]["auxiliary_artifacts"] for role in roles
            },
            "total_uncompressed_bytes": total_bundle,
        }
        if "profile" in selected:
            result["profile"] = selected["profile"]
        with receipt.open("x") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
        return result
    except (
        zipfile.BadZipFile,
        zipfile.LargeZipFile,
        zlib.error,
        NotImplementedError,
        EOFError,
        RuntimeError,
    ) as error:
        with (output_root / "failure.log").open("x") as stream:
            stream.write(f"invalid ZIP archive: {type(error).__name__}: {error}\n")
        raise ValueError(f"invalid ZIP archive: {error}") from error
    except BaseException as error:
        with (output_root / "failure.log").open("x") as stream:
            stream.write(f"{type(error).__name__}: {error}\n")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("materialize")
    command.add_argument("--selected-runs", type=Path, required=True)
    command.add_argument("--output-root", type=Path, required=True)
    command.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    try:
        materialize(read_receipt(args.selected_runs), args.output_root, args.receipt)
    except (
        OSError,
        ValueError,
        zipfile.BadZipFile,
        subprocess.SubprocessError,
    ) as error:
        parser.exit(1, f"payload materialization failed: {error}\n")


if __name__ == "__main__":
    main()
