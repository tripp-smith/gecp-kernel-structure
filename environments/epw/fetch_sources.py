"""Fetch and verify the exact QE sources and native MgB2 smoke inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Any


def _download(item: dict[str, Any], output: Path) -> None:
    digest = hashlib.sha256()
    size = 0
    with (
        urllib.request.urlopen(item["url"], timeout=120) as response,
        output.open("wb") as stream,
    ):
        while block := response.read(1024 * 1024):
            digest.update(block)
            size += len(block)
            stream.write(block)
    if digest.hexdigest() != item["sha256"] or size != item["bytes"]:
        raise ValueError(f"source identity mismatch: {item['name']}")


def fetch_sources(root: Path, lock_path: Path) -> None:
    """Keep upstream source code unchanged and reject unpinned downloads."""
    lock = json.loads(lock_path.read_text())
    if lock["schema"] != "kernelgecp.epw-sources.v1":
        raise ValueError("unsupported EPW source lock")
    root.parent.mkdir(parents=True, exist_ok=True)
    for item in lock["sources"]:
        print(f"Fetching pinned {item['name']} {item['commit']}", flush=True)
        with tempfile.TemporaryDirectory(prefix="epw-source-") as temporary:
            temporary_path = Path(temporary)
            archive = temporary_path / "source.tar.gz"
            _download(item, archive)
            extracted = temporary_path / "extracted"
            extracted.mkdir()
            with tarfile.open(archive) as source:
                source.extractall(extracted, filter="data")
            top = list(extracted.iterdir())
            if len(top) != 1 or not top[0].is_dir():
                raise ValueError("source archive must contain one source directory")
            destination = root / item["destination"]
            if destination.exists():
                # codeload may include an empty gitlink directory.
                destination.rmdir()
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(top[0]), destination)
            if item["name"] != "qe":
                # QE's tarball build otherwise fetches submodules implicitly.
                # The archive above was already checked by SHA256 and byte size.
                subprocess.run(["git", "init", "--quiet", destination], check=True)
    for item in lock["pseudopotentials"]:
        destination = root / "pseudo" / item["name"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        _download(item, destination)
    version_source = (root / "EPW/src/epw.f90").read_text()
    if "version_number = '6.0'" not in version_source:
        raise ValueError("pinned source does not declare EPW 6.0")
    metadata = {
        "schema": lock["schema"],
        "lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "sources": lock["sources"],
        "pseudopotentials": lock["pseudopotentials"],
        "native_input_sha256": {
            name: hashlib.sha256(
                (root / lock["calculation"]["upstream_directory"] / name).read_bytes()
            ).hexdigest()
            for name in lock["calculation"]["inputs"]
        },
    }
    (root / "kernelgecp-source-manifest.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    arguments = parser.parse_args()
    fetch_sources(arguments.root, arguments.lock)


if __name__ == "__main__":
    main()
