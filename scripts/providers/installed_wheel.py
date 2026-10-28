"""Bind imported source modules to an actual non-editable wheel distribution."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from importlib import metadata
from pathlib import Path
from typing import Any


def installed_wheel_identity(
    module_paths: dict[str, Path], package_files: dict[str, str]
) -> dict[str, Any]:
    try:
        distribution = metadata.distribution("kernelgecp")
    except metadata.PackageNotFoundError as error:
        raise RuntimeError(
            "installed kernelgecp distribution metadata is required"
        ) from error
    if set(module_paths) != set(package_files) or not package_files:
        raise RuntimeError("installed module identity is incomplete")
    files, record = distribution.files, distribution.read_text("RECORD")
    if files is None or not isinstance(record, str) or not record.strip():
        raise RuntimeError("an actual installed wheel RECORD is required")
    direct_url = distribution.read_text("direct_url.json")
    if direct_url is not None:
        try:
            direct = json.loads(direct_url)
        except (TypeError, ValueError) as error:
            raise RuntimeError("invalid wheel installation metadata") from error
        if not isinstance(direct, dict) or not isinstance(
            direct.get("dir_info", {}), dict
        ):
            raise RuntimeError("invalid wheel installation metadata")
        directory = direct.get("dir_info", {})
        if "editable" in directory and type(directory["editable"]) is not bool:
            raise RuntimeError("invalid editable installation metadata")
        if directory.get("editable") is True:
            raise RuntimeError(
                "editable installations are not installed-wheel evidence"
            )
    try:
        rows = list(csv.reader(io.StringIO(record), strict=True))
    except csv.Error as error:
        raise RuntimeError("invalid wheel RECORD") from error
    if any(len(row) != 3 for row in rows):
        raise RuntimeError("invalid wheel RECORD")
    listed = {str(item) for item in files}
    for name, relative in package_files.items():
        if relative not in listed or sum(row[0] == relative for row in rows) != 1:
            raise RuntimeError("imported module is absent from wheel file metadata")
        if (
            module_paths[name].resolve()
            != Path(distribution.locate_file(relative)).resolve()
        ):
            raise RuntimeError(
                "imported module does not belong to the installed distribution"
            )
    return {
        "distribution": "kernelgecp",
        "version": distribution.version,
        "record_sha256": hashlib.sha256(record.encode()).hexdigest(),
        "module_count": len(package_files),
        "editable": False,
    }
