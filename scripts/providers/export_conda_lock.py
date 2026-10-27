#!/usr/bin/env python3
"""Export replayable, hash-pinned public conda-forge installed-package locks.

This reads conda-meta in the executing interpreter's prefix. It does not infer
provider execution from installed packages; the numerical receipt is separate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

SCHEMA = "kernelgecp.conda-runtime-lock.v1"


class CondaLockError(ValueError):
    """Installed metadata cannot establish a safe replayable exact lock."""


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise CondaLockError("duplicate metadata key")
        result[key] = value
    return result


def _constant(value: str) -> Any:
    raise CondaLockError("nonfinite metadata number")


def _float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise CondaLockError("nonfinite metadata number")
    return result


def _metadata(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 16 * 1024**2:
        raise CondaLockError("metadata must be a bounded regular JSON file")
    try:
        result = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_object,
            parse_constant=_constant,
            parse_float=_float,
        )
    except CondaLockError:
        raise
    except (UnicodeError, ValueError) as error:
        raise CondaLockError("invalid metadata JSON") from error
    if not isinstance(result, dict):
        raise CondaLockError("metadata must be an object")
    return result


def _identifier(record: dict[str, Any], key: str) -> str:
    value = record.get(key)
    if (
        not isinstance(value, str)
        or len(value) > 256
        or re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.+!*-]*", value) is None
    ):
        raise CondaLockError(f"missing or invalid package {key}")
    return value


def _public_text(value: Any, field: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 4096
        or any(ord(character) < 32 for character in value)
        or "/" in value
        or "\\" in value
        or any(word in value.lower() for word in ("token=", "password=", "api_key="))
    ):
        raise CondaLockError(f"missing or unsafe package {field}")
    return value


def _package(record: dict[str, Any]) -> dict[str, Any]:
    package: dict[str, Any] = {
        key: _identifier(record, key) for key in ("name", "version", "build", "subdir")
    }
    digest = record.get("sha256")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise CondaLockError(
            "missing or invalid artifact SHA256; no unverified fallback"
        )
    url = record.get("url")
    if (
        not isinstance(url, str)
        or len(url) > 4096
        or any(ord(character) < 33 for character in url)
    ):
        raise CondaLockError("missing artifact URL")
    try:
        parsed = urlsplit(url)
    except ValueError as error:
        raise CondaLockError("invalid artifact URL") from error
    if (
        parsed.scheme != "https"
        or parsed.netloc != "conda.anaconda.org"
        or parsed.query
        or parsed.fragment
        or url != parsed.geturl()
    ):
        raise CondaLockError("only credential-free public conda-forge URLs are allowed")
    stem = "{name}-{version}-{build}".format(**package)
    allowed = {
        f"/conda-forge/{package['subdir']}/{stem}{suffix}"
        for suffix in (".conda", ".tar.bz2")
    }
    if parsed.path not in allowed:
        raise CondaLockError("artifact URL does not match the public package identity")
    filename = parsed.path.rsplit("/", 1)[1]
    if "fn" in record and record["fn"] != filename:
        raise CondaLockError("artifact filename disagrees with installed metadata")
    dependencies = record.get("depends")
    if not isinstance(dependencies, list):
        raise CondaLockError("missing package dependency list")
    package.update(
        {
            "url": url,
            "sha256": digest,
            "depends": sorted(
                _public_text(item, "dependency") for item in dependencies
            ),
            "license": _public_text(record.get("license"), "license"),
        }
    )
    return package


def build_lock(prefix: Path) -> tuple[str, dict[str, Any]]:
    """Return exact URL#SHA256 lines and a deterministic installed inventory."""
    directory = prefix / "conda-meta"
    if not directory.is_dir() or directory.is_symlink():
        raise CondaLockError("executing prefix has no regular conda-meta directory")
    paths = sorted(directory.glob("*.json"))
    if not paths:
        raise CondaLockError("conda-meta contains no package records")
    packages = [_package(_metadata(path)) for path in paths]
    packages.sort(key=lambda package: package["name"])
    if len({package["name"] for package in packages}) != len(packages):
        raise CondaLockError("duplicate installed package name")
    platforms = {package["subdir"] for package in packages} - {"noarch"}
    if len(platforms) != 1:
        raise CondaLockError(
            "one native platform plus optional noarch packages is required"
        )
    platform = platforms.pop()
    if re.fullmatch(r"(?:linux|osx|win)-[A-Za-z0-9_]+", platform) is None:
        raise CondaLockError("invalid native conda platform")
    explicit = "@EXPLICIT\n" + "".join(
        f"{package['url']}#{package['sha256']}\n" for package in packages
    )
    inventory = {
        "schema": SCHEMA,
        "provenance_scope": "installed_conda_metadata_not_numerical_execution",
        "platform": platform,
        "package_count": len(packages),
        "explicit_sha256": hashlib.sha256(explicit.encode()).hexdigest(),
        "packages": packages,
    }
    return explicit, inventory


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, default=Path(sys.prefix))
    parser.add_argument("--explicit", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.explicit.resolve() == args.inventory.resolve():
            raise CondaLockError("explicit lock and inventory must use different files")
        if args.explicit.is_symlink() or args.inventory.is_symlink():
            raise CondaLockError("output files must not be symlinks")
        explicit, inventory = build_lock(args.prefix)
        rendered = (
            json.dumps(inventory, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        for path, payload in ((args.explicit, explicit), (args.inventory, rendered)):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
    except CondaLockError as error:
        parser.exit(1, f"conda lock export rejected: {error}\n")
    except OSError as error:
        parser.exit(1, f"conda lock export failed: {type(error).__name__}\n")
    print(
        json.dumps(
            {
                key: inventory[key]
                for key in ("platform", "package_count", "explicit_sha256")
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
