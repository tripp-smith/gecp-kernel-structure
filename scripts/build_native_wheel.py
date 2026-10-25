"""Build an optional native wheel from an inspected corresponding-source asset.

Hatchling keeps producing the universal wheel. This command uses an isolated
copy of the exact source asset and gives the native distribution its own license
metadata, complete notices and sanitized dependency inventory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from pathlib import Path
from typing import Any


class NativeWheelError(ValueError):
    """A source asset or native build does not match the declared release."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise NativeWheelError(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def read_source_report(path: Path) -> dict[str, Any]:
    _require(path.stat().st_size <= 1024**2, "source report exceeds size budget")
    report = json.loads(path.read_text())
    _require(isinstance(report, dict), "source report must be an object")
    _require(
        report.get("schema") == "kernelgecp.corresponding-source.receipt.v1",
        "unsupported source report schema",
    )
    for name in ("source_revision", "source_tree"):
        _require(
            isinstance(report.get(name), str)
            and re.fullmatch(r"[0-9a-f]{40}", report[name]) is not None,
            f"invalid source {name}",
        )
    for name, hash_name in (
        ("archive", "sha256"),
        ("license_bundle", "license_bundle_sha256"),
        ("sbom", "sbom_sha256"),
    ):
        filename, digest = report.get(name), report.get(hash_name)
        _require(
            isinstance(filename, str)
            and Path(filename).name == filename
            and filename not in ("", ".", "..")
            and "/" not in filename
            and "\\" not in filename,
            f"{name} must be a basename",
        )
        _require(
            isinstance(digest, str)
            and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            f"invalid {hash_name}",
        )
        artifact = path.parent / filename
        _require(artifact.is_file() and not artifact.is_symlink(), f"missing {name}")
        limit = 256 * 1024**2 if name == "archive" else 16 * 1024**2
        _require(artifact.stat().st_size <= limit, f"{name} exceeds size budget")
        _require(_sha256(artifact) == digest, f"{name} hash mismatch")
    expression = report.get("combined_native_license_expression")
    _require(
        isinstance(expression, str)
        and len(expression) <= 16_384
        and re.fullmatch(r"[A-Za-z0-9 .()+\-]+", expression) is not None
        and "Apache-2.0" in expression
        and "LGPL-3.0" in expression,
        "missing or invalid combined native license expression",
    )
    return report


def overlay_project_metadata(text: str, expression: str) -> str:
    lines = text.splitlines()
    license_line = 'license = { text = "Apache-2.0" }'
    _require(lines.count(license_line) == 1, "unexpected base project license metadata")
    index = lines.index(license_line)
    lines[index : index + 1] = [
        "license = " + json.dumps(expression),
        'license-files = ["LICENSE", "native/_wheel_licenses.txt"]',
    ]
    section = "[tool.maturin.sbom]"
    _require(lines.count(section) == 1, "missing controlled Maturin SBOM configuration")
    index = lines.index(section) + 1
    lines.insert(index, 'include = ["native/_wheel_dependencies.cyclonedx.json"]')
    result = "\n".join(lines) + "\n"
    parsed = tomllib.loads(result)
    _require(parsed["project"]["license"] == expression, "license overlay mismatch")
    _require(
        parsed["tool"]["maturin"]["sbom"].get("rust") is False
        and parsed["tool"]["maturin"]["sbom"].get("auditwheel") is False,
        "automatic host-specific SBOM generation must be disabled",
    )
    return result


def native_build_environment(
    working: Path,
    source_revision: str,
    *,
    original: dict[str, str] | None = None,
    platform: str | None = None,
) -> dict[str, str]:
    """Bind source provenance and remove compiler/linker host paths."""
    environment = dict(os.environ if original is None else original)
    canonical = working.resolve()
    environment["KERNELGECP_SOURCE_REVISION"] = source_revision
    # A cached C library may have been compiled without path remapping.
    environment["GMP_MPFR_SYS_CACHE"] = "_"
    flags = (
        environment.get("CARGO_ENCODED_RUSTFLAGS", "").split("\x1f")
        if environment.get("CARGO_ENCODED_RUSTFLAGS")
        else shlex.split(environment.get("RUSTFLAGS", ""))
    )
    flags.append(f"--remap-path-prefix={canonical}=/kernelgecp-source")
    if (sys.platform if platform is None else platform) == "darwin":
        # Rust's default Mach-O LC_ID_DYLIB embeds the absolute output path.
        # Compiler prefix maps do not rewrite linker load commands.
        flags += ["-C", "link-arg=-Wl,-install_name,@rpath/_native.abi3.so"]
    environment.pop("RUSTFLAGS", None)
    environment["CARGO_ENCODED_RUSTFLAGS"] = "\x1f".join(flags)
    for name in ("CFLAGS", "CXXFLAGS", "CPPFLAGS"):
        c_flags = shlex.split(environment.get(name, ""))
        c_flags.append(f"-ffile-prefix-map={canonical}=/kernelgecp-source")
        environment[name] = shlex.join(c_flags)
    return environment


def build_native_wheel(
    *, project_root: Path, source_report: Path, output: Path, compatibility: str | None
) -> dict[str, Any]:
    # Import here so report/metadata validation remains a small independent gate.
    if __package__:
        from .prepare_corresponding_source import inspect_corresponding_source
    else:
        from prepare_corresponding_source import inspect_corresponding_source

    root = project_root.resolve()
    report_path = source_report.resolve()
    report = read_source_report(report_path)
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    _require(project["name"] == "kernelgecp", "unexpected project")
    _require(
        project["version"] == report.get("version"), "source/project version mismatch"
    )
    archive = report_path.parent / report["archive"]
    inspected = inspect_corresponding_source(
        archive,
        expected_sha256=report["sha256"],
        project_root=root,
        revision=report["source_revision"],
    )
    for key in (
        "version",
        "source_revision",
        "source_tree",
        "sha256",
        "license_bundle_sha256",
        "sbom_sha256",
        "combined_native_license_expression",
    ):
        _require(
            inspected.get(key) == report.get(key), f"source report disagrees: {key}"
        )
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _require(not list(output.glob("*.whl")), "output already contains wheels")
    with tempfile.TemporaryDirectory(prefix="kernelgecp-native-wheel-") as temporary:
        stage = Path(temporary).resolve()
        with tarfile.open(archive, "r:gz") as source:
            source.extractall(stage, filter="data")
        roots = list(stage.iterdir())
        _require(len(roots) == 1 and roots[0].is_dir(), "source archive root mismatch")
        working = roots[0]
        pyproject = working / "pyproject.toml"
        pyproject.write_text(
            overlay_project_metadata(
                pyproject.read_text(), report["combined_native_license_expression"]
            )
        )
        for source_key, destination in (
            ("license_bundle", "native/_wheel_licenses.txt"),
            ("sbom", "native/_wheel_dependencies.cyclonedx.json"),
        ):
            (working / destination).write_bytes(
                (report_path.parent / report[source_key]).read_bytes()
            )
        inventory = {
            "schema": "kernelgecp.native-wheel.v1",
            "version": report["version"],
            "source_revision": report["source_revision"],
            "source_tree": report["source_tree"],
            "source_asset": {
                "filename": report["archive"],
                "sha256": report["sha256"],
                "url": "https://github.com/tripp-smith/gecp-kernel-structure/releases/download/"
                f"v{report['version']}/{report['archive']}",
            },
            "license_bundle_sha256": report["license_bundle_sha256"],
            "sbom_sha256": report["sbom_sha256"],
            "license_expression": report["combined_native_license_expression"],
        }
        (working / "kernelgecp/_native_inventory.json").write_text(
            json.dumps(inventory, sort_keys=True, indent=2) + "\n"
        )
        environment = native_build_environment(working, report["source_revision"])
        command = [
            sys.executable,
            "-m",
            "maturin",
            "build",
            "--release",
            "--locked",
            "--offline",
            "--out",
            str(output),
        ]
        if compatibility is not None:
            command += ["--compatibility", compatibility]
        subprocess.run(command, cwd=working, env=environment, check=True)
    wheels = list(output.glob("*.whl"))
    _require(len(wheels) == 1, "native build must produce exactly one wheel")
    return {
        "wheel": wheels[0].name,
        "sha256": _sha256(wheels[0]),
        "source_revision": report["source_revision"],
        "source_asset_sha256": report["sha256"],
        "license_expression": report["combined_native_license_expression"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--source-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compatibility")
    args = parser.parse_args()
    result = build_native_wheel(
        project_root=args.project_root,
        source_report=args.source_report,
        output=args.output,
        compatibility=args.compatibility,
    )
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
