"""Build exact Green sources inside the isolated Linux provider image."""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path


def run(arguments: list[str], *, cwd: Path | None = None) -> None:
    subprocess.run(arguments, cwd=cwd, check=True)


def main() -> None:
    root = Path("/opt/green-src")
    root.mkdir()
    lock = json.loads(Path("/opt/green-lock/sources.json").read_text())
    source_overrides = []
    for source in lock["sources"]:
        name, revision = source["name"], source["revision"]
        if not re.fullmatch(r"green-[a-z0-9-]+", name) or not re.fullmatch(
            r"[0-9a-f]{40}", revision
        ):
            raise ValueError("invalid source lock entry")
        checkout = root / name
        run(["git", "init", str(checkout)])
        run(["git", "remote", "add", "origin", source["repository"]], cwd=checkout)
        run(["git", "fetch", "--depth=1", "origin", revision], cwd=checkout)
        run(["git", "checkout", "--detach", "FETCH_HEAD"], cwd=checkout)
        actual = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=checkout, text=True
        ).strip()
        if actual != revision:
            raise RuntimeError(f"source identity mismatch for {name}")
        source_overrides.append(f"-DFETCHCONTENT_SOURCE_DIR_{name.upper()}={checkout}")
    build = Path("/opt/green-build")
    run(
        [
            "cmake",
            "-S",
            str(root / "green-mbpt"),
            "-B",
            str(build),
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DGREEN_RELEASE=v1.0.0",
            "-DBuild_Tests=OFF",
            "-DFETCHCONTENT_FULLY_DISCONNECTED=ON",
            "-DFETCHCONTENT_UPDATES_DISCONNECTED=ON",
            *source_overrides,
        ]
    )
    run(["cmake", "--build", str(build), "--target", "mbpt.exe", "--parallel", "2"])
    destination = Path("/opt/green/bin")
    destination.mkdir(parents=True)
    run(["cp", str(build / "mbpt.exe"), str(destination / "mbpt.exe")])
    provenance = {
        "sources": lock,
        "native_executable_sha256": subprocess.check_output(
            ["sha256sum", str(destination / "mbpt.exe")], text=True
        ).split()[0],
        "system_packages": subprocess.check_output(
            ["dpkg-query", "-W", "-f=${Package}=${Version}\n"], text=True
        ).splitlines(),
        "python_packages": subprocess.check_output(
            ["/opt/venv/bin/python", "-m", "pip", "freeze", "--all"], text=True
        ).splitlines(),
    }
    Path("/opt/green-lock/build-provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
