"""Replay frozen spectroscopy auditors over hash-pinned evidence; never sample."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
import stat
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
FROZEN = "feaac669eb49d700f0059f5be755a9e4403a1719"
RECORDED_WORKSPACE = Path(
    "/home/runner/work/gecp-kernel-structure/gecp-kernel-structure"
)
ARCHIVES = {
    "spectroscopy-linux-feaac66-passed.tar.gz": (
        "63ac446677a9c36dfb17736e9823582a92ce49c59d85c3c0f393fd35c4ee2576"
    ),
    "cthyb-linux-7200ec0-passed.zip": (
        "5fba88be6a6796d29947b5e3bc5bf360e155813dec18d15e506ebe1f51832e77"
    ),
    "cthyb_9d28053_campaign.tar.gz": (
        "c50e2db4464001e6bc30a948e14a0da9176d3d005d3ae3c7cb43d13c3267bfd3"
    ),
}
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def validate_name(name):
    require(
        isinstance(name, str)
        and name
        and not name.startswith("/")
        and "\\" not in name
        and all(x not in ("", ".", "..") for x in name.rstrip("/").split("/")),
        "unsafe archive path",
    )
    return PurePosixPath(name.rstrip("/"))


def extract_tar(archive, destination, *, count, total, prefix=None, apple=False):
    """Validate the entire inventory before extracting, including skipped metadata."""
    with tarfile.open(archive, "r:gz") as z:
        entries = z.getmembers()
        require(len(entries) <= 4000, "too many archive members")
        seen, files = set(), []
        for entry in entries:
            path = validate_name(entry.name)
            key = str(path)
            require(
                key not in seen and (entry.isfile() or entry.isdir()),
                "duplicate or nonregular archive member",
            )
            seen.add(key)
            require(entry.size <= 16_000_000, "archive member too large")
            if not entry.isfile():
                continue
            if apple and any(part.startswith("._") for part in path.parts):
                continue
            if prefix:
                require(
                    path.parts[0] == prefix and len(path.parts) > 1,
                    "unexpected archive prefix",
                )
            files.append((entry, path))
        require(
            len(files) == count and sum(e.size for e, _ in files) == total,
            "archive size or count mismatch",
        )
        for entry, path in files:
            target = destination.joinpath(*path.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            source = z.extractfile(entry)
            require(source is not None, "missing archive data")
            with source, target.open("xb") as out:
                while chunk := source.read(65536):
                    out.write(chunk)


def extract_zip(archive, destination):
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        require(
            len(entries) == 467 and sum(e.file_size for e in entries) == 21621432,
            "parent inventory mismatch",
        )
        seen = set()
        for e in entries:
            path = validate_name(e.filename)
            mode = e.external_attr >> 16
            require(
                str(path) not in seen
                and not e.is_dir()
                and not stat.S_ISLNK(mode)
                and stat.S_IFMT(mode) in (0, stat.S_IFREG)
                and e.file_size <= 1_000_000,
                "unsafe ZIP member",
            )
            seen.add(str(path))
        for e in entries:
            target = destination.joinpath(*PurePosixPath(e.filename).parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(e) as source, target.open("xb") as out:
                while chunk := source.read(65536):
                    out.write(chunk)


def runtime_versions():
    """Use each provider's actual version interface, including Conda namespaces."""
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir")
    }
    for name in ("triqs", "triqs_cthyb", "triqs_maxent"):
        versions[name] = str(importlib.import_module(name + ".version").version)
    return versions


def require_runtime(evidence):
    require(sys.platform.startswith("linux"), "Linux required for scientific replay")
    require(
        all(os.environ.get(x) == "1" for x in THREADS), "thread environment mismatch"
    )
    for name in (
        "triqs",
        "triqs_cthyb",
        "triqs_maxent",
        "pydlr",
        "sparse_ir",
        "numpy",
        "scipy",
    ):
        importlib.import_module(name)
    from triqs_maxent import TauMaxEnt  # noqa: F401

    snapshot = json.loads(
        (evidence / "native-calibration/calibration.json").read_text()
    )["snapshot"]
    require(
        sys.version == snapshot["python"], "Python runtime differs from calibration"
    )
    versions = runtime_versions()
    require(
        set(versions) == set(snapshot["packages"]), "runtime package inventory mismatch"
    )
    for name, expected in snapshot["packages"].items():
        require(versions[name] == expected, "provider version mismatch: " + name)
    for name, record in snapshot["provider_modules"].items():
        module = importlib.import_module(name)
        require(
            digest(Path(module.__file__)) == record["sha256"],
            "provider module mismatch: " + name,
        )


def verified_archives(root):
    archives = {name: root / "experiments/data" / name for name in ARCHIVES}
    for name, path in archives.items():
        require(
            path.is_file()
            and not path.is_symlink()
            and path.stat().st_size <= 64 * 2**20,
            "missing, unsafe or oversized archive: " + name,
        )
        require(digest(path) == ARCHIVES[name], "archive hash mismatch: " + name)
    return archives


def validate_workspace(root):
    """Frozen worker receipts bind absolute paths in the original Linux runner."""
    require(
        root == root.resolve() and root == RECORDED_WORKSPACE,
        "scientific replay requires the recorded Linux CI workspace: "
        + str(RECORDED_WORKSPACE),
    )


def prepare_evidence(root):
    validate_workspace(root)
    evidence = root / "artifacts"
    require(not evidence.is_symlink(), "recorded evidence workspace is a symlink")
    evidence.mkdir(exist_ok=True)
    require(not any(evidence.iterdir()), "recorded evidence workspace must be empty")
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        require(
            sys.platform.startswith("linux"), "Linux required for scientific replay"
        )
        archives = verified_archives(ROOT)
        candidate = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
        with tempfile.TemporaryDirectory(prefix="gecp-eh-replay-") as temporary:
            base = Path(temporary)
            # Immutable worker receipts contain absolute original CI paths.
            # Restore that layout; never rewrite path-bearing evidence or auditors.
            evidence = prepare_evidence(ROOT)
            extract_tar(
                archives["spectroscopy-linux-feaac66-passed.tar.gz"],
                evidence,
                count=1451,
                total=60859249,
            )
            extract_zip(archives["cthyb-linux-7200ec0-passed.zip"], evidence / "parent")
            extract_tar(
                archives["cthyb_9d28053_campaign.tar.gz"],
                evidence / "cthyb-input",
                count=423,
                total=20361064,
                prefix="gecp-ee-odd-9d28053",
                apple=True,
            )
            require_runtime(evidence)
            frozen = base / "frozen"
            subprocess.run(
                ["git", "clone", "--shared", "--no-checkout", str(ROOT), str(frozen)],
                check=True,
                stdout=subprocess.DEVNULL,
            )
            subprocess.run(
                ["git", "checkout", "--detach", FROZEN],
                cwd=frozen,
                check=True,
                stdout=subprocess.DEVNULL,
            )
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts/check_spectroscopy_release.py"),
                    str(evidence),
                    "--frozen-checkout",
                    str(frozen),
                    "--candidate",
                    candidate,
                    "--input-repository",
                    str(ROOT),
                    "--output",
                    str(output),
                ],
                check=True,
            )
            result = json.loads(output.read_text())
            require(result.get("status") == "passed", "independent replay failed")
            result["archive_sha256"] = dict(ARCHIVES)
            result["checker_sha256"] = digest(
                ROOT / "scripts/check_spectroscopy_release.py"
            )
            result["replay_script_sha256"] = digest(Path(__file__))
            output.write_text(
                json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
            )
        return 0
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        ImportError,
        tarfile.TarError,
        zipfile.BadZipFile,
        subprocess.SubprocessError,
    ) as exc:
        failure = {}
        if output.is_file():
            try:
                value = json.loads(output.read_text())
                if isinstance(value, dict):
                    failure = value
            except (OSError, ValueError):
                pass
        failure["status"] = "failed"
        failure.setdefault("error", str(exc))
        failure["wrapper_error"] = str(exc)
        output.write_text(json.dumps(failure, indent=2, allow_nan=False) + "\n")
        print("spectroscopy replay failed: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
