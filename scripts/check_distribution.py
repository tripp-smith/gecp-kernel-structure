"""Inspect release archives and exercise an installed wheel outside its checkout.

This gate does not publish anything.  Archive inspection uses only the standard
library; the installed check needs the package's runtime dependencies and mypy.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import re
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import zipfile
from email.message import Message
from email.parser import BytesParser
from pathlib import Path, PurePosixPath
from typing import Any

SDIST_ROOT_FILES = {
    ".gitignore",
    "APPLICATION.md",
    "LICENSE",
    "PKG-INFO",
    "README.md",
    "pyproject.toml",
}
FORBIDDEN_PARTS = {
    ".agents",
    ".cache",
    ".git",
    ".github",
    ".lake",
    ".venv",
    "__pycache__",
    "build",
    "data",
    "dist",
    "experiments",
    "target",
    "tests",
}
FORBIDDEN_NAMES = {
    ".env",
    ".pypirc",
    "FINDINGS 2.md",
    "MATHLIB 2.md",
    "autonomous-implementation.md",
    "id_rsa",
    "id_ed25519",
}


class DistributionError(ValueError):
    """The distribution is incomplete or violates the release boundary."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise DistributionError(message)


def _safe_name(name: str) -> PurePosixPath:
    path = PurePosixPath(name)
    _require(
        bool(name)
        and not path.is_absolute()
        and "\\" not in name
        and not any(part in {"", ".", ".."} for part in name.rstrip("/").split("/"))
        and not any(part in FORBIDDEN_PARTS for part in path.parts)
        and path.name not in FORBIDDEN_NAMES,
        f"unsafe or excluded archive member: {name!r}",
    )
    return path


def _is_package_source(name: str) -> bool:
    path = PurePosixPath(name)
    return (
        path.parts[0] == "kernelgecp"
        and len(path.parts) >= 2
        and (path.suffix in {".py", ".pyi"} or path.name == "py.typed")
    )


def _metadata(
    payload: bytes,
    project: dict[str, Any],
    location: str,
    license_expression: str = "Apache-2.0",
) -> Message:
    metadata = BytesParser().parsebytes(payload)
    for field in ("Metadata-Version", "Name", "Version", "Requires-Python"):
        _require(
            len(metadata.get_all(field, [])) == 1,
            f"{location}: require exactly one {field}",
        )

    def normalize(name: object) -> str:
        return re.sub(r"[-_.]+", "-", str(name)).lower()

    _require(
        normalize(metadata["Name"]) == normalize(project["name"]),
        f"{location}: distribution name disagrees with pyproject.toml",
    )
    _require(
        metadata["Version"] == project["version"],
        f"{location}: version disagrees with pyproject.toml",
    )
    _require(
        metadata["Requires-Python"] == project["requires-python"],
        f"{location}: Python requirement disagrees with pyproject.toml",
    )
    _require(
        metadata.get("License-Expression") == license_expression
        or (
            license_expression == "Apache-2.0"
            and metadata.get("License", "").strip() == license_expression
        ),
        f"{location}: required {license_expression} license metadata is missing",
    )
    _require(
        not (metadata.get("License-Expression") and metadata.get("License")),
        f"{location}: mutually exclusive license fields",
    )
    return metadata


def _native_source(name: str) -> bool:
    path = PurePosixPath(name)
    return (
        name == "rust-toolchain.toml"
        or name
        in {
            "native/Cargo.toml",
            "native/Cargo.lock",
            "native/build.rs",
            "native/README.md",
            "native/SOURCE_BUILD.md",
            "native/THIRD_PARTY_NOTICES.md",
        }
        or (path.parts[:2] == ("native", "src") and path.suffix == ".rs")
        or (path.parts[:2] == ("native", "licenses") and path.suffix == ".txt")
    )


def _json_unique(payload: bytes) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in items:
            _require(key not in result, "duplicate inventory JSON key")
            result[key] = value
        return result

    result = json.loads(payload, object_pairs_hook=pairs)
    _require(isinstance(result, dict), "inventory must be a JSON object")
    return result


def _no_host_paths(value: Any) -> None:
    if isinstance(value, str):
        _require(
            "file://" not in value
            and re.search(r"/(?:Users|home|tmp|private)/", value) is None,
            "native metadata contains a host path",
        )
    elif isinstance(value, dict):
        for item in value.values():
            _no_host_paths(item)
    elif isinstance(value, list):
        for item in value:
            _no_host_paths(item)


def _native_inventory(payload: bytes, project: dict[str, Any]) -> dict[str, Any]:
    inventory = _json_unique(payload)
    _require(
        inventory.get("schema") == "kernelgecp.native-wheel.v1"
        and inventory.get("version") == project["version"],
        "native inventory schema/version mismatch",
    )
    for key in ("source_revision", "source_tree"):
        _require(
            isinstance(inventory.get(key), str)
            and re.fullmatch(r"[0-9a-f]{40}", inventory[key]) is not None,
            "invalid native source identity",
        )
    source = inventory.get("source_asset")
    _require(isinstance(source, dict), "missing native corresponding source")
    filename = f"kernelgecp-{project['version']}-native-source.tar.gz"
    _require(
        source.get("filename") == filename
        and source.get("url")
        == (
            "https://github.com/tripp-smith/gecp-kernel-structure/releases/download/"
            f"v{project['version']}/{filename}"
        ),
        "native source URL/name mismatch",
    )
    for digest in (
        source.get("sha256"),
        inventory.get("license_bundle_sha256"),
        inventory.get("sbom_sha256"),
    ):
        _require(
            isinstance(digest, str)
            and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
            "invalid native artifact hash",
        )
    expression = inventory.get("license_expression")
    _require(
        isinstance(expression, str)
        and "Apache-2.0" in expression
        and "LGPL-3.0" in expression,
        "native composite license is missing",
    )
    _no_host_paths(inventory)
    return inventory


def check_archives(
    *, project_root: Path, wheel: Path, sdist: Path, kind: str = "pure"
) -> dict[str, Any]:
    """Validate both archive payloads against the actual package source tree."""
    project_root = project_root.resolve()
    project = tomllib.loads(
        (project_root / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    _require(project["name"] == "kernelgecp", "unexpected distribution project")
    source_files = {
        path.relative_to(project_root).as_posix(): path.read_bytes()
        for path in (project_root / "kernelgecp").rglob("*")
        if path.is_file()
        and _is_package_source(path.relative_to(project_root).as_posix())
    }
    _require("kernelgecp/__init__.py" in source_files, "missing package __init__.py")
    _require("kernelgecp/py.typed" in source_files, "missing source typing marker")
    license_bytes = (project_root / "LICENSE").read_bytes()
    _require(bool(license_bytes.strip()), "empty source license")
    _require(kind in {"pure", "native"}, "unknown wheel kind")
    native_sources = {
        path.relative_to(project_root).as_posix(): path.read_bytes()
        for path in (project_root / "native").rglob("*")
        if path.is_file() and _native_source(path.relative_to(project_root).as_posix())
    }
    if (project_root / "rust-toolchain.toml").is_file():
        native_sources["rust-toolchain.toml"] = (
            project_root / "rust-toolchain.toml"
        ).read_bytes()

    platform = ""
    if kind == "pure":
        _require(
            wheel.name.endswith("-py3-none-any.whl"), "wheel must be universal Python 3"
        )
    else:
        match = re.fullmatch(
            rf"kernelgecp-{re.escape(project['version'])}-cp312-abi3-([A-Za-z0-9_.]+)[.]whl",
            wheel.name,
        )
        _require(match is not None, "native wheel must use cp312-abi3 platform tags")
        assert match is not None
        platform = match.group(1)
        _require(
            all(
                re.fullmatch(
                    r"(?:manylinux[_0-9]+|macosx_[0-9]+_[0-9]+)_(?:x86_64|aarch64|arm64|universal2)",
                    p,
                )
                for p in platform.split(".")
            ),
            "unsupported native platform tag",
        )
    with zipfile.ZipFile(wheel) as archive:
        members = archive.infolist()
        names = [member.filename for member in members]
        _require(len(names) == len(set(names)), "duplicate wheel archive member")
        for member in members:
            _safe_name(member.filename)
            mode = member.external_attr >> 16
            _require(mode & 0o170000 != 0o120000, "wheel contains a symbolic link")
        wheel_files = {
            member.filename: archive.read(member)
            for member in members
            if not member.is_dir()
        }

    metadata_names = [
        name for name in wheel_files if name.endswith(".dist-info/METADATA")
    ]
    _require(len(metadata_names) == 1, "wheel must contain exactly one METADATA")
    dist_info = metadata_names[0].split("/", 1)[0]
    _require(
        dist_info == f"kernelgecp-{project['version']}.dist-info",
        "wheel dist-info directory disagrees with project version",
    )
    dist_info_files = {
        f"{dist_info}/{name}"
        for name in (
            "METADATA",
            "WHEEL",
            "RECORD",
            "entry_points.txt",
            "top_level.txt",
            "licenses/LICENSE",
        )
    }
    inventory: dict[str, Any] | None = None
    if kind == "native":
        inventory_name = "kernelgecp/_native_inventory.json"
        _require(inventory_name in wheel_files, "native wheel missing inventory")
        inventory = _native_inventory(wheel_files[inventory_name], project)
        dist_info_files.update(
            {
                "kernelgecp/_native.abi3.so",
                inventory_name,
                f"{dist_info}/licenses/native/_wheel_licenses.txt",
                f"{dist_info}/sboms/_wheel_dependencies.cyclonedx.json",
            }
        )
        for name, hash_key in (
            (
                f"{dist_info}/licenses/native/_wheel_licenses.txt",
                "license_bundle_sha256",
            ),
            (f"{dist_info}/sboms/_wheel_dependencies.cyclonedx.json", "sbom_sha256"),
        ):
            _require(name in wheel_files, "native wheel missing notices or SBOM")
            _require(
                hashlib.sha256(wheel_files[name]).hexdigest() == inventory[hash_key],
                "native license/SBOM hash mismatch",
            )
        sbom = _json_unique(
            wheel_files[f"{dist_info}/sboms/_wheel_dependencies.cyclonedx.json"]
        )
        _no_host_paths(sbom)
        binary = wheel_files.get("kernelgecp/_native.abi3.so", b"")
        _require(
            len(binary) > 16
            and binary[:4] in {b"\x7fELF", b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe"},
            "native extension is missing or malformed",
        )
        _require(
            not any(
                path in binary
                for path in (b"/Users/", b"/home/", b"/private/var/", b"/tmp/")
            ),
            "native extension contains a host build path",
        )
    for name in wheel_files:
        _require(
            _is_package_source(name) or name in dist_info_files,
            f"unexpected wheel payload: {name}",
        )
    parsed_metadata = _metadata(
        wheel_files[metadata_names[0]],
        project,
        "wheel METADATA",
        inventory["license_expression"] if inventory is not None else "Apache-2.0",
    )
    if inventory is not None:
        _require(
            set(parsed_metadata.get_all("License-File", []))
            == {"LICENSE", "native/_wheel_licenses.txt"},
            "native license-file metadata is incomplete",
        )
    for required in ("WHEEL", "RECORD", "licenses/LICENSE"):
        _require(f"{dist_info}/{required}" in wheel_files, f"wheel missing {required}")
    _require(
        wheel_files[f"{dist_info}/licenses/LICENSE"] == license_bytes,
        "wheel license differs from source license",
    )
    wheel_metadata = BytesParser().parsebytes(wheel_files[f"{dist_info}/WHEEL"])
    tags = set(wheel_metadata.get_all("Tag", []))
    expected_tags = {f"cp312-abi3-{p}" for p in platform.split(".")}
    _require(
        wheel_metadata.get("Wheel-Version") == "1.0",
        "unsupported pure universal/native wheel metadata",
    )
    _require(
        (
            kind == "pure"
            and wheel_metadata.get("Root-Is-Purelib", "").lower() == "true"
            and tags == {"py3-none-any"}
        )
        or (
            kind == "native"
            and wheel_metadata.get("Root-Is-Purelib", "").lower() == "false"
            and tags in (expected_tags, {f"cp312-abi3-{platform}"})
        ),
        "wheel metadata does not match its pure/native platform tags",
    )
    packaged_sources = {
        name: payload
        for name, payload in wheel_files.items()
        if _is_package_source(name)
    }
    _require(
        packaged_sources == source_files, "wheel package sources differ from checkout"
    )

    _require(sdist.name.endswith(".tar.gz"), "sdist must be a .tar.gz archive")
    with tarfile.open(sdist, "r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        _require(len(names) == len(set(names)), "duplicate sdist archive member")
        roots: set[str] = set()
        sdist_files: dict[str, bytes] = {}
        for member in members:
            path = _safe_name(member.name)
            roots.add(path.parts[0])
            _require(
                member.isdir() or member.isfile(), "sdist contains a nonregular member"
            )
            if member.isfile():
                _require(
                    len(path.parts) > 1, "sdist file must have a project-root prefix"
                )
                relative = PurePosixPath(*path.parts[1:]).as_posix()
                _require(
                    relative in SDIST_ROOT_FILES
                    or _is_package_source(relative)
                    or _native_source(relative),
                    f"unexpected sdist payload: {relative}",
                )
                stream = archive.extractfile(member)
                _require(stream is not None, f"unreadable sdist member: {member.name}")
                assert stream is not None
                sdist_files[relative] = stream.read()
    _require(len(roots) == 1, "sdist must contain exactly one project-root directory")
    _require(
        roots == {f"kernelgecp-{project['version']}"},
        "sdist root directory disagrees with project version",
    )
    for required in ("pyproject.toml", "APPLICATION.md", "LICENSE", "PKG-INFO"):
        _require(required in sdist_files, f"sdist missing {required}")
    for name in (
        "pyproject.toml",
        "APPLICATION.md",
        "LICENSE",
        "README.md",
        ".gitignore",
    ):
        if name in sdist_files:
            _require(
                sdist_files[name] == (project_root / name).read_bytes(),
                f"sdist {name} differs from checkout",
            )
    _metadata(sdist_files["PKG-INFO"], project, "sdist PKG-INFO")
    _require(
        {name: data for name, data in sdist_files.items() if _is_package_source(name)}
        == source_files,
        "sdist package sources differ from checkout",
    )
    _require(
        {name: data for name, data in sdist_files.items() if _native_source(name)}
        == native_sources,
        "sdist native sources differ from checkout",
    )
    return {
        "name": project["name"],
        "version": project["version"],
        "kind": kind,
        "package_source_files": len(source_files),
        "archives": [
            {
                "filename": path.name,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (wheel, sdist)
        ],
    }


TYPING_CONSUMER = """\
import numpy as np
import numpy.typing as npt
from typing import assert_type
from kernelgecp import (
    FermionicBasis, FermionicFit, LeastSquaresPlan, FitPrecisionCertificate,
    GridPrecisionCertificate, certify_fit_precision, certify_grid_pivot_precision,
)

design: npt.NDArray[np.float64] = np.eye(3)
values: npt.NDArray[np.float64] = np.ones(3)
plan = LeastSquaresPlan(design)
assert_type(plan.fit(values), npt.NDArray[np.float64] | npt.NDArray[np.complex128])
basis = FermionicBasis(
    beta=2.0, wmax=3.0, tol=1e-8, reference_backend="python",
)
assert_type(
    FermionicBasis.load("type-only.npz", reference_backend="rust"), FermionicBasis,
)
assert_type(
    FermionicFit.load("type-only.npz", reference_backend="rust"), FermionicFit,
)
fit = basis.fit_tau(values, tau=np.linspace(0.0, 2.0, 3), axis=0)
tau_values = fit.evaluate_tau(np.array([0.2, 0.8]))
iw_values = fit.evaluate_matsubara(np.array([-1, 0, 2]))
precision = certify_fit_precision(
    fit, values, np.linspace(0.0, 2.0, 3), query_tau=[0.2],
)
assert_type(precision, FitPrecisionCertificate)
grid = certify_grid_pivot_precision(basis, [0.2, 0.8], [-1.0, 1.0])
assert_type(grid, GridPrecisionCertificate)
"""


def check_installed(*, project_root: Path, backend: str = "python") -> dict[str, Any]:
    """Check actual installed imports, physical conventions and distributed typing."""
    project_root = project_root.resolve()
    _require(backend in {"python", "rust"}, "unknown installed backend")
    _require(
        not Path.cwd().resolve().is_relative_to(project_root),
        "installed checks must run outside the checkout",
    )
    package = importlib.import_module("kernelgecp")
    origin = Path(package.__file__).resolve()
    _require(
        not origin.is_relative_to(project_root), "import resolved to checkout source"
    )
    _require(
        origin.is_relative_to(Path(sys.prefix).resolve()),
        "import is outside the active environment",
    )
    distribution = importlib.metadata.distribution("kernelgecp")
    _require(
        distribution.locate_file("kernelgecp/py.typed").is_file(),
        "installed package lacks py.typed",
    )
    project = tomllib.loads(
        (project_root / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]
    _require(
        distribution.version == package.__version__ == project["version"],
        "installed, runtime and project versions disagree",
    )

    import numpy as np

    basis = package.FermionicBasis(
        beta=2.0, wmax=3.0, tol=1e-8, precision="reference", reference_backend=backend
    )
    tau = np.linspace(0.0, 2.0, 65)
    pole = 0.7
    values = -np.exp(-tau * pole) / (1.0 + np.exp(-2.0 * pole))
    fitted = basis.fit_tau(values, tau=tau, axis=0)
    heldout = (np.arange(33) + 0.5) * 2.0 / 33
    expected_tau = -np.exp(-heldout * pole) / (1.0 + np.exp(-2.0 * pole))
    actual_tau = np.asarray(fitted.evaluate_tau(heldout))
    indices = np.array([-3, -1, 0, 2])
    expected_iw = 1.0 / (1j * (2 * indices + 1) * np.pi / 2.0 - pole)
    actual_iw = np.asarray(fitted.evaluate_matsubara(indices))
    for actual, expected, domain in (
        (actual_tau, expected_tau, "tau"),
        (actual_iw, expected_iw, "Matsubara"),
    ):
        _require(actual.shape == expected.shape, f"{domain} smoke output shape changed")
        _require(
            bool(np.all(np.isfinite(actual)))
            and float(np.max(np.abs(actual - expected))) < 1e-5,
            f"{domain} physical one-pole smoke failed",
        )
    certificate = package.certify_fit_precision(
        fitted,
        values,
        tau,
        query_tau=heldout[[0, 15, 32]],
        matsubara_indices=indices,
    )
    _require(
        certificate.certified
        and certificate.time is not None
        and certificate.matsubara is not None
        and certificate.time.max_error < 1e-5
        and certificate.matsubara.max_error < 1e-5,
        "installed finite-query coefficient/solve certificate failed",
    )
    grid_certificate = package.certify_grid_pivot_precision(
        basis, heldout[[0, 15, 32]], [-2.3, -0.3, 1.7]
    )
    _require(
        grid_certificate.certified
        and grid_certificate.eta is not None
        and grid_certificate.eta > 0,
        "installed finite-grid pivot certificate failed",
    )
    _require(
        basis.diagnostics.requested_reference_backend == backend
        and "numpy" in basis.diagnostics.fit_backend,
        "installed backend diagnostics disagree with execution policy",
    )
    if backend == "rust":
        inventory = _native_inventory(
            (origin.parent / "_native_inventory.json").read_bytes(), project
        )
        native = importlib.import_module("kernelgecp._native")
        _require(
            Path(native.__file__).resolve().parent == origin.parent,
            "native module is outside the installed package",
        )
        info = native.build_info()
        _require(
            info.get("package_version") == distribution.version
            and info.get("source_revision") == inventory["source_revision"]
            and info.get("working_precision_bits") == 128
            and info.get("protocol_version") == 1
            and info.get("mpfr_tls") is True,
            "native build provenance or arithmetic protocol mismatch",
        )
        _require(
            basis.diagnostics.time_evaluation_backend == "rust_mpfr_128"
            and basis.diagnostics.frequency_evaluation_backend
            == "python_mpmath_reference",
            "native component scope is not time-only 128-bit evaluation",
        )
        with tempfile.TemporaryDirectory(
            prefix="kernelgecp-installed-archive-"
        ) as saved:
            path = Path(saved) / "fit.npz"
            fitted.save(path)
            portable = package.FermionicFit.load(path)
            native_restored = package.FermionicFit.load(path, reference_backend="rust")
            _require(
                portable.basis.diagnostics.requested_reference_backend == "python"
                and native_restored.basis.diagnostics.requested_reference_backend
                == "rust",
                "archive restoration activated an unexpected backend",
            )
            for restored in (portable, native_restored):
                _require(
                    bool(
                        np.allclose(
                            restored.evaluate_tau(heldout),
                            actual_tau,
                            rtol=1e-12,
                            atol=1e-12,
                        )
                    ),
                    "installed native/portable archive round trip changed predictions",
                )
    design = np.array([[1.0, 0.0], [1.0, 1.0], [1.0, 2.0]])
    expected_coefficients = np.array([[2.0 + 1j, 0.5 - 2j], [1j, -3.0]])
    fitted_coefficients = package.LeastSquaresPlan(design).fit(
        expected_coefficients @ design.T, axis=1
    )
    _require(
        np.asarray(fitted_coefficients).shape == expected_coefficients.shape
        and bool(
            np.allclose(
                fitted_coefficients, expected_coefficients, rtol=1e-12, atol=1e-12
            )
        ),
        "installed least-squares complex/batch/axis smoke failed",
    )
    with tempfile.TemporaryDirectory(prefix="kernelgecp-typing-") as temporary:
        consumer = Path(temporary) / "consumer.py"
        typing_consumer = TYPING_CONSUMER
        if backend == "rust":
            typing_consumer += """
from kernelgecp import _native
assert_type(_native.build_info(), dict[str, str | int | bool])
def native_consumer(handle: _native.ReferenceCardinal128) -> None:
    assert_type(handle.rank, int)
    assert_type(handle.design_tau([0.25, 0.75]), bytes)
"""
        consumer.write_text(typing_consumer, encoding="utf-8")
        mypy = [sys.executable, "-I", "-m", "mypy", "--strict", "--no-incremental"]
        subprocess.run(
            [*mypy, str(consumer)],
            cwd=temporary,
            check=True,
        )
        consumer.write_text(
            typing_consumer + '\nplan.fit(values, axis="invalid-axis")\n',
            encoding="utf-8",
        )
        rejected = subprocess.run(
            [*mypy, str(consumer)],
            cwd=temporary,
            capture_output=True,
            text=True,
            check=False,
        )
        _require(
            rejected.returncode != 0 and "[arg-type]" in rejected.stdout,
            "installed typing consumer did not reject an invalid public argument",
        )
    return {
        "name": "kernelgecp",
        "version": distribution.version,
        "imported_from": str(origin),
        "tau_smoke_error": float(np.max(np.abs(actual_tau - expected_tau))),
        "matsubara_smoke_error": float(np.max(np.abs(actual_iw - expected_iw))),
        "typing_consumer": "passed",
        "typing_rejection_consumer": "passed",
        "scope": "installed_api_smoke_not_a_uniform_accuracy_certificate",
        "requested_backend": backend,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    inspect = subparsers.add_parser("inspect", help="validate wheel and sdist contents")
    inspect.add_argument("--project-root", type=Path, required=True)
    inspect.add_argument("--wheel", type=Path, required=True)
    inspect.add_argument("--sdist", type=Path, required=True)
    inspect.add_argument("--kind", choices=("pure", "native"), default="pure")
    installed = subparsers.add_parser(
        "installed", help="check an installed wheel and typing"
    )
    installed.add_argument("--project-root", type=Path, required=True)
    installed.add_argument("--backend", choices=("python", "rust"), default="python")
    args = parser.parse_args()
    try:
        if args.command == "inspect":
            report = check_archives(
                project_root=args.project_root,
                wheel=args.wheel,
                sdist=args.sdist,
                kind=args.kind,
            )
        else:
            report = check_installed(
                project_root=args.project_root, backend=args.backend
            )
    except (OSError, ValueError, zipfile.BadZipFile, tarfile.TarError) as exc:
        parser.exit(1, f"distribution check failed: {exc}\n")
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
