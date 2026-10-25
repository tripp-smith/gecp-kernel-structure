"""Prepare and inspect complete, Git-bound native corresponding-source assets.

No package is published. Cargo vendoring is an explicit prior step. The build
command prints a JSON receipt; --report additionally saves that receipt.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import subprocess
import tarfile
import tomllib
from pathlib import Path, PurePosixPath
from typing import Any

MAX_FILES = 20_000
MAX_TOTAL_BYTES = 256 * 2**20
MAX_FILE_BYTES = 16 * 2**20
MAX_ARCHIVE_BYTES = 128 * 2**20
_MANIFEST = "SOURCE_MANIFEST.json"
_LICENSE_BUNDLE = "NATIVE_THIRD_PARTY_LICENSES.txt"
_SBOM = "NATIVE_DEPENDENCIES.cyclonedx.json"
_CONFIG = (
    b'[source.crates-io]\nreplace-with = "vendored-sources"\n'
    b'[source.vendored-sources]\ndirectory = "vendor"\n'
)
_REQUIRED = {
    "pyproject.toml",
    "uv.lock",
    "LICENSE",
    "README.md",
    "APPLICATION.md",
    "rust-toolchain.toml",
    "scripts/prepare_corresponding_source.py",
    "scripts/offline_native_build.sh",
    "native/Cargo.toml",
    "native/Cargo.lock",
    "native/build.rs",
    "native/src/lib.rs",
    "native/README.md",
    "native/SOURCE_BUILD.md",
    "native/THIRD_PARTY_NOTICES.md",
    "native/licenses/LGPL-3.0.txt",
    "native/licenses/GPL-3.0.txt",
    "kernelgecp/__init__.py",
    "kernelgecp/py.typed",
}
_RELEASE_DOCUMENTS = {
    "PACKAGE_README.md",
    "METHODS.md",
    "CITATION.cff",
    "RELEASE_NOTES.md",
}
_SCHEMA = "kernelgecp.corresponding-source.v1"
_RECEIPT_SCHEMA = "kernelgecp.corresponding-source.receipt.v1"


class SourceAssetError(ValueError):
    """A source asset is incomplete, unbound, unsafe or inconsistent."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceAssetError(message)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _hash(value: object, length: int = 64) -> bool:
    return (
        isinstance(value, str)
        and re.fullmatch(rf"[0-9a-f]{{{length}}}", value) is not None
    )


def _path(value: object) -> str:
    _require(isinstance(value, str) and 0 < len(value) <= 4096, "invalid source path")
    assert isinstance(value, str)
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and "\\" not in value
        and ":" not in value
        and all(part not in {"", ".", ".."} for part in value.split("/")),
        "unsafe source path",
    )
    return value


def _git(root: Path, *args: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, check=False
    )
    _require(result.returncode == 0, "Git source binding failed")
    return result.stdout


def _tree(root: Path, revision: str) -> tuple[str, dict[str, tuple[str, int]]]:
    _require(_hash(revision, 40), "revision must be a full lowercase Git commit SHA")
    actual = (
        _git(root, "rev-parse", "--verify", f"{revision}^{{commit}}").decode().strip()
    )
    _require(actual == revision, "source revision is not the requested commit")
    tree = _git(root, "rev-parse", f"{revision}^{{tree}}").decode().strip()
    selected: dict[str, tuple[str, int]] = {}
    for raw in _git(root, "ls-tree", "-r", "-z", revision).split(b"\0"):
        if not raw:
            continue
        header, raw_name = raw.split(b"\t", 1)
        mode, kind, blob = header.decode().split()
        name = raw_name.decode("utf-8")
        include = name in _REQUIRED or name in _RELEASE_DOCUMENTS
        if name.startswith("kernelgecp/"):
            include = name.endswith((".py", ".pyi")) or name.endswith("/py.typed")
        if name.startswith("native/"):
            include = True
            _require(
                not set(PurePosixPath(name).parts)
                & {"target", ".git", ".venv", "__pycache__", "dist"},
                "generated native source in Git tree",
            )
            _require(
                PurePosixPath(name).suffix
                in {".rs", ".toml", ".lock", ".md", ".txt", ".sh", ".json"},
                "unexpected native source type",
            )
        if name.startswith("scripts/") and name.endswith((".py", ".sh")):
            include = True
        if include:
            _path(name)
            _require(
                kind == "blob" and mode in {"100644", "100755"},
                "source must be a regular Git file",
            )
            selected[name] = (blob, 0o755 if mode == "100755" else 0o644)
    _require(selected.keys() >= _REQUIRED, "required source file missing from revision")
    return tree, selected


def _version(payloads: dict[str, bytes]) -> str:
    project = tomllib.loads(payloads["pyproject.toml"].decode())["project"]
    native = tomllib.loads(payloads["native/Cargo.toml"].decode())["package"]
    _require(
        isinstance(project, dict) and isinstance(native, dict), "invalid project tables"
    )
    version = project.get("version")
    _require(
        project.get("name") == "kernelgecp"
        and native.get("name") == "kernelgecp-native",
        "unexpected project identity",
    )
    _require(
        isinstance(version, str)
        and re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(?:[A-Za-z0-9.+-]*)?", version)
        is not None,
        "invalid project version",
    )
    _require(native.get("version") == version, "native/project version mismatch")
    if "readme" in project:
        readme = project["readme"]
        _require(
            isinstance(readme, str) and readme in payloads,
            "configured package readme missing from corresponding source",
        )
    return str(version)


def _normalize_license(value: object) -> str:
    _require(
        isinstance(value, str) and 0 < len(value) <= 512,
        "crate license declaration missing",
    )
    expression = " ".join(str(value).split())
    # az 1.3.0 retains Cargo's legacy slash spelling; its bundled README and
    # LICENSE-APACHE/LICENSE-MIT explicitly state the alternative licenses.
    if expression == "MIT/Apache-2.0":
        expression = "MIT OR Apache-2.0"
    for old, new in {
        "LGPL-3.0+": "LGPL-3.0-or-later",
        "LGPL-2.1+": "LGPL-2.1-or-later",
        "GPL-3.0+": "GPL-3.0-or-later",
        "GPL-2.0+": "GPL-2.0-or-later",
        "AGPL-3.0+": "AGPL-3.0-or-later",
    }.items():
        expression = expression.replace(old, new)
    tokens = re.findall(r"\(|\)|[A-Za-z0-9][A-Za-z0-9.-]*", expression)
    _require(
        "".join(tokens) == re.sub(r"\s", "", expression), "invalid SPDX license syntax"
    )
    position = 0

    def atom() -> None:
        nonlocal position
        _require(position < len(tokens), "incomplete SPDX expression")
        if tokens[position] == "(":
            position += 1
            expr()
            _require(
                position < len(tokens) and tokens[position] == ")",
                "unbalanced SPDX expression",
            )
            position += 1
        else:
            _require(
                tokens[position] not in {"AND", "OR", "WITH", ")"},
                "invalid SPDX operand",
            )
            position += 1
            if position < len(tokens) and tokens[position] == "WITH":
                position += 1
                _require(
                    position < len(tokens)
                    and tokens[position] not in {"AND", "OR", "WITH", "(", ")"},
                    "invalid SPDX exception",
                )
                position += 1

    def expr() -> None:
        nonlocal position
        atom()
        while position < len(tokens) and tokens[position] in {"AND", "OR"}:
            position += 1
            atom()

    expr()
    _require(position == len(tokens), "invalid SPDX expression tail")
    return expression


def _lock_packages(payload: bytes) -> dict[tuple[str, str], dict[str, Any]]:
    lock = tomllib.loads(payload.decode())
    result = {}
    for item in lock.get("package", []):
        if "source" not in item:
            continue
        _require(
            str(item["source"]).startswith("registry+"),
            "only locked registry dependencies are supported",
        )
        _require(
            isinstance(item.get("name"), str)
            and isinstance(item.get("version"), str)
            and _hash(item.get("checksum")),
            "invalid locked crate identity",
        )
        key = item["name"], item["version"]
        _require(key not in result, "duplicate locked crate identity")
        result[key] = item
    _require(
        {name for name, _ in result} >= {"rug", "gmp-mpfr-sys", "pyo3"},
        "native dependency lock is incomplete",
    )
    return result


def _vendor_metadata(
    payloads: dict[str, bytes], lock_bytes: bytes
) -> list[dict[str, Any]]:
    expected = _lock_packages(lock_bytes)
    found: dict[tuple[str, str], dict[str, Any]] = {}
    directories = sorted(
        {
            PurePosixPath(name).parts[1]
            for name in payloads
            if name.startswith("vendor/") and len(PurePosixPath(name).parts) > 2
        }
    )
    for directory in directories:
        prefix = f"vendor/{directory}/"
        _require(
            prefix + "Cargo.toml" in payloads
            and prefix + ".cargo-checksum.json" in payloads,
            "vendored crate metadata missing",
        )
        package = tomllib.loads(payloads[prefix + "Cargo.toml"].decode())["package"]
        _require(isinstance(package, dict), "vendored crate package table malformed")
        key = package.get("name"), package.get("version")
        _require(
            key in expected and key not in found,
            "unexpected or duplicate vendored crate",
        )
        checksums = json.loads(payloads[prefix + ".cargo-checksum.json"])
        _require(
            isinstance(checksums, dict) and set(checksums) == {"files", "package"},
            "vendored checksum envelope malformed",
        )
        _require(
            checksums.get("package") == expected[key]["checksum"]
            and isinstance(checksums.get("files"), dict),
            "vendor package checksum mismatch",
        )
        actual_files = {
            name[len(prefix) :] for name in payloads if name.startswith(prefix)
        } - {".cargo-checksum.json"}
        _require(
            actual_files == set(checksums["files"]), "vendor file coverage mismatch"
        )
        for relative, checksum in checksums["files"].items():
            _path(relative)
            _require(
                _hash(checksum) and _sha(payloads[prefix + relative]) == checksum,
                "vendor file checksum mismatch",
            )
        licenses = sorted(
            prefix + name
            for name in actual_files
            if PurePosixPath(name)
            .name.upper()
            .startswith(("LICENSE", "COPYING", "NOTICE", "UNLICENSE"))
        )
        _require(
            bool(licenses) and all(payloads[name].strip() for name in licenses),
            "vendor license coverage missing",
        )
        found[key] = dict(
            name=key[0],
            version=key[1],
            registry_checksum=expected[key]["checksum"],
            license_declared=package["license"],
            license=_normalize_license(package.get("license")),
            license_files=licenses,
        )
    _require(found.keys() == expected.keys(), "locked vendor dependencies missing")
    return [found[key] for key in sorted(found)]


def _generated(
    packages: list[dict[str, Any]], payloads: dict[str, bytes], version: str
) -> tuple[bytes, bytes, str]:
    expressions = sorted({"Apache-2.0", *(package["license"] for package in packages)})
    combined = " AND ".join(f"({expression})" for expression in expressions)
    lines = [f"kernelgecp {version}: bundled native dependency notices\n\n"]
    seen = set()
    for package in packages:
        for name in package["license_files"]:
            data = payloads[name]
            checksum = _sha(data)
            lines.append(
                f"\n--- {package['name']}@{package['version']}: {name} "
                f"(SHA256 {checksum}) ---\n"
            )
            if checksum in seen:
                lines.append(
                    "Identical license text is included above under the same SHA256.\n"
                )
            else:
                lines.append(data.decode("utf-8") + "\n")
                seen.add(checksum)
    licenses = "".join(lines).encode("utf-8")
    _require(len(licenses) <= MAX_FILE_BYTES, "license bundle exceeds file budget")
    sbom = dict(
        bomFormat="CycloneDX",
        specVersion="1.5",
        version=1,
        metadata=dict(
            component=dict(type="library", name="kernelgecp", version=version)
        ),
        components=[
            dict(
                type="library",
                name=p["name"],
                version=p["version"],
                purl=f"pkg:cargo/{p['name']}@{p['version']}",
                hashes=[dict(alg="SHA-256", content=p["registry_checksum"])],
                licenses=[dict(expression=p["license"])],
            )
            for p in packages
        ],
    )
    encoded = (json.dumps(sbom, indent=2, sort_keys=True) + "\n").encode()
    _require(
        not any(
            marker in encoded
            for marker in (b"/Users/", b"/private/", b"/tmp/", b"file://", b"C:\\")
        ),
        "host path in dependency inventory",
    )
    return licenses, encoded, combined


def _limits(payloads: dict[str, bytes]) -> None:
    _require(
        len(payloads) <= MAX_FILES
        and sum(map(len, payloads.values())) <= MAX_TOTAL_BYTES,
        "source asset exceeds resource budget",
    )
    _require(
        all(len(data) <= MAX_FILE_BYTES for data in payloads.values()),
        "source file exceeds size budget",
    )
    for name in payloads:
        _path(name)


def _receipt(manifest: dict[str, Any], container: bytes) -> dict[str, Any]:
    version = manifest["version"]
    return dict(
        schema=_RECEIPT_SCHEMA,
        version=version,
        source_revision=manifest["source_revision"],
        source_tree=manifest["source_tree"],
        archive=f"kernelgecp-{version}-native-source.tar.gz",
        sha256=_sha(container),
        archive_bytes=len(container),
        license_bundle=f"kernelgecp-{version}-native-licenses.txt",
        license_bundle_sha256=manifest["license_bundle_sha256"],
        sbom=f"kernelgecp-{version}-native-dependencies.cyclonedx.json",
        sbom_sha256=manifest["sbom_sha256"],
        combined_native_license_expression=manifest[
            "combined_native_license_expression"
        ],
        vendored_crates=len(manifest["packages"]),
        payload_files=len(manifest["files"]),
    )


def prepare_corresponding_source(
    *, project_root: Path, revision: str, vendor: Path, output: Path
) -> dict[str, Any]:
    """Build an immutable-source asset and its wheel-staging license/inventory files."""
    try:
        root, output = project_root.resolve(), output.resolve()
        _require(
            not vendor.is_symlink() and vendor.is_dir(),
            "vendor must be a real directory",
        )
        tree, selected = _tree(root, revision)
        payloads: dict[str, bytes] = {}
        modes: dict[str, int] = {}
        blobs = []
        for name, (blob, mode) in sorted(selected.items()):
            source = root / name
            _require(
                source.resolve().is_relative_to(root)
                and source.is_file()
                and source.stat().st_size <= MAX_FILE_BYTES,
                "selected source exceeds file budget or escapes root",
            )
            data = _git(root, "cat-file", "blob", blob)
            _require(
                source.is_file()
                and not source.is_symlink()
                and source.read_bytes() == data,
                "selected source is dirty or missing",
            )
            _require(
                name.endswith("/py.typed") or bool(data.strip()),
                "selected source is empty",
            )
            payloads[name], modes[name] = data, mode
            blobs.append(dict(path=name, git_blob=blob, sha256=_sha(data), mode=mode))
        version = _version(payloads)
        vendor_total = sum(map(len, payloads.values()))
        for source in sorted(vendor.rglob("*")):
            _require(not source.is_symlink(), "vendor symbolic link rejected")
            if source.is_file():
                name = "vendor/" + _path(source.relative_to(vendor).as_posix())
                _require(
                    source.stat().st_size <= MAX_FILE_BYTES,
                    "vendor file exceeds size budget",
                )
                vendor_total += source.stat().st_size
                _require(
                    vendor_total <= MAX_TOTAL_BYTES and len(payloads) < MAX_FILES,
                    "source asset exceeds resource budget",
                )
                payloads[name] = source.read_bytes()
                modes[name] = 0o755 if source.stat().st_mode & 0o111 else 0o644
        packages = _vendor_metadata(payloads, payloads["native/Cargo.lock"])
        licenses, sbom, combined = _generated(packages, payloads, version)
        payloads[_LICENSE_BUNDLE], payloads[_SBOM], payloads[".cargo/config.toml"] = (
            licenses,
            sbom,
            _CONFIG,
        )
        modes.update(
            {_LICENSE_BUNDLE: 0o644, _SBOM: 0o644, ".cargo/config.toml": 0o644}
        )
        _limits(payloads)
        manifest = dict(
            schema=_SCHEMA,
            version=version,
            source_revision=revision,
            source_tree=tree,
            root=f"kernelgecp-{version}-native-source",
            native_precision_bits=128,
            native_protocol_version=1,
            source_blobs=blobs,
            lock_sha256=_sha(payloads["native/Cargo.lock"]),
            packages=packages,
            license_bundle_sha256=_sha(licenses),
            sbom_sha256=_sha(sbom),
            combined_native_license_expression=combined,
            files=[
                dict(path=name, sha256=_sha(data), size=len(data), mode=modes[name])
                for name, data in sorted(payloads.items())
            ],
        )
        payloads[_MANIFEST] = (
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        ).encode()
        modes[_MANIFEST] = 0o644
        _limits(payloads)
        output.mkdir(parents=True, exist_ok=True)
        names = [
            f"kernelgecp-{version}-native-source.tar.gz",
            f"kernelgecp-{version}-native-licenses.txt",
            f"kernelgecp-{version}-native-dependencies.cyclonedx.json",
        ]
        _require(
            not any((output / name).exists() for name in names),
            "refusing to replace an existing output artifact",
        )
        buffer = io.BytesIO()
        with (
            gzip.GzipFile(
                filename="", mode="wb", fileobj=buffer, mtime=0
            ) as compressed,
            tarfile.open(
                mode="w|", fileobj=compressed, format=tarfile.PAX_FORMAT
            ) as archive,
        ):
            for name, data in sorted(payloads.items()):
                info = tarfile.TarInfo(f"{manifest['root']}/{name}")
                info.size, info.mode = len(data), modes[name]
                info.mtime = info.uid = info.gid = 0
                info.uname = info.gname = ""
                archive.addfile(info, io.BytesIO(data))
        container = buffer.getvalue()
        _require(
            len(container) <= MAX_ARCHIVE_BYTES,
            "compressed source asset exceeds budget",
        )
        for name, data in zip(names, (container, licenses, sbom), strict=True):
            with (output / name).open("xb") as target:
                target.write(data)
        receipt = inspect_corresponding_source(
            output / names[0],
            expected_sha256=_sha(container),
            project_root=root,
            revision=revision,
        )
        return receipt
    except SourceAssetError:
        raise
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        UnicodeError,
        tarfile.TarError,
    ) as error:
        raise SourceAssetError(
            f"could not prepare corresponding sources: {error}"
        ) from error


def inspect_corresponding_source(
    archive: Path,
    *,
    expected_sha256: str,
    project_root: Path | None = None,
    revision: str | None = None,
) -> dict[str, Any]:
    """Inspect without extraction/execution, including exact Git/blob binding."""
    try:
        _require(_hash(expected_sha256), "expected source SHA256 is invalid")
        _require(
            archive.is_file() and archive.stat().st_size <= MAX_ARCHIVE_BYTES,
            "source archive exceeds size budget",
        )
        container = archive.read_bytes()
        _require(_sha(container) == expected_sha256, "source archive SHA256 mismatch")
        payloads: dict[str, bytes] = {}
        modes: dict[str, int] = {}
        roots = set()
        total = 0
        with tarfile.open(fileobj=io.BytesIO(container), mode="r:gz") as source:
            for item in source:
                name = _path(item.name)
                parts = PurePosixPath(name).parts
                _require(
                    len(parts) >= 2
                    and item.isfile()
                    and not item.issym()
                    and not item.islnk(),
                    "source archive has a nonregular member",
                )
                roots.add(parts[0])
                relative = "/".join(parts[1:])
                _require(
                    relative not in payloads and 0 <= item.size <= MAX_FILE_BYTES,
                    "duplicate or oversized source member",
                )
                total += item.size
                _require(
                    total <= MAX_TOTAL_BYTES and len(payloads) < MAX_FILES,
                    "source archive exceeds expanded budget",
                )
                _require(
                    item.mode in {0o644, 0o755}
                    and item.uid == item.gid == item.mtime == 0
                    and item.uname == item.gname == "",
                    "source archive is not normalized",
                )
                stream = source.extractfile(item)
                _require(stream is not None, "source member unreadable")
                assert stream is not None
                data = stream.read(MAX_FILE_BYTES + 1)
                _require(len(data) == item.size, "truncated source member")
                payloads[relative], modes[relative] = data, item.mode
        _require(
            _MANIFEST in payloads and payloads.keys() >= _REQUIRED,
            "source manifest or required source missing",
        )
        manifest = json.loads(payloads[_MANIFEST])
        _require(
            isinstance(manifest, dict) and manifest.get("schema") == _SCHEMA,
            "source manifest schema mismatch",
        )
        _require(
            set(manifest)
            == {
                "schema",
                "version",
                "source_revision",
                "source_tree",
                "root",
                "native_precision_bits",
                "native_protocol_version",
                "source_blobs",
                "lock_sha256",
                "packages",
                "license_bundle_sha256",
                "sbom_sha256",
                "combined_native_license_expression",
                "files",
            },
            "source manifest envelope mismatch",
        )
        version = _version(payloads)
        _require(
            manifest.get("version") == version
            and roots == {f"kernelgecp-{version}-native-source"}
            and manifest.get("root") in roots,
            "source version/root mismatch",
        )
        _require(
            _hash(manifest.get("source_revision"), 40)
            and _hash(manifest.get("source_tree"), 40),
            "source Git provenance malformed",
        )
        _require(
            manifest.get("native_precision_bits") == 128
            and type(manifest.get("native_precision_bits")) is int
            and manifest.get("native_protocol_version") == 1
            and type(manifest.get("native_protocol_version")) is int,
            "native source contract mismatch",
        )
        files = manifest.get("files")
        _require(
            isinstance(files, list) and len(files) <= MAX_FILES,
            "source file inventory missing",
        )
        indexed = {}
        for item in files:
            _require(
                isinstance(item, dict)
                and set(item) == {"path", "sha256", "size", "mode"},
                "malformed source file record",
            )
            name = _path(item.get("path"))
            _require(
                name not in indexed and name != _MANIFEST and name in payloads,
                "source file inventory coverage mismatch",
            )
            _require(
                type(item.get("size")) is int
                and item["size"] == len(payloads[name])
                and _hash(item.get("sha256"))
                and item["sha256"] == _sha(payloads[name])
                and item.get("mode") == modes[name],
                "source payload hash/size/mode mismatch",
            )
            indexed[name] = item
        _require(
            indexed.keys() == payloads.keys() - {_MANIFEST},
            "source payload inventory coverage mismatch",
        )
        _require(
            manifest.get("lock_sha256") == _sha(payloads["native/Cargo.lock"]),
            "source lock hash mismatch",
        )
        packages = _vendor_metadata(payloads, payloads["native/Cargo.lock"])
        _require(
            packages == manifest.get("packages"),
            "source package/license inventory mismatch",
        )
        licenses, sbom, combined = _generated(packages, payloads, version)
        _require(
            payloads.get(_LICENSE_BUNDLE) == licenses
            and payloads.get(_SBOM) == sbom
            and payloads.get(".cargo/config.toml") == _CONFIG,
            "generated source license/SBOM/config mismatch",
        )
        _require(
            manifest.get("license_bundle_sha256") == _sha(licenses)
            and manifest.get("sbom_sha256") == _sha(sbom)
            and manifest.get("combined_native_license_expression") == combined,
            "generated source receipt mismatch",
        )
        blobs = manifest.get("source_blobs")
        _require(isinstance(blobs, list), "source Git blobs missing")
        bound = {}
        for entry in blobs:
            _require(
                isinstance(entry, dict)
                and set(entry) == {"path", "git_blob", "sha256", "mode"},
                "source Git blob malformed",
            )
            name = _path(entry.get("path"))
            _require(
                name not in bound
                and name in payloads
                and not name.startswith("vendor/"),
                "source Git blob coverage mismatch",
            )
            data = payloads[name]
            blob_hash = hashlib.sha1(
                b"blob " + str(len(data)).encode() + b"\0" + data
            ).hexdigest()
            _require(
                entry.get("git_blob") == blob_hash
                and entry.get("sha256") == _sha(data)
                and entry.get("mode") == modes[name],
                "source Git blob hash mismatch",
            )
            bound[name] = (blob_hash, modes[name])
        _require(bound.keys() >= _REQUIRED, "required Git source binding missing")
        _require(
            set(bound)
            | {name for name in payloads if name.startswith("vendor/")}
            | {_MANIFEST, _LICENSE_BUNDLE, _SBOM, ".cargo/config.toml"}
            == payloads.keys(),
            "unbound source payload",
        )
        if revision is not None:
            _require(
                revision == manifest["source_revision"], "source revision mismatch"
            )
        if project_root is not None:
            tree, selected = _tree(project_root.resolve(), manifest["source_revision"])
            _require(
                tree == manifest["source_tree"] and selected == bound,
                "source Git tree/binding mismatch",
            )
            for name, (blob, _) in selected.items():
                _require(
                    _git(project_root.resolve(), "cat-file", "blob", blob)
                    == payloads[name],
                    "source bytes disagree with Git",
                )
        return _receipt(manifest, container)
    except SourceAssetError:
        raise
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        UnicodeError,
        tarfile.TarError,
    ) as error:
        raise SourceAssetError(
            f"could not inspect corresponding sources: {error}"
        ) from error


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build")
    build.add_argument("--project-root", type=Path, required=True)
    build.add_argument("--revision", required=True)
    build.add_argument("--vendor", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--report", type=Path)
    inspect = commands.add_parser("inspect")
    inspect.add_argument("--archive", type=Path, required=True)
    inspect.add_argument("--expected-sha256", required=True)
    inspect.add_argument("--project-root", type=Path)
    inspect.add_argument("--revision")
    args = parser.parse_args()
    if args.command == "build":
        report = prepare_corresponding_source(
            project_root=args.project_root,
            revision=args.revision,
            vendor=args.vendor,
            output=args.output,
        )
        if args.report is not None:
            with args.report.open("x", encoding="utf-8") as target:
                json.dump(report, target, indent=2, sort_keys=True)
                target.write("\n")
    else:
        report = inspect_corresponding_source(
            args.archive,
            expected_sha256=args.expected_sha256,
            project_root=args.project_root,
            revision=args.revision,
        )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
