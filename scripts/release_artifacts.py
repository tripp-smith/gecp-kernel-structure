"""Promote the exact tested main artifacts; never build during publication."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import time
import tomllib
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

if __package__:
    from .build_native_wheel import read_source_report
    from .check_distribution import check_archives
    from .prepare_corresponding_source import inspect_corresponding_source
else:
    from build_native_wheel import read_source_report
    from check_distribution import check_archives
    from prepare_corresponding_source import inspect_corresponding_source

REPOSITORY = "tripp-smith/gecp-kernel-structure"
ARTIFACTS = {
    "kernelgecp-release-source",
    "kernelgecp-release-native-linux",
    "kernelgecp-release-native-macos",
}


# Verifier-owned 1.2 delivery inventory. Presence is an asset gate only; it
# does not establish scientific eligibility or any provider-domain gate.
# Missing future guides deliberately remain requirements until authored.
REQUIRED_12_ASSETS = {
    "core": ("METHODS.md", "CITATION.cff", "RELEASE_NOTES.md", "FINAL_HANDOFF.md"),
    "domain-guides": (
        "INTERACTING_DMFT.md",
        "CTHYB_DMFT.md",
        "MOLECULAR_GREEN.md",
        "EPW_INTEROPERABILITY.md",
        "SPECTROSCOPY.md",
        "COVARIANCE.md",
    ),
    "user-guides": (
        "INSTALLATION_COMPATIBILITY.md",
        "TROUBLESHOOTING.md",
        "MIGRATION.md",
    ),
    "examples": (
        "examples/interacting_dmft.py",
        "examples/molecular_green.py",
        "examples/triqs_interoperability.py",
        "examples/epw_restart_archive.py",
        "examples/spectral_continuation.py",
        "examples/covariance_operator.py",
        "experiments/covariance_benchmark.py",
    ),
    "reports": (
        "data/cthyb/linux-audit-binding-7200ec0.json",
        "data/cthyb/linux-independent-audit-7200ec0.json",
        "data/epw/native-archive-linux-audit.json",
        "data/covariance/audit.json",
    ),
    "source-notices": (
        "native/SOURCE_BUILD.md",
        "native/THIRD_PARTY_NOTICES.md",
        "PROVIDER_NOTICES.md",
        "environments/provider-licenses.json",
        "environments/locks/triqs-linux-64-conda-inventory.json",
        "environments/locks/triqs-osx-arm64-conda-inventory.json",
        "environments/locks/dlr-physics-inventory.json",
        "environments/locks/green-linux-build.json",
    ),
}


def _is_12(version: str) -> bool:
    return tuple(int(part) for part in version.split(".")) >= (1, 2, 0)


def _required_12_records() -> dict[str, tuple[str, str]]:
    records = {}
    for group, paths in REQUIRED_12_ASSETS.items():
        for source in paths:
            name = Path(source).name
            _require(name not in records, "duplicate required 1.2 asset basename")
            records[name] = (source, group)
    return records


def _check_12_content(path: Path, group: str) -> None:
    _require(
        path.is_file() and not path.is_symlink(),
        f"missing required 1.2 {group} asset: {path.name}",
    )
    data = path.read_bytes()
    _require(
        bool(data.strip()),
        f"empty required 1.2 {group} asset: {path.name}",
    )
    try:
        data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ReleaseError(
            f"invalid text required 1.2 {group} asset: {path.name}"
        ) from exc
    if path.suffix == ".json":
        try:
            parsed = _json(path)
        except (ValueError, UnicodeDecodeError) as exc:
            raise ReleaseError(
                f"invalid required 1.2 {group} JSON: {path.name}"
            ) from exc
        _require(
            isinstance(parsed, dict) and bool(parsed),
            f"empty required 1.2 {group} JSON: {path.name}",
        )


class ReleaseError(ValueError):
    """A release does not satisfy the exact-source/artifact contract."""


def _require(value: bool, message: str) -> None:
    if not value:
        raise ReleaseError(message)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    _require(
        path.is_file() and not path.is_symlink() and path.stat().st_size <= 2**20,
        "invalid JSON receipt",
    )

    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    result = json.loads(path.read_text(), object_pairs_hook=unique)
    _require(isinstance(result, dict), "receipt must be an object")
    return result


def _git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, text=True, check=False
    )
    _require(result.returncode == 0, "Git release identity/ancestry check failed")
    return result.stdout.strip()


def _gh(path: str, *query: str) -> dict[str, Any]:
    result = subprocess.run(
        ["gh", "api", path, *query], capture_output=True, text=True, check=False
    )
    _require(result.returncode == 0, "GitHub release metadata request failed")
    value = json.loads(result.stdout)
    _require(isinstance(value, dict), "invalid GitHub metadata response")
    return value


def select_runs(
    repository: str,
    revision: str,
    run_items: dict[str, list[dict[str, Any]]],
    artifact_items: list[dict[str, Any]],
) -> dict[str, Any]:
    """Require successful main push producers; PR merge jobs cannot publish."""
    _require(repository == REPOSITORY, "unexpected release repository")
    _require(re.fullmatch(r"[0-9a-f]{40}", revision) is not None, "invalid source SHA")
    chosen = {}
    for role, workflow in (("ci", "ci.yml"), ("native", "native.yml")):
        eligible = [
            run
            for run in run_items.get(role, [])
            if run.get("head_sha") == revision
            and run.get("head_branch") == "main"
            and run.get("event") == "push"
            and run.get("path") == f".github/workflows/{workflow}"
            and run.get("repository", {}).get("full_name") == repository
            and run.get("head_repository", {}).get("full_name") == repository
            and type(run.get("id")) is int
            and run["id"] > 0
            and type(run.get("run_attempt")) is int
            and run["run_attempt"] > 0
        ]
        _require(bool(eligible), f"no successful exact-main {role} producer")
        run = max(eligible, key=lambda item: item["id"])
        _require(
            run.get("status") == "completed" and run.get("conclusion") == "success",
            f"latest exact-main {role} producer has not succeeded",
        )
        chosen[f"{role}_run"] = {
            "id": run["id"],
            "attempt": run["run_attempt"],
            "url": run["html_url"],
        }
    _require(len(artifact_items) == len(ARTIFACTS), "wrong producer artifact count")
    _require(
        {a.get("name") for a in artifact_items} == ARTIFACTS, "wrong artifact names"
    )
    for artifact in artifact_items:
        run = artifact.get("workflow_run", {})
        _require(
            artifact.get("expired") is False
            and type(artifact.get("id")) is int
            and artifact["id"] > 0
            and run.get("id") == chosen["native_run"]["id"]
            and run.get("head_sha") == revision
            and run.get("head_branch") == "main",
            "expired or mismatched producer artifact",
        )
    return {**chosen, "artifact_ids": {a["name"]: a["id"] for a in artifact_items}}


def select_release(
    root: Path, ref: str, event_sha: str, repository: str
) -> dict[str, Any]:
    match = re.fullmatch(
        r"refs/tags/v((?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*))", ref
    )
    _require(match is not None, "release requires a stable version tag")
    assert match is not None
    version = match[1]
    revision = _git(root, "rev-parse", f"{ref}^{{commit}}")
    _require(
        revision == event_sha == _git(root, "rev-parse", "HEAD"),
        "tag/event/HEAD mismatch",
    )
    _git(root, "merge-base", "--is-ancestor", revision, "origin/main")
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    native = tomllib.loads((root / "native/Cargo.toml").read_text())["package"]
    runtime = re.search(
        r'^__version__ = "([^"]+)"$',
        (root / "kernelgecp/__init__.py").read_text(),
        re.M,
    )
    _require(
        project["name"] == "kernelgecp"
        and native["name"] == "kernelgecp-native"
        and project["version"] == native["version"] == version
        and runtime is not None
        and runtime[1] == version,
        "tag/Python/Cargo/runtime version mismatch",
    )
    runs = {
        role: _gh(
            f"repos/{repository}/actions/workflows/{workflow}/runs",
            "--method",
            "GET",
            "-f",
            f"head_sha={revision}",
            "-f",
            "branch=main",
            "-f",
            "event=push",
            "-f",
            "per_page=100",
        )["workflow_runs"]
        for role, workflow in (("ci", "ci.yml"), ("native", "native.yml"))
    }
    # Select the same latest successful candidate whose artifacts are queried.
    native_candidates = [
        r
        for r in runs["native"]
        if r.get("status") == "completed" and r.get("conclusion") == "success"
    ]
    _require(bool(native_candidates), "native main producer is not complete")
    producer = max(native_candidates, key=lambda r: r["id"])
    artifacts = _gh(
        f"repos/{repository}/actions/runs/{producer['id']}/artifacts",
        "-f",
        "per_page=100",
        "--method",
        "GET",
    )["artifacts"]
    selected = select_runs(repository, revision, runs, artifacts)
    _require(
        selected["native_run"]["id"] == producer["id"], "producer selection changed"
    )
    return {
        "schema": "kernelgecp.selected-release.v1",
        "repository": repository,
        "version": version,
        "tag": f"v{version}",
        "source_revision": revision,
        "source_tree": _git(root, "rev-parse", "HEAD^{tree}"),
        **selected,
    }


def _files(root: Path) -> set[Path]:
    _require(
        root.is_dir() and not root.is_symlink(),
        "artifact directory missing or symlinked",
    )
    result = set()
    for path in root.rglob("*"):
        _require(not path.is_symlink(), "symlinked artifact payload")
        _require(path.is_dir() or path.is_file(), "nonregular artifact payload")
        if path.is_file():
            _require(path.stat().st_size <= 256 * 2**20, "artifact exceeds byte budget")
            result.add(path.relative_to(root))
    return result


def collect_artifacts(
    *, project_root: Path, selected: dict[str, Any], artifacts: Path, output: Path
) -> dict[str, Any]:
    root = project_root.resolve()
    _require(
        selected.get("schema") == "kernelgecp.selected-release.v1"
        and selected.get("repository") == REPOSITORY,
        "invalid selected release",
    )
    version, revision = selected["version"], selected["source_revision"]
    _files(artifacts)
    _require(
        _git(root, "rev-parse", "HEAD") == revision
        and _git(root, "rev-parse", "HEAD^{tree}") == selected["source_tree"],
        "collector source mismatch",
    )
    _require(
        {p.name for p in artifacts.iterdir()} == ARTIFACTS,
        "unexpected downloaded artifacts",
    )
    base = artifacts / "kernelgecp-release-source"
    report_path = base / "source/source-report.json"
    report = read_source_report(report_path)
    _require(
        report["version"] == version
        and report["source_revision"] == revision
        and report["source_tree"] == selected["source_tree"],
        "source report identity mismatch",
    )
    source_names = {
        "source-report.json",
        report["archive"],
        report["license_bundle"],
        report["sbom"],
    }
    pure_name, sdist_name = (
        f"kernelgecp-{version}-py3-none-any.whl",
        f"kernelgecp-{version}.tar.gz",
    )
    expected = {Path("source") / name for name in source_names} | {
        Path("dist") / pure_name,
        Path("dist") / sdist_name,
        Path("validation/sdist-wheel") / pure_name,
    }
    _require(
        _files(base) == expected, "missing or extra source/portable/validation payload"
    )
    inspected = inspect_corresponding_source(
        base / "source" / report["archive"],
        expected_sha256=report["sha256"],
        project_root=root,
        revision=revision,
    )
    _require(inspected == report, "source receipt differs from inspected asset")
    sdist = base / "dist" / sdist_name
    for pure in (
        base / "dist" / pure_name,
        base / "validation/sdist-wheel" / pure_name,
    ):
        check_archives(project_root=root, wheel=pure, sdist=sdist)
    distributions = [(base / "dist" / pure_name, "pure"), (sdist, "sdist")]
    for role, suffix in (
        ("linux", "manylinux_2_28_x86_64"),
        ("macos", "macosx_11_0_arm64"),
    ):
        directory = artifacts / f"kernelgecp-release-native-{role}"
        filename = f"kernelgecp-{version}-cp312-abi3-{suffix}.whl"
        _require(
            _files(directory) == {Path(filename)}, "unexpected native target or payload"
        )
        wheel = directory / filename
        check_archives(project_root=root, wheel=wheel, sdist=sdist, kind="native")
        with zipfile.ZipFile(wheel) as archive:
            inventory = json.loads(archive.read("kernelgecp/_native_inventory.json"))
        _require(
            inventory["source_revision"] == revision
            and inventory["source_tree"] == selected["source_tree"]
            and inventory["source_asset"]["sha256"] == report["sha256"]
            and inventory["source_asset"]["filename"] == report["archive"]
            and inventory["license_bundle_sha256"] == report["license_bundle_sha256"]
            and inventory["sbom_sha256"] == report["sbom_sha256"],
            "native/source inventory mismatch",
        )
        distributions.append((wheel, f"native-{role}"))
    _require(not output.exists(), "release output already exists")
    documents = [
        (root / name, "document")
        for name in ("METHODS.md", "CITATION.cff", "RELEASE_NOTES.md")
    ]
    if _is_12(version):
        documents = []
        for source, group in _required_12_records().values():
            path = root / source
            _check_12_content(path, group)
            blob = subprocess.run(
                ["git", "-C", str(root), "show", f"{revision}:{source}"],
                capture_output=True,
                check=False,
            )
            _require(
                blob.returncode == 0 and blob.stdout == path.read_bytes(),
                f"required 1.2 asset differs from Git source: {source}",
            )
            documents.append((path, group))
    groups = {
        "dist": distributions,
        "source": [(base / "source" / name, "source") for name in sorted(source_names)],
        "documents": documents,
    }
    manifest = {"schema": "kernelgecp.release.v1", "project": "kernelgecp", **selected}
    manifest["schema"] = "kernelgecp.release.v1"
    for group, records in groups.items():
        (output / group).mkdir(parents=True)
        manifest[group] = []
        for path, role in records:
            shutil.copyfile(path, output / group / path.name)
            manifest[group].append(
                {"filename": path.name, "role": role, "sha256": _hash(path)}
            )
    (output / "release-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    verify_payload(output)
    return manifest


def verify_payload(root: Path) -> dict[str, Any]:
    manifest = _json(root / "release-manifest.json")
    _require(
        manifest.get("schema") == "kernelgecp.release.v1"
        and manifest.get("project") == "kernelgecp",
        "invalid release manifest",
    )
    version = manifest.get("version")
    _require(
        manifest.get("repository") == REPOSITORY
        and isinstance(version, str)
        and re.fullmatch(
            r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)", version
        )
        is not None
        and manifest.get("tag") == f"v{version}"
        and all(
            isinstance(manifest.get(key), str)
            and re.fullmatch(r"[0-9a-f]{40}", manifest[key]) is not None
            for key in ("source_revision", "source_tree")
        ),
        "invalid release source/version identity",
    )
    expected = {Path("release-manifest.json")}
    all_names = {"release-manifest.json"}
    for group in ("dist", "source", "documents"):
        entries = manifest.get(group)
        _require(
            isinstance(entries, list) and bool(entries),
            f"missing manifest files: {group}",
        )
        names = set()
        for entry in entries:
            name, digest = entry["filename"], entry["sha256"]
            _require(
                isinstance(name, str)
                and name not in ("", ".", "..")
                and Path(name).name == name
                and "/" not in name
                and "\\" not in name
                and name not in names,
                "invalid or duplicate manifest filename",
            )
            _require(name not in all_names, "duplicate release asset basename")
            all_names.add(name)
            _require(
                isinstance(digest, str)
                and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
                "invalid manifest hash",
            )
            names.add(name)
            path = root / group / name
            _require(
                path.is_file() and not path.is_symlink() and _hash(path) == digest,
                "release file hash mismatch",
            )
            expected.add(Path(group) / name)
    if _is_12(version):
        required = _required_12_records()
        found = {
            entry["filename"]: entry.get("role") for entry in manifest["documents"]
        }
        for name, (_, group) in required.items():
            _require(
                found.get(name) == group, f"missing required 1.2 {group} asset: {name}"
            )
            _check_12_content(root / "documents" / name, group)
    _require(
        len(manifest["dist"]) == 4
        and {x["role"] for x in manifest["dist"]}
        == {"pure", "sdist", "native-linux", "native-macos"},
        "wrong distribution roles",
    )
    _require(_files(root) == expected, "unexpected release payload files")
    return manifest


def verify_pypi(root: Path) -> dict[str, Any]:
    manifest = verify_payload(root)
    expected = {item["filename"]: item["sha256"] for item in manifest["dist"]}
    for attempt in range(6):
        try:
            with urllib.request.urlopen(
                f"https://pypi.org/pypi/kernelgecp/{manifest['version']}/json",
                timeout=30,
            ) as response:
                data = json.load(response)
            _require(
                data["info"]["name"] == "kernelgecp"
                and data["info"]["version"] == manifest["version"],
                "public project/version mismatch",
            )
            observed = {
                item["filename"]: item["digests"]["sha256"] for item in data["urls"]
            }
            _require(
                observed == expected and len(data["urls"]) == 4,
                "public PyPI files/hashes disagree",
            )
            return {
                "version": manifest["version"],
                "files": observed,
                "status": "verified",
            }
        except urllib.error.HTTPError as exc:
            if exc.code != 404 or attempt == 5:
                raise ReleaseError("public PyPI metadata unavailable") from exc
        # A newly created project's JSON endpoint may lag its public file index.
        # Keep six strict attempts, but allow bounded propagation backoff.
        time.sleep(min(60, 5 * 2**attempt))
    raise ReleaseError("PyPI publication not visible")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    select = sub.add_parser("select")
    select.add_argument("--project-root", type=Path, required=True)
    select.add_argument("--ref", required=True)
    select.add_argument("--event-sha", required=True)
    select.add_argument("--repository", default=REPOSITORY)
    select.add_argument("--output", type=Path, required=True)
    select.add_argument("--github-output", type=Path)
    collect = sub.add_parser("collect")
    collect.add_argument("--project-root", type=Path, required=True)
    collect.add_argument("--selected-runs", type=Path, required=True)
    collect.add_argument("--artifacts", type=Path, required=True)
    collect.add_argument("--output", type=Path, required=True)
    for command in ("verify", "verify-pypi"):
        sub.add_parser(command).add_argument("--payload", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "select":
            result = select_release(
                args.project_root, args.ref, args.event_sha, args.repository
            )
            with args.output.open("x") as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write("\n")
            if args.github_output:
                with args.github_output.open("a") as stream:
                    stream.write(
                        f"native_run_id={result['native_run']['id']}\nversion={result['version']}\nsource_sha={result['source_revision']}\n"
                    )
        elif args.command == "collect":
            result = collect_artifacts(
                project_root=args.project_root,
                selected=_json(args.selected_runs),
                artifacts=args.artifacts,
                output=args.output,
            )
        else:
            result = (
                verify_payload(args.payload)
                if args.command == "verify"
                else verify_pypi(args.payload)
            )
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(1, f"release check failed: {exc}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
