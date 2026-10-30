"""Replay the accepted EE audit on immutable native measurements, without sampling."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import stat
import subprocess
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
GENERATION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
ZIP_SHA = "5fba88be6a6796d29947b5e3bc5bf360e155813dec18d15e506ebe1f51832e77"
BINDING_SHA = "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
AUDITOR_SHA = "7e86c5bf1ce81ca54069c285d02f069033b47f65db3d2f50a0f13ca7f473cc9c"
SOURCES = {
    "environments/locks/dlr-physics-pip.txt",
    "experiments/cthyb_dmft.py",
    "experiments/cthyb_exact_state.toml",
    "experiments/cthyb_loop_state.py",
    "experiments/dlr_adapters.py",
    "kernelgecp/basis.py",
    "kernelgecp/dmft.py",
    "kernelgecp/fitting.py",
    "kernelgecp/integrations/cthyb_dmft.py",
    "kernelgecp/integrations/triqs.py",
}
GATES = {
    "loop_controllers",
    "standard_errors",
    "physical_observable_agreement",
    "native_and_represented_self_consistency",
    "all_provider_replay_and_cutoff",
    "independent_native_refinement",
    "independent_low_frequency_quadrature",
}


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(archive: Path, destination: Path) -> None:
    """Reject unsafe, duplicate, oversized or non-regular ZIP members."""
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        require(len(entries) == 467, "archive inventory count mismatch")
        require(sum(x.file_size for x in entries) <= 32_000_000, "archive too large")
        seen: set[str] = set()
        for entry in entries:
            name = entry.filename
            path = PurePosixPath(name)
            mode = entry.external_attr >> 16
            require(
                name
                and not name.startswith("/")
                and "\\" not in name
                and all(part not in (".", "..") for part in name.split("/"))
                and not entry.is_dir()
                and not stat.S_ISLNK(mode)
                and (stat.S_IFMT(mode) in (0, stat.S_IFREG))
                and entry.file_size <= 1_000_000
                and name not in seen,
                "unsafe archive member",
            )
            seen.add(name)
            target = destination.joinpath(*path.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(entry) as src, target.open("xb") as out:
                while chunk := src.read(65536):
                    out.write(chunk)


def check_gates(audit: dict[str, Any]) -> None:
    require(audit.get("status") == "passed", "audit status is not passed")
    gates = audit.get("gates")
    require(isinstance(gates, dict) and set(gates) == GATES, "missing audit gates")
    require(all(value is True for value in gates.values()), "audit gate failed")
    require(audit.get("generation_revision") == GENERATION, "generation mismatch")
    require(audit.get("audit_sha256") == AUDITOR_SHA, "auditor hash mismatch")


def require_providers() -> None:
    for name in ("triqs", "triqs_cthyb", "pydlr", "sparse_ir"):
        importlib.import_module(name)


def verify_tree(
    folder: Path, binding: dict[str, Any], root: Path, *, check_sources: bool = True
) -> None:
    require(
        binding.get("generation_revision") == GENERATION, "binding generation mismatch"
    )
    for key, name in (
        ("manifest_sha256", "manifest.json"),
        ("freeze_sha256", "freeze.json"),
        ("result_sha256", "result.json"),
        ("audit_receipt_sha256", "independent-audit.json"),
    ):
        require(digest(folder / name) == binding[key], f"{name} digest mismatch")
    manifest = json.loads((folder / "manifest.json").read_text())
    require(manifest["source_revision"] == GENERATION, "manifest generation mismatch")
    sources = manifest["sources"]
    require(set(sources) == SOURCES, "scientific source inventory mismatch")
    if check_sources:
        for name, expected in sources.items():
            require(
                digest(root / name) == expected, f"scientific source changed: {name}"
            )
            blob = subprocess.check_output(["git", "show", f"HEAD:{name}"], cwd=root)
            require(
                hashlib.sha256(blob).hexdigest() == expected,
                f"HEAD source changed: {name}",
            )
        auditor = root / "experiments/audit_cthyb_dmft.py"
        require(digest(auditor) == AUDITOR_SHA, "auditor changed")
        blob = subprocess.check_output(
            ["git", "show", "HEAD:experiments/audit_cthyb_dmft.py"], cwd=root
        )
        require(hashlib.sha256(blob).hexdigest() == AUDITOR_SHA, "HEAD auditor changed")
    archives = binding["measurement_archives"]
    require(len(archives) == 208, "measurement inventory count mismatch")
    for name, expected in archives.items():
        path = PurePosixPath(name)
        require(
            not path.is_absolute()
            and ".." not in path.parts
            and path.suffix == ".json",
            "unsafe binding path",
        )
        receipt = folder.joinpath(*path.parts)
        require(digest(receipt) == expected, f"measurement changed: {name}")
        data = json.loads(receipt.read_text())
        require(
            digest(receipt.with_suffix(".npz")) == data["arrays_sha256"],
            f"measurement arrays changed: {name}",
        )
    check_gates(json.loads((folder / "independent-audit.json").read_text()))


# Additional inputs govern whether the native campaign must be regenerated.
# Workflow-only edits do not change the scientific computation.
CAMPAIGN_PATHS = (
    "experiments/cthyb*",
    "experiments/audit_cthyb_dmft.py",
    "experiments/dlr_adapters.py",
    "environments/locks",
    "scripts/providers/cthyb_checkpoint.py",
    "pyproject.toml",
)


def sources_unchanged(root: Path, sources: dict[str, str]) -> bool:
    expected = dict(sources)
    expected["experiments/audit_cthyb_dmft.py"] = AUDITOR_SHA
    try:
        for name, value in expected.items():
            if digest(root / name) != value:
                return False
            blob = subprocess.check_output(["git", "show", f"HEAD:{name}"], cwd=root)
            if hashlib.sha256(blob).hexdigest() != value:
                return False
        for revisions in ((GENERATION, "HEAD"), ("HEAD",)):
            command = ["git", "diff", "--quiet", *revisions, "--", *CAMPAIGN_PATHS]
            if subprocess.run(command, cwd=root, check=False).returncode != 0:
                return False
    except (OSError, subprocess.CalledProcessError):
        return False
    return True


def campaign_needed(root: Path, *, force: bool = False) -> bool:
    """Use pinned metadata only; this does not perform the independent replay."""
    if force:
        return True
    archive = root / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    binding_file = root / "data/cthyb/linux-audit-binding-7200ec0.json"
    if not archive.exists() or not binding_file.exists():
        return True
    require(digest(archive) == ZIP_SHA, "archive digest mismatch")
    require(digest(binding_file) == BINDING_SHA, "binding digest mismatch")
    binding = json.loads(binding_file.read_text())
    with tempfile.TemporaryDirectory(prefix="cthyb-preflight-") as temp:
        folder = Path(temp)
        extract(archive, folder)
        verify_tree(folder, binding, root, check_sources=False)
        manifest = json.loads((folder / "manifest.json").read_text())
        return not sources_unchanged(root, manifest["sources"])


def replay(root: Path, output: Path) -> dict[str, Any]:
    archive = root / "experiments/data/cthyb-linux-7200ec0-passed.zip"
    binding_file = root / "data/cthyb/linux-audit-binding-7200ec0.json"
    require(digest(archive) == ZIP_SHA, "archive digest mismatch")
    require(digest(binding_file) == BINDING_SHA, "binding digest mismatch")
    binding = json.loads(binding_file.read_text())
    require_providers()
    with tempfile.TemporaryDirectory(prefix="cthyb-replay-") as temp:
        folder = Path(temp)
        extract(archive, folder)
        verify_tree(folder, binding, root)
        from experiments.audit_cthyb_dmft import _plain, audit

        result = _plain(audit(folder))
    check_gates(result)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()
    require(result["audit_revision"] == revision, "replay revision mismatch")
    receipt = {
        "generation_revision": GENERATION,
        "replay_revision": revision,
        "archive_sha256": ZIP_SHA,
        "binding_sha256": BINDING_SHA,
        "recomputed_audit": result,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        print(
            json.dumps({"requires_campaign": campaign_needed(ROOT, force=args.force)})
        )
    else:
        require(
            args.output is not None and not args.force, "--output required for replay"
        )
        replay(ROOT, args.output)
        print("accepted frozen evidence independently replayed")


if __name__ == "__main__":
    main()
