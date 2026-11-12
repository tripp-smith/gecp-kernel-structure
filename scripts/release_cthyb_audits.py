"""Replay authenticated CT-HYB measurements; not release eligibility."""

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from scripts.providers.verify_receipt import read_receipt
from scripts.release_epw_audits import as_text, bounded_runner
from scripts.release_epw_native import checked_path, inventory
from scripts.release_provider_audits import digest, require, source_snapshot

ROOT = Path(__file__).resolve().parents[1]
GENERATION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
ZIP_SHA = "5fba88be6a6796d29947b5e3bc5bf360e155813dec18d15e506ebe1f51832e77"
BINDING_SHA = "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
AUDITOR_SHA = "7e86c5bf1ce81ca54069c285d02f069033b47f65db3d2f50a0f13ca7f473cc9c"
REPLAYER_SHA = "4752a544d39b33267ec5773081407875e5ff6db733702219abf8be53db5dbcc7"
REFERENCE_SHA = "2ab18e8d19f2d62b9cd6654991db30e03059a1e0402512f0d9b497ec3084443a"
GATES = {
    "loop_controllers",
    "standard_errors",
    "physical_observable_agreement",
    "native_and_represented_self_consistency",
    "all_provider_replay_and_cutoff",
    "independent_native_refinement",
    "independent_low_frequency_quadrature",
}
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
IDENTITY_PATHS = (
    "scripts/release_cthyb_audits.py",
    "scripts/release_epw_audits.py",
    "scripts/release_epw_native.py",
    "scripts/release_provider_audits.py",
    "scripts/providers/verify_receipt.py",
)
INPUTS = (
    "experiments/data/cthyb-linux-7200ec0-passed.zip",
    "data/cthyb/linux-audit-binding-7200ec0.json",
    "environments/locks/triqs-linux-64-conda-explicit.txt",
    "environments/locks/dlr-physics-pip.txt",
    "experiments/replay_cthyb_acceptance.py",
    "experiments/audit_cthyb_dmft.py",
    "data/cthyb/linux-independent-audit-7200ec0.json",
    "experiments/cthyb_dmft.py",
    "experiments/cthyb_exact_state.toml",
    "experiments/cthyb_loop_state.py",
    "experiments/dlr_adapters.py",
    "kernelgecp/basis.py",
    "kernelgecp/dmft.py",
    "kernelgecp/fitting.py",
    "kernelgecp/integrations/cthyb_dmft.py",
    "kernelgecp/integrations/triqs.py",
)
TIMEOUT = 3300
MAX_RECEIPT = 2_000_000
BOOTSTRAP = """import importlib,importlib.metadata as md,json,os,platform,runpy,sys
assert sys.platform == "linux" and platform.machine() in ("x86_64", "AMD64")
assert platform.python_version() == "3.12.16"
names=("triqs","triqs_cthyb","pydlr","sparse_ir","numpy","scipy","pylibsparseir")
modules={n:importlib.import_module(n) for n in names}
versions={n:md.version(n) for n in
          ("numpy","scipy","pydlr","sparse-ir","pylibsparseir")}
versions.update(triqs=importlib.import_module("triqs.version").version,cthyb=importlib.import_module("triqs_cthyb.version").version)
info={"platform":sys.platform,"machine":platform.machine(),
      "python":platform.python_version(),"executable":sys.executable,
      "prefix":sys.prefix,"versions":versions,
      "modules":{n:os.path.realpath(m.__file__) for n,m in modules.items()},
      "threads":{n:os.environ.get(n) for n in
                 ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS","BLIS_NUM_THREADS")}}
with open(sys.argv.pop(1),"x") as f: json.dump(info,f,sort_keys=True)
expected={"triqs":"4.0.2","cthyb":"4.0.1","numpy":"2.2.6","scipy":"1.18.1",
          "pydlr":"1.0.1","sparse-ir":"2.1.6","pylibsparseir":"0.12.0"}
assert versions == expected and all(v == "1" for v in info["threads"].values())
prefix=os.path.realpath(sys.prefix)
assert all(os.path.commonpath((prefix,p)) == prefix
           for p in info["modules"].values())
source=sys.argv.pop(1);sys.path.insert(0,source)
runpy.run_module("experiments.replay_cthyb_acceptance",run_name="__main__",alter_sys=True)
"""


def current_identity(repository):
    repository = Path(repository)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True, timeout=60
    ).strip()
    require(
        re.fullmatch(r"[0-9a-f]{40}", revision) is not None,
        "invalid controller revision",
    )
    status = subprocess.check_output(
        [
            "git",
            "status",
            "--porcelain",
            "--untracked-files=all",
            "--",
            *IDENTITY_PATHS,
        ],
        cwd=repository,
        text=True,
        timeout=60,
    )
    require(not status.strip(), "controller source dirty")
    hashes = {}
    for name in IDENTITY_PATHS:
        path = repository / name
        require(path.is_file() and not path.is_symlink(), "unsafe controller source")
        blob = subprocess.check_output(
            ["git", "show", revision + ":" + name], cwd=repository, timeout=60
        )
        hashes[name] = digest(path)
        require(
            hashes[name] == hashlib.sha256(blob).hexdigest(), "controller blob mismatch"
        )
    return {"revision": revision, "source_sha256": hashes}


def shape(value, reference, path="audit"):
    """Enforce complete historical shape; floating results need not be identical."""
    if isinstance(reference, dict):
        require(
            isinstance(value, dict) and set(value) == set(reference), path + " keys"
        )
        for k in reference:
            shape(value[k], reference[k], path + "." + k)
    elif isinstance(reference, list):
        require(
            isinstance(value, list) and len(value) == len(reference), path + " length"
        )
        for v, r in zip(value, reference, strict=True):
            shape(v, r, path)
    elif isinstance(reference, bool):
        require(type(value) is bool and value is reference, path + " boolean")
    elif isinstance(reference, int):
        require(type(value) is int and value == reference, path + " integer")
    elif isinstance(reference, float):
        require(type(value) in (int, float) and math.isfinite(value), path + " numeric")
    else:
        require(value == reference, path + " identity")


def validate_replay(value, revision, reference):
    require(
        isinstance(value, dict)
        and set(value)
        == {
            "generation_revision",
            "replay_revision",
            "archive_sha256",
            "binding_sha256",
            "recomputed_audit",
        },
        "replay fields",
    )
    require(
        value["generation_revision"] == GENERATION
        and value["replay_revision"] == revision
        and value["archive_sha256"] == ZIP_SHA
        and value["binding_sha256"] == BINDING_SHA,
        "replay identity",
    )
    result = value["recomputed_audit"]
    require(
        isinstance(result, dict)
        and result.get("schema") == "kernelgecp.cthyb-audit.v1"
        and result.get("status") == "passed"
        and result.get("generation_revision") == GENERATION
        and result.get("audit_revision") == revision
        and result.get("audit_sha256") == AUDITOR_SHA
        and result.get("independent_chain_count") == 208,
        "audit identity",
    )
    gates = result.get("gates")
    require(
        isinstance(gates, dict)
        and set(gates) == GATES
        and all(v is True for v in gates.values()),
        "audit gates",
    )
    expected = dict(reference)
    expected["audit_revision"] = revision
    shape(result, expected)


def validate_runtime(runtime, python):
    versions = {
        "python": "3.12.16",
        "triqs": "4.0.2",
        "cthyb": "4.0.1",
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "pydlr": "1.0.1",
        "sparse-ir": "2.1.6",
        "pylibsparseir": "0.12.0",
    }
    require(
        isinstance(runtime, dict)
        and runtime.get("platform") == "linux"
        and runtime.get("machine") in ("x86_64", "AMD64")
        and runtime.get("python") == versions.pop("python")
        and runtime.get("versions") == versions
        and runtime.get("threads") == dict.fromkeys(THREADS, "1"),
        "runtime mismatch",
    )
    prefix = runtime.get("prefix")
    modules = runtime.get("modules")
    require(
        isinstance(prefix, str)
        and Path(prefix).is_absolute()
        and isinstance(runtime.get("executable"), str)
        and Path(runtime["executable"]).is_absolute()
        and Path(runtime["executable"]).resolve() == Path(python).resolve()
        and Path(prefix).resolve() in Path(runtime["executable"]).resolve().parents
        and isinstance(modules, dict)
        and set(modules)
        == {
            "triqs",
            "triqs_cthyb",
            "pydlr",
            "sparse_ir",
            "numpy",
            "scipy",
            "pylibsparseir",
        }
        and all(
            isinstance(p, str)
            and Path(p).is_absolute()
            and Path(prefix) in Path(p).parents
            for p in modules.values()
        ),
        "runtime module paths",
    )


def read_json(path):
    require(
        path.is_file()
        and not path.is_symlink()
        and path.stat().st_nlink == 1
        and path.stat().st_size <= MAX_RECEIPT,
        "invalid receipt file",
    )
    return read_receipt(path)


def payload_inventory(root):
    require(
        {p.name for p in root.iterdir()} == {"cthyb-evidence-replay.json"},
        "invalid payload entries",
    )
    return inventory(root)


def source_inventory(source):
    hashes = {}
    for name in INPUTS:
        path = source / name
        require(
            path.is_file() and not path.is_symlink() and path.stat().st_nlink == 1,
            "unsafe scientific input",
        )
        hashes[name] = digest(path)
    return hashes


def source_blob(source, revision, name):
    return subprocess.check_output(
        ["git", "show", revision + ":" + name], cwd=source, timeout=60
    )


def audit(
    payload,
    revision,
    python,
    output,
    *,
    repository=ROOT,
    runner=bounded_runner,
    snapshot=source_snapshot,
    identity_probe=current_identity,
    read_blob=source_blob,
):
    require(
        ".." not in Path(payload).parts and ".." not in Path(output).parts,
        "noncanonical path",
    )
    root = checked_path(payload).resolve(strict=True)
    out_raw = checked_path(output)
    out = out_raw.resolve(strict=False)
    require(
        root.is_dir()
        and out_raw.parent.resolve(strict=True).is_dir()
        and out != root
        and root not in out.parents
        and not out.exists(),
        "unsafe output or payload",
    )
    before = None
    out.mkdir(exist_ok=False)
    receipt = {
        "schema": "kernelgecp.cthyb-payload-audit.v1",
        "status": "failed",
        "source_revision": revision,
        "scope": "saved measurement audit only; not release eligibility",
    }
    deadline = time.monotonic() + TIMEOUT
    try:
        before = payload_inventory(root)
        require(
            set(before) == {"cthyb-evidence-replay.json"}
            and (root / "cthyb-evidence-replay.json").stat().st_size <= MAX_RECEIPT,
            "invalid payload inventory",
        )
        require(
            isinstance(revision, str) and re.fullmatch(r"[0-9a-f]{40}", revision),
            "invalid source revision",
        )
        identity = identity_probe(ROOT)
        env = dict(os.environ)
        env.update(dict.fromkeys(THREADS, "1"))
        env.pop("PYTHONPATH", None)
        with snapshot(repository, revision) as source:
            source_hashes = source_inventory(source)
            require(
                source_hashes[INPUTS[0]] == ZIP_SHA
                and source_hashes[INPUTS[1]] == BINDING_SHA
                and source_hashes[INPUTS[5]] == AUDITOR_SHA
                and source_hashes[INPUTS[4]] == REPLAYER_SHA
                and source_hashes[INPUTS[6]] == REFERENCE_SHA,
                "frozen source/input mismatch",
            )
            for name in (INPUTS[4], INPUTS[5], INPUTS[6]):
                require(
                    hashlib.sha256(read_blob(source, revision, name)).hexdigest()
                    == source_hashes[name],
                    "scientific Git blob mismatch",
                )
            reference = read_json(source / INPUTS[6])
            validate_replay(
                read_json(root / "cthyb-evidence-replay.json"), revision, reference
            )
            cmd = [
                str(Path(python).absolute()),
                "-I",
                "-c",
                BOOTSTRAP,
                str(out / "runtime.json"),
                str(source),
                "--output",
                str(out / "independent.json"),
            ]
            try:
                remaining = deadline - time.monotonic()
                require(remaining > 0, "audit deadline exhausted")
                process = runner(
                    cmd,
                    cwd=source,
                    env=env,
                    text=True,
                    capture_output=True,
                    check=False,
                    timeout=remaining,
                )
            except subprocess.TimeoutExpired as error:
                (out / "auditor.log").write_text(
                    as_text(error.stdout)
                    + "\n"
                    + as_text(error.stderr)
                    + "\nTimeoutExpired\n"
                )
                raise ValueError("audit timed out") from error
            (out / "auditor.log").write_text(
                as_text(process.stdout) + "\n" + as_text(process.stderr)
            )
            require(process.returncode == 0, "independent auditor failed")
            runtime = read_json(out / "runtime.json")
            validate_runtime(runtime, python)
            validate_replay(read_json(out / "independent.json"), revision, reference)
            require(
                source_inventory(source) == source_hashes,
                "scientific source/input changed",
            )
            require(time.monotonic() <= deadline, "audit deadline exhausted")
        require(identity_probe(ROOT) == identity, "controller source changed")
        receipt.update(
            status="passed",
            controller_identity=identity,
            source_sha256=source_hashes,
            runtime=runtime,
            payload_sha256=before,
            independent_sha256=digest(out / "independent.json"),
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        subprocess.SubprocessError,
    ) as error:
        receipt["error"] = f"{type(error).__name__}: {error}"
        log = out / "auditor.log"
        if not log.exists():
            log.write_text(receipt["error"] + "\n")
    target = out / "receipt.json"
    receipt["payload_unchanged"] = before is not None
    target.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    try:
        if before is not None:
            require(payload_inventory(root) == before, "payload changed during audit")
    except (OSError, ValueError) as error:
        receipt.update(
            status="failed",
            payload_unchanged=False,
            error=f"final inventory failed: {error}",
        )
        target.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--payload", required=True, type=Path)
    p.add_argument("--source-revision", required=True)
    p.add_argument("--triqs-python", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    try:
        result = audit(
            args.payload, args.source_revision, args.triqs_python, args.output
        )
    except (OSError, ValueError) as error:
        print(f"CT-HYB payload audit failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
