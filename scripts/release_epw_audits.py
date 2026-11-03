"""Compose independent EPW archive and saved native-output audits.

Caller authenticates and safely materializes a static, exclusively owned payload.
No executable bytes are supplied; this does not establish release eligibility.
"""

import argparse
import contextlib
import hashlib
import json
import math
import os
import re
import signal
import subprocess
import sys
from pathlib import Path

from scripts.providers.verify_receipt import read_receipt
from scripts.release_epw_native import checked_path, inventory
from scripts.release_provider_audits import digest, require, source_snapshot

ROOT = Path(__file__).resolve().parents[1]
TIMEOUT = 360
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
CHANNELS = {
    "Z_real",
    "Z_imag",
    "gap_real_ev",
    "gap_imag_ev",
    "shift_real_ev",
    "shift_imag_ev",
}

# Set the address-space limit before numerical imports. The executable supplied
# by the caller must be the isolated pinned EPW Python runtime.
BOOTSTRAP = """import sys,os,platform,resource,json,runpy
assert (
    sys.platform == "linux" and platform.machine() in ("x86_64", "AMD64")
), "EPW replay requires Linux x86_64"
resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
import numpy,scipy,mpmath
info={"platform":sys.platform,"machine":platform.machine(),"python":platform.python_version(),"numpy":numpy.__version__,"scipy":scipy.__version__,"mpmath":mpmath.__version__,"executable":sys.executable}
assert sys.version_info[:2] == (3,12), "EPW replay requires Python 3.12"
versions=(numpy.__version__,scipy.__version__,mpmath.__version__)
assert versions == ("2.2.6","1.18.1","1.3.0"), (
    "EPW replay requires the pinned numerical runtime"
)
with open(sys.argv.pop(1),"x") as f: json.dump(info,f)
source=sys.argv.pop(1);sys.path.insert(0,source)
runpy.run_module("experiments.audit_epw_restart_roundtrip",run_name="__main__",alter_sys=True)
"""
NATIVE_BOOTSTRAP = """import sys,platform,resource,runpy
assert (
    sys.platform == "linux" and platform.machine() in ("x86_64", "AMD64")
), "EPW replay requires Linux x86_64"
assert sys.version_info[:2] == (3,12), "EPW replay requires Python 3.12"
resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
runpy.run_path(sys.argv.pop(1),run_name="__main__")
"""


def as_text(value):
    return (
        value.decode("utf-8", errors="replace")
        if isinstance(value, bytes)
        else (value or "")
    )


def nonnegative(value):
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value) and value >= 0
    except OverflowError:
        return False


def validate_archive(value):
    require(
        isinstance(value, dict)
        and set(value)
        == {
            "schema",
            "status",
            "cases",
            "scope",
            "interpreter_identity_scope",
            "native_acceptance",
            "public_decoder",
        }
        and value.get("scope")
        == (
            "finite 20 K fit-only holdout and all-field truth-assisted exact "
            "reconstruction; native solver comparator required"
        )
        and value.get("interpreter_identity_scope")
        == ("approved producer-recorded identity only; executable bytes not provided")
        and value.get("schema") == "kernelgecp.epw-restart-roundtrip-audit.v1"
        and value.get("status") == "passed"
        and value.get("native_acceptance") == "not_tested"
        and value.get("public_decoder") == "consistency only",
        "archive audit identity",
    )
    cases = value.get("cases")
    require(isinstance(cases, list) and len(cases) == 3, "archive audit cohort")
    metric_keys = {
        "green_scaled",
        "green_global_scaled",
        "z",
        "gap_ev",
        "chi_ev",
        "imaginary_defect",
    }
    budgets = {
        "green_scaled": 1e-6,
        "z": 1e-5,
        "gap_ev": 1e-6,
        "chi_ev": 1e-6,
        "imaginary_defect": 1e-8,
    }
    for i, case in enumerate(cases):
        require(
            isinstance(case, dict)
            and set(case)
            == {
                "repeat",
                "results",
                "bytes",
                "training_design_condition_max",
                "xor_correction_bytes",
                "xor_correction_values",
            }
            and type(case["repeat"]) is int
            and case["repeat"] == i,
            "archive case identity",
        )
        sizes = case["bytes"]
        require(
            isinstance(sizes, dict)
            and set(sizes)
            == {
                "gzip_ascii",
                "model_only",
                "model_plus_xor_bundle",
                "raw_ascii",
                "raw_binary_arrays",
            }
            and all(type(v) is int and v > 0 for v in sizes.values())
            and type(case["xor_correction_values"]) is int
            and case["xor_correction_values"] > 0
            and type(case["xor_correction_bytes"]) is int
            and case["xor_correction_bytes"] == 8 * case["xor_correction_values"]
            and nonnegative(case["training_design_condition_max"])
            and 1 <= case["training_design_condition_max"] <= 1e4,
            "archive size or conditioning",
        )
        splits = case["results"]
        require(
            isinstance(splits, dict) and set(splits) == {"heldout", "training", "full"},
            "archive split inventory",
        )
        for name, split in splits.items():
            require(
                isinstance(split, dict)
                and set(split)
                == {
                    "metrics",
                    "cardinal_condition_max",
                    "recovery_condition_max",
                    "z_derivative_gain_max",
                }
                and nonnegative(split["cardinal_condition_max"])
                and 1 <= split["cardinal_condition_max"] <= 1e4
                and nonnegative(split["recovery_condition_max"])
                and split["recovery_condition_max"] >= 1
                and nonnegative(split["z_derivative_gain_max"]),
                "archive split conditioning",
            )
            metrics = split["metrics"]
            require(
                isinstance(metrics, dict)
                and set(metrics) == metric_keys
                and all(nonnegative(v) for v in metrics.values()),
                "archive metrics",
            )
            if name == "heldout":
                require(
                    all(metrics[k] <= v for k, v in budgets.items()),
                    "archive heldout budget",
                )


def validate_native(value, revision, expected_inventory):
    require(
        isinstance(value, dict)
        and value.get("schema") == "kernelgecp.epw-native-payload-audit.v1"
        and value.get("status") == "passed"
        and value.get("source_revision") == revision
        and value.get("payload_unchanged") is True
        and value.get("payload_sha256") == expected_inventory
        and value.get("cases") == ["raw", "reconstructed", "wrong-frequency"]
        and value.get("comparison_scope")
        == "exact finite printed-number budgets, not continuous certification",
        "native audit identity",
    )
    comp = value.get("comparison")
    require(
        isinstance(comp, dict)
        and set(comp)
        == {
            "state_count",
            "frequency_count",
            "row_count",
            "atol",
            "rtol",
            "budget_reference",
            "channels",
            "passed",
        }
        and comp.get("budget_reference") == "absolute raw value"
        and comp.get("passed") is True
        and type(comp.get("row_count")) is int
        and comp["row_count"] == 630000
        and type(comp.get("state_count")) is int
        and comp["state_count"] == 140
        and type(comp.get("frequency_count")) is int
        and comp["frequency_count"] == 4500
        and type(comp.get("atol")) is float
        and comp["atol"] == 1e-5
        and type(comp.get("rtol")) is float
        and comp["rtol"] == 1e-3
        and isinstance(comp.get("channels"), dict)
        and set(comp["channels"]) == CHANNELS,
        "native comparison cohort",
    )
    for name, item in comp["channels"].items():
        require(
            isinstance(item, dict)
            and set(item)
            == {"failures", "max_absolute_error", "max_budget_ratio", "unit"}
            and type(item["failures"]) is int
            and item["failures"] == 0
            and nonnegative(item["max_absolute_error"])
            and nonnegative(item["max_budget_ratio"])
            and item["max_budget_ratio"] <= 1
            and item["unit"] == ("dimensionless" if name.startswith("Z_") else "eV"),
            "native comparison channel",
        )


IDENTITY_PATHS = (
    "scripts/release_epw_audits.py",
    "scripts/release_epw_native.py",
    "scripts/release_provider_audits.py",
    "scripts/providers/verify_receipt.py",
)


def current_identity(repository):
    """Bind executed helpers to clean tracked bytes at K, separate from producer P."""
    repo = Path(repository)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True, timeout=60
    ).strip()
    require(
        re.fullmatch(r"[0-9a-f]{40}", revision) is not None, "invalid auditor revision"
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
        cwd=repo,
        text=True,
        timeout=60,
    )
    require(not status.strip(), "auditor source is dirty or untracked")
    hashes = {}
    for path in IDENTITY_PATHS:
        source = ROOT / path
        require(source.is_file() and not source.is_symlink(), "unsafe auditor source")
        blob = subprocess.check_output(
            ["git", "show", f"{revision}:{path}"], cwd=repo, timeout=60
        )
        observed = digest(source)
        require(
            hashlib.sha256(blob).hexdigest() == observed,
            "auditor source differs from K",
        )
        hashes[path] = observed
    return {"revision": revision, "source_sha256": hashes}


def bounded_runner(cmd, *, cwd, env, text, capture_output, check, timeout):
    """Timeout the complete child process group, including spawned descendants."""
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        text=text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as error:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=2)
        except subprocess.TimeoutExpired:
            stdout = stderr = None
        # A descendant can survive TERM even after its output pipes close.
        with contextlib.suppress(ProcessLookupError):
            os.killpg(process.pid, signal.SIGKILL)
        if stdout is None and stderr is None:
            try:
                stdout, stderr = process.communicate(timeout=2)
            except subprocess.TimeoutExpired:
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()
                with contextlib.suppress(subprocess.TimeoutExpired):
                    process.wait(timeout=2)
        raise subprocess.TimeoutExpired(
            cmd, timeout, output=stdout or error.output, stderr=stderr or error.stderr
        ) from error
    return subprocess.CompletedProcess(cmd, process.returncode, stdout, stderr)


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
):
    root = checked_path(payload).resolve(strict=True)
    out_raw = checked_path(output)
    out = out_raw.resolve(strict=False)
    require(
        root.is_dir() and out_raw.parent.resolve(strict=True).is_dir(),
        "invalid directory",
    )
    require(
        out != root and root not in out.parents and not out.exists(), "unsafe output"
    )
    before = inventory(root)
    out.mkdir(exist_ok=False)
    receipt = {
        "schema": "kernelgecp.epw-payload-audit.v1",
        "status": "failed",
        "source_revision": revision,
        "scope": (
            "static authenticated payload; recorded native binary lineage only; "
            "not release eligibility"
        ),
    }
    try:
        require(
            isinstance(revision, str) and re.fullmatch("[0-9a-f]{40}", revision),
            "invalid revision",
        )
        freeze = read_receipt(root / "restart20/freeze.json")
        native = read_receipt(root / "epw-anisotropic.json")
        require(
            freeze["environment"]["revision"] == revision
            and native["source_revision"] == revision,
            "producer revision mismatch",
        )
        env = dict(os.environ)
        env.update(dict.fromkeys(THREADS, "1"))
        env.pop("PYTHONPATH", None)
        auditor_identity = identity_probe(ROOT)
        helper = ROOT / "scripts/release_epw_native.py"
        helper_sha256 = digest(helper)
        with snapshot(repository, revision) as source:
            commands = [
                (
                    "archive",
                    [
                        str(Path(python).absolute()),
                        "-I",
                        "-c",
                        BOOTSTRAP,
                        str(out / "runtime.json"),
                        str(source),
                        "--manifest",
                        str(root / "restart20/campaign.json"),
                        "--native-root-override",
                        str(root),
                        "--output",
                        str(out / "archive-independent.json"),
                    ],
                ),
                (
                    "native",
                    [
                        str(Path(python).absolute()),
                        "-I",
                        "-c",
                        NATIVE_BOOTSTRAP,
                        str(helper),
                        "--payload",
                        str(root),
                        "--source-revision",
                        revision,
                        "--repository",
                        str(source),
                        "--output",
                        str(out / "native"),
                    ],
                ),
            ]
            for name, cmd in commands:
                try:
                    process = runner(
                        cmd,
                        cwd=source,
                        env=env,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=TIMEOUT,
                    )
                except subprocess.TimeoutExpired as error:
                    (out / (name + ".log")).write_text(
                        as_text(error.stdout)
                        + "\n"
                        + as_text(error.stderr)
                        + "\nTimeoutExpired\n"
                    )
                    raise ValueError(name + " audit timed out") from error
                (out / (name + ".log")).write_text(
                    as_text(process.stdout) + "\n" + as_text(process.stderr)
                )
                require(process.returncode == 0, name + " auditor failed")
                if name == "archive":
                    runtime = read_receipt(out / "runtime.json")
                    require(
                        isinstance(runtime, dict)
                        and runtime.get("platform") == "linux"
                        and runtime.get("machine") in ("x86_64", "AMD64")
                        and isinstance(runtime.get("python"), str)
                        and re.fullmatch(r"3\.12\.\d+", runtime["python"])
                        and (
                            runtime.get("numpy"),
                            runtime.get("scipy"),
                            runtime.get("mpmath"),
                        )
                        == ("2.2.6", "1.18.1", "1.3.0")
                        and isinstance(runtime.get("executable"), str)
                        and runtime["executable"],
                        "invalid runtime identity",
                    )
                    validate_archive(read_receipt(out / "archive-independent.json"))
                else:
                    validate_native(
                        read_receipt(out / "native/receipt.json"), revision, before
                    )
        require(digest(helper) == helper_sha256, "native auditor changed during replay")
        require(
            identity_probe(ROOT) == auditor_identity,
            "auditor source changed during replay",
        )
        receipt.update(
            status="passed",
            auditor_identity=auditor_identity,
            native_auditor_sha256=helper_sha256,
            runtime=runtime,
            archive_sha256=digest(out / "archive-independent.json"),
            native_sha256=digest(out / "native/receipt.json"),
            payload_sha256=before,
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
    receipt["payload_unchanged"] = True
    target = out / "receipt.json"
    target.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    try:
        require(inventory(root) == before, "payload changed during audit")
    except (OSError, ValueError) as error:
        receipt.update(
            status="failed",
            payload_unchanged=False,
            error=f"final inventory failed: {error}",
        )
        target.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        with contextlib.suppress(OSError, ValueError):
            inventory(root)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--payload", required=True, type=Path)
    p.add_argument("--source-revision", required=True)
    p.add_argument("--epw-python", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    try:
        result = audit(a.payload, a.source_revision, a.epw_python, a.output)
    except (OSError, ValueError) as error:
        print(f"EPW payload audit failed: {error}", file=sys.stderr)
        return 1
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
