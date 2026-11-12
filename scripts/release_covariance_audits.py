"""Audit saved covariance payloads and replay one bounded provider consistency check."""

import argparse
import ast
import hashlib
import io
import json
import math
import os
import re
import struct
import subprocess
import sys
import time
import zipfile
from pathlib import Path

from scripts.providers.verify_receipt import read_receipt
from scripts.release_epw_audits import as_text, bounded_runner
from scripts.release_epw_native import checked_path, inventory
from scripts.release_provider_audits import digest, require, source_snapshot

ROOT = Path(__file__).resolve().parents[1]
SIZES = (512, 1024, 2048, 4096, 8192, 16384)
SOURCE = (
    "kernelgecp/kernel_cholesky.py",
    "experiments/covariance_benchmark.py",
    "scripts/audit_covariance.py",
    "environments/gp/uv.lock",
    "environments/gp/pyproject.toml",
)
IDENTITY = (
    "scripts/release_covariance_audits.py",
    "scripts/release_covariance_provider.py",
    "scripts/release_epw_audits.py",
    "scripts/release_epw_native.py",
    "scripts/release_provider_audits.py",
    "scripts/providers/verify_receipt.py",
)
OUTER_SHA = "b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4"
CSV_SHA = "a94fbcfc71708b6ffdc033360163d91efb7f675e496f6b5860929aa96b96351b"
AUDIT_SHA = "9c73af4585397dac260aac04f38305f87eee3fa5c4828cad151ab635e677d72c"
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
CASE_KEYS = set(
    [
        "schema",
        "status",
        "canonical",
        "size",
        "repeat",
        "rank_budget",
        "actual_rank",
        "checks",
        "provenance",
        "dataset",
        "frozen_calibration_sha256",
        "normalization",
        "split",
        "accuracy",
        "standard_deviation_definition",
        "solve_relative_residual",
        "cg_relative_residual",
        "cg_info",
        "cg_iterations",
        "retained_numerical_bytes",
        "retained_components_bytes",
        "peak_rss_bytes",
        "factor_metadata",
        "timing_scopes",
        "timing_seconds",
        "predictive_rmse_original_units",
        "artifact",
        "scope",
    ]
)
PROVIDER_KEYS = set(
    [
        "schema",
        "status",
        "size",
        "provider",
        "versions",
        "provenance",
        "dataset",
        "metrics",
        "standard_deviation_definition",
        "checks",
        "frozen_calibration_sha256",
        "artifact",
    ]
)
SPLIT = set(
    [
        "train_ordinal_sha256",
        "query_ordinal_sha256",
        "train_first",
        "train_last",
        "query_first",
        "query_last",
    ]
)
ARRAYS = set(
    [
        "root",
        "mean",
        "std",
        "expected_mean",
        "expected_std",
        "query_y",
        "cg_solution",
        "cg_residual",
        "cg_rhs",
        "pivots",
        "approximate_solution",
    ]
)
PROVIDER_ARRAYS = set(
    ["root", "mean", "std", "expected_mean", "expected_std", "solution", "rhs"]
)
TIMEOUT = 4 * 60 * 60
MAX_PAYLOAD = 4 * 1024**3
MAX_FILE = 160 * 1024**2
MAX_DATA = 20 * 1024**2
BOOTSTRAP = """import importlib,importlib.metadata as md,json,os,platform,runpy,sys
import urllib.request
assert sys.platform == "linux" and platform.machine() in ("x86_64","AMD64")
assert sys.version_info[:2] == (3,12)
mods={n:importlib.import_module(n) for n in
      ("numpy","scipy","torch","gpytorch","linear_operator")}
versions={n:md.version(n) for n in
          ("numpy","scipy","torch","gpytorch","linear-operator")}
expected={"numpy":"2.2.6","scipy":"1.18.1","torch":"2.14.1+cpu",
          "gpytorch":"1.15.2","linear-operator":"0.6.1"}
assert versions == expected and mods["torch"].version.cuda is None
threads={n:os.environ.get(n) for n in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS",
         "MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")}
assert all(x == "1" for x in threads.values())
prefix=os.path.realpath(sys.prefix)
paths={n:os.path.realpath(m.__file__) for n,m in mods.items()}
assert all(os.path.commonpath((prefix,p)) == prefix for p in paths.values())
info={"platform":sys.platform,"machine":platform.machine(),"python":platform.python_version(),
      "executable":sys.executable,"prefix":sys.prefix,"versions":versions,
      "modules":paths,"threads":threads}
with open(sys.argv.pop(1),"x") as f: json.dump(info,f,sort_keys=True)
source=sys.argv.pop(1);sys.path.insert(0,source)
for name in ("experiments.covariance_benchmark", "scripts.audit_covariance",
             "kernelgecp.kernel_cholesky"):
    module=importlib.import_module(name)
    assert os.path.commonpath((source,os.path.realpath(module.__file__))) == source
def no_network(*args,**kwargs): raise RuntimeError("network disabled in saved audit")
urllib.request.urlopen=no_network
kind=sys.argv.pop(1)
if kind == "science":
    runpy.run_module("scripts.audit_covariance",run_name="__main__",alter_sys=True)
else: runpy.run_path(sys.argv.pop(1),run_name="__main__")
"""


def identity(repository):
    repo = Path(repository)
    rev = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True, timeout=60
    ).strip()
    require(re.fullmatch(r"[0-9a-f]{40}", rev), "bad controller revision")
    require(
        not subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=all", "--", *IDENTITY],
            cwd=repo,
            text=True,
            timeout=60,
        ).strip(),
        "controller source dirty",
    )
    hashes = {}
    for name in IDENTITY:
        path = repo / name
        require(path.is_file() and not path.is_symlink(), "unsafe controller source")
        blob = subprocess.check_output(
            ["git", "show", rev + ":" + name], cwd=repo, timeout=60
        )
        hashes[name] = digest(path)
        require(
            hashes[name] == hashlib.sha256(blob).hexdigest(), "controller blob mismatch"
        )
    return {"revision": rev, "source_sha256": hashes}


def regular(path, limit):
    require(
        path.is_file()
        and not path.is_symlink()
        and path.stat().st_nlink == 1
        and path.stat().st_size <= limit,
        "unsafe or oversized file",
    )
    return path


def obj(path):
    value = read_receipt(regular(path, 2_000_000))
    require(type(value) is dict, "JSON object required")

    def check(item):
        if isinstance(item, dict):
            for v in item.values():
                check(v)
        elif isinstance(item, list):
            for v in item:
                check(v)
        elif isinstance(item, float):
            require(math.isfinite(item), "nonfinite JSON value")

    check(value)
    return value


def keyset(value, keys, label):
    require(type(value) is dict and set(value) == set(keys), label + " keyset")


def finite(value):
    require(type(value) in (int, float), "nonfinite numerical value")
    try:
        require(math.isfinite(value), "nonfinite numerical value")
    except OverflowError as error:
        raise ValueError("nonfinite numerical value") from error
    return value


def dataset_check(path):
    regular(path, MAX_DATA)
    require(digest(path) == OUTER_SHA, "UCI outer hash")
    with zipfile.ZipFile(path) as outer:
        members = outer.infolist()
        names = [x.filename for x in members]
        require(len(names) == len(set(names)), "duplicate outer ZIP member")
        name = "PRSA2017_Data_20130301-20170228.zip"
        require(
            names.count(name) == 1 and outer.getinfo(name).file_size <= MAX_DATA,
            "nested archive",
        )
        nested = outer.read(name)
    with zipfile.ZipFile(io.BytesIO(nested)) as archive:
        members = archive.infolist()
        names = [x.filename for x in members]
        require(len(names) == len(set(names)), "duplicate nested ZIP member")
        station = [
            x
            for x in names
            if x.endswith("/PRSA_Data_Aotizhongxin_20130301-20170228.csv")
        ]
        require(
            len(station) == 1 and archive.getinfo(station[0]).file_size <= MAX_DATA,
            "station member",
        )
        require(
            hashlib.sha256(archive.read(station[0])).hexdigest() == CSV_SHA,
            "station hash",
        )
    return OUTER_SHA


def npz_headers(path, keys, size, rank, *, provider=False):
    regular(path, MAX_FILE)
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        names = [x.filename for x in entries]
        require(
            len(names) == len(set(names)) and set(names) == {x + ".npy" for x in keys},
            "NPZ keyset",
        )
        require(
            sum(x.file_size for x in entries) <= 140 * 1024**2, "NPZ expanded budget"
        )
        for item in entries:
            with archive.open(item) as f:
                require(f.read(6) == b"\x93NUMPY", "NPY magic")
                version = f.read(2)
                require(version in (b"\x01\x00", b"\x02\x00"), "NPY version")
                n = struct.unpack(
                    "<H" if version[0] == 1 else "<I",
                    f.read(2 if version[0] == 1 else 4),
                )[0]
                require(n <= 4096, "NPY header size")
                header = ast.literal_eval(f.read(n).decode("latin1"))
                require(
                    type(header) is dict
                    and set(header) == {"descr", "fortran_order", "shape"}
                    and type(header["fortran_order"]) is bool,
                    "NPY header",
                )
                name = item.filename[:-4]
                if name == "root":
                    shape = (size, rank)
                elif name == "pivots":
                    shape = (rank,)
                elif name in ("mean", "std", "query_y"):
                    shape = (512,)
                elif name in ("expected_mean", "expected_std"):
                    shape = (512,) if provider or size <= 2048 else (0,)
                else:
                    shape = (size,)
                require(
                    header["shape"] == shape
                    and header["descr"]
                    in (("<i8", "|i8") if name == "pivots" else ("<f8", "|f8")),
                    "NPY dtype/shape",
                )
                require(
                    item.file_size
                    == 6 + 2 + (2 if version[0] == 1 else 4) + n + 8 * math.prod(shape),
                    "NPY byte size",
                )


def payload_inventory(root):
    require(root.is_dir() and not root.is_symlink(), "payload directory")
    cases = {
        f"n{size}-repeat{rep}.{ext}"
        for size in SIZES
        for rep in range(3)
        for ext in ("json", "npz", "log")
    }
    require(
        {p.name for p in root.iterdir()}
        == {
            "n512-repeat0.json",
            "n512-repeat0.npz",
            "provider.json",
            "provider.npz",
            "small-audit.json",
            "full-audit.json",
            "full",
        },
        "payload root inventory",
    )
    folder = root / "full"
    require(
        folder.is_dir()
        and not folder.is_symlink()
        and {p.name for p in folder.iterdir()} == cases | {"campaign.json"},
        "full inventory",
    )
    values = inventory(root)
    require(
        len(values) == 61
        and sum((root / n).stat().st_size for n in values) <= MAX_PAYLOAD
        and all((root / n).stat().st_size <= MAX_FILE for n in values),
        "payload budget",
    )
    return values


def validate_case(case, size, repeat, frozen, revision):
    keyset(case, CASE_KEYS, "case")
    require(
        case["schema"] == "kernelgecp.covariance-case.v1"
        and case["status"] == "passed"
        and case["canonical"] is True
        and case["size"] == size
        and case["repeat"] == repeat
        and type(case["actual_rank"]) is int
        and 0
        < case["actual_rank"]
        <= case["rank_budget"]
        == frozen["selected_ranks"][str(size)],
        "case identity",
    )
    keyset(
        case["provenance"],
        {
            "source_revision",
            "source_hashes",
            "source_clean_and_tracked",
            "source_binding_scope",
            "development_lock_sha256_context_only",
            "versions",
            "python",
            "python_executable",
            "scientific_runtime_lock",
            "platform",
            "threads",
        },
        "provenance",
    )
    keyset(
        case["factor_metadata"],
        {
            "algorithm",
            "arithmetic",
            "certification_scope",
            "certification_status",
            "conditioning_history",
            "conditioning_method",
            "conditioning_status",
            "configuration",
            "converged",
            "convergence_scope",
            "precision_bits",
            "precision_status",
            "refinement",
            "requested_precision_bits",
            "schema_version",
            "stop_reason",
            "tie_counts",
            "tie_policy",
            "tie_rtol",
            "tie_status",
            "tie_tolerance_rule",
        },
        "factor metadata",
    )
    meta = case["factor_metadata"]
    keyset(
        meta["configuration"],
        {
            "column_evaluations",
            "retained_numerical_bytes",
            "returned_diagonal_error",
            "roundoff_slack",
            "threshold",
        },
        "factor configuration",
    )
    keyset(
        meta["refinement"],
        {
            "accepted",
            "evaluations",
            "max_evaluations",
            "message",
            "optimizer_status",
            "optimizer_success",
            "status",
        },
        "refinement",
    )
    require(type(meta["conditioning_history"]) is list, "conditioning history")
    for row in meta["conditioning_history"]:
        keyset(
            row,
            {
                "columns",
                "condition_number",
                "effective_rank",
                "sigma_min",
                "sigma_min_decimal",
                "stage",
                "status",
            },
            "conditioning record",
        )
    if size <= 2048:
        keyset(
            case["accuracy"], {"mean_rms", "mean_max", "std_rms", "std_max"}, "accuracy"
        )
    else:
        require(case["accuracy"] is None, "large dense accuracy")
    keyset(case["split"], SPLIT, "split")
    keyset(
        case["normalization"],
        {"features_mean", "features_std", "target_mean", "target_std", "fit_rows"},
        "normalization",
    )
    require(
        case["normalization"]["fit_rows"] == size
        and len(case["normalization"]["features_mean"])
        == len(case["normalization"]["features_std"])
        == 5,
        "normalization",
    )
    require(
        case["dataset"] == frozen["dataset"]
        and case["provenance"]["source_revision"] == revision
        and case["provenance"]["source_hashes"]
        == frozen["provenance"]["source_hashes"],
        "case provenance",
    )
    keyset(
        case["checks"],
        {
            "dense_accuracy",
            "approximate_solve",
            "original_kernel_cg",
            "retained_state",
            "rss",
            "source_clean",
        },
        "case gates",
    )
    require(all(v is True for v in case["checks"].values()), "case gates")
    keyset(case["artifact"], {"path", "sha256"}, "case artifact")
    require(case["artifact"]["path"] == f"n{size}-repeat{repeat}.npz", "case filename")
    keyset(
        case["timing_seconds"],
        {
            "application_workflow",
            "construction",
            "prediction",
            "dense_reference",
            "original_kernel_cg",
            "full_workflow",
        },
        "timings",
    )
    for k, v in case["timing_seconds"].items():
        require(
            v is None if k == "dense_reference" and size > 2048 else finite(v) > 0,
            "case timing",
        )
    require(
        finite(case["peak_rss_bytes"]) > 0
        and finite(case["retained_numerical_bytes"]) > 0,
        "case resources",
    )
    keyset(
        case["retained_components_bytes"],
        {
            "shifted_factor",
            "source_values",
            "source_ordinals",
            "normalized_split",
            "prediction_arrays",
            "solve_arrays",
        },
        "resource components",
    )
    require(
        all(
            type(v) is int and v > 0 for v in case["retained_components_bytes"].values()
        )
        and sum(case["retained_components_bytes"].values())
        == case["retained_numerical_bytes"],
        "retained resource values",
    )


def validate_audit(value, revision, *, small=False):
    if small:
        keyset(
            value,
            {"schema", "status", "source_revision", "case", "scope"},
            "small audit",
        )
        require(
            value["schema"] == "kernelgecp.covariance-small-audit.v1"
            and value["status"] == "passed"
            and value["source_revision"] == revision
            and value["scope"] == "small_real_CI_case_not_full_EI_campaign_acceptance",
            "small audit identity",
        )
        rows = [value["case"]]
        expected = {(512, 0)}
    else:
        keyset(
            value,
            {
                "schema",
                "status",
                "source_revision",
                "cases",
                "real_data",
                "actual_provider",
                "independent_full_kernel_cg_recomputed_once_per_size",
            },
            "full audit",
        )
        require(
            value["schema"] == "kernelgecp.covariance-audit.v1"
            and value["status"] == "passed"
            and value["source_revision"] == revision
            and value["real_data"] is True
            and value["actual_provider"] is True
            and value["independent_full_kernel_cg_recomputed_once_per_size"] is True,
            "full audit identity",
        )
        rows = value["cases"]
        expected = {(n, r) for n in SIZES for r in range(3)}
    require(type(rows) is list and len(rows) == len(expected), "audit cohort count")
    found = set()
    for row in rows:
        keyset(row, {"size", "repeat", "solve_residual", "cg_residual"}, "audit row")
        identity = (row["size"], row["repeat"])
        require(
            type(row["size"]) is int
            and type(row["repeat"]) is int
            and identity in expected
            and identity not in found
            and 0 <= finite(row["solve_residual"]) <= 1e-10
            and 0 <= finite(row["cg_residual"]) <= 1e-8,
            "audit row values",
        )
        found.add(identity)
    require(found == expected, "audit cohort")


def validate_payload(root, frozen, revision):
    manifest = obj(root / "full/campaign.json")
    keyset(manifest, {"schema", "provenance", "frozen_path", "cases"}, "campaign")
    require(
        manifest["schema"] == "kernelgecp.covariance-campaign.v1"
        and manifest["provenance"]["source_revision"] == revision
        and manifest["provenance"]["source_hashes"]
        == frozen["provenance"]["source_hashes"],
        "campaign provenance",
    )
    rows = manifest["cases"]
    require(type(rows) is list and len(rows) == 18, "campaign count")
    expected = {(n, r) for n in SIZES for r in range(3)}
    seen = set()
    for item in rows:
        keyset(
            item,
            {"path", "exit_code", "process_wall_seconds", "sha256"},
            "campaign item",
        )
        match = re.fullmatch(
            r"n(512|1024|2048|4096|8192|16384)-repeat([012])\.json", item["path"]
        )
        require(
            match is not None
            and type(item["exit_code"]) is int
            and item["exit_code"] == 0
            and finite(item["process_wall_seconds"]) > 0,
            "campaign child",
        )
        size, repeat = map(int, match.groups())
        require((size, repeat) not in seen, "duplicate case")
        seen.add((size, repeat))
        path = root / "full" / item["path"]
        require(digest(path) == item["sha256"], "case receipt hash")
        case = obj(path)
        validate_case(case, size, repeat, frozen, revision)
        arrays = path.with_suffix(".npz")
        require(digest(arrays) == case["artifact"]["sha256"], "case arrays hash")
        npz_headers(arrays, ARRAYS, size, case["actual_rank"])
    require(seen == expected, "cohort")
    small = obj(root / "n512-repeat0.json")
    validate_case(small, 512, 0, frozen, revision)
    require(
        digest(root / "n512-repeat0.npz") == small["artifact"]["sha256"],
        "small arrays hash",
    )
    npz_headers(root / "n512-repeat0.npz", ARRAYS, 512, small["actual_rank"])
    provider = obj(root / "provider.json")
    keyset(provider, PROVIDER_KEYS, "provider")
    keyset(
        provider["checks"],
        {"cpu_float64", "frozen_factor", "numerical"},
        "provider checks",
    )
    require(all(v is True for v in provider["checks"].values()), "provider checks")
    keyset(
        provider["metrics"],
        {
            "mean_rms",
            "mean_max",
            "std_rms",
            "std_max",
            "solve_max",
            "log_probability_error",
        },
        "provider metrics",
    )
    require(
        all(0 <= finite(v) <= 1e-9 for v in provider["metrics"].values()),
        "provider metrics",
    )
    require(
        provider["schema"] == "kernelgecp.covariance-provider.v1"
        and provider["status"] == "passed"
        and provider["provider"] == "gpytorch"
        and provider["size"] == 512
        and provider["dataset"] == frozen["dataset"]
        and provider["provenance"]["source_revision"] == revision
        and provider["provenance"]["source_hashes"]
        == frozen["provenance"]["source_hashes"],
        "provider identity",
    )
    keyset(provider["artifact"], {"path", "sha256"}, "provider artifact")
    require(
        provider["artifact"]["path"] == "provider.npz"
        and provider["artifact"]["sha256"] == digest(root / "provider.npz"),
        "provider hash",
    )
    npz_headers(
        root / "provider.npz", PROVIDER_ARRAYS, 512, small["actual_rank"], provider=True
    )
    validate_audit(obj(root / "full-audit.json"), revision)
    validate_audit(obj(root / "small-audit.json"), revision, small=True)


def validate_runtime(value, python):
    require(
        type(value) is dict
        and value.get("platform") == "linux"
        and value.get("machine") in ("x86_64", "AMD64")
        and re.fullmatch(r"3\.12\.\d+", value.get("python", ""))
        and value.get("versions")
        == {
            "numpy": "2.2.6",
            "scipy": "1.18.1",
            "torch": "2.14.1+cpu",
            "gpytorch": "1.15.2",
            "linear-operator": "0.6.1",
        }
        and value.get("threads") == dict.fromkeys(THREADS, "1"),
        "runtime version",
    )
    prefix = Path(value["prefix"]).resolve()
    executable = Path(value["executable"]).absolute()
    require(
        executable.resolve() == Path(python).resolve()
        and prefix in executable.parents
        and set(value["modules"])
        == {"numpy", "scipy", "torch", "gpytorch", "linear_operator"}
        and all(prefix in Path(p).resolve().parents for p in value["modules"].values()),
        "runtime origin",
    )


def source_hashes(source):
    return {
        name: digest(regular(source / name, MAX_FILE))
        for name in (*SOURCE, "data/covariance/calibration.json")
    }


def audit(
    payload,
    revision,
    python,
    data,
    output,
    *,
    repository=ROOT,
    runner=bounded_runner,
    snapshot=source_snapshot,
    identity_probe=identity,
):
    require(
        all(".." not in Path(p).parts for p in (payload, data, output)),
        "noncanonical path",
    )
    root, data = (
        checked_path(payload).resolve(strict=True),
        checked_path(data).resolve(strict=False),
    )
    rawout = checked_path(output)
    out = rawout.resolve(strict=False)
    require(
        root.is_dir()
        and rawout.parent.resolve(strict=True).is_dir()
        and out != root
        and root not in out.parents
        and not out.exists()
        and data != out
        and out not in data.parents,
        "unsafe output",
    )
    out.mkdir(exist_ok=False)
    receipt = {
        "schema": "kernelgecp.covariance-payload-audit.v1",
        "status": "failed",
        "source_revision": revision,
        "scope": (
            "saved arrays and fresh n512 provider consistency; "
            "historical performance; not release eligibility"
        ),
    }
    before = None
    data_hash = None
    deadline = time.monotonic() + TIMEOUT
    try:
        require(re.fullmatch(r"[0-9a-f]{40}", revision) is not None, "revision")
        before = payload_inventory(root)
        data_hash = dataset_check(data)
        controller = identity_probe(ROOT)
        env = dict(os.environ)
        env.update(dict.fromkeys(THREADS, "1"))
        env.pop("PYTHONPATH", None)
        with snapshot(repository, revision) as source:
            hashes = source_hashes(source)
            require(hashes["scripts/audit_covariance.py"] == AUDIT_SHA, "auditor pin")
            for name in (*SOURCE, "data/covariance/calibration.json"):
                blob = subprocess.check_output(
                    ["git", "show", revision + ":" + name], cwd=source, timeout=60
                )
                require(
                    hashlib.sha256(blob).hexdigest() == hashes[name], "source Git blob"
                )
            frozen = obj(source / "data/covariance/calibration.json")
            keyset(
                frozen,
                {
                    "binding_reproduction",
                    "calibration_only",
                    "dataset",
                    "exhausted_sizes",
                    "hyperparameter_candidates",
                    "hyperparameter_policy",
                    "lengthscale",
                    "noise",
                    "provenance",
                    "random_seed",
                    "rank_calibration",
                    "rank_policy",
                    "ranks",
                    "schema",
                    "selected_ranks",
                    "sizes",
                    "standard_deviation_definition",
                    "status",
                    "stochastic_operations",
                    "test_evaluated",
                    "timing_scopes",
                },
                "calibration",
            )
            require(
                frozen["schema"] == "kernelgecp.covariance-calibration.v1"
                and frozen["status"] == "calibrated"
                and frozen["calibration_only"] is True
                and frozen["test_evaluated"] is False
                and frozen["sizes"] == list(SIZES)
                and len(frozen["rank_calibration"]) == 30,
                "calibration identity",
            )
            for row in frozen["rank_calibration"]:
                keyset(
                    row,
                    {
                        "accuracy",
                        "actual_rank",
                        "eligible_small_dense_margin",
                        "predictive_rmse",
                        "rank_budget",
                        "seconds_pilot_only",
                        "size",
                        "solve_relative_residual",
                    },
                    "calibration row",
                )
            for row in frozen["hyperparameter_candidates"]:
                keyset(
                    row,
                    {
                        "lengthscale",
                        "noise",
                        "nll",
                        "rank512_accuracy",
                        "rank512_admissible",
                    },
                    "hyperparameter row",
                )
            old_revision = frozen["provenance"]["source_revision"]
            require(
                isinstance(old_revision, str)
                and re.fullmatch(r"[0-9a-f]{40}", old_revision),
                "calibration revision",
            )
            for name in SOURCE:
                blob = subprocess.check_output(
                    ["git", "show", old_revision + ":" + name], cwd=source, timeout=60
                )
                require(
                    hashlib.sha256(blob).hexdigest() == hashes[name],
                    "calibration source Git blob",
                )
            require(
                set(frozen["provenance"]["source_hashes"]) == set(SOURCE)
                and all(
                    hashes[n] == frozen["provenance"]["source_hashes"][n]
                    for n in SOURCE
                ),
                "frozen scientific source",
            )
            require(
                frozen["dataset"]["outer_sha256"] == OUTER_SHA
                and frozen["dataset"]["csv_sha256"] == CSV_SHA
                and frozen["dataset"]["doi"] == "10.24432/C5RK5G"
                and frozen["dataset"]["license"] == "CC BY 4.0",
                "dataset provenance",
            )
            validate_payload(root, frozen, revision)
            commands = [
                (
                    "science",
                    [
                        "--manifest",
                        str(root / "full/campaign.json"),
                        "--provider",
                        str(root / "provider.json"),
                        "--data",
                        str(data),
                        "--frozen",
                        str(source / "data/covariance/calibration.json"),
                        "--source-revision",
                        revision,
                        "--output",
                        str(out / "independent.json"),
                    ],
                ),
                (
                    "provider",
                    [
                        str(ROOT / "scripts/release_covariance_provider.py"),
                        "--payload",
                        str(root),
                        "--source-root",
                        str(source),
                        "--data",
                        str(data),
                        "--output",
                        str(out / "provider.json"),
                    ],
                ),
            ]
            runtimes = {}
            for kind, arguments in commands:
                regular(data, MAX_DATA)
                remain = deadline - time.monotonic()
                require(remain > 0, "audit deadline exhausted")
                cmd = [
                    str(Path(python).absolute()),
                    "-I",
                    "-c",
                    BOOTSTRAP,
                    str(out / (kind + "-runtime.json")),
                    str(source),
                    kind,
                    *arguments,
                ]
                try:
                    process = runner(
                        cmd,
                        cwd=source,
                        env=env,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=remain,
                    )
                except subprocess.TimeoutExpired as error:
                    (out / (kind + ".log")).write_text(
                        as_text(error.stdout)
                        + "\n"
                        + as_text(error.stderr)
                        + "\nTimeoutExpired\n"
                    )
                    raise ValueError(kind + " timed out") from error
                (out / (kind + ".log")).write_text(
                    as_text(process.stdout) + "\n" + as_text(process.stderr)
                )
                require(process.returncode == 0, kind + " auditor failed")
                runtimes[kind] = obj(out / (kind + "-runtime.json"))
                validate_runtime(runtimes[kind], python)
            checked = obj(out / "independent.json")
            validate_audit(checked, revision)
            provider = obj(out / "provider.json")
            keyset(
                provider,
                {
                    "schema",
                    "status",
                    "scope",
                    "historical_reported_log_probability_error",
                    "historical_log_probability_value_available",
                    "scipy_reference_log_probability",
                    "actual_replay",
                    "input_sha256",
                },
                "provider replay",
            )
            require(
                provider["schema"] == "kernelgecp.covariance-provider-replay.v1"
                and provider["status"] == "passed"
                and provider["historical_log_probability_value_available"] is False,
                "provider replay result",
            )
            actual = provider["actual_replay"]
            keyset(
                actual,
                {
                    "errors",
                    "log_probability",
                    "device",
                    "dtype",
                    "finite_counts",
                    "versions",
                    "module_paths",
                },
                "provider actual",
            )
            keyset(
                actual["errors"],
                {
                    "solution_max_abs_error",
                    "mean_max_abs_error",
                    "std_max_abs_error",
                    "log_probability_abs_error",
                },
                "provider errors",
            )
            require(
                actual["device"] == "cpu"
                and actual["dtype"] == "torch.float64"
                and all(0 <= finite(v) <= 1e-9 for v in actual["errors"].values())
                and finite(actual["log_probability"]) == actual["log_probability"]
                and finite(provider["scipy_reference_log_probability"])
                == provider["scipy_reference_log_probability"],
                "provider arithmetic",
            )
            require(
                actual["finite_counts"] == {"solution": 512, "mean": 512, "std": 512}
                and actual["versions"]
                == {
                    "numpy": "2.2.6",
                    "scipy": "1.18.1",
                    "torch": "2.14.1+cpu",
                    "gpytorch": "1.15.2",
                    "linear_operator": "0.6.1",
                }
                and actual["module_paths"] == runtimes["provider"]["modules"],
                "provider runtime binding",
            )
            expected_inputs = {
                "provider_json": root / "provider.json",
                "provider_npz": root / "provider.npz",
                "case_json": root / "full/n512-repeat0.json",
                "case_npz": root / "full/n512-repeat0.npz",
                "manifest": root / "full/campaign.json",
                "calibration": source / "data/covariance/calibration.json",
                "benchmark": source / "experiments/covariance_benchmark.py",
                "data": data,
            }
            require(
                provider["input_sha256"]
                == {k: digest(v) for k, v in expected_inputs.items()},
                "provider replay input binding",
            )
            require(
                source_hashes(source) == hashes and time.monotonic() <= deadline,
                "scientific source changed or timeout",
            )
        require(identity_probe(ROOT) == controller, "controller changed")
        receipt.update(
            status="passed",
            controller_identity=controller,
            source_sha256=hashes,
            dataset_sha256=data_hash,
            payload_sha256=before,
            runtime=runtimes,
            independent_sha256=digest(out / "independent.json"),
            provider_sha256=digest(out / "provider.json"),
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        subprocess.SubprocessError,
        zipfile.BadZipFile,
        SyntaxError,
        struct.error,
        OverflowError,
    ) as error:
        receipt["error"] = f"{type(error).__name__}: {error}"
        if not any(out.glob("*.log")):
            (out / "audit.log").write_text(receipt["error"] + "\n")
    receipt["inputs_unchanged"] = before is not None and data_hash is not None
    target = out / "receipt.json"
    target.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    try:
        if before is not None:
            require(payload_inventory(root) == before, "payload mutated")
        if "controller_identity" in receipt:
            require(
                identity_probe(ROOT) == receipt["controller_identity"],
                "controller mutated",
            )
        if data_hash is not None:
            require(dataset_check(data) == data_hash, "dataset mutated")
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        receipt.update(
            status="failed", inputs_unchanged=False, error=f"final input check: {error}"
        )
        target.write_text(
            json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--gp-python", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(
            args.payload, args.source_revision, args.gp_python, args.data, args.output
        )
    except (OSError, ValueError) as error:
        print(f"covariance audit failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
