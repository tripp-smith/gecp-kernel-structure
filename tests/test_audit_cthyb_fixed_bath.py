"""Failure-first checks for the independent fixed-bath evidence boundary."""

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from experiments import audit_cthyb_fixed_bath as audit

PARENT = Path("/tmp/gecp-ee-linux-7200-passed")
BINDING = Path("/tmp/gecp-ee-linux-7200-audit-binding.json")


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def test_json_rejects_duplicates_and_nonfinite(tmp_path):
    p = tmp_path / "x.json"
    for content in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":1e999}'):
        p.write_text(content)
        with pytest.raises(ValueError):
            audit.read_json(p)


def test_npz_rejects_unbounded_or_wrong_header(tmp_path):
    p = tmp_path / "x.npz"
    np.savez_compressed(p, bath_iw=np.zeros(512))
    with pytest.raises(ValueError, match="inventory"):
        audit.arrays(p)
    values = {
        k: np.zeros(shape, dtype=np.float64 if k == "occupancies" else np.complex128)
        for k, shape in audit.SHAPES.items()
    }
    values["g_iw"] = np.zeros((1, 2), dtype=complex)
    np.savez_compressed(p, **values)
    with pytest.raises(ValueError, match="header"):
        audit.arrays(p)


@pytest.mark.skipif(
    not BINDING.exists() or not PARENT.exists(),
    reason="external accepted parent not present",
)
def test_parent_binding_and_source_blob_are_valid(monkeypatch):
    monkeypatch.setattr(audit, "ROOT", Path.cwd())
    binding, bath, seeds = audit.parent_bath(PARENT, BINDING)
    assert binding["status"] == "passed"
    assert bath.shape == (512,)
    assert len(seeds) == 208
    assert not set(range(32100001, 32100513, 2)) & seeds


@pytest.mark.skipif(
    not BINDING.exists() or not PARENT.exists(),
    reason="external accepted parent not present",
)
def test_parent_binding_tamper_fails(tmp_path):
    p = tmp_path / "binding.json"
    shutil.copy(BINDING, p)
    p.write_bytes(p.read_bytes() + b" ")
    with pytest.raises(ValueError, match="binding hash"):
        audit.parent_bath(PARENT, p)


def test_orphan_missing_and_source_mismatch_fail_closed(tmp_path, monkeypatch):
    folder = tmp_path / "cohort"
    folder.mkdir()
    (folder / "chains").mkdir()
    config = {"beta": 10.0}
    bath = np.zeros(512, dtype=complex)
    parent = {
        "audit_receipt_sha256": "a",
        "manifest_sha256": "b",
        "freeze_sha256": "c",
        "result_sha256": "d",
    }
    monkeypatch.setattr(audit, "parent_bath", lambda *_: (parent, bath, set()))
    parent_dir = tmp_path / "parent"
    (parent_dir / "dense-final").mkdir(parents=True)
    write(parent_dir / "dense-final/replica-00.json", {"config": config})
    (parent_dir / "dense-final/replica-00.npz").write_bytes(b"bath")
    ref = {
        "binding_sha256": audit.PARENT_BINDING_SHA256,
        "audit_sha256": "a",
        "manifest_sha256": "b",
        "freeze_sha256": "c",
        "result_sha256": "d",
        "bath_receipt_sha256": audit.sha(parent_dir / "dense-final/replica-00.json"),
        "bath_archive_sha256": audit.sha(parent_dir / "dense-final/replica-00.npz"),
        "bath_bytes_sha256": hashlib.sha256(bath.tobytes()).hexdigest(),
    }
    manifest = {
        "schema": "kernelgecp.cthyb-fixed-bath.v1",
        "protocol": {
            "seeds": list(range(32100001, 32100513, 2)),
            "train_count": 128,
            "cycles": 1600000,
            "segment_seconds": 14400,
            "memory_limit_bytes": 4294967296,
            "config": config,
        },
        "parent": ref,
        "provenance": {
            "revision": "a" * 40,
            "sources": {"x.py": hashlib.sha256(b"good").hexdigest()},
            "runtime": {
                "python": "test",
                "platform": "Linux",
                "versions": {"numpy": "test"},
                "provider_files": {"x": "a" * 64},
            },
            "threads": {
                n: "1"
                for n in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                    "BLIS_NUM_THREADS",
                )
            },
        },
    }
    write(folder / "manifest.json", manifest)
    write(
        folder / "result.json",
        {
            "schema": "kernelgecp.cthyb-fixed-bath-result.v1",
            "status": "complete",
            "manifest_sha256": audit.sha(folder / "manifest.json"),
            "chains": [{}] * 256,
        },
    )
    monkeypatch.setattr(audit, "validate_calibration", lambda *_: None)
    monkeypatch.setattr(audit.subprocess, "check_output", lambda *a, **kw: "x.py\n")
    monkeypatch.setattr(audit, "git_blob", lambda *_: b"bad")
    manifest["calibration"] = {}
    write(folder / "manifest.json", manifest)
    result = audit.read_json(folder / "result.json")
    result["manifest_sha256"] = audit.sha(folder / "manifest.json")
    write(folder / "result.json", result)
    args = (
        folder,
        parent_dir,
        BINDING,
        tmp_path / "receipt",
        tmp_path / "binding",
        tmp_path / "cal",
        tmp_path / "cal-audit",
    )
    with pytest.raises(ValueError, match="generation source"):
        audit.audit_cohort(*args)
    monkeypatch.setattr(audit, "git_blob", lambda *_: b"good")
    with pytest.raises(ValueError, match="orphan or missing"):
        audit.audit_cohort(*args)
    (folder / "chains/orphan.json").write_text("{}")
    with pytest.raises(ValueError, match="orphan or missing"):
        audit.audit_cohort(*args)


@pytest.mark.skipif(
    not BINDING.exists() or not PARENT.exists(),
    reason="external accepted parent not present",
)
def test_full_synthetic_inventory_and_tamper(tmp_path, monkeypatch):
    # Reuse one real native sample solely to exercise parsing and provenance;
    # the fabricated repeated payloads do not represent independent chains.
    folder = tmp_path / "cohort"
    (folder / "chains").mkdir(parents=True)
    binding = audit.read_json(BINDING)
    original = audit.read_json(PARENT / "dense-final/replica-00.json")
    source = PARENT / "dense-final/replica-00.npz"
    bath = audit.arrays(source)["bath_iw"]
    seeds = list(range(32100001, 32100513, 2))
    ref = {
        "binding_sha256": audit.PARENT_BINDING_SHA256,
        "audit_sha256": binding["audit_receipt_sha256"],
        "manifest_sha256": binding["manifest_sha256"],
        "freeze_sha256": binding["freeze_sha256"],
        "result_sha256": binding["result_sha256"],
        "bath_receipt_sha256": audit.sha(PARENT / "dense-final/replica-00.json"),
        "bath_archive_sha256": audit.sha(source),
        "bath_bytes_sha256": hashlib.sha256(bath.tobytes()).hexdigest(),
    }
    threads = {
        n: "1"
        for n in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "BLIS_NUM_THREADS",
        )
    }
    manifest = {
        "schema": "kernelgecp.cthyb-fixed-bath.v1",
        "protocol": {
            "seeds": seeds,
            "train_count": 128,
            "cycles": 1600000,
            "segment_seconds": 14400,
            "memory_limit_bytes": 4294967296,
            "config": original["config"],
        },
        "parent": ref,
        "provenance": {
            "revision": "a" * 40,
            "sources": {"x.py": hashlib.sha256(b"good").hexdigest()},
            "runtime": {
                "python": "test",
                "platform": "Linux",
                "versions": {"numpy": "test"},
                "provider_files": {"x": "a" * 64},
            },
            "threads": threads,
        },
    }
    manifest["calibration"] = {}
    write(folder / "manifest.json", manifest)
    manifest_hash = audit.sha(folder / "manifest.json")
    chains = []
    for i, seed in enumerate(seeds):
        role = "train" if i < 128 else "heldout"
        stem = folder / "chains" / f"chain-{i:03d}"
        stem.with_suffix(".npz").hardlink_to(source)
        raw = {
            **original,
            "seed": seed,
            "identity": {"manifest_sha256": manifest_hash, "index": i, "role": role},
        }
        write(stem.with_suffix(".json"), raw)
        entry = {
            "index": i,
            "seed": seed,
            "role": role,
            "receipt_sha256": audit.sha(stem.with_suffix(".json")),
            "arrays_sha256": audit.sha(source),
        }
        timing = {
            "schema": "kernelgecp.cthyb-fixed-bath-timing.v1",
            "seed": seed,
            "elapsed_seconds": 1,
            "peak_rss_bytes": 1,
            "solve_status": 0,
            "average_sign": raw["average_sign"],
            "receipt_sha256": entry["receipt_sha256"],
            "arrays_sha256": entry["arrays_sha256"],
        }
        write(stem.with_suffix(".timing.json"), timing)
        entry["timing_sha256"] = audit.sha(stem.with_suffix(".timing.json"))
        chains.append(entry)
    write(
        folder / "result.json",
        {
            "schema": "kernelgecp.cthyb-fixed-bath-result.v1",
            "manifest_sha256": manifest_hash,
            "status": "complete",
            "chains": chains,
        },
    )
    monkeypatch.setattr(audit, "parent_bath", lambda *_: (binding, bath, set()))
    monkeypatch.setattr(audit, "validate_calibration", lambda *_: None)
    monkeypatch.setattr(
        audit.subprocess,
        "check_output",
        lambda cmd, **kw: "x.py\n" if cmd[1] == "ls-tree" else "b" * 40 + "\n",
    )
    monkeypatch.setattr(
        audit,
        "git_blob",
        lambda revision, path: (
            b"good" if path == "x.py" else Path(audit.__file__).read_bytes()
        ),
    )
    args = (
        folder,
        PARENT,
        BINDING,
        tmp_path / "receipt",
        tmp_path / "binding",
        tmp_path / "cal",
        tmp_path / "cal-audit",
    )
    receipt, _ = audit.audit_cohort(*args)
    assert receipt["claim_scope"] == "fixed_bath_raw_cohort_only_not_EE_campaign"
    (tmp_path / "receipt").unlink()
    (tmp_path / "binding").unlink()
    actual_arrays = audit.arrays

    def invalid_moment(path):
        values = actual_arrays(path)
        values["g_moments"][1, 0] = 0.5
        return values

    with monkeypatch.context() as attack:
        attack.setattr(audit, "arrays", invalid_moment)
        with pytest.raises(ValueError, match="zeroth/unit moment"):
            audit.audit_cohort(*args)
    with (folder / "chains/chain-000.json").open("a") as f:
        f.write(" ")
    with pytest.raises(ValueError, match="chain inventory"):
        audit.audit_cohort(*args)


def test_calibration_freeze_and_adversarial_changes(tmp_path, monkeypatch):
    revision = "a" * 40
    source_bytes = b"independent audit"
    source_hash = hashlib.sha256(source_bytes).hexdigest()
    monkeypatch.setattr(audit, "git_blob", lambda *_: source_bytes)
    sources = {"experiments/audit_cthyb_spectroscopy.py": source_hash}
    threads = {"OMP_NUM_THREADS": "1"}
    packages = (
        "numpy",
        "scipy",
        "triqs",
        "triqs_cthyb",
        "pydlr",
        "sparse_ir",
        "pylibsparseir",
    )
    versions = {k: "version" for k in packages}
    provider_files = {}
    trees = {}
    for name in (*packages, "triqs_maxent", "mpmath"):
        root = f"/fake/{name}"
        file = root + "/module.py"
        file_hash = hashlib.sha256(file.encode()).hexdigest()
        digest = hashlib.sha256(b"module.py\0" + bytes.fromhex(file_hash)).hexdigest()
        trees[name] = {"root": root, "files": 1, "sha256": digest}
        if name in packages:
            provider_files[file] = file_hash
    runtime = {
        "python": "python",
        "platform": "Linux",
        "versions": versions,
        "provider_files": provider_files,
    }
    provenance = {
        "revision": revision,
        "sources": sources,
        "threads": threads,
        "runtime": runtime,
    }
    snapshot = {
        "revision": revision,
        "sources": sources,
        "threads": threads,
        "python": "python",
        "platform": "Linux",
        "packages": {
            ("sparse-ir" if k == "sparse_ir" else k): v for k, v in versions.items()
        },
        "runtime_trees": trees,
    }
    methods = ["raw", "kernelgecp", "pydlr", "sparse_ir"]
    defaults = ["flat", "broad_gaussian"]
    cases = {}
    for method in methods:
        for model in defaults:
            path = tmp_path / f"{method}-{model}-r0.json"
            archive = path.with_suffix(".npz")
            archive.write_bytes(b"case data")
            write(
                path,
                {
                    "status": "passed",
                    "eps": 1e-10,
                    "rcond": 1e-12,
                    "array_sha256": audit.sha(archive),
                },
            )
            cases[path.name] = {
                "result_sha256": audit.sha(path),
                "array_sha256": audit.sha(archive),
            }
    cal = {
        "schema": "kernelgecp.cthyb-spectroscopy-calibration.v1",
        "status": "passed",
        "binding_sha256": audit.CAL_BINDING,
        "cohort": {
            "stage": "calibration",
            "role": "kernelgecp-final",
            "generation_revision": audit.CAL_REVISION,
            "seed_base": 29100000,
            "seeds": list(range(29140001, 29140129, 2)),
            "binding_sha256": audit.CAL_BINDING,
        },
        "protocol": {
            "methods": methods,
            "default_models": defaults,
            "eps": 1e-10,
            "rcond": 1e-12,
            "n_omega": 401,
            "n_alpha": 33,
            "minimizer_mu0": 1e-8,
            "covariance_rtol": 1e-4,
            "bootstrap_count": 8,
            "bootstrap_repeat": 0,
            "final_repeats": 3,
            "alpha_selection": "normalized_training_only_discrepancy",
        },
        "snapshot": snapshot,
        "cases": cases,
    }
    calibration = tmp_path / "calibration.json"
    audit_path = tmp_path / "calibration-audit.json"
    write(calibration, cal)
    cal_hash = audit.sha(calibration)
    receipt = {
        "schema": "kernelgecp.cthyb-spectroscopy-independent-audit.v1",
        "status": "passed",
        "kind": "calibration",
        "cases": 8,
        "campaign_sha256": cal_hash,
        "calibration_sha256": cal_hash,
        "binding_sha256": audit.CAL_BINDING,
        "auditor_revision": revision,
        "auditor_sha256": source_hash,
        "gates": dict.fromkeys(audit.CAL_GATES, True),
    }
    write(audit_path, receipt)
    ref = {
        "receipt_sha256": cal_hash,
        "audit_sha256": audit.sha(audit_path),
        "binding_sha256": audit.CAL_BINDING,
        "auditor_revision": revision,
        "auditor_sha256": source_hash,
    }
    audit.validate_calibration(calibration, audit_path, ref, provenance)
    # Rebinding a failed audit must not convert it into a pass.
    receipt["gates"]["heldout_rms"] = False
    write(audit_path, receipt)
    ref["audit_sha256"] = audit.sha(audit_path)
    with pytest.raises(ValueError, match="independent audit did not pass"):
        audit.validate_calibration(calibration, audit_path, ref, provenance)
    receipt["gates"]["heldout_rms"] = True
    write(audit_path, receipt)
    ref["audit_sha256"] = audit.sha(audit_path)
    runtime["provider_files"]["/fake/triqs/module.py"] = "0" * 64
    with pytest.raises(ValueError, match="provider files mismatch"):
        audit.validate_calibration(calibration, audit_path, ref, provenance)
    runtime["provider_files"]["/fake/triqs/module.py"] = hashlib.sha256(
        b"/fake/triqs/module.py"
    ).hexdigest()
    (tmp_path / "raw-flat-r0.npz").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="case hash mismatch"):
        audit.validate_calibration(calibration, audit_path, ref, provenance)
