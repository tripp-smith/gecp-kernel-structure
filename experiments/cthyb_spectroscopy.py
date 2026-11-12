"""Actual fixed-bath CT-HYB spectroscopy with independent chain validation."""

from __future__ import annotations

import argparse
import copy
import dataclasses
import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import threading
import time
from contextlib import suppress
from pathlib import Path

import numpy as np

from experiments.cthyb_dmft import load_measurement
from experiments.spectroscopy import METHODS, SOURCE_FILES, THREADS, basis, write_json
from kernelgecp.integrations.maxent import (
    MaxEntConfig,
    SpectralDataset,
    continue_tau_maxent,
    covariance_modes,
    prepare_spectral_projection,
)

ROOT = Path(__file__).resolve().parents[1]
CALIBRATION_REVISION = "9d28053060329cd20647dc03e038a383c351adf1"
FINAL_REVISION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
CALIBRATION_SEED_BASE = 29100000
FINAL_SEED_BASE = 30100000
OBS = np.arange(0, 3073, 48)
QUERY = 24 + 192 * np.arange(16)
RTOL = 1e-4
IMAG_TOL = 1e-10


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def design_hash(array):
    return hashlib.sha256(
        np.ascontiguousarray(array, dtype="<f8").tobytes()
    ).hexdigest()


def identity_arithmetic(design):
    digest = design_hash(design)
    return {
        "method": "identity",
        "precision_bits": 53,
        "source_design_sha256": digest,
        "decomposed_design_sha256": digest,
        "source_design_shape": list(design.shape),
        "decomposed_design_shape": list(design.shape),
        "hash_scope": "C-order little-endian binary64 array bytes with recorded shapes",
        "rank_guard": None,
        "transformed_svd": None,
    }


def validate_case_projection(value, archive, method, *, profile="v2"):
    """Bind v2 receipts to the serialized source and decomposed designs."""
    if (
        value.get("schema") != schema_for("", profile)
        or value.get("method") != method
        or value.get("array_sha256") != sha(archive)
    ):
        raise ValueError("case schema, method or archive binding mismatch")
    arithmetic = value.get("projection_arithmetic")
    bits = 53 if method == "raw" else PROTOCOL["projection_precision_bits"]
    if (
        not isinstance(arithmetic, dict)
        or type(value.get("precision_bits")) is not int
        or value["precision_bits"] != bits
        or arithmetic.get("precision_bits") != bits
    ):
        raise ValueError("case projection precision mismatch")
    with np.load(archive, allow_pickle=False) as arrays:
        source = arrays["observation_design"]
        decomposed = arrays["decomposed_design"]
        if any(
            a.ndim != 2 or a.dtype != np.dtype("float64") or not np.isfinite(a).all()
            for a in (source, decomposed)
        ):
            raise ValueError("case design arrays invalid")
        if any(
            arithmetic.get(name + "_sha256") != design_hash(a)
            or arithmetic.get(name + "_shape") != list(a.shape)
            for name, a in (
                ("source_design", source),
                ("decomposed_design", decomposed),
            )
        ):
            raise ValueError("case design hash or shape mismatch")
        if method == "raw":
            if (
                arithmetic != identity_arithmetic(source)
                or not np.array_equal(source, np.eye(len(source)))
                or not np.array_equal(source, decomposed)
                or value.get("projection_svd") is not None
                or value.get("transformed_projection_svd") is not None
            ):
                raise ValueError("case raw projection identity mismatch")
        elif (
            arithmetic.get("method") != "mpmath_qr_then_binary64_cached_svd"
            or not isinstance(value.get("projection_svd"), dict)
            or not isinstance(value.get("transformed_projection_svd"), dict)
            or arithmetic.get("transformed_svd") != value["transformed_projection_svd"]
        ):
            raise ValueError("case projection diagnostics mismatch")


def mean_cov(samples):
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 2 or samples.shape[0] < 2 or not np.isfinite(samples).all():
        raise ValueError(
            "finite independent sample matrix with at least two rows required"
        )
    mean = samples.mean(axis=0)
    residual = samples - mean
    return mean, residual.T @ residual / (len(samples) * (len(samples) - 1))


PARENT_BINDING = "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
FIXED_SEEDS = [33100001 + 2 * i for i in range(256)]


def load_fixed_cohort(
    folder, receipt_path, binding_path, binding_sha256, *, profile="v2"
):
    fixed_seeds = [protocol_for(profile)["final_seed_base"] + 2 * i for i in range(256)]
    """Consume only the independently audited, content-bound fixed-bath cohort."""
    if len(binding_sha256) != 64 or sha(binding_path) != binding_sha256:
        raise ValueError("fixed-bath audit binding digest mismatch")
    binding = json.loads(binding_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    manifest = json.loads((folder / "manifest.json").read_text())
    result = json.loads((folder / "result.json").read_text())
    manifest_hash = sha(folder / "manifest.json")
    if (
        binding.get("schema") != f"kernelgecp.cthyb-fixed-bath-audit-binding.{profile}"
        or binding.get("status") != "passed"
        or binding.get("receipt_sha256") != sha(receipt_path)
        or binding.get("manifest_sha256") != manifest_hash
        or binding.get("result_sha256") != sha(folder / "result.json")
        or binding.get("parent_binding_sha256") != PARENT_BINDING
        or receipt.get("schema") != f"kernelgecp.cthyb-fixed-bath-audit.{profile}"
        or receipt.get("status") != "passed"
        or receipt.get("parent_binding_sha256") != PARENT_BINDING
        or receipt.get("manifest_sha256") != manifest_hash
        or receipt.get("result_sha256") != sha(folder / "result.json")
        or receipt.get("chain_count") != 256
        or receipt.get("train_count") != 128
        or receipt.get("heldout_count") != 128
        or receipt.get("audit_revision") != binding.get("audit_revision")
        or receipt.get("audit_source_sha256") != binding.get("audit_source_sha256")
    ):
        raise ValueError("fixed-bath audit or parent binding mismatch")
    rev = binding.get("audit_revision")
    if (
        not isinstance(rev, str)
        or len(rev) != 40
        or any(c not in "0123456789abcdef" for c in rev)
    ):
        raise ValueError("invalid fixed-bath auditor revision")
    blob = subprocess.check_output(
        ["git", "show", rev + ":experiments/audit_cthyb_fixed_bath.py"], cwd=ROOT
    )
    if hashlib.sha256(blob).hexdigest() != binding.get("audit_source_sha256"):
        raise ValueError("fixed-bath auditor source mismatch")
    protocol = manifest.get("protocol", {})
    provenance = manifest.get("provenance", {})
    calibration = manifest.get("calibration")
    keys = {
        "receipt_sha256",
        "audit_sha256",
        "binding_sha256",
        "auditor_revision",
        "auditor_sha256",
    }
    if (
        not isinstance(calibration, dict)
        or set(calibration) != keys
        or receipt.get("calibration") != calibration
        or binding.get("calibration") != calibration
        or calibration.get("auditor_revision") != provenance.get("revision")
        or calibration.get("auditor_sha256")
        != provenance.get("sources", {}).get("experiments/audit_cthyb_spectroscopy.py")
        or any(
            not isinstance(calibration[k], str)
            or len(calibration[k]) != (40 if k == "auditor_revision" else 64)
            or any(c not in "0123456789abcdef" for c in calibration[k])
            for k in keys
        )
    ):
        raise ValueError("fixed-bath frozen calibration binding mismatch")
    if (
        manifest.get("schema") != f"kernelgecp.cthyb-fixed-bath.{profile}"
        or manifest.get("parent", {}).get("binding_sha256") != PARENT_BINDING
        or protocol.get("seeds") != fixed_seeds
        or protocol.get("train_count") != 128
        or protocol.get("cycles") != 1600000
        or result.get("schema") != f"kernelgecp.cthyb-fixed-bath-result.{profile}"
        or result.get("status") != "complete"
        or result.get("manifest_sha256") != manifest_hash
        or len(result.get("chains", [])) != 256
        or not provenance.get("revision")
    ):
        raise ValueError("fixed-bath manifest, split or result mismatch")
    expected = {
        f"chains/chain-{i:03d}{ext}"
        for i in range(256)
        for ext in (".json", ".npz", ".timing.json")
    }
    inventory = binding.get("inventory")
    if not isinstance(inventory, dict) or set(inventory) != expected:
        raise ValueError("fixed-bath inventory mismatch")
    if any(
        (folder / name).is_symlink() or sha(folder / name) != digest
        for name, digest in inventory.items()
    ):
        raise ValueError("fixed-bath inventory hash mismatch")
    if {p.name for p in (folder / "chains").iterdir()} != {
        Path(x).name for x in expected
    }:
        raise ValueError("fixed-bath orphan or missing file")
    samples, bath = [], None
    for i, entry in enumerate(result["chains"]):
        role = "train" if i < 128 else "heldout"
        stem = folder / "chains" / f"chain-{i:03d}"
        if entry != {
            "index": i,
            "seed": fixed_seeds[i],
            "role": role,
            "receipt_sha256": sha(stem.with_suffix(".json")),
            "arrays_sha256": sha(stem.with_suffix(".npz")),
            "timing_sha256": sha(stem.with_suffix(".timing.json")),
        }:
            raise ValueError("fixed-bath chain entry mismatch")
        m = load_measurement(
            stem, {"manifest_sha256": manifest_hash, "index": i, "role": role}
        )
        if (
            m.seed != fixed_seeds[i]
            or m.cycles != 1600000
            or m.solve_status != 0
            or not 0.99 < m.average_sign <= 1.000000001
            or dataclasses.asdict(m.config) != protocol.get("config")
        ):
            raise ValueError("fixed-bath chain configuration or solve mismatch")
        if bath is None:
            bath = m.bath_iw.copy()
        if not np.array_equal(bath, m.bath_iw) or hashlib.sha256(
            m.bath_iw.tobytes()
        ).hexdigest() != manifest["parent"].get("bath_bytes_sha256"):
            raise ValueError("fixed-bath chain bath mismatch")
        x = m.g_tau
        if (
            x.shape != (3073, 2)
            or not np.isfinite(x).all()
            or np.max(np.abs(x.imag)) > IMAG_TOL
        ):
            raise ValueError("fixed-bath invalid Green function")
        samples.append(x.real.mean(axis=1))
    return np.array(samples), {
        "stage": "final",
        "role": "fixed-bath-256",
        "generation_revision": provenance["revision"],
        "seeds": fixed_seeds,
        "manifest_sha256": manifest_hash,
        "result_sha256": sha(folder / "result.json"),
        "receipt_sha256": sha(receipt_path),
        "binding_sha256": binding_sha256,
        "parent_binding_sha256": PARENT_BINDING,
        "calibration": calibration,
        "bath_bytes_sha256": manifest["parent"]["bath_bytes_sha256"],
        "measurement_receipt_sha256": [
            inventory[f"chains/chain-{i:03d}.json"] for i in range(256)
        ],
        "spin_operation": "paramagnetic_arithmetic_mean_up_down",
        "imaginary_absolute_tolerance": IMAG_TOL,
    }


def load_cohort(
    folder, receipt_path, binding_path, binding_sha256, stage, *, profile="v2"
):
    protocol_for(profile)
    """Validate audit, raw receipts, hashes, seeds, shape and one common bath."""
    if stage not in ("calibration", "final"):
        raise ValueError("unsupported stage")
    if len(binding_sha256) != 64 or sha(binding_path) != binding_sha256:
        raise ValueError("audit binding digest mismatch")
    binding = json.loads(binding_path.read_text())
    if stage == "final":
        if (
            binding.get("schema")
            != f"kernelgecp.cthyb-fixed-bath-audit-binding.{profile}"
        ):
            raise ValueError("final requires a fresh v2 fixed-bath audit binding")
        return load_fixed_cohort(
            folder, receipt_path, binding_path, binding_sha256, profile=profile
        )
    revision = CALIBRATION_REVISION
    seed_base = CALIBRATION_SEED_BASE
    if (
        binding.get("schema") != "kernelgecp.cthyb-audit-binding.v1"
        or binding.get("status") != "passed"
        or binding.get("generation_revision") != revision
        or binding.get("audit_receipt_sha256") != sha(receipt_path)
        or binding.get("result_sha256") != sha(folder / "result.json")
        or binding.get("manifest_sha256") != sha(folder / "manifest.json")
        or binding.get("freeze_sha256") != sha(folder / "freeze.json")
    ):
        raise ValueError("audit binding does not match source")
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("schema") != "kernelgecp.cthyb-audit.v1"
        or receipt.get("status") != "passed"
        or receipt.get("generation_revision") != revision
        or not receipt.get("gates")
        or not all(v is True for v in receipt["gates"].values())
    ):
        raise ValueError("independent EE audit did not pass for expected revision")
    audit_blob = subprocess.check_output(
        ["git", "show", receipt["audit_revision"] + ":experiments/audit_cthyb_dmft.py"],
        cwd=ROOT,
    )
    if hashlib.sha256(audit_blob).hexdigest() != receipt.get("audit_sha256"):
        raise ValueError("audit source hash mismatch")
    identity = json.loads((folder / "manifest.json").read_text())
    result = json.loads((folder / "result.json").read_text())
    config = identity["config"]
    if (
        identity.get("source_revision") != revision
        or result.get("identity") != identity
        or result.get("status") != "passed"
        or result.get("freeze_sha256") != sha(folder / "freeze.json")
        or config.get("seed_base") != seed_base
        or config.get("seed_mapping") != "odd_stride_2_v1"
        or config.get("final_replicas") != 64
        or config.get("n_tau") != 3073
        or config.get("beta") != 10
        or config.get("interaction") != 2
    ):
        raise ValueError("campaign identity, result, or configuration mismatch")
    role, offset = "kernelgecp-final", 40000
    paths = sorted((folder / role).glob("replica-*.json"))
    if len(paths) != 64 or paths != [
        folder / role / f"replica-{i:02d}.json" for i in range(64)
    ]:
        raise ValueError("cohort must have exactly 64 ordered replicas")
    samples, seeds, hashes, bath = [], set(), [], None
    for i, path in enumerate(paths):
        digest = sha(path)
        if (
            result["measurement_archives"].get(str(path.relative_to(folder))) != digest
            or binding.get("measurement_archives", {}).get(
                str(path.relative_to(folder))
            )
            != digest
        ):
            raise ValueError("measurement receipt hash mismatch")
        measurement = load_measurement(path.with_suffix(""), identity)
        if (
            measurement.seed != config["seed_base"] + offset + 1 + 2 * i
            or measurement.seed in seeds
        ):
            raise ValueError("duplicate or incorrect odd-stride seed")
        seeds.add(measurement.seed)
        if (
            measurement.config.beta != 10
            or measurement.config.interaction != 2
            or measurement.config.n_tau != 3073
            or measurement.cycles != 1600000
        ):
            raise ValueError("measurement physics or cycle count mismatch")
        if measurement.solve_status != 0 or measurement.average_sign <= 0.99:
            raise ValueError("native chain solve or sign failed")
        if bath is None:
            bath = measurement.bath_iw.copy()
        elif not np.array_equal(bath, measurement.bath_iw):
            raise ValueError("cohort contains different baths")
        values = measurement.g_tau
        if values.shape != (3073, 2) or not np.isfinite(values).all():
            raise ValueError("invalid G(tau) array")
        if np.max(np.abs(values.imag)) > IMAG_TOL:
            raise ValueError("imaginary G(tau) exceeds declared tolerance")
        samples.append(values.real.mean(axis=1))
        hashes.append(digest)
    return np.array(samples), {
        "stage": stage,
        "role": role,
        "generation_revision": revision,
        "seed_base": seed_base,
        "manifest_sha256": sha(folder / "manifest.json"),
        "freeze_sha256": sha(folder / "freeze.json"),
        "seeds": sorted(seeds),
        "receipt_sha256": sha(receipt_path),
        "binding_sha256": binding_sha256,
        "result_sha256": sha(folder / "result.json"),
        "measurement_receipt_sha256": hashes,
        "spin_operation": "paramagnetic_arithmetic_mean_up_down",
        "imaginary_absolute_tolerance": IMAG_TOL,
    }


def dataset_and_validation(samples, metadata):
    if samples.shape not in ((64, 3073), (256, 3073)) or not np.isfinite(samples).all():
        raise ValueError("cohort must contain 64 or 256 finite 3073-point chains")
    split = len(samples) // 2
    train, held = samples[:split], samples[split:]
    obs_mean, obs_cov = mean_cov(train[:, OBS])
    query_mean, held_cov = mean_cov(held[:, QUERY])
    tau = np.linspace(0, 10, 3073)
    data = SpectralDataset(
        10, tau[OBS], obs_mean, obs_cov, covariance_rtol=RTOL, provenance=metadata
    )
    return data, tau[QUERY], query_mean, held_cov


def validation_score(covariance, residual):
    _, modes = covariance_modes(covariance, rtol=RTOL)
    rms = modes.whitened_rms(residual)
    return rms, modes, bool(rms <= 2)


def weights(omega, density):
    cdf = np.r_[0.0, np.cumsum(np.diff(omega) * (density[:-1] + density[1:]) / 2)]
    return np.diff(np.interp([-4, -1, 1, 4], omega, cdf))


def mode_record(modes):
    return {
        "eigenvalues": modes.eigenvalues.tolist(),
        "threshold": modes.threshold,
        "relative_cutoff": modes.relative_cutoff,
        "absolute_cutoff": modes.absolute_cutoff,
        "effective_rank": modes.effective_rank,
        "discarded_eigenvalues": modes.eigenvalues[
            modes.eigenvalues <= modes.threshold
        ].tolist(),
    }


# The protocol is fixed before calibration and final validation. Repeats are
# fresh solver invocations on the same cohort, not independent data realizations.
PROTOCOL = {
    "profile": "corrective_projection_v2",
    "projection_precision_bits": 256,
    "methods": list(METHODS),
    "default_models": ["flat", "broad_gaussian"],
    # One common native-data candidate for every method and model; the ladder
    # is a declared history, not an automatic or method-dependent fallback.
    "basis_control_ladder": [
        {
            "name": "baseline",
            "eps": 1e-8,
            "rcond": 1e-10,
            "outcome": "failed_partial_calibration_pydlr_flat",
            "source_revision": "d5070509b675706e3e18164de96660530ac0489a",
            "failed_whitened_rms": 2.0585741570745255,
            "failed_projection_rank": 20,
            "failed_projected_covariance_rank": 19,
            "progress_sha256": (
                "c84df6e1dc8d2b2afce0f5ea7954ebdc8ee444cfac569d2fb803e582e7de1747"
            ),
            "failed_result_sha256": (
                "03c506bc6c864b8acf7c7f88ba9b08b16c3a17c13c67de97cace5ff28bd5739b"
            ),
            "failed_array_sha256": (
                "ebe689065afd27d657047af795710adb373e96bdb3319d05d55273feb339c523"
            ),
        },
        {"name": "candidate1", "eps": 1e-10, "rcond": 1e-12},
    ],
    "selected_basis_control": "candidate1",
    "eps": 1e-10,
    "rcond": 1e-12,
    "n_omega": 401,
    "n_alpha": 33,
    "minimizer_mu0": 1e-8,
    "covariance_rtol": RTOL,
    "bootstrap_count": 8,
    "bootstrap_repeat": 0,
    "bootstrap_seed": 193,
    "calibration_repeats": 1,
    "final_repeats": 3,
    "validation": "independent_heldout_chain_mean_conditional_on_training_fit",
    "alpha_selection": "normalized_training_only_discrepancy",
    "consumed_final_failure_array_sha256": (
        "407c85405e9f56135442db29bb2a4c75df8ef77b37e3f831a58c754e562817ba"
    ),
    "calibration_generation_revision": CALIBRATION_REVISION,
    "parent_generation_revision": FINAL_REVISION,
    "final_generation_revision_policy": "equal_campaign_snapshot_revision",
    "final_seed_base": 33100001,
    "final_seed_stride": 2,
    "final_chain_count": 256,
    "final_training_count": 128,
    "final_heldout_count": 128,
}


def protocol_for(profile="v2"):
    if profile not in ("v2", "v3", "v4"):
        raise ValueError("unsupported spectroscopy profile")
    protocol = copy.deepcopy(PROTOCOL)
    if profile in ("v3", "v4"):
        protocol.update(
            profile="native_interior_v3",
            final_seed_base=34100001,
            mass_mode="constrained",
            mass_curvature="projected",
            mass_interior_fraction=1e-12,
            minimizer_marquardt=True,
            maxiter=1000,
            development_scope="9d28053_cohort_consumed_development",
        )
    if profile == "v4":
        protocol.update(
            profile="native_shifted_curvature_v4",
            mass_curvature_floor_relative=1e-10,
        )
    return protocol


def schema_for(name, profile="v2"):
    protocol_for(profile)
    return f"kernelgecp.cthyb-spectroscopy{name}.{profile}"


def validate_basis_controls(*, profile="v2"):
    protocol = protocol_for(profile)
    """Fail closed if the declared candidate and executed controls diverge."""
    selected = [
        item
        for item in protocol["basis_control_ladder"]
        if item["name"] == protocol["selected_basis_control"]
    ]
    if len(selected) != 1 or any(
        selected[0][key] != protocol[key] for key in ("eps", "rcond")
    ):
        raise ValueError("selected native basis controls do not match protocol")


def snapshot(*, profile="v2"):
    protocol_for(profile)
    """Require all scientific source blobs to match the committed revision."""
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    tracked = subprocess.check_output(
        ["git", "ls-files", "kernelgecp", "experiments"], cwd=ROOT, text=True
    ).splitlines()
    files = sorted(
        set(SOURCE_FILES)
        | {
            "experiments/cthyb_spectroscopy.py",
            "experiments/cthyb_dmft.py",
            "experiments/audit_cthyb_spectroscopy.py",
        }
        | {p for p in tracked if p.startswith("kernelgecp/") and p.endswith(".py")}
    )
    if profile == "v4" and "kernelgecp/integrations/_maxent_curvature.py" not in files:
        raise ValueError("v4 curvature source must be tracked")
    hashes = {}
    for name in files:
        blob = subprocess.check_output(["git", "show", revision + ":" + name], cwd=ROOT)
        digest = sha(ROOT / name)
        if digest != hashlib.sha256(blob).hexdigest():
            raise ValueError("scientific source differs from committed blob: " + name)
        hashes[name] = digest
    threads = {name: os.environ.get(name) for name in (*THREADS, "BLIS_NUM_THREADS")}
    if any(value != "1" for value in threads.values()):
        raise ValueError("campaign requires declared single-thread environment")
    modules = [
        "triqs_maxent",
        "triqs_maxent.kernels",
        "triqs_maxent.version",
        "triqs",
        "triqs_cthyb",
    ]
    imported = {}
    for name in modules:
        module = importlib.import_module(name)
        path = Path(module.__file__).resolve()
        imported[name] = {"path": str(path), "sha256": sha(path)}
    # Hash installed code and native extensions, excluding generated caches.
    runtime_trees = {}
    for package in (
        "triqs_maxent",
        "triqs",
        "triqs_cthyb",
        "numpy",
        "scipy",
        "mpmath",
        "pydlr",
        "sparse_ir",
        "pylibsparseir",
    ):
        module = importlib.import_module(package)
        root = Path(module.__file__).resolve().parent
        digest = hashlib.sha256()
        count = 0
        for file in sorted(root.rglob("*")):
            if (
                file.is_file()
                and file.suffix in (".py", ".so", ".dylib")
                and "__pycache__" not in file.parts
            ):
                digest.update(str(file.relative_to(root)).encode() + b"\0")
                digest.update(bytes.fromhex(sha(file)))
                count += 1
        if not count:
            raise ValueError(
                "runtime package contains no source or native modules: " + package
            )
        runtime_trees[package] = {
            "root": str(root),
            "files": count,
            "sha256": digest.hexdigest(),
        }
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir")
    }
    version = importlib.import_module("triqs_maxent.version").version
    if str(version) != "4.0.1":
        raise ValueError("requires TRIQS MaxEnt 4.0.1")
    versions["triqs_maxent"] = str(version)
    versions["triqs"] = str(importlib.import_module("triqs.version").version)
    versions["triqs_cthyb"] = str(
        importlib.import_module("triqs_cthyb.version").version
    )
    if versions["triqs"] != "4.0.2" or versions["triqs_cthyb"] != "4.0.1":
        raise ValueError("unsupported TRIQS/CT-HYB provider versions")
    return {
        "revision": revision,
        "sources": hashes,
        "threads": threads,
        "python": sys.version,
        "executable": str(Path(sys.executable).resolve()),
        "executable_sha256": sha(Path(sys.executable).resolve()),
        "platform": platform.platform(),
        "packages": versions,
        "provider_modules": imported,
        "runtime_trees": runtime_trees,
    }


def verify_acquisition_freeze(folder, frozen):
    """Match acquisition provenance to the already calibrated source/runtime freeze."""
    provenance = json.loads((folder / "manifest.json").read_text())["provenance"]
    runtime = provenance.get("runtime", {})
    if (
        provenance.get("revision") != frozen.get("revision")
        or runtime.get("python") != frozen.get("python")
        or runtime.get("platform") != frozen.get("platform")
        or provenance.get("threads") != frozen.get("threads")
        or any(
            provenance.get("sources", {}).get(k) != v
            for k, v in frozen["sources"].items()
        )
    ):
        raise ValueError(
            "acquisition source or runtime differs from calibration freeze"
        )
    names = {"sparse_ir": "sparse-ir"}
    for name, value in runtime.get("versions", {}).items():
        if frozen["packages"].get(names.get(name, name)) != value:
            raise ValueError(
                "acquisition provider version differs from calibration freeze"
            )
    files = runtime.get("provider_files", {})
    if not files:
        raise ValueError("acquisition provider files missing")
    for tree in frozen["runtime_trees"].values():
        root = Path(tree["root"])
        selected = sorted(
            (Path(p), h)
            for p, h in files.items()
            if Path(p).is_relative_to(root)
            and Path(p).suffix in (".py", ".so", ".dylib")
            and "__pycache__" not in Path(p).parts
        )
        # MaxEnt and mpmath are analysis-only, absent from acquisition.
        if not selected and root.name in ("triqs_maxent", "mpmath"):
            continue
        digest = hashlib.sha256()
        for path, value in selected:
            digest.update(str(path.relative_to(root)).encode() + b"\0")
            digest.update(bytes.fromhex(value))
        if len(selected) != tree["files"] or digest.hexdigest() != tree["sha256"]:
            raise ValueError(
                "acquisition provider tree differs from calibration freeze"
            )


def verify_calibration(path, digest, snapshot_value, *, profile="v2"):
    protocol = protocol_for(profile)
    validate_basis_controls(profile=profile)
    if not path or not digest or len(digest) != 64 or sha(path) != digest:
        raise ValueError("final requires pinned passed calibration receipt")
    receipt = json.loads(path.read_text())
    if (
        receipt.get("schema") != schema_for("-calibration", profile)
        or receipt.get("status") != "passed"
        or receipt.get("protocol") != protocol
        or receipt.get("snapshot") != snapshot_value
    ):
        raise ValueError(
            "calibration freeze does not match protocol, source, runtime or EE binding"
        )
    cohort = receipt.get("cohort", {})
    if (
        cohort.get("stage") != "calibration"
        or cohort.get("role") != "kernelgecp-final"
        or cohort.get("generation_revision") != CALIBRATION_REVISION
        or cohort.get("seed_base") != CALIBRATION_SEED_BASE
        or cohort.get("seeds")
        != [CALIBRATION_SEED_BASE + 40001 + 2 * i for i in range(64)]
        or cohort.get("binding_sha256") != receipt.get("binding_sha256")
        or not isinstance(receipt.get("binding_sha256"), str)
        or len(receipt["binding_sha256"]) != 64
    ):
        raise ValueError("calibration cohort binding or identity mismatch")
    expected = {
        f"{method}-{model}-r0.json"
        for method in METHODS
        for model in protocol["default_models"]
    }
    entries = receipt.get("cases", {})
    if set(entries) != expected:
        raise ValueError("calibration case inventory mismatch")
    for name, hashes in entries.items():
        case = path.parent / name
        archive = case.with_suffix(".npz")
        if (
            sha(case) != hashes["result_sha256"]
            or sha(archive) != hashes["array_sha256"]
        ):
            raise ValueError("calibration case hash mismatch")
        value = json.loads(case.read_text())
        if (
            value.get("status") != "passed"
            or value.get("array_sha256") != sha(archive)
            or value.get("eps") != protocol["eps"]
            or value.get("rcond") != protocol["rcond"]
            or value.get("config", {}).get("minimizer_mu0") != protocol["minimizer_mu0"]
        ):
            raise ValueError("calibration case failed")
        if profile == "v2" and "mass_curvature_floor_relative" in value.get(
            "config", {}
        ):
            raise ValueError("v2 calibration case curvature floor mismatch")
        if profile in ("v3", "v4"):
            expected_config = MaxEntConfig(
                default_model=value.get("config", {}).get("default_model"),
                n_omega=401,
                n_alpha=33,
                minimizer_mu0=1e-8,
                mass_mode="constrained",
                mass_curvature="projected",
                mass_interior_fraction=1e-12,
                minimizer_marquardt=True,
                **({"mass_curvature_floor_relative": 1e-10} if profile == "v4" else {}),
            ).as_record()
            if value.get("config") != expected_config:
                raise ValueError("v3 calibration case controls mismatch")
        validate_case_projection(value, archive, name.split("-", 1)[0], profile=profile)
    return receipt


def check_fixed_calibration(cohort, calibration_path, audit_path, calibration, audit):
    expected = {
        "receipt_sha256": sha(calibration_path),
        "audit_sha256": sha(audit_path),
        "binding_sha256": calibration["binding_sha256"],
        "auditor_revision": audit["auditor_revision"],
        "auditor_sha256": audit["auditor_sha256"],
    }
    if cohort.get("calibration") != expected:
        raise ValueError("final calibration differs from pre-acquisition freeze")


def validate_cohort_pair(calibration, final):
    if (
        calibration["generation_revision"] == final["generation_revision"]
        or calibration["binding_sha256"] == final["binding_sha256"]
        or set(calibration["seeds"]) & set(final["seeds"])
    ):
        raise ValueError("calibration and final cohort bindings or seeds overlap")


def verify_calibration_audit(
    path, digest, calibration_digest, binding_digest, *, profile="v2"
):
    if not path or not digest or len(digest) != 64 or sha(path) != digest:
        raise ValueError("final requires a pinned independent calibration audit")
    audit = json.loads(path.read_text())
    if (
        audit.get("schema") != schema_for("-independent-audit", profile)
        or audit.get("status") != "passed"
        or audit.get("kind") != "calibration"
        or audit.get("calibration_sha256") != calibration_digest
        or audit.get("binding_sha256") != binding_digest
        or not audit.get("gates")
        or not all(value is True for value in audit["gates"].values())
    ):
        raise ValueError(
            "independent calibration audit did not pass or binding mismatch"
        )
    blob = subprocess.check_output(
        [
            "git",
            "show",
            audit["auditor_revision"] + ":experiments/audit_cthyb_spectroscopy.py",
        ],
        cwd=ROOT,
    )
    if hashlib.sha256(blob).hexdigest() != audit.get("auditor_sha256"):
        raise ValueError("independent calibration auditor hash mismatch")
    return audit


def run(
    folder,
    receipt,
    binding,
    digest,
    stage,
    method,
    model,
    output,
    *,
    bootstrap=0,
    seed=193,
    eps=1e-10,
    rcond=1e-12,
    n_omega=401,
    n_alpha=33,
    minimizer_mu0=1e-8,
    precision_bits=53,
    profile="v2",
):
    protocol = protocol_for(profile)
    if profile in ("v3", "v4") and (
        eps != protocol["eps"]
        or rcond != protocol["rcond"]
        or n_omega != protocol["n_omega"]
        or n_alpha != protocol["n_alpha"]
        or minimizer_mu0 != protocol["minimizer_mu0"]
        or precision_bits != protocol["projection_precision_bits"]
    ):
        raise ValueError("v3 fixed scientific controls mismatch")
    if method not in METHODS or model not in ("flat", "broad_gaussian"):
        raise ValueError("unsupported method or default model")
    if bootstrap not in (0, *range(2, 65)):
        raise ValueError("bootstrap must be zero or between 2 and 64")
    workflow_start = time.perf_counter()
    samples, metadata = load_cohort(
        folder, receipt, binding, digest, stage, profile=profile
    )
    data, query, reference, validation_cov = dataset_and_validation(samples, metadata)
    construction_start = time.perf_counter()
    projection = None
    design = np.eye(len(data.tau))
    decomposed = design
    if method != "raw":
        b = basis(method, data.beta, eps)
        projection = prepare_spectral_projection(
            b,
            data.tau,
            data.tau,
            rcond=rcond,
            metadata={"method": method, "eps": eps},
            precision_bits=precision_bits,
        )
        design = projection.source_design
        decomposed = projection.decomposed_design
    construction_seconds = time.perf_counter() - construction_start
    projection_start = time.perf_counter()
    fitted = data if projection is None else projection.apply(data)
    projection_seconds = time.perf_counter() - projection_start
    config = MaxEntConfig(
        n_omega=n_omega,
        n_alpha=n_alpha,
        default_model=model,
        minimizer_mu0=minimizer_mu0,
        **(
            {
                "mass_mode": "constrained",
                "mass_curvature": "projected",
                "mass_interior_fraction": 1e-12,
                "minimizer_marquardt": True,
            }
            if profile in ("v3", "v4")
            else {}
        ),
        **({"mass_curvature_floor_relative": 1e-10} if profile == "v4" else {}),
    )
    continuation_start = time.perf_counter()
    spectrum = continue_tau_maxent(fitted, config)
    continuation_seconds = time.perf_counter() - continuation_start
    prediction_start = time.perf_counter()
    predicted = spectrum.evaluate_tau(data.beta, query)
    prediction_seconds = time.perf_counter() - prediction_start
    rms, modes, heldout_pass = validation_score(validation_cov, predicted - reference)
    p = np.eye(len(data.tau)) if projection is None else projection.matrix
    from triqs_maxent import LinearOmegaMesh
    from triqs_maxent.kernels import TauKernel

    native_kernel = TauKernel(
        data.tau.copy(),
        LinearOmegaMesh(config.omega_min, config.omega_max, config.n_omega),
        beta=data.beta,
    ).K
    boot_weights, boot_spectra, boot_indices, boot_records = [], [], [], []
    rng = np.random.default_rng(seed)
    bootstrap_start = time.perf_counter()
    for _ in range(bootstrap):
        train_count = len(samples) // 2
        indices = rng.integers(0, train_count, size=train_count)
        boot_indices.append(indices)
        sample_mean, sample_cov = mean_cov(samples[indices][:, OBS])
        boot_data = SpectralDataset(
            10,
            data.tau,
            sample_mean,
            sample_cov,
            covariance_rtol=RTOL,
            provenance={"bootstrap": "independent_chain_resampling"},
        )
        boot_data = boot_data if projection is None else projection.apply(boot_data)
        trial = continue_tau_maxent(boot_data, config)
        boot_weights.append(weights(trial.omega, trial.density))
        boot_spectra.append(trial.density)
        boot_records.append(
            {
                "selected_index": trial.selected_index,
                "converged": trial.converged,
                "iterations": trial.iterations,
                "effective_rank": trial.effective_rank,
            }
        )
        if profile in ("v3", "v4"):
            boot_records[-1].update(
                raw_mass=trial.raw_mass,
                provider_diagnostics=trial.diagnostics,
            )
    bootstrap_seconds = time.perf_counter() - bootstrap_start
    gates = {
        "positive": bool(np.min(spectrum.density) >= -1e-6),
        "normalized": bool(abs(spectrum.mass - 1) <= 1e-6),
        "heldout_whitened_rms": heldout_pass,
        "provider_selected_converged": bool(
            spectrum.converged[spectrum.selected_index]
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    archive = output.with_suffix(".npz")
    if output.exists() or archive.exists():
        raise FileExistsError("refusing to overwrite output")
    np.savez_compressed(
        archive,
        omega=spectrum.omega,
        density=spectrum.density,
        observation_tau=data.tau,
        observation_values=data.values,
        covariance=data.covariance,
        projected_values=fitted.values,
        projected_covariance=fitted.covariance,
        projection=p,
        observation_design=design,
        decomposed_design=decomposed,
        native_kernel=native_kernel,
        projected_kernel=p @ native_kernel,
        query=query,
        query_mean=reference,
        validation_covariance=validation_cov,
        prediction=predicted,
        alphas=spectrum.alphas,
        chi2=spectrum.chi2,
        trial_spectra=spectrum.trial_spectra,
        bootstrap_weights=np.array(boot_weights),
        bootstrap_spectra=np.array(boot_spectra),
        bootstrap_indices=np.array(boot_indices),
    )
    result = {
        "schema": schema_for("", profile),
        "status": "passed" if all(gates.values()) else "failed",
        "gates": gates,
        "source": metadata,
        "method": method,
        "config": config.as_record(),
        "eps": eps,
        "rcond": rcond,
        "observation_indices": OBS.tolist(),
        "query_indices": QUERY.tolist(),
        "covariance_estimator": "sum_residual_outer_over_N_times_N_minus_1",
        "covariance_rtol": RTOL,
        "validation_covariance": (
            "heldout_chain_mean_covariance_conditional_on_training_fit; "
            "excludes_nonlinear_prediction_uncertainty"
        ),
        "training_covariance_modes": mode_record(data.modes),
        "projected_covariance_modes": mode_record(fitted.modes),
        "validation_covariance_modes": mode_record(modes),
        "precision_bits": 53 if projection is None else precision_bits,
        "projection_arithmetic": identity_arithmetic(design)
        if projection is None
        else json.loads(projection.projection_provenance_json),
        "transformed_projection_svd": None
        if projection is None or projection.transformed_diagnostics is None
        else json.loads(projection.projection_provenance_json)["transformed_svd"],
        "projection_svd": None
        if projection is None
        else dataclasses.asdict(projection.diagnostics),
        "provider_diagnostics": spectrum.diagnostics,
        "normalized_training_whitened_rms": spectrum.normalized_training_whitened_rms,
        "units": {"tau": "inverse_energy", "omega": "energy", "G_tau": "dimensionless"},
        "training_rank": data.modes.effective_rank,
        "projected_rank": fitted.modes.effective_rank,
        "validation_rank": modes.effective_rank,
        "whitened_rms": rms,
        "mass": spectrum.mass,
        "raw_mass": spectrum.raw_mass,
        "selected_index": spectrum.selected_index,
        "converged": spectrum.converged,
        "iterations": spectrum.iterations,
        "window_weights": weights(spectrum.omega, spectrum.density).tolist(),
        "bootstrap_count": bootstrap,
        "bootstrap_seed": seed,
        "bootstrap_records": boot_records,
        "bootstrap_window_standard_deviation": np.std(
            boot_weights, axis=0, ddof=1
        ).tolist()
        if bootstrap
        else None,
        "claim_scope": (
            "finite_grid_fixed_bath_empirical_diagnostics_no_native_truth_"
            "no_unique_spectrum_or_coverage_claim"
        ),
        "array_sha256": sha(archive),
        "timing_seconds": {
            "construction": construction_seconds,
            "projection": projection_seconds,
            "continuation": continuation_seconds,
            "bootstrap": bootstrap_seconds,
            "physical_prediction": prediction_seconds,
            "workflow_before_receipt": time.perf_counter() - workflow_start,
        },
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        * (1 if sys.platform == "darwin" else 1024),
        "archive_bytes": archive.stat().st_size,
    }
    write_json(output, result)
    return result


GLOBAL_SECONDS = 330 * 60
CASE_SECONDS = 60 * 60


def launch_worker(command, log_path, timeout, on_started=None):
    """Own a fresh process group, reap it on every exit, and retain its log."""
    started = time.monotonic()
    with log_path.open("xb") as log:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            if on_started is not None:
                on_started(process.pid)
            try:
                code = process.wait(timeout=timeout)
                status = "passed" if code == 0 else "failed"
            except subprocess.TimeoutExpired:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                code, status = process.returncode, "timeout"
        except BaseException:
            if process.poll() is None:
                with suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    return {
        "pid": process.pid,
        "command": command,
        "exit_code": code,
        "elapsed_seconds": time.monotonic() - started,
        "status": status,
        "log_file": log_path.name,
        "log_sha256": sha(log_path),
    }


def worker(job_path, job_digest, *, profile="v2"):
    protocol = protocol_for(profile)
    validate_basis_controls(profile=profile)
    if len(job_digest) != 64 or sha(job_path) != job_digest:
        raise ValueError("worker job hash mismatch")
    job = json.loads(job_path.read_text())
    if (
        job.get("schema") != schema_for("-worker-job", profile)
        or job.get("protocol") != protocol
    ):
        raise ValueError("invalid worker job protocol")
    method, model, repeat, stage = (
        job["method"],
        job["model"],
        job["repeat"],
        job["stage"],
    )
    if (
        stage not in ("calibration", "final")
        or method not in METHODS
        or model not in protocol["default_models"]
        or type(repeat) is not int
        or not 0 <= repeat < protocol[stage + "_repeats"]
    ):
        raise ValueError("worker case outside frozen inventory")
    if job.get("parent_pid") != os.getppid():
        raise ValueError("worker parent process mismatch")
    output = Path(job["output"])
    if output.name != f"{method}-{model}-r{repeat}.json":
        raise ValueError("worker output does not match case")
    frozen = snapshot(profile=profile)
    if frozen != job["snapshot"]:
        raise ValueError("worker source/runtime snapshot mismatch")
    if os.getpgrp() != os.getpid():
        raise ValueError("worker requires its own process group")
    parent_pid = job["parent_pid"]

    def watch_parent():
        while True:
            time.sleep(1)
            if os.getppid() != parent_pid:
                os.killpg(os.getpgrp(), signal.SIGKILL)

    threading.Thread(target=watch_parent, daemon=True).start()
    if stage == "final":
        calibration_value = verify_calibration(
            Path(job["calibration"]), job["calibration_sha256"], frozen, profile=profile
        )
        calibration_audit_value = verify_calibration_audit(
            Path(job["calibration_audit"]),
            job["calibration_audit_sha256"],
            job["calibration_sha256"],
            calibration_value["binding_sha256"],
            profile=profile,
        )
    if stage == "final":
        _, cohort = load_cohort(
            Path(job["folder"]),
            Path(job["audit_receipt"]),
            Path(job["audit_binding"]),
            job["binding_sha256"],
            stage,
            profile=profile,
        )
        check_fixed_calibration(
            cohort,
            Path(job["calibration"]),
            Path(job["calibration_audit"]),
            calibration_value,
            calibration_audit_value,
        )
    result = run(
        Path(job["folder"]),
        Path(job["audit_receipt"]),
        Path(job["audit_binding"]),
        job["binding_sha256"],
        stage,
        method,
        model,
        output,
        bootstrap=protocol["bootstrap_count"] if repeat == 0 else 0,
        seed=protocol["bootstrap_seed"],
        eps=protocol["eps"],
        rcond=protocol["rcond"],
        n_omega=protocol["n_omega"],
        n_alpha=protocol["n_alpha"],
        minimizer_mu0=protocol["minimizer_mu0"],
        precision_bits=protocol["projection_precision_bits"],
        profile=profile,
    )
    # Revalidate every bound raw receipt/archive after the calculation.
    _, cohort = load_cohort(
        Path(job["folder"]),
        Path(job["audit_receipt"]),
        Path(job["audit_binding"]),
        job["binding_sha256"],
        stage,
        profile=profile,
    )
    if snapshot(profile=profile) != frozen or sha(job_path) != job_digest:
        raise RuntimeError("worker source/runtime or job changed during execution")
    if stage == "final":
        calibration_value = verify_calibration(
            Path(job["calibration"]), job["calibration_sha256"], frozen, profile=profile
        )
        calibration_audit_value = verify_calibration_audit(
            Path(job["calibration_audit"]),
            job["calibration_audit_sha256"],
            job["calibration_sha256"],
            calibration_value["binding_sha256"],
            profile=profile,
        )
    if stage == "final":
        check_fixed_calibration(
            cohort,
            Path(job["calibration"]),
            Path(job["calibration_audit"]),
            calibration_value,
            calibration_audit_value,
        )
    if result["status"] != "passed":
        raise RuntimeError("worker case gates failed")
    validate_case_projection(
        result, output.with_suffix(".npz"), method, profile=profile
    )
    write_json(
        output.with_name(output.stem + "-worker.json"),
        {
            "schema": schema_for("-worker", profile),
            "status": "passed",
            "pid": os.getpid(),
            "job_sha256": job_digest,
            "result_sha256": sha(output),
            "array_sha256": sha(output.with_suffix(".npz")),
        },
    )


def write_progress(path, result):
    temporary = path.with_suffix(".tmp")
    write_json(temporary, result)
    os.replace(temporary, path)


def campaign(
    folder,
    receipt,
    binding,
    binding_digest,
    stage,
    output,
    *,
    calibration=None,
    calibration_digest=None,
    calibration_audit=None,
    calibration_audit_digest=None,
    profile="v2",
):
    protocol = protocol_for(profile)
    started = time.monotonic()
    validate_basis_controls(profile=profile)
    if stage not in ("calibration", "final"):
        raise ValueError("invalid campaign stage")
    frozen = snapshot(profile=profile)
    if stage == "final":
        calibration_value = verify_calibration(
            calibration, calibration_digest, frozen, profile=profile
        )
        calibration_audit_value = verify_calibration_audit(
            calibration_audit,
            calibration_audit_digest,
            calibration_digest,
            calibration_value["binding_sha256"],
            profile=profile,
        )
    _, cohort = load_cohort(
        folder, receipt, binding, binding_digest, stage, profile=profile
    )
    if stage == "final" and (
        cohort.get("role") != "fixed-bath-256"
        or cohort.get("generation_revision") != frozen["revision"]
    ):
        raise ValueError(
            "final requires fresh fixed-bath cohort from frozen source revision"
        )
    if stage == "final":
        verify_acquisition_freeze(folder, frozen)
        check_fixed_calibration(
            cohort,
            calibration,
            calibration_audit,
            calibration_value,
            calibration_audit_value,
        )
        validate_cohort_pair(calibration_value["cohort"], cohort)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError("campaign output directory must be empty")
    output = output.resolve()
    target = output / ("calibration.json" if stage == "calibration" else "final.json")
    progress_path = output / "progress.json"
    cases, records, workers = {}, {}, {}
    progress = {
        "schema": schema_for("-progress", profile),
        "status": "running",
        "stage": stage,
        "parent_pid": os.getpid(),
        "protocol": protocol,
        "snapshot": frozen,
        "binding_sha256": binding_digest,
        "cohort": cohort,
        "calibration_cohort": calibration_value["cohort"] if stage == "final" else None,
        "cases": cases,
        "workers": workers,
        "global_timeout_seconds": GLOBAL_SECONDS,
        "case_timeout_seconds": CASE_SECONDS,
    }
    write_progress(progress_path, progress)
    repeats = protocol[stage + "_repeats"]
    try:
        for repeat in range(repeats):
            for method in METHODS:
                for model in protocol["default_models"]:
                    name = f"{method}-{model}-r{repeat}.json"
                    case = output / name
                    job_path = output / (case.stem + "-job.json")
                    log_path = output / (case.stem + "-worker.log")
                    job = {
                        "schema": schema_for("-worker-job", profile),
                        "stage": stage,
                        "parent_pid": os.getpid(),
                        "method": method,
                        "model": model,
                        "repeat": repeat,
                        "protocol": protocol,
                        "snapshot": frozen,
                        "folder": str(folder.resolve()),
                        "audit_receipt": str(receipt.resolve()),
                        "audit_binding": str(binding.resolve()),
                        "binding_sha256": binding_digest,
                        "output": str(case),
                        "calibration": str(calibration.resolve())
                        if calibration
                        else None,
                        "calibration_sha256": calibration_digest,
                        "calibration_audit": str(calibration_audit.resolve())
                        if calibration_audit
                        else None,
                        "calibration_audit_sha256": calibration_audit_digest,
                    }
                    # Empty output and exclusive job creation forbid retries.
                    with job_path.open("x") as stream:
                        json.dump(
                            job, stream, sort_keys=True, indent=2, allow_nan=False
                        )
                        stream.write("\n")
                    job_digest = sha(job_path)
                    command = [
                        sys.executable,
                        "-m",
                        "experiments.cthyb_spectroscopy",
                        "--worker",
                        str(job_path),
                        "--job-sha256",
                        job_digest,
                    ]
                    if profile in ("v3", "v4"):
                        command.extend(["--profile", profile])
                    remaining = GLOBAL_SECONDS - (time.monotonic() - started)
                    if remaining <= 0:
                        raise TimeoutError("campaign global time budget exhausted")
                    workers[name] = {
                        "status": "starting",
                        "command": command,
                        "job_file": job_path.name,
                        "job_sha256": job_digest,
                        "log_file": log_path.name,
                    }
                    write_progress(progress_path, progress)

                    def mark_started(pid, case_name=name):
                        workers[case_name].update(pid=pid, status="running")
                        write_progress(progress_path, progress)

                    child = launch_worker(
                        command, log_path, min(CASE_SECONDS, remaining), mark_started
                    )
                    workers[name].update(child)
                    write_progress(progress_path, progress)
                    if child["status"] != "passed":
                        raise RuntimeError(
                            f"worker {name} {child['status']}, "
                            f"exit {child['exit_code']}"
                        )
                    if any(
                        other.get("pid") == child["pid"]
                        for key, other in workers.items()
                        if key != name
                    ):
                        raise RuntimeError("worker process ID unexpectedly reused")
                    worker_receipt = case.with_name(case.stem + "-worker.json")
                    finished = json.loads(worker_receipt.read_text())
                    if (
                        finished.get("schema") != schema_for("-worker", profile)
                        or finished.get("status") != "passed"
                        or finished.get("pid") != child["pid"]
                        or finished.get("job_sha256") != job_digest
                        or finished.get("result_sha256") != sha(case)
                        or finished.get("array_sha256") != sha(case.with_suffix(".npz"))
                    ):
                        raise ValueError("worker receipt hash or process mismatch")
                    workers[name]["worker_receipt_file"] = worker_receipt.name
                    workers[name]["worker_receipt_sha256"] = sha(worker_receipt)
                    records[name] = json.loads(case.read_text())
                    validate_case_projection(
                        records[name], case.with_suffix(".npz"), method, profile=profile
                    )
                    cases[name] = {
                        "result_sha256": sha(case),
                        "array_sha256": sha(case.with_suffix(".npz")),
                    }
                    write_progress(progress_path, progress)
    except BaseException as error:
        progress.update(status="failed", error=f"{type(error).__name__}: {error}")
        write_progress(progress_path, progress)
        raise
    try:
        if snapshot(profile=profile) != frozen:
            raise RuntimeError("scientific source or runtime changed during campaign")
        if time.monotonic() - started > GLOBAL_SECONDS:
            raise TimeoutError("campaign global time budget exhausted")
    except BaseException as error:
        progress.update(status="failed", error=f"{type(error).__name__}: {error}")
        write_progress(progress_path, progress)
        raise
    sensitivities = {}
    for repeat in range(repeats):
        for method in METHODS:
            a = records[f"{method}-flat-r{repeat}.json"]
            b = records[f"{method}-broad_gaussian-r{repeat}.json"]
            sensitivities[f"{method}-r{repeat}"] = {
                "broad_minus_flat_window_weights": (
                    np.array(b["window_weights"]) - a["window_weights"]
                ).tolist(),
                "flat_bootstrap_window_std": a["bootstrap_window_standard_deviation"],
                "broad_bootstrap_window_std": b["bootstrap_window_standard_deviation"],
            }
    result = {
        "schema": schema_for("-" + stage, profile),
        "status": "passed"
        if all(x["status"] == "passed" for x in records.values())
        else "failed",
        "protocol": protocol,
        "snapshot": frozen,
        "binding_sha256": binding_digest,
        "cohort": cohort,
        "calibration_cohort": calibration_value["cohort"] if stage == "final" else None,
        "calibration_sha256": calibration_digest if stage == "final" else None,
        "calibration_audit_sha256": calibration_audit_digest
        if stage == "final"
        else None,
        "cases": cases,
        "workers": workers,
        "parent_pid": os.getpid(),
        "elapsed_seconds": time.monotonic() - started,
        "global_timeout_seconds": GLOBAL_SECONDS,
        "case_timeout_seconds": CASE_SECONDS,
        "sensitivity": sensitivities,
        "repeat_scope": (
            "fresh_solver_invocations_same_fixed_bath_cohort_not_independent_data"
        ),
        "claim_scope": (
            "conditional_heldout_measurement_noise_not_total_prediction_"
            "uncertainty_or_spectral_truth"
        ),
    }
    write_json(target, result)
    progress["status"] = result["status"]
    progress["campaign_sha256"] = sha(target)
    write_progress(progress_path, progress)
    return result


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        parser = argparse.ArgumentParser(
            description="Execute one frozen campaign worker"
        )
        parser.add_argument("--worker", type=Path, required=True)
        parser.add_argument("--job-sha256", required=True)
        parser.add_argument("--profile", choices=["v2", "v3", "v4"], default="v2")
        args = parser.parse_args()
        worker(args.worker, args.job_sha256, profile=args.profile)
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--audit-receipt", type=Path, required=True)
    p.add_argument("--audit-binding", type=Path, required=True)
    p.add_argument("--binding-sha256", required=True)
    p.add_argument("--stage", choices=["calibration", "final"], required=True)
    p.add_argument("--calibration", type=Path)
    p.add_argument("--calibration-sha256")
    p.add_argument("--calibration-audit", type=Path)
    p.add_argument("--calibration-audit-sha256")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--profile", choices=["v2", "v3", "v4"], default="v2")
    args = p.parse_args()
    result = campaign(
        args.folder,
        args.audit_receipt,
        args.audit_binding,
        args.binding_sha256,
        args.stage,
        args.output_dir,
        calibration=args.calibration,
        calibration_digest=args.calibration_sha256,
        calibration_audit=args.calibration_audit,
        calibration_audit_digest=args.calibration_audit_sha256,
        profile=args.profile,
    )
    if result["status"] != "passed":
        raise SystemExit("campaign gates failed; artifacts retained")


if __name__ == "__main__":
    main()
