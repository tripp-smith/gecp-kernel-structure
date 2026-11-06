"""Independent audit of frozen CT-HYB spectroscopy; never runs MaxEnt."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import mpmath
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import expit

OBS = np.arange(0, 3073, 48)
QUERY = 24 + 192 * np.arange(16)
METHODS = ("raw", "kernelgecp", "pydlr", "sparse_ir")
MODELS = ("flat", "broad_gaussian")
CAL_REV = "9d28053060329cd20647dc03e038a383c351adf1"
FINAL_REV = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
CAL_SEED = 29100000
FINAL_SEED = 30100000
BINDING = "bdc256787c761d95d00ecd59b66ba3763a85fd8da4571e64306ae8c9f5ecef88"
ROOT = Path(__file__).resolve().parents[1]


def require(x, msg):
    if not x:
        raise ValueError(msg)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    x = json.loads(Path(p).read_text())
    require(isinstance(x, dict) and bool(x), "invalid JSON")
    return x


def close(a, b, msg, rtol=1e-8, atol=1e-11):
    a, b = np.asarray(a), np.asarray(b)
    require(
        a.shape == b.shape
        and np.isfinite(a).all()
        and np.isfinite(b).all()
        and np.allclose(a, b, rtol=rtol, atol=atol),
        msg,
    )


def moments(a):
    a = np.asarray(a, float)
    require(a.ndim == 2 and len(a) > 1 and np.isfinite(a).all(), "sample matrix")
    # Separate two-pass accumulation, without calling producer covariance code.
    mu = np.sum(a, axis=0) / len(a)
    d = a - mu
    return mu, np.einsum("ni,nj->ij", d, d) / (len(a) * (len(a) - 1))


def modes(c):
    c = np.asarray(c)
    require(
        c.ndim == 2 and c.shape[0] == c.shape[1] and np.isfinite(c).all(), "covariance"
    )
    scale = np.max(np.abs(c))
    slack = 128 * np.finfo(float).eps * len(c) * scale
    require(np.max(np.abs(c - c.T)) <= slack, "covariance asymmetry")
    e, v = np.linalg.eigh((c + c.T) / 2)
    require(e[0] >= -slack, "negative covariance")
    keep = e > max(1e-18, 1e-4 * e[-1])
    require(np.any(keep), "no covariance modes")
    return e[keep], v[:, keep]


def weights(w, d):
    cum = np.r_[0, np.cumsum(np.diff(w) * (d[1:] + d[:-1]) / 2)]
    return np.diff(np.interp([-4, -1, 1, 4], w, cum))


def forward(t, w, d):
    t = np.asarray(t)

    def integrand(x):
        k = -np.exp(np.where(x >= 0, -t * x, (10 - t) * x)) * expit(10 * abs(x))
        return k * np.interp(x, w, d)

    return quad_vec(integrand, -4, 4, points=w[1:-1], epsabs=1e-11, epsrel=1e-11)[0]


PARENT_BINDING = "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
FIXED_SEEDS = [33100001 + 2 * i for i in range(256)]


def fixed_raw_cohort(folder, receipt_path, binding_path, *, profile="v2"):
    fixed_seeds = [
        (33100001 if profile == "v2" else 34100001) + 2 * i for i in range(256)
    ]
    """Independently recompute fixed-bath chain moments."""
    folder = Path(folder)
    receipt_path = Path(receipt_path)
    binding_path = Path(binding_path)
    binding_sha256 = sha(binding_path)
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
        require(
            entry
            == {
                "index": i,
                "seed": fixed_seeds[i],
                "role": role,
                "receipt_sha256": sha(stem.with_suffix(".json")),
                "arrays_sha256": sha(stem.with_suffix(".npz")),
                "timing_sha256": sha(stem.with_suffix(".timing.json")),
            },
            "chain entry",
        )
        d = read(stem.with_suffix(".json"))
        require(
            d.get("schema") == "kernelgecp.cthyb-measurement.v1"
            and d.get("identity")
            == {"manifest_sha256": manifest_hash, "index": i, "role": role}
            and d.get("config") == protocol.get("config")
            and d.get("seed") == fixed_seeds[i]
            and d.get("arrays_sha256") == sha(stem.with_suffix(".npz"))
            and d.get("cycles") == 1600000
            and d.get("solve_status") == 0
            and 0.99 < d.get("average_sign", 0) <= 1.000000001,
            "chain identity",
        )
        with np.load(stem.with_suffix(".npz"), allow_pickle=False) as z:
            x, t = z["g_tau"], z["bath_iw"]
            require(
                x.shape == (3073, 2)
                and np.isfinite(x).all()
                and np.max(np.abs(x.imag)) <= 1e-10,
                "chain data",
            )
            if bath is None:
                bath = t.copy()
            require(
                np.array_equal(bath, t)
                and hashlib.sha256(t.tobytes()).hexdigest()
                == manifest["parent"].get("bath_bytes_sha256"),
                "bath",
            )
            samples.append((x[:, 0].real + x[:, 1].real) / 2)
    samples = np.array(samples)
    om, oc = moments(samples[:128, OBS])
    qm, qc = moments(samples[128:, QUERY])
    cohort = {
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
        "imaginary_absolute_tolerance": 1e-10,
    }
    return (samples, om, oc, qm, qc), cohort


def audit_acquisition_freeze(folder, frozen):
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


def raw_cohort(root, audit, binding, stage):
    root = Path(root)
    require(stage in ("calibration", "final"), "EE stage")
    binding_hash = sha(binding)
    if stage == "calibration":
        require(binding_hash == BINDING, "EE binding hash")
    revision = CAL_REV if stage == "calibration" else FINAL_REV
    seed_base = CAL_SEED if stage == "calibration" else FINAL_SEED
    b = read(binding)
    a = read(audit)
    manifest = read(root / "manifest.json")
    result = read(root / "result.json")
    require(
        b.get("schema") == "kernelgecp.cthyb-audit-binding.v1"
        and b.get("status") == "passed"
        and b.get("generation_revision") == revision
        and b.get("audit_receipt_sha256") == sha(audit),
        "EE binding",
    )
    require(
        a.get("schema") == "kernelgecp.cthyb-audit.v1"
        and a.get("status") == "passed"
        and a.get("generation_revision") == revision
        and a.get("gates")
        and all(v is True for v in a["gates"].values()),
        "EE independent audit",
    )
    audit_blob = subprocess.check_output(
        ["git", "show", a["audit_revision"] + ":experiments/audit_cthyb_dmft.py"],
        cwd=ROOT,
    )
    require(
        hashlib.sha256(audit_blob).hexdigest() == a.get("audit_sha256"),
        "EE auditor source",
    )
    for key, name in [
        ("manifest_sha256", "manifest.json"),
        ("result_sha256", "result.json"),
        ("freeze_sha256", "freeze.json"),
    ]:
        require(b.get(key) == sha(root / name), "EE " + name)
    require(
        manifest.get("source_revision") == revision
        and result.get("identity") == manifest
        and result.get("status") == "passed"
        and result.get("freeze_sha256") == sha(root / "freeze.json"),
        "EE provenance",
    )
    cfg = manifest["config"]
    require(
        cfg.get("seed_base") == seed_base
        and cfg.get("beta") == 10
        and cfg.get("interaction") == 2
        and cfg.get("seed_mapping") == "odd_stride_2_v1"
        and cfg.get("n_tau") == 3073
        and cfg.get("final_replicas") == 64,
        "EE configuration",
    )
    role, offset = (
        ("kernelgecp-final", 40000)
        if stage == "calibration"
        else ("dense-final", 30000)
    )
    samples = []
    fingerprints = set()
    bath = None
    for i in range(64):
        p = root / role / f"replica-{i:02d}.json"
        rel = p.relative_to(root).as_posix()
        d = read(p)
        require(
            b["measurement_archives"].get(rel)
            == sha(p)
            == result["measurement_archives"].get(rel),
            "chain receipt hash",
        )
        require(
            d.get("seed") == cfg["seed_base"] + offset + 1 + 2 * i
            and d.get("cycles") == 1600000
            and d.get("solve_status") == 0
            and d.get("identity") == manifest
            and d.get("average_sign", 0) > 0.99,
            "chain seed/cycles/sign",
        )
        require(
            d.get("arrays_sha256") == sha(p.with_suffix(".npz")), "chain array hash"
        )
        with np.load(p.with_suffix(".npz"), allow_pickle=False) as z:
            x = z["g_tau"]
            t = z["bath_iw"]
            require(
                x.shape == (3073, 2)
                and np.isfinite(x).all()
                and np.max(np.abs(x.imag)) <= 1e-10,
                "spin data",
            )
            if bath is None:
                bath = t.copy()
            else:
                require(np.array_equal(bath, t), "cohort bath")
            h = hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
            require(h not in fingerprints, "duplicate chain data")
            fingerprints.add(h)
            samples.append((x[:, 0].real + x[:, 1].real) / 2)
    samples = np.array(samples)
    om, oc = moments(samples[:32, OBS])
    qm, qc = moments(samples[32:, QUERY])
    cohort = {
        "stage": stage,
        "role": role,
        "generation_revision": revision,
        "seed_base": seed_base,
        "manifest_sha256": sha(root / "manifest.json"),
        "freeze_sha256": sha(root / "freeze.json"),
        "seeds": [seed_base + offset + 1 + 2 * i for i in range(64)],
        "receipt_sha256": sha(audit),
        "binding_sha256": binding_hash,
        "result_sha256": sha(root / "result.json"),
        "measurement_receipt_sha256": [
            sha(root / role / f"replica-{i:02d}.json") for i in range(64)
        ],
        "spin_operation": "paramagnetic_arithmetic_mean_up_down",
        "imaginary_absolute_tolerance": 1e-10,
    }
    return (samples, om, oc, qm, qc), cohort


def audit_svd_record(record, design):
    """Check recorded binary64 diagnostics without using its vectors as truth."""
    require(isinstance(record, dict), "projection SVD record")
    singular = np.linalg.svd(design, compute_uv=False)
    threshold = 1e-12 * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    require(
        design.shape[0] >= design.shape[1]
        and rank == design.shape[1]
        and singular[-1]
        > max(2 * threshold, 8 * np.finfo(float).eps * max(design.shape) * singular[0]),
        "projection SVD rank guard",
    )
    require(
        record.get("design_shape") == list(design.shape)
        and record.get("design_dtype") == "float64"
        and record.get("effective_rank") == rank
        and record.get("rcond") == 1e-12
        and record.get("conditioning_status") == "full_rank"
        and record.get("factorization") == "cached_thin_svd"
        and record.get("arithmetic") == "numpy_float64_or_complex128"
        and record.get("precision_bits") == 53
        and record.get("converged") is True
        and record.get("stop_reason") == "svd_factorization_completed"
        and record.get("convergence_scope")
        == "numerical_factorization_of_supplied_design"
        and record.get("solve_scope")
        == "truncated_svd_minimum_norm_fit_on_supplied_samples"
        and record.get("residual_scope")
        == "design_times_coefficients_minus_supplied_training_values"
        and record.get("certification_status") == "not_certified"
        and record.get("certification_scope") == "no_analytic_or_rounding_certificate",
        "projection SVD diagnostics",
    )
    close(
        record.get("singular_values"),
        singular,
        "projection singular values",
        rtol=1e-5,
        atol=0,
    )
    close(
        record.get("rank_threshold"),
        threshold,
        "projection rank threshold",
        rtol=1e-10,
        atol=0,
    )
    condition = singular[0] / singular[-1]
    for key in ("condition_number", "retained_condition_number"):
        close(record.get(key), condition, "projection condition", rtol=1e-5, atol=0)


def audit_projection(record, arrays, method):
    design, q, p = (
        arrays[k] for k in ("observation_design", "decomposed_design", "projection")
    )
    require(
        design.ndim == 2 and design.shape[0] == 65 and 0 < design.shape[1] <= 65,
        "design",
    )
    require(
        q.shape == design.shape
        and p.shape == (65, 65)
        and all(
            x.dtype == np.dtype("float64") and np.isfinite(x).all()
            for x in (design, q, p)
        ),
        "decomposed design shape or dtype",
    )
    info = record.get("projection_arithmetic")
    require(isinstance(info, dict), "projection arithmetic")
    for label, value in (("source_design", design), ("decomposed_design", q)):
        digest = hashlib.sha256(
            np.ascontiguousarray(value, dtype="<f8").tobytes()
        ).hexdigest()
        require(
            info.get(label + "_sha256") == digest
            and info.get(label + "_shape") == list(value.shape),
            "projection design hash or shape",
        )
    require(
        info.get("hash_scope")
        == "C-order little-endian binary64 array bytes with recorded shapes",
        "projection hash scope",
    )
    if method == "raw":
        require(
            record.get("precision_bits") == 53
            and info.get("precision_bits") == 53
            and info.get("method") == "identity"
            and info.get("rank_guard") is None
            and info.get("transformed_svd") is None
            and record.get("projection_svd") is None
            and record.get("transformed_projection_svd") is None
            and np.array_equal(design, np.eye(65))
            and np.array_equal(q, np.eye(65))
            and np.array_equal(p, np.eye(65)),
            "raw identity projection",
        )
        return
    require(
        record.get("precision_bits") == 256
        and info.get("precision_bits") == 256
        and info.get("method") == "mpmath_qr_then_binary64_cached_svd"
        and info.get("rank_guard")
        == "smallest > max(2 * threshold, 8 * eps64 * max(shape) * largest); heuristic"
        and info.get("transformed_svd") == record.get("transformed_projection_svd"),
        "projection precision or arithmetic",
    )
    audit_svd_record(record.get("projection_svd"), design)
    audit_svd_record(record.get("transformed_projection_svd"), q)
    close(q.T @ q, np.eye(q.shape[1]), "orthogonal coordinates", rtol=0, atol=2e-12)
    close(p, q @ q.T, "decomposition projection", rtol=0, atol=2e-12)
    ctx = mpmath.mp.clone()
    ctx.dps = 80
    matrix = ctx.matrix([[ctx.mpf(float(x)) for x in row] for row in design])
    u, singular, _ = ctx.svd(matrix, full_matrices=False)
    rank = sum(
        singular[i] > ctx.mpf("1e-12") * singular[0] for i in range(len(singular))
    )
    require(rank == design.shape[1], "independent projection rank")
    oracle = np.array(
        [
            [
                float(sum(u[i, k] * u[j, k] for k in range(rank)))
                for j in range(len(design))
            ]
            for i in range(len(design))
        ]
    )
    close(p, oracle, "independent projection", atol=2e-8, rtol=1e-8)


def audit_curvature_metric(parameter, config):
    metric = parameter.get("curvature_metric")
    if config.get("mass_curvature_floor_relative") != 1e-10:
        require(metric is None, "unexpected curvature metric")
        return
    require(
        isinstance(metric, dict)
        and metric.get("certified") is False
        and type(metric.get("evaluation_count")) is int
        and metric["evaluation_count"] > 0,
        "curvature metric evaluation",
    )
    value = metric.get("last_evaluation")
    require(
        isinstance(value, dict)
        and value.get("certified") is False
        and value.get("floor_relative") == 1e-10
        and type(value.get("dimension")) is int
        and 1 <= value["dimension"] <= 128,
        "curvature metric identity",
    )
    keys = (
        "scale",
        "normalized_min_eigenvalue",
        "normalized_max_eigenvalue",
        "normalized_shift",
        "shift",
        "shifted_normalized_min_eigenvalue",
        "shifted_normalized_max_eigenvalue",
        "condition_estimate",
    )
    require(
        all(type(value.get(k)) in (int, float) and np.isfinite(value[k]) for k in keys),
        "curvature metric finite estimates",
    )
    lo, hi = value["normalized_min_eigenvalue"], value["normalized_max_eigenvalue"]
    shifted_lo = value["shifted_normalized_min_eigenvalue"]
    shifted_hi = value["shifted_normalized_max_eigenvalue"]
    shift = value["normalized_shift"]
    require(
        value["scale"] > 0
        and lo <= hi
        and shift >= 0
        and np.isclose(
            shift, max(0.0, 1e-10 * value["dimension"] - lo), rtol=1e-12, atol=1e-15
        )
        and np.isclose(value["shift"], shift * value["scale"], rtol=1e-12, atol=0)
        and 0 < shifted_lo <= shifted_hi
        and np.isclose(
            value["condition_estimate"], shifted_hi / shifted_lo, rtol=1e-12
        ),
        "curvature metric inconsistent estimates",
    )


def audit_v3_native_metadata(diagnostics, config):
    param = diagnostics.get("mass_parameterization", {})
    require(
        diagnostics.get("mass_operation") == "native_mass_constraint"
        and param.get("curvature") == "projected"
        and param.get("interior_fraction") == 1e-12
        and param.get("d_dv") is False
        and param.get("dA_projection") == 2
        and param.get("chi2_factor") == 1.0
        and param.get("floor_policy") == "reject_native_entropy_floor_or_overflow"
        and diagnostics.get("config", {}) == config
        and diagnostics.get("provider") == "TRIQS/MaxEnt"
        and diagnostics.get("version") == "4.0.1"
        and diagnostics.get("native_minimizer", {}).get("class") == "LevenbergMinimizer"
        and diagnostics.get("native_minimizer", {}).get("marquardt") is True
        and diagnostics.get("native_minimizer", {}).get("J_squared") is False
        and diagnostics.get("native_minimizer", {}).get("mu0") == 1e-8
        and diagnostics.get("native_minimizer", {}).get("maxiter") == 1000,
        "native mass and damping metadata",
    )
    audit_curvature_metric(param, config or {})
    for key in ("min_observed_H", "min_observed_mixture_floor"):
        value = param.get(key)
        require(
            type(value) in (float, int) and np.isfinite(value) and value > 1e-100,
            "native mass floor metadata",
        )


def audit_case(path, method, model, repeat, stats, *, profile="v2"):
    require(profile in ("v2", "v3", "v4"), "profile")
    d = read(path)
    npz = path.with_suffix(".npz")
    require(
        d.get("schema") == f"kernelgecp.cthyb-spectroscopy.{profile}"
        and d.get("status") == "passed"
        and d.get("method") == method
        and d.get("config", {}).get("default_model") == model
        and d.get("array_sha256") == sha(npz),
        "case identity",
    )
    require(
        d.get("covariance_rtol") == 1e-4
        and d.get("observation_indices") == OBS.tolist()
        and d.get("query_indices") == QUERY.tolist(),
        "partition",
    )
    require(
        d.get("eps") == 1e-10
        and d.get("rcond") == 1e-12
        and d.get("config")
        == {
            "omega_min": -4.0,
            "omega_max": 4.0,
            "n_omega": 401,
            "alpha_min": 1e-4,
            "alpha_max": 1e4,
            "n_alpha": 33,
            "maxiter": 1000,
            "default_model": model,
            "normalize_mass": 1.0,
            "minimizer_mu0": 1e-8,
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
        },
        "solver controls",
    )
    count = 8 if repeat == 0 else 0
    require(
        d.get("bootstrap_count") == count and d.get("bootstrap_seed") == 193,
        "bootstrap count/seed",
    )
    with np.load(npz, allow_pickle=False) as z:
        a = {k: z[k] for k in z.files}
    required = {
        "omega",
        "density",
        "observation_tau",
        "observation_values",
        "covariance",
        "projected_values",
        "projected_covariance",
        "projection",
        "observation_design",
        "decomposed_design",
        "native_kernel",
        "projected_kernel",
        "query",
        "query_mean",
        "validation_covariance",
        "prediction",
        "alphas",
        "chi2",
        "trial_spectra",
        "bootstrap_weights",
        "bootstrap_spectra",
        "bootstrap_indices",
    }
    require(
        set(a) == required and all(np.isfinite(x).all() for x in a.values()),
        "array inventory",
    )
    samples, om, oc, qm, qc = stats
    t = np.linspace(0, 10, 3073)
    close(a["observation_tau"], t[OBS], "obs tau")
    close(a["query"], t[QUERY], "query tau")
    require(not np.intersect1d(OBS, QUERY).size, "overlap")
    close(a["observation_values"], om, "train mean")
    close(a["covariance"], oc, "train covariance", atol=1e-18)
    close(a["query_mean"], qm, "heldout mean")
    close(a["validation_covariance"], qc, "heldout covariance", atol=1e-18)
    p = a["projection"]
    require(p.shape == (65, 65), "projection shape")
    close(p, p.T, "projection symmetry", atol=2e-10)
    close(p @ p, p, "projection idempotence", atol=2e-10)
    audit_projection(d, a, method)
    native_times = t[OBS]
    w0 = np.linspace(-4, 4, 401)
    widths = np.r_[np.diff(w0)[0] / 2, (w0[2:] - w0[:-2]) / 2, np.diff(w0)[-1] / 2]
    native = (
        -np.exp(
            np.where(
                w0[None, :] >= 0,
                -native_times[:, None] * w0[None, :],
                (10 - native_times[:, None]) * w0[None, :],
            )
        )
        * expit(10 * np.abs(w0))[None, :]
    )
    close(a["native_kernel"], native, "native kernel", rtol=1e-10, atol=1e-12)
    close(a["projected_values"], p @ om, "projected values")
    close(a["projected_covariance"], p @ oc @ p.T, "projected covariance", atol=2e-18)
    close(a["projected_kernel"], p @ a["native_kernel"], "projected kernel")
    for key, cov in [
        ("training_covariance_modes", oc),
        ("projected_covariance_modes", p @ oc @ p.T),
        ("validation_covariance_modes", qc),
    ]:
        eig = np.linalg.eigvalsh((cov + cov.T) / 2)
        record = d[key]
        close(
            record["eigenvalues"],
            eig,
            "covariance eigenspectrum",
            rtol=1e-6,
            atol=1e-18,
        )
        close(
            record["threshold"],
            max(1e-18, 1e-4 * eig[-1]),
            "covariance threshold",
            rtol=1e-6,
            atol=1e-18,
        )
        require(
            record["relative_cutoff"] == 1e-4
            and record["absolute_cutoff"] == 1e-18
            and record["effective_rank"] == len(modes(cov)[0]),
            "covariance mode record",
        )
    e, v = modes(p @ oc @ p.T)
    qe, qv = modes(qc)
    require(
        d["training_rank"] == len(modes(oc)[0])
        and d["projected_rank"] == len(e)
        and d["validation_rank"] == len(qe),
        "covariance rank",
    )
    w = a["omega"]
    den = a["density"]
    close(w, np.linspace(-4, 4, 401), "omega")
    require(
        den.shape == (401,)
        and np.min(den) >= -1e-6
        and abs(np.trapezoid(den, w) - 1) <= 1e-6,
        "spectrum positivity or mass",
    )
    trials = a["trial_spectra"]
    chi = a["chi2"]
    alph = a["alphas"]
    conv = d["converged"]
    require(
        trials.shape == (33, 401)
        and chi.shape == (33,)
        and len(conv) == len(d["iterations"]) == 33
        and all(type(x) is bool for x in conv),
        "trials",
    )
    close(alph, np.geomspace(1e4, 1e-4, 33), "alphas")
    calc = np.sum(
        (
            (
                v.T
                @ (
                    p @ a["native_kernel"] @ (widths[:, None] * trials.T)
                    - a["projected_values"][:, None]
                )
            )
            / np.sqrt(e)[:, None]
        )
        ** 2,
        axis=0,
    )
    close(chi, calc, "native chi2", rtol=3e-5, atol=2e-5)
    masses = np.trapezoid(trials, w, axis=1)
    eligible = (
        np.array(conv)
        & (chi >= 0)
        & np.all(trials >= 0, axis=1)
        & np.isfinite(masses)
        & (masses > 0)
    )
    require(eligible.any(), "convergence")
    # Integrate all eligible piecewise-linear spectra with adaptive SciPy quadrature.
    indices = np.flatnonzero(eligible)
    if profile in ("v3", "v4"):
        require(
            np.all(trials > 0) and np.all(np.abs(masses - 1.0) <= 1e-10),
            "native trial masses or positivity",
        )
        normalized = trials[indices]
    else:
        normalized = trials[indices] / masses[indices, None]

    def integrand(x):
        j = min(np.searchsorted(w, x, side="right") - 1, len(w) - 2)
        j = max(j, 0)
        f = (x - w[j]) / (w[j + 1] - w[j])
        density = normalized[:, j] * (1 - f) + normalized[:, j + 1] * f
        kernel = -np.exp(
            np.where(x >= 0, -native_times * x, (10 - native_times) * x)
        ) * expit(10 * abs(x))
        return kernel[:, None] * density[None, :]

    physical = quad_vec(integrand, -4, 4, points=w[1:-1], epsabs=1e-11, epsrel=1e-11)[0]
    residual = p @ physical - a["projected_values"][:, None]
    selection = np.sum((v.T @ residual / np.sqrt(e)[:, None]) ** 2, axis=0)
    diagnostics = d["provider_diagnostics"]
    require(
        diagnostics.get("alpha_selection")
        == "normalized_physical_forward_training_discrepancy"
        and diagnostics.get("selection_rule")
        == "closest_abs_log_chi2_per_retained_rank_to_one",
        "selection rule",
    )
    close(diagnostics.get("trial_masses"), masses, "trial masses")
    reported = diagnostics.get("selection_chi2")
    require(
        isinstance(reported, list)
        and len(reported) == 33
        and all(reported[i] is None for i in range(33) if not eligible[i]),
        "selection chi2 inventory",
    )
    close(
        [reported[i] for i in indices],
        selection,
        "selection chi2",
        rtol=3e-5,
        atol=2e-5,
    )
    score = np.abs(np.log(np.maximum(selection, np.finfo(float).tiny) / len(e)))
    idx = int(indices[np.argmin(score)])
    require(d["selected_index"] == idx, "alpha selection")
    raw = masses[idx]
    close(
        den,
        trials[idx] if profile in ("v3", "v4") else trials[idx] / raw,
        "selected spectrum",
    )
    if profile in ("v3", "v4"):
        require(np.array_equal(den, trials[idx]), "post-fit constrained rescale")
        audit_v3_native_metadata(diagnostics, d["config"])
    prediction = forward(a["query"], w, den)
    close(prediction, a["prediction"], "independent quadrature", rtol=2e-6, atol=2e-7)
    rms = float(
        np.linalg.norm((qv.T @ (prediction - qm)) / np.sqrt(qe)) / np.sqrt(len(qe))
    )
    require(rms <= 2, "heldout whitened RMS")
    close(d["whitened_rms"], rms, "reported RMS", rtol=2e-4, atol=2e-4)
    close(d["window_weights"], weights(w, den), "window weights")
    require(
        a["bootstrap_spectra"].shape == ((8, 401) if count else (0,))
        and a["bootstrap_indices"].shape == ((8, len(samples) // 2) if count else (0,)),
        "bootstrap shape",
    )
    if count:
        ix = np.random.default_rng(d["bootstrap_seed"]).integers(
            0, len(samples) // 2, size=(8, len(samples) // 2)
        )
        require(np.array_equal(ix, a["bootstrap_indices"]), "bootstrap indices")
        bw = []
        for i, x in enumerate(a["bootstrap_spectra"]):
            if profile in ("v3", "v4"):
                rec = d["bootstrap_records"][i]
                require(
                    abs(rec["raw_mass"] - 1) <= 1e-10 and np.all(x > 0),
                    "bootstrap native mass",
                )
                audit_v3_native_metadata(
                    rec.get("provider_diagnostics", {}), d["config"]
                )
            require(
                np.min(x) >= -1e-6 and abs(np.trapezoid(x, w) - 1) <= 1e-6,
                "bootstrap spectrum",
            )
            bw.append(weights(w, x))
        close(a["bootstrap_weights"], bw, "bootstrap weights")
        close(
            d["bootstrap_window_standard_deviation"],
            np.std(bw, axis=0, ddof=1),
            "bootstrap variation",
        )
        require(
            len(d["bootstrap_records"]) == 8
            and all(
                x["converged"][x["selected_index"]] for x in d["bootstrap_records"]
            ),
            "bootstrap convergence",
        )
    return d, a


def check_snapshot(snapshot, profile="v2"):
    require(
        isinstance(snapshot, dict)
        and snapshot.get("revision")
        and snapshot.get("sources")
        and snapshot.get("provider_modules")
        and snapshot.get("packages", {}).get("triqs_maxent") == "4.0.1",
        "runtime snapshot",
    )
    require(
        all(v == "1" for v in snapshot.get("threads", {}).values())
        and len(snapshot.get("threads", {})) >= 5,
        "thread snapshot",
    )
    rev = snapshot["revision"]
    require(
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", rev, "HEAD"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0,
        "snapshot revision",
    )
    tracked = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", rev, "kernelgecp"], cwd=ROOT, text=True
    ).splitlines()
    required = {p for p in tracked if p.endswith(".py")}
    required |= {
        "experiments/cthyb_spectroscopy.py",
        "experiments/audit_cthyb_spectroscopy.py",
        "experiments/cthyb_dmft.py",
        "experiments/spectroscopy_campaign.py",
        "experiments/audit_spectroscopy.py",
        "experiments/spectroscopy.py",
        "experiments/dlr_adapters.py",
        "environments/locks/triqs-linux-64-conda-explicit.txt",
        "environments/locks/triqs-maxent-pip.txt",
        "environments/locks/dlr-physics-pip.txt",
    }
    if profile == "v4":
        required.add("kernelgecp/integrations/_maxent_curvature.py")
    require(
        required.issubset(snapshot["sources"]), "incomplete scientific source inventory"
    )
    for name, h in snapshot["sources"].items():
        require(
            not Path(name).is_absolute() and ".." not in Path(name).parts, "source path"
        )
        blob = subprocess.check_output(["git", "show", rev + ":" + name], cwd=ROOT)
        require(hashlib.sha256(blob).hexdigest() == h, "source hash " + name)
    require(
        "experiments/cthyb_spectroscopy.py" in snapshot["sources"]
        and "experiments/audit_cthyb_spectroscopy.py" in snapshot["sources"],
        "source inventory",
    )


def check_workers(
    manifest,
    campaign,
    root,
    ee_audit,
    binding,
    stage,
    calibration,
    calibration_audit,
    *,
    profile="v2",
):
    """Verify the complete pinned child-job and process receipt inventory."""
    workers = manifest.get("workers")
    require(
        isinstance(workers, dict) and set(workers) == set(manifest["cases"]),
        "worker inventory",
    )
    require(
        manifest.get("global_timeout_seconds") == 19800
        and manifest.get("case_timeout_seconds") == 3600,
        "worker deadline controls",
    )
    elapsed_total = manifest.get("elapsed_seconds")
    require(
        type(elapsed_total) in (int, float)
        and np.isfinite(elapsed_total)
        and 0 < elapsed_total <= 19800,
        "campaign duration",
    )
    seen = set()
    parent_pid = manifest.get("parent_pid")
    require(type(parent_pid) is int and parent_pid > 0, "campaign parent PID")
    campaign_dir = Path(campaign).resolve().parent
    for name, record in workers.items():
        method, model, repeat = name[:-5].rsplit("-", 2)
        require(
            record.get("status") == "passed" and record.get("exit_code") == 0,
            "worker exit",
        )
        pid = record.get("pid")
        require(type(pid) is int and pid > 0 and pid not in seen, "worker PID")
        seen.add(pid)
        elapsed = record.get("elapsed_seconds")
        require(
            type(elapsed) in (int, float)
            and np.isfinite(elapsed)
            and 0 < elapsed <= 3600,
            "worker duration",
        )

        def local(key, current_record=record):
            p = Path(current_record[key])
            p = p if p.is_absolute() else campaign_dir / p
            p = p.resolve()
            require(p.is_relative_to(campaign_dir) and p.is_file(), "worker file path")
            return p

        job_path, log_path, worker_path = (
            local(k) for k in ("job_file", "log_file", "worker_receipt_file")
        )
        require(
            sha(job_path) == record.get("job_sha256")
            and sha(log_path) == record.get("log_sha256")
            and sha(worker_path) == record.get("worker_receipt_sha256"),
            "worker artifact hash",
        )
        job = read(job_path)
        parent = job.get("parent_pid")
        require(
            type(parent) is int and parent > 0 and parent != pid, "worker parent PID"
        )
        require(parent == parent_pid, "worker parent mismatch")
        require(
            job.get("schema") == f"kernelgecp.cthyb-spectroscopy-worker-job.{profile}"
            and job.get("stage") == stage
            and job.get("method") == method
            and job.get("model") == model
            and job.get("repeat") == int(repeat[1:])
            and job.get("protocol") == manifest["protocol"]
            and job.get("snapshot") == manifest["snapshot"],
            "worker job controls",
        )
        require(
            Path(job["folder"]).resolve() == Path(root).resolve()
            and Path(job["audit_receipt"]).resolve() == Path(ee_audit).resolve()
            and Path(job["audit_binding"]).resolve() == Path(binding).resolve()
            and job.get("binding_sha256") == sha(binding)
            and Path(job["output"]).resolve() == campaign_dir / name,
            "worker input paths",
        )
        if stage == "final":
            require(
                Path(job["calibration"]).resolve() == Path(calibration).resolve()
                and job.get("calibration_sha256") == sha(calibration)
                and Path(job["calibration_audit"]).resolve()
                == Path(calibration_audit).resolve()
                and job.get("calibration_audit_sha256") == sha(calibration_audit),
                "worker calibration binding",
            )
        else:
            require(
                all(
                    job.get(k) is None
                    for k in (
                        "calibration",
                        "calibration_sha256",
                        "calibration_audit",
                        "calibration_audit_sha256",
                    )
                ),
                "calibration worker controls",
            )
        command = record.get("command")
        require(
            isinstance(command, list)
            and len(command) == (7 if profile == "v2" else 9)
            and Path(command[0]).resolve()
            == Path(manifest["snapshot"]["executable"]).resolve()
            and command[1:4] == ["-m", "experiments.cthyb_spectroscopy", "--worker"]
            and Path(command[4]).resolve() == job_path
            and command[5:]
            == ["--job-sha256", sha(job_path)]
            + ([] if profile == "v2" else ["--profile", profile]),
            "worker command",
        )
        receipt = read(worker_path)
        require(
            receipt.get("schema") == f"kernelgecp.cthyb-spectroscopy-worker.{profile}"
            and receipt.get("status") == "passed"
            and receipt.get("pid") == pid
            and receipt.get("job_sha256") == sha(job_path)
            and receipt.get("result_sha256") == manifest["cases"][name]["result_sha256"]
            and receipt.get("array_sha256") == manifest["cases"][name]["array_sha256"],
            "worker receipt",
        )
    require(len(seen) == len(workers), "worker PID inventory")


def check_calibration_binding(cal, ca, calibration, *, profile="v2"):
    require(
        ca.get("schema") == f"kernelgecp.cthyb-spectroscopy-independent-audit.{profile}"
        and ca.get("status") == "passed"
        and ca.get("kind") == "calibration"
        and ca.get("calibration_sha256") == sha(calibration)
        and ca.get("binding_sha256") == cal.get("binding_sha256")
        and cal.get("binding_sha256") == BINDING
        and ca.get("gates")
        and all(v is True for v in ca["gates"].values()),
        "calibration audit binding",
    )


def campaign_audit(
    campaign,
    root,
    ee_audit,
    binding,
    calibration=None,
    calibration_audit=None,
    *,
    profile="v2",
):
    require(profile in ("v2", "v3", "v4"), "profile")
    campaign = Path(campaign)
    m = read(campaign)
    stage = (
        "calibration"
        if m.get("schema") == f"kernelgecp.cthyb-spectroscopy-calibration.{profile}"
        else "final"
    )
    require(
        m.get("schema") == f"kernelgecp.cthyb-spectroscopy-{stage}.{profile}"
        and m.get("status") == "passed"
        and m.get("binding_sha256") == sha(binding),
        "campaign schema",
    )
    protocol = m["protocol"]
    expected_protocol = {
        "profile": "corrective_projection_v2",
        "projection_precision_bits": 256,
        "methods": list(METHODS),
        "default_models": list(MODELS),
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
        "covariance_rtol": 1e-4,
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
        "calibration_generation_revision": CAL_REV,
        "parent_generation_revision": FINAL_REV,
        "final_generation_revision_policy": "equal_campaign_snapshot_revision",
        "final_seed_base": 33100001,
        "final_seed_stride": 2,
        "final_chain_count": 256,
        "final_training_count": 128,
        "final_heldout_count": 128,
    }
    if profile in ("v3", "v4"):
        expected_protocol.update(
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
        expected_protocol.update(
            profile="native_shifted_curvature_v4", mass_curvature_floor_relative=1e-10
        )
    require(protocol == expected_protocol, "protocol")
    check_snapshot(m["snapshot"], profile)
    if stage == "final":
        require(
            calibration is not None and calibration_audit is not None,
            "independent calibration audit required",
        )
        ca = read(calibration_audit)
        cal = read(calibration)
        require(m.get("calibration_sha256") == sha(calibration), "calibration hash")
        check_calibration_binding(cal, ca, calibration, profile=profile)
        audit_blob = subprocess.check_output(
            [
                "git",
                "show",
                ca["auditor_revision"] + ":experiments/audit_cthyb_spectroscopy.py",
            ],
            cwd=ROOT,
        )
        require(
            hashlib.sha256(audit_blob).hexdigest() == ca.get("auditor_sha256"),
            "calibration auditor source",
        )
        require(
            cal.get("protocol") == protocol and cal.get("snapshot") == m["snapshot"],
            "calibration freeze",
        )
    if (
        stage == "final"
        and read(binding).get("schema")
        == f"kernelgecp.cthyb-fixed-bath-audit-binding.{profile}"
    ):
        samples, cohort = fixed_raw_cohort(root, ee_audit, binding, profile=profile)
        require(
            cohort["generation_revision"] == m["snapshot"]["revision"],
            "fixed-bath source freeze",
        )
        audit_acquisition_freeze(Path(root), m["snapshot"])
        expected_calibration = {
            "receipt_sha256": sha(calibration),
            "audit_sha256": sha(calibration_audit),
            "binding_sha256": cal["binding_sha256"],
            "auditor_revision": ca["auditor_revision"],
            "auditor_sha256": ca["auditor_sha256"],
        }
        require(
            cohort.get("calibration") == expected_calibration
            and m.get("calibration_audit_sha256") == sha(calibration_audit),
            "final calibration differs from pre-acquisition freeze",
        )
    else:
        require(stage == "calibration", "v2 fixed-bath binding required")
        samples, cohort = raw_cohort(root, ee_audit, binding, stage)
    require(m.get("cohort") == cohort, "campaign cohort metadata")
    if stage == "calibration":
        require(m.get("calibration_cohort") is None, "calibration cohort metadata")
    else:
        cc = cal.get("cohort", {})
        require(
            cal.get("schema") == f"kernelgecp.cthyb-spectroscopy-calibration.{profile}"
            and cal.get("status") == "passed"
            and cc.get("stage") == "calibration"
            and cc.get("role") == "kernelgecp-final"
            and cc.get("generation_revision") == CAL_REV
            and cc.get("seed_base") == CAL_SEED
            and cc.get("binding_sha256") == BINDING
            and cc.get("seeds") == [CAL_SEED + 40001 + 2 * i for i in range(64)]
            and m.get("calibration_cohort") == cc,
            "calibration cohort metadata",
        )
        require(
            cc["generation_revision"] != cohort["generation_revision"]
            and cc["binding_sha256"] != cohort["binding_sha256"]
            and not set(cc["seeds"]) & set(cohort["seeds"]),
            "cohort overlap",
        )
    repeats = 1 if stage == "calibration" else 3
    expected = {
        f"{method}-{model}-r{repeat}.json"
        for method in METHODS
        for model in MODELS
        for repeat in range(repeats)
    }
    require(set(m["cases"]) == expected, "child inventory")
    check_workers(
        m,
        campaign,
        root,
        ee_audit,
        binding,
        stage,
        calibration,
        calibration_audit,
        profile=profile,
    )
    records = {}
    arrays = {}
    for name in sorted(expected):
        method, model, r = name[:-5].rsplit("-", 2)
        repeat = int(r[1:])
        p = campaign.parent / name
        h = m["cases"][name]
        require(
            sha(p) == h["result_sha256"]
            and sha(p.with_suffix(".npz")) == h["array_sha256"],
            "child hashes",
        )
        records[name], arrays[name] = audit_case(
            p, method, model, repeat, samples, profile=profile
        )
        require(records[name].get("source") == cohort, "case cohort metadata")
    if stage == "final":
        for method in METHODS:
            for model in MODELS:
                first = arrays[f"{method}-{model}-r0.json"]
                for rep in (1, 2):
                    other = arrays[f"{method}-{model}-r{rep}.json"]
                    for k in first:
                        if not k.startswith("bootstrap_"):
                            require(
                                np.array_equal(first[k], other[k]),
                                "repeated scientific arrays " + k,
                            )
    sensitivity = {}
    for repeat in range(repeats):
        for method in METHODS:
            a = records[f"{method}-flat-r{repeat}.json"]
            b = records[f"{method}-broad_gaussian-r{repeat}.json"]
            delta = np.array(b["window_weights"]) - a["window_weights"]
            key = f"{method}-r{repeat}"
            close(
                m["sensitivity"][key]["broad_minus_flat_window_weights"],
                delta,
                "default sensitivity",
            )
            sensitivity[key] = delta.tolist()
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    source = ROOT / "experiments/audit_cthyb_spectroscopy.py"
    blob = subprocess.check_output(
        ["git", "show", revision + ":experiments/audit_cthyb_spectroscopy.py"], cwd=ROOT
    )
    require(
        sha(source) == hashlib.sha256(blob).hexdigest(),
        "auditor source must be committed",
    )
    return {
        "schema": f"kernelgecp.cthyb-spectroscopy-independent-audit.{profile}",
        "status": "passed",
        "kind": stage,
        "calibration_sha256": sha(campaign)
        if stage == "calibration"
        else sha(calibration),
        "campaign_sha256": sha(campaign),
        "binding_sha256": sha(binding),
        "auditor_revision": revision,
        "auditor_sha256": sha(source),
        "gates": {
            "ee_source_and_raw_chains": True,
            "frozen_protocol_and_hashes": True,
            "covariance_and_projection": True,
            "native_trials_and_quadrature": True,
            "heldout_rms": True,
            "bootstrap_and_defaults": True,
            "repeat_arrays": True,
            "fresh_process_receipts": True,
        },
        "cases": len(records),
        "default_sensitivity": sensitivity,
        "scope": (
            "heldout chain mean only; finite native model, no truth CDF, "
            "nonlinear coverage or unique spectrum claim"
        ),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("--ee-root", type=Path, required=True)
    p.add_argument("--ee-audit", type=Path, required=True)
    p.add_argument("--ee-binding", type=Path, required=True)
    p.add_argument("--calibration", type=Path)
    p.add_argument("--calibration-audit", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--profile", choices=["v2", "v3", "v4"], default="v2")
    a = p.parse_args()
    try:
        r = campaign_audit(
            a.campaign,
            a.ee_root,
            a.ee_audit,
            a.ee_binding,
            a.calibration,
            a.calibration_audit,
            profile=a.profile,
        )
        code = 0
    except Exception as e:
        r = {
            "schema": f"kernelgecp.cthyb-spectroscopy-independent-audit.{a.profile}",
            "status": "failed",
            "error": f"{type(e).__name__}: {e}",
        }
        code = 1
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(r, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
