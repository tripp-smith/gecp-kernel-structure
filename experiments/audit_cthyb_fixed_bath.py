"""Independent audit of the bounded EH fixed-bath chain acquisition.

A passed receipt concerns the raw cohort only; it is not an EE campaign pass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PARENT_REVISION = "7200ec0e912954c568bedb2092d9ff1c053f9ba7"
PARENT_BINDING_SHA256 = (
    "6e345920a20f700d16dbad82ad1ff240c28cf827321d3fac7e04ba7904348fcd"
)
SHAPES = {
    "bath_iw": (512,),
    "g_tau": (3073, 2),
    "g_iw": (512, 2),
    "g_iw_raw": (512, 2),
    "sigma_iw": (512, 2),
    "g_moments": (3, 2),
    "sigma_moments": (2, 2),
    "occupancies": (2,),
}
LIMIT = 64 * 1024 * 1024


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    require(
        path.is_file() and not path.is_symlink() and path.stat().st_size <= LIMIT,
        f"unsafe file: {path}",
    )
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    require(
        path.is_file() and not path.is_symlink() and path.stat().st_size <= LIMIT,
        f"unsafe JSON: {path}",
    )

    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    def bad_constant(_):
        raise ValueError("nonfinite JSON number")

    value = json.loads(
        path.read_bytes(), object_pairs_hook=pairs, parse_constant=bad_constant
    )

    def finite_values(item):
        if isinstance(item, float):
            return math.isfinite(item)
        if isinstance(item, list):
            return all(finite_values(x) for x in item)
        if isinstance(item, dict):
            return all(finite_values(x) for x in item.values())
        return True

    require(finite_values(value), "nonfinite JSON number")
    return value


def arrays(path):
    sha(path)
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        require(
            len(entries) == len(SHAPES)
            and {e.filename for e in entries} == {k + ".npy" for k in SHAPES},
            "NPZ inventory mismatch",
        )
        require(sum(e.file_size for e in entries) <= LIMIT, "NPZ allocation too large")
        for key, shape in SHAPES.items():
            entry = archive.getinfo(key + ".npy")
            with archive.open(entry) as stream:
                version = np.lib.format.read_magic(stream)
                require(version in ((1, 0), (2, 0)), "unsupported NPY version")
                if version == (1, 0):
                    actual, order, dtype = np.lib.format.read_array_header_1_0(
                        stream, max_header_size=10000
                    )
                else:
                    actual, order, dtype = np.lib.format.read_array_header_2_0(
                        stream, max_header_size=10000
                    )
                expected = np.dtype("float64" if key == "occupancies" else "complex128")
                require(
                    actual == shape and not order and dtype == expected,
                    "NPY header mismatch",
                )
                require(
                    stream.tell() + math.prod(shape) * dtype.itemsize
                    == entry.file_size,
                    "NPY size mismatch",
                )
    with np.load(path, allow_pickle=False) as data:
        result = {key: data[key] for key in SHAPES}
    require(all(np.isfinite(a).all() for a in result.values()), "nonfinite array")
    return result


def git_blob(revision, path):
    require(
        isinstance(revision, str)
        and len(revision) == 40
        and all(c in "0123456789abcdef" for c in revision),
        "invalid revision",
    )
    require(
        isinstance(path, str)
        and path
        and not path.startswith("/")
        and ".." not in Path(path).parts,
        "invalid source path",
    )
    return subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=ROOT)


CAL_BINDING = "bdc256787c761d95d00ecd59b66ba3763a85fd8da4571e64306ae8c9f5ecef88"
CAL_REVISION = "9d28053060329cd20647dc03e038a383c351adf1"
CAL_GATES = {
    "ee_source_and_raw_chains",
    "frozen_protocol_and_hashes",
    "covariance_and_projection",
    "native_trials_and_quadrature",
    "heldout_rms",
    "bootstrap_and_defaults",
    "repeat_arrays",
    "fresh_process_receipts",
}


def validate_calibration(calibration_path, audit_path, ref, provenance):
    """Independently check the frozen calibration and its independent audit."""
    calibration_path, audit_path = Path(calibration_path), Path(audit_path)
    cal_hash, audit_hash = sha(calibration_path), sha(audit_path)
    cal, ca = read_json(calibration_path), read_json(audit_path)
    require(
        ref
        == {
            "receipt_sha256": cal_hash,
            "audit_sha256": audit_hash,
            "binding_sha256": CAL_BINDING,
            "auditor_revision": ca.get("auditor_revision"),
            "auditor_sha256": ca.get("auditor_sha256"),
        },
        "calibration freeze hash mismatch",
    )
    require(
        cal.get("schema") == "kernelgecp.cthyb-spectroscopy-calibration.v1"
        and cal.get("status") == "passed"
        and cal.get("binding_sha256") == CAL_BINDING,
        "calibration receipt did not pass",
    )
    gates = ca.get("gates")
    require(
        ca.get("schema") == "kernelgecp.cthyb-spectroscopy-independent-audit.v1"
        and ca.get("status") == "passed"
        and ca.get("kind") == "calibration"
        and ca.get("cases") == 8
        and ca.get("campaign_sha256") == cal_hash
        and ca.get("calibration_sha256") == cal_hash
        and ca.get("binding_sha256") == CAL_BINDING
        and isinstance(gates, dict)
        and set(gates) == CAL_GATES
        and all(v is True for v in gates.values()),
        "calibration independent audit did not pass",
    )
    require(
        ca["auditor_revision"] == provenance["revision"]
        and hashlib.sha256(
            git_blob(ca["auditor_revision"], "experiments/audit_cthyb_spectroscopy.py")
        ).hexdigest()
        == ca["auditor_sha256"]
        and provenance["sources"].get("experiments/audit_cthyb_spectroscopy.py")
        == ca["auditor_sha256"],
        "calibration auditor source mismatch",
    )
    cohort = cal.get("cohort", {})
    expected_seeds = [29140001 + 2 * i for i in range(64)]
    require(
        cohort.get("stage") == "calibration"
        and cohort.get("role") == "kernelgecp-final"
        and cohort.get("generation_revision") == CAL_REVISION
        and cohort.get("seed_base") == 29100000
        and cohort.get("seeds") == expected_seeds
        and cohort.get("binding_sha256") == CAL_BINDING
        and not set(expected_seeds) & set(range(32100001, 32100513, 2)),
        "calibration cohort mismatch",
    )
    protocol = cal.get("protocol", {})
    require(
        protocol.get("methods") == ["raw", "kernelgecp", "pydlr", "sparse_ir"]
        and protocol.get("default_models") == ["flat", "broad_gaussian"]
        and protocol.get("eps") == 1e-10
        and protocol.get("rcond") == 1e-12
        and protocol.get("n_omega") == 401
        and protocol.get("n_alpha") == 33
        and protocol.get("minimizer_mu0") == 1e-8
        and protocol.get("covariance_rtol") == 1e-4
        and protocol.get("bootstrap_count") == 8
        and protocol.get("bootstrap_repeat") == 0
        and protocol.get("final_repeats") == 3
        and protocol.get("alpha_selection") == "normalized_training_only_discrepancy",
        "calibration controls mismatch",
    )
    snapshot = cal.get("snapshot", {})
    runtime = provenance["runtime"]
    require(
        snapshot.get("revision") == provenance["revision"]
        and snapshot.get("python") == runtime.get("python")
        and snapshot.get("platform") == runtime.get("platform")
        and snapshot.get("threads") == provenance["threads"]
        and isinstance(snapshot.get("sources"), dict)
        and snapshot["sources"]
        and all(
            provenance["sources"].get(p) == h for p, h in snapshot["sources"].items()
        ),
        "calibration source/runtime mismatch",
    )
    names = {"sparse_ir": "sparse-ir"}
    require(
        all(
            snapshot.get("packages", {}).get(names.get(k, k)) == v
            for k, v in runtime["versions"].items()
        ),
        "calibration provider version mismatch",
    )
    require(
        set(snapshot.get("runtime_trees", {}))
        == {
            "triqs_maxent",
            "triqs",
            "triqs_cthyb",
            "numpy",
            "scipy",
            "mpmath",
            "pydlr",
            "sparse_ir",
            "pylibsparseir",
        },
        "calibration runtime tree inventory mismatch",
    )
    require(
        set(runtime["versions"])
        == {
            "numpy",
            "scipy",
            "triqs",
            "triqs_cthyb",
            "pydlr",
            "sparse_ir",
            "pylibsparseir",
        },
        "provider version inventory mismatch",
    )
    for name, tree in snapshot["runtime_trees"].items():
        if name in ("triqs_maxent", "mpmath"):
            continue
        root = Path(tree["root"])
        selected = sorted(
            (Path(p), h)
            for p, h in runtime["provider_files"].items()
            if Path(p).is_relative_to(root)
            and Path(p).suffix in (".py", ".so", ".dylib")
            and "__pycache__" not in Path(p).parts
        )
        digest = hashlib.sha256()
        for path, value in selected:
            digest.update(str(path.relative_to(root)).encode() + b"\0")
            digest.update(bytes.fromhex(value))
        require(
            len(selected) == tree["files"] and digest.hexdigest() == tree["sha256"],
            "calibration provider files mismatch",
        )
    expected = {
        f"{m}-{d}-r0.json"
        for m in protocol["methods"]
        for d in protocol["default_models"]
    }
    cases = cal.get("cases", {})
    require(set(cases) == expected, "calibration case inventory mismatch")
    for name, record in cases.items():
        path = calibration_path.parent / name
        archive = path.with_suffix(".npz")
        require(
            sha(path) == record["result_sha256"]
            and sha(archive) == record["array_sha256"],
            "calibration case hash mismatch",
        )
        case = read_json(path)
        require(
            case.get("status") == "passed"
            and case.get("array_sha256") == record["array_sha256"]
            and case.get("eps") == 1e-10
            and case.get("rcond") == 1e-12,
            "calibration case failed",
        )
    return ref


def parent_bath(folder, binding_path):
    require(sha(binding_path) == PARENT_BINDING_SHA256, "parent binding hash mismatch")
    binding = read_json(binding_path)
    require(
        binding.get("schema") == "kernelgecp.cthyb-audit-binding.v1"
        and binding.get("status") == "passed"
        and binding.get("generation_revision") == PARENT_REVISION,
        "parent binding identity mismatch",
    )
    for name, key in [
        ("independent-audit.json", "audit_receipt_sha256"),
        ("manifest.json", "manifest_sha256"),
        ("freeze.json", "freeze_sha256"),
        ("result.json", "result_sha256"),
    ]:
        require(sha(folder / name) == binding[key], "parent file hash mismatch")
    audit, manifest, result = (
        read_json(folder / n)
        for n in ("independent-audit.json", "manifest.json", "result.json")
    )
    require(
        audit.get("schema") == "kernelgecp.cthyb-audit.v1"
        and audit.get("status") == "passed"
        and audit.get("generation_revision") == PARENT_REVISION
        and audit.get("gates")
        and all(x is True for x in audit["gates"].values()),
        "parent audit failed",
    )
    require(
        hashlib.sha256(
            git_blob(audit["audit_revision"], "experiments/audit_cthyb_dmft.py")
        ).hexdigest()
        == audit["audit_sha256"],
        "parent audit source mismatch",
    )
    require(
        manifest.get("source_revision") == PARENT_REVISION
        and result.get("identity") == manifest
        and result.get("status") == "passed"
        and result.get("freeze_sha256") == binding["freeze_sha256"],
        "parent result mismatch",
    )
    cfg = manifest["config"]
    require(
        cfg["seed_base"] == 30100000
        and cfg["final_replicas"] == 64
        and cfg["beta"] == 10
        and cfg["interaction"] == 2
        and cfg["n_iw"] == 256
        and cfg["n_tau"] == 3073,
        "parent config mismatch",
    )
    seeds = set()
    inventory = binding["measurement_archives"]
    require(inventory == result["measurement_archives"], "parent inventory mismatch")
    for name, digest in inventory.items():
        path = Path(name)
        require(
            not path.is_absolute()
            and ".." not in path.parts
            and path.suffix == ".json",
            "unsafe parent path",
        )
        p = folder / path
        require(sha(p) == digest, "parent chain hash mismatch")
        receipt = read_json(p)
        require(
            receipt.get("identity") == manifest
            and sha(p.with_suffix(".npz")) == receipt.get("arrays_sha256"),
            "parent chain identity mismatch",
        )
        seed = receipt["seed"]
        require(type(seed) is int and seed not in seeds, "parent duplicate seed")
        seeds.add(seed)
    stem = folder / "dense-final/replica-00"
    receipt = read_json(stem.with_suffix(".json"))
    require(
        receipt["seed"] == 30130001
        and receipt["cycles"] == 1600000
        and receipt["solve_status"] == 0
        and receipt["average_sign"] > 0.99,
        "parent bath chain invalid",
    )
    bath = arrays(stem.with_suffix(".npz"))["bath_iw"]
    return binding, bath, seeds


def audit_cohort(
    folder,
    parent_dir,
    parent_binding,
    receipt_path,
    binding_path,
    calibration_path,
    calibration_audit_path,
):
    folder, parent_dir, parent_binding = map(Path, (folder, parent_dir, parent_binding))
    parent, bath, parent_seeds = parent_bath(parent_dir, parent_binding)
    manifest = read_json(folder / "manifest.json")
    result = read_json(folder / "result.json")
    manifest_hash = sha(folder / "manifest.json")
    require(
        manifest.get("schema") == "kernelgecp.cthyb-fixed-bath.v1",
        "manifest schema mismatch",
    )
    protocol, ref, provenance = (
        manifest["protocol"],
        manifest["parent"],
        manifest["provenance"],
    )
    seeds = [32100001 + 2 * i for i in range(256)]
    require(
        protocol["seeds"] == seeds
        and protocol["train_count"] == 128
        and protocol["cycles"] == 1600000
        and protocol["segment_seconds"] == 14400
        and protocol["memory_limit_bytes"] == 4294967296
        and not set(seeds) & parent_seeds,
        "protocol or seed mismatch",
    )
    expected_ref = {
        "binding_sha256": PARENT_BINDING_SHA256,
        "audit_sha256": parent["audit_receipt_sha256"],
        "manifest_sha256": parent["manifest_sha256"],
        "freeze_sha256": parent["freeze_sha256"],
        "result_sha256": parent["result_sha256"],
        "bath_receipt_sha256": sha(parent_dir / "dense-final/replica-00.json"),
        "bath_archive_sha256": sha(parent_dir / "dense-final/replica-00.npz"),
        "bath_bytes_sha256": hashlib.sha256(bath.tobytes()).hexdigest(),
    }
    require(ref == expected_ref, "parent bath provenance mismatch")
    parent_config = read_json(parent_dir / "dense-final/replica-00.json")["config"]
    require(protocol["config"] == parent_config, "native config mismatch")
    require(
        provenance.get("sources")
        and provenance.get("runtime")
        and provenance.get("threads")
        and all(v == "1" for v in provenance["threads"].values()),
        "runtime identity missing",
    )
    runtime = provenance["runtime"]
    require(
        isinstance(runtime.get("python"), str)
        and isinstance(runtime.get("platform"), str)
        and "Linux" in runtime["platform"]
        and isinstance(runtime.get("versions"), dict)
        and bool(runtime["versions"])
        and isinstance(runtime.get("provider_files"), dict)
        and bool(runtime["provider_files"])
        and all(
            isinstance(d, str)
            and len(d) == 64
            and all(c in "0123456789abcdef" for c in d)
            for d in runtime["provider_files"].values()
        ),
        "runtime/provider identity invalid",
    )
    require(
        set(provenance["threads"])
        == {
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "BLIS_NUM_THREADS",
        },
        "thread environment incomplete",
    )
    tracked = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", provenance["revision"]],
        cwd=ROOT,
        text=True,
    ).splitlines()
    expected_sources = {
        p for p in tracked if p.endswith((".py", ".toml", ".lock", ".txt"))
    }
    require(
        set(provenance["sources"]) == expected_sources and expected_sources,
        "generation source inventory mismatch",
    )
    validate_calibration(
        calibration_path, calibration_audit_path, manifest["calibration"], provenance
    )
    for path, digest in provenance["sources"].items():
        require(
            hashlib.sha256(git_blob(provenance["revision"], path)).hexdigest()
            == digest,
            "generation source mismatch",
        )
    require(
        result.get("schema") == "kernelgecp.cthyb-fixed-bath-result.v1"
        and result.get("status") == "complete"
        and result.get("manifest_sha256") == manifest_hash
        and len(result.get("chains", [])) == 256,
        "incomplete acquisition",
    )
    expected_files = {
        f"chain-{i:03d}{suffix}"
        for i in range(256)
        for suffix in (".json", ".npz", ".timing.json")
    }
    actual_files = {p.name for p in (folder / "chains").iterdir()}
    require(actual_files == expected_files, "orphan or missing chain file")
    inventory = {}
    for i, entry in enumerate(result["chains"]):
        role = "train" if i < 128 else "heldout"
        stem = folder / "chains" / f"chain-{i:03d}"
        paths = [
            stem.with_suffix(".json"),
            stem.with_suffix(".npz"),
            stem.with_suffix(".timing.json"),
        ]
        digests = [sha(p) for p in paths]
        require(
            entry
            == {
                "index": i,
                "seed": seeds[i],
                "role": role,
                "receipt_sha256": digests[0],
                "arrays_sha256": digests[1],
                "timing_sha256": digests[2],
            },
            "chain inventory mismatch",
        )
        inventory.update(
            {str(p.relative_to(folder)): h for p, h in zip(paths, digests, strict=True)}
        )
        raw, timing = read_json(paths[0]), read_json(paths[2])
        require(
            raw.get("schema") == "kernelgecp.cthyb-measurement.v1"
            and raw.get("identity")
            == {"manifest_sha256": manifest_hash, "index": i, "role": role}
            and raw.get("arrays_sha256") == digests[1]
            and raw.get("seed") == seeds[i]
            and raw.get("cycles") == 1600000
            and raw.get("config") == parent_config,
            "chain identity/config mismatch",
        )
        require(
            raw.get("spin_blocks") == ["up", "down"]
            and raw.get("ordinary_indices") == "-n_iw..n_iw-1"
            and raw.get("density_source") == "native_reduced_density_matrix"
            and raw.get("moment_source") == "native_density_matrix_commutators"
            and raw.get("thermalization") == "fresh_solver_and_warmup_for_each_archive",
            "native metadata mismatch",
        )
        sign, double = raw["average_sign"], raw["double_occupancy"]
        require(
            type(raw["solve_status"]) is int
            and raw["solve_status"] == 0
            and type(sign) in (int, float)
            and math.isfinite(sign)
            and 0.99 < sign <= 1.000000001
            and type(double) in (int, float)
            and math.isfinite(double)
            and 0 <= double <= 1,
            "native solve/sign/density invalid",
        )
        require(
            timing.get("schema") == "kernelgecp.cthyb-fixed-bath-timing.v1"
            and timing.get("seed") == seeds[i]
            and type(timing.get("solve_status")) is int
            and timing.get("solve_status") == 0
            and timing.get("average_sign") == sign
            and timing.get("arrays_sha256") == digests[1]
            and timing.get("receipt_sha256") == digests[0]
            and 0 < timing.get("elapsed_seconds", 0) <= 14400
            and 0 < timing.get("peak_rss_bytes", 0) <= 4294967296,
            "timing mismatch",
        )
        a = arrays(paths[1])
        require(np.array_equal(a["bath_iw"], bath), "bath bits mismatch")
        require(
            np.max(np.abs(a["g_moments"][:2] - np.array([[0, 0], [1, 1]]))) <= 1e-12,
            "native Green zeroth/unit moment mismatch",
        )
        occ = a["occupancies"]
        require(
            np.all((occ >= 0) & (occ <= 1))
            and double <= min(occ) + 1e-10
            and double >= max(0, sum(occ) - 1) - 1e-10,
            "density consistency mismatch",
        )
        require(
            np.max(np.abs(a["sigma_moments"][0] - 2 * occ[::-1])) < 1e-10
            and np.max(np.abs(a["g_moments"][2] - (a["sigma_moments"][0] - 1))) < 1e-10,
            "native moment mismatch",
        )
        require(
            np.max(np.abs(a["g_iw"] - a["g_iw"][::-1].conj())) < 1e-10
            and np.max(a["g_iw"][256:].imag) < 1e-8,
            "Green function consistency mismatch",
        )
        mask = np.any(a["g_iw"] != a["g_iw_raw"], axis=1)
        require(np.any(mask) and not np.any(mask[252:260]), "native tail mask mismatch")
        iw = 1j * (2 * np.arange(-256, 256) + 1) * np.pi / 10
        tail = (
            a["g_moments"][1][None, :] / iw[:, None]
            + a["g_moments"][2][None, :] / iw[:, None] ** 2
        )
        require(
            np.max(np.abs(a["g_iw"][mask] - tail[mask])) < 1e-12,
            "native tail moment mismatch",
        )
    source = Path(__file__).resolve()
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    source_hash = sha(source)
    require(
        hashlib.sha256(
            git_blob(revision, "experiments/audit_cthyb_fixed_bath.py")
        ).hexdigest()
        == source_hash,
        "audit source must match committed revision",
    )
    receipt = {
        "schema": "kernelgecp.cthyb-fixed-bath-audit.v1",
        "status": "passed",
        "claim_scope": "fixed_bath_raw_cohort_only_not_EE_campaign",
        "chain_count": 256,
        "train_count": 128,
        "heldout_count": 128,
        "manifest_sha256": manifest_hash,
        "result_sha256": sha(folder / "result.json"),
        "parent_binding_sha256": PARENT_BINDING_SHA256,
        "calibration": manifest["calibration"],
        "audit_revision": revision,
        "audit_source_sha256": source_hash,
    }
    receipt_path, binding_path = Path(receipt_path), Path(binding_path)
    require(
        not receipt_path.exists() and not binding_path.exists(),
        "refusing to overwrite audit outputs",
    )
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    binding = {
        "schema": "kernelgecp.cthyb-fixed-bath-audit-binding.v1",
        "status": "passed",
        "receipt_sha256": sha(receipt_path),
        "manifest_sha256": manifest_hash,
        "result_sha256": receipt["result_sha256"],
        "inventory": inventory,
        "parent_binding_sha256": PARENT_BINDING_SHA256,
        "calibration": manifest["calibration"],
        "audit_revision": revision,
        "audit_source_sha256": source_hash,
    }
    binding_path.write_text(json.dumps(binding, indent=2, sort_keys=True) + "\n")
    return receipt, binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "folder",
        "parent-dir",
        "parent-binding",
        "receipt",
        "binding",
        "calibration",
        "calibration-audit",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit_cohort(
        args.folder,
        args.parent_dir,
        args.parent_binding,
        args.receipt,
        args.binding,
        args.calibration,
        args.calibration_audit,
    )


if __name__ == "__main__":
    main()
