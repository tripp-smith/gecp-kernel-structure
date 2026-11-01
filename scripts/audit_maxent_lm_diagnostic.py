#!/usr/bin/env python3
"""Audit six diagnostic reports; no production acceptance or solver rerun."""

import argparse
import hashlib
import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_exact(a):
    rows = [[Fraction(float(x)) for x in row] for row in a]
    m, n = len(rows), len(rows[0])
    rank = 0
    for c in range(n):
        p = next((i for i in range(rank, m) if rows[i][c]), None)
        if p is None:
            continue
        rows[rank], rows[p] = rows[p], rows[rank]
        pivot = rows[rank][c]
        rows[rank] = [x / pivot for x in rows[rank]]
        for i in range(rank + 1, m):
            k = rows[i][c]
            if k:
                rows[i] = [x - k * y for x, y in zip(rows[i], rows[rank], strict=True)]
        rank += 1
        if rank == m:
            break
    return rank


def fail_constant(value):
    raise ValueError(f"Nonstandard JSON numeric constant: {value}")


def check_hex(value, n):
    return (
        isinstance(value, str) and re.fullmatch(rf"[0-9a-f]{{{n}}}", value) is not None
    )


def audit(root):
    paths = sorted(root.iterdir())
    if len(paths) != 6 or any(not p.is_dir() for p in paths):
        raise ValueError("reports root must contain exactly six report directories")
    pattern = re.compile(r"maxent-lm-(unit|calibration)-([123])-(.+)")
    labels = set()
    suffixes = set()
    common = None
    configs = {}
    runs = []
    for folder in paths:
        match = pattern.fullmatch(folder.name)
        if match is None:
            raise ValueError(f"unexpected report directory {folder.name}")
        kind, replica, suffix = match.groups()
        labels.add((kind, replica))
        suffixes.add(suffix)
        report = folder / "result.json"
        data = json.loads(report.read_text(), parse_constant=fail_constant)
        if (
            data.get("schema") != "kernelgecp.maxent-lm-diagnostic.v1"
            or data.get("scope")
            != "Consumed calibration fixtures only; not production acceptance."
        ):
            raise ValueError(f"unexpected schema or scope in {report}")
        identity = {
            k: data[k]
            for k in (
                "repo_head",
                "script_sha256",
                "integration_sha256",
                "provider",
                "python",
                "numpy",
                "scipy",
            )
        }
        identity["minimizer_sha256"] = data["provider_source_hashes"][
            "minimizers/levenberg_minimizer.py"
        ]
        if (
            not check_hex(identity["repo_head"], 40)
            or not all(
                check_hex(identity[k], 64)
                for k in ("script_sha256", "integration_sha256", "minimizer_sha256")
            )
            or not isinstance(identity["provider"], str)
            or not identity["provider"]
            or not all(
                isinstance(identity[k], str) and identity[k]
                for k in ("python", "numpy", "scipy")
            )
            or not str(data.get("platform", "")).startswith("Linux-")
            or not check_hex(data.get("fixture_fingerprint"), 64)
            or (kind == "calibration" and not check_hex(data.get("dataset_sha256"), 64))
        ):
            raise ValueError("invalid source identity")
        if common is None:
            common = identity
        elif common != identity:
            raise ValueError("inconsistent source/provider identity")
        config = data["config"]
        if not isinstance(config, dict) or config.get("n_alpha") != (
            9 if kind == "unit" else 33
        ):
            raise ValueError("unexpected fixture config")
        if kind in configs and configs[kind] != config:
            raise ValueError("inconsistent per-fixture configuration")
        configs[kind] = config
        trials = data["trials"]
        if [t.get("name") for t in trials] != ["baseline", "mu0_1e-12", "mu0_1e-8"]:
            raise ValueError("unexpected or incomplete trial list")
        result = {
            "directory": folder.name,
            "result_sha256": digest(report),
            "fixture_fingerprint": data.get("fixture_fingerprint"),
            "dataset_sha256": data.get("dataset_sha256"),
            "trials": [],
        }
        captures = []
        for t, mu in zip(trials, (1e-18, 1e-12, 1e-8), strict=True):
            if (
                t.get("effective_mu0") != mu
                or t.get("status") not in ("completed", "exception")
                or not isinstance(t.get("history"), list)
            ):
                raise ValueError("invalid trial metadata")
            row = {
                "name": t["name"],
                "status": t["status"],
                "effective_mu0": mu,
                "completed_history": sum(
                    h.get("converged") is True for h in t["history"]
                ),
            }
            if t["status"] == "completed":
                if (
                    len(t.get("alphas", [])) != config["n_alpha"]
                    or len(t.get("converged", [])) != config["n_alpha"]
                ):
                    raise ValueError("incomplete completed trial")
                row.update(
                    {
                        "alpha_count": len(t["alphas"]),
                        "all_converged": all(t["converged"]),
                        "selected_index": t.get("selected_index"),
                        "normalized_training_whitened_rms": t.get(
                            "normalized_training_whitened_rms"
                        ),
                    }
                )
            else:
                if (
                    not isinstance(t.get("exception_chain"), str)
                    or not t["exception_chain"]
                ):
                    raise ValueError("exception record missing")
                row["original_linux_exception_chain"] = t["exception_chain"]
                row["failure_metadata"] = t.get("failure", {})
                row["capture_recorded"] = bool(t.get("failure"))
                if t.get("failure"):
                    captures.append(t["failure"])
            result["trials"].append(row)
        npz = folder / "result.failure.npz"
        if len(captures) > 1 or (bool(captures) != npz.is_file()):
            raise ValueError("missing or ambiguous failure capture")
        if captures:
            meta = captures[0]
            if (
                not check_hex(meta.get("npz_sha256"), 64)
                or digest(npz) != meta["npz_sha256"]
            ):
                raise ValueError("failure NPZ hash mismatch")
            with np.load(npz, allow_pickle=False) as archive:
                if set(archive.files) != {"matrix", "rhs", "mu", "J", "f", "Id"}:
                    raise ValueError("invalid capture keys")
                z = {k: archive[k] for k in archive.files}
            if any(
                x.dtype != np.dtype("float64") or not np.isfinite(x).all()
                for x in z.values()
            ):
                raise ValueError("capture must be finite float64")
            a, j, b, ident = z["matrix"], z["J"], z["rhs"], z["Id"]
            n = len(b)
            if (
                n < 1
                or n > 100
                or b.shape != (n,)
                or z["f"].shape != (n,)
                or any(x.shape != (n, n) for x in (a, j, ident))
                or z["mu"].shape != ()
                or list(a.shape) != meta.get("matrix_shape")
            ):
                raise ValueError("invalid capture shapes")
            mu = float(z["mu"])
            if (
                mu != meta.get("mu")
                or not np.array_equal(a, j + mu * ident)
                or not np.array_equal(b, z["f"])
                or not np.array_equal(ident, np.eye(n))
            ):
                raise ValueError("inconsistent captured linear system")
            sv = np.linalg.svd(a, compute_uv=False)
            eig = np.linalg.eigvals(a)
            x, _, numerical_rank, _ = np.linalg.lstsq(a, b, rcond=None)
            bnorm = float(np.linalg.norm(b))
            entry = {
                "npz_sha256": digest(npz),
                "shape": list(a.shape),
                "mu": mu,
                "matrix_equals_rounded_j_plus_mu_id": True,
                "rhs_equals_f": True,
                "id_equals_identity": True,
                "unchanged_diagonal_count": int(np.sum(np.diag(a) == np.diag(j))),
                "diagonal_ulp": np.abs(np.spacing(np.diag(j))).tolist(),
                "singular_values": sv.tolist(),
                "eigenvalues_real_imag": [[float(v.real), float(v.imag)] for v in eig],
                "exact_binary_matrix_rank": rank_exact(a),
                "exact_binary_augmented_rank": rank_exact(np.column_stack((a, b))),
                "lstsq_rank_default_tolerance": int(numerical_rank),
                "lstsq_relative_residual": float(np.linalg.norm(a @ x - b) / bnorm)
                if bnorm
                else None,
                "fixed_j_hypotheticals": {},
            }
            for hypothetical in (1e-18, 1e-12, 1e-8):
                shifted = j + hypothetical * ident
                s = np.linalg.svd(shifted, compute_uv=False)
                cond = float(s[0] / s[-1]) if s[-1] else None
                entry["fixed_j_hypotheticals"][str(hypothetical)] = {
                    "unchanged_diagonal_count": int(
                        np.sum(np.diag(shifted) == np.diag(j))
                    ),
                    "singular_values": s.tolist(),
                    "condition_estimate": cond
                    if cond is None or math.isfinite(cond)
                    else None,
                    "exact_binary_rank": rank_exact(shifted),
                }
            result["captured_system"] = entry
        runs.append(result)
    if (
        labels
        != {(kind, str(i)) for kind in ("unit", "calibration") for i in (1, 2, 3)}
        or len(suffixes) != 1
    ):
        raise ValueError(
            "expected one common campaign and unit/calibration replicas 1, 2, 3"
        )
    return {
        "schema": "gecp.lm-independent-audit.v1",
        "scope": (
            "Diagnostic numerical evidence only; original Linux exceptions preserved; "
            "not production acceptance or a local solver rerun."
        ),
        "python": sys.version,
        "numpy": np.__version__,
        "source_identity": common,
        "configs": configs,
        "runs": runs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports_root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(args.reports_root)
        output = json.dumps(result, indent=2, allow_nan=False) + "\n"
        args.output.write_text(output)
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        IndexError,
        np.linalg.LinAlgError,
    ) as error:
        parser.exit(1, f"audit failed: {error}\n")


if __name__ == "__main__":
    main()
