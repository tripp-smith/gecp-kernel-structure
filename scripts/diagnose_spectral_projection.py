#!/usr/bin/env python3
"""Diagnose saved projections from the pinned archive without provider execution."""

import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import sys
import tarfile
from pathlib import Path, PurePosixPath

import numpy as np
import scipy
from scipy.linalg import svd

ARCHIVE_SHA256 = "63ac446677a9c36dfb17736e9823582a92ce49c59d85c3c0f393fd35c4ee2576"
ARCHIVE_BYTES = 55_453_529
WIDTHS = {"raw": 65, "kernelgecp": 24, "pydlr": 23, "sparse_ir": 23}
EXPECTED = {
    f"native-calibration/{m}-{d}-r0.npz"
    for m in WIDTHS
    for d in ("flat", "broad_gaussian")
}


def load_archive(path):
    if path.stat().st_size != ARCHIVE_BYTES:
        raise ValueError("archive byte count mismatch")
    with path.open("rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError("archive SHA-256 mismatch")
    with tarfile.open(path, "r:gz") as tar:
        entries = tar.getmembers()
        if len(entries) > 4000 or sum(e.size for e in entries) > 150_000_000:
            raise ValueError("archive resource cap exceeded")
        seen, selected = set(), {}
        for e in entries:
            name = e.name
            if (
                not name
                or name.startswith("/")
                or "\\" in name
                or any(x in ("", ".", "..") for x in name.rstrip("/").split("/"))
                or not (e.isfile() or e.isdir())
                or e.size > 16_000_000
            ):
                raise ValueError("unsafe archive member")
            name = str(PurePosixPath(name))
            if name in seen:
                raise ValueError("duplicate archive member")
            seen.add(name)
            if name in EXPECTED:
                if not e.isfile():
                    raise ValueError("expected member is not regular")
                selected[name] = e
        if set(selected) != EXPECTED:
            raise ValueError("calibration array inventory mismatch")
        result = {}
        for name, e in sorted(selected.items()):
            with tar.extractfile(e) as f:
                data = f.read()
            result[name] = read_matrices(
                data, WIDTHS[PurePosixPath(name).name.split("-")[0]]
            )
        return result


def read_matrices(data, width):
    with np.load(io.BytesIO(data), allow_pickle=False) as z:
        design, projection = z["observation_design"], z["projection"]
    if (
        design.shape != (65, width)
        or projection.shape != (65, 65)
        or not np.isrealobj(design)
        or not np.isrealobj(projection)
        or not np.isfinite(design).all()
        or not np.isfinite(projection).all()
    ):
        raise ValueError("invalid design or projection matrices")
    return design, projection


def compare(p, q):
    diff = np.abs(p - q)
    tol = 2e-8 + 1e-8 * np.abs(q)
    return {
        "max_abs": float(diff.max()),
        "max_tolerance_ratio": float((diff / tol).max()),
        "within_tolerance": bool(np.all(diff <= tol)),
    }


def diagnostic(a, p):
    eig = np.linalg.eigvalsh((p + p.T) / 2)
    result = {
        "shape": list(a.shape),
        "projection_trace": float(np.trace(p)),
        "projection_symmetry_max_abs": float(np.max(np.abs(p - p.T))),
        "projection_idempotence_max_abs": float(np.max(np.abs(p @ p - p))),
        "projection_eigenvalues_below_half": int((eig < 0.5).sum()),
        "projection_eigenvalues_at_least_half": int((eig >= 0.5).sum()),
        "methods": {},
    }
    for method in ("numpy", "gesdd", "gesvd"):
        u, s, _ = (
            np.linalg.svd(a, full_matrices=False)
            if method == "numpy"
            else svd(a, full_matrices=False, lapack_driver=method)
        )
        keep = s > 1e-12 * s[0]
        rank = int(keep.sum())
        if rank == 0:
            raise ValueError("no retained design modes")
        q = u[:, keep] @ u[:, keep].T
        result["methods"][method] = {
            **compare(p, q),
            "rank": rank,
            "smallest_retained_singular_ratio": float(s[rank - 1] / s[0]),
            "largest_discarded_singular_ratio": float(s[rank] / s[0])
            if rank < len(s)
            else None,
            "retained_condition": float(s[0] / s[rank - 1]),
        }
    return result


def high_precision(a, p):
    import mpmath as mp

    with mp.workdps(50):
        matrix = mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
        u, s, _ = mp.svd(matrix, full_matrices=False)
        rank = sum(x > mp.mpf("1e-12") * s[0] for x in s)
        if rank == 0:
            raise ValueError("no retained high-precision modes")
        u64 = np.array(
            [[float(u[i, j]) for j in range(rank)] for i in range(a.shape[0])]
        )
        return {
            **compare(p, u64 @ u64.T),
            "rank": rank,
            "decimal_digits": 50,
            "smallest_retained_singular_ratio": float(s[rank - 1] / s[0]),
        }


def runtime():
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        np.show_config()
        np.show_runtime()
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numpy_config_and_runtime": stream.getvalue(),
        "thread_environment": {
            x: os.environ.get(x)
            for x in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--high-precision", action="store_true")
    args = parser.parse_args(argv)
    if args.output.resolve() == args.archive.resolve():
        print("output must not overwrite the archive", file=sys.stderr)
        return 1
    try:
        arrays = load_archive(args.archive)
        cases = {name: diagnostic(*v) for name, v in arrays.items()}
        worst = max(
            (name for name in cases if not Path(name).name.startswith("raw-")),
            key=lambda name: cases[name]["methods"]["numpy"]["max_tolerance_ratio"],
        )
        result = {
            "status": "diagnostic_complete",
            "archive_sha256": ARCHIVE_SHA256,
            "runtime": runtime(),
            "relative_svd_cutoff": 1e-12,
            "comparison_atol": 2e-8,
            "comparison_rtol": 1e-8,
            "cases": cases,
            "worst_nonraw_numpy_case": worst,
            "high_precision": high_precision(*arrays[worst])
            if args.high_precision
            else None,
            "scope": "Diagnostic only; no scientific acceptance or provider execution.",
        }
        code = 0
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as exc:
        result, code = {"status": "diagnostic_failed", "reason": str(exc)}, 1
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
