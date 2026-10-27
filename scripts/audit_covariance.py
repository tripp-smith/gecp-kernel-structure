"""Independently recompute the complete frozen EI receipt and real-data gates."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from experiments import covariance_benchmark as benchmark


class CovarianceAuditError(ValueError):
    """Missing, stale, malformed, or numerically unsuccessful EI evidence."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CovarianceAuditError(message)


def _pairs(entries):
    result = {}
    for key, value in entries:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: Path) -> dict[str, Any]:
    def invalid(value):
        raise CovarianceAuditError(f"nonfinite JSON constant: {value}")

    result = json.loads(
        path.read_text(), object_pairs_hook=_pairs, parse_constant=invalid
    )
    require(type(result) is dict, "JSON object required")
    return result


def artifact(parent: Path, name: object, expected_hash: object) -> Path:
    require(
        isinstance(name, str) and Path(name).name == name,
        "artifact must be a local filename",
    )
    path = (parent / str(name)).resolve()
    require(
        path.parent == parent.resolve() and path.is_file(), "missing/escaped artifact"
    )
    require(
        isinstance(expected_hash, str)
        and re.fullmatch(r"[0-9a-f]{64}", expected_hash) is not None,
        "invalid artifact hash",
    )
    require(
        benchmark.sha256(path.read_bytes()) == expected_hash, "artifact hash mismatch"
    )
    return path


def finite_number(value: object, name: str) -> float:
    require(type(value) in (float, int), f"{name} must be a numeric scalar")
    result = float(value)
    require(math.isfinite(result), f"{name} must be finite")
    return result


def independent_rbf(left, right, lengthscale):
    """Independent 2D sum of coordinate differences, without cdist/3D tensors."""
    distances = np.zeros((len(left), len(right)), dtype=np.float64)
    for dimension in range(left.shape[1]):
        difference = left[:, dimension, None] - right[None, :, dimension]
        distances += (difference / lengthscale) ** 2
    return np.exp(-0.5 * distances)


def original_action(x, solution, lengthscale, noise):
    result = noise * solution
    for start in range(0, len(x), 64):
        end = min(start + 64, len(x))
        result[start:end] += independent_rbf(x[start:end], x, lengthscale) @ solution
    return result


def _arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        result = {name: archive[name] for name in archive.files}
    for name, value in result.items():
        require(
            value.dtype.kind in "fiu" and np.all(np.isfinite(value)),
            f"invalid {name} array",
        )
    return result


def _close(actual, expected, name, *, atol=1e-10):
    require(np.shape(actual) == np.shape(expected), f"wrong {name} shape")
    require(
        bool(np.allclose(actual, expected, rtol=1e-10, atol=atol)),
        f"inconsistent {name}",
    )


def _source(receipt, frozen, revision):
    require(
        receipt.get("standard_deviation_definition") == benchmark.STD_DEFINITION
        and frozen.get("standard_deviation_definition") == benchmark.STD_DEFINITION,
        "ambiguous predictive standard-deviation convention",
    )
    current = receipt["provenance"]
    require(current["source_revision"] == revision, "wrong/stale source revision")
    require(
        current["source_clean_and_tracked"] is True, "uncommitted scientific source"
    )
    require(
        current["source_hashes"] == frozen["provenance"]["source_hashes"],
        "wrong source hashes",
    )
    require(
        current["versions"] == frozen["provenance"]["versions"], "dependency mismatch"
    )
    require(
        all(current["threads"].get(key) == "1" for key in benchmark.THREADS),
        "wrong thread policy",
    )
    require(receipt["dataset"] == frozen["dataset"], "wrong data provenance")
    digest = benchmark.sha256(json.dumps(frozen, sort_keys=True).encode())
    require(receipt["frozen_calibration_sha256"] == digest, "wrong frozen calibration")


def bind_source_to_git(record: dict[str, Any], revision: str | None = None) -> None:
    """Bind every declared scientific source hash to its actual immutable blob."""
    declared = record.get("source_revision")
    require(
        isinstance(declared, str)
        and re.fullmatch(r"[0-9a-f]{40}", declared) is not None,
        "full scientific source revision required",
    )
    if revision is not None:
        require(declared == revision, "declared generation revision mismatch")
    hashes = record.get("source_hashes")
    require(
        type(hashes) is dict and set(hashes) == set(benchmark.SOURCE_FILES),
        "unexpected scientific source keyset",
    )
    try:
        resolved = subprocess.check_output(
            ["git", "rev-parse", "--verify", f"{declared}^{{commit}}"],
            cwd=benchmark.ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        require(
            resolved == declared, "scientific source must identify the commit itself"
        )
        for name in benchmark.SOURCE_FILES:
            digest = hashes[name]
            require(
                isinstance(digest, str)
                and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
                "invalid scientific source hash",
            )
            content = subprocess.check_output(
                ["git", "show", f"{declared}:{name}"],
                cwd=benchmark.ROOT,
                stderr=subprocess.DEVNULL,
            )
            require(
                benchmark.sha256(content) == digest,
                f"scientific source differs from Git blob: {name}",
            )
    except subprocess.CalledProcessError as error:
        raise CovarianceAuditError("scientific revision/blob is unavailable") from error


def audit_case(receipt, arrays, frozen, data, revision, *, previous=None):
    require(
        receipt.get("schema") == "kernelgecp.covariance-case.v1", "wrong case schema"
    )
    require(
        receipt.get("canonical") is True and receipt.get("status") == "passed",
        "case did not pass canonically",
    )
    _source(receipt, frozen, revision)
    size, repeat = receipt["size"], receipt["repeat"]
    require(type(size) is int and size in benchmark.SIZES, "wrong size")
    require(type(repeat) is int and repeat in (0, 1, 2), "wrong repeat")
    rank = receipt["actual_rank"]
    require(type(rank) is int and 0 <= rank <= receipt["rank_budget"], "invalid rank")
    require(
        receipt["rank_budget"] == frozen["selected_ranks"][str(size)],
        "rank changed after freeze",
    )
    checks = receipt["checks"]
    require(
        set(checks)
        == {
            "dense_accuracy",
            "approximate_solve",
            "original_kernel_cg",
            "retained_state",
            "rss",
            "source_clean",
        },
        "missing case gates",
    )
    require(all(value is True for value in checks.values()), "a case gate did not pass")
    expected_keys = {
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
    }
    require(set(arrays) == expected_keys, "missing/extra numerical arrays")
    for key in ("mean", "std", "query_y"):
        require(arrays[key].shape == (512,), f"wrong {key} shape")
    for key in ("cg_solution", "cg_residual", "cg_rhs", "approximate_solution"):
        require(arrays[key].shape == (size,), f"wrong {key} shape")
    root, pivots = arrays["root"], arrays["pivots"]
    require(root.shape == (size, rank), "wrong retained root shape")
    require(pivots.dtype.kind in "iu" and pivots.shape == (rank,), "wrong pivot array")
    require(
        len(set(pivots.tolist())) == rank and np.all((pivots >= 0) & (pivots < size)),
        "repeated/out-of-range pivot",
    )
    split = benchmark.split_data(data, size, "test")
    require(
        receipt["normalization"] == split["normalization"],
        "training normalization mismatch",
    )
    for name in receipt["split"]:
        require(receipt["split"][name] == split[name], "chronological split mismatch")
    _close(arrays["cg_rhs"], split["y"], "real training target", atol=0)
    _close(arrays["query_y"], split["query_y"], "real test target", atol=0)
    lengthscale, noise = frozen["lengthscale"], frozen["noise"]
    residual = (
        root @ (root.T @ arrays["approximate_solution"])
        + noise * arrays["approximate_solution"]
        - split["y"]
    )
    solve_error = float(np.linalg.norm(residual) / np.linalg.norm(split["y"]))
    require(solve_error <= 1e-10, "approximate solve residual exceeds acceptance")
    require(
        0
        <= finite_number(receipt["solve_relative_residual"], "solve residual")
        <= 1e-10,
        "bad reported solve residual",
    )
    # Later repetitions must reproduce actual arrays bit for bit. Reuse only
    # previously independently checked exact numerical values, never metrics.
    if previous is not None:
        for name in arrays:
            require(
                np.array_equal(arrays[name], previous[name]),
                "fresh repeats disagree numerically",
            )
    else:
        cg_residual = (
            original_action(split["x"], arrays["cg_solution"], lengthscale, noise)
            - split["y"]
        )
        cg_error = float(np.linalg.norm(cg_residual) / np.linalg.norm(split["y"]))
        require(
            cg_error <= 1e-8,
            "original-kernel independent CG residual exceeds acceptance",
        )
        _close(
            arrays["cg_residual"], cg_residual, "original-kernel residual", atol=5e-11
        )
        selected_columns = independent_rbf(split["x"], split["x"][pivots], lengthscale)
        _close(root @ root[pivots].T, selected_columns, "selected true-kernel columns")
        selected_cross = independent_rbf(
            split["query_x"], split["x"][pivots], lengthscale
        )
        features = solve_triangular(root[pivots], selected_cross.T, lower=True).T
        mean = features @ (root.T @ arrays["approximate_solution"])
        chol = cho_factor(noise * np.eye(rank) + root.T @ root, lower=True)
        residual_prior = 1 - np.sum(features**2, axis=1)
        std = np.sqrt(
            noise
            + residual_prior
            + noise * np.sum(features * cho_solve(chol, features.T).T, axis=1)
        )
        _close(arrays["mean"], mean, "approximate prediction mean")
        _close(arrays["std"], std, "approximate prediction std")
        if size <= 2048:
            matrix = independent_rbf(
                split["x"], split["x"], lengthscale
            ) + noise * np.eye(size)
            full_chol = cho_factor(matrix, lower=True)
            cross = independent_rbf(split["query_x"], split["x"], lengthscale)
            expected_mean = cross @ cho_solve(full_chol, split["y"])
            expected_std = np.sqrt(
                1 + noise - np.sum(cross * cho_solve(full_chol, cross.T).T, axis=1)
            )
            _close(arrays["expected_mean"], expected_mean, "independent dense mean")
            _close(arrays["expected_std"], expected_std, "independent dense std")
    cg_error = float(np.linalg.norm(arrays["cg_residual"]) / np.linalg.norm(split["y"]))
    require(
        cg_error <= 1e-8
        and type(receipt["cg_info"]) is int
        and receipt["cg_info"] == 0,
        "CG did not converge",
    )
    require(
        math.isclose(
            cg_error,
            finite_number(receipt["cg_relative_residual"], "CG residual"),
            rel_tol=1e-8,
            abs_tol=1e-14,
        ),
        "reported CG residual mismatch",
    )
    require(np.all(arrays["std"] > 0), "nonpositive prediction uncertainty")
    if size <= 2048:
        metrics = benchmark.accuracy(
            arrays["mean"],
            arrays["std"],
            arrays["expected_mean"],
            arrays["expected_std"],
        )
        require(
            max(metrics["mean_rms"], metrics["std_rms"]) <= 1e-3
            and max(metrics["mean_max"], metrics["std_max"]) <= 1e-2,
            "dense prediction accuracy gate failed",
        )
        for name, value in metrics.items():
            require(
                math.isclose(
                    value,
                    finite_number(receipt["accuracy"][name], name),
                    rel_tol=1e-8,
                    abs_tol=1e-14,
                ),
                "reported accuracy mismatch",
            )
    else:
        require(
            receipt["accuracy"] is None
            and arrays["expected_mean"].size == arrays["expected_std"].size == 0,
            "unsupported large dense reference claim",
        )
    retained = finite_number(receipt["retained_numerical_bytes"], "retained state")
    require(
        retained == sum(receipt["retained_components_bytes"].values()),
        "incomplete retained-state sum",
    )
    require(
        root.nbytes <= retained <= 128 * 1024**2,
        "retained numerical memory gate failed",
    )
    require(
        0 < finite_number(receipt["peak_rss_bytes"], "RSS") <= 4 * 1024**3,
        "RSS gate failed",
    )
    for name, value in receipt["timing_seconds"].items():
        if name == "dense_reference" and size > 2048:
            require(value is None, "unexpected large dense timing")
        else:
            require(finite_number(value, name) > 0, "invalid timing")
    return {
        "size": size,
        "repeat": repeat,
        "solve_residual": solve_error,
        "cg_residual": cg_error,
    }


def audit_campaign(
    manifest_path: Path,
    frozen_path: Path,
    data_path: Path,
    provider_path: Path,
    revision: str,
):
    require(
        re.fullmatch(r"[0-9a-f]{40}", revision) is not None, "full revision required"
    )
    manifest, frozen = read_json(manifest_path), read_json(frozen_path)
    benchmark.validate_frozen(frozen, canonical=False)
    bind_source_to_git(frozen["provenance"])
    bind_source_to_git(manifest["provenance"], revision)
    require(
        manifest["provenance"]["source_hashes"]
        == frozen["provenance"]["source_hashes"],
        "generation Git blobs differ from frozen calibration Git blobs",
    )
    require(
        frozen["sizes"] == list(benchmark.SIZES)
        and frozen["ranks"] == list(benchmark.RANKS),
        "incomplete calibration grid",
    )
    require(len(frozen["rank_calibration"]) == 30, "missing rank calibration cases")
    calibration_grid = {
        (row["size"], row["rank_budget"]) for row in frozen["rank_calibration"]
    }
    require(
        calibration_grid
        == {(size, rank) for size in benchmark.SIZES for rank in benchmark.RANKS},
        "duplicate/incomplete rank calibration grid",
    )
    require(
        manifest["schema"] == "kernelgecp.covariance-campaign.v1",
        "wrong campaign schema",
    )
    require(
        manifest["provenance"]["source_revision"] == revision, "stale campaign source"
    )
    expected = {(size, repeat) for size in benchmark.SIZES for repeat in range(3)}
    require(len(manifest["cases"]) == len(expected), "missing fresh-process repeats")
    data = benchmark.load_data(data_path)
    require(
        frozen["dataset"] == benchmark.dataset_identity(data),
        "frozen dataset identity mismatch",
    )
    checked, previous, summaries = set(), {}, []
    for item in manifest["cases"]:
        require(
            type(item["exit_code"]) is int and item["exit_code"] == 0,
            "failed child process",
        )
        path = artifact(manifest_path.parent, item["path"], item["sha256"])
        receipt = read_json(path)
        identity = (receipt["size"], receipt["repeat"])
        require(
            identity in expected and identity not in checked,
            "duplicate/unexpected case",
        )
        checked.add(identity)
        arrays = _arrays(
            artifact(
                path.parent, receipt["artifact"]["path"], receipt["artifact"]["sha256"]
            )
        )
        summaries.append(
            audit_case(
                receipt,
                arrays,
                frozen,
                data,
                revision,
                previous=previous.get(identity[0]),
            )
        )
        if identity[0] not in previous:
            previous[identity[0]] = arrays
    require(checked == expected, "incomplete canonical cases")
    provider = read_json(provider_path)
    require(
        provider.get("schema") == "kernelgecp.covariance-provider.v1"
        and provider.get("status") == "passed",
        "actual real-data provider missing/failed",
    )
    _source(provider, frozen, revision)
    require(
        provider["provider"] == "gpytorch" and provider["size"] == 512,
        "wrong actual provider",
    )
    require(
        set(provider["checks"]) == {"cpu_float64", "frozen_factor", "numerical"}
        and all(value is True for value in provider["checks"].values()),
        "actual provider gate failed",
    )
    require(
        {name: value.split("+")[0] for name, value in provider["versions"].items()}
        == {
            "numpy": "2.2.6",
            "scipy": "1.18.1",
            "torch": "2.14.1",
            "gpytorch": "1.15.2",
            "linear_operator": "0.6.1",
        },
        "wrong actual provider versions",
    )
    values = _arrays(
        artifact(
            provider_path.parent,
            provider["artifact"]["path"],
            provider["artifact"]["sha256"],
        )
    )
    _close(
        values["root"],
        previous[512]["root"],
        "identical frozen provider factor",
        atol=0,
    )
    _close(values["mean"], previous[512]["mean"], "actual provider mean", atol=1e-9)
    _close(values["std"], previous[512]["std"], "actual provider std", atol=1e-9)
    require(
        set(provider["metrics"])
        == {
            "mean_rms",
            "mean_max",
            "std_rms",
            "std_max",
            "solve_max",
            "log_probability_error",
        }
        and all(
            0 <= finite_number(value, name) <= 1e-9
            for name, value in provider["metrics"].items()
        ),
        "actual provider numerical gate failed",
    )
    return {
        "schema": "kernelgecp.covariance-audit.v1",
        "status": "passed",
        "source_revision": revision,
        "cases": summaries,
        "real_data": True,
        "actual_provider": True,
        "independent_full_kernel_cg_recomputed_once_per_size": True,
    }


def audit_single_case(
    case_path: Path, frozen_path: Path, data_path: Path, revision: str
):
    """Small real CI evidence, explicitly insufficient for full EI acceptance."""
    receipt, frozen = read_json(case_path), read_json(frozen_path)
    benchmark.validate_frozen(frozen, canonical=False)
    bind_source_to_git(frozen["provenance"])
    bind_source_to_git(receipt["provenance"], revision)
    require(
        receipt["provenance"]["source_hashes"] == frozen["provenance"]["source_hashes"],
        "generation Git blobs differ from calibration",
    )
    arrays = _arrays(
        artifact(
            case_path.parent, receipt["artifact"]["path"], receipt["artifact"]["sha256"]
        )
    )
    data = benchmark.load_data(data_path)
    require(
        receipt["size"] <= 2048, "single-case CI audit is bounded to small dense cases"
    )
    summary = audit_case(receipt, arrays, frozen, data, revision)
    return {
        "schema": "kernelgecp.covariance-small-audit.v1",
        "status": "passed",
        "source_revision": revision,
        "case": summary,
        "scope": "small_real_CI_case_not_full_EI_campaign_acceptance",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--manifest", type=Path)
    target.add_argument("--single-case", type=Path)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--provider", type=Path)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.single_case is not None:
            result = audit_single_case(
                args.single_case, args.frozen, args.data, args.source_revision
            )
        else:
            if args.provider is None:
                parser.error("--provider is required for a full campaign audit")
            result = audit_campaign(
                args.manifest,
                args.frozen,
                args.data,
                args.provider,
                args.source_revision,
            )
    except Exception as error:
        benchmark.write_json(
            args.output,
            {
                "schema": "kernelgecp.covariance-audit.v1",
                "status": "failed",
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise SystemExit(1) from error
    benchmark.write_json(args.output, result)


if __name__ == "__main__":
    main()
