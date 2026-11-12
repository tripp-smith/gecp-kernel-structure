"""Bounded saved-factor replay using SciPy and actual CPU GPyTorch."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import math
import struct
import sys
import zipfile
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

OUTER_SHA256 = "b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4"
PROVIDER_KEYS = {
    "root",
    "mean",
    "std",
    "expected_mean",
    "expected_std",
    "solution",
    "rhs",
}
CASE_KEYS = {
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
LIMIT = 32 * 1024**2


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe(path, maximum):
    path = Path(path).absolute()
    require(".." not in path.parts, "noncanonical path")
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        require(not current.is_symlink(), "symlink path")
    require(
        path.is_file()
        and path.stat().st_nlink == 1
        and 0 < path.stat().st_size <= maximum,
        "unsafe or oversized input",
    )
    return path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    def pairs(entries):
        result = {}
        for key, value in entries:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    def invalid(value):
        raise ValueError("nonfinite JSON: " + value)

    result = json.loads(
        safe(path, 2_000_000).read_text(),
        object_pairs_hook=pairs,
        parse_constant=invalid,
    )
    require(type(result) is dict, "JSON object required")
    return result


def read_arrays(path, keys):
    safe(path, LIMIT)
    with zipfile.ZipFile(path) as archive:
        names = [item.filename for item in archive.infolist()]
        require(
            len(names) == len(set(names)) and set(names) == {k + ".npy" for k in keys},
            "NPZ keyset mismatch",
        )
        require(
            sum(item.file_size for item in archive.infolist()) <= LIMIT,
            "NPZ expanded size limit",
        )
        for item in archive.infolist():
            with archive.open(item) as stream:
                version = np.lib.format.read_magic(stream)
                require(version in ((1, 0), (2, 0)), "unsupported NPY version")
                length_bytes = 2 if version == (1, 0) else 4
                length_position = stream.tell()
                encoded_length = stream.read(length_bytes)
                require(
                    len(encoded_length) == length_bytes, "truncated NPY header length"
                )
                header_length = struct.unpack(
                    "<H" if length_bytes == 2 else "<I", encoded_length
                )[0]
                require(header_length <= 10000, "NPY header size limit")
                stream.seek(length_position)
                header = (
                    np.lib.format.read_array_header_1_0
                    if version == (1, 0)
                    else np.lib.format.read_array_header_2_0
                )
                shape, _, dtype = header(stream, max_header_size=10000)
                require(
                    dtype.kind in "fiu" and 0 < dtype.itemsize <= LIMIT,
                    "invalid NPY dtype",
                )
                require(
                    len(shape) <= 2
                    and all(type(n) is int and 0 <= n <= LIMIT for n in shape),
                    "invalid NPY dimensions",
                )
                count = math.prod(shape)
                require(count <= LIMIT // dtype.itemsize, "NPY allocation limit")
                require(
                    stream.tell() + count * dtype.itemsize == item.file_size,
                    "NPY payload size mismatch",
                )
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in keys}
    for value in arrays.values():
        require(
            value.dtype.kind in "fiu" and np.isfinite(value).all(), "invalid NPZ array"
        )
    return arrays


def close(a, b, message, *, atol=1e-9):
    require(
        np.shape(a) == np.shape(b) and bool(np.allclose(a, b, atol=atol, rtol=0)),
        message,
    )


def calculate_reference(provider, case, split, noise, lengthscale):
    """Independently evaluate the already-selected, bounded saved factor."""
    root = case["root"]
    require(root.ndim == 2, "root dimensions")
    n, rank = root.shape
    x, q, y = split["x"], split["query_x"], split["y"]
    require(1 <= n <= 512 and 1 <= rank <= n and len(q) <= 512, "matrix budget")
    require(
        x.ndim == q.ndim == 2
        and x.shape[0] == n
        and x.shape[1] == q.shape[1]
        and y.shape == (n,),
        "split dimensions",
    )
    require(
        np.isfinite(noise)
        and noise > 0
        and np.isfinite(lengthscale)
        and lengthscale > 0,
        "invalid frozen controls",
    )
    pivots = case["pivots"]
    require(
        pivots.dtype.kind in "iu"
        and pivots.shape == (rank,)
        and len(np.unique(pivots)) == rank
        and np.all((pivots >= 0) & (pivots < n)),
        "invalid pivots",
    )
    core = root[pivots]
    require(
        np.all(np.diag(core) > 0)
        and np.allclose(core, np.tril(core), atol=1e-12, rtol=0),
        "invalid pivot core",
    )
    close(provider["root"], root, "provider root mismatch", atol=0)
    close(provider["rhs"], y, "provider rhs mismatch", atol=0)
    close(case["cg_rhs"], y, "case rhs mismatch", atol=0)
    diff = (q[:, None, :] - x[pivots][None, :, :]) / lengthscale
    kernel = np.exp(-0.5 * np.sum(diff**2, axis=2))
    features = solve_triangular(core, kernel.T, lower=True).T
    cross = features @ root.T
    matrix = root @ root.T + noise * np.eye(n)
    chol = cho_factor(matrix, lower=True)
    solution = cho_solve(chol, y)
    mean = cross @ solution
    variance = 1 + noise - np.sum(cross * cho_solve(chol, cross.T).T, axis=1)
    require(
        np.isfinite(variance).all() and np.all(variance > 0),
        "invalid predictive variance",
    )
    std = np.sqrt(variance)
    logprob = float(
        -0.5
        * (y @ solution + 2 * np.log(np.diag(chol[0])).sum() + n * np.log(2 * np.pi))
    )
    require(math.isfinite(logprob), "nonfinite log probability")
    for name, expected in (
        ("solution", solution),
        ("mean", mean),
        ("std", std),
        ("expected_mean", case["mean"]),
        ("expected_std", case["std"]),
    ):
        close(provider[name], expected, "provider " + name + " mismatch")
    close(case["approximate_solution"], solution, "case approximate solve mismatch")
    close(case["mean"], mean, "case predicted mean mismatch")
    close(case["std"], std, "case predicted std mismatch")
    return {
        "root": root,
        "rhs": y,
        "cross": cross,
        "solution": solution,
        "mean": mean,
        "std": std,
        "log_probability": logprob,
    }


def actual_provider(reference, noise):
    # No optional provider import occurs until an actual replay is requested.
    torch = importlib.import_module("torch")
    gpytorch = importlib.import_module("gpytorch")
    operators = importlib.import_module("linear_operator.operators")
    torch.set_num_threads(1)
    if torch.get_num_interop_threads() != 1:
        torch.set_num_interop_threads(1)
    root = torch.tensor(reference["root"], dtype=torch.float64, device="cpu")
    rhs = torch.tensor(reference["rhs"], dtype=torch.float64, device="cpu")
    cross = torch.tensor(reference["cross"], dtype=torch.float64, device="cpu")
    operator = operators.RootLinearOperator(root) + operators.DiagLinearOperator(
        torch.full((len(rhs),), noise, dtype=torch.float64, device="cpu")
    )
    require(
        operator.device.type == "cpu"
        and operator.dtype == torch.float64
        and not operator.requires_grad,
        "provider runtime mismatch",
    )
    with (
        torch.no_grad(),
        gpytorch.settings.fast_computations(
            covar_root_decomposition=False, log_prob=False, solves=False
        ),
        gpytorch.settings.max_cholesky_size(1024),
    ):
        solution = operator.solve(rhs)
        mean = cross @ solution
        variance = 1 + noise - (cross * operator.solve(cross.T).T).sum(dim=1)
        require(
            bool(torch.isfinite(variance).all()) and bool((variance > 0).all()),
            "provider variance invalid",
        )
        std = variance.sqrt()
        distribution = gpytorch.distributions.MultivariateNormal(
            torch.zeros_like(rhs), operator
        )
        logprob = float(distribution.log_prob(rhs))
    values = {"solution": solution.numpy(), "mean": mean.numpy(), "std": std.numpy()}
    errors = {
        key + "_max_abs_error": float(np.max(np.abs(value - reference[key])))
        for key, value in values.items()
    }
    errors["log_probability_abs_error"] = abs(logprob - reference["log_probability"])
    require(
        all(math.isfinite(v) and v <= 1e-9 for v in errors.values()),
        "actual provider numerical mismatch",
    )
    return {
        "errors": errors,
        "log_probability": logprob,
        "device": str(operator.device),
        "dtype": str(operator.dtype),
        "finite_counts": {k: int(np.isfinite(v).sum()) for k, v in values.items()},
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "torch", "gpytorch", "linear_operator")
        },
        "module_paths": {
            name: str(Path(importlib.import_module(name).__file__).resolve())
            for name in ("numpy", "scipy", "torch", "gpytorch", "linear_operator")
        },
    }


def replay(payload, source_root, data_path, output):
    payload, source_root, data_path, output = map(
        lambda p: Path(p).absolute(), (payload, source_root, data_path, output)
    )
    require(
        not output.exists()
        and output.parent.is_dir()
        and not output.parent.is_symlink(),
        "output must be new",
    )
    require(
        output != payload
        and payload not in output.parents
        and output != source_root
        and source_root not in output.parents
        and output != data_path,
        "output must be external",
    )
    paths = {
        "provider_json": safe(payload / "provider.json", 2_000_000),
        "provider_npz": safe(payload / "provider.npz", LIMIT),
        "case_json": safe(payload / "full/n512-repeat0.json", 2_000_000),
        "case_npz": safe(payload / "full/n512-repeat0.npz", LIMIT),
        "manifest": safe(payload / "full/campaign.json", 2_000_000),
        "calibration": safe(
            source_root / "data/covariance/calibration.json", 2_000_000
        ),
        "benchmark": safe(
            source_root / "experiments/covariance_benchmark.py", 2_000_000
        ),
        "data": safe(data_path, 20 * 1024**2),
    }
    before = {name: digest(path) for name, path in paths.items()}
    require(before["data"] == OUTER_SHA256, "UCI archive hash mismatch")
    provider, case, manifest, frozen = (
        read_json(paths[k])
        for k in ("provider_json", "case_json", "manifest", "calibration")
    )
    require(
        set(provider)
        == {
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
        },
        "provider receipt keyset",
    )
    require(
        provider["schema"] == "kernelgecp.covariance-provider.v1"
        and provider["status"] == "passed"
        and type(provider["size"]) is int
        and provider["size"] == 512
        and provider["provider"] == "gpytorch",
        "provider identity",
    )
    require(
        type(provider["checks"]) is dict
        and set(provider["checks"]) == {"cpu_float64", "frozen_factor", "numerical"}
        and all(value is True for value in provider["checks"].values()),
        "provider historical checks",
    )
    require(
        case.get("schema") == "kernelgecp.covariance-case.v1"
        and case.get("status") == "passed"
        and case.get("canonical") is True
        and case.get("size") == 512
        and case.get("repeat") == 0,
        "case identity",
    )
    require(
        provider["artifact"]
        == {"path": "provider.npz", "sha256": before["provider_npz"]}
        and case["artifact"]
        == {"path": "n512-repeat0.npz", "sha256": before["case_npz"]},
        "NPZ hash binding",
    )
    require(
        manifest.get("schema") == "kernelgecp.covariance-campaign.v1"
        and type(manifest.get("cases")) is list,
        "campaign identity",
    )
    rows = [
        item for item in manifest["cases"] if item.get("path") == "n512-repeat0.json"
    ]
    require(
        len(rows) == 1
        and rows[0].get("sha256") == before["case_json"]
        and type(rows[0].get("exit_code")) is int
        and rows[0]["exit_code"] == 0,
        "case manifest binding",
    )
    frozen_hash = hashlib.sha256(
        json.dumps(frozen, sort_keys=True).encode()
    ).hexdigest()
    require(
        provider["frozen_calibration_sha256"]
        == case["frozen_calibration_sha256"]
        == frozen_hash,
        "frozen calibration binding",
    )
    require(
        provider["dataset"] == case["dataset"] == frozen["dataset"],
        "dataset metadata mismatch",
    )
    require(
        provider["provenance"] == case["provenance"] == manifest["provenance"],
        "generation provenance mismatch",
    )
    metrics = provider["metrics"]
    require(
        type(metrics) is dict
        and set(metrics)
        == {
            "mean_rms",
            "mean_max",
            "std_rms",
            "std_max",
            "solve_max",
            "log_probability_error",
        }
        and all(
            type(v) in (int, float) and math.isfinite(v) and 0 <= v <= 1e-9
            for v in metrics.values()
        ),
        "historical provider metrics",
    )
    sys.path.insert(0, str(source_root))
    benchmark = importlib.import_module("experiments.covariance_benchmark")
    require(
        Path(benchmark.__file__).resolve() == paths["benchmark"],
        "benchmark imported outside source root",
    )

    # File existence and outer hash are checked before load_data; no download path.
    def no_network(*args, **kwargs):
        raise ValueError("network access forbidden during replay")

    benchmark.urllib.request.urlopen = no_network
    data = benchmark.load_data(paths["data"])
    require(
        benchmark.dataset_identity(data) == frozen["dataset"],
        "actual dataset identity mismatch",
    )
    split = benchmark.split_data(data, 512, "test")
    require(
        case["normalization"] == split["normalization"],
        "training normalization mismatch",
    )
    split_keys = {
        "train_ordinal_sha256",
        "query_ordinal_sha256",
        "train_first",
        "train_last",
        "query_first",
        "query_last",
    }
    require(
        type(case["split"]) is dict
        and set(case["split"]) == split_keys
        and all(case["split"][key] == split[key] for key in split_keys),
        "chronological split mismatch",
    )
    provider_arrays, case_arrays = (
        read_arrays(paths["provider_npz"], PROVIDER_KEYS),
        read_arrays(paths["case_npz"], CASE_KEYS),
    )
    require(
        case_arrays["root"].shape[0] == 512
        and case_arrays["root"].shape[1] == case["actual_rank"]
        and case["rank_budget"] == frozen["selected_ranks"]["512"],
        "case rank mismatch",
    )
    reference = calculate_reference(
        provider_arrays, case_arrays, split, frozen["noise"], frozen["lengthscale"]
    )
    actual = actual_provider(reference, frozen["noise"])
    receipt = {
        "schema": "kernelgecp.covariance-provider-replay.v1",
        "status": "passed",
        "scope": (
            "new saved-factor provider replay; no fitting, campaign, "
            "timing or release eligibility"
        ),
        "historical_reported_log_probability_error": metrics["log_probability_error"],
        "historical_log_probability_value_available": False,
        "scipy_reference_log_probability": reference["log_probability"],
        "actual_replay": actual,
        "input_sha256": before,
    }
    with output.open("x") as file:
        json.dump(receipt, file, sort_keys=True, indent=2, allow_nan=False)
        file.write("\n")
    if {
        name: digest(safe(path, 20 * 1024**2 if name == "data" else LIMIT))
        for name, path in paths.items()
    } != before:
        output.unlink()
        raise ValueError("input changed during replay")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        replay(args.payload, args.source_root, args.data, args.output)
    except Exception as error:
        print(
            f"provider replay failed: {type(error).__name__}: {error}", file=sys.stderr
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
