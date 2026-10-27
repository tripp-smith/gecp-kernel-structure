"""Run an actual isolated CPU GPyTorch operator and probability calculation."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import re
from pathlib import Path


def run(source_revision: str) -> dict[str, object]:
    """Compare real provider solves/log probabilities with independent SciPy."""
    import gpytorch
    import numpy as np
    import torch
    from linear_operator.operators import LowRankRootLinearOperator
    from scipy.linalg import cho_factor, cho_solve

    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "torch", "gpytorch", "linear-operator")
    }
    versions["python"] = platform.python_version()
    expected = {
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "mpmath": "1.3.0",
        "torch": "2.14.1",
        "gpytorch": "1.15.2",
        "linear-operator": "0.6.1",
    }
    if any(
        versions[name].split("+")[0] != version for name, version in expected.items()
    ):
        raise RuntimeError(
            f"provider versions differ from the supported lock: {versions}"
        )
    if platform.python_version_tuple()[:2] != ("3", "12"):
        raise RuntimeError("the locked GP runtime requires Python 3.12")
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    root = np.array([[1.0, 0.2], [0.4, -0.3], [-0.2, 0.8], [0.7, 0.1]])
    shift = np.array([0.1, 0.2, 0.15, 0.3])
    rhs = np.array([[0.3, 0.4], [-0.8, 0.2], [1.1, -0.6], [0.7, 0.5]])
    observed = rhs[:, 0]
    dense = root @ root.T + np.diag(shift)
    chol = cho_factor(dense, lower=True)
    expected_solve = cho_solve(chol, rhs)
    expected_logdet = float(2 * np.log(np.diag(chol[0])).sum())
    expected_logprob = float(
        -0.5
        * (observed @ expected_solve[:, 0] + expected_logdet + 4 * np.log(2 * np.pi))
    )
    operator = LowRankRootLinearOperator(torch.from_numpy(root)).add_diagonal(
        torch.from_numpy(shift)
    )
    solved = operator.solve(torch.from_numpy(rhs)).numpy()
    logdet = float(operator.logdet())
    distribution = gpytorch.distributions.MultivariateNormal(
        torch.zeros(4, dtype=torch.float64), operator
    )
    logprob = float(distribution.log_prob(torch.from_numpy(observed)))
    metrics = {
        "solve_max_error": float(np.max(np.abs(solved - expected_solve))),
        "logdet_error": abs(logdet - expected_logdet),
        "log_probability_error": abs(logprob - expected_logprob),
    }
    if not all(np.isfinite(value) and value <= 1e-10 for value in metrics.values()):
        raise RuntimeError(f"CPU provider numerical smoke failed: {metrics}")
    if operator.device.type != "cpu" or operator.dtype != torch.float64:
        raise RuntimeError("GP provider smoke must execute CPU float64 arithmetic")
    return {
        "schema": "kernelgecp.provider-smoke.v1",
        "provider": "gp",
        "status": "passed",
        "source_revision": source_revision,
        "checks": {"imports": True, "numerical": True, "cpu_float64": True},
        "versions": versions,
        "metrics": metrics,
        "platform": {"system": platform.system(), "machine": platform.machine()},
        "scope": "actual_provider_small_operator_not_large_data_or_GP_phase_acceptance",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    if re.fullmatch(r"[0-9a-f]{40}", args.source_revision) is None:
        parser.error("--source-revision must be a full lowercase Git commit SHA")
    status = 0
    try:
        result = run(args.source_revision)
    except Exception as error:
        result = {
            "schema": "kernelgecp.provider-smoke.v1",
            "provider": "gp",
            "status": "failed",
            "source_revision": args.source_revision,
            "error": f"{type(error).__name__}: {error}",
        }
        status = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(status)


if __name__ == "__main__":
    main()
