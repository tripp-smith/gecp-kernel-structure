"""Build a column-only covariance operator and solve its shifted approximation."""

from __future__ import annotations

import argparse
import json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gpytorch", action="store_true", help="check optional CPU parity"
    )
    args = parser.parse_args()
    import numpy as np

    from kernelgecp.kernel_cholesky import KernelCholeskyFactor

    x = np.linspace(-2, 2, 256)
    factor = KernelCholeskyFactor.from_columns(
        np.ones(len(x)),
        lambda j: np.exp(-(((x - x[j]) / 0.6) ** 2) / 2),
        max_rank=64,
        atol=1e-10,
    )
    rhs = np.sin(x)
    noise = 0.05
    result = factor.shifted(noise).solve_with_diagnostics(rhs)
    operator = factor.as_linear_operator(noise)
    residual = float(
        np.linalg.norm(operator @ result.solution - rhs) / np.linalg.norm(rhs)
    )
    if residual > 1e-10:
        raise RuntimeError(
            "approximate shifted-system residual exceeds declared budget"
        )
    report = {
        "rank": factor.rank,
        "stop_reason": factor.stop_reason,
        "retained_factor_bytes": factor.storage_nbytes,
        "approximate_system_relative_residual": residual,
        "scope": "finite RBF example; not original-kernel error or GP calibration",
    }
    if args.gpytorch:
        import torch

        actual = factor.to_gpytorch(noise).solve(torch.tensor(rhs, dtype=torch.float64))
        error = float(np.max(np.abs(actual.numpy() - result.solution)))
        if error > 1e-10:
            raise RuntimeError("identical-factor CPU GPyTorch parity failed")
        report["gpytorch_max_difference"] = error
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
