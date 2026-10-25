"""Standalone upstream-recipe candidate; needs only NumPy and sparse-ir.

Run ``python -m examples.sparse_ir_physical_recipe`` from the checkout with
sparse-ir installed. This is not a calibrated comparative performance test.
"""

from __future__ import annotations

import importlib.metadata
import json

import numpy as np


def run_recipe() -> dict[str, str | float | int | bool]:
    import sparse_ir

    beta = 10.0
    basis = sparse_ir.FiniteTempBasis("F", beta=beta, wmax=1.0, eps=1e-10)
    tau = np.linspace(0.0, beta, 121)
    amplitude = np.array([[1.0, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])
    # A normalized spectral atom at zero: G(tau)=-1/2, G(iw)=1/iw.
    values = -0.5 * np.ones((tau.size, 1, 1)) * amplitude
    sampling = sparse_ir.TauSampling(basis, sampling_points=tau)
    coefficients = sampling.fit(values, axis=0)
    heldout_tau = (np.arange(73) + 0.5) * beta / 73
    predicted_tau = np.tensordot(basis.u(heldout_tau).T, coefficients, axes=(1, 0))
    n = np.arange(-32, 32)  # Ordinary indices: sparse-ir instead expects 2*n+1.
    predicted_iw = np.tensordot(basis.uhat(2 * n + 1).T, coefficients, axes=(1, 0))
    expected_iw = amplitude / (1j * (2 * n[:, None, None] + 1) * np.pi / beta)
    time_error = float(np.max(np.abs(predicted_tau + 0.5 * amplitude)))
    frequency_error = float(np.max(np.abs(predicted_iw - expected_iw)))
    return {
        "sparse_ir_version": importlib.metadata.version("sparse-ir"),
        "rank": int(basis.size),
        "time_error": time_error,
        "frequency_error": frequency_error,
        "accuracy_passed": max(time_error, frequency_error) <= 1e-8,
        "scope": "sampled_physical_check_not_a_uniform_or_rounding_certificate",
    }


if __name__ == "__main__":
    result = run_recipe()
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    if not result["accuracy_passed"]:
        raise SystemExit("physical accuracy check failed")
