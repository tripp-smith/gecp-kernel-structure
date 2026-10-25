"""Standalone upstream-recipe candidate; needs only NumPy and sparse-ir.

Run ``python -m examples.sparse_ir_physical_recipe`` from the checkout with
sparse-ir installed. This is not a calibrated comparative performance test.
"""

from __future__ import annotations

import importlib.metadata
import json

import numpy as np
import numpy.typing as npt

type RealArray = npt.NDArray[np.float64]
type ComplexArray = npt.NDArray[np.complex128]


def _spectral_fixture() -> tuple[float, float, RealArray, ComplexArray]:
    """Two real poles with positive Hermitian residues summing to identity."""
    residue = np.array([[0.65, 0.15 + 0.2j], [0.15 - 0.2j, 0.35]])
    return (
        10.0,
        1.0,
        np.array([-0.3, 0.4]),
        np.array([residue, np.eye(2) - residue]),
    )


def _exact_tau(
    tau: RealArray, beta: float, energies: RealArray, residues: ComplexArray
) -> ComplexArray:
    kernel = -np.exp(-tau[:, None] * energies) / (1 + np.exp(-beta * energies))
    return np.einsum("tp,pab->tab", kernel, residues)


def _exact_iw(
    reduced_n: npt.NDArray[np.int64],
    beta: float,
    energies: RealArray,
    residues: ComplexArray,
) -> ComplexArray:
    z = 1j * reduced_n * np.pi / beta
    return np.einsum("np,pab->nab", 1 / (z[:, None] - energies), residues)


def run_recipe() -> dict[str, str | float | int | bool]:
    import sparse_ir

    beta, wmax, energies, residues = _spectral_fixture()
    basis = sparse_ir.FiniteTempBasis("F", beta=beta, wmax=wmax, eps=1e-10)
    tau_fit = np.linspace(0.0, beta, 121)
    values = _exact_tau(tau_fit, beta, energies, residues)
    sampling = sparse_ir.TauSampling(basis, sampling_points=tau_fit)
    coefficients = sampling.fit(values, axis=0)
    tau_check = beta * (np.arange(73) + 0.37) / 73
    if np.intersect1d(tau_fit, tau_check).size:
        raise RuntimeError("physical validation points must be disjoint from fitting")
    # Follow upstream notation: ordinary m maps to reduced fermionic n=2*m+1.
    m = np.arange(-32, 32)
    n = 2 * m + 1
    predicted_tau = np.einsum("lt,lab->tab", basis.u(tau_check), coefficients)
    predicted_iw = np.einsum("ln,lab->nab", basis.uhat(n), coefficients)
    expected_tau = _exact_tau(tau_check, beta, energies, residues)
    expected_iw = _exact_iw(n, beta, energies, residues)
    time_error = float(np.max(np.abs(predicted_tau - expected_tau)))
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
