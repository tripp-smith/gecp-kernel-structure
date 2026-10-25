"""Run ``uv run python -m examples.precision_bridge`` for a finite-query audit.

QR below independently checks the interval certificate's algebraic reference.
Physical pole errors are reported separately; they are sampled evidence.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any

import mpmath as mp
import numpy as np

from kernelgecp import FermionicBasis
from kernelgecp.precision import certify_fit_precision
from kernelgecp.precision_grid import certify_grid_pivot_precision


def audited_case(
    *,
    beta: float = 3.7,
    wmax: float = 2.4,
    tolerance: float = 1e-9,
    training: str = "uniform",
    precision_bits: int = 256,
    max_work: int = 2**22,
) -> dict[str, Any]:
    """Compare actual complex matrix fits to a 384-bit independent QR solve."""
    basis = FermionicBasis(
        beta=beta, wmax=wmax, tol=tolerance, grid_order=6, precision="reference"
    )
    if training == "uniform":
        tau = np.linspace(0, beta, max(40, basis.rank + 12))
    elif training == "selected_nodes":
        tau = basis.tau_nodes
    elif training == "resolved":
        from experiments.dlr_physics import physical_tau_meshes

        tau = physical_tau_meshes(beta=beta, wmax=wmax, order=4)[0]
    else:
        raise ValueError("training must be uniform, selected_nodes or resolved")
    energies = wmax * np.array([-0.31, 0.23, 0.71])
    weights = np.array([0.25, 0.5, 0.25])
    amplitudes = np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])
    # Stable scalar pole kernel; these actual binary64 observations define
    # the exact-data LS target. The physical model is a separate comparison.
    exponents = np.where(
        energies >= 0, -tau[:, None] * energies, (beta - tau[:, None]) * energies
    )
    positive_kernel = np.exp(exponents) / (1 + np.exp(-beta * np.abs(energies)))
    values = -(positive_kernel @ weights)[:, None, None] * amplitudes
    fit = basis.fit_tau(values, tau)
    queries = beta * np.array([0.071, 0.23, 0.41, 0.63, 0.91])
    modes = np.array([-17, -2, 0, 4, 31])
    certificate = certify_fit_precision(
        fit,
        values,
        tau,
        query_tau=queries,
        matsubara_indices=modes,
        precision_bits=precision_bits,
        max_work=max_work,
    )
    grid = certify_grid_pivot_precision(
        basis,
        queries,
        basis.cutoff * np.array([-0.91, -0.41, 0.03, 0.37, 0.83]),
        precision_bits=precision_bits,
        max_work=max_work,
    )
    metadata = {
        "beta": beta,
        "wmax": wmax,
        "rank": basis.rank,
        "training": training,
        "selection_evaluation_bits": basis.diagnostics.evaluation_precision_bits,
        "actual_least_squares": asdict(basis.prepare_tau(tau).diagnostics),
        "fit_certificate": certificate.to_dict(),
        "grid_certificate": asdict(grid),
        "scope": "finite_query_arithmetic_reference;physical_errors_are_sampled",
    }
    if not certificate.certified:
        return {**metadata, "independent_qr_performed": False, "queries": {}}
    assert certificate.coefficients is not None
    assert certificate.time is not None and certificate.matsubara is not None
    context = mp.mp.clone()
    context.prec = 384
    exact_beta = context.mpf(beta)

    def scalar(value: complex | float) -> Any:
        return context.mpc(float(np.real(value)), float(np.imag(value)))

    def kernel(time: float, omega: float) -> Any:
        t, w = context.mpf(float(time)) / exact_beta, context.mpf(float(omega))
        return context.exp(-t * w) / (1 + context.exp(-w))

    core = context.matrix(
        [[kernel(t, w) for w in basis.omega_nodes] for t in basis.tau_nodes]
    )
    core_q, core_r = context.qr(core.T)

    def qr_solve(q: Any, r: Any, rhs: Any) -> Any:
        projected = q.T * rhs
        result = context.matrix(r.cols, 1)
        for i in range(r.cols - 1, -1, -1):
            result[i] = (
                projected[i]
                - context.fsum(r[i, j] * result[j] for j in range(i + 1, r.cols))
            ) / r[i, i]
        return result

    def row(time: float) -> Any:
        return qr_solve(
            core_q, core_r, context.matrix([kernel(time, w) for w in basis.omega_nodes])
        )

    training_design = context.matrix([list(row(float(t))) for t in tau])
    training_q, training_r = context.qr(training_design)
    flattened = values.reshape(len(tau), -1)
    coefficients = [
        qr_solve(
            training_q,
            training_r,
            context.matrix([scalar(v) for v in flattened[:, column]]),
        )
        for column in range(flattened.shape[1])
    ]
    actual_coefficients = fit.coefficients.reshape(basis.rank, -1)
    coefficient_bound = certificate.coefficients.actual_error_bound.reshape(
        basis.rank, -1
    )
    max_coefficient_error = context.zero
    for column, solution in enumerate(coefficients):
        for i in range(basis.rank):
            error = abs(solution[i] - scalar(actual_coefficients[i, column]))
            if error > context.mpf(float(coefficient_bound[i, column])):
                raise AssertionError("independent QR exceeds coefficient enclosure")
            max_coefficient_error = max(max_coefficient_error, error)

    query_records = {}
    for name, points, bounds in (
        ("time", queries, certificate.time),
        ("matsubara", modes, certificate.matsubara),
    ):
        if name == "time":
            designs = [row(float(t)) for t in points]
            actual = fit.evaluate_tau(points).reshape(len(points), -1)
        else:
            designs = [
                qr_solve(
                    core_q,
                    core_r,
                    context.matrix(
                        [
                            -exact_beta
                            / (context.j * context.pi * (2 * int(n) + 1) - float(w))
                            for w in basis.omega_nodes
                        ]
                    ),
                )
                for n in points
            ]
            actual = fit.evaluate_matsubara(points).reshape(len(points), -1)
        direct = bounds.direct_error_bound.reshape(len(points), -1)
        decomposed = bounds.total_error_bound.reshape(len(points), -1)
        lower_real = bounds.ideal_prediction.real_lower.reshape(len(points), -1)
        upper_real = bounds.ideal_prediction.real_upper.reshape(len(points), -1)
        lower_imag = bounds.ideal_prediction.imag_lower.reshape(len(points), -1)
        upper_imag = bounds.ideal_prediction.imag_upper.reshape(len(points), -1)
        max_error, max_physical_error = context.zero, context.zero
        for i, design in enumerate(designs):
            for j, solution in enumerate(coefficients):
                ideal = context.fsum(design[k] * solution[k] for k in range(basis.rank))
                error = abs(ideal - scalar(actual[i, j]))
                if not (
                    error <= context.mpf(float(direct[i, j]))
                    and error <= context.mpf(float(decomposed[i, j]))
                    and context.mpf(float(lower_real[i, j]))
                    <= ideal.real
                    <= context.mpf(float(upper_real[i, j]))
                    and context.mpf(float(lower_imag[i, j]))
                    <= ideal.imag
                    <= context.mpf(float(upper_imag[i, j]))
                ):
                    raise AssertionError(
                        "independent QR exceeds query enclosure/budget"
                    )
                if name == "time":
                    truth = -context.fsum(
                        context.mpf(float(weight))
                        * context.exp(-context.mpf(float(points[i])) * float(energy))
                        / (1 + context.exp(-exact_beta * float(energy)))
                        for weight, energy in zip(weights, energies, strict=True)
                    )
                else:
                    truth = context.fsum(
                        context.mpf(float(weight))
                        / (
                            context.j
                            * context.pi
                            * (2 * int(points[i]) + 1)
                            / exact_beta
                            - float(energy)
                        )
                        for weight, energy in zip(weights, energies, strict=True)
                    )
                physical_error = abs(
                    scalar(actual[i, j]) - truth * scalar(amplitudes.flat[j])
                )
                max_error, max_physical_error = (
                    max(max_error, error),
                    max(max_physical_error, physical_error),
                )
        query_records[name] = {
            "independent_qr_error": context.nstr(max_error, 35),
            "direct_error_upper": float(np.max(direct)),
            "decomposed_error_upper": float(np.max(decomposed)),
            "sampled_physical_pole_error": context.nstr(max_physical_error, 35),
        }
    return {
        **metadata,
        "independent_qr_performed": True,
        "independent_qr_bits": context.prec,
        "independent_coefficient_error": context.nstr(max_coefficient_error, 35),
        "queries": query_records,
    }


def main() -> None:
    result = audited_case()
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "fit_certificate"},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
