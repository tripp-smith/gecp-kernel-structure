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
from kernelgecp.kernels import FermionicKernel
from kernelgecp.precision import _certificate_hash, certify_fit_precision
from kernelgecp.precision_grid import _grid_content_hash, certify_grid_pivot_precision


def audited_case(
    *,
    beta: float = 3.7,
    wmax: float = 2.4,
    tolerance: float = 1e-9,
    training: str = "uniform",
    precision_bits: int = 256,
    max_work: int = 2**22,
    audit_precision_bits: int = 1024,
) -> dict[str, Any]:
    """Compare actual complex matrix fits to an independent QR solve."""
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
    arrays = {
        "basis_parameters": np.array([basis.beta, basis.wmax, basis.cutoff]),
        "basis_tau_nodes": basis.tau_nodes,
        "basis_omega_nodes": basis.omega_nodes,
        "training_tau": tau,
        "observations": values,
        "coefficients": fit.coefficients,
        "query_tau": queries,
        "matsubara_indices": modes,
        "actual_training_design": basis.design_tau(tau),
        "actual_time_design": basis.design_tau(queries),
        "actual_frequency_design": basis.design_matsubara(modes),
        "actual_time_output": fit.evaluate_tau(queries),
        "actual_frequency_output": fit.evaluate_matsubara(modes),
    }
    grid_omega = basis.cutoff * np.array([-0.91, -0.41, 0.03, 0.37, 0.83])
    kernel_binary64 = FermionicKernel(basis.cutoff)
    with np.errstate(all="ignore"):
        arrays["reconstructed_binary64_core"] = kernel_binary64(
            basis.tau_nodes[:, None] / basis.beta, basis.omega_nodes[None, :]
        )
        arrays["actual_grid_residual"] = kernel_binary64(
            queries[:, None] / basis.beta, grid_omega[None, :]
        ) - arrays["actual_time_design"] @ kernel_binary64(
            basis.tau_nodes[:, None] / basis.beta, grid_omega[None, :]
        )
    if (
        _grid_content_hash(basis, queries, grid_omega, arrays["actual_grid_residual"])
        != grid.content_sha256
    ):
        raise AssertionError("actual rounded grid changed after certification")
    recorded_hashes = dict(certificate.input_hashes)
    for name, array in arrays.items():
        if (
            name in recorded_hashes
            and _certificate_hash(array) != recorded_hashes[name]
        ):
            raise AssertionError(f"frozen actual input/output changed: {name}")

    def packed(array: Any) -> dict[str, Any]:
        value = np.asarray(array)
        return {
            "dtype": value.dtype.str,
            "shape": value.shape,
            "real": value.real.tolist(),
            "imag": value.imag.tolist() if np.iscomplexobj(value) else None,
            "sha256": _certificate_hash(value),
        }

    metadata = {
        "beta": beta,
        "wmax": wmax,
        "rank": basis.rank,
        "training": training,
        "selection_evaluation_bits": basis.diagnostics.evaluation_precision_bits,
        "actual_least_squares": asdict(basis.prepare_tau(tau).diagnostics),
        "fit_certificate": certificate.to_dict(),
        "grid_certificate": asdict(grid),
        "inputs": {name: packed(array) for name, array in arrays.items()},
        "scope": "finite_query_arithmetic_reference;physical_errors_are_sampled",
    }
    if not certificate.certified:
        return {**metadata, "independent_qr_performed": False, "queries": {}}
    assert certificate.coefficients is not None
    assert certificate.time is not None and certificate.matsubara is not None
    context = mp.mp.clone()
    context.prec = audit_precision_bits
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

    flattened = values.reshape(len(tau), -1)
    if training == "selected_nodes":
        # Once the core is nonsingular its exact cardinal design at the same
        # nodes is I. This exact algebraic oracle avoids interpreting QR's
        # tiny nonzero roundoff as an error in an exactly reproduced datum.
        coefficients = [
            context.matrix([scalar(v) for v in flattened[:, column]])
            for column in range(flattened.shape[1])
        ]
    else:
        training_design = context.matrix([list(row(float(t))) for t in tau])
        training_q, training_r = context.qr(training_design)
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
                raise AssertionError(
                    f"independent QR exceeds coefficient enclosure at {(i, column)}: "
                    f"error={context.nstr(error, 25)}, "
                    f"bound={coefficient_bound[i, column]}"
                )
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
        "independent_coefficient_oracle": (
            "exact_selected_node_reproduction" if training == "selected_nodes" else "QR"
        ),
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
