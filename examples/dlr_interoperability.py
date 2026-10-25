"""Same physical workload and solver using GECP, pydlr and sparse-ir.

Run from the checkout with optional providers installed, for example::

    uv run --with pydlr==1.0.1 --with sparse-ir==2.1.6 \
        python -m examples.dlr_interoperability

This is an executable tutorial, not the frozen matched-accuracy benchmark.
It deliberately uses one fixed native eps; those eps meanings differ between
providers. See experiments.dlr_application_benchmark for calibration and timing.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Any

import numpy as np

from experiments.dlr_adapters import build_adapter
from experiments.dlr_physics import (
    BetheSolverPlan,
    frequency_index_meshes,
    matrix_amplitudes,
    matsubara_green_function,
    physical_tau_meshes,
    semicircle_matsubara,
    semicircle_tau,
    solver_validation_index_mesh,
)
from experiments.dlr_workloads import green_function, make_pole_workloads
from kernelgecp import FermionicBasis, LeastSquaresPlan

METHODS = ("kernelgecp", "pydlr", "sparse_ir")


def run_demo(methods: tuple[str, ...] = METHODS) -> dict[str, Any]:
    """Fit identical complex matrix observations and run the same solver loop."""
    if (
        not methods
        or len(set(methods)) != len(methods)
        or any(name not in METHODS for name in methods)
    ):
        raise ValueError("methods must be distinct supported provider names")
    beta, wmax, native_eps = 10.0, 2.0, 1e-10
    training, _, heldout = physical_tau_meshes(beta=beta, wmax=wmax, order=8)
    _, heldout_n = frequency_index_meshes()
    solver_heldout_n = solver_validation_index_mesh()
    workload = make_pole_workloads(beta=beta, wmax=wmax)[1]
    amplitude = matrix_amplitudes()
    observations = (
        green_function(workload, training, beta=beta, wmax=wmax)[:, None, None]
        * amplitude
    )
    expected_tau = (
        green_function(workload, heldout, beta=beta, wmax=wmax)[:, None, None]
        * amplitude
    )
    expected_iw = (
        matsubara_green_function(workload, heldout_n, beta=beta, wmax=wmax)[
            :, None, None
        ]
        * amplitude
    )
    records = []
    for method in methods:
        basis = (
            FermionicBasis(beta, wmax, native_eps, grid_order=8, max_rank=128)
            if method == "kernelgecp"
            else build_adapter(method, beta, wmax, native_eps)
        )
        # The common fitter consumes provider design matrices. It does not
        # silently apply a Hermitian transpose/conjugation to complex batches.
        fit_plan = LeastSquaresPlan(basis.design_tau(training), rcond=1e-14)
        coefficients = fit_plan.fit(observations)
        predicted_tau = np.tensordot(
            basis.design_tau(heldout), coefficients, axes=(1, 0)
        )
        predicted_iw = np.tensordot(
            basis.design_matsubara(heldout_n), coefficients, axes=(1, 0)
        )
        solver = BetheSolverPlan(basis, hopping=1.0)
        solved = solver.solve()
        solver_tau = basis.design_tau(heldout) @ solved.coefficients
        solver_iw = basis.design_matsubara(solver_heldout_n) @ solved.coefficients
        errors = {
            "matrix_time": float(np.max(np.abs(predicted_tau - expected_tau))),
            "matrix_frequency": float(np.max(np.abs(predicted_iw - expected_iw))),
            "solver_time": float(
                np.max(
                    np.abs(solver_tau - semicircle_tau(heldout, beta=beta, hopping=1))
                )
            ),
            "solver_frequency": float(
                np.max(
                    np.abs(
                        solver_iw
                        - semicircle_matsubara(solver_heldout_n, beta=beta, hopping=1)
                    )
                )
            ),
        }
        errors["solver_observed_grid_independent_truth"] = float(
            np.max(
                np.abs(
                    solved.frequency_values
                    - semicircle_matsubara(solved.indices, beta=beta, hopping=1)
                )
            )
        )
        records.append(
            {
                "method": method,
                "rank": basis.rank,
                "errors": errors,
                "accuracy_passed": all(error <= 1e-8 for error in errors.values()),
                "solver_converged": solved.converged,
                "solver_stop_reason": solved.stop_reason,
                "solver_iterations": len(solved.history),
                "last_iteration": asdict(solved.history[-1]),
            }
        )
    return {
        "beta": beta,
        "wmax": wmax,
        "native_eps": native_eps,
        "comparison_scope": (
            "fixed_eps_tutorial_not_matched_accuracy_performance_evidence"
        ),
        "certification_status": "not_certified_sampled_physical_checks",
        "solver_scope": "noninteracting_half_filled_bethe_not_interacting_DMFT",
        "observed_grid_error_scope": "independent_truth_after_solve_not_a_heldout_mesh",
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", choices=METHODS, nargs="+", default=list(METHODS))
    args = parser.parse_args()
    result = run_demo(tuple(args.methods))
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    if not all(
        record["accuracy_passed"] and record["solver_converged"]
        for record in result["records"]
    ):
        raise SystemExit("a physical accuracy or convergence check failed")


if __name__ == "__main__":
    main()
