"""Calibrate, freeze and validate interacting Bethe IPT representations.

Run inside the optional pinned TRIQS/DLR/IR environment. Final physical truth
is loaded only after solving; calibration and final reference files are separate.
This benchmark compares an approximate IPT model, not exact Hubbard physics.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
import tomllib
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("kernelgecp", "pydlr", "sparse_ir")
SOURCE_FILES = (
    "kernelgecp/dmft.py",
    "kernelgecp/lehmann.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "kernelgecp/integrations/triqs_dmft.py",
    "experiments/interacting_dmft.py",
    "experiments/dlr_adapters.py",
    "environments/locks/dlr-physics-pip.txt",
)
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def _normal(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _normal(dataclasses.asdict(value))
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return _normal(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {key: _normal(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_normal(item) for item in value]
    return value


def _write(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(_normal(value), sort_keys=True, indent=2, allow_nan=False) + "\n"
    )


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    data = tomllib.loads(path.read_text())
    if data.get("schema_version") != 1 or data.get("methods") != list(METHODS):
        raise ValueError("schema 1 requires all three declared comparison providers")
    for name in (
        "n_iw",
        "n_tau",
        "reference_n_iw",
        "reference_n_tau",
        "max_iterations",
        "grid_order",
        "max_rank",
        "precision_bits",
        "repeats",
    ):
        if type(data.get(name)) is not int or data[name] < 1:
            raise ValueError(f"{name} must be a positive integer")
    for name in (
        "beta",
        "hopping",
        "target",
        "calibration_target",
        "reference_refinement_target",
        "iteration_tolerance",
        "fit_rcond",
        "mixing",
        "selfenergy_cutoff_multiplier",
        "worker_timeout_seconds",
    ):
        value = data.get(name)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"{name} must be finite and positive")
    if (
        not 0 < data["calibration_target"] < data["target"] < 1
        or data["reference_refinement_target"] >= data["target"]
    ):
        raise ValueError(
            "calibration and reference budgets must be stricter than validation"
        )
    if (
        data["mixing"] > 1
        or data["reference_n_iw"] != 2 * data["n_iw"]
        or data["reference_n_tau"] != 2 * data["n_tau"] - 1
    ):
        raise ValueError("invalid mixing or incompatible doubled reference grids")
    for name in ("calibration_indices", "validation_indices"):
        values = data[name]
        limit = (
            data["n_iw"] if name == "calibration_indices" else data["reference_n_iw"]
        )
        if (
            not values
            or any(type(x) is not int or not -limit <= x < limit for x in values)
            or len(set(values)) != len(values)
        ):
            raise ValueError(
                "holdout indices must be distinct ordinary indices "
                "inside the declared reference grid"
            )
        if any(-18 <= x < 18 for x in values):
            raise ValueError("essential low frequencies must remain observed")
        if name == "validation_indices" and any(
            -data["n_iw"] <= x < data["n_iw"] for x in values
        ):
            raise ValueError(
                "final frequencies must be outside all primary observations"
            )
    if set(data["calibration_indices"]) & set(data["validation_indices"]):
        raise ValueError("calibration and final frequencies must be disjoint")
    if (
        data["calibration_tau_modulo"] != 1
        or data["validation_tau_rule"] != "gauss_legendre_64"
    ):
        raise ValueError("time validation must use the frozen off-grid rule")
    for name in ("green_cutoffs", "native_eps_candidates", "interactions"):
        values = data.get(name)
        if (
            not isinstance(values, list)
            or not values
            or any(
                isinstance(x, bool)
                or not isinstance(x, (int, float))
                or not math.isfinite(x)
                for x in values
            )
        ):
            raise ValueError(f"invalid numeric candidate list: {name}")
    if (
        any(x < 0 for x in data["interactions"])
        or any(x <= 0 or data["beta"] * x < 1 for x in data["green_cutoffs"])
        or any(not 0 < x < 1 for x in data["native_eps_candidates"])
    ):
        raise ValueError("invalid physical cutoff, interaction or native eps")
    if data["reference_backend"] not in ("python", "rust"):
        raise ValueError("reference backend must be explicit")
    if data["repeats"] != 3:
        raise ValueError(
            "canonical performance requires exactly three fresh repetitions"
        )
    if len(set(data["interactions"])) != len(data["interactions"]):
        raise ValueError("interactions must be unique")
    return data


def _config(data: dict[str, Any], interaction: float, *, refined: bool = False) -> Any:
    from kernelgecp.dmft import IPTConfig

    return IPTConfig(
        beta=data["beta"],
        hopping=data["hopping"],
        interaction=interaction,
        n_iw=data["reference_n_iw" if refined else "n_iw"],
        n_tau=data["reference_n_tau" if refined else "n_tau"],
        mixing=data["mixing"],
        max_iterations=data["max_iterations"],
        tolerance=data["iteration_tolerance"],
        rcond=data["fit_rcond"],
    )


def _fit_indices(data: dict[str, Any]) -> np.ndarray:
    excluded = data["calibration_indices"] + data["validation_indices"]
    return np.setdiff1d(
        np.arange(-data["n_iw"], data["n_iw"], dtype=np.int64), excluded
    )


def _reference(
    data: dict[str, Any], interaction: float, output: Path
) -> dict[str, Any]:
    from kernelgecp.integrations.triqs_dmft import DenseIPTPlan

    primary = DenseIPTPlan(_config(data, interaction)).solve(record_trajectory=True)
    refined = DenseIPTPlan(_config(data, interaction, refined=True)).solve()
    if not primary.converged or not refined.converged:
        raise RuntimeError(
            "dense native reference failed its independent stopping gate"
        )
    offset = data["reference_n_iw"] - data["n_iw"]
    differences = {
        "g_iw": float(
            np.max(
                np.abs(primary.g_iw - refined.g_iw[offset : offset + 2 * data["n_iw"]])
            )
        ),
        "sigma_iw": float(
            np.max(
                np.abs(
                    primary.sigma_iw
                    - refined.sigma_iw[offset : offset + 2 * data["n_iw"]]
                )
            )
        ),
        "g_tau": float(np.max(np.abs(primary.g_tau - refined.g_tau[::2]))),
        "sigma_tau": float(np.max(np.abs(primary.sigma_tau - refined.sigma_tau[::2]))),
    }
    if max(differences.values()) >= data["reference_refinement_target"]:
        raise RuntimeError(f"dense reference refinement failed: {differences}")
    np.savez_compressed(output / "trajectory.npz", bath=np.stack(primary.trajectory))
    quadrature_nodes = np.polynomial.legendre.leggauss(64)[0]
    final_tau = data["beta"] * (quadrature_nodes + 1) / 2
    if np.any(np.isin(final_tau, refined.tau)):
        raise RuntimeError("final quadrature unexpectedly overlaps the native mesh")
    interpolation_errors = {
        field: float(
            np.max(
                np.abs(
                    primary.evaluate_tau(final_tau, field=field)
                    - refined.evaluate_tau(final_tau, field=field)
                )
            )
        )
        for field in ("green", "sigma")
    }
    if max(interpolation_errors.values()) >= data["reference_refinement_target"]:
        raise RuntimeError("native off-grid interpolation did not refine")
    for role in ("calibration", "validation"):
        indices = np.asarray(data[role + "_indices"], dtype=np.int64)
        slots = indices + data["reference_n_iw"]
        time_slots = np.arange(1, data["reference_n_tau"], 4)
        query_tau = (
            np.asarray(refined.tau)[time_slots] if role == "calibration" else final_tau
        )
        np.savez_compressed(
            output / (role + ".npz"),
            indices=indices,
            tau=query_tau,
            g_iw=np.asarray(refined.g_iw)[slots],
            sigma_iw=np.asarray(refined.sigma_iw)[slots],
            g_tau=refined.evaluate_tau(query_tau, field="green"),
            sigma_tau=refined.evaluate_tau(query_tau, field="sigma"),
        )
    return {
        "status": "completed",
        "refinement_errors": differences,
        "off_grid_interpolation_refinement_errors": interpolation_errors,
        "primary_history": _normal(primary.history),
        "refined_history": _normal(refined.history),
        "scope": "dense_native_IPT_reference_not_exact_Hubbard_physics",
    }


def _basis(method: str, data: dict[str, Any], cutoff: float, eps: float) -> Any:
    if method == "kernelgecp":
        from kernelgecp import FermionicBasis

        return FermionicBasis(
            data["beta"],
            cutoff,
            eps,
            grid_order=data["grid_order"],
            max_rank=data["max_rank"],
            precision="reference",
            precision_bits=data["precision_bits"],
            reference_backend=data["reference_backend"],
        )
    from experiments.dlr_adapters import build_adapter

    return build_adapter(method, data["beta"], cutoff, eps)


def _predictions(plan: Any, update: Any, reference: Any) -> dict[str, np.ndarray]:
    return {
        "g_iw": plan.green_basis.design_matsubara(reference["indices"])
        @ update.g_coefficients,
        "sigma_iw": plan.selfenergy_basis.design_matsubara(reference["indices"])
        @ update.sigma_coefficients,
        "g_tau": plan.green_basis.design_tau(reference["tau"]) @ update.g_coefficients,
        "sigma_tau": plan.selfenergy_basis.design_tau(reference["tau"])
        @ update.sigma_coefficients,
    }


def _errors(plan: Any, update: Any, reference: Any) -> dict[str, float]:
    return {
        name: float(np.max(np.abs(values - reference[name])))
        for name, values in _predictions(plan, update, reference).items()
    }


def _fixture_errors(plan: Any, interaction: float) -> dict[str, float]:
    """Independent finite-model truth, separate from DMFT validation data."""
    from kernelgecp.lehmann import (
        FermionicPoleModel,
        hubbard_atom,
        ipt_self_energy_poles,
    )

    indices = plan.indices
    z = 1j * (2 * indices + 1) * np.pi / plan.config.beta
    weiss = FermionicPoleModel(plan.config.beta, [-1.3, 0.2, 1.1], [0.2, 0.5, 0.3])
    update = plan.update(z - 1 / weiss.evaluate_matsubara(indices))
    expected = ipt_self_energy_poles(weiss, interaction)
    atom = plan.update(np.zeros_like(z))
    atom_g = hubbard_atom(plan.config.beta, interaction)
    return {
        "asymmetric_sigma_iw": float(
            np.max(np.abs(update.sigma_iw - expected.evaluate_matsubara(indices)))
        ),
        "asymmetric_sigma_tau": float(
            np.max(np.abs(update.sigma_tau - expected.evaluate_tau(plan.tau)))
        ),
        "atomic_green_iw": float(
            np.max(np.abs(atom.g_iw - atom_g.evaluate_matsubara(indices)))
        ),
        "atomic_green_tau": float(
            np.max(np.abs(atom.g_tau - atom_g.evaluate_tau(plan.tau)))
        ),
    }


def _physics(
    data: dict[str, Any],
    interaction: float,
    indices: np.ndarray,
    green_iw: np.ndarray,
    sigma_iw: np.ndarray,
    green_endpoints: np.ndarray,
) -> dict[str, Any]:
    """Report sampled diagnostics without repairing fitted data or moments.

    Tail fits use scaled powers in two declared positive-frequency windows.
    Their difference is a sensitivity diagnostic, not a remainder bound.
    """
    positive = indices >= 0
    tails: dict[str, Any] = {}
    for name, values, expected in (
        ("green", green_iw, [1, 0, data["hopping"] ** 2 + interaction**2 / 4]),
        ("sigma", sigma_iw, [interaction**2 / 4, 0]),
    ):
        estimates = []
        for fraction in (0.5, 0.75):
            slots = indices >= int(data["n_iw"] * fraction)
            z = 1j * (2 * indices[slots] + 1) * np.pi / data["beta"]
            scale = float(np.min(np.abs(z)))
            powers = np.arange(1, 6)
            design = (scale / z[:, None]) ** powers
            coefficients, _, rank, singular = np.linalg.lstsq(
                design, values[slots], rcond=data["fit_rcond"]
            )
            moments = coefficients * scale**powers
            estimates.append(
                {
                    "first_ordinary_index": int(indices[slots][0]),
                    "moments": _normal(moments[: len(expected)]),
                    "effective_rank": int(rank),
                    "condition": float(singular[0] / singular[-1]),
                    "residual": float(
                        np.max(np.abs(design @ coefficients - values[slots]))
                    ),
                    "expected_moment_errors": np.abs(
                        moments[: len(expected)] - expected
                    ).tolist(),
                }
            )
        tails[name] = {
            "model_expected_moments": expected,
            "window_estimates": estimates,
        }
    return {
        "per_spin_density_from_beta_minus": _normal(-green_endpoints[-1]),
        "endpoint_mass": _normal(-sum(green_endpoints)),
        "density_error_from_half": float(abs(-green_endpoints[-1] - 0.5)),
        "endpoint_mass_error": float(abs(-sum(green_endpoints) - 1)),
        "positive_imag_green_max": float(np.max(green_iw[positive].imag)),
        "positive_imag_sigma_max": float(np.max(sigma_iw[positive].imag)),
        "particle_hole_green_real_max": float(np.max(np.abs(green_iw.real))),
        "particle_hole_sigma_real_max": float(np.max(np.abs(sigma_iw.real))),
        "tail_fits": tails,
        "scope": "sampled_density_causality_and_unconstrained_tail_diagnostics",
        "certification_status": "not_certified_no_tail_remainder_bound",
    }


def _peak_rss() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _dense_case(
    data: dict[str, Any], interaction: float, reference: Path, arrays_output: Path
) -> dict[str, Any]:
    from kernelgecp.integrations.triqs_dmft import DenseIPTPlan

    importlib.import_module("triqs.gfs")
    importlib.import_module("triqs.mesh")
    started = time.perf_counter()
    plan = DenseIPTPlan(_config(data, interaction))
    prepared = time.perf_counter()
    result = plan.solve()
    solved = time.perf_counter()
    loop_peak = _peak_rss()
    with np.load(reference, allow_pickle=False) as truth:
        errors = {}
        arrays = {
            name: getattr(result, name)
            for name in (
                "indices",
                "tau",
                "bath_iw",
                "g_iw",
                "sigma_iw",
                "weiss_iw",
                "weiss_tau",
                "sigma_tau",
                "g_tau",
            )
        }
        arrays.update(validation_indices=truth["indices"], validation_tau=truth["tau"])
        for name, field in (("g", "green"), ("sigma", "sigma")):
            for domain, coordinates in (("iw", "indices"), ("tau", "tau")):
                values = (
                    result.evaluate_matsubara(truth[coordinates], field=field)
                    if domain == "iw"
                    else result.evaluate_tau(truth[coordinates], field=field)
                )
                errors[name + "_" + domain] = float(
                    np.max(np.abs(values - truth[name + "_" + domain]))
                )
                arrays["validation_" + name + "_" + domain] = values
    np.savez_compressed(arrays_output, **arrays)
    return {
        "status": "completed",
        "method": "dense_triqs",
        "interaction": interaction,
        "converged": result.converged,
        "stop_reason": result.stop_reason,
        "errors": errors,
        "history": _normal(result.history),
        "native_tail": _normal(result.native_tail),
        "physics": _physics(
            data,
            interaction,
            result.indices,
            result.g_iw,
            result.sigma_iw,
            result.g_tau[[0, -1]],
        ),
        "measurements": {
            "preparation_seconds": prepared - started,
            "loop_seconds": solved - prepared,
            "application_seconds": solved - started,
            "evaluation_seconds": time.perf_counter() - solved,
            "peak_rss_bytes": _peak_rss(),
            "loop_peak_rss_bytes": loop_peak,
            "output_array_bytes": sum(
                getattr(result, x).nbytes
                for x in (
                    "g_iw",
                    "sigma_iw",
                    "weiss_iw",
                    "g_tau",
                    "sigma_tau",
                    "weiss_tau",
                )
            ),
        },
        "numerical_fingerprints": {
            name: hashlib.sha256(
                np.ascontiguousarray(getattr(result, name)).tobytes()
            ).hexdigest()
            for name in ("g_iw", "sigma_iw", "g_tau", "sigma_tau")
        },
    }


def _run_case(
    data: dict[str, Any],
    interaction: float,
    method: str,
    cutoff: float,
    eps: float,
    reference_path: Path,
    *,
    replay: Path | None = None,
    arrays_output: Path | None = None,
) -> dict[str, Any]:
    from kernelgecp.dmft import IPTPlan

    # Warm provider imports for every constructor timing. Native import is
    # warmed in _dense_case too; process wall includes all validation overhead.
    if method != "kernelgecp":
        importlib.import_module("sparse_ir" if method == "sparse_ir" else method)
    started = time.perf_counter()
    green = _basis(method, data, cutoff, eps)
    sigma = _basis(method, data, data["selfenergy_cutoff_multiplier"] * cutoff, eps)
    constructed = time.perf_counter()
    plan = IPTPlan(
        green,
        sigma,
        config=_config(data, interaction),
        frequency_fit_indices=_fit_indices(data),
    )
    prepared = time.perf_counter()
    result = plan.solve()
    solved = time.perf_counter()
    loop_peak = _peak_rss()
    if result.final_update is None:
        raise RuntimeError("no valid interacting update: " + str(result.stop_reason))
    # Deliberately load the role-specific truth only after the solve.
    with np.load(reference_path, allow_pickle=False) as truth:
        errors = _errors(plan, result.final_update, truth)
        validation_arrays = {
            "validation_" + name: values
            for name, values in _predictions(plan, result.final_update, truth).items()
        }
        validation_arrays.update(
            validation_indices=truth["indices"], validation_tau=truth["tau"]
        )
    evaluated = time.perf_counter()
    represented_g = green.design_matsubara(plan.indices) @ result.g_coefficients
    physics = _physics(
        data,
        interaction,
        plan.indices,
        represented_g,
        result.sigma_iw,
        green.design_tau([0, data["beta"]]) @ result.g_coefficients,
    )
    fixture_errors = _fixture_errors(plan, interaction)
    replay_record = None
    if replay is not None:
        from kernelgecp.integrations.triqs_dmft import DenseIPTPlan

        native = DenseIPTPlan(_config(data, interaction))
        maxima = {"g_iw": 0.0, "sigma_iw": 0.0, "g_tau": 0.0, "sigma_tau": 0.0}
        with np.load(replay, allow_pickle=False) as records:
            baths = records["bath"]
            if (
                baths.ndim != 2
                or baths.shape[1] != plan.indices.size
                or baths.shape[0] < 2
            ):
                raise ValueError(
                    "replay requires a complete nonempty native trajectory"
                )
            for bath in baths:
                candidate = plan.update(bath)
                expected = native.update(bath)
                for name, value in {
                    "g_iw": np.max(
                        np.abs(
                            plan.green_basis.design_matsubara(plan.indices)
                            @ candidate.g_coefficients
                            - expected.g_iw
                        )
                    ),
                    "sigma_iw": np.max(np.abs(candidate.sigma_iw - expected.sigma_iw)),
                    "g_tau": np.max(np.abs(candidate.g_tau - expected.g_tau)),
                    "sigma_tau": np.max(
                        np.abs(candidate.sigma_tau - expected.sigma_tau)
                    ),
                }.items():
                    maxima[name] = max(maxima[name], float(value))
        replay_record = {
            "scope": "common_native_bath_trajectory_not_independent_closed_loop",
            "errors": maxima,
            "updates_executed": int(baths.shape[0]),
        }
    if arrays_output is not None:
        np.savez_compressed(
            arrays_output,
            indices=plan.indices,
            frequency_fit_indices=plan.frequency_fit_indices,
            tau=plan.tau,
            bath_iw=result.bath_iw,
            g_iw=result.g_iw,
            represented_g_iw=represented_g,
            sigma_iw=result.sigma_iw,
            weiss_iw=result.weiss_iw,
            weiss_tau=result.weiss_tau,
            weiss_reflected_tau=result.final_update.weiss_reflected_tau,
            sigma_tau_product=result.sigma_tau_product,
            sigma_tau=result.sigma_tau,
            g_tau=result.g_tau,
            g_coefficients=result.g_coefficients,
            sigma_coefficients=result.sigma_coefficients,
            **validation_arrays,
        )
        if method == "kernelgecp":
            from kernelgecp import FermionicFit

            FermionicFit(
                green, result.g_coefficients, observation_domain="matsubara"
            ).save(arrays_output.parent / "green-fit.npz")
            FermionicFit(
                sigma, result.sigma_coefficients, observation_domain="tau"
            ).save(arrays_output.parent / "sigma-fit.npz")
    ended = time.perf_counter()
    return {
        "status": "completed",
        "method": method,
        "interaction": interaction,
        "green_cutoff": cutoff,
        "sigma_cutoff": data["selfenergy_cutoff_multiplier"] * cutoff,
        "native_eps": eps,
        "green_rank": green.rank,
        "sigma_rank": sigma.rank,
        "errors": errors,
        "converged": result.converged,
        "stop_reason": result.stop_reason,
        "history": _normal(result.history),
        "diagnostics": _normal(result.diagnostics),
        "physics": physics,
        "finite_fixture_errors": fixture_errors,
        "basis_metadata": {
            name: _normal(
                basis.diagnostics if method == "kernelgecp" else basis.metadata
            )
            for name, basis in (("green", green), ("sigma", sigma))
        },
        "replay": replay_record,
        "numerical_fingerprints": {
            name: hashlib.sha256(
                np.ascontiguousarray(getattr(result.final_update, name)).tobytes()
            ).hexdigest()
            for name in ("g_iw", "sigma_iw", "g_coefficients", "sigma_coefficients")
        },
        "measurements": {
            "construction_seconds": constructed - started,
            "preparation_seconds": prepared - constructed,
            "loop_seconds": solved - prepared,
            "evaluation_seconds": evaluated - solved,
            "validation_replay_seconds": ended - evaluated,
            "application_seconds": solved - started,
            "peak_rss_bytes": _peak_rss(),
            "loop_peak_rss_bytes": loop_peak,
            "coefficient_bytes": sum(
                getattr(result, name).nbytes
                for name in (
                    "g_coefficients",
                    "sigma_coefficients",
                    "weiss_coefficients",
                )
            ),
        },
        "convergence_scope": "unmixed_equation_defect_on_frequency_observations_only",
        "certification_status": "not_certified_sampled_IPT_model_comparison",
    }


def _worker(request: dict[str, Any], output: Path) -> None:
    data = request["config"]
    if request["operation"] == "reference":
        result = _reference(data, request["interaction"], output.parent)
    elif request["operation"] == "dense":
        result = _dense_case(
            data,
            request["interaction"],
            Path(request["reference"]),
            output.with_suffix(".npz"),
        )
    else:
        result = _run_case(
            data,
            request["interaction"],
            request["method"],
            request["cutoff"],
            request["eps"],
            Path(request["reference"]),
            replay=Path(request["replay"]) if request.get("replay") else None,
            arrays_output=output.with_suffix(".npz"),
        )
    _write(output, result)


def _child(request: dict[str, Any], directory: Path) -> dict[str, Any]:
    directory.mkdir(parents=True)
    request_path = directory / "request.json"
    output = directory / "result.json"
    _write(request_path, request)
    started = time.perf_counter()
    try:
        with (directory / "worker.log").open("w") as stream:
            process = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.interacting_dmft",
                    "--worker",
                    str(request_path),
                    "--worker-output",
                    str(output),
                ],
                cwd=ROOT,
                env={**os.environ, **dict.fromkeys(THREADS, "1")},
                stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=request["config"]["worker_timeout_seconds"],
                check=False,
            )
        if process.returncode != 0 or not output.is_file():
            return {
                "status": "failed",
                "returncode": process.returncode,
                "log_sha256": _sha(directory / "worker.log"),
            }
        result = json.loads(output.read_text())
        result.setdefault("measurements", {})["process_seconds"] = (
            time.perf_counter() - started
        )
        result["measurements"]["process_wall_scope"] = (
            "full_worker_including_validation_and_optional_replay"
        )
        result["measurements"]["loop_peak_rss_scope"] = (
            "process_high_water_after_solve_before_truth_and_replay"
        )
        result["measurements"]["peak_rss_scope"] = (
            "full_worker_including_validation_and_optional_native_replay"
        )
        result["measurements"]["application_time_scope"] = (
            "warm_provider_imports_construction_preparation_and_loop"
        )
        result["artifacts"] = {
            path.name: _sha(path) for path in directory.iterdir() if path.is_file()
        }
        return result
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "log_sha256": _sha(directory / "worker.log")}


def run(config_path: Path, output: Path, *, pilot: bool) -> dict[str, Any]:
    data = load_config(config_path)
    config_hash = _sha(config_path)
    source_hashes = {name: _sha(ROOT / name) for name in SOURCE_FILES}
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    if not pilot:
        subprocess.run(["git", "diff", "--quiet", "HEAD"], cwd=ROOT, check=True)
        for name in (
            "experiments/interacting_dmft.py",
            "kernelgecp/dmft.py",
            "kernelgecp/lehmann.py",
            "kernelgecp/integrations/triqs_dmft.py",
        ):
            subprocess.run(
                ["git", "ls-files", "--error-unmatch", name],
                cwd=ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
            )
    output.mkdir(parents=True, exist_ok=False)
    report: dict[str, Any] = {
        "schema": "kernelgecp.interacting-dmft.v1",
        "manifest_version": 1,
        "units": {
            "energy": "hopping_t",
            "beta": "inverse_hopping_t",
            "tau": "inverse_hopping_t",
            "g_tau": "dimensionless",
            "g_iw": "inverse_hopping_t",
            "sigma_iw": "hopping_t",
            "sigma_tau": "hopping_t_squared",
        },
        "state_ordering": "scalar_paramagnetic_per_spin_no_block_reordering",
        "chemical_potential_convention": "mu_U_over_2_cancelled_with_static_Hartree",
        "source_revision": revision,
        "config_sha256": _sha(config_path),
        "pilot": pilot,
        "config": data,
        "references": {},
        "calibration": [],
        "selection": [],
        "validation": [],
        "cutoff_refinement": [],
        "dense_validation": [],
    }
    selections = []
    for interaction in data["interactions"]:
        case = "U" + str(interaction).replace(".", "_")
        reference_dir = output / case / "reference"
        reference = _child(
            {"operation": "reference", "config": data, "interaction": interaction},
            reference_dir,
        )
        report["references"][case] = reference
        if reference.get("status") in ("failed", "timeout"):
            continue
        for method in METHODS:
            chosen = None
            for cutoff in data["green_cutoffs"]:
                for eps in data["native_eps_candidates"]:
                    request = {
                        "operation": "case",
                        "config": data,
                        "interaction": interaction,
                        "method": method,
                        "cutoff": cutoff,
                        "eps": eps,
                        "reference": str(reference_dir / "calibration.npz"),
                    }
                    directory = output / case / method / f"calibration-{cutoff}-{eps}"
                    result = _child(request, directory)
                    report["calibration"].append(result)
                    _write(output / "report.json", report)
                    if (
                        result.get("converged") is True
                        and max(result.get("errors", {"missing": math.inf}).values())
                        <= data["calibration_target"]
                        and max(
                            result.get(
                                "finite_fixture_errors", {"missing": math.inf}
                            ).values()
                        )
                        <= data["target"]
                    ):
                        chosen = request
                        break
                if chosen is not None:
                    break
            if chosen is not None:
                selections.append(
                    {
                        **chosen,
                        "reference": str(reference_dir / "validation.npz"),
                        "replay": str(reference_dir / "trajectory.npz"),
                    }
                )
                report["selection"].append(
                    {
                        "case": case,
                        "method": method,
                        "cutoff": chosen["cutoff"],
                        "eps": chosen["eps"],
                    }
                )
                refined = {**chosen, "cutoff": 2 * chosen["cutoff"]}
                refinement = _child(
                    refined, output / case / method / "cutoff-refinement"
                )
                report["cutoff_refinement"].append(refinement)
    # Freeze all choices before opening any final reference in a selected worker.
    _write(output / "selection.json", report["selection"])
    for repeat in range(1 if pilot else data["repeats"]):
        ordered = (
            selections[repeat % max(1, len(selections)) :]
            + selections[: repeat % max(1, len(selections))]
        )
        for request in ordered:
            result = _child(
                {
                    **request,
                    "reference": request["reference"].replace(
                        "validation.npz", "calibration.npz"
                    ),
                    "replay": None,
                }
                if pilot
                else request,
                output
                / "validation"
                / f"repeat-{repeat}-{request['method']}-U{request['interaction']}",
            )
            result["repeat"] = repeat
            result["accuracy_passed"] = (
                result.get("converged") is True
                and max(result.get("errors", {"missing": math.inf}).values())
                <= data["target"]
                and max(
                    result.get("finite_fixture_errors", {"missing": math.inf}).values()
                )
                <= data["target"]
                and (
                    pilot
                    or (
                        result.get("replay") is not None
                        and max(result["replay"]["errors"].values()) <= data["target"]
                    )
                )
                and result.get("physics", {}).get("density_error_from_half", math.inf)
                <= data["target"]
                and result.get("physics", {}).get("positive_imag_green_max", math.inf)
                <= data["target"]
                and result.get("physics", {}).get("positive_imag_sigma_max", math.inf)
                <= data["target"]
            )
            report["validation"].append(result)
            _write(output / "report.json", report)
        for interaction in data["interactions"]:
            case = "U" + str(interaction).replace(".", "_")
            result = _child(
                {
                    "operation": "dense",
                    "config": data,
                    "interaction": interaction,
                    "reference": str(
                        output
                        / case
                        / "reference"
                        / ("calibration.npz" if pilot else "validation.npz")
                    ),
                },
                output / "dense-validation" / f"repeat-{repeat}-U{interaction}",
            )
            result["repeat"] = repeat
            result["accuracy_passed"] = (
                result.get("converged") is True
                and max(result.get("errors", {"missing": math.inf}).values())
                <= data["target"]
            )
            report["dense_validation"].append(result)
    report["all_mandatory_passed"] = (
        len(report["selection"]) == len(METHODS) * len(data["interactions"])
        and len(report["validation"])
        == (1 if pilot else 3) * len(METHODS) * len(data["interactions"])
        and len(report["dense_validation"])
        == (1 if pilot else 3) * len(data["interactions"])
        and all(x["accuracy_passed"] for x in report["validation"])
        and len(report["cutoff_refinement"]) == len(report["selection"])
        and all(
            x.get("converged") is True
            and max(x.get("errors", {"missing": math.inf}).values())
            <= data["calibration_target"]
            for x in report["cutoff_refinement"]
        )
        and all(x["accuracy_passed"] for x in report["dense_validation"])
    )
    report["versions"] = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "mpmath", "pydlr", "sparse-ir", "pylibsparseir")
    }
    report["platform"] = {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "threads": dict.fromkeys(THREADS, "1"),
    }
    report["source_hashes"] = source_hashes
    report["triqs_version"] = importlib.import_module("triqs.version").version
    report["environment_freeze"] = subprocess.check_output(
        [sys.executable, "-m", "pip", "freeze", "--all"], text=True
    ).splitlines()
    if not pilot:
        subprocess.run(["git", "diff", "--quiet", "HEAD"], cwd=ROOT, check=True)
        if (
            revision
            != subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip()
        ):
            raise RuntimeError("source revision changed during canonical generation")
        if source_hashes != {
            name: _sha(ROOT / name) for name in SOURCE_FILES
        } or config_hash != _sha(config_path):
            raise RuntimeError(
                "source/config bytes changed during canonical generation"
            )
    fingerprints: dict[str, list[Any]] = {}
    for item in report["validation"] + report["dense_validation"]:
        fingerprints.setdefault(
            f"{item.get('method')}-U{item.get('interaction')}", []
        ).append(item.get("numerical_fingerprints"))
    report["deterministic_repetitions"] = {
        name: len(items) == 3 and bool(items[0]) and all(x == items[0] for x in items)
        for name, items in fingerprints.items()
    }
    report["determinism_scope"] = (
        "final_field_arrays_and_coefficients_only_timing_excluded"
    )
    report["campaign_checks_passed"] = report["all_mandatory_passed"]
    report["all_mandatory_passed"] &= not pilot and all(
        report["deterministic_repetitions"].values()
    )
    _write(output / "report.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=ROOT / "experiments/interacting_dmft.toml"
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker is not None:
        _worker(json.loads(args.worker.read_text()), args.worker_output)
        return
    if args.output_dir is None:
        parser.error("--output-dir is required and must not already exist")
    result = run(args.config, args.output_dir.resolve(), pilot=args.pilot)
    print(
        json.dumps(
            {
                "all_mandatory_passed": result["all_mandatory_passed"],
                "output": str(args.output_dir),
            }
        )
    )
    if not result["campaign_checks_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
