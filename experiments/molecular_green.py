"""Calibrate and validate matrix compression of four actual native GF2 runs.

Native model data, its independent IR reconstruction and exact-ERI/DF finite
Hamiltonians remain separate from the compressed representations. This runner
does not replace a Green solver transform or claim a native solver speedup.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib
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
THREADS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)
SOURCES = (
    "experiments/molecular_green.py",
    "experiments/molecular_green.toml",
    "experiments/molecular_native.py",
    "experiments/molecular_reference.py",
    "experiments/dlr_adapters.py",
    "kernelgecp/integrations/green.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "environments/green-science/uv.lock",
)


def normal(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return normal(dataclasses.asdict(value))
    if isinstance(value, np.ndarray):
        return (
            {"real": value.real.tolist(), "imag": value.imag.tolist()}
            if np.iscomplexobj(value)
            else value.tolist()
        )
    if isinstance(value, np.generic):
        return normal(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {key: normal(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [normal(item) for item in value]
    return value


def write(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(normal(value), sort_keys=True, indent=2, allow_nan=False) + "\n"
    )


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    data = tomllib.loads(path.read_text())
    if (
        type(data.get("schema_version")) is not int
        or data.get("schema_version") != 1
        or data.get("methods") != list(METHODS)
        or data.get("bonds_angstrom") != [0.74, 1.5]
        or data.get("native_lambdas") != [1e4, 1e5]
        or data.get("beta") != 20
        or data.get("repeats") != 3
        or data.get("validation_tau_rule") != "native_mod4_3_plus_gauss64"
    ):
        raise ValueError("the complete declared EF campaign is required")
    for name in (
        "target",
        "calibration_target",
        "native_refinement_target",
        "native_dyson_target",
        "fit_rcond",
        "selfenergy_cutoff_multiplier",
        "worker_timeout_seconds",
    ):
        value = data.get(name)
        if (
            isinstance(value, bool)
            or not isinstance(value, (float, int))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError("invalid positive finite control: " + name)
    if (
        not data["calibration_target"] < data["target"]
        or data["selfenergy_cutoff_multiplier"] < 3
    ):
        raise ValueError(
            "calibration must be stricter; Sigma support must cover "
            "the three-factor scale"
        )
    for name, maximum in (
        ("target", 1e-6),
        ("calibration_target", 1e-7),
        ("native_refinement_target", 1e-7),
        ("native_dyson_target", 1e-6),
    ):
        if data[name] > maximum:
            raise ValueError("accepted EF accuracy budget cannot be weakened: " + name)
    for name in ("grid_order", "max_rank", "reference_precision_bits"):
        if type(data.get(name)) is not int or data[name] < 1:
            raise ValueError("invalid integer control: " + name)
    if data.get("construction_precision") not in ("float64", "reference"):
        raise ValueError("construction precision must be explicit")
    for name in ("green_cutoffs_hartree", "native_eps_candidates"):
        values = data.get(name)
        if (
            not isinstance(values, list)
            or not values
            or any(
                type(x) not in (float, int) or not math.isfinite(x) or x <= 0
                for x in values
            )
        ):
            raise ValueError("invalid candidate ladder: " + name)
        if len(set(values)) != len(values):
            raise ValueError("candidate ladders must have no duplicates")
    for name in ("calibration_indices", "validation_indices"):
        values = data.get(name)
        if (
            not isinstance(values, list)
            or not values
            or any(type(x) is not int for x in values)
            or len(set(values)) != len(values)
        ):
            raise ValueError("invalid ordinary-index query list: " + name)
    if set(data["calibration_indices"]) & set(data["validation_indices"]):
        raise ValueError("final and calibration frequency queries must be disjoint")
    if any(eps >= 1 for eps in data["native_eps_candidates"]):
        raise ValueError("native representation tolerances must be below one")
    return data


def split_time_rows(size: int) -> tuple[np.ndarray, np.ndarray]:
    training = np.unique(np.r_[np.arange(0, size, 2), size - 1])
    calibration = np.setdiff1d(np.arange(1, size, 4), training)
    return training, calibration


def query_tau(data: Any, role: str) -> np.ndarray:
    if role == "calibration":
        return data.tau[split_time_rows(len(data.tau))[1]]
    if role != "validation":
        raise ValueError("role must be calibration or validation")
    training, calibration = split_time_rows(len(data.tau))
    final_rows = np.setdiff1d(np.arange(3, len(data.tau), 4), training)
    gauss = data.beta * (np.polynomial.legendre.leggauss(64)[0] + 1) / 2
    final = np.unique(np.r_[data.tau[final_rows], gauss])
    if np.any(np.isin(final, data.tau[np.r_[training, calibration]])):
        raise ValueError("final time queries overlap fitting or calibration rows")
    return final


def _basis(method: str, config: dict[str, Any], cutoff: float, eps: float) -> Any:
    if method == "kernelgecp":
        from kernelgecp import FermionicBasis

        return FermionicBasis(
            config["beta"],
            cutoff,
            eps,
            grid_order=config["grid_order"],
            max_rank=config["max_rank"],
            precision=config["construction_precision"],
            precision_bits=config["reference_precision_bits"],
            reference_backend=config["reference_backend"],
        )
    from experiments.dlr_adapters import build_adapter

    return build_adapter(method, config["beta"], cutoff, eps)


def _peak() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _read(folder: Path) -> Any:
    from kernelgecp.integrations.green import read_green_tau

    return read_green_tau(
        folder / "input.h5", folder / "sim.h5", folder / "bound-manifest.json"
    )


def _evaluate(
    basis: Any, coefficients: np.ndarray, query: np.ndarray, domain: str
) -> np.ndarray:
    design = (
        basis.design_tau(query) if domain == "tau" else basis.design_matsubara(query)
    )
    return np.tensordot(design, coefficients, axes=(1, 0))


def _worker(request: dict[str, Any], output: Path) -> dict[str, Any]:
    from kernelgecp.fitting import LeastSquaresPlan

    config, role = request["config"], request["role"]
    for module in ("kernelgecp", "pydlr", "sparse_ir"):
        importlib.import_module(module)
    read_start = time.perf_counter()
    data = _read(Path(request["native_directory"]))
    if data.manifest_sha256 != request["native_manifest_sha256"]:
        raise ValueError("native input provenance changed before fitting")
    read_end = time.perf_counter()
    training, _ = split_time_rows(len(data.tau))
    times = query_tau(data, role)
    indices = np.asarray(config[role + "_indices"], dtype=np.int64)
    started = time.perf_counter()
    green = _basis(request["method"], config, request["cutoff"], request["eps"])
    sigma = _basis(
        request["method"],
        config,
        config["selfenergy_cutoff_multiplier"] * request["cutoff"],
        request["eps"],
    )
    constructed = time.perf_counter()
    green_plan = LeastSquaresPlan(
        green.design_tau(data.tau[training]), rcond=config["fit_rcond"]
    )
    sigma_plan = LeastSquaresPlan(
        sigma.design_tau(data.tau[training]), rcond=config["fit_rcond"]
    )
    prepared = time.perf_counter()
    cg = green_plan.fit(data.g_tau[training])
    cs = sigma_plan.fit(data.sigma_tau[training])
    fitted = time.perf_counter()
    predictions = {
        field + "_" + domain: _evaluate(basis, coefficients, query, domain)
        for field, basis, coefficients in (("green", green, cg), ("sigma", sigma, cs))
        for domain, query in (("tau", times), ("iw", indices))
    }
    endpoints = _evaluate(green, cg, np.array([0.0, data.beta]), "tau")
    evaluated = time.perf_counter()
    peak = _peak()
    # No role-specific reference values are read until both fits/predictions
    # are complete. They cannot influence rank/coefficients inside a worker.
    with np.load(request["reference"], allow_pickle=False) as reference:
        if not np.array_equal(times, reference["tau"]) or not np.array_equal(
            indices, reference["indices"]
        ):
            raise ValueError("native reference query identity changed")
        errors = {}
        for name, values in predictions.items():
            truth = reference[name]
            if truth.shape != values.shape or not np.all(np.isfinite(truth)):
                raise ValueError("invalid native reference shape/values: " + name)
            errors[name] = float(np.max(abs(values - truth)))
    native_density = complex(-2 * np.trace(data.g_tau[-1, 0, 0]))
    density = complex(-2 * np.trace(endpoints[-1, 0, 0]))
    positive = indices >= 0
    causality = {}
    for field in ("green", "sigma"):
        matrix = predictions[field + "_iw"][positive]
        imaginary_part = (matrix - matrix.swapaxes(-1, -2).conj()) / (2j)
        causality[field] = float(np.max(np.linalg.eigvalsh(imaginary_part)))
    arrays = {
        "training_rows": training,
        "tau": times,
        "indices": indices,
        "g_coefficients": cg,
        "sigma_coefficients": cs,
        "green_endpoints": endpoints,
        **predictions,
    }
    np.savez_compressed(output.with_suffix(".npz"), **arrays)
    result = {
        "status": "completed",
        "method": request["method"],
        "case": request["case"],
        "role": role,
        "green_cutoff_hartree": request["cutoff"],
        "sigma_cutoff_hartree": config["selfenergy_cutoff_multiplier"]
        * request["cutoff"],
        "native_eps": request["eps"],
        "green_rank": green.rank,
        "sigma_rank": sigma.rank,
        "errors": errors,
        "density": density,
        "native_density": native_density,
        "density_error": abs(density - native_density),
        "causality_imag_max": causality,
        "plans": {"green": green_plan.diagnostics, "sigma": sigma_plan.diagnostics},
        "basis_metadata": {
            name: basis.diagnostics
            if request["method"] == "kernelgecp"
            else basis.metadata
            for name, basis in (("green", green), ("sigma", sigma))
        },
        "training_residuals": {
            "green": float(np.max(abs(green_plan.residual(data.g_tau[training], cg)))),
            "sigma": float(
                np.max(abs(sigma_plan.residual(data.sigma_tau[training], cs)))
            ),
        },
        "native_iteration": data.iteration,
        "native_manifest_sha256": data.manifest_sha256,
        "numerical_arrays_sha256": {
            key: hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
            for key, value in arrays.items()
        },
        "measurements": {
            "read_seconds": read_end - read_start,
            "construction_seconds": constructed - started,
            "preparation_seconds": prepared - constructed,
            "fit_seconds": fitted - prepared,
            "evaluation_seconds": evaluated - fitted,
            "application_seconds": evaluated - started,
            "pre_validation_peak_rss_bytes": peak,
            "timing_scope": (
                "matched_warm_provider_imports_matrix_construction_SVD_fit_evaluation"
            ),
            "coefficient_bytes": cg.nbytes + cs.nbytes,
        },
        "scope": (
            "same_native_GF2_solution_representation_not_native_solver_acceleration"
        ),
        "certification_status": "not_certified_sampled_representation",
    }
    write(output, result)
    return result


def _child(request: dict[str, Any], folder: Path) -> dict[str, Any]:
    folder.mkdir(parents=True)
    write(folder / "request.json", request)
    started = time.perf_counter()
    try:
        with (folder / "worker.log").open("w") as log:
            process = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "experiments.molecular_green",
                    "--worker",
                    str(folder / "request.json"),
                    "--worker-output",
                    str(folder / "result.json"),
                ],
                cwd=ROOT,
                env={**os.environ, **dict.fromkeys(THREADS, "1")},
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=request["config"]["worker_timeout_seconds"],
            )
        if process.returncode:
            result = {"status": "failed", "returncode": process.returncode}
        else:
            result = json.loads((folder / "result.json").read_text())
    except subprocess.TimeoutExpired:
        result = {"status": "timeout"}
    result.update(case=request["case"], method=request["method"])
    result.setdefault("measurements", {})["process_seconds"] = (
        time.perf_counter() - started
    )
    result["artifacts"] = {
        path.name: sha(path) for path in folder.iterdir() if path.is_file()
    }
    return result


def _passed(result: dict[str, Any], budget: float) -> bool:
    errors = result.get("errors")
    return (
        result.get("status") == "completed"
        and isinstance(errors, dict)
        and set(errors) == {"green_tau", "green_iw", "sigma_tau", "sigma_iw"}
        and all(
            type(value) in (float, int)
            and math.isfinite(value)
            and 0 <= value <= budget
            for value in errors.values()
        )
    )


def _df_reference(
    folder: Path, data: Any, reference: Any
) -> tuple[Any, dict[str, Any]]:
    import h5py

    from experiments.molecular_reference import reference_from_integrals

    with h5py.File(folder / "df_hf_int/meta.h5", "r") as source:
        if (
            not np.array_equal(source["chunk_indices"][:], [0])
            or source["chunk_size"][()] != 1
        ):
            raise ValueError("unexpected bounded H2 density-fitting chunk layout")
    with h5py.File(folder / "df_hf_int/VQ_0.h5", "r") as source:
        raw = np.asarray(source["0"], dtype=np.float64)
    if raw.ndim != 4 or raw.shape[0] != 1 or raw.shape[2:] != (2, 4):
        raise ValueError("unexpected native molecular density-fitting vertex layout")
    vertices = np.ascontiguousarray(raw).view(np.complex128).reshape(raw.shape[1], 2, 2)
    eri = np.einsum("Qpq,Qrs->pqrs", vertices, vertices.conj())
    if np.max(abs(eri.imag)) > 1e-12 or np.max(abs(data.h_core.imag)) > 1e-12:
        raise ValueError(
            "the supported real H2 Hamiltonian has complex native integrals"
        )
    # Independent restricted-HF identities check the native vertex layout and
    # chemists' ERI index convention before using that Hamiltonian in FCI.
    h_core = data.h_core[0, 0].real
    _, orbitals = np.linalg.eigh(data.hf_fock[0, 0].real)
    density = 2 * np.outer(orbitals[:, 0], orbitals[:, 0])
    fock = h_core + np.einsum("rs,pqrs->pq", density, eri.real)
    fock -= 0.5 * np.einsum("rs,prqs->pq", density, eri.real)
    fock_error = float(np.max(abs(fock - data.hf_fock[0, 0])))
    energy = data.energies.nuclear_hartree + 0.5 * np.einsum(
        "pq,qp", density, h_core + fock
    )
    energy_error = float(abs(energy - data.energies.initial_hf_hartree))
    if max(fock_error, energy_error) > 3e-9:
        raise ValueError(
            "native DF vertices do not reproduce the independent HF identities"
        )
    model = reference_from_integrals(
        data.h_core[0, 0].real,
        eri.real,
        beta=data.beta,
        mu=data.mu,
        nuclear_energy=data.energies.nuclear_hartree,
    )
    return model, {
        "core_integral_error_vs_exact_hartree": float(
            np.max(abs(data.h_core[0, 0] - reference.h_core_lowdin))
        ),
        "eri_error_vs_exact_hartree": float(np.max(abs(eri - reference.eri_lowdin))),
        "initial_HF_energy_difference_hartree": data.energies.initial_hf_hartree
        - reference.exact_hf_energy_hartree,
        "native_DF_HF_fock_reconstruction_error_hartree": fock_error,
        "native_DF_HF_energy_reconstruction_error_hartree": energy_error,
        "scope": (
            "density_fitting_Hamiltonian_difference_separate_from_GF2_approximation"
        ),
    }


def run(
    native_root: Path, output: Path, config_path: Path, *, pilot: bool
) -> dict[str, Any]:
    from experiments.molecular_native import NativeIRSolution

    config = load_config(config_path)
    hashes = {name: sha(ROOT / name) for name in SOURCES}
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    if not pilot:
        if (
            config_path.resolve()
            != (ROOT / "experiments/molecular_green.toml").resolve()
        ):
            raise ValueError("canonical runs require the committed configuration path")
        subprocess.run(["git", "diff", "--quiet", "HEAD"], cwd=ROOT, check=True)
        for name in SOURCES:
            subprocess.run(
                ["git", "ls-files", "--error-unmatch", name],
                cwd=ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
            )
    native_report = json.loads((native_root / "native-campaign.json").read_text())
    if (
        native_report.get("schema") != "kernelgecp.green-molecular-native.v1"
        or native_report.get("status") != "passed"
        or len(native_report.get("cases", [])) != 4
    ):
        raise ValueError("all four actual converged native runs are required")
    output.mkdir(parents=True, exist_ok=False)
    report: dict[str, Any] = {
        "schema": "kernelgecp.molecular-green.v1",
        "pilot": pilot,
        "source_revision": revision,
        "source_hashes": hashes,
        "config": config,
        "config_sha256": sha(config_path),
        "native_source_revision": native_report["source_revision"],
        "native_campaign_sha256": sha(native_root / "native-campaign.json"),
        "native": [],
        "calibration": [],
        "selection": [],
        "cutoff_refinement": [],
        "validation": [],
        "model_differences": [],
        "reference_artifacts": {},
    }
    solutions, requests = {}, []
    for row in native_report["cases"]:
        bond, native_lambda = row["bond_length_angstrom"], row["native_grid_lambda"]
        if (
            bond not in config["bonds_angstrom"]
            or native_lambda not in config["native_lambdas"]
        ):
            raise ValueError("unexpected native model/grid case")
        relative = Path(row["directory"])
        if relative.name != str(relative):
            raise ValueError("native directories must be single relative names")
        folder = native_root / relative
        if sha(folder / "bound-manifest.json") != row["manifest_sha256"]:
            raise ValueError(
                "native manifest identity differs from the actual-run receipt"
            )
        data = _read(folder)
        if not data.converged or data.beta != config["beta"]:
            raise ValueError("native run is not the specified converged physical model")
        if (
            data.convergence.energy_threshold_hartree > 1e-10
            or data.convergence.green_threshold > 1e-6
        ):
            raise ValueError("native convergence thresholds exceed the EF contract")
        if (
            data.molecule != "H2"
            or data.orbital_basis.lower() != "sto-3g"
            or data.electrons != 2
            or data.atoms_angstrom != (("H", 0, 0, 0), ("H", 0, 0, bond))
        ):
            raise ValueError(
                "native geometry/basis differs from its declared campaign case"
            )
        native = NativeIRSolution(data, folder / "native-grid.h5", native_lambda)
        key = relative.name
        if (bond, native_lambda) in solutions:
            raise ValueError("duplicate native case")
        solutions[bond, native_lambda] = native
        reference_folder = output / key
        reference_folder.mkdir()
        times = query_tau(data, "calibration")
        indices = np.asarray(config["calibration_indices"], dtype=np.int64)
        values = {
            field + "_" + domain: getattr(
                native, "evaluate_" + ("matsubara" if domain == "iw" else "tau")
            )(query, field=field)
            for field in ("green", "sigma")
            for domain, query in (("tau", times), ("iw", indices))
        }
        np.savez_compressed(
            reference_folder / "calibration.npz", tau=times, indices=indices, **values
        )
        report["reference_artifacts"][key] = {
            "calibration.npz": sha(reference_folder / "calibration.npz")
        }
        native_dyson = native.dyson_residual(indices)
        independent_errors = native.independent_field_reconstruction_errors
        roundtrip_error = native.physical_tau_roundtrip_matrix_error
        report["native"].append(
            {
                "case": key,
                "bond_length_angstrom": bond,
                "native_lambda": native_lambda,
                "manifest_sha256": data.manifest_sha256,
                "table_errors": native.table_errors,
                "endpoint_reconstruction_errors": native.endpoint_reconstruction_errors,
                "independent_field_reconstruction_errors": independent_errors,
                "physical_tau_roundtrip_matrix_error": roundtrip_error,
                "native_dyson_residual": native_dyson,
                "iteration": data.iteration,
                "condition_number": native.condition_number,
            }
        )
        for method in METHODS:
            chosen = None
            for cutoff in config["green_cutoffs_hartree"]:
                for eps in config["native_eps_candidates"]:
                    request = {
                        "case": key,
                        "method": method,
                        "cutoff": cutoff,
                        "eps": eps,
                        "config": config,
                        "native_directory": str(folder.resolve()),
                        "native_manifest_sha256": data.manifest_sha256,
                        "reference": str(
                            (reference_folder / "calibration.npz").resolve()
                        ),
                        "role": "calibration",
                    }
                    result = _child(
                        request,
                        reference_folder / method / f"calibration-{cutoff}-{eps}",
                    )
                    report["calibration"].append(result)
                    write(output / "report.json", report)
                    if _passed(result, config["calibration_target"]):
                        chosen = request
                        break
                if chosen is not None:
                    break
            if chosen is None:
                continue
            requests.append(chosen)
            report["selection"].append(
                {key: chosen[key] for key in ("case", "method", "cutoff", "eps")}
            )
            report["cutoff_refinement"].append(
                _child(
                    {**chosen, "cutoff": 2 * chosen["cutoff"]},
                    reference_folder / method / "cutoff-refinement",
                )
            )
    if len(solutions) != 4:
        raise ValueError("the complete two-bond/two-grid cohort is missing")
    # No final queries have been evaluated before this immutable selection file.
    write(output / "selection.json", report["selection"])
    report["selection_sha256"] = sha(output / "selection.json")
    report["selection_frozen_before_final_queries"] = True
    report["grid_refinement"] = []
    for bond in config["bonds_angstrom"]:
        coarse, refined = solutions[bond, 1e4], solutions[bond, 1e5]
        times = query_tau(refined.data, "calibration" if pilot else "validation")
        indices = np.asarray(
            config["calibration_indices" if pilot else "validation_indices"],
            dtype=np.int64,
        )
        errors = {
            field + "_" + domain: float(
                np.max(
                    abs(
                        getattr(
                            coarse,
                            "evaluate_" + ("matsubara" if domain == "iw" else "tau"),
                        )(query, field=field)
                        - getattr(
                            refined,
                            "evaluate_" + ("matsubara" if domain == "iw" else "tau"),
                        )(query, field=field)
                    )
                )
            )
            for field in ("green", "sigma")
            for domain, query in (("tau", times), ("iw", indices))
        }
        report["grid_refinement"].append(
            {
                "bond_length_angstrom": bond,
                "errors": errors,
                "energy_difference_hartree": abs(
                    coarse.data.energies.total_hartree
                    - refined.data.energies.total_hartree
                ),
                "mu_difference_hartree": abs(coarse.data.mu - refined.data.mu),
            }
        )
    if not pilot:
        from experiments.molecular_reference import build_h2_reference

        for row in native_report["cases"]:
            bond, native_lambda, key = (
                row["bond_length_angstrom"],
                row["native_grid_lambda"],
                row["directory"],
            )
            native = solutions[bond, native_lambda]
            times, indices = (
                query_tau(native.data, "validation"),
                np.asarray(config["validation_indices"], dtype=np.int64),
            )
            values = {
                field + "_" + domain: getattr(
                    native, "evaluate_" + ("matsubara" if domain == "iw" else "tau")
                )(query, field=field)
                for field in ("green", "sigma")
                for domain, query in (("tau", times), ("iw", indices))
            }
            np.savez_compressed(
                output / key / "validation.npz", tau=times, indices=indices, **values
            )
            report["reference_artifacts"][key]["validation.npz"] = sha(
                output / key / "validation.npz"
            )
            for native_record in report["native"]:
                if native_record["case"] == key:
                    native_record["final_native_dyson_residual"] = (
                        native.dyson_residual(indices)
                    )
            reference = build_h2_reference(bond, beta=config["beta"], mu=native.data.mu)
            df_reference, df_differences = _df_reference(
                native_root / key, native.data, reference
            )
            model = {
                "case": key,
                "mu_hartree": native.data.mu,
                "density_fitting": df_differences,
                "exact_FCI_energy_hartree": reference.thermal_energy_hartree,
                "DF_FCI_energy_hartree": df_reference.thermal_energy_hartree,
                "GF2_energy_hartree": native.data.energies.total_hartree,
                "GF2_minus_DF_FCI_energy_hartree": native.data.energies.total_hartree
                - df_reference.thermal_energy_hartree,
                "DF_minus_exact_FCI_energy_hartree": df_reference.thermal_energy_hartree
                - reference.thermal_energy_hartree,
                "exact_FCI_density": reference.electron_count,
                "DF_FCI_density": df_reference.electron_count,
                "sector_energy_error_max": max(
                    x[2] for x in reference.sector_energy_errors
                ),
            }
            for label, other in (("exact_FCI", reference), ("DF_FCI", df_reference)):
                model[label + "_green_errors"] = {
                    domain: float(
                        np.max(
                            abs(
                                values["green_" + domain]
                                - getattr(
                                    other,
                                    "evaluate_"
                                    + ("matsubara" if domain == "iw" else "tau"),
                                )(query)
                            )
                        )
                    )
                    for domain, query in (("tau", times), ("iw", indices))
                }
            report["model_differences"].append(model)
        for repeat in range(3):
            ordered = requests[repeat:] + requests[:repeat]
            for request in ordered:
                key = request["case"]
                final_request = {
                    **request,
                    "role": "validation",
                    "reference": str((output / key / "validation.npz").resolve()),
                }
                result = _child(
                    final_request,
                    output
                    / "validation"
                    / f"repeat-{repeat}-{request['method']}-{key}",
                )
                result["repeat"] = repeat
                report["validation"].append(result)
                write(output / "report.json", report)
    grouped: dict[str, list[Any]] = {}
    for row in report["validation"]:
        grouped.setdefault(row["case"] + ":" + row["method"], []).append(
            row.get("numerical_arrays_sha256")
        )
    report["deterministic_repetitions"] = {
        key: len(values) == 3
        and bool(values[0])
        and all(x == values[0] for x in values)
        for key, values in grouped.items()
    }
    report["all_mandatory_passed"] = (
        not pilot
        and len(report["selection"]) == 12
        and len(report["validation"]) == 36
        and len(report["cutoff_refinement"]) == 12
        and len(report["model_differences"]) == 4
        and all(
            _passed(row, config["target"])
            and row["density_error"] <= config["target"]
            and max(row["causality_imag_max"].values()) <= config["target"]
            for row in report["validation"]
        )
        and all(
            _passed(row, config["calibration_target"])
            for row in report["cutoff_refinement"]
        )
        and all(
            row["native_dyson_residual"] <= config["native_dyson_target"]
            and row["final_native_dyson_residual"] <= config["native_dyson_target"]
            for row in report["native"]
        )
        and all(
            max(row["errors"].values()) <= config["native_refinement_target"]
            and row["energy_difference_hartree"] <= config["native_refinement_target"]
            for row in report["grid_refinement"]
        )
        and len(grouped) == 12
        and all(report["deterministic_repetitions"].values())
    )
    report["calibration_gates_passed"] = (
        len(report["selection"]) == len(report["cutoff_refinement"]) == 12
        and all(
            _passed(row, config["calibration_target"])
            for row in report["cutoff_refinement"]
        )
        and all(
            row["native_dyson_residual"] <= config["native_dyson_target"]
            for row in report["native"]
        )
        and all(
            max(row["errors"].values()) <= config["native_refinement_target"]
            for row in report["grid_refinement"]
        )
    )
    report["versions"] = {
        name: importlib.metadata.version(name)
        for name in (
            "numpy",
            "scipy",
            "mpmath",
            "h5py",
            "pyscf",
            "pydlr",
            "sparse-ir",
            "pylibsparseir",
        )
    }
    report["platform"] = {
        "python": sys.version,
        "platform": platform.platform(),
        "threads": dict.fromkeys(THREADS, "1"),
    }
    report["environment_freeze"] = sorted(
        f"{distribution.metadata['Name']}=={distribution.version}"
        for distribution in importlib.metadata.distributions()
    )
    report["environment_freeze_scope"] = (
        "actual_python_distribution_inventory_with_committed_uv_lock"
    )
    report["source_unchanged_during_generation"] = (
        hashes == {name: sha(ROOT / name) for name in SOURCES}
        and revision
        == subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    )
    if not pilot:
        subprocess.run(["git", "diff", "--quiet", "HEAD"], cwd=ROOT, check=True)
        if not report["source_unchanged_during_generation"]:
            raise RuntimeError("source changed during canonical generation")
        if report["native_campaign_sha256"] != sha(
            native_root / "native-campaign.json"
        ):
            raise RuntimeError("native campaign provenance changed during generation")
    write(output / "report.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--config", type=Path, default=ROOT / "experiments/molecular_green.toml"
    )
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker is not None:
        _worker(json.loads(args.worker.read_text()), args.worker_output)
        return
    if args.native_root is None or args.output_dir is None:
        parser.error("--native-root and --output-dir are required")
    report = run(
        args.native_root.resolve(),
        args.output_dir.resolve(),
        args.config,
        pilot=args.pilot,
    )
    print(
        json.dumps(
            {
                "all_mandatory_passed": report["all_mandatory_passed"],
                "calibration_gates_passed": report["calibration_gates_passed"],
                "pilot": args.pilot,
            }
        )
    )
    if not report["calibration_gates_passed"] or (
        not args.pilot and not report["all_mandatory_passed"]
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
