"""Resumable actual CT-HYB calibration, frozen feedback and replica campaign.

Run under the locked optional TRIQS/DLR/IR environment, with all thread
variables set to one. Fast smoke tests do not satisfy canonical uncertainty
acceptance. Every native measurement is retained for independent EE/EH work.
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
import tempfile
import time
import tomllib
import zipfile
from pathlib import Path
from typing import Any

import numpy as np

from kernelgecp.integrations.cthyb_dmft import (
    CTHYBConfig,
    CTHYBMeasurement,
    CTHYBPlan,
    estimate_replicas,
)

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("kernelgecp", "pydlr", "sparse_ir")
ARRAYS = (
    "bath_iw",
    "g_tau",
    "g_iw",
    "g_iw_raw",
    "sigma_iw",
    "g_moments",
    "sigma_moments",
    "occupancies",
)
SOURCES = (
    "kernelgecp/integrations/cthyb_dmft.py",
    "kernelgecp/integrations/triqs.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "kernelgecp/dmft.py",
    "experiments/cthyb_dmft.py",
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


def _plain(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _plain(dataclasses.asdict(value))
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    return value


def _bytes(value: Any) -> bytes:
    return (
        json.dumps(_plain(value), sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _verify_sources(identity: dict[str, Any]) -> None:
    for relative, digest in identity.get("sources", {}).items():
        if _sha(ROOT / relative) != digest:
            raise ValueError("source changed during campaign execution")


def _write(path: Path, value: Any) -> None:
    data = _bytes(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
    os.replace(temporary, path)


def _read(path: Path) -> dict[str, Any]:
    if path.stat().st_size > 16 * 1024**2:
        raise ValueError("JSON receipt exceeds its bounded allocation")

    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def reject(value: str) -> Any:
        raise ValueError(f"nonfinite JSON number: {value}")

    value = json.loads(
        path.read_text(), object_pairs_hook=unique, parse_constant=reject
    )
    # Encoding also detects an otherwise finite-looking exponent such as 1e999.
    _bytes(value)
    if not isinstance(value, dict):
        raise ValueError("receipt must be a JSON object")
    return value


def load_config(path: Path) -> dict[str, Any]:
    data = tomllib.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("CT-HYB campaign requires schema_version=1")
    CTHYBConfig(
        **{field.name: data[field.name] for field in dataclasses.fields(CTHYBConfig)}
    )
    for name in (
        "calibration_replicas",
        "final_replicas",
        "refinement_replicas",
        "calibration_loop_cycles",
        "calibration_loop_iterations",
        "max_iterations",
        "min_iterations",
        "consecutive_passes",
        "final_bath_average",
        "grid_order",
        "max_rank",
        "precision_bits",
        "fit_holdout_modulus",
        "seed_base",
        "memory_limit_bytes",
    ):
        if type(data.get(name)) is not int or data[name] < 1:
            raise ValueError(f"{name} must be a positive integer")
    for name in (
        "loop_tolerance",
        "max_standard_error",
        "agreement_absolute",
        "agreement_sigma",
        "self_consistency_absolute",
        "replay_selected_target",
        "native_eps",
    ):
        value = data.get(name)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ValueError(f"{name} must be positive and finite")
    if (
        min(
            data["calibration_replicas"],
            data["final_replicas"],
            data["refinement_replicas"],
        )
        < 2
    ):
        raise ValueError("independent uncertainty requires at least two replicas")
    if (
        max(
            data["calibration_replicas"],
            data["final_replicas"],
            data["refinement_replicas"],
        )
        > 64
    ):
        raise ValueError("replicas exceed the declared bounded seed allocation")
    if (
        not data["final_bath_average"]
        <= data["min_iterations"]
        <= data["max_iterations"]
        <= 64
    ):
        raise ValueError("invalid bounded loop/averaging allocation")
    if data["calibration_loop_iterations"] > 64 or data["seed_base"] + 51000 >= 2**31:
        raise ValueError("seed allocation exceeds the signed native seed range")
    if (
        type(data.get("fit_holdout_remainder")) is not int
        or not 0 <= data["fit_holdout_remainder"] < data["fit_holdout_modulus"]
    ):
        raise ValueError("invalid excluded measured-time row rule")
    for role in ("calibration_indices", "selected_indices"):
        indices = data.get(role)
        if (
            not isinstance(indices, list)
            or not indices
            or any(type(x) is not int or not 0 <= x < data["n_iw"] for x in indices)
            or len(set(indices)) != len(indices)
        ):
            raise ValueError(
                "Green query indices must be distinct nonnegative mesh indices"
            )
    if set(data["calibration_indices"]) & set(data["selected_indices"]):
        raise ValueError("calibration and final Green queries must be disjoint")
    budgets = data.get("cycle_budgets")
    if (
        not isinstance(budgets, list)
        or not budgets
        or any(type(x) is not int or x < 1 for x in budgets)
        or budgets != sorted(set(budgets))
    ):
        raise ValueError("cycle ladder must be sorted distinct positive integers")
    cutoffs = data.get("cutoffs")
    if (
        not isinstance(cutoffs, list)
        or not cutoffs
        or any(
            isinstance(x, bool)
            or not isinstance(x, (int, float))
            or not math.isfinite(x)
            or x <= 0
            for x in cutoffs
        )
    ):
        raise ValueError("cutoffs must be finite positive values")
    if (
        data["native_eps"] >= 1
        or not 2 <= data["grid_order"] <= 256
        or not 1 <= data["max_rank"] <= 256
        or not 53 <= data["precision_bits"] <= 1024
        or not 64 * 1024**2 <= data["memory_limit_bytes"] <= 16 * 1024**3
        or any(
            not math.isfinite(data["beta"] * x) or data["beta"] * x < 1 for x in cutoffs
        )
    ):
        raise ValueError(
            "basis/CPU-memory parameters exceed the bounded runtime contract"
        )
    if data.get("canonical") is not True and data.get("canonical") is not False:
        raise ValueError("canonical must be an explicit boolean")
    if data["canonical"] and (
        budgets != [100000, 400000, 1600000]
        or data["beta"] != 10
        or data["interaction"] != 2
        or data["hopping"] != 1
        or data["max_standard_error"] > 3e-4
        or data["agreement_absolute"] > 1e-3
        or data["agreement_sigma"] != 3
    ):
        raise ValueError(
            "canonical EE physical, budget or uncertainty gates were weakened"
        )
    return data


def _native_config(data: dict[str, Any], *, refined: bool = False) -> CTHYBConfig:
    values = {field.name: data[field.name] for field in dataclasses.fields(CTHYBConfig)}
    if refined:
        values["n_iw"] *= 2
        values["n_tau"] = 2 * values["n_tau"] - 1
    return CTHYBConfig(**values)


def _identity(data: dict[str, Any]) -> dict[str, Any]:
    for variable in THREADS:
        if os.environ.get(variable) != "1":
            raise ValueError(f"{variable}=1 is required")
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    source_hashes = {name: _sha(ROOT / name) for name in SOURCES}
    if data["canonical"]:
        for name, digest in source_hashes.items():
            committed = subprocess.check_output(
                ["git", "show", revision + ":" + name], cwd=ROOT
            )
            if hashlib.sha256(committed).hexdigest() != digest:
                raise ValueError(
                    "canonical generation requires committed frozen source bytes"
                )
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "scipy", "pydlr", "sparse-ir", "pylibsparseir")
    }
    import triqs.version
    import triqs_cthyb.version

    versions.update(triqs=triqs.version.version, cthyb=triqs_cthyb.version.version)
    if versions["triqs"] != "4.0.2" or versions["cthyb"] != "4.0.1":
        raise ValueError("the pinned actual provider versions are required")
    return {
        "schema": "kernelgecp.cthyb-campaign.v1",
        "source_revision": revision,
        "sources": source_hashes,
        "config": data,
        "versions": versions,
        "python": platform.python_version(),
        "platform": platform.system() + "-" + platform.machine(),
        "threads": {name: "1" for name in THREADS},
        "claim_scope": "actual_CT-HYB_finite_grid_stochastic_DMFT",
        "certification_status": "not_certified",
    }


def save_measurement(
    stem: Path, measurement: CTHYBMeasurement, identity: dict[str, Any]
) -> None:
    """Atomic numeric archive plus a hash-bound, independently readable receipt."""
    array_path, receipt_path = stem.with_suffix(".npz"), stem.with_suffix(".json")
    if array_path.exists() or receipt_path.exists():
        raise FileExistsError("refusing to overwrite a retained native measurement")
    stem.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=stem.parent, suffix=".npz", delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        np.savez_compressed(
            temporary, **{name: getattr(measurement, name) for name in ARRAYS}
        )
        os.replace(temporary, array_path)
    finally:
        temporary.unlink(missing_ok=True)
    receipt = {
        "schema": "kernelgecp.cthyb-measurement.v1",
        "identity": identity,
        "arrays_sha256": _sha(array_path),
        "config": dataclasses.asdict(measurement.config),
        **{
            name: getattr(measurement, name)
            for name in (
                "seed",
                "cycles",
                "double_occupancy",
                "average_sign",
                "solve_status",
            )
        },
        "spin_blocks": ["up", "down"],
        "ordinary_indices": "-n_iw..n_iw-1",
        "density_source": "native_reduced_density_matrix",
        "moment_source": "native_density_matrix_commutators",
        "thermalization": "fresh_solver_and_warmup_for_each_archive",
    }
    _write(receipt_path, receipt)


def load_measurement(stem: Path, identity: dict[str, Any]) -> CTHYBMeasurement:
    array_path, receipt_path = stem.with_suffix(".npz"), stem.with_suffix(".json")
    if array_path.is_symlink() or receipt_path.is_symlink():
        raise ValueError("measurement archives cannot be symlinks")
    receipt = _read(receipt_path)
    if (
        receipt.get("schema") != "kernelgecp.cthyb-measurement.v1"
        or receipt.get("identity") != identity
    ):
        raise ValueError("measurement source/config/role identity mismatch")
    if receipt.get("arrays_sha256") != _sha(array_path):
        raise ValueError("measurement array hash mismatch")
    if receipt.get("spin_blocks") != ["up", "down"]:
        raise ValueError("measurement spin ordering mismatch")
    config = CTHYBConfig(**receipt["config"])
    shapes = {
        "bath_iw": (2 * config.n_iw,),
        "g_tau": (config.n_tau, 2),
        "g_iw": (2 * config.n_iw, 2),
        "g_iw_raw": (2 * config.n_iw, 2),
        "sigma_iw": (2 * config.n_iw, 2),
        "g_moments": (3, 2),
        "sigma_moments": (2, 2),
        "occupancies": (2,),
    }
    if array_path.stat().st_size > 64 * 1024**2:
        raise ValueError("measurement ZIP exceeds its bounded allocation")
    with zipfile.ZipFile(array_path) as zipped:
        if (
            len(zipped.infolist()) != len(ARRAYS)
            or sum(item.file_size for item in zipped.infolist()) > 64 * 1024**2
        ):
            raise ValueError("measurement ZIP allocation exceeds its bounded schema")
        if {entry.filename for entry in zipped.infolist()} != {
            name + ".npy" for name in ARRAYS
        }:
            raise ValueError("measurement ZIP member inventory mismatch")
        for name in ARRAYS:
            with zipped.open(name + ".npy") as stream:
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
                elif version == (2, 0):
                    shape, fortran, dtype = np.lib.format.read_array_header_2_0(stream)
                else:
                    raise ValueError("unsupported measurement NPY version")
                expected = np.dtype(
                    np.float64 if name == "occupancies" else np.complex128
                )
                if shape != shapes[name] or fortran or dtype != expected:
                    raise ValueError(
                        "measurement NPY header violates its bounded schema"
                    )
                required = stream.tell() + math.prod(shape) * dtype.itemsize
                if required != zipped.getinfo(name + ".npy").file_size:
                    raise ValueError("measurement NPY payload size mismatch")
    with np.load(array_path, allow_pickle=False) as archive:
        if set(archive.files) != set(ARRAYS):
            raise ValueError("measurement array inventory mismatch")
        arrays = {name: archive[name] for name in ARRAYS}
    return CTHYBMeasurement(
        config=config,
        **arrays,
        **{
            name: receipt[name]
            for name in (
                "seed",
                "cycles",
                "double_occupancy",
                "average_sign",
                "solve_status",
            )
        },
    )


def _sample(
    plan: CTHYBPlan,
    bath: np.ndarray,
    cycles: int,
    seed: int,
    stem: Path,
    identity: dict[str, Any],
) -> CTHYBMeasurement:
    _verify_sources(identity)
    if stem.with_suffix(".json").exists():
        measurement = load_measurement(stem, identity)
        if (
            measurement.config != plan.config
            or measurement.cycles != cycles
            or measurement.seed != seed
            or not np.array_equal(measurement.bath_iw, bath)
        ):
            raise ValueError("resumed chain does not match its requested input")
        return measurement
    if stem.with_suffix(".npz").exists():
        raise ValueError(
            "incomplete measurement checkpoint requires explicit inspection"
        )
    start = time.perf_counter()
    measurement = plan.sample(bath, cycles=cycles, seed=seed)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = int(peak if platform.system() == "Darwin" else 1024 * peak)
    if peak_bytes > identity["config"]["memory_limit_bytes"]:
        raise MemoryError("actual campaign RSS exceeded the frozen CPU memory budget")
    save_measurement(stem, measurement, identity)
    print(
        json.dumps(
            {
                "measurement": str(stem.relative_to(stem.parents[1])),
                "cycles": cycles,
                "seed": seed,
                "seconds": time.perf_counter() - start,
            }
        ),
        flush=True,
    )
    return measurement


def _basis(method: str, data: dict[str, Any], cutoff: float) -> Any:
    if method == "kernelgecp":
        from kernelgecp import FermionicBasis

        return FermionicBasis(
            data["beta"],
            cutoff,
            data["native_eps"],
            grid_order=data["grid_order"],
            max_rank=data["max_rank"],
            precision="reference",
            precision_bits=data["precision_bits"],
        )
    from experiments.dlr_adapters import build_adapter

    return build_adapter(method, data["beta"], cutoff, data["native_eps"])


def _plans(data: dict[str, Any], cutoff: float) -> dict[str, CTHYBPlan]:
    config = _native_config(data)
    observed = np.flatnonzero(
        np.arange(config.n_tau) % data["fit_holdout_modulus"]
        != data["fit_holdout_remainder"]
    )
    return {
        name: CTHYBPlan(config, _basis(name, data, cutoff), fit_tau_indices=observed)
        for name in METHODS
    }


def _replicas(
    plan: CTHYBPlan,
    bath: np.ndarray,
    cycles: int,
    count: int,
    seed: int,
    folder: Path,
    identity: dict[str, Any],
    indices: list[int],
    *,
    represented: bool = False,
) -> tuple[list[CTHYBMeasurement], Any]:
    measurements = [
        _sample(plan, bath, cycles, seed + i, folder / f"replica-{i:02d}", identity)
        for i in range(count)
    ]
    estimate = estimate_replicas(
        np.array(
            [
                plan.observables(
                    m,
                    indices,
                    representation=plan.represent(m) if represented else None,
                )
                for m in measurements
            ]
        )
    )
    return measurements, estimate


def _loop(
    plan: CTHYBPlan,
    data: dict[str, Any],
    cycles: int,
    seed: int,
    folder: Path,
    identity: dict[str, Any],
    *,
    compressed: bool,
    pilot: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    bath = plan.initial_bath()
    history = []
    baths = []
    consecutive = 0
    maximum = data["calibration_loop_iterations"] if pilot else data["max_iterations"]
    for iteration in range(maximum):
        measurement = _sample(
            plan,
            bath,
            cycles,
            seed + iteration,
            folder / f"iteration-{iteration:02d}",
            identity,
        )
        representation = plan.represent(measurement) if compressed else None
        updated, defect = plan.bath_update(measurement, representation)
        baths.append(bath)
        history.append(
            {
                "iteration": iteration + 1,
                "unmixed_defect": defect,
                "mixed_step": float(np.max(np.abs(updated - bath))),
                "density": measurement.density,
                "double_occupancy": measurement.double_occupancy,
                "representation_discrepancy": None
                if representation is None
                else representation.native_iw_discrepancy,
            }
        )
        consecutive = consecutive + 1 if defect < data["loop_tolerance"] else 0
        if (
            not pilot
            and iteration + 1 >= data["min_iterations"]
            and consecutive >= data["consecutive_passes"]
        ):
            break
        bath = updated
    converged = (
        not pilot
        and len(history) >= data["min_iterations"]
        and consecutive >= data["consecutive_passes"]
    )
    # Average actual input baths, never include an unmeasured proposed iterate.
    final = np.mean(baths[-data["final_bath_average"] :], axis=0)
    report = {
        "history": history,
        "converged": converged,
        "stop_reason": "calibration_pilot_budget"
        if pilot
        else "stochastic_controller"
        if converged
        else "max_iterations_inconclusive",
        "final_bath_rule": "mean_of_last_declared_number_of_measured_input_baths",
        "averaged_baths": min(len(baths), data["final_bath_average"]),
    }
    _write(folder / "loop.json", report)
    return final, report


def calibrate(
    data: dict[str, Any], output: Path, identity: dict[str, Any]
) -> dict[str, Any]:
    freeze_path = output / "freeze.json"
    if freeze_path.exists():
        frozen = _read(freeze_path)
        if frozen.get("identity") != identity:
            raise ValueError("frozen calibration differs from current source/config")
        return frozen
    dense = CTHYBPlan(_native_config(data))
    bath, pilot = _loop(
        dense,
        data,
        data["calibration_loop_cycles"],
        data["seed_base"],
        output / "calibration-loop",
        identity,
        compressed=False,
        pilot=True,
    )
    records = []
    candidate_measurements = []
    for level, cycles in enumerate(data["cycle_budgets"]):
        measurements, estimate = _replicas(
            dense,
            bath,
            cycles,
            data["calibration_replicas"],
            data["seed_base"] + 1000 * (level + 1),
            output / f"calibration-{cycles}",
            identity,
            data["calibration_indices"],
        )
        maximum = float(np.max(estimate.standard_error))
        records.append(
            {
                "cycles": cycles,
                "estimate": _plain(estimate),
                "max_standard_error": maximum,
                "uncertainty_pass": maximum <= data["max_standard_error"],
            }
        )
        candidate_measurements.append(measurements)
        _write(
            output / "calibration-progress.json",
            {"identity": identity, "budgets": records},
        )
    passing = [i for i, record in enumerate(records) if record["uncertainty_pass"]]
    if not passing:
        result = {
            "identity": identity,
            "status": "inconclusive_uncertainty_budget_exhausted",
            "budgets": records,
        }
        _write(output / "calibration-result.json", result)
        return result
    selected = passing[0]
    measurements = candidate_measurements[selected]
    replay = []
    selected_cutoff = None
    offset = np.array(data["calibration_indices"]) + data["n_iw"]
    for cutoff in data["cutoffs"]:
        plans = _plans(data, cutoff)
        providers = {}
        for method, plan in plans.items():
            reps = [plan.represent(m) for m in measurements]
            native = np.mean([m.g_iw for m in measurements], axis=0)
            represented = np.mean([r.g_iw for r in reps], axis=0)
            error = float(np.max(np.abs(represented[offset] - native[offset])))
            providers[method] = {
                "rank": reps[0].rank,
                "retained_rank": reps[0].retained_rank,
                "condition_estimate": reps[0].condition_estimate,
                "selected_mean_iw_discrepancy": error,
                "maximum_noisy_time_holdout_residual": max(
                    r.validation_residual for r in reps
                ),
                "pass": error <= data["replay_selected_target"],
            }
        replay.append({"cutoff": cutoff, "providers": providers})
        if selected_cutoff is None and all(p["pass"] for p in providers.values()):
            selected_cutoff = cutoff
    if selected_cutoff is None:
        result = {
            "identity": identity,
            "status": "inconclusive_representation_calibration",
            "budgets": records,
            "replay": replay,
        }
        _write(output / "calibration-result.json", result)
        return result
    frozen = {
        "identity": identity,
        "status": "frozen_calibration_not_final_acceptance",
        "cycles": records[selected]["cycles"],
        "cutoff": selected_cutoff,
        "budgets": records,
        "replay": replay,
        "pilot": pilot,
        "independent_final_seeds": True,
        "calibration_archives": {
            str(p.relative_to(output)): _sha(p)
            for p in sorted(output.glob("calibration-*/*.*"))
            if p.is_file()
        },
    }
    _verify_sources(identity)
    _write(freeze_path, frozen)
    return frozen


def run_final(
    data: dict[str, Any], output: Path, identity: dict[str, Any], frozen: dict[str, Any]
) -> dict[str, Any]:
    if (
        frozen.get("status") != "frozen_calibration_not_final_acceptance"
        or frozen.get("identity") != identity
    ):
        raise ValueError("final generation requires a matching frozen calibration")
    for relative, digest in frozen["calibration_archives"].items():
        if _sha(output / relative) != digest:
            raise ValueError("frozen calibration evidence changed")
    plans = _plans(data, frozen["cutoff"])
    dense = CTHYBPlan(_native_config(data))
    cycles = frozen["cycles"]
    loops = {}
    measurements = {}
    estimates = {}
    native_estimates = {}
    residuals = {}
    final_baths = {}
    for name, plan, seed in (
        ("dense", dense, 10000),
        ("kernelgecp", plans["kernelgecp"], 20000),
    ):
        bath, loop = _loop(
            plan,
            data,
            cycles,
            data["seed_base"] + seed,
            output / (name + "-loop"),
            identity,
            compressed=name != "dense",
        )
        loops[name] = loop
        final_baths[name] = bath
        values, estimate = _replicas(
            plan,
            bath,
            cycles,
            data["final_replicas"],
            data["seed_base"] + seed + 20000,
            output / (name + "-final"),
            identity,
            data["selected_indices"],
            represented=name != "dense",
        )
        measurements[name] = values
        estimates[name] = estimate
        native_estimate = estimate_replicas(
            np.array([dense.observables(m, data["selected_indices"]) for m in values])
        )
        native_estimates[name] = native_estimate
        mean = np.mean([m.g_iw.mean(axis=1) for m in values], axis=0)
        positions = np.array(data["selected_indices"]) + data["n_iw"]
        expected = bath[positions] / data["hopping"] ** 2
        selected = mean[positions]
        count = len(positions)
        errors = native_estimate.standard_error[2:]
        defects = np.r_[(selected - expected).real, (selected - expected).imag]
        limit = data["self_consistency_absolute"] + 3 * errors
        residuals[name] = {
            "selected_real_imag_defects": defects,
            "standard_errors": errors,
            "limits": limit,
            "pass": bool(np.all(np.abs(defects) <= limit)),
            "full_grid_mean_defect": float(
                np.max(np.abs(data["hopping"] ** 2 * mean - bath))
            ),
        }
        algorithm_g = estimate.mean[2 : 2 + count] + 1j * estimate.mean[2 + count :]
        algorithm_defect = np.r_[
            (algorithm_g - expected).real, (algorithm_g - expected).imag
        ]
        algorithm_limit = (
            data["self_consistency_absolute"] + 3 * estimate.standard_error[2:]
        )
        residuals[name]["represented_selected_defects"] = algorithm_defect
        residuals[name]["represented_limits"] = algorithm_limit
        residuals[name]["represented_pass"] = bool(
            np.all(np.abs(algorithm_defect) <= algorithm_limit)
        )
        assert errors.size == 2 * count
    first, second = estimates["dense"], estimates["kernelgecp"]
    combined = np.sqrt(first.standard_error**2 + second.standard_error**2)
    difference = np.abs(first.mean - second.mean)
    limits = data["agreement_absolute"] + data["agreement_sigma"] * combined
    agreement = bool(np.all(difference <= limits))
    uncertainty = bool(
        max(np.max(first.standard_error), np.max(second.standard_error))
        <= data["max_standard_error"]
    )
    replay = {}
    positions = np.array(data["selected_indices"]) + data["n_iw"]
    for trajectory, values in measurements.items():
        replay[trajectory] = {}
        for method, plan in plans.items():
            reps = [plan.represent(m) for m in values]
            native = np.mean([m.g_iw for m in values], axis=0)
            represented = np.mean([r.g_iw for r in reps], axis=0)
            error = float(np.max(np.abs(represented[positions] - native[positions])))
            replay[trajectory][method] = {
                "rank": reps[0].rank,
                "retained_rank": reps[0].retained_rank,
                "condition_estimate": reps[0].condition_estimate,
                "selected_mean_iw_discrepancy": error,
                "pass": error <= data["replay_selected_target"],
                "native_tail_change_max": max(r.native_tail_change for r in reps),
                "noisy_time_holdout_residual_max": max(
                    r.validation_residual for r in reps
                ),
            }
    # Preserve the exact primary bath rows and extend only new frequencies
    # with its measured moment tail. This is a declared finite-grid extension,
    # not an analyticity or compact-support certificate.
    refined = CTHYBPlan(_native_config(data, refined=True))
    old_bath = final_baths["dense"]
    iw = 1j * (2 * dense.indices + 1) * np.pi / data["beta"]
    tail_rows = np.abs(dense.indices + 0.5) >= 0.8 * data["n_iw"]
    second = float(
        np.mean(((old_bath - data["hopping"] ** 2 / iw) * iw**2)[tail_rows]).real
    )
    extended_iw = 1j * (2 * refined.indices + 1) * np.pi / data["beta"]
    extended = data["hopping"] ** 2 / extended_iw + second / extended_iw**2
    extended[data["n_iw"] : 3 * data["n_iw"]] = old_bath
    _, refined_estimate = _replicas(
        refined,
        extended,
        cycles,
        data["refinement_replicas"],
        data["seed_base"] + 50000,
        output / "dense-refinement",
        identity,
        data["selected_indices"],
    )
    refinement_difference = np.abs(refined_estimate.mean - first.mean)
    refinement_se = np.sqrt(
        refined_estimate.standard_error**2 + first.standard_error**2
    )
    refinement_limits = (
        data["agreement_absolute"] + data["agreement_sigma"] * refinement_se
    )
    refinement = {
        "primary_grid": [data["n_iw"], data["n_tau"]],
        "refined_grid": [2 * data["n_iw"], 2 * data["n_tau"] - 1],
        "bath_extension": "exact_primary_rows_plus_first_two_moment_tail",
        "bath_second_moment": second,
        "estimate": refined_estimate,
        "absolute_difference": refinement_difference,
        "combined_standard_error": refinement_se,
        "limits": refinement_limits,
        "pass": bool(np.all(refinement_difference <= refinement_limits)),
    }
    cutoff_replay = {}
    for cutoff in data["cutoffs"]:
        cutoff_plans = _plans(data, cutoff)
        cutoff_replay[str(cutoff)] = {}
        for method, plan in cutoff_plans.items():
            represented = np.mean(
                [plan.represent(m).g_iw for m in measurements["dense"]], axis=0
            )
            native = np.mean([m.g_iw for m in measurements["dense"]], axis=0)
            error = float(np.max(np.abs(represented[positions] - native[positions])))
            cutoff_replay[str(cutoff)][method] = {
                "selected_mean_iw_discrepancy": error,
                "pass": error <= data["replay_selected_target"],
            }
    checks = {
        "both_stochastic_controllers": all(v["converged"] for v in loops.values()),
        "independent_final_uncertainty": uncertainty,
        "physical_observable_agreement": agreement,
        "fresh_fixed_bath_self_consistency": all(
            v["pass"] and v["represented_pass"] for v in residuals.values()
        ),
        "all_three_frozen_replays": all(
            v["pass"] for row in replay.values() for v in row.values()
        ),
        "native_grid_refinement": refinement["pass"],
        "all_provider_cutoff_sensitivity": all(
            v["pass"] for row in cutoff_replay.values() for v in row.values()
        ),
    }
    result = {
        "schema": "kernelgecp.cthyb-result.v1",
        "identity": identity,
        "freeze_sha256": _sha(output / "freeze.json"),
        "canonical": data["canonical"],
        "status": "passed" if all(checks.values()) else "inconclusive",
        "checks": checks,
        "loops": loops,
        "estimates": estimates,
        "native_estimates": native_estimates,
        "difference": difference,
        "combined_standard_error": combined,
        "agreement_limits": limits,
        "fixed_bath_residuals": residuals,
        "frozen_replay": replay,
        "native_grid_refinement": refinement,
        "cutoff_replay": cutoff_replay,
        "measurement_archives": {
            str(p.relative_to(output)): _sha(p)
            for p in sorted(output.glob("*/*.json"))
            if p.name.startswith(("replica-", "iteration-"))
        },
        "observables": ["density_total", "double_occupancy"]
        + [f"G_real_n{x}" for x in data["selected_indices"]]
        + [f"G_imag_n{x}" for x in data["selected_indices"]],
        "seed_policy": "disjoint_calibration_loop_final_and_replica_families",
        "uncertainty_scope": "independent_fixed_bath_chains_not_DMFT_trajectory_error",
        "native_grid_refinement_status": (
            "independent_fixed_bath_replicas_with_declared_tail_extension"
        ),
    }
    _verify_sources(identity)
    _write(output / "result.json", result)
    return result


def _execute(data: dict[str, Any], args: argparse.Namespace) -> int:
    identity = _identity(data)
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = args.output / "manifest.json"
    if manifest.exists():
        previous = _read(manifest)
        # A documentation-only commit may follow generation. Resume retains
        # the original generation revision while rechecking every source byte.
        old = {
            key: value for key, value in previous.items() if key != "source_revision"
        }
        current = {
            key: value for key, value in identity.items() if key != "source_revision"
        }
        if old != current:
            raise ValueError(
                "campaign source/config/runtime changed; use a new output directory"
            )
        identity = previous
    _write(manifest, identity)
    if args.stage == "final" and not (args.output / "freeze.json").exists():
        raise ValueError("final stage cannot calibrate or tune missing frozen evidence")
    frozen = calibrate(data, args.output, identity)
    if frozen.get("status") != "frozen_calibration_not_final_acceptance":
        print(json.dumps({"status": frozen["status"]}), flush=True)
        return 2
    if args.stage != "calibrate":
        result = run_final(data, args.output, identity, frozen)
        print(
            json.dumps({"status": result["status"], "checks": result["checks"]}),
            flush=True,
        )
        return 0 if result["status"] == "passed" else 2
    print(
        json.dumps(
            {
                "status": frozen["status"],
                "cycles": frozen["cycles"],
                "cutoff": frozen["cutoff"],
            }
        ),
        flush=True,
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("calibrate", "final", "all"))
    parser.add_argument(
        "--config", type=Path, default=ROOT / "experiments/cthyb_dmft.toml"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = load_config(args.config)
    args.output.mkdir(parents=True, exist_ok=True)
    # The supported native runtimes are Linux/macOS. An OS advisory lock is
    # released even after process death, so a crashed campaign can resume
    # without double-running or overwriting an active chain.
    import fcntl

    with (args.output / ".campaign.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("another campaign owns this output directory") from error
        return _execute(data, args)


if __name__ == "__main__":
    raise SystemExit(main())
