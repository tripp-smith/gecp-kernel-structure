"""Independently recompute EE saved-chain statistics, equations and replay.

No campaign acceptance boolean, stored estimate, or estimate_replicas call is
used to decide acceptance. Scalar two-pass covariance and pivoted QR fitting
provide separate implementations of the statistical and representation gates.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import lstsq

from experiments.cthyb_dmft import METHODS, ROOT, _read, _sha, load_measurement


def _require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _statistics(rows: np.ndarray) -> dict[str, Any]:
    """Scalar independent-chain sample covariance, without campaign helpers."""
    _require(
        rows.ndim == 2 and rows.shape[0] >= 2 and np.all(np.isfinite(rows)),
        "independent finite replica matrix required",
    )
    count, width = rows.shape
    means = np.array(
        [math.fsum(float(x) for x in rows[:, j]) / count for j in range(width)]
    )
    covariance = np.empty((width, width))
    for i in range(width):
        for j in range(width):
            covariance[i, j] = math.fsum(
                float((row[i] - means[i]) * (row[j] - means[j])) for row in rows
            ) / (count - 1)
    _require(
        np.all(np.isfinite(covariance)) and np.all(np.diag(covariance) >= 0),
        "replica covariance is nonfinite or has a negative variance",
    )
    errors = np.array([math.sqrt(covariance[j, j] / count) for j in range(width)])
    return {
        "mean": means,
        "covariance": covariance,
        "covariance_of_mean": covariance / count,
        "standard_error": errors,
        "replicas": count,
    }


def _observables(
    m: Any, indices: list[int], green: np.ndarray | None = None
) -> np.ndarray:
    values = m.g_iw if green is None else green
    selected = values[np.array(indices) + m.config.n_iw].mean(axis=1)
    return np.r_[
        math.fsum(m.occupancies), m.double_occupancy, selected.real, selected.imag
    ]


class _QRReplay:
    def __init__(self, method: str, data: dict[str, Any], cutoff: float) -> None:
        if method == "kernelgecp":
            from kernelgecp import FermionicBasis

            basis = FermionicBasis(
                data["beta"],
                cutoff,
                data["native_eps"],
                grid_order=data["grid_order"],
                max_rank=data["max_rank"],
                precision="reference",
                precision_bits=data["precision_bits"],
            )
        else:
            from experiments.dlr_adapters import build_adapter

            basis = build_adapter(method, data["beta"], cutoff, data["native_eps"])
        self.tau = np.linspace(0, data["beta"], data["n_tau"])
        self.indices = np.arange(-data["n_iw"], data["n_iw"])
        self.observed = np.flatnonzero(
            np.arange(data["n_tau"]) % data["fit_holdout_modulus"]
            != data["fit_holdout_remainder"]
        )
        self.time = np.asarray(basis.design_tau(self.tau), complex)
        self.frequency = np.asarray(basis.design_matsubara(self.indices), complex)
        self.rcond = data["rcond"]

    def apply(self, m: Any) -> np.ndarray:
        coefficients, _, rank, _ = lstsq(
            self.time[self.observed],
            m.g_tau[self.observed],
            cond=self.rcond,
            lapack_driver="gelsy",
        )
        _require(rank > 0, "independent QR retained zero rank")
        green = self.frequency @ coefficients
        # The native tail replacement is fixed by the common native moments,
        # not a noisy fitted tail. Identify its actual changed rows directly
        # from the independently retained pre/post native Fourier arrays.
        mask = np.any(m.g_iw != m.g_iw_raw, axis=1)
        _require(
            np.any(mask) and not np.any(mask[m.config.n_iw - 4 : m.config.n_iw + 4]),
            "native tail window is absent or invades essential low frequencies",
        )
        iw = 1j * (2 * self.indices + 1) * np.pi / m.config.beta
        tail = (
            m.g_moments[1][None, :] / iw[:, None]
            + m.g_moments[2][None, :] / iw[:, None] ** 2
        )
        _require(
            np.max(np.abs(m.g_iw[mask] - tail[mask])) < 1e-12,
            "stored native tail differs from the density-matrix moment expansion",
        )
        green[mask] = tail[mask]
        return green


class _Archive:
    def __init__(self, folder: Path, identity: dict[str, Any]) -> None:
        self.folder, self.identity = folder, identity
        self.seeds: set[int] = set()
        self.digests: dict[str, str] = {}
        self.quadrature_max = 0.0

    def get(
        self, role: str, index: int, seed: int, cycles: int, *, iteration: bool = False
    ) -> Any:
        name = f"iteration-{index:02d}" if iteration else f"replica-{index:02d}"
        stem = self.folder / role / name
        m = load_measurement(stem, self.identity)
        expected_config = {
            name: self.identity["config"][name] for name in dataclasses.asdict(m.config)
        }
        if role == "dense-refinement":
            expected_config["n_iw"] *= 2
            expected_config["n_tau"] = 2 * expected_config["n_tau"] - 1
        _require(
            dataclasses.asdict(m.config) == expected_config,
            "native measurement config differs from its role",
        )
        _require(
            m.seed == seed and m.cycles == cycles,
            "chain seed/cycle allocation mismatch",
        )
        _require(
            seed not in self.seeds, "seed reused across supposedly independent chains"
        )
        self.seeds.add(seed)
        self.digests[str(stem.with_suffix(".json").relative_to(self.folder))] = _sha(
            stem.with_suffix(".json")
        )
        _require(
            m.solve_status == 0 and m.average_sign > 0.99,
            "native solve/sign gate failed",
        )
        _require(
            np.max(np.abs(m.g_iw - m.g_iw[::-1].conj())) < 1e-10,
            "measured scalar Green function violates frequency adjoint symmetry",
        )
        _require(
            np.max(m.g_iw[m.config.n_iw :].imag) < 1e-8,
            "measured Green function violates causality",
        )
        _require(
            np.max(
                np.abs(m.sigma_moments[0] - m.config.interaction * m.occupancies[::-1])
            )
            < 1e-10,
            "native opposite-spin Hartree moment mismatch",
        )
        _require(
            np.max(
                np.abs(m.g_moments[2] - (m.sigma_moments[0] - m.config.interaction / 2))
            )
            < 1e-10,
            "native Green/static-mu moment mismatch",
        )
        # Direct composite trapezoid has a different formulation from the
        # native known-moment transform. Its finite-grid discrepancy is kept
        # separate from stochastic observable agreement.
        labels = np.array(self.identity["config"]["selected_indices"])
        tau = np.linspace(0, m.config.beta, m.config.n_tau)
        phase = np.exp(
            1j * (2 * labels[:, None] + 1) * np.pi / m.config.beta * tau[None, :]
        )
        quadrature = np.trapezoid(phase[:, :, None] * m.g_tau[None, :, :], tau, axis=1)
        error = float(np.max(np.abs(quadrature - m.g_iw_raw[labels + m.config.n_iw])))
        self.quadrature_max = max(self.quadrature_max, error)
        return m


def _assert_bath(actual: np.ndarray, expected: np.ndarray) -> None:
    _require(
        actual.shape == expected.shape and np.max(np.abs(actual - expected)) < 2e-11,
        "saved input bath does not match independent update/averaging",
    )


def _initial(data: dict[str, Any]) -> np.ndarray:
    n = np.arange(-data["n_iw"], data["n_iw"])
    w = (2 * n + 1) * np.pi / data["beta"]
    # Independent positive-frequency algebraic semicircle formula, continued
    # with the explicitly physical scalar adjoint relation for this model.
    return (
        -2j
        * data["hopping"] ** 2
        * np.sign(w)
        / (np.abs(w) + np.sqrt(w * w + 4 * data["hopping"] ** 2))
    )


def _loop(
    archive: _Archive,
    data: dict[str, Any],
    role: str,
    seed: int,
    cycles: int,
    replay: _QRReplay | None,
    *,
    pilot: bool = False,
) -> tuple[np.ndarray, bool, list[Any]]:
    paths = sorted((archive.folder / role).glob("iteration-*.json"))
    maximum = data["calibration_loop_iterations"] if pilot else data["max_iterations"]
    _require(
        1 <= len(paths) <= maximum,
        "loop archive count exceeds or misses its frozen allocation",
    )
    if pilot:
        _require(len(paths) == maximum, "calibration pilot has incomplete iterations")
    expected = _initial(data)
    baths = []
    samples = []
    consecutive = 0
    converged = False
    for i in range(len(paths)):
        m = archive.get(role, i, seed + 2 * i + 1, cycles, iteration=True)
        _assert_bath(m.bath_iw, expected)
        green = m.g_iw if replay is None else replay.apply(m)
        candidate = data["hopping"] ** 2 * green.mean(axis=1)
        defect = float(np.max(np.abs(candidate - m.bath_iw)))
        consecutive = consecutive + 1 if defect < data["loop_tolerance"] else 0
        reached = (
            not pilot
            and i + 1 >= data["min_iterations"]
            and consecutive >= data["consecutive_passes"]
        )
        _require(
            not reached or i == len(paths) - 1,
            "controller continued after its frozen stopping rule",
        )
        converged = reached
        expected = (1 - data["mixing"]) * m.bath_iw + data["mixing"] * candidate
        baths.append(m.bath_iw)
        samples.append(m)
    if not pilot and not converged:
        _require(len(paths) == maximum, "loop stopped early without convergence")
    return np.mean(baths[-data["final_bath_average"] :], axis=0), converged, samples


def _replicas(
    archive: _Archive,
    role: str,
    seed: int,
    cycles: int,
    count: int,
    bath: np.ndarray,
    indices: list[int],
    replay: _QRReplay | None = None,
) -> tuple[list[Any], dict[str, Any]]:
    _require(
        len(list((archive.folder / role).glob("replica-*.json"))) == count,
        "replica population differs from frozen allocation",
    )
    samples = []
    rows = []
    fingerprints = set()
    for i in range(count):
        fingerprint = _read(archive.folder / role / f"replica-{i:02d}.json")[
            "arrays_sha256"
        ]
        _require(
            fingerprint not in fingerprints,
            "duplicate native measurement arrays invalidate independent replicas",
        )
        fingerprints.add(fingerprint)
        m = archive.get(role, i, seed + 2 * i + 1, cycles)
        _assert_bath(m.bath_iw, bath)
        samples.append(m)
        rows.append(
            _observables(m, indices, None if replay is None else replay.apply(m))
        )
    return samples, _statistics(np.array(rows))


def _agreement(
    a: dict[str, Any], b: dict[str, Any], data: dict[str, Any]
) -> dict[str, Any]:
    difference = np.abs(a["mean"] - b["mean"])
    combined = np.sqrt(a["standard_error"] ** 2 + b["standard_error"] ** 2)
    limits = data["agreement_absolute"] + data["agreement_sigma"] * combined
    return {
        "absolute_difference": difference,
        "combined_standard_error": combined,
        "limits": limits,
        "passed": bool(np.all(difference <= limits)),
    }


def audit(folder: Path) -> dict[str, Any]:
    identity = _read(folder / "manifest.json")
    data = identity["config"]
    _require(
        identity.get("schema") == "kernelgecp.cthyb-campaign.v1",
        "wrong campaign schema",
    )
    _require(
        data["canonical"] is True
        and data["cycle_budgets"] == [100000, 400000, 1600000],
        "canonical EE audit refuses a tuning pilot",
    )
    _require(
        data["beta"] == 10 and data["interaction"] == 2 and data["hopping"] == 1,
        "canonical EE physical model mismatch",
    )
    _require(
        data["max_standard_error"] <= 3e-4
        and data["agreement_absolute"] <= 1e-3
        and data["agreement_sigma"] == 3,
        "canonical statistical gates weakened",
    )
    _require(
        not set(data["calibration_indices"]) & set(data["selected_indices"]),
        "calibration and final Green queries overlap",
    )
    _require(
        all(
            identity["versions"].get(key) == value
            for key, value in {
                "triqs": "4.0.2",
                "cthyb": "4.0.1",
                "pydlr": "1.0.1",
                "sparse-ir": "2.1.6",
            }.items()
        ),
        "unsupported generation provider tuple",
    )
    for name, digest in identity["sources"].items():
        blob = subprocess.check_output(
            ["git", "show", identity["source_revision"] + ":" + name], cwd=ROOT
        )
        _require(
            hashlib.sha256(blob).hexdigest() == digest,
            "generation source blob hash mismatch",
        )
    frozen = _read(folder / "freeze.json")
    _require(frozen["identity"] == identity, "frozen identity mismatch")
    for relative, digest in frozen["calibration_archives"].items():
        path = Path(relative)
        _require(
            not path.is_absolute() and ".." not in path.parts, "unsafe evidence path"
        )
        _require(_sha(folder / path) == digest, "frozen calibration file changed")
    saved = _read(folder / "result.json")
    _require(
        saved["identity"] == identity
        and saved["freeze_sha256"] == _sha(folder / "freeze.json"),
        "final result identity/freeze mismatch",
    )
    archive = _Archive(folder, identity)
    audit_revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    audit_relative = "experiments/audit_cthyb_dmft.py"
    audit_hash = _sha(ROOT / audit_relative)
    audit_blob = subprocess.check_output(
        ["git", "show", audit_revision + ":" + audit_relative], cwd=ROOT
    )
    _require(
        hashlib.sha256(audit_blob).hexdigest() == audit_hash,
        "independent audit code must be committed before canonical acceptance",
    )
    _require(
        data.get("seed_mapping") == "odd_stride_2_v1" and data["seed_base"] % 2 == 0,
        "native seed mapping is not validated; nominal independence is invalid",
    )
    seed = data["seed_base"]
    if "revalidation" in identity:
        parent = audit(folder / "parent")
        prior = _read(folder / "parent-audit.json")
        proof = frozen["parent_proof"]
        _require(
            parent["status"] == prior["status"] == "inconclusive",
            "failed parent must remain inconclusive",
        )
        _require(
            parent["generation_revision"]
            == identity["revalidation"]["parent_generation_revision"],
            "parent revision mismatch",
        )
        _require(
            _sha(folder / "parent-audit.json") == proof["independent_audit_sha256"],
            "parent audit hash mismatch",
        )
        budget_records = parent["calibration_budgets"]
        calibration_replay = parent["calibration_QR_replay"]
        chosen = next((row["cycles"] for row in budget_records if row["passed"]), None)
        first_cutoff = next(
            (
                float(cutoff)
                for cutoff, errors in calibration_replay.items()
                if all(
                    value <= data["replay_selected_target"] for value in errors.values()
                )
            ),
            None,
        )
        _require(
            chosen == frozen["cycles"] and first_cutoff == frozen["cutoff"],
            "parent calibration selection mismatch",
        )
        replay_cache = {
            cutoff: {method: _QRReplay(method, data, cutoff) for method in METHODS}
            for cutoff in data["cutoffs"]
        }
    else:
        calibration_bath, _, _ = _loop(
            archive,
            data,
            "calibration-loop",
            seed,
            data["calibration_loop_cycles"],
            None,
            pilot=True,
        )
        budget_records = []
        chosen = None
        calibration_samples = []
        for level, cycles in enumerate(data["cycle_budgets"]):
            samples, estimate = _replicas(
                archive,
                f"calibration-{cycles}",
                seed + 1000 * (level + 1),
                cycles,
                data["calibration_replicas"],
                calibration_bath,
                data["calibration_indices"],
            )
            maximum = float(np.max(estimate["standard_error"]))
            passed = maximum <= data["max_standard_error"]
            budget_records.append(
                {"cycles": cycles, "max_standard_error": maximum, "passed": passed}
            )
            if passed and chosen is None:
                chosen = cycles
                calibration_samples = samples
        _require(
            chosen == frozen["cycles"],
            "frozen budget is not the first independently passing budget",
        )
        replay_cache = {
            cutoff: {method: _QRReplay(method, data, cutoff) for method in METHODS}
            for cutoff in data["cutoffs"]
        }
        first_cutoff = None
        calibration_replay = {}
        positions = np.array(data["calibration_indices"]) + data["n_iw"]
        for cutoff, providers in replay_cache.items():
            errors = {}
            for method, replay in providers.items():
                difference = np.mean(
                    [replay.apply(m) - m.g_iw for m in calibration_samples], axis=0
                )
                errors[method] = float(np.max(np.abs(difference[positions])))
            calibration_replay[str(cutoff)] = errors
            if first_cutoff is None and all(
                x <= data["replay_selected_target"] for x in errors.values()
            ):
                first_cutoff = cutoff
        _require(
            first_cutoff == frozen["cutoff"],
            "frozen cutoff is not the first independently passing cutoff",
        )
    providers = replay_cache[frozen["cutoff"]]
    final_samples = {}
    final_estimates = {}
    final_baths = {}
    residuals = {}
    loop_pass = {}
    for name, offset in (("dense", 10000), ("kernelgecp", 20000)):
        replay = None if name == "dense" else providers["kernelgecp"]
        bath, converged, _ = _loop(
            archive, data, name + "-loop", seed + offset, chosen, replay
        )
        final_baths[name] = bath
        loop_pass[name] = converged
        samples, estimate = _replicas(
            archive,
            name + "-final",
            seed + offset + 20000,
            chosen,
            data["final_replicas"],
            bath,
            data["selected_indices"],
            replay,
        )
        final_samples[name], final_estimates[name] = samples, estimate
        for field in ("mean", "standard_error", "covariance"):
            stored = np.asarray(saved["estimates"][name][field], float)
            _require(
                stored.shape == estimate[field].shape
                and np.max(np.abs(stored - estimate[field])) < 2e-10,
                "reported estimate differs from independent chain recomputation",
            )
        native_estimate = _statistics(
            np.array([_observables(m, data["selected_indices"]) for m in samples])
        )
        positions = np.array(data["selected_indices"]) + data["n_iw"]
        expected = bath[positions] / data["hopping"] ** 2
        target = np.r_[expected.real, expected.imag]
        checks = []
        defect_details = []
        for stats in (native_estimate, estimate):
            difference = np.abs(stats["mean"][2:] - target)
            limits = data["self_consistency_absolute"] + 3 * stats["standard_error"][2:]
            checks.append(bool(np.all(difference <= limits)))
            defect_details.append(
                {
                    "absolute_equation_defect": difference,
                    "monte_carlo_standard_error": stats["standard_error"][2:],
                    "acceptance_limit": limits,
                }
            )
        residuals[name] = {
            "native_and_represented_pass": all(checks),
            "native": defect_details[0],
            "represented": defect_details[1],
        }
    agreement = _agreement(
        final_estimates["dense"], final_estimates["kernelgecp"], data
    )
    maximum_se = max(
        float(np.max(stats["standard_error"])) for stats in final_estimates.values()
    )
    replay_errors = {}
    positions = np.array(data["selected_indices"]) + data["n_iw"]
    for cutoff, methods in replay_cache.items():
        replay_errors[str(cutoff)] = {}
        for trajectory, samples in final_samples.items():
            replay_errors[str(cutoff)][trajectory] = {
                method: float(
                    np.max(
                        np.abs(
                            np.mean(
                                [replay.apply(m) - m.g_iw for m in samples], axis=0
                            )[positions]
                        )
                    )
                )
                for method, replay in methods.items()
            }
    old = final_baths["dense"]
    n = np.arange(-data["n_iw"], data["n_iw"])
    iw = 1j * (2 * n + 1) * np.pi / data["beta"]
    tail_rows = np.abs(n + 0.5) >= 0.8 * data["n_iw"]
    second = float(np.mean(((old - data["hopping"] ** 2 / iw) * iw**2)[tail_rows]).real)
    new_n = np.arange(-2 * data["n_iw"], 2 * data["n_iw"])
    new_iw = 1j * (2 * new_n + 1) * np.pi / data["beta"]
    extended = data["hopping"] ** 2 / new_iw + second / new_iw**2
    extended[data["n_iw"] : 3 * data["n_iw"]] = old
    _, refined_estimate = _replicas(
        archive,
        "dense-refinement",
        seed + 50000,
        chosen,
        data["refinement_replicas"],
        extended,
        data["selected_indices"],
    )
    refinement = _agreement(final_estimates["dense"], refined_estimate, data)
    # Every retained final chain must be present in the saved receipt inventory;
    # acceptance itself remains independently recomputed from its arrays.
    _require(
        saved["measurement_archives"] == archive.digests,
        "saved chain inventory incomplete or changed",
    )
    gates = {
        "loop_controllers": all(loop_pass.values()),
        "standard_errors": maximum_se <= 3e-4,
        "physical_observable_agreement": agreement["passed"],
        "native_and_represented_self_consistency": all(
            row["native_and_represented_pass"] for row in residuals.values()
        ),
        "all_provider_replay_and_cutoff": all(
            error <= data["replay_selected_target"]
            for cut in replay_errors.values()
            for trajectory in cut.values()
            for error in trajectory.values()
        ),
        "independent_native_refinement": refinement["passed"],
        "independent_low_frequency_quadrature": archive.quadrature_max <= 1e-3,
    }
    return {
        "schema": "kernelgecp.cthyb-audit.v1",
        "status": "passed" if all(gates.values()) else "inconclusive",
        "generation_revision": identity["source_revision"],
        "audit_revision": audit_revision,
        "audit_sha256": audit_hash,
        "gates": gates,
        "calibration_budgets": budget_records,
        "calibration_QR_replay": calibration_replay,
        "final_estimates": final_estimates,
        "agreement": agreement,
        "max_standard_error": maximum_se,
        "fixed_bath_residuals": residuals,
        "QR_replay_discrepancies": replay_errors,
        "native_grid_refinement": refinement,
        "quadrature_max": archive.quadrature_max,
        "independent_chain_count": len(archive.seeds),
        "independence_scope": (
            "disjoint_seeds_fresh_thermalized_fixed_bath_chains_not_paired_DMFT_paths"
        ),
        "statistical_formulation": "scalar_two_pass_sample_covariance_of_chain_means",
        "representation_formulation": "independent_scipy_pivoted_QR_gelsy",
    }


def _plain(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_plain(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(args.folder)
    except Exception as error:
        result = {
            "schema": "kernelgecp.cthyb-audit.v1",
            "status": "failed",
            "error_type": type(error).__name__,
            "error": str(error),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(_plain(result), sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
