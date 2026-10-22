"""Independent selected-core and held-out cross checks for schema-2 census data.

These are higher-precision numerical checks, not interval certificates. The
held-out check evaluates the *stored* selected cross on representative interior
points of every order-32 dyadic/endpoint panel; it does not rerun denser GECP.
An optional resolution dataset compares separately completed full-grid GECP runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import subprocess
from collections.abc import Mapping, Sequence
from decimal import Decimal, InvalidOperation
from functools import cache
from importlib.metadata import version
from pathlib import Path
from typing import Any

import mpmath as mp

from kernelgecp import grids as imported_grids
from kernelgecp.grids import chebyshev_lobatto, composite_chebyshev_grids

CORE_TOLERANCE = "1e-18"
SIGMA_TOLERANCE = "1e-15"
REPRESENTATIVE_CUTOFFS = ("1", "1000", "1000000")
HISTORIES = (
    "pivots",
    "t_nodes",
    "omega_nodes",
    "core_det_history",
    "core_sigma_min",
)
GECP_CONTROLS = {
    "tol",
    "max_rank",
    "pivot",
    "precision_bits",
    "grid_order",
    "certificate_abs_tol",
    "certificate_rel_tol",
    "max_cells",
    "tie_rtol",
}
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _source_configuration(record: Mapping[str, Any]) -> Mapping[str, Any]:
    metadata = record.get("metadata")
    if not isinstance(metadata, dict) or not isinstance(
        metadata.get("configuration"), dict
    ):
        raise ValueError("checks require source configuration metadata")
    return metadata["configuration"]


def _validate_sampled_grid(
    record: Mapping[str, Any], *, require_explicit: bool = False
) -> None:
    """Reject contradictory labels without requiring optional fixture metadata."""
    metadata = record.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    configuration = metadata.get("configuration", {})
    if not isinstance(configuration, dict):
        configuration = {}
    labels = (
        (record.get("algorithm"), "grid"),
        (metadata.get("algorithm"), "gecp_grid"),
        (configuration.get("pivot"), "grid"),
        (metadata.get("convergence_scope"), "sampled_grid"),
    )
    if any(value is not None and value != expected for value, expected in labels):
        raise ValueError("finite-grid checks require sampled grid GECP, not adaptive")
    if require_explicit and not any(value is not None for value, _ in labels[:3]):
        raise ValueError("resolution comparison requires explicit grid algorithm")


@cache
def _revision_grid_sha256(revision: str) -> str | None:
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        return None
    completed = subprocess.run(
        ["git", "show", f"{revision}:kernelgecp/grids.py"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    return hashlib.sha256(completed.stdout).hexdigest()


def _record_provenance(record: Mapping[str, Any]) -> dict[str, Any]:
    revision = record.get("git_commit")
    grid_hash = _revision_grid_sha256(revision) if isinstance(revision, str) else None
    return {
        "git_commit": revision,
        "git_dirty": record.get("git_dirty"),
        "config_hash": record.get("config_hash"),
        "revision_grid_source_sha256": grid_hash,
        "revision_grid_source_status": "recorded_revision_blob"
        if grid_hash is not None
        else "unavailable_recorded_revision_blob",
        "recorded_environment": record.get("environment"),
    }


def _decimal(value: Any, name: str) -> Decimal:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a decimal string")
    try:
        number = Decimal(value)
    except InvalidOperation:
        raise ValueError(f"{name} must be a finite decimal string") from None
    if not number.is_finite():
        raise ValueError(f"{name} must be finite")
    return number


def _integer(value: Any, name: str, minimum: int) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _validate_record(record: Mapping[str, Any]) -> None:
    if record.get("schema_version") != 2:
        raise ValueError("independent checks require schema_version 2")
    rank = _integer(record.get("rank"), "rank", 0)
    _integer(record.get("precision_bits"), "precision_bits", 53)
    cutoff = _decimal(record.get("cutoff"), "cutoff")
    tolerance = _decimal(record.get("tolerance"), "tolerance")
    if cutoff < 1 or not 0 < tolerance < 1:
        raise ValueError("require cutoff >= 1 and 0 < tolerance < 1")
    if not math.isfinite(float(cutoff)):
        raise ValueError("cutoff must be representable as a finite float")
    for name in (*HISTORIES, "residual_history"):
        values = record.get(name)
        expected = rank + (name == "residual_history")
        if not isinstance(values, list) or len(values) != expected:
            raise ValueError(f"misaligned {name}")
        numbers = [_decimal(value, name) for value in values]
        if name == "t_nodes" and any(not 0 <= value <= 1 for value in numbers):
            raise ValueError("t_nodes must lie in [0, 1]")
        if name == "omega_nodes" and any(abs(value) > cutoff for value in numbers):
            raise ValueError("omega_nodes lie outside the cutoff")
        if name == "residual_history" and any(value < 0 for value in numbers):
            raise ValueError("residual_history must be nonnegative")
        if name == "core_sigma_min" and any(value <= 0 for value in numbers):
            raise ValueError("selected-core singular values must be positive")
        if name in ("pivots", "core_det_history") and any(
            value == 0 for value in numbers
        ):
            raise ValueError("selected pivots and determinants must be nonzero")
    if _decimal(record.get("residual"), "residual") != Decimal(
        record["residual_history"][-1]
    ):
        raise ValueError("terminal residual disagrees with residual_history")
    if (
        len({Decimal(value) for value in record["t_nodes"]}) != rank
        or len({Decimal(value) for value in record["omega_nodes"]}) != rank
    ):
        raise ValueError("selected coordinates must be distinct")


def strictest_records(
    records: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    """Validate all records and select the deepest tolerance per cutoff."""
    if not records:
        raise ValueError("census must contain at least one record")
    grouped: dict[Decimal, list[Mapping[str, Any]]] = {}
    for record in records:
        _validate_record(record)
        grouped.setdefault(Decimal(record["cutoff"]), []).append(record)
    selected = []
    for cutoff in sorted(grouped):
        cases = grouped[cutoff]
        if len({Decimal(case["tolerance"]) for case in cases}) != len(cases):
            raise ValueError("duplicate cutoff/tolerance pair")
        deepest = min(cases, key=lambda case: Decimal(case["tolerance"]))
        for case in cases:
            if case["rank"] > deepest["rank"]:
                raise ValueError("a looser prefix is deeper than the strictest case")
            for name in (*HISTORIES, "residual_history"):
                if case[name] != deepest[name][: len(case[name])]:
                    raise ValueError(f"incompatible prefix {name}")
        selected.append(deepest)
    return selected


def _kernel(t: Any, omega: Any) -> Any:
    """Independent sign-stable formula, with no call to GECP/kernel code."""
    if omega >= 0:
        return mp.exp(-t * omega) / (1 + mp.exp(-omega))
    return mp.exp((1 - t) * omega) / (1 + mp.exp(omega))


def _core(record: Mapping[str, Any], rank: int) -> Any:
    times = [mp.mpf(value) for value in record["t_nodes"][:rank]]
    frequencies = [mp.mpf(value) for value in record["omega_nodes"][:rank]]
    return mp.matrix([[_kernel(t, omega) for omega in frequencies] for t in times])


def _encode(value: Any, bits: int) -> str:
    if not mp.isfinite(value):
        raise ArithmeticError("nonfinite independent check")
    return str(mp.nstr(value, math.ceil(bits * math.log10(2)) + 3))


def _relative_error(observed: Any, reference: Any) -> Any:
    if reference == 0:
        raise ArithmeticError("zero independent reference")
    return abs(observed - reference) / abs(reference)


def check_selected_cores(
    record: Mapping[str, Any], precision_bits: int = 192
) -> dict[str, Any]:
    """Re-evaluate determinant/SVD and signed pivot products at fixed precision."""
    _validate_record(record)
    _integer(precision_bits, "check precision_bits", 54)
    if precision_bits <= record["precision_bits"]:
        raise ValueError("check precision must exceed source precision")
    rank = record["rank"]
    ranks = sorted({1, max(1, rank // 2), rank}) if rank else []
    checks = []
    with mp.workprec(precision_bits):
        for prefix in ranks:
            core = _core(record, prefix)
            determinant = mp.det(core)
            singular_values = mp.svd(core, compute_uv=False)
            sigma = singular_values[prefix - 1]
            if determinant == 0 or sigma <= 0:
                raise ArithmeticError("singular independently evaluated selected core")
            pivot_product = mp.fprod(mp.mpf(v) for v in record["pivots"][:prefix])
            det_error = _relative_error(
                mp.mpf(record["core_det_history"][prefix - 1]), determinant
            )
            product_error = _relative_error(pivot_product, determinant)
            sigma_error = _relative_error(
                mp.mpf(record["core_sigma_min"][prefix - 1]), sigma
            )
            condition = singular_values[0] / sigma
            passed = (
                det_error <= mp.mpf(CORE_TOLERANCE)
                and product_error <= mp.mpf(CORE_TOLERANCE)
                and sigma_error <= mp.mpf(SIGMA_TOLERANCE)
            )
            checks.append(
                {
                    "rank": prefix,
                    "determinant": _encode(determinant, precision_bits),
                    "sigma_min": _encode(sigma, precision_bits),
                    "condition_number_2": _encode(condition, precision_bits),
                    "stored_determinant_relative_error": _encode(
                        det_error, precision_bits
                    ),
                    "signed_pivot_product_relative_error": _encode(
                        product_error, precision_bits
                    ),
                    "stored_sigma_relative_error": _encode(sigma_error, precision_bits),
                    "source_roundoff_indicator": _encode(
                        mp.mpf(2) ** (1 - record["precision_bits"])
                        * prefix**2
                        * condition,
                        precision_bits,
                    ),
                    "status": "passed" if passed else "disagreement",
                }
            )
    return {
        "cutoff": record["cutoff"],
        "tolerance": record["tolerance"],
        "source_rank": rank,
        "checks": checks,
        "status": "not_applicable_rank_zero"
        if rank == 0
        else "passed"
        if all(check["status"] == "passed" for check in checks)
        else "disagreement",
    }


def _panels(cutoff: float) -> tuple[list[tuple[float, float]], ...]:
    first = min(0.5, 1 / cutoff)
    time = [(0.0, first)]
    left = first
    while left < 0.5:
        right = min(0.5, 2 * left)
        time.append((left, right))
        left = right
    time += [(1 - right, 1 - left) for left, right in reversed(time)]
    positive = [(0.0, 1.0)]
    left = 1.0
    while left < cutoff:
        right = min(cutoff, 2 * left)
        positive.append((left, right))
        left = right
    frequency = [(-right, -left) for left, right in reversed(positive)] + positive
    return time, frequency


def panel_samples(
    cutoff: float, source_order: int, denser_order: int = 32
) -> tuple[list[str], list[str], dict[str, Any]]:
    """Take up to three held-out interior quantiles from every denser panel."""
    if not math.isfinite(cutoff) or cutoff < 1:
        raise ValueError("cutoff must be finite and >= 1")
    _integer(source_order, "source_order", 2)
    _integer(denser_order, "denser_order", source_order + 1)
    source_grids = composite_chebyshev_grids(cutoff, source_order)
    denser_grids = composite_chebyshev_grids(cutoff, denser_order)
    panels = _panels(cutoff)
    samples = []
    panel_counts = []
    for intervals, source in zip(panels, source_grids, strict=True):
        excluded = set(source)
        values: set[float] = set()
        counts = []
        for left, right in intervals:
            candidates = [
                float(value)
                for value in chebyshev_lobatto(left, right, denser_order)[1:-1]
                if float(value) not in excluded
            ]
            if not candidates:
                raise ValueError("panel has no representable held-out interior node")
            count = min(3, len(candidates))
            chosen = {
                candidates[(index + 1) * len(candidates) // (count + 1)]
                for index in range(count)
            }
            values.update(chosen)
            counts.append(len(chosen))
        samples.append([repr(value) for value in sorted(values)])
        panel_counts.append(counts)
    return (
        samples[0],
        samples[1],
        {
            "source_order": source_order,
            "denser_order": denser_order,
            "source_grid_shape": [len(grid) for grid in source_grids],
            "full_denser_grid_shape": [len(grid) for grid in denser_grids],
            "panel_counts": [len(axis) for axis in panels],
            "samples_per_panel": panel_counts,
            "sample_shape": [len(axis) for axis in samples],
            "sample_count": len(samples[0]) * len(samples[1]),
            "policy": (
                "Each canonical dyadic/endpoint panel: order-N Chebyshev-Lobatto "
                "interiors, remove exact float64 source-axis coordinates, choose "
                "up to 3 remaining ascending-index quantiles at floor((j+1)*n/(m+1)); "
                "take the Cartesian product of the unique axis samples."
            ),
            "coordinate_precision_bits": 53,
        },
    )


def _batched_lu_solve(core: Any, rhs: Any) -> Any:
    """Factor once and solve all columns by triangular substitution, no inverse."""
    permutation, lower, upper = mp.lu(core)
    size = core.rows
    forward = mp.matrix(size, rhs.cols)
    solution = mp.matrix(size, rhs.cols)
    rows = [next(j for j in range(size) if permutation[i, j]) for i in range(size)]
    for i in range(size):
        for column in range(rhs.cols):
            forward[i, column] = (
                rhs[rows[i], column]
                - mp.fsum(lower[i, j] * forward[j, column] for j in range(i))
            ) / lower[i, i]
    for i in range(size - 1, -1, -1):
        for column in range(rhs.cols):
            solution[i, column] = (
                forward[i, column]
                - mp.fsum(upper[i, j] * solution[j, column] for j in range(i + 1, size))
            ) / upper[i, i]
    return solution


def check_held_out_cross(
    record: Mapping[str, Any], precision_bits: int = 192, denser_order: int = 32
) -> dict[str, Any]:
    """Evaluate the final stored selected cross on a panel-spanning held-out set."""
    _validate_record(record)
    _validate_sampled_grid(record)
    _integer(precision_bits, "check precision_bits", 54)
    if precision_bits <= record["precision_bits"]:
        raise ValueError("check precision must exceed source precision")
    configuration = _source_configuration(record)
    source_order = _integer(configuration.get("grid_order"), "source grid_order", 2)
    times, frequencies, sampling = panel_samples(
        float(record["cutoff"]), source_order, denser_order
    )
    with mp.workprec(precision_bits):
        t_values = [mp.mpf(value) for value in times]
        w_values = [mp.mpf(value) for value in frequencies]
        rank = record["rank"]
        solve_error = None
        if rank:
            core = _core(record, rank)
            selected_t = [mp.mpf(value) for value in record["t_nodes"]]
            selected_w = [mp.mpf(value) for value in record["omega_nodes"]]
            rhs = mp.matrix([[_kernel(t, w) for w in w_values] for t in selected_t])
            solution = _batched_lu_solve(core, rhs)
            left = mp.matrix([[_kernel(t, w) for w in selected_w] for t in t_values])
            approximation = left * solution
            solve_error = max(abs(value) for value in core * solution - rhs) / max(
                abs(value) for value in rhs
            )
        else:
            approximation = mp.matrix(len(t_values), len(w_values))
        peak = mp.mpf(0)
        worst = (0, 0)
        for i, t in enumerate(t_values):
            for j, omega in enumerate(w_values):
                error = abs(_kernel(t, omega) - approximation[i, j])
                if error > peak:
                    peak, worst = error, (i, j)
        source = mp.mpf(record["residual"])
        tolerance = mp.mpf(record["tolerance"])
        return {
            "cutoff": record["cutoff"],
            "tolerance": record["tolerance"],
            "rank": rank,
            "sampling": sampling,
            "t_samples": times,
            "omega_samples": frequencies,
            "held_out_max_abs_residual": _encode(peak, precision_bits),
            "worst_coordinates": [times[worst[0]], frequencies[worst[1]]],
            "source_terminal_residual": record["residual"],
            "ratio_to_source_terminal_residual": _encode(peak / source, precision_bits)
            if source
            else None,
            "ratio_status": "defined" if source else "undefined_source_zero",
            "exceeds_source_terminal_residual": bool(peak > source),
            "within_requested_tolerance_on_sample": bool(peak <= tolerance),
            "solve_relative_max_entry_residual": _encode(solve_error, precision_bits)
            if solve_error is not None
            else None,
            "solve_method": "mpmath_lu_once_batched_forward_backward_substitution"
            if rank
            else "not_applicable_rank_zero",
            "status": "observed_held_out_sample",
            "certification_status": "not_certified",
            "scope": "stored_selected_cross_on_panel_spanning_held_out_sample",
            "not_claimed": [
                "denser_GECP_run",
                "full_denser_grid_max",
                "continuous_bound",
            ],
        }


def _read_records(data: bytes) -> list[Mapping[str, Any]]:
    records = []
    for line in data.decode("utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            if not isinstance(record, dict):
                raise ValueError("census rows must be JSON objects")
            records.append(record)
    return records


def _resolution_case(record: Mapping[str, Any]) -> dict[str, Any]:
    configuration = _source_configuration(record)
    order = _integer(configuration.get("grid_order"), "source grid_order", 2)
    shape = record.get("grid_shape")
    expected_shape = [
        len(axis) for axis in composite_chebyshev_grids(float(record["cutoff"]), order)
    ]
    if shape != expected_shape:
        raise ValueError("recorded grid_shape disagrees with effective grid controls")
    ties = record.get("tie_counts")
    if (
        not isinstance(ties, list)
        or len(ties) != record["rank"]
        or any(type(value) is not int or value < 1 for value in ties)
    ):
        raise ValueError("resolution comparison requires aligned positive tie_counts")
    if type(record.get("converged")) is not bool or not isinstance(
        record.get("stop_reason"), str
    ):
        raise ValueError("resolution comparison requires recorded convergence and stop")
    metadata = record["metadata"]
    return {
        "provenance": _record_provenance(record),
        "effective_configuration": dict(configuration),
        "grid_order": order,
        "grid_shape": shape,
        "grid_coordinate_precision_bits": record.get("grid_coordinate_precision_bits"),
        "grid_generation": record.get("grid_generation"),
        "rank": record["rank"],
        "residual": record["residual"],
        "tolerance": record["tolerance"],
        "precision_bits": record["precision_bits"],
        "max_tie_count": max(ties) if ties else None,
        "converged": record["converged"],
        "stop_reason": record["stop_reason"],
        "recorded_convergence_scope": metadata.get("convergence_scope"),
        "recorded_certification_scope": metadata.get("certification_scope"),
        "comparison_sampling_scope": "each_run_full_finite_training_grid_only",
    }


def compare_resolution_records(
    source_trajectories: Sequence[Mapping[str, Any]],
    resolution_trajectories: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Compare actual GECP runs at changed resolution, without rank monotonicity."""
    source = {Decimal(record["cutoff"]): record for record in source_trajectories}
    comparisons = []
    for denser in resolution_trajectories:
        _validate_record(denser)
        cutoff = Decimal(denser["cutoff"])
        if cutoff not in source:
            raise ValueError("resolution case has no canonical source counterpart")
        original = source[cutoff]
        _validate_record(original)
        for record in (original, denser):
            _validate_sampled_grid(record, require_explicit=True)
            configuration = _source_configuration(record)
            if not GECP_CONTROLS.issubset(configuration):
                raise ValueError(
                    "resolution comparison requires all effective GECP controls"
                )
            if (
                record["precision_bits"] != 128
                or configuration["precision_bits"] != 128
            ):
                raise ValueError("resolution comparison requires 128-bit source runs")
            if Decimal(str(configuration["tol"])) != Decimal(record["tolerance"]):
                raise ValueError(
                    "effective tolerance disagrees with recorded tolerance"
                )
        original_configuration = dict(_source_configuration(original))
        denser_configuration = dict(_source_configuration(denser))
        original_order = _integer(
            original_configuration.pop("grid_order"), "source grid_order", 2
        )
        denser_order = _integer(
            denser_configuration.pop("grid_order"), "resolution grid_order", 2
        )
        if denser_order <= original_order:
            raise ValueError("resolution grid_order must exceed source grid_order")
        if (
            Decimal(original["tolerance"]) != Decimal(denser["tolerance"])
            or original_configuration != denser_configuration
        ):
            raise ValueError("effective GECP controls differ beyond grid_order")
        original_case, denser_case = (
            _resolution_case(record) for record in (original, denser)
        )
        common_prefix = 0
        for t, w, other_t, other_w in zip(
            original["t_nodes"],
            original["omega_nodes"],
            denser["t_nodes"],
            denser["omega_nodes"],
            strict=False,
        ):
            if (Decimal(t), Decimal(w)) != (Decimal(other_t), Decimal(other_w)):
                break
            common_prefix += 1
        original_grid_hash = original_case["provenance"]["revision_grid_source_sha256"]
        denser_grid_hash = denser_case["provenance"]["revision_grid_source_sha256"]
        comparisons.append(
            {
                "cutoff": original["cutoff"],
                "source": original_case,
                "resolution": denser_case,
                "controls_compatible_except_grid_order": True,
                "source_grid_implementation_hashes_equal": (
                    original_grid_hash == denser_grid_hash
                    if original_grid_hash is not None and denser_grid_hash is not None
                    else None
                ),
                "rank_difference_resolution_minus_source": denser["rank"]
                - original["rank"],
                "residual_comparison": (
                    "resolution_smaller"
                    if Decimal(denser["residual"]) < Decimal(original["residual"])
                    else "equal"
                    if Decimal(denser["residual"]) == Decimal(original["residual"])
                    else "resolution_larger"
                ),
                "common_selected_prefix_length": common_prefix,
                "common_prefix_definition": "exact_decimal_node_pairs_in_pivot_order",
                "status": "observed_resolution_comparison",
                "certification_status": "not_certified",
                "scope": "separate_full_finite_grid_GECP_runs",
                "not_claimed": [
                    "monotone_rank_under_grid_refinement",
                    "nested_grids",
                    "continuous_residual_bound",
                    "grid_independent_pivot_trajectory",
                ],
            }
        )
    return comparisons


def check_census(
    input_path: Path,
    *,
    precision_bits: int = 192,
    denser_order: int = 32,
    representative_cutoffs: Sequence[str] = REPRESENTATIVE_CUTOFFS,
    resolution_input_path: Path | None = None,
) -> dict[str, Any]:
    """Return deterministic independent diagnostics, including any disagreement."""
    data = input_path.read_bytes()
    records = _read_records(data)
    trajectories = strictest_records(records)
    resolution = None
    resolution_trajectories = []
    if resolution_input_path is not None:
        resolution_data = resolution_input_path.read_bytes()
        resolution_records = _read_records(resolution_data)
        resolution_trajectories = strictest_records(resolution_records)
        resolution = {
            "input_sha256": hashlib.sha256(resolution_data).hexdigest(),
            "source_record_count": len(resolution_records),
            "strictest_trajectory_count": len(resolution_trajectories),
            "comparisons": compare_resolution_records(
                trajectories, resolution_trajectories
            ),
            "scope": "actual_GECP_resolution_ladder_not_held_out_cross_re_evaluation",
            "certification_status": "not_certified",
        }
        for comparison in resolution["comparisons"]:
            comparison["source"]["provenance"]["input_sha256"] = hashlib.sha256(
                data
            ).hexdigest()
            comparison["resolution"]["provenance"]["input_sha256"] = resolution[
                "input_sha256"
            ]
    representatives = {
        _decimal(value, "representative cutoff") for value in representative_cutoffs
    }
    if any(value < 1 for value in representatives):
        raise ValueError("representative cutoffs must be >= 1")
    core_checks = [
        check_selected_cores(record, precision_bits) for record in trajectories
    ]
    held_out = [
        check_held_out_cross(record, precision_bits, denser_order)
        for record in trajectories
        if Decimal(record["cutoff"]) in representatives
    ]
    resolution_core_checks = [
        check_selected_cores(record, precision_bits)
        for record in resolution_trajectories
    ]
    if resolution is not None:
        resolution["selected_core_checks"] = resolution_core_checks
        resolution["status"] = (
            "passed_numerical_core_checks"
            if all(case["status"] != "disagreement" for case in resolution_core_checks)
            else "numerical_core_disagreement"
        )
    return {
        "schema_version": 1,
        "input_sha256": hashlib.sha256(data).hexdigest(),
        "checker_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "imported_grid_source_sha256": hashlib.sha256(
            Path(imported_grids.__file__).read_bytes()
        ).hexdigest(),
        "environment": {
            "python": platform.python_version(),
            "mpmath": version("mpmath"),
            "numpy": version("numpy"),
            "mpmath_backend": str(mp.libmp.BACKEND),
            "uv_lock_sha256": hashlib.sha256(
                (PROJECT_ROOT / "uv.lock").read_bytes()
            ).hexdigest(),
        },
        "precision_bits": precision_bits,
        "source_record_count": len(records),
        "strictest_trajectory_count": len(trajectories),
        "checked_cutoffs": [record["cutoff"] for record in trajectories],
        "core_check_ranks": (
            "unique sorted {1, max(1, floor(rank/2)), rank}; empty at rank 0"
        ),
        "relative_error_thresholds": {
            "stored_determinant": CORE_TOLERANCE,
            "signed_pivot_product": CORE_TOLERANCE,
            "stored_sigma_min": SIGMA_TOLERANCE,
        },
        "threshold_scope": (
            "Fixed numerical acceptance thresholds, never adapted to observed error. "
            "For context only, "
            "source_roundoff_indicator=2^(1-source_bits)*rank^2*kappa2; "
            "it is a heuristic sensitivity indicator, not a proved error bound."
        ),
        "selected_core_checks": core_checks,
        "held_out_cross_checks": held_out,
        "resolution_ladder": resolution,
        "unavailable_representative_cutoffs": [
            str(value)
            for value in sorted(representatives)
            if value not in {Decimal(record["cutoff"]) for record in trajectories}
        ],
        "status": "passed_numerical_core_checks"
        if all(
            case["status"] != "disagreement"
            for case in (*core_checks, *resolution_core_checks)
        )
        else "numerical_core_disagreement",
        "certification_status": "not_certified",
        "scope": "independent_point_arithmetic_and_held_out_selected_cross_evidence",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--precision-bits", type=int, default=192)
    parser.add_argument("--denser-order", type=int, default=32)
    parser.add_argument("--resolution-input", type=Path)
    args = parser.parse_args()
    report = check_census(
        args.input,
        precision_bits=args.precision_bits,
        denser_order=args.denser_order,
        resolution_input_path=args.resolution_input,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(
        json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(args.output)
    if report["status"] == "numerical_core_disagreement":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
