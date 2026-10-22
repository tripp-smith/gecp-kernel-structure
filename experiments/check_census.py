"""Independent selected-core and held-out cross checks for schema-2 census data.

These are higher-precision numerical checks, not interval certificates. The
held-out check evaluates the *stored* selected cross on representative interior
points of every order-32 dyadic/endpoint panel; it does not rerun denser GECP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from collections.abc import Mapping, Sequence
from decimal import Decimal, InvalidOperation
from importlib.metadata import version
from pathlib import Path
from typing import Any

import mpmath as mp

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
    _integer(precision_bits, "check precision_bits", 54)
    if precision_bits <= record["precision_bits"]:
        raise ValueError("check precision must exceed source precision")
    metadata = record.get("metadata")
    if not isinstance(metadata, dict) or not isinstance(
        metadata.get("configuration"), dict
    ):
        raise ValueError("held-out checks require source configuration metadata")
    configuration = metadata["configuration"]
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


def check_census(
    input_path: Path,
    *,
    precision_bits: int = 192,
    denser_order: int = 32,
    representative_cutoffs: Sequence[str] = REPRESENTATIVE_CUTOFFS,
) -> dict[str, Any]:
    """Return deterministic independent diagnostics, including any disagreement."""
    data = input_path.read_bytes()
    records = []
    for line in data.decode("utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            if not isinstance(record, dict):
                raise ValueError("census rows must be JSON objects")
            records.append(record)
    trajectories = strictest_records(records)
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
    return {
        "schema_version": 1,
        "input_sha256": hashlib.sha256(data).hexdigest(),
        "checker_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "environment": {
            "python": platform.python_version(),
            "mpmath": version("mpmath"),
            "numpy": version("numpy"),
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
        "unavailable_representative_cutoffs": [
            str(value)
            for value in sorted(representatives)
            if value not in {Decimal(record["cutoff"]) for record in trajectories}
        ],
        "status": "passed_numerical_core_checks"
        if all(case["status"] != "disagreement" for case in core_checks)
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
    args = parser.parse_args()
    report = check_census(
        args.input, precision_bits=args.precision_bits, denser_order=args.denser_order
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
