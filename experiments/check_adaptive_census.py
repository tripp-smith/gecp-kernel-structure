"""Independent checks of explicitly recorded residual-adaptive sampled grids.

The kernel, panel probes and selected-cross solves below deliberately do not
call the adaptive controller.  All maxima remain finite-sample observations;
neither a successful run nor this checker supplies a continuum certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from collections.abc import Callable, Mapping, Sequence
from decimal import (
    MAX_EMAX,
    MIN_EMIN,
    ROUND_HALF_EVEN,
    Context,
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    localcontext,
)
from functools import wraps
from pathlib import Path
from typing import Any

import mpmath as mp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_CUTOFFS = tuple(Decimal(10) ** index for index in range(7))
CANONICAL_TOLERANCES = tuple(Decimal(10) ** (-index) for index in (6, 8, 10, 12))
HISTORIES = (
    "pivots",
    "t_nodes",
    "omega_nodes",
    "core_det_history",
    "core_sigma_min",
)


def _stable_decimal[**P, T](function: Callable[P, T]) -> Callable[P, T]:
    """Keep checker decisions independent of a caller's decimal settings."""

    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> T:
        context = Context(
            prec=120,
            rounding=ROUND_HALF_EVEN,
            Emin=MIN_EMIN,
            Emax=MAX_EMAX,
            traps=[InvalidOperation, DivisionByZero, Overflow],
        )
        with localcontext(context):
            return function(*args, **kwargs)

    return wrapped


def _decimal(value: Any, name: str, *, nonnegative: bool = False) -> Decimal:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a decimal string")
    try:
        result = Decimal(value)
    except InvalidOperation:
        raise ValueError(f"{name} must be a finite decimal string") from None
    if not result.is_finite() or (nonnegative and result < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def _integer(value: Any, name: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _object(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _list(value: Any, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a list")
    return value


def _encode(value: Any, bits: int) -> str:
    if not mp.isfinite(value):
        raise ArithmeticError("nonfinite independent check")
    return str(mp.nstr(value, math.ceil(bits * math.log10(2)) + 3))


def _kernel(t: Any, omega: Any) -> Any:
    """Independent sign-stable original fermionic kernel."""
    if omega >= 0:
        return mp.exp(-t * omega) / (1 + mp.exp(-omega))
    return mp.exp((1 - t) * omega) / (1 + mp.exp(omega))


@_stable_decimal
def _validate_record(record: Mapping[str, Any]) -> None:
    if record.get("schema_version") != 3:
        raise ValueError("adaptive checks require schema_version 3")
    if record.get("algorithm") != "residual_adaptive_grid":
        raise ValueError("adaptive census trajectories require sampled grid GECP")
    rank = _integer(record.get("rank"), "rank")
    _integer(record.get("precision_bits"), "precision_bits", 54)
    cutoff = _decimal(record.get("cutoff"), "cutoff")
    tolerance = _decimal(record.get("tolerance"), "tolerance")
    if cutoff < 1 or not 0 < tolerance < 1:
        raise ValueError("require cutoff >= 1 and 0 < tolerance < 1")
    for name in (*HISTORIES, "residual_history"):
        values = _list(record.get(name), name)
        if len(values) != rank + (name == "residual_history"):
            raise ValueError(f"misaligned {name}")
        numbers = [_decimal(value, name) for value in values]
        if name == "t_nodes" and any(not 0 <= value <= 1 for value in numbers):
            raise ValueError("selected times outside physical domain")
        coordinate_roundoff = Decimal(2) ** (1 - record["precision_bits"]) * max(
            Decimal(1), cutoff
        )
        if name == "omega_nodes" and any(
            abs(value) > cutoff + coordinate_roundoff for value in numbers
        ):
            raise ValueError("selected frequencies outside physical domain")
        if name == "residual_history" and any(value < 0 for value in numbers):
            raise ValueError("residual_history must be nonnegative")
        if name == "core_sigma_min" and any(value <= 0 for value in numbers):
            raise ValueError("selected-core singular values must be positive")
        if name in ("pivots", "core_det_history") and any(
            value == 0 for value in numbers
        ):
            raise ValueError("pivots and determinants must be nonzero")
    if _decimal(record.get("residual"), "residual") != Decimal(
        record["residual_history"][-1]
    ):
        raise ValueError("terminal residual disagrees with residual_history")
    if any(
        len({Decimal(value) for value in record[name]}) != rank
        for name in ("t_nodes", "omega_nodes")
    ):
        raise ValueError("selected coordinates must be distinct")
    ties = _list(record.get("tie_counts"), "tie_counts")
    if len(ties) != rank or any(type(value) is not int or value < 1 for value in ties):
        raise ValueError("tie_counts must align and be positive")
    if type(record.get("converged")) is not bool:
        raise ValueError("converged must be boolean")
    if not isinstance(record.get("stop_reason"), str) or not record["stop_reason"]:
        raise ValueError("stop_reason must be a nonempty string")


def _core(record: Mapping[str, Any], rank: int) -> Any:
    times = [mp.mpf(value) for value in record["t_nodes"][:rank]]
    frequencies = [mp.mpf(value) for value in record["omega_nodes"][:rank]]
    return mp.matrix([[_kernel(t, omega) for omega in frequencies] for t in times])


def check_selected_cores(
    record: Mapping[str, Any], precision_bits: int = 192
) -> dict[str, Any]:
    """Recompute original selected cores at first, middle and terminal rank."""
    _validate_record(record)
    _integer(precision_bits, "check precision_bits", 192)
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
                raise ArithmeticError("singular independent selected core")
            product = mp.fprod(mp.mpf(v) for v in record["pivots"][:prefix])
            errors = {
                "stored_determinant_relative_error": abs(
                    mp.mpf(record["core_det_history"][prefix - 1]) / determinant - 1
                ),
                "signed_pivot_product_relative_error": abs(product / determinant - 1),
                "stored_sigma_relative_error": abs(
                    mp.mpf(record["core_sigma_min"][prefix - 1]) / sigma - 1
                ),
            }
            passed = all(
                value
                <= mp.mpf("1e-15" if name == "stored_sigma_relative_error" else "1e-18")
                for name, value in errors.items()
            )
            checks.append(
                {
                    "rank": prefix,
                    "determinant": _encode(determinant, precision_bits),
                    "sigma_min": _encode(sigma, precision_bits),
                    "condition_number_2": _encode(
                        singular_values[0] / sigma, precision_bits
                    ),
                    **{
                        name: _encode(value, precision_bits)
                        for name, value in errors.items()
                    },
                    "status": "passed" if passed else "disagreement",
                }
            )
    return {
        "cutoff": record["cutoff"],
        "source_rank": rank,
        "checks": checks,
        "status": "passed"
        if all(c["status"] == "passed" for c in checks)
        else "disagreement",
        "scope": "original_kernel_selected_cores_first_middle_final",
        "certification_status": "not_certified",
    }


def _batched_lu_solve(core: Any, rhs: Any) -> Any:
    """One independent LU factorization followed by vector triangular solves."""
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


def evaluate_cross(
    record: Mapping[str, Any],
    t_samples: Sequence[str],
    omega_samples: Sequence[str],
    precision_bits: int = 192,
) -> dict[str, Any]:
    """Independently evaluate a stored selected cross on a supplied product set."""
    _validate_record(record)
    _integer(precision_bits, "check precision_bits", 192)
    if precision_bits <= record["precision_bits"]:
        raise ValueError("check precision must exceed source precision")
    if not t_samples or not omega_samples:
        raise ValueError("cross check requires nonempty sample axes")
    cutoff = _decimal(record["cutoff"], "cutoff")
    if any(not 0 <= _decimal(t, "sample time") <= 1 for t in t_samples) or any(
        abs(_decimal(w, "sample frequency")) > cutoff for w in omega_samples
    ):
        raise ValueError("cross check samples outside physical domain")
    with mp.workprec(precision_bits):
        ts = [mp.mpf(value) for value in t_samples]
        ws = [mp.mpf(value) for value in omega_samples]
        rank = record["rank"]
        solve_error = None
        if rank:
            core = _core(record, rank)
            selected_t = [mp.mpf(value) for value in record["t_nodes"]]
            selected_w = [mp.mpf(value) for value in record["omega_nodes"]]
            rhs = mp.matrix([[_kernel(t, w) for w in ws] for t in selected_t])
            solution = _batched_lu_solve(core, rhs)
            left = mp.matrix([[_kernel(t, w) for w in selected_w] for t in ts])
            approximation = left * solution
            solve_error = max(abs(value) for value in core * solution - rhs) / max(
                abs(value) for value in rhs
            )
        else:
            approximation = mp.matrix(len(ts), len(ws))
        peak = mp.mpf(0)
        worst = (0, 0)
        for i, t in enumerate(ts):
            for j, omega in enumerate(ws):
                error = abs(_kernel(t, omega) - approximation[i, j])
                if error > peak:
                    peak, worst = error, (i, j)
        return {
            "rank": rank,
            "sample_shape": [len(ts), len(ws)],
            "max_abs_residual": _encode(peak, precision_bits),
            "worst_coordinates": [t_samples[worst[0]], omega_samples[worst[1]]],
            "within_requested_tolerance_on_sample": bool(
                peak <= mp.mpf(record["tolerance"])
            ),
            "solve_relative_max_entry_residual": _encode(solve_error, precision_bits)
            if solve_error is not None
            else None,
            "solve_method": "mpmath_lu_once_batched_forward_backward_substitution"
            if rank
            else "not_applicable_rank_zero",
            "precision_bits": precision_bits,
            "scope": "stored_selected_cross_on_independent_finite_sample",
            "certification_status": "not_certified",
        }


def _panel_axis(raw_panels: Any, cutoff: Decimal, axis: str) -> list[Mapping[str, Any]]:
    panels = [
        _object(value, f"{axis} panel") for value in _list(raw_panels, f"{axis} panels")
    ]
    if not panels:
        raise ValueError("panel axes must be nonempty")
    identifiers = []
    previous = Decimal(0) if axis == "t" else -cutoff
    for panel in panels:
        identifier = panel.get("identifier")
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("panel identifier must be a nonempty string")
        identifiers.append(identifier)
        left, right = panel.get("left"), panel.get("right")
        if (
            isinstance(left, bool)
            or isinstance(right, bool)
            or not isinstance(left, (int, float))
            or not isinstance(right, (int, float))
            or not math.isfinite(left)
            or not math.isfinite(right)
        ):
            raise ValueError("panel bounds must be finite physical numbers")
        a, b = Decimal(str(left)), Decimal(str(right))
        if a != previous or a >= b:
            raise ValueError("panels must tile the physical axis in ascending order")
        previous = b
        depth = _integer(panel.get("depth"), "panel depth")
        if (depth == 0 and panel.get("parent") is not None) or (
            depth > 0 and not isinstance(panel.get("parent"), str)
        ):
            raise ValueError("panel parent must agree with its depth")
    if previous != (Decimal(1) if axis == "t" else cutoff):
        raise ValueError("panels do not reach the physical endpoint")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("duplicate panel identifiers")
    return panels


def _panel_points(panels: Sequence[Mapping[str, Any]], denominator: int) -> list[str]:
    """Recreate the recorded float lattice without adaptive-grid imports."""
    fractions = range(1, denominator)
    if denominator == 8:
        fractions = range(1, denominator, 2)
    values = {
        float(
            (1 - j / denominator) * panel["left"] + (j / denominator) * panel["right"]
        )
        for panel in panels
        for j in fractions
    }
    return [repr(value) for value in sorted(values)]


def _training_points(panels: Sequence[Mapping[str, Any]]) -> list[str]:
    values = {
        float((1 - j / 4) * panel["left"] + (j / 4) * panel["right"])
        for panel in panels
        for j in range(5)
    }
    return [repr(value) for value in sorted(values)]


def _axis(raw: Any, cutoff: Decimal, axis: str) -> list[str]:
    values = _list(raw, f"{axis} explicit axis")
    numbers = [_decimal(value, f"{axis} explicit axis") for value in values]
    low, high = (Decimal(0), Decimal(1)) if axis == "t" else (-cutoff, cutoff)
    if len(numbers) < 2 or numbers[0] != low or numbers[-1] != high:
        raise ValueError("explicit axes must contain both physical endpoints")
    if any(a >= b for a, b in zip(numbers, numbers[1:], strict=False)):
        raise ValueError("explicit axes must be strictly increasing")
    if any(
        repr(float(value)) != text for value, text in zip(numbers, values, strict=True)
    ):
        raise ValueError("explicit axes must preserve repr-float provenance")
    return values


@_stable_decimal
def _mirrored_marks(
    panels: Sequence[Mapping[str, Any]],
    scores: Mapping[str, Any],
    axis: str,
    max_marked: int,
    relative_threshold: Decimal,
) -> list[str]:
    """Independently sort influential scores and include unique reflection mates."""
    identifiers = [panel["identifier"] for panel in panels]
    if set(scores) != set(identifiers):
        raise ValueError("panel scores must cover every current leaf panel")
    values = {
        identifier: _decimal(
            _object(scores[identifier], "panel score").get("score"),
            "panel score",
            nonnegative=True,
        )
        for identifier in identifiers
    }
    peak = max(values.values())
    with localcontext() as context:
        context.prec = len(peak.as_tuple().digits) + len(
            relative_threshold.as_tuple().digits
        )
        context.Emax, context.Emin = MAX_EMAX, MIN_EMIN
        threshold = peak * relative_threshold
    chosen = sorted(
        identifier
        for identifier in identifiers
        if values[identifier] > 0 and values[identifier] >= threshold
    )
    chosen.sort(key=values.__getitem__, reverse=True)
    chosen = chosen[:max_marked]
    marked = set(chosen)
    scale = (
        1.0 if axis == "t" else max(abs(panels[0]["left"]), abs(panels[-1]["right"]))
    )
    for identifier in chosen:
        panel = next(panel for panel in panels if panel["identifier"] == identifier)
        target = (
            (1 - panel["right"], 1 - panel["left"])
            if axis == "t"
            else (-panel["right"], -panel["left"])
        )
        matches = [
            other["identifier"]
            for other in panels
            if all(
                abs(other[bound] - endpoint)
                <= 4 * max(math.ulp(other[bound]), math.ulp(endpoint), math.ulp(scale))
                for bound, endpoint in zip(("left", "right"), target, strict=True)
            )
        ]
        if len(matches) != 1:
            raise ValueError("marked panel must have exactly one reflection mate")
        marked.add(matches[0])
    return [identifier for identifier in identifiers if identifier in marked]


def _refined_panels(
    panels: Sequence[Mapping[str, Any]], marks: Sequence[str]
) -> list[dict[str, Any]]:
    result = []
    for panel in panels:
        if panel["identifier"] not in marks:
            result.append(dict(panel))
            continue
        middle = (panel["left"] + panel["right"]) / 2
        for side, left, right in (
            ("L", panel["left"], middle),
            ("R", middle, panel["right"]),
        ):
            result.append(
                {
                    "identifier": panel["identifier"] + side,
                    "left": left,
                    "right": right,
                    "depth": panel["depth"] + 1,
                    "parent": panel["identifier"],
                }
            )
    return result


def check_final_cross(
    record: Mapping[str, Any],
    final_panels: Mapping[str, Any],
    validation: Mapping[str, Any],
    precision_bits: int = 192,
) -> dict[str, Any]:
    """Check all six seventh-fraction nodes of every final leaf panel."""
    cutoff = _decimal(record.get("cutoff"), "cutoff")
    panels = {
        axis: _panel_axis(final_panels.get(axis), cutoff, axis)
        for axis in ("t", "omega")
    }
    points = {axis: _panel_points(panels[axis], 7) for axis in panels}
    if (
        validation.get("t_nodes") != points["t"]
        or validation.get("omega_nodes") != points["omega"]
    ):
        raise ValueError(
            "final validation must contain exactly the untouched seventh-fraction axes"
        )
    validation_bits = _integer(
        validation.get("precision_bits"), "validation precision", 55
    )
    if (
        validation.get("rank") != record["rank"]
        or not record["precision_bits"] < validation_bits <= precision_bits
    ):
        raise ValueError("final validation rank or precision disagrees")
    for axis in panels:
        training = set(_training_points(panels[axis]))
        probes = set(_panel_points(panels[axis], 8))
        if set(points[axis]).intersection(training | probes):
            raise ValueError(
                "final validation coordinates overlap training or stopping probes"
            )
    report = evaluate_cross(record, points["t"], points["omega"], precision_bits)
    with mp.workprec(precision_bits):
        reference = mp.mpf(report["max_abs_residual"])
        stored = mp.mpf(str(validation.get("max_residual")))
        disagreement = abs(stored - reference)
        consistent = disagreement <= max(
            mp.mpf("1e-45"), abs(reference) * mp.mpf("1e-25")
        )
    passed = report["within_requested_tolerance_on_sample"]
    if type(validation.get("passed")) is not bool or validation["passed"] != passed:
        raise ValueError(
            "final validation pass flag disagrees with the independent cross"
        )
    report.update(
        {
            "cutoff": record["cutoff"],
            "tolerance": record["tolerance"],
            "stored_max_abs_residual": validation.get("max_residual"),
            "stored_absolute_disagreement": _encode(disagreement, precision_bits),
            "status": "passed" if consistent else "disagreement",
            "scope": "stored_selected_cross_on_all_final_leaf_seventh_fraction_points",
            "used_for_refinement_or_stopping": False,
        }
    )
    return report


def _provenance(
    payload: Mapping[str, Any], *, check_revision: bool, require_clean: bool = True
) -> dict[str, Any]:
    configuration = _object(payload.get("configuration"), "configuration")
    canonical = json.dumps(
        configuration, sort_keys=True, separators=(",", ":")
    ).encode()
    if payload.get("config_hash") != hashlib.sha256(canonical).hexdigest():
        raise ValueError("configuration hash disagrees with explicit controls")
    revision = payload.get("git_commit")
    if not isinstance(revision, str) or re.fullmatch("[0-9a-f]{40}", revision) is None:
        raise ValueError("git_commit must identify an exact source revision")
    if type(payload.get("git_dirty")) is not bool or (
        require_clean and payload["git_dirty"]
    ):
        raise ValueError("canonical adaptive census requires a clean source revision")
    sources = _object(payload.get("algorithm_source_sha256"), "algorithm_source_sha256")
    if not sources:
        raise ValueError("algorithm source provenance must be nonempty")
    checked = {}
    for path, expected in {**sources, "uv.lock": payload.get("lock_sha256")}.items():
        if (
            not isinstance(path, str)
            or Path(path).is_absolute()
            or ".." in Path(path).parts
        ):
            raise ValueError("source provenance paths must be repository-relative")
        if (
            not isinstance(expected, str)
            or re.fullmatch("[0-9a-f]{64}", expected) is None
        ):
            raise ValueError("source provenance requires SHA-256 values")
        if check_revision:
            completed = subprocess.run(
                ["git", "show", f"{revision}:{path}"],
                cwd=PROJECT_ROOT,
                capture_output=True,
                check=False,
            )
            if completed.returncode != 0:
                raise ValueError(f"recorded source revision blob unavailable: {path}")
            if hashlib.sha256(completed.stdout).hexdigest() != expected:
                raise ValueError(
                    f"recorded source hash disagrees with revision: {path}"
                )
        checked[path] = expected
    environment = _object(payload.get("environment"), "environment")
    if not all(
        isinstance(environment.get(name), str) and environment[name]
        for name in ("python", "numpy", "scipy", "mpmath")
    ):
        raise ValueError("numerical environment versions must be explicit")
    return {
        "git_commit": revision,
        "config_hash": payload["config_hash"],
        "source_hashes": checked,
        "revision_blobs_checked": check_revision,
        "recorded_environment": dict(environment),
    }


def _initial_panels(cutoff: float, axis: str) -> list[dict[str, Any]]:
    if axis == "t":
        left = [0.0, min(0.5, 1 / cutoff)]
        while left[-1] < 0.5:
            left.append(min(0.5, 2 * left[-1]))
        bounds = left + [1 - value for value in left[-2::-1]]
    else:
        positive = [0.0, 1.0]
        while positive[-1] < cutoff:
            positive.append(min(cutoff, 2 * positive[-1]))
        bounds = [-value for value in positive[:0:-1]] + positive
    prefix = "t" if axis == "t" else "omega"
    return [
        {
            "identifier": f"{prefix}{index:04d}",
            "left": left,
            "right": right,
            "depth": 0,
            "parent": None,
        }
        for index, (left, right) in enumerate(zip(bounds, bounds[1:], strict=False))
    ]


@_stable_decimal
def _classification(
    check: Mapping[str, Any], threshold: Decimal, *, joint: str
) -> bool:
    values = {
        name: _decimal(check.get(name), name, nonnegative=True)
        for name in ("training_max", "probe_max", "surplus_max")
    }
    flags = {
        "training_passed": values["training_max"] <= threshold,
        "probe_passed": values["probe_max"] <= threshold,
        "surplus_passed": values["surplus_max"] <= threshold / 4,
    }
    for name, value in {**flags, joint: all(flags.values())}.items():
        if type(check.get(name)) is not bool or check[name] != value:
            raise ValueError(f"recorded {name} disagrees with sampled errors")
    return all(flags.values())


@_stable_decimal
def _probe(
    check: Mapping[str, Any],
    panels: Mapping[str, Sequence[Mapping[str, Any]]],
    threshold: Decimal,
    check_bits: int,
) -> bool:
    if check.get("precision_bits") != check_bits:
        raise ValueError("probe precision disagrees with configuration")
    if _decimal(check.get("tolerance"), "probe tolerance") != threshold:
        raise ValueError("probe tolerance disagrees with requested scope")
    for axis in ("t", "omega"):
        if check.get(f"{axis}_nodes") != _panel_points(panels[axis], 8):
            raise ValueError("stopping probes must be exactly odd-eighth panel points")
        scores = _object(check.get(f"{axis}_scores"), "panel scores")
        if set(scores) != {panel["identifier"] for panel in panels[axis]}:
            raise ValueError("panel scores must cover every leaf")
        peaks = []
        for raw in scores.values():
            score = _object(raw, "panel score")
            peak = _decimal(score.get("probe_max"), "panel probe_max", nonnegative=True)
            surplus = _decimal(score.get("surplus"), "panel surplus", nonnegative=True)
            with localcontext() as context:
                context.prec = max(80, math.ceil(check_bits * math.log10(2)) + 20)
                context.Emax, context.Emin = MAX_EMAX, MIN_EMIN
                expected = max(peak / threshold, 4 * surplus / threshold)
            if (
                _decimal(score.get("score"), "panel score", nonnegative=True)
                != expected
            ):
                raise ValueError(
                    "panel score disagrees with recorded residual and surplus"
                )
            peaks.append(peak)
        if max(peaks) != _decimal(
            check.get("probe_max"), "probe_max", nonnegative=True
        ):
            raise ValueError("panel residual maxima disagree with probe maximum")
    surpluses = [
        _decimal(raw["surplus"], "panel surplus", nonnegative=True)
        for axis in ("t", "omega")
        for raw in check[f"{axis}_scores"].values()
    ]
    if max(surpluses) != _decimal(
        check.get("surplus_max"), "surplus_max", nonnegative=True
    ):
        raise ValueError("panel surpluses disagree with surplus maximum")
    return _classification(check, threshold, joint="passed")


def _selected_on_axis(
    record: Mapping[str, Any], axes: Mapping[str, Sequence[str]]
) -> None:
    with mp.workprec(record["precision_bits"] + 32):
        for key, axis in (("t_nodes", "t"), ("omega_nodes", "omega")):
            grid = [mp.mpf(value) for value in axes[axis]]
            for value in record[key]:
                number = mp.mpf(value)
                slack = mp.mpf(2) ** (5 - record["precision_bits"]) * max(
                    1, abs(number)
                )
                if not any(abs(number - coordinate) <= slack for coordinate in grid):
                    raise ValueError(
                        "selected coordinates absent from the explicit training axes"
                    )


def _geometry(record: Mapping[str, Any]) -> None:
    """Recompute physical distances/scales and nearest reflected selected nodes."""
    geometry = _list(record.get("pivot_geometry"), "pivot_geometry")
    if len(geometry) != record["rank"]:
        raise ValueError("pivot geometry must align with the selected prefix")
    with mp.workprec(record["precision_bits"] + 64):
        cutoff = mp.mpf(record["cutoff"])
        match = mp.mpf(record.get("reflection_match_tolerance", "1e-12"))
        nodes = [
            (mp.mpf(t), mp.mpf(w))
            for t, w in zip(record["t_nodes"], record["omega_nodes"], strict=True)
        ]
        slack = mp.mpf(2) ** (5 - record["precision_bits"])
        for (t, w), raw in zip(nodes, geometry, strict=True):
            item = _object(raw, "pivot geometry")
            for name, values in (
                ("t_boundary_distances", (t, 1 - t)),
                ("omega_boundary_distances", (w + cutoff, cutoff - w)),
            ):
                stored = _list(item.get(name), name)
                if len(stored) != 2 or any(
                    abs(mp.mpf(_decimal(value, name)) - expected)
                    > slack * max(1, abs(expected))
                    for value, expected in zip(stored, values, strict=True)
                ):
                    raise ValueError(
                        "pivot boundary distances disagree with physical coordinates"
                    )
            scale = 0
            while abs(w) > mp.mpf(2) ** scale:
                scale += 1
            if item.get("dyadic_scale") != scale or item.get("dyadic_panel") != (
                "central" if scale == 0 else "positive" if w > 0 else "negative"
            ):
                raise ValueError("pivot dyadic scale disagrees with frequency")
            distances = [max(abs(s - (1 - t)), abs(v + w) / cutoff) for s, v in nodes]
            partner = _integer(
                item.get("reflection_nearest_index"), "reflection partner"
            )
            if partner >= len(nodes) or distances[partner] > min(distances) + slack:
                raise ValueError("reflection partner is not nearest in this prefix")
            for name, value in (
                ("reflection_t", 1 - t),
                ("reflection_omega", -w),
                ("reflection_distance", distances[partner]),
            ):
                if abs(mp.mpf(_decimal(item.get(name), name)) - value) > slack * max(
                    1, abs(value)
                ):
                    raise ValueError("reflected coordinate or distance disagrees")
            if type(item.get("reflection_matched")) is not bool or item[
                "reflection_matched"
            ] != (distances[partner] <= match):
                raise ValueError(
                    "reflection match flag disagrees with recorded tolerance"
                )


@_stable_decimal
def validate_payload(
    payload: Mapping[str, Any],
    *,
    require_full: bool = True,
    check_revision: bool = True,
) -> dict[str, Any]:
    """Check controls, explicit axes, genuine refinement and all recorded decisions.

    Unsuccessful trajectories are valid evidence when their failure is retained.
    ``require_full=False`` labels a pilot; it cannot establish seven-cutoff scope.
    """
    if (
        payload.get("schema_version") != 3
        or payload.get("kind") != "residual_adaptive_sampled_grid"
    ):
        raise ValueError(
            "adaptive census requires schema-3 residual-adaptive sampled-grid data"
        )
    config = _object(payload.get("configuration"), "configuration")
    tolerances = [
        Decimal(str(value)) for value in _list(config.get("tolerances"), "tolerances")
    ]
    cutoffs = [Decimal(str(value)) for value in _list(config.get("cutoffs"), "cutoffs")]
    if not tolerances or any(
        not value.is_finite() or not 0 < value < 1 for value in tolerances
    ):
        raise ValueError("tolerances must be finite and in (0,1)")
    if tolerances != sorted(set(tolerances), reverse=True):
        raise ValueError("tolerances must be distinct and decreasing")
    if (
        not cutoffs
        or len(set(cutoffs)) != len(cutoffs)
        or any(not value.is_finite() or value < 1 for value in cutoffs)
    ):
        raise ValueError("cutoffs must be distinct, finite and at least one")
    if require_full and (
        set(cutoffs) != set(CANONICAL_CUTOFFS)
        or set(tolerances) != set(CANONICAL_TOLERANCES)
    ):
        raise ValueError("full artifact requires all seven cutoffs and four tolerances")
    controls = {
        name: _integer(config.get(name), name, minimum)
        for name, minimum in (
            ("precision_bits", 54),
            ("check_precision_bits", 55),
            ("max_rank", 0),
            ("max_stages", 1),
            ("max_nodes", 2),
            ("max_marked", 1),
        )
    }
    if controls["check_precision_bits"] <= controls["precision_bits"]:
        raise ValueError("check precision must exceed elimination precision")
    if require_full and (
        controls["precision_bits"] != 128 or controls["check_precision_bits"] != 192
    ):
        raise ValueError("full artifact requires MP128 elimination and MP192 checks")
    margin = Decimal(str(config.get("generation_margin")))
    if not margin.is_finite() or not 0 < margin < 1:
        raise ValueError("generation margin must be finite and in (0,1)")
    provenance = _provenance(
        payload, check_revision=check_revision, require_clean=require_full
    )
    required_sources = {
        "kernelgecp/adaptive_grid.py",
        "kernelgecp/adaptive_grids.py",
        "kernelgecp/gecp.py",
        "kernelgecp/census.py",
        "kernelgecp/kernels.py",
        "kernelgecp/types.py",
        "experiments/adaptive_census.py",
    }
    if not required_sources.issubset(provenance["source_hashes"]):
        raise ValueError("source provenance omits an algorithm input")
    cases = [_object(case, "case") for case in _list(payload.get("cases"), "cases")]
    if len(cases) != len(cutoffs) or {
        _decimal(case.get("cutoff"), "cutoff") for case in cases
    } != set(cutoffs):
        raise ValueError("case coverage disagrees with configured cutoffs")
    summaries = []
    for case in cases:
        cutoff = _decimal(case["cutoff"], "cutoff")
        stages = [
            _object(stage, "stage") for stage in _list(case.get("stages"), "stages")
        ]
        if not stages or len(stages) > controls["max_stages"]:
            raise ValueError("stage count exceeds recorded budget")
        previous_axes: dict[str, list[str]] | None = None
        previous_panels: dict[str, Sequence[Mapping[str, Any]]] | None = None
        previous_marks: dict[str, list[str]] | None = None
        used = {"t": set(), "omega": set()}
        consecutive = 0
        refinements = 0
        for index, stage in enumerate(stages):
            if stage.get("index") != index:
                raise ValueError("stage indices must start at zero and be consecutive")
            budgets = _object(stage.get("budgets"), "stage budgets")
            expected_budgets = {
                **controls,
                "tolerances": config["tolerances"],
                "generation_margin": config["generation_margin"],
            }
            if dict(budgets) != expected_budgets:
                raise ValueError(
                    "stage budget snapshot disagrees with frozen configuration"
                )
            if _decimal(
                stage.get("generation_tolerance"), "generation tolerance"
            ) != Decimal(str(float(min(tolerances) * margin))):
                raise ValueError("generation tolerance disagrees with frozen margin")
            panels = {
                axis: _panel_axis(stage.get(f"{axis}_panels"), cutoff, axis)
                for axis in ("t", "omega")
            }
            axes = {
                axis: _axis(stage.get(f"{axis}_grid"), cutoff, axis) for axis in panels
            }
            if stage.get("grid_shape") != [len(axes["t"]), len(axes["omega"])]:
                raise ValueError("stage shape disagrees with explicit axes")
            if (
                max(map(len, axes.values())) > controls["max_nodes"]
                or stage.get("nested") is not True
            ):
                raise ValueError("stage violates node budget or nesting metadata")
            for axis in panels:
                expected_panels = (
                    _initial_panels(float(cutoff), axis)
                    if previous_panels is None
                    else _refined_panels(previous_panels[axis], previous_marks[axis])
                )
                if panels[axis] != expected_panels:
                    raise ValueError(
                        "leaf panels disagree with actual recorded refinement"
                    )
                expected_axis = set(_training_points(panels[axis]))
                if previous_axes is not None:
                    expected_axis.update(previous_axes[axis])
                    if not set(previous_axes[axis]).issubset(axes[axis]):
                        raise ValueError("training axes are not nested")
                if set(axes[axis]) != expected_axis:
                    raise ValueError(
                        "explicit axis disagrees with retained quarter-lattice "
                        "coordinates"
                    )
                used[axis].update(axes[axis])
            gecp = _object(stage.get("gecp"), "stage gecp")
            rank = _integer(gecp.get("rank"), "stage rank")
            if rank > controls["max_rank"] or rank > min(map(len, axes.values())):
                raise ValueError("stage rank exceeds its budget or axis shape")
            high = _object(gecp.get("high_precision"), "stage high_precision")
            stage_metadata = _object(gecp.get("metadata"), "stage GECP metadata")
            if (
                stage_metadata.get("precision_bits") != controls["precision_bits"]
                or stage_metadata.get("convergence_scope") != "sampled_grid"
            ):
                raise ValueError(
                    "stage metadata disagrees with sampled source precision/scope"
                )
            trajectory = {
                **high,
                "core_sigma_min": high.get("core_sigma_min_history"),
                "schema_version": 3,
                "algorithm": "residual_adaptive_grid",
                "precision_bits": controls["precision_bits"],
                "cutoff": case["cutoff"],
                "tolerance": str(min(tolerances) * margin),
                "rank": rank,
                "residual": high.get("residual_history", [None])[-1],
                "tie_counts": gecp.get("tie_counts"),
                "converged": gecp.get("converged"),
                "stop_reason": gecp.get("stop_reason"),
            }
            _validate_record(trajectory)
            _selected_on_axis(trajectory, axes)
            if gecp["converged"] != (
                Decimal(trajectory["residual"]) <= min(tolerances) * margin
            ):
                raise ValueError(
                    "GECP convergence disagrees with its generation margin"
                )
            probe = _object(stage.get("probe"), "stage probe")
            if (
                probe.get("rank") != rank
                or probe.get("training_max") != trajectory["residual"]
            ):
                raise ValueError("stage probe disagrees with actual trajectory prefix")
            jointly = _probe(
                probe, panels, min(tolerances), controls["check_precision_bits"]
            )
            for axis in panels:
                probe_nodes = probe[f"{axis}_nodes"]
                if set(probe_nodes).intersection(axes[axis]):
                    raise ValueError(
                        "stopping probes overlap explicit training coordinates"
                    )
                used[axis].update(probe_nodes)
            passed = gecp["converged"] and jointly
            consecutive = consecutive + 1 if passed else 0
            if (
                stage.get("passed") != passed
                or stage.get("consecutive_passes") != consecutive
            ):
                raise ValueError("stage pass count disagrees with joint sampled checks")
            marks = {
                axis: _list(stage.get(f"marked_{axis}"), "marked panels")
                for axis in panels
            }
            expected_marks = {
                axis: _mirrored_marks(
                    panels[axis],
                    probe[f"{axis}_scores"],
                    axis,
                    controls["max_marked"],
                    Decimal("0.5"),
                )
                for axis in panels
            }
            applied = stage.get("refinement_applied")
            if type(applied) is not bool:
                raise ValueError("refinement_applied must be boolean")
            if index < len(stages) - 1:
                if not applied or not any(marks.values()) or marks != expected_marks:
                    raise ValueError(
                        "actual refinement must follow residual scores and "
                        "reflection marks"
                    )
                refinements += 1
                if consecutive >= 2:
                    raise ValueError(
                        "controller continued after two successive passing stages"
                    )
            elif applied:
                raise ValueError("last stage cannot apply an unrecorded refinement")
            elif any(marks.values()) and marks != expected_marks:
                raise ValueError(
                    "rejected refinement marks disagree with residual scores"
                )
            if applied:
                proposed = {
                    axis: set(axes[axis])
                    | set(_training_points(_refined_panels(panels[axis], marks[axis])))
                    for axis in panels
                }
                if not any(len(proposed[axis]) > len(axes[axis]) for axis in panels):
                    raise ValueError("refinement did not add real training coordinates")
                if stage.get("proposed_grid_shape") != [
                    len(proposed["t"]),
                    len(proposed["omega"]),
                ]:
                    raise ValueError(
                        "proposed grid shape disagrees with actual panel splits"
                    )
            previous_axes, previous_panels, previous_marks = axes, panels, marks
        final_axes = _object(case.get("final_axes"), "final_axes")
        final_panels = _object(case.get("final_panels"), "final_panels")
        if final_axes != previous_axes or final_panels != previous_panels:
            raise ValueError("final geometry disagrees with the last actual stage")
        controller_passed = len(stages) >= 2 and consecutive >= 2 and refinements >= 1
        records = [
            _object(record, "record")
            for record in _list(case.get("records"), "records")
        ]
        prefixes = [
            _object(prefix, "prefix")
            for prefix in _list(case.get("prefixes"), "prefixes")
        ]
        if len(records) != len(tolerances) or len(prefixes) != len(tolerances):
            raise ValueError("each case must retain all tolerance records and prefixes")
        for record, prefix, tolerance in zip(
            records, prefixes, tolerances, strict=True
        ):
            _validate_record(record)
            if record["precision_bits"] != controls["precision_bits"]:
                raise ValueError(
                    "record precision disagrees with frozen source precision"
                )
            if (
                _decimal(record["cutoff"], "cutoff") != cutoff
                or _decimal(record["tolerance"], "tolerance") != tolerance
                or _decimal(prefix.get("tolerance"), "prefix tolerance") != tolerance
            ):
                raise ValueError("cutoff/tolerance record alignment disagrees")
            if (
                record.get("git_commit") != payload["git_commit"]
                or record.get("config_hash") != payload["config_hash"]
                or record.get("git_dirty") != payload["git_dirty"]
            ):
                raise ValueError("record provenance disagrees with artifact provenance")
            if record.get("grid_shape") != [
                len(final_axes["t"]),
                len(final_axes["omega"]),
            ]:
                raise ValueError("record shape must use explicit final axes")
            metadata = _object(record.get("metadata"), "record metadata")
            if (
                metadata.get("convergence_scope")
                != "sampled_training_probe_and_independent_validation"
                or metadata.get("certification_status") != "not_certified"
            ):
                raise ValueError(
                    "record must label the finite sampled convergence scope"
                )
            _selected_on_axis(record, final_axes)
            _geometry(record)
            rank = record["rank"]
            if prefix.get("chosen_rank") != rank or rank > trajectory["rank"]:
                raise ValueError("chosen prefix rank disagrees with actual trajectory")
            for name in (*HISTORIES, "residual_history"):
                if (
                    record[name]
                    != trajectory[name][: rank + (name == "residual_history")]
                ):
                    raise ValueError(
                        "record histories disagree with final trajectory prefix"
                    )
            training_rank = next(
                (
                    i
                    for i, value in enumerate(trajectory["residual_history"])
                    if Decimal(value) <= tolerance
                ),
                None,
            )
            if prefix.get("training_rank") != training_rank:
                raise ValueError(
                    "first training-passing rank disagrees with full history"
                )
            start = training_rank if training_rank is not None else trajectory["rank"]
            if prefix.get("tested_ranks") != list(range(start, rank + 1)):
                raise ValueError("tested prefix ranks must retain every attempted rank")
            if prefix.get("training_max") != record["residual"]:
                raise ValueError(
                    "prefix training maximum disagrees with selected history"
                )
            joint = _classification(prefix, tolerance, joint="joint_passed")
            if "probe" in prefix:
                _probe(
                    _object(prefix["probe"], "prefix probe"),
                    previous_panels,
                    tolerance,
                    controls["check_precision_bits"],
                )
            validation = _object(prefix.get("validation"), "final validation")
            if (
                validation.get("rank") != rank
                or validation.get("precision_bits") != controls["check_precision_bits"]
            ):
                raise ValueError("final validation rank or precision disagrees")
            for axis in ("t", "omega"):
                expected = _panel_points(previous_panels[axis], 7)
                if validation.get(f"{axis}_nodes") != expected or used[
                    axis
                ].intersection(expected):
                    raise ValueError(
                        "final validation must be untouched seventh-fraction "
                        "samples, disjoint from every stage"
                    )
            validation_passed = (
                _decimal(
                    validation.get("max_residual"),
                    "final validation max",
                    nonnegative=True,
                )
                <= tolerance
            )
            if (
                type(validation.get("passed")) is not bool
                or validation["passed"] != validation_passed
            ):
                raise ValueError("final validation pass flag disagrees with error")
            fit_passed = joint and validation_passed
            success = controller_passed and fit_passed
            if (
                prefix.get("fit_passed") != fit_passed
                or prefix.get("controller_passed") != controller_passed
            ):
                raise ValueError(
                    "prefix fit/controller flags disagree with recorded evidence"
                )
            if prefix.get("converged") != success or record["converged"] != success:
                raise ValueError(
                    "prefix convergence disagrees with joint and final "
                    "validation checks"
                )
        success = controller_passed and all(prefix["converged"] for prefix in prefixes)
        if type(case.get("converged")) is not bool or case["converged"] != success:
            raise ValueError(
                "case convergence requires two actual passing stages and "
                "all final prefixes"
            )
        if success and case.get("stop_reason") != "sampled_tolerances":
            raise ValueError("successful case must record sampled_tolerances stop")
        if not success and case.get("stop_reason") == "sampled_tolerances":
            raise ValueError("unsuccessful case cannot claim sampled convergence")
        stop = stages[-1].get("stop_reason")
        if stop == "max_stages" and len(stages) != controls["max_stages"]:
            raise ValueError("max_stages stop does not exhaust its recorded budget")
        if stop == "max_nodes":
            proposed = stages[-1].get("proposed_grid_shape")
            if (
                not isinstance(proposed, list)
                or len(proposed) != 2
                or max(proposed) <= controls["max_nodes"]
            ):
                raise ValueError("max_nodes stop does not exhaust its recorded budget")
        if stop == "strictest_sampled_stages_passed" and not controller_passed:
            raise ValueError("controller stop lacks two genuine passing stages")
        summaries.append(
            {
                "cutoff": case["cutoff"],
                "stage_count": len(stages),
                "real_refinement_count": refinements,
                "converged": success,
                "stop_reason": case.get("stop_reason"),
            }
        )
    return {
        "scope": "finite_training_probe_and_independent_validation_sets_only",
        "full_requested_coverage": require_full,
        "case_count": len(cases),
        "record_count": sum(len(case["records"]) for case in cases),
        "provenance": provenance,
        "cases": summaries,
        "status": "passed",
        "all_cases_converged": all(case["converged"] for case in cases),
        "certification_status": "not_certified",
    }


def check_adaptive_census(
    input_path: Path,
    *,
    precision_bits: int = 192,
    require_full: bool = True,
    check_revision: bool = True,
) -> dict[str, Any]:
    """Return deterministic independent checks, retaining numerical disagreements."""
    data = input_path.read_bytes()
    payload = _object(json.loads(data), "adaptive census")
    structure = validate_payload(
        payload, require_full=require_full, check_revision=check_revision
    )
    core_checks, cross_checks = [], []
    for case in payload["cases"]:
        deepest = max(case["records"], key=lambda record: record["rank"])
        core_checks.append(check_selected_cores(deepest, precision_bits))
        for record, prefix in zip(case["records"], case["prefixes"], strict=True):
            cross_checks.append(
                check_final_cross(
                    record, case["final_panels"], prefix["validation"], precision_bits
                )
            )
    passed = all(check["status"] == "passed" for check in (*core_checks, *cross_checks))
    return {
        "schema_version": 1,
        "kind": "independent_residual_adaptive_sampled_grid_checks",
        "input_sha256": hashlib.sha256(data).hexdigest(),
        "checker_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "precision_bits": precision_bits,
        "structure": structure,
        "selected_cores": core_checks,
        "final_crosses": cross_checks,
        "status": "passed" if passed else "disagreement",
        "certification_status": "not_certified",
        "scope": (
            "original_selected_cores_and_final_untouched_seventh_fraction_cross_samples"
        ),
        "not_claimed": [
            "continuous_supremum_bound",
            "interval_certificate",
            "exact_inexact_pivot_rate",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path("experiments/data/adaptive_census.json")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/data/adaptive_census_checks.json"),
    )
    parser.add_argument("--precision-bits", type=int, default=192)
    parser.add_argument(
        "--allow-partial", action="store_true", help="labelled pilot data only"
    )
    args = parser.parse_args(argv)
    try:
        report = check_adaptive_census(
            args.input,
            precision_bits=args.precision_bits,
            require_full=not args.allow_partial,
        )
    except (ValueError, ArithmeticError) as error:
        parser.error(str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
