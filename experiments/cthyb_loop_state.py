"""Content-bound exact state for newly generated CT-HYB feedback loops.

These hashes detect accidental changes; they are not signatures. Old loops
without state are rejected rather than interpreted as compatible checkpoints.
"""

from __future__ import annotations

import base64
import dataclasses
from collections.abc import Callable
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np


def _encode(array: np.ndarray) -> str:
    values = np.asarray(array)
    if (
        values.dtype != np.dtype("complex128")
        or values.ndim != 1
        or not np.isfinite(values).all()
    ):
        raise ValueError("invalid loop state bath")
    return base64.b64encode(
        np.ascontiguousarray(values, dtype="<c16").tobytes()
    ).decode("ascii")


def _decode(value: str, size: int) -> np.ndarray:
    try:
        raw = base64.b64decode(value, validate=True)
    except (ValueError, TypeError) as error:
        raise ValueError("invalid loop state bath encoding") from error
    if len(raw) != size * 16:
        raise ValueError("loop state bath size mismatch")
    values = np.frombuffer(raw, dtype="<c16").astype(np.complex128, copy=True)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite loop state bath")
    return values


def _validate_final_mean(final: np.ndarray, baths: list[np.ndarray]) -> None:
    """Check a floating mean against exact sums with an IEEE binary64 bound.

    For n terms, k=n+4 bounds the n-1 reduction additions and the scalar
    division implementation (including a rounded reciprocal and component
    arithmetic). For each real/imag component, gamma_k=k*u/(1-k*u),
    u=2^-53; one minimum subnormal per operation bounds gradual-underflow
    error. Exact rational arithmetic evaluates the bound without rounding it
    downward. This assumes finite IEEE round-to-nearest binary64 operations.
    """
    n = len(baths)
    k = n + 4
    u = Fraction(1, 2**53)
    eta = Fraction(1, 2**1074)
    if not n or k * u >= 1:
        raise ValueError("invalid final bath mean count")
    gamma = k * u / (1 - k * u)
    subnormal = k * eta / (1 - k * u)
    for component in ("real", "imag"):
        result = getattr(final, component)
        values = [getattr(bath, component) for bath in baths]
        for index in range(len(final)):
            terms = [Fraction.from_float(float(a[index])) for a in values]
            exact = sum(terms, Fraction()) / n
            magnitude = sum((abs(x) for x in terms), Fraction()) / n
            error = abs(Fraction.from_float(float(result[index])) - exact)
            if error > gamma * magnitude + subnormal:
                raise ValueError("final bath violates recorded mean rounding bound")


def restore(
    folder: Path,
    plan: Any,
    cycles: int,
    seed: int,
    identity: dict[str, Any],
    *,
    maximum: int,
    pilot: bool,
    compressed: bool,
    data: dict[str, Any],
    read: Callable,
    sha: Callable,
    load: Callable,
) -> dict[str, Any] | None:
    """Verify all state links and recorded measurement inputs, without replay."""
    head_path = folder / "state-head.json"
    states = sorted(folder.glob("state-[0-9][0-9].json"))
    archives = sorted(folder.glob("iteration-*.*"))
    if not head_path.exists():
        if states or archives or (folder / "loop.json").exists():
            raise ValueError("legacy or orphan loop checkpoint lacks state head")
        return None
    if head_path.is_symlink():
        raise ValueError("loop state head cannot be a symlink")
    head = read(head_path)
    count = head.get("next_index")
    if (
        type(count) is not int
        or not 1 <= count <= maximum
        or head.get("identity") != identity
        or head.get("model_config") != dataclasses.asdict(plan.config)
    ):
        raise ValueError("invalid loop state head")
    expected_states = [folder / f"state-{i:02d}.json" for i in range(count)]
    expected_archives = sorted(
        [
            folder / f"iteration-{i:02d}{suffix}"
            for i in range(count)
            for suffix in (".json", ".npz")
        ]
    )
    if states != expected_states or archives != expected_archives:
        raise ValueError("orphan or missing loop state or measurement")
    if any(p.is_symlink() for p in states + archives):
        raise ValueError("loop state paths cannot be symlinks")
    initial_bath = _decode(head["initial_bath"], 2 * plan.config.n_iw)
    previous_sha = None
    previous_bath = initial_bath
    previous_history: list[dict[str, Any]] = []
    baths = []
    for i, path in enumerate(states):
        state = read(path)
        stem = folder / f"iteration-{i:02d}"
        if (
            state.get("schema") != "kernelgecp.cthyb-loop-state.v1"
            or state.get("identity") != identity
            or state.get("initial_bath") != head.get("initial_bath")
            or state.get("model_config") != head.get("model_config")
            or state.get("previous_state_sha256") != previous_sha
            or state.get("next_index") != i + 1
            or state.get("pilot") != pilot
            or state.get("compressed") != compressed
            or state.get("measurement_json_sha256") != sha(stem.with_suffix(".json"))
            or state.get("measurement_npz_sha256") != sha(stem.with_suffix(".npz"))
        ):
            raise ValueError("loop state chain or measurement hash mismatch")
        measurement = load(stem, identity)
        if (
            measurement.config != plan.config
            or measurement.cycles != cycles
            or measurement.seed != seed + 2 * i + 1
        ):
            raise ValueError("loop state measurement config, seed or cycles mismatch")
        if previous_bath is not None and not np.array_equal(
            measurement.bath_iw, previous_bath
        ):
            raise ValueError("loop state measurement bath mismatch")
        baths.append(measurement.bath_iw)
        history = state.get("history")
        if (
            not isinstance(history, list)
            or len(history) != i + 1
            or history[:-1] != previous_history
            or history[-1].get("iteration") != i + 1
        ):
            raise ValueError("loop state history mismatch")
        consecutive = 0
        for h in history:
            consecutive = (
                consecutive + 1 if h["unmixed_defect"] < data["loop_tolerance"] else 0
            )
        if state.get("consecutive") != consecutive:
            raise ValueError("loop state controller mismatch")
        if type(state.get("completed")) is not bool:
            raise ValueError("invalid loop completion flag")
        if not state["completed"] and (
            state.get("report") is not None or state.get("final_bath") is not None
        ):
            raise ValueError("partial loop contains a final result")
        should_stop = (
            not pilot
            and i + 1 >= data["min_iterations"]
            and consecutive >= data["consecutive_passes"]
        ) or i + 1 == maximum
        if bool(state.get("completed")) != should_stop or (
            should_stop and i != count - 1
        ):
            raise ValueError("loop state completion mismatch")
        previous_bath = _decode(state["next_bath"], 2 * plan.config.n_iw)
        previous_history = history
        previous_sha = sha(path)
    if head.get("state_sha256") != previous_sha:
        raise ValueError("loop state head hash mismatch")
    final_bath = None
    report = None
    if state["completed"]:
        final_bath = _decode(state["final_bath"], 2 * plan.config.n_iw)
        _validate_final_mean(final_bath, baths[-data["final_bath_average"] :])
        report = state.get("report")
        converged = (
            not pilot
            and consecutive >= data["consecutive_passes"]
            and count >= data["min_iterations"]
        )
        expected_report = {
            "history": history,
            "converged": converged,
            "stop_reason": "calibration_pilot_budget"
            if pilot
            else "stochastic_controller"
            if converged
            else "max_iterations_inconclusive",
            "final_bath_rule": "mean_of_last_declared_number_of_measured_input_baths",
            "averaged_baths": min(count, data["final_bath_average"]),
        }
        if report != expected_report:
            raise ValueError("loop state final report mismatch")
    return {
        "next_index": count,
        "initial_bath": initial_bath,
        "next_bath": previous_bath,
        "history": history,
        "consecutive": consecutive,
        "baths": baths,
        "previous_state_sha256": previous_sha,
        "final_bath": final_bath,
        "report": report,
    }


def commit(
    folder: Path,
    *,
    index: int,
    identity: dict[str, Any],
    previous_state_sha256: str | None,
    next_bath: np.ndarray,
    initial_bath: np.ndarray,
    model_config: dict[str, Any],
    history: list[dict[str, Any]],
    consecutive: int,
    pilot: bool,
    compressed: bool,
    report: dict[str, Any] | None,
    final_bath: np.ndarray | None,
    sha: Callable,
    write: Callable,
) -> None:
    """Write step first, then atomically advance the external state head."""
    path = folder / f"state-{index:02d}.json"
    stem = folder / f"iteration-{index:02d}"
    if path.exists():
        raise ValueError("refusing to overwrite loop state")
    state = {
        "schema": "kernelgecp.cthyb-loop-state.v1",
        "identity": identity,
        "initial_bath": _encode(initial_bath),
        "model_config": model_config,
        "previous_state_sha256": previous_state_sha256,
        "next_index": index + 1,
        "measurement_json_sha256": sha(stem.with_suffix(".json")),
        "measurement_npz_sha256": sha(stem.with_suffix(".npz")),
        "next_bath": _encode(next_bath),
        "history": history,
        "consecutive": consecutive,
        "pilot": pilot,
        "compressed": compressed,
        "completed": report is not None,
        "report": report,
        "final_bath": None if final_bath is None else _encode(final_bath),
    }
    write(path, state)
    write(
        folder / "state-head.json",
        {
            "identity": identity,
            "initial_bath": _encode(initial_bath),
            "model_config": model_config,
            "next_index": index + 1,
            "state_sha256": sha(path),
        },
    )
