"""A native exit-zero message cannot stand in for Green convergence evidence."""

from __future__ import annotations

import pytest

from scripts.providers.green_smoke import (
    ENERGY_THRESHOLD,
    GREEN_CHANGE_THRESHOLD,
    MAX_ITERATIONS,
    validate_convergence,
)


def test_native_marker_and_independent_stored_deltas_are_both_required():
    validate_convergence(
        "============== Simulation Converged ==============\nCompleted.",
        iteration=7,
        energy_change_sum=1e-12,
        green_change=1e-8,
    )


@pytest.mark.parametrize(
    "log",
    [
        "",
        "Completed.",
        "====== Reached Maximum number of iterations ======\nCompleted.",
        "Simulation Converged\nReached Maximum number of iterations\nCompleted.",
    ],
)
def test_exit_zero_completed_and_max_iteration_logs_do_not_pass(log):
    with pytest.raises(RuntimeError, match="native self-consistency convergence"):
        validate_convergence(
            log, iteration=7, energy_change_sum=1e-12, green_change=1e-8
        )


@pytest.mark.parametrize("iteration", [True, False, 0, 1, -1, 3.0, MAX_ITERATIONS + 1])
def test_incomplete_or_invalid_iteration_pair_cannot_pass(iteration):
    with pytest.raises(RuntimeError, match="converged pair"):
        validate_convergence(
            "Simulation Converged",
            iteration=iteration,
            energy_change_sum=1e-12,
            green_change=1e-8,
        )


@pytest.mark.parametrize("field", ["energy_change_sum", "green_change"])
@pytest.mark.parametrize(
    "value", [True, False, -1e-12, float("nan"), float("inf"), -float("inf")]
)
def test_bad_stored_delta_cannot_be_promoted_to_convergence(field, value):
    deltas = {"energy_change_sum": 1e-12, "green_change": 1e-8, field: value}
    with pytest.raises(RuntimeError):
        validate_convergence("Simulation Converged", iteration=7, **deltas)


@pytest.mark.parametrize(
    "energy_change_sum,green_change",
    [(ENERGY_THRESHOLD, 1e-8), (1e-12, GREEN_CHANGE_THRESHOLD)],
)
def test_native_marker_does_not_override_failed_independent_delta(
    energy_change_sum, green_change
):
    with pytest.raises(RuntimeError, match="threshold|too large"):
        validate_convergence(
            "Simulation Converged",
            iteration=7,
            energy_change_sum=energy_change_sum,
            green_change=green_change,
        )
