"""Controller budgets, independent crosses, deterministic choices and isolation."""

import importlib
import json
from collections import Counter
from dataclasses import replace
from decimal import Decimal, Inexact, localcontext

import mpmath as mp
import numpy as np
import pytest

from kernelgecp.adaptive_grid import (
    AdaptiveGridCase,
    AdaptiveGridConfig,
    _classify_probe,
    _lu_solve,
    _residual_matrix,
    run_adaptive_grid,
)

MODULE = importlib.import_module("kernelgecp.adaptive_grid")


def _config(**changes: object) -> AdaptiveGridConfig:
    return replace(
        AdaptiveGridConfig(
            tolerances=(1e-4, 1e-6),
            precision_bits=80,
            check_precision_bits=128,
            max_stages=3,
        ),
        **changes,
    )


def _canonical(case: AdaptiveGridCase) -> str:
    return json.dumps(
        {
            "stages": case.stages,
            "prefixes": case.prefixes,
            "converged": case.converged,
            "stop_reason": case.stop_reason,
        },
        sort_keys=True,
        allow_nan=False,
    )


@pytest.fixture(scope="module")
def case() -> AdaptiveGridCase:
    return run_adaptive_grid(1.0, config=_config())


@pytest.mark.parametrize(
    "changes",
    [
        {"precision_bits": 53},
        {"precision_bits": True},
        {"check_precision_bits": 80},
        {"max_rank": -1},
        {"max_stages": 0},
        {"max_nodes": 1},
        {"max_marked": False},
        {"generation_margin": 1.0},
        {"generation_margin": "bad"},
        {"generation_margin": float("nan")},
        {"tolerances": ()},
        {"tolerances": (1e-6, 1e-4)},
        {"tolerances": (1e-6, 1e-6)},
        {"tolerances": (float("nan"),)},
        {"tolerances": (True,)},
        {"tolerances": ("bad",)},
    ],
)
def test_invalid_controls_are_rejected(changes: dict) -> None:
    with pytest.raises(ValueError):
        _config(**changes)


def test_two_actual_successive_stages_and_complete_sampled_metadata(
    case: AdaptiveGridCase,
) -> None:
    assert case.converged and case.stop_reason == "sampled_tolerances"
    assert len(case.stages) >= 2
    assert case.stages[-1]["consecutive_passes"] >= 2
    previous_t: set[float] = set()
    previous_w: set[float] = set()
    for index, stage in enumerate(case.stages):
        assert stage["index"] == index and stage["nested"]
        times = set(map(float, stage["t_grid"]))
        frequencies = set(map(float, stage["omega_grid"]))
        assert previous_t.issubset(times) and previous_w.issubset(frequencies)
        if index:
            assert times != previous_t or frequencies != previous_w
        previous_t, previous_w = times, frequencies
        assert stage["coordinate_precision_bits"] == 53
        assert stage["coordinate_semantics"] == "repr_float64_decimal"
        assert stage["certification_status"] == "not_certified"
        rank = stage["gecp"]["rank"]
        high = stage["gecp"]["high_precision"]
        assert len(high["residual_history"]) == rank + 1
        assert len(high["core_sigma_min_history"]) == rank
        assert stage["probe"]["precision_bits"] == 128
        assert stage["probe"]["factorizations"] == 1
        assert all(
            Decimal(values["score"]).is_finite()
            for axis in ("t", "omega")
            for values in stage["probe"][f"{axis}_scores"].values()
        )
    for prefix in case.prefixes:
        assert prefix["controller_passed"] and prefix["fit_passed"]
        assert prefix["converged"] and prefix["joint_passed"]
        assert (
            prefix["training_rank"] <= prefix["chosen_rank"] <= case.final_result.rank
        )
        assert Decimal(prefix["validation"]["max_residual"]) <= Decimal(
            prefix["tolerance"]
        )
        assert "after_refinement_stops" in prefix["validation"]["scope"]


def test_repeat_outputs_are_identical(case: AdaptiveGridCase) -> None:
    repeated = run_adaptive_grid(1.0, config=_config())
    assert _canonical(case) == _canonical(repeated)


def test_actual_selective_refinement_does_not_relabel_a_uniform_order_change() -> None:
    result = run_adaptive_grid(10.0, config=_config(tolerances=(1e-5,), max_marked=1))
    first, second = result.stages[:2]
    assert 0 < len(first["marked_t"]) < len(first["t_panels"])
    assert 0 < len(first["marked_omega"]) < len(first["omega_panels"])
    assert first["refinement_applied"]
    assert second["grid_shape"][0] > first["grid_shape"][0]
    assert second["grid_shape"][1] > first["grid_shape"][1]


@pytest.mark.parametrize(
    "changes,reason",
    [
        ({"max_stages": 1}, "max_stages"),
        ({"max_nodes": 9}, "max_nodes"),
        ({"max_rank": 1}, "max_rank"),
    ],
)
def test_budget_exhaustion_never_promotes_a_local_fit_to_protocol_success(
    changes: dict,
    reason: str,
) -> None:
    result = run_adaptive_grid(1.0, config=_config(**changes))
    assert not result.converged and result.stop_reason == reason
    assert result.stages[-1]["stop_reason"] == reason
    assert all(not prefix["controller_passed"] for prefix in result.prefixes)
    assert all(not prefix["converged"] for prefix in result.prefixes)
    assert all(prefix["stop_reason"] == reason for prefix in result.prefixes)
    if reason == "max_nodes":
        assert result.stages[-1]["marked_t"]
        assert not result.stages[-1]["refinement_applied"]
        assert max(result.stages[-1]["proposed_grid_shape"]) > 9


def test_initial_node_budget_failure_cannot_create_an_empty_success() -> None:
    with pytest.raises(ValueError, match="initial grid exceeds max_nodes"):
        run_adaptive_grid(1.0, config=_config(max_nodes=3))


def test_noninvolutive_lu_permutation_matches_independent_vector_solves() -> None:
    with mp.workprec(160):
        matrix = mp.matrix([[0, 2, 0], [0, 0, 3], [5, 0, 0]])
        right = mp.matrix([[1, 4], [2, 5], [3, 6]])
        observed = _lu_solve(matrix, right)
        for column in range(right.cols):
            reference = mp.lu_solve(matrix, right[:, column])
            assert max(
                abs(observed[i, column] - reference[i]) for i in range(3)
            ) < mp.mpf("1e-40")


def test_cross_probe_matches_independent_direct_formula_and_solve(
    case: AdaptiveGridCase,
) -> None:
    times, frequencies = np.array([0.123, 0.789]), np.array([-0.73, 0.91])
    rank = case.final_result.rank
    observed, diagnostics = _residual_matrix(
        case.final_result, rank, times, frequencies, 160
    )
    high = case.final_result.high_precision
    with mp.workprec(256):

        def direct(t, w):
            return mp.exp(-t * w) / (1 + mp.exp(-w))

        selected_t = [mp.mpf(value) for value in high.t_nodes]
        selected_w = [mp.mpf(value) for value in high.omega_nodes]
        core = mp.matrix([[direct(t, w) for w in selected_w] for t in selected_t])
        for j, w in enumerate(frequencies):
            omega = mp.mpf(repr(float(w)))
            coefficients = mp.lu_solve(
                core, mp.matrix([direct(t, omega) for t in selected_t])
            )
            for i, t in enumerate(times):
                time = mp.mpf(repr(float(t)))
                reference = direct(time, omega) - mp.fsum(
                    direct(time, node) * coefficients[k]
                    for k, node in enumerate(selected_w)
                )
                assert abs(observed[i, j] - reference) < mp.mpf("1e-35")
    assert diagnostics["factorizations"] == 1
    assert diagnostics["solve_residual_rhs_indices"] == [0, 1]


def test_validation_is_only_after_all_marking_and_prefix_choices(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    events = []
    original_probe = MODULE._probe_metrics
    original_validation_nodes = MODULE.validation_nodes
    original_mark = MODULE.mark_panels

    def probe(*args, **kwargs):
        events.append("probe")
        return original_probe(*args, **kwargs)

    def validation(*args, **kwargs):
        events.append("validation")
        return original_validation_nodes(*args, **kwargs)

    def mark(*args, **kwargs):
        events.append("mark")
        return original_mark(*args, **kwargs)

    monkeypatch.setattr(MODULE, "_probe_metrics", probe)
    monkeypatch.setattr(MODULE, "validation_nodes", validation)
    monkeypatch.setattr(MODULE, "mark_panels", mark)
    result = run_adaptive_grid(1.0, config=_config())
    assert result.converged
    first_validation = events.index("validation")
    assert "mark" in events[:first_validation]
    assert all(event == "validation" for event in events[first_validation:])


def test_failed_final_validation_never_changes_controller_or_selected_prefixes(
    monkeypatch: pytest.MonkeyPatch,
    case: AdaptiveGridCase,
) -> None:
    original = MODULE._validation_metrics

    def failed(*args, **kwargs):
        values = original(*args, **kwargs)
        return {**values, "max_residual": "1.0"}

    monkeypatch.setattr(MODULE, "_validation_metrics", failed)
    failed_case = run_adaptive_grid(1.0, config=_config())
    assert failed_case.stages == case.stages
    assert [p["chosen_rank"] for p in failed_case.prefixes] == [
        p["chosen_rank"] for p in case.prefixes
    ]
    assert not failed_case.converged and failed_case.stop_reason == "validation_failed"
    assert all(
        p["joint_passed"] and not p["validation"]["passed"]
        for p in failed_case.prefixes
    )


def test_failed_early_prefix_advances_and_probe_solves_are_cached(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = MODULE._probe_metrics
    calls = Counter()

    def distorted(result, rank, times, frequencies, *args, **kwargs):
        calls[(tuple(times), tuple(frequencies), rank)] += 1
        metrics = original(result, rank, times, frequencies, *args, **kwargs)
        if rank < result.rank:
            metrics["probe_max"] = "1.0"
        return metrics

    monkeypatch.setattr(MODULE, "_probe_metrics", distorted)
    result = run_adaptive_grid(1.0, config=_config())
    assert result.converged
    assert all(
        prefix["chosen_rank"] == result.final_result.rank for prefix in result.prefixes
    )
    assert any(len(prefix["tested_ranks"]) > 1 for prefix in result.prefixes)
    assert max(calls.values()) == 1


def test_decimal_context_cannot_change_generation_threshold_or_scores() -> None:
    config = _config(tolerances=(1.23456789e-10,), generation_margin=0.123456789)
    metrics = {
        "precision_bits": 192,
        "probe_max": "1.23456789123456789e-11",
        "surplus_max": "3.086419725e-11",
        "t_scores": {
            "t0": {"probe_max": "1.23456789123456789e-11", "surplus": "3.086419725e-11"}
        },
        "omega_scores": {},
    }
    outputs = []
    for precision in (2, 80):
        with localcontext() as context:
            context.prec = precision
            context.Emax = 9
            context.Emin = -9
            context.traps[Inexact] = True
            outputs.append(
                (
                    config.generation_tolerance,
                    _classify_probe(metrics, config.tolerances[0], "1e-12"),
                )
            )
    assert outputs[0] == outputs[1]
    assert outputs[0][0] == 1.5241578750190522e-11
