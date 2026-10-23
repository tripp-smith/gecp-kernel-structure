"""Adaptive plots retain stage changes, nonmonotonic samples and failures."""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any

import pytest

from experiments import plot_adaptive_census as plotting


def _payload() -> dict[str, Any]:
    history = ["1", "0.1", "0.2", "0"]
    records, prefixes = [], []
    for rank, tolerance, passed, error in (
        (2, "1e-6", False, "0.3"),
        (3, "1e-12", True, "0"),
    ):
        records.append(
            {
                "schema_version": 3,
                "algorithm": "residual_adaptive_grid",
                "cutoff": "4",
                "tolerance": tolerance,
                "precision_bits": 128,
                "rank": rank,
                "residual": history[rank],
                "residual_history": history[: rank + 1],
                "t_nodes": ["0", "0.5", "1"][:rank],
                "omega_nodes": ["-1", "0", "1"][:rank],
                "pivots": ["1"] * rank,
                "core_det_history": ["1"] * rank,
                "core_sigma_min": ["1"] * rank,
                "tie_counts": [1] * rank,
                "converged": passed,
                "stop_reason": "sampled_tolerance" if passed else "validation_failed",
            }
        )
        prefixes.append(
            {
                "tolerance": tolerance,
                "chosen_rank": rank,
                "training_rank": rank,
                "converged": passed,
                "validation": {
                    "max_residual": error,
                    "passed": passed,
                    "precision_bits": 192,
                },
            }
        )
    return {
        "schema_version": 3,
        "kind": "residual_adaptive_sampled_grid",
        "cases": [
            {
                "cutoff": "4",
                "converged": False,
                "stop_reason": "validation_failed",
                "final_axes": {
                    "t": ["0.0", "0.5", "1.0"],
                    "omega": ["-4.0", "0.0", "4.0"],
                },
                "records": records,
                "prefixes": prefixes,
                "stages": [
                    {
                        "t_grid": ["0", "1"],
                        "omega_grid": ["-4", "4"],
                        "probe": {
                            "training_max": "1",
                            "probe_max": "0.8",
                            "surplus_max": "0.5",
                        },
                    },
                    {
                        "t_grid": ["0", ".5", "1"],
                        "omega_grid": ["-4", "0", "4"],
                        "probe": {
                            "training_max": "0",
                            "probe_max": "0.01",
                            "surplus_max": "0",
                        },
                    },
                ],
            }
        ],
    }


def test_curves_use_log1p_and_retain_nonmonotonic_zeros_and_failures() -> None:
    curve = plotting.normalized_curves(_payload(), require_full=False)[0]
    assert curve.normalized_ranks == pytest.approx(
        tuple(k / math.log1p(4) for k in range(4))
    )
    assert curve.residuals == (1, 0.1, 0.2, 0)
    assert curve.zero_count == 2
    assert curve.nonmonotonic
    assert not curve.converged
    assert curve.validation_passed == (False, True)
    assert curve.stage_count == 2
    with pytest.raises(ValueError, match="all four"):
        plotting.normalized_curves(_payload())


@pytest.mark.parametrize(
    "alteration", ["scope", "history", "axis", "alignment", "zero_precision"]
)
def test_invalid_adaptive_plot_data_rejected(alteration: str) -> None:
    payload = _payload()
    case = payload["cases"][0]
    if alteration == "scope":
        payload["kind"] = "continuous_interval_certificate"
    elif alteration == "history":
        case["records"][0]["residual_history"][1] = "0.9"
    elif alteration == "axis":
        case["final_axes"]["t"][-1] = "0.9"
    elif alteration == "alignment":
        case["prefixes"][0]["chosen_rank"] += 1
    else:
        case["records"][0]["precision_bits"] = 0
    with pytest.raises(ValueError):
        plotting.normalized_curves(payload, require_full=False)


def test_plot_repeat_pngs_mask_zeros_and_disclose_scope(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed: dict[str, Any] = {}
    original = plotting.plt.Figure.savefig

    def inspect(figure: Any, *args: Any, **kwargs: Any) -> None:
        observed["titles"] = [axis.get_title() for axis in figure.axes]
        observed["texts"] = [text.get_text() for text in figure.texts]
        observed["ys"] = figure.axes[0].lines[0].get_ydata()
        observed["legends"] = [
            text.get_text() for text in figure.axes[0].get_legend().get_texts()
        ]
        original(figure, *args, **kwargs)

    monkeypatch.setattr(plotting.plt.Figure, "savefig", inspect)
    before = plotting.plt.get_fignums()
    first, second = tmp_path / "a.png", tmp_path / "b.png"
    summary = plotting.plot_adaptive_census(_payload(), first, require_full=False)
    plotting.plot_adaptive_census(_payload(), second, require_full=False)
    assert first.read_bytes() == second.read_bytes()
    assert summary.curve_count == 1
    assert summary.omitted_zero_count == 4
    assert summary.unsuccessful_case_count == 1
    assert summary.nonmonotonic_trajectory_count == 1
    assert math.isnan(observed["ys"][-1])
    assert "untouched validation" in observed["titles"][0]
    assert "chosen joint rank" in observed["titles"][1]
    assert "final validation" in observed["legends"]
    assert any("no continuous error certificate" in text for text in observed["texts"])
    assert any("no positive display floor" in text for text in observed["texts"])
    assert plotting.plt.get_fignums() == before


def test_failed_save_closes_all_plot_figures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*args: Any, **kwargs: Any) -> None:
        raise OSError("test save failure")

    monkeypatch.setattr(plotting.plt.Figure, "savefig", fail)
    before = plotting.plt.get_fignums()
    with pytest.raises(OSError, match="test save failure"):
        plotting.plot_adaptive_census(
            _payload(), tmp_path / "fail.png", require_full=False
        )
    assert plotting.plt.get_fignums() == before


def test_all_zero_and_missing_training_pass_rank_remain_explicit(
    tmp_path: Path,
) -> None:
    payload = copy.deepcopy(_payload())
    for record in payload["cases"][0]["records"]:
        record["residual_history"] = ["0"] * (record["rank"] + 1)
        record["residual"] = "0"
    payload["cases"][0]["prefixes"][0]["training_rank"] = None
    report = plotting.plot_adaptive_census(
        payload, tmp_path / "zero.png", require_full=False
    )
    assert report.omitted_zero_count == 7
