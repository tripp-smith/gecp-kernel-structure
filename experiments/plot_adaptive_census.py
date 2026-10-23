"""Plot residual-adaptive sampled-grid histories, ranks and actual stages."""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from experiments.check_adaptive_census import (
    CANONICAL_CUTOFFS,
    CANONICAL_TOLERANCES,
    _axis,
    _decimal,
    _list,
    _object,
    _validate_record,
)


@dataclass(frozen=True)
class AdaptiveCurve:
    """Deepest recorded training trajectory and all final validation samples."""

    cutoff: float
    ranks: tuple[int, ...]
    normalized_ranks: tuple[float, ...]
    residuals: tuple[float, ...]
    validation_ranks: tuple[int, ...]
    validation_errors: tuple[float, ...]
    validation_passed: tuple[bool, ...]
    converged: bool
    stage_count: int

    @property
    def zero_count(self) -> int:
        return self.residuals.count(0) + self.validation_errors.count(0)

    @property
    def nonmonotonic(self) -> bool:
        return any(
            a < b for a, b in zip(self.residuals, self.residuals[1:], strict=False)
        )


@dataclass(frozen=True)
class AdaptivePlotSummary:
    curve_count: int
    omitted_zero_count: int
    unsuccessful_case_count: int
    nonmonotonic_trajectory_count: int
    stage_count: int


def _float(value: Any, name: str) -> float:
    number = _decimal(value, name, nonnegative=True)
    result = float(number)
    if not math.isfinite(result) or (number != 0 and result == 0):
        raise ValueError(f"{name} outside the finite positive plotting range")
    return result


def read_adaptive_census(path: Path) -> Mapping[str, Any]:
    try:
        result = json.loads(path.read_bytes())
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid adaptive census JSON: {error}") from error
    return _object(result, "adaptive census")


def normalized_curves(
    payload: Mapping[str, Any], *, require_full: bool = True
) -> list[AdaptiveCurve]:
    """Retain failures and nonmonotonic samples instead of filtering them away."""
    if (
        payload.get("schema_version") != 3
        or payload.get("kind") != "residual_adaptive_sampled_grid"
    ):
        raise ValueError("plot requires schema-3 residual-adaptive sampled-grid data")
    cases = _list(payload.get("cases"), "cases")
    if not cases:
        raise ValueError("adaptive census must contain at least one case")
    curves = []
    seen = set()
    for case in sorted(cases, key=lambda item: _decimal(item.get("cutoff"), "cutoff")):
        case = _object(case, "case")
        cutoff = _decimal(case.get("cutoff"), "cutoff")
        if cutoff < 1 or cutoff in seen:
            raise ValueError("cutoffs must be distinct and at least one")
        seen.add(cutoff)
        axes = _object(case.get("final_axes"), "final_axes")
        for axis in ("t", "omega"):
            _axis(axes.get(axis), cutoff, axis)
        records = _list(case.get("records"), "records")
        prefixes = _list(case.get("prefixes"), "prefixes")
        if not records or len(records) != len(prefixes):
            raise ValueError("records and prefixes must be nonempty and aligned")
        for record in records:
            _validate_record(_object(record, "record"))
            if _decimal(record["cutoff"], "cutoff") != cutoff:
                raise ValueError("record cutoff disagrees with case")
        deepest = max(records, key=lambda record: record["rank"])
        for record in records:
            if (
                record["residual_history"]
                != deepest["residual_history"][: record["rank"] + 1]
            ):
                raise ValueError("incompatible final-grid training histories")
        tolerances = {_decimal(record["tolerance"], "tolerance") for record in records}
        if len(tolerances) != len(records):
            raise ValueError("duplicate tolerance records")
        if require_full and tolerances != set(CANONICAL_TOLERANCES):
            raise ValueError(
                "full adaptive plot requires all four requested tolerances"
            )
        ranks, errors, passed = [], [], []
        for record, prefix in zip(records, prefixes, strict=True):
            prefix = _object(prefix, "prefix")
            validation = _object(prefix.get("validation"), "validation")
            if (
                prefix.get("chosen_rank") != record["rank"]
                or prefix.get("tolerance") != record["tolerance"]
            ):
                raise ValueError("prefix rank/tolerance disagrees with record")
            if (
                type(prefix.get("converged")) is not bool
                or type(validation.get("passed")) is not bool
            ):
                raise ValueError(
                    "prefix convergence and validation pass must be boolean"
                )
            ranks.append(record["rank"])
            errors.append(_float(validation.get("max_residual"), "validation max"))
            passed.append(prefix["converged"])
        if type(case.get("converged")) is not bool:
            raise ValueError("case convergence must be boolean")
        stages = _list(case.get("stages"), "stages")
        if not stages:
            raise ValueError("adaptive plot requires recorded stages")
        cutoff_float = _float(str(cutoff), "cutoff")
        denominator = math.log1p(cutoff_float)
        rank_axis = tuple(range(deepest["rank"] + 1))
        curves.append(
            AdaptiveCurve(
                cutoff_float,
                rank_axis,
                tuple(rank / denominator for rank in rank_axis),
                tuple(
                    _float(value, "training residual")
                    for value in deepest["residual_history"]
                ),
                tuple(ranks),
                tuple(errors),
                tuple(passed),
                case["converged"],
                len(stages),
            )
        )
    if require_full and seen != set(CANONICAL_CUTOFFS):
        raise ValueError("full adaptive plot requires all seven cutoffs")
    return curves


def _positive(values: list[float] | tuple[float, ...]) -> list[float]:
    return [value if value > 0 else math.nan for value in values]


def plot_adaptive_census(
    payload: Mapping[str, Any], output: Path, *, require_full: bool = True
) -> AdaptivePlotSummary:
    """Save deterministic stage-aware plots; numerical zeros have no fake floor."""
    curves = normalized_curves(payload, require_full=require_full)
    cases = sorted(payload["cases"], key=lambda case: Decimal(case["cutoff"]))
    zero_count = sum(curve.zero_count for curve in curves)
    failures = sum(not curve.converged for curve in curves)
    nonmonotonic = sum(curve.nonmonotonic for curve in curves)
    figure, axes = plt.subplots(2, 2, figsize=(14, 9))
    colors = plt.get_cmap("tab10")
    try:
        residual_axis, rank_axis, stage_axis, nodes_axis = axes.flat
        residual_axis.set_yscale("log")
        stage_axis.set_yscale("log")
        for index, (curve, case) in enumerate(zip(curves, cases, strict=True)):
            color = colors(index)
            label = f"Λ={curve.cutoff:g}" + (" (failed)" if not curve.converged else "")
            residual_axis.plot(
                curve.normalized_ranks,
                _positive(curve.residuals),
                color=color,
                label=label,
            )
            for rank, error, passed in zip(
                curve.validation_ranks,
                curve.validation_errors,
                curve.validation_passed,
                strict=True,
            ):
                residual_axis.plot(
                    [rank / math.log1p(curve.cutoff)],
                    _positive([error]),
                    linestyle="none",
                    marker="D" if passed else "x",
                    markersize=5,
                    color=color,
                    markerfacecolor="none",
                )
            stages = case["stages"]
            indices = list(range(1, len(stages) + 1))
            for name, style in (
                ("training_max", "-"),
                ("probe_max", "--"),
                ("surplus_max", ":"),
            ):
                values = [_float(stage["probe"][name], name) for stage in stages]
                zero_count += values.count(0)
                stage_axis.plot(
                    indices, _positive(values), style, color=color, marker="."
                )
            nodes_axis.plot(
                indices,
                [len(stage["t_grid"]) for stage in stages],
                "-o",
                color=color,
                markersize=3,
                label=label,
            )
            nodes_axis.plot(
                indices,
                [len(stage["omega_grid"]) for stage in stages],
                "--s",
                color=color,
                markersize=3,
            )
        tolerances = sorted(
            {
                Decimal(prefix["tolerance"])
                for case in cases
                for prefix in case["prefixes"]
            },
            reverse=True,
        )
        for index, tolerance in enumerate(tolerances):
            coordinates, training, chosen, success = [], [], [], []
            for case in cases:
                prefix = next(
                    (
                        p
                        for p in case["prefixes"]
                        if Decimal(p["tolerance"]) == tolerance
                    ),
                    None,
                )
                if prefix is not None:
                    coordinates.append(float(case["cutoff"]))
                    training.append(prefix["training_rank"])
                    chosen.append(prefix["chosen_rank"])
                    success.append(prefix["converged"])
            color = colors(index)
            rank_axis.plot(
                coordinates,
                chosen,
                "-o",
                color=color,
                markersize=4,
                label=f"tol={float(tolerance):.0e}",
            )
            rank_axis.plot(coordinates, training, "--", color=color)
            for cutoff, rank, passed in zip(coordinates, chosen, success, strict=True):
                if not passed:
                    rank_axis.plot(cutoff, rank, "x", color="black", markersize=7)
        for axis in (residual_axis, stage_axis):
            if not any(
                math.isfinite(y) and y > 0
                for line in axis.lines
                for y in line.get_ydata()
            ):
                axis.set_ylim(0.1, 1)
            axis.grid(alpha=0.2)
        residual_axis.set(
            xlabel="k / log(1 + cutoff)",
            ylabel="sampled |Rₖ|",
            title="Final training histories and untouched validation",
        )
        cutoff_legend = residual_axis.legend(fontsize=8, ncol=2)
        residual_axis.add_artist(cutoff_legend)
        residual_axis.add_artist(
            residual_axis.legend(
                handles=[
                    Line2D([], [], color="black", label="training", linestyle="-"),
                    Line2D(
                        [],
                        [],
                        color="black",
                        label="final validation",
                        marker="D",
                        markerfacecolor="none",
                        linestyle="none",
                    ),
                    Line2D(
                        [],
                        [],
                        color="black",
                        label="unsuccessful prefix",
                        marker="x",
                        linestyle="none",
                    ),
                ],
                loc="lower left",
                fontsize=7,
            )
        )
        rank_axis.set(
            xscale="log",
            xlabel="cutoff Λ",
            ylabel="sampled rank",
            title="Solid: chosen joint rank; dashed: first training rank",
        )
        rank_axis.legend(fontsize=8)
        rank_axis.grid(alpha=0.2)
        stage_axis.set(
            xlabel="actual stage",
            ylabel="sampled maximum / surplus",
            title="Solid: training; dashed: probe; dotted: surplus",
        )
        nodes_axis.set(
            xlabel="actual stage",
            ylabel="explicit axis node count",
            title="Solid: time nodes; dashed: frequency nodes",
        )
        nodes_axis.grid(alpha=0.2)
        for axis in (stage_axis, nodes_axis):
            axis.set_xticks(range(1, max(curve.stage_count for curve in curves) + 1))
        training_bits = sorted(
            {record["precision_bits"] for case in cases for record in case["records"]}
        )
        checking_bits = sorted(
            {
                prefix["validation"]["precision_bits"]
                for case in cases
                for prefix in case["prefixes"]
            }
        )
        figure.suptitle(
            "Residual-adaptive sampled GECP — "
            f"MP{'/'.join(map(str, training_bits))} training / "
            f"MP{'/'.join(map(str, checking_bits))} checks\n"
            "Finite samples only; no continuous error certificate",
            fontsize=13,
        )
        figure.text(
            0.02,
            0.014,
            f"{zero_count} zero samples omitted from log axes; "
            "no positive display floor. "
            f"{nonmonotonic} nonmonotonic within-run histories retained. "
            f"{failures} unsuccessful cases retained (× marks unsuccessful prefixes).",
            fontsize=8,
        )
        figure.tight_layout(rect=(0, 0.035, 1, 0.94))
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(
            output, dpi=160, metadata={"Software": "kernelgecp adaptive census"}
        )
    finally:
        plt.close(figure)
    return AdaptivePlotSummary(
        len(curves),
        zero_count,
        failures,
        nonmonotonic,
        sum(curve.stage_count for curve in curves),
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path("experiments/data/adaptive_census.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("experiments/data/adaptive_census.png")
    )
    parser.add_argument(
        "--allow-partial", action="store_true", help="labelled pilot data only"
    )
    args = parser.parse_args(argv)
    try:
        plot_adaptive_census(
            read_adaptive_census(args.input),
            args.output,
            require_full=not args.allow_partial,
        )
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
