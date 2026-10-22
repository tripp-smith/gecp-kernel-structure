"""Plot sampled GECP ranks or residuals normalized by log(1 + cutoff)."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class NormalizedCurve:
    """One deepest-tolerance trajectory, including its rank-zero residual."""

    cutoff: float
    tolerance: float
    precision_bits: int
    ranks: tuple[int, ...]
    normalized_ranks: tuple[float, ...]
    residuals: tuple[float, ...]

    @property
    def zero_indices(self) -> tuple[int, ...]:
        return tuple(index for index, value in enumerate(self.residuals) if value == 0)


@dataclass(frozen=True)
class NormalizedPlotSummary:
    """Display metadata; omitted zero samples are not positive error estimates."""

    curve_count: int
    omitted_zero_count: int
    precision_bits: tuple[int, ...]


def _decimal(value: Any, name: str, *, positive: bool = False) -> Decimal:
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError):
        raise ValueError(f"{name} must be a finite decimal") from None
    if not result.is_finite() or result < 0 or (positive and result == 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return result


def _as_float(value: Decimal, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or (value != 0 and result == 0):
        raise ValueError(f"{name} is outside the positive float plotting range")
    return result


def _integer(value: Any, name: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if value < (1 if positive else 0):
        raise ValueError(f"{name} must be {'positive' if positive else 'nonnegative'}")
    return value


def _validate_sampled_record(record: Mapping[str, Any]) -> None:
    """Do not label interval bounds or unsuccessful runs as grid evidence.

    Missing fields remain accepted for old, minimal census records. Explicit
    algorithm, scope and convergence metadata must agree with the plot labels.
    """

    if "algorithm" in record and record["algorithm"] != "grid":
        raise ValueError("census plots require the grid algorithm, not adaptive bounds")
    metadata = record.get("metadata", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("census plot metadata must be an object")
    if (
        "convergence_scope" in metadata
        and metadata["convergence_scope"] != "sampled_grid"
    ):
        raise ValueError("census plots require sampled_grid convergence scope")
    for name, source in (("record", record), ("metadata", metadata)):
        if "converged" in source and source["converged"] is not True:
            raise ValueError(
                "census plots require converged records; "
                f"{name} does not establish the rank required for its tolerance"
            )


def read_census(path: Path) -> list[dict[str, Any]]:
    """Read JSONL without depending on optional, evolving census metadata."""

    records = []
    with path.open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"invalid census JSON at line {number}: {error}"
                ) from error
            if not isinstance(record, dict):
                raise ValueError(f"census line {number} must be an object")
            records.append(record)
    if not records:
        raise ValueError("census must contain at least one record")
    return records


def normalized_curves(records: Sequence[Mapping[str, Any]]) -> list[NormalizedCurve]:
    """Select one compatible strictest-tolerance prefix for each numeric cutoff.

    All histories are checked before deduplication, so incompatible or truncated
    prefixes cannot silently disappear behind a valid deepest trajectory.
    """

    if not records:
        raise ValueError("census must contain at least one record")
    grouped: dict[Decimal, list[tuple[Decimal, int, int, tuple[Decimal, ...]]]] = (
        defaultdict(list)
    )
    for record in records:
        _validate_sampled_record(record)
        if record.get("schema_version") != 2:
            raise ValueError(
                "normalized mode requires schema_version 2 residual histories"
            )
        try:
            cutoff = _decimal(record["cutoff"], "cutoff", positive=True)
            tolerance = _decimal(record["tolerance"], "tolerance", positive=True)
            rank = _integer(record["rank"], "rank")
            precision = _integer(
                record["precision_bits"], "precision_bits", positive=True
            )
            raw_history = record["residual_history"]
        except KeyError as error:
            raise ValueError(f"normalized mode requires {error.args[0]}") from error
        if not isinstance(raw_history, list) or len(raw_history) != rank + 1:
            raise ValueError("residual_history must contain exactly rank + 1 values")
        if not all(isinstance(value, str) for value in raw_history):
            raise ValueError("schema-2 residual_history values must be decimal strings")
        history = tuple(_decimal(value, "residual_history") for value in raw_history)
        grouped[cutoff].append((tolerance, rank, precision, history))

    curves = []
    for cutoff, candidates in sorted(grouped.items()):
        candidates.sort(key=lambda item: (item[0], -item[1]))
        tolerance, rank, precision, history = candidates[0]
        for _, prefix_rank, prefix_precision, prefix_history in candidates:
            if (
                prefix_precision != precision
                or prefix_rank > rank
                or prefix_history != history[: prefix_rank + 1]
            ):
                raise ValueError(
                    f"incompatible residual-history prefixes for cutoff {cutoff}"
                )
        cutoff_float = _as_float(cutoff, "cutoff")
        denominator = math.log1p(cutoff_float)
        ranks = tuple(range(rank + 1))
        normalized = tuple(index / denominator for index in ranks)
        if not all(math.isfinite(value) for value in normalized):
            raise ValueError(
                "normalized rank is outside the finite float plotting range"
            )
        curves.append(
            NormalizedCurve(
                cutoff=cutoff_float,
                tolerance=_as_float(tolerance, "tolerance"),
                precision_bits=precision,
                ranks=ranks,
                normalized_ranks=normalized,
                residuals=tuple(
                    _as_float(value, "residual_history") for value in history
                ),
            )
        )
    return curves


def _precision_label(precisions: Sequence[int]) -> str:
    return ", ".join(str(value) for value in sorted(set(precisions))) + " bits"


def plot_normalized(
    records: Sequence[Mapping[str, Any]], output: Path
) -> NormalizedPlotSummary:
    """Save log residual curves; numerical zeros are masked, never floored."""

    curves = normalized_curves(records)
    precisions = tuple(sorted({curve.precision_bits for curve in curves}))
    zero_count = sum(len(curve.zero_indices) for curve in curves)
    figure, axis = plt.subplots()
    try:
        axis.set_yscale("log")
        for curve in curves:
            axis.plot(
                curve.normalized_ranks,
                [value if value > 0 else math.nan for value in curve.residuals],
                marker=".",
                label=f"Λ={curve.cutoff:g}, tol={curve.tolerance:.0e}",
            )
        if not any(value > 0 for curve in curves for value in curve.residuals):
            axis.set_ylim(0.1, 1)
        if zero_count:
            axis.text(
                0.02,
                0.02,
                f"{zero_count} numerically zero sample(s) omitted from log scale;\n"
                "no positive display floor was substituted",
                transform=axis.transAxes,
                fontsize="small",
            )
        axis.set_xlabel("k / log(1 + cutoff)")
        axis.set_ylabel("sampled finite-grid residual maximum |Rₖ|")
        axis.set_title(
            "Finite-grid GECP residuals (not continuum certificates)\n"
            f"observed precision: {_precision_label(precisions)}"
        )
        axis.legend()
        figure.tight_layout()
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=180)
    finally:
        plt.close(figure)
    return NormalizedPlotSummary(len(curves), zero_count, precisions)


def plot_ranks(records: Sequence[Mapping[str, Any]], output: Path) -> None:
    """Save the original rank-versus-cutoff plot, accepting legacy schema 1."""

    if not records:
        raise ValueError("census must contain at least one record")
    grouped: dict[Decimal, list[tuple[float, int]]] = defaultdict(list)
    precisions = []
    for record in records:
        _validate_sampled_record(record)
        try:
            cutoff = _as_float(
                _decimal(record["cutoff"], "cutoff", positive=True), "cutoff"
            )
            tolerance = _decimal(record["tolerance"], "tolerance", positive=True)
            rank = _integer(record["rank"], "rank")
            precisions.append(
                _integer(record["precision_bits"], "precision_bits", positive=True)
            )
        except KeyError as error:
            raise ValueError(f"rank mode requires {error.args[0]}") from error
        grouped[tolerance].append((cutoff, rank))
    figure, axis = plt.subplots()
    try:
        for tolerance, points in sorted(grouped.items()):
            points.sort()
            axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                marker="o",
                label=f"tol={float(tolerance):.0e}",
            )
        axis.set_xscale("log")
        axis.set_xlabel("cutoff")
        axis.set_ylabel("finite-grid GECP rank")
        axis.set_title(
            "Sampled finite-grid fermionic GECP census\n"
            f"observed precision: {_precision_label(precisions)}"
        )
        axis.legend()
        figure.tight_layout()
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=180)
    finally:
        plt.close(figure)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("ranks", "normalized"), default="ranks")
    arguments = parser.parse_args(argv)
    try:
        records = read_census(arguments.input)
        if arguments.mode == "normalized":
            plot_normalized(records, arguments.output)
        else:
            plot_ranks(records, arguments.output)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
