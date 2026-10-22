"""Plot every retained fit, including failures; never floor numerical zeros."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def plot(data: dict[str, Any], output: Path) -> None:
    if data.get("schema_version") != 1 or len(data.get("fixtures", [])) != 4:
        raise ValueError("expected the four-fixture schema-1 comparison")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True)
    omitted = 0
    try:
        for axis, fixture in zip(axes.flat, data["fixtures"], strict=True):
            for method in data["methods"]:
                records = [r for r in fixture["results"] if r["method"] == method]
                valid = [r for r in records if r["held_out_error"] > 0]
                omitted += len(records) - len(valid)
                axis.semilogy(
                    [r["atom_count"] for r in valid],
                    [r["held_out_error"] for r in valid],
                    "o-",
                    markersize=3,
                    label=method.replace("_", " "),
                )
            for tolerance in data["config"]["tolerances"]:
                axis.axhline(tolerance, color="0.65", linestyle=":", linewidth=0.8)
            axis.set_title(fixture["name"].replace("_", " "))
            axis.set_xlabel("Returned spectral atoms (tested budgets)")
            axis.set_ylabel("Noiseless held-out max error")
            axis.grid(True, alpha=0.2)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
        fig.suptitle("Six methods, shared time observations; failures retained")
        fig.text(
            0.5,
            0.09,
            "Sampled errors, not continuous/support certificates. "
            f"Zero-error samples omitted: {omitted}.",
            ha="center",
            fontsize=9,
        )
        fig.tight_layout(rect=(0, 0.12, 1, 0.95))
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            output, dpi=140, metadata={"Software": "kernelgecp spectral comparison"}
        )
    finally:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path("experiments/data/spectral_comparison.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("experiments/data/spectral_comparison.png")
    )
    args = parser.parse_args()
    plot(json.loads(args.input.read_text(encoding="utf-8")), args.output)


if __name__ == "__main__":
    main()
