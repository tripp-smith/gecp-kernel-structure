"""Plot rank against cutoff from the canonical census JSONL."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    grouped: dict[str, list[tuple[float, int]]] = defaultdict(list)
    with arguments.input.open(encoding="utf-8") as stream:
        for line in stream:
            record = json.loads(line)
            grouped[record["tolerance"]].append(
                (float(record["cutoff"]), int(record["rank"]))
            )
    for tolerance, points in sorted(grouped.items()):
        points.sort()
        plt.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            marker="o",
            label=f"tol={tolerance}",
        )
    plt.xscale("log")
    plt.xlabel("cutoff")
    plt.ylabel("finite-grid GECP rank")
    plt.legend()
    plt.tight_layout()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(arguments.output, dpi=180)


if __name__ == "__main__":
    main()
