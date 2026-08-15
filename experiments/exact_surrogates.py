"""Generate canonical exact-surrogate evidence for sizes two through eight."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from kernelgecp.surrogates import exact_surrogate_census


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/data/exact_surrogates.json"),
    )
    arguments = parser.parse_args()
    records = [record.as_dict() for record in exact_surrogate_census()]
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
