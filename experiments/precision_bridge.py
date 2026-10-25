"""Freeze-aware deterministic finite-query arithmetic bridge evidence.

Run as a module from a committed, clean tracked source revision. This is an
arithmetic certificate demonstration plus independent QR consistency audit;
physical pole errors remain finite-sample observations.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import tomllib
from pathlib import Path

from examples.precision_bridge import audited_case
from experiments.dlr_benchmark import source_snapshot


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    config = args.config.resolve()
    raw = tomllib.loads(config.read_text())
    if (
        set(raw) != {"schema_version", "cases"}
        or type(raw["schema_version"]) is not int
        or raw["schema_version"] != 1
    ):
        raise ValueError("unsupported arithmetic bridge config")
    if not isinstance(raw["cases"], list) or not raw["cases"]:
        raise ValueError("cases must be a nonempty list")
    for case in raw["cases"]:
        if not {"beta", "wmax", "tolerance"} <= set(case):
            raise ValueError("case must specify beta, wmax and tolerance")

    def snapshot() -> dict:
        result = source_snapshot(root, config)
        path = root / "examples/precision_bridge.py"
        result["source_sha256"]["examples/precision_bridge.py"] = hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        return result

    before = snapshot()
    records = [audited_case(**case) for case in raw["cases"]]
    if snapshot() != before:
        raise RuntimeError("source/config changed during certificate generation")
    result = {
        "schema_version": 1,
        "source": before,
        "config": raw,
        "python": platform.python_version(),
        "versions": {
            package: importlib.metadata.version(package)
            for package in ("kernelgecp", "numpy", "scipy", "mpmath")
        },
        "records": records,
        "scope": "finite_query_arithmetic;independent_QR;physical_errors_sampled",
    }
    text = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as output:
        output.write(text)
    print(f"Wrote {len(records)} independently checked cases to {args.output}")


if __name__ == "__main__":
    main()
