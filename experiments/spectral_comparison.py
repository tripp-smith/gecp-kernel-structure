"""Reproducible six-method spectral comparison, not an inverse-problem proof."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import tomllib
from pathlib import Path

from kernelgecp.comparison import ComparisonConfig, run_spectral_comparison


def run(config_path: Path, output_path: Path, *, require_clean: bool = True) -> None:
    root = Path(__file__).resolve().parents[1]
    raw = tomllib.loads(config_path.read_text(encoding="utf-8"))
    for key in ("atom_budgets", "tolerances"):
        if key in raw:
            raw[key] = tuple(raw[key])
    config = ComparisonConfig(**raw)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()
    dirty = bool(
        subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=root,
            text=True,
        ).strip()
    )
    if require_clean and dirty:
        raise ValueError("canonical comparison requires a clean tracked revision")
    result = run_spectral_comparison(config, git_commit=revision)
    result["git_dirty"] = dirty
    result["environment"] = {
        "python": platform.python_version(),
        **{
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "mpmath")
        },
    }
    result["lock_sha256"] = hashlib.sha256((root / "uv.lock").read_bytes()).hexdigest()
    result["algorithm_source_sha256"] = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in (
            "kernelgecp/comparison.py",
            "kernelgecp/rational.py",
            "kernelgecp/prony.py",
            "kernelgecp/vector_fitting.py",
            "kernelgecp/spectral_common.py",
            "kernelgecp/sparse.py",
            "kernelgecp/gecp.py",
            "experiments/spectral_comparison.py",
        )
    }
    encoded = (
        json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(encoded, encoding="utf-8", newline="\n")
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("experiments/spectral_comparison.toml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("experiments/data/spectral_comparison.json")
    )
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="exploratory only; dirty status is recorded",
    )
    args = parser.parse_args()
    run(args.config, args.output, require_clean=not args.allow_dirty)


if __name__ == "__main__":
    main()
