"""Deterministic command-line census for the fermionic GECP geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tomllib
from dataclasses import asdict
from decimal import Decimal
from pathlib import Path
from typing import Any

from .gecp import gecp
from .kernels import FermionicKernel
from .types import CensusRecord, GECPConfig


def _decimal(value: float) -> str:
    return format(float(value), ".17g")


def _git_revision(root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _config_hash(config: dict[str, Any]) -> str:
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def load_census_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate a census TOML file."""

    with path.open("rb") as stream:
        raw: dict[str, Any] = tomllib.load(stream)
    if raw.get("schema_version") != 1:
        raise ValueError("unsupported census schema_version")
    if not raw.get("cutoffs") or not raw.get("tolerances"):
        raise ValueError("census requires nonempty cutoffs and tolerances")
    return raw


def run_census(
    config_path: Path, output_path: Path, *, root: Path | None = None
) -> None:
    """Run every configured case and atomically replace canonical JSONL output."""

    project_root = config_path.resolve().parents[1] if root is None else root.resolve()
    raw = load_census_config(config_path)
    revision = _git_revision(project_root)
    config_hash = _config_hash(raw)
    records: list[CensusRecord] = []
    for cutoff in raw["cutoffs"]:
        tolerances = [float(value) for value in raw["tolerances"]]
        strictest_tolerance = min(tolerances)
        config = GECPConfig(
            tol=strictest_tolerance,
            max_rank=int(raw["max_rank"]),
            pivot=str(raw["pivot"]),  # type: ignore[arg-type]
            precision_bits=int(raw["precision_bits"]),
            grid_order=int(raw["grid_order"]),
            certificate_abs_tol=float(raw.get("certificate_abs_tol", 1e-12)),
            certificate_rel_tol=float(raw.get("certificate_rel_tol", 1e-8)),
            max_cells=int(raw.get("max_cells", 100_000)),
        )
        result = gecp(FermionicKernel(float(cutoff)), config=config)
        high_precision = result.high_precision
        for tolerance in raw["tolerances"]:
            requested_tolerance = Decimal(str(tolerance))
            if high_precision is None:
                rank = next(
                    (
                        index
                        for index, value in enumerate(result.residual_history)
                        if value <= float(tolerance)
                    ),
                    result.rank,
                )
                residual = _decimal(result.residual_history[rank])
                pivots = [_decimal(value) for value in result.pivots[:rank]]
                t_nodes = [_decimal(value) for value in result.t_nodes[:rank]]
                omega_nodes = [_decimal(value) for value in result.omega_nodes[:rank]]
                sigma_min = [
                    _decimal(value) for value in result.core_sigma_min_history[:rank]
                ]
            else:
                rank = next(
                    (
                        index
                        for index, value in enumerate(high_precision.residual_history)
                        if Decimal(value) <= requested_tolerance
                    ),
                    result.rank,
                )
                residual = high_precision.residual_history[rank]
                pivots = high_precision.pivots[:rank]
                t_nodes = high_precision.t_nodes[:rank]
                omega_nodes = high_precision.omega_nodes[:rank]
                sigma_min = high_precision.core_sigma_min_history[:rank]
            converged = Decimal(residual) <= requested_tolerance
            records.append(
                CensusRecord(
                    schema_version=1,
                    git_commit=revision,
                    config_hash=config_hash,
                    cutoff=_decimal(float(cutoff)),
                    tolerance=_decimal(float(tolerance)),
                    precision_bits=result.precision_bits,
                    algorithm=config.pivot,
                    rank=rank,
                    residual=residual,
                    pivots=pivots,
                    t_nodes=t_nodes,
                    omega_nodes=omega_nodes,
                    core_sigma_min=sigma_min,
                    tie_counts=result.tie_counts[:rank],
                    converged=converged,
                    stop_reason="tolerance" if converged else result.stop_reason,
                )
            )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        for record in records:
            stream.write(
                json.dumps(asdict(record), sort_keys=True, separators=(",", ":"))
            )
            stream.write("\n")
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("experiments/census.toml"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/data/gecp_census.jsonl"),
    )
    arguments = parser.parse_args()
    run_census(arguments.config, arguments.output)


if __name__ == "__main__":
    main()
