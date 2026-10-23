"""Generate explicit-axis, residual-adaptive sampled GECP evidence, not bounds."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tomllib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from decimal import Decimal
from itertools import repeat
from pathlib import Path
from typing import Any

from kernelgecp.adaptive_grid import (
    AdaptiveGridCase,
    AdaptiveGridConfig,
    run_adaptive_grid,
)
from kernelgecp.census import _config_hash, _git_revision, _pivot_geometry, _provenance
from kernelgecp.types import CensusRecord

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    "kernelgecp/adaptive_grid.py",
    "kernelgecp/adaptive_grids.py",
    "kernelgecp/gecp.py",
    "kernelgecp/census.py",
    "kernelgecp/kernels.py",
    "kernelgecp/types.py",
    "experiments/adaptive_census.py",
)


def load_config(path: Path) -> tuple[dict[str, Any], AdaptiveGridConfig]:
    with path.open("rb") as stream:
        raw: dict[str, Any] = tomllib.load(stream)
    if type(raw.get("schema_version")) is not int or raw["schema_version"] != 3:
        raise ValueError("adaptive census requires schema_version 3")
    cutoffs = raw.get("cutoffs")
    if not isinstance(cutoffs, list) or not cutoffs:
        raise ValueError("cutoffs must be a nonempty list")
    for value in cutoffs:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not Decimal(str(value)).is_finite()
            or value < 1
        ):
            raise ValueError("cutoffs must be finite and at least one")
    if len(set(cutoffs)) != len(cutoffs):
        raise ValueError("duplicate cutoffs")
    if type(raw.get("workers", 1)) is not int or raw.get("workers", 1) < 1:
        raise ValueError("workers must be a positive integer")
    if type(raw.get("require_clean_revision", True)) is not bool:
        raise ValueError("require_clean_revision must be boolean")
    controls = {
        name: value
        for name, value in raw.items()
        if name
        not in ("schema_version", "cutoffs", "workers", "require_clean_revision")
    }
    if "tolerances" in controls:
        if not isinstance(controls["tolerances"], list):
            raise ValueError("tolerances must be a list")
        controls["tolerances"] = tuple(controls["tolerances"])
    return raw, AdaptiveGridConfig(**controls)


def _record(
    case: AdaptiveGridCase,
    prefix: dict[str, Any],
    revision: str,
    config_hash: str,
    provenance: dict[str, Any],
) -> dict[str, Any]:
    """Serialize the actual jointly sampled prefix, never infer a fixed grid."""
    result = case.final_result
    high = result.high_precision
    if high is None:
        raise ValueError("adaptive census requires canonical MP histories")
    rank = prefix["chosen_rank"]
    if type(rank) is not int or not 0 <= rank <= result.rank:
        raise ValueError("invalid sampled prefix rank")
    for name in (
        "t_nodes",
        "omega_nodes",
        "pivots",
        "core_det_history",
        "core_sigma_min_history",
    ):
        if len(getattr(high, name)) != result.rank:
            raise ValueError(f"misaligned {name}")
    if len(high.residual_history) != result.rank + 1:
        raise ValueError("misaligned residual history")
    source = result.metadata
    conditioning = source.conditioning_history[:rank]
    converged = bool(prefix["converged"])
    stop = str(prefix["stop_reason"])
    metadata = replace(
        source,
        algorithm="gecp_residual_adaptive_grid",
        conditioning_history=conditioning,
        conditioning_status=(
            "not_representable"
            if any(item.status == "not_representable" for item in conditioning)
            else "estimated"
            if conditioning
            else "not_applicable"
        ),
        tie_counts=source.tie_counts[:rank],
        convergence_scope="sampled_training_probe_and_independent_validation",
        converged=converged,
        stop_reason=stop,
        certification_status="not_certified",
        certification_scope="sampled_only_no_continuous_bound",
        configuration={
            **source.configuration,
            "tol": float(prefix["tolerance"]),
            "source_trajectory_tolerance": result.config["tol"],
            "prefix_extraction": True,
            "prefix_rule": "first_tested_training_and_probe_and_surplus_passing",
            "final_validation_used_for_refinement": False,
            "coordinate_precision_bits": 53,
            "check_precision_bits": prefix["validation"]["precision_bits"],
        },
    )
    times, frequencies, pivots = (
        high.t_nodes[:rank],
        high.omega_nodes[:rank],
        high.pivots[:rank],
    )
    record = CensusRecord(
        schema_version=3,
        git_commit=revision,
        config_hash=config_hash,
        cutoff=str(case.cutoff),
        tolerance=str(prefix["tolerance"]),
        precision_bits=result.precision_bits,
        algorithm="residual_adaptive_grid",
        rank=rank,
        residual=high.residual_history[rank],
        pivots=pivots,
        t_nodes=times,
        omega_nodes=frequencies,
        core_sigma_min=high.core_sigma_min_history[:rank],
        tie_counts=result.tie_counts[:rank],
        converged=converged,
        stop_reason=stop,
        residual_history=high.residual_history[: rank + 1],
        core_det_history=high.core_det_history[:rank],
        pivot_magnitudes=[value.lstrip("-") for value in pivots],
        pivot_geometry=_pivot_geometry(
            times, frequencies, str(case.cutoff), result.precision_bits, "1e-12"
        ),
        grid_shape=(len(case.t_grid), len(case.omega_grid)),
        grid_coordinate_precision_bits=53,
        grid_generation="float64_nested_dyadic_panel_lattice_residual_refinement",
        git_dirty=bool(provenance["git_dirty"]),
        environment=provenance["environment"],
        metadata=metadata,
    )
    return asdict(record)


def _run_case(
    cutoff: float,
    config: AdaptiveGridConfig,
    revision: str,
    config_hash: str,
    provenance: dict[str, Any],
) -> dict[str, Any]:
    print(f"adaptive census: cutoff={cutoff:g} started", file=sys.stderr, flush=True)
    case = run_adaptive_grid(cutoff, config=config)
    print(
        f"adaptive census: cutoff={cutoff:g} stages={len(case.stages)} "
        f"converged={case.converged} stop={case.stop_reason}",
        file=sys.stderr,
        flush=True,
    )
    return {
        "cutoff": str(cutoff),
        "stages": case.stages,
        "final_axes": {
            "t": [repr(float(value)) for value in case.t_grid],
            "omega": [repr(float(value)) for value in case.omega_grid],
        },
        "final_panels": {
            "t": [asdict(panel) for panel in case.t_panels],
            "omega": [asdict(panel) for panel in case.omega_panels],
        },
        "prefixes": case.prefixes,
        "records": [
            _record(case, prefix, revision, config_hash, provenance)
            for prefix in case.prefixes
        ],
        "converged": case.converged,
        "stop_reason": case.stop_reason,
    }


def _untracked_runtime_code(root: Path) -> bool:
    """Include experiment bootstrap code in both source transaction guards."""
    untracked = subprocess.check_output(
        [
            "git",
            "ls-files",
            "--others",
            "--exclude-standard",
            "-z",
            "--",
            "kernelgecp",
            "experiments",
        ],
        cwd=root,
        text=True,
    ).split("\0")
    return any(
        path
        and (
            path.startswith("kernelgecp/")
            or Path(path).suffix in (".py", ".so", ".pyd")
        )
        for path in untracked
    )


def run(config_path: Path, output_path: Path, *, root: Path = ROOT) -> None:
    """Atomically write deterministic JSON, retaining unsuccessful cases.

    Canonical generation uses a fresh CLI process from this checkout. Disk
    guards cannot attest previously imported code in a long-lived process.
    """
    raw, config = load_config(config_path)
    if raw.get("require_clean_revision", True) and root.resolve() != ROOT.resolve():
        raise ValueError("canonical source checkout must match the executing runner")
    revision, provenance = _git_revision(root), _provenance(root)
    source_hashes = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in SOURCES
    }
    lock_hash = hashlib.sha256((root / "uv.lock").read_bytes()).hexdigest()
    provenance["git_dirty"] = bool(
        provenance["git_dirty"] or _untracked_runtime_code(root)
    )
    if raw.get("require_clean_revision", True) and provenance["git_dirty"]:
        raise ValueError("canonical adaptive census requires clean tracked source")
    config_hash = _config_hash(raw)
    arguments = (
        [float(value) for value in raw["cutoffs"]],
        repeat(config),
        repeat(revision),
        repeat(config_hash),
        repeat(provenance),
    )
    if raw.get("workers", 1) == 1:
        cases = list(map(_run_case, *arguments))
    else:
        with ProcessPoolExecutor(max_workers=raw["workers"]) as executor:
            cases = list(executor.map(_run_case, *arguments))
    final_provenance = _provenance(root)
    if (
        _git_revision(root) != revision
        or hashlib.sha256((root / "uv.lock").read_bytes()).hexdigest() != lock_hash
        or any(
            hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
            for name, digest in source_hashes.items()
        )
        or final_provenance["environment"] != provenance["environment"]
    ):
        raise RuntimeError(
            "source, revision, lock or environment changed during census"
        )
    provenance["git_dirty"] = bool(
        provenance["git_dirty"]
        or final_provenance["git_dirty"]
        or _untracked_runtime_code(root)
    )
    if raw.get("require_clean_revision", True) and provenance["git_dirty"]:
        raise ValueError("canonical source became dirty during census")
    for case in cases:
        for record in case["records"]:
            record["git_dirty"] = provenance["git_dirty"]
    payload = {
        "schema_version": 3,
        "kind": "residual_adaptive_sampled_grid",
        "configuration": raw,
        "config_hash": config_hash,
        "git_commit": revision,
        **provenance,
        "lock_sha256": lock_hash,
        "algorithm_source_sha256": source_hashes,
        "certification_status": "not_certified",
        "error_scope": "finite_training_probe_and_independent_validation_sets_only",
        "cases": cases,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary.write_text(encoded + "\n", encoding="utf-8", newline="\n")
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("experiments/adaptive_census.toml")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("experiments/data/adaptive_census.json")
    )
    args = parser.parse_args()
    run(args.config, args.output)


if __name__ == "__main__":
    main()
