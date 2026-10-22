"""Deterministic command-line census for the fermionic GECP geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import tomllib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from decimal import Decimal
from importlib.metadata import version
from itertools import repeat
from pathlib import Path
from typing import Any

import mpmath as mp

from .gecp import gecp
from .grids import composite_chebyshev_grids
from .kernels import FermionicKernel
from .types import CensusRecord, GECPConfig, GECPResult, PivotGeometry


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
    """Load a finite, deterministic census configuration (legacy v1 accepted)."""

    with path.open("rb") as stream:
        raw: dict[str, Any] = tomllib.load(stream)
    if type(raw.get("schema_version")) is not int or raw["schema_version"] not in (
        1,
        2,
    ):
        raise ValueError("unsupported census schema_version")
    if not raw.get("cutoffs") or not raw.get("tolerances"):
        raise ValueError("census requires nonempty cutoffs and tolerances")
    for name, lower, upper in (("cutoffs", 1, math.inf), ("tolerances", 0, 1)):
        values = raw[name]
        if not isinstance(values, list) or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or (value < lower if name == "cutoffs" else value <= lower)
            or value >= upper
            for value in values
        ):
            raise ValueError(f"invalid finite {name}")
        if len(set(values)) != len(values):
            raise ValueError(f"duplicate {name}")
    for name, default, minimum in (
        ("max_rank", 256, 0),
        ("precision_bits", 128, 53),
        ("grid_order", 24, 2),
        ("workers", 1, 1),
        ("max_cells", 100_000, 1),
    ):
        value = raw.get(name, default)
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if raw.get("pivot", "grid") not in ("grid", "adaptive"):
        raise ValueError("unsupported pivot mode")
    for name in (
        "reflection_match_tolerance",
        "certificate_abs_tol",
        "certificate_rel_tol",
        "tie_rtol",
    ):
        value = raw.get(name, 0)
        if (
            isinstance(value, bool)
            or not math.isfinite(float(value))
            or float(value) < 0
        ):
            raise ValueError(f"{name} must be finite and nonnegative")
    if type(raw.get("require_clean_revision", False)) is not bool:
        raise ValueError("require_clean_revision must be boolean")
    if "tie_rtol" in raw:
        raw["tie_rtol"] = float(raw["tie_rtol"])
    return raw


def _provenance(root: Path) -> dict[str, Any]:
    dirty = subprocess.run(
        ["git", "diff", "HEAD", "--quiet"], cwd=root, check=False
    ).returncode
    if dirty not in (0, 1):
        raise RuntimeError("cannot inspect tracked revision cleanliness")
    # Unknown untracked runtime code matters; unrelated user notes do not.
    untracked = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard", "--", "kernelgecp"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {
        "git_dirty": bool(dirty or untracked),
        "environment": {
            "python": platform.python_version(),
            "mpmath_backend": str(mp.libmp.BACKEND),
            **{name: version(name) for name in ("numpy", "mpmath", "scipy")},
            "uv_lock_sha256": hashlib.sha256(
                (root / "uv.lock").read_bytes()
            ).hexdigest(),
        },
    }


def _pivot_geometry(
    t_nodes: list[str],
    omega_nodes: list[str],
    cutoff: str,
    precision_bits: int,
    match_tolerance: str,
) -> list[PivotGeometry]:
    """Recompute nearest partners within this prefix, using a stated metric."""
    if len(t_nodes) != len(omega_nodes):
        raise ValueError("node histories must align")
    with mp.workprec(precision_bits):
        limit, tolerance = mp.mpf(cutoff), mp.mpf(match_tolerance)
        if (
            not mp.isfinite(limit)
            or limit < 1
            or not mp.isfinite(tolerance)
            or tolerance < 0
        ):
            raise ValueError("invalid geometry controls")
        nodes = [
            (mp.mpf(t), mp.mpf(w)) for t, w in zip(t_nodes, omega_nodes, strict=True)
        ]
        if any(
            not (mp.isfinite(t) and mp.isfinite(w) and 0 <= t <= 1 and abs(w) <= limit)
            for t, w in nodes
        ):
            raise ValueError("nodes must be finite and in the physical domain")
        digits = math.ceil(precision_bits * math.log10(2)) + 3

        def encode(value: Any) -> str:
            return str(mp.nstr(value, digits))

        geometry = []
        for t, omega in nodes:
            scale = 0
            while abs(omega) > mp.mpf(2) ** scale:
                scale += 1
            reflected_t, reflected_omega = 1 - t, -omega
            distances = [
                max(abs(s - reflected_t), abs(w - reflected_omega) / limit)
                for s, w in nodes
            ]
            index = min(range(len(nodes)), key=distances.__getitem__)
            geometry.append(
                PivotGeometry(
                    t_boundary_distances=(encode(t), encode(1 - t)),
                    omega_boundary_distances=(
                        encode(omega + limit),
                        encode(limit - omega),
                    ),
                    dyadic_scale=scale,
                    dyadic_panel="central"
                    if scale == 0
                    else "positive"
                    if omega > 0
                    else "negative",
                    reflection_t=encode(reflected_t),
                    reflection_omega=encode(reflected_omega),
                    reflection_nearest_index=index,
                    reflection_distance=encode(distances[index]),
                    reflection_matched=bool(distances[index] <= tolerance),
                )
            )
    return geometry


def _records_from_result(
    raw: dict[str, Any],
    cutoff: float,
    revision: str,
    config_hash: str,
    result: GECPResult,
    provenance: dict[str, Any],
) -> list[CensusRecord]:
    """Extract aligned prefixes without promoting an uncertified bound to evidence."""
    high = result.high_precision
    histories = {
        name: list(getattr(high, name))
        if high is not None
        else [_decimal(value) for value in getattr(result, name)]
        for name in (
            "pivots",
            "t_nodes",
            "omega_nodes",
            "residual_history",
            "core_det_history",
            "core_sigma_min_history",
        )
    }
    for name, values in histories.items():
        expected = result.rank + (name == "residual_history")
        if len(values) != expected or any(
            not Decimal(value).is_finite() for value in values
        ):
            raise ValueError(f"invalid or misaligned {name}")
    if any(Decimal(value) < 0 for value in histories["residual_history"]):
        raise ValueError("residual bounds must be nonnegative")
    cutoff_text = str(cutoff)
    reflection_tolerance = str(raw.get("reflection_match_tolerance", 1e-12))
    grids = composite_chebyshev_grids(cutoff, int(raw.get("grid_order", 24)))
    source_metadata = result.metadata
    adaptive = source_metadata.convergence_scope == "continuous_interval_bound"
    records = []
    for tolerance in raw["tolerances"]:
        threshold = Decimal(str(tolerance))

        def certified_at(index: int) -> bool:
            if not adaptive:
                return True
            if index < len(result.pivot_certificates):
                return result.pivot_certificates[index].certified
            return result.converged and result.stop_reason == "certified_tolerance"

        rank = next(
            (
                i
                for i, value in enumerate(histories["residual_history"])
                if Decimal(value) <= threshold and certified_at(i)
            ),
            result.rank,
        )
        residual_history = histories["residual_history"][: rank + 1]
        residual = residual_history[-1]
        converged = Decimal(residual) <= threshold and certified_at(rank)
        stop = (
            ("certified_tolerance" if adaptive else "tolerance")
            if converged
            else result.stop_reason
        )
        conditioning = source_metadata.conditioning_history[:rank]
        metadata = replace(
            source_metadata,
            conditioning_history=conditioning,
            conditioning_status=(
                "not_representable"
                if any(item.status == "not_representable" for item in conditioning)
                else "estimated"
                if conditioning
                else "not_applicable"
            ),
            tie_counts=source_metadata.tie_counts[:rank],
            converged=bool(converged),
            stop_reason=stop,
            certification_status="certified"
            if adaptive and converged
            else "not_certified",
            configuration={
                **result.config,
                "tol": float(tolerance),
                "source_trajectory_tolerance": str(min(raw["tolerances"])),
                "prefix_extraction": True,
            },
        )
        pivots = histories["pivots"][:rank]
        t_nodes, omega_nodes = (
            histories["t_nodes"][:rank],
            histories["omega_nodes"][:rank],
        )
        records.append(
            CensusRecord(
                schema_version=2,
                git_commit=revision,
                config_hash=config_hash,
                cutoff=cutoff_text,
                tolerance=str(tolerance),
                precision_bits=result.precision_bits,
                algorithm=str(raw.get("pivot", "grid")),
                rank=rank,
                residual=residual,
                pivots=pivots,
                t_nodes=t_nodes,
                omega_nodes=omega_nodes,
                core_sigma_min=histories["core_sigma_min_history"][:rank],
                tie_counts=[] if adaptive else result.tie_counts[:rank],
                converged=bool(converged),
                stop_reason=stop,
                residual_history=residual_history,
                core_det_history=histories["core_det_history"][:rank],
                pivot_magnitudes=[value.lstrip("-") for value in pivots],
                pivot_geometry=_pivot_geometry(
                    t_nodes,
                    omega_nodes,
                    cutoff_text,
                    result.precision_bits,
                    reflection_tolerance,
                ),
                grid_shape=(0, 0) if adaptive else (len(grids[0]), len(grids[1])),
                grid_coordinate_precision_bits=None if adaptive else 53,
                grid_generation="not_applicable_interval_adaptive"
                if adaptive
                else "float64_composite_chebyshev_lobatto_dyadic_endpoint_panels",
                reflection_match_tolerance=reflection_tolerance,
                git_dirty=bool(provenance["git_dirty"]),
                environment=provenance["environment"],
                metadata=metadata,
            )
        )
    return records


def _run_cutoff(
    raw: dict[str, Any],
    cutoff: float,
    revision: str,
    config_hash: str,
    provenance: dict[str, Any],
) -> list[CensusRecord]:
    """Run one strict trajectory and extract every requested tolerance prefix."""

    tolerances = [float(value) for value in raw["tolerances"]]
    config = GECPConfig(
        tol=min(tolerances),
        max_rank=int(raw.get("max_rank", 256)),
        pivot=str(raw.get("pivot", "grid")),  # type: ignore[arg-type]
        precision_bits=int(raw.get("precision_bits", 128)),
        grid_order=int(raw.get("grid_order", 24)),
        certificate_abs_tol=float(raw.get("certificate_abs_tol", 1e-12)),
        certificate_rel_tol=float(raw.get("certificate_rel_tol", 1e-8)),
        max_cells=int(raw.get("max_cells", 100_000)),
        tie_rtol=raw.get("tie_rtol"),
    )
    result = gecp(FermionicKernel(float(cutoff)), config=config)
    return _records_from_result(raw, cutoff, revision, config_hash, result, provenance)


def run_census(
    config_path: Path, output_path: Path, *, root: Path | None = None
) -> None:
    """Run every configured case and atomically replace canonical JSONL output."""

    project_root = config_path.resolve().parents[1] if root is None else root.resolve()
    raw = load_census_config(config_path)
    revision = _git_revision(project_root)
    provenance = _provenance(project_root)
    if raw.get("require_clean_revision", False) and provenance["git_dirty"]:
        raise ValueError("canonical census requires a clean tracked revision")
    config_hash = _config_hash(raw)
    workers = int(raw.get("workers", 1))
    if workers < 1:
        raise ValueError("census workers must be positive")
    cutoffs = [float(value) for value in raw["cutoffs"]]
    if workers == 1:
        grouped_records = [
            _run_cutoff(raw, cutoff, revision, config_hash, provenance)
            for cutoff in cutoffs
        ]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            grouped_records = list(
                executor.map(
                    _run_cutoff,
                    repeat(raw),
                    cutoffs,
                    repeat(revision),
                    repeat(config_hash),
                    repeat(provenance),
                )
            )
    records = [record for group in grouped_records for record in group]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        for record in records:
            stream.write(
                json.dumps(
                    asdict(record),
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
            )
            stream.write("\n")
    temporary.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("experiments/census.toml"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("experiments/data/gecp_census_v2.jsonl"),
    )
    arguments = parser.parse_args()
    run_census(arguments.config, arguments.output)


if __name__ == "__main__":
    main()
