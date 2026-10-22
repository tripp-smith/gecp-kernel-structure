"""Physical endpoint and historical-grid preservation regressions."""

from __future__ import annotations

import json
import math
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

from kernelgecp.census import run_census
from kernelgecp.grids import (
    chebyshev_lobatto,
    composite_chebyshev_grids,
    frequency_grid,
    time_grid,
)
from kernelgecp.types import FloatArray


def _legacy_frequency_grid(cutoff: float, order: int) -> FloatArray:
    """Reproduce the historical arithmetic, without the endpoint repair."""
    parts = [chebyshev_lobatto(0, 1, order)]
    left = 1.0
    for _ in range(max(0, math.ceil(math.log2(cutoff)))):
        right = min(cutoff, 2 * left)
        if right > left:
            parts.append(chebyshev_lobatto(left, right, order))
        left = right
        if left >= cutoff:
            break
    positive = np.unique(np.concatenate(parts))
    return np.concatenate((-positive[:0:-1], positive))


@pytest.mark.parametrize("cutoff", [1, 10, 100, 1000, 10000, 100000, 1000000])
def test_canonical_integer_grids_preserve_every_historical_byte(cutoff: float) -> None:
    times, frequencies = composite_chebyshev_grids(cutoff, 24)
    assert times.tobytes() == time_grid(24, cutoff).tobytes()
    assert frequencies.tobytes() == _legacy_frequency_grid(cutoff, 24).tobytes()


@pytest.mark.parametrize("cutoff", [1.2, 1.7, 1000000.1])
@pytest.mark.parametrize("order", [2, 6, 24])
def test_noninteger_frequency_endpoints_equal_physical_cutoff(
    cutoff: float, order: int
) -> None:
    frequencies = frequency_grid(cutoff, order)
    legacy = _legacy_frequency_grid(cutoff, order)
    assert frequencies[0] == -cutoff
    assert frequencies[-1] == cutoff
    assert np.all(np.abs(frequencies) <= cutoff)
    assert np.array_equal(frequencies[1:-1], legacy[1:-1])
    assert np.array_equal(frequencies, -frequencies[::-1])
    assert np.all(np.diff(frequencies) > 0)


def test_noninteger_census_records_remain_inside_physical_domain(
    tmp_path: Path,
) -> None:
    config = tmp_path / "census.toml"
    config.write_text(
        """schema_version = 2
cutoffs = [1.2, 1.7]
tolerances = [1e-6]
precision_bits = 128
max_rank = 8
grid_order = 6
pivot = "grid"
require_clean_revision = false
""",
        encoding="utf-8",
    )
    output = tmp_path / "census.jsonl"
    run_census(config, output, root=Path(__file__).parents[1])
    records = [json.loads(line) for line in output.read_text().splitlines()]
    assert len(records) == 2
    for record in records:
        cutoff = float(record["cutoff"])
        assert record["converged"]
        # Grid coordinates originate in float64; their 128-bit decimal
        # encodings may differ from the exact short decimal by MP roundoff.
        assert all(abs(float(value)) <= cutoff for value in record["omega_nodes"])
        assert float(record["omega_nodes"][0]) in (-cutoff, cutoff)
        assert all(
            Decimal(distance) >= 0
            for geometry in record["pivot_geometry"]
            for distance in geometry["omega_boundary_distances"]
        )
