"""Cheap cross-platform checks of the stored, authenticated numeric outputs."""

import gzip
import json
from pathlib import Path

from scripts.audit_precision_bridge import validate_record


def test_canonical_v2_records_validate_without_replaying_platform_blas() -> None:
    root = Path(__file__).resolve().parents[1]
    artifact = json.loads(
        gzip.decompress(
            (root / "experiments/data/precision_bridge.json.gz").read_bytes()
        )
    )
    assert artifact["schema_version"] == 2
    assert len(artifact["records"]) == len(artifact["config"]["cases"]) == 6
    enclosed = 0
    for record, config in zip(
        artifact["records"], artifact["config"]["cases"], strict=True
    ):
        arrays, grid, _ = validate_record(record, config)
        assert len(arrays) == 15 and grid[2].shape == (5, 5)
        enclosed += record["fit_certificate"]["certified"]
    assert enclosed == 5
    assert artifact["records"][2]["fit_certificate"]["termination_reason"] == (
        "normal_pivot_contains_zero"
    )
