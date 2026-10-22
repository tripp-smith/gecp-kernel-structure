import json
from pathlib import Path
from typing import Any

import pytest

from kernelgecp import census
from kernelgecp.types import CensusRecord


def _config(tmp_path: Path, **changes: str) -> Path:
    fields = {
        "schema_version": "2",
        "cutoffs": "[1, 4]",
        "tolerances": "[1e-6, 1e-12]",
        "precision_bits": "128",
        "max_rank": "12",
        "grid_order": "6",
        "workers": "1",
        "pivot": '"grid"',
    }
    fields.update(changes)
    path = tmp_path / "census.toml"
    path.write_text(
        "\n".join(f"{name} = {value}" for name, value in fields.items()) + "\n",
        encoding="utf-8",
    )
    return path


@pytest.mark.parametrize("field", ["cutoffs", "tolerances"])
@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "true"])
def test_census_rejects_nonfinite_or_boolean_values(
    tmp_path: Path, field: str, value: str
) -> None:
    with pytest.raises(ValueError, match=field):
        census.load_census_config(_config(tmp_path, **{field: f"[{value}]"}))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("cutoffs", "[0.5]"),
        ("tolerances", "[0]"),
        ("tolerances", "[1]"),
        ("cutoffs", "[]"),
        ("tolerances", "[]"),
        ("cutoffs", "1"),
        ("tolerances", '"1e-6"'),
    ],
)
def test_census_rejects_invalid_ranges_and_container_types(
    tmp_path: Path, field: str, value: str
) -> None:
    with pytest.raises(ValueError):
        census.load_census_config(_config(tmp_path, **{field: value}))


@pytest.mark.parametrize(
    ("field", "minimum"),
    [
        ("workers", 1),
        ("grid_order", 2),
        ("max_rank", 0),
        ("precision_bits", 53),
        ("max_cells", 1),
    ],
)
@pytest.mark.parametrize("value", ["true", "1.5", "nan", "below_minimum"])
def test_census_integer_controls_are_strict(
    tmp_path: Path, field: str, minimum: int, value: str
) -> None:
    encoded = str(minimum - 1) if value == "below_minimum" else value
    with pytest.raises(ValueError, match=field):
        census.load_census_config(_config(tmp_path, **{field: encoded}))


@pytest.mark.parametrize(
    ("field", "value"),
    [("cutoffs", "[1, 1.0]"), ("tolerances", "[1e-6, 0.000001]")],
)
def test_census_rejects_numeric_duplicates(
    tmp_path: Path, field: str, value: str
) -> None:
    with pytest.raises(ValueError, match=f"duplicate {field}"):
        census.load_census_config(_config(tmp_path, **{field: value}))


@pytest.mark.parametrize("value", ['"complete"', '"gridd"', "true", "1"])
def test_census_rejects_unknown_pivot_modes(tmp_path: Path, value: str) -> None:
    with pytest.raises(ValueError, match="pivot"):
        census.load_census_config(_config(tmp_path, pivot=value))


@pytest.mark.parametrize(
    "field",
    [
        "reflection_match_tolerance",
        "certificate_abs_tol",
        "certificate_rel_tol",
        "tie_rtol",
    ],
)
@pytest.mark.parametrize("value", ["nan", "inf", "-1e-6", "true"])
def test_census_rejects_invalid_diagnostic_tolerances(
    tmp_path: Path, field: str, value: str
) -> None:
    with pytest.raises(ValueError, match=field):
        census.load_census_config(_config(tmp_path, **{field: value}))


@pytest.mark.parametrize("schema", ["1", "2"])
def test_census_schema_one_and_two_configs_remain_supported(
    tmp_path: Path, schema: str
) -> None:
    raw = census.load_census_config(
        _config(
            tmp_path,
            schema_version=schema,
            max_rank="0",
            reflection_match_tolerance="0",
            tie_rtol="0",
        )
    )
    assert raw["schema_version"] == int(schema)
    assert raw["cutoffs"] == [1, 4]
    assert raw["max_rank"] == 0


@pytest.mark.parametrize("value", ["0", '"false"'])
def test_require_clean_revision_must_be_boolean(tmp_path: Path, value: str) -> None:
    with pytest.raises(ValueError, match="require_clean_revision"):
        census.load_census_config(_config(tmp_path, require_clean_revision=value))


def _record(cutoff: float = 1) -> CensusRecord:
    return CensusRecord(
        schema_version=2,
        git_commit="a" * 40,
        config_hash="b" * 64,
        cutoff=str(cutoff),
        tolerance="1e-6",
        precision_bits=128,
        algorithm="grid",
        rank=0,
        residual="1",
        pivots=[],
        t_nodes=[],
        omega_nodes=[],
        core_sigma_min=[],
        tie_counts=[],
        converged=False,
        stop_reason="max_rank",
        residual_history=["1"],
    )


def _mock_provenance(monkeypatch: pytest.MonkeyPatch, *, dirty: bool) -> None:
    monkeypatch.setattr(census, "_git_revision", lambda root: "a" * 40)
    monkeypatch.setattr(
        census,
        "_provenance",
        lambda root: {"git_dirty": dirty, "environment": {"python": "test-runtime"}},
    )


@pytest.mark.parametrize("existing_output", [False, True])
def test_dirty_revision_gate_prevents_work_or_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, existing_output: bool
) -> None:
    config = _config(tmp_path, require_clean_revision="true")
    output = tmp_path / "new-output" / "census.jsonl"
    if existing_output:
        output.parent.mkdir()
        output.write_bytes(b"existing canonical evidence\n")
    _mock_provenance(monkeypatch, dirty=True)

    def forbidden_run(*args: Any) -> list[CensusRecord]:
        pytest.fail("dirty canonical census must not start a trajectory")

    monkeypatch.setattr(census, "_run_cutoff", forbidden_run)
    with pytest.raises(ValueError, match="clean tracked revision"):
        census.run_census(config, output, root=tmp_path)
    if existing_output:
        assert output.read_bytes() == b"existing canonical evidence\n"
    else:
        assert not output.parent.exists()
    assert not output.with_suffix(".jsonl.tmp").exists()


def test_failed_trajectory_does_not_replace_existing_census(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path)
    output = tmp_path / "census.jsonl"
    output.write_bytes(b"existing canonical evidence\n")
    _mock_provenance(monkeypatch, dirty=False)
    cutoffs_seen = []

    def fail_second(
        raw: dict[str, Any],
        cutoff: float,
        revision: str,
        config_hash: str,
        provenance: dict[str, Any],
    ) -> list[CensusRecord]:
        cutoffs_seen.append(cutoff)
        assert revision == "a" * 40
        assert len(config_hash) == 64
        assert provenance["git_dirty"] is False
        if cutoff == 4:
            raise RuntimeError("injected trajectory failure")
        return [_record(cutoff)]

    monkeypatch.setattr(census, "_run_cutoff", fail_second)
    with pytest.raises(RuntimeError, match="trajectory failure"):
        census.run_census(config, output, root=tmp_path)
    assert cutoffs_seen == [1, 4]
    assert output.read_bytes() == b"existing canonical evidence\n"
    assert not output.with_suffix(".jsonl.tmp").exists()


def test_mid_serialization_failure_preserves_existing_census(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path, cutoffs="[1]")
    output = tmp_path / "census.jsonl"
    output.write_bytes(b"existing canonical evidence\n")
    _mock_provenance(monkeypatch, dirty=False)
    valid = _record()
    invalid = _record()
    invalid.environment["nonfinite"] = float("nan")  # type: ignore[assignment]
    monkeypatch.setattr(
        census,
        "_run_cutoff",
        lambda raw, cutoff, revision, digest, source: [valid, invalid],
    )
    with pytest.raises(ValueError, match="Out of range float"):
        census.run_census(config, output, root=tmp_path)
    assert output.read_bytes() == b"existing canonical evidence\n"


def test_clean_revision_allows_atomic_schema_two_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path, schema_version="1", require_clean_revision="true")
    output = tmp_path / "census.jsonl"
    _mock_provenance(monkeypatch, dirty=False)
    monkeypatch.setattr(
        census,
        "_run_cutoff",
        lambda raw, cutoff, revision, digest, source: [_record(cutoff)],
    )
    census.run_census(config, output, root=tmp_path)
    records = [json.loads(line) for line in output.read_text().splitlines()]
    assert [record["cutoff"] for record in records] == ["1.0", "4.0"]
    assert all(record["schema_version"] == 2 for record in records)
    assert not output.with_suffix(".jsonl.tmp").exists()
