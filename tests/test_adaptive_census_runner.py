"""Explicit-axis census serialization and source-cleanliness failure paths."""

import json
from pathlib import Path
from typing import Any

import pytest

from experiments import adaptive_census as runner


def _config(tmp_path: Path, **changes: str) -> Path:
    fields = {
        "schema_version": "3",
        "cutoffs": "[1]",
        "tolerances": "[1e-5, 1e-8]",
        "precision_bits": "80",
        "check_precision_bits": "96",
        "generation_margin": "0.25",
        "max_rank": "24",
        "max_stages": "3",
        "max_nodes": "128",
        "max_marked": "2",
        "workers": "1",
        "require_clean_revision": "false",
    }
    fields.update(changes)
    path = tmp_path / "adaptive.toml"
    path.write_text("\n".join(f"{k} = {v}" for k, v in fields.items()) + "\n")
    return path


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("schema_version", "2"),
        ("schema_version", "true"),
        ("cutoffs", "[]"),
        ("cutoffs", "[1, 1.0]"),
        ("cutoffs", "[nan]"),
        ("cutoffs", "[true]"),
        ("cutoffs", "[0.5]"),
        ("workers", "0"),
        ("workers", "true"),
        ("require_clean_revision", "1"),
        ("tolerances", "1e-8"),
    ],
)
def test_invalid_census_configuration_fails_before_work(
    tmp_path: Path, field: str, value: str
) -> None:
    with pytest.raises(ValueError):
        runner.load_config(_config(tmp_path, **{field: value}))


@pytest.mark.parametrize(
    "unknown_code", ["", "kernelgecp/adaptive_grid.py", "experiments/__init__.py"]
)
def test_clean_source_gate_preserves_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, unknown_code: str
) -> None:
    config = _config(tmp_path, require_clean_revision="true")
    output = tmp_path / "previous.json"
    output.write_text("previous evidence\n")
    monkeypatch.setattr(runner, "_git_revision", lambda root: "a" * 40)
    monkeypatch.setattr(
        runner,
        "_provenance",
        lambda root: {"git_dirty": not unknown_code, "environment": {}},
    )
    monkeypatch.setattr(
        runner.subprocess,
        "check_output",
        lambda *a, **k: unknown_code,
    )

    def forbidden(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("work must not start on an unclean source")

    monkeypatch.setattr(runner, "_run_case", forbidden)
    with pytest.raises(ValueError, match="clean tracked source"):
        runner.run(config, output)
    assert output.read_text() == "previous evidence\n"
    assert not output.with_suffix(".json.tmp").exists()


def test_canonical_checkout_must_match_executing_source(tmp_path: Path) -> None:
    config = _config(tmp_path, require_clean_revision="true")
    with pytest.raises(ValueError, match="checkout must match"):
        runner.run(config, tmp_path / "rejected.json", root=tmp_path)


def test_untracked_generated_evidence_does_not_dirty_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path, require_clean_revision="true")
    monkeypatch.setattr(runner, "_git_revision", lambda root: "a" * 40)
    monkeypatch.setattr(
        runner, "_provenance", lambda root: {"git_dirty": False, "environment": {}}
    )
    monkeypatch.setattr(
        runner.subprocess,
        "check_output",
        lambda *a, **k: "experiments/data/adaptive_census.json\0",
    )
    monkeypatch.setattr(runner, "_run_case", lambda *a: {"pilot": "unused"})
    output = tmp_path / "repeat.json"
    runner.run(config, output)
    assert json.loads(output.read_text())["git_dirty"] is False


def test_small_explicit_axis_census_repeats_and_keeps_joint_prefix_metadata(
    tmp_path: Path,
) -> None:
    config = _config(tmp_path)
    first, second = tmp_path / "first.json", tmp_path / "second.json"
    runner.run(config, first)
    runner.run(config, second)
    assert first.read_bytes() == second.read_bytes()
    payload = json.loads(first.read_text())
    assert payload["schema_version"] == 3
    assert payload["kind"] == "residual_adaptive_sampled_grid"
    assert payload["certification_status"] == "not_certified"
    assert set(payload["algorithm_source_sha256"]) == set(runner.SOURCES)
    case = payload["cases"][0]
    assert len(case["stages"]) >= 2
    axes = case["final_axes"]
    assert case["stages"][-1]["t_grid"] == axes["t"]
    assert case["stages"][-1]["omega_grid"] == axes["omega"]
    assert len(case["records"]) == 2
    for record, prefix in zip(case["records"], case["prefixes"], strict=True):
        assert record["rank"] == prefix["chosen_rank"]
        assert record["tolerance"] == prefix["tolerance"]
        assert record["converged"] is prefix["converged"]
        assert record["metadata"]["converged"] is prefix["converged"]
        assert record["metadata"]["certification_status"] == "not_certified"
        assert (
            record["metadata"]["configuration"]["final_validation_used_for_refinement"]
            is False
        )
        assert record["grid_shape"] == [len(axes["t"]), len(axes["omega"])]
        assert record["algorithm"] == "residual_adaptive_grid"
        assert len(record["residual_history"]) == record["rank"] + 1
        assert all(float(x) in set(map(float, axes["t"])) for x in record["t_nodes"])
        assert all(
            float(x) in set(map(float, axes["omega"])) for x in record["omega_nodes"]
        )
        assert len(record["pivot_geometry"]) == record["rank"]


def test_stage_budget_failure_is_not_discarded(tmp_path: Path) -> None:
    config = _config(tmp_path, max_stages="1")
    output = tmp_path / "failed.json"
    runner.run(config, output)
    payload = json.loads(output.read_text())
    case = payload["cases"][0]
    assert not case["converged"]
    assert case["stop_reason"] == "max_stages"
    assert len(case["stages"]) == 1
    assert len(case["records"]) == 2
    assert all(not prefix["converged"] for prefix in case["prefixes"])
    assert all(not record["converged"] for record in case["records"])
    assert all(not record["metadata"]["converged"] for record in case["records"])
    assert payload["certification_status"] == "not_certified"


def test_mid_generation_source_change_preserves_previous_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "source"
    root.mkdir()
    for source in runner.SOURCES:
        path = root / source
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("source before generation\n")
    (root / "uv.lock").write_text("frozen lock\n")
    monkeypatch.setattr(runner, "_git_revision", lambda root: "a" * 40)
    monkeypatch.setattr(
        runner, "_provenance", lambda root: {"git_dirty": False, "environment": {}}
    )
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k: "")

    def changed(*args: Any, **kwargs: Any) -> dict[str, Any]:
        (root / runner.SOURCES[0]).write_text("source changed during generation\n")
        return {}

    monkeypatch.setattr(runner, "_run_case", changed)
    output = tmp_path / "previous.json"
    output.write_text("previous evidence\n")
    with pytest.raises(RuntimeError, match="changed during census"):
        runner.run(_config(tmp_path), output, root=root)
    assert output.read_text() == "previous evidence\n"
