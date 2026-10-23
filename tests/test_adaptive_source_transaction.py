"""Canonical evidence cannot replace prior output after new runtime code appears."""

import json
from pathlib import Path

import pytest

from experiments import adaptive_census as runner


@pytest.mark.parametrize("require_clean", [True, False])
def test_new_experiment_bootstrap_during_generation_is_not_clean_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, require_clean: bool
) -> None:
    config = tmp_path / "controls.toml"
    config.write_text(
        "schema_version = 3\ncutoffs = [1]\ntolerances = [1e-6]\n"
        f"workers = 1\nrequire_clean_revision = {str(require_clean).lower()}\n"
    )
    previous = tmp_path / "evidence.json"
    previous.write_text("previous evidence\n")
    bootstrap_created = False
    monkeypatch.setattr(runner, "_git_revision", lambda root: "a" * 40)
    monkeypatch.setattr(
        runner, "_provenance", lambda root: {"git_dirty": False, "environment": {}}
    )

    def untracked(*args: object, **kwargs: object) -> str:
        return "experiments/__init__.py\0" if bootstrap_created else ""

    def generate(*args: object) -> dict[str, object]:
        nonlocal bootstrap_created
        bootstrap_created = True
        return {
            "fixture": "no numerical work",
            "records": [{"git_dirty": False}, {"git_dirty": False}],
        }

    monkeypatch.setattr(runner.subprocess, "check_output", untracked)
    monkeypatch.setattr(runner, "_run_case", generate)
    if require_clean:
        with pytest.raises(ValueError, match="became dirty"):
            runner.run(config, previous)
        assert previous.read_text() == "previous evidence\n"
    else:
        runner.run(config, previous)
        payload = json.loads(previous.read_text())
        assert payload["git_dirty"] is True
        records = payload["cases"][0]["records"]
        assert len(records) == 2
        assert all(record["git_dirty"] is True for record in records)
    assert not previous.with_suffix(".json.tmp").exists()
