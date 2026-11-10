"""Offline tests for bounded authenticated metadata collection."""

import runpy
from pathlib import Path
from unittest.mock import patch

import pytest

from scripts.release_science_fetch import _gh, _pages, fetch_science_runs
from scripts.release_science_runs import REPOSITORY, ROLES

fixture = runpy.run_path(str(Path(__file__).with_name("test_release_science_runs.py")))
metadata, SHA = fixture["metadata"], fixture["SHA"]


def adapter(data, calls):
    runs, jobs, artifacts = data

    def get(path, page):
        calls.append((path, page))
        parts = path.split("/")
        if "workflows" in parts:
            workflow = parts[parts.index("workflows") + 1]
            role = next(r for r, v in ROLES.items() if v[0] == workflow)
            key, values = "workflow_runs", runs[role]
        else:
            run_id = int(parts[parts.index("runs") + 1])
            key = "jobs" if parts[-1] == "jobs" else "artifacts"
            values = (jobs if key == "jobs" else artifacts)[run_id]
        return {"total_count": len(values), key: values[(page - 1) * 100 : page * 100]}

    return get


def test_full_collection_and_pagination():
    data = metadata()
    template = data[0]["molecular"][0]
    # Historical unrelated runs force a second page without changing selection.
    data[0]["molecular"] += [
        dict(template, id=1000 + i, head_sha="b" * 40) for i in range(101)
    ]
    calls = []
    result = fetch_science_runs(REPOSITORY, SHA, adapter(data, calls))
    assert result["runs"]["molecular"]["id"] == 2
    assert result["provenance"]["payloads_downloaded"] is False
    assert any("molecular-green.yml/runs" in path and page == 2 for path, page in calls)
    assert any("/runs/2/attempts/1/jobs" in path for path, _ in calls)
    workflow_paths = [path for path, _ in calls if "/workflows/" in path]
    assert all(f"?branch=main&head_sha={SHA}" in path for path in workflow_paths)


def test_late_failed_run_is_not_bypassed():
    data = metadata()
    old = data[0]["molecular"][0]
    data[0]["molecular"] += [
        dict(old, id=1000 + i, head_sha="b" * 40) for i in range(100)
    ]
    data[0]["molecular"].append(dict(old, id=3000, conclusion="failure"))
    with pytest.raises(ValueError, match="latest molecular run failed"):
        fetch_science_runs(REPOSITORY, SHA, adapter(data, []))


def test_stale_attempt_rejected():
    data = metadata()
    data[1][2][0]["run_attempt"] = 2
    with pytest.raises(ValueError, match="invalid molecular job"):
        fetch_science_runs(REPOSITORY, SHA, adapter(data, []))


@pytest.mark.parametrize(
    "repository,revision", [("wrong/repo", SHA), (REPOSITORY, "A" * 40)]
)
def test_binding_rejected_before_api(repository, revision):
    with pytest.raises(ValueError):
        fetch_science_runs(repository, revision, lambda *_: pytest.fail("API called"))


@pytest.mark.parametrize(
    "response",
    [
        None,
        {},
        {"total_count": True, "jobs": []},
        {"total_count": 2, "jobs": [{"id": 1}]},
        {"total_count": 2, "jobs": [{"id": 1}, {"id": 1}]},
        {"total_count": 10001, "jobs": []},
    ],
)
def test_bad_pages(response):
    with pytest.raises(ValueError):
        _pages("path", "jobs", lambda *_: response)


def test_changed_total_rejected():
    def get(_, page):
        return {
            "total_count": 101 if page == 1 else 102,
            "jobs": [{"id": i} for i in range(1, 101)] if page == 1 else [{"id": 101}],
        }

    with pytest.raises(ValueError, match="total changed"):
        _pages("path", "jobs", get)


def test_gh_errors_and_get_arguments():
    class Result:
        returncode = 1
        stdout = ""

    with patch(
        "scripts.release_science_fetch.subprocess.run", return_value=Result()
    ) as run:
        with pytest.raises(ValueError, match="request failed"):
            _gh("endpoint", 2)
        assert run.call_args.args[0] == [
            "gh",
            "api",
            "--method",
            "GET",
            "endpoint",
            "-f",
            "per_page=100",
            "-f",
            "page=2",
        ]
    Result.returncode, Result.stdout = 0, "not json"
    with (
        patch("scripts.release_science_fetch.subprocess.run", return_value=Result()),
        pytest.raises(ValueError, match="JSON"),
    ):
        _gh("endpoint", 1)
