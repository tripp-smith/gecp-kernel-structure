import runpy
from copy import deepcopy
from pathlib import Path

import pytest

from scripts import release_science_payloads as payloads
from scripts.release_science_fetch import fetch_spectroscopy_origin
from scripts.release_science_runs import (
    REPOSITORY,
    V4_JOBS,
    select_science_runs,
    select_spectroscopy_origin,
)

SHA = "a" * 40
RUN_ID = 123


def origin():
    run = dict(
        id=RUN_ID,
        run_attempt=1,
        head_sha=SHA,
        head_branch="main",
        event="workflow_dispatch",
        path=".github/workflows/maxent-science.yml",
        repository={"full_name": REPOSITORY},
        head_repository={"full_name": REPOSITORY},
        status="completed",
        conclusion="success",
    )
    jobs = [
        dict(
            id=10 + i,
            run_id=RUN_ID,
            run_attempt=1,
            head_sha=SHA,
            name=name,
            status="completed",
            conclusion="success",
            steps=[
                dict(name=step, status="completed", conclusion="success")
                for step in steps
            ],
        )
        for i, (name, steps) in enumerate(V4_JOBS.items())
    ]
    artifacts = [
        dict(
            id=555,
            name=f"maxent-fixed-bath-final-{RUN_ID}-1",
            expired=False,
            size_in_bytes=20,
            digest="sha256:" + "b" * 64,
            workflow_run=dict(
                id=RUN_ID, head_sha=SHA, head_branch="main", run_attempt=1
            ),
        )
    ]
    return run, jobs, artifacts


def select(data):
    return select_spectroscopy_origin(REPOSITORY, SHA, RUN_ID, *data)


def test_origin_metadata_is_not_scientific_acceptance():
    result = select(origin())
    assert result["artifact"]["id"] == 555
    assert "replay required" in result["scope"]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d[0].update(run_attempt=2),
        lambda d: d[0].update(conclusion="failure"),
        lambda d: d[0].update(head_sha="c" * 40),
        lambda d: d[0].update(event="push"),
        lambda d: d[1].pop(),
        lambda d: d[1][1]["steps"][1].update(conclusion="failure"),
        lambda d: d[2].clear(),
        lambda d: d[2][0].update(expired=True),
        lambda d: d[2][0].update(digest="sha256:bad"),
        lambda d: d[2][0]["workflow_run"].update(head_sha="c" * 40),
        lambda d: d[2].append(deepcopy(d[2][0])),
    ],
)
def test_origin_rejects_incomplete_or_mismatched(mutate):
    data = origin()
    mutate(data)
    with pytest.raises(ValueError):
        select(data)


def test_fetch_authenticates_current_attempt_and_pages():
    run, jobs, artifacts = origin()

    def get(path, page):
        if path.endswith("/runs/123"):
            return run
        key, values = (
            ("jobs", jobs) if path.endswith("/jobs") else ("artifacts", artifacts)
        )
        return {"total_count": len(values), key: values}

    assert (
        fetch_spectroscopy_origin(REPOSITORY, SHA, RUN_ID, get)["artifact"]["id"] == 555
    )
    run["run_attempt"] = 2
    with pytest.raises(ValueError, match="attempt 1"):
        fetch_spectroscopy_origin(REPOSITORY, SHA, RUN_ID, get)


def test_production_profile_exact_roles_and_legacy_shape(tmp_path):
    fixtures = runpy.run_path(
        str(Path(__file__).with_name("test_release_science_runs.py"))
    )
    data = fixtures["metadata"]()
    assert "profile" not in select_science_runs(REPOSITORY, SHA, *data)
    data[0].pop("spectroscopy")
    result = select_science_runs(
        REPOSITORY, SHA, *data, profile="production-components"
    )
    assert result["profile"] == "production-components"
    assert set(result["runs"]) == {
        "providers",
        "molecular",
        "epw",
        "covariance",
        "cthyb",
    }
    data[0]["spectroscopy"] = []
    with pytest.raises(ValueError, match="inventory"):
        select_science_runs(REPOSITORY, SHA, *data, profile="production-components")
    fixture = runpy.run_path(
        str(Path(__file__).with_name("test_release_science_payloads.py"))
    )
    blob = fixture["zip_bytes"]([("one", "data")])
    selected = fixture["selection"](blob)
    selected["profile"] = "production-components"
    with pytest.raises(ValueError, match="roles"):
        payloads.validate_selection(selected)
    selected["runs"].pop("spectroscopy")
    out = payloads.materialize(
        selected, tmp_path / "out", tmp_path / "receipt", fixture["fake"](blob)
    )
    assert out["profile"] == "production-components"
    assert set(out["auxiliary_artifacts"]) == set(selected["runs"])
