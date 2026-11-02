"""Metadata selection fixtures; these contain no scientific measurements."""

from copy import deepcopy

import pytest

from scripts.release_science_runs import REPOSITORY, ROLES, SCOPE, select_science_runs

SHA = "a" * 40


def metadata():
    runs, jobs, artifacts = {}, {}, {}
    for number, (role, (workflow, requirements, names)) in enumerate(ROLES.items(), 1):
        runs[role] = [
            {
                "id": number,
                "run_attempt": 1,
                "head_sha": SHA,
                "head_branch": "main",
                "event": "push",
                "path": f".github/workflows/{workflow}",
                "repository": {"full_name": REPOSITORY},
                "head_repository": {"full_name": REPOSITORY},
                "status": "completed",
                "conclusion": "success",
            }
        ]
        jobs[number] = [
            {
                "id": number * 100 + i,
                "run_id": number,
                "run_attempt": 1,
                "head_sha": SHA,
                "name": name,
                "status": "completed",
                "conclusion": "success",
                "steps": [
                    {"name": step, "status": "completed", "conclusion": "success"}
                    for step in steps
                ],
            }
            for i, (name, steps) in enumerate(requirements.items(), 1)
        ]
        artifacts[number] = [
            {
                "id": number * 100 + i,
                "name": name.format(id=number, attempt=1),
                "expired": False,
                "size_in_bytes": 20,
                "digest": "sha256:" + "b" * 64,
                "workflow_run": {
                    "id": number,
                    "head_sha": SHA,
                    "head_branch": "main",
                    "run_attempt": 1,
                },
            }
            for i, name in enumerate(names, 1)
        ]
    return runs, jobs, artifacts


def select(data):
    return select_science_runs(REPOSITORY, SHA, *data)


def test_selection_is_metadata_only():
    result = select(metadata())
    assert result["scope"] == SCOPE
    assert len(result["runs"]["providers"]["jobs"]) == 7


@pytest.mark.parametrize(
    "field,value",
    [
        ("head_sha", "b" * 40),
        ("head_branch", "other"),
        ("event", "pull_request"),
        ("repository", {"full_name": "wrong/repo"}),
        ("head_repository", {"full_name": "wrong/repo"}),
        ("path", "other"),
    ],
)
def test_wrong_run_binding(field, value):
    data = metadata()
    data[0]["molecular"][0][field] = value
    with pytest.raises(ValueError):
        select(data)


def test_failed_latest_is_not_bypassed():
    data = metadata()
    latest = deepcopy(data[0]["molecular"][0])
    latest.update(id=100, conclusion="failure")
    data[0]["molecular"].append(latest)
    with pytest.raises(ValueError):
        select(data)


@pytest.mark.parametrize(
    "role,index", [("molecular", 0), ("covariance", 1), ("epw", 2)]
)
def test_skipped_essential_step(role, index):
    data = metadata()
    run_id = data[0][role][0]["id"]
    data[1][run_id][0]["steps"][index]["conclusion"] = "skipped"
    with pytest.raises(ValueError):
        select(data)


def test_missing_provider_lane():
    data = metadata()
    data[1][1].pop()
    with pytest.raises(ValueError):
        select(data)


@pytest.mark.parametrize(
    "field,value",
    [
        ("digest", None),
        ("size_in_bytes", float("nan")),
        ("size_in_bytes", 0),
        ("expired", True),
        ("id", 0),
    ],
)
def test_invalid_artifact(field, value):
    data = metadata()
    data[2][2][0][field] = value
    with pytest.raises(ValueError):
        select(data)


@pytest.mark.parametrize(
    "field,value",
    [("id", 99), ("head_sha", "c" * 40), ("head_branch", "other"), ("run_attempt", 2)],
)
def test_invalid_artifact_binding(field, value):
    data = metadata()
    data[2][2][0]["workflow_run"][field] = value
    with pytest.raises(ValueError):
        select(data)


def test_stale_job_attempt():
    data = metadata()
    data[1][2][0]["run_attempt"] = 2
    with pytest.raises(ValueError):
        select(data)


def docker_record(data, role, suffix="JMEZHQ"):
    run_id = data[0][role][0]["id"]
    artifact = deepcopy(data[2][run_id][0])
    artifact["id"] += 10000
    artifact["name"] = f"tripp-smith~gecp-kernel-structure~{suffix}.dockerbuild"
    return artifact


@pytest.mark.parametrize("role", ["providers", "molecular", "epw"])
def test_auxiliary_build_record_is_retained(role):
    data = metadata()
    run_id = data[0][role][0]["id"]
    artifact = docker_record(data, role)
    data[2][run_id].append(artifact)
    result = select(data)["runs"][role]
    assert result["auxiliary_artifacts"][0]["id"] == artifact["id"]
    assert len(result["artifacts"]) == len(ROLES[role][2])


def test_two_provider_build_records_allowed():
    data = metadata()
    first = docker_record(data, "providers")
    second = docker_record(data, "providers", "IIR8U4")
    second["id"] += 1
    data[2][1].extend([first, second])
    assert len(select(data)["runs"]["providers"]["auxiliary_artifacts"]) == 2


@pytest.mark.parametrize(
    "case",
    [
        "wrongname",
        "excess",
        "badbinding",
        "duplicateid",
        "duplicatename",
        "bad_digest",
        "wrong_role",
    ],
)
def test_reject_bad_auxiliary(case):
    data = metadata()
    artifact = docker_record(data, "molecular")
    if case == "wrongname":
        artifact["name"] = "unrelated.dockerbuild"
    if case == "badbinding":
        artifact["workflow_run"]["head_sha"] = "c" * 40
    if case == "duplicateid":
        artifact["id"] = data[2][2][0]["id"]
    if case == "bad_digest":
        artifact["digest"] = None
    if case == "wrong_role":
        artifact = docker_record(data, "covariance")
        data[2][4].append(artifact)
    else:
        data[2][2].append(artifact)
    if case in ("excess", "duplicatename"):
        second = deepcopy(artifact)
        second["id"] += 1
        if case == "excess":
            second["name"] = "tripp-smith~gecp-kernel-structure~second.dockerbuild"
        data[2][2].append(second)
    with pytest.raises(ValueError):
        select(data)


def test_missing_provider_payload_audit_job():
    data = metadata()
    data[1][1] = [job for job in data[1][1] if job["name"] != "payload-audit"]
    with pytest.raises(ValueError, match="providers jobs incomplete"):
        select(data)


def test_skipped_provider_payload_audit_step():
    data = metadata()
    job = next(job for job in data[1][1] if job["name"] == "payload-audit")
    job["steps"][0]["conclusion"] = "skipped"
    with pytest.raises(ValueError, match="providers required step failed or skipped"):
        select(data)
