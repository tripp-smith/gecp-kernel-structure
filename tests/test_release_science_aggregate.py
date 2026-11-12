import json
from pathlib import Path

import pytest
from test_release_science_gate import fixtures

from scripts import release_science_aggregate as gate
from scripts import release_science_gate as roles

P = "1" * 40
K = "2" * 40


@pytest.fixture
def replays(tmp_path):
    _, selection, origin, _, _ = fixtures(tmp_path)
    selection["revision"] = P
    plan_dir = tmp_path / "plan"
    replay_dir = tmp_path / "replays"
    plan_dir.mkdir()
    replay_dir.mkdir()
    selected_path = plan_dir / "selected-runs.json"
    origin_path = plan_dir / "spectroscopy-origin.json"
    gate.write(selected_path, selection)
    gate.write(origin_path, origin)
    plan = {
        "schema": gate.PLAN_SCHEMA,
        "repository": gate.REPOSITORY,
        "run_id": "500",
        "run_attempt": "1",
        "producer_revision": P,
        "candidate_revision": K,
        "science_revision": origin["revision"],
        "mode": "tag",
        "selection_sha256": gate.digest(selected_path),
        "origin_sha256": gate.digest(origin_path),
    }
    gate.write(plan_dir / "plan.json", plan)
    for role in gate.ROLES:
        folder = replay_dir / role
        folder.mkdir()
        source = origin["revision"] if role == "spectroscopy" else P
        wrapper = {key: {"present": True} for key in roles.FIELDS[role]}
        wrapper.update(
            schema=roles.SCHEMAS[role],
            source_revision=source,
            status="passed",
            payload_unchanged=True,
            inputs_unchanged=True,
        )
        if role == "providers":
            from scripts.release_provider_audits import LANES

            wrapper["lanes"] = {
                name: {
                    "receipt_valid": True,
                    "receipt_sha256": "a" * 64,
                    "bound_payload_sha256": {},
                    "independent_dmft": "audit.json",
                }
                for name in LANES
            }
        if role == "spectroscopy":
            wrapper.update(
                counts={
                    "synthetic_calibration": 24,
                    "synthetic_final": 72,
                    "measured_calibration": 8,
                    "measured_final": 24,
                    "raw_chains": 256,
                },
                receipt_sha256={
                    name: "b" * 64
                    for name in (
                        "synthetic-cal",
                        "synthetic-final",
                        "measured-cal",
                        "raw",
                        "measured-final",
                    )
                },
                archive_sha256=origin["artifact"]["digest"][7:],
                archive_size=origin["artifact"]["size_in_bytes"],
                run_id=str(origin["id"]),
            )
        wrapper_path = folder / "wrapper-receipt.json"
        gate.write(wrapper_path, wrapper)
        gate.write(
            folder / "role-receipt.json",
            {
                "schema": "kernelgecp.science-role-replay.v1",
                "role": role,
                "status": "passed",
                "source_revision": source,
                "producer_revision": P,
                "selection_sha256": plan["selection_sha256"],
                "origin_sha256": plan["origin_sha256"],
                "wrapper_receipt_sha256": gate.digest(wrapper_path),
                "command_exit": 0,
                "artifacts": [origin["artifact"]]
                if role == "spectroscopy"
                else selection["runs"][role]["artifacts"],
            },
        )
    return plan_dir, replay_dir


def aggregate(replays, out, **fields):
    kwargs = dict(run_id="500", run_attempt="1", candidate=K)
    kwargs.update(fields)
    gate.aggregate(*replays, out, **kwargs)


def test_accept_six_bound_receipts(tmp_path, replays):
    out = tmp_path / "accepted.json"
    aggregate(replays, out)
    value = json.loads(out.read_text())
    assert value["schema"] == gate.RECEIPT_SCHEMA
    assert value["status"] == "passed" and set(value["replays"]) == gate.ROLES


@pytest.mark.parametrize(
    "change",
    [
        "unknown",
        "missing",
        "wrong_attempt",
        "wrong_commit",
        "bad_hash",
        "failed",
        "forged_wrapper",
        "wrong_artifact",
    ],
)
def test_reject_unbound_or_incomplete_receipts(tmp_path, replays, change):
    _, replay_dir = replays
    args = {}
    if change == "unknown":
        (replay_dir / "extra").mkdir()
    elif change == "missing":
        (replay_dir / "epw").rename(tmp_path / "removed")
    elif change == "wrong_attempt":
        args["run_attempt"] = "2"
    elif change == "wrong_commit":
        args["candidate"] = P
    else:
        receipt = replay_dir / "spectroscopy/role-receipt.json"
        data = json.loads(receipt.read_text())
        if change == "bad_hash":
            data["wrapper_receipt_sha256"] = "0" * 64
        elif change == "failed":
            data["command_exit"] = 1
        elif change == "wrong_artifact":
            data["artifacts"][0]["id"] += 1
        else:
            wrapper = replay_dir / "spectroscopy/wrapper-receipt.json"
            wrapper.write_text(
                json.dumps(
                    {
                        "status": "passed",
                        "schema": roles.SCHEMAS["spectroscopy"],
                        "source_revision": data["source_revision"],
                    }
                )
            )
            data["wrapper_receipt_sha256"] = gate.digest(wrapper)
        receipt.write_text(json.dumps(data))
    out = tmp_path / "accept.json"
    with pytest.raises(ValueError):
        aggregate(replays, out, **args)
    assert not out.exists()


def test_only_known_legacy_version_skips_replays(tmp_path, monkeypatch):
    env = dict(
        GITHUB_EVENT_NAME="push",
        GITHUB_REF="refs/tags/v1.1.0",
        GITHUB_SHA=K,
        GITHUB_REPOSITORY=gate.REPOSITORY,
    )
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(gate, "git", lambda *args: K if args[0] == "rev-parse" else "")
    assert gate.prepare(tmp_path / "missing", "1.1.0", K) == ("legacy", K)
    with pytest.raises(ValueError, match="unsupported publication"):
        gate.prepare(tmp_path / "missing", "1.3.0", K)


def test_release_jobs_depend_on_scientific_acceptance():
    # Workflow YAML uses GitHub expressions and does not require a runtime
    # YAML parser in the project's verification dependencies.
    root = Path(__file__).resolve().parents[1]
    workflow = (root / ".github/workflows/publish.yml").read_text()
    assert "uses: ./.github/workflows/release-science.yml" in workflow
    assert "  release-sources:\n    needs: [collect, scientific-acceptance]" in workflow
    assert "  publish-pypi:\n    needs: [collect, release-sources]" in workflow
