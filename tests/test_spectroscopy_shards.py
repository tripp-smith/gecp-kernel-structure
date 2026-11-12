import shutil

import pytest

from experiments import spectroscopy_campaign as driver


def test_partition_keeps_families_together():
    groups = [driver.shard_jobs(i) for i in range(6)]
    assert [len(g) for g in groups] == [12] * 6
    assert [j for group in groups for j in group] == list(driver.jobs(False, "v4"))
    for group in groups:
        assert all(
            [j["repeat"] for j in group[n : n + 3]] == [0, 1, 2]
            for n in range(0, 12, 3)
        )
        assert all(
            group[n]["bootstrap"] == 8
            and group[n + 1]["bootstrap"] == group[n + 2]["bootstrap"] == 0
            for n in range(0, 12, 3)
        )
    with pytest.raises(ValueError):
        driver.shard_jobs(6)


def test_assembly_preserves_and_rejects_tampering(tmp_path, monkeypatch):
    env = {"source_revision": "frozen", "python_executable": "python"}
    monkeypatch.setattr(driver, "environment", lambda *args: env)
    monkeypatch.setattr(
        driver,
        "compare_repeats",
        lambda d, m, *a: m.update(repeat_comparison={"status": "passed"}),
    )
    monkeypatch.setattr(
        driver,
        "validate_child",
        lambda j, r, e, d: {"json_sha256": driver.sha(d / (driver.stem(j) + ".json"))},
    )
    freeze = tmp_path / "freeze.json"
    driver.write(
        freeze,
        {
            "schema": "kernelgecp.spectroscopy-freeze.v4",
            "status": "passed",
            "controls": driver.controls("v4"),
            "provenance": env,
        },
    )
    root = tmp_path / "shards"
    root.mkdir()
    for i in range(6):
        folder = root / f"shard-{i}"
        folder.mkdir()
        selected = driver.shard_jobs(i)
        for name in driver.expected_files(selected):
            (folder / name).write_text('{"fake": true}')
        children = {}
        for job in selected:
            name = driver.stem(job)
            children[name] = {
                "job": job,
                "status": "passed",
                "exit_code": 0,
                "hashes": driver.validate_child(job, {}, env, folder),
                "log_sha256": driver.sha(folder / (name + ".log")),
            }
        driver.write(
            folder / "campaign.json",
            {
                "schema": "kernelgecp.spectroscopy-campaign-shard.v4",
                "kind": "campaign",
                "ordinal": i,
                "status": "passed",
                "controls": driver.controls("v4"),
                "provenance": env,
                "freeze_sha256": driver.sha(freeze),
                "children": children,
                "repeat_comparison": {"status": "passed"},
                "inventory": driver.inventory(folder, selected),
            },
        )
    result = driver.assemble_shards(tmp_path / "assembled", freeze, root)
    assert result["status"] == "passed" and len(result["children"]) == 72
    for i in range(6):
        assert (
            tmp_path / "assembled" / "shard-receipts" / f"shard-{i}.json"
        ).read_bytes() == (root / f"shard-{i}" / "campaign.json").read_bytes()
        for name in driver.expected_files(driver.shard_jobs(i)):
            assert (tmp_path / "assembled" / name).read_bytes() == (
                root / f"shard-{i}" / name
            ).read_bytes()
    victim = root / "shard-0" / next(iter(driver.expected_files(driver.shard_jobs(0))))
    victim.write_text("tampered")
    with pytest.raises(ValueError):
        driver.assemble_shards(tmp_path / "bad", freeze, root)
    victim.unlink()
    with pytest.raises(ValueError):
        driver.assemble_shards(tmp_path / "bad", freeze, root)
    shutil.rmtree(root / "shard-5")
    with pytest.raises(ValueError):
        driver.assemble_shards(tmp_path / "bad", freeze, root)
