import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from experiments import spectroscopy_campaign as driver


@pytest.fixture
def setup(monkeypatch, tmp_path):
    monkeypatch.setattr(driver.science, "CASES", ("one",))
    monkeypatch.setattr(driver.science, "METHODS", ("raw",))
    monkeypatch.setattr(driver, "DEFAULTS", ("flat",))
    env = {
        "source_revision": "a" * 40,
        "source_hashes": {"a": "b"},
        "tracked_clean": True,
        "threads": {k: "1" for k in driver.science.THREADS},
        "python": "test",
        "python_executable": "test",
    }
    monkeypatch.setattr(driver, "environment", lambda _: env.copy())
    state = {"calls": 0, "mode": "ok"}

    def child(cmd, stdout, stderr, timeout, cwd):
        state["calls"] += 1

        def option(flag):
            return cmd[cmd.index(flag) + 1]

        out = Path(option("--output"))
        npz = out.with_suffix(".npz")
        mode = state["mode"]
        if mode == "exit":
            return SimpleNamespace(returncode=3)
        if mode == "missing":
            return SimpleNamespace(returncode=0)
        x = 2 if mode == "different" and "-r1." in str(out) else 1
        np.savez(
            npz,
            density=np.array([x]),
            bootstrap_spectra=np.array([8] if option("--bootstrap") == "8" else []),
        )
        np.savez(out.with_name(out.stem + "-dataset.npz"), x=[1])
        if mode != "manifest":
            out.with_name(out.stem + "-dataset.json").write_text("{}\n")
        result = {
            "status": "passed",
            "gates": {"positive": True},
            "case": option("--case"),
            "method": option("--method"),
            "seed": int(option("--seed")),
            "eps": float(option("--eps")),
            "pilot": "--pilot" in cmd,
            "query_order": int(option("--query-order")),
            "bootstrap_count": int(option("--bootstrap")),
            "config": {
                "default_model": option("--default-model"),
                "n_omega": int(option("--n-omega")),
                "n_alpha": int(option("--n-alpha")),
            },
            "provenance": {k: v for k, v in env.items() if k != "python_executable"},
            "array_sha256": driver.sha(npz),
            "coefficient_fit_archive": None,
        }
        if mode == "gates":
            result["gates"] = {}
        out.write_text(json.dumps(result))
        if mode == "corrupt":
            npz.write_bytes(b"bad")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(driver.subprocess, "run", child)
    return tmp_path, env, state


def calibrate(ctx):
    tmp, _, _ = ctx
    return driver.execute("calibration", tmp / "cal", tmp / "freeze.json", "test")


def test_success_freeze_and_repeat_resume(setup):
    tmp, _, state = setup
    assert calibrate(setup)["status"] == "passed"
    freeze = driver.read(tmp / "freeze.json")
    assert freeze["controls"]["final_seeds"] == [2001, 2002, 2003]
    result = driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")
    assert result["status"] == "passed" and len(result["children"]) == 3
    assert len(result["repeat_comparison"]["skipped"]) == 2
    count = state["calls"]
    assert (
        driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")["status"]
        == "passed"
    )
    assert state["calls"] == count


@pytest.mark.parametrize("mode", ["exit", "missing", "gates", "corrupt", "manifest"])
def test_failed_calibration_never_freezes(setup, mode):
    tmp, _, state = setup
    state["mode"] = mode
    result = calibrate(setup)
    assert result["status"] == "failed" and not (tmp / "freeze.json").exists()
    assert (tmp / "cal" / "campaign.json").exists()


def test_orphan_fails_without_retry(setup):
    tmp, _, state = setup
    calibrate(setup)
    out = tmp / "final"
    out.mkdir()
    (out / (driver.stem(next(driver.jobs(False))) + ".npz")).write_bytes(b"orphan")
    before = state["calls"]
    result = driver.execute("campaign", out, tmp / "freeze.json", "test")
    assert (
        result["status"] == "failed"
        and "orphan" in result["error"]
        and state["calls"] == before
    )


def test_corrupt_resume_fails(setup):
    tmp, _, state = setup
    calibrate(setup)
    assert (
        driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")["status"]
        == "passed"
    )
    path = tmp / "final" / (driver.stem(next(driver.jobs(False))) + ".npz")
    path.write_bytes(b"bad")
    before = state["calls"]
    result = driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")
    assert result["status"] == "failed" and state["calls"] == before


def test_source_change_rejected(setup):
    tmp, env, _ = setup
    calibrate(setup)
    env["source_hashes"] = {"changed": "yes"}
    with pytest.raises(ValueError, match="mismatch"):
        driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")


def test_repeat_difference_is_failure(setup):
    tmp, _, state = setup
    calibrate(setup)
    state["mode"] = "different"
    result = driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")
    assert result["status"] == "failed" and "repeat array differs" in result["error"]


def test_deadline_is_inconclusive(setup, monkeypatch):
    tmp, _, _ = setup
    monkeypatch.setattr(driver, "LIMIT_SECONDS", 0)
    result = calibrate(setup)
    assert result["status"] == "inconclusive" and not (tmp / "freeze.json").exists()


def test_changed_frozen_config_rejected(setup):
    tmp, _, _ = setup
    calibrate(setup)
    path = tmp / "freeze.json"
    frozen = driver.read(path)
    frozen["controls"]["n_alpha"] = 99
    path.write_text(json.dumps(frozen))
    with pytest.raises(ValueError, match="mismatch"):
        driver.execute("campaign", tmp / "final", path, "test")


def test_failed_child_is_not_retried(setup):
    tmp, _, state = setup
    calibrate(setup)
    state["mode"] = "exit"
    first = driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")
    assert first["status"] == "failed"
    before = state["calls"]
    state["mode"] = "ok"
    second = driver.execute("campaign", tmp / "final", tmp / "freeze.json", "test")
    assert second["status"] == "failed" and state["calls"] == before


def test_failed_provider_retains_projected_input_and_failed_receipt(
    monkeypatch, tmp_path
):
    import sys

    from experiments import spectroscopy
    from kernelgecp.integrations.maxent import load_spectral_dataset
    from kernelgecp.types import NumericalBreakdownError

    class Basis:
        beta = 10

        def design_tau(self, tau):
            return np.polynomial.legendre.legvander(2 * np.asarray(tau) / 10 - 1, 3)

    monkeypatch.setattr(spectroscopy, "basis", lambda *args: Basis())
    monkeypatch.setattr(spectroscopy, "provenance", lambda: {})
    output = tmp_path / "nested" / "case.json"
    dataset_path = output.with_name("case-dataset.npz")

    def fail(fitted, config):
        restored = load_spectral_dataset(dataset_path)
        np.testing.assert_array_equal(restored.covariance, fitted.covariance)
        np.testing.assert_array_equal(restored.projection, fitted.projection)
        np.testing.assert_array_equal(
            restored.projection_source_tau, fitted.projection_source_tau
        )
        raise NumericalBreakdownError("stub provider failure")

    monkeypatch.setattr(spectroscopy, "continue_tau_maxent", fail)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "spectroscopy",
            "--case",
            "offgrid_atoms_correlated",
            "--method",
            "pydlr",
            "--seed",
            "102",
            "--pilot",
            "--output",
            str(output),
        ],
    )
    with pytest.raises(SystemExit) as caught:
        spectroscopy.main()
    assert caught.value.code == 1
    assert json.loads(output.read_text())["status"] == "failed"
    assert dataset_path.exists() and dataset_path.with_suffix(".json").exists()
    assert not output.with_suffix(".npz").exists()
