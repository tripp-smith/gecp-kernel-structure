from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from experiments.cthyb_dmft import _loop
from kernelgecp.integrations.cthyb_dmft import CTHYBConfig, CTHYBMeasurement

IDENTITY = {"config": {"memory_limit_bytes": 10**12}}
DATA = {
    "max_iterations": 3,
    "calibration_loop_iterations": 3,
    "min_iterations": 3,
    "consecutive_passes": 1,
    "loop_tolerance": 0.5,
    "final_bath_average": 2,
}


class FakePlan:
    config = CTHYBConfig(n_iw=8, n_tau=97)

    def initial_bath(self):
        if os.environ.get("FAKE_MODE") == "verify":
            raise AssertionError("initial bath recomputed")
        return np.zeros(16, dtype=complex)

    def sample(self, bath, *, cycles, seed):
        if os.environ.get("FAKE_MODE") == "stop" and seed == 103:
            raise RuntimeError("simulated stop")
        return CTHYBMeasurement(
            self.config,
            seed,
            cycles,
            bath,
            np.zeros((97, 2), complex),
            np.zeros((16, 2), complex),
            np.zeros((16, 2), complex),
            np.zeros((16, 2), complex),
            np.array([[0, 0], [1, 1], [0, 0]], complex),
            np.ones((2, 2), complex),
            np.array([0.5, 0.5]),
            0.2,
            1.0,
            0,
        )

    def represent(self, measurement):
        if os.environ.get("FAKE_MODE") == "verify" or (
            os.environ.get("FAKE_MODE") == "resume" and measurement.seed == 101
        ):
            raise AssertionError("historical representation replayed")

        class Representation:
            native_iw_discrepancy = 0.0

        return Representation()

    def bath_update(self, measurement, representation):
        if os.environ.get("FAKE_MODE") == "verify":
            raise AssertionError("historical update replayed")
        delta = (
            0.125
            if os.environ.get("FAKE_MODE") == "stop"
            else np.nextafter(0.125, np.inf)
        )
        return measurement.bath_iw + delta, 0.1 if measurement.seed == 105 else 1.0


DATA.update(dataclasses.asdict(FakePlan.config))
IDENTITY["config"].update(dataclasses.asdict(FakePlan.config))


def run(folder):
    return _loop(FakePlan(), DATA, 100, 100, folder, IDENTITY, compressed=True)


def test_fresh_process_resume_and_completed_state(tmp_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
    env["FAKE_MODE"] = "stop"
    first = subprocess.run(
        [sys.executable, __file__, str(tmp_path)],
        env=env,
        capture_output=True,
        text=True,
    )
    assert first.returncode != 0 and "simulated stop" in first.stderr
    env["FAKE_MODE"] = "resume"
    subprocess.run([sys.executable, __file__, str(tmp_path)], env=env, check=True)
    env["FAKE_MODE"] = "verify"
    subprocess.run([sys.executable, __file__, str(tmp_path)], env=env, check=True)
    with np.load(tmp_path / "iteration-01.npz") as z:
        assert np.array_equal(z["bath_iw"], np.full(16, 0.125, complex))
    with (
        np.load(tmp_path / "iteration-01.npz") as a,
        np.load(tmp_path / "iteration-02.npz") as b,
    ):
        expected = np.mean([a["bath_iw"], b["bath_iw"]], axis=0)
    # The exact completed bath is read from state even if mean cannot be called.
    os.environ["FAKE_MODE"] = "verify"
    original = np.mean
    try:
        np.mean = lambda *a, **k: (_ for _ in ()).throw(AssertionError("mean replayed"))
        bath, report = run(tmp_path)
    finally:
        np.mean = original
        os.environ.pop("FAKE_MODE")
    assert report["converged"] and np.array_equal(bath, expected)


def test_orphan_legacy_identity_and_hash_fail(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_MODE", "stop")
    with pytest.raises(RuntimeError, match="simulated stop"):
        run(tmp_path)
    monkeypatch.setenv("FAKE_MODE", "resume")
    with pytest.raises(ValueError, match="model config|head"):
        _loop(
            FakePlan(),
            DATA,
            100,
            100,
            tmp_path,
            {"config": {"memory_limit_bytes": 1}},
            compressed=True,
        )
    with pytest.raises(ValueError, match="seed or cycles"):
        _loop(FakePlan(), DATA, 101, 100, tmp_path, IDENTITY, compressed=True)
    with pytest.raises(ValueError, match="seed or cycles"):
        _loop(FakePlan(), DATA, 100, 102, tmp_path, IDENTITY, compressed=True)
    wrong = dataclasses.replace(FakePlan.config, n_tau=99)
    plan = FakePlan()
    plan.config = wrong
    with pytest.raises(ValueError, match="model config|config, seed or cycles"):
        _loop(plan, DATA, 100, 100, tmp_path, IDENTITY, compressed=True)
    npz = tmp_path / "iteration-00.npz"
    npz.write_bytes(npz.read_bytes() + b"x")
    with pytest.raises(ValueError, match="measurement hash"):
        run(tmp_path)
    (tmp_path / "state-head.json").unlink()
    with pytest.raises(ValueError, match="lacks state head"):
        run(tmp_path)


if __name__ == "__main__":
    run(Path(sys.argv[1]))


def test_wrong_bath_rejected_even_with_rebound_file_hashes(tmp_path, monkeypatch):
    import hashlib
    import json

    monkeypatch.setenv("FAKE_MODE", "normal")
    run(tmp_path)
    npz = tmp_path / "iteration-01.npz"
    with np.load(npz) as archive:
        arrays = {name: archive[name] for name in archive.files}
    arrays["bath_iw"][0] = np.nextafter(arrays["bath_iw"][0].real, np.inf)
    np.savez_compressed(npz, **arrays)
    receipt = tmp_path / "iteration-01.json"
    record = json.loads(receipt.read_text())

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    record["arrays_sha256"] = digest(npz)
    receipt.write_text(json.dumps(record))
    state1 = tmp_path / "state-01.json"
    record = json.loads(state1.read_text())
    record["measurement_npz_sha256"] = digest(npz)
    record["measurement_json_sha256"] = digest(receipt)
    state1.write_text(json.dumps(record))
    state2 = tmp_path / "state-02.json"
    record = json.loads(state2.read_text())
    record["previous_state_sha256"] = digest(state1)
    state2.write_text(json.dumps(record))
    head = tmp_path / "state-head.json"
    record = json.loads(head.read_text())
    record["state_sha256"] = digest(state2)
    head.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="measurement bath mismatch"):
        run(tmp_path)


def test_changed_source_and_legacy_report_rejected(tmp_path, monkeypatch):
    import hashlib

    import experiments.cthyb_dmft as campaign

    source = tmp_path / "frozen.py"
    source.write_text("original")
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    identity = {
        "config": {"memory_limit_bytes": 10**12},
        "sources": {"frozen.py": hashlib.sha256(source.read_bytes()).hexdigest()},
    }
    source.write_text("changed")
    with pytest.raises(ValueError, match="source changed"):
        _loop(FakePlan(), DATA, 100, 100, tmp_path / "loop", identity, compressed=True)
    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "loop.json").write_text("{}")
    with pytest.raises(ValueError, match="lacks state head"):
        run(legacy)


def test_checkpoint_interruption_boundaries_and_notification(
    tmp_path, monkeypatch, capsys
):
    import hashlib
    import json

    import experiments.cthyb_dmft as campaign

    original = campaign._write
    monkeypatch.setenv("FAKE_MODE", "normal")
    for stage in ("before_state", "before_head"):
        folder = tmp_path / stage

        def injected(path, value, stage=stage):
            if stage == "before_state" and path.name == "state-00.json":
                raise RuntimeError(stage)
            if stage == "before_head" and path.name == "state-head.json":
                raise RuntimeError(stage)
            original(path, value)

        with monkeypatch.context() as m:
            m.setattr(campaign, "_write", injected)
            with pytest.raises(RuntimeError, match=stage):
                run(folder)
        assert "measurement" not in capsys.readouterr().out
        with pytest.raises(ValueError, match="lacks state head"):
            run(folder)
    normal = tmp_path / "normal"
    original_print = __import__("builtins").print

    def checked_print(*args, **kwargs):
        event = json.loads(args[0])
        assert "measurement" in event
        head = json.loads((normal / "state-head.json").read_text())
        state = normal / f"state-{head['next_index'] - 1:02d}.json"
        assert hashlib.sha256(state.read_bytes()).hexdigest() == head["state_sha256"]
        original_print(*args, **kwargs)

    monkeypatch.setattr("builtins.print", checked_print)
    run(normal)


@pytest.mark.parametrize(
    "mode,defect,minimum,expected",
    [
        ("pilot", 0.1, 1, "calibration_pilot_budget"),
        ("final", 0.1, 1, "stochastic_controller"),
        ("final", 0.5, 1, "max_iterations_inconclusive"),
    ],
)
def test_controller_boundaries(tmp_path, monkeypatch, mode, defect, minimum, expected):
    monkeypatch.setenv("FAKE_MODE", "normal")

    class Plan(FakePlan):
        def bath_update(self, measurement, representation):
            return measurement.bath_iw + 0.125, defect

    data = dict(DATA, min_iterations=minimum)
    bath, report = _loop(
        Plan(),
        data,
        100,
        100,
        tmp_path,
        IDENTITY,
        compressed=True,
        pilot=mode == "pilot",
    )
    assert report["stop_reason"] == expected
    assert len(report["history"]) == (1 if expected == "stochastic_controller" else 3)
    monkeypatch.setenv("FAKE_MODE", "verify")
    again, restored = _loop(
        Plan(),
        data,
        100,
        100,
        tmp_path,
        IDENTITY,
        compressed=True,
        pilot=mode == "pilot",
    )
    assert np.array_equal(bath, again) and report == restored


def test_wrong_first_bath_and_final_mean_rejected_after_hash_rebinding(
    tmp_path, monkeypatch
):
    import hashlib
    import json

    from experiments.cthyb_loop_state import _decode, _encode

    monkeypatch.setenv("FAKE_MODE", "normal")
    for case in ("first", "final"):
        folder = tmp_path / case
        run(folder)

        def digest(path):
            return hashlib.sha256(path.read_bytes()).hexdigest()

        if case == "first":
            npz = folder / "iteration-00.npz"
            with np.load(npz) as archive:
                arrays = {name: archive[name] for name in archive.files}
            arrays["bath_iw"][0] = np.nextafter(0.0, 1.0)
            np.savez_compressed(npz, **arrays)
            receipt = folder / "iteration-00.json"
            r = json.loads(receipt.read_text())
            r["arrays_sha256"] = digest(npz)
            receipt.write_text(json.dumps(r))
            state = folder / "state-00.json"
            r = json.loads(state.read_text())
            r["measurement_npz_sha256"] = digest(npz)
            r["measurement_json_sha256"] = digest(receipt)
            state.write_text(json.dumps(r))
            for i in (1, 2):
                previous = folder / f"state-{i - 1:02d}.json"
                path = folder / f"state-{i:02d}.json"
                r = json.loads(path.read_text())
                r["previous_state_sha256"] = digest(previous)
                path.write_text(json.dumps(r))
        else:
            path = folder / "state-02.json"
            r = json.loads(path.read_text())
            bath = _decode(r["final_bath"], 16)
            bath[0] += 0.01
            r["final_bath"] = _encode(bath)
            path.write_text(json.dumps(r))
        head = folder / "state-head.json"
        r = json.loads(head.read_text())
        r["state_sha256"] = digest(folder / "state-02.json")
        head.write_text(json.dumps(r))
        with pytest.raises(
            ValueError,
            match="measurement bath mismatch"
            if case == "first"
            else "mean rounding bound",
        ):
            run(folder)


def test_final_mean_rounding_bound_extremes():
    from experiments.cthyb_loop_state import _validate_final_mean

    samples = [
        np.array([complex(1e200, 1e-310), complex(-1e200, -1e-310)], complex),
        np.array([complex(-1e200, -1e-310), complex(1e200, 1e-310)], complex),
        np.array([complex(3.0, 5e-324), complex(5.0, -5e-324)], complex),
    ]
    mean = np.mean(samples, axis=0)
    _validate_final_mean(mean, samples)
    small = [np.array([1.0 + 0j]), np.array([2.0 + 0j])]
    bad = np.mean(small, axis=0) + 1e-6
    with pytest.raises(ValueError, match="mean rounding bound"):
        _validate_final_mean(bad, small)


def test_missing_or_bad_predecessor_rejected(tmp_path, monkeypatch):
    import hashlib
    import json

    monkeypatch.setenv("FAKE_MODE", "normal")
    for mode in ("missing", "bad_link"):
        folder = tmp_path / mode
        run(folder)
        if mode == "missing":
            (folder / "state-00.json").unlink()
        else:
            state = folder / "state-02.json"
            value = json.loads(state.read_text())
            value["previous_state_sha256"] = "0" * 64
            state.write_text(json.dumps(value))
            head = folder / "state-head.json"
            value = json.loads(head.read_text())
            value["state_sha256"] = hashlib.sha256(state.read_bytes()).hexdigest()
            head.write_text(json.dumps(value))
        with pytest.raises(ValueError, match="missing loop state|chain or measurement"):
            run(folder)
