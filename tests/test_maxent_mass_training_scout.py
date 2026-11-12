"""Training access and provenance checks; not native-provider acceptance."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).parents[1] / "scripts/providers/maxent_mass_training_scout.py"
spec = importlib.util.spec_from_file_location("training_scout", PATH)
scout = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scout)


def test_load_training_accesses_only_three_members():
    seen = []

    class Archive:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def __getitem__(self, name):
            seen.append(name)
            assert name in {"observation_tau", "observation_values", "covariance"}
            return np.array([1.0])

    class FakeNumpy:
        def load(self, *args, **kwargs):
            assert kwargs == {"allow_pickle": False}
            return Archive()

    scout.load_training(FakeNumpy(), "archive.npz")
    assert seen == ["observation_tau", "observation_values", "covariance"]


def test_wrong_source_fails_binding(monkeypatch, tmp_path):
    monkeypatch.setattr(scout, "TRACKED", ["source.py"])
    (tmp_path / "source.py").write_bytes(b"different")
    monkeypatch.setattr(scout.subprocess, "check_output", lambda cmd: b"git source")
    with pytest.raises(RuntimeError, match="working source differs"):
        scout.bind_source(tmp_path, "revision")


def test_wrong_fixture_fails_binding(monkeypatch, tmp_path):
    fixture = scout.FIXTURE + ".npz"
    monkeypatch.setattr(scout, "TRACKED", [fixture])
    path = tmp_path / fixture
    path.parent.mkdir(parents=True)
    path.write_bytes(b"wrong fixture")
    monkeypatch.setattr(scout.subprocess, "check_output", lambda cmd: b"wrong fixture")
    with pytest.raises(RuntimeError, match="fixture archive hash mismatch"):
        scout.bind_source(tmp_path, "revision")


def test_forward_quadrature_mass_endpoints():
    omega = np.linspace(-4, 4, 101)
    spectra = np.ones((1, len(omega))) / 8
    result = scout.forward12(np, 8.0, np.array([0.0, 8.0]), omega, spectra)
    assert np.allclose(result, [[-0.5, -0.5]], atol=1e-14)


def test_fixture_json_hash_is_bound(monkeypatch, tmp_path):
    import hashlib

    archive = scout.FIXTURE + ".npz"
    metadata = scout.FIXTURE + ".json"
    monkeypatch.setattr(scout, "TRACKED", [archive, metadata])
    payloads = {archive: b"archive", metadata: b"wrong metadata"}
    monkeypatch.setattr(scout, "ARRAY_SHA256", hashlib.sha256(b"archive").hexdigest())
    for name, payload in payloads.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    monkeypatch.setattr(
        scout.subprocess,
        "check_output",
        lambda cmd: payloads[cmd[-1].split(":", 1)[1]],
    )
    with pytest.raises(RuntimeError, match="fixture JSON hash mismatch"):
        scout.bind_source(tmp_path, "revision")


def test_capture_native_descriptors_and_bind_outputs(tmp_path):
    import hashlib
    import os

    log_path = tmp_path / "worker.log"
    with log_path.open("x") as log, scout.capture_fds(log):
        os.write(1, b"native stdout\n")
        os.write(2, b"native stderr\n")
        print("python stdout")
    data = log_path.read_bytes()
    assert b"native stdout" in data and b"native stderr" in data
    assert b"python stdout" in data
    (tmp_path / "identity-flat-trials.npz").write_bytes(b"trial")
    (tmp_path / "marquardt-flat-failed-solve.npz").write_bytes(b"solve")
    hashes = scout.output_hashes(tmp_path)
    assert hashes == {
        "worker.log": hashlib.sha256(data).hexdigest(),
        "identity-flat-trials.npz": hashlib.sha256(b"trial").hexdigest(),
        "marquardt-flat-failed-solve.npz": hashlib.sha256(b"solve").hexdigest(),
    }


def test_invalid_training_fails_before_continuation():
    tau = np.linspace(0, 10, 65)
    covariance = np.eye(65)
    covariance[0, 0] = np.nan
    with pytest.raises(RuntimeError, match="invalid finite"):
        scout.validate_training(np, tau, np.zeros(65), covariance)


def test_candidate_completion_never_hides_control_failure():
    complete = {
        f"{damping}/{model}": {"status": "completed"}
        for damping in ("identity", "marquardt")
        for model in ("flat", "broad_gaussian")
    }
    assert scout.summarize_cases(complete) == (True, "diagnostic_completed")
    complete["identity/flat"] = {"status": "failed"}
    assert scout.summarize_cases(complete) == (True, "failed")
    complete.pop("marquardt/flat")
    assert scout.summarize_cases(complete) == (False, "failed")
