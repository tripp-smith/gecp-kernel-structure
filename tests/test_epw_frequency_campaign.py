import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from experiments import epw_frequency as exp

NATIVE = Path("/tmp/gecp-eg-native-13915c8")


def test_partitions_exclude_all_reserved():
    parts = exp.partition(109)
    assert set(parts["train"]).isdisjoint(parts["calibration"] + parts["final"])
    assert parts["calibration"] == [36, 40, 44] and parts["final"] == [45, 49, 53]
    assert sorted(parts["train"] + parts["calibration"] + parts["final"]) == list(
        range(109)
    )
    with pytest.raises(ValueError):
        exp.partition(53)


def test_budget_and_matrix_normalization():
    assert exp.passed(
        {
            "green_scaled": 1e-7,
            "z": 1e-6,
            "gap_ev": 1e-7,
            "chi_ev": 1e-7,
            "imaginary_defect": 0,
        },
        0.1,
    )
    assert not exp.passed(
        {
            "green_scaled": 1.01e-7,
            "z": 0,
            "gap_ev": 0,
            "chi_ev": 0,
            "imaginary_defect": 0,
        },
        0.1,
    )


@pytest.mark.skipif(not NATIVE.exists(), reason="actual native source not available")
def test_actual_native_context_and_explicit_offset():
    data, context, path = exp.native_context(NATIVE, 0.0)
    assert (
        data.z.shape == (109, 140)
        and context["reference_fermi_energy_ev"] == 7.6644748416
    )
    assert context["source_provenance"]["corpus_sha256"][path.name] == exp.sha(path)
    with pytest.raises(ValueError, match="offset"):
        exp.native_context(NATIVE, 0.01)


def test_calibration_never_constructs_reserved_rows(tmp_path, monkeypatch):
    env = {"source_hashes": {"x": "y"}, "tracked_clean": True}
    monkeypatch.setattr(exp, "source_environment", lambda **kw: env)
    monkeypatch.setattr(
        exp,
        "native_context",
        lambda root, offset: (
            SimpleNamespace(indices=np.arange(109)),
            {},
            tmp_path / "native",
        ),
    )
    (tmp_path / "native").write_text("source")
    seen = []

    def fit(data, rows, rank):
        assert set(rows).isdisjoint(exp.CAL + exp.FINAL)

        class FakeModel:
            def save(self, path):
                path.write_text("fit manifest")
                path.with_suffix(".npz").write_bytes(b"fit archive")

            def evaluate_matsubara(self, indices):
                return np.asarray(indices)

        model = FakeModel()
        model.diagnostics = {
            "rank": rank,
            "timing_seconds": {
                k: 0.0
                for k in (
                    "basis_construction",
                    "nambu_observation_construction",
                    "least_squares_preparation",
                    "coefficient_fit",
                    "training_replay",
                )
            },
        }
        model.rank = rank
        return model

    monkeypatch.setattr(exp.EPWFrequencyCrossFit, "load", lambda path: fit(None, [], 0))

    def measure(data, model, rows):
        seen.extend(rows)
        assert set(rows).isdisjoint(exp.FINAL)
        metrics = {
            "green_scaled": 0,
            "z": 0,
            "gap_ev": 0,
            "chi_ev": 0,
            "imaginary_defect": 0,
        }
        if model.rank == 48:
            metrics["z"] = 1
        return metrics, {}

    monkeypatch.setattr(exp, "fit", fit)
    monkeypatch.setattr(exp, "measure", measure)
    report = exp.calibrate(tmp_path, tmp_path / "out", tmp_path / "freeze", 0)
    assert report["selected_rank"] == 60 and [
        a["passed"] for a in report["attempts"]
    ] == [False, True]
    for attempt in report["attempts"]:
        rank = attempt["rank"]
        assert attempt["fit_json_sha256"] == exp.sha(
            tmp_path / "out" / f"rank-{rank}-fit.json"
        )
        assert attempt["fit_npz_sha256"] == exp.sha(
            tmp_path / "out" / f"rank-{rank}-fit.npz"
        )
    assert (
        set(seen).isdisjoint(exp.FINAL)
        and exp.load(tmp_path / "freeze")["selected_rank"] == 60
    )


def test_no_pass_no_freeze(tmp_path, monkeypatch):
    monkeypatch.setattr(exp, "source_environment", lambda **kw: {"tracked_clean": True})
    monkeypatch.setattr(
        exp,
        "native_context",
        lambda *a: (SimpleNamespace(indices=np.arange(109)), {}, tmp_path / "native"),
    )
    monkeypatch.setattr(
        exp,
        "fit",
        lambda *a: (_ for _ in ()).throw(exp.NumericalBreakdownError("conditioning")),
    )
    with pytest.raises(ValueError, match="no rank"):
        exp.calibrate(tmp_path, tmp_path / "out", tmp_path / "freeze", 0)
    assert (
        not (tmp_path / "freeze").exists()
        and len(exp.load(tmp_path / "out/calibration.json")["attempts"]) == 3
    )


def test_bad_freeze_rejected(monkeypatch, tmp_path):
    f = tmp_path / "freeze"
    f.write_text(json.dumps({"schema": "bad", "status": "passed"}))
    monkeypatch.setattr(exp, "source_environment", lambda: {})
    with pytest.raises(ValueError, match="freeze"):
        exp.checked_freeze(f, tmp_path, 0)


def test_campaign_rejects_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        exp.campaign(tmp_path, tmp_path / "freeze", 0, tmp_path)


def test_campaign_child_failure_records_receipt(tmp_path, monkeypatch):
    freeze = tmp_path / "freeze"
    freeze.write_text("{}")
    monkeypatch.setattr(
        exp, "checked_freeze", lambda *a: ({}, {"revision": "a"}, None, None)
    )
    monkeypatch.setattr(
        exp.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=3)
    )
    result = exp.campaign(tmp_path, freeze, 0, tmp_path / "out")
    assert result["status"] == "failed" and result["children"][0]["exit_code"] == 3
    assert exp.load(tmp_path / "out/campaign.json")["status"] == "failed"


def test_measure_uses_per_matrix_scale_and_independent_recovery(monkeypatch):
    truth = np.zeros((2, 1, 2, 2), complex)
    truth[:, 0, 0, 0] = -1j
    truth[:, 0, 1, 1] = -1j
    predicted = truth.copy()
    predicted[0] *= 1.00000001
    predicted[1] *= 1.00000002
    data = SimpleNamespace(
        indices=np.arange(2),
        frequencies_ev=np.ones(2),
        state_energies_ev=np.zeros(1),
        context={"chemical_potential_offset_ev": 0},
        z=np.ones((2, 1)),
        delta_ev=np.zeros((2, 1)),
        chi_ev=np.zeros((2, 1)),
    )
    monkeypatch.setattr(exp, "nambu_green", lambda data, rows: truth)
    model = SimpleNamespace(evaluate_matsubara=lambda rows: predicted)
    metrics, arrays = exp.measure(data, model, [0, 1])
    assert metrics["green_scaled"] == pytest.approx(2e-8, rel=1e-7)
    assert metrics["z"] == pytest.approx(2e-8, rel=1e-6)
    assert np.array_equal(arrays["truth_green"], truth)


def test_repeat_array_disagreement_fails_campaign(tmp_path, monkeypatch):
    freeze = tmp_path / "freeze"
    freeze.write_text("{}")
    env = {"revision": "a"}
    monkeypatch.setattr(exp, "checked_freeze", lambda *a: ({}, env, None, None))

    def child(cmd, **kwargs):
        folder = Path(cmd[cmd.index("--output") + 1])
        folder.mkdir()
        repeat = int(cmd[cmd.index("--repeat") + 1])
        array = folder / "scientific-arrays.npz"
        np.savez(array, x=[repeat])
        (folder / "case.json").write_text(
            json.dumps(
                {
                    "status": "passed",
                    "freeze_sha256": exp.sha(freeze),
                    "repeat": repeat,
                    "provenance": env,
                    "files_sha256": {array.name: exp.sha(array)},
                }
            )
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(exp.subprocess, "run", child)
    result = exp.campaign(tmp_path, freeze, 0, tmp_path / "out")
    assert (
        result["status"] == "failed"
        and "arrays differ" in result["error"]
        and len(result["children"]) == 3
    )


def test_missing_candidate_archive_is_fatal(tmp_path, monkeypatch):
    monkeypatch.setattr(exp, "source_environment", lambda **kw: {"tracked_clean": True})
    monkeypatch.setattr(
        exp,
        "native_context",
        lambda *args: (
            SimpleNamespace(indices=np.arange(109)),
            {},
            tmp_path / "native",
        ),
    )
    timings = {
        k: 0.0
        for k in (
            "basis_construction",
            "nambu_observation_construction",
            "least_squares_preparation",
            "coefficient_fit",
            "training_replay",
        )
    }

    class Model:
        diagnostics = {"timing_seconds": timings}

        def save(self, path):
            path.write_text("missing archive")

        def evaluate_matsubara(self, indices):
            return np.asarray(indices)

    monkeypatch.setattr(exp, "fit", lambda *args: Model())
    monkeypatch.setattr(exp.EPWFrequencyCrossFit, "load", lambda path: Model())
    monkeypatch.setattr(
        exp,
        "measure",
        lambda *args: (
            {
                "green_scaled": 0,
                "z": 0,
                "gap_ev": 0,
                "chi_ev": 0,
                "imaginary_defect": 0,
            },
            {},
        ),
    )
    with pytest.raises(FileNotFoundError):
        exp.calibrate(tmp_path, tmp_path / "out", tmp_path / "freeze", 0)
    report = exp.load(tmp_path / "out" / "calibration.json")
    assert report["status"] == "failed" and len(report["attempts"]) == 1
    assert (
        report["attempts"][0]["status"] == "failed"
        and not (tmp_path / "freeze").exists()
    )


def test_native_validation_failure_records_receipt(tmp_path, monkeypatch):
    monkeypatch.setattr(exp, "source_environment", lambda **kw: {"tracked_clean": True})
    monkeypatch.setattr(
        exp,
        "native_context",
        lambda *args: (_ for _ in ()).throw(ValueError("bad input hash")),
    )
    with pytest.raises(ValueError, match="bad input hash"):
        exp.calibrate(tmp_path, tmp_path / "out", tmp_path / "freeze", 0)
    report = exp.load(tmp_path / "out" / "calibration.json")
    assert report["status"] == "failed" and "bad input hash" in report["error"]
