"""Negative tests for the independent spectroscopy evidence auditor."""

import hashlib
import json

import numpy as np
import pytest

from experiments import audit_spectroscopy as audit


def _hash(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def _dataset(tmp_path):
    t = np.linspace(0, 10, 65)
    y = np.sin(t)
    c = np.eye(65) * 1e-8
    p = np.eye(65)
    p[-1, -1] = 0.0
    values = p @ y
    cov = p @ c @ p.T
    arrays = {
        "tau": t,
        "values": values,
        "covariance": cov,
        "projection": p,
        "projection_source_tau": t,
    }
    path = tmp_path / "sample-dataset.npz"
    np.savez_compressed(path, **arrays)
    config = {"default_model": "flat"}
    meta = {
        "schema": "kernelgecp.spectral-dataset.v1",
        "observation_domain": "imaginary_time",
        "beta": 10.0,
        "units": {"tau": "inverse_energy", "omega": "energy", "G_tau": "dimensionless"},
        "state_ordering": ["scalar"],
        "provider": {"name": "kernelgecp", "continuation": "TRIQS/MaxEnt4.0.1"},
        "calculation_settings": config,
        "file_sha256": audit.sha(path),
        "covariance_rtol": 1e-10,
        "covariance_atol": 1e-18,
        "effective_rank": 64,
        "fingerprint": audit.fingerprint(10, 1e-10, 1e-18, t, values, cov, p, t),
        "source_fingerprint": audit.fingerprint(10, 1e-10, 1e-18, t, y, c),
        "arrays": {
            k: {"shape": list(v.shape), "sha256": _hash(v)} for k, v in arrays.items()
        },
    }
    path.with_suffix(".json").write_text(json.dumps(meta))
    a = {
        "tau": t,
        "values": y,
        "covariance": c,
        "projected_values": values,
        "projected_covariance": cov,
        "projection": p,
    }
    result = {"method": "kernelgecp", "config": config, "effective_rank": 64}
    return path, meta, a, result


def test_covariance_cutoff_and_corruption():
    e, v = audit.modes(np.diag([1.0, 1e-12, 0.0]))
    assert e.size == 1
    assert audit.rms(np.array([2.0, 99.0, 8.0]), e, v) == pytest.approx(2.0)
    with pytest.raises(ValueError, match="indefinite"):
        audit.modes(np.diag([1.0, -0.1]))
    with pytest.raises(ValueError, match="asymmetric"):
        audit.modes(np.array([[1.0, 0.1], [0.0, 1.0]]))
    with pytest.raises(ValueError, match="no retained"):
        audit.modes(np.zeros((2, 2)))


def test_dataset_checks_full_hashes_and_source_fingerprint(tmp_path):
    path, meta, a, result = _dataset(tmp_path)
    assert len(audit.dataset_check(path, a, result)[0]) == 64
    meta["source_fingerprint"] = "0" * 64
    path.with_suffix(".json").write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="source fingerprint"):
        audit.dataset_check(path, a, result)


def test_dataset_rejects_rehashed_changed_covariance(tmp_path):
    path, meta, a, result = _dataset(tmp_path)
    a["projected_covariance"] = a["projected_covariance"].copy()
    a["projected_covariance"][0, 0] *= 2
    with pytest.raises(ValueError, match="dataset covariance"):
        audit.dataset_check(path, a, result)


def test_dataset_rejects_tampered_array_even_if_file_hash_replaced(tmp_path):
    path, meta, a, result = _dataset(tmp_path)
    with np.load(path, allow_pickle=False) as z:
        arr = {k: z[k] for k in z.files}
    arr["values"][0] += 1
    np.savez_compressed(path, **arr)
    meta["file_sha256"] = audit.sha(path)
    path.with_suffix(".json").write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="array hash"):
        audit.dataset_check(path, a, result)


def test_inventory_and_failure_output(tmp_path):
    assert len(dict(audit.jobs(True))) == 24 and len(dict(audit.jobs(False))) == 72
    assert all(
        j["bootstrap"] == (8 if j["repeat"] == 0 else 0) for _, j in audit.jobs(False)
    )
    manifest = tmp_path / "campaign.json"
    freeze = tmp_path / "freeze.json"
    out = tmp_path / "audit.json"
    manifest.write_text("{}")
    freeze.write_text("{}")
    import subprocess
    import sys

    run = subprocess.run(
        [
            sys.executable,
            "-m",
            "experiments.audit_spectroscopy",
            "--manifest",
            str(manifest),
            "--freeze",
            str(freeze),
            "--output",
            str(out),
        ],
        check=False,
    )
    assert run.returncode != 0 and json.loads(out.read_text())["status"] == "failed"


def test_analytic_truth_and_cdf():
    t = np.array([0.0, 5.0, 10.0])
    actual = audit.truth("offgrid_atoms_independent", t)
    assert np.all(actual < 0) and actual.shape == (3,)
    assert audit.true_cdf(
        "offgrid_atoms_correlated", np.array([-4.0, -1.0, 1.0, 4.0])
    ).tolist() == pytest.approx([0, 0.2, 0.75, 1])
    assert audit.true_cdf(
        "smooth_correlated", np.array([-4.0, 4.0])
    ).tolist() == pytest.approx([0, 1])


def _campaign(tmp_path):
    control = {
        "cases": list(audit.CASES),
        "methods": list(audit.METHODS),
        "defaults": list(audit.DEFAULTS),
        "eps": 1e-8,
        "n_omega": 401,
        "n_alpha": 33,
        "minimizer_mu0": 1e-8,
        "calibration_seeds": [101, 102, 103],
        "calibration_query_order": 24,
        "final_seeds": [3001, 3002, 3003],
        "final_query_order": 64,
        "bootstrap_first_repeat": 8,
        "repeats": 3,
    }
    p = {
        "source_revision": "fixture",
        "source_hashes": {},
        "tracked_clean": True,
        "python_executable": "fixture",
    }
    freeze = tmp_path / "freeze.json"
    freeze.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.spectroscopy-freeze.v1",
                "status": "passed",
                "controls": control,
                "provenance": p,
            }
        )
    )
    m = {
        "schema": "kernelgecp.spectroscopy-campaign.v1",
        "kind": "campaign",
        "status": "passed",
        "controls": control,
        "provenance": p,
        "freeze_sha256": audit.sha(freeze),
        "children": {},
    }
    for name, job in audit.jobs(False):
        fit = None
        if job["method"] == "kernelgecp":
            fit = {"path": name + "-gecp-fit.npz"}
            (tmp_path / fit["path"]).write_bytes(b"fit")
        meta = {
            "provenance": {k: v for k, v in p.items() if k != "python_executable"},
            "coefficient_fit_archive": fit,
        }
        path = tmp_path / (name + ".json")
        path.write_text(json.dumps(meta))
        for suffix in (".npz", "-dataset.npz", "-dataset.json", ".log"):
            (tmp_path / (name + suffix)).write_bytes(b"fixture")
        hashes = {
            "json_sha256": audit.sha(path),
            "npz_sha256": audit.sha(path.with_suffix(".npz")),
            "dataset_sha256": audit.sha(tmp_path / (name + "-dataset.npz")),
            "dataset_json_sha256": audit.sha(tmp_path / (name + "-dataset.json")),
            "fit_sha256": audit.sha(tmp_path / fit["path"]) if fit else None,
        }
        m["children"][name] = {
            "job": job,
            "status": "passed",
            "exit_code": 0,
            "hashes": hashes,
            "log_sha256": audit.sha(tmp_path / (name + ".log")),
        }
    path = tmp_path / "campaign.json"
    path.write_text(json.dumps(m))
    return path, freeze, m


def test_rejects_incomplete_campaign_and_nonzero_child(tmp_path, monkeypatch):
    path, freeze, m = _campaign(tmp_path)
    monkeypatch.setattr(audit, "check_provenance", lambda p: None)
    removed = next(iter(m["children"]))
    record = m["children"].pop(removed)
    path.write_text(json.dumps(m))
    with pytest.raises(ValueError, match="incomplete"):
        audit.audit_campaign(path, freeze)
    record["exit_code"] = 1
    m["children"][removed] = record
    path.write_text(json.dumps(m))
    with pytest.raises(ValueError, match="job/exit"):
        audit.audit_campaign(path, freeze)


def test_rejects_changed_repeat_even_with_valid_hashes(tmp_path, monkeypatch):
    path, freeze, _ = _campaign(tmp_path)
    monkeypatch.setattr(audit, "check_provenance", lambda p: None)

    def fake(path, job, full):
        arrays = {k: np.array([1.0]) for k in audit.ARRAYS}
        if job["repeat"] == 1:
            arrays["density"] = np.array([2.0])
        return arrays

    monkeypatch.setattr(audit, "audit_case", fake)
    with pytest.raises(ValueError, match="repeat mismatch"):
        audit.audit_campaign(path, freeze)
