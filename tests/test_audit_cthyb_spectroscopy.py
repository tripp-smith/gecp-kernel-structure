import json

import numpy as np
import pytest

from experiments import audit_cthyb_spectroscopy as a


def test_covariance_two_pass():
    x = np.random.default_rng(4).normal(size=(32, 7))
    m, c = a.moments(x)
    np.testing.assert_allclose(m, x.mean(axis=0))
    np.testing.assert_allclose(c, np.cov(x, rowvar=False) / 32)


def test_covariance_modes_reject_indefinite():
    with pytest.raises(ValueError, match="negative covariance"):
        a.modes(np.diag([1.0, -0.1]))


def test_covariance_modes_reject_asymmetry():
    with pytest.raises(ValueError, match="asymmetry"):
        a.modes(np.array([[1.0, 0.1], [0.0, 1.0]]))


def test_adaptive_quadrature_against_symmetric_spectrum():
    w = np.linspace(-4, 4, 401)
    d = np.full(401, 1 / 8)
    pred = a.forward(np.array([0.0, 5.0, 10.0]), w, d)
    assert np.all(np.isfinite(pred))
    assert pred[0] == pytest.approx(-0.5, abs=1e-10)
    assert pred[2] == pytest.approx(-0.5, abs=1e-10)


def test_window_weights_constant_density():
    w = np.linspace(-4, 4, 401)
    np.testing.assert_allclose(a.weights(w, np.full(401, 1 / 8)), [0.375, 0.25, 0.375])


def test_bad_binding_fails_before_chain_access(tmp_path):
    p = tmp_path / "binding"
    p.write_text("{}")
    with pytest.raises(ValueError, match="binding hash"):
        a.raw_cohort(tmp_path, p, p, "calibration")


def test_observation_and_query_disjoint():
    assert (
        len(a.OBS) == 65
        and len(a.QUERY) == 16
        and not np.intersect1d(a.OBS, a.QUERY).size
    )


def test_close_rejects_nan():
    with pytest.raises(ValueError):
        a.close([np.nan], [np.nan], "bad")


def test_fake_provider_case_end_to_end(tmp_path, monkeypatch):
    """Producer writes a finite fake-provider case; the real auditor checks it."""
    import sys
    from types import ModuleType, SimpleNamespace

    from scipy.special import expit

    from experiments import cthyb_spectroscopy as producer

    full_tau = np.linspace(0, 10, 3073)
    omega = np.linspace(-4, 4, 401)
    density = np.full(401, 0.125)
    # This test-only spectrum defines the sample mean; no native cohort is used.
    reference = a.forward(full_tau, omega, density)
    samples = reference + np.random.default_rng(5).normal(0, 0.01, (64, 3073))
    monkeypatch.setattr(producer, "load_cohort", lambda *args: (samples, {}))

    def kernel(t):
        return (
            -np.exp(
                np.where(
                    omega[None, :] >= 0, -t[:, None] * omega, (10 - t[:, None]) * omega
                )
            )
            * expit(10 * np.abs(omega))[None, :]
        )

    widths = np.r_[0.01, np.full(399, 0.02), 0.01]

    class Spectrum:
        def __init__(self, dataset):
            self.omega = omega
            self.density = density
            self.alphas = np.geomspace(1e4, 1e-4, 33)
            self.trial_spectra = np.tile(density, (33, 1))
            residual = kernel(dataset.tau) @ (widths * density) - dataset.values
            transformed = dataset.modes.whiten(residual)
            self.chi2 = np.full(33, transformed @ transformed)
            self.selected_index = 0
            self.converged = (True,) * 33
            self.iterations = (1,) * 33
            self.effective_rank = dataset.modes.effective_rank
            self.raw_mass = self.mass = 1.0
            self.normalized_training_whitened_rms = 0.0
            self.diagnostics = {"provider": "fake test fixture"}

        def evaluate_tau(self, beta, tau):
            return a.forward(tau, omega, density)

    monkeypatch.setattr(
        producer, "continue_tau_maxent", lambda dataset, config: Spectrum(dataset)
    )
    provider = ModuleType("triqs_maxent")
    provider.LinearOmegaMesh = lambda *args: omega
    kernels = ModuleType("triqs_maxent.kernels")
    kernels.TauKernel = lambda tau, mesh, beta: SimpleNamespace(K=kernel(tau))
    monkeypatch.setitem(sys.modules, "triqs_maxent", provider)
    monkeypatch.setitem(sys.modules, "triqs_maxent.kernels", kernels)
    case = tmp_path / "raw-flat-r0.json"
    producer.run(
        tmp_path,
        tmp_path,
        tmp_path,
        "test",
        "calibration",
        "raw",
        "flat",
        case,
        bootstrap=8,
    )
    train, train_cov = a.moments(samples[:32, a.OBS])
    held, held_cov = a.moments(samples[32:, a.QUERY])
    stats = samples, train, train_cov, held, held_cov
    result, _ = a.audit_case(case, "raw", "flat", 0, stats)
    assert result["status"] == "passed"
    altered = {**result, "rcond": 1e-10}
    case.write_text(json.dumps(altered))
    with pytest.raises(ValueError, match="solver controls"):
        a.audit_case(case, "raw", "flat", 0, stats)
    case.write_text(json.dumps(result))
    with np.load(case.with_suffix(".npz"), allow_pickle=False) as archive:
        arrays = {k: archive[k] for k in archive.files}
    arrays["validation_covariance"] += train_cov[:16, :16]
    np.savez_compressed(case.with_suffix(".npz"), **arrays)

    result["array_sha256"] = a.sha(case.with_suffix(".npz"))
    case.write_text(json.dumps(result))
    with pytest.raises(ValueError, match="heldout covariance"):
        a.audit_case(case, "raw", "flat", 0, stats)


def test_24_test_only_workers_and_tampering(tmp_path):
    """Fake worker subprocesses exercise receipt and frozen-input verification."""
    import subprocess
    import sys

    root = tmp_path / "ee"
    root.mkdir()
    ee_audit = tmp_path / "ee-audit.json"
    ee_binding = tmp_path / "ee-binding.json"
    cal = tmp_path / "calibration.json"
    cal_audit = tmp_path / "calibration-audit.json"
    for path in (ee_audit, ee_binding, cal, cal_audit):
        path.write_text("{}")
    snapshot = {"executable": sys.executable}
    protocol = {"fixture": True}
    manifest = {
        "snapshot": snapshot,
        "protocol": protocol,
        "cases": {},
        "workers": {},
        "parent_pid": 9000000,
        "elapsed_seconds": 10.0,
        "global_timeout_seconds": 19800,
        "case_timeout_seconds": 3600,
    }
    for repeat in range(3):
        for method in a.METHODS:
            for model in a.MODELS:
                name = f"{method}-{model}-r{repeat}.json"
                result = tmp_path / name
                result.write_text("{}")
                result.with_suffix(".npz").write_bytes(b"test fixture")
                jobfile = tmp_path / (name + ".job")
                logfile = tmp_path / (name + ".log")
                workerfile = tmp_path / (name + ".worker")
                job = {
                    "schema": "kernelgecp.cthyb-spectroscopy-worker-job.v1",
                    "stage": "final",
                    "parent_pid": 9000000,
                    "method": method,
                    "model": model,
                    "repeat": repeat,
                    "protocol": protocol,
                    "snapshot": snapshot,
                    "folder": str(root),
                    "audit_receipt": str(ee_audit),
                    "audit_binding": str(ee_binding),
                    "binding_sha256": a.BINDING,
                    "output": str(result),
                    "calibration": str(cal),
                    "calibration_sha256": a.sha(cal),
                    "calibration_audit": str(cal_audit),
                    "calibration_audit_sha256": a.sha(cal_audit),
                }
                jobfile.write_text(json.dumps(job))
                hashes = {
                    "result_sha256": a.sha(result),
                    "array_sha256": a.sha(result.with_suffix(".npz")),
                }
                manifest["cases"][name] = hashes
                # Real fresh OS process; fixture receipt, no scientific provider.
                cmd = [sys.executable, "-c", "import os; print(os.getpid())"]
                child = subprocess.run(cmd, check=True, capture_output=True, text=True)
                pid = int(child.stdout.strip())
                logfile.write_text(child.stdout)
                workerfile.write_text(
                    json.dumps(
                        {
                            "schema": "kernelgecp.cthyb-spectroscopy-worker.v1",
                            "status": "passed",
                            "pid": pid,
                            "job_sha256": a.sha(jobfile),
                            **hashes,
                        }
                    )
                )
                manifest["workers"][name] = {
                    "status": "passed",
                    "exit_code": 0,
                    "pid": pid,
                    "elapsed_seconds": 0.1,
                    "job_file": jobfile.name,
                    "job_sha256": a.sha(jobfile),
                    "log_file": logfile.name,
                    "log_sha256": a.sha(logfile),
                    "worker_receipt_file": workerfile.name,
                    "worker_receipt_sha256": a.sha(workerfile),
                    "command": [
                        sys.executable,
                        "-m",
                        "experiments.cthyb_spectroscopy",
                        "--worker",
                        str(jobfile),
                        "--job-sha256",
                        a.sha(jobfile),
                    ],
                }
    campaign = tmp_path / "final.json"
    a.check_workers(
        manifest, campaign, root, ee_audit, ee_binding, "final", cal, cal_audit
    )
    assert len(manifest["workers"]) == 24
    records = list(manifest["workers"].values())
    original_pid = records[1]["pid"]
    records[1]["pid"] = records[0]["pid"]
    with pytest.raises(ValueError, match="PID"):
        a.check_workers(
            manifest, campaign, root, ee_audit, ee_binding, "final", cal, cal_audit
        )
    records[1]["pid"] = original_pid
    records[0]["exit_code"] = 1
    with pytest.raises(ValueError, match="exit"):
        a.check_workers(
            manifest, campaign, root, ee_audit, ee_binding, "final", cal, cal_audit
        )
    records[0]["exit_code"] = 0
    (tmp_path / records[0]["log_file"]).write_text("tampered")
    with pytest.raises(ValueError, match="artifact hash"):
        a.check_workers(
            manifest, campaign, root, ee_audit, ee_binding, "final", cal, cal_audit
        )
