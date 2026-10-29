import json
from types import SimpleNamespace

import numpy as np
import pytest

from experiments import cthyb_spectroscopy as m


def test_independent_mean_and_covariance():
    rng = np.random.default_rng(8)
    x = rng.normal(size=(64, 3073))
    data, query, held, cov = m.dataset_and_validation(x, {})
    assert np.array_equal(m.OBS, np.arange(0, 3073, 48))
    assert not set(m.OBS) & set(m.QUERY)
    assert len(query) == 16 and data.modes.relative_cutoff == 1e-4
    np.testing.assert_allclose(data.values, x[:32, m.OBS].mean(axis=0))
    np.testing.assert_allclose(
        data.covariance, np.cov(x[:32, m.OBS], rowvar=False) / 32, atol=1e-15
    )
    np.testing.assert_allclose(held, x[32:, m.QUERY].mean(axis=0))
    np.testing.assert_allclose(
        cov,
        np.cov(x[32:, m.QUERY], rowvar=False) / 32,
        atol=1e-15,
    )
    assert data.modes.effective_rank <= 31
    with pytest.raises(ValueError):
        m.dataset_and_validation(x[:63], {})


def fixture(tmp_path, monkeypatch):
    role = tmp_path / "dense-final"
    role.mkdir()
    identity = {
        "source_revision": m.REVISION,
        "config": {
            "seed_base": 29100000,
            "seed_mapping": "odd_stride_2_v1",
            "final_replicas": 64,
            "n_tau": 3073,
            "beta": 10,
            "interaction": 2,
        },
    }
    (tmp_path / "manifest.json").write_text(json.dumps(identity))
    (tmp_path / "freeze.json").write_text("{}")
    archives = {}
    for i in range(64):
        path = role / f"replica-{i:02d}.json"
        path.write_text(json.dumps({"index": i}))
        archives[str(path.relative_to(tmp_path))] = m.sha(path)
    result = {
        "identity": identity,
        "status": "passed",
        "freeze_sha256": m.sha(tmp_path / "freeze.json"),
        "measurement_archives": archives,
    }
    (tmp_path / "result.json").write_text(json.dumps(result))
    audit = {
        "schema": "kernelgecp.cthyb-audit.v1",
        "status": "passed",
        "generation_revision": m.REVISION,
        "gates": {"gate": True},
        "audit_revision": "revision",
        "audit_sha256": m.hashlib.sha256(b"audit").hexdigest(),
    }
    receipt = tmp_path / "audit.json"
    receipt.write_text(json.dumps(audit))
    monkeypatch.setattr(m.subprocess, "check_output", lambda *a, **k: b"audit")

    def load(stem, identity):
        i = int(stem.name.split("-")[-1])
        return SimpleNamespace(
            seed=29130001 + 2 * i,
            config=SimpleNamespace(beta=10, interaction=2, n_tau=3073),
            cycles=1600000,
            solve_status=0,
            average_sign=1.0,
            bath_iw=np.ones(4),
            g_tau=np.full((3073, 2), i + 0j),
        )

    monkeypatch.setattr(m, "load_measurement", load)
    binding = tmp_path / "binding.json"
    binding.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.cthyb-audit-binding.v1",
                "status": "passed",
                "generation_revision": m.REVISION,
                "audit_receipt_sha256": m.sha(receipt),
                "result_sha256": m.sha(tmp_path / "result.json"),
                "manifest_sha256": m.sha(tmp_path / "manifest.json"),
                "freeze_sha256": m.sha(tmp_path / "freeze.json"),
                "measurement_archives": archives,
            }
        )
    )
    return receipt, binding


def test_raw_cohort_and_failure_paths(tmp_path, monkeypatch):
    receipt, binding = fixture(tmp_path, monkeypatch)
    data, meta = m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")
    assert data.shape == (64, 3073) and meta["seeds"][0] == 29130001
    with pytest.raises(ValueError, match="digest"):
        m.load_cohort(tmp_path, receipt, binding, "0" * 64, "final")
    path = tmp_path / "dense-final" / "replica-01.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="receipt hash"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")


def test_duplicate_seed_and_imaginary_rejected(tmp_path, monkeypatch):
    receipt, binding = fixture(tmp_path, monkeypatch)

    def bad(stem, identity):
        return SimpleNamespace(
            seed=29130001,
            config=SimpleNamespace(beta=10, interaction=2, n_tau=3073),
            cycles=1600000,
            solve_status=0,
            average_sign=1.0,
            bath_iw=np.ones(4),
            g_tau=np.ones((3073, 2), complex),
        )

    monkeypatch.setattr(m, "load_measurement", bad)
    with pytest.raises(ValueError, match="seed"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")

    def imag(stem, identity):
        i = int(stem.name.split("-")[-1])
        return SimpleNamespace(
            seed=29130001 + 2 * i,
            config=SimpleNamespace(beta=10, interaction=2, n_tau=3073),
            cycles=1600000,
            solve_status=0,
            average_sign=1.0,
            bath_iw=np.ones(4),
            g_tau=np.full((3073, 2), 1 + 1e-4j),
        )

    monkeypatch.setattr(m, "load_measurement", imag)
    with pytest.raises(ValueError, match="imaginary"):
        m.load_cohort(tmp_path, receipt, binding, m.sha(binding), "final")


def test_heldout_noise_gate():
    cov = np.eye(16) * 0.01
    rms, modes, passed = m.validation_score(cov, np.full(16, 0.3))
    assert rms == pytest.approx(3) and modes.effective_rank == 16 and not passed
    with pytest.raises(ValueError):
        m.validation_score(np.zeros((16, 16)), np.ones(16))


def test_calibration_freeze_must_match_sources_and_cases(tmp_path):
    snapshot = {"sources": {"core": "original"}, "runtime": "pinned"}
    cases = {}
    for method in m.METHODS:
        for model in m.PROTOCOL["default_models"]:
            path = tmp_path / f"{method}-{model}-r0.json"
            archive = path.with_suffix(".npz")
            archive.write_bytes(b"finite fixture")
            path.write_text(
                json.dumps({"status": "passed", "array_sha256": m.sha(archive)})
            )
            cases[path.name] = {
                "result_sha256": m.sha(path),
                "array_sha256": m.sha(archive),
            }
    receipt = tmp_path / "calibration.json"
    receipt.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.cthyb-spectroscopy-calibration.v1",
                "status": "passed",
                "protocol": m.PROTOCOL,
                "snapshot": snapshot,
                "binding_sha256": "binding",
                "cases": cases,
            }
        )
    )
    m.verify_calibration(receipt, m.sha(receipt), snapshot, "binding")
    with pytest.raises(ValueError, match="requires"):
        m.verify_calibration(None, None, snapshot, "binding")
    with pytest.raises(ValueError, match="source"):
        m.verify_calibration(
            receipt, m.sha(receipt), {"sources": {"core": "changed"}}, "binding"
        )
    (tmp_path / "raw-flat-r0.npz").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="hash"):
        m.verify_calibration(receipt, m.sha(receipt), snapshot, "binding")


def test_real_measurement_archive_tampering(tmp_path):
    from experiments.cthyb_dmft import load_measurement, save_measurement
    from kernelgecp.integrations.cthyb_dmft import CTHYBConfig, CTHYBMeasurement

    c = CTHYBConfig(n_iw=8, n_tau=49)

    def z(shape):
        return np.zeros(shape, dtype=np.complex128)

    moments = z((3, 2))
    moments[1] = 1
    measurement = CTHYBMeasurement(
        config=c,
        seed=1,
        cycles=2,
        bath_iw=z((16,)),
        g_tau=z((49, 2)),
        g_iw=z((16, 2)),
        g_iw_raw=z((16, 2)),
        sigma_iw=z((16, 2)),
        g_moments=moments,
        sigma_moments=z((2, 2)),
        occupancies=np.array([0.5, 0.5]),
        double_occupancy=0.2,
        average_sign=1.0,
        solve_status=0,
    )
    stem = tmp_path / "measurement"
    save_measurement(stem, measurement, {"fixture": True})
    assert load_measurement(stem, {"fixture": True}).seed == 1
    with stem.with_suffix(".npz").open("ab") as output:
        output.write(b"tampered")
    with pytest.raises(ValueError, match="hash"):
        load_measurement(stem, {"fixture": True})


def test_independent_audit_is_required_and_bound(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="requires"):
        m.verify_calibration_audit(None, None, "cal", "binding")
    path = tmp_path / "independent.json"
    value = {
        "schema": "kernelgecp.cthyb-spectroscopy-independent-audit.v1",
        "status": "passed",
        "kind": "calibration",
        "calibration_sha256": "cal",
        "binding_sha256": "binding",
        "gates": {"recomputed": True},
        "auditor_revision": "rev",
        "auditor_sha256": m.hashlib.sha256(b"auditor").hexdigest(),
    }
    path.write_text(json.dumps(value))
    monkeypatch.setattr(m.subprocess, "check_output", lambda *a, **k: b"auditor")
    m.verify_calibration_audit(path, m.sha(path), "cal", "binding")
    with pytest.raises(ValueError, match="binding"):
        m.verify_calibration_audit(path, m.sha(path), "changed", "binding")
    value["gates"] = {"recomputed": False}
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="did not pass"):
        m.verify_calibration_audit(path, m.sha(path), "cal", "binding")


def test_snapshot_rejects_changed_committed_scientific_source(tmp_path, monkeypatch):
    names = [
        "experiments/cthyb_spectroscopy.py",
        "experiments/cthyb_dmft.py",
        "experiments/audit_cthyb_spectroscopy.py",
        "kernelgecp/model.py",
    ]
    for name in names:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"committed")
    monkeypatch.setattr(m, "ROOT", tmp_path)
    monkeypatch.setattr(m, "SOURCE_FILES", ())

    def git(args, **kwargs):
        if args[1] == "rev-parse":
            return "revision\n"
        if args[1] == "ls-files":
            return "\n".join(names)
        return b"committed"

    monkeypatch.setattr(m.subprocess, "check_output", git)
    (tmp_path / "kernelgecp/model.py").write_bytes(b"changed")
    with pytest.raises(ValueError, match="scientific source differs"):
        m.snapshot()


def test_case_propagates_covariance_kernel_and_bootstrap(tmp_path, monkeypatch):
    import sys
    from types import ModuleType

    samples = np.random.default_rng(7).normal(size=(64, 3073)) * 1e-3
    monkeypatch.setattr(m, "load_cohort", lambda *args: (samples, {}))

    class Basis:
        beta = 10.0

        def design_tau(self, tau):
            return np.column_stack([np.ones(len(tau)), tau, tau**2])

    monkeypatch.setattr(m, "basis", lambda *args: Basis())
    omega = np.linspace(-4, 4, 5)

    class Spectrum:
        density = np.ones(5) / 8
        alphas = np.array([1.0, 2.0])
        chi2 = np.array([1.0, 2.0])
        trial_spectra = np.stack([density, density])
        selected_index = 0
        converged = (True, True)
        iterations = (1, 1)
        effective_rank = 3
        raw_mass = 1.0
        mass = 1.0
        normalized_training_whitened_rms = 0.1

        def __init__(self):
            self.omega = omega
            self.diagnostics = {"provider": "fixture"}

        def evaluate_tau(self, beta, tau):
            return samples[32:, m.QUERY].mean(axis=0)

    monkeypatch.setattr(m, "continue_tau_maxent", lambda *args: Spectrum())
    provider = ModuleType("triqs_maxent")
    provider.LinearOmegaMesh = lambda *args: omega
    kernels = ModuleType("triqs_maxent.kernels")
    kernels.TauKernel = lambda tau, mesh, beta: SimpleNamespace(
        K=np.outer(tau + 1, mesh + 5)
    )
    monkeypatch.setitem(sys.modules, "triqs_maxent", provider)
    monkeypatch.setitem(sys.modules, "triqs_maxent.kernels", kernels)
    output = tmp_path / "case.json"
    result = m.run(
        tmp_path,
        tmp_path,
        tmp_path,
        "digest",
        "calibration",
        "kernelgecp",
        "flat",
        output,
        bootstrap=2,
        n_omega=5,
        n_alpha=2,
    )
    assert result["status"] == "passed" and result["bootstrap_count"] == 2
    with np.load(output.with_suffix(".npz")) as a:
        np.testing.assert_allclose(
            a["projected_covariance"],
            a["projection"] @ a["covariance"] @ a["projection"].T,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            a["projected_kernel"], a["projection"] @ a["native_kernel"]
        )
        np.testing.assert_allclose(
            a["validation_covariance"], np.cov(samples[32:, m.QUERY], rowvar=False) / 32
        )
        assert a["bootstrap_indices"].shape == (2, 32)


def test_fresh_worker_processes_exit_and_timeout(tmp_path):
    import sys

    pids = []
    for i in range(3):
        log = tmp_path / f"worker-{i}.log"
        record = m.launch_worker(
            [sys.executable, "-c", "print('fixture worker')"],
            log,
            5,
            on_started=pids.append,
        )
        assert record["status"] == "passed" and record["exit_code"] == 0
        assert record["pid"] == pids[-1] and record["log_sha256"] == m.sha(log)
    assert len(set(pids)) == 3
    failed = m.launch_worker(
        [sys.executable, "-c", "raise SystemExit(7)"], tmp_path / "fail.log", 5
    )
    assert failed["status"] == "failed" and failed["exit_code"] == 7
    timed = m.launch_worker(
        [sys.executable, "-c", "import time; time.sleep(5)"],
        tmp_path / "timeout.log",
        0.1,
    )
    assert timed["status"] == "timeout" and timed["exit_code"] != 0
    with pytest.raises(FileExistsError):
        m.launch_worker([sys.executable, "-c", "pass"], tmp_path / "fail.log", 5)


def test_worker_rejects_changed_freeze_before_case(tmp_path, monkeypatch):
    job = {
        "schema": "kernelgecp.cthyb-spectroscopy-worker-job.v1",
        "protocol": m.PROTOCOL,
        "method": "raw",
        "model": "flat",
        "repeat": 0,
        "stage": "final",
        "output": str(tmp_path / "raw-flat-r0.json"),
        "snapshot": {"source": "before"},
        "parent_pid": __import__("os").getppid(),
    }
    path = tmp_path / "job.json"
    path.write_text(json.dumps(job))
    monkeypatch.setattr(m, "snapshot", lambda: {"source": "changed"})
    with pytest.raises(ValueError, match="snapshot mismatch"):
        m.worker(path, m.sha(path))
    with pytest.raises(ValueError, match="job hash"):
        m.worker(path, "0" * 64)


def test_final_requires_freezes_before_any_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "snapshot", lambda: {"source": "fixture"})
    monkeypatch.setattr(
        m, "launch_worker", lambda *a, **kw: pytest.fail("worker must not start")
    )
    with pytest.raises(ValueError, match="requires"):
        m.campaign(tmp_path, tmp_path, tmp_path, "0" * 64, "final", tmp_path / "output")
    assert not (tmp_path / "output").exists()
