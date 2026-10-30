"""Fast integrity tests; no native acquisition or external provider required."""

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from experiments import cthyb_fixed_bath as driver
from kernelgecp.integrations.cthyb_dmft import CTHYBConfig, CTHYBMeasurement


@pytest.mark.skipif(
    sys.platform != "linux" or importlib.util.find_spec("triqs_cthyb") is None,
    reason="canonical Linux acquisition-runtime preflight requires actual CT-HYB",
)
def test_linux_native_runtime_freeze_survives_first_solve():
    config = CTHYBConfig(n_iw=32, n_tau=385, warmup_cycles=500)
    plan = driver.CTHYBPlan(config)
    before = driver.provenance()
    # A separate disposable wiring probe, outside every declared final seed family.
    measured = plan.sample(plan.initial_bath(), cycles=2000, seed=32290001)
    assert measured.solve_status == 0
    assert driver.provenance() == before


def test_parent_uses_original_audit_receipt_binding_key(monkeypatch, tmp_path):
    from types import SimpleNamespace

    parent = tmp_path / "parent"
    (parent / "dense-final").mkdir(parents=True)
    config = CTHYBConfig()
    manifest = {"config": driver.dataclasses.asdict(config)}
    files = {
        "manifest.json": manifest,
        "independent-audit.json": {"status": "passed"},
        "freeze.json": {},
        "result.json": {},
    }
    for name, data in files.items():
        (parent / name).write_text(json.dumps(data))
    stem = parent / "dense-final/replica-00"
    stem.with_suffix(".npz").write_bytes(b"hash-bound native fixture placeholder")
    stem.with_suffix(".json").write_text(
        json.dumps(
            {
                "identity": manifest,
                "seed": 30130001,
                "arrays_sha256": driver._sha(stem.with_suffix(".npz")),
            }
        )
    )
    binding = {
        "schema": "kernelgecp.cthyb-audit-binding.v1",
        "status": "passed",
        "generation_revision": driver.REVISION,
        "audit_receipt_sha256": driver._sha(parent / "independent-audit.json"),
        "manifest_sha256": driver._sha(parent / "manifest.json"),
        "freeze_sha256": driver._sha(parent / "freeze.json"),
        "result_sha256": driver._sha(parent / "result.json"),
        "measurement_archives": {
            "dense-final/replica-00.json": driver._sha(stem.with_suffix(".json"))
        },
    }
    path = tmp_path / "binding.json"
    path.write_text(json.dumps(binding))
    monkeypatch.setattr(driver, "BINDING_SHA", driver._sha(path))
    monkeypatch.setattr(
        driver,
        "load_measurement",
        lambda *args: SimpleNamespace(
            config=config,
            cycles=1600000,
            bath_iw=np.zeros(2 * config.n_iw, dtype=complex),
        ),
    )
    info, _, actual = driver.parent_data(parent, path)
    assert actual == config
    assert info["audit_sha256"] == binding["audit_receipt_sha256"]


def calargs():
    return (Path("cal.json"), "a" * 64, Path("audit.json"), "b" * 64)


def fixture_setup(monkeypatch, tmp_path):
    config = CTHYBConfig(n_iw=8, n_tau=49)
    bath = np.zeros(16, dtype=complex)
    parent = tmp_path / "parent"
    binding = tmp_path / "binding"
    output = tmp_path / "cohort"
    monkeypatch.setattr(
        driver, "parent_data", lambda p, b: ({"binding_sha256": "parent"}, bath, config)
    )
    monkeypatch.setattr(driver, "provenance", lambda: {"snapshot": "frozen"})
    monkeypatch.setattr(
        driver,
        "calibration_binding",
        lambda *args: {"receipt_sha256": "frozen-calibration"},
    )

    class Plan:
        def __init__(self, config):
            self.config = config

        def sample(self, b, *, cycles, seed):
            return CTHYBMeasurement(
                self.config,
                seed,
                cycles,
                b,
                np.zeros((49, 2), complex),
                np.zeros((16, 2), complex),
                np.zeros((16, 2), complex),
                np.zeros((16, 2), complex),
                np.array([[0, 0], [1, 1], [0, 0]], complex),
                np.zeros((2, 2), complex),
                np.array([0.5, 0.5]),
                0.2,
                1.0,
                0,
            )

    monkeypatch.setattr(driver, "CTHYBPlan", Plan)
    driver.freeze(parent, binding, output, *calargs())
    return parent, binding, output


def run_worker(parent, binding, output):
    lock = (output / ".acquire.lock").open("a")
    try:
        driver.worker(output, parent, binding, *calargs(), 0, lock.fileno())
    finally:
        lock.close()


def test_freeze_and_resume(monkeypatch, tmp_path):
    parent, binding, output = fixture_setup(monkeypatch, tmp_path)
    run_worker(parent, binding, output)
    manifest, result, _, _ = driver.validate(output, parent, binding, *calargs())
    assert manifest["protocol"]["seeds"] == [32100001 + 2 * i for i in range(256)]
    assert result["chains"][0]["role"] == "train"
    with pytest.raises(FileExistsError):
        driver.freeze(parent, binding, output, *calargs())
    with pytest.raises(ValueError):
        run_worker(parent, binding, output)


@pytest.mark.parametrize(
    "name", ["chain-001.npz", "chain-001.timing.json", "chain-001.attempt"]
)
def test_orphan_rejected(monkeypatch, tmp_path, name):
    parent, binding, output = fixture_setup(monkeypatch, tmp_path)
    (output / "chains").mkdir()
    (output / "chains" / name).touch()
    with pytest.raises(ValueError, match="orphan"):
        driver.validate(output, parent, binding, *calargs())


def test_tampered_archive_and_source_rejected(monkeypatch, tmp_path):
    parent, binding, output = fixture_setup(monkeypatch, tmp_path)
    run_worker(parent, binding, output)
    archive = output / "chains/chain-000.npz"
    archive.write_bytes(archive.read_bytes() + b"bad")
    with pytest.raises(ValueError, match="hash"):
        driver.validate(output, parent, binding, *calargs())
    monkeypatch.setattr(driver, "provenance", lambda: {"snapshot": "changed"})
    with pytest.raises(ValueError, match="changed"):
        driver.validate(output, parent, binding, *calargs())


def test_failed_chain_retains_attempt(monkeypatch, tmp_path):
    parent, binding, output = fixture_setup(monkeypatch, tmp_path)

    class Failing:
        def __init__(self, config):
            pass

        def sample(self, *args, **kwargs):
            raise RuntimeError("native failure")

    monkeypatch.setattr(driver, "CTHYBPlan", Failing)
    with pytest.raises(RuntimeError, match="native failure"):
        run_worker(parent, binding, output)
    with pytest.raises(ValueError, match="orphan"):
        driver.validate(output, parent, binding, *calargs())


def calibration_fixture(monkeypatch, tmp_path):
    spec = driver.spectroscopy
    auditor_blob = b"independent auditor"
    auditor_hash = hashlib.sha256(auditor_blob).hexdigest()
    snapshot = {
        "revision": "revision",
        "sources": {"experiments/audit_cthyb_spectroscopy.py": auditor_hash},
    }
    monkeypatch.setattr(spec, "snapshot", lambda: snapshot)
    monkeypatch.setattr(spec, "validate_basis_controls", lambda: None)
    monkeypatch.setattr(
        spec.subprocess, "check_output", lambda *args, **kwargs: auditor_blob
    )
    binding = "b" * 64
    entries = {}
    for method in spec.METHODS:
        for model in spec.PROTOCOL["default_models"]:
            name = f"{method}-{model}-r0.json"
            case = tmp_path / name
            archive = case.with_suffix(".npz")
            archive.write_bytes(b"fake case archive")
            case.write_text(
                json.dumps(
                    {
                        "status": "passed",
                        "array_sha256": driver._sha(archive),
                        "eps": spec.PROTOCOL["eps"],
                        "rcond": spec.PROTOCOL["rcond"],
                        "config": {"minimizer_mu0": spec.PROTOCOL["minimizer_mu0"]},
                    }
                )
            )
            entries[name] = {
                "result_sha256": driver._sha(case),
                "array_sha256": driver._sha(archive),
            }
    cohort = {
        "stage": "calibration",
        "role": "kernelgecp-final",
        "generation_revision": spec.CALIBRATION_REVISION,
        "seed_base": spec.CALIBRATION_SEED_BASE,
        "seeds": [spec.CALIBRATION_SEED_BASE + 40001 + 2 * i for i in range(64)],
        "binding_sha256": binding,
    }
    receipt = tmp_path / "calibration.json"
    receipt.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.cthyb-spectroscopy-calibration.v1",
                "status": "passed",
                "protocol": spec.PROTOCOL,
                "snapshot": snapshot,
                "binding_sha256": binding,
                "cohort": cohort,
                "cases": entries,
            }
        )
    )
    digest = driver._sha(receipt)
    audit = tmp_path / "audit.json"
    audit.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.cthyb-spectroscopy-independent-audit.v1",
                "status": "passed",
                "kind": "calibration",
                "calibration_sha256": digest,
                "campaign_sha256": digest,
                "binding_sha256": binding,
                "auditor_revision": "revision",
                "auditor_sha256": auditor_hash,
                "gates": {g: True for g in driver.CALIBRATION_GATES},
                "cases": 8,
            }
        )
    )
    return receipt, digest, audit, driver._sha(audit)


def test_calibration_binding_and_stale_receipt(monkeypatch, tmp_path):
    args = calibration_fixture(monkeypatch, tmp_path)
    result = driver.calibration_binding(*args)
    assert result["receipt_sha256"] == args[1]
    args[0].write_text(args[0].read_text() + " ")
    with pytest.raises(ValueError, match="pinned"):
        driver.calibration_binding(*args)


@pytest.mark.parametrize(
    "mutation", ["failed", "missing_gate", "stale_auditor", "wrong_count"]
)
def test_calibration_audit_rejected(monkeypatch, tmp_path, mutation):
    receipt, digest, audit, _ = calibration_fixture(monkeypatch, tmp_path)
    data = json.loads(audit.read_text())
    if mutation == "failed":
        data["status"] = "failed"
    if mutation == "missing_gate":
        data["gates"].pop("heldout_rms")
    if mutation == "stale_auditor":
        data["auditor_revision"] = "old"
    if mutation == "wrong_count":
        data["cases"] = 7
    audit.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        driver.calibration_binding(receipt, digest, audit, driver._sha(audit))


def test_calibration_missing(monkeypatch, tmp_path):
    args = calibration_fixture(monkeypatch, tmp_path)
    args[0].unlink()
    with pytest.raises((FileNotFoundError, ValueError)):
        driver.calibration_binding(*args)
