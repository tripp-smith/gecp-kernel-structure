"""Synthetic and fake-process protocol checks; no native EPW acceptance."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from test_epw_frequency import _static_data

from experiments import epw_restart_roundtrip as driver
from kernelgecp.integrations.epw import EPWAnisotropicData


def test_bitwise_native_comparison_detects_signed_zero():
    data = _static_data()
    driver.same_native(data, data)
    chi = data.chi_ev.copy()
    chi[0, 0] = -0.0 if chi[0, 0] == 0 else chi[0, 0] + 1
    changed = EPWAnisotropicData(
        data.frequencies_ev,
        data.state_energies_ev,
        data.z,
        data.delta_ev,
        chi,
        context=data.context,
    )
    with pytest.raises(ValueError, match="bits mismatch"):
        driver.same_native(data, changed)


def test_controls_holdout_and_explicit_boundary():
    assert driver.HOLDOUT == (45, 49, 53)
    assert driver.CONTROLS["rank"] == 60
    assert driver.CONTROLS["real_nambu_parity"] is True
    assert driver.CONTROLS["policy"] == "lossless_native_xor_v1"
    assert all(
        i not in [j for j in range(54) if j not in driver.HOLDOUT]
        for i in driver.HOLDOUT
    )


def test_fake_process_campaign_frozen_first_and_repeat_comparison(
    tmp_path, monkeypatch
):
    env = {"source": "synthetic"}
    context = {"synthetic": True}
    monkeypatch.setattr(driver, "source_environment", lambda: env)
    monkeypatch.setattr(driver, "native_context", lambda *args: (context, None, None))

    def fake_run(cmd, **kwargs):
        out = __import__("pathlib").Path(cmd[cmd.index("--output") + 1])
        frozen = __import__("pathlib").Path(cmd[cmd.index("--freeze") + 1])
        assert frozen.exists()
        freeze = json.loads(frozen.read_text())
        assert freeze["controls"] == driver.CONTROLS and freeze["environment"] == env
        out.mkdir()
        np.savez_compressed(
            out / "scientific-arrays.npz", synthetic=np.array([1.0, 2.0])
        )
        driver.write(
            out / "case.json",
            {
                "status": "passed",
                "repeat": int(cmd[-1]),
                "freeze_sha256": driver.sha(frozen),
                "files_sha256": {
                    "scientific-arrays.npz": driver.sha(out / "scientific-arrays.npz")
                },
            },
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(driver.subprocess, "run", fake_run)
    result = driver.campaign(tmp_path / "synthetic", 20, 0.0, tmp_path / "out")
    assert result["status"] == "passed"
    assert result["repeat_arrays_identical"] is True
    assert result["native_acceptance"] == "not_tested"
    assert len(result["children"]) == 3
    with pytest.raises(ValueError, match="20 K"):
        driver.campaign(tmp_path, 17, 0.0, tmp_path / "other")


def test_synthetic_receipt_extracts_second_temperature_and_rejects_wrong_offset(
    tmp_path, monkeypatch
):
    work = tmp_path / "epw-native/test-suite/epw_super"
    work.mkdir(parents=True)
    (work / "epw6.in").write_text("synthetic input")
    (work / "epw6.out").write_text(
        "Fermi level (eV) = 1.000E+00\nChemical potential (itemp = 2) = 1.125E+00\n"
        + "Convergence was reached in nsiter\n" * 3
    )
    monkeypatch.setattr(driver, "UPSTREAM_INPUT_SHA256", "input-pin")
    monkeypatch.setattr(driver, "UPSTREAM_IO_SHA256", "io-pin")
    monkeypatch.setattr(driver, "SCOUT_INPUT_SHA256", driver.sha(work / "epw6.in"))
    files = []
    for temperature in (17, 20, 30):
        name = f"MgB2.imag_aniso_{temperature:06.2f}"
        (work / name).write_text("synthetic placeholder; no provider data")
        files.append(
            {
                "temperature_kelvin": temperature,
                "filename": name,
                "sha256": driver.sha(work / name),
                "bytes": (work / name).stat().st_size,
                "columns": [
                    "frequency_ev",
                    "energy_minus_ef0_ev",
                    "Z",
                    "gap_ev",
                    "shift_ev",
                ],
                "positive_only": True,
                "state_ordering": "native_frequency_major_then_ik_then_ibnd_ordinal",
                "state_count": 1,
                "frequency_count": 1,
            }
        )
    digest = driver.sha(work / "epw6.in")
    receipt = {
        "schema": "kernelgecp.epw-native-anisotropic-fixture.v1",
        "status": "passed",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "wscut_ev": 1.0,
        "script_sha256": driver.sha(
            driver.ROOT / "scripts/providers/epw_anisotropic_fixture.py"
        ),
        "baseline": {
            "checks": {"synthetic": True},
            "native_executions": [
                {"input": "epw1.in", "executable_sha256": "synthetic-binary"}
            ],
        },
        "native_anisotropic_execution": {
            "returncode": 0,
            "executable_sha256": "synthetic-binary",
            "input_sha256": digest,
            "log_sha256": driver.sha(work / "epw6.out"),
        },
        "executed_input_sha256": digest,
        "upstream_input_sha256": "input-pin",
        "upstream_io_sha256": "io-pin",
        "source_revision": "synthetic",
        "files": files,
    }
    driver.write(tmp_path / "epw-anisotropic.json", receipt)
    context, path, _ = driver.native_context(tmp_path, 0.125, 20)
    assert context["temperature_kelvin"] == 20
    assert context["reference_fermi_energy_ev"] == 1.0
    assert path.name == "MgB2.imag_aniso_020.00"
    with pytest.raises(ValueError, match="chemical offset"):
        driver.native_context(tmp_path, 0.0, 20)
    (work / files[1]["filename"]).write_text("tampered")
    with pytest.raises(ValueError, match="hash mismatch"):
        driver.native_context(tmp_path, 0.125, 20)


def test_pure_failure_is_recorded_before_any_correction(tmp_path, monkeypatch):
    env, context = {"synthetic": "source"}, {"synthetic": "context"}
    monkeypatch.setattr(driver, "source_environment", lambda: env)
    monkeypatch.setattr(
        driver, "native_context", lambda *args: (context, tmp_path / "synthetic", {})
    )
    monkeypatch.setattr(
        driver, "read_checked", lambda *args: SimpleNamespace(indices=np.arange(54))
    )
    seen = {}

    def fake_fit(data, rows, rank):
        seen["rows"] = rows
        return SimpleNamespace(diagnostics={"synthetic": True})

    monkeypatch.setattr(driver, "fit", fake_fit)
    monkeypatch.setattr(
        driver,
        "measure",
        lambda *args: (
            {
                "green_scaled": 1.0,
                "z": 0.0,
                "gap_ev": 0.0,
                "chi_ev": 0.0,
                "imaginary_defect": 0.0,
            },
            {},
        ),
    )
    monkeypatch.setattr(
        driver.EPWRestartArchive,
        "from_data",
        lambda *args, **kwargs: pytest.fail("corrections created after pure failure"),
    )
    freeze = tmp_path / "freeze.json"
    driver.write(
        freeze, {"environment": env, "context": context, "controls": driver.CONTROLS}
    )
    out = tmp_path / "case"
    with pytest.raises(ValueError, match="pure-model heldout"):
        driver._case(tmp_path, 0.0, out, 0, freeze)
    result = json.loads((out / "case.json").read_text())
    assert result["status"] == "failed" and result["pure_model_heldout_passed"] is False
    assert not set(driver.HOLDOUT).intersection(seen["rows"])
    assert not (out / "bundle.json").exists()
