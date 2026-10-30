"""Deterministic tests with a fake solver; these do not run native science."""

import json

import numpy as np
import pytest

from kernelgecp.integrations.epw import (
    KB_EV_PER_K,
    EPWAnisotropicData,
    write_epw_anisotropic,
)
from scripts.providers import epw_pade_calibration as calibration
from scripts.providers import epw_restart as epw
from scripts.providers import epw_restart_sensitivity as sensitivity

INPUT = """&inputepw
 iverbosity = 0
 limag = .true.
 wscut = 1.0 ! eV
 temps = 17 20 30 ! K
 /
"""


def fixture(tmp_path, monkeypatch):
    artifact, qe = tmp_path / "artifact", tmp_path / "qe"
    campaign = artifact / "campaign/case-0"
    workspace = artifact / "epw-native/test-suite/epw_super"
    old = artifact / "restart-comparison-restart-cases/raw" / epw.OUTPUT_NAME
    for directory in (
        campaign,
        workspace,
        old.parent,
        qe / "EPW/bin",
        qe / "EPW/src/io",
        qe / "test-suite/epw_super",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    context = {
        "schema": "kernelgecp.epw-anisotropic-input.v1",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "temperature_kelvin": 17,
        "units": "eV",
        "state_order": "frequency_major_ik_ibnd_ordinal",
        "state_count": 2,
        "matsubara_indices": [0, 1],
        "positive_only": True,
        "full_bandwidth": True,
        "chi_convention": "native_fbw_column",
        "reference_fermi_energy_ev": 3.2,
        "chemical_potential_offset_ev": 0.03,
    }
    frequencies = (2 * np.arange(2) + 1) * np.pi * KB_EV_PER_K * 17
    a = EPWAnisotropicData(
        frequencies,
        [-1, 1],
        np.ones((2, 2)),
        np.full((2, 2), 0.02),
        np.full((2, 2), 0.01),
        context=context,
    )
    b = a.replace_fields(a.z + 0.01, a.delta_ev + 0.001, a.chi_ev + 0.0001)
    raw, fit = campaign / "original-replay.dat", campaign / "fitted-native.dat"
    write_epw_anisotropic(a, raw)
    write_epw_anisotropic(b, fit)
    (campaign / "case.json").write_text(json.dumps({"fitted_context": context}))
    original, official, io, binary = (
        workspace / epw.INPUT_NAME,
        qe / "test-suite/epw_super/epw6.in",
        qe / "EPW/src/io/io_supercond.f90",
        qe / "EPW/bin/epw.x",
    )
    for path, content in (
        (original, INPUT),
        (official, INPUT),
        (io, "source"),
        (binary, "fake"),
    ):
        path.write_text(content)
    a2f = workspace / "MgB2.a2f"
    a2f.write_text("w[meV]\n" + "".join(f"{2 * i} 0\n" for i in range(1, 501)))
    header = (
        "# "
        + " ".join(
            (
                "Re[znorm(w)]",
                "Im[znorm(w)]",
                "Re[delta(w)]",
                "Im[delta(w)]",
                "Re[shift(w)]",
                "Im[shift(w)]",
            )
        )
        + "\n"
    )
    old.write_text(
        header
        + "".join(
            "".join(
                f"{x:20.10E}"
                for x in ((i + 1) * 0.002, state, 1, 0.01, 0.02, 0.001, 0.03, 0.002)
            )
            + "\n"
            for state in (-1, 1)
            for i in range(500)
        )
    )
    for module, name, path in (
        (calibration, "RAW", raw),
        (calibration, "FITTED", fit),
        (calibration, "BINARY", binary),
        (epw, "SCOUT_INPUT_SHA256", original),
        (epw, "UPSTREAM_INPUT_SHA256", official),
        (epw, "UPSTREAM_IO_SHA256", io),
    ):
        monkeypatch.setattr(module, name, epw.sha(path))
    saved = {
        "provider": {"qe": "7.5", "epw": "6.0"},
        "source": {"binary_sha256": epw.sha(binary)},
        "inputs": {
            "raw": {"sha256": epw.sha(raw)},
            "reconstructed": {"sha256": epw.sha(fit)},
        },
        "cases": {"raw": {"output": {"sha256": epw.sha(old)}}},
        "expected_real_grid": {"a2f_sha256": epw.sha(a2f)},
    }
    (artifact / "restart-comparison.json").write_text(json.dumps(saved))

    def fake_execute(_binary, work, timeout):
        assert 1 <= timeout <= 1200
        log = work / "restart.out"
        log.write_text(
            "EPW v.6.0\nPade approximant of full-bandwidth anisotropic\n"
            "Convergence was reached for N = 98\nraxis_pade\n"
            "Total program execution\n"
        )
        (work / epw.OUTPUT_NAME).write_bytes(old.read_bytes())
        return {
            "log": str(log),
            "log_bytes": log.stat().st_size,
            "returncode": 0,
            "timed_out": False,
            "wall_seconds": 0,
        }

    monkeypatch.setattr(epw, "execute", fake_execute)
    return artifact, qe, raw, fit, context


def test_formatted_perturbation_and_zero_control(tmp_path, monkeypatch):
    _, _, raw, fit, context = fixture(tmp_path, monkeypatch)
    zero = tmp_path / "zero"
    result = sensitivity.prepare(raw, fit, context, zero, 0)
    assert zero.read_bytes() == raw.read_bytes()
    assert result["fields"]["Z"]["equal_text_count"] == 4
    result = sensitivity.prepare(raw, fit, context, tmp_path / "scaled", 0.1)
    assert result["frequency_and_state_text_unchanged"]
    assert result["fields"]["gap_ev"]["max_absolute_error"] == pytest.approx(0.0001)
    assert result["fields"]["gap_ev"]["equal_value_count"] == 0
    with pytest.raises(ValueError, match="predeclared"):
        sensitivity.prepare(raw, fit, context, tmp_path / "invalid", 0.2)


def test_full_mocked_run_and_receipt(tmp_path, monkeypatch):
    artifact, qe, *_ = fixture(tmp_path, monkeypatch)
    out = tmp_path / "run"
    result = sensitivity.run(artifact, qe, out, "a" * 40)
    assert result["status"] == "completed", result["errors"]
    assert len(result["cases"]) == 6
    assert all(
        result["cases"][f"scale-{x:g}"]["comparison"]["passed"]
        for x in sensitivity.SCALES
    )
    assert len(result["cases"]["scale-0"]["comparison"]["channels"]) == 6
    assert json.loads((out / "sensitivity-progress.json").read_text()) == result


def test_failure_receipt_and_budgets(tmp_path, monkeypatch):
    artifact, qe, raw, _, _ = fixture(tmp_path, monkeypatch)
    raw.write_text("corrupt")
    out = tmp_path / "bad"
    result = sensitivity.run(artifact, qe, out, "a" * 40)
    assert (
        result["status"] == "inconclusive"
        and "input hash mismatch" in result["errors"][0]
    )
    assert (out / "sensitivity-progress.json").exists()
    result = sensitivity.run(
        artifact, qe, tmp_path / "budget", "a" * 40, budget_seconds=18001
    )
    assert "budget exceeds" in result["errors"][0]
    with pytest.raises(TimeoutError):
        sensitivity.remaining(0)


def test_malformed_input_fails_closed(tmp_path, monkeypatch):
    _, _, raw, fit, context = fixture(tmp_path, monkeypatch)
    fit.write_text(fit.read_text().replace("E+00", "E+xx", 1))
    with pytest.raises(ValueError):
        sensitivity.prepare(raw, fit, context, tmp_path / "out", 1)


def test_nonzero_solver_exit_is_recorded(tmp_path, monkeypatch):
    artifact, qe, *_ = fixture(tmp_path, monkeypatch)
    execute = epw.execute

    def failed(binary, work, timeout):
        status = execute(binary, work, timeout)
        status["returncode"] = 7
        return status

    monkeypatch.setattr(epw, "execute", failed)
    result = sensitivity.run(artifact, qe, tmp_path / "failed", "a" * 40)
    assert result["status"] == "inconclusive"
    assert result["cases"]["raw"]["execution"]["returncode"] == 7
    assert result["cases"]["scale-0"]["status"] == "not_run"


def test_scale_one_preserves_axes_and_matches_fitted_fields(tmp_path, monkeypatch):
    _, _, raw, fit, context = fixture(tmp_path, monkeypatch)
    target = tmp_path / "scale-one"
    sensitivity.prepare(raw, fit, context, target, 1)
    raw_lines = raw.read_bytes().splitlines()
    fit_lines = fit.read_bytes().splitlines()
    result_lines = target.read_bytes().splitlines()
    assert result_lines[0] == raw_lines[0]
    assert all(
        a[:40] == b[:40] for a, b in zip(raw_lines[1:], result_lines[1:], strict=True)
    )
    assert all(
        a[40:] == b[40:] for a, b in zip(fit_lines[1:], result_lines[1:], strict=True)
    )


def test_axis_mismatch_fails_before_writing(tmp_path, monkeypatch):
    _, _, raw, fit, context = fixture(tmp_path, monkeypatch)
    lines = fit.read_text().splitlines(keepends=True)
    for i in (1, 3):
        lines[i] = lines[i][:20] + f"{-2.0:20.10E}" + lines[i][40:]
    fit.write_text("".join(lines))
    target = tmp_path / "bad-axes"
    with pytest.raises(ValueError, match="axes mismatch"):
        sensitivity.prepare(raw, fit, context, target, 0.1)
    assert not target.exists()


def test_deadline_caps_timeout(monkeypatch):
    monkeypatch.setattr(sensitivity.time, "monotonic", lambda: 10.0)
    assert sensitivity.remaining(12.9) == 2
    assert sensitivity.remaining(5010) == 1200
    with pytest.raises(TimeoutError):
        sensitivity.remaining(10.9)
