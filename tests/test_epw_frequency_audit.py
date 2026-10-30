"""Independent arithmetic and fail-closed EPW frequency audit tests."""

import json
import subprocess
import sys

import numpy as np
import pytest

from experiments import audit_epw_frequency as audit


def test_direct_product_cardinal_and_condition():
    nodes = np.array([1.0, 2.0])
    poles = np.array([0.0, -1.0])
    weights = np.array([-2.0, 6.0])
    design, condition = audit.direct_design(
        np.array([1.0, 1.5, 2.0]), nodes, poles, weights
    )
    np.testing.assert_array_equal(design[[0, 2]], np.eye(2))
    expected = np.array([4 / 15, 4 / 5])
    np.testing.assert_allclose(design[1], expected)
    assert condition >= 1 and np.all(np.isfinite(design))
    with pytest.raises(ValueError, match="cardinal dimensions"):
        audit.direct_design([1.0], nodes, poles[:1], weights)


def test_native_truth_and_determinant_recovery():
    data = {
        "frequency": np.array([0.1, 0.2]),
        "energy": np.array([-0.3, 0.4]),
        "fields": {
            0: [[1.2, 0.01, 0.03], [1.3, 0.02, -0.01]],
            1: [[1.4, 0.03, 0.02], [1.5, 0.04, 0.01]],
        },
        "context": {"chemical_potential_offset_ev": 0.0},
    }
    g, z, gap, chi = audit.truth(data, [0, 1])
    zz, dd, cc, cond, gain = audit.recover(g, data["frequency"], data["energy"], 0.0)
    np.testing.assert_allclose(zz, z, atol=1e-13)
    np.testing.assert_allclose(dd, gap, atol=1e-13)
    np.testing.assert_allclose(cc, chi, atol=1e-13)
    assert cond >= 1 and gain > 0
    with pytest.raises(ValueError, match="singular"):
        audit.recover(np.zeros_like(g), data["frequency"], data["energy"], 0.0)


def test_calibration_does_not_read_reserved_fields(tmp_path):
    path = tmp_path / "native"
    lines = ["# shift(w)\n"]
    for i in range(54):
        row = "".join(f"{v:20.10E}" for v in (i + 0.5, -1.0, 1.0, 0.01, 0.02))
        if i in audit.FINAL:
            row = row[:40] + ("not a number".rjust(20)) * 3
        lines.append(row + "\n")
    path.write_text("".join(lines))
    allowed = [x for x in range(54) if x not in audit.FINAL]
    frequencies, energies, fields = audit.parse_native_file(path, 54, 1, allowed)
    assert len(frequencies) == 54 and energies == [-1.0] and set(fields) == set(allowed)
    with pytest.raises(ValueError):
        audit.parse_native_file(path, 54, 1, [45])
    path.write_text("".join(lines[:-1]))
    with pytest.raises(ValueError, match="truncated"):
        audit.parse_native_file(path, 54, 1, allowed)


def test_native_state_order_and_nonfinite_fail(tmp_path):
    path = tmp_path / "native"

    def row(*v):
        return "".join(f"{x:20.10E}" for x in v) + "\n"

    path.write_text(
        "# shift(w)\n" + row(1.0, 0.0, 1.0, 0.0, 0.0) + row(2.0, 1.0, 1.0, 0.0, 0.0)
    )
    with pytest.raises(ValueError, match="ordering"):
        audit.parse_native_file(path, 2, 1, [0, 1])
    path.write_text("# shift(w)\n" + row(1.0, 0.0, float("nan"), 0.0, 0.0))
    with pytest.raises(ValueError, match="field invalid"):
        audit.parse_native_file(path, 1, 1, [0])


def test_missing_and_invalid_archive_fail(tmp_path):
    p = tmp_path / "empty.npz"
    np.savez_compressed(p)
    with pytest.raises(ValueError, match="empty"):
        audit.npz(p)
    np.savez_compressed(p, a=np.array([np.nan]))
    with pytest.raises(ValueError, match="nonfinite"):
        audit.npz(p)
    p.write_bytes(b"not an archive")
    with pytest.raises(ValueError):
        audit.npz(p)


def test_reject_noncanonical_scout(tmp_path):
    cal = tmp_path / "cal.json"
    freeze = tmp_path / "freeze.json"
    cal.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.epw-frequency-calibration.v1",
                "status": "passed",
                "pilot_noncanonical": True,
            }
        )
    )
    freeze.write_text(
        json.dumps(
            {
                "schema": "kernelgecp.epw-frequency-freeze.v1",
                "status": "passed",
                "pilot_noncanonical": True,
            }
        )
    )
    with pytest.raises(ValueError, match="noncanonical"):
        audit.audit_calibration(cal, freeze, tmp_path)


def test_failed_cli_receipt(tmp_path):
    freeze = tmp_path / "freeze.json"
    freeze.write_text("{}")
    out = tmp_path / "audit.json"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "experiments.audit_epw_frequency",
            "--calibration",
            str(tmp_path / "missing.json"),
            "--freeze",
            str(freeze),
            "--native-root",
            str(tmp_path),
            "--output",
            str(out),
        ],
        check=False,
    )
    assert proc.returncode != 0 and json.loads(out.read_text())["status"] == "failed"


def test_candidate_inventory_is_ordered_and_failure_closed(tmp_path, monkeypatch):
    cal = tmp_path / "cal.json"
    freeze = tmp_path / "freeze.json"
    parts = audit.partitions()
    p = {"revision": "fake"}
    ctx = {}
    calibration = {
        "schema": "kernelgecp.epw-frequency-calibration.v1",
        "status": "passed",
        "pilot_noncanonical": False,
        "provenance": p,
        "context": ctx,
        "partitions": parts,
        "selected_rank": 60,
        "attempts": [{"rank": 60, "passed": True}],
    }
    cal.write_text(json.dumps(calibration))
    native = tmp_path / "native"
    native.write_bytes(b"fixture")
    f = {
        "schema": "kernelgecp.epw-frequency-freeze.v1",
        "status": "passed",
        "pilot_noncanonical": False,
        "provenance": p,
        "context": ctx,
        "partitions": parts,
        "selected_rank": 60,
        "calibration_sha256": audit.sha(cal),
        "rank_ladder": audit.RANKS,
        "budget": audit.BUDGET,
        "controls": audit.CONTROLS,
        "calibration_margin": 0.1,
        "native_sha256": audit.sha(native),
    }
    freeze.write_text(json.dumps(f))
    monkeypatch.setattr(audit, "provenance", lambda p: None)
    monkeypatch.setattr(audit, "native", lambda *x: {"path": native})
    with pytest.raises(ValueError, match="candidate inventory"):
        audit.audit_calibration(cal, freeze, tmp_path)
    calibration["attempts"] = [
        {"rank": 48, "passed": False, "error": "unexplained"},
        {"rank": 60, "passed": True},
    ]
    cal.write_text(json.dumps(calibration))
    f["calibration_sha256"] = audit.sha(cal)
    freeze.write_text(json.dumps(f))
    with pytest.raises(ValueError, match="unresolved candidate"):
        audit.audit_calibration(cal, freeze, tmp_path)
