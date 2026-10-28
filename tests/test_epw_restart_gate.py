"""Fail-closed tests using a fake process, never substituting for native EPW."""

from pathlib import Path

import pytest

from scripts.providers import epw_restart as restart

INPUT = """&inputepw
 prefix = 'MgB2'
 iverbosity  = 0
 limag = .true.
 lpade = .true.
 wscut = 1.0 ! eV
 temps = 17 20 30 ! K
 /
"""

FAKE = """#!/usr/bin/env python3
import pathlib, sys
p=pathlib.Path('.')
mode=(p/'mode').read_text()
if mode == 'timeout':
    import time
    time.sleep(2)
lines=(p/'MgB2.imag_aniso_017.00').read_text().splitlines()
first=float(lines[1][:20]); last=float(lines[2][:20])
if first != last and mode != 'badcontrol':
    print('Error in routine eliashberg_read_aniso_iaxis: '
          'temperature not the same with the input')
    sys.exit(2)
if mode == 'nonzero':
    print('injected native failure')
    sys.exit(4)
print('Program EPW v.6.0')
print('Pade approximant of full-bandwidth anisotropic')
print('Convergence was reached for N = 5 Pade approximants')
print('raxis_pade   : done')
print('Total program execution')
if mode == 'missing': sys.exit(0)
with (p/'MgB2.pade_aniso_017.00').open('w') as out:
    if mode == 'empty': sys.exit(0)
    out.write('# w [eV] Enk-Ef [eV] Re[znorm(w)] Im[znorm(w)] '
              'Re[delta(w)] Im[delta(w)] Re[shift(w)] Im[shift(w)]\\n')
    for energy in (-1., 1.):
        for iw in range(500):
            if mode == 'truncated' and energy == 1. and iw == 499: break
            vals=[(iw+1)*.002,energy,1.1,.01,.003,.0001,.02,.001]
            if mode == 'budget' and p.resolve().name == 'reconstructed': vals[4]+=.1
            if mode == 'malformed' and iw == 1: out.write('bad\\n'); continue
            out.write(''.join(f'{x:20.10E}' for x in vals)+'\\n')
"""


def fixture(tmp_path, monkeypatch, mode="ok"):
    work = tmp_path / "baseline"
    work.mkdir()
    (work / "epw6.in").write_text(INPUT)
    (work / "mode").write_text(mode)
    (work / "MgB2.a2f").write_text(
        "w[meV] fake\n" + "".join(f"{2 * i}.0 0\n" for i in range(1, 501))
    )
    rows = ["# w [eV] Enk-Ef [eV] znorm(w) delta(w) [eV] shift(w) [eV]\n"]
    for freq in (0.01, 0.03):
        for energy in (-1.0, 1.0):
            rows.append(
                "".join(f"{x:20.10E}" for x in (freq, energy, 1.1, 0.002, 0.001)) + "\n"
            )
    raw = tmp_path / "raw"
    raw.write_text("".join(rows))
    fitted = tmp_path / "fit"
    fitted.write_text("".join(rows))
    qe = tmp_path / "qe"
    (qe / "EPW/bin").mkdir(parents=True)
    (qe / "EPW/src/io").mkdir(parents=True)
    (qe / "test-suite/epw_super").mkdir(parents=True)
    executable = qe / "EPW/bin/epw.x"
    executable.write_text(FAKE)
    executable.chmod(0o755)
    official = qe / "test-suite/epw_super/epw6.in"
    official.write_text(INPUT)
    io = qe / "EPW/src/io/io_supercond.f90"
    io.write_text("pinned source fixture")
    monkeypatch.setattr(restart, "SCOUT_INPUT_SHA256", restart.sha(work / "epw6.in"))
    monkeypatch.setattr(restart, "UPSTREAM_INPUT_SHA256", restart.sha(official))
    monkeypatch.setattr(restart, "UPSTREAM_IO_SHA256", restart.sha(io))
    return work, raw, fitted, qe


def run_fixture(tmp_path, monkeypatch, mode="ok"):
    work, raw, fitted, qe = fixture(tmp_path, monkeypatch, mode)
    result = restart.run(work, raw, fitted, "a" * 40, qe, tmp_path / "receipt.json", 10)
    return result, work, raw


def test_input_changes_exactly_three():
    original = INPUT
    modified = restart.restart_input(original)
    assert modified.count("imag_read = .true.") == 1
    assert "temps = 17 ! K" in modified and "iverbosity  = 2" in modified
    assert (
        modified.replace("temps = 17 ! K", "temps = 17 20 30 ! K")
        .replace("iverbosity  = 2", "iverbosity  = 0")
        .replace("  imag_read = .true.\n", "")
        == original
    )
    with pytest.raises(ValueError):
        restart.restart_input(original + " imag_read = .false.\n")


def test_fake_process_happy_path_and_corruption(tmp_path, monkeypatch):
    result, work, raw = run_fixture(tmp_path, monkeypatch)
    assert result["status"] == "passed", result["errors"]
    assert result["comparison"]["state_count"] == 2
    assert result["comparison"]["frequency_count"] == 500
    assert result["cases"]["wrong-frequency"]["execution"]["returncode"] == 2
    assert result["input"]["original_sha256"] == restart.sha(work / "epw6.in")
    assert len(result["input"]["changes"]) == 3
    for case in result["cases"].values():
        assert case["input_sha256"] == result["input"]["modified_sha256"]
    wrong = Path(result["cases"]["wrong-frequency"]["workdir"]) / restart.RAW_NAME
    a = raw.read_text().splitlines()
    b = wrong.read_text().splitlines()
    assert [i for i, (x, y) in enumerate(zip(a, b, strict=True)) if x != y] == [2]
    assert a[2][20:] == b[2][20:]


@pytest.mark.parametrize(
    "mode",
    ["missing", "empty", "nonzero", "malformed", "truncated", "badcontrol", "budget"],
)
def test_fake_process_fail_closed(tmp_path, monkeypatch, mode):
    result, _, _ = run_fixture(tmp_path, monkeypatch, mode)
    assert result["status"] == "failed" and result["errors"]
    if mode == "budget":
        assert result["comparison"]["channels"]["gap_real_ev"]["failures"] > 0


def test_wrong_source_or_timeout_has_failed_receipt(tmp_path, monkeypatch):
    work, raw, fitted, qe = fixture(tmp_path, monkeypatch)
    result = restart.run(
        work, raw, fitted, "a" * 40, qe, tmp_path / "receipt.json", 1201
    )
    assert result["status"] == "failed" and "timeout" in result["errors"][0]
    result = restart.run(
        work, raw, fitted, "not-a-sha", qe, tmp_path / "receipt.json", 10
    )
    assert result["status"] == "failed" and "revision" in result["errors"][0]


def test_log_success_without_completion_or_bad_negative_fails(tmp_path):
    log = tmp_path / "log"
    log.write_text("Program EPW v.6.0\n")
    execution = {
        "log": str(log),
        "log_bytes": log.stat().st_size,
        "timed_out": False,
        "returncode": 0,
    }
    with pytest.raises(ValueError, match="completion"):
        restart.check_execution(execution, negative=False)
    execution["returncode"] = 1
    with pytest.raises(ValueError, match="error not observed"):
        restart.check_execution(execution, negative=True)


def test_timeout_is_failure_with_durable_log(tmp_path, monkeypatch):
    work, raw, fitted, qe = fixture(tmp_path, monkeypatch, "timeout")
    result = restart.run(work, raw, fitted, "a" * 40, qe, tmp_path / "receipt.json", 1)
    assert result["status"] == "failed"
    assert all(case["execution"]["timed_out"] for case in result["cases"].values())
    assert all(
        Path(case["execution"]["log"]).is_file() for case in result["cases"].values()
    )
