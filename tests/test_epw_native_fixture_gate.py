"""Fail-closed checks for the first authentic EPW fixture-generation gate."""

from pathlib import Path

import pytest

from scripts.providers.epw_anisotropic_fixture import inspect_native_file


def _write(path: Path, rows: list[list[float]]) -> Path:
    header = "# w [eV] Enk-Ef [eV] znorm(w) delta(w) [eV] shift(w) [eV]\n"
    path.write_text(
        header + "".join("".join(f"{x:20.10E}" for x in row) + "\n" for row in rows)
    )
    return path


def test_native_gate_preserves_degenerate_state_ordinals(tmp_path):
    rows = [[w, -0.2, 1.5, gap, 0.01] for w in (0.01, 0.03) for gap in (0.002, 0.004)]
    result = inspect_native_file(_write(tmp_path / "native", rows))
    assert result["frequency_count"] == result["state_count"] == 2
    assert result["state_energy_minus_ef0_ev"] == [-0.2, -0.2]
    assert result["positive_only"] is True
    assert len(result["sha256"]) == 64


@pytest.mark.parametrize(
    "bad",
    ["nan", "zero_frequency", "zero_Z", "repeated_block", "truncated", "state_order"],
)
def test_malformed_or_unsupported_native_file_fails(tmp_path, bad):
    rows = [[w, e, 1.5, 0.004, 0.01] for w in (0.01, 0.03) for e in (-0.2, 0.3)]
    if bad == "nan":
        rows[1][3] = float("nan")
    elif bad == "zero_frequency":
        rows[0][0] = 0
    elif bad == "zero_Z":
        rows[0][2] = 0
    elif bad == "repeated_block":
        rows.extend(rows[:2])
    elif bad == "truncated":
        rows.pop()
    elif bad == "state_order":
        rows[-1][1], rows[-2][1] = rows[-2][1], rows[-1][1]
    with pytest.raises(ValueError):
        inspect_native_file(_write(tmp_path / "bad", rows))


def test_unrounded_or_short_fields_and_wrong_header_fail(tmp_path):
    source = tmp_path / "bad"
    for text in ("native\n0.01 0.0 1.0 0.1 0.0\n", "shift(w)\n0.01 0 1 0.1 0\n", ""):
        source.write_text(text)
        with pytest.raises(ValueError):
            inspect_native_file(source)
