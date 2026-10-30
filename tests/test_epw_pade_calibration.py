"""Cheap fail-closed tests; no fake process is native acceptance."""

import pytest

from scripts.providers import epw_pade_calibration as m


def test_order_mapping_and_nontrivial_orders():
    assert tuple(map(m.order, m.PERCENTAGES)) == m.ORDERS
    assert min(m.ORDERS) == 16


def test_resource_exhaustion_is_explicit(monkeypatch):
    monkeypatch.setattr(m, "_deadline", m.time.monotonic() - 1)
    with pytest.raises(TimeoutError, match="budget exhausted"):
        m.remaining_timeout()


def test_extract_requires_exactly_one_of_each():
    text = (
        "SUBROUTINE pade_coeff(N)\nEND SUBROUTINE pade_coeff\n"
        "SUBROUTINE pade_eval(N)\nEND SUBROUTINE pade_eval\n"
    )
    assert "pade_eval" in m.extract(text)
    with pytest.raises(ValueError):
        m.extract(text + "\n" + text)
    with pytest.raises(ValueError):
        m.extract("SUBROUTINE pade_eval(N)\nEND SUBROUTINE pade_eval\n")


def test_bad_identity_fails_closed(tmp_path):
    prior = tmp_path / "prior"
    prior.mkdir()
    qe = tmp_path / "qe"
    (qe / "EPW/bin").mkdir(parents=True)
    (qe / "EPW/bin/epw.x").write_bytes(b"not pinned")
    (prior / "restart-comparison.json").write_text(
        '{"source":{"binary_sha256":"wrong"}}'
    )
    result = m.run(prior, qe, tmp_path / "out", "a" * 40)
    assert (
        result["status"] == "inconclusive" and "binary mismatch" in result["errors"][0]
    )


def test_padded_header_keeps_native_checks(tmp_path):
    required = (
        "# w Re[znorm(w)] Im[znorm(w)] Re[delta(w)] Im[delta(w)] "
        "Re[shift(w)] Im[shift(w)]\n"
    )
    row = (
        "".join(f"{v:20.10E}" for v in (0.2, 1.0, 1.0, 0.0, 0.1, 0.0, 0.2, 0.0)) + "\n"
    )
    row2 = (
        "".join(f"{v:20.10E}" for v in (0.4, 1.0, 1.0, 0.0, 0.1, 0.0, 0.2, 0.0)) + "\n"
    )
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.write_text("     " + required + row + row2)
    b.write_text("     " + required + row + row2)
    assert m.compare(a, b, [1.0], 2, 0.2)["passed"]
    b.write_text(" bad\n" + row + row2)
    with pytest.raises(ValueError, match="header"):
        m.compare(a, b, [1.0], 2, 0.2)
