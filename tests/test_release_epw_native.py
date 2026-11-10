import gzip
import hashlib
import time

import pytest

from scripts import release_epw_native as draft

HEADER = (
    b"# Re[znorm(w)] Im[znorm(w)] Re[delta(w)] Im[delta(w)] Re[shift(w)] Im[shift(w)]\n"
)


def data():
    return HEADER + b"".join(
        (
            "".join(f"{x:20.10E}" for x in (freq, -1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
            + "\n"
        ).encode()
        for freq in (0.1, 0.2)
    )


def inputs(tmp_path):
    content = data()
    files = []
    for name in ("a.gz", "b.gz"):
        path = tmp_path / name
        with gzip.open(path, "wb") as f:
            f.write(content)
        files.append(path)
    meta = {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
    return files, [meta, meta]


def test_small_positive(tmp_path):
    paths, meta = inputs(tmp_path)
    result = draft.table_streams(*paths, [-1.0], 2, 0.1, meta, time.monotonic() + 5)
    assert result["passed"] and result["row_count"] == 2


def test_tamper_and_timeout(tmp_path):
    paths, meta = inputs(tmp_path)
    meta[0] = {**meta[0], "sha256": "0" * 64}
    with pytest.raises(ValueError, match="hash/size"):
        draft.table_streams(*paths, [-1.0], 2, 0.1, meta, time.monotonic() + 5)
    with pytest.raises(ValueError, match="timeout"):
        draft.table_streams(*paths, [-1.0], 2, 0.1, meta, time.monotonic() - 1)


def test_expansion_limit(tmp_path, monkeypatch):
    paths, meta = inputs(tmp_path)
    monkeypatch.setattr(draft, "MAX_FILE", len(HEADER) + 10)
    with pytest.raises(ValueError, match="expansion limit"):
        draft.table_streams(*paths, [-1.0], 2, 0.1, meta, time.monotonic() + 5)


def test_restart_axes_validation(tmp_path):
    path = tmp_path / "native"

    def row(freq, energy):
        return (
            "".join(f"{x:20.10E}" for x in (freq, energy, 1.0, 0.1, 0.2)) + "\n"
        ).encode()

    path.write_bytes(b"# shift(w)\n" + row(0.1, -1) + row(0.2, -1))
    assert draft.restart_axes(path, 2, 1) == ([0.1, 0.2], [-1.0])
    path.write_bytes(b"# shift(w)\n" + row(0.2, -1) + row(0.1, -1))
    with pytest.raises(ValueError, match="frequency ordering"):
        draft.restart_axes(path, 2, 1)


def test_modified_input_exact_transform():
    original = (
        "&inputepw\n temps = 17 20 30 ! K\n iverbosity = 0\n"
        " limag = .true.\n wscut = 1.0 ! eV\n/\n"
    )
    out = draft.expected_input(original)
    assert "temps = 20 ! K" in out and "imag_read = .true." in out
    with pytest.raises(ValueError):
        draft.expected_input(original.replace("17 20 30", "17 30"))


def test_paths_and_inventory_reject_symlinks_hardlinks_and_containment(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    (root / "data").write_text("x")
    (tmp_path / "alias").symlink_to(root, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        draft.checked_path(tmp_path / "alias" / "x")
    with pytest.raises(ValueError, match="unsafe output"):
        draft.audit_native_comparator(root, "a" * 40, root / "x")
    assert not (root / "x").exists()
    (root / "other").hardlink_to(root / "data")
    with pytest.raises(ValueError, match="unsafe payload"):
        draft.inventory(root)


def test_whole_deadline_preserves_failed_receipt(tmp_path, monkeypatch):
    root = tmp_path / "input"
    root.mkdir()
    (root / "restart20-native-comparison.json").write_text("{}")
    monkeypatch.setattr(draft, "LIMIT", 0.03)

    def stall(_):
        time.sleep(0.1)

    monkeypatch.setattr(draft, "read_json", stall)
    out = tmp_path / "out"
    result = draft.audit_native_comparator(root, "a" * 40, out)
    assert result["status"] == "failed" and "TimeoutError" in result["error"]
    assert (out / "audit.log").exists() and (out / "receipt.json").exists()


def test_bad_gzip_is_reported(tmp_path):
    paths, meta = inputs(tmp_path)
    paths[0].write_bytes(b"not gzip")
    with pytest.raises(OSError):
        draft.table_streams(*paths, [-1.0], 2, 0.1, meta, time.monotonic() + 5)


def test_complete_saved_payload_immutable_on_success_and_failure(tmp_path):
    import os
    from pathlib import Path

    root_name = os.environ.get("EPW_DRAFT_FIXTURE")
    if not root_name:
        pytest.skip("no saved authenticated payload fixture")
    root = Path(root_name)
    before = draft.inventory(root)
    good = draft.audit_native_comparator(
        root, "72eae3353cab7a440e4154f9d906f204ffe099b8", tmp_path / "good"
    )
    assert good["status"] == "passed" and draft.inventory(root) == before
    bad = draft.audit_native_comparator(root, "a" * 40, tmp_path / "bad")
    assert bad["status"] == "failed" and draft.inventory(root) == before
    assert (tmp_path / "good" / "audit.log").read_text().endswith("passed \n")
    assert "failed" in (tmp_path / "bad" / "audit.log").read_text()


def test_coherent_axis_metadata_tampering_still_fails(tmp_path):
    first = tmp_path / "raw"
    second = tmp_path / "restored"

    def make(path, delta):
        def row(w, e):
            return (
                "".join(f"{x:20.10E}" for x in (w, e, 1.0, 0.1, 0.2)) + "\n"
            ).encode()

        path.write_bytes(
            b"# shift(w)\n" + row(0.1 + delta, -1.0) + row(0.2 + delta, -1.0)
        )

    make(first, 0.0)
    make(second, 0.01)
    raw = draft.restart_axes(first, 2, 1)
    restored = draft.restart_axes(second, 2, 1)
    # Each parsed file agrees with its own (coherently changed) metadata.
    assert restored == ([0.11, 0.21], [-1.0])
    with pytest.raises(ValueError, match="axes differ"):
        draft.compare_native_axes(raw, restored)


def test_active_timer_rejected_without_output(tmp_path):
    import signal

    root = tmp_path / "input"
    root.mkdir()
    previous = signal.getsignal(signal.SIGALRM)
    signal.setitimer(signal.ITIMER_REAL, 10)
    try:
        with pytest.raises(ValueError, match="active caller timer"):
            draft.audit_native_comparator(root, "a" * 40, tmp_path / "out")
        assert (
            signal.getitimer(signal.ITIMER_REAL)[0] > 0
            and not (tmp_path / "out").exists()
        )
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def test_hardlink_rejected_before_output_creation(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    (root / "a").write_text("x")
    (root / "b").hardlink_to(root / "a")
    with pytest.raises(ValueError, match="unsafe payload"):
        draft.audit_native_comparator(root, "a" * 40, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_aggregate_inventory_budget(tmp_path, monkeypatch):
    root = tmp_path / "input"
    root.mkdir()
    (root / "a").write_bytes(b"123")
    (root / "b").write_bytes(b"456")
    monkeypatch.setattr(draft, "MAX_INVENTORY", 5)
    with pytest.raises(ValueError, match="inventory budget"):
        draft.inventory(root)


def test_final_inventory_after_writes_cannot_report_pass(tmp_path, monkeypatch):
    import os
    from pathlib import Path

    root_name = os.environ.get("EPW_DRAFT_FIXTURE")
    if not root_name:
        pytest.skip("no saved authenticated payload fixture")
    root = Path(root_name)
    original = draft.inventory
    calls = []
    output = tmp_path / "result"

    def altered(root_path):
        data = original(root_path)
        calls.append(1)
        if len(calls) == 2:
            assert (output / "receipt.json").exists() and (
                output / "audit.log"
            ).exists()
            data["simulated-changed-file"] = "0" * 64
        return data

    monkeypatch.setattr(draft, "inventory", altered)
    result = draft.audit_native_comparator(
        root, "72eae3353cab7a440e4154f9d906f204ffe099b8", output
    )
    assert result["status"] == "failed" and result["payload_unchanged"] is False
    assert (
        __import__("json").loads((output / "receipt.json").read_text())["status"]
        == "failed"
    )


def test_exact_printed_decimal_catches_binary64_boundary():
    raw = b"6.1933334092759E-2"
    fitted = b"6.200526742685176E-2"
    assert len(raw) <= 20 and len(fitted) <= 20
    assert abs(float(fitted) - float(raw)) <= 1e-5 + 1e-3 * abs(float(raw))
    diff, ratio = draft.exact_field_budget(raw, fitted)
    assert diff > 0 and ratio > 1


def test_equal_numeric_different_printed_formats():
    diff, ratio = draft.exact_field_budget(b"1.0000000E+00", b"1.0E+00")
    assert diff == 0 and ratio == 0
