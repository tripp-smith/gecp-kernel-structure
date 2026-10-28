"""Boundary stopping must preserve a resumable, complete chain pair."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from scripts.providers.cthyb_checkpoint import complete_pairs, segment


def test_segment_stops_only_after_completed_chain(tmp_path: Path) -> None:
    chain = tmp_path / "dense-final"
    chain.mkdir()
    command = [
        sys.executable,
        "-c",
        "import json,pathlib,sys,time; p=pathlib.Path(sys.argv[1]); "
        "time.sleep(.1); (p/'replica-00.npz').write_bytes(b'x'); "
        "(p/'replica-00.json').write_text('{}'); "
        "print(json.dumps({'measurement':'dense-final/replica-00'}),flush=True); "
        "time.sleep(20)",
        str(chain),
    ]
    assert segment(command, tmp_path, 0.02, 5) == "checkpointed"
    assert complete_pairs(tmp_path)
    (chain / "replica-00.json").unlink()
    assert not complete_pairs(tmp_path)


def test_unexpected_exit_cannot_look_like_checkpoint(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="unexpectedly"):
        segment([sys.executable, "-c", "raise SystemExit(3)"], tmp_path, 1, 1)
