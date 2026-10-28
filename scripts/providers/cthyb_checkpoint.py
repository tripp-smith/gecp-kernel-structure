"""Bound a CT-HYB validation segment at a completed measurement boundary.

A zero exit means either a complete campaign or a clean resumable checkpoint;
the final job must still run the scientific auditor and check its result.
"""

from __future__ import annotations

import argparse
import json
import os
import selectors
import subprocess
import sys
import time
from pathlib import Path


def complete_pairs(folder: Path) -> bool:
    for pattern in ("replica-*", "iteration-*"):
        stems: dict[Path, set[str]] = {}
        for item in folder.glob(f"*/{pattern}"):
            if item.suffix in (".json", ".npz"):
                stems.setdefault(item.with_suffix(""), set()).add(item.suffix)
        if any(pair != {".json", ".npz"} for pair in stems.values()):
            return False
    return True


def segment(command: list[str], folder: Path, seconds: float, grace: float) -> str:
    """Stop only after a complete chain is durably logged; never remove data."""
    if seconds <= 0 or grace <= 0:
        raise ValueError("segment and grace durations must be positive")
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    assert process.stdout is not None
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    start = time.monotonic()
    pending = b""
    try:
        while True:
            ready = selector.select(timeout=1)
            if ready:
                block = os.read(process.stdout.fileno(), 65536)
                if block:
                    pending += block
                    while b"\n" in pending:
                        line, pending = pending.split(b"\n", 1)
                        print(line.decode("utf-8", errors="replace"), flush=True)
                        try:
                            message = json.loads(line)
                        except (ValueError, UnicodeDecodeError):
                            message = None
                        if (
                            time.monotonic() - start >= seconds
                            and isinstance(message, dict)
                            and isinstance(message.get("measurement"), str)
                        ):
                            process.terminate()
                            process.wait(timeout=30)
                            if not complete_pairs(folder):
                                raise RuntimeError(
                                    "incomplete checkpoint: manual inspection required"
                                )
                            return "checkpointed"
                else:
                    selector.unregister(process.stdout)
            if process.poll() is not None:
                if pending:
                    print(pending.decode("utf-8", errors="replace"), flush=True)
                if not complete_pairs(folder):
                    raise RuntimeError(
                        "incomplete checkpoint: manual inspection required"
                    )
                if process.returncode in (0, 2) and (folder / "result.json").is_file():
                    return "completed"
                raise RuntimeError(
                    f"campaign exited unexpectedly: {process.returncode}"
                )
            if time.monotonic() - start > seconds + grace:
                process.kill()
                process.wait(timeout=30)
                raise RuntimeError(
                    "no completed measurement within checkpoint grace period"
                )
    finally:
        selector.close()
        if process.poll() is None:
            process.kill()
            process.wait(timeout=30)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--minutes", type=float, default=240)
    parser.add_argument("--grace-minutes", type=float, default=30)
    args = parser.parse_args()
    if (args.folder / "result.json").exists():
        if not complete_pairs(args.folder):
            raise RuntimeError("completed campaign has incomplete chain files")
        print(json.dumps({"segment": "already_completed"}))
        return 0
    status = segment(
        [
            sys.executable,
            "-m",
            "experiments.cthyb_dmft",
            "final",
            "--config",
            "experiments/cthyb_fresh_odd.toml",
            "--output",
            str(args.folder),
        ],
        args.folder,
        args.minutes * 60,
        args.grace_minutes * 60,
    )
    print(json.dumps({"segment": status}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
