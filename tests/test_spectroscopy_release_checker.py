"""Orchestration rejection tests only; these do not create scientific evidence."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "gate",
    Path(__file__).resolve().parents[1] / "scripts/check_spectroscopy_release.py",
)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


class RejectionTests(unittest.TestCase):
    def test_json_rejects_ambiguity_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "e.json"
            for value in ('{"x": 1, "x": 2}', '{"x": NaN}', '{"x": 1e999}', "{"):
                p.write_text(value)
                with self.subTest(value=value), self.assertRaises(ValueError):
                    gate.load(p)

    def test_missing_future_evidence_rejected(self):
        with (
            tempfile.TemporaryDirectory() as d,
            patch.object(gate, "source_check"),
            self.assertRaises(ValueError),
        ):
            gate.check(Path(d), gate.FROZEN)

    def test_synthetic_requires_24_independent_and_72_children(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            manifest, freeze, audit = [
                root / x for x in ("manifest.json", "freeze.json", "audit.json")
            ]
            freeze.write_text("{}")
            manifest.write_text(
                json.dumps({"children": {str(i): {} for i in range(72)}})
            )
            base = {
                "status": "passed",
                "kind": "campaign",
                "cases": 72,
                "independently_audited_cases": 24,
                "manifest_sha256": gate.sha(manifest),
                "freeze_sha256": gate.sha(freeze),
            }
            for change in (
                {"status": "failed"},
                {"cases": 24},
                {"independently_audited_cases": 72},
                {"manifest_sha256": "0" * 64},
            ):
                audit.write_text(json.dumps({**base, **change}))
                with self.subTest(change=change), self.assertRaises(ValueError):
                    gate.synthetic_summary(audit, manifest, freeze, "campaign", 72)
            audit.write_text(json.dumps(base))
            manifest.write_text(
                json.dumps({"children": {str(i): {} for i in range(71)}})
            )
            base["manifest_sha256"] = gate.sha(manifest)
            audit.write_text(json.dumps(base))
            with self.assertRaises(ValueError):
                gate.synthetic_summary(audit, manifest, freeze, "campaign", 72)

    def test_full_candidate_and_clean_frozen_checkout_required(self):
        with (
            patch.object(gate, "git", side_effect=[(gate.FROZEN + "\n").encode(), b""]),
            self.assertRaises(ValueError),
        ):
            gate.source_check("short")
        with (
            patch.object(
                gate, "git", side_effect=[(gate.FROZEN + "\n").encode(), b" M file"]
            ),
            self.assertRaises(ValueError),
        ):
            gate.source_check(gate.FROZEN)


if __name__ == "__main__":
    unittest.main()
