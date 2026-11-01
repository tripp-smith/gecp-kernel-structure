"""Focused checks of replay infrastructure, not a scientific acceptance fixture."""

import importlib.util
import io
import subprocess
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


replay = module("replay", ROOT / "experiments/replay_spectroscopy_acceptance.py")
gate = module("gate", ROOT / "scripts/check_spectroscopy_release.py")


class ReplayTests(unittest.TestCase):
    def test_unsafe_paths_and_links_rejected(self):
        for name, kind in [
            ("../escape", "file"),
            ("/absolute", "file"),
            ("normal", "link"),
            ("a/../b", "file"),
        ]:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                archive = Path(directory) / "bad.tar.gz"
                with tarfile.open(archive, "w:gz") as z:
                    item = tarfile.TarInfo(name)
                    if kind == "link":
                        item.type = tarfile.SYMTYPE
                        item.linkname = "elsewhere"
                        z.addfile(item)
                    else:
                        item.size = 1
                        z.addfile(item, io.BytesIO(b"x"))
                with self.assertRaises(ValueError):
                    replay.extract_tar(
                        archive, Path(directory) / "out", count=1, total=1
                    )

    def test_duplicate_and_incorrect_count_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "bad.tar.gz"
            with tarfile.open(archive, "w:gz") as z:
                for _ in range(2):
                    entry = tarfile.TarInfo("x")
                    entry.size = 1
                    z.addfile(entry, io.BytesIO(b"x"))
            with self.assertRaises(ValueError):
                replay.extract_tar(archive, Path(directory) / "out", count=2, total=2)

    def test_corrupt_archive_fails_before_clone_or_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / "experiments/data"
            folder.mkdir(parents=True)
            archive = folder / "spectroscopy-linux-feaac66-passed.tar.gz"
            archive.write_bytes(b"corrupt")
            output = root / "result.json"
            with (
                patch.object(replay, "ROOT", root),
                patch.object(replay.sys, "platform", "linux"),
                patch.object(replay.sys, "argv", ["replay", "--output", str(output)]),
                patch.object(replay.subprocess, "run") as invoke,
                patch.object(replay.sys, "stderr", io.StringIO()),
            ):
                self.assertEqual(replay.main(), 1)
                invoke.assert_not_called()
            self.assertIn("archive hash mismatch", output.read_text())

    def test_native_provider_versions_do_not_require_distribution_metadata(self):
        from types import SimpleNamespace

        native = {"triqs": "4.0.2", "triqs_cthyb": "4.0.1", "triqs_maxent": "4.0.1"}

        def distribution(name):
            if name in native:
                raise AssertionError("Conda native namespaces have no wheel metadata")
            return "declared-python-distribution"

        with (
            patch.object(
                replay.importlib.metadata, "version", side_effect=distribution
            ),
            patch.object(
                replay.importlib,
                "import_module",
                side_effect=lambda name: SimpleNamespace(
                    version=native[name.removesuffix(".version")]
                ),
            ),
        ):
            actual = replay.runtime_versions()
        self.assertEqual({name: actual[name] for name in native}, native)
        self.assertEqual(len(actual), 9)

    def test_native_source_change_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            checkout = Path(directory) / "frozen"
            subprocess.run(
                ["git", "clone", "--shared", "--no-checkout", str(ROOT), str(checkout)],
                check=True,
                capture_output=True,
            )

            def git_command(*args):
                return (
                    subprocess.check_output(
                        [
                            "git",
                            "-c",
                            "core.hooksPath=/dev/null",
                            "-C",
                            str(checkout),
                            *args,
                        ],
                        stderr=subprocess.PIPE,
                    )
                    .decode()
                    .strip()
                )

            git_command("checkout", "--detach", gate.FROZEN)
            native = checkout / "native/src/lib.rs"
            native.write_text(native.read_text() + "\n// Regression mutation\n")
            git_command("add", "native/src/lib.rs")
            git_command(
                "-c",
                "user.name=Fixture",
                "-c",
                "user.email=fixture@example.invalid",
                "-c",
                "commit.gpgsign=false",
                "commit",
                "-qm",
                "source mutation",
            )
            changed = git_command("rev-parse", "HEAD")
            git_command("checkout", "--detach", gate.FROZEN)
            with (
                patch.object(gate, "ROOT", checkout),
                self.assertRaisesRegex(ValueError, "scientific source or lock changed"),
            ):
                gate.source_check(changed)


if __name__ == "__main__":
    unittest.main()
