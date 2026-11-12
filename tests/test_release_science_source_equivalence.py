import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).parents[1] / "scripts" / "release_science_source_equivalence.py"
spec = importlib.util.spec_from_file_location("guard", MODULE)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class ScienceSourceEquivalenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        self.git("init", "-q")
        self.git("config", "user.name", "Test")
        self.git("config", "user.email", "test@example.org")
        for path in (
            "kernelgecp/core.py",
            "experiments/run.py",
            "scripts/providers/run.py",
            "tests/data/example.json",
            "uv.lock",
            ".github/workflows/science.yml",
        ):
            self.write(path, "original\n")
        self.science = self.save()

    def git(self, *args):
        return (
            subprocess.check_output(
                ["git", "-C", str(self.repo), *args], stderr=subprocess.PIPE
            )
            .decode()
            .strip()
        )

    def write(self, path, content):
        dest = self.repo / path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(content)

    def save(self):
        self.git("add", "-A")
        self.git("commit", "-qm", "fixture")
        return self.git("rev-parse", "HEAD")

    def compare(self, producer):
        return guard.compare(self.repo, self.science, producer)

    def test_same_and_allowlisted_new_helpers(self):
        self.assertEqual(self.compare(self.science)["protected_paths"], 6)
        self.write("scripts/release_spectroscopy_v4.py", "helper\n")
        self.write("README.md", "documentation\n")
        result = self.compare(self.save())
        self.assertEqual(
            result["allowed_additions"], ["scripts/release_spectroscopy_v4.py"]
        )

    def test_modification_or_deletion_of_each_protected_kind_rejected(self):
        for path in (
            "kernelgecp/core.py",
            "experiments/run.py",
            "scripts/providers/run.py",
            "tests/data/example.json",
            "uv.lock",
            ".github/workflows/science.yml",
        ):
            with self.subTest(path=path):
                self.git("reset", "--hard", self.science)
                self.write(path, "changed\n")
                with self.assertRaises(ValueError):
                    self.compare(self.save())
                self.git("reset", "--hard", self.science)
                (self.repo / path).unlink()
                with self.assertRaises(ValueError):
                    self.compare(self.save())

    def test_mode_change_and_unknown_additions_rejected(self):
        (self.repo / "kernelgecp/core.py").chmod(0o755)
        with self.assertRaises(ValueError):
            self.compare(self.save())
        self.git("reset", "--hard", self.science)
        self.write("scripts/release_unsafe.py", "new\n")
        with self.assertRaises(ValueError):
            self.compare(self.save())

    def test_existing_release_path_may_change_but_not_be_removed(self):
        self.write("scripts/release_spectroscopy_v4.py", "first\n")
        origin = self.save()
        self.write("scripts/release_spectroscopy_v4.py", "changed\n")
        target = self.save()
        self.assertEqual(
            guard.compare(self.repo, origin, target)["changed_release_only"],
            ["scripts/release_spectroscopy_v4.py"],
        )
        (self.repo / "scripts/release_spectroscopy_v4.py").unlink()
        with self.assertRaises(ValueError):
            guard.compare(self.repo, origin, self.save())

    def test_missing_invalid_and_nonancestor_commits_rejected(self):
        for identity in ("a" * 40, "HEAD", self.science.upper()):
            with (
                self.subTest(identity=identity),
                self.assertRaises((ValueError, subprocess.CalledProcessError)),
            ):
                self.compare(identity)
        self.write("kernelgecp/core.py", "later\n")
        later = self.save()
        with self.assertRaises(subprocess.CalledProcessError):
            guard.compare(self.repo, later, self.science)

    def test_symlink_rejected(self):
        (self.repo / "kernelgecp/core.py").unlink()
        (self.repo / "kernelgecp/core.py").symlink_to("../uv.lock")
        with self.assertRaises(ValueError):
            self.compare(self.save())

    def test_special_mode_rejected_for_excluded_and_release_paths(self):
        for path in ("README.md", "scripts/release_spectroscopy_v4.py"):
            with self.subTest(path=path):
                self.git("reset", "--hard", self.science)
                dest = self.repo / path
                dest.symlink_to("uv.lock")
                with self.assertRaises(ValueError):
                    self.compare(self.save())

    def test_exact_release_documents_can_change_but_neighbor_is_protected(self):
        for path in (
            "PACKAGE_README.md",
            "INSTALLATION_COMPATIBILITY.md",
            "MIGRATION.md",
            "PROVIDER_NOTICES.md",
            "RELEASE_NOTES.md",
            "CITATION.cff",
            "METHODS.md",
            "INTERACTING_DMFT.md",
            "CTHYB_DMFT.md",
            "MOLECULAR_GREEN.md",
            "EPW_INTEROPERABILITY.md",
            "COVARIANCE.md",
        ):
            self.write(path, "initial\n")
        origin = self.save()
        for path in (
            "PACKAGE_README.md",
            "INSTALLATION_COMPATIBILITY.md",
            "MIGRATION.md",
            "PROVIDER_NOTICES.md",
            "RELEASE_NOTES.md",
            "CITATION.cff",
            "METHODS.md",
            "INTERACTING_DMFT.md",
            "CTHYB_DMFT.md",
            "MOLECULAR_GREEN.md",
            "EPW_INTEROPERABILITY.md",
            "COVARIANCE.md",
        ):
            self.write(path, "updated\n")
        guard.compare(self.repo, origin, self.save())
        self.write("OTHER_RELEASE_NOTES.md", "unknown\n")
        with self.assertRaises(ValueError):
            guard.compare(self.repo, origin, self.save())

    def test_cthyb_scope_only_and_exact_helper_additions(self):
        path = ".github/workflows/cthyb-science.yml"
        self.write(
            path,
            (
                "name: science\njobs:\n  scope:\n    runs-on: old\n"
                "  calibration:\n    runs-on: fixed\n"
            ),
        )
        origin = self.save()
        self.write(
            path,
            (
                "name: science\njobs:\n  scope:\n    runs-on: new\n"
                "  calibration:\n    runs-on: fixed\n"
            ),
        )
        self.write("scripts/release_cthyb_campaign_scope.py", "helper\n")
        self.write("tests/test_release_cthyb_campaign_scope.py", "tests\n")
        result = guard.compare(self.repo, origin, self.save())
        self.assertEqual(
            result["allowed_additions"],
            [
                "scripts/release_cthyb_campaign_scope.py",
                "tests/test_release_cthyb_campaign_scope.py",
            ],
        )

    def test_cthyb_changes_outside_scope_or_mode_rejected(self):
        path = ".github/workflows/cthyb-science.yml"
        original = (
            "name: science\njobs:\n  scope:\n    runs-on: old\n"
            "  calibration:\n    runs-on: fixed\n"
        )
        self.write(path, original)
        origin = self.save()
        for changed in (
            original.replace("name: science", "name: changed"),
            original.replace("runs-on: fixed", "runs-on: changed"),
            original.replace(
                "  calibration:", "  other:\n    runs-on: changed\n  calibration:"
            ),
        ):
            with self.subTest(changed=changed):
                self.git("reset", "--hard", origin)
                self.write(path, changed)
                with self.assertRaises(ValueError):
                    guard.compare(self.repo, origin, self.save())
        self.git("reset", "--hard", origin)
        (self.repo / path).chmod(0o755)
        with self.assertRaises(ValueError):
            guard.compare(self.repo, origin, self.save())
        self.git("reset", "--hard", origin)
        (self.repo / path).unlink()
        with self.assertRaises(ValueError):
            guard.compare(self.repo, origin, self.save())

    def test_cthyb_malformed_scope_delimiters_rejected(self):
        path = ".github/workflows/cthyb-science.yml"
        original = (
            "jobs:\n  scope:\n    runs-on: old\n  calibration:\n    runs-on: fixed\n"
        )
        self.write(path, original)
        origin = self.save()
        for changed in (
            original.replace("  scope:", "  missing:"),
            original.replace("  calibration:", "  missing:"),
            original.replace(
                "  calibration:", "  scope:\n    runs-on: duplicate\n  calibration:"
            ),
            original.replace("  calibration:", "  'calibration':"),
            original.replace("jobs:", "jobs: # malformed"),
            original.replace("  calibration:", "env:\n  X: unsafe\n  calibration:"),
        ):
            with self.subTest(changed=changed):
                self.git("reset", "--hard", origin)
                self.write(path, changed)
                with self.assertRaises(ValueError):
                    guard.compare(self.repo, origin, self.save())

    def test_license_is_protected(self):
        self.write("LICENSE", "license\n")
        with self.assertRaises(ValueError):
            self.compare(self.save())


if __name__ == "__main__":
    unittest.main()
