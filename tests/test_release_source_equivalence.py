import importlib.util
import pathlib
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "eq", ROOT / "scripts" / "release_source_equivalence.py"
)
eq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eq)
SOURCE = str(ROOT)


def run(repo, *args):
    return (
        subprocess.check_output(["git", "-C", str(repo), *args], stderr=subprocess.PIPE)
        .decode()
        .strip()
    )


class EquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="gecp-equiv-fixture-")
        cls.repo = pathlib.Path(cls.tmp.name)
        run(cls.repo, "init", "-q")
        run(cls.repo, "config", "user.name", "Fixture")
        run(cls.repo, "config", "user.email", "fixture@example.invalid")
        paths = sorted(
            eq.VERSION_FILES
            | {
                "kernelgecp/kernels.py",
                "native/src/lib.rs",
                "rust-toolchain.toml",
                "environments/epw/sources.lock.json",
                "experiments/cthyb_dmft.py",
            }
        )
        for p in paths:
            dest = cls.repo / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(eq.git(SOURCE, "show", f"{eq.BASELINE}:{p}"))
        run(cls.repo, "add", ".")
        run(cls.repo, "commit", "-qm", "baseline")
        cls.base = run(cls.repo, "rev-parse", "HEAD")
        for p in eq.VERSION_FILES:
            f = cls.repo / p
            data = f.read_bytes()
            # Replace only the selected field by using the checker's strict locator.
            f.write_bytes(eq.normalized(p, data, eq.OLD).replace(b"<VERSION>", eq.NEW))
        run(cls.repo, "add", ".")
        run(cls.repo, "commit", "-qm", "controlled release")
        cls.good = run(cls.repo, "rev-parse", "HEAD")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        run(self.repo, "reset", "--hard", "-q", self.good)
        run(self.repo, "clean", "-fdq")

    def compare(self, candidate, baseline=None):
        with patch.object(eq, "BASELINE", baseline or self.base):
            return eq.compare(str(self.repo), candidate)

    def rejected(self, path, transform):
        f = self.repo / path
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(transform(f.read_bytes() if f.exists() else b""))
        run(self.repo, "add", ".")
        run(self.repo, "commit", "-qm", "mutation")
        with self.assertRaises((ValueError, subprocess.CalledProcessError)):
            self.compare(run(self.repo, "rev-parse", "HEAD"))

    def test_accept_baseline_identity(self):
        result = self.compare(self.base)
        self.assertEqual(result["version_transition"], "1.1.0 -> 1.1.0")

    def test_accept_controlled_bump(self):
        r = self.compare(self.good)
        self.assertEqual(r["candidate_commit"], self.good)

    def test_exact_fixture_diff(self):
        targets = {
            "kernelgecp/__init__.py": b'__version__ = "1.1.0"',
            "pyproject.toml": b'name = "kernelgecp"\nversion = "1.1.0"',
            "uv.lock": b'name = "kernelgecp"\nversion = "1.1.0"',
            "native/Cargo.toml": b'name = "kernelgecp-native"\nversion = "1.1.0"',
            "native/Cargo.lock": b'name = "kernelgecp-native"\nversion = "1.1.0"',
        }
        for p, old in targets.items():
            before = eq.git(str(self.repo), "show", f"{self.base}:{p}")
            after = eq.git(str(self.repo), "show", f"{self.good}:{p}")
            self.assertEqual(before.count(old), 1)
            self.assertEqual(
                after, before.replace(old, old.replace(b"1.1.0", b"1.2.0"), 1)
            )

    def test_dependency_change(self):
        self.rejected(
            "uv.lock",
            lambda s: s.replace(
                b'name = "mypy-extensions"\nversion = "1.1.0"',
                b'name = "mypy-extensions"\nversion = "1.1.1"',
            ),
        )

    def test_feature_flag(self):
        self.rejected(
            "native/Cargo.toml",
            lambda s: s.replace(b"[features]", b"[features]\nnew_feature = []"),
        )

    def test_numerical_function(self):
        self.rejected(
            "kernelgecp/kernels.py",
            lambda s: s + b"\ndef injected_numeric_function(): return 0\n",
        )

    def test_ambiguous_stanza(self):
        self.rejected(
            "uv.lock",
            lambda s: s + b'\n[[package]]\nname = "kernelgecp"\nversion = "1.2.0"\n',
        )

    def test_wrong_version(self):
        self.rejected(
            "pyproject.toml",
            lambda s: s.replace(b'version = "1.2.0"', b'version = "1.2.1"', 1),
        )

    def test_extra_import(self):
        self.rejected(
            "kernelgecp/__init__.py",
            lambda s: s.replace(b"__version__ =", b"import decimal\n__version__ ="),
        )

    def test_new_module(self):
        self.rejected("kernelgecp/new_numeric.py", lambda s: b"def f(): return 1\n")

    def test_native_addition(self):
        self.rejected("native/src/new.rs", lambda s: b"pub fn f() {}\n")

    def test_release_additions_allowed(self):
        for p in eq.ALLOWED_ADDITIONS:
            f = self.repo / p
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text("# audit tool\n")
        run(self.repo, "add", ".")
        run(self.repo, "commit", "-qm", "audit additions")
        self.assertEqual(
            set(self.compare(run(self.repo, "rev-parse", "HEAD"))["allowed_additions"]),
            eq.ALLOWED_ADDITIONS,
        )

    def test_new_numeric_experiment(self):
        self.rejected("experiments/new_numeric.py", lambda s: b"def f(): return 1\n")

    def test_new_script(self):
        self.rejected("scripts/new_numeric.py", lambda s: b"def f(): return 1\n")

    def test_new_lock(self):
        self.rejected("other.lock", lambda s: b"dependency = 1\n")

    def test_changed_lock(self):
        self.rejected("environments/epw/sources.lock.json", lambda s: s + b" ")

    def test_changed_build_source(self):
        self.rejected("native/src/lib.rs", lambda s: s + b"// changed\n")

    def test_partial_bump(self):
        self.rejected(
            "pyproject.toml",
            lambda s: s.replace(b'version = "1.2.0"', b'version = "1.1.0"', 1),
        )

    def test_actual_baseline_inventory(self):
        tree = eq.tree(SOURCE, eq.BASELINE)
        for p in (
            "experiments/audit_cthyb_fixed_bath.py",
            "experiments/replay_cthyb_acceptance.py",
            "experiments/cthyb_loop_state.py",
            "scripts/providers/cthyb_checkpoint.py",
            "environments/triqs.yml",
            "native/build.rs",
            "tests/data/dmft_native_ipt_v1.npz",
        ):
            self.assertIn(p, tree)
        for p in eq.VERSION_FILES:
            eq.normalized(p, eq.blob(SOURCE, tree[p][1]), eq.OLD)

    def test_bad_identity_and_ancestry(self):
        with self.assertRaises(ValueError):
            self.compare(self.good[:7])
        with self.assertRaises(subprocess.CalledProcessError):
            self.compare(self.base, baseline=self.good)


if __name__ == "__main__":
    unittest.main(verbosity=2)
