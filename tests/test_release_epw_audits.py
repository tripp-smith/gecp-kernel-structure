"""Synthetic controller receipts; not actual scientific evidence."""

import json
import subprocess
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import release_epw_audits as audit

SHA = "a" * 40


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def payload(root):
    write(root / "restart20/freeze.json", {"environment": {"revision": SHA}})
    write(root / "epw-anisotropic.json", {"source_revision": SHA})


@contextmanager
def snapshot(repository, revision):
    assert revision == SHA
    yield Path(repository)


def runner(cmd, **kwargs):
    assert (
        kwargs["timeout"] == 360
        and cmd[1] == "-I"
        and "PYTHONPATH" not in kwargs["env"]
    )
    if "--manifest" in cmd:
        out = Path(cmd[-1])
        runtime = out.parent / "runtime.json"
        write(
            runtime,
            {
                "platform": "linux",
                "machine": "x86_64",
                "python": "3.12.12",
                "numpy": "2.2.6",
                "scipy": "1.18.1",
                "mpmath": "1.3.0",
                "executable": "/pinned/python",
            },
        )
        metrics = dict.fromkeys(
            (
                "green_scaled",
                "green_global_scaled",
                "z",
                "gap_ev",
                "chi_ev",
                "imaginary_defect",
            ),
            0.0,
        )
        split = {
            "metrics": metrics,
            "cardinal_condition_max": 1.0,
            "recovery_condition_max": 1.0,
            "z_derivative_gain_max": 0.0,
        }
        cases = [
            {
                "repeat": i,
                "results": {
                    name: dict(split) for name in ("heldout", "training", "full")
                },
                "bytes": dict.fromkeys(
                    (
                        "gzip_ascii",
                        "model_only",
                        "model_plus_xor_bundle",
                        "raw_ascii",
                        "raw_binary_arrays",
                    ),
                    1,
                ),
                "training_design_condition_max": 1.0,
                "xor_correction_bytes": 8,
                "xor_correction_values": 1,
            }
            for i in range(3)
        ]
        write(
            out,
            {
                "schema": "kernelgecp.epw-restart-roundtrip-audit.v1",
                "status": "passed",
                "native_acceptance": "not_tested",
                "public_decoder": "consistency only",
                "scope": (
                    "finite 20 K fit-only holdout and all-field truth-assisted exact "
                    "reconstruction; native solver comparator required"
                ),
                "interpreter_identity_scope": (
                    "approved producer-recorded identity only; "
                    "executable bytes not provided"
                ),
                "cases": cases,
            },
        )
    else:
        out = Path(cmd[-1]) / "receipt.json"
        comparison = {
            "passed": True,
            "row_count": 630000,
            "state_count": 140,
            "frequency_count": 4500,
            "atol": 1e-5,
            "rtol": 1e-3,
            "budget_reference": "absolute raw value",
            "channels": {
                k: {
                    "failures": 0,
                    "max_budget_ratio": 0.0,
                    "max_absolute_error": 0.0,
                    "unit": "dimensionless" if k.startswith("Z_") else "eV",
                }
                for k in audit.CHANNELS
            },
        }
        write(
            out,
            {
                "schema": "kernelgecp.epw-native-payload-audit.v1",
                "status": "passed",
                "source_revision": SHA,
                "payload_unchanged": True,
                "payload_sha256": audit.inventory(
                    Path(cmd[cmd.index("--payload") + 1])
                ),
                "cases": ["raw", "reconstructed", "wrong-frequency"],
                "comparison_scope": (
                    "exact finite printed-number budgets, not continuous certification"
                ),
                "comparison": comparison,
            },
        )
    return SimpleNamespace(returncode=0, stdout="ok", stderr="")


def run(tmp_path, execute=runner):
    return audit.audit(
        tmp_path / "payload",
        SHA,
        "/runtime/python",
        tmp_path / "output",
        repository=tmp_path,
        runner=execute,
        snapshot=snapshot,
        identity_probe=lambda _: {"revision": "b" * 40, "source_sha256": {}},
    )


def test_success_preserves_inputs_and_is_exclusive(tmp_path):
    payload(tmp_path / "payload")
    before = audit.inventory(tmp_path / "payload")
    result = run(tmp_path)
    assert (
        result["status"] == "passed" and audit.inventory(tmp_path / "payload") == before
    )
    with pytest.raises(ValueError):
        run(tmp_path)


@pytest.mark.parametrize("mode", ["timeout", "missing", "bad_schema", "failed"])
def test_child_failure_and_partial_logs_are_preserved(tmp_path, mode):
    payload(tmp_path / "payload")

    def bad(cmd, **kwargs):
        if mode == "timeout":
            raise subprocess.TimeoutExpired(
                cmd, 360, output=b"partial out", stderr=b"partial err"
            )
        if mode == "missing":
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if mode == "failed":
            return SimpleNamespace(returncode=1, stdout="failed", stderr="detail")
        response = runner(cmd, **kwargs)
        p = Path(cmd[-1])
        doc = json.loads(p.read_text())
        doc["status"] = "failed"
        p.write_text(json.dumps(doc))
        return response

    result = run(tmp_path, bad)
    assert result["status"] == "failed" and (tmp_path / "output/receipt.json").exists()
    if mode == "timeout":
        assert "partial out" in (tmp_path / "output/archive.log").read_text()
        assert "partial err" in (tmp_path / "output/archive.log").read_text()


def test_rejects_contained_output_before_creation(tmp_path):
    payload(tmp_path / "payload")
    with pytest.raises(ValueError, match="unsafe output"):
        audit.audit(
            tmp_path / "payload",
            SHA,
            "/python",
            tmp_path / "payload/out",
            snapshot=snapshot,
        )
    assert not (tmp_path / "payload/out").exists()


@pytest.mark.parametrize(
    "mode", ["missing_native", "native_timeout", "runtime", "mutation"]
)
def test_native_runtime_and_final_immutability_fail_closed(tmp_path, mode):
    payload(tmp_path / "payload")

    def controlled(cmd, **kwargs):
        is_archive = "--manifest" in cmd
        if mode == "native_timeout" and not is_archive:
            raise subprocess.TimeoutExpired(
                cmd, 360, output=b"partial-native", stderr=b"detail"
            )
        if mode == "missing_native" and not is_archive:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        response = runner(cmd, **kwargs)
        if mode == "runtime" and is_archive:
            runtime = Path(cmd[-1]).parent / "runtime.json"
            data = json.loads(runtime.read_text())
            data["machine"] = "arm64"
            write(runtime, data)
        if mode == "mutation" and not is_archive:
            (tmp_path / "payload" / "tamper").write_text("changed")
        return response

    result = run(tmp_path, controlled)
    assert result["status"] == "failed"
    if mode == "native_timeout":
        assert "partial-native" in (tmp_path / "output/native.log").read_text()
    if mode == "mutation":
        assert result["payload_unchanged"] is False


@pytest.mark.parametrize(
    "bad",
    [
        "empty_split",
        "empty_bytes",
        "negative",
        "boolean",
        "nan",
        "missing_native_inventory",
        "native_boolean",
        "nonzero_passed_json",
    ],
)
def test_rejects_incomplete_or_untrusted_receipts(tmp_path, bad):
    payload(tmp_path / "payload")

    def invalid(cmd, **kwargs):
        response = runner(cmd, **kwargs)
        archive = "--manifest" in cmd
        if bad == "nonzero_passed_json" and archive:
            return SimpleNamespace(returncode=1, stdout="failed", stderr="failed")
        if archive and bad in (
            "empty_split",
            "empty_bytes",
            "negative",
            "boolean",
            "nan",
        ):
            path = Path(cmd[-1])
            data = json.loads(path.read_text())
            case = data["cases"][0]
            if bad == "empty_split":
                case["results"]["heldout"] = {}
            elif bad == "empty_bytes":
                case["bytes"] = {}
            else:
                case["results"]["heldout"]["metrics"]["green_scaled"] = {
                    "negative": -1.0,
                    "boolean": True,
                    "nan": float("nan"),
                }[bad]
            write(path, data)
        if not archive and bad in ("missing_native_inventory", "native_boolean"):
            path = Path(cmd[-1]) / "receipt.json"
            data = json.loads(path.read_text())
            if bad == "missing_native_inventory":
                data["payload_sha256"] = {}
            else:
                data["comparison"]["channels"]["Z_real"]["max_budget_ratio"] = True
            write(path, data)
        return response

    assert run(tmp_path, invalid)["status"] == "failed"


def test_bounded_runner_terminates_and_kills_process_group(monkeypatch):
    signals = []

    class FakeProcess:
        pid = 12345
        returncode = -9
        calls = 0

        def communicate(self, timeout=None):
            self.calls += 1
            if self.calls <= 2:
                raise subprocess.TimeoutExpired(
                    "cmd", timeout, output=b"partial", stderr=b"err"
                )
            return (b"saved", b"detail")

    process = FakeProcess()

    def spawn(*args, **kwargs):
        assert kwargs["start_new_session"] is True
        return process

    monkeypatch.setattr(audit.subprocess, "Popen", spawn)
    monkeypatch.setattr(audit.os, "killpg", lambda pid, sig: signals.append((pid, sig)))
    with pytest.raises(subprocess.TimeoutExpired) as result:
        audit.bounded_runner(
            ["cmd"],
            cwd="/",
            env={},
            text=True,
            capture_output=True,
            check=False,
            timeout=360,
        )
    assert signals == [(12345, audit.signal.SIGTERM), (12345, audit.signal.SIGKILL)]
    assert result.value.output == b"saved" and result.value.stderr == b"detail"


@pytest.mark.parametrize(
    "mode", ["wrong_policy", "missing_policy", "extra_field", "archive_scope"]
)
def test_strict_child_receipt_policy_and_schema(tmp_path, mode):
    payload(tmp_path / "payload")

    def altered(cmd, **kwargs):
        response = runner(cmd, **kwargs)
        archive = "--manifest" in cmd
        if archive and mode == "archive_scope":
            path = Path(cmd[-1])
            data = json.loads(path.read_text())
            data.pop("interpreter_identity_scope")
            write(path, data)
        elif not archive and mode != "archive_scope":
            path = Path(cmd[-1]) / "receipt.json"
            data = json.loads(path.read_text())
            if mode == "wrong_policy":
                data["comparison"]["budget_reference"] = "fitted value"
            elif mode == "missing_policy":
                data["comparison"].pop("budget_reference")
            else:
                data["comparison"]["new_policy"] = "unsupported"
            write(path, data)
        return response

    assert run(tmp_path, altered)["status"] == "failed"
