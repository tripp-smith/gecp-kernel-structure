import copy
import json
import subprocess
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import release_cthyb_audits as gate

P = "a" * 40


def example():
    ref = json.loads((gate.ROOT / gate.INPUTS[6]).read_text())
    ref["audit_revision"] = P
    return {
        "generation_revision": gate.GENERATION,
        "replay_revision": P,
        "archive_sha256": gate.ZIP_SHA,
        "binding_sha256": gate.BINDING_SHA,
        "recomputed_audit": ref,
    }


def runtime():
    return {
        "platform": "linux",
        "machine": "x86_64",
        "python": "3.12.16",
        "executable": "/env/bin/python",
        "prefix": "/env",
        "versions": {
            "triqs": "4.0.2",
            "cthyb": "4.0.1",
            "numpy": "2.2.6",
            "scipy": "1.18.1",
            "pydlr": "1.0.1",
            "sparse-ir": "2.1.6",
            "pylibsparseir": "0.12.0",
        },
        "threads": dict.fromkeys(gate.THREADS, "1"),
        "modules": {
            x: f"/env/lib/{x}.py"
            for x in (
                "triqs",
                "triqs_cthyb",
                "pydlr",
                "sparse_ir",
                "numpy",
                "scipy",
                "pylibsparseir",
            )
        },
    }


def test_real_receipt_shape_and_numeric_tolerance():
    value = example()
    ref = copy.deepcopy(value["recomputed_audit"])
    value["recomputed_audit"]["max_standard_error"] *= 1.00001
    gate.validate_replay(value, P, ref)


@pytest.mark.parametrize(
    "change",
    [
        lambda x: x.update(recomputed_audit={"status": "passed"}),
        lambda x: x.update(archive_sha256="0" * 64),
        lambda x: x.update(replay_revision="b" * 40),
        lambda x: x["recomputed_audit"].update(schema="fake"),
        lambda x: x["recomputed_audit"].update(independent_chain_count=207),
        lambda x: x["recomputed_audit"]["gates"].pop("standard_errors"),
        lambda x: x["recomputed_audit"]["gates"].update(extra=True),
        lambda x: x["recomputed_audit"]["gates"].update(standard_errors=1),
        lambda x: x["recomputed_audit"].pop("final_estimates"),
        lambda x: x["recomputed_audit"]["final_estimates"]["dense"]["mean"].pop(),
        lambda x: x["recomputed_audit"].update(max_standard_error=float("nan")),
    ],
)
def test_reject_fake_or_partial_receipt(change):
    value = example()
    ref = copy.deepcopy(value["recomputed_audit"])
    change(value)
    with pytest.raises((ValueError, KeyError)):
        gate.validate_replay(value, P, ref)


@pytest.mark.parametrize(
    "change",
    [
        lambda x: x.update(platform="darwin"),
        lambda x: x["versions"].update(triqs="0"),
        lambda x: x["modules"].update(triqs="/tmp/triqs.py"),
        lambda x: x["threads"].update(OMP_NUM_THREADS="2"),
    ],
)
def test_runtime_checks(change):
    value = runtime()
    change(value)
    with pytest.raises(ValueError):
        gate.validate_runtime(value, "/env/bin/python")


@contextmanager
def source(*args):
    yield gate.ROOT


def setup(tmp_path):
    root = tmp_path / "payload"
    root.mkdir()
    (root / "cthyb-evidence-replay.json").write_text(json.dumps(example()))
    return root, tmp_path / "out"


def fake_run(cmd, **kwargs):
    Path(cmd[4]).write_text(json.dumps(runtime()))
    Path(cmd[-1]).write_text(json.dumps(example()))
    return SimpleNamespace(returncode=0, stdout="ok", stderr="")


def run(root, out, runner=fake_run, identity=lambda _: {"revision": "c" * 40}):
    return gate.audit(
        root,
        P,
        "/env/bin/python",
        out,
        runner=runner,
        snapshot=source,
        identity_probe=identity,
        read_blob=lambda source, revision, name: (source / name).read_bytes(),
    )


def test_success_uses_separate_receipt_and_logs(tmp_path):
    root, out = setup(tmp_path)
    result = run(root, out)
    assert result["status"] == "passed"
    assert result["payload_sha256"]["cthyb-evidence-replay.json"] == gate.digest(
        root / "cthyb-evidence-replay.json"
    )
    assert (out / "auditor.log").is_file()


def test_nonzero_even_with_valid_output_and_timeout(tmp_path):
    root, out = setup(tmp_path)

    def failed(*args, **kwargs):
        fake_run(*args, **kwargs)
        return SimpleNamespace(returncode=1, stdout="", stderr="error")

    assert run(root, out, failed)["status"] == "failed"

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1, output="partial")

    result = run(root, tmp_path / "out2", timeout)
    assert result["status"] == "failed" and "timed out" in result["error"]
    assert "partial" in (tmp_path / "out2/auditor.log").read_text()


def test_reject_unsafe_or_extra_payload(tmp_path):
    root, out = setup(tmp_path)
    (root / "extra").write_text("bad")
    assert run(root, out)["status"] == "failed"
    (root / "extra").unlink()
    receipt = root / "cthyb-evidence-replay.json"
    contents = receipt.read_text()
    receipt.unlink()
    outside = tmp_path / "outside"
    outside.write_text(contents)
    receipt.symlink_to(outside)
    assert run(root, tmp_path / "out2")["status"] == "failed"


def test_mutated_payload_and_controller_fail(tmp_path):
    root, out = setup(tmp_path)

    def mutate(*args, **kwargs):
        result = fake_run(*args, **kwargs)
        (root / "cthyb-evidence-replay.json").write_text("{}")
        return result

    result = run(root, out, mutate)
    assert result["status"] == "failed" and result["payload_unchanged"] is False
    (root / "cthyb-evidence-replay.json").write_text(json.dumps(example()))
    values = iter([{"revision": "a"}, {"revision": "b"}])
    assert (
        run(root, tmp_path / "out2", identity=lambda _: next(values))["status"]
        == "failed"
    )


def test_bad_runtime_and_wrong_input_hash_fail(tmp_path, monkeypatch):
    root, out = setup(tmp_path)

    def bad_runtime(*args, **kwargs):
        result = fake_run(*args, **kwargs)
        path = Path(args[0][4])
        data = runtime()
        data["versions"]["triqs"] = "invalid"
        path.write_text(json.dumps(data))
        return result

    assert run(root, out, bad_runtime)["status"] == "failed"
    old_digest = gate.digest
    monkeypatch.setattr(
        gate,
        "digest",
        lambda p: "bad" if str(p).endswith(gate.INPUTS[0]) else old_digest(p),
    )
    result = run(root, tmp_path / "out2")
    assert result["status"] == "failed" and "source/input" in result["error"]


def test_hardlinked_receipt_and_nested_output_rejected(tmp_path):
    root, out = setup(tmp_path)
    import os

    os.link(root / "cthyb-evidence-replay.json", tmp_path / "linked")
    assert run(root, out)["status"] == "failed"
    (tmp_path / "linked").unlink()
    with pytest.raises(ValueError):
        run(root, root / "output")


def test_duplicate_key_and_empty_directory_rejected(tmp_path):
    root, out = setup(tmp_path)
    path = root / "cthyb-evidence-replay.json"
    path.write_text('{"generation_revision":"x","generation_revision":"y"}')
    assert run(root, out)["status"] == "failed"
    path.write_text(json.dumps(example()))
    (root / "empty").mkdir()
    assert run(root, tmp_path / "out2")["status"] == "failed"


def test_dot_dot_path_rejected(tmp_path):
    root, out = setup(tmp_path)
    with pytest.raises(ValueError, match="noncanonical"):
        run(root / ".." / "payload", out)


def test_runtime_rejects_unrelated_executable():
    value = runtime()
    value["executable"] = "/elsewhere/bin/python"
    with pytest.raises(ValueError, match="runtime module paths"):
        gate.validate_runtime(value, "/env/bin/python")


def test_scientific_blob_mismatch_rejected(tmp_path):
    root, out = setup(tmp_path)
    result = gate.audit(
        root,
        P,
        "/env/bin/python",
        out,
        runner=fake_run,
        snapshot=source,
        identity_probe=lambda _: {"revision": "c" * 40},
        read_blob=lambda *_: b"wrong",
    )
    assert result["status"] == "failed" and "Git blob" in result["error"]
