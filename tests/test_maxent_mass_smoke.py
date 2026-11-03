"""Failure-path and fixture checks; these do not claim native acceptance."""

import contextlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/providers/maxent_mass_smoke.py"
spec = importlib.util.spec_from_file_location("maxent_mass_smoke", PATH)
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


def test_fixture_is_deterministic_and_physical():
    beta, tau, values = smoke.fixture(np)
    assert beta == 8 and tau.shape == values.shape == (17,)
    assert np.isfinite(values).all() and np.all(values < 0)
    assert np.array_equal(values, smoke.fixture(np)[2])
    assert abs(values[0] + values[-1] + 1) < 1e-14


def test_missing_provider_fails_after_bound_source(monkeypatch):
    helper = SimpleNamespace(
        source_identity=lambda root: {"source_revision": "test-revision"},
        environment_identity=lambda prefix: {},
        provider_output_to_stderr=contextlib.nullcontext,
    )
    monkeypatch.setattr(smoke, "load_helper", lambda *args: helper)
    monkeypatch.setattr(
        smoke.subprocess,
        "check_output",
        lambda cmd: (ROOT / cmd[-1].split(":", 1)[1]).read_bytes(),
    )
    original = smoke.importlib.import_module

    def import_module(name):
        if name == "threadpoolctl":
            return SimpleNamespace()
        if name == "triqs_maxent":
            raise ImportError("provider unavailable")
        return original(name)

    monkeypatch.setattr(smoke.importlib, "import_module", import_module)
    result = smoke.run_probe(ROOT, "test-revision", False)
    assert result["status"] == "failed"
    assert result["error"]["stage"] == "imports"
    assert result["error"]["type"] == "ImportError"


def test_record_hash_mismatch_fails(monkeypatch, tmp_path):
    files = {}
    for key in smoke.FILES:
        path = tmp_path / key
        path.write_bytes(b"source")
        files[key] = path
    record = "".join(f"{name},sha256=wrong,6\n" for name in smoke.FILES.values())
    distribution = SimpleNamespace(read_text=lambda name: record)
    monkeypatch.setattr(
        smoke.importlib.metadata, "distribution", lambda name: distribution
    )
    with pytest.raises(RuntimeError, match="hash or size mismatch"):
        smoke.verify_record_hashes(files)


def test_cli_reports_source_mismatch(tmp_path):
    output = tmp_path / "receipt.json"
    run = subprocess.run(
        [
            sys.executable,
            str(PATH),
            "--source-revision",
            "not-the-revision",
            "--output",
            str(output),
        ],
        capture_output=True,
        check=False,
    )
    assert run.returncode == 1
    receipt = json.loads(output.read_text())
    assert receipt["status"] == "failed"
    assert receipt["error"]["stage"] == "source"


def test_invalid_native_curvature_rejected_before_provider_import(monkeypatch):
    from kernelgecp.integrations import _maxent_native_mass as native

    def fail_import(name):
        raise AssertionError("provider import was attempted")

    monkeypatch.setattr(native.importlib, "import_module", fail_import)
    with pytest.raises(ValueError, match="curvature"):
        native.make_native_mass_cost(1.0, curvature="invalid")


@pytest.mark.parametrize("fraction", [-1.0, 1.0, float("nan"), True, 1e-20])
def test_invalid_interior_fraction_rejected_before_provider_import(
    monkeypatch, fraction
):
    from kernelgecp.integrations import _maxent_native_mass as native

    def fail_import(name):
        raise AssertionError("provider import was attempted")

    monkeypatch.setattr(native.importlib, "import_module", fail_import)
    with pytest.raises(ValueError, match="interior_fraction"):
        native.make_native_mass_cost(1.0, interior_fraction=fraction)
