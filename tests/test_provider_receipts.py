"""Provider gates reject setup-only, stale and incompletely checked runtimes."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.providers.verify_receipt import (
    EXPECTED,
    REQUIRED_CHECKS,
    ProviderReceiptError,
    read_receipt,
    validate_receipt,
)

SHA = "a" * 40


def _receipt(provider: str = "gp") -> dict:
    receipt = {
        "schema": "kernelgecp.provider-smoke.v1",
        "provider": provider,
        "source_revision": SHA,
        "status": "passed",
        "checks": {
            name: True for name in {"imports", "numerical"} | REQUIRED_CHECKS[provider]
        },
        "versions": {"python": "3.12.12", **EXPECTED[provider]},
    }
    receipt["metrics"] = {
        "solve_max_error": 4e-15,
        "logdet_error": 4e-16,
        "log_probability_error": 9e-16,
        "triqs_resolvent_max_error": 3e-17,
        "cthyb_cycles": 200,
        "cthyb_average_sign_real": 1.0,
        "maxent_alpha_count": 3,
        "maxent_min_chi2": 15.5,
    }
    receipt["measurements"] = {
        "iteration": 20,
        "energy_change_sum": 1e-12,
        "green_change_max": 1e-8,
        "correlation_energy_hartree": -0.02,
        "dynamic_selfenergy_max": 0.1,
        "endpoint_normalization_error": 1e-9,
        "electron_count": 2.0,
    }
    receipt["artifacts"] = dict.fromkeys(
        [
            "green/input.log",
            "green/gf2.log",
            "green/input.h5",
            "green/sim.h5",
            "green/run-manifest.json",
            "green/build-provenance.json",
        ],
        "b" * 64,
    )
    receipt["calculation"] = {
        "electron_phonon_coupling": 0.7,
        "native_frequency_count": 30,
        "lowest_frequency_gap_ev": 0.002,
        "temperature_kelvin": 15,
        "phonon_irreducible_qpoints": 4,
    }
    receipt["native_executions"] = [
        {
            "returncode": 0,
            "mpi_ranks": 1,
            "omp_threads": 1,
            "input_sha256": "b" * 64,
            "log_sha256": "b" * 64,
            "executable_sha256": "c" * 64,
        }
        for _ in range(5)
    ]
    return receipt


@pytest.mark.parametrize("provider", EXPECTED)
def test_provider_receipt_checks_actual_versions_and_numerical_status(
    provider: str,
) -> None:
    receipt = _receipt(provider)
    assert validate_receipt(receipt, provider, SHA) is receipt


@pytest.mark.parametrize(
    "key,value",
    [
        ("schema", "unknown"),
        ("provider", "triqs"),
        ("source_revision", "b" * 40),
        ("status", "failed"),
        ("status", "skipped"),
        ("status", "resolved"),
        ("checks", {}),
        ("checks", {"imports": True}),
        ("checks", {"imports": True, "numerical": False}),
        ("checks", {"imports": True, "numerical": 1}),
        ("checks", {"imports": True, "numerical": True, "extra": False}),
        ("versions", {}),
    ],
)
def test_provider_receipt_rejects_missing_failed_or_stale_checks(
    key: str, value: object
) -> None:
    receipt = _receipt()
    receipt[key] = value
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, "gp", SHA)


@pytest.mark.parametrize("package", ["python", *EXPECTED["gp"]])
def test_provider_receipt_rejects_unverified_version(package: str) -> None:
    receipt = copy.deepcopy(_receipt())
    receipt["versions"][package] = "999.0"
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, "gp", SHA)


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"metric":NaN}', "[] trailing"])
def test_provider_receipt_rejects_ambiguous_json(tmp_path: Path, text: str) -> None:
    path = tmp_path / "receipt.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_receipt(path)


def test_provider_receipt_cli_failure_is_nonzero(tmp_path: Path) -> None:
    path = tmp_path / "failed.json"
    receipt = _receipt()
    receipt["status"] = "failed"
    path.write_text(json.dumps(receipt))
    result = subprocess.run(
        [
            sys.executable,
            "scripts/providers/verify_receipt.py",
            str(path),
            "--provider",
            "gp",
            "--source-revision",
            SHA,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "rejected" in result.stderr


def test_optional_gp_import_failure_cannot_emit_success(tmp_path: Path) -> None:
    path = tmp_path / "gp.json"
    result = subprocess.run(
        [
            sys.executable,
            "-S",  # Deliberately disable optional site packages in this child.
            "scripts/providers/gp_smoke.py",
            "--output",
            str(path),
            "--source-revision",
            SHA,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    receipt = json.loads(path.read_text())
    assert receipt["status"] == "failed"
    assert "ModuleNotFoundError" in receipt["error"]
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, "gp", SHA)


def test_only_torch_accepts_the_supported_cpu_version_suffix() -> None:
    receipt = _receipt()
    receipt["versions"]["torch"] += "+cpu"
    validate_receipt(receipt, "gp", SHA)
    receipt["versions"]["numpy"] += "+cpu"
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, "gp", SHA)


@pytest.mark.parametrize("value", [None, -1, 99, float("inf"), float("nan"), True])
def test_gp_receipt_rejects_missing_or_failed_numeric_evidence(value: object) -> None:
    receipt = _receipt()
    receipt["metrics"]["solve_max_error"] = value
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, "gp", SHA)


@pytest.mark.parametrize("provider", EXPECTED)
def test_domain_checks_cannot_be_replaced_by_import_success(provider: str) -> None:
    receipt = _receipt(provider)
    receipt["checks"] = {"imports": True, "numerical": True}
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, provider, SHA)


@pytest.mark.parametrize("provider", EXPECTED)
def test_domain_numeric_evidence_cannot_be_omitted(provider: str) -> None:
    receipt = _receipt(provider)
    for field in ("metrics", "measurements", "calculation", "native_executions"):
        receipt.pop(field)
    with pytest.raises(ProviderReceiptError):
        validate_receipt(receipt, provider, SHA)
