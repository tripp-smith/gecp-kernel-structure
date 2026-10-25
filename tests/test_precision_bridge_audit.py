"""Independent precision audit identity, rejection gates and tiny math oracle."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np
import pytest

from scripts import audit_precision_bridge as audit

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def canonical() -> dict[str, Any]:
    return json.loads(
        gzip.decompress(
            (ROOT / "experiments/data/precision_bridge.json.gz").read_bytes()
        )
    )


def packed(value: Any) -> dict[str, Any]:
    array = np.asarray(value)
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "real": array.real.tolist(),
        "imag": array.imag.tolist() if np.iscomplexobj(array) else None,
        "sha256": audit.array_hash(array),
    }


@pytest.mark.parametrize("compressed", [False, True])
def test_reader_binds_both_container_and_payload(
    tmp_path: Path, compressed: bool
) -> None:
    payload = b'{"schema_version":1,"records":[]}\n'
    container = gzip.compress(payload, mtime=0) if compressed else payload
    path = tmp_path / "artifact"
    path.write_bytes(container)
    result, compressed_hash, payload_hash = audit.read_artifact(
        path, audit.sha256(container), audit.sha256(payload)
    )
    assert result == {"schema_version": 1, "records": []}
    assert compressed_hash == audit.sha256(container)
    assert payload_hash == audit.sha256(payload)
    for expected_container, expected_payload in (
        ("0" * 64, audit.sha256(payload)),
        (audit.sha256(container), "0" * 64),
        ("bad", audit.sha256(payload)),
    ):
        with pytest.raises(audit.AuditError):
            audit.read_artifact(path, expected_container, expected_payload)


@pytest.mark.parametrize(
    "payload",
    [b"{broken", b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":1e999}'],
)
def test_reader_rejects_malformed_duplicate_and_nonfinite_json(
    tmp_path: Path, payload: bytes
) -> None:
    path = tmp_path / "bad.json"
    path.write_bytes(payload)
    with pytest.raises(audit.AuditError):
        audit.read_artifact(path, audit.sha256(payload), audit.sha256(payload))


@pytest.mark.parametrize("compressed", [False, True])
def test_reader_has_real_container_and_decompression_budgets(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, compressed: bool
) -> None:
    monkeypatch.setattr(audit, "MAX_ARTIFACT_BYTES", 100)
    payload = json.dumps({"data": "0" * 1000}).encode()
    container = gzip.compress(payload, mtime=0) if compressed else payload
    path = tmp_path / "oversized"
    path.write_bytes(container)
    with pytest.raises(audit.AuditError, match="size budget"):
        audit.read_artifact(path, audit.sha256(container), audit.sha256(payload))


def test_complete_complex_arrays_are_preserved_and_imaginary_bytes_matter() -> None:
    value = np.array([[1 + 2j, -0.5j], [3, 4 - 0.125j]])
    payload = packed(value)
    np.testing.assert_array_equal(audit.decode_array(payload), value)
    payload["imag"][0][0] = -2
    with pytest.raises(audit.AuditError, match="hash"):
        audit.decode_array(payload)


@pytest.mark.parametrize(
    "change",
    [
        {"dtype": "<f4"},
        {"dtype": "O"},
        {"shape": [4]},
        {"shape": [True, 2]},
        {"imag": [0]},
        {"real": [[True, 1], [1, 1]]},
        {"real": [["1", "1"], ["1", "1"]]},
        {"sha256": "0" * 64},
    ],
)
def test_array_decoder_rejects_dtype_shape_component_and_hash_forgery(
    change: dict[str, Any],
) -> None:
    payload = packed(np.ones((2, 2), dtype=np.complex128))
    payload.update(change)
    with pytest.raises(audit.AuditError):
        audit.decode_array(payload)


def test_all_canonical_record_bindings_and_honest_inconclusive_status(
    canonical: dict[str, Any],
) -> None:
    count = 0
    for record, config in zip(
        canonical["records"], canonical["config"]["cases"], strict=True
    ):
        arrays, grid, _ = audit.validate_record(record, config)
        assert len(arrays) == 13 and grid[2].shape == (5, 5)
        count += record["fit_certificate"]["certified"]
    assert count == 5
    failure = canonical["records"][2]
    audit.validate_fit_status(failure)
    assert (
        failure["fit_certificate"]["termination_reason"] == "normal_pivot_contains_zero"
    )
    assert failure["fit_certificate"]["coefficients"] is None


@pytest.mark.parametrize(
    "field,value",
    [
        ("certified", False),
        ("certified", 1),
        ("termination_reason", "normal_pivot_contains_zero"),
        ("continuous_certificate", True),
        ("physical_truth_certificate", True),
        ("objective_scope", "continuous_physical_truth"),
        ("coefficients", None),
    ],
)
def test_forged_success_flags_and_scopes_are_rejected(
    canonical: dict[str, Any], field: str, value: Any
) -> None:
    record = copy.deepcopy(canonical["records"][0])
    record["fit_certificate"][field] = value
    with pytest.raises(audit.AuditError):
        audit.validate_fit_status(record)


def test_unknown_failure_and_failure_with_success_data_are_rejected(
    canonical: dict[str, Any],
) -> None:
    record = copy.deepcopy(canonical["records"][2])
    record["fit_certificate"]["termination_reason"] = "all_done"
    with pytest.raises(audit.AuditError, match="failure status"):
        audit.validate_fit_status(record)
    record = copy.deepcopy(canonical["records"][2])
    record["queries"] = {"time": {"direct_error_upper": 0}}
    with pytest.raises(audit.AuditError, match="query success"):
        audit.validate_fit_status(record)


@pytest.mark.parametrize(
    "change",
    [
        {"eta": 1.0},
        {"eta": 0.0},
        {"certified": False},
        {"continuous_certificate": True},
        {"candidate_index": [0, 0]},
        {"content_sha256": "0" * 64},
        {"true_global_upper": 0.0},
        {"termination_reason": "exact_zero_grid_residual"},
    ],
)
def test_grid_eta_flags_and_identity_forgery_are_rejected(
    canonical: dict[str, Any], change: dict[str, Any]
) -> None:
    record = copy.deepcopy(canonical["records"][0])
    record["grid_certificate"].update(change)
    arrays = {k: audit.decode_array(v) for k, v in record["inputs"].items()}
    with pytest.raises(audit.AuditError):
        audit.validate_grid_binding(record, arrays)


def test_generation_provenance_uses_historical_git_blobs(
    canonical: dict[str, Any],
) -> None:
    count = audit.validate_provenance(
        canonical, ROOT, canonical["source"]["git_commit"]
    )
    assert count == len(canonical["source"]["source_sha256"])
    assert count >= 40


@pytest.mark.parametrize("mutation", ["commit", "tree", "blob", "coverage", "config"])
def test_generation_source_and_config_forgery_are_rejected(
    canonical: dict[str, Any], mutation: str
) -> None:
    bad = copy.deepcopy(canonical)
    if mutation == "commit":
        bad["source"]["git_commit"] = "0" * 40
    elif mutation == "tree":
        bad["source"]["git_tree"] = "0" * 40
    elif mutation == "blob":
        bad["source"]["source_sha256"]["kernelgecp/precision.py"] = "0" * 64
    elif mutation == "coverage":
        bad["source"]["source_sha256"].pop("kernelgecp/precision.py")
    else:
        bad["config"]["cases"][0]["beta"] = 3.8
    with pytest.raises(audit.AuditError):
        audit.validate_provenance(bad, ROOT, canonical["source"]["git_commit"])


def test_wrong_root_and_noncanonical_commit_are_rejected(
    canonical: dict[str, Any], tmp_path: Path
) -> None:
    with pytest.raises(audit.AuditError):
        audit.validate_provenance(
            canonical, tmp_path, canonical["source"]["git_commit"]
        )
    with pytest.raises(audit.AuditError, match="full lowercase SHA"):
        audit.validate_provenance(canonical, ROOT, "HEAD")


def tiny_grid() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    beta, wmax = 2.0, 1.0
    times = beta * np.array([0.071, 0.23, 0.41, 0.63, 0.91])
    omega = 2 * np.array([-0.91, -0.41, 0.03, 0.37, 0.83])
    arrays = {
        "basis_parameters": np.array([beta, wmax, beta * wmax]),
        "basis_tau_nodes": np.array([0.0]),
        "basis_omega_nodes": np.array([0.0]),
        "query_tau": times,
        "actual_time_design": np.ones((5, 1)),
    }
    rounded = audit.kernel_binary64(
        times[:, None] / beta, omega[None, :]
    ) - audit.kernel_binary64(0.0, omega[None, :])
    location = np.unravel_index(np.argmax(np.abs(rounded)), rounded.shape)
    magnitude = float(abs(rounded[location]))
    discrepancy = 1e-12
    lower = np.nextafter(magnitude - discrepancy, -np.inf)
    upper = np.nextafter(magnitude + discrepancy, np.inf)
    eta = np.nextafter(
        float(
            (Fraction(magnitude) - Fraction(discrepancy))
            / (Fraction(magnitude) + Fraction(discrepancy))
        ),
        0.0,
    )
    digest = hashlib.sha256(b"kernelgecp-grid-precision-v1")
    for array in (
        np.array([beta, wmax]),
        arrays["basis_tau_nodes"],
        arrays["basis_omega_nodes"],
        times,
        omega,
    ):
        digest.update(str(array.shape).encode())
        digest.update(array.dtype.str.encode())
        digest.update(array.tobytes())
    digest.update(rounded.dtype.str.encode())
    digest.update(rounded.tobytes())
    ties = int(np.count_nonzero(abs(rounded) == magnitude))
    record = {
        "beta": beta,
        "grid_certificate": {
            "certified": True,
            "continuous_certificate": False,
            "objective_scope": "exact_physical_selected_cross_on_supplied_finite_grid",
            "content_sha256": digest.hexdigest(),
            "grid_shape": [5, 5],
            "candidate_index": list(location),
            "candidate_tau": float(times[location[0]]),
            "candidate_omega": float(omega[location[1]]),
            "rounded_candidate_magnitude": magnitude,
            "rounded_global_upper": magnitude,
            "discrepancy_upper": discrepancy,
            "true_candidate_lower": lower,
            "true_global_upper": upper,
            "eta": eta,
            "termination_reason": "positive_finite_grid_quality",
            "exact_tie_count": ties,
            "near_tie_count": ties,
            "tie_rtol": 0,
        },
    }
    return record, arrays


def test_tiny_independent_bordered_determinant_success_and_understated_error() -> None:
    record, arrays = tiny_grid()
    grid = audit.validate_grid_binding(record, arrays)
    context = mp.mp.clone()
    context.prec = audit.AUDIT_BITS
    result, count = audit.audit_grid_math(record, arrays, grid, context, {})
    assert count == 25
    assert float(result["independent_quality"]) >= record["grid_certificate"]["eta"]
    bad = copy.deepcopy(record)
    bad["grid_certificate"].update(discrepancy_upper=0.0, eta=1.0)
    with pytest.raises(audit.AuditError, match="discrepancy"):
        audit.audit_grid_math(bad, arrays, grid, context, {})


def test_cli_hash_failure_exits_nonzero_without_receipt(tmp_path: Path) -> None:
    source, output = tmp_path / "bad.json", tmp_path / "receipt.json"
    source.write_text("{}")
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_precision_bridge.py"),
            "--source",
            str(source),
            "--output",
            str(output),
            "--project-root",
            str(ROOT),
            "--generation-commit",
            "0" * 40,
            "--expected-sha256",
            "0" * 64,
            "--expected-payload-sha256",
            "0" * 64,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1 and "container hash mismatch" in result.stderr
    assert not output.exists()
