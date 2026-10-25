"""Independent application-audit arithmetic, identity, and failure handling."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import subprocess
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np
import pytest

from scripts import audit_dlr_applications as audit


def _payload(value: Any) -> dict[str, Any]:
    array = np.asarray(value)
    array = np.asarray(array, dtype=array.dtype.newbyteorder("<"), order="C")
    header = json.dumps(
        {"shape": array.shape, "dtype": array.dtype.str}, sort_keys=True
    ).encode()
    result: dict[str, Any] = {
        "shape": list(array.shape),
        "dtype": array.dtype.str,
        "sha256": hashlib.sha256(header + array.tobytes()).hexdigest(),
    }
    if np.iscomplexobj(array):
        result.update(real=array.real.tolist(), imag=array.imag.tolist())
    else:
        result["values"] = array.tolist()
    return result


def _error_record(predicted: Any, truth: Any, target: float) -> dict[str, Any]:
    predicted_array = np.asarray(predicted)
    truth_array = np.asarray(truth)
    absolute = float(np.max(np.abs(predicted_array - truth_array)))
    scale = float(np.max(np.abs(truth_array)))
    return {
        "predicted": _payload(predicted_array),
        "truth": _payload(truth_array),
        "absolute_error": absolute,
        "relative_error": absolute / scale if scale else None,
        "truth_scale": scale,
        "accuracy_passed": absolute <= target,
    }


def test_lift_binary64_preserves_actual_float_instead_of_short_decimal() -> None:
    with mp.workdps(95):
        value = audit.lift_binary64(Decimal("0.1"))
        numerator, denominator = (0.1).as_integer_ratio()
        assert value == mp.mpf(numerator) / denominator
        assert value != mp.mpf("0.1")
        assert value == audit.lift_binary64(0.1)


@pytest.mark.parametrize("compressed", [False, True])
def test_artifact_reader_checks_container_and_uncompressed_identity(
    tmp_path: Path, compressed: bool
) -> None:
    raw = b'{"schema_version": 1, "records": []}\n'
    container = gzip.compress(raw, mtime=0) if compressed else raw
    path = tmp_path / ("audit.json.gz" if compressed else "audit.json")
    path.write_bytes(container)
    decoded, container_hash, payload_hash = audit.read_artifact(
        path, hashlib.sha256(container).hexdigest()
    )
    assert decoded == {"schema_version": 1, "records": []}
    assert container_hash == hashlib.sha256(container).hexdigest()
    assert payload_hash == hashlib.sha256(raw).hexdigest()
    if compressed:
        assert container_hash != payload_hash


def test_artifact_reader_rejects_wrong_hash_and_malformed_json(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_bytes(b"{broken")
    with pytest.raises(audit.AuditError):
        audit.read_artifact(path, "0" * 64)
    with pytest.raises(audit.AuditError):
        audit.read_artifact(path, hashlib.sha256(path.read_bytes()).hexdigest())


@pytest.mark.parametrize(
    "values",
    [np.array([1.0, -0.25]), np.array([[1 + 2j, -0.5j]], dtype=np.complex128)],
)
def test_array_decoder_preserves_complex_values_and_validates_content(
    values: Any,
) -> None:
    payload = _payload(values)
    np.testing.assert_array_equal(audit.decode_array(payload), values)
    bad_hash = copy.deepcopy(payload)
    bad_hash["sha256"] = "0" * 64
    with pytest.raises(audit.AuditError):
        audit.decode_array(bad_hash)
    bad_shape = copy.deepcopy(payload)
    bad_shape["shape"] = [np.asarray(values).size + 1]
    with pytest.raises(audit.AuditError):
        audit.decode_array(bad_shape)


@pytest.mark.parametrize("dtype", ["O", "bool", "<U4"])
def test_array_decoder_rejects_non_numeric_payload_types(dtype: str) -> None:
    payload = _payload(np.array([1.0]))
    payload["dtype"] = dtype
    with pytest.raises(audit.AuditError):
        audit.decode_array(payload)


def test_error_audit_retains_honest_scientific_failure() -> None:
    record = _error_record([0.875], [1.0], target=0.01)
    report: list[dict[str, Any]] = []
    with mp.workdps(95):
        assert not audit.verify_error(
            record,
            [mp.mpf(1)],
            mp.mpf("0.01"),
            "honest_failure",
            report,
            expected_shape=(1,),
        )
    assert len(report) == 1
    assert report[0]["accuracy_passed"] is False


@pytest.mark.parametrize("corruption", ["flag", "shape", "reported_error", "truth"])
def test_error_audit_rejects_inconsistent_flags_shapes_and_physical_truth(
    corruption: str,
) -> None:
    record = _error_record([1.0], [1.0], target=1e-8)
    expected_shape = (1,)
    if corruption == "flag":
        record["accuracy_passed"] = False
    elif corruption == "shape":
        expected_shape = (1, 1)
    elif corruption == "reported_error":
        record["absolute_error"] = 0.1
    else:
        record["truth"] = _payload([0.5])
    with mp.workdps(95), pytest.raises(audit.AuditError):
        audit.verify_error(
            record,
            [mp.mpf(1)],
            mp.mpf("1e-8"),
            "corrupted",
            [],
            expected_shape=expected_shape,
        )


def test_expected_workloads_scale_fixed_fixtures_by_binary64_wmax() -> None:
    case = {"beta": 3.7, "wmax": 0.1, "target": 1e-8, "solver": False}
    expected = audit.expected_workloads(case)
    assert [item["name"] for item in expected] == [
        "central_pole",
        "asymmetric_satellites",
        "gapped",
    ]
    assert expected[0]["physical_poles"] == [0.0]
    assert expected[0]["spectral_weights"] == [1.0]
    assert expected[1]["physical_poles"] == [
        0.1 * fraction for fraction in (-0.75, 0.08, 0.55)
    ]
    assert expected[1]["spectral_weights"] == [0.2, 0.5, 0.3]
    assert expected[2]["physical_poles"] == [
        0.1 * fraction for fraction in (-0.8, -0.3, 0.35, 0.9)
    ]
    assert expected[2]["spectral_weights"] == [0.15, 0.35, 0.3, 0.2]


def test_semicircle_root_is_causal_and_satisfies_equation_for_both_signs() -> None:
    with mp.workdps(95):
        beta, hopping = mp.mpf("2.5"), mp.mpf("0.7")
        for index in (-4, -1, 0, 3):
            z = mp.j * (2 * index + 1) * mp.pi / beta
            result = audit.semicircle_resolvent(index, beta, hopping)
            assert mp.im(result) * mp.im(z) < 0
            assert abs(result - 1 / (z - hopping**2 * result)) < mp.mpf("1e-85")
        assert abs(
            audit.semicircle_resolvent(-4, beta, hopping)
            - mp.conj(audit.semicircle_resolvent(3, beta, hopping))
        ) < mp.mpf("1e-85")


def test_semicircle_energy_integral_agrees_with_independent_angle_formulation() -> None:
    with mp.workdps(95):
        beta, hopping = mp.mpf("0.5"), mp.mpf("0.7")
        tau = beta / 2
        actual = audit.semicircle_energy_integral(tau, beta, hopping)

        def angle_integrand(theta: Any) -> Any:
            energy = 2 * hopping * mp.cos(theta)
            return (
                -2
                * mp.sin(theta) ** 2
                / mp.pi
                * mp.exp(-tau * energy)
                / (1 + mp.exp(-beta * energy))
            )

        expected = mp.quad(angle_integrand, [0, mp.pi / 2, mp.pi])
        assert abs(actual - expected) < mp.mpf("1e-80")


def test_cli_requires_explicit_source_and_rejects_bad_container_hash(
    tmp_path: Path,
) -> None:
    path = tmp_path / "input.json"
    path.write_text("{}\n")
    output = tmp_path / "report.json"
    script = Path(audit.__file__).resolve()
    command = [
        sys.executable,
        str(script),
        "--input",
        str(path),
        "--output",
        str(output),
        "--root",
        str(tmp_path),
        "--expected-sha256",
        "0" * 64,
    ]
    missing_source = subprocess.run(
        command, capture_output=True, text=True, check=False
    )
    assert missing_source.returncode != 0
    assert "--source" in missing_source.stderr
    assert not output.exists()
    bad_hash = subprocess.run(
        [*command, "--source", "1" * 40],
        capture_output=True,
        text=True,
        check=False,
    )
    assert bad_hash.returncode != 0
    assert not output.exists()
