"""Independently audit frozen finite-query arithmetic bridge evidence.

No kernelgecp, benchmark or certificate-engine module is imported. Complete
arrays and generation-era Git blobs are checked before 1280-bit centered-kernel
bordered determinants and exact selected-node coefficient checks. These are
independent consistency checks, not interval proofs or verified runtime code.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import io
import json
import math
import re
import subprocess
import sys
import tomllib
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

AUDIT_BITS = 1280
MAX_ARTIFACT_BYTES = 64 * 1024**2
HEX256 = re.compile(r"[0-9a-f]{64}\Z")
HEX_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
INPUT_NAMES = (
    "basis_parameters",
    "basis_tau_nodes",
    "basis_omega_nodes",
    "training_tau",
    "observations",
    "coefficients",
    "query_tau",
    "matsubara_indices",
    "actual_training_design",
    "actual_time_design",
    "actual_frequency_design",
    "actual_time_output",
    "actual_frequency_output",
)
FIT_FAILURES = {
    "rank_budget",
    "work_budget",
    "underdetermined_training_design",
    "core_pivot_contains_zero",
    "normal_pivot_contains_zero",
    "bound_export_overflow",
    "arithmetic_inconclusive:ArithmeticError",
    "arithmetic_inconclusive:ValueError",
    "arithmetic_inconclusive:NumericalBreakdownError",
    "arithmetic_inconclusive:FloatingPointError",
    "arithmetic_inconclusive:OverflowError",
    "arithmetic_inconclusive:ZeroDivisionError",
}


class AuditError(ValueError):
    """Evidence is malformed, inconsistent or cannot be independently checked."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _object(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def _nonfinite(value: str) -> None:
    raise AuditError(f"nonfinite JSON constant: {value}")


def _finite_tree(value: Any) -> None:
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON number")
    elif isinstance(value, dict):
        for item in value.values():
            _finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            _finite_tree(item)


def read_artifact(
    path: Path, expected_sha256: str, expected_payload_sha256: str
) -> tuple[dict[str, Any], str, str]:
    """Check both immutable identities before accepting any JSON evidence."""
    for expected in (expected_sha256, expected_payload_sha256):
        require(
            isinstance(expected, str) and HEX256.fullmatch(expected) is not None,
            "expected SHA256 must be 64 lowercase hexadecimal characters",
        )
    with path.open("rb") as stream:
        container = stream.read(MAX_ARTIFACT_BYTES + 1)
    require(len(container) <= MAX_ARTIFACT_BYTES, "artifact exceeds size budget")
    require(sha256(container) == expected_sha256, "artifact container hash mismatch")
    try:
        if container[:2] == b"\x1f\x8b":
            with gzip.GzipFile(fileobj=io.BytesIO(container)) as stream:
                payload = stream.read(MAX_ARTIFACT_BYTES + 1)
        else:
            payload = container
        require(len(payload) <= MAX_ARTIFACT_BYTES, "payload exceeds size budget")
        require(
            sha256(payload) == expected_payload_sha256, "artifact payload hash mismatch"
        )
        result = json.loads(
            payload, object_pairs_hook=_object, parse_constant=_nonfinite
        )
    except (OSError, EOFError, UnicodeError, json.JSONDecodeError) as exc:
        raise AuditError("malformed JSON/gzip artifact") from exc
    require(isinstance(result, dict), "artifact must be a JSON object")
    _finite_tree(result)
    return result, sha256(container), sha256(payload)


def array_hash(array: np.ndarray) -> str:
    array = np.asarray(array, dtype=array.dtype.newbyteorder("<"), order="C")
    header = json.dumps(
        {"dtype": array.dtype.str, "shape": array.shape}, sort_keys=True
    ).encode()
    return sha256(header + array.tobytes())


def _contains_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return True
    return isinstance(value, list) and any(_contains_bool(item) for item in value)


def decode_array(payload: Any) -> np.ndarray:
    """Validate every component, dtype, shape and byte of a packed array."""
    require(isinstance(payload, dict), "array payload must be an object")
    require(
        set(payload) == {"dtype", "shape", "real", "imag", "sha256"},
        "unexpected array payload fields",
    )
    dtype, shape = payload["dtype"], payload["shape"]
    require(dtype in ("<f8", "<c16", "<i8"), "unsupported array dtype")
    require(
        isinstance(shape, list)
        and all(type(size) is int and size >= 0 for size in shape),
        "invalid array shape",
    )
    try:
        require(not _contains_bool(payload["real"]), "boolean array component")
        raw = np.asarray(payload["real"])
        require(raw.dtype.kind in "iuf", "nonnumeric real array component")
        if dtype == "<c16":
            require(not _contains_bool(payload["imag"]), "boolean array component")
            imag = np.asarray(payload["imag"])
            require(imag.dtype.kind in "iuf", "nonnumeric imaginary array component")
            require(
                raw.shape == imag.shape == tuple(shape), "complex array shape mismatch"
            )
            array = np.empty(tuple(shape), dtype="<c16")
            array.real, array.imag = raw, imag
        else:
            require(payload["imag"] is None, "real array has imaginary components")
            if dtype == "<i8":
                require(raw.dtype.kind in "iu", "integer array has noninteger values")
                require(
                    np.all(raw >= np.iinfo(np.int64).min)
                    and np.all(raw <= np.iinfo(np.int64).max),
                    "integer array exceeds int64",
                )
            array = np.asarray(raw, dtype=dtype)
        require(array.shape == tuple(shape), "array shape mismatch")
        require(bool(np.all(np.isfinite(array))), "nonfinite array component")
        require(array_hash(array) == payload["sha256"], "array content hash mismatch")
    except (TypeError, OverflowError) as exc:
        raise AuditError("invalid numeric array") from exc
    return array


def _git(root: Path, *args: str) -> bytes:
    try:
        return subprocess.run(
            ["git", *args], cwd=root, check=True, capture_output=True
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise AuditError("generation Git object unavailable") from exc


def validate_provenance(
    artifact: dict[str, Any], root: Path, generation_commit: str
) -> int:
    """Read historical blobs, never substitute the current checkout/lockfile."""
    require(
        HEX_COMMIT.fullmatch(generation_commit) is not None,
        "generation commit must be a full lowercase SHA",
    )
    require(root.is_dir(), "project root must be an existing directory")
    actual_root = Path(_git(root, "rev-parse", "--show-toplevel").decode().strip())
    require(actual_root.resolve() == root.resolve(), "project root is not the Git root")
    source = artifact.get("source")
    require(isinstance(source, dict), "missing generation provenance")
    require(source.get("git_commit") == generation_commit, "generation commit mismatch")
    tree = _git(root, "rev-parse", generation_commit + "^{tree}").decode().strip()
    require(source.get("git_tree") == tree, "generation tree mismatch")
    paths = (
        _git(root, "ls-tree", "-r", "--name-only", generation_commit)
        .decode()
        .splitlines()
    )
    expected = {
        name
        for name in paths
        if re.fullmatch(r"(?:kernelgecp|experiments)/[^/]+\.py", name)
    } | {"examples/precision_bridge.py", "pyproject.toml", "uv.lock"}
    manifest = source.get("source_sha256")
    require(isinstance(manifest, dict), "missing source hash manifest")
    require(set(manifest) == expected, "source hash manifest coverage mismatch")
    for name, expected_hash in manifest.items():
        require(
            isinstance(expected_hash, str)
            and HEX256.fullmatch(expected_hash) is not None,
            "invalid source hash",
        )
        require(
            sha256(_git(root, "show", generation_commit + ":" + name)) == expected_hash,
            "generation source blob hash mismatch: " + name,
        )
    config = _git(
        root, "show", generation_commit + ":experiments/precision_bridge.toml"
    )
    require(sha256(config) == source.get("config_sha256"), "config blob hash mismatch")
    try:
        require(
            tomllib.loads(config.decode()) == artifact.get("config"), "config mismatch"
        )
        project = tomllib.loads(
            _git(root, "show", generation_commit + ":pyproject.toml").decode()
        )
        lock = tomllib.loads(
            _git(root, "show", generation_commit + ":uv.lock").decode()
        )
    except (UnicodeError, tomllib.TOMLDecodeError) as exc:
        raise AuditError("malformed generation TOML") from exc
    versions = artifact.get("versions")
    require(
        isinstance(versions, dict)
        and set(versions) == {"kernelgecp", "mpmath", "numpy", "scipy"},
        "invalid generation versions",
    )
    require(
        versions["kernelgecp"] == project["project"]["version"],
        "project version mismatch",
    )
    packages = lock.get("package", [])
    for name in ("mpmath", "numpy", "scipy"):
        require(
            any(
                p.get("name") == name and p.get("version") == versions[name]
                for p in packages
            ),
            "generation dependency version absent from frozen lock: " + name,
        )
    return len(manifest)


def kernel_binary64(time: Any, omega: Any) -> np.ndarray:
    """Reproduce the frozen rounded input solely to check its content binding."""
    time, omega = np.broadcast_arrays(
        np.asarray(time, dtype=np.float64), np.asarray(omega, dtype=np.float64)
    )
    result = np.empty_like(omega)
    positive = omega >= 0
    with np.errstate(all="ignore"):
        result[positive] = np.exp(-time[positive] * omega[positive]) / (
            1 + np.exp(-omega[positive])
        )
        result[~positive] = np.exp((1 - time[~positive]) * omega[~positive]) / (
            1 + np.exp(omega[~positive])
        )
    return result


def validate_fit_status(record: dict[str, Any]) -> None:
    certificate = record.get("fit_certificate")
    require(isinstance(certificate, dict), "missing fit certificate")
    require(type(certificate.get("certified")) is bool, "invalid fit certified flag")
    reason = certificate.get("termination_reason")
    require(
        certificate.get("continuous_certificate") is False
        and certificate.get("physical_truth_certificate") is False,
        "invalid fit scope flags",
    )
    require(
        certificate.get("objective_scope")
        == "finite_query_exact_physical_full_column_rank_least_squares",
        "unsupported fit objective scope",
    )
    if certificate["certified"]:
        require(reason == "enclosed", "certified fit has a failure status")
        require(
            all(
                isinstance(certificate.get(k), dict)
                for k in (
                    "core",
                    "training_design",
                    "coefficients",
                    "time",
                    "matsubara",
                )
            ),
            "certified fit lacks completed enclosures",
        )
        require(
            record.get("independent_qr_performed") is True,
            "missing recorded independent check",
        )
        require(
            type(record.get("independent_qr_bits")) is int
            and record["independent_qr_bits"] >= 1024,
            "recorded independent oracle precision is insufficient",
        )
        require(
            set(record.get("queries", {})) == {"time", "matsubara"},
            "missing query results",
        )
    else:
        require(reason in FIT_FAILURES, "unknown fit failure status")
        require(
            all(
                certificate.get(k) is None
                for k in (
                    "core",
                    "training_design",
                    "coefficients",
                    "time",
                    "matsubara",
                )
            ),
            "inconclusive fit contains completed bound claims",
        )
        require(
            record.get("independent_qr_performed") is False,
            "inconclusive fit claims an oracle check",
        )
        require(record.get("queries") == {}, "inconclusive fit contains query success")


def _numeric_bound(value: Any, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value)
    require(array.dtype.kind in "iuf" and array.shape == shape, "invalid " + name)
    require(
        bool(np.all(np.isfinite(array))) and bool(np.all(array >= 0)), "invalid " + name
    )
    return np.asarray(array, dtype=np.float64)


def validate_grid_binding(
    record: dict[str, Any], arrays: dict[str, np.ndarray]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Check finite-grid identity and exact outward transfer arithmetic."""
    grid = record.get("grid_certificate")
    require(isinstance(grid, dict), "missing grid certificate")
    require(
        grid.get("certified") is True,
        "this protocol requires a completed grid certificate",
    )
    require(grid.get("continuous_certificate") is False, "invalid grid scope flag")
    require(
        grid.get("objective_scope")
        == "exact_physical_selected_cross_on_supplied_finite_grid",
        "unsupported grid objective scope",
    )
    beta, wmax, cutoff = arrays["basis_parameters"]
    times = arrays["query_tau"]
    omega = cutoff * np.array([-0.91, -0.41, 0.03, 0.37, 0.83])
    rounded = kernel_binary64(times[:, None] / beta, omega[None, :]) - arrays[
        "actual_time_design"
    ] @ kernel_binary64(arrays["basis_tau_nodes"][:, None] / beta, omega[None, :])
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
    require(
        digest.hexdigest() == grid.get("content_sha256"), "grid content hash mismatch"
    )
    require(grid.get("grid_shape") == [5, 5], "grid shape mismatch")
    location = np.unravel_index(int(np.argmax(np.abs(rounded))), rounded.shape)
    require(grid.get("candidate_index") == list(location), "grid candidate mismatch")
    require(
        grid.get("candidate_tau") == float(times[location[0]]), "candidate tau mismatch"
    )
    require(
        grid.get("candidate_omega") == float(omega[location[1]]),
        "candidate omega mismatch",
    )
    magnitude = float(abs(rounded[location]))
    require(
        grid.get("rounded_candidate_magnitude")
        == grid.get("rounded_global_upper")
        == magnitude,
        "rounded grid magnitude mismatch",
    )
    delta = grid.get("discrepancy_upper")
    require(
        isinstance(delta, (float, int)) and not isinstance(delta, bool) and delta >= 0,
        "invalid grid discrepancy",
    )
    lower = max(Fraction(0), Fraction(magnitude) - Fraction(delta))
    upper = Fraction(magnitude) + Fraction(delta)
    require(
        Fraction(grid["true_candidate_lower"]) <= lower,
        "grid lower bound rounds inward",
    )
    require(
        Fraction(grid["true_global_upper"]) >= upper, "grid upper bound rounds inward"
    )
    if upper:
        require(
            grid.get("termination_reason") == "positive_finite_grid_quality",
            "invalid positive grid status",
        )
        require(
            0 < Fraction(grid["eta"]) <= lower / upper,
            "grid eta rounds inward or is nonpositive",
        )
    else:
        require(
            grid.get("termination_reason") == "exact_zero_grid_residual"
            and grid.get("eta") == 1,
            "invalid zero grid status",
        )
    ties = int(np.count_nonzero(np.abs(rounded) == magnitude))
    require(
        grid.get("exact_tie_count") == grid.get("near_tie_count") == ties,
        "grid tie count mismatch",
    )
    require(grid.get("tie_rtol") == 0, "unsupported grid tie tolerance")
    return times, omega, rounded


def validate_record(
    record: dict[str, Any], config: dict[str, Any]
) -> tuple[dict[str, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray], int]:
    require(
        isinstance(record, dict) and isinstance(config, dict),
        "invalid case record/config",
    )
    rank = record.get("rank")
    require(type(rank) is int and 1 <= rank <= 256, "invalid rank")
    packed = record.get("inputs")
    require(
        isinstance(packed, dict) and set(packed) == set(INPUT_NAMES),
        "input manifest mismatch",
    )
    arrays = {name: decode_array(value) for name, value in packed.items()}
    beta, wmax, cutoff = arrays["basis_parameters"]
    require(
        beta == config.get("beta") == record.get("beta")
        and wmax == config.get("wmax") == record.get("wmax"),
        "case physical parameters mismatch",
    )
    require(
        beta > 0 and wmax > 0 and cutoff == float(beta * wmax) >= 1,
        "invalid physical parameters",
    )
    m = arrays["training_tau"].size
    shapes = {
        "basis_parameters": (3,),
        "basis_tau_nodes": (rank,),
        "basis_omega_nodes": (rank,),
        "training_tau": (m,),
        "observations": (m, 2, 2),
        "coefficients": (rank, 2, 2),
        "query_tau": (5,),
        "matsubara_indices": (5,),
        "actual_training_design": (m, rank),
        "actual_time_design": (5, rank),
        "actual_frequency_design": (5, rank),
        "actual_time_output": (5, 2, 2),
        "actual_frequency_output": (5, 2, 2),
    }
    complex_names = {
        "observations",
        "coefficients",
        "actual_frequency_design",
        "actual_time_output",
        "actual_frequency_output",
    }
    for name, shape in shapes.items():
        require(arrays[name].shape == shape, "input shape mismatch: " + name)
        expected_dtype = (
            "<c16"
            if name in complex_names
            else ("<i8" if name == "matsubara_indices" else "<f8")
        )
        require(
            arrays[name].dtype.str == expected_dtype, "input dtype mismatch: " + name
        )
    require(
        np.array_equal(
            arrays["query_tau"], beta * np.array([0.071, 0.23, 0.41, 0.63, 0.91])
        ),
        "query times mismatch",
    )
    require(
        np.array_equal(arrays["matsubara_indices"], [-17, -2, 0, 4, 31]),
        "Matsubara convention mismatch",
    )
    require(
        np.unique(arrays["basis_tau_nodes"]).size == rank
        and np.unique(arrays["basis_omega_nodes"]).size == rank,
        "repeated basis nodes",
    )
    require(
        bool(
            np.all(
                (arrays["basis_tau_nodes"] >= 0) & (arrays["basis_tau_nodes"] <= beta)
            )
        ),
        "basis time outside domain",
    )
    require(
        bool(np.all(np.abs(arrays["basis_omega_nodes"]) <= cutoff)),
        "basis frequency outside domain",
    )
    training = config.get("training", "uniform")
    require(
        record.get("training") == training
        and training in {"uniform", "resolved", "selected_nodes"},
        "training protocol mismatch",
    )
    require(
        m > 0
        and bool(
            np.all((arrays["training_tau"] >= 0) & (arrays["training_tau"] <= beta))
        ),
        "invalid training domain",
    )
    if training == "uniform":
        require(
            np.array_equal(
                arrays["training_tau"], np.linspace(0, beta, max(40, rank + 12))
            ),
            "uniform training mesh mismatch",
        )
    if training == "selected_nodes":
        require(
            np.array_equal(arrays["training_tau"], arrays["basis_tau_nodes"]),
            "selected training mesh mismatch",
        )
    energies = wmax * np.array([-0.31, 0.23, 0.71])
    exponents = np.where(
        energies >= 0,
        -arrays["training_tau"][:, None] * energies,
        (beta - arrays["training_tau"][:, None]) * energies,
    )
    with np.errstate(all="ignore"):
        observed = -(
            np.exp(exponents)
            / (1 + np.exp(-beta * np.abs(energies)))
            @ np.array([0.25, 0.5, 0.25])
        )[:, None, None] * np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])
    require(
        np.array_equal(observed, arrays["observations"]),
        "physical observation fixture mismatch",
    )
    validate_fit_status(record)
    certificate = record["fit_certificate"]
    hashes = certificate.get("input_hashes")
    require(
        isinstance(hashes, list)
        and all(isinstance(item, list) and len(item) == 2 for item in hashes),
        "invalid certificate input hashes",
    )
    recorded = dict(hashes)
    require(len(recorded) == len(hashes), "duplicate certificate input hash")
    require(
        set(INPUT_NAMES[:8])
        <= set(recorded)
        <= set(INPUT_NAMES) | {"reconstructed_binary64_core"},
        "certificate input coverage mismatch",
    )
    for name, expected in recorded.items():
        if name != "reconstructed_binary64_core":
            require(
                array_hash(arrays[name]) == expected,
                "certificate input hash mismatch: " + name,
            )
    core_count = 0
    if "reconstructed_binary64_core" in recorded:
        core = kernel_binary64(
            arrays["basis_tau_nodes"][:, None] / beta,
            arrays["basis_omega_nodes"][None, :],
        )
        require(
            array_hash(core) == recorded["reconstructed_binary64_core"],
            "rounded core hash mismatch",
        )
        core_count = 1
    if certificate["certified"]:
        require(
            set(recorded) == set(INPUT_NAMES) | {"reconstructed_binary64_core"},
            "completed certificate input coverage mismatch",
        )
        for name in ("time", "matsubara"):
            for field in (
                "direct_error_bound",
                "coefficient_contribution",
                "design_contribution",
                "dot_rounding_bound",
                "total_error_bound",
            ):
                _numeric_bound(certificate[name][field], (5, 2, 2), name + " " + field)
    metadata = {
        key: certificate[key]
        for key in (
            "input_hashes",
            "rank",
            "precision_bits",
            "observation_axis",
            "estimated_work",
        )
    }
    require(
        sha256(json.dumps(metadata, sort_keys=True).encode())
        == certificate.get("content_sha256"),
        "fit certificate content hash mismatch",
    )
    require(
        certificate.get("rank") == rank and certificate.get("observation_axis") == 0,
        "fit certificate rank/axis mismatch",
    )
    require(
        certificate.get("precision_bits") == config.get("precision_bits", 256),
        "certificate precision mismatch",
    )
    return arrays, validate_grid_binding(record, arrays), core_count


def centered_kernel(ctx: Any, time: float, omega: float, beta: float) -> Any:
    time, omega, beta = (
        ctx.mpf(float(time)),
        ctx.mpf(float(omega)),
        ctx.mpf(float(beta)),
    )
    return ctx.exp((ctx.mpf(".5") - time / beta) * omega) / (2 * ctx.cosh(omega / 2))


def audit_grid_math(
    record: dict[str, Any],
    arrays: dict[str, np.ndarray],
    grid_arrays: tuple[np.ndarray, np.ndarray, np.ndarray],
    ctx: Any,
    cache: dict[Any, Any],
) -> tuple[dict[str, Any], int]:
    """Check all true grid values by the bordered-determinant identity."""
    times, omega, rounded = grid_arrays
    beta = record["beta"]
    tau_nodes, nodes = arrays["basis_tau_nodes"], arrays["basis_omega_nodes"]
    key = (tau_nodes.tobytes(), nodes.tobytes(), times.tobytes(), omega.tobytes(), beta)
    evaluated = 0
    if key not in cache:
        core = ctx.matrix(
            [[centered_kernel(ctx, t, w, beta) for w in nodes] for t in tau_nodes]
        )
        determinant = ctx.det(core)
        require(determinant != 0, "independent core determinant unresolved")
        truth = []
        for time in times:
            row = []
            for frequency in omega:
                border = ctx.matrix(
                    [
                        [centered_kernel(ctx, t, w, beta) for w in [*nodes, frequency]]
                        for t in [*tau_nodes, time]
                    ]
                )
                row.append(ctx.det(border) / determinant)
                evaluated += 1
            truth.append(row)
        cache[key] = truth
    truth = cache[key]
    grid = record["grid_certificate"]
    maximum = max(abs(value) for row in truth for value in row)
    discrepancy = max(
        abs(truth[i][j] - ctx.mpf(float(rounded[i, j])))
        for i in range(5)
        for j in range(5)
    )
    i, j = grid["candidate_index"]
    candidate = abs(truth[i][j])
    quality = candidate / maximum if maximum else ctx.one
    require(
        discrepancy <= ctx.mpf(grid["discrepancy_upper"]),
        "independent true-grid discrepancy exceeds enclosure",
    )
    require(
        candidate >= ctx.mpf(grid["true_candidate_lower"]),
        "independent candidate below lower bound",
    )
    require(
        maximum <= ctx.mpf(grid["true_global_upper"]),
        "independent grid exceeds global upper bound",
    )
    require(quality >= ctx.mpf(grid["eta"]), "independent grid quality below eta")
    return {
        "candidate": grid["candidate_index"],
        "reported_eta": grid["eta"],
        "independent_quality": ctx.nstr(quality, 40),
        "independent_maximum": ctx.nstr(maximum, 40),
        "independent_discrepancy": ctx.nstr(discrepancy, 40),
    }, evaluated


def _scalar(ctx: Any, value: Any) -> Any:
    return ctx.mpc(float(np.real(value)), float(np.imag(value)))


def audit_selected_nodes(
    record: dict[str, Any], arrays: dict[str, np.ndarray], ctx: Any
) -> dict[str, Any]:
    """Exact coefficient oracle; independent centered-core LU query checks."""
    require(
        record["training"] == "selected_nodes",
        "selected-node oracle requires selected nodes",
    )
    certificate = record["fit_certificate"]
    require(
        certificate["certified"], "selected-node oracle requires completed enclosures"
    )
    coef, exact_coef = arrays["coefficients"], arrays["observations"]
    bounds = certificate["coefficients"]
    errors = _numeric_bound(
        bounds["actual_error_bound"], coef.shape, "coefficient error"
    )
    max_coefficient = ctx.zero
    for index in np.ndindex(coef.shape):
        ideal = _scalar(ctx, exact_coef[index])
        error = abs(ideal - _scalar(ctx, coef[index]))
        require(
            error <= ctx.mpf(float(errors[index])),
            "exact selected-node coefficient exceeds bound",
        )
        for part in ("real", "imag"):
            lower = np.asarray(bounds["ideal"][part + "_lower"])
            upper = np.asarray(bounds["ideal"][part + "_upper"])
            require(
                ctx.mpf(float(lower[index]))
                <= getattr(ideal, part)
                <= ctx.mpf(float(upper[index])),
                "exact coefficient outside enclosure",
            )
        max_coefficient = max(max_coefficient, error)
    beta, nodes, tau_nodes = (
        record["beta"],
        arrays["basis_omega_nodes"],
        arrays["basis_tau_nodes"],
    )
    core = ctx.matrix(
        [[centered_kernel(ctx, t, w, beta) for w in nodes] for t in tau_nodes]
    )
    query_records = {}
    for name, points, actual_name in (
        ("time", arrays["query_tau"], "actual_time_output"),
        ("matsubara", arrays["matsubara_indices"], "actual_frequency_output"),
    ):
        query = certificate[name]
        maximum = ctx.zero
        for row, point in enumerate(points):
            rhs = ctx.matrix(
                [
                    centered_kernel(ctx, point, w, beta)
                    if name == "time"
                    else -ctx.mpf(beta)
                    / (ctx.j * ctx.pi * (2 * int(point) + 1) - ctx.mpf(float(w)))
                    for w in nodes
                ]
            )
            design = ctx.lu_solve(core.T, rhs)
            for batch in np.ndindex(2, 2):
                index = (row, *batch)
                ideal = ctx.fsum(
                    design[k] * _scalar(ctx, exact_coef[(k, *batch)])
                    for k in range(len(nodes))
                )
                error = abs(ideal - _scalar(ctx, arrays[actual_name][index]))
                for field in ("direct_error_bound", "total_error_bound"):
                    require(
                        error <= ctx.mpf(float(np.asarray(query[field])[index])),
                        "independent selected-node query exceeds budget",
                    )
                for part in ("real", "imag"):
                    lower = np.asarray(query["ideal_prediction"][part + "_lower"])
                    upper = np.asarray(query["ideal_prediction"][part + "_upper"])
                    require(
                        ctx.mpf(float(lower[index]))
                        <= getattr(ideal, part)
                        <= ctx.mpf(float(upper[index])),
                        "independent selected-node query outside enclosure",
                    )
                maximum = max(maximum, error)
        query_records[name] = {
            "checked_values": 20,
            "independent_max_error": ctx.nstr(maximum, 40),
        }
    return {
        "coefficient_components": int(coef.size),
        "max_coefficient_error": ctx.nstr(max_coefficient, 40),
        "query_values_checked": 40,
        "queries": query_records,
    }


def audit_physical_errors(
    record: dict[str, Any], arrays: dict[str, np.ndarray], ctx: Any
) -> dict[str, str]:
    """Independent finite-pole truth; never infer a physical accuracy pass."""
    beta = record["beta"]
    energies = record["wmax"] * np.array([-0.31, 0.23, 0.71])
    weights = [0.25, 0.5, 0.25]
    amplitude = np.array([[1, 0.25 + 0.5j], [0.25 - 0.5j, 0.75]])
    result = {}
    for name, points, actual_name in (
        ("time", arrays["query_tau"], "actual_time_output"),
        ("matsubara", arrays["matsubara_indices"], "actual_frequency_output"),
    ):
        error = ctx.zero
        for row, point in enumerate(points):
            if name == "time":
                truth = -ctx.fsum(
                    ctx.mpf(w)
                    * ctx.exp(-ctx.mpf(float(point)) * ctx.mpf(float(e)))
                    / (1 + ctx.exp(-ctx.mpf(beta) * ctx.mpf(float(e))))
                    for w, e in zip(weights, energies, strict=True)
                )
            else:
                truth = ctx.fsum(
                    ctx.mpf(w)
                    / (
                        ctx.j * ctx.pi * (2 * int(point) + 1) / ctx.mpf(beta)
                        - ctx.mpf(float(e))
                    )
                    for w, e in zip(weights, energies, strict=True)
                )
            for batch in np.ndindex(2, 2):
                error = max(
                    error,
                    abs(
                        _scalar(ctx, arrays[actual_name][(row, *batch)])
                        - truth * _scalar(ctx, amplitude[batch])
                    ),
                )
        reported = (
            record.get("queries", {}).get(name, {}).get("sampled_physical_pole_error")
        )
        if reported is not None:
            require(
                abs(error - ctx.mpf(reported))
                <= max(ctx.mpf("1e-45"), abs(error) * ctx.mpf("1e-33")),
                "reported physical error mismatch",
            )
        result[name] = ctx.nstr(error, 40)
    return result


def audit_artifact(
    source: Path,
    root: Path,
    generation_commit: str,
    expected_sha256: str,
    expected_payload_sha256: str,
) -> dict[str, Any]:
    artifact, container_hash, payload_hash = read_artifact(
        source, expected_sha256, expected_payload_sha256
    )
    require(
        type(artifact.get("schema_version")) is int and artifact["schema_version"] == 1,
        "unsupported artifact schema",
    )
    blob_count = validate_provenance(artifact, root, generation_commit)
    config = artifact.get("config")
    require(
        isinstance(config, dict)
        and type(config.get("schema_version")) is int
        and config["schema_version"] == 1,
        "unsupported config schema",
    )
    records, cases = artifact.get("records"), config.get("cases")
    require(
        isinstance(records, list)
        and isinstance(cases, list)
        and 0 < len(records) == len(cases),
        "case count mismatch",
    )
    decoded = [
        validate_record(record, case)
        for record, case in zip(records, cases, strict=True)
    ]
    ctx = mp.mp.clone()
    ctx.prec = AUDIT_BITS
    cache: dict[Any, Any] = {}
    grid_records, selected, physical = [], [], []
    determinant_count = 0
    for index, (record, (arrays, grid, _)) in enumerate(
        zip(records, decoded, strict=True)
    ):
        math_record, evaluated = audit_grid_math(record, arrays, grid, ctx, cache)
        grid_records.append({"record": index, **math_record})
        determinant_count += evaluated
        if (
            record["training"] == "selected_nodes"
            and record["fit_certificate"]["certified"]
        ):
            selected.append(
                {"record": index, **audit_selected_nodes(record, arrays, ctx)}
            )
        physical.append({"record": index, **audit_physical_errors(record, arrays, ctx)})
    return {
        "status": "passed",
        "generation_commit": generation_commit,
        "generation_tree": artifact["source"]["git_tree"],
        "artifact_container_sha256": container_hash,
        "artifact_payload_sha256": payload_hash,
        "audit_script_sha256": sha256(Path(__file__).read_bytes()),
        "source_blob_count": blob_count,
        "packed_array_count": len(INPUT_NAMES) * len(records),
        "reconstructed_core_hashes": sum(item[2] for item in decoded),
        "config_git_blob_verified": True,
        "generation_versions": artifact["versions"],
        "audit_versions": {
            name: importlib.metadata.version(name) for name in ("mpmath", "numpy")
        },
        "audit_precision_bits": ctx.prec,
        "unique_grid_computations": len(cache),
        "bordered_determinants_evaluated": determinant_count,
        "grid_records": grid_records,
        "selected_node_records": selected,
        "physical_pole_values_checked": len(records) * 40,
        "physical_errors": physical,
        "enclosed_fit_records": sum(
            record["fit_certificate"]["certified"] for record in records
        ),
        "inconclusive_fit_records": [
            {"record": i, "reason": r["fit_certificate"]["termination_reason"]}
            for i, r in enumerate(records)
            if not r["fit_certificate"]["certified"]
        ],
        "limitations": [
            "No kernelgecp module, interval engine or benchmark helper was imported.",
            "Independent high-precision computations are consistency evidence, "
            "not interval proofs or Lean verification of runtime.",
            "True grids use centered-kernel bordered determinants; selected-node "
            "coefficients use exact reproduction and queries use direct LU.",
            "No independent LS recomputation of other enclosed fits is claimed; "
            "their recorded QR checks remain generation evidence.",
            "Inconclusive interval pivot reports are retained, not reclassified "
            "as singularity proofs or successful certificates.",
            "Physical pole errors are sampled observations; small arithmetic "
            "bounds do not imply physical accuracy.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--project-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--generation-commit", required=True)
    parser.add_argument(
        "--expected-sha256", "--expected-container-sha256", required=True
    )
    parser.add_argument("--expected-payload-sha256", required=True)
    args = parser.parse_args()
    try:
        result = audit_artifact(
            args.source,
            args.project_root.resolve(),
            args.generation_commit,
            args.expected_sha256,
            args.expected_payload_sha256,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    except (AuditError, KeyError, TypeError, IndexError, OverflowError, OSError) as exc:
        print(f"Precision bridge audit failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(f"Independent precision bridge audit passed: {args.output}")


if __name__ == "__main__":
    main()
