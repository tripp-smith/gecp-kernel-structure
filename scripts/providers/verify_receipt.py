"""Fail closed on missing, failed, stale or wrong-version provider smoke receipts."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

EXPECTED = {
    "triqs": {"triqs": "4.0.2", "cthyb": "4.0.1", "maxent": "4.0.1"},
    "green": {"green": "1.0.0", "pyscf": "2.14.0"},
    "epw": {"qe": "7.5", "epw": "6.0"},
    "gp": {
        "numpy": "2.2.6",
        "scipy": "1.18.1",
        "mpmath": "1.3.0",
        "torch": "2.14.1",
        "gpytorch": "1.15.2",
        "linear-operator": "0.6.1",
    },
}
REQUIRED_CHECKS = {
    "gp": {"cpu_float64"},
    "triqs": {
        "source_revision",
        "pinned_versions",
        "one_mpi_rank",
        "fermionic_frequency_labels",
        "triqs_resolvent",
        "cthyb_sampling",
        "maxent_minimization",
        "single_thread",
    },
    "green": {
        "native_molecular_input_generation",
        "native_gf2_self_consistency",
        "stored_final_energy_and_green_deltas",
        "nonzero_interacting_selfenergy",
        "endpoint_normalization_and_particle_count",
    },
    "epw": {
        "source_lock_matches_image",
        "native_qe_7_5",
        "native_epw_6_0",
        "scf_converged",
        "phonon_calculation_completed",
        "official_inputs_unmodified",
        "native_eliashberg_converged",
        "finite_native_gap_mesh",
        "positive_mgb2_gap_at_15_kelvin",
    },
}


class ProviderReceiptError(ValueError):
    """The supplied receipt does not establish the promised execution gate."""


def _number(record: Any, name: str) -> float:
    value = record.get(name) if isinstance(record, dict) else None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProviderReceiptError(f"missing numeric execution evidence: {name}")
    try:
        result = float(value)
    except OverflowError as error:
        raise ProviderReceiptError(f"overflowed execution evidence: {name}") from error
    if not math.isfinite(result):
        raise ProviderReceiptError(f"nonfinite execution evidence: {name}")
    return result


def _evidence(data: dict[str, Any], provider: str) -> None:
    metrics = data.get("metrics")
    if provider == "gp":
        for name in ("solve_max_error", "logdet_error", "log_probability_error"):
            if not 0 <= _number(metrics, name) <= 1e-10:
                raise ProviderReceiptError(f"failed GP accuracy evidence: {name}")
    elif provider == "triqs":
        if (
            not 0 <= _number(metrics, "triqs_resolvent_max_error") < 1e-12
            or _number(metrics, "cthyb_cycles") != 200
            or not 0 < _number(metrics, "cthyb_average_sign_real") <= 1 + 1e-12
            or _number(metrics, "maxent_alpha_count") != 3
            or _number(metrics, "maxent_min_chi2") < 0
        ):
            raise ProviderReceiptError("failed TRIQS/CT-HYB/MaxEnt numerical evidence")
    elif provider == "green":
        values = data.get("measurements")
        if (
            not 2 <= _number(values, "iteration") <= 100
            or not 0 <= _number(values, "energy_change_sum") < 1e-10
            or not 0 <= _number(values, "green_change_max") < 1e-6
            or abs(_number(values, "correlation_energy_hartree")) <= 1e-8
            or _number(values, "dynamic_selfenergy_max") <= 1e-8
            or not 0 <= _number(values, "endpoint_normalization_error") <= 1e-5
            or abs(_number(values, "electron_count") - 2) > 1e-5
        ):
            raise ProviderReceiptError(
                "failed native Green convergence/numerical evidence"
            )
        artifacts = data.get("artifacts")
        if (
            not isinstance(artifacts, dict)
            or not {
                "green/input.log",
                "green/gf2.log",
                "green/input.h5",
                "green/sim.h5",
                "green/run-manifest.json",
                "green/build-provenance.json",
            }.issubset(artifacts)
            or any(
                not isinstance(value, str)
                or re.fullmatch(r"[0-9a-f]{64}", value) is None
                for value in artifacts.values()
            )
        ):
            raise ProviderReceiptError("missing native Green evidence identities")
    else:
        values = data.get("calculation")
        if (
            not 0.1 < _number(values, "electron_phonon_coupling") < 3
            or _number(values, "native_frequency_count") < 8
            or _number(values, "lowest_frequency_gap_ev") <= 0
            or _number(values, "temperature_kelvin") != 15
            or _number(values, "phonon_irreducible_qpoints") < 1
        ):
            raise ProviderReceiptError("failed native EPW numerical evidence")
        executions = data.get("native_executions")
        if not isinstance(executions, list) or len(executions) != 5:
            raise ProviderReceiptError("missing actual native EPW pipeline executions")
        for execution in executions:
            if (
                _number(execution, "returncode") != 0
                or _number(execution, "mpi_ranks") != 1
                or _number(execution, "omp_threads") != 1
            ):
                raise ProviderReceiptError("failed or unbounded native EPW execution")
            for name in ("input_sha256", "log_sha256", "executable_sha256"):
                value = execution.get(name)
                if (
                    not isinstance(value, str)
                    or re.fullmatch(r"[0-9a-f]{64}", value) is None
                ):
                    raise ProviderReceiptError(
                        f"missing native EPW evidence hash: {name}"
                    )


def validate_receipt(data: Any, provider: str, source_revision: str) -> dict[str, Any]:
    """Validate identity and successful checks, never infer execution from a lock."""
    if (
        provider not in EXPECTED
        or re.fullmatch(r"[0-9a-f]{40}", source_revision) is None
    ):
        raise ProviderReceiptError("invalid expected provider or source revision")
    if not isinstance(data, dict):
        raise ProviderReceiptError("provider receipt must be a JSON object")
    if (
        data.get("schema") != "kernelgecp.provider-smoke.v1"
        or data.get("provider") != provider
        or data.get("source_revision") != source_revision
        or data.get("status") != "passed"
    ):
        raise ProviderReceiptError("wrong, failed, missing or stale provider identity")
    checks = data.get("checks")
    if (
        not isinstance(checks, dict)
        or not ({"imports", "numerical"} | REQUIRED_CHECKS[provider]).issubset(checks)
        or any(not isinstance(key, str) or not key for key in checks)
        or any(value is not True for value in checks.values())
    ):
        raise ProviderReceiptError(
            "all required imports and numerical checks must pass"
        )
    versions = data.get("versions")
    if not isinstance(versions, dict):
        raise ProviderReceiptError("missing actual provider versions")
    python = versions.get("python")
    if not isinstance(python, str) or re.fullmatch(r"3\.12\.[0-9]+", python) is None:
        raise ProviderReceiptError(
            "initial provider runtimes require actual Python 3.12"
        )
    for package, expected in EXPECTED[provider].items():
        actual = versions.get(package)
        supported = {expected, expected + "+cpu"} if package == "torch" else {expected}
        if not isinstance(actual, str) or actual not in supported:
            raise ProviderReceiptError(
                f"wrong or missing {package} version: {actual!r}"
            )
    _evidence(data, provider)
    return data


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ProviderReceiptError(f"duplicate receipt key: {key}")
        result[key] = value
    return result


def _constant(value: str) -> Any:
    raise ProviderReceiptError(f"nonfinite JSON constant: {value}")


def read_receipt(path: Path) -> Any:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 4 * 1024 * 1024:
        raise ProviderReceiptError("receipt must be a regular file of at most 4 MiB")
    return json.loads(
        path.read_text(), object_pairs_hook=_object, parse_constant=_constant
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--provider", choices=sorted(EXPECTED), required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    try:
        validate_receipt(
            read_receipt(args.receipt), args.provider, args.source_revision
        )
    except (OSError, ValueError) as error:
        parser.exit(1, f"provider receipt rejected: {error}\n")
    print(
        f"verified {args.provider} numerical runtime receipt at {args.source_revision}"
    )


if __name__ == "__main__":
    main()
