"""Independently audit saved ED equations, residuals and finite-model tails.

This auditor does not import the campaign runner, its adapters or its truth
functions. Saved native reference data remain numerical IPT evidence; neither
those data nor the finite-pole remainder test certifies interacting-model tails.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import tomllib
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("kernelgecp", "pydlr", "sparse_ir")
ERROR_FIELDS = ("g_iw", "sigma_iw", "g_tau", "sigma_tau")
SOURCE_FILES = (
    "kernelgecp/dmft.py",
    "kernelgecp/lehmann.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "kernelgecp/integrations/triqs_dmft.py",
    "experiments/interacting_dmft.py",
    "experiments/dlr_adapters.py",
    "environments/locks/dlr-physics-pip.txt",
)


class AuditError(ValueError):
    """A missing or inconsistent artifact cannot establish a passing audit."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _number(value: Any, name: str, *, nonnegative: bool = True) -> float:
    _require(
        type(value) in (int, float) and math.isfinite(value),
        f"{name} must be a finite real scalar",
    )
    result = float(value)
    _require(not nonnegative or result >= 0, f"{name} must be nonnegative")
    return result


def _integer(value: Any, name: str, minimum: int = 0) -> int:
    _require(type(value) is int and value >= minimum, f"invalid integer {name}")
    return value


def _json(path: Path) -> dict[str, Any]:
    def invalid_constant(value: str) -> None:
        raise AuditError(f"nonfinite JSON constant {value}")

    value = json.loads(path.read_text(), parse_constant=invalid_constant)
    _require(isinstance(value, dict), f"{path.name} must contain an object")
    return value


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _blob(repository: Path, revision: str, name: str) -> bytes:
    path = Path(name)
    _require(
        not path.is_absolute() and ".." not in path.parts and ":" not in name,
        "unsafe source-relative path",
    )
    return subprocess.check_output(
        ["git", "show", f"{revision}:{name}"], cwd=repository, stderr=subprocess.PIPE
    )


def _source_bindings(report: dict[str, Any], repository: Path) -> None:
    revision = report.get("source_revision")
    _require(
        isinstance(revision, str)
        and re.fullmatch("[0-9a-f]{40}", revision) is not None,
        "source_revision must be a full Git object identity",
    )
    hashes = report.get("source_hashes")
    _require(isinstance(hashes, dict), "source_hashes must be present")
    _require(set(SOURCE_FILES) <= set(hashes), "source bindings are incomplete")
    for name, expected in hashes.items():
        _require(
            isinstance(expected, str)
            and re.fullmatch("[0-9a-f]{64}", expected) is not None,
            f"invalid source hash: {name}",
        )
        actual = hashlib.sha256(_blob(repository, revision, name)).hexdigest()
        _require(actual == expected, f"Git source blob mismatch: {name}")
    config_blob = _blob(repository, revision, "experiments/interacting_dmft.toml")
    _require(
        hashlib.sha256(config_blob).hexdigest() == report.get("config_sha256"),
        "configuration bytes do not match the committed generation configuration",
    )
    _require(
        tomllib.loads(config_blob.decode()) == report.get("config"),
        "reported configuration differs from the committed configuration",
    )


def _config(data: Any) -> dict[str, Any]:
    _require(isinstance(data, dict), "configuration is missing")
    # This auditor closes the declared ED campaign, not an arbitrary easier run.
    expected = {
        "beta": 10.0,
        "hopping": 1.0,
        "interactions": [0.0, 1.0, 2.0],
        "n_iw": 512,
        "n_tau": 6145,
        "reference_n_iw": 1024,
        "reference_n_tau": 12289,
        "mixing": 0.5,
        "max_iterations": 200,
        "iteration_tolerance": 1e-8,
        "target": 1e-6,
        "calibration_target": 1e-7,
        "reference_refinement_target": 1e-7,
        "fit_rcond": 1e-13,
        "selfenergy_cutoff_multiplier": 3.0,
        "repeats": 3,
        "methods": list(METHODS),
    }
    for name, value in expected.items():
        _require(data.get(name) == value, f"ED contract configuration mismatch: {name}")
        if type(value) in (int, float):
            _number(data.get(name), name)
    _require(
        all(type(value) in (int, float) for value in data["interactions"]),
        "interaction values must be real numbers, not booleans",
    )
    for name in ("calibration_indices", "validation_indices"):
        values = data.get(name)
        _require(isinstance(values, list) and bool(values), f"missing {name}")
        _require(
            all(type(index) is int for index in values)
            and len(values) == len(set(values)),
            f"invalid {name}",
        )
    _require(
        not set(data["calibration_indices"]) & set(data["validation_indices"]),
        "calibration and final frequencies overlap",
    )
    _require(
        all(-512 <= n < 512 and not -18 <= n < 18 for n in data["calibration_indices"]),
        "calibration frequencies must leave essential low modes observed",
    )
    _require(
        all(512 <= abs(n + 0.5) < 1024 for n in data["validation_indices"]),
        "final frequencies must be outside the primary mesh and inside the reference",
    )
    return data


def _selection_records(report: dict[str, Any], data: dict[str, Any]) -> None:
    rows = report.get("selection")
    refinements = report.get("cutoff_refinement")
    calibration = report.get("calibration")
    _require(isinstance(rows, list) and len(rows) == 9, "missing frozen selections")
    _require(
        isinstance(refinements, list) and len(refinements) == 9,
        "missing doubled-cutoff controls",
    )
    _require(
        isinstance(calibration, list) and bool(calibration),
        "missing calibration records",
    )
    selected = {}
    for row in rows:
        identity = row.get("method"), row.get("case")
        _require(
            identity[0] in METHODS
            and identity[1] in ("U0_0", "U1_0", "U2_0")
            and identity not in selected,
            "duplicate or invalid frozen selection",
        )
        _require(
            row.get("cutoff") in data["green_cutoffs"]
            and row.get("eps") in data["native_eps_candidates"],
            "selection is outside the frozen calibration ladder",
        )
        selected[identity] = row
        matched = [
            record
            for record in calibration
            if record.get("method") == identity[0]
            and "U" + str(record.get("interaction")).replace(".", "_") == identity[1]
            and record.get("green_cutoff") == row["cutoff"]
            and record.get("native_eps") == row["eps"]
            and record.get("converged") is True
        ]
        _require(
            len(matched) == 1, "selection lacks exactly one passing calibration record"
        )
        _errors(
            matched[0].get("errors"),
            data["calibration_target"],
            "selected calibration errors",
        )
    seen = set()
    for row in refinements:
        identity = (
            row.get("method"),
            "U" + str(row.get("interaction")).replace(".", "_"),
        )
        _require(
            identity in selected and identity not in seen, "duplicate cutoff refinement"
        )
        seen.add(identity)
        choice = selected[identity]
        _require(
            row.get("green_cutoff") == 2 * choice["cutoff"]
            and row.get("sigma_cutoff") == 6 * choice["cutoff"]
            and row.get("native_eps") == choice["eps"]
            and row.get("converged") is True
            and row.get("status") == "completed",
            "cutoff control changed more than the declared support bounds",
        )
        _errors(row.get("errors"), data["calibration_target"], "doubled-cutoff errors")


def _errors(value: Any, budget: float, name: str) -> dict[str, float]:
    _require(isinstance(value, dict), f"missing {name}")
    _require(set(value) == set(ERROR_FIELDS), f"incomplete {name}")
    result = {
        field: _number(value[field], name + "." + field) for field in ERROR_FIELDS
    }
    _require(max(result.values()) <= budget, f"{name} exceeds its accuracy budget")
    return result


def _history(value: Any, data: dict[str, Any], name: str) -> list[dict[str, Any]]:
    _require(
        isinstance(value, list) and 2 <= len(value) <= data["max_iterations"],
        f"{name} requires at least two completed updates",
    )
    for i, row in enumerate(value, 1):
        _require(
            isinstance(row, dict)
            and type(row.get("iteration")) is int
            and row["iteration"] == i,
            f"{name} iteration sequence is invalid",
        )
        _number(row.get("unmixed_defect"), name + ".unmixed_defect")
        _number(row.get("mixed_step"), name + ".mixed_step")
    _require(
        all(row["unmixed_defect"] < data["iteration_tolerance"] for row in value[-2:]),
        f"{name} did not pass the strict unmixed defect twice",
    )
    return value


def _artifacts(directory: Path, record: dict[str, Any], required: set[str]) -> None:
    hashes = record.get("artifacts")
    _require(
        isinstance(hashes, dict) and required <= set(hashes),
        "required artifact bindings are missing",
    )
    for name, expected in hashes.items():
        _require(Path(name).name == name, "artifact paths must be local filenames")
        _require(
            isinstance(expected, str) and _sha(directory / name) == expected,
            f"artifact hash mismatch: {directory.name}/{name}",
        )


def _arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        result = {name: np.array(archive[name]) for name in archive.files}
    _require(bool(result), "empty numerical artifact")
    for name, value in result.items():
        _require(
            value.dtype.kind in "iufc" and value.size > 0,
            f"invalid numeric array {name}",
        )
        _require(bool(np.all(np.isfinite(value))), f"nonfinite numeric array {name}")
    return result


def _same(actual: Any, expected: Any, name: str, *, atol: float = 2e-11) -> float:
    a, b = np.asarray(actual), np.asarray(expected)
    _require(a.shape == b.shape, f"shape mismatch: {name}")
    error = float(np.max(np.abs(a - b), initial=0))
    scale = max(1.0, float(np.max(np.abs(b), initial=0)))
    _require(
        math.isfinite(error) and error <= atol * scale,
        f"numerical identity failed: {name} ({error})",
    )
    return error


def _physics(
    arrays: dict[str, np.ndarray], record: dict[str, Any], budget: float
) -> None:
    physical = record.get("physics")
    _require(isinstance(physical, dict), "missing physical diagnostics")
    g = arrays.get("represented_g_iw", arrays["g_iw"])
    sigma = arrays["sigma_iw"]
    positive = arrays["indices"] >= 0
    endpoints = arrays["g_tau"][[0, -1]]
    diagnostics = {
        "density_error_from_half": float(abs(-endpoints[-1] - 0.5)),
        "endpoint_mass_error": float(abs(-sum(endpoints) - 1)),
        "positive_imag_green_max": float(np.max(g[positive].imag)),
        "positive_imag_sigma_max": float(np.max(sigma[positive].imag)),
        "particle_hole_green_real_max": float(np.max(np.abs(g.real))),
        "particle_hole_sigma_real_max": float(np.max(np.abs(sigma.real))),
    }
    for name, actual in diagnostics.items():
        recorded = _number(physical.get(name), "physics." + name, nonnegative=False)
        _same(recorded, actual, "reported physics." + name)
        _require(actual <= budget, f"physical diagnostic exceeds budget: {name}")
    _require(
        physical.get("certification_status") == "not_certified_no_tail_remainder_bound",
        "interacting tail diagnostics must retain their uncertified scope",
    )
    tails = physical.get("tail_fits")
    _require(
        isinstance(tails, dict) and set(tails) == {"green", "sigma"},
        "missing unconstrained tail diagnostics",
    )
    for tail in tails.values():
        windows = tail.get("window_estimates")
        _require(
            isinstance(windows, list) and len(windows) == 2,
            "two independent tail windows must be recorded",
        )
        for window in windows:
            _number(window.get("residual"), "tail residual")
            _number(window.get("condition"), "tail condition")
            _integer(window.get("effective_rank"), "tail rank", 1)


def _numerics(
    arrays: dict[str, np.ndarray],
    record: dict[str, Any],
    data: dict[str, Any],
    *,
    compressed: bool,
) -> dict[str, float]:
    frequency_fields = {"bath_iw", "g_iw", "sigma_iw", "weiss_iw"}
    time_fields = {"weiss_tau", "sigma_tau", "g_tau"}
    required = frequency_fields | time_fields | {"indices", "tau"}
    if compressed:
        required |= {
            "represented_g_iw",
            "sigma_tau_product",
            "weiss_reflected_tau",
            "frequency_fit_indices",
            "g_coefficients",
            "sigma_coefficients",
        }
        frequency_fields.add("represented_g_iw")
        time_fields |= {"sigma_tau_product", "weiss_reflected_tau"}
    _require(required <= set(arrays), "missing saved numerical fields")
    indices, tau = arrays["indices"], arrays["tau"]
    _require(
        indices.dtype.kind in "iu"
        and np.array_equal(indices, np.arange(-data["n_iw"], data["n_iw"])),
        "wrong primary frequency mesh",
    )
    _same(tau, np.linspace(0, data["beta"], data["n_tau"]), "uniform primary time mesh")
    for name in frequency_fields:
        _require(arrays[name].shape == indices.shape, f"wrong shape: {name}")
    for name in time_fields:
        _require(arrays[name].shape == tau.shape, f"wrong shape: {name}")
    interaction = _number(record.get("interaction"), "interaction")
    z = 1j * (2 * indices.astype(float) + 1) * (np.pi / data["beta"])
    # Residual form avoids reproducing the solver's scalar division.
    metrics = {
        "weiss_equation": _same(
            (z - arrays["bath_iw"]) * arrays["weiss_iw"],
            np.ones(indices.size),
            "Weiss equation",
        ),
        "dyson_equation": _same(
            (z - arrays["bath_iw"] - arrays["sigma_iw"]) * arrays["g_iw"],
            np.ones(indices.size),
            "Dyson equation",
        ),
    }
    observed = np.arange(indices.size)
    product = interaction**2 * arrays["weiss_tau"] ** 2 * arrays["weiss_tau"][::-1]
    if compressed:
        fit_indices = arrays["frequency_fit_indices"]
        expected_fit = np.setdiff1d(
            indices, data["calibration_indices"] + data["validation_indices"]
        )
        _require(
            fit_indices.dtype.kind in "iu"
            and np.array_equal(fit_indices, expected_fit),
            "saved observation subset differs from frozen permitted observations",
        )
        observed = fit_indices + data["n_iw"]
        _same(
            arrays["weiss_reflected_tau"],
            arrays["weiss_tau"][::-1],
            "reflected Weiss samples",
        )
        metrics["reflected_product"] = _same(
            arrays["sigma_tau_product"], product, "reflected IPT product"
        )
        for field, rank_field in (
            ("g_coefficients", "green_rank"),
            ("sigma_coefficients", "sigma_rank"),
        ):
            rank = _integer(record.get(rank_field), rank_field, 1)
            _require(arrays[field].shape == (rank,), f"wrong {field} shape")
        diagnostic = record.get("diagnostics")
        _require(isinstance(diagnostic, dict), "missing projection diagnostics")
        _same(
            diagnostic.get("sigma_projection_residual"),
            np.max(abs(arrays["sigma_tau"] - arrays["sigma_tau_product"])),
            "self-energy projection residual",
        )
        _same(
            diagnostic.get("green_projection_residual"),
            np.max(
                abs(arrays["represented_g_iw"][observed] - arrays["g_iw"][observed])
            ),
            "Green projection residual",
        )
        _require(
            diagnostic.get("certification_status") == "not_certified",
            "projected fits cannot claim certification",
        )
    else:
        metrics["reflected_product"] = _same(
            arrays["sigma_tau"], product, "native reflected IPT product"
        )
    defect = data["hopping"] ** 2 * arrays["g_iw"] - arrays["bath_iw"]
    observed_defect = float(np.max(abs(defect[observed])))
    history = _history(record.get("history"), data, "history")
    _same(
        history[-1]["unmixed_defect"],
        observed_defect,
        "final observed defect",
        atol=5e-13,
    )
    _require(
        observed_defect < data["iteration_tolerance"], "saved bath is not converged"
    )
    _same(
        history[-1]["mixed_step"],
        data["mixing"] * np.max(abs(defect)),
        "final proposed mixed step",
        atol=5e-13,
    )
    if compressed:
        _same(
            history[-1].get("full_grid_unmixed_defect"),
            np.max(abs(defect)),
            "reported full-grid defect",
            atol=5e-13,
        )
    _physics(arrays, record, data["target"])
    return metrics


def _validation_errors(
    arrays: dict[str, np.ndarray],
    truth: dict[str, np.ndarray],
    record: dict[str, Any],
    data: dict[str, Any],
) -> None:
    for name in ("indices", "tau"):
        _require("validation_" + name in arrays, "missing final query arrays")
        _same(arrays["validation_" + name], truth[name], "final query " + name, atol=0)
    errors = _errors(record.get("errors"), data["target"], "final errors")
    for field in ERROR_FIELDS:
        _require("validation_" + field in arrays, "missing final prediction: " + field)
        _require(
            arrays["validation_" + field].shape == truth[field].shape,
            "wrong final prediction shape: " + field,
        )
        actual = float(np.max(abs(arrays["validation_" + field] - truth[field])))
        _require(actual <= data["target"], "saved final prediction misses: " + field)
        _same(errors[field], actual, "recomputed final error " + field, atol=5e-14)


def _fingerprint(arrays: dict[str, np.ndarray]) -> str:
    digest = hashlib.sha256()
    for name, value in sorted(arrays.items()):
        digest.update(name.encode())
        digest.update(value.dtype.str.encode())
        digest.update(str(value.shape).encode())
        digest.update(np.ascontiguousarray(value).tobytes())
    return digest.hexdigest()


def finite_pole_audit() -> dict[str, Any]:
    """Exp quadrature and a bounded geometric tail for one finite fixture only."""
    from kernelgecp.lehmann import FermionicPoleModel, ipt_self_energy_poles

    beta, interaction = 10.0, 2.0
    energies, weights = [-1.3, 0.2, 2.0], [0.2, 0.5, 0.3]
    model = FermionicPoleModel(beta, energies, weights)
    sigma = ipt_self_energy_poles(model, interaction)
    modes = [-17, -4, 0, 5, 19]
    actual = sigma.evaluate_matsubara(modes)
    with mp.workdps(70):
        b, u = mp.mpf(10), mp.mpf(2)
        nodes = [mp.mpf(str(x)) for x in energies]
        masses = [mp.mpf(str(x)) for x in weights]

        def g(t: Any) -> Any:
            return -mp.fsum(
                w * mp.exp(-t * e) / (1 + mp.exp(-b * e))
                for e, w in zip(nodes, masses, strict=True)
            )

        errors = []
        for n, observed in zip(modes, actual, strict=True):
            z = mp.j * (2 * n + 1) * mp.pi / b
            truth = mp.quad(
                lambda t, z=z: mp.exp(z * t) * u**2 * g(t) ** 2 * g(b - t),
                [0, b / 4, b / 2, 3 * b / 4, b],
            )
            errors.append(float(abs(mp.mpc(complex(observed)) - truth)))
        density = mp.fsum(
            w / (1 + mp.exp(b * e)) for e, w in zip(nodes, masses, strict=True)
        )
        mass = u**2 * density * (1 - density)
        mass_error = float(abs(mp.mpf(float(sigma.mass)) - mass))
        _require(
            max(errors) < 3e-12 and mass_error < 3e-13,
            "independent finite-pole Fourier/sign/normalization check failed",
        )
        # The finite measure has nonnegative mass <= U^2/4 and |Omega|<=3 max|E|.
        # For p retained Laurent terms, summing the exact geometric remainder
        # gives mass/|z| * (W/|z|)^p / (1-W/|z|), when |z|>W.
        p, n = 6, 64
        z = mp.j * (2 * n + 1) * mp.pi / b
        support = 3 * max(abs(e) for e in nodes)
        moments = [mp.mpf(0) for _ in range(p)]
        for ea, wa in zip(nodes, masses, strict=True):
            for eb, wb in zip(nodes, masses, strict=True):
                for ec, wc in zip(nodes, masses, strict=True):
                    fa, fb, fc = [1 / (1 + mp.exp(b * e)) for e in (ea, eb, ec)]
                    omega = ea + eb - ec
                    weight = (
                        u**2
                        * wa
                        * wb
                        * wc
                        * ((1 - fa) * (1 - fb) * fc + fa * fb * (1 - fc))
                    )
                    for k in range(p):
                        moments[k] += weight * omega**k
        polynomial = mp.fsum(moment / z ** (k + 1) for k, moment in enumerate(moments))
        exact = mp.quad(
            lambda t: mp.exp(z * t) * u**2 * g(t) ** 2 * g(b - t),
            [b * k / 32 for k in range(33)],
        )
        remainder = abs(exact - polynomial)
        bound = mass / abs(z) * (support / abs(z)) ** p / (1 - support / abs(z))
        _require(remainder <= bound, "finite geometric tail remainder bound failed")
    return {
        "passed": True,
        "quadrature_decimal_digits": 70,
        "frequency_error_max": max(errors),
        "spectral_mass_error": mass_error,
        "tail_terms": p,
        "tail_index": n,
        "tail_remainder": float(remainder),
        "tail_remainder_bound": float(bound),
        "scope": "positive_finite_pole_fixture_only_not_interacting_DMFT_tail_bound",
    }


def audit_campaign(directory: Path, *, repository: Path = ROOT) -> dict[str, Any]:
    """Return a fail-closed audit receipt for the frozen ED campaign."""
    result: dict[str, Any] = {
        "schema": "kernelgecp.interacting-dmft-audit.v1",
        "passed": False,
        "checks": {},
        "failures": [],
        "certification_status": "not_certified",
        "scope": "saved_equations_sampled_native_reference_and_finite_fixture_bound",
        "auditor_source_sha256": _sha(Path(__file__)),
    }
    stage = "report"
    try:
        report = _json(directory / "report.json")
        report_hash = _sha(directory / "report.json")
        result["campaign_report_sha256"] = report_hash
        _require(
            report.get("schema") == "kernelgecp.interacting-dmft.v1",
            "wrong report schema",
        )
        _require(
            report.get("pilot") is False, "pilot data cannot close canonical validation"
        )
        _require(
            report.get("all_mandatory_passed") is True,
            "campaign did not pass its own gates",
        )
        data = _config(report.get("config"))
        stage = "source_bindings"
        _source_bindings(report, repository)
        result["source_revision"] = report["source_revision"]
        result["checks"][stage] = True
        stage = "selection_and_environment"
        _selection_records(report, data)
        versions = report.get("versions")
        _require(isinstance(versions, dict), "missing actual dependency versions")
        for name, expected in (
            ("pydlr", "1.0.1"),
            ("sparse-ir", "2.1.6"),
            ("pylibsparseir", "0.12.0"),
        ):
            _require(
                versions.get(name) == expected, f"provider version mismatch: {name}"
            )
        _require(
            report.get("triqs_version") == "4.0.2", "native provider version mismatch"
        )
        freeze = report.get("environment_freeze")
        _require(
            isinstance(freeze, list)
            and bool(freeze)
            and all(isinstance(line, str) for line in freeze),
            "missing actual environment freeze",
        )
        result["checks"][stage] = True
        stage = "reference"
        references = report.get("references")
        _require(
            isinstance(references, dict) and len(references) == 3,
            "missing native references",
        )
        truths, replay_counts = {}, {}
        for interaction in data["interactions"]:
            case = "U" + str(interaction).replace(".", "_")
            reference = references[case]
            _require(reference.get("status") == "completed", "native reference failed")
            _errors(
                reference.get("refinement_errors"),
                data["reference_refinement_target"],
                "reference refinement",
            )
            refinement = reference.get("off_grid_interpolation_refinement_errors")
            _require(
                isinstance(refinement, dict) and set(refinement) == {"green", "sigma"},
                "missing off-grid native interpolation refinement",
            )
            for error in refinement.values():
                _require(
                    _number(error, "interpolation refinement")
                    < data["reference_refinement_target"],
                    "native interpolation failed refinement",
                )
            history = _history(
                reference.get("primary_history"), data, "primary reference history"
            )
            _history(
                reference.get("refined_history"), data, "refined reference history"
            )
            folder = directory / case / "reference"
            _artifacts(folder, reference, {"trajectory.npz", "validation.npz"})
            trajectory = _arrays(folder / "trajectory.npz")
            _require(
                set(trajectory) == {"bath"}
                and trajectory["bath"].shape == (len(history), 2 * data["n_iw"]),
                "native replay trajectory is incomplete",
            )
            replay_counts[interaction] = len(history)
            truth = _arrays(folder / "validation.npz")
            _require(
                set(truth) == {"indices", "tau", *ERROR_FIELDS},
                "malformed native final reference",
            )
            _require(
                np.array_equal(truth["indices"], data["validation_indices"]),
                "wrong final frequency reference",
            )
            tau = data["beta"] * (np.polynomial.legendre.leggauss(64)[0] + 1) / 2
            # Different NumPy/LAPACK versions round these independently
            # generated coordinates by a few ulps. The actual query arrays
            # used for prediction remain byte-bound to the saved reference.
            _same(
                truth["tau"],
                tau,
                "frozen Gauss final queries",
                atol=8 * np.finfo(float).eps,
            )
            _require(
                not np.any(
                    np.isin(tau, np.linspace(0, data["beta"], data["reference_n_tau"]))
                ),
                "final time queries overlap training/refinement rows",
            )
            for name in ERROR_FIELDS:
                size = len(truth["indices"] if name.endswith("iw") else tau)
                _require(
                    truth[name].shape == (size,), "wrong native final reference shape"
                )
            truths[interaction] = truth
        result["checks"][stage] = True
        stage = "records"
        fingerprints: dict[tuple[str, float], list[str]] = {}
        counts = {"compressed": 0, "dense": 0}
        equation_errors: dict[str, float] = {}
        for compressed, field in ((True, "validation"), (False, "dense_validation")):
            rows = report.get(field)
            methods = METHODS if compressed else ("dense_triqs",)
            _require(
                isinstance(rows, list) and len(rows) == 3 * len(methods) * 3,
                f"missing or extra {field} records",
            )
            seen = set()
            for record in rows:
                method = record.get("method")
                interaction = _number(record.get("interaction"), "interaction")
                repeat = _integer(record.get("repeat"), "repeat")
                identity = method, interaction, repeat
                _require(
                    method in methods
                    and interaction in data["interactions"]
                    and repeat < 3
                    and identity not in seen,
                    "duplicate or unexpected numerical repetition",
                )
                seen.add(identity)
                _require(
                    record.get("status") == "completed"
                    and record.get("converged") is True
                    and record.get("accuracy_passed") is True
                    and record.get("stop_reason")
                    == (
                        "self_consistency_tolerance_reached"
                        if compressed
                        else "tolerance"
                    ),
                    "worker did not pass actual convergence and validation",
                )
                folder = directory / (
                    "validation" if compressed else "dense-validation"
                )
                folder /= (
                    f"repeat-{repeat}-{method}-U{interaction}"
                    if compressed
                    else f"repeat-{repeat}-U{interaction}"
                )
                _artifacts(folder, record, {"result.npz", "result.json"})
                saved = _json(folder / "result.json")
                for name in ("method", "interaction", "history", "errors", "physics"):
                    _require(
                        saved.get(name) == record.get(name),
                        f"report/worker mismatch: {name}",
                    )
                arrays = _arrays(folder / "result.npz")
                numerical = _numerics(arrays, record, data, compressed=compressed)
                for name, error in numerical.items():
                    equation_errors[name] = max(equation_errors.get(name, 0.0), error)
                _validation_errors(arrays, truths[interaction], record, data)
                if compressed:
                    replay = record.get("replay")
                    _require(
                        isinstance(replay, dict), "missing native trajectory replay"
                    )
                    _errors(replay.get("errors"), data["target"], "replay errors")
                    _require(
                        replay.get("updates_executed") == replay_counts[interaction],
                        "replay count differs from the complete native trajectory",
                    )
                    _require(
                        replay.get("scope")
                        == "common_native_bath_trajectory_not_independent_closed_loop",
                        "native replay scope is missing",
                    )
                fingerprints.setdefault((method, interaction), []).append(
                    _fingerprint(arrays)
                )
                counts["compressed" if compressed else "dense"] += 1
        result["checks"][stage] = True
        stage = "repetitions"
        _require(
            len(fingerprints) == 12
            and all(
                len(values) == 3 and len(set(values)) == 1
                for values in fingerprints.values()
            ),
            "saved numerical arrays differ across the three fresh repetitions",
        )
        result["checks"][stage] = True
        result["audited_records"] = counts
        result["saved_equation_residual_maxima"] = equation_errors
        result["repetition_array_sha256"] = {
            f"{method}-U{interaction}": values[0]
            for (method, interaction), values in sorted(fingerprints.items())
        }
        stage = "finite_pole_oracle"
        result[stage] = finite_pole_audit()
        result["checks"][stage] = True
        _require(
            _sha(directory / "report.json") == report_hash,
            "campaign report changed during its independent audit",
        )
        result["passed"] = True
    except (
        AuditError,
        AttributeError,
        IndexError,
        OverflowError,
        ValueError,
        TypeError,
        KeyError,
        OSError,
        subprocess.SubprocessError,
    ) as error:
        result["failures"].append(f"{stage}: {type(error).__name__}: {error}")
        result["checks"][stage] = False
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repository", type=Path, default=ROOT)
    args = parser.parse_args()
    receipt = audit_campaign(args.campaign, repository=args.repository)
    output = args.output or args.campaign / "audit.json"
    output.write_text(
        json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"passed": receipt["passed"], "output": str(output)}))
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
