"""Artifact auditor attacks; synthetic fixtures test equations, not IPT physics."""

from __future__ import annotations

import copy
import hashlib
import json
import shutil
import subprocess
import sys
import tomllib

import numpy as np
import pytest

from experiments import audit_interacting_dmft as audit


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False))


def bind(folder, record):
    record["artifacts"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in folder.iterdir()
        if path.is_file()
    }


def history():
    return [
        {
            "iteration": i,
            "unmixed_defect": 0.0,
            "mixed_step": 0.0,
            "full_grid_unmixed_defect": 0.0,
        }
        for i in (1, 2)
    ]


def physical(arrays):
    g = arrays.get("represented_g_iw", arrays["g_iw"])
    positive = arrays["indices"] >= 0
    return {
        "density_error_from_half": 0.0,
        "endpoint_mass_error": 0.0,
        "positive_imag_green_max": float(np.max(g[positive].imag)),
        "positive_imag_sigma_max": float(np.max(arrays["sigma_iw"][positive].imag)),
        "particle_hole_green_real_max": 0.0,
        "particle_hole_sigma_real_max": 0.0,
        "certification_status": "not_certified_no_tail_remainder_bound",
        "tail_fits": {
            name: {
                "window_estimates": [
                    {"residual": 0.0, "condition": 1.0, "effective_rank": 5}
                    for _ in (0, 1)
                ]
            }
            for name in ("green", "sigma")
        },
    }


def numerical_fixture(data, interaction):
    indices = np.arange(-data["n_iw"], data["n_iw"], dtype=np.int64)
    tau = np.linspace(0, data["beta"], data["n_tau"])
    z = 1j * (2 * indices + 1) * np.pi / data["beta"]
    sigma = interaction**2 / (4 * z)
    # This closed scalar Dyson fixture is deliberately not presented as a
    # complete physical IPT loop. In particular its constant Weiss time rows
    # are an algebraic product fixture, not a transform-accuracy oracle.
    a = (z - sigma).imag
    g = -2j * np.sign(a) / (abs(a) + np.sqrt(a**2 + 4))
    weiss_tau = np.full(tau.size, -0.5, complex)
    product = np.full(tau.size, -(interaction**2) / 8, complex)
    return {
        "indices": indices,
        "tau": tau,
        "frequency_fit_indices": np.setdiff1d(
            indices, data["calibration_indices"] + data["validation_indices"]
        ),
        "bath_iw": g.copy(),
        "g_iw": g,
        "represented_g_iw": g.copy(),
        "sigma_iw": sigma,
        "weiss_iw": 1 / (z - g),
        "weiss_tau": weiss_tau,
        "weiss_reflected_tau": weiss_tau.copy(),
        "sigma_tau_product": product,
        "sigma_tau": product.copy(),
        "g_tau": np.full(tau.size, -0.5, complex),
        "g_coefficients": np.array([1.0], complex),
        "sigma_coefficients": np.array([interaction**2 / 4], complex),
    }


def full_fixture(tmp_path):
    config_blob = (audit.ROOT / "experiments/interacting_dmft.toml").read_bytes()
    data = tomllib.loads(config_blob.decode())
    repository = tmp_path / "source"
    repository.mkdir()
    for name in audit.SOURCE_FILES:
        path = repository / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic source identity fixture for " + name)
    path = repository / "experiments/interacting_dmft.toml"
    path.write_bytes(config_blob)
    for command in (
        ["git", "init", "-q"],
        ["git", "add", "."],
        [
            "git",
            "-c",
            "user.name=Artifact Test",
            "-c",
            "user.email=artifact@example.invalid",
            "commit",
            "-qm",
            "Synthetic artifact provenance",
        ],
    ):
        subprocess.run(command, cwd=repository, check=True, capture_output=True)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True
    ).strip()
    campaign = tmp_path / "campaign"
    campaign.mkdir()
    report = {
        "schema": "kernelgecp.interacting-dmft.v1",
        "source_revision": revision,
        "source_hashes": {
            name: audit._sha(repository / name) for name in audit.SOURCE_FILES
        },
        "config_sha256": hashlib.sha256(config_blob).hexdigest(),
        "config": data,
        "pilot": False,
        "all_mandatory_passed": True,
        "references": {},
        "selection": [],
        "calibration": [],
        "cutoff_refinement": [],
        "validation": [],
        "dense_validation": [],
        "versions": {"pydlr": "1.0.1", "sparse-ir": "2.1.6", "pylibsparseir": "0.12.0"},
        "triqs_version": "4.0.2",
        "environment_freeze": ["pydlr==1.0.1"],
    }
    for interaction in data["interactions"]:
        case = "U" + str(interaction).replace(".", "_")
        arrays = numerical_fixture(data, interaction)
        folder = campaign / case / "reference"
        folder.mkdir(parents=True)
        final_indices = np.asarray(data["validation_indices"], dtype=np.int64)
        final_tau = data["beta"] * (np.polynomial.legendre.leggauss(64)[0] + 1) / 2
        z = 1j * (2 * final_indices + 1) * np.pi / data["beta"]
        truth = {
            "indices": final_indices,
            "tau": final_tau,
            "g_iw": 1 / z,
            "sigma_iw": interaction**2 / (4 * z),
            "g_tau": np.full(64, -0.5, complex),
            "sigma_tau": np.full(64, -(interaction**2) / 8, complex),
        }
        np.savez_compressed(folder / "validation.npz", **truth)
        np.savez_compressed(
            folder / "trajectory.npz", bath=np.stack([arrays["bath_iw"]] * 2)
        )
        reference = {
            "status": "completed",
            "refinement_errors": dict.fromkeys(audit.ERROR_FIELDS, 0.0),
            "off_grid_interpolation_refinement_errors": {"green": 0.0, "sigma": 0.0},
            "primary_history": history(),
            "refined_history": history(),
        }
        bind(folder, reference)
        report["references"][case] = reference
        arrays.update({"validation_" + name: values for name, values in truth.items()})
        for method in (*audit.METHODS, "dense_triqs"):
            compressed = method != "dense_triqs"
            if compressed:
                cutoff, eps = data["green_cutoffs"][0], data["native_eps_candidates"][0]
                report["selection"].append(
                    {"case": case, "method": method, "cutoff": cutoff, "eps": eps}
                )
                report["calibration"].append(
                    {
                        "method": method,
                        "interaction": interaction,
                        "green_cutoff": cutoff,
                        "native_eps": eps,
                        "converged": True,
                        "errors": dict.fromkeys(audit.ERROR_FIELDS, 0.0),
                    }
                )
                report["cutoff_refinement"].append(
                    {
                        "method": method,
                        "interaction": interaction,
                        "green_cutoff": 2 * cutoff,
                        "sigma_cutoff": 6 * cutoff,
                        "native_eps": eps,
                        "converged": True,
                        "status": "completed",
                        "errors": dict.fromkeys(audit.ERROR_FIELDS, 0.0),
                    }
                )
            for repeat in range(3):
                name = (
                    f"repeat-{repeat}-{method}-U{interaction}"
                    if compressed
                    else f"repeat-{repeat}-U{interaction}"
                )
                folder = (
                    campaign
                    / ("validation" if compressed else "dense-validation")
                    / name
                )
                folder.mkdir(parents=True)
                record = {
                    "method": method,
                    "interaction": interaction,
                    "repeat": repeat,
                    "status": "completed",
                    "converged": True,
                    "accuracy_passed": True,
                    "stop_reason": (
                        "self_consistency_tolerance_reached"
                        if compressed
                        else "tolerance"
                    ),
                    "history": history(),
                    "errors": dict.fromkeys(audit.ERROR_FIELDS, 0.0),
                    "physics": physical(arrays),
                    "green_rank": 1,
                    "sigma_rank": 1,
                    "diagnostics": {
                        "sigma_projection_residual": 0.0,
                        "green_projection_residual": 0.0,
                        "certification_status": "not_certified",
                    },
                    "replay": {
                        "updates_executed": 2,
                        "errors": dict.fromkeys(audit.ERROR_FIELDS, 0.0),
                        "scope": (
                            "common_native_bath_trajectory_not_independent_closed_loop"
                        ),
                    },
                }
                np.savez_compressed(folder / "result.npz", **arrays)
                write_json(folder / "result.json", record)
                bind(folder, record)
                report["validation" if compressed else "dense_validation"].append(
                    record
                )
    write_json(campaign / "report.json", report)
    return campaign, repository, report


@pytest.fixture(scope="module")
def complete_fixture(tmp_path_factory):
    return full_fixture(tmp_path_factory.mktemp("dmft-audit"))


def test_independent_finite_pole_quadrature_mass_and_geometric_remainder():
    receipt = audit.finite_pole_audit()
    assert receipt["passed"]
    assert receipt["frequency_error_max"] < 3e-12
    assert 0 < receipt["tail_remainder"] <= receipt["tail_remainder_bound"]
    assert "not_interacting" in receipt["scope"]


def test_complete_artifact_audit_independently_checks_all_saved_records(
    complete_fixture,
):
    campaign, repository, _ = complete_fixture
    receipt = audit.audit_campaign(campaign, repository=repository)
    assert receipt["passed"], receipt
    assert receipt["audited_records"] == {"compressed": 27, "dense": 9}
    assert all(receipt["checks"].values())


@pytest.mark.parametrize(
    "attack",
    [
        "pilot",
        "source",
        "config",
        "missing_repetition",
        "duplicate",
        "missing_selection",
        "empty_refinement",
        "missing_environment",
        "empty_history",
        "threshold",
        "replay_zero",
        "errors_missing",
        "density",
        "missing_artifact",
    ],
)
def test_campaign_rejects_missing_malformed_or_false_passing_metadata(
    complete_fixture, monkeypatch, attack
):
    campaign, repository, original = complete_fixture
    report = copy.deepcopy(original)
    row = report["validation"][0]
    if attack == "pilot":
        report["pilot"] = True
    elif attack == "source":
        report["source_hashes"][audit.SOURCE_FILES[0]] = "0" * 64
    elif attack == "config":
        report["config"]["target"] = 0.1
    elif attack == "missing_repetition":
        report["validation"].pop()
    elif attack == "duplicate":
        report["validation"][1] = copy.deepcopy(row)
    elif attack == "missing_selection":
        report["selection"].pop()
    elif attack == "empty_refinement":
        report["cutoff_refinement"] = []
    elif attack == "missing_environment":
        report.pop("environment_freeze")
    elif attack == "empty_history":
        row["history"] = []
    elif attack == "threshold":
        row["history"][-1]["unmixed_defect"] = report["config"]["iteration_tolerance"]
    elif attack == "replay_zero":
        row["replay"]["updates_executed"] = 0
    elif attack == "errors_missing":
        row["errors"] = {}
    elif attack == "density":
        row["physics"]["density_error_from_half"] = 1
    elif attack == "missing_artifact":
        row["artifacts"].pop("result.npz")
    actual_json = audit._json
    monkeypatch.setattr(
        audit,
        "_json",
        lambda path: report if path.name == "report.json" else actual_json(path),
    )
    receipt = audit.audit_campaign(campaign, repository=repository)
    assert not receipt["passed"]
    assert receipt["failures"]


def test_equation_audit_is_not_only_an_artifact_hash_check(complete_fixture):
    campaign, _, report = complete_fixture
    path = campaign / "validation/repeat-0-kernelgecp-U0.0/result.npz"
    arrays = audit._arrays(path)
    arrays["weiss_iw"][0] *= 1.01
    with pytest.raises(audit.AuditError, match="Weiss equation"):
        audit._numerics(
            arrays, report["validation"][0], report["config"], compressed=True
        )


@pytest.mark.parametrize(
    "attack",
    [
        "wrong_product",
        "wrong_observations",
        "empty_coefficients",
        "false_projection",
        "final_defect",
        "wrong_mesh",
        "false_causality",
    ],
)
def test_independent_saved_equation_and_projection_attacks(complete_fixture, attack):
    campaign, _, original = complete_fixture
    row = copy.deepcopy(original["validation"][9])
    path = (
        campaign
        / "validation"
        / f"repeat-{row['repeat']}-{row['method']}-U{row['interaction']}"
        / "result.npz"
    )
    arrays = audit._arrays(path)
    if attack == "wrong_product":
        arrays["sigma_tau_product"] += 0.1
    elif attack == "wrong_observations":
        arrays["frequency_fit_indices"] = arrays["indices"]
    elif attack == "empty_coefficients":
        arrays["g_coefficients"] = np.array([])
    elif attack == "false_projection":
        row["diagnostics"]["sigma_projection_residual"] = 0.01
    elif attack == "final_defect":
        row["history"][-1]["unmixed_defect"] = 5e-9
    elif attack == "wrong_mesh":
        arrays["indices"] = arrays["indices"] + 1
    elif attack == "false_causality":
        row["physics"]["positive_imag_green_max"] = 0
    with pytest.raises(audit.AuditError):
        audit._numerics(arrays, row, original["config"], compressed=True)


def test_final_error_is_recomputed_instead_of_trusting_reported_pass(complete_fixture):
    campaign, _, report = complete_fixture
    arrays = audit._arrays(campaign / "validation/repeat-0-kernelgecp-U0.0/result.npz")
    truth = audit._arrays(campaign / "U0_0/reference/validation.npz")
    arrays["validation_g_iw"][0] += 1e-3
    with pytest.raises(audit.AuditError, match="prediction misses"):
        audit._validation_errors(
            arrays, truth, report["validation"][0], report["config"]
        )


def test_short_zero_prediction_cannot_broadcast_across_final_holdouts(complete_fixture):
    campaign, _, report = complete_fixture
    arrays = audit._arrays(campaign / "validation/repeat-0-kernelgecp-U0.0/result.npz")
    truth = audit._arrays(campaign / "U0_0/reference/validation.npz")
    arrays["validation_sigma_iw"] = np.array([0j])
    with pytest.raises(audit.AuditError, match="prediction shape"):
        audit._validation_errors(
            arrays, truth, report["validation"][0], report["config"]
        )


def test_repetition_fingerprint_checks_arrays_not_only_recorded_fingerprints(
    complete_fixture,
):
    campaign, _, _ = complete_fixture
    arrays = audit._arrays(campaign / "validation/repeat-0-kernelgecp-U0.0/result.npz")
    before = audit._fingerprint(arrays)
    arrays["g_coefficients"][0] += 1e-13
    assert before != audit._fingerprint(arrays)


def test_changed_saved_array_cannot_pass_with_unchanged_reported_fingerprints(
    complete_fixture, tmp_path
):
    campaign, repository, original = complete_fixture
    copied = tmp_path / "copied"
    shutil.copytree(campaign, copied)
    report = copy.deepcopy(original)
    row = report["validation"][1]
    folder = copied / "validation/repeat-1-kernelgecp-U0.0"
    arrays = audit._arrays(folder / "result.npz")
    arrays["g_coefficients"][0] += 1e-13
    np.savez_compressed(folder / "result.npz", **arrays)
    row["artifacts"]["result.npz"] = audit._sha(folder / "result.npz")
    write_json(copied / "report.json", report)
    receipt = audit.audit_campaign(copied, repository=repository)
    assert not receipt["passed"]
    assert receipt["failures"][0].startswith("repetitions:")


@pytest.mark.parametrize("attack", ["sign", "normalization"])
def test_independent_finite_oracle_rejects_wrong_sign_or_unit_normalization(
    monkeypatch, attack
):
    from kernelgecp import lehmann

    real = lehmann.ipt_self_energy_poles

    def wrong(weiss, interaction):
        sigma = real(weiss, interaction)
        if attack == "normalization":
            return lehmann.FermionicPoleModel(
                sigma.beta, sigma.energies, sigma.weights / sigma.mass
            )

        class WrongSign:
            mass = sigma.mass

            def evaluate_matsubara(self, indices):
                return -sigma.evaluate_matsubara(indices)

        return WrongSign()

    monkeypatch.setattr(lehmann, "ipt_self_energy_poles", wrong)
    with pytest.raises(audit.AuditError, match="Fourier/sign/normalization"):
        audit.finite_pole_audit()


@pytest.mark.parametrize("value", [True, "0", None, float("nan"), float("inf")])
def test_numeric_metadata_does_not_coerce_invalid_values(value):
    with pytest.raises(audit.AuditError):
        audit._number(value, "test")


def test_cli_writes_failed_receipt_and_exits_nonzero(tmp_path):
    write_json(tmp_path / "report.json", {"schema": "wrong"})
    process = subprocess.run(
        [sys.executable, "-m", "experiments.audit_interacting_dmft", str(tmp_path)],
        cwd=audit.ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert process.returncode == 1
    receipt = json.loads((tmp_path / "audit.json").read_text())
    assert receipt["passed"] is False
    assert receipt["failures"]
