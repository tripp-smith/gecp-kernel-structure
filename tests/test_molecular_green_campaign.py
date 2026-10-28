"""Independent matrix-pole worker oracle and fail-closed campaign controls."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from experiments import molecular_green as campaign


def configuration():
    return campaign.load_config(campaign.ROOT / "experiments/molecular_green.toml")


class PoleBasis:
    rank = 1
    diagnostics = {"scope": "independent_analytic_test_basis"}

    def __init__(self, beta, energy):
        self.beta, self.energy = beta, energy

    def design_tau(self, tau):
        return (
            -np.exp(-self.energy * np.asarray(tau))
            / (1 + np.exp(-self.beta * self.energy))
        )[:, None]

    def design_matsubara(self, indices):
        return (
            1 / (1j * (2 * np.asarray(indices) + 1) * np.pi / self.beta - self.energy)
        )[:, None]


def worker_fixture(tmp_path, monkeypatch):
    config = configuration()
    beta, energy = config["beta"], 0.2
    tau = np.linspace(0, beta, 37)
    amplitude = np.array([[0.6, 0.1j], [-0.1j, 0.4]])[None, None]

    def exact_tau(points):
        return (-np.exp(-energy * np.asarray(points)) / (1 + np.exp(-beta * energy)))[
            :, None, None, None, None
        ] * amplitude

    def exact_iw(points):
        return (1 / (1j * (2 * np.asarray(points) + 1) * np.pi / beta - energy))[
            :, None, None, None, None
        ] * amplitude

    data = SimpleNamespace(
        beta=beta,
        tau=tau,
        g_tau=exact_tau(tau),
        sigma_tau=0.3 * exact_tau(tau),
        manifest_sha256="a" * 64,
        iteration=7,
    )
    times = campaign.query_tau(data, "calibration")
    indices = np.asarray(config["calibration_indices"], dtype=np.int64)
    reference = {
        "tau": times,
        "indices": indices,
        "green_tau": exact_tau(times),
        "green_iw": exact_iw(indices),
        "sigma_tau": 0.3 * exact_tau(times),
        "sigma_iw": 0.3 * exact_iw(indices),
    }
    path = tmp_path / "truth.npz"
    np.savez_compressed(path, **reference)
    monkeypatch.setattr(campaign, "_read", lambda folder: data)
    monkeypatch.setattr(campaign, "_basis", lambda *args: PoleBasis(beta, energy))
    actual_import = campaign.importlib.import_module

    def optional_import(name, *args, **kwargs):
        if name in ("pydlr", "sparse_ir"):
            return SimpleNamespace()
        return actual_import(name, *args, **kwargs)

    monkeypatch.setattr(campaign.importlib, "import_module", optional_import)
    request = {
        "config": config,
        "role": "calibration",
        "native_directory": str(tmp_path),
        "native_manifest_sha256": data.manifest_sha256,
        "case": "matrix-pole-oracle",
        "method": "kernelgecp",
        "cutoff": 4.0,
        "eps": 1e-8,
        "reference": str(path),
    }
    return request, data, reference, amplitude


def test_matrix_worker_matches_analytic_complex_pole_and_preserves_all_axes(
    tmp_path, monkeypatch
):
    request, data, reference, amplitude = worker_fixture(tmp_path, monkeypatch)
    output = tmp_path / "result.json"
    result = campaign._worker(request, output)
    assert campaign._passed(result, 1e-12)
    assert max(result["errors"].values()) < 2e-15
    assert result["density_error"] < 1e-14
    assert result["causality_imag_max"]["green"] < 0
    assert result["native_manifest_sha256"] == data.manifest_sha256
    with np.load(output.with_suffix(".npz"), allow_pickle=False) as arrays:
        assert arrays["g_coefficients"].shape == (1, 1, 1, 2, 2)
        np.testing.assert_allclose(arrays["g_coefficients"][0], amplitude)
        np.testing.assert_allclose(arrays["sigma_coefficients"][0], 0.3 * amplitude)
        for name in ("green_tau", "green_iw", "sigma_tau", "sigma_iw"):
            np.testing.assert_allclose(arrays[name], reference[name], atol=2e-15)
    assert result["scope"].endswith("not_native_solver_acceleration")


def test_role_truth_is_unavailable_until_fitting_and_predictions_complete(
    tmp_path, monkeypatch
):
    request, _, _, _ = worker_fixture(tmp_path, monkeypatch)
    # A missing truth artifact is not used to choose coefficients or silently
    # converted into zero residual. Its absence fails the completed worker.
    request["reference"] = str(tmp_path / "missing.npz")
    with pytest.raises(FileNotFoundError):
        campaign._worker(request, tmp_path / "result.json")
    assert not (tmp_path / "result.json").exists()


def test_worker_rejects_changed_native_provenance(tmp_path, monkeypatch):
    request, _, _, _ = worker_fixture(tmp_path, monkeypatch)
    request["native_manifest_sha256"] = "b" * 64
    with pytest.raises(ValueError, match="provenance changed"):
        campaign._worker(request, tmp_path / "result.json")


def test_short_reference_cannot_broadcast_as_a_complete_matrix_reference(
    tmp_path, monkeypatch
):
    request, _, reference, _ = worker_fixture(tmp_path, monkeypatch)
    reference["sigma_iw"] = np.zeros((1,))
    np.savez_compressed(request["reference"], **reference)
    with pytest.raises(ValueError, match="shape/values"):
        campaign._worker(request, tmp_path / "result.json")


@pytest.mark.parametrize("size", [37, 110, 144])
def test_observations_calibration_and_final_time_queries_are_disjoint(size):
    training, calibration = campaign.split_time_rows(size)
    assert 0 in training and size - 1 in training
    assert not np.intersect1d(training, calibration).size
    data = SimpleNamespace(tau=np.linspace(0, 20, size), beta=20)
    final = campaign.query_tau(data, "validation")
    assert final.ndim == 1 and len(final) > 64
    assert not np.any(np.isin(final, data.tau[training]))
    assert not np.any(np.isin(final, data.tau[calibration]))
    assert np.count_nonzero(~np.isin(final, data.tau)) == 64
    with pytest.raises(ValueError):
        campaign.query_tau(data, "unknown")


@pytest.mark.parametrize(
    "errors",
    [
        {},
        {"green_iw": 0},
        dict.fromkeys(("green_tau", "green_iw", "sigma_tau", "sigma_iw"), np.nan),
        dict.fromkeys(("green_tau", "green_iw", "sigma_tau", "sigma_iw"), -1),
        dict.fromkeys(("green_tau", "green_iw", "sigma_tau", "sigma_iw"), True),
    ],
)
def test_incomplete_nonfinite_negative_or_boolean_errors_cannot_pass(errors):
    assert not campaign._passed({"status": "completed", "errors": errors}, 1e-6)


def test_valid_worker_numbers_alone_do_not_mask_failed_status():
    errors = dict.fromkeys(("green_tau", "green_iw", "sigma_tau", "sigma_iw"), 0.0)
    assert campaign._passed({"status": "completed", "errors": errors}, 1e-6)
    assert not campaign._passed({"status": "timeout", "errors": errors}, 1e-6)


@pytest.mark.parametrize(
    "name",
    ["target", "calibration_target", "native_refinement_target", "native_dyson_target"],
)
def test_committed_controls_cannot_weaken_the_accepted_scientific_budgets(
    monkeypatch, name
):
    data = configuration()
    data[name] *= 2
    monkeypatch.setattr(campaign.tomllib, "loads", lambda _: data)
    with pytest.raises(ValueError, match="cannot be weakened"):
        campaign.load_config(Path(__file__))


def test_complex_json_roundtrip_and_nonfinite_rejection(tmp_path):
    path = tmp_path / "record.json"
    campaign.write(path, {"value": np.complex128(1 + 2j), "array": np.array([1j])})
    value = json.loads(path.read_text())
    assert value == {
        "value": {"real": 1, "imag": 2},
        "array": {"real": [0], "imag": [1]},
    }
    with pytest.raises(ValueError):
        campaign.write(path, {"invalid": np.float64(np.nan)})


@pytest.mark.parametrize(
    "pilot,calibration,mandatory,exit_code",
    [
        (False, True, False, 1),
        (False, True, True, None),
        (True, False, False, 1),
        (True, True, False, None),
    ],
)
def test_cli_exit_status_does_not_treat_incomplete_science_as_a_pass(
    tmp_path, monkeypatch, pilot, calibration, mandatory, exit_code
):
    monkeypatch.setattr(
        campaign,
        "run",
        lambda *args, **kwargs: {
            "calibration_gates_passed": calibration,
            "all_mandatory_passed": mandatory,
        },
    )
    argv = [
        "molecular_green",
        "--native-root",
        str(tmp_path),
        "--output-dir",
        str(tmp_path / "output"),
    ]
    if pilot:
        argv.append("--pilot")
    monkeypatch.setattr(campaign.sys, "argv", argv)
    if exit_code is None:
        campaign.main()
    else:
        with pytest.raises(SystemExit) as failure:
            campaign.main()
        assert failure.value.code == exit_code


@pytest.mark.parametrize(
    "attack",
    [
        "methods",
        "bonds",
        "grids",
        "repeat",
        "budget",
        "target",
        "duplicate",
        "overlap",
        "indices",
    ],
)
def test_campaign_controls_require_the_full_scope(monkeypatch, attack):
    data = copy.deepcopy(configuration())
    if attack == "methods":
        data["methods"] = ["kernelgecp"]
    elif attack == "bonds":
        data["bonds_angstrom"] = [0.74]
    elif attack == "grids":
        data["native_lambdas"] = [1e4]
    elif attack == "repeat":
        data["repeats"] = 1
    elif attack == "budget":
        data["worker_timeout_seconds"] = np.inf
    elif attack == "target":
        data["calibration_target"] = data["target"]
    elif attack == "duplicate":
        data["green_cutoffs_hartree"] = [4, 4]
    elif attack == "overlap":
        data["validation_indices"] = data["calibration_indices"]
    elif attack == "indices":
        data["validation_indices"] = [False]
    monkeypatch.setattr(campaign.tomllib, "loads", lambda _: data)
    with pytest.raises(ValueError):
        campaign.load_config(Path(__file__))
