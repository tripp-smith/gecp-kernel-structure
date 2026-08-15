import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from kernelgecp import (
    SyntheticApplicationConfig,
    fermionic_green_from_density,
    quadrature_weights,
    run_synthetic_application_suite,
    synthetic_spectral_density,
)

ROOT = Path(__file__).parents[1]


def quick_config() -> SyntheticApplicationConfig:
    return SyntheticApplicationConfig(
        cutoff=6.0,
        gecp_tolerance=1e-5,
        precision_bits=80,
        grid_order=6,
        max_rank=30,
        quadrature_order=401,
        validation_time_order=41,
        sparse_time_order=61,
        sparse_candidate_order=401,
        sparse_max_atoms=8,
        sparse_tolerance=1e-7,
        noise_sigma=1e-8,
        psd_point_count=32,
        psd_max_rank=30,
        psd_tolerance=1e-6,
        seed=20260815,
    )


def test_spectral_fixtures_are_physical_and_integrate() -> None:
    omega = np.linspace(-8, 8, 2_001)
    weights = quadrature_weights(omega)
    for fixture in ("hubbard_three_peak", "gapped_two_band"):
        density = synthetic_spectral_density(fixture, omega)
        assert np.all(density >= 0)
        assert float(weights @ density) == pytest.approx(1.0, rel=3e-4)
        green = fermionic_green_from_density(
            np.linspace(0, 1, 17), omega, density, cutoff=8
        )
        assert green.shape == (17,)
        assert np.all(np.isfinite(green))
    gapped = synthetic_spectral_density("gapped_two_band", omega)
    assert np.max(gapped[np.abs(omega) < 2]) == 0


def test_application_suite_verifies_transfer_recovery_and_psd_control() -> None:
    suite = run_synthetic_application_suite(quick_config(), git_commit="test")
    assert suite.schema_version == 1
    assert len(suite.config_hash) == 64
    assert all(record.gecp_converged for record in suite.green_compression)
    assert all(record.transfer_verified for record in suite.green_compression)
    assert all(
        record.green_validation_error <= record.transfer_bound + 1e-13
        for record in suite.green_compression
    )

    clean, noisy = suite.sparse_recovery
    assert clean.converged
    assert clean.candidate_strategy == "known_transition_library"
    assert clean.recovered_atom_count == clean.true_atom_count == 4
    assert clean.max_frequency_error == pytest.approx(0.0, abs=1e-10)
    assert clean.held_out_error < 1e-9
    assert noisy.candidate_strategy == "dense_blind_scan"
    assert noisy.held_out_error < 1e-5

    landmarks = suite.psd_landmarks
    assert landmarks.pivot_paths_agree
    assert landmarks.selected_rank < landmarks.point_count
    assert landmarks.cross_max_error < 2e-5
    assert landmarks.residual_min_eigenvalue > -1e-10

    first = suite.canonical_json()
    second = suite.canonical_json()
    assert first == second
    assert json.loads(first)["git_commit"] == "test"


def test_application_config_rejects_underresolved_inputs() -> None:
    with pytest.raises(ValueError):
        SyntheticApplicationConfig(quadrature_order=20)


def test_canonical_application_artifacts_match_documented_delivery() -> None:
    data_path = ROOT / "experiments/data/synthetic_applications.json"
    plot_path = ROOT / "experiments/data/synthetic_applications.png"
    assert hashlib.sha256(data_path.read_bytes()).hexdigest() == (
        "f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2"
    )
    assert hashlib.sha256(plot_path.read_bytes()).hexdigest() == (
        "0529263affd2fe10773f71784f06007650bff833374fbbebbb16f2507b2ab169"
    )
    delivered = json.loads(data_path.read_text(encoding="utf-8"))
    assert delivered["git_commit"] == "6404311e8c864d0b9ad9f655e8997f119194acff"
    assert delivered["config_hash"] == (
        "b41f4c17e76ccb108f2999b0e07af1003eb011bf8d96f2c7456e92b9d36f8d35"
    )
    assert [record["gecp_rank"] for record in delivered["green_compression"]] == [
        12,
        12,
    ]
    assert delivered["psd_landmarks"]["pivot_paths_agree"]
