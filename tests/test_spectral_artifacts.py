"""Independent reconstruction of every committed comparison's measured error."""

import hashlib
import json
import struct
import subprocess
from pathlib import Path

import numpy as np

from kernelgecp.applications import synthetic_spectral_density

ROOT = Path(__file__).parents[1]
DATA = ROOT / "experiments/data"
GENERATION_REVISION = "072ce1004de11422f05ef1fbe3e738093d6f2f13"
GENERATION_LOCK = "ec9204011e82b7d4bc7e1f626d4449ee7ed86d0c3e3e6375b90eea7fc026969c"
# Immutable generation manifest, not a requirement that all future source stays
# identical to DN. Shallow CI/archive checkouts may lack the historical blobs.
GENERATION_SOURCES = {
    "experiments/spectral_comparison.py": (
        "64549c6ca229f7192ba720a0b873ee7a2f7acec670d8e53b45a162c166b2dd36"
    ),
    "kernelgecp/comparison.py": (
        "47b5c43795ec987080200d311b7168de4df20d6b16f31c213c97c9456094fbb1"
    ),
    "kernelgecp/gecp.py": (
        "ea7b15a45c8f0f74eb85f1bc0cdcc043b282428810008e86e1eac3c9981a69a5"
    ),
    "kernelgecp/prony.py": (
        "70df18eecc1ad9fcceb721e6432d801b8c65abd9d2d0cc1c1d09fb42c88169c1"
    ),
    "kernelgecp/rational.py": (
        "bf1d00b00ad9dbbc241649ac1a0c13c1d89453efa27f5811fb49cc3d60a222f2"
    ),
    "kernelgecp/sparse.py": (
        "3cb62cb4612eba112f53949d1c9b2708bff731e023e8956203fbae9f091c41b5"
    ),
    "kernelgecp/spectral_common.py": (
        "ead7c4ddcc3ebe73fe60de5ba2f299040b97a6186c3d33b0be9027816ffbc8fb"
    ),
    "kernelgecp/vector_fitting.py": (
        "e6ed5fec88e85a311ed2fc3649f7fe8dd5b5f5f87174dbdd6d5365d46079b4cf"
    ),
}


def _generation_blob(path: str) -> bytes | None:
    try:
        result = subprocess.run(
            ["git", "show", f"{GENERATION_REVISION}:{path}"],
            cwd=ROOT,
            capture_output=True,
            check=False,
        )
    except FileNotFoundError:
        return None
    return result.stdout if result.returncode == 0 else None


def _load() -> dict:
    return json.loads((DATA / "spectral_comparison.json").read_text())


def _kernel(t: np.ndarray, omega: np.ndarray) -> np.ndarray:
    # Independent direct formula is safe on this fixed cutoff-eight cohort.
    return np.exp(-t[:, None] * omega) / (1 + np.exp(-omega))


def _truth(fixture: dict, t: np.ndarray, config: dict) -> np.ndarray:
    if fixture["true_frequencies"] is not None:
        return _kernel(t, np.array(fixture["true_frequencies"])) @ np.array(
            fixture["true_weights"]
        )
    omega = np.linspace(-config["cutoff"], config["cutoff"], config["density_count"])
    density = synthetic_spectral_density(fixture["name"], omega)
    masses = np.asarray(density).copy()
    masses[[0, -1]] /= 2
    masses /= np.sum(masses)  # Uniform spacing cancels in the normalization.
    return _kernel(t, omega) @ masses


def test_committed_comparison_complete_coverage_and_provenance() -> None:
    d = _load()
    assert d["schema_version"] == 1 and not d["git_dirty"]
    assert d["git_commit"] == GENERATION_REVISION
    assert (
        d["config_hash"]
        == "665e158dd17b5722ecd1e34b2d549cbe0dcd413613dd7f5d9c33a3af54cdeeb7"
    )
    assert len(d["methods"]) == 6 and len(d["fixtures"]) == 4
    assert d["validation_count_actual"] == 1024
    assert [r["rank"] for r in d["universal_ranks"]] == [11, 12]
    assert all(r["converged"] for r in d["universal_ranks"])
    assert d["algorithm_source_sha256"] == GENERATION_SOURCES
    for path, digest in GENERATION_SOURCES.items():
        blob = _generation_blob(path)
        if blob is not None:
            assert hashlib.sha256(blob).hexdigest() == digest
    assert d["lock_sha256"] == GENERATION_LOCK
    generation_lock = _generation_blob("uv.lock")
    if generation_lock is not None:
        assert hashlib.sha256(generation_lock).hexdigest() == GENERATION_LOCK
    for fixture in d["fixtures"]:
        assert len(fixture["results"]) == 48
        assert {
            (r["method"], r["requested_atom_budget"]) for r in fixture["results"]
        } == {
            (method, budget)
            for method in d["methods"]
            for budget in d["config"]["atom_budgets"]
        }
        assert all(
            r["training_input_hash"] == fixture["training_input_hash"]
            for r in fixture["results"]
        )
    assert "not_blind" in d["known_library_control"]["candidate_strategy"]


def test_every_reported_time_error_and_smallest_tested_support_recomputes() -> None:
    d = _load()
    training_t = np.array(d["training_times"])
    validation_t = np.array(d["validation_times"])
    assert not np.intersect1d(training_t, validation_t).size
    for fixture in d["fixtures"]:
        observed = np.array(fixture["observations"])
        truth = _truth(fixture, validation_t, d["config"])
        for r in fixture["results"]:
            omega, weights = np.array(r["frequencies"]), np.array(r["weights"])
            assert len(omega) == r["atom_count"] <= r["requested_atom_budget"]
            assert np.all(np.abs(omega) <= d["config"]["cutoff"])
            error = np.max(np.abs(_kernel(validation_t, omega) @ weights - truth))
            training_error = np.max(
                np.abs(_kernel(training_t, omega) @ weights - observed)
            )
            assert abs(error - r["held_out_error"]) < 5e-14
            assert abs(training_error - r["training_error"]) < 5e-14
            assert (
                abs(error / np.max(np.abs(truth)) - r["held_out_relative_error"])
                < 1e-12
            )
            assert r["negative_weight_count"] == int(np.count_nonzero(weights < 0))
            assert all(
                p["passed"] == (r["held_out_error"] <= p["tolerance"])
                for p in r["passes"]
            )
            assert not r["certified"] and r["precision_bits"] == 53
        for summary in fixture["smallest_tested"]:
            passing = [
                r["atom_count"]
                for r in fixture["results"]
                if r["method"] == summary["method"]
                and r["held_out_error"] <= summary["tolerance"]
            ]
            assert summary["smallest_tested_passing_atoms"] == min(
                passing, default=None
            )


def test_failed_aaa_fits_and_transform_error_are_not_hidden() -> None:
    d = _load()
    assert d["fixtures"][0]["true_atom_count"] == 4
    for fixture in d["fixtures"]:
        transform = fixture["transform"]
        assert transform["sample_count"] == 513
        assert transform["source_ladder_max_abs_difference"] > 0
        assert transform["oracle_scope"].startswith("evaluation_only")
        assert transform["certification_status"] == "not_certified"
        if fixture["name"] != "clean_off_grid":
            aaa = [r for r in fixture["results"] if r["method"] == "aaa"]
            assert all(not p["passed"] for r in aaa for p in r["passes"])
            assert any(r["native_diagnostics"]["rejected_pole_count"] > 0 for r in aaa)


def test_canonical_json_and_reviewed_png_hashes() -> None:
    assert (
        hashlib.sha256((DATA / "spectral_comparison.json").read_bytes()).hexdigest()
        == "9f85509c929bc0511b0f34f25af5188420902a165ddd3697631607d16fd5ec31"
    )
    png = (DATA / "spectral_comparison.png").read_bytes()
    assert (
        hashlib.sha256(png).hexdigest()
        == "0368a73f27aa4b15ed9c313cb73a7326157d706b6aea22ca5c798cd75fe834d2"
    )
    assert struct.unpack(">II", png[16:24]) == (1680, 1120)
