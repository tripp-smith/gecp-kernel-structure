"""Actual Green output, independent convergence checks and malformed inputs."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest

from kernelgecp.integrations.green import GreenFormatError, read_green_tau

h5py = pytest.importorskip("h5py")
FIXTURE = Path(__file__).parent / "data/green_h2_074"


def read(folder=FIXTURE, **kwargs):
    return read_green_tau(
        folder / "input.h5", folder / "sim.h5", folder / "bound-manifest.json", **kwargs
    )


def rewrite_manifest(folder, mutate):
    path = folder / "bound-manifest.json"
    data = json.loads(path.read_text())
    mutate(data)
    path.write_text(json.dumps(data, allow_nan=False))


def rehash(folder, key):
    def update(data):
        path = folder / data["files"][key]["path"]
        data["files"][key]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

    rewrite_manifest(folder, update)


@pytest.fixture
def copied(tmp_path):
    folder = tmp_path / "run"
    shutil.copytree(FIXTURE, folder)
    return folder


def test_actual_native_data_axes_static_terms_energies_and_provenance():
    data = read()
    assert data.iteration == data.source_iteration == 83
    assert data.converged
    assert data.beta == 20
    assert data.axes == ("tau", "spin", "k", "orbital_row", "orbital_column")
    assert data.g_tau.shape == data.sigma_tau.shape == (110, 1, 1, 2, 2)
    assert data.sigma_static.shape == data.overlap.shape == (1, 1, 2, 2)
    assert data.g_tau.dtype == data.sigma_tau.dtype == np.complex128
    np.testing.assert_array_equal(data.overlap[0, 0], np.eye(2))
    assert data.provider.green == data.provider.green_mbtools == "1.0.0"
    assert data.provider.grids == "0.3.0"
    assert data.provider.pyscf == "2.14.0"
    assert data.atoms_angstrom == (("H", 0, 0, 0), ("H", 0, 0, 0.74))
    assert data.electrons == 2
    assert data.energies.total_hartree == pytest.approx(-1.129672596129894)
    assert data.energies.one_body_hartree == pytest.approx(-2.4988280219409544)
    assert data.energies.nuclear_hartree == pytest.approx(0.7151043390810812)
    assert data.mu == pytest.approx(0.0462183952331543)
    assert abs(data.sigma_static[0, 0, 0, 0] - 0.91042814) < 1e-8
    assert np.max(abs(data.sigma_tau)) > 1e-3
    assert data.units.dynamic_selfenergy_tau == "hartree^2"
    assert data.convergence.energy_change_sum_hartree == pytest.approx(
        9.792607350012084e-11
    )
    assert data.convergence.green_change_max == pytest.approx(4.46208903159849e-11)
    assert data.certification_status == "not_certified_native_data_import"
    assert len(data.manifest_sha256) == 64
    assert (
        dict(data.file_sha256)["input"]
        == hashlib.sha256((FIXTURE / "input.h5").read_bytes()).hexdigest()
    )
    assert data.green_tau is data.g_tau
    assert data.dynamic_selfenergy_tau is data.sigma_tau


def test_all_returned_arrays_are_owned_and_immutable(copied):
    data = read(copied)
    expected = data.g_tau.copy()
    for field in data.__dataclass_fields__:
        value = getattr(data, field)
        if isinstance(value, np.ndarray):
            assert not value.flags.writeable
            with pytest.raises(ValueError):
                value.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        data.mu = 0
    with h5py.File(copied / "sim.h5", "r+") as f:
        f["iter83/G_tau/data"][:] = 0
    np.testing.assert_array_equal(data.g_tau, expected)
    with pytest.raises(GreenFormatError, match="hash mismatch"):
        read(copied)


def test_historical_iteration_is_explicit_and_not_labeled_converged():
    final, previous = read(), read(iteration=82)
    assert previous.iteration == 82 and previous.source_iteration == 83
    assert not previous.converged
    assert previous.convergence == final.convergence
    assert np.max(abs(previous.g_tau - final.g_tau)) > 0


def test_cell_serialized_object_is_never_read(monkeypatch):
    actual = h5py.Dataset.__getitem__

    def forbidden(dataset, key):
        if dataset.name == "/Cell":
            raise AssertionError("serialized Cell was accessed")
        return actual(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", forbidden)
    assert read().iteration == 83


@pytest.mark.parametrize("value", [True, 83.0, 0, 84, -1])
def test_invalid_iteration_rejected(value):
    with pytest.raises(GreenFormatError, match="iteration"):
        read(iteration=value)


@pytest.mark.parametrize(
    "attack",
    [
        "schema",
        "units",
        "version",
        "sources",
        "beta",
        "geometry",
        "basis",
        "method",
        "binary",
        "bad_hash",
        "path_escape",
    ],
)
def test_manifest_missing_or_inconsistent_metadata_fails_closed(copied, attack):
    def mutate(data):
        if attack == "schema":
            data["schema"] = "old"
        elif attack == "units":
            data["units"]["energy"] = "eV"
        elif attack == "version":
            data["provider"]["pyscf"] = "2.13.0"
        elif attack == "sources":
            data["provider"]["source_revisions"]["green-mbpt"] = "0" * 40
        elif attack == "beta":
            data["beta_inverse_hartree"] = 30
        elif attack == "geometry":
            data.pop("geometry")
        elif attack == "basis":
            data["storage_basis"] = "unknown"
        elif attack == "method":
            data["method"] = "GW"
        elif attack == "binary":
            data["provider"]["native_executable_sha256"] = "0" * 64
        elif attack == "bad_hash":
            data["files"]["results"]["sha256"] = "0" * 64
        elif attack == "path_escape":
            data["files"]["native_log"]["path"] = "../gf2.log"

    rewrite_manifest(copied, mutate)
    with pytest.raises(GreenFormatError):
        read(copied)


@pytest.mark.parametrize(
    "attack",
    [
        "incomplete",
        "moved_iter",
        "missing_previous",
        "float32",
        "nan",
        "mesh",
        "sigma_mesh",
        "shape",
        "soft_link",
    ],
)
def test_incomplete_or_malformed_output_cannot_fall_back(copied, attack):
    with h5py.File(copied / "sim.h5", "r+") as f:
        group = f["iter83"]
        if attack == "incomplete":
            del group["Sigma1"]
        elif attack == "moved_iter":
            f["iter"][()] = 84
        elif attack == "missing_previous":
            del f["iter82"]
        elif attack == "float32":
            value = group["mu"][()]
            del group["mu"]
            group.create_dataset("mu", data=np.float32(value))
        elif attack == "nan":
            group["G_tau/data"][0, 0, 0, 0, 0] = np.nan
        elif attack == "mesh":
            group["G_tau/mesh"][2] = 0
        elif attack == "sigma_mesh":
            group["Selfenergy/mesh"][2] += 1e-5
        elif attack == "shape":
            value = group["G_tau/data"][:, 0, 0]
            del group["G_tau/data"]
            group.create_dataset("G_tau/data", data=value)
        elif attack == "soft_link":
            group.move("G_tau/data", "actual_green")
            group["G_tau/data"] = h5py.SoftLink("/iter83/actual_green")
    rehash(copied, "results")
    with pytest.raises(GreenFormatError):
        read(copied)
    if attack == "moved_iter":
        with pytest.raises(GreenFormatError):
            read(copied, iteration=82)


@pytest.mark.parametrize(
    "attack", ["complex_marker", "overlap", "not_lowdin", "hermitian", "version"]
)
def test_native_input_contract_is_validated_after_hash_binding(copied, attack):
    with h5py.File(copied / "input.h5", "r+") as f:
        if attack == "complex_marker":
            del f["HF/S-k"].attrs["__complex__"]
        elif attack == "overlap":
            f["HF/S-k"][0, 0, 0, 0, 0] = -1
        elif attack == "not_lowdin":
            f["HF/S-k"][0, 0, 0, 0, 0] = 1.1
        elif attack == "hermitian":
            f["HF/H-k"][0, 0, 0, 1, 0] = 99
        elif attack == "version":
            f.attrs["__green_version__"] = "2.0.0"
    rehash(copied, "input")
    with pytest.raises(GreenFormatError):
        read(copied)


@pytest.mark.parametrize(
    "attack", ["exit_zero_only", "iteration_budget", "stored_deltas"]
)
def test_native_log_and_stored_pair_both_establish_convergence(copied, attack):
    path = copied / "gf2.log"
    if attack == "exit_zero_only":
        path.write_text("Completed.\n")
    elif attack == "iteration_budget":
        path.write_text(path.read_text() + "Reached Maximum number of iterations\n")
    else:
        with h5py.File(copied / "sim.h5", "r+") as f:
            f["iter82/Energy_HF"][()] += 0.01
        rehash(copied, "results")
    rehash(copied, "native_log")
    with pytest.raises(GreenFormatError, match="converg"):
        read(copied)


def test_module_import_is_lazy_and_missing_optional_provider_has_actionable_error():
    code = """
import builtins,sys
actual=builtins.__import__
def guarded(name,*args,**kwargs):
    if name.split('.')[0] in ('h5py','pyscf'):
        raise AssertionError('optional provider imported eagerly')
    return actual(name,*args,**kwargs)
builtins.__import__=guarded
import kernelgecp
import kernelgecp.integrations.green as adapter
assert 'h5py' not in sys.modules and 'pyscf' not in sys.modules
def unavailable(name):
    raise ImportError('provider intentionally absent')
adapter.importlib.import_module=unavailable
try:
    adapter.read_green_tau('input.h5','sim.h5','manifest.json')
except ImportError as error:
    assert 'kernelgecp[green]' in str(error)
else:
    raise AssertionError('missing provider accepted')
"""
    result = subprocess.run(
        [sys.executable, "-c", code], check=False, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
