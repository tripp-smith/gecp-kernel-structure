"""Independent checks of the explicit real-parity EPW frequency adapter."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import kernelgecp.integrations.epw_frequency as epw_frequency
from kernelgecp._basis_io import read_archive, write_archive
from kernelgecp.integrations.epw import (
    KB_EV_PER_K,
    EPWAnisotropicData,
    read_epw_anisotropic,
)
from kernelgecp.integrations.epw_frequency import (
    EPWFrequencyCrossFit,
    fit_epw_frequency_cross,
)
from kernelgecp.types import NumericalBreakdownError

FIXTURES = Path(__file__).parent / "data/epw"


def _static_data():
    indices = np.arange(24)
    context = {
        "schema": "kernelgecp.epw-anisotropic-input.v1",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "temperature_kelvin": 500,
        "units": "eV",
        "state_order": "frequency_major_ik_ibnd_ordinal",
        "state_count": 3,
        "matsubara_indices": indices.tolist(),
        "positive_only": True,
        "full_bandwidth": True,
        "chi_convention": "native_fbw_column",
        "reference_fermi_energy_ev": 3.2,
        "chemical_potential_offset_ev": 0.03,
    }
    omega = (2 * indices + 1) * np.pi * KB_EV_PER_K * 500
    energies = np.array([0.63, 0.03, -0.77])
    chi = np.broadcast_to([0.01, 0.0, 0.0], (24, 3)).copy()
    gap = np.broadcast_to([np.sqrt(1 - 0.61**2), 1.0, 0.6], (24, 3)).copy()
    return EPWAnisotropicData(
        omega, energies, np.ones((24, 3)), gap, chi, context=context
    )


def _fit_static(data, **kwargs):
    return fit_epw_frequency_cross(
        data,
        real_nambu_parity=True,
        fit_rows=np.arange(0, 24, 2),
        rank=3,
        dictionary_size=101,
        pole_cutoff_ev=1,
        rcond=1e-12,
        **kwargs,
    )


def _analytic_green(data, indices):
    omega = (2 * np.asarray(indices) + 1) * np.pi / data.beta
    output = np.empty((len(omega), len(data.state_energies_ev), 2, 2), complex)
    for j, energy in enumerate(data.state_energies_ev):
        xi = energy - data.context["chemical_potential_offset_ev"] + data.chi_ev[0, j]
        gap = data.delta_ev[0, j]
        eigenvalues, eigenvectors = np.linalg.eigh([[xi, gap], [gap, -xi]])
        for i, w in enumerate(omega):
            output[i, j] = (eigenvectors / (1j * w - eigenvalues)) @ eigenvectors.T
    return output


def test_static_resolvent_heldout_adjoint_and_recovery():
    data = _static_data()
    fit = _fit_static(data)
    modes = np.arange(1, 24, 2)
    np.testing.assert_allclose(
        fit.evaluate_matsubara(modes), _analytic_green(data, modes), atol=1e-11
    )
    negative = -modes - 1
    np.testing.assert_allclose(
        fit.evaluate_matsubara(negative),
        fit.evaluate_matsubara(modes).conj().swapaxes(-1, -2),
        atol=1e-13,
    )
    recovered = fit.reconstruct_native()
    np.testing.assert_allclose(recovered.data.z, data.z, atol=1e-10)
    np.testing.assert_allclose(recovered.data.delta_ev, data.delta_ev, atol=1e-10)
    np.testing.assert_allclose(recovered.data.chi_ev, data.chi_ev, atol=1e-10)


def test_rows_isolation_and_ownership(monkeypatch):
    data = _static_data()
    calls = []
    original = epw_frequency.nambu_green

    def track(data, *, rows=None):
        calls.append(np.array(rows, copy=True))
        return original(data, rows=rows)

    monkeypatch.setattr(epw_frequency, "nambu_green", track)
    fit = _fit_static(data)
    np.testing.assert_array_equal(calls[0], np.arange(0, 24, 2))
    assert len(calls) == 1
    context = fit.source_context
    frequencies = data.frequencies_ev.copy()
    energies = data.state_energies_ev.copy()
    coefficients = fit.coefficients.copy()
    rows = fit.fit_rows.copy()
    manual = EPWFrequencyCrossFit(
        fit.plan,
        coefficients,
        frequencies_ev=frequencies,
        state_energies_ev=energies,
        context=context,
        fit_rows=rows,
        diagnostics=fit.diagnostics,
    )
    coefficients[:] = 0
    frequencies[:] = 0
    energies[:] = 0
    rows[:] = 0
    context["temperature_kelvin"] = 1
    np.testing.assert_allclose(
        manual.evaluate_matsubara([1]), fit.evaluate_matsubara([1])
    )
    for arr in (
        manual.coefficients,
        manual.fit_rows,
        manual.indices,
        manual.evaluate_matsubara([1]),
    ):
        with pytest.raises(ValueError):
            arr.setflags(write=True)
    with pytest.raises(AttributeError):
        manual.plan.rank = 1
    assert manual.plan.prepare(rcond=1e-12) is manual.plan.prepare(rcond=1e-12)


def test_input_rejections():
    data = _static_data()
    with pytest.raises(TypeError):
        fit_epw_frequency_cross(data)
    for flag in (False, 1, None):
        with pytest.raises(ValueError):
            fit_epw_frequency_cross(data, real_nambu_parity=flag)
    for rows in ([True], [-1], [24], [0, 0]):
        with pytest.raises(ValueError):
            fit_epw_frequency_cross(data, real_nambu_parity=True, fit_rows=rows, rank=1)
    with pytest.raises(ValueError):
        fit_epw_frequency_cross(data, real_nambu_parity=True, rank=25)
    with pytest.raises(NumericalBreakdownError):
        _fit_static(data, max_retained_condition=1)
    fit = _fit_static(data)
    kwargs = dict(
        frequencies_ev=data.frequencies_ev,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=fit.fit_rows,
        diagnostics=fit.diagnostics,
    )
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit(fit.plan, fit.coefficients[:, :, :2], **kwargs)
    bad = dict(kwargs)
    bad["context"] = {**data.context, "temperature_kelvin": 501}
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit(fit.plan, fit.coefficients, **bad)
    bad["context"] = {k: v for k, v in data.context.items() if k != "chi_convention"}
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit(fit.plan, fit.coefficients, **bad)
    with pytest.raises(ValueError):
        fit.evaluate_matsubara([2**51])


def _rehash_manifest(path, record):
    record.pop("record_sha256", None)
    canonical = json.dumps(
        record, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    record["record_sha256"] = hashlib.sha256(canonical).hexdigest()
    path.write_text(json.dumps(record))


def test_archive_roundtrip_and_attacks(tmp_path):
    fit = _fit_static(_static_data())
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    for directory in (a, b):
        fit.save(directory / "fit.json")
    assert (a / "fit.json").read_bytes() == (b / "fit.json").read_bytes()
    assert (a / "fit.npz").read_bytes() == (b / "fit.npz").read_bytes()
    loaded = EPWFrequencyCrossFit.load(a / "fit.json")
    record = json.loads((a / "fit.json").read_text())
    assert record["representation_units"]["channels_and_coefficients"] == [
        "eV^-1", "eV^-2", "eV^-1"
    ]
    np.testing.assert_array_equal(loaded.coefficients, fit.coefficients)
    np.testing.assert_array_equal(
        loaded.evaluate_matsubara([1, 3]), fit.evaluate_matsubara([1, 3])
    )
    with pytest.raises(FileExistsError):
        fit.save(a / "fit.json")
    record["representation_units"]["green_function"] = "eV"
    _rehash_manifest(a / "fit.json", record)
    with pytest.raises(ValueError, match="unsupported"):
        EPWFrequencyCrossFit.load(a / "fit.json")
    fit.save(a / "fit.json", overwrite=True)
    raw = (a / "fit.npz").read_bytes()
    (a / "fit.npz").write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit.load(a / "fit.json")
    (a / "fit.npz").write_bytes(raw)
    record = json.loads((a / "fit.json").read_text())
    record["archive"]["file"] = "../fit.npz"
    _rehash_manifest(a / "fit.json", record)
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit.load(a / "fit.json")
    fit.save(a / "fit.json", overwrite=True)
    metadata, arrays = read_archive(a / "fit.npz")
    metadata["completion_policy"] = "implicit"
    write_archive(a / "fit.npz", metadata, arrays)
    record = json.loads((a / "fit.json").read_text())
    record["archive"]["sha256"] = hashlib.sha256(
        (a / "fit.npz").read_bytes()
    ).hexdigest()
    _rehash_manifest(a / "fit.json", record)
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit.load(a / "fit.json")
    fit.save(a / "fit.json", overwrite=True)
    metadata, arrays = read_archive(a / "fit.npz")
    arrays["weights"][0] *= 2
    write_archive(a / "fit.npz", metadata, arrays)
    record = json.loads((a / "fit.json").read_text())
    record["archive"]["sha256"] = hashlib.sha256(
        (a / "fit.npz").read_bytes()
    ).hexdigest()
    _rehash_manifest(a / "fit.json", record)
    with pytest.raises(ValueError):
        EPWFrequencyCrossFit.load(a / "fit.json")


def test_actual_17k_calibration_only_and_total_storage(tmp_path, monkeypatch):
    entry = json.loads((FIXTURES / "manifest-scout-wscut1.json").read_text())["files"][
        0
    ]
    path = tmp_path / entry["context"]["native_filename"]
    path.write_bytes(
        gzip.decompress((FIXTURES / entry["compressed_file"]).read_bytes())
    )
    data = read_epw_anisotropic(path, entry["context"])
    calibration = np.array([36, 40, 44])
    reserved = {45, 49, 53}
    train = np.array(
        [i for i in range(len(data.indices)) if i not in set(calibration) | reserved]
    )
    original = epw_frequency.nambu_green

    def track(data, *, rows=None):
        np.testing.assert_array_equal(rows, train)
        return original(data, rows=rows)

    monkeypatch.setattr(epw_frequency, "nambu_green", track)
    fit = fit_epw_frequency_cross(
        data,
        real_nambu_parity=True,
        fit_rows=train,
        rank=60,
        pole_cutoff_ev=13,
        precision_bits=256,
        rcond=1e-12,
    )
    fit.save(tmp_path / "calibration.json")
    fit = EPWFrequencyCrossFit.load(tmp_path / "calibration.json")
    predicted = fit.evaluate_matsubara(data.indices[calibration])
    w = data.frequencies_ev[calibration, None]
    z, chi, gap = (
        data.z[calibration],
        data.chi_ev[calibration],
        data.delta_ev[calibration],
    )
    xi = data.state_energies_ev[None, :] - data.context["chemical_potential_offset_ev"]
    matrix = np.empty(predicted.shape, complex)
    matrix[..., 0, 0], matrix[..., 1, 1] = 1j * w * z - xi - chi, 1j * w * z + xi + chi
    matrix[..., 0, 1] = matrix[..., 1, 0] = -z * gap
    truth = np.linalg.solve(matrix, np.broadcast_to(np.eye(2), matrix.shape))
    inverse = np.linalg.solve(predicted, np.broadcast_to(np.eye(2), predicted.shape))
    recovered_z = (inverse[..., 0, 0] + inverse[..., 1, 1]) / (2j * w)
    recovered_chi = (inverse[..., 1, 1] - inverse[..., 0, 0]) / 2 - xi
    recovered_gap = -(inverse[..., 0, 1] + inverse[..., 1, 0]) / (2 * recovered_z)
    errors = (
        np.max(abs(predicted - truth)) / np.max(abs(truth)),
        np.max(abs(recovered_z - z)),
        np.max(abs(recovered_gap - gap)),
        np.max(abs(recovered_chi - chi)),
    )
    assert np.all(np.array(errors) < [1e-6, 1e-5, 1e-6, 1e-6])
    assert np.max(abs(recovered_z.imag)) < 1e-8
    assert np.max(abs(recovered_gap.imag)) < 1e-8
    assert np.max(abs(recovered_chi.imag)) < 1e-8
    assert np.ptp(z, axis=1).max() > 1e-4
    assert np.ptp(gap, axis=1).max() > 1e-4
    assert np.ptp(chi, axis=1).max() > 1e-4
    archive_bytes = (tmp_path / "calibration.json").stat().st_size + (
        tmp_path / "calibration.npz"
    ).stat().st_size
    native_bytes = sum(
        a.nbytes
        for a in (
            data.z,
            data.chi_ev,
            data.delta_ev,
            data.frequencies_ev,
            data.state_energies_ev,
        )
    )
    print(
        "EPW calibration errors",
        errors,
        "archive bytes",
        archive_bytes,
        "native bytes",
        native_bytes,
    )
    assert archive_bytes < native_bytes
