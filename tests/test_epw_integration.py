"""Independent physical, native-order and ownership checks for EPW interchange."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from kernelgecp import FermionicBasis
from kernelgecp.integrations.epw import (
    KB_EV_PER_K,
    EPWAnisotropicData,
    EPWNambuFit,
    fit_epw_nambu,
    nambu_green,
    read_epw_anisotropic,
    recover_epw_nambu,
    write_epw_anisotropic,
)
from kernelgecp.types import NumericalBreakdownError

FIXTURES = Path(__file__).parent / "data/epw"


def _context(temperature=500, states=3, modes=None, full=True):
    indices = list(range(24)) if modes is None else list(modes)
    return {
        "schema": "kernelgecp.epw-anisotropic-input.v1",
        "provider": {"qe": "7.5", "epw": "6.0"},
        "temperature_kelvin": temperature,
        "units": "eV",
        "state_order": "frequency_major_ik_ibnd_ordinal",
        "state_count": states,
        "matsubara_indices": indices,
        "positive_only": min(indices) >= 0,
        "full_bandwidth": full,
        "chi_convention": "native_fbw_column" if full else "absent_chi_is_zero_non_fbw",
        "reference_fermi_energy_ev": 3.2,
        "chemical_potential_offset_ev": 0.03,
    }


def _data(*, full=True, signed=False):
    indices = np.arange(-12, 12) if signed else np.arange(24)
    context = _context(modes=indices.tolist(), full=full)
    frequencies = (
        (2 * indices + 1) * np.pi * KB_EV_PER_K * context["temperature_kelvin"]
    )
    z = np.ones((24, 3))
    gap = np.broadcast_to([0.02, -0.03, 0.01], z.shape).copy()
    chi = np.broadcast_to([0.01, 0.02, -0.015], z.shape).copy() if full else None
    return EPWAnisotropicData(
        frequencies, [-0.4, 0.2, 0.2], z, gap, chi, context=context
    )


@pytest.fixture(scope="module")
def basis():
    return FermionicBasis(
        beta=1 / (KB_EV_PER_K * 500),
        wmax=1.5,
        tol=1e-10,
        precision="reference",
        precision_bits=128,
    )


def test_bcs_nambu_resolvent_matches_independent_hermitian_eigensystem():
    data = _data()
    observed = nambu_green(data)
    for state in range(3):
        xi = data.state_energies_ev[state] - 0.03 + data.chi_ev[0, state]
        delta = data.delta_ev[0, state]
        h = np.array([[xi, delta], [delta, -xi]])
        energies, vectors = np.linalg.eigh(h)
        for row, frequency in enumerate(data.frequencies_ev):
            expected = (
                vectors * (1 / (1j * frequency - energies))[None, :]
            ) @ vectors.T
            np.testing.assert_allclose(observed[row, state], expected, atol=3e-15)
    recovered = recover_epw_nambu(data, observed)
    np.testing.assert_allclose(recovered.data.z, data.z, atol=1e-14)
    np.testing.assert_allclose(recovered.data.delta_ev, data.delta_ev, atol=1e-14)
    np.testing.assert_allclose(recovered.data.chi_ev, data.chi_ev, atol=1e-14)
    assert recovered.diagnostics.solve_residual_max < 1e-14


def test_nambu_scalar_complex_pole_formula_with_dynamic_fields():
    data = _data()
    z = 1.2 + np.arange(24)[:, None] / 80 + np.zeros((24, 3))
    delta = data.delta_ev * (1 + np.arange(24)[:, None] / 50)
    chi = data.chi_ev + np.arange(24)[:, None] / 200
    data = data.replace_fields(z, delta, chi)
    iwz = 1j * data.frequencies_ev[:, None] * z
    xi = data.state_energies_ev[None, :] - 0.03 + chi
    determinant = iwz**2 - xi**2 - (z * delta) ** 2
    expected = np.empty((*z.shape, 2, 2), dtype=complex)
    expected[..., 0, 0] = (iwz + xi) / determinant
    expected[..., 1, 1] = (iwz - xi) / determinant
    expected[..., 0, 1] = expected[..., 1, 0] = z * delta / determinant
    np.testing.assert_allclose(nambu_green(data), expected, rtol=1e-14, atol=1e-14)


def test_signed_mesh_retains_actual_order_and_adjoint_relation():
    data = _data(signed=True)
    values = nambu_green(data)
    np.testing.assert_allclose(
        values[:12], values[12:][::-1].swapaxes(-1, -2).conj(), atol=1e-14
    )
    assert data.indices.tolist() == list(range(-12, 12))


@pytest.mark.parametrize("temperature", [17, 20, 30])
def test_authentic_native_fixture_preserves_bytes_fields_and_all_state_ordinals(
    tmp_path, temperature
):
    record = json.loads((FIXTURES / "manifest.json").read_text())
    entry = next(
        value
        for value in record["files"]
        if value["context"]["temperature_kelvin"] == temperature
    )
    context = entry["context"]
    raw = gzip.decompress((FIXTURES / entry["compressed_file"]).read_bytes())
    source = tmp_path / context["native_filename"]
    source.write_bytes(raw)
    data = read_epw_anisotropic(source, context)
    assert data.source_bytes_verified
    assert data.z.shape == (len(context["matsubara_indices"]), 140)
    assert data.state_energies_ev[0] == -12.526663257
    target = tmp_path / "replayed"
    manifest = write_epw_anisotropic(data, target)
    assert target.read_bytes() == raw
    assert manifest["native_sha256"] == hashlib.sha256(raw).hexdigest()
    recovered = recover_epw_nambu(data, nambu_green(data))
    np.testing.assert_allclose(recovered.data.z, data.z, atol=5e-13)
    np.testing.assert_allclose(recovered.data.delta_ev, data.delta_ev, atol=1e-14)
    np.testing.assert_allclose(recovered.data.chi_ev, data.chi_ev, atol=1e-13)
    assert recovered.diagnostics.condition_max < 1 + 1e-12
    assert recovered.diagnostics.z_derivative_gain_max > 1e4
    assert not recovered.data.source_bytes_verified


def test_authentic_denser_scout_file_retains_source_order_and_bytes(tmp_path):
    record = json.loads((FIXTURES / "manifest-scout-wscut1.json").read_text())
    entry = record["files"][0]
    context = entry["context"]
    raw = gzip.decompress((FIXTURES / entry["compressed_file"]).read_bytes())
    source = tmp_path / context["native_filename"]
    source.write_bytes(raw)
    data = read_epw_anisotropic(source, context)
    assert data.z.shape == (109, 140)
    assert data.context["chemical_potential_offset_ev"] == 0
    target = tmp_path / "replayed"
    write_epw_anisotropic(data, target)
    assert target.read_bytes() == raw


@pytest.mark.parametrize("full", [False, True])
def test_native_writer_rounding_degenerate_energies_and_explicit_non_fbw_shift(
    tmp_path, full
):
    data = _data(full=full)
    source = tmp_path / "native"
    context = write_epw_anisotropic(data, source)
    parsed = read_epw_anisotropic(source, context)
    assert parsed.state_energies_ev.tolist() == [-0.4, 0.2, 0.2]
    assert parsed.z.shape == (24, 3)
    assert parsed.delta_ev[0, 1] != parsed.delta_ev[0, 2]
    with pytest.raises(FileExistsError):
        write_epw_anisotropic(data, source)
    if not full:
        assert np.all(parsed.chi_ev == 0)
        recovered = recover_epw_nambu(parsed, nambu_green(parsed))
        np.testing.assert_array_equal(recovered.data.chi_ev, 0)
        assert recovered.diagnostics.non_fbw_shift_zero_defect_ev < 1e-14


def test_owned_arrays_context_and_return_values_are_deeply_immutable():
    context = _context()
    frequencies = (2 * np.arange(24) + 1) * np.pi * KB_EV_PER_K * 500
    z = np.ones((24, 3))
    gap = np.zeros((24, 3))
    data = EPWAnisotropicData(
        frequencies, [-0.4, 0.2, 0.2], z, gap, gap, context=context
    )
    frequencies[:] = 0
    z[:] = 0
    context["provider"]["epw"] = "wrong"
    assert np.all(data.frequencies_ev > 0) and np.all(data.z == 1)
    assert data.context["provider"]["epw"] == "6.0"
    for array in (
        data.frequencies_ev,
        data.state_energies_ev,
        data.z,
        data.indices,
        nambu_green(data),
    ):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(AttributeError):
        data.z = np.zeros((24, 3))


def test_shared_complex_frequency_fit_completion_and_safe_archive(basis, tmp_path):
    data = _data()
    rows = np.array([j for j in range(24) if j % 4 != 1])
    result = fit_epw_nambu(
        data,
        basis,
        fit_rows=rows,
        rcond=1e-13,
        complete_adjoint=True,
        analysis_provenance={"case": "independent BCS"},
    )
    assert result.fit.observation_domain == "matsubara"
    assert result.fit.coefficients.shape == (basis.rank, 3, 2, 2)
    assert result.completion_policy == "explicit_adjoint_real_eliashberg"
    np.testing.assert_allclose(
        result.evaluate_matsubara(data.indices[rows]),
        nambu_green(data)[rows],
        atol=2e-9,
    )
    positive = fit_epw_nambu(data, basis, rcond=1e-10)
    assert positive.completion_policy == "none"
    manifest = tmp_path / "nambu.json"
    result.save(manifest)
    loaded = EPWNambuFit.load(manifest)
    assert loaded.source_context == result.source_context
    np.testing.assert_array_equal(loaded.fit.coefficients, result.fit.coefficients)
    np.testing.assert_allclose(
        loaded.evaluate_matsubara(data.indices),
        result.evaluate_matsubara(data.indices),
        atol=1e-14,
    )
    with pytest.raises(FileExistsError):
        result.save(manifest)
    original = manifest.read_bytes()
    record = json.loads(original)
    record["layout"]["context"]["chemical_potential_offset_ev"] = 0.04
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="checksum"):
        EPWNambuFit.load(manifest)
    manifest.write_bytes(original)
    archive = manifest.with_suffix(".npz")
    archive.write_bytes(archive.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="hash"):
        EPWNambuFit.load(manifest)


@pytest.mark.parametrize(
    "attack",
    [
        "missing_mu",
        "wrong_units",
        "wrong_provider",
        "boolean_temperature",
        "wrong_temperature",
        "boolean_index",
        "incomplete_signed",
        "missing_chi_convention",
        "unknown_field",
    ],
)
def test_invalid_context_cannot_silently_change_native_conventions(attack):
    data = _data()
    context = data.context
    if attack == "missing_mu":
        context.pop("chemical_potential_offset_ev")
    elif attack == "wrong_units":
        context["units"] = "meV"
    elif attack == "wrong_provider":
        context["provider"]["epw"] = "5"
    elif attack == "boolean_temperature":
        context["temperature_kelvin"] = True
    elif attack == "wrong_temperature":
        context["temperature_kelvin"] = 510
    elif attack == "boolean_index":
        context["matsubara_indices"][0] = True
    elif attack == "incomplete_signed":
        context["positive_only"] = False
    elif attack == "missing_chi_convention":
        context.pop("chi_convention")
    elif attack == "unknown_field":
        context["guess_missing_fields"] = True
    with pytest.raises(ValueError):
        EPWAnisotropicData(
            data.frequencies_ev,
            data.state_energies_ev,
            data.z,
            data.delta_ev,
            data.chi_ev,
            context=context,
        )


@pytest.mark.parametrize(
    "attack",
    [
        "wrong_hash",
        "short_field",
        "truncated",
        "last_state_frequency",
        "state_order",
        "zero_z",
        "nonfinite",
    ],
)
def test_native_reader_rejects_malformed_rows_including_last_state_frequency(
    tmp_path, attack
):
    source = tmp_path / "native"
    context = write_epw_anisotropic(_data(), source)
    lines = source.read_text().splitlines()
    if attack == "short_field":
        lines[1] = lines[1][1:]
    elif attack == "truncated":
        lines.pop()
    elif attack == "last_state_frequency":
        lines[3] = f"{0.02:20.10E}" + lines[3][20:]
    elif attack == "state_order":
        lines[4], lines[5] = lines[5], lines[4]
    elif attack == "zero_z":
        lines[1] = lines[1][:40] + f"{0:20.10E}" + lines[1][60:]
    elif attack == "nonfinite":
        lines[1] = lines[1][:40] + f"{float('nan'):20.10E}" + lines[1][60:]
    raw = ("\n".join(lines) + "\n").encode()
    source.write_bytes(raw)
    if attack == "wrong_hash":
        context["native_sha256"] = "0" * 64
    else:
        context["native_sha256"] = hashlib.sha256(raw).hexdigest()
    with pytest.raises(ValueError):
        read_epw_anisotropic(source, context)


def test_recovery_rejects_singular_nonreal_and_large_absolute_amplification():
    data = _data()
    values = nambu_green(data)
    for bad in (np.zeros_like(values), values * 1e-10):
        with pytest.raises(NumericalBreakdownError, match="conditioning"):
            recover_epw_nambu(data, bad)
    changed = values.copy()
    changed[:, :, 0, 1] += 0.01j
    with pytest.raises(NumericalBreakdownError, match="non-real/asymmetric"):
        recover_epw_nambu(data, changed)
    with pytest.raises(NumericalBreakdownError, match="amplification"):
        recover_epw_nambu(data, values, max_z_amplification=1)


def test_completion_is_explicit_and_fit_rows_cannot_escape_source(basis):
    data = _data()
    for rows in ([True], [-1], [24], [0, 0]):
        with pytest.raises(ValueError):
            fit_epw_nambu(data, basis, fit_rows=rows)
    with pytest.raises(ValueError, match="positive-only"):
        fit_epw_nambu(_data(signed=True), basis, complete_adjoint=True)


def test_fit_rejects_undersampling_and_conditioning(basis):
    data = _data()
    with pytest.raises(NumericalBreakdownError, match="undersampled"):
        fit_epw_nambu(data, basis, fit_rows=[0])
    with pytest.raises(NumericalBreakdownError, match="conditioning"):
        fit_epw_nambu(data, basis, complete_adjoint=True, max_retained_condition=1)
