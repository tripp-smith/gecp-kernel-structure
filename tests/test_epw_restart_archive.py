"""Exact restart fields are additional information, not pure-fit evidence."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from test_epw_frequency import _fit_static, _static_data

from kernelgecp._basis_io import read_archive, write_archive
from kernelgecp.integrations.epw import EPWAnisotropicData, read_epw_anisotropic
from kernelgecp.integrations.epw_frequency import EPWFrequencyCrossFit
from kernelgecp.integrations.epw_restart_archive import EPWRestartArchive, _baseline
from kernelgecp.types import NumericalBreakdownError


def test_actual_macos_fixture_decodes_on_the_current_platform(tmp_path):
    fixture = Path(__file__).parent / "data/epw/restart20-a9c2c39"
    provenance = json.loads((fixture / "provenance.json").read_text())
    for name, key in (
        ("bundle.json", "bundle_json_sha256"),
        ("bundle.npz", "bundle_npz_sha256"),
    ):
        assert (
            hashlib.sha256((fixture / name).read_bytes()).hexdigest() == provenance[key]
        )
    archive = EPWRestartArchive.load(fixture / "bundle.json")
    native = archive.reconstruct_native()
    assert native.z.shape == (provenance["frequency_count"], provenance["state_count"])
    assert (
        native.context["source_provenance"]["source_revision"]
        == (provenance["native_generation_revision"])
    )
    assert native.temperature_kelvin == 20
    assert native.context["chemical_potential_offset_ev"] == 0
    output = tmp_path / "restored-native.dat"
    archive.write_native_restart(output)
    assert (
        hashlib.sha256(output.read_bytes()).hexdigest()
        == (provenance["original_native_sha256"])
    )
    # Byte equality verifies the original native row order and rounding as well.


def _archive():
    data = _static_data()
    # Include negative zero as an exact bit-level test.
    gap = data.delta_ev.copy()
    chi = data.chi_ev.copy()
    chi[0, 1] = -0.0
    data = EPWAnisotropicData(
        data.frequencies_ev,
        data.state_energies_ev,
        data.z,
        gap,
        chi,
        context=data.context,
    )
    return data, EPWRestartArchive.from_data(
        data, _fit_static(data), policy="lossless_native_xor_v1"
    )


def test_independent_bcs_baseline_and_exact_roundtrip(tmp_path):
    data, archive = _archive()
    # The independent closed-form BCS fields are the synthetic input constants.
    baseline = _baseline(archive.fit)
    np.testing.assert_allclose(baseline[..., 0], data.z, atol=1e-10)
    np.testing.assert_allclose(baseline[..., 1], data.delta_ev, atol=1e-10)
    np.testing.assert_allclose(baseline[..., 2], data.chi_ev, atol=1e-10)
    manifest = tmp_path / "restart.json"
    archive.save(manifest)
    loaded = EPWRestartArchive.load(manifest)
    reconstructed = loaded.reconstruct_native()
    for name in ("z", "delta_ev", "chi_ev"):
        np.testing.assert_array_equal(
            getattr(reconstructed, name).view("u8"), getattr(data, name).view("u8")
        )
    assert np.signbit(reconstructed.chi_ev[0, 1])
    context = loaded.write_native_restart(tmp_path / "native.dat")
    reread = read_epw_anisotropic(tmp_path / "native.dat", context)
    np.testing.assert_array_equal(reread.chi_ev.view("u8"), data.chi_ev.view("u8"))
    with pytest.raises(ValueError):
        loaded._corrections.setflags(write=True)
    with pytest.raises(ValueError):
        loaded.fit.coefficients.setflags(write=True)
    assert loaded.fit.plan is not archive.fit.plan


def test_explicit_policy_layout_and_recovery_failures(monkeypatch):
    data = _static_data()
    fit = _fit_static(data)
    with pytest.raises(TypeError):
        EPWRestartArchive.from_data(data, fit)
    with pytest.raises(ValueError, match="policy"):
        EPWRestartArchive.from_data(data, fit, policy="lossy")
    changed = EPWAnisotropicData(
        data.frequencies_ev,
        data.state_energies_ev,
        data.z,
        data.delta_ev,
        data.chi_ev,
        context={**data.context, "reference_fermi_energy_ev": 4.0},
    )
    with pytest.raises(ValueError, match="layouts"):
        EPWRestartArchive.from_data(changed, fit, policy="lossless_native_xor_v1")
    broken = EPWFrequencyCrossFit(
        fit.plan,
        np.zeros_like(fit.coefficients),
        frequencies_ev=data.frequencies_ev,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=fit.fit_rows,
        diagnostics=fit.diagnostics,
    )
    with pytest.raises(NumericalBreakdownError):
        EPWRestartArchive.from_data(data, broken, policy="lossless_native_xor_v1")
    monkeypatch.setattr("kernelgecp.integrations.epw_restart_archive._MAX_WORK", 1)
    with pytest.raises(ValueError, match="resource"):
        EPWRestartArchive.from_data(data, fit, policy="lossless_native_xor_v1")


def test_hash_failures_and_no_reselection(tmp_path, monkeypatch):
    _, archive = _archive()
    manifest = tmp_path / "a.json"
    archive.save(manifest)
    monkeypatch.setattr(
        "kernelgecp.integrations._frequency_cross.CauchyFrequencyPlan.__init__",
        lambda *a, **k: pytest.fail("reselected"),
    )
    loaded = EPWRestartArchive.load(manifest)
    object.__setattr__(loaded, "_baseline_sha256", "0" * 64)
    with pytest.raises(ValueError, match="baseline SHA256"):
        loaded.reconstruct_native()
    object.__setattr__(loaded, "_baseline_sha256", archive._baseline_sha256)
    object.__setattr__(loaded, "_target_sha256", "0" * 64)
    with pytest.raises(ValueError, match="target native"):
        loaded.reconstruct_native()
    component = manifest.with_suffix(".npz")
    component.write_bytes(component.read_bytes() + b"x")
    with pytest.raises(ValueError, match="hash"):
        EPWRestartArchive.load(manifest)


def test_bad_dtype_and_path(tmp_path):
    from kernelgecp.integrations.epw import _digest, _json_object

    _, archive = _archive()
    manifest = tmp_path / "a.json"
    archive.save(manifest)
    component = manifest.with_suffix(".npz")
    meta, arrays = read_archive(component)
    arrays["corrections"] = arrays["corrections"].astype("float64")
    write_archive(component, meta, arrays)
    record = json.loads(manifest.read_text())
    record["archive"]["sha256"] = _digest(component.read_bytes())
    record.pop("record_sha256")
    record["record_sha256"] = _digest(_json_object(record, "record").encode())
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="dtype"):
        EPWRestartArchive.load(manifest)
    record["archive"]["file"] = "../a.npz"
    record.pop("record_sha256")
    record["record_sha256"] = _digest(_json_object(record, "record").encode())
    manifest.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="local filename"):
        EPWRestartArchive.load(manifest)


def test_native_frequency_recovery_matches_independent_matrix_solve():
    data = _static_data()
    native = data.frequencies_ev * (1 + 5e-11)
    fit0 = _fit_static(data)
    fit = EPWFrequencyCrossFit(
        fit0.plan,
        fit0.coefficients,
        frequencies_ev=native,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=fit0.fit_rows,
        diagnostics=fit0.diagnostics,
    )
    green = fit.evaluate_matsubara(data.indices)
    inverse = np.linalg.solve(green, np.broadcast_to(np.eye(2), green.shape))
    z = ((inverse[..., 0, 0] + inverse[..., 1, 1]) / (2j * native[:, None])).real
    gap = (-(inverse[..., 0, 1] + inverse[..., 1, 0]) / (2 * z)).real
    chi = ((inverse[..., 1, 1] - inverse[..., 0, 0]) / 2).real
    chi = (
        chi
        - data.state_energies_ev[None, :]
        + data.context["chemical_potential_offset_ev"]
    )
    baseline = _baseline(fit)
    np.testing.assert_allclose(baseline[..., 0], z, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(baseline[..., 1], gap, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(baseline[..., 2], chi, rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize("value", [1e-15, 1e-200, 1e200])
def test_gap_sensitivity_rejects_unsafe_channels(value):
    data = _static_data()
    fit = _fit_static(data)
    coefficients = fit.coefficients.copy()
    coefficients[..., 1] = value
    bad = EPWFrequencyCrossFit(
        fit.plan,
        coefficients,
        frequencies_ev=data.frequencies_ev,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=fit.fit_rows,
        diagnostics=fit.diagnostics,
    )
    with pytest.raises(NumericalBreakdownError):
        EPWRestartArchive.from_data(data, bad, policy="lossless_native_xor_v1")


def test_save_failure_preserves_existing_pair_and_cleans_new_files(
    tmp_path, monkeypatch
):
    import kernelgecp.integrations.epw_restart_archive as module

    _, archive = _archive()
    existing = tmp_path / "existing.json"
    archive.save(existing)
    before = (existing.read_bytes(), existing.with_suffix(".npz").read_bytes())
    with pytest.raises(ValueError, match="overwrite"):
        archive.save(existing, overwrite=True)
    with pytest.raises(FileExistsError):
        archive.save(existing)
    assert before == (existing.read_bytes(), existing.with_suffix(".npz").read_bytes())
    assert EPWRestartArchive.load(existing).reconstruct_native().z.size
    target = tmp_path / "new.json"
    link = module.os.link
    calls = 0

    def fail_second(*args):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected manifest publication failure")
        return link(*args)

    monkeypatch.setattr(module.os, "link", fail_second)
    with pytest.raises(OSError, match="injected"):
        archive.save(target)
    assert not target.exists() and not target.with_suffix(".npz").exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "existing.json",
        "existing.npz",
    ]
    assert before == (existing.read_bytes(), existing.with_suffix(".npz").read_bytes())


def test_consumed_17k_native_printed_roundtrip(tmp_path):
    import gzip
    from pathlib import Path

    from kernelgecp.integrations.epw_frequency import fit_epw_frequency_cross

    fixture = Path(__file__).parent / "data" / "epw"
    entry = json.loads((fixture / "manifest-scout-wscut1.json").read_text())["files"][0]
    path = tmp_path / "source"
    path.write_bytes(gzip.decompress((fixture / entry["compressed_file"]).read_bytes()))
    data = read_epw_anisotropic(path, entry["context"])
    rows = np.array(
        [i for i in range(len(data.indices)) if i not in {36, 40, 44, 45, 49, 53}]
    )
    fit = fit_epw_frequency_cross(
        data,
        real_nambu_parity=True,
        fit_rows=rows,
        rank=60,
        pole_cutoff_ev=13,
        precision_bits=256,
        rcond=1e-12,
    )
    archive = EPWRestartArchive.from_data(data, fit, policy="lossless_native_xor_v1")
    archive.save(tmp_path / "full.json")
    restored = EPWRestartArchive.load(tmp_path / "full.json")
    reconstructed = restored.reconstruct_native()
    for field in ("z", "delta_ev", "chi_ev"):
        np.testing.assert_array_equal(
            getattr(data, field).view("u8"), getattr(reconstructed, field).view("u8")
        )
    context = restored.write_native_restart(tmp_path / "written")
    parsed = read_epw_anisotropic(tmp_path / "written", context)
    for field in (
        "frequencies_ev",
        "state_energies_ev",
        "z",
        "delta_ev",
        "chi_ev",
        "indices",
    ):
        np.testing.assert_array_equal(getattr(data, field), getattr(parsed, field))
    assert parsed.context["state_order"] == data.context["state_order"]
    assert parsed.context["matsubara_indices"] == data.context["matsubara_indices"]


def test_native_text_rounding_is_not_binary64_roundtrip(tmp_path):
    data = _static_data()
    gap = data.delta_ev.copy()
    gap[0, 0] = np.nextafter(gap[0, 0], np.inf)
    changed = EPWAnisotropicData(
        data.frequencies_ev,
        data.state_energies_ev,
        data.z,
        gap,
        data.chi_ev,
        context=data.context,
    )
    archive = EPWRestartArchive.from_data(
        changed, _fit_static(changed), policy="lossless_native_xor_v1"
    )
    assert archive.reconstruct_native().delta_ev[0, 0].view("u8") == gap[0, 0].view(
        "u8"
    )
    context = archive.write_native_restart(tmp_path / "native")
    parsed = read_epw_anisotropic(tmp_path / "native", context)
    assert parsed.delta_ev[0, 0] == float(f"{gap[0, 0]:.10E}")
    assert parsed.delta_ev[0, 0].view("u8") != gap[0, 0].view("u8")
