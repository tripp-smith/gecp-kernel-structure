"""Pure API/uncertainty checks and explicit optional actual native execution."""

from __future__ import annotations

import builtins
import dataclasses
import hashlib
import importlib.util
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from experiments.cthyb_dmft import (
    ARRAYS,
    _read,
    _sample,
    load_config,
    load_measurement,
    save_measurement,
)
from kernelgecp.integrations.cthyb_dmft import (
    CTHYBConfig,
    CTHYBMeasurement,
    CTHYBPlan,
    estimate_replicas,
)


@pytest.mark.parametrize(
    "change",
    [
        {"beta": True},
        {"hopping": 0},
        {"interaction": -1},
        {"n_iw": 7},
        {"n_tau": 1536},
        {"warmup_cycles": 0},
        {"length_cycle": 1.5},
        {"mixing": 1.01},
        {"rcond": float("nan")},
        {"tail_fraction": 1},
        {"n_iw": 8192},
        {"n_tau": 100000},
        {"beta": float("inf")},
    ],
)
def test_config_rejects_invalid(change: dict[str, object]) -> None:
    with pytest.raises((ValueError, TypeError)):
        CTHYBConfig(**change)


def test_import_is_lazy(monkeypatch: pytest.MonkeyPatch) -> None:
    original = builtins.__import__

    def guarded(name: str, *args: object, **kwargs: object) -> object:
        if name in ("triqs", "triqs_cthyb") and not kwargs.get("level", 0):
            raise AssertionError("optional native runtime imported eagerly")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    import kernelgecp.integrations.cthyb_dmft

    assert kernelgecp.integrations.cthyb_dmft.CTHYBConfig().beta == 10


def test_replica_estimator_independent_closed_form() -> None:
    values = np.array([[1, 4], [2, 2], [3, 6], [4, 8]], dtype=float)
    result = estimate_replicas(values)
    # Independent scalar two-pass sample covariance calculation.
    means = [sum(values[:, i]) / 4 for i in range(2)]
    covariance = np.array(
        [
            [
                sum((row[i] - means[i]) * (row[j] - means[j]) for row in values) / 3
                for j in range(2)
            ]
            for i in range(2)
        ]
    )
    np.testing.assert_allclose(result.mean, means)
    np.testing.assert_allclose(result.covariance, covariance)
    np.testing.assert_allclose(result.standard_error, np.sqrt(np.diag(covariance) / 4))
    values[:] = 0
    assert result.mean[0] == 2.5
    for array in (result.mean, result.covariance, result.standard_error):
        with pytest.raises(ValueError):
            array.setflags(write=True)


@pytest.mark.parametrize(
    "values",
    [
        [[1.0, 2.0]],
        [[1.0, float("nan")], [2.0, 3.0]],
        [1, 2],
        [[1j, 2], [2, 3]],
        [[True], [False]],
    ],
)
def test_replica_estimator_rejects_bad(values: object) -> None:
    with pytest.raises(ValueError):
        estimate_replicas(values)


def fixture_measurement() -> CTHYBMeasurement:
    config = CTHYBConfig(n_iw=8, n_tau=97)
    moments = np.array([[0, 0], [1, 1], [0, 0]], dtype=complex)
    return CTHYBMeasurement(
        config,
        123,
        100,
        np.zeros(16, complex),
        np.zeros((97, 2), complex),
        np.zeros((16, 2), complex),
        np.zeros((16, 2), complex),
        np.zeros((16, 2), complex),
        moments,
        np.array([[1, 1], [1, 1]], complex),
        np.array([0.5, 0.5]),
        0.2,
        1.0,
        0,
    )


def test_measurement_owns_arrays() -> None:
    measurement = fixture_measurement()
    for name in ARRAYS:
        with pytest.raises(ValueError):
            getattr(measurement, name).setflags(write=True)
    assert measurement.density == 1
    with pytest.raises(ValueError):
        dataclasses.replace(measurement, cycles=True)
    with pytest.raises(ValueError):
        dataclasses.replace(measurement, seed=2**31)
    with pytest.raises(ValueError):
        dataclasses.replace(measurement, solve_status=1)
    with pytest.raises(ValueError):
        dataclasses.replace(measurement, g_tau=np.full((97, 2), np.nan))
    with pytest.raises(ValueError):
        dataclasses.replace(measurement, g_moments=np.zeros((3, 2)))


def test_archive_roundtrip_and_mutation_rejection(tmp_path: Path) -> None:
    measurement = fixture_measurement()
    identity = {"source": "test fixture, not provider evidence"}
    stem = tmp_path / "sample"
    save_measurement(stem, measurement, identity)
    result = load_measurement(stem, identity)
    for name in ARRAYS:
        np.testing.assert_array_equal(getattr(result, name), getattr(measurement, name))
    with pytest.raises(FileExistsError):
        save_measurement(stem, measurement, identity)
    with pytest.raises(ValueError, match="identity"):
        load_measurement(stem, {})
    stem.with_suffix(".npz").write_bytes(
        stem.with_suffix(".npz").read_bytes() + b"corruption"
    )
    with pytest.raises(ValueError, match="hash"):
        load_measurement(stem, identity)


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":1e999}', "[]"])
def test_receipt_json_rejects_bad_numbers_and_duplicate_keys(
    tmp_path: Path, text: str
) -> None:
    path = tmp_path / "test.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        _read(path)


def test_canonical_config_is_not_weakened(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    original = (root / "experiments/cthyb_dmft.toml").read_text()
    data = load_config(root / "experiments/cthyb_dmft.toml")
    assert data["cycle_budgets"] == [100000, 400000, 1600000]
    assert data["final_replicas"] >= 8
    for before, after in [
        ("max_standard_error = 0.0003", "max_standard_error = 0.003"),
        ("agreement_absolute = 0.001", "agreement_absolute = 0.01"),
        ("400000, 1600000", "400000, 3200000"),
    ]:
        path = tmp_path / "bad.toml"
        path.write_text(original.replace(before, after))
        with pytest.raises(ValueError, match="weakened"):
            load_config(path)


native = pytest.mark.skipif(
    importlib.util.find_spec("triqs_cthyb") is None,
    reason="optional actual CT-HYB runtime; mandatory smoke never skips",
)


@pytest.fixture(scope="module")
def actual() -> tuple[CTHYBPlan, CTHYBMeasurement]:
    from kernelgecp import FermionicBasis

    config = CTHYBConfig(n_iw=32, n_tau=385, warmup_cycles=500)
    basis = FermionicBasis(10, 8, 1e-6, grid_order=24, max_rank=32)
    plan = CTHYBPlan(config, basis, fit_tau_indices=np.arange(0, config.n_tau, 2))
    return plan, plan.sample(plan.initial_bath(), cycles=2000, seed=27199921)


@native
def test_actual_measurement_has_density_matrix_moments(
    actual: tuple[CTHYBPlan, CTHYBMeasurement],
) -> None:
    plan, m = actual
    assert m.solve_status == 0 and m.average_sign > 0.99
    assert 0.7 < m.density < 1.3 and 0 < m.double_occupancy < 0.4
    # Native commutator moments use the opposite spin's measured occupation.
    np.testing.assert_allclose(m.sigma_moments[0], 2 * m.occupancies[::-1], atol=1e-12)
    np.testing.assert_allclose(m.g_moments[2], m.sigma_moments[0] - 1, atol=1e-12)
    np.testing.assert_allclose(m.g_iw, m.g_iw[::-1].conj(), atol=1e-12)
    np.testing.assert_allclose(plan.indices, np.arange(-32, 32))


@native
def test_actual_measured_feedback_and_unobserved_rows(
    actual: tuple[CTHYBPlan, CTHYBMeasurement],
) -> None:
    plan, m = actual
    rep = plan.represent(m)
    bath, defect = plan.bath_update(m, rep)
    expected = m.bath_iw + 0.5 * (rep.g_iw.mean(axis=1) - m.bath_iw)
    np.testing.assert_allclose(bath, expected, atol=1e-14)
    assert defect == pytest.approx(np.max(np.abs(rep.g_iw.mean(axis=1) - m.bath_iw)))
    changed = m.g_tau.copy()
    changed[1::2] += 100
    modified = dataclasses.replace(m, g_tau=changed)
    second = plan.represent(modified)
    np.testing.assert_array_equal(rep.coefficients, second.coefficients)
    assert second.validation_residual > 99
    for array in (rep.g_tau, rep.g_iw, rep.coefficients, bath, plan.indices, plan.tau):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(ValueError):
        plan.sample(m.bath_iw, cycles=0, seed=3)
    with pytest.raises(ValueError):
        plan.observables(m, [32])
    with pytest.raises(AttributeError):
        plan.config = CTHYBConfig()


@native
def test_actual_fresh_chain_repeatability(
    actual: tuple[CTHYBPlan, CTHYBMeasurement],
) -> None:
    plan, m = actual
    repeated = plan.sample(m.bath_iw, cycles=m.cycles, seed=m.seed)
    np.testing.assert_array_equal(repeated.g_tau, m.g_tau)
    assert repeated.double_occupancy == m.double_occupancy


def test_symlink_archive_is_rejected(tmp_path: Path) -> None:
    measurement = fixture_measurement()
    save_measurement(tmp_path / "original", measurement, {})
    (tmp_path / "alias.npz").symlink_to(tmp_path / "original.npz")
    (tmp_path / "alias.json").symlink_to(tmp_path / "original.json")
    with pytest.raises(ValueError, match="symlink"):
        load_measurement(tmp_path / "alias", {})


@pytest.mark.parametrize(
    "mode", ["object", "huge_header", "wrong_shape", "missing_member"]
)
def test_array_headers_are_validated_before_allocation(
    tmp_path: Path, mode: str
) -> None:
    m = fixture_measurement()
    stem = tmp_path / "sample"
    save_measurement(stem, m, {})
    contents = {name + ".npy": None for name in ARRAYS}
    with zipfile.ZipFile(stem.with_suffix(".npz")) as archive:
        contents = {name: archive.read(name) for name in contents}
    if mode == "missing_member":
        del contents["g_tau.npy"]
    else:
        stream = io.BytesIO()
        if mode == "huge_header":
            np.lib.format.write_array_header_1_0(
                stream, {"descr": "<c16", "fortran_order": False, "shape": (2**50, 2)}
            )
        else:
            values = (
                np.array([object()], dtype=object)
                if mode == "object"
                else np.zeros((1, 2), complex)
            )
            np.lib.format.write_array(stream, values, allow_pickle=True)
        contents["g_tau.npy"] = stream.getvalue()
    with zipfile.ZipFile(stem.with_suffix(".npz"), "w") as archive:
        for name, value in contents.items():
            archive.writestr(name, value)
    receipt = json.loads(stem.with_suffix(".json").read_text())
    receipt["arrays_sha256"] = hashlib.sha256(
        stem.with_suffix(".npz").read_bytes()
    ).hexdigest()
    stem.with_suffix(".json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError):
        load_measurement(stem, {})


def test_resume_checks_seed_cycles_and_actual_bath(tmp_path: Path) -> None:
    m = fixture_measurement()
    stem = tmp_path / "sample"
    save_measurement(stem, m, {})

    class Plan:
        config = m.config

    for bath, cycles, seed in [
        (m.bath_iw + 1, 100, 123),
        (m.bath_iw, 101, 123),
        (m.bath_iw, 100, 124),
    ]:
        with pytest.raises(ValueError, match="requested input"):
            _sample(Plan(), bath, cycles, seed, stem, {})
    result = _sample(Plan(), m.bath_iw, 100, 123, stem, {})
    assert result.seed == 123
