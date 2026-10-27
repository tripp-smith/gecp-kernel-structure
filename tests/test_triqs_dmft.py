"""Actual native IPT reference checks; skipped providers are not acceptance."""

from __future__ import annotations

import importlib
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.special import expit

from kernelgecp.dmft import IPTConfig
from kernelgecp.integrations import triqs_dmft as native
from kernelgecp.lehmann import (
    FermionicPoleModel,
    hubbard_atom,
    impurity_one_bath,
    ipt_self_energy_poles,
)


@pytest.fixture(scope="module")
def provider() -> tuple[Any, Any]:
    gfs = pytest.importorskip("triqs.gfs", exc_type=ModuleNotFoundError)
    assert importlib.import_module("triqs.version").version == "4.0.2"
    return gfs, importlib.import_module("triqs.mesh")


def test_native_reference_import_is_lazy() -> None:
    code = """
import builtins,sys
original=builtins.__import__
def guard(name,globals=None,locals=None,fromlist=(),level=0):
    if level==0 and (name=='triqs' or name.startswith('triqs.')):
        raise RuntimeError('eager optional import')
    return original(name,globals,locals,fromlist,level)
builtins.__import__=guard
import kernelgecp.integrations.triqs_dmft
assert not any(n=='triqs' or n.startswith('triqs.') for n in sys.modules)
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_wrong_config_or_missing_provider_never_falls_back(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(TypeError):
        native.DenseIPTPlan(object())

    def absent() -> Any:
        raise ImportError("deliberately absent actual TRIQS")

    monkeypatch.setattr(native, "_provider", absent)
    with pytest.raises(ImportError, match="actual TRIQS"):
        native.solve_dense_ipt(IPTConfig())


@pytest.mark.parametrize("bad", [[True, False], [0, np.nan], [[0, 1]], ["0", "1"]])
def test_invalid_bath_vector_rejected_without_native_calls(bad: Any) -> None:
    with pytest.raises(ValueError):
        native._complex_vector(bad, 2, "bath")


@pytest.mark.parametrize("interaction,iterations", [(0, 2), (1, 19), (2, 21)])
def test_canonical_native_loops_and_unmixed_equations(
    provider: tuple[Any, Any], interaction: float, iterations: int
) -> None:
    config = IPTConfig(interaction=interaction)
    result = native.solve_dense_ipt(config, record_trajectory=True)
    assert result.converged and result.stop_reason == "tolerance"
    assert len(result.history) == iterations
    assert len(result.trajectory) == iterations
    assert all(row.unmixed_defect < config.tolerance for row in result.history[-2:])
    assert abs(result.density - 0.5) < 1e-10
    assert result.native_tail.positive_frequency_imag_g_max <= 0
    assert result.native_tail.positive_frequency_imag_sigma_max <= 0
    assert result.native_tail.particle_hole_g_error < 1e-10
    iw = 1j * (2 * result.indices + 1) * np.pi / config.beta
    np.testing.assert_allclose(
        result.g_iw, 1 / (iw - result.bath_iw - result.sigma_iw), atol=1e-14
    )
    np.testing.assert_allclose(result.weiss_iw, 1 / (iw - result.bath_iw), atol=1e-14)
    np.testing.assert_allclose(
        result.sigma_tau,
        config.interaction**2 * result.weiss_tau**2 * result.weiss_tau[::-1],
        atol=1e-14,
    )
    assert result.history[-1].unmixed_defect == pytest.approx(
        np.max(np.abs(config.hopping**2 * result.g_iw - result.bath_iw))
    )
    assert result.history[-1].mixed_step == pytest.approx(
        config.mixing * result.history[-1].unmixed_defect
    )
    if interaction:
        assert result.native_tail.sigma_first_moment_error < 1e-6
    else:
        np.testing.assert_array_equal(result.sigma_iw, 0)
        exact = (
            -2j
            * np.sign(iw.imag)
            / (np.abs(iw.imag) + np.hypot(iw.imag, 2 * config.hopping))
        )
        np.testing.assert_allclose(result.g_iw, exact, atol=1e-14)


def test_one_small_defect_is_not_convergence_and_budget_is_retained(
    provider: tuple[Any, Any],
) -> None:
    config = IPTConfig(interaction=0, n_iw=32, n_tau=385, max_iterations=1)
    result = native.solve_dense_ipt(config, record_trajectory=True)
    assert result.history[0].unmixed_defect < config.tolerance
    assert not result.converged and result.stop_reason == "max_iterations"
    assert len(result.history) == len(result.trajectory) == 1
    interacting = native.solve_dense_ipt(replace(config, interaction=2))
    assert not interacting.converged
    assert interacting.history[-1].unmixed_defect > config.tolerance
    assert interacting.trajectory == ()


def test_plan_and_result_own_immutable_state(provider: tuple[Any, Any]) -> None:
    plan = native.DenseIPTPlan(IPTConfig(n_iw=32, n_tau=385))
    bath = plan.initial_bath()
    first = plan.update(bath)
    expected = first.g_iw.copy()
    bath[:] = 0
    plan.update(bath)
    np.testing.assert_array_equal(first.g_iw, expected)
    for array in (
        first.g_iw,
        first.sigma_iw,
        first.weiss_tau,
        first.indices,
        first.tau,
    ):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(AttributeError):
        plan.config = IPTConfig()
    with pytest.raises(AttributeError):
        plan.indices = np.arange(2)
    result = plan.solve(record_trajectory=True)
    with pytest.raises(ValueError):
        result.trajectory[0].setflags(write=True)
    with pytest.raises(TypeError):
        plan.solve(record_trajectory=1)


def test_native_refinement_disjoint_time_and_unseen_frequency(
    provider: tuple[Any, Any],
) -> None:
    coarse = native.solve_dense_ipt(IPTConfig())
    fine = native.solve_dense_ipt(IPTConfig(n_iw=1024, n_tau=12289))
    assert np.max(np.abs(coarse.g_iw - fine.g_iw[512:1536])) < 1e-7
    assert np.max(np.abs(coarse.sigma_iw - fine.sigma_iw[512:1536])) < 1e-7
    assert np.max(np.abs(coarse.g_tau - fine.g_tau[::2])) < 1e-7
    assert np.max(np.abs(coarse.sigma_tau - fine.sigma_tau[::2])) < 1e-7
    heldout = np.arange(3, fine.tau.size - 1, 4)
    for field, values in (("green", fine.g_tau), ("sigma", fine.sigma_tau)):
        assert (
            np.max(
                np.abs(
                    coarse.evaluate_tau(fine.tau[heldout], field=field)
                    - values[heldout]
                )
            )
            < 1e-7
        )
    indices = np.array([-701, -601, -513, 512, 600, 700])
    assert not np.any(np.isin(indices, coarse.indices))
    for field in ("green", "sigma"):
        assert (
            np.max(
                np.abs(
                    coarse.evaluate_matsubara(indices, field=field)
                    - fine.evaluate_matsubara(indices, field=field)
                )
            )
            < 1e-7
        )
    with pytest.raises(ValueError, match="refine"):
        coarse.evaluate_matsubara([2000])
    nodes, weights = np.polynomial.legendre.leggauss(256)
    tau = coarse.config.beta * (nodes + 1) / 2
    labels = np.array([-8, -2, 0, 6])
    iw = 1j * (2 * labels + 1) * np.pi / coarse.config.beta
    for field in ("green", "sigma"):
        integral = (
            np.exp(iw[:, None] * tau) * (coarse.config.beta * weights / 2)
        ) @ coarse.evaluate_tau(tau, field=field)
        assert (
            np.max(np.abs(integral - coarse.evaluate_matsubara(labels, field=field)))
            < 1e-7
        )


def _native_pole_update(
    provider: tuple[Any, Any], model: FermionicPoleModel, interaction: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gfs, meshes = provider
    indices = np.arange(-512, 512)
    mesh_iw = meshes.MeshImFreq(beta=model.beta, statistic="Fermion", n_iw=512)
    mesh_tau = meshes.MeshImTime(beta=model.beta, statistic="Fermion", n_tau=6145)
    giw = gfs.Gf(mesh=mesh_iw.copy(), target_shape=[])
    gtau = gfs.Gf(mesh=mesh_tau.copy(), target_shape=[])
    siw = gfs.Gf(mesh=mesh_iw.copy(), target_shape=[])
    stau = gfs.Gf(mesh=mesh_tau.copy(), target_shape=[])
    giw.data[:] = model.evaluate_matsubara(indices)
    gtau.set_from_fourier(giw, known_moments=np.r_[0, model.moments(2)].astype(complex))
    stau.data[:] = interaction**2 * gtau.data**2 * gtau.data[::-1]
    density = model.density()
    left, right = -(1 - density), -density
    derivative0 = np.sum(
        model.weights * model.energies * expit(model.beta * model.energies)
    )
    derivativeb = np.sum(
        model.weights * model.energies * expit(-model.beta * model.energies)
    )
    first = interaction**2 * density * (1 - density)
    second = interaction**2 * (
        2 * left * derivative0 * right
        - left**2 * derivativeb
        + 2 * right * derivativeb * left
        - right**2 * derivative0
    )
    siw.set_from_fourier(stau, known_moments=np.array([0, first, second], complex))
    return np.array([p.value for p in mesh_tau]), stau.data.copy(), siw.data.copy()


@pytest.mark.parametrize("asymmetric", [False, True])
def test_native_reflected_product_against_independent_pole_oracle(
    provider: tuple[Any, Any], asymmetric: bool
) -> None:
    model = (
        FermionicPoleModel(10, [-1.3, 0.2, 1.1], [0.2, 0.5, 0.3])
        if asymmetric
        else FermionicPoleModel(10, [-1.7, -0.5, 0.5, 1.7], [0.2, 0.3, 0.3, 0.2])
    )
    oracle = ipt_self_energy_poles(model, 2)
    tau, native_tau, native_iw = _native_pole_update(provider, model, 2)
    assert np.max(np.abs(native_tau - oracle.evaluate_tau(tau))) < 1e-8
    assert (
        np.max(np.abs(native_iw - oracle.evaluate_matsubara(np.arange(-512, 512))))
        < 1e-8
    )
    if asymmetric:
        wrong_cube = 4 * model.evaluate_tau(tau) ** 3
        assert np.max(np.abs(wrong_cube - oracle.evaluate_tau(tau))) > 0.1


def test_atomic_dyson_and_finite_bath_native_fourier_match_lehmann(
    provider: tuple[Any, Any],
) -> None:
    tau, _, sigma = _native_pole_update(provider, FermionicPoleModel(10, [0], [1]), 2)
    indices = np.arange(-512, 512)
    iw = 1j * (2 * indices + 1) * np.pi / 10
    assert (
        np.max(
            np.abs(1 / (iw - sigma) - hubbard_atom(10, 2).evaluate_matsubara(indices))
        )
        < 1e-10
    )
    gfs, meshes = provider
    model = impurity_one_bath(10, 2, 0.4)
    giw = gfs.Gf(
        mesh=meshes.MeshImFreq(beta=10, statistic="Fermion", n_iw=512), target_shape=[]
    )
    gtau = gfs.Gf(
        mesh=meshes.MeshImTime(beta=10, statistic="Fermion", n_tau=tau.size),
        target_shape=[],
    )
    giw.data[:] = model.evaluate_matsubara(indices)
    gtau.set_from_fourier(giw, known_moments=np.r_[0, model.moments(2)].astype(complex))
    assert np.max(np.abs(gtau.data - model.evaluate_tau(tau))) < 1e-8


@pytest.fixture
def smoke_script(monkeypatch: pytest.MonkeyPatch) -> Any:
    directory = Path(__file__).resolve().parents[1] / "scripts" / "providers"
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location(
        "dmft_smoke_test", directory / "dmft_smoke.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in module.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    return module


def test_smoke_rejects_stale_revision_before_any_import(
    smoke_script: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        smoke_script, "source_identity", lambda root: {"source_revision": "a" * 40}
    )
    result = smoke_script.run_probe("b" * 40)
    assert result["status"] == "failed"
    assert result["checks"]["source_revision"] is False
    assert result["checks"]["imports"] is False
    assert result["checks"]["numerical"] is False


def test_smoke_rejects_absent_actual_runtime_without_skipping(
    smoke_script: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        smoke_script, "source_identity", lambda root: {"source_revision": "a" * 40}
    )
    monkeypatch.setattr(smoke_script, "environment_identity", lambda prefix: {})
    original = importlib.import_module

    def missing(name: str) -> Any:
        if name.startswith("triqs"):
            raise ModuleNotFoundError("deliberately missing native provider")
        return original(name)

    monkeypatch.setattr(smoke_script.importlib, "import_module", missing)
    result = smoke_script.run_probe("a" * 40)
    assert result["status"] == "failed"
    assert result["error"]["stage"] == "imports"
    assert result["checks"]["numerical"] is False


def test_extended_wheel_guard_rejects_a_shadowed_native_reference(
    smoke_script: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    paths = {
        name: tmp_path / relative
        for name, relative in smoke_script.PACKAGE_FILES.items()
    }

    class Distribution:
        files = list(smoke_script.PACKAGE_FILES.values())

        def read_text(self, name: str) -> str:
            return "\n".join(relative + ",," for relative in self.files)

        def locate_file(self, relative: str) -> Path:
            return tmp_path / relative

    monkeypatch.setattr(smoke_script, "installed_wheel_identity", lambda value: {})
    monkeypatch.setattr(
        smoke_script.metadata, "distribution", lambda name: Distribution()
    )
    paths["reference"] = tmp_path / "shadow.py"
    with pytest.raises(RuntimeError, match="not owned by the wheel"):
        smoke_script._wheel_identity(paths)


def test_failed_smoke_cli_remains_nonzero_with_a_failed_json_receipt(
    smoke_script: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    receipt = {"status": "failed", "checks": {"imports": False, "numerical": False}}
    monkeypatch.setattr(smoke_script, "run_probe", lambda *args, **kwargs: receipt)
    path = tmp_path / "receipt.json"
    monkeypatch.setattr(
        smoke_script.sys,
        "argv",
        [
            "dmft_smoke",
            "--output",
            str(path),
            "--source-revision",
            "a" * 40,
        ],
    )
    assert smoke_script.main() == 1
    assert json.loads(path.read_text()) == receipt
