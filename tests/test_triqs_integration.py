"""Optional actual-provider interoperability plus always-run boundary checks.

The fail-closed scripts/providers/triqs_adapter_smoke.py supplies mandatory
provider acceptance; skips in the base environment are not that acceptance.
"""

from __future__ import annotations

import importlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from kernelgecp import FermionicBasis
from kernelgecp.integrations import triqs as adapter


@pytest.fixture(scope="module")
def provider() -> tuple[Any, Any]:
    gfs = pytest.importorskip("triqs.gfs", exc_type=ModuleNotFoundError)
    assert importlib.import_module("triqs.version").version == "4.0.2"
    return gfs, importlib.import_module("triqs.mesh")


@pytest.fixture(scope="module")
def basis() -> FermionicBasis:
    return FermionicBasis(
        beta=4, wmax=2, tol=1e-12, precision="reference", grid_order=8, max_rank=64
    )


def _time_values(tau: np.ndarray, residues: np.ndarray) -> np.ndarray:
    energies = np.array([-0.8, 0.6])
    scalar = -np.exp(-tau[:, None] * energies) / (1 + np.exp(-4 * energies))
    return np.tensordot(scalar, residues, axes=(1, 0))


def _iw_values(indices: np.ndarray, residues: np.ndarray) -> np.ndarray:
    energies = np.array([-0.8, 0.6])
    iw = 1j * (2 * indices + 1) * np.pi / 4
    return np.tensordot(1 / (iw[:, None] - energies), residues, axes=(1, 0))


def _residues() -> np.ndarray:
    first = np.array([[0.7, 0.15 + 0.2j], [0.15 - 0.2j, 0.3]])
    return np.array([first, np.eye(2) - first])


def _source(provider: tuple[Any, Any], residues: np.ndarray) -> Any:
    gfs, meshes = provider
    mesh = meshes.MeshImTime(beta=4, statistic="Fermion", n_tau=97)
    source = gfs.Gf(mesh=mesh, target_shape=list(residues.shape[1:]))
    source.data[...] = _time_values(np.array([p.value for p in mesh]), residues)
    return source


def test_importing_integration_never_imports_optional_provider() -> None:
    code = """
import builtins
import sys
original = builtins.__import__
def guarded(name, *args, **kwargs):
    if name == 'triqs' or name.startswith('triqs.'):
        raise RuntimeError('optional provider imported eagerly')
    return original(name, *args, **kwargs)
builtins.__import__ = guarded
import kernelgecp
import kernelgecp.integrations
import kernelgecp.integrations.triqs
assert not any(n == 'triqs' or n.startswith('triqs.') for n in sys.modules)
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        env={**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1])},
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_missing_or_wrong_provider_fails_clearly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def absent(name: str) -> Any:
        raise ModuleNotFoundError(name)

    monkeypatch.setattr(adapter.importlib, "import_module", absent)
    with pytest.raises(ImportError, match="TRIQS 4.0.2"):
        adapter._provider()
    monkeypatch.setattr(
        adapter.importlib,
        "import_module",
        lambda name: SimpleNamespace(version="3.3.3"),
    )
    with pytest.raises(ImportError, match="pinned TRIQS 4.0.2"):
        adapter._provider()


@pytest.mark.parametrize("matrix", [False, True])
def test_scalar_and_complex_matrix_physical_fits(
    provider: tuple[Any, Any], basis: FermionicBasis, matrix: bool
) -> None:
    _, meshes = provider
    residues = _residues() if matrix else np.array([0.4, 0.6])
    source = _source(provider, residues)
    fitted = adapter.fit_triqs_tau(basis, source)
    query_mesh = meshes.MeshImTime(beta=4, statistic="Fermion", n_tau=98)
    queried = adapter.to_triqs_tau(fitted, query_mesh)
    tau = np.array([point.value for point in query_mesh])
    assert queried.data.shape == (98, *residues.shape[1:])
    assert np.max(np.abs(queried.data - _time_values(tau, residues))) < 1e-8
    assert not np.allclose(queried.data[0], queried.data[-1])
    normalization = np.eye(2) if matrix else 1.0
    np.testing.assert_allclose(
        queried.data[0] + queried.data[-1], -normalization, atol=1e-8
    )
    for positive in (False, True):
        mesh = meshes.MeshImFreq(
            beta=4,
            statistic="Fermion",
            n_iw=19,
            opt="positive_frequencies_only" if positive else "all_frequencies",
        )
        output = adapter.to_triqs_iw(fitted, mesh)
        indices = np.array([point.index for point in mesh])
        assert output.data.shape == (19 if positive else 38, *residues.shape[1:])
        assert np.max(np.abs(output.data - _iw_values(indices, residues))) < 1e-8
        assert indices[0] == (0 if positive else -19)
        if matrix and not positive:
            np.testing.assert_allclose(
                output.data[18], output.data[19].conj().T, atol=1e-8
            )
            assert np.max(np.abs(output.data[18] - output.data[19].conj())) > 1e-3


def test_independent_fourier_sign_and_complex_linearity(
    provider: tuple[Any, Any], basis: FermionicBasis
) -> None:
    _, meshes = provider
    residues = np.array([[[1 + 2j, 3 + 4j], [-0.6j, 2 - 0.5j]], np.zeros((2, 2))])
    fitted = adapter.fit_triqs_tau(basis, _source(provider, residues))
    mesh = meshes.MeshImFreq(beta=4, statistic="Fermion", n_iw=4)
    output = adapter.to_triqs_iw(fitted, mesh)
    indices = np.array([point.index for point in mesh])
    assert np.max(np.abs(output.data - _iw_values(indices, residues))) < 1e-8
    nodes, weights = np.polynomial.legendre.leggauss(96)
    tau = 2 * (nodes + 1)
    iw = 1j * (2 * indices + 1) * np.pi / 4
    fourier = np.tensordot(
        np.exp(iw[:, None] * tau) * (2 * weights),
        fitted.evaluate_tau(tau),
        axes=(1, 0),
    )
    assert np.max(np.abs(output.data - fourier)) < 1e-8


def test_ownership_and_unusual_observation_axis(
    provider: tuple[Any, Any], basis: FermionicBasis
) -> None:
    _, meshes = provider
    source = _source(provider, _residues())
    before = source.data.copy()
    fitted = adapter.fit_triqs_tau(basis, source)
    tau = np.array([point.value for point in source.mesh])
    alternate = basis.fit_tau(np.moveaxis(before, 0, 1), tau, axis=1)
    mesh = meshes.MeshImTime(beta=4, statistic="Fermion", n_tau=98)
    output = adapter.to_triqs_tau(fitted, mesh)
    other = adapter.to_triqs_tau(alternate, mesh)
    assert other.target_shape == (2, 2)
    np.testing.assert_allclose(output.data, other.data, atol=1e-12)
    expected = output.data.copy()
    assert not np.shares_memory(source.data, output.data)
    assert not np.shares_memory(fitted.coefficients, output.data)
    source.data[...] = 0
    np.testing.assert_allclose(
        adapter.to_triqs_tau(fitted, mesh).data, expected, atol=0
    )
    output.data[...] = 9 + 4j
    np.testing.assert_allclose(
        adapter.to_triqs_tau(fitted, mesh).data, expected, atol=0
    )


@pytest.mark.parametrize("domain", ["tau", "iw"])
def test_export_owns_mutable_mesh_in_both_directions(
    provider: tuple[Any, Any], basis: FermionicBasis, domain: str
) -> None:
    """TRIQS Gf(mesh=mesh) aliases copy_from-capable meshes unless copied."""
    _, meshes = provider
    fitted = adapter.fit_triqs_tau(basis, _source(provider, _residues()))

    def make_mesh(beta: float, count: int) -> Any:
        if domain == "tau":
            return meshes.MeshImTime(beta=beta, statistic="Fermion", n_tau=count)
        return meshes.MeshImFreq(beta=beta, statistic="Fermion", n_iw=count)

    def coordinates(mesh: Any) -> np.ndarray:
        return np.array([complex(point.value) for point in mesh])

    original = make_mesh(4, 9)
    converter = adapter.to_triqs_tau if domain == "tau" else adapter.to_triqs_iw
    exported = converter(fitted, original)
    before = exported.data.copy()
    expected_coordinates = coordinates(original)
    expected_count = len(original)
    # Previously this changed the output to beta=5 and a different mesh size
    # while its data still contained the beta=4, original-sized evaluation.
    original.copy_from(make_mesh(5, 11))
    assert exported.mesh.beta == 4
    assert len(exported.mesh) == expected_count
    np.testing.assert_array_equal(coordinates(exported.mesh), expected_coordinates)
    np.testing.assert_array_equal(exported.data, before)
    caller_coordinates = coordinates(original)
    caller_count = len(original)
    # The owning output is mutable too, but its mutation cannot reach the
    # caller's mesh. Direct mesh mutation does not resize a Gf's data buffer.
    exported.mesh.copy_from(make_mesh(6, 13))
    assert original.beta == 5
    assert len(original) == caller_count
    np.testing.assert_array_equal(coordinates(original), caller_coordinates)
    np.testing.assert_array_equal(exported.data, before)


def test_ordered_unequal_block_recipe(
    provider: tuple[Any, Any], basis: FermionicBasis
) -> None:
    gfs, meshes = provider
    names = ["z_single", "a_matrix"]
    blocks = [
        _source(provider, np.array([0.4, 0.6])[:, None, None]),
        _source(provider, _residues()),
    ]
    original = gfs.BlockGf(name_list=names, block_list=blocks, make_copies=True)
    fits = [(name, adapter.fit_triqs_tau(basis, block)) for name, block in original]
    mesh = meshes.MeshImFreq(beta=4, statistic="Fermion", n_iw=7)
    output = gfs.BlockGf(
        name_list=[name for name, _ in fits],
        block_list=[adapter.to_triqs_iw(fit, mesh) for _, fit in fits],
        make_copies=True,
    )
    assert [name for name, _ in output] == names
    assert output[names[0]].data.shape == (14, 1, 1)
    assert output[names[1]].data.shape == (14, 2, 2)
    with pytest.raises(TypeError, match="BlockGf blocks in order"):
        adapter.fit_triqs_tau(basis, original)


@pytest.mark.parametrize("statistic,beta", [("Boson", 4), ("Fermion", 5)])
def test_actual_mesh_statistics_and_temperature_are_checked(
    provider: tuple[Any, Any], basis: FermionicBasis, statistic: str, beta: float
) -> None:
    gfs, meshes = provider
    tau = meshes.MeshImTime(beta=beta, statistic=statistic, n_tau=9)
    source = gfs.Gf(mesh=tau, target_shape=[])
    with pytest.raises(ValueError):
        adapter.fit_triqs_tau(basis, source)
    fitted = adapter.fit_triqs_tau(basis, _source(provider, np.array([0.4, 0.6])))
    with pytest.raises(ValueError):
        adapter.to_triqs_tau(fitted, tau)
    iw = meshes.MeshImFreq(beta=beta, statistic=statistic, n_iw=3)
    with pytest.raises(ValueError):
        adapter.to_triqs_iw(fitted, iw)


def test_wrong_mesh_empty_mesh_nonfinite_data_are_rejected(
    provider: tuple[Any, Any], basis: FermionicBasis
) -> None:
    _, meshes = provider
    source = _source(provider, np.array([0.4, 0.6]))
    fitted = adapter.fit_triqs_tau(basis, source)
    with pytest.raises(TypeError):
        adapter.to_triqs_iw(fitted, source.mesh)
    with pytest.raises(TypeError):
        adapter.fit_triqs_tau(
            basis,
            adapter.to_triqs_iw(
                fitted, meshes.MeshImFreq(beta=4, statistic="Fermion", n_iw=3)
            ),
        )
    with pytest.raises(ValueError):
        adapter.to_triqs_tau(fitted, meshes.MeshImTime())
    with pytest.raises(ValueError):
        adapter.to_triqs_iw(fitted, meshes.MeshImFreq())
    source.data[1] = np.nan
    with pytest.raises(ValueError, match="finite target shape"):
        adapter.fit_triqs_tau(basis, source)


class _FrequencyMesh:
    beta = 4.0
    statistic = "Fermion"
    positive_only = False
    first_index = -1
    last_index = 0

    def __init__(self, indices: list[Any], values: list[complex]) -> None:
        self.points = [
            SimpleNamespace(index=n, value=w)
            for n, w in zip(indices, values, strict=True)
        ]

    def __len__(self) -> int:
        return len(self.points)

    def __iter__(self) -> Any:
        return iter(self.points)


@pytest.mark.parametrize(
    "indices,frequencies",
    [
        ([-1, 0], [1j, 2j]),
        ([0, -1], [1j, -1j]),
        ([-1, -1], [-1j, -1j]),
        ([True, 0], [-1j, 1j]),
        ([-1.5, 0], [-1j, 1j]),
        ([-1, 0], [complex(np.nan, 0), 1j]),
    ],
)
def test_inconsistent_frequency_coordinate_metadata_is_rejected(
    indices: list[Any], frequencies: list[complex]
) -> None:
    mesh = _FrequencyMesh(indices, frequencies)
    with pytest.raises(ValueError):
        adapter._iw_indices(mesh, 4, SimpleNamespace(MeshImFreq=_FrequencyMesh))


class _TimeMesh:
    beta = 4.0
    statistic = "Fermion"

    def __init__(self, coordinates: list[float]) -> None:
        self.points = [SimpleNamespace(value=value) for value in coordinates]

    def __len__(self) -> int:
        return len(self.points)

    def __iter__(self) -> Any:
        return iter(self.points)


@pytest.mark.parametrize(
    "coordinates",
    [
        [],
        [0.0],
        [0.0, 4.0, np.nan],
        [0.0, 3.0, 2.0, 4.0],
        [0.0, 0.0, 4.0],
        [-1.0, 0.0, 4.0],
        [0.0, 3.9],
    ],
)
def test_invalid_time_order_or_boundary_metadata_is_rejected(
    coordinates: list[float],
) -> None:
    with pytest.raises(ValueError):
        adapter._tau_coordinates(
            _TimeMesh(coordinates), 4, SimpleNamespace(MeshImTime=_TimeMesh)
        )


@pytest.fixture
def adapter_probe(monkeypatch: pytest.MonkeyPatch) -> Any:
    directory = Path(__file__).resolve().parents[1] / "scripts/providers"
    monkeypatch.syspath_prepend(str(directory))
    specification = importlib.util.spec_from_file_location(
        "_kernelgecp_adapter_probe_test", directory / "triqs_adapter_smoke.py"
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class _WheelDistribution:
    version = "1.2.0-test"

    def __init__(self, root: Path) -> None:
        self.root = root
        self.files: list[str] | None = [
            "kernelgecp/basis.py",
            "kernelgecp/integrations/triqs.py",
        ]
        self.record: str | None = "\n".join(f"{name},," for name in self.files)
        self.direct_url: str | None = None

    def read_text(self, name: str) -> str | None:
        return self.record if name == "RECORD" else self.direct_url

    def locate_file(self, relative: str) -> Path:
        return self.root / relative


def test_installed_wheel_guard_binds_distribution_files_without_private_paths(
    adapter_probe: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    distribution = _WheelDistribution(tmp_path)
    monkeypatch.setattr(
        adapter_probe.package_metadata, "distribution", lambda _: distribution
    )
    modules = {
        name: tmp_path / relative
        for name, relative in adapter_probe.PACKAGE_FILES.items()
    }
    identity = adapter_probe.installed_wheel_identity(modules)
    assert identity["module_count"] == 2
    assert identity["editable"] is False
    assert identity["version"] == "1.2.0-test"
    rendered = json.dumps(identity)
    assert str(tmp_path) not in rendered
    assert "file://" not in rendered


@pytest.mark.parametrize(
    "corruption",
    [
        "missing_files",
        "missing_record",
        "missing_member",
        "record_omits_member",
        "duplicate_record_member",
        "wrong_import_path",
        "editable",
        "malformed_direct_url",
    ],
)
def test_installed_wheel_guard_rejects_shadowed_or_nonwheel_imports(
    adapter_probe: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    corruption: str,
) -> None:
    distribution = _WheelDistribution(tmp_path / "installed")
    modules = {
        name: distribution.root / relative
        for name, relative in adapter_probe.PACKAGE_FILES.items()
    }
    if corruption == "missing_files":
        distribution.files = None
    elif corruption == "missing_record":
        distribution.record = None
    elif corruption == "missing_member":
        distribution.files = ["kernelgecp/basis.py"]
    elif corruption == "record_omits_member":
        distribution.record = "kernelgecp/basis.py,,\n"
    elif corruption == "duplicate_record_member":
        distribution.record += "\nkernelgecp/basis.py,,\n"
    elif corruption == "wrong_import_path":
        # Another checkout may have matching module bytes and lie outside the
        # current checkout, but that does not make it the installed wheel.
        modules["basis"] = tmp_path / "other-checkout/kernelgecp/basis.py"
    elif corruption == "editable":
        distribution.direct_url = json.dumps({"dir_info": {"editable": True}})
    else:
        distribution.direct_url = "not-json"
    monkeypatch.setattr(
        adapter_probe.package_metadata, "distribution", lambda _: distribution
    )
    with pytest.raises(RuntimeError):
        adapter_probe.installed_wheel_identity(modules)


def test_installed_wheel_guard_rejects_missing_distribution_metadata(
    adapter_probe: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    def absent(name: str) -> Any:
        raise adapter_probe.package_metadata.PackageNotFoundError(name)

    monkeypatch.setattr(adapter_probe.package_metadata, "distribution", absent)
    with pytest.raises(RuntimeError, match="distribution metadata is required"):
        adapter_probe.installed_wheel_identity({})
