"""Owned molecular Green 1.0 GF2 time data with bound external provenance.

HDF5 alone does not identify units, method or convergence. This interface
requires a versioned JSON manifest binding those claims to the actual native
files/log/build receipt. It never reads or unpickles Green's serialized Cell.
Optional h5py is imported only by read_green_tau; PySCF is not needed to read.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import operator
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import numpy.typing as npt

from ..fitting import _immutable_copy
from ..types import FloatArray

type ComplexArray = npt.NDArray[np.complex128]
_MAX_ARRAY_ELEMENTS = 2**25
_SHA256 = re.compile(r"[0-9a-f]{64}")
_REVISION = re.compile(r"[0-9a-f]{40}")
_ENERGIES = ("Energy_1b", "Energy_HF", "Energy_2b")
_REQUIRED = (
    "G_tau/data",
    "G_tau/mesh",
    "Selfenergy/data",
    "Selfenergy/mesh",
    "Sigma1",
    "mu",
    *_ENERGIES,
)


class GreenFormatError(ValueError):
    """A molecular checkpoint or its provenance violates the supported format."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GreenFormatError(message)


def _number(value: Any, name: str, *, positive: bool = False) -> float:
    _require(type(value) in (int, float), f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except OverflowError as error:
        raise GreenFormatError(f"{name} must be finite in binary64") from error
    _require(
        np.isfinite(result) and (not positive or result > 0),
        f"{name} must be finite" + (" and positive" if positive else ""),
    )
    return result


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    def invalid_constant(value: str) -> None:
        raise GreenFormatError(f"nonfinite JSON constant {value}")

    try:
        value = json.loads(path.read_text(), parse_constant=invalid_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise GreenFormatError("cannot read the required JSON provenance") from error
    _require(isinstance(value, dict), "JSON provenance must be an object")
    return cast(dict[str, Any], value)


def _string(value: Any, name: str) -> str:
    _require(isinstance(value, str) and bool(value), f"missing {name}")
    return cast(str, value)


def _version(value: Any, name: str) -> str:
    text = _string(value, name)
    _require(
        re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", text) is not None,
        f"{name} must identify an exact release",
    )
    return text


def _frozen_complex(value: npt.ArrayLike) -> ComplexArray:
    return cast(ComplexArray, _immutable_copy(np.asarray(value, dtype=np.complex128)))


@dataclass(frozen=True, slots=True)
class GreenUnits:
    energy: str = "hartree"
    imaginary_time: str = "hartree^-1"
    geometry: str = "angstrom"
    green_tau: str = "dimensionless"
    dynamic_selfenergy_tau: str = "hartree^2"


@dataclass(frozen=True, slots=True)
class GreenEnergies:
    one_body_hartree: float
    hartree_fock_hartree: float
    correlation_hartree: float
    nuclear_hartree: float
    initial_hf_hartree: float

    @property
    def total_hartree(self) -> float:
        """Energy_HF already includes nuclear energy; do not add Energy_1b."""
        return self.hartree_fock_hartree + self.correlation_hartree


@dataclass(frozen=True, slots=True)
class GreenConvergence:
    source_iteration: int
    energy_change_sum_hartree: float
    green_change_max: float
    energy_threshold_hartree: float
    green_threshold: float
    native_log_reports_convergence: bool = True
    scope: str = "native_final_iteration_pair_and_declared_thresholds"


@dataclass(frozen=True, slots=True)
class GreenProvider:
    green: str
    grids: str
    green_mbtools: str
    pyscf: str
    writer_numpy: str
    writer_h5py: str
    native_executable_sha256: str
    source_revisions: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class GreenTauData:
    beta: float
    tau: FloatArray
    g_tau: ComplexArray
    sigma_tau: ComplexArray
    overlap: ComplexArray
    h_core: ComplexArray
    hf_fock: ComplexArray
    sigma_static: ComplexArray
    mu: float
    iteration: int
    source_iteration: int
    converged: bool
    convergence: GreenConvergence
    energies: GreenEnergies
    molecule: str
    orbital_basis: str
    storage_basis: str
    atoms_angstrom: tuple[tuple[str, float, float, float], ...]
    electrons: int
    provider: GreenProvider
    file_sha256: tuple[tuple[str, str], ...]
    manifest_sha256: str
    overlap_min_eigenvalue: float
    units: GreenUnits = GreenUnits()
    axes: tuple[str, ...] = ("tau", "spin", "k", "orbital_row", "orbital_column")
    certification_status: str = "not_certified_native_data_import"

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "tau", _immutable_copy(np.asarray(self.tau, dtype=float))
        )
        for name in (
            "g_tau",
            "sigma_tau",
            "overlap",
            "h_core",
            "hf_fock",
            "sigma_static",
        ):
            object.__setattr__(self, name, _frozen_complex(getattr(self, name)))

    @property
    def green_tau(self) -> ComplexArray:
        return self.g_tau

    @property
    def dynamic_selfenergy_tau(self) -> ComplexArray:
        return self.sigma_tau


def _manifest(
    path: Path,
    input_path: Path,
    results_path: Path,
) -> tuple[dict[str, Any], GreenProvider, tuple[tuple[str, str], ...], dict[str, Path]]:
    manifest = _json(path)
    _require(
        manifest.get("schema") == "kernelgecp.green-molecular-run.v1",
        "a bound kernelgecp.green-molecular-run.v1 manifest is required",
    )
    _require(
        manifest.get("method") == "GF2" and manifest.get("restricted") is True,
        "the supported molecular interface requires restricted GF2",
    )
    _require(
        manifest.get("storage_basis") == "symmetric_lowdin",
        "the supported storage basis is symmetric_lowdin",
    )
    units = manifest.get("units")
    expected_units = {
        "energy": "hartree",
        "imaginary_time": "hartree^-1",
        "geometry": "angstrom",
        "green_tau": "dimensionless",
        "dynamic_selfenergy_tau": "hartree^2",
    }
    _require(units == expected_units, "explicit atomic units are required")
    files = manifest.get("files")
    _require(
        isinstance(files, dict)
        and {"input", "results", "native_log", "build_provenance", "provider_receipt"}
        <= set(files),
        "required native file bindings are missing",
    )
    files = cast(dict[str, Any], files)
    paths: dict[str, Path] = {}
    hashes: list[tuple[str, str]] = []
    for name, entry in files.items():
        _require(isinstance(entry, dict), "invalid native file binding")
        relative = Path(_string(entry.get("path"), "artifact path"))
        _require(
            not relative.is_absolute() and ".." not in relative.parts,
            "manifest artifact paths must stay below the manifest directory",
        )
        expected_hash = entry.get("sha256")
        _require(
            isinstance(expected_hash, str)
            and _SHA256.fullmatch(expected_hash) is not None,
            "native artifacts require full SHA256 identities",
        )
        target = {"input": input_path, "results": results_path}.get(
            name, path.parent / relative
        )
        try:
            actual_hash = _hash(target)
        except OSError as error:
            raise GreenFormatError(
                f"bound native artifact is unavailable: {name}"
            ) from error
        _require(actual_hash == expected_hash, f"native artifact hash mismatch: {name}")
        paths[name] = target
        hashes.append((name, actual_hash))
    provider = manifest.get("provider")
    _require(isinstance(provider, dict), "provider versions are required")
    provider = cast(dict[str, Any], provider)
    _require(
        provider.get("green") == "1.0.0" and provider.get("green-mbtools") == "1.0.0",
        "this reader supports the Green/green-mbtools 1.0.0 layout",
    )
    revisions = provider.get("source_revisions")
    _require(
        isinstance(revisions, dict) and bool(revisions), "source revisions are required"
    )
    revisions = cast(dict[str, str], revisions)
    for name, revision in revisions.items():
        _string(name, "source name")
        _require(
            isinstance(revision, str) and _REVISION.fullmatch(revision) is not None,
            "source revisions must be full Git object identities",
        )
    binary_hash = provider.get("native_executable_sha256")
    _require(
        isinstance(binary_hash, str) and _SHA256.fullmatch(binary_hash) is not None,
        "native executable identity is required",
    )
    build = _json(paths["build_provenance"])
    _require(
        build.get("native_executable_sha256") == binary_hash,
        "native executable identity disagrees with build provenance",
    )
    build_sources = build.get("sources")
    _require(isinstance(build_sources, dict), "missing build source identities")
    build_sources = cast(dict[str, Any], build_sources)
    source_rows = build_sources.get("sources")
    _require(
        isinstance(source_rows, list) and bool(source_rows),
        "missing native source revisions",
    )
    source_rows = cast(list[dict[str, str]], source_rows)
    _require(
        all(
            isinstance(row, dict)
            and isinstance(row.get("name"), str)
            and isinstance(row.get("revision"), str)
            for row in source_rows
        ),
        "malformed build source identity list",
    )
    _require(
        len(source_rows) == len(revisions)
        and {row["name"]: row["revision"] for row in source_rows} == revisions,
        "source revisions disagree with build provenance",
    )
    receipt = _json(paths["provider_receipt"])
    _require(
        receipt.get("schema") == "kernelgecp.provider-smoke.v1"
        and receipt.get("provider") == "green"
        and receipt.get("status") == "passed",
        "a passed actual native provider receipt is required",
    )
    checks = receipt.get("checks")
    _require(
        isinstance(checks, dict)
        and checks.get("numerical") is True
        and checks.get("native_gf2_self_consistency") is True,
        "native provider receipt does not contain a real GF2 convergence check",
    )
    versions = receipt.get("versions")
    _require(isinstance(versions, dict), "actual writer versions are missing")
    versions = cast(dict[str, Any], versions)
    for name in ("green", "green-mbtools", "pyscf", "numpy", "h5py"):
        _require(
            versions.get(name) == provider.get(name),
            "writer version/receipt mismatch: " + name,
        )
    identity = GreenProvider(
        _version(provider.get("green"), "Green"),
        _version(provider.get("green-grids"), "Green grids"),
        _version(provider.get("green-mbtools"), "green-mbtools"),
        _version(provider.get("pyscf"), "PySCF"),
        _version(provider.get("numpy"), "writer NumPy"),
        _version(provider.get("h5py"), "writer h5py"),
        cast(str, binary_hash),
        tuple(sorted(revisions.items())),
    )
    return manifest, identity, tuple(sorted(hashes)), paths


def _dataset(group: Any, path: str, h5py: Any) -> Any:
    current = group
    for component in path.split("/"):
        _require(
            isinstance(current, h5py.Group) and component in current,
            "missing required native dataset: " + path,
        )
        _require(
            isinstance(current.get(component, getlink=True), h5py.HardLink),
            "native required datasets cannot use external or symbolic links",
        )
        current = current[component]
    _require(
        isinstance(current, h5py.Dataset), "required path is not a dataset: " + path
    )
    return current


def _scalar(group: Any, path: str, h5py: Any, *, integer: bool = False) -> float | int:
    dataset = _dataset(group, path, h5py)
    _require(
        dataset.shape == ()
        and (
            dataset.dtype.kind in "iu"
            if integer
            else dataset.dtype == np.dtype(np.float64)
        ),
        "invalid scalar dataset: " + path,
    )
    value = dataset[()]
    _require(bool(np.isfinite(value)), "nonfinite scalar dataset: " + path)
    return int(value) if integer else float(value)


def _array(
    group: Any,
    path: str,
    shape: tuple[int, ...],
    h5py: Any,
    *,
    complex_output: bool = True,
) -> np.ndarray:
    dataset = _dataset(group, path, h5py)
    expected_dtype = np.dtype(np.complex128 if complex_output else np.float64)
    _require(
        dataset.shape == shape and dataset.dtype == expected_dtype,
        "unexpected shape/dtype for native dataset: " + path,
    )
    _require(
        math_product(shape) <= _MAX_ARRAY_ELEMENTS,
        "native data exceeds the reader allocation budget",
    )
    value = np.asarray(dataset[()])
    _require(bool(np.all(np.isfinite(value))), "nonfinite native data: " + path)
    return value


def math_product(shape: tuple[int, ...]) -> int:
    result = 1
    for size in shape:
        result *= size
    return result


def _input_matrix(
    source: Any, name: str, shape: tuple[int, ...], h5py: Any
) -> ComplexArray:
    dataset = _dataset(source, "HF/" + name, h5py)
    marker = dataset.attrs.get("__complex__")
    _require(
        isinstance(marker, (int, np.integer))
        and not isinstance(marker, (bool, np.bool_))
        and marker == 1,
        "input HF matrices require __complex__=1",
    )
    rows = _array(source, "HF/" + name, (*shape, 2), h5py, complex_output=False)
    return np.asarray(rows[..., 0] + 1j * rows[..., 1], dtype=np.complex128)


def _checkpoint(
    source: Any, iteration: int, shape: tuple[int, ...], beta: float, h5py: Any
) -> dict[str, Any]:
    name = f"iter{iteration}"
    _require(
        name in source and isinstance(source.get(name, getlink=True), h5py.HardLink),
        f"selected native checkpoint {name} is unavailable",
    )
    group = source[name]
    for path in _REQUIRED:
        _dataset(group, path, h5py)
    mesh = _dataset(group, "G_tau/mesh", h5py)
    _require(
        len(mesh.shape) == 1 and mesh.shape[0] >= 2,
        "native time mesh must contain ordered one-sided endpoints",
    )
    tau = _array(group, "G_tau/mesh", mesh.shape, h5py, complex_output=False)
    _require(
        tau[0] == 0
        and abs(tau[-1] - beta) <= 1e-13 * max(1, beta)
        and bool(np.all(np.diff(tau) > 0)),
        "native imaginary-time mesh/beta mismatch",
    )
    sigma_mesh = _array(group, "Selfenergy/mesh", tau.shape, h5py, complex_output=False)
    _require(np.array_equal(tau, sigma_mesh), "Green and self-energy meshes differ")
    return {
        "tau": tau,
        "g_tau": _array(group, "G_tau/data", (len(tau), *shape), h5py),
        "sigma_tau": _array(group, "Selfenergy/data", (len(tau), *shape), h5py),
        "sigma_static": _array(group, "Sigma1", shape, h5py),
        "mu": _scalar(group, "mu", h5py),
        **{field: _scalar(group, field, h5py) for field in _ENERGIES},
    }


def read_green_tau(
    input_path: str | Path,
    results_path: str | Path,
    run_manifest: str | Path,
    *,
    iteration: int | None = None,
) -> GreenTauData:
    """Read a bound complete molecular checkpoint without altering native axes.

    None selects the actual /iter checkpoint, which must be complete and agree
    with the converged run manifest. Explicit complete earlier checkpoints are
    returned as historical data with converged=False; no older fallback occurs.
    """
    try:
        h5py = importlib.import_module("h5py")
    except ImportError as error:
        raise ImportError("read_green_tau requires kernelgecp[green] (h5py)") from error
    input_file, results_file, manifest_file = map(
        Path, (input_path, results_path, run_manifest)
    )
    manifest_hash = _hash(manifest_file)
    manifest, provider, hashes, paths = _manifest(
        manifest_file, input_file, results_file
    )
    beta = _number(manifest.get("beta_inverse_hartree"), "physical beta", positive=True)
    geometry = manifest.get("geometry")
    _require(
        isinstance(geometry, list) and bool(geometry),
        "explicit geometry is required; Cell is never deserialized",
    )
    geometry = cast(list[dict[str, Any]], geometry)
    atoms: list[tuple[str, float, float, float]] = []
    for atom in geometry:
        _require(isinstance(atom, dict), "geometry entries must be objects")
        symbol = _string(atom.get("element"), "element")
        xyz = atom.get("xyz")
        _require(
            isinstance(xyz, list) and len(xyz) == 3,
            "geometry requires three Cartesian coordinates",
        )
        xyz = cast(list[Any], xyz)
        atoms.append(
            (
                symbol,
                _number(xyz[0], "coordinate"),
                _number(xyz[1], "coordinate"),
                _number(xyz[2], "coordinate"),
            )
        )
    convergence = manifest.get("convergence")
    _require(isinstance(convergence, dict), "actual convergence provenance is required")
    convergence = cast(dict[str, Any], convergence)
    final = convergence.get("iteration")
    _require(
        type(final) is int and final >= 2, "converged final iteration pair is required"
    )
    final = cast(int, final)
    energy_threshold = _number(
        convergence.get("energy_change_threshold_hartree"),
        "energy threshold",
        positive=True,
    )
    green_threshold = _number(
        convergence.get("green_change_threshold"), "Green threshold", positive=True
    )
    log = paths["native_log"].read_text()
    _require(
        "Simulation Converged" in log
        and "Reached Maximum number of iterations" not in log,
        "native log does not establish convergence",
    )
    if iteration is not None:
        _require(
            not isinstance(iteration, (bool, np.bool_)), "iteration must be an integer"
        )
        try:
            selected = operator.index(iteration)
        except TypeError as error:
            raise GreenFormatError("iteration must be an integer") from error
        _require(
            1 <= selected <= final, "requested iteration is outside the recorded run"
        )
    else:
        selected = final
    try:
        with (
            h5py.File(input_file, "r") as source,
            h5py.File(results_file, "r") as results,
        ):
            _require(
                source.attrs.get("__green_version__") == provider.green,
                "input provider version disagrees with manifest",
            )
            _require(
                results.attrs.get("__grids_version__") == provider.grids,
                "output grids version disagrees with manifest",
            )
            ns = int(_scalar(source, "params/ns", h5py, integer=True))
            nk = int(_scalar(source, "params/nk", h5py, integer=True))
            nso = int(_scalar(source, "params/nso", h5py, integer=True))
            electrons = int(_scalar(source, "params/nel_cell", h5py, integer=True))
            _require(
                ns == nk == 1 and nso >= 1 and 0 <= electrons <= 2 * nso,
                "unsupported restricted molecular spin/k/orbital dimensions",
            )
            shape = (ns, nk, nso, nso)
            overlap = _input_matrix(source, "S-k", shape, h5py)
            h_core = _input_matrix(source, "H-k", shape, h5py)
            hf_fock = _input_matrix(source, "Fock-k", shape, h5py)
            for name, matrix in (
                ("overlap", overlap),
                ("one-body", h_core),
                ("HF Fock", hf_fock),
            ):
                _require(
                    float(np.max(abs(matrix - matrix.swapaxes(-1, -2).conj())))
                    <= 1e-10,
                    f"non-Hermitian native {name} matrix",
                )
            minimum = float(np.min(np.linalg.eigvalsh(overlap)))
            _require(
                minimum > 0 and float(np.max(abs(overlap - np.eye(nso)))) <= 1e-10,
                "declared Lowdin storage requires positive identity overlap",
            )
            source_iteration = int(_scalar(results, "iter", h5py, integer=True))
            _require(
                source_iteration == final,
                "/iter differs from the converged manifest; no older fallback",
            )
            current = _checkpoint(results, final, shape, beta, h5py)
            previous = _checkpoint(results, final - 1, shape, beta, h5py)
            _require(
                np.array_equal(current["tau"], previous["tau"]),
                "final iteration pair meshes differ",
            )
            energy_change = sum(
                abs(current[field] - previous[field]) for field in _ENERGIES
            )
            green_change = float(np.max(abs(current["g_tau"] - previous["g_tau"])))
            _require(
                np.isfinite(energy_change)
                and energy_change < energy_threshold
                and green_change < green_threshold,
                "stored final iteration pair does not converge",
            )
            data = (
                current
                if selected == final
                else _checkpoint(results, selected, shape, beta, h5py)
            )
            nuclear = float(_scalar(source, "HF/Energy_nuc", h5py))
            initial_hf = float(_scalar(source, "HF/Energy", h5py))
            _require(
                np.isfinite(data["Energy_HF"] + data["Energy_2b"]),
                "native total energy is nonfinite",
            )
    except (OSError, KeyError) as error:
        raise GreenFormatError(
            "cannot read a complete supported native checkpoint"
        ) from error
    _require(
        all(_hash(paths[name]) == expected for name, expected in hashes)
        and _hash(manifest_file) == manifest_hash,
        "native files changed while reading",
    )
    return GreenTauData(
        beta,
        data["tau"],
        data["g_tau"],
        data["sigma_tau"],
        overlap,
        h_core,
        hf_fock,
        data["sigma_static"],
        data["mu"],
        selected,
        source_iteration,
        selected == final,
        GreenConvergence(
            final, energy_change, green_change, energy_threshold, green_threshold
        ),
        GreenEnergies(
            data["Energy_1b"], data["Energy_HF"], data["Energy_2b"], nuclear, initial_hf
        ),
        _string(manifest.get("molecule"), "molecule"),
        _string(manifest.get("orbital_basis"), "orbital basis"),
        "symmetric_lowdin",
        tuple(atoms),
        electrons,
        provider,
        hashes,
        manifest_hash,
        minimum,
    )
