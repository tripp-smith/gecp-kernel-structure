"""Independent saved-array and raw-HDF5 audit of the four-case EF campaign.

The campaign runner and NativeIRSolution are never imported. Native fields
are reconstructed by overdetermined pivoted QR and direct exponential Fourier
quadrature. Finite-model discrepancy remains separate from representation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np
from scipy.linalg import expm, lstsq

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("kernelgecp", "pydlr", "sparse_ir")
SOURCE_FILES = {
    "experiments/molecular_green.py",
    "experiments/molecular_green.toml",
    "experiments/molecular_native.py",
    "experiments/molecular_reference.py",
    "experiments/dlr_adapters.py",
    "kernelgecp/integrations/green.py",
    "kernelgecp/basis.py",
    "kernelgecp/fitting.py",
    "environments/green-science/uv.lock",
}
NATIVE_MECHANISMS = (
    "environments/green/Dockerfile",
    "environments/green/build_green.py",
    "environments/green/sources.json",
    "environments/green/requirements.lock",
    "scripts/providers/green_molecular_runs.py",
    "scripts/providers/bind_green_manifest.py",
)
FIELDS = ("green_tau", "green_iw", "sigma_tau", "sigma_iw")
WORKER_KEYS = {
    "training_rows",
    "tau",
    "indices",
    "g_coefficients",
    "sigma_coefficients",
    "green_endpoints",
    *FIELDS,
}


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> dict[str, Any]:
    require(
        not path.is_symlink() and path.stat().st_size < 16 * 1024**2,
        "unsafe JSON receipt",
    )

    def unique(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result

    def invalid(value: str) -> Any:
        raise ValueError("nonfinite JSON number: " + value)

    result = json.loads(
        path.read_text(), object_pairs_hook=unique, parse_constant=invalid
    )
    json.dumps(result, allow_nan=False)
    require(isinstance(result, dict), "receipt must be an object")
    return result


def relative(base: Path, name: str) -> Path:
    part = PurePosixPath(name)
    require(
        not part.is_absolute() and ".." not in part.parts and "\\" not in name,
        "unsafe relative artifact path",
    )
    path = base.joinpath(*part.parts)
    require(not path.is_symlink(), "artifact cannot be a symlink")
    return path


def numeric_archive(path: Path, keys: set[str]) -> dict[str, np.ndarray]:
    require(
        not path.is_symlink() and path.stat().st_size <= 64 * 1024**2,
        "oversized numeric archive",
    )
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        require(
            len(entries) == len(keys)
            and {item.filename for item in entries} == {key + ".npy" for key in keys},
            "numeric archive member identity mismatch",
        )
        require(
            sum(item.file_size for item in entries) <= 64 * 1024**2,
            "numeric archive expansion exceeds budget",
        )
        for entry in entries:
            with archive.open(entry) as stream:
                version = np.lib.format.read_magic(stream)
                require(version in ((1, 0), (2, 0)), "unsupported numeric array format")
                header = (
                    np.lib.format.read_array_header_1_0
                    if version == (1, 0)
                    else np.lib.format.read_array_header_2_0
                )(stream)
                shape, _, dtype = header
                require(
                    0 < len(shape) <= 6
                    and all(type(x) is int and x > 0 for x in shape)
                    and math.prod(shape) <= 2**23
                    and dtype.kind in "fciu"
                    and not dtype.hasobject,
                    "unsafe numeric array header",
                )
                require(
                    stream.tell() + math.prod(shape) * dtype.itemsize
                    == entry.file_size,
                    "numeric payload size mismatch",
                )
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in keys}
    require(
        all(np.all(np.isfinite(value)) for value in arrays.values()),
        "nonfinite numeric artifact",
    )
    return arrays


def maximum(value: np.ndarray) -> float:
    require(
        value.size > 0 and np.all(np.isfinite(value)),
        "empty or nonfinite numerical residual",
    )
    return float(np.max(np.abs(value)))


def complex_number(value: Any) -> complex:
    if isinstance(value, dict):
        require(set(value) == {"real", "imag"}, "invalid complex scalar encoding")
        result = complex(value["real"], value["imag"])
    else:
        result = complex(value)
    require(np.isfinite(result), "nonfinite complex scalar")
    return result


def causality(values: np.ndarray, indices: np.ndarray) -> float:
    require(
        values.ndim == 5
        and values.shape[1:] == (1, 1, 2, 2)
        and values.shape[0] == len(indices),
        "unexpected molecular matrix axes",
    )
    positive = values[indices >= 0]
    require(len(positive) > 0, "causality requires positive physical frequencies")
    anti = (positive - positive.conj().swapaxes(-1, -2)) / (2j)
    return float(np.max(np.linalg.eigvalsh(anti)))


class RawNative:
    """Read native matrices directly; reconstruct with QR rather than runner LU."""

    def __init__(
        self,
        folder: Path,
        bond: float,
        native_lambda: float,
        expected_revision: str | None = None,
    ) -> None:
        import h5py
        import sparse_ir

        self.folder = folder
        self.manifest = read(folder / "bound-manifest.json")
        m = self.manifest
        require(
            m["schema"] == "kernelgecp.green-molecular-run.v1"
            and re.fullmatch("[0-9a-f]{40}", m["source_revision"]) is not None
            and (
                expected_revision is None or m["source_revision"] == expected_revision
            ),
            "unexpected native generation identity",
        )
        require(
            m["beta_inverse_hartree"] == 20
            and m["storage_basis"] == "symmetric_lowdin"
            and m["method"] == "GF2"
            and m["restricted"] is True,
            "wrong native physical model",
        )
        require(
            m["molecule"] == "H2"
            and m["orbital_basis"].lower() == "sto-3g"
            and m["geometry"]
            == [
                {"element": "H", "xyz": [0.0, 0.0, 0.0]},
                {"element": "H", "xyz": [0.0, 0.0, bond]},
            ],
            "native geometry/basis mismatch",
        )
        require(
            m["native_grid_lambda"] == native_lambda, "native cutoff label mismatch"
        )
        require(
            m["convergence"]["energy_change_threshold_hartree"] <= 1e-10
            and m["convergence"]["green_change_threshold"] <= 1e-6,
            "native convergence budget weakened",
        )
        for record in m["files"].values():
            require(
                sha(relative(folder, record["path"])) == record["sha256"],
                "native file hash mismatch",
            )
        log = relative(folder, m["files"]["native_log"]["path"]).read_text()
        require(
            "Simulation Converged" in log
            and "Reached Maximum number of iterations" not in log,
            "native log does not establish convergence",
        )
        self.beta = 20.0
        self.bond = bond
        self.native_lambda = native_lambda
        with (
            h5py.File(folder / "input.h5") as inputs,
            h5py.File(folder / "sim.h5") as results,
        ):
            require(
                inputs.attrs["__green_version__"] == "1.0.0"
                and results.attrs["__grids_version__"] == "0.3.0",
                "native writer version mismatch",
            )
            require(
                [
                    int(inputs["params/" + name][()])
                    for name in ("ns", "nk", "nso", "nel_cell")
                ]
                == [1, 1, 2, 2],
                "wrong native restricted dimensions",
            )

            def matrix(name: str) -> np.ndarray:
                raw = np.asarray(inputs["HF/" + name])
                require(
                    raw.shape == (1, 1, 2, 2, 2)
                    and inputs["HF/" + name].attrs["__complex__"] == 1,
                    "wrong native input matrix layout",
                )
                return raw[..., 0] + 1j * raw[..., 1]

            self.overlap, self.core, self.fock = (
                matrix(name) for name in ("S-k", "H-k", "Fock-k")
            )
            require(
                maximum(self.overlap - np.eye(2)) <= 1e-10,
                "native storage is not declared orthonormal Lowdin",
            )
            self.nuclear = float(inputs["HF/Energy_nuc"][()])
            self.initial_hf = float(inputs["HF/Energy"][()])
            self.iteration = int(results["iter"][()])
            require(
                self.iteration == m["convergence"]["iteration"] and self.iteration >= 2,
                "native final iteration mismatch",
            )
            group = results[f"iter{self.iteration}"]
            previous = results[f"iter{self.iteration - 1}"]
            self.tau = np.asarray(group["G_tau/mesh"])
            self.fields = {
                "green": np.asarray(group["G_tau/data"]),
                "sigma": np.asarray(group["Selfenergy/data"]),
            }
            self.static = np.asarray(group["Sigma1"])
            self.mu = float(group["mu"][()])
            self.energy = float(group["Energy_HF"][()] + group["Energy_2b"][()])
            energy_change = sum(
                abs(float(group[name][()] - previous[name][()]))
                for name in ("Energy_1b", "Energy_HF", "Energy_2b")
            )
            green_change = maximum(
                self.fields["green"] - np.asarray(previous["G_tau/data"])
            )
            require(
                energy_change < m["convergence"]["energy_change_threshold_hartree"]
                and green_change < m["convergence"]["green_change_threshold"],
                "native final pair is not converged",
            )
            require(
                np.array_equal(self.tau, np.asarray(group["Selfenergy/mesh"]))
                and np.array_equal(self.tau, np.asarray(previous["G_tau/mesh"])),
                "native meshes differ",
            )
        require(
            self.tau[0] == 0 and self.tau[-1] == 20 and np.all(np.diff(self.tau) > 0),
            "native physical tau/beta mismatch",
        )
        require(
            all(
                value.shape == (len(self.tau), 1, 1, 2, 2)
                and np.all(np.isfinite(value))
                for value in self.fields.values()
            ),
            "native field axes/values invalid",
        )
        with h5py.File(folder / "native-grid.h5") as grid:
            require(
                grid.attrs["__grids_version__"] == "0.3.0"
                and float(grid["fermi/metadata/lambda"][()]) == native_lambda,
                "wrong native IR grid identity",
            )
            x = np.asarray(grid["fermi/xgrid"])
            interior = np.asarray(grid["fermi/uxl"])
            table = np.vstack(
                [
                    np.asarray(grid["fermi/u1l_neg"]),
                    interior,
                    np.asarray(grid["fermi/u1l_pos"]),
                ]
            ) * math.sqrt(2 / self.beta)
            stored_frequency = np.asarray(grid["fermi/uwl"]) * math.sqrt(self.beta)
            stored_indices = np.asarray(grid["fermi/ngrid"], np.int64)
        self.dimension = table.shape[1]
        require(
            table.shape == (len(self.tau), self.dimension)
            and maximum(self.tau - np.r_[0.0, (x + 1) * self.beta / 2, self.beta])
            <= 1e-12,
            "native IR physical grid shape mismatch",
        )
        self.basis = sparse_ir.FiniteTempBasis("F", beta=1.0, wmax=native_lambda)
        defining = np.r_[0.0, (x + 1) / 2, 1.0]
        normalized = np.asarray(self.basis.u[: self.dimension](defining)).T / math.sqrt(
            self.beta
        )
        self.table_error = maximum(normalized - table)
        transformed = np.asarray(
            self.basis.uhat[: self.dimension](2 * stored_indices + 1)
        ).T * math.sqrt(self.beta)
        require(
            max(self.table_error, maximum(transformed - stored_frequency)) <= 2e-10,
            "native IR normalization differs independently",
        )
        self.coefficients = {}
        self.reconstruction_errors = {}
        for field, values in self.fields.items():
            coefficients, _, rank, _ = lstsq(
                table, values.reshape(len(table), 4), cond=1e-13, lapack_driver="gelsy"
            )
            require(rank == self.dimension, "independent native QR is rank deficient")
            self.coefficients[field] = coefficients
            self.reconstruction_errors[field] = maximum(
                (table @ coefficients).reshape(values.shape) - values
            )
        require(
            max(self.reconstruction_errors.values()) <= 2e-9,
            "independent native QR fails its stored fields",
        )
        x, w = np.polynomial.legendre.leggauss(512)
        self.quadrature_tau = self.beta * (x + 1) / 2
        self.quadrature_weights = self.beta * w / 2
        self.quadrature_fields = {
            name: self.time(self.quadrature_tau, name) for name in self.fields
        }
        self.density = complex(-2 * np.trace(self.fields["green"][-1, 0, 0]))

    def time(self, tau: np.ndarray, field: str) -> np.ndarray:
        design = np.asarray(
            self.basis.u[: self.dimension](np.asarray(tau) / self.beta)
        ).T / math.sqrt(self.beta)
        return (design @ self.coefficients[field]).reshape(len(tau), 1, 1, 2, 2)

    def frequency(self, indices: np.ndarray, field: str) -> np.ndarray:
        frequency = (2 * np.asarray(indices) + 1) * np.pi / self.beta
        phase = (
            np.exp(1j * frequency[:, None] * self.quadrature_tau)
            * self.quadrature_weights
        )
        return (phase @ self.quadrature_fields[field].reshape(512, 4)).reshape(
            len(indices), 1, 1, 2, 2
        )

    def basis_frequency_error(self, indices: np.ndarray, field: str) -> float:
        design = np.asarray(
            self.basis.uhat[: self.dimension](2 * np.asarray(indices) + 1)
        ).T * math.sqrt(self.beta)
        return maximum(
            self.frequency(indices, field)
            - (design @ self.coefficients[field]).reshape(len(indices), 1, 1, 2, 2)
        )

    def dyson(self, indices: np.ndarray) -> float:
        z = 1j * (2 * np.asarray(indices) + 1) * np.pi / self.beta
        green, sigma = (self.frequency(indices, name) for name in ("green", "sigma"))
        operator = (
            (z + self.mu)[:, None, None, None, None] * self.overlap
            - self.core
            - self.static
            - sigma
        )
        return maximum(operator @ green - np.eye(2))


def query_rows(tau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions = np.arange(len(tau))
    training = np.unique(np.r_[positions[::2], len(tau) - 1])
    calibration = np.setdiff1d(positions[1::4], training)
    validation = np.setdiff1d(positions[3::4], training)
    return training, calibration, validation


def expected_time(native: RawNative, role: str) -> np.ndarray:
    _, calibration, validation = query_rows(native.tau)
    if role == "calibration":
        return native.tau[calibration]
    gauss = 10 * (np.polynomial.legendre.leggauss(64)[0] + 1)
    if role == "regression":
        return np.unique(np.r_[native.tau[validation], gauss])
    return 10 * (np.polynomial.legendre.leggauss(96)[0] + 1)


def build_basis(method: str, config: dict[str, Any], cutoff: float, eps: float) -> Any:
    if method == "kernelgecp":
        from kernelgecp import FermionicBasis

        return FermionicBasis(
            20,
            cutoff,
            eps,
            grid_order=config["grid_order"],
            max_rank=config["max_rank"],
            precision=config["construction_precision"],
            precision_bits=config["reference_precision_bits"],
            reference_backend=config["reference_backend"],
        )
    from experiments.dlr_adapters import build_adapter

    return build_adapter(method, 20, cutoff, eps)


def validate_arrays(
    arrays: dict[str, np.ndarray],
    native: RawNative,
    config: dict[str, Any],
    role: str,
    green_rank: int,
    sigma_rank: int,
) -> None:
    training, calibration, _ = query_rows(native.tau)
    require(
        np.array_equal(arrays["training_rows"], training),
        "worker training rows changed",
    )
    require(
        np.array_equal(arrays["indices"], config[role + "_indices"]),
        "worker frequency query identity changed",
    )
    expected = expected_time(native, role)
    require(
        arrays["tau"].shape == expected.shape
        and maximum(arrays["tau"] - expected) <= 32 * np.finfo(float).eps * 20,
        "worker physical time query rule changed",
    )
    require(
        not np.intersect1d(arrays["tau"], native.tau[training]).size,
        "validation observes fitted time rows",
    )
    if role == "validation":
        require(
            not np.intersect1d(arrays["tau"], expected_time(native, "regression")).size,
            "final times overlap retired queries",
        )
        require(
            not np.intersect1d(arrays["tau"], native.tau[calibration]).size,
            "final time observations overlap calibration",
        )
    shapes = {
        "g_coefficients": (green_rank, 1, 1, 2, 2),
        "sigma_coefficients": (sigma_rank, 1, 1, 2, 2),
        "green_endpoints": (2, 1, 1, 2, 2),
    }
    for name in FIELDS:
        shapes[name] = (
            len(arrays["tau"] if name.endswith("tau") else arrays["indices"]),
            1,
            1,
            2,
            2,
        )
    require(
        all(arrays[name].shape == shape for name, shape in shapes.items()),
        "worker matrix/coefficients shape mismatch",
    )


def worker_metrics(arrays: dict[str, np.ndarray], native: RawNative) -> dict[str, Any]:
    errors = {
        name: maximum(
            arrays[name]
            - (
                native.time(arrays["tau"], name.split("_")[0])
                if name.endswith("tau")
                else native.frequency(arrays["indices"], name.split("_")[0])
            )
        )
        for name in FIELDS
    }
    density = complex(-2 * np.trace(arrays["green_endpoints"][-1, 0, 0]))
    return {
        "errors": errors,
        "density": density,
        "native_density": native.density,
        "density_error": float(abs(density - native.density)),
        "causality_imag_max": {
            name: causality(arrays[name + "_iw"], arrays["indices"])
            for name in ("green", "sigma")
        },
    }


class FiniteTrace:
    """Four-mode Fock space built by bit parity, evaluated by thermal expm."""

    def __init__(
        self, core: np.ndarray, eri: np.ndarray, mu: float, nuclear: float
    ) -> None:
        self.operators = []
        for mode in range(4):
            annihilator = np.zeros((16, 16))
            for state in range(16):
                if state & (1 << mode):
                    sign = -1 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1
                    annihilator[state ^ (1 << mode), state] = sign
            self.operators.append(annihilator)
        hamiltonian = np.zeros((16, 16))
        for p, q in np.ndindex(2, 2):
            for spin in range(2):
                hamiltonian += (
                    core[p, q]
                    * self.operators[2 * p + spin].T
                    @ self.operators[2 * q + spin]
                )
        for p, q, r, s in np.ndindex(2, 2, 2, 2):
            for alpha, beta in np.ndindex(2, 2):
                hamiltonian += (
                    0.5
                    * eri[p, q, r, s]
                    * (
                        self.operators[2 * p + alpha].T
                        @ self.operators[2 * r + beta].T
                        @ self.operators[2 * s + beta]
                        @ self.operators[2 * q + alpha]
                    )
                )
        require(
            maximum(hamiltonian - hamiltonian.T) <= 1e-10,
            "finite molecular Hamiltonian is not Hermitian",
        )
        self.hamiltonian = (hamiltonian + hamiltonian.T) / 2
        number = np.diag([state.bit_count() for state in range(16)])
        grand = self.hamiltonian - mu * number
        self.shifted = grand - np.linalg.eigvalsh(grand)[0] * np.eye(16)
        thermal = expm(-20 * self.shifted)
        self.partition = float(np.trace(thermal))
        self.energy = float(
            np.trace(thermal @ self.hamiltonian) / self.partition + nuclear
        )
        self.density = float(np.trace(thermal @ number) / self.partition)
        self.cache = {}

    def time(self, tau: np.ndarray) -> np.ndarray:
        output = []
        for value in tau:
            key = float(value)
            if key not in self.cache:
                left, right = (
                    expm(-(20 - key) * self.shifted),
                    expm(-key * self.shifted),
                )
                matrix = np.zeros((2, 2))
                for p, q in np.ndindex(2, 2):
                    matrix[p, q] = (
                        -0.5
                        * sum(
                            np.trace(
                                left
                                @ self.operators[2 * p + spin]
                                @ right
                                @ self.operators[2 * q + spin].T
                            )
                            for spin in range(2)
                        )
                        / self.partition
                    )
                self.cache[key] = matrix
            output.append(self.cache[key])
        return np.asarray(output, dtype=complex)[:, None, None]

    def frequency(self, indices: np.ndarray) -> np.ndarray:
        x, w = np.polynomial.legendre.leggauss(512)
        times = 10 * (x + 1)
        phase = (
            np.exp(1j * (2 * np.asarray(indices)[:, None] + 1) * np.pi / 20 * times)
            * w
            * 10
        )
        return (phase @ self.time(times).reshape(512, 4)).reshape(
            len(indices), 1, 1, 2, 2
        )


def model_accounting(
    native: RawNative, tau: np.ndarray, indices: np.ndarray
) -> dict[str, Any]:
    import h5py
    from pyscf import fci, gto, scf

    molecule = gto.M(
        atom=[("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, native.bond))],
        basis="sto-3g",
        unit="Angstrom",
        spin=0,
        charge=0,
        verbose=0,
    )
    molecule.max_memory = 512
    overlap = molecule.intor("int1e_ovlp")
    eigenvalues, vectors = np.linalg.eigh(overlap)
    transform = (vectors / np.sqrt(eigenvalues)) @ vectors.T
    core = (
        transform.T
        @ (molecule.intor("int1e_kin") + molecule.intor("int1e_nuc"))
        @ transform
    )
    eri = np.einsum(
        "ap,bq,cr,ds,abcd->pqrs",
        transform,
        transform,
        transform,
        transform,
        molecule.intor("int2e"),
        optimize=True,
    )
    with h5py.File(native.folder / "df_hf_int/meta.h5") as source:
        require(
            np.array_equal(source["chunk_indices"][:], [0])
            and source["chunk_size"][()] == 1,
            "unknown native DF chunk layout",
        )
    with h5py.File(native.folder / "df_hf_int/VQ_0.h5") as source:
        raw = np.asarray(source["0"])
    require(
        raw.ndim == 4 and raw.shape[0] == 1 and raw.shape[2:] == (2, 4),
        "unknown native DF vertex axes",
    )
    pairs = raw[0].reshape(raw.shape[1], 2, 2, 2)
    vertex = pairs[..., 0] + 1j * pairs[..., 1]
    df = np.empty((2, 2, 2, 2), complex)
    for p, q, r, s in np.ndindex(2, 2, 2, 2):
        df[p, q, r, s] = sum(
            vertex[k, p, q] * vertex[k, r, s].conjugate() for k in range(len(vertex))
        )
    require(
        maximum(df.imag) <= 1e-12 and maximum(native.core.imag) <= 1e-12,
        "bounded H2 reference expects real integrals",
    )
    df = df.real
    _, orbitals = np.linalg.eigh(native.fock[0, 0].real)
    density = 2 * np.outer(orbitals[:, 0], orbitals[:, 0])
    fock = native.core[0, 0].real.copy()
    for p, q in np.ndindex(2, 2):
        for r, s in np.ndindex(2, 2):
            fock[p, q] += density[r, s] * (df[p, q, r, s] - 0.5 * df[p, r, q, s])
    fock_error = maximum(fock - native.fock[0, 0])
    hf_energy = native.nuclear + 0.5 * sum(
        density[p, q] * (native.core[0, 0, q, p].real + fock[q, p])
        for p, q in np.ndindex(2, 2)
    )
    hf_error = abs(hf_energy - native.initial_hf)
    require(max(fock_error, hf_error) <= 3e-9, "independent DF HF reconstruction fails")
    exact = FiniteTrace(core, eri, native.mu, native.nuclear)
    fitted = FiniteTrace(native.core[0, 0].real, df, native.mu, native.nuclear)
    sector_errors = []
    for na, nb in np.ndindex(3, 3):
        states = [
            state
            for state in range(16)
            if sum((state >> (2 * p)) & 1 for p in range(2)) == na
            and sum((state >> (2 * p + 1)) & 1 for p in range(2)) == nb
        ]
        lowest = float(np.linalg.eigvalsh(exact.hamiltonian[np.ix_(states, states)])[0])
        if na + nb:
            solver = fci.direct_spin1.FCI()
            solver.verbose, solver.conv_tol = 0, 1e-13
            expected = float(solver.kernel(core, eri, 2, (na, nb))[0])
        else:
            expected = 0.0
        sector_errors.append(abs(lowest - expected))
    require(
        max(sector_errors) <= 3e-9,
        "independent finite model differs from a PySCF FCI sector",
    )
    hf = scf.RHF(molecule)
    hf.verbose, hf.conv_tol = 0, 1e-13
    exact_hf = float(hf.kernel())
    require(hf.converged, "independent exact-ERI restricted HF did not converge")
    result = {
        "mu_hartree": native.mu,
        "exact_FCI_energy_hartree": exact.energy,
        "DF_FCI_energy_hartree": fitted.energy,
        "GF2_energy_hartree": native.energy,
        "GF2_minus_DF_FCI_energy_hartree": native.energy - fitted.energy,
        "DF_minus_exact_FCI_energy_hartree": fitted.energy - exact.energy,
        "exact_FCI_density": exact.density,
        "DF_FCI_density": fitted.density,
        "sector_energy_error_max": max(sector_errors),
        "density_fitting": {
            "core_integral_error_vs_exact_hartree": maximum(native.core[0, 0] - core),
            "eri_error_vs_exact_hartree": maximum(df - eri),
            "initial_HF_energy_difference_hartree": native.initial_hf - exact_hf,
            "native_DF_HF_fock_reconstruction_error_hartree": fock_error,
            "native_DF_HF_energy_reconstruction_error_hartree": hf_error,
        },
    }
    for label, model in (("exact_FCI", exact), ("DF_FCI", fitted)):
        result[label + "_green_errors"] = {
            "tau": maximum(native.time(tau, "green") - model.time(tau)),
            "iw": maximum(
                native.frequency(indices, "green") - model.frequency(indices)
            ),
        }
    return result


def audit(folder: Path, native_root: Path) -> dict[str, Any]:
    report = read(folder / "report.json")
    require(
        report.get("schema") == "kernelgecp.molecular-green.v2"
        and report.get("pilot") is False,
        "canonical audit refuses a pilot or wrong schema",
    )
    config = report["config"]
    require(
        config["schema_version"] == 2
        and config["validation_tau_rule"] == "gauss_legendre_96"
        and config["regression_tau_rule"] == "native_mod4_3_plus_gauss64"
        and config["retired_validation_source"]
        == "6050ddd33af30deafe0d99057ea5d8a277c411a9",
        "wrong fresh validation protocol",
    )
    require(
        not set(config["validation_indices"])
        & (set(config["calibration_indices"]) | set(config["regression_indices"])),
        "final frequency query overlap",
    )
    require(
        config["target"] <= 1e-6
        and config["calibration_target"] <= 1e-7
        and config["native_refinement_target"] <= 1e-7
        and config["native_dyson_target"] <= 1e-6,
        "EF acceptance budgets were weakened",
    )
    require(
        config["methods"] == list(METHODS)
        and config["bonds_angstrom"] == [0.74, 1.5]
        and config["native_lambdas"] == [1e4, 1e5]
        and config["beta"] == 20
        and config["repeats"] == 3,
        "incomplete EF protocol",
    )
    require(
        not set(config["calibration_indices"]) & set(config["validation_indices"]),
        "calibration/final frequency queries overlap",
    )
    require(
        set(report["source_hashes"]) == SOURCE_FILES
        and re.fullmatch("[0-9a-f]{40}", report["source_revision"]) is not None,
        "generation source inventory/revision invalid",
    )
    for name, digest in report["source_hashes"].items():
        blob = subprocess.check_output(
            ["git", "show", report["source_revision"] + ":" + name], cwd=ROOT
        )
        require(
            hashlib.sha256(blob).hexdigest() == digest,
            "generation source hash mismatch",
        )
    config_blob = subprocess.check_output(
        [
            "git",
            "show",
            report["source_revision"] + ":experiments/molecular_green.toml",
        ],
        cwd=ROOT,
    )
    require(
        hashlib.sha256(config_blob).hexdigest() == report["config_sha256"],
        "committed configuration hash mismatch",
    )
    import tomllib

    require(
        tomllib.loads(config_blob.decode()) == config,
        "reported configuration differs from its frozen source",
    )
    require(
        re.fullmatch("[0-9a-f]{40}", report["native_source_revision"]) is not None
        and sha(native_root / "native-campaign.json")
        == report["native_campaign_sha256"],
        "actual native campaign identity changed",
    )
    native_report = read(native_root / "native-campaign.json")
    require(
        native_report["schema"] == "kernelgecp.green-molecular-native.v1"
        and native_report["source_revision"] == report["native_source_revision"],
        "native campaign/wrapper revision mismatch",
    )
    native_mechanism_hashes = {}
    for name in NATIVE_MECHANISMS:
        blob = subprocess.check_output(
            ["git", "show", report["native_source_revision"] + ":" + name], cwd=ROOT
        )
        native_mechanism_hashes[name] = hashlib.sha256(blob).hexdigest()
    require(len(native_report["cases"]) == 4, "native campaign case count changed")
    natives = {}
    for item in native_report["cases"]:
        case = item["directory"]
        require(case not in natives, "duplicate native case")
        native_folder = relative(native_root, case)
        require(
            sha(native_folder / "bound-manifest.json") == item["manifest_sha256"],
            "native manifest hash changed",
        )
        natives[case] = RawNative(
            native_folder,
            item["bond_length_angstrom"],
            item["native_grid_lambda"],
            report["native_source_revision"],
        )
    require(
        {(n.bond, n.native_lambda) for n in natives.values()}
        == {(bond, cutoff) for bond in (0.74, 1.5) for cutoff in (1e4, 1e5)},
        "native cohort is not complete",
    )
    require(
        sha(folder / "selection.json") == report["selection_sha256"],
        "frozen selection file changed",
    )
    selection = json.loads((folder / "selection.json").read_text())
    require(
        selection == report["selection"] and len(selection) == 12,
        "selection identity/count changed",
    )
    chosen = {(row["case"], row["method"]): row for row in selection}
    require(
        set(chosen) == {(case, method) for case in natives for method in METHODS},
        "selection Cartesian cohort incomplete",
    )
    require(
        len(report["validation"]) == 36
        and len(report["cutoff_refinement"]) == 12
        and len(report["model_differences"]) == 4,
        "required artifacts omitted",
    )
    control_metrics = []
    first_passes = []
    for case, native in natives.items():
        path = folder / case / "calibration.npz"
        require(
            sha(path) == report["reference_artifacts"][case]["calibration.npz"],
            "calibration reference hash changed",
        )
        reference = numeric_archive(path, {"tau", "indices", *FIELDS})
        require(
            np.array_equal(reference["indices"], config["calibration_indices"]),
            "calibration frequency query identity changed",
        )
        for method in METHODS:
            selected = chosen[case, method]
            ladder = [
                (cutoff, eps)
                for cutoff in config["green_cutoffs_hartree"]
                for eps in config["native_eps_candidates"]
            ]
            require(
                (selected["cutoff"], selected["eps"]) in ladder,
                "selection is outside frozen candidates",
            )
            selected_index = ladder.index((selected["cutoff"], selected["eps"]))
            attempts = [
                row
                for row in report["calibration"]
                if row["case"] == case and row["method"] == method
            ]
            require(
                len(attempts) == selected_index + 1,
                "calibration attempt ledger omits a candidate",
            )
            controls = [
                row
                for row in report["cutoff_refinement_attempts"]
                if row["case"] == case and row["method"] == method
            ]
            control_index = 0
            for index, (cutoff, eps) in enumerate(ladder[: selected_index + 1]):
                row = attempts[index]
                expected_dir = f"{case}/{method}/calibration-{cutoff}-{eps}"
                require(
                    row["worker_directory"] == expected_dir,
                    "calibration directory mismatch",
                )

                def check_worker(
                    row: dict[str, Any], directory: str, actual_cutoff: float
                ) -> dict[str, Any] | None:
                    worker = relative(folder, directory)
                    for name, digest in row["artifacts"].items():
                        require(
                            sha(relative(worker, name)) == digest,
                            "calibration artifact hash changed",
                        )
                    request = read(worker / "request.json")
                    require(
                        request["case"] == case
                        and request["method"] == method
                        and request["cutoff"] == actual_cutoff
                        and request["eps"] == eps
                        and request["config"] == config
                        and request["role"] == "calibration",
                        "calibration request identity changed",
                    )
                    if row["status"] != "completed":
                        require(
                            row["status"] == "timeout" or row.get("returncode", 0) != 0,
                            "invalid failed worker status",
                        )
                        return None
                    arrays = numeric_archive(worker / "result.npz", WORKER_KEYS)
                    validate_arrays(
                        arrays,
                        native,
                        config,
                        "calibration",
                        row["green_rank"],
                        row["sigma_rank"],
                    )
                    return worker_metrics(arrays, native)

                calculated = check_worker(row, expected_dir, cutoff)
                primary_pass = (
                    calculated is not None
                    and max(calculated["errors"].values())
                    <= config["calibration_target"]
                )
                refinement_pass = False
                if primary_pass:
                    require(
                        control_index < len(controls), "missing doubled-cutoff attempt"
                    )
                    control = controls[control_index]
                    control_index += 1
                    ref_dir = f"{case}/{method}/refinement-{cutoff}-{eps}"
                    require(
                        control["worker_directory"] == ref_dir,
                        "refinement directory mismatch",
                    )
                    refined = check_worker(control, ref_dir, 2 * cutoff)
                    refinement_pass = (
                        refined is not None
                        and max(refined["errors"].values())
                        <= config["calibration_target"]
                    )
                    if index == selected_index:
                        require(
                            selected["calibration_directory"] == expected_dir
                            and selected["refinement_directory"] == ref_dir,
                            "selected worker directories changed",
                        )
                        selected_controls = [
                            item
                            for item in report["cutoff_refinement"]
                            if item["case"] == case and item["method"] == method
                        ]
                        require(
                            len(selected_controls) == 1
                            and selected_controls[0] == control,
                            "selected control does not match attempt",
                        )
                        control_metrics.append(
                            {"case": case, "method": method, **refined}
                        )
                require(
                    (primary_pass and refinement_pass) == (index == selected_index),
                    "selection is not first passing both calibration gates",
                )
            require(control_index == len(controls), "unexpected doubled-cutoff attempt")
            first_passes.append(
                {"case": case, "method": method, "candidate_index": selected_index}
            )
    seen = set()
    fingerprints = {}
    metrics = []
    bases = {}
    reconstruction_max = 0.0
    for row in report["validation"]:
        case, method, repeat = row["case"], row["method"], row["repeat"]
        require(
            row["status"] == "completed" and row["role"] == "validation",
            "final worker did not complete its declared role",
        )
        require(
            case in natives
            and method in METHODS
            and type(repeat) is int
            and 0 <= repeat < 3,
            "unknown final trial identity",
        )
        key = (case, method, repeat)
        require(key not in seen, "duplicate final trial identity")
        seen.add(key)
        worker = folder / "validation" / f"repeat-{repeat}-{method}-{case}"
        for name, digest in row["artifacts"].items():
            require(
                sha(relative(worker, name)) == digest, "worker artifact hash changed"
            )
        arrays = numeric_archive(worker / "result.npz", WORKER_KEYS)
        native = natives[case]
        picked = chosen[case, method]
        require(
            row["green_cutoff_hartree"] == picked["cutoff"]
            and row["native_eps"] == picked["eps"]
            and row["sigma_cutoff_hartree"]
            == config["selfenergy_cutoff_multiplier"] * picked["cutoff"],
            "final provider differs from frozen selection",
        )
        require(
            row["native_manifest_sha256"] == sha(native.folder / "bound-manifest.json")
            and row["native_iteration"] == native.iteration,
            "final native identity changed",
        )
        validate_arrays(
            arrays, native, config, "validation", row["green_rank"], row["sigma_rank"]
        )
        reference_path = folder / case / "validation.npz"
        require(
            sha(reference_path)
            == report["reference_artifacts"][case]["validation.npz"],
            "saved native reference hash changed",
        )
        reference = numeric_archive(reference_path, {"tau", "indices", *FIELDS})
        require(
            np.array_equal(arrays["tau"], reference["tau"])
            and np.array_equal(arrays["indices"], reference["indices"]),
            "actual final query identities differ",
        )
        calculated = worker_metrics(arrays, native)
        for field in FIELDS:
            require(
                abs(calculated["errors"][field] - row["errors"][field]) < 2e-8,
                "reported error differs from independent native reconstruction",
            )
            require(
                maximum(
                    reference[field]
                    - (
                        native.time(arrays["tau"], field.split("_")[0])
                        if field.endswith("tau")
                        else native.frequency(arrays["indices"], field.split("_")[0])
                    )
                )
                <= 2e-8,
                "saved truth differs from independent native reconstruction",
            )
        require(
            abs(calculated["density"] - complex_number(row["density"])) < 1e-10
            and abs(
                calculated["native_density"] - complex_number(row["native_density"])
            )
            < 1e-10
            and abs(calculated["density_error"] - row["density_error"]) < 1e-10,
            "reported density mismatch",
        )
        require(
            all(
                abs(
                    calculated["causality_imag_max"][field]
                    - row["causality_imag_max"][field]
                )
                < 1e-10
                for field in ("green", "sigma")
            ),
            "reported matrix causality differs",
        )
        actual_fingerprints = {
            name: hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
            for name, value in arrays.items()
        }
        require(
            actual_fingerprints == row["numerical_arrays_sha256"],
            "numerical array fingerprints changed",
        )
        previous = fingerprints.setdefault((case, method), actual_fingerprints)
        require(
            previous == actual_fingerprints,
            "three fresh processes do not reproduce exact arrays",
        )
        if (case, method) not in bases:
            bases[case, method] = {
                field: build_basis(
                    method,
                    config,
                    picked["cutoff"]
                    * (
                        1
                        if field == "green"
                        else config["selfenergy_cutoff_multiplier"]
                    ),
                    picked["eps"],
                )
                for field in ("green", "sigma")
            }
        for field, basis in bases[case, method].items():
            coefficients = arrays[
                "g_coefficients" if field == "green" else "sigma_coefficients"
            ].reshape(basis.rank, 4)
            for domain, query in (("tau", arrays["tau"]), ("iw", arrays["indices"])):
                design = (
                    basis.design_tau(query)
                    if domain == "tau"
                    else basis.design_matsubara(query)
                )
                represented = (np.asarray(design) @ coefficients).reshape(
                    len(query), 1, 1, 2, 2
                )
                reconstruction_max = max(
                    reconstruction_max,
                    maximum(represented - arrays[field + "_" + domain]),
                )
        endpoint = (
            bases[case, method]["green"].design_tau(np.array([0.0, 20.0]))
            @ arrays["g_coefficients"].reshape(row["green_rank"], 4)
        ).reshape(2, 1, 1, 2, 2)
        reconstruction_max = max(
            reconstruction_max, maximum(endpoint - arrays["green_endpoints"])
        )
        metrics.append({"case": case, "method": method, "repeat": repeat, **calculated})
    require(len(seen) == 36, "complete three-repeat cohort missing")
    native_checks = []
    for case, native in natives.items():
        indices = np.asarray(config["validation_indices"], np.int64)
        native_checks.append(
            {
                "case": case,
                "dyson": native.dyson(indices),
                "direct_fourier_discrepancy": max(
                    native.basis_frequency_error(indices, field)
                    for field in ("green", "sigma")
                ),
                "QR_reconstruction_max": max(native.reconstruction_errors.values()),
            }
        )
    grid = []
    for bond in (0.74, 1.5):
        coarse = next(
            n for n in natives.values() if n.bond == bond and n.native_lambda == 1e4
        )
        refined = next(
            n for n in natives.values() if n.bond == bond and n.native_lambda == 1e5
        )
        tau = expected_time(refined, "validation")
        indices = np.asarray(config["validation_indices"], np.int64)
        errors = {
            field + "_" + domain: maximum(
                (coarse.time(tau, field) - refined.time(tau, field))
                if domain == "tau"
                else (
                    coarse.frequency(indices, field) - refined.frequency(indices, field)
                )
            )
            for field in ("green", "sigma")
            for domain in ("tau", "iw")
        }
        grid.append(
            {
                "bond_length_angstrom": bond,
                "errors": errors,
                "energy_difference_hartree": abs(coarse.energy - refined.energy),
                "mu_difference_hartree": abs(coarse.mu - refined.mu),
            }
        )
    models = []
    for case, native in natives.items():
        rows = [row for row in report["model_differences"] if row["case"] == case]
        require(len(rows) == 1, "finite-model identity missing/duplicated")
        actual = model_accounting(
            native,
            expected_time(native, "validation"),
            np.asarray(config["validation_indices"], np.int64),
        )
        reported = rows[0]
        for name, value in actual.items():
            if isinstance(value, dict):
                for key, number in value.items():
                    require(
                        abs(number - reported[name][key]) <= 2e-8,
                        "reported finite-model/DF accounting differs independently",
                    )
            else:
                require(
                    abs(value - reported[name]) <= 2e-8,
                    "reported finite-model energy/density differs independently",
                )
        models.append({"case": case, **actual})
    gates = {
        "representation": all(
            max(row["errors"].values()) <= config["target"] for row in metrics
        ),
        "density": all(row["density_error"] <= config["target"] for row in metrics),
        "matrix_causality": all(
            max(row["causality_imag_max"].values()) <= config["target"]
            for row in metrics
        ),
        "coefficient_prediction_consistency": reconstruction_max <= 1e-9,
        "native_Dyson_and_independent_Fourier": all(
            row["dyson"] <= config["native_dyson_target"]
            and row["direct_fourier_discrepancy"] <= 2e-8
            for row in native_checks
        ),
        "native_grid_refinement": all(
            max(row["errors"].values()) <= config["native_refinement_target"]
            and row["energy_difference_hartree"] <= config["native_refinement_target"]
            for row in grid
        ),
        "doubled_cutoff_controls": len(control_metrics) == 12
        and all(
            max(row["errors"].values()) <= config["calibration_target"]
            for row in control_metrics
        ),
        "finite_model_and_DF_accounting": len(models) == 4,
    }
    return {
        "schema": "kernelgecp.molecular-green-audit.v1",
        "status": "passed" if all(gates.values()) else "failed",
        "generation_revision": report["source_revision"],
        "native_source_revision": report["native_source_revision"],
        "native_mechanism_hashes": native_mechanism_hashes,
        "gates": gates,
        "final_trials": metrics,
        "native": native_checks,
        "native_grid_refinement": grid,
        "cutoff_controls": control_metrics,
        "first_passing_calibration": first_passes,
        "model_differences": models,
        "prediction_reconstruction_max": reconstruction_max,
        "validated_identity_count": len(seen),
        "independent_formulations": [
            "raw_HDF5",
            "overdetermined_pivoted_QR",
            "direct_exp_iw_tau_512_quadrature",
        ],
        "model_scope": "finite_basis_GF2_model_errors_are_separate_no_1e-6_model_claim",
    }


def normal(value: Any) -> Any:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, np.generic):
        return normal(value.item())
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: normal(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [normal(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(args.folder, args.native_root)
    except Exception as error:
        result = {
            "schema": "kernelgecp.molecular-green-audit.v1",
            "status": "failed",
            "error_type": type(error).__name__,
            "error": str(error),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(normal(result), sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
