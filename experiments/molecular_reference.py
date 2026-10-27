"""Bounded four-mode finite-temperature molecular reference for the EF campaign.

Exact-ERI FCI is a finite-basis model reference, evaluated numerically in
binary64. It is not equivalent to Green's density-fitted GF2 approximation.
PySCF is imported only by the explicit H2 factory; core matrix operations need
only the ordinary package dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
from typing import Any

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from kernelgecp.basis import _indices_array, _tau_array
from kernelgecp.fitting import _numeric_array
from kernelgecp.lehmann import fock_annihilation_operators
from kernelgecp.types import NumericalBreakdownError

MODE_ORDER = ("orbital0_up", "orbital0_down", "orbital1_up", "orbital1_down")


def _positive(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} is not representable in binary64") from error
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive real")
    return result


def _owned(value: Any) -> np.ndarray:
    array = np.ascontiguousarray(value)
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


def _real_array(value: Any, shape: tuple[int, ...], name: str) -> np.ndarray:
    original = _numeric_array(value, name)
    if np.iscomplexobj(original) or original.shape != shape:
        raise ValueError(f"{name} must be a real array of shape {shape}")
    result = np.asarray(original, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values")
    return result


def molecular_hamiltonian(h_core: Any, eri: Any) -> tuple[np.ndarray, np.ndarray]:
    """Build the chemists-ERI Hamiltonian with interleaved spin ordering.

    H = h[p,q] c†[p,s] c[q,s] + (pq|rs)/2 c†[p,a] c†[r,b] c[s,b] c[q,a].
    Only roundoff-sized Hermitian defects are symmetrized, explicitly here.
    No density fitting, frozen-core terms or chemical-potential shift is added.
    """
    one = _real_array(h_core, (2, 2), "h_core")
    two = _real_array(eri, (2, 2, 2, 2), "eri")
    if not np.allclose(one, one.T, rtol=0, atol=1e-12):
        raise ValueError("one-body integrals must be symmetric")
    for transpose in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        if not np.allclose(two, two.transpose(transpose), rtol=0, atol=1e-12):
            raise ValueError("ERIs must obey real chemists' permutation symmetries")
    operators = fock_annihilation_operators(4)
    hamiltonian = np.zeros((16, 16))
    number = sum(c.T @ c for c in operators)
    for p, q in np.ndindex(2, 2):
        for spin in range(2):
            hamiltonian += (
                one[p, q] * operators[2 * p + spin].T @ operators[2 * q + spin]
            )
    for p, q, r, s in np.ndindex(2, 2, 2, 2):
        for alpha, beta in np.ndindex(2, 2):
            hamiltonian += (
                0.5
                * two[p, q, r, s]
                * (
                    operators[2 * p + alpha].T
                    @ operators[2 * r + beta].T
                    @ operators[2 * s + beta]
                    @ operators[2 * q + alpha]
                )
            )
    defect = float(np.max(np.abs(hamiltonian - hamiltonian.T)))
    if not np.all(np.isfinite(hamiltonian)) or defect > 1e-11:
        raise NumericalBreakdownError("finite molecular Hamiltonian is unresolved")
    return _owned((hamiltonian + hamiltonian.T) / 2), _owned(number)


def _thermal(hamiltonian: np.ndarray, number: np.ndarray, beta: float, mu: float):
    with np.errstate(all="ignore"):
        grand = hamiltonian - mu * number
    if not np.all(np.isfinite(grand)):
        raise NumericalBreakdownError("grand-canonical Hamiltonian is nonfinite")
    try:
        values, vectors = np.linalg.eigh(grand)
    except np.linalg.LinAlgError as error:
        raise NumericalBreakdownError(
            "grand-canonical eigensolve did not converge"
        ) from error
    with np.errstate(all="ignore"):
        exponents = -beta * (values - values[0])
        raw = np.exp(exponents)
        weights = raw / raw.sum()
    if not np.all(np.isfinite(weights)):
        raise NumericalBreakdownError("finite grand-canonical partition is unresolved")
    counts = np.real(np.diag(vectors.T @ number @ vectors))
    return values, vectors, weights, counts


def _neutral_mu(hamiltonian: np.ndarray, number: np.ndarray, beta: float) -> float:
    scale = max(1.0, float(np.linalg.norm(hamiltonian, ord=2)))

    def equation(mu):
        _, _, weights, counts = _thermal(hamiltonian, number, beta, mu)
        return float(weights @ counts - 2)

    return float(brentq(equation, -4 * scale, 4 * scale, xtol=1e-14, rtol=1e-14))


@dataclass(frozen=True, slots=True)
class MolecularReference:
    beta: float
    mu: float
    bond_length_angstrom: float | None
    h_core_lowdin: np.ndarray
    eri_lowdin: np.ndarray
    overlap_ao: np.ndarray
    lowdin_x: np.ndarray
    hamiltonian: np.ndarray
    number_operator: np.ndarray
    grand_energies: np.ndarray
    thermal_weights: np.ndarray
    spectral_energies: np.ndarray
    spectral_weights: np.ndarray
    nuclear_energy_hartree: float
    exact_hf_energy_hartree: float | None
    thermal_energy_hartree: float
    electron_count: float
    sector_energy_errors: tuple[tuple[int, int, float], ...]
    provider: str = "numpy_eigh_four_mode_exact_ERI_model"
    certification_status: str = "not_certified_binary64_finite_model"

    def __post_init__(self) -> None:
        for name in (
            "h_core_lowdin",
            "eri_lowdin",
            "overlap_ao",
            "lowdin_x",
            "hamiltonian",
            "number_operator",
            "grand_energies",
            "thermal_weights",
            "spectral_energies",
            "spectral_weights",
        ):
            object.__setattr__(self, name, _owned(getattr(self, name)))
        object.__setattr__(
            self, "sector_energy_errors", tuple(self.sector_energy_errors)
        )

    def evaluate_spin_tau(self, tau: Any) -> np.ndarray:
        times = _tau_array(tau, self.beta)
        energies = self.spectral_energies
        positive = energies >= 0
        kernel = np.empty((len(times), len(energies)))
        with np.errstate(all="ignore"):
            kernel[:, positive] = -np.exp(-times[:, None] * energies[positive]) * expit(
                self.beta * energies[positive]
            )
            kernel[:, ~positive] = -np.exp(
                (self.beta - times[:, None]) * energies[~positive]
            ) * expit(-self.beta * energies[~positive])
            result = np.einsum("tw,wpq->tpq", kernel, self.spectral_weights)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("nonfinite molecular imaginary-time values")
        return _owned(result)

    def evaluate_spin_matsubara(self, indices: Any) -> np.ndarray:
        labels = _indices_array(indices)
        z = 1j * (2.0 * labels + 1) * (np.pi / self.beta)
        with np.errstate(all="ignore"):
            design = 1 / (z[:, None] - self.spectral_energies)
            result = np.einsum("nw,wpq->npq", design, self.spectral_weights)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("nonfinite molecular Matsubara values")
        return _owned(result)

    @staticmethod
    def _restricted(values: np.ndarray) -> np.ndarray:
        up, down = values[:, ::2, ::2], values[:, 1::2, 1::2]
        if (
            max(
                np.max(np.abs(up - down), initial=0),
                np.max(np.abs(values[:, ::2, 1::2]), initial=0),
                np.max(np.abs(values[:, 1::2, ::2]), initial=0),
            )
            > 1e-11
        ):
            raise NumericalBreakdownError("reference is not a restricted spin scalar")
        return _owned(((up + down) / 2)[:, None, None])

    def evaluate_tau(self, tau: Any) -> np.ndarray:
        """Explicit spin average; output axes are time, spin=1, k=1, p, q."""
        return self._restricted(self.evaluate_spin_tau(tau))

    def evaluate_matsubara(self, indices: Any) -> np.ndarray:
        """Explicit spin average; output uses ordinary integer Matsubara labels."""
        return self._restricted(self.evaluate_spin_matsubara(indices))


def reference_from_integrals(
    h_core: Any,
    eri: Any,
    *,
    beta: float = 20.0,
    mu: float | None = None,
    nuclear_energy: float = 0.0,
) -> MolecularReference:
    temperature = _positive(beta, "beta")
    for value, name in ((nuclear_energy, "nuclear_energy"), (mu, "mu")):
        if value is None and name == "mu":
            continue
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise ValueError(f"{name} must be a finite real")
        try:
            valid = np.isfinite(float(value))
        except (ValueError, OverflowError) as error:
            raise ValueError(f"{name} must be finite in binary64") from error
        if not valid:
            raise ValueError(f"{name} must be finite")
    hamiltonian, number = molecular_hamiltonian(h_core, eri)
    if mu is None:
        chemical = _neutral_mu(hamiltonian, number, temperature)
    elif (
        isinstance(mu, (bool, np.bool_))
        or not isinstance(mu, Real)
        or not np.isfinite(mu)
    ):
        raise ValueError("mu must be a finite real or None")
    else:
        chemical = float(mu)
    energies, vectors, weights, counts = _thermal(
        hamiltonian, number, temperature, chemical
    )
    transformed = np.stack(
        [vectors.T @ c @ vectors for c in fock_annihilation_operators(4)]
    )
    amplitudes = transformed.transpose(1, 2, 0)
    measure = (weights[:, None] + weights[None, :])[:, :, None, None] * (
        amplitudes[:, :, :, None] * amplitudes[:, :, None, :].conj()
    )
    poles = energies[None, :] - energies[:, None]
    if not np.all(np.isfinite(poles)) or not np.all(np.isfinite(measure)):
        raise NumericalBreakdownError("molecular spectral measure is nonfinite")
    electronic = float(weights @ (energies + chemical * counts))
    if not np.isfinite(electronic + nuclear_energy):
        raise NumericalBreakdownError("nonfinite molecular thermal energy")
    return MolecularReference(
        beta=temperature,
        mu=chemical,
        bond_length_angstrom=None,
        h_core_lowdin=np.asarray(h_core, float),
        eri_lowdin=np.asarray(eri, float),
        overlap_ao=np.eye(2),
        lowdin_x=np.eye(2),
        hamiltonian=hamiltonian,
        number_operator=number,
        grand_energies=energies,
        thermal_weights=weights,
        spectral_energies=poles.ravel(),
        spectral_weights=measure.reshape(-1, 4, 4),
        nuclear_energy_hartree=float(nuclear_energy),
        exact_hf_energy_hartree=None,
        thermal_energy_hartree=electronic + nuclear_energy,
        electron_count=float(weights @ counts),
        sector_energy_errors=(),
    )


def build_h2_reference(
    bond_length_angstrom: float, *, beta: float = 20.0, mu: float | None = None
) -> MolecularReference:
    """Actual PySCF 2.14.0 H2/STO-3G integrals and all nine FCI spin sectors.

    Exact four-index ERIs are transformed by symmetric Lowdin, without density
    fitting. None mu chooses neutral grand-canonical expectation N=2; a supplied
    native GF2 mu is retained exactly, with its own FCI particle count reported.
    """
    from dataclasses import replace

    distance, temperature = (
        _positive(bond_length_angstrom, "bond length"),
        _positive(beta, "beta"),
    )
    try:
        import pyscf
        from pyscf import fci, gto, scf
    except ImportError as error:
        raise ImportError(
            "the molecular FCI reference requires optional PySCF 2.14.0"
        ) from error
    if pyscf.__version__ != "2.14.0":
        raise RuntimeError("the supported molecular reference requires PySCF 2.14.0")
    mol = gto.M(
        atom=f"H 0 0 0; H 0 0 {distance:.17g}",
        basis="sto-3g",
        unit="Angstrom",
        spin=0,
        charge=0,
        verbose=0,
        max_memory=512,
    )
    mean_field = scf.RHF(mol)
    mean_field.conv_tol = 1e-12
    mean_field.max_cycle = 200
    hf_energy = float(mean_field.kernel())
    if not mean_field.converged or not np.isfinite(hf_energy):
        raise NumericalBreakdownError("PySCF exact-ERI restricted HF did not converge")
    overlap = mol.intor("int1e_ovlp")
    eigenvalues, vectors = np.linalg.eigh(overlap)
    if eigenvalues[0] <= 0:
        raise NumericalBreakdownError("AO overlap is not positive definite")
    lowdin = (vectors / np.sqrt(eigenvalues)) @ vectors.T
    one = lowdin.T @ mean_field.get_hcore() @ lowdin
    two = np.einsum(
        "ap,bq,cr,ds,abcd->pqrs",
        lowdin,
        lowdin,
        lowdin,
        lowdin,
        mol.intor("int2e", aosym="s1"),
        optimize=True,
    )
    result = reference_from_integrals(
        one, two, beta=temperature, mu=mu, nuclear_energy=mol.energy_nuc()
    )
    errors = []
    for up, down in np.ndindex(3, 3):
        states = [
            state
            for state in range(16)
            if sum((state >> mode) & 1 for mode in (0, 2)) == up
            and sum((state >> mode) & 1 for mode in (1, 3)) == down
        ]
        expected = (
            np.linalg.eigvalsh(result.hamiltonian[np.ix_(states, states)])
            + mol.energy_nuc()
        )
        actual, _ = fci.direct_spin1.kernel(
            one,
            two,
            2,
            (up, down),
            nroots=len(states),
            ecore=mol.energy_nuc(),
            verbose=0,
        )
        error = float(np.max(np.abs(np.sort(np.atleast_1d(actual)) - expected)))
        if not np.isfinite(error) or error > 3e-10:
            raise NumericalBreakdownError(
                f"PySCF FCI sector mismatch at {(up, down)}: {error}"
            )
        errors.append((up, down, error))
    return replace(
        result,
        bond_length_angstrom=distance,
        overlap_ao=overlap,
        lowdin_x=lowdin,
        exact_hf_energy_hartree=hf_energy,
        sector_energy_errors=tuple(errors),
        provider="PySCF2.14.0_exact_ERI_all_FCI_sectors",
    )
