"""Finite spectral and Fock-space references, independent of an IPT solver.

These are binary64 closed finite-model calculations, not a certificate for an
interacting continuum bath. Spectral weights are never silently normalized:
canonical Green functions have unit mass, whereas a dynamic self-energy does not.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from numbers import Real
from typing import cast

import numpy as np
import numpy.typing as npt
from scipy.linalg import eigh, eigvalsh, expm
from scipy.special import expit, logsumexp

from .fitting import (
    NumericArray,
    _contains_boolean_sequence,
    _immutable_copy,
    _numeric_array,
)
from .types import FloatArray, NumericalBreakdownError

type ComplexArray = npt.NDArray[np.complex128]
_MAX_FOCK_MODES = 4
_MAX_IPT_TERMS = 2**24


def _real(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be finite in binary64") from error
    if not math.isfinite(result) or (positive and result <= 0):
        qualifier = " and positive" if positive else ""
        raise ValueError(f"{name} must be finite{qualifier}")
    return result


def _nonnegative(value: object, name: str) -> float:
    result = _real(value, name)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _natural(value: object, name: str, upper: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer in [0, {upper}]")
    try:
        result = operator.index(value)  # type: ignore[arg-type]
    except TypeError as error:
        raise ValueError(f"{name} must be an integer in [0, {upper}]") from error
    if not 0 <= result <= upper:
        raise ValueError(f"{name} must be an integer in [0, {upper}]")
    return result


def _vector(value: npt.ArrayLike, name: str) -> FloatArray:
    array = _numeric_array(value, name)
    if array.ndim != 1 or np.iscomplexobj(array):
        raise ValueError(f"{name} must be a one-dimensional real array")
    return cast(FloatArray, array)


def _times(tau: npt.ArrayLike, beta: float) -> FloatArray:
    values = _vector(tau, "tau")
    if np.any(values < 0) or np.any(values > beta):
        raise ValueError("tau must lie in [0, beta]")
    return values


def _indices(value: npt.ArrayLike) -> npt.NDArray[np.int64]:
    if _contains_boolean_sequence(value):
        raise ValueError("Matsubara indices must not contain booleans")
    array = np.asarray(value)
    if array.ndim == 1 and array.size == 0:
        return np.array([], dtype=np.int64)
    if array.ndim != 1 or array.dtype.kind not in "iu":
        raise ValueError("Matsubara indices must be a one-dimensional integer array")
    limit = 2**51 - 1
    if np.any(array > limit) or (array.dtype.kind == "i" and np.any(array < -limit)):
        raise ValueError("Matsubara index magnitude must be smaller than 2**51")
    return np.asarray(array, dtype=np.int64)


def _factors(beta: float, energies: FloatArray) -> tuple[FloatArray, FloatArray]:
    # Compute both sides directly: 1 - expit(x) loses small thermal populations.
    with np.errstate(over="ignore", under="ignore"):
        scaled = beta * energies
        return np.asarray(expit(-scaled)), np.asarray(expit(scaled))


@dataclass(frozen=True, slots=True, init=False, eq=False)
class FermionicPoleModel:
    """Owned nonnegative scalar spectral measure in physical energy units.

    ``G(tau)=-sum(w*K(tau,E))`` and ``G(iw)=sum(w/(iw-E))``. Energies are
    stored in stable ascending order, including repeated poles. Zero weights
    and zero total mass are allowed, for example for a zero self-energy.
    ``cutoff`` optionally validates a supplied support bound; it never clips
    poles. This model does not impose normalization or infer physical truth.
    """

    beta: float
    energies: FloatArray
    weights: FloatArray
    declared_cutoff: float | None
    _mass: float

    def __init__(
        self,
        beta: float,
        energies: npt.ArrayLike,
        weights: npt.ArrayLike,
        *,
        cutoff: float | None = None,
    ) -> None:
        beta_value = _real(beta, "beta", positive=True)
        nodes = _vector(energies, "energies")
        masses = _vector(weights, "weights")
        if not nodes.size or nodes.shape != masses.shape or np.any(masses < 0):
            raise ValueError("matching nonempty poles and nonnegative weights required")
        support = None if cutoff is None else _nonnegative(cutoff, "cutoff")
        if support is not None and np.any(np.abs(nodes) > support):
            raise ValueError("a pole lies outside the declared cutoff")
        try:
            mass = math.fsum(float(weight) for weight in masses)
        except OverflowError as error:
            raise ValueError("total pole mass must be finite in binary64") from error
        if not math.isfinite(mass):
            raise ValueError("total pole mass must be finite in binary64")
        order = np.argsort(nodes, kind="stable")
        object.__setattr__(self, "beta", beta_value)
        object.__setattr__(self, "energies", _immutable_copy(nodes[order]))
        object.__setattr__(self, "weights", _immutable_copy(masses[order]))
        object.__setattr__(self, "declared_cutoff", support)
        object.__setattr__(self, "_mass", mass)

    @property
    def mass(self) -> float:
        """Spectral zeroth moment; never implicitly replaced by one."""
        return self._mass

    @property
    def spectral_cutoff(self) -> float:
        """Largest supplied absolute pole energy, including zero-weight poles."""
        return float(np.max(np.abs(self.energies)))

    def evaluate_tau(self, tau: npt.ArrayLike) -> FloatArray:
        """Evaluate both 0+ and beta- directly with nonpositive exponents."""
        times = _times(tau, self.beta)
        occupied, empty = _factors(self.beta, self.energies)
        positive = self.energies >= 0
        rows = np.empty((times.size, self.energies.size), dtype=np.float64)
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            rows[:, positive] = (
                np.exp(-times[:, None] * self.energies[positive]) * empty[positive]
            )
            rows[:, ~positive] = (
                np.exp((self.beta - times[:, None]) * self.energies[~positive])
                * occupied[~positive]
            )
            result = -(rows @ self.weights)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("finite-pole time evaluation is nonfinite")
        return cast(FloatArray, _immutable_copy(result))

    def evaluate_matsubara(self, indices: npt.ArrayLike) -> ComplexArray:
        """Evaluate ordinary signed integer n, with wn=(2*n+1)*pi/beta."""
        modes = _indices(indices)
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            frequencies = (2 * modes + 1) * np.pi / self.beta
            if not np.all(np.isfinite(frequencies)):
                raise NumericalBreakdownError("Matsubara frequencies are nonfinite")
            result = (1 / (1j * frequencies[:, None] - self.energies)) @ self.weights
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("pole frequency evaluation is nonfinite")
        return cast(ComplexArray, _immutable_copy(result))

    def density(self) -> float:
        """Occupied mass; a particle count only for a canonical Green function."""
        occupied, _ = _factors(self.beta, self.energies)
        return math.fsum(
            float(w) * float(f) for w, f in zip(self.weights, occupied, strict=True)
        )

    def moments(self, max_order: int = 2) -> FloatArray:
        """Return immutable moments M0 through M[max_order], for orders 0..16."""
        order = _natural(max_order, "max_order", 16)
        values = [self.mass]
        try:
            for power in range(1, order + 1):
                values.append(
                    math.fsum(
                        float(w) * float(e) ** power
                        for e, w in zip(self.energies, self.weights, strict=True)
                        if w != 0
                    )
                )
        except (OverflowError, ValueError) as error:
            raise NumericalBreakdownError("spectral moments exceed binary64") from error
        result = np.asarray(values, dtype=np.float64)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("spectral moments exceed binary64")
        return cast(FloatArray, _immutable_copy(result))


def fock_annihilation_operators(modes: int) -> tuple[FloatArray, ...]:
    """Canonical annihilation matrices, low-bit-first ordering, for 1..4 modes."""
    count = _natural(modes, "modes", _MAX_FOCK_MODES)
    if count == 0:
        raise ValueError("at least one Fock mode is required")
    dimension = 2**count
    operators = []
    for mode in range(count):
        matrix = np.zeros((dimension, dimension), dtype=np.float64)
        for state in range(dimension):
            if state & (1 << mode):
                sign = -1 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1
                matrix[state ^ (1 << mode), state] = sign
        operators.append(cast(FloatArray, _immutable_copy(matrix)))
    return tuple(operators)


def _fock_inputs(
    hamiltonian: npt.ArrayLike, annihilation: npt.ArrayLike
) -> tuple[NumericArray, NumericArray]:
    h = _numeric_array(hamiltonian, "hamiltonian")
    c = _numeric_array(annihilation, "annihilation")
    if h.ndim != 2 or h.shape[0] != h.shape[1] or h.shape[0] not in (2, 4, 8, 16):
        raise ValueError("Fock Hamiltonian must be square with dimension 2, 4, 8 or 16")
    if c.shape != h.shape or not np.array_equal(h, h.conj().T):
        raise ValueError("matching annihilator and Hermitian Hamiltonian required")
    with np.errstate(all="ignore"):
        canonical = np.allclose(
            c @ c.conj().T + c.conj().T @ c, np.eye(len(h)), rtol=1e-12, atol=1e-12
        )
        nilpotent = np.allclose(c @ c, 0, rtol=0, atol=1e-12)
    if not canonical or not nilpotent:
        raise ValueError("annihilation must satisfy canonical anticommutation")
    return _immutable_copy(h), _immutable_copy(c)


def finite_fock_lehmann(
    beta: float, hamiltonian: npt.ArrayLike, annihilation: npt.ArrayLike
) -> FermionicPoleModel:
    """Build an untruncated scalar Lehmann measure from a <=16-state Hamiltonian."""
    temperature = _real(beta, "beta", positive=True)
    h, c = _fock_inputs(hamiltonian, annihilation)
    energies, vectors = eigh(h)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        exponents = -temperature * (energies - energies[0])
        probabilities = np.exp(exponents - logsumexp(exponents))
        amplitudes = vectors.conj().T @ c @ vectors
        poles = energies[None, :] - energies[:, None]
        weights = (probabilities[:, None] + probabilities[None, :]) * (
            np.abs(amplitudes) ** 2
        )
    if not np.all(np.isfinite(poles)) or not np.all(np.isfinite(weights)):
        raise NumericalBreakdownError("Fock Lehmann construction is nonfinite")
    selected = weights > 0
    if not np.any(selected):
        raise NumericalBreakdownError("canonical Fock spectrum has no positive mass")
    return FermionicPoleModel(temperature, poles[selected], weights[selected])


def thermal_trace_tau(
    beta: float,
    hamiltonian: npt.ArrayLike,
    annihilation: npt.ArrayLike,
    tau: npt.ArrayLike,
) -> FloatArray:
    """Independent matrix-exponential trace, without Lehmann eigenvectors/weights."""
    temperature = _real(beta, "beta", positive=True)
    h, c = _fock_inputs(hamiltonian, annihilation)
    times = _times(tau, temperature)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        shifted = h - float(eigvalsh(h)[0]) * np.eye(len(h))
        if not np.all(np.isfinite(shifted)):
            raise NumericalBreakdownError("ground-state energy shift is nonfinite")
        thermal = expm(-temperature * shifted)
        partition = complex(np.trace(thermal))
        if (
            not math.isfinite(partition.real)
            or partition.real <= 0
            or not math.isfinite(partition.imag)
            or abs(partition.imag) > 1e-10 * partition.real
        ):
            raise NumericalBreakdownError("thermal partition trace is inconclusive")
        values = []
        for time in times:
            left = expm(-(temperature - float(time)) * shifted)
            right = expm(-float(time) * shifted)
            values.append(-complex(np.trace(left @ c @ right @ c.conj().T)) / partition)
    complex_values = np.asarray(values, dtype=np.complex128)
    if not np.all(np.isfinite(complex_values)) or np.any(
        np.abs(complex_values.imag) > 1e-10 * (1 + np.abs(complex_values.real))
    ):
        raise NumericalBreakdownError("thermal matrix trace is inconclusive")
    return cast(FloatArray, _immutable_copy(complex_values.real))


def hubbard_atom(beta: float, interaction: float) -> FermionicPoleModel:
    """Exact PH atom for H=U*(n_up-1/2)*(n_down-1/2), from all four Fock states."""
    u = _nonnegative(interaction, "interaction")
    up, down = fock_annihilation_operators(2)
    identity = np.eye(4)
    h = u * ((up.T @ up - 0.5 * identity) @ (down.T @ down - 0.5 * identity))
    return finite_fock_lehmann(beta, h, up)


def impurity_one_bath(
    beta: float, interaction: float, hybridization: float
) -> FermionicPoleModel:
    """Exact PH impurity plus one zero-energy bath site per spin (four modes).

    Mode order is impurity up/down, bath up/down. The hybridization is real;
    its sign is retained, not interpreted as a change of bath model.
    """
    u = _nonnegative(interaction, "interaction")
    v = _real(hybridization, "hybridization")
    up, down, bath_up, bath_down = fock_annihilation_operators(4)
    identity = np.eye(16)
    h = u * ((up.T @ up - 0.5 * identity) @ (down.T @ down - 0.5 * identity))
    h += v * (up.T @ bath_up + bath_up.T @ up + down.T @ bath_down + bath_down.T @ down)
    return finite_fock_lehmann(beta, h, up)


def ipt_self_energy_poles(
    weiss: FermionicPoleModel, interaction: float
) -> FermionicPoleModel:
    """Independent reflected second-order triple sum for a normalized Weiss G0.

    Pole energies are Ea+Eb-Ec. Positive weights use occupied/unoccupied
    factors separately; their sum is U^2*n0*(1-n0), rather than one. This is
    perturbation theory of the supplied Weiss function, not the exact finite-
    bath interacting self-energy. The finite oracle permits <=2^24 triples.
    """
    if not isinstance(weiss, FermionicPoleModel):
        raise TypeError("weiss must be a FermionicPoleModel")
    u = _nonnegative(interaction, "interaction")
    if not math.isclose(weiss.mass, 1, rel_tol=0, abs_tol=1e-12):
        raise ValueError("a normalized unit-mass Weiss Green function is required")
    if weiss.energies.size**3 > _MAX_IPT_TERMS:
        raise ValueError("finite IPT oracle exceeds its 2^24 triple budget")
    if u == 0:
        return FermionicPoleModel(weiss.beta, [0.0], [0.0])
    occupied, empty = _factors(weiss.beta, weiss.energies)
    e, w = weiss.energies, weiss.weights
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        poles = e[:, None, None] + e[None, :, None] - e[None, None, :]
        thermal = (
            empty[:, None, None] * empty[None, :, None] * occupied[None, None, :]
            + occupied[:, None, None] * occupied[None, :, None] * empty[None, None, :]
        )
        weights = (
            (u * u) * w[:, None, None] * w[None, :, None] * w[None, None, :] * thermal
        )
    if not np.all(np.isfinite(poles)) or not np.all(np.isfinite(weights)):
        raise NumericalBreakdownError("finite-pole IPT construction is nonfinite")
    selected = weights > 0
    if not np.any(selected):
        return FermionicPoleModel(weiss.beta, [0.0], [0.0])
    return FermionicPoleModel(weiss.beta, poles[selected], weights[selected])
