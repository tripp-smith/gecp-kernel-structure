"""Actual optional CT-HYB measurements and measured-data Bethe bath updates.

The explicitly paramagnetic, single-band problem uses mu=U/2. Measured spin
blocks are retained; only bath feedback is spin averaged. Native density-matrix
moments and TRIQS tail replacement are used identically for native and fitted
frequency values. Monte Carlo uncertainty is not a representation certificate.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import numpy.typing as npt

from ..dmft import PhysicalBasis, _integer, _real
from ..fitting import LeastSquaresPlan, _immutable_copy
from ..types import NumericalBreakdownError
from .triqs import _provider

type ComplexArray = npt.NDArray[np.complex128]
type FloatArray = npt.NDArray[np.float64]
type IndexArray = npt.NDArray[np.int64]


def _owned(value: npt.ArrayLike, dtype: Any = np.complex128) -> Any:
    return _immutable_copy(np.asarray(value, dtype=dtype))


def _array(value: npt.ArrayLike, name: str, shape: tuple[int, ...]) -> ComplexArray:
    raw = np.asarray(value)
    if raw.shape != shape or raw.dtype.kind not in "fci":
        raise ValueError(f"{name} must have finite numeric shape {shape}")
    array = np.asarray(raw, dtype=np.complex128)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite numbers")
    return array


@dataclass(frozen=True, slots=True)
class CTHYBConfig:
    """Bounded serial half-filled paramagnetic CT-HYB Bethe parameters."""

    beta: float = 10.0
    hopping: float = 1.0
    interaction: float = 2.0
    n_iw: int = 256
    n_tau: int = 3073
    mixing: float = 0.5
    warmup_cycles: int = 10000
    length_cycle: int = 50
    tail_fraction: float = 0.2
    rcond: float = 1e-10

    def __post_init__(self) -> None:
        for name in ("beta", "hopping", "mixing", "tail_fraction"):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        for name in ("interaction", "rcond"):
            object.__setattr__(
                self, name, _real(getattr(self, name), name, zero_allowed=True)
            )
        for name, minimum in (
            ("n_iw", 8),
            ("n_tau", 2),
            ("warmup_cycles", 1),
            ("length_cycle", 1),
        ):
            object.__setattr__(self, name, _integer(getattr(self, name), name, minimum))
        if self.n_tau < 6 * self.n_iw + 1:
            raise ValueError("n_tau must be >= 6*n_iw+1 for the native transform")
        if self.n_iw > 4096 or self.n_tau > 49153:
            raise ValueError("CT-HYB grid exceeds the supported bounded allocation")
        if self.mixing > 1 or self.tail_fraction >= 1 or self.rcond >= 1:
            raise ValueError("mixing must be <=1; tail_fraction and rcond must be <1")
        if not np.isfinite(self.beta * self.interaction):
            raise ValueError("beta*interaction must be finite")
        if (
            not np.isfinite(self.hopping * self.hopping)
            or self.hopping * self.hopping == 0
        ):
            raise ValueError("hopping squared must be finite and positive")


@dataclass(frozen=True, slots=True)
class CTHYBMeasurement:
    """One independent native Markov chain; no standard error is inferred."""

    config: CTHYBConfig
    seed: int
    cycles: int
    bath_iw: ComplexArray
    g_tau: ComplexArray
    g_iw: ComplexArray
    g_iw_raw: ComplexArray
    sigma_iw: ComplexArray
    g_moments: ComplexArray
    sigma_moments: ComplexArray
    occupancies: FloatArray
    double_occupancy: float
    average_sign: float
    solve_status: int

    def __post_init__(self) -> None:
        if not isinstance(self.config, CTHYBConfig):
            raise TypeError("config must be CTHYBConfig")
        object.__setattr__(self, "seed", _integer(self.seed, "seed", 1))
        object.__setattr__(self, "cycles", _integer(self.cycles, "cycles", 1))
        if self.seed >= 2**31:
            raise ValueError("seed must fit a positive signed 32-bit integer")
        config = self.config
        for name, shape in (
            ("bath_iw", (2 * config.n_iw,)),
            ("g_tau", (config.n_tau, 2)),
            ("g_iw", (2 * config.n_iw, 2)),
            ("g_iw_raw", (2 * config.n_iw, 2)),
            ("sigma_iw", (2 * config.n_iw, 2)),
            ("g_moments", (3, 2)),
            ("sigma_moments", (2, 2)),
        ):
            object.__setattr__(
                self, name, _owned(_array(getattr(self, name), name, shape))
            )
        occ = _array(self.occupancies, "occupancies", (2,))
        if (
            np.max(np.abs(occ.imag)) > 1e-12
            or np.any(occ.real < 0)
            or np.any(occ.real > 1)
        ):
            raise ValueError("occupancies must be real probabilities")
        object.__setattr__(self, "occupancies", _owned(occ.real, np.float64))
        for name in ("double_occupancy", "average_sign"):
            value = _real(getattr(self, name), name, zero_allowed=True)
            if value > 1:
                raise ValueError(f"{name} must lie in [0,1]")
            object.__setattr__(self, name, value)
        if type(self.solve_status) is not int or self.solve_status != 0:
            raise ValueError(
                "only successfully completed native solves are measurements"
            )
        if self.double_occupancy > min(self.occupancies) + 1e-12:
            raise ValueError("double occupancy exceeds a single-spin occupation")
        if not np.allclose(self.g_moments[:2], [[0, 0], [1, 1]], atol=1e-12, rtol=0):
            raise ValueError("native Green moments must start with zero and unit mass")

    @property
    def density(self) -> float:
        """Total impurity occupation, obtained from the native density matrix."""
        return float(np.sum(self.occupancies))


@dataclass(frozen=True, slots=True)
class CTHYBRepresentation:
    """A fit to one frozen measurement, with explicit native tail handling."""

    g_tau: ComplexArray
    g_iw: ComplexArray
    coefficients: ComplexArray
    rank: int
    retained_rank: int
    condition_estimate: float | None
    training_residual: float
    validation_residual: float
    native_iw_discrepancy: float
    native_tail_change: float


@dataclass(frozen=True, slots=True)
class ReplicaEstimate:
    """Independent-chain mean/covariance and standard error of that mean."""

    mean: FloatArray
    covariance: FloatArray
    standard_error: FloatArray
    replicas: int


def estimate_replicas(observations: npt.ArrayLike) -> ReplicaEstimate:
    """Sample covariance across independent chain means, divided by N for SE.

    This function cannot establish independence from the values themselves.
    Callers must use disjoint seeds and separately thermalized Markov chains.
    It deliberately refuses a single chain and complex-valued observable axes.
    """
    raw = np.asarray(observations)
    if raw.dtype.kind not in "fi" or raw.ndim != 2 or min(raw.shape) < 1:
        raise ValueError("replica observations must be a real two-dimensional array")
    array = np.asarray(raw, dtype=np.float64)
    if array.shape[0] < 2 or not np.all(np.isfinite(array)):
        raise ValueError("at least two finite independent replica means are required")
    mean = array.mean(axis=0)
    centered = array - mean
    with np.errstate(over="ignore", invalid="ignore"):
        covariance = centered.T @ centered / (array.shape[0] - 1)
        error = np.sqrt(np.diag(covariance) / array.shape[0])
    if not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(error)):
        raise NumericalBreakdownError("replica covariance overflow")
    return ReplicaEstimate(
        _owned(mean, np.float64),
        _owned(covariance, np.float64),
        _owned(error, np.float64),
        array.shape[0],
    )


_DEFAULT_CONFIG = CTHYBConfig()


class CTHYBPlan:
    """Cached actual-provider meshes and optional common time-design fit.

    Constructing the plan requires the pinned native runtime. Calling sample
    always creates a fresh Solver and thermalizes it; previous chain states
    cannot silently couple fixed-bath replicas.
    """

    def __init__(
        self,
        config: CTHYBConfig = _DEFAULT_CONFIG,
        basis: PhysicalBasis | None = None,
        *,
        fit_tau_indices: npt.ArrayLike | None = None,
    ) -> None:
        if not isinstance(config, CTHYBConfig):
            raise TypeError("config must be CTHYBConfig")
        self._config = config
        self._gfs, self._meshes = _provider()
        try:
            native = importlib.import_module("triqs_cthyb")
            version = importlib.import_module("triqs_cthyb.version")
            self._operators = importlib.import_module("triqs.operators")
            self._atom = importlib.import_module("triqs.atom_diag")
            mpi = importlib.import_module("triqs.utility.mpi")
        except ImportError as error:
            raise ImportError(
                "CT-HYB 4.0.1 is required for this optional integration"
            ) from error
        if getattr(version, "version", None) != "4.0.1":
            raise ImportError(
                "this integration supports the pinned CT-HYB 4.0.1 runtime"
            )
        if mpi.size != 1:
            raise ValueError("this bounded integration requires one MPI rank")
        self._solver_class = native.Solver
        self._indices = _owned(np.arange(-config.n_iw, config.n_iw), np.int64)
        self._tau = _owned(np.linspace(0, config.beta, config.n_tau), np.float64)
        self._mesh = self._meshes.MeshImFreq(
            beta=config.beta, statistic="Fermion", n_iw=config.n_iw
        )
        self._mesh.set_tail_fit_parameters(config.tail_fraction)
        self._basis = basis
        self._fit: LeastSquaresPlan | None = None
        if fit_tau_indices is None:
            observed = np.arange(config.n_tau, dtype=np.int64)
        else:
            raw = np.asarray(fit_tau_indices)
            if raw.ndim != 1 or raw.dtype.kind not in "iu" or not raw.size:
                raise ValueError("fit_tau_indices must be nonempty integer positions")
            if np.any(raw < 0) or np.any(raw >= config.n_tau):
                raise ValueError("fit_tau_indices lie outside the native time mesh")
            observed = np.asarray(raw, dtype=np.int64)
            if np.unique(observed).size != observed.size:
                raise ValueError("fit_tau_indices must be distinct")
        self._observed = _owned(observed, np.int64)
        self._heldout = _owned(
            np.setdiff1d(np.arange(config.n_tau), observed), np.int64
        )
        self._tau_design: ComplexArray | None = None
        self._iw_design: ComplexArray | None = None
        if basis is not None:
            if basis.beta != config.beta or basis.rank < 1:
                raise ValueError("basis must have matching beta and a positive rank")
            tau_design = _array(
                basis.design_tau(self._tau), "time design", (config.n_tau, basis.rank)
            )
            iw_design = _array(
                basis.design_matsubara(self._indices),
                "frequency design",
                (2 * config.n_iw, basis.rank),
            )
            if observed.size < basis.rank:
                raise ValueError("fitting requires at least rank measured time rows")
            self._tau_design = _owned(tau_design)
            self._iw_design = _owned(iw_design)
            self._fit = LeastSquaresPlan(tau_design[observed], rcond=config.rcond)
            if self._fit.diagnostics.effective_rank == 0:
                raise ValueError(
                    "the measured-data representation has zero retained rank"
                )

    @property
    def config(self) -> CTHYBConfig:
        return self._config

    @property
    def indices(self) -> IndexArray:
        return cast(IndexArray, self._indices)

    @property
    def tau(self) -> FloatArray:
        return cast(FloatArray, self._tau)

    def initial_bath(self) -> ComplexArray:
        """Independent analytic noninteracting semicircular Bethe bath."""
        iw = 1j * (2 * self.indices + 1) * np.pi / self.config.beta
        root = np.sqrt(iw * iw - 4 * self.config.hopping**2)
        root = np.where(np.imag(root) * np.imag(iw) < 0, -root, root)
        return cast(ComplexArray, _owned(2 * self.config.hopping**2 / (iw + root)))

    def _bath(self, bath_iw: npt.ArrayLike) -> ComplexArray:
        bath = _array(bath_iw, "bath_iw", (2 * self.config.n_iw,))
        if np.max(np.abs(bath - bath[::-1].conj())) > 1e-10:
            raise ValueError(
                "the physical scalar bath must satisfy frequency adjoint symmetry"
            )
        if np.max(bath[self.config.n_iw :].imag) > 1e-10:
            raise ValueError("the physical bath must be causal at positive frequencies")
        return bath

    def sample(
        self, bath_iw: npt.ArrayLike, *, cycles: int, seed: int
    ) -> CTHYBMeasurement:
        """Run a fresh actual native chain with density-matrix moments."""
        bath = self._bath(bath_iw)
        cycles, seed = _integer(cycles, "cycles", 1), _integer(seed, "seed", 1)
        if seed >= 2**31:
            raise ValueError("seed must fit a positive signed 32-bit integer")
        config = self.config
        solver = self._solver_class(
            beta=config.beta,
            gf_struct=[("up", 1), ("down", 1)],
            n_iw=config.n_iw,
            n_tau=config.n_tau,
        )
        iw = 1j * (2 * self.indices + 1) * np.pi / config.beta
        for _, g0 in solver.G0_iw:
            g0.data[:, 0, 0] = 1 / (iw + config.interaction / 2 - bath)
        up, down = self._operators.n("up", 0), self._operators.n("down", 0)
        solver.solve(
            h_int=config.interaction * up * down,
            n_cycles=cycles,
            n_warmup_cycles=config.warmup_cycles,
            length_cycle=config.length_cycle,
            random_seed=seed,
            measure_G_tau=True,
            measure_density_matrix=True,
            use_norm_as_weight=True,
            perform_post_proc=True,
            perform_tail_fit=False,
            verbosity=0,
        )
        if int(solver.solve_status) != 0:
            raise NumericalBreakdownError("native CT-HYB did not complete successfully")
        if solver.G_moments is None or solver.Sigma_moments is None:
            raise NumericalBreakdownError("native density-matrix moments are missing")
        blocks = ("up", "down")
        g_tau = np.column_stack([solver.G_tau[b].data[:, 0, 0] for b in blocks])
        raw = np.column_stack([solver.G_iw_raw[b].data[:, 0, 0] for b in blocks])
        moments = np.column_stack([solver.G_moments[b][:, 0, 0] for b in blocks])
        sigma_moments = np.column_stack(
            [solver.Sigma_moments[b][:, 0, 0] for b in blocks]
        )
        density_matrix, diagonalization = (
            solver.density_matrix,
            solver.h_loc_diagonalization,
        )
        occ = [
            self._atom.trace_rho_op(density_matrix, op, diagonalization)
            for op in (up, down)
        ]
        double = self._atom.trace_rho_op(density_matrix, up * down, diagonalization)
        if (
            max(abs(complex(x).imag) for x in (*occ, double, solver.average_sign))
            > 1e-12
        ):
            raise NumericalBreakdownError("nonreal static native observables")
        return CTHYBMeasurement(
            config,
            seed,
            cycles,
            bath,
            g_tau,
            self._native_tail(raw, moments),
            raw,
            np.column_stack([solver.Sigma_iw[b].data[:, 0, 0] for b in blocks]),
            moments,
            sigma_moments,
            np.array([complex(x).real for x in occ]),
            float(complex(double).real),
            float(complex(solver.average_sign).real),
            int(solver.solve_status),
        )

    def _native_tail(self, values: ComplexArray, moments: ComplexArray) -> ComplexArray:
        output = np.empty_like(values)
        for spin in range(2):
            gf = self._gfs.Gf(mesh=self._mesh.copy(), target_shape=[1, 1])
            gf.data[:, 0, 0] = values[:, spin]
            gf.replace_by_tail_in_fit_window(moments[:, spin, None, None])
            output[:, spin] = gf.data[:, 0, 0]
        if not np.all(np.isfinite(output)):
            raise NumericalBreakdownError("nonfinite native tail replacement")
        return cast(ComplexArray, _owned(output))

    def represent(self, measurement: CTHYBMeasurement) -> CTHYBRepresentation:
        """Fit only the declared measured rows; validate on excluded rows."""
        if measurement.config != self.config:
            raise ValueError("measurement config differs from this plan")
        if self._fit is None or self._tau_design is None or self._iw_design is None:
            raise ValueError("representation requires a basis supplied to the plan")
        coefficients = np.asarray(
            self._fit.fit(measurement.g_tau[self._observed]), dtype=complex
        )
        tau = self._tau_design @ coefficients
        raw = self._iw_design @ coefficients
        iw = self._native_tail(raw, measurement.g_moments)
        train = float(
            np.max(np.abs(tau[self._observed] - measurement.g_tau[self._observed]))
        )
        heldout = (
            float(np.max(np.abs(tau[self._heldout] - measurement.g_tau[self._heldout])))
            if self._heldout.size
            else 0.0
        )
        diagnostics = self._fit.diagnostics
        return CTHYBRepresentation(
            _owned(tau),
            _owned(iw),
            _owned(coefficients),
            diagnostics.design_shape[1],
            diagnostics.effective_rank,
            diagnostics.retained_condition_number,
            train,
            heldout,
            float(np.max(np.abs(iw - measurement.g_iw))),
            float(np.max(np.abs(iw - raw))),
        )

    def bath_update(
        self,
        measurement: CTHYBMeasurement,
        representation: CTHYBRepresentation | None = None,
    ) -> tuple[ComplexArray, float]:
        """Return the proposed mixed bath and unmixed full-grid defect.

        A defect from one noisy chain is a controller diagnostic, not a
        self-consistency or physical-accuracy certificate.
        """
        if measurement.config != self.config:
            raise ValueError("measurement config differs from this plan")
        values = measurement.g_iw if representation is None else representation.g_iw
        values = _array(values, "represented g_iw", (2 * self.config.n_iw, 2))
        candidate = self.config.hopping**2 * values.mean(axis=1)
        defect = float(np.max(np.abs(candidate - measurement.bath_iw)))
        result = measurement.bath_iw + self.config.mixing * (
            candidate - measurement.bath_iw
        )
        self._bath(result)
        return cast(ComplexArray, _owned(result)), defect

    def observables(
        self,
        measurement: CTHYBMeasurement,
        indices: npt.ArrayLike,
        *,
        representation: CTHYBRepresentation | None = None,
    ) -> FloatArray:
        """Total density, double occupancy, then real/imag spin-mean G(iw)."""
        if measurement.config != self.config:
            raise ValueError("measurement config differs from this plan")
        labels = np.asarray(indices)
        if labels.ndim != 1 or labels.dtype.kind not in "iu" or not labels.size:
            raise ValueError("selected indices must be a nonempty integer vector")
        if np.any(labels < -self.config.n_iw) or np.any(labels >= self.config.n_iw):
            raise ValueError("selected indices lie outside the measured grid")
        values = measurement.g_iw if representation is None else representation.g_iw
        values = _array(values, "observable Green values", (2 * self.config.n_iw, 2))
        selected = values[labels + self.config.n_iw].mean(axis=1)
        return cast(
            FloatArray,
            _owned(
                np.r_[
                    measurement.density,
                    measurement.double_occupancy,
                    selected.real,
                    selected.imag,
                ],
                np.float64,
            ),
        )
