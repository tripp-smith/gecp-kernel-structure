"""Optional dense TRIQS Fourier reference for half-filled Bethe IPT.

IPT is an approximate impurity model. Native Fourier grids, interpolation and
fitted tail diagnostics are numerical evidence, not exact Hubbard physics or
continuous error certificates. TRIQS 4.0.2 is loaded only when a plan is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, cast

import numpy as np
import numpy.typing as npt
from scipy.interpolate import CubicSpline

from ..basis import _indices_array, _tau_array
from ..dmft import IPTConfig, IPTIteration
from ..types import FloatArray, NumericalBreakdownError
from .triqs import _iw_indices, _provider, _tau_coordinates

type ComplexArray = npt.NDArray[np.complex128]
type IntArray = npt.NDArray[np.int64]
type DenseField = Literal["green", "sigma", "weiss"]


def _frozen(values: npt.ArrayLike, dtype: npt.DTypeLike) -> Any:
    array = np.ascontiguousarray(values, dtype=dtype)
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


def _complex_vector(values: npt.ArrayLike, size: int, name: str) -> ComplexArray:
    original = np.asarray(values)
    if original.dtype.kind not in "fci" or original.shape != (size,):
        raise ValueError(
            f"{name} must be a numeric vector matching the full native mesh"
        )
    result = np.array(original, dtype=np.complex128, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values")
    return result


def _moments(config: IPTConfig, field: DenseField) -> tuple[float, ...]:
    if field == "green":
        return (0.0, 1.0, 0.0, config.hopping**2 + config.interaction**2 / 4)
    if field == "weiss":
        return (0.0, 1.0, 0.0, config.hopping**2)
    if field == "sigma":
        return (0.0, config.interaction**2 / 4, 0.0)
    raise ValueError("field must be green, sigma or weiss")


@dataclass(frozen=True, slots=True)
class DenseIPTUpdate:
    """Owning native update at one supplied bath, before any bath mixing."""

    config: IPTConfig
    indices: IntArray
    tau: FloatArray
    bath_iw: ComplexArray
    weiss_iw: ComplexArray
    weiss_tau: ComplexArray
    sigma_iw: ComplexArray
    sigma_tau: ComplexArray
    g_iw: ComplexArray
    g_tau: ComplexArray

    def __post_init__(self) -> None:
        if not isinstance(self.config, IPTConfig):
            raise TypeError("config must be IPTConfig")
        indices = np.asarray(self.indices)
        if indices.dtype.kind not in "iu" or not np.array_equal(
            indices, np.arange(-self.config.n_iw, self.config.n_iw)
        ):
            raise ValueError("indices must be the complete signed ordinary native mesh")
        object.__setattr__(self, "indices", _frozen(indices, np.int64))
        tau = _tau_array(self.tau, self.config.beta)
        if (
            tau.size != self.config.n_tau
            or tau[0] != 0
            or tau[-1] != self.config.beta
            or not np.all(np.diff(tau) > 0)
        ):
            raise ValueError("tau must be the ordered complete native time mesh")
        object.__setattr__(self, "tau", _frozen(tau, np.float64))
        for name in ("bath_iw", "weiss_iw", "sigma_iw", "g_iw"):
            array = _complex_vector(getattr(self, name), indices.size, name)
            object.__setattr__(self, name, _frozen(array, np.complex128))
        for name in ("weiss_tau", "sigma_tau", "g_tau"):
            array = _complex_vector(getattr(self, name), tau.size, name)
            object.__setattr__(self, name, _frozen(array, np.complex128))

    @property
    def sigma_tau_product(self) -> ComplexArray:
        """The native reference retains the unprojected reflected product."""
        return self.sigma_tau

    def evaluate_tau(
        self, tau: npt.ArrayLike, *, field: DenseField = "green"
    ) -> ComplexArray:
        """Interpolate native time samples with a not-a-knot cubic spline.

        The caller must verify reference refinement for its queries. Endpoints
        remain 0+ and beta-; interpolation does not create a continuous bound.
        """
        query = _tau_array(tau, self.config.beta)
        _moments(self.config, field)
        source = {
            "green": self.g_tau,
            "sigma": self.sigma_tau,
            "weiss": self.weiss_tau,
        }[field]
        return cast(ComplexArray, np.asarray(CubicSpline(self.tau, source)(query)))

    def evaluate_matsubara(
        self, indices: npt.ArrayLike, *, field: DenseField = "green"
    ) -> ComplexArray:
        """Read native samples, or Fourier-extend using declared model moments.

        Extension is limited by the existing time mesh; it never silently
        refines the reference. No values are filled by conjugation.
        """
        query = _indices_array(indices)
        moments = _moments(self.config, field)
        source_iw = {
            "green": self.g_iw,
            "sigma": self.sigma_iw,
            "weiss": self.weiss_iw,
        }[field]
        if np.all((query >= -self.config.n_iw) & (query < self.config.n_iw)):
            return np.array(source_iw[query + self.config.n_iw], copy=True)
        count = max(self.config.n_iw, int(np.max(query)) + 1, -int(np.min(query)))
        if self.config.n_tau < 6 * count + 1:
            raise ValueError("refine the native time mesh before extending this far")
        gfs, meshes = _provider()
        output = gfs.Gf(
            mesh=meshes.MeshImFreq(
                beta=self.config.beta, statistic="Fermion", n_iw=count
            ),
            target_shape=[],
        )
        source = gfs.Gf(
            mesh=meshes.MeshImTime(
                beta=self.config.beta, statistic="Fermion", n_tau=self.config.n_tau
            ),
            target_shape=[],
        )
        source.data[:] = {
            "green": self.g_tau,
            "sigma": self.sigma_tau,
            "weiss": self.weiss_tau,
        }[field]
        output.set_from_fourier(source, known_moments=np.array(moments, complex))
        result = np.array(output.data[query + count], dtype=np.complex128, copy=True)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("native Fourier extension is nonfinite")
        return result


@dataclass(frozen=True, slots=True)
class NativeTailDiagnostics:
    """Applied PH model moments and separately fitted numerical diagnostics."""

    weiss_known_moments: tuple[float, ...]
    green_known_moments: tuple[float, ...]
    sigma_known_moments: tuple[float, ...]
    sigma_fitted_first_moment: float | None
    sigma_first_moment_error: float | None
    sigma_tail_fit_error: float | None
    sigma_tail_fit_status: str
    density_imaginary: float
    positive_frequency_imag_g_max: float
    positive_frequency_imag_sigma_max: float
    particle_hole_g_error: float
    particle_hole_sigma_error: float
    time_grid_recommendation_met: bool
    arithmetic: str = "triqs_native_complex128_fourier"
    moment_scope: str = "half_filled_paramagnetic_bethe_IPT_model_assumptions"
    density_scope: str = "native_moment_constrained_per_spin_density"
    certification_status: str = "not_certified"


@dataclass(frozen=True, slots=True)
class DenseIPTResult(DenseIPTUpdate):
    """Native candidate and actual input bath at the final tested iteration."""

    history: tuple[IPTIteration, ...]
    converged: bool
    stop_reason: str
    density: float
    native_tail: NativeTailDiagnostics
    trajectory: tuple[ComplexArray, ...] = ()

    def __post_init__(self) -> None:
        DenseIPTUpdate.__post_init__(self)
        object.__setattr__(self, "history", tuple(self.history))
        object.__setattr__(
            self,
            "trajectory",
            tuple(
                _frozen(
                    _complex_vector(bath, self.indices.size, "trajectory bath"),
                    np.complex128,
                )
                for bath in self.trajectory
            ),
        )


class DenseIPTPlan:
    """Cached native meshes and Fourier buffers for the explicit PH model."""

    @property
    def config(self) -> IPTConfig:
        return self._config

    @property
    def indices(self) -> IntArray:
        return self._indices

    @property
    def tau(self) -> FloatArray:
        return self._tau

    def __init__(self, config: IPTConfig) -> None:
        if not isinstance(config, IPTConfig):
            raise TypeError("config must be IPTConfig")
        self._config = config
        self._gfs, self._meshes = _provider()
        iw_mesh = self._meshes.MeshImFreq(
            beta=config.beta, statistic="Fermion", n_iw=config.n_iw
        )
        tau_mesh = self._meshes.MeshImTime(
            beta=config.beta, statistic="Fermion", n_tau=config.n_tau
        )
        self._indices: IntArray = _frozen(
            _iw_indices(iw_mesh, config.beta, self._meshes), np.int64
        )
        self._tau: FloatArray = _frozen(
            _tau_coordinates(tau_mesh, config.beta, self._meshes), np.float64
        )
        self._iw = 1j * (2 * self.indices + 1) * np.pi / config.beta
        self._weiss_iw = self._gfs.Gf(mesh=iw_mesh.copy(), target_shape=[])
        self._sigma_iw = self._gfs.Gf(mesh=iw_mesh.copy(), target_shape=[])
        self._green_iw = self._gfs.Gf(mesh=iw_mesh.copy(), target_shape=[])
        self._weiss_tau = self._gfs.Gf(mesh=tau_mesh.copy(), target_shape=[])
        self._sigma_tau = self._gfs.Gf(mesh=tau_mesh.copy(), target_shape=[])
        self._green_tau = self._gfs.Gf(mesh=tau_mesh.copy(), target_shape=[])
        self._weiss_moments = np.array(_moments(config, "weiss"), complex)
        self._sigma_moments = np.array(_moments(config, "sigma"), complex)
        self._green_moments = np.array(_moments(config, "green"), complex)
        if not all(
            np.all(np.isfinite(array))
            for array in (
                self._iw,
                self._weiss_moments,
                self._sigma_moments,
                self._green_moments,
            )
        ):
            raise ValueError(
                "native frequencies and model moments must be finite in binary64"
            )

    def initial_bath(self) -> ComplexArray:
        """Independent causal noninteracting semicircle, Delta=t^2 G."""
        w = self._iw.imag
        scale = np.maximum(np.abs(w), 2 * self.config.hopping)
        normalized = np.abs(w) / scale
        magnitude = (
            2 / (normalized + np.hypot(normalized, 2 * self.config.hopping / scale))
        ) / scale
        result = np.asarray(-1j * np.sign(w) * self.config.hopping**2 * magnitude)
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("native initial bath is nonfinite")
        return result

    def update(self, bath_iw: npt.ArrayLike) -> DenseIPTUpdate:
        """Compute one native update; do not mix or alter the supplied bath."""
        bath = _complex_vector(bath_iw, self.indices.size, "bath_iw")
        with np.errstate(all="ignore"):
            self._weiss_iw.data[:] = 1 / (self._iw - bath)
        if not np.all(np.isfinite(self._weiss_iw.data)):
            raise NumericalBreakdownError("native Weiss inversion is nonfinite")
        self._weiss_tau.set_from_fourier(
            self._weiss_iw, known_moments=self._weiss_moments
        )
        self._sigma_tau.data[:] = (
            self.config.interaction**2
            * self._weiss_tau.data**2
            * self._weiss_tau.data[::-1]
        )
        if self.config.interaction == 0:
            self._sigma_iw.data[:] = 0
        else:
            self._sigma_iw.set_from_fourier(
                self._sigma_tau, known_moments=self._sigma_moments
            )
        with np.errstate(all="ignore"):
            self._green_iw.data[:] = 1 / (self._iw - bath - self._sigma_iw.data)
        if not all(
            np.all(np.isfinite(gf.data))
            for gf in (self._weiss_tau, self._sigma_tau, self._sigma_iw, self._green_iw)
        ):
            raise NumericalBreakdownError("native IPT update is nonfinite")
        self._green_tau.set_from_fourier(
            self._green_iw, known_moments=self._green_moments
        )
        return DenseIPTUpdate(
            config=self.config,
            indices=self.indices,
            tau=self.tau,
            bath_iw=bath,
            weiss_iw=self._weiss_iw.data,
            weiss_tau=self._weiss_tau.data,
            sigma_iw=self._sigma_iw.data,
            sigma_tau=self._sigma_tau.data,
            g_iw=self._green_iw.data,
            g_tau=self._green_tau.data,
        )

    def _diagnostics(
        self, update: DenseIPTUpdate
    ) -> tuple[float, NativeTailDiagnostics]:
        self._green_iw.data[:] = update.g_iw
        density = complex(
            self._gfs.density(self._green_iw, known_moments=self._green_moments)
        )
        if not np.isfinite(density):
            raise NumericalBreakdownError("native density is nonfinite")
        fitted = error = moment_error = None
        status = "zero_interaction_exact_zero"
        if self.config.interaction:
            self._sigma_iw.data[:] = update.sigma_iw
            try:
                moments, native_error = self._gfs.fit_hermitian_tail(
                    self._sigma_iw, known_moments=np.array([0], complex)
                )
                first = complex(moments[1])
                if not np.isfinite(first) or not np.isfinite(native_error):
                    raise ValueError("nonfinite native tail diagnostic")
                fitted, error = float(first.real), float(native_error)
                moment_error = abs(first - self.config.interaction**2 / 4)
                status = "native_unconstrained_first_moment_fit"
            except (RuntimeError, ValueError, IndexError):
                status = "native_fit_unavailable"
        native_tail = NativeTailDiagnostics(
            weiss_known_moments=_moments(self.config, "weiss"),
            green_known_moments=_moments(self.config, "green"),
            sigma_known_moments=_moments(self.config, "sigma"),
            sigma_fitted_first_moment=fitted,
            sigma_first_moment_error=moment_error,
            sigma_tail_fit_error=error,
            sigma_tail_fit_status=status,
            density_imaginary=float(density.imag),
            positive_frequency_imag_g_max=float(
                np.max(update.g_iw[self.config.n_iw :].imag)
            ),
            positive_frequency_imag_sigma_max=float(
                np.max(update.sigma_iw[self.config.n_iw :].imag)
            ),
            particle_hole_g_error=float(
                np.max(np.abs(update.g_iw + update.g_iw[::-1]))
            ),
            particle_hole_sigma_error=float(
                np.max(np.abs(update.sigma_iw + update.sigma_iw[::-1]))
            ),
            time_grid_recommendation_met=self.config.n_tau >= 6 * self.config.n_iw + 1,
        )
        return float(density.real), native_tail

    def solve(self, *, record_trajectory: bool = False) -> DenseIPTResult:
        """Stop only after the unmixed full-grid defect is small twice.

        The last proposed mixed step is recorded but not applied on either
        convergence or the iteration budget. bath_iw is the actual input to
        the returned candidate, so its Dyson equation remains explicit.
        """
        if not isinstance(record_trajectory, bool):
            raise TypeError("record_trajectory must be boolean")
        bath = self.initial_bath()
        history = []
        trajectory = []
        consecutive = 0
        for iteration in range(1, self.config.max_iterations + 1):
            if record_trajectory:
                trajectory.append(_frozen(bath, np.complex128))
            update = self.update(bath)
            step = self.config.hopping**2 * update.g_iw - bath
            defect = float(np.max(np.abs(step)))
            history.append(
                IPTIteration(
                    iteration=iteration,
                    unmixed_defect=defect,
                    mixed_step=float(np.max(np.abs(self.config.mixing * step))),
                    full_grid_unmixed_defect=defect,
                )
            )
            consecutive = consecutive + 1 if defect < self.config.tolerance else 0
            if consecutive >= 2 or iteration == self.config.max_iterations:
                break
            bath = bath + self.config.mixing * step
        density, native_tail = self._diagnostics(update)
        return DenseIPTResult(
            config=update.config,
            indices=update.indices,
            tau=update.tau,
            bath_iw=update.bath_iw,
            weiss_iw=update.weiss_iw,
            weiss_tau=update.weiss_tau,
            sigma_iw=update.sigma_iw,
            sigma_tau=update.sigma_tau,
            g_iw=update.g_iw,
            g_tau=update.g_tau,
            history=tuple(history),
            converged=consecutive >= 2,
            stop_reason="tolerance" if consecutive >= 2 else "max_iterations",
            density=density,
            native_tail=native_tail,
            trajectory=tuple(trajectory),
        )


def solve_dense_ipt(
    config: IPTConfig, *, record_trajectory: bool = False
) -> DenseIPTResult:
    """Run the optional native reference without importing TRIQS eagerly."""
    return DenseIPTPlan(config).solve(record_trajectory=record_trajectory)
