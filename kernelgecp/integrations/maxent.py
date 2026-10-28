"""Covariance-preserving positive-spectrum data and optional actual TauMaxEnt.

Projection, covariance truncation, scalar marginal continuation and spectral
mass normalization are explicit operations. No interpolated point becomes an
independent measurement. Optional providers are imported only when invoked.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import math
import re
import zipfile
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path
from typing import Any, Protocol, cast

import numpy as np
import numpy.typing as npt
from scipy.special import expit

from ..basis import _positive_real, _tau_array
from ..fitting import (
    LeastSquaresDiagnostics,
    LeastSquaresPlan,
    _immutable_copy,
    _numeric_array,
)
from ..types import FloatArray, NumericalBreakdownError


def _real_array(value: npt.ArrayLike, name: str) -> FloatArray:
    result = _numeric_array(value, name)
    if np.iscomplexobj(result):
        raise ValueError(f"{name} must be real; no implicit real-part projection")
    return cast(FloatArray, result)


def _owned(value: npt.ArrayLike) -> FloatArray:
    return cast(FloatArray, _immutable_copy(np.asarray(value, dtype=np.float64)))


def _jsonable(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _jsonable(asdict(value))
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _document(value: Any) -> str:
    return json.dumps(_jsonable(value), sort_keys=True, allow_nan=False)


@dataclass(frozen=True, slots=True)
class CovarianceModes:
    """Declared retained noise modes; discarded modes carry no coverage claim."""

    eigenvalues: FloatArray
    retained_values: FloatArray
    retained_vectors: FloatArray
    threshold: float
    relative_cutoff: float
    absolute_cutoff: float
    symmetry_correction: float

    def __post_init__(self) -> None:
        for name in ("eigenvalues", "retained_values", "retained_vectors"):
            object.__setattr__(self, name, _owned(getattr(self, name)))

    @property
    def effective_rank(self) -> int:
        return int(self.retained_values.size)

    @property
    def retained_covariance(self) -> FloatArray:
        value = (self.retained_vectors * self.retained_values) @ self.retained_vectors.T
        return _owned((value + value.T) / 2)

    def whiten(self, residual: npt.ArrayLike) -> FloatArray:
        values = _real_array(residual, "residual")
        if values.shape != (self.retained_vectors.shape[0],):
            raise ValueError("residual must match the original observation coordinates")
        transformed = self.retained_vectors.T @ values / np.sqrt(self.retained_values)
        if not np.all(np.isfinite(transformed)):
            raise NumericalBreakdownError("whitened residual is nonfinite")
        return _owned(transformed)

    def whitened_rms(self, residual: npt.ArrayLike) -> float:
        return float(
            np.linalg.norm(self.whiten(residual)) / math.sqrt(self.effective_rank)
        )


def covariance_modes(
    covariance: npt.ArrayLike, *, rtol: float = 1e-10, atol: float = 1e-18
) -> tuple[FloatArray, CovarianceModes]:
    """Check/symmetrize roundoff, then retain eigenvalues strictly above cutoff."""
    matrix = _real_array(covariance, "covariance")
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or not 1 <= len(matrix) <= 2048
    ):
        raise ValueError(
            "covariance must be a nonempty square matrix of at most2048 rows"
        )
    relative, absolute = (
        _positive_real(rtol, "covariance rtol"),
        _positive_real(atol, "covariance atol"),
    )
    if relative >= 1:
        raise ValueError("covariance rtol must lie in (0,1)")
    scale = float(np.max(np.abs(matrix), initial=0))
    correction = float(np.max(np.abs(matrix - matrix.T), initial=0))
    slack = 128 * np.finfo(float).eps * len(matrix) * scale
    if correction > slack:
        raise ValueError("covariance must be symmetric within roundoff")
    symmetric = matrix / 2 + matrix.T / 2
    try:
        values, vectors = np.linalg.eigh(symmetric)
    except np.linalg.LinAlgError as error:
        raise NumericalBreakdownError(
            "covariance eigensolve did not converge"
        ) from error
    if not np.all(np.isfinite(values)) or values[0] < -slack:
        raise ValueError("covariance is materially indefinite or nonfinite")
    threshold = max(absolute, relative * float(values[-1]))
    retained = values > threshold
    if not np.any(retained):
        raise ValueError("no covariance eigenmodes survive the declared cutoffs")
    modes = CovarianceModes(
        values,
        values[retained],
        vectors[:, retained],
        threshold,
        relative,
        absolute,
        correction,
    )
    return _owned(symmetric), modes


@dataclass(frozen=True, slots=True, init=False, eq=False)
class SpectralDataset:
    """Real scalar G(tau) observations with the full covariance of their mean.

    Tau includes distinct0+ and beta- endpoints. Empirical covariance must
    already be divided by the number of independent chain means. The caller
    records that provenance; independence cannot be inferred from values.
    """

    beta: float
    tau: FloatArray
    values: FloatArray
    covariance: FloatArray
    modes: CovarianceModes
    provenance_json: str
    projection: FloatArray | None
    source_fingerprint: str | None

    def __init__(
        self,
        beta: float,
        tau: npt.ArrayLike,
        values: npt.ArrayLike,
        covariance: npt.ArrayLike,
        *,
        covariance_rtol: float = 1e-10,
        covariance_atol: float = 1e-18,
        provenance: Any = None,
        projection: npt.ArrayLike | None = None,
        source_fingerprint: str | None = None,
    ) -> None:
        temperature = _positive_real(beta, "beta")
        points = _tau_array(tau, temperature)
        observed = _real_array(values, "observations")
        if (
            points.size < 3
            or points[0] != 0
            or points[-1] != temperature
            or np.any(np.diff(points) <= 0)
            or observed.shape != points.shape
        ):
            raise ValueError(
                "ordered observations must include distinct0+ and beta- endpoints"
            )
        matrix, modes = covariance_modes(
            covariance, rtol=covariance_rtol, atol=covariance_atol
        )
        if matrix.shape != (points.size, points.size):
            raise ValueError("covariance must match all observation coordinates")
        mapping = None if projection is None else _real_array(projection, "projection")
        if mapping is not None and (
            mapping.ndim != 2 or mapping.shape[0] != points.size
        ):
            raise ValueError("projection rows must match output observations")
        object.__setattr__(self, "beta", temperature)
        object.__setattr__(self, "tau", _owned(points))
        object.__setattr__(self, "values", _owned(observed))
        object.__setattr__(self, "covariance", matrix)
        object.__setattr__(self, "modes", modes)
        if provenance is not None and not isinstance(provenance, dict):
            raise ValueError("provenance must be a JSON object")
        if source_fingerprint is not None and not re.fullmatch(
            r"[0-9a-f]{64}", source_fingerprint
        ):
            raise ValueError("source_fingerprint must be a SHA-256 digest")
        object.__setattr__(
            self, "provenance_json", _document({} if provenance is None else provenance)
        )
        object.__setattr__(
            self, "projection", None if mapping is None else _owned(mapping)
        )
        object.__setattr__(self, "source_fingerprint", source_fingerprint)

    @property
    def provenance(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self.provenance_json))

    @property
    def fingerprint(self) -> str:
        digest = hashlib.sha256()
        for array in (
            np.array(
                [self.beta, self.modes.relative_cutoff, self.modes.absolute_cutoff]
            ),
            self.tau,
            self.values,
            self.covariance,
        ):
            digest.update(str(array.shape).encode())
            digest.update(np.ascontiguousarray(array, dtype="<f8").tobytes())
        return digest.hexdigest()


def _array_hash(array: FloatArray) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(array, dtype="<f8").tobytes()
    ).hexdigest()


def save_spectral_dataset(
    data: SpectralDataset,
    path: str | Path,
    *,
    units: dict[str, str],
    provider: dict[str, Any],
    state_ordering: tuple[str, ...] = ("scalar",),
    calculation_settings: dict[str, Any] | None = None,
) -> Path:
    """Write owning numeric NPZ plus versioned hash-bound application JSON.

    Both destinations must be new. Units and provider identity are explicit
    caller declarations; the archive does not infer physical conventions.
    """
    target = Path(path)
    if target.suffix != ".npz":
        raise ValueError("spectral archive path must end in .npz")
    manifest_path = target.with_suffix(".json")
    if target.exists() or manifest_path.exists():
        raise FileExistsError("spectral archive or application manifest already exists")
    if (
        not units
        or any(not isinstance(v, str) or not v for v in units.values())
        or not provider
        or not state_ordering
        or len(set(state_ordering)) != len(state_ordering)
    ):
        raise ValueError("archive requires units, provider identity and ordered states")
    arrays = {"tau": data.tau, "values": data.values, "covariance": data.covariance}
    if data.projection is not None:
        arrays["projection"] = data.projection
    manifest = {
        "schema": "kernelgecp.spectral-dataset.v1",
        "observation_domain": "imaginary_time",
        "beta": data.beta,
        "covariance_rtol": data.modes.relative_cutoff,
        "covariance_atol": data.modes.absolute_cutoff,
        "effective_rank": data.modes.effective_rank,
        "fingerprint": data.fingerprint,
        "source_fingerprint": data.source_fingerprint,
        "provenance": data.provenance,
        "units": units,
        "provider": provider,
        "state_ordering": state_ordering,
        "calculation_settings": {}
        if calculation_settings is None
        else calculation_settings,
        "arrays": {
            name: {"shape": array.shape, "sha256": _array_hash(array)}
            for name, array in arrays.items()
        },
    }
    _document(manifest)  # Validate JSON before writing either destination.
    with target.open("xb") as stream:
        np.savez_compressed(stream, **cast(dict[str, Any], arrays))
    manifest["file_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    with manifest_path.open("x") as stream:
        stream.write(_document(manifest) + "\n")
    return manifest_path


def _unique_json(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON manifest key")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> Any:
    raise ValueError(f"nonfinite JSON constant {value}")


def load_spectral_dataset(path: str | Path) -> SpectralDataset:
    """Read only bounded, hash-bound, pickle-free spectral application data."""
    target = Path(path)
    if target.suffix != ".npz" or target.stat().st_size > 128 * 1024**2:
        raise ValueError("unsupported or oversized spectral archive")
    manifest_path = target.with_suffix(".json")
    if manifest_path.stat().st_size > 1024**2:
        raise ValueError("oversized spectral application manifest")
    manifest = json.loads(
        manifest_path.read_text(),
        object_pairs_hook=_unique_json,
        parse_constant=_reject_json_constant,
    )
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema") != "kernelgecp.spectral-dataset.v1"
        or manifest.get("observation_domain") != "imaginary_time"
        or manifest.get("file_sha256")
        != hashlib.sha256(target.read_bytes()).hexdigest()
    ):
        raise ValueError("unsupported spectral manifest or archive hash mismatch")
    expected = {"tau", "values", "covariance"}
    if "projection" in manifest.get("arrays", {}):
        expected.add("projection")
    if set(manifest.get("arrays", {})) != expected:
        raise ValueError("spectral manifest has unsupported array members")
    with zipfile.ZipFile(target) as archive:
        info = archive.infolist()
        if (
            len(info) != len(expected)
            or {item.filename for item in info} != {name + ".npy" for name in expected}
            or sum(item.file_size for item in info) > 128 * 1024**2
        ):
            raise ValueError(
                "spectral archive has duplicate, oversized or unknown members"
            )
    with np.load(target, allow_pickle=False) as archive:
        arrays = {name: _real_array(archive[name], name) for name in expected}
        for name, array in arrays.items():
            record = manifest["arrays"][name]
            if record.get("shape") != list(array.shape) or record.get(
                "sha256"
            ) != _array_hash(array):
                raise ValueError("spectral array shape or hash mismatch")
        data = SpectralDataset(
            manifest["beta"],
            arrays["tau"],
            arrays["values"],
            arrays["covariance"],
            covariance_rtol=manifest["covariance_rtol"],
            covariance_atol=manifest["covariance_atol"],
            provenance=manifest["provenance"],
            projection=arrays.get("projection"),
            source_fingerprint=manifest["source_fingerprint"],
        )
    if (
        data.fingerprint != manifest["fingerprint"]
        or data.modes.effective_rank != manifest["effective_rank"]
    ):
        raise ValueError("spectral fingerprint or retained covariance rank mismatch")
    return data


class TauDesign(Protocol):
    """Real time-domain design, with physical rather than reduced times."""

    @property
    def beta(self) -> float: ...

    def design_tau(self, tau: npt.ArrayLike) -> npt.ArrayLike: ...


@dataclass(frozen=True, slots=True, eq=False)
class SpectralProjection:
    beta: float
    source_tau: FloatArray
    target_tau: FloatArray
    matrix: FloatArray
    diagnostics: LeastSquaresDiagnostics
    metadata_json: str

    def __post_init__(self) -> None:
        for name in ("source_tau", "target_tau", "matrix"):
            object.__setattr__(self, name, _owned(getattr(self, name)))

    def apply(
        self, data: SpectralDataset, *, covariance_rtol: float | None = None
    ) -> SpectralDataset:
        if data.beta != self.beta or not np.array_equal(data.tau, self.source_tau):
            raise ValueError("projection requires the exact prepared observation mesh")
        with np.errstate(over="ignore", invalid="ignore"):
            values = self.matrix @ data.values
            covariance = self.matrix @ data.covariance @ self.matrix.T
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(covariance)):
            raise NumericalBreakdownError("projected data or covariance is nonfinite")
        return SpectralDataset(
            self.beta,
            self.target_tau,
            values,
            covariance,
            covariance_rtol=(
                data.modes.relative_cutoff
                if covariance_rtol is None
                else covariance_rtol
            ),
            covariance_atol=data.modes.absolute_cutoff,
            provenance={
                "operation": "complete_linear_projection",
                "source": data.provenance,
                "basis": json.loads(self.metadata_json),
                "svd": asdict(self.diagnostics),
                "independent_interpolated_observations": False,
            },
            projection=self.matrix,
            source_fingerprint=data.fingerprint,
        )


def prepare_spectral_projection(
    basis: TauDesign,
    source_tau: npt.ArrayLike,
    target_tau: npt.ArrayLike,
    *,
    rcond: float | None = None,
    metadata: Any = None,
) -> SpectralProjection:
    """Prepare P=A_target A_source^+ through a cached transposed SVD solve.

    No explicit inverse or observation-sized identity is constructed. This is
    ordinary least squares, not generalized least squares; covariance is then
    propagated through exactly the same complete linear map.
    """
    beta = _positive_real(basis.beta, "basis beta")
    source, target = _tau_array(source_tau, beta), _tau_array(target_tau, beta)
    for points in (source, target):
        if (
            not 3 <= points.size <= 2048
            or points[0] != 0
            or points[-1] != beta
            or np.any(np.diff(points) <= 0)
        ):
            raise ValueError("projection meshes must be ordered and include endpoints")
    observed = _real_array(basis.design_tau(source), "source design")
    evaluated = _real_array(basis.design_tau(target), "target design")
    if (
        observed.ndim != 2
        or evaluated.ndim != 2
        or observed.shape[0] != source.size
        or evaluated.shape[0] != target.size
        or observed.shape[1] != evaluated.shape[1]
    ):
        raise ValueError("basis designs must have matching trailing coefficient axes")
    plan = LeastSquaresPlan(observed.T, rcond=rcond)
    if plan.diagnostics.effective_rank == 0:
        raise NumericalBreakdownError("projection design has no retained singular mode")
    matrix = plan.fit(evaluated.T).T
    return SpectralProjection(
        beta,
        source,
        target,
        _real_array(matrix, "projection"),
        plan.diagnostics,
        _document({} if metadata is None else metadata),
    )


def diagonal_spectral_datasets(
    beta: float,
    tau: npt.ArrayLike,
    values: npt.ArrayLike,
    covariance: npt.ArrayLike,
    labels: tuple[str, ...],
    **kwargs: Any,
) -> tuple[SpectralDataset, ...]:
    """Extract ordered scalar marginals from an explicit time-major joint C.

    Values have shape (time, diagonal_channel); C uses flattened time-major
    ordering. Cross-channel covariances are checked but separate scalar MaxEnt
    runs make no joint posterior or cross-channel uncertainty claim.
    """
    observed = _real_array(values, "diagonal observations")
    if observed.ndim != 2 or observed.shape[1] != len(labels) or not labels:
        raise ValueError("diagonal values must match nonempty ordered labels")
    if len(set(labels)) != len(labels) or any(not label for label in labels):
        raise ValueError("diagonal labels must be unique nonempty strings")
    matrix, _ = covariance_modes(covariance)
    if matrix.shape != (observed.size, observed.size):
        raise ValueError("joint covariance must use complete time-major ordering")
    results = []
    for column, label in enumerate(labels):
        indices = np.arange(column, observed.size, len(labels))
        results.append(
            SpectralDataset(
                beta,
                tau,
                observed[:, column],
                matrix[np.ix_(indices, indices)],
                provenance={
                    "channel": label,
                    "ordering": labels,
                    "operation": "scalar_marginal_no_joint_inference",
                },
                **kwargs,
            )
        )
    return tuple(results)


@dataclass(frozen=True, slots=True)
class MaxEntConfig:
    """Explicit finite grids and training-only discrepancy alpha selection."""

    omega_min: float = -4.0
    omega_max: float = 4.0
    n_omega: int = 401
    alpha_min: float = 1e-4
    alpha_max: float = 1e4
    n_alpha: int = 25
    maxiter: int = 1000
    default_model: str = "flat"
    normalize_mass: float | None = 1.0

    def __post_init__(self) -> None:
        for name in ("omega_min", "omega_max"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not math.isfinite(value):
                raise ValueError("spectral window must be finite")
        if not self.omega_min < 0 < self.omega_max:
            raise ValueError("spectral window must straddle zero")
        for name, lower, upper in (
            ("n_omega", 5, 2048),
            ("n_alpha", 2, 128),
            ("maxiter", 1, 10000),
        ):
            value = getattr(self, name)
            if type(value) is not int or not lower <= value <= upper:
                raise ValueError(f"{name} must be an integer in [{lower},{upper}]")
        if not _positive_real(self.alpha_min, "alpha_min") < _positive_real(
            self.alpha_max, "alpha_max"
        ):
            raise ValueError("alpha range must be increasing")
        if self.default_model not in {"flat", "broad_gaussian"}:
            raise ValueError("default_model must be flat or broad_gaussian")
        if self.normalize_mass is not None:
            _positive_real(self.normalize_mass, "normalize_mass")


@dataclass(frozen=True, slots=True, eq=False)
class MaxEntSpectrum:
    omega: FloatArray
    density: FloatArray
    alphas: FloatArray
    chi2: FloatArray
    selected_index: int
    converged: tuple[bool, ...]
    iterations: tuple[int, ...]
    effective_rank: int
    raw_mass: float
    mass: float
    normalized_training_whitened_rms: float
    diagnostics_json: str

    def __post_init__(self) -> None:
        for name in ("omega", "density", "alphas", "chi2"):
            object.__setattr__(self, name, _owned(getattr(self, name)))

    @property
    def diagnostics(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self.diagnostics_json))

    def evaluate_tau(self, beta: float, tau: npt.ArrayLike) -> FloatArray:
        """Gauss quadrature of the returned piecewise-linear spectral density."""
        return spectral_forward(beta, tau, self.omega, self.density)


def spectral_forward(
    beta: float, tau: npt.ArrayLike, omega: npt.ArrayLike, density: npt.ArrayLike
) -> FloatArray:
    """Independent eight-point quadrature on each linear spectral interval."""
    temperature = _positive_real(beta, "beta")
    times = _tau_array(tau, temperature)
    energies, spectral = _real_array(omega, "omega"), _real_array(density, "density")
    if (
        energies.ndim != 1
        or not 2 <= energies.size <= 2048
        or spectral.shape != energies.shape
        or np.any(np.diff(energies) <= 0)
    ):
        raise ValueError("spectral coordinates must be strictly increasing and matched")
    nodes, weights = np.polynomial.legendre.leggauss(8)
    widths = np.diff(energies)
    points = (
        (energies[:-1, None] + energies[1:, None]) / 2 + widths[:, None] * nodes / 2
    ).ravel()
    fractions = (nodes + 1) / 2
    amplitudes = spectral[:-1, None] * (1 - fractions) + spectral[1:, None] * fractions
    quadrature = (amplitudes * widths[:, None] * weights / 2).ravel()
    output = np.empty(times.size)
    # exp(-tau*w)*expit(beta*w), split by sign to avoid inf*0.
    for start in range(0, times.size, 64):
        chunk = times[start : start + 64, None]
        exponent = np.where(
            points >= 0, -chunk * points, (temperature - chunk) * points
        )
        kernel = -np.exp(exponent) * expit(temperature * np.abs(points))
        output[start : start + 64] = kernel @ quadrature
    if not np.all(np.isfinite(output)):
        raise NumericalBreakdownError("spectral forward quadrature is nonfinite")
    return _owned(output)


def continue_tau_maxent(
    data: SpectralDataset, config: MaxEntConfig | None = None
) -> MaxEntSpectrum:
    """Run actual TRIQS/MaxEnt4.0.1; unavailable or exhausted providers fail.

    Alpha minimizes |log(chi2 / retained noise rank)| over converged provider
    trials. Optional mass rescaling is explicit and its forward defect is
    recomputed. This is a regularized estimate, not unique peak recovery.
    """
    settings = MaxEntConfig() if config is None else config
    try:
        provider = importlib.import_module("triqs_maxent")
        minimizers = importlib.import_module("triqs_maxent.minimizers")
        version_module = importlib.import_module("triqs_maxent.version")
    except ImportError as error:
        raise ImportError(
            "continuation requires the optional TRIQS/MaxEnt4.0.1 runtime"
        ) from error
    version = str(getattr(version_module, "version", "unavailable"))
    if version != "4.0.1":
        raise ValueError(f"unsupported TRIQS/MaxEnt version {version}; require4.0.1")

    history: list[tuple[bool, int]] = []

    base_minimizer: Any = minimizers.LevenbergMinimizer

    class RecordingMinimizer(base_minimizer):  # type: ignore[misc]
        def minimize(self, *args: Any, **kwargs: Any) -> Any:
            result = super().minimize(*args, **kwargs)
            history.append((bool(self.converged), int(self.n_iter_last)))
            return result

    solver = provider.TauMaxEnt(
        cost_function="bryan",
        minimizer=RecordingMinimizer(maxiter=settings.maxiter),
        analyzers=[],
        interactive=False,
        scale_alpha=None,
        cov_threshold=data.modes.threshold,
    )
    solver.omega = provider.LinearOmegaMesh(
        settings.omega_min, settings.omega_max, settings.n_omega
    )
    solver.alpha_mesh = provider.LogAlphaMesh(
        settings.alpha_min, settings.alpha_max, settings.n_alpha
    )
    omega = np.asarray(solver.omega, dtype=float)
    default = np.ones(omega.size)
    if settings.default_model == "broad_gaussian":
        sigma = (settings.omega_max - settings.omega_min) / 4
        default = np.exp(-0.5 * (omega / sigma) ** 2)
    default /= np.trapezoid(default, omega)
    solver.D = provider.DataDefaultModel(default, solver.omega)
    solver.set_G_tau_data(data.tau.copy(), data.values.copy())
    solver.set_cov(data.modes.retained_covariance.copy())
    if len(solver.err) != data.modes.effective_rank:
        raise NumericalBreakdownError(
            "provider covariance truncation changed declared rank"
        )
    try:
        result = solver.run()
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        raise NumericalBreakdownError(
            "actual TauMaxEnt failed; no substitute method"
        ) from error
    spectra = _real_array(result.A, "provider spectra")
    alphas, chi2 = (
        _real_array(result.alpha, "provider alphas"),
        _real_array(result.chi2, "provider chi2"),
    )
    if (
        spectra.shape != (settings.n_alpha, settings.n_omega)
        or alphas.shape != (settings.n_alpha,)
        or chi2.shape != alphas.shape
        or len(history) != settings.n_alpha
    ):
        raise NumericalBreakdownError("provider returned incomplete alpha trials")
    converged = np.array([entry[0] for entry in history])
    eligible = converged & (chi2 >= 0) & np.all(spectra >= 0, axis=1)
    if not np.any(eligible):
        raise NumericalBreakdownError("no positive converged TauMaxEnt alpha trial")
    scores = np.full(chi2.size, np.inf)
    scores[eligible] = np.abs(
        np.log(
            np.maximum(chi2[eligible], np.finfo(float).tiny) / data.modes.effective_rank
        )
    )
    selected = int(np.argmin(scores))
    density = spectra[selected].copy()
    raw_mass = float(np.trapezoid(density, omega))
    if not math.isfinite(raw_mass) or raw_mass <= 0:
        raise NumericalBreakdownError(
            "selected spectral mass is nonpositive or nonfinite"
        )
    if settings.normalize_mass is not None:
        density *= settings.normalize_mass / raw_mass
    mass = float(np.trapezoid(density, omega))
    rms = data.modes.whitened_rms(
        spectral_forward(data.beta, data.tau, omega, density) - data.values
    )
    return MaxEntSpectrum(
        omega,
        density,
        alphas,
        chi2,
        selected,
        tuple(entry[0] for entry in history),
        tuple(entry[1] for entry in history),
        data.modes.effective_rank,
        raw_mass,
        mass,
        rms,
        _document(
            {
                "provider": "TRIQS/MaxEnt",
                "version": version,
                "config": asdict(settings),
                "dataset": data.fingerprint,
                "alpha_selection": (
                    "training_discrepancy_closest_chi2_per_retained_rank_to_one"
                ),
                "covariance_threshold": data.modes.threshold,
                "effective_rank": data.modes.effective_rank,
                "mass_operation": "none"
                if settings.normalize_mass is None
                else "explicit_rescale",
                "certification_status": "not_certified",
                "uncertainty_scope": (
                    "retained_covariance_modes_no_unique_peak_or_coverage_claim"
                ),
            }
        ),
    )
