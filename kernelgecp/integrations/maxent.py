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
    NumericArray,
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


class _SpectralLeastSquaresPlan(LeastSquaresPlan):
    """Integration-only observation projector over the inherited cached SVD."""

    __slots__ = ()

    def observation_projection(self) -> NumericArray:
        """Return the retained observation projector U_r U_r^H.

        This avoids forming A A^+ through a badly scaled coefficient basis.
        The full observation-sized matrix is explicitly bounded to2048 rows;
        it describes this finite SVD, not held-out prediction or certification.
        """
        if self._design.shape[0] > 2048:
            raise ValueError("observation projection matrix is limited to2048 rows")
        rank = self._diagnostics.effective_rank
        left = self._u[:, :rank]
        result = left @ left.conj().T
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("observation projection is nonfinite")
        return _immutable_copy(result)


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
    projection_source_tau: FloatArray | None
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
        projection_source_tau: npt.ArrayLike | None = None,
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
        original_tau = (
            None
            if projection_source_tau is None
            else _tau_array(projection_source_tau, temperature)
        )
        if (mapping is None) != (original_tau is None):
            raise ValueError("projection requires its original source time coordinates")
        if original_tau is not None and (
            not 3 <= original_tau.size <= 2048
            or original_tau[0] != 0
            or original_tau[-1] != temperature
            or np.any(np.diff(original_tau) <= 0)
            or mapping is None
            or mapping.shape[1] != original_tau.size
        ):
            raise ValueError("projection columns must match the ordered source mesh")
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
        object.__setattr__(
            self,
            "projection_source_tau",
            None if original_tau is None else _owned(original_tau),
        )
        object.__setattr__(self, "source_fingerprint", source_fingerprint)

    @property
    def provenance(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self.provenance_json))

    @property
    def fingerprint(self) -> str:
        digest = hashlib.sha256()
        arrays = [
            np.array(
                [self.beta, self.modes.relative_cutoff, self.modes.absolute_cutoff]
            ),
            self.tau,
            self.values,
            self.covariance,
        ]
        if self.projection is not None and self.projection_source_tau is not None:
            arrays.extend((self.projection, self.projection_source_tau))
        for array in arrays:
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
        assert data.projection_source_tau is not None
        arrays["projection_source_tau"] = data.projection_source_tau
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
        expected.update(("projection", "projection_source_tau"))
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
            projection_source_tau=arrays.get("projection_source_tau"),
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
    source_design: FloatArray | None = None
    decomposed_design: FloatArray | None = None
    transformed_diagnostics: LeastSquaresDiagnostics | None = None
    projection_provenance_json: str = "{}"

    def __post_init__(self) -> None:
        for name in (
            "source_tau",
            "target_tau",
            "matrix",
            "source_design",
            "decomposed_design",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _owned(value))

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
                "projection_arithmetic": json.loads(self.projection_provenance_json),
                "independent_interpolated_observations": False,
            },
            projection=(
                self.matrix
                if data.projection is None
                else self.matrix @ data.projection
            ),
            projection_source_tau=(
                data.tau
                if data.projection_source_tau is None
                else data.projection_source_tau
            ),
            source_fingerprint=data.source_fingerprint or data.fingerprint,
        )


def prepare_spectral_projection(
    basis: TauDesign,
    source_tau: npt.ArrayLike,
    target_tau: npt.ArrayLike,
    *,
    rcond: float | None = None,
    metadata: Any = None,
    precision_bits: int = 53,
) -> SpectralProjection:
    """Prepare the complete linear observation map through a cached SVD.

    On the same mesh, form U_r U_r^T directly for a stable projector. Otherwise
    use the transposed solve P=A_target A_source^+. No explicit inverse or
    observation-sized identity is constructed. This is
    ordinary least squares, not generalized least squares; covariance is then
    propagated through exactly the same complete linear map.

    ``precision_bits=53`` retains the binary64 path. An explicit integer from
    128 through 1024 enables mpmath QR followed by a cached binary64 SVD of the
    rounded orthogonal coordinates. That option supports only same-mesh designs
    with full column rank, at most 256 rows and 128 columns. Truncated,
    underdetermined and insufficiently separated rank decisions are rejected;
    the separation guard is a heuristic, not a certificate. Original design
    singular values, cutoff and conditioning remain in ``diagnostics``.
    Transformed diagnostics, immutable exact decomposition inputs and arithmetic
    provenance are exposed separately. No automatic precision fallback occurs.
    """
    if type(precision_bits) is not int or not (
        precision_bits == 53 or 128 <= precision_bits <= 1024
    ):
        raise ValueError("precision_bits must be 53 or an integer from 128 to 1024")
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
    same_mesh = np.array_equal(source, target)
    if precision_bits != 53 and not same_mesh:
        raise ValueError("higher precision projection requires the same mesh")
    observed = _owned(_real_array(basis.design_tau(source), "source design"))
    evaluated = (
        observed
        if same_mesh
        else _real_array(basis.design_tau(target), "target design")
    )
    if (
        observed.ndim != 2
        or evaluated.ndim != 2
        or observed.shape[0] != source.size
        or evaluated.shape[0] != target.size
        or observed.shape[1] != evaluated.shape[1]
    ):
        raise ValueError("basis designs must have matching trailing coefficient axes")
    if precision_bits != 53 and (observed.shape[0] > 256 or observed.shape[1] > 128):
        raise ValueError("higher precision projection exceeds resource limits")
    plan = _SpectralLeastSquaresPlan(observed if same_mesh else observed.T, rcond=rcond)
    if plan.diagnostics.effective_rank == 0:
        raise NumericalBreakdownError("projection design has no retained singular mode")
    decomposed = _real_array(plan.design, "decomposed design")
    transformed_diagnostics = None
    method = "binary64_cached_svd"
    if precision_bits != 53:
        rows, columns = observed.shape
        diagnostics = plan.diagnostics
        if rows < columns or diagnostics.effective_rank != columns:
            raise ValueError("higher precision projection requires full column rank")
        # Both margins are heuristics, not certified binary64 separation bounds.
        threshold = diagnostics.rank_threshold
        rounding_floor = (
            8
            * np.finfo(np.float64).eps
            * max(rows, columns)
            * diagnostics.singular_values[0]
        )
        if threshold is None or diagnostics.singular_values[-1] <= max(
            2 * threshold, rounding_floor
        ):
            raise ValueError("higher precision projection has ambiguous rank boundary")
        try:
            import mpmath

            context = mpmath.mp.clone()
            context.prec = precision_bits
            design_mp = context.matrix(
                [[context.mpf(float(value)) for value in row] for row in observed]
            )
            q, _ = context.qr(design_mp)
            decomposed = _owned(
                np.array(
                    [[float(q[i, j]) for j in range(columns)] for i in range(rows)]
                )
            )
        except (ArithmeticError, ValueError, OverflowError, MemoryError) as error:
            raise NumericalBreakdownError("higher precision QR failed") from error
        transformed = _SpectralLeastSquaresPlan(decomposed, rcond=rcond)
        if transformed.diagnostics.effective_rank != columns:
            raise NumericalBreakdownError("rounded orthogonal design lost rank")
        transformed_diagnostics = transformed.diagnostics
        matrix = transformed.observation_projection()
        method = "mpmath_qr_then_binary64_cached_svd"
    else:
        matrix = plan.observation_projection() if same_mesh else plan.fit(evaluated.T).T
    provenance = {
        "method": method,
        "precision_bits": precision_bits,
        "source_design_sha256": _array_hash(observed),
        "decomposed_design_sha256": _array_hash(decomposed),
        "source_design_shape": list(observed.shape),
        "decomposed_design_shape": list(decomposed.shape),
        "hash_scope": "C-order little-endian binary64 array bytes with recorded shapes",
        "rank_guard": (
            "smallest > max(2 * threshold, 8 * eps64 * max(shape) * largest); heuristic"
        )
        if precision_bits != 53
        else None,
        "transformed_svd": asdict(transformed_diagnostics)
        if transformed_diagnostics
        else None,
    }
    return SpectralProjection(
        beta,
        source,
        target,
        _real_array(matrix, "projection"),
        plan.diagnostics,
        _document({} if metadata is None else metadata),
        observed,
        decomposed,
        transformed_diagnostics,
        _document(provenance),
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
    caller_provenance = kwargs.pop("provenance", {})
    if not isinstance(caller_provenance, dict):
        raise ValueError("provenance must be a JSON object")
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
                    "source": caller_provenance,
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
    """Explicit finite grids and training-only discrepancy alpha selection.

    minimizer_mu0 is the native minimizer's initial damping, not a floor;
    the native algorithm can reduce it during a solve. The default preserves
    the provider setting; scientific campaigns declare their candidate explicitly.
    The rescale default applies mass after fitting. Constrained mode instead
    uses a normalized native parameterization and full chain-rule derivatives,
    with explicit numerical rank and resource limits. Convergence is not assured.
    """

    omega_min: float = -4.0
    omega_max: float = 4.0
    n_omega: int = 401
    alpha_min: float = 1e-4
    alpha_max: float = 1e4
    n_alpha: int = 25
    maxiter: int = 1000
    default_model: str = "flat"
    normalize_mass: float | None = 1.0
    minimizer_mu0: float = 1e-18
    mass_mode: str = "rescale"

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
        _positive_real(self.minimizer_mu0, "minimizer_mu0")
        if self.default_model not in {"flat", "broad_gaussian"}:
            raise ValueError("default_model must be flat or broad_gaussian")
        if self.normalize_mass is not None:
            _positive_real(self.normalize_mass, "normalize_mass")
        if not isinstance(self.mass_mode, str) or self.mass_mode not in {
            "rescale",
            "constrained",
        }:
            raise ValueError("mass_mode must be rescale or constrained")
        if self.mass_mode == "constrained" and self.normalize_mass is None:
            raise ValueError("constrained mass requires an explicit normalize_mass")

    def as_record(self) -> dict[str, Any]:
        """Preserve legacy rescale receipts; record an explicit nondefault mode."""
        result = asdict(self)
        if self.mass_mode == "rescale":
            result.pop("mass_mode")
        return result


@dataclass(frozen=True, slots=True, eq=False)
class MaxEntSpectrum:
    omega: FloatArray
    density: FloatArray
    alphas: FloatArray
    chi2: FloatArray
    trial_spectra: FloatArray
    selected_index: int
    converged: tuple[bool, ...]
    iterations: tuple[int, ...]
    effective_rank: int
    raw_mass: float
    mass: float
    normalized_training_whitened_rms: float
    diagnostics_json: str

    def __post_init__(self) -> None:
        for name in ("omega", "density", "alphas", "chi2", "trial_spectra"):
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


def _selection_metrics(
    data: SpectralDataset,
    omega: FloatArray,
    spectra: FloatArray,
    chi2: FloatArray,
    converged: npt.NDArray[np.bool_],
    normalize_mass: float | None,
    *,
    mass_mode: str = "rescale",
) -> tuple[int, list[float | None], list[float | None]]:
    """Score training fits after the delivered mass operation."""
    masses = np.asarray(np.trapezoid(spectra, omega, axis=1), dtype=np.float64)
    eligible = (
        converged
        & (chi2 >= 0)
        & np.all(spectra >= 0, axis=1)
        & np.isfinite(masses)
        & (masses > 0)
    )
    if mass_mode == "constrained":
        if normalize_mass is None:
            raise ValueError("constrained mass requires a target")
        eligible &= np.isclose(masses, normalize_mass, rtol=1e-10, atol=0)
    mass_records = [float(m) if np.isfinite(m) else None for m in masses]
    metrics: list[float | None] = [None] * len(masses)
    for i in np.flatnonzero(eligible):
        if normalize_mass is None:
            value = float(chi2[i])
        else:
            density = (
                spectra[i]
                if mass_mode == "constrained"
                else spectra[i] * (normalize_mass / masses[i])
            )
            prediction = spectral_forward(
                data.beta,
                data.tau
                if data.projection_source_tau is None
                else data.projection_source_tau,
                omega,
                density,
            )
            if data.projection is not None:
                prediction = data.projection @ prediction
            whitened = data.modes.whiten(prediction - data.values)
            value = float(np.dot(whitened, whitened))
        if math.isfinite(value):
            metrics[i] = value
    candidates = [i for i, value in enumerate(metrics) if value is not None]
    if not candidates:
        raise NumericalBreakdownError(
            "no positive finite-mass converged TauMaxEnt alpha trial"
        )
    selected = min(
        candidates,
        key=lambda i: abs(
            math.log(
                max(cast(float, metrics[i]), np.finfo(float).tiny)
                / data.modes.effective_rank
            )
        ),
    )
    return selected, mass_records, metrics


def continue_tau_maxent(
    data: SpectralDataset, config: MaxEntConfig | None = None
) -> MaxEntSpectrum:
    """Run actual TRIQS/MaxEnt4.0.1; unavailable or exhausted providers fail.

    Alpha minimizes |log(chi2 / retained noise rank)| over eligible converged
    trials. With mass rescaling, chi2 is recomputed from the normalized physical
    forward model and training data; without it, provider chi2 is used. The
    default generates the family without a mass constraint in optimization.
    Explicit constrained mode uses a normalized positive native parameterization
    with full cost derivatives and no posthoc mass scaling. Its numerical rank
    and resource limits are recorded; no alternative optimizer is substituted.
    This is a regularized estimate, not unique peak recovery.
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

    cost_function: Any = "bryan"
    mass_parameterization: Any = None
    if settings.mass_mode == "constrained":
        from ._maxent_native_mass import make_native_mass_cost

        assert settings.normalize_mass is not None
        cost_function, mass_parameterization = make_native_mass_cost(
            settings.normalize_mass
        )
    solver = provider.TauMaxEnt(
        cost_function=cost_function,
        minimizer=RecordingMinimizer(
            maxiter=settings.maxiter, mu0=settings.minimizer_mu0
        ),
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
    if data.projection is not None:
        # Apply the same complete observation map to the native physical
        # forward kernel, before the provider rotates covariance eigenmodes.
        # Interpolation changes observations, covariance AND the model.
        kernels = importlib.import_module("triqs_maxent.kernels")
        assert data.projection_source_tau is not None
        original_kernel = kernels.TauKernel(
            data.projection_source_tau.copy(), solver.omega, beta=data.beta
        )
        solver.K = kernels.DataKernel(
            data.tau.copy(), solver.omega, data.projection @ original_kernel.K
        )
    solver.set_cov(data.modes.retained_covariance.copy())
    if len(solver.err) != data.modes.effective_rank:
        raise NumericalBreakdownError(
            "provider covariance truncation changed declared rank"
        )
    if settings.mass_mode == "constrained":
        # The provider multiplies A_init by mesh weights before H.inv().
        # Supply a density, so its initial hidden image is the weighted default.
        assert settings.normalize_mass is not None
        weighted_default = np.asarray(solver.D.D, dtype=float)
        solver.A_init = (
            settings.normalize_mass
            * weighted_default
            / np.sum(weighted_default)
            / np.asarray(solver.omega.delta, dtype=float)
        )
    try:
        # The native minimizer changes NumPy's error policy internally. Keep
        # the caller's policy intact even if a native solve raises midway.
        with np.errstate(**np.geterr()):
            result = solver.run()
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        raise NumericalBreakdownError(
            f"actual TauMaxEnt failed; no substitute method; "
            f"{type(error).__name__}: {error}; completed trials: {len(history)}; "
            f"history (converged, iterations): {history}"
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
    selected, trial_masses, selection_chi2 = _selection_metrics(
        data,
        omega,
        spectra,
        chi2,
        converged,
        settings.normalize_mass,
        mass_mode=settings.mass_mode,
    )
    density = spectra[selected].copy()
    raw_mass = float(np.trapezoid(density, omega))
    if not math.isfinite(raw_mass) or raw_mass <= 0:
        raise NumericalBreakdownError(
            "selected spectral mass is nonpositive or nonfinite"
        )
    if settings.normalize_mass is not None and settings.mass_mode == "rescale":
        density *= settings.normalize_mass / raw_mass
    mass = float(np.trapezoid(density, omega))
    prediction = spectral_forward(
        data.beta,
        data.tau if data.projection_source_tau is None else data.projection_source_tau,
        omega,
        density,
    )
    if data.projection is not None:
        prediction = data.projection @ prediction
    rms = data.modes.whitened_rms(prediction - data.values)
    return MaxEntSpectrum(
        omega,
        density,
        alphas,
        chi2,
        spectra,
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
                "config": settings.as_record(),
                "dataset": data.fingerprint,
                "alpha_selection": (
                    "provider_chi2_training_discrepancy"
                    if settings.normalize_mass is None
                    else "normalized_physical_forward_training_discrepancy"
                ),
                "selection_rule": "closest_abs_log_chi2_per_retained_rank_to_one",
                "trial_masses": trial_masses,
                "selection_chi2": selection_chi2,
                "covariance_threshold": data.modes.threshold,
                "effective_rank": data.modes.effective_rank,
                "forward_operator": "native_tau_kernel"
                if data.projection is None
                else "complete_projection_times_native_source_tau_kernel",
                "mass_operation": "none"
                if settings.normalize_mass is None
                else (
                    "native_mass_constraint"
                    if settings.mass_mode == "constrained"
                    else "explicit_rescale"
                ),
                "mass_parameterization": (
                    mass_parameterization.diagnostics
                    if mass_parameterization is not None
                    else None
                ),
                "certification_status": "not_certified",
                "uncertainty_scope": (
                    "retained_covariance_modes_no_unique_peak_or_coverage_claim"
                ),
            }
        ),
    )
