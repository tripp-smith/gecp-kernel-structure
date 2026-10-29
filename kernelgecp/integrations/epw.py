"""EPW 6 real anisotropic Matsubara/Nambu interchange, in eV units.

This is a downstream representation of declared native observations. It does
not implement EPW's specialized IR objects or replace its internal solver.
Native state ordinals, including degenerate energies, are never deduplicated.
"""

from __future__ import annotations

import hashlib
import json
import math
import operator
import re
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from typing import Any, cast

import numpy as np
import numpy.typing as npt

from ..basis import FermionicBasis, FermionicFit
from ..fitting import _contains_boolean_sequence, _immutable_copy, _numeric_array
from ..types import FloatArray, NumericalBreakdownError

__all__ = [
    "EPWAnisotropicData",
    "EPWNambuFit",
    "EPWRecoveryDiagnostics",
    "EPWRecoveryResult",
    "fit_epw_nambu",
    "nambu_green",
    "read_epw_anisotropic",
    "recover_epw_nambu",
    "write_epw_anisotropic",
]

type ComplexArray = npt.NDArray[np.complex128]
type IntArray = npt.NDArray[np.int64]
KB_EV_PER_K = 1.380649e-23 / 1.602176634e-19
_INPUT_SCHEMA = "kernelgecp.epw-anisotropic-input.v1"
_ARCHIVE_SCHEMA = "kernelgecp.epw-nambu-fit.v1"
_ORDER = "frequency_major_ik_ibnd_ordinal"
_FIELD = re.compile(r" *[+-]?\d\.\d{10}[EeDd][+-]\d{2}")
_HEADER_FIELDS = (
    "#        w [eV]",
    "Enk-Ef [eV]",
    "znorm(w)",
    "delta(w) [eV]",
    "shift(w) [eV]",
)
_REQUIRED = {
    "schema",
    "provider",
    "temperature_kelvin",
    "units",
    "state_order",
    "state_count",
    "matsubara_indices",
    "positive_only",
    "full_bandwidth",
    "chi_convention",
    "reference_fermi_energy_ev",
    "chemical_potential_offset_ev",
}
_OPTIONAL = {
    "native_sha256",
    "native_filename",
    "source_provenance",
    "source_native_sha256",
}


def _real(value: object, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} is unrepresentable in float64") from error
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite{' and positive' if positive else ''}")
    return result


def _integer(value: object, name: str, *, positive: bool = False) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        result = operator.index(value)  # type: ignore[arg-type]
    except TypeError as error:
        raise ValueError(f"{name} must be an integer") from error
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _array(values: npt.ArrayLike, name: str, ndim: int) -> FloatArray:
    result = _numeric_array(values, name)
    if np.iscomplexobj(result) or result.ndim != ndim:
        raise ValueError(f"{name} must be a real {ndim}-dimensional array")
    return cast(FloatArray, result)


def _own(values: FloatArray) -> FloatArray:
    return cast(FloatArray, _immutable_copy(values))


def _indices(values: npt.ArrayLike, name: str = "matsubara_indices") -> IntArray:
    if _contains_boolean_sequence(values):
        raise ValueError(f"{name} must not contain booleans")
    result = np.asarray(values)
    if result.ndim != 1 or not result.size or result.dtype.kind not in "iu":
        raise ValueError(f"{name} must be a nonempty integer vector")
    if np.any(result > 2**50) or (
        result.dtype.kind == "i" and np.any(result < -(2**50))
    ):
        raise ValueError(f"{name} magnitude exceeds float64 frequency resolution")
    converted = np.asarray(result, dtype=np.int64)
    if np.unique(converted).size != converted.size:
        raise ValueError(f"{name} must be unique and retain explicit order")
    return np.frombuffer(converted.tobytes(), dtype=np.int64)


def _digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json_object(value: object, name: str) -> str:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    try:
        result = json.dumps(
            dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False
        )
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must contain finite JSON values") from error
    return result


def _load_json(path: Path) -> dict[str, Any]:
    if path.stat().st_size > 1024**2:
        raise ValueError("EPW application manifest exceeds1 MiB")

    def pairs(entries: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in entries:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value: str) -> None:
        raise ValueError(f"nonfinite JSON constant: {value}")

    result = json.loads(
        path.read_text(), object_pairs_hook=pairs, parse_constant=constant
    )
    if not isinstance(result, dict):
        raise ValueError("manifest must be a JSON object")
    return cast(dict[str, Any], result)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class _Layout:
    frequencies_ev: FloatArray
    state_energies_ev: FloatArray
    indices: IntArray
    temperature_kelvin: float
    beta: float
    delta_mu_ev: float
    full_bandwidth: bool
    context_json: str

    def __init__(
        self,
        frequencies: npt.ArrayLike,
        energies: npt.ArrayLike,
        context: Mapping[str, Any],
    ) -> None:
        encoded = _json_object(context, "context")
        metadata = json.loads(encoded)
        if not metadata.keys() >= _REQUIRED or metadata.keys() - _REQUIRED - _OPTIONAL:
            raise ValueError("missing or unsupported EPW context fields")
        if metadata["schema"] != _INPUT_SCHEMA or metadata["provider"] != {
            "qe": "7.5",
            "epw": "6.0",
        }:
            raise ValueError("this native format supports declared qe-7.5/EPW 6.0")
        if metadata["units"] != "eV" or metadata["state_order"] != _ORDER:
            raise ValueError(
                "explicit eV units and native state ordinal order required"
            )
        for key in ("full_bandwidth", "positive_only"):
            if type(metadata[key]) is not bool:
                raise ValueError(f"{key} must be an explicit boolean")
        convention = (
            "native_fbw_column"
            if metadata["full_bandwidth"]
            else "absent_chi_is_zero_non_fbw"
        )
        if metadata["chi_convention"] != convention:
            raise ValueError("explicit matching native chi convention required")
        temperature = _real(
            metadata["temperature_kelvin"], "temperature_kelvin", positive=True
        )
        _real(metadata["reference_fermi_energy_ev"], "reference_fermi_energy_ev")
        offset = _real(
            metadata["chemical_potential_offset_ev"], "chemical_potential_offset_ev"
        )
        thermal_energy = KB_EV_PER_K * temperature
        if not math.isfinite(thermal_energy) or thermal_energy <= 0:
            raise ValueError("temperature is unrepresentable in physical eV units")
        beta = 1 / thermal_energy
        if not math.isfinite(beta) or beta <= 0:
            raise ValueError(
                "temperature is unrepresentable in physical inverse-eV units"
            )
        modes = _indices(metadata["matsubara_indices"])
        freq, energy = (
            _array(frequencies, "frequencies_ev", 1),
            _array(energies, "state_energies_ev", 1),
        )
        count = _integer(metadata["state_count"], "state_count", positive=True)
        if freq.shape != modes.shape or energy.shape != (count,):
            raise ValueError(
                "context frequency/state shapes do not match native arrays"
            )
        if metadata["positive_only"]:
            if not np.array_equal(modes, np.arange(len(modes))):
                raise ValueError(
                    "positive native mesh must contain ordered n=0,...,N-1"
                )
        elif set(modes.tolist()) != set((-modes - 1).tolist()):
            raise ValueError(
                "signed native mesh must contain both declared fermionic signs"
            )
        with np.errstate(all="ignore"):
            expected = (2 * modes + 1) * np.pi / beta
            rounding = (6e-11 + 32 * np.finfo(float).eps) * np.abs(expected)
        if (
            not np.all(np.isfinite(expected))
            or np.any(expected == 0)
            or np.any(np.abs(freq - expected) > rounding)
        ):
            raise ValueError(
                "native frequencies disagree with temperature/ordinary indices"
            )
        for key in ("native_sha256", "source_native_sha256"):
            if key in metadata and (
                not isinstance(metadata[key], str)
                or re.fullmatch(r"[0-9a-f]{64}", metadata[key]) is None
            ):
                raise ValueError(f"{key} must be a lowercase SHA256")
        object.__setattr__(self, "frequencies_ev", _own(freq))
        object.__setattr__(self, "state_energies_ev", _own(energy))
        object.__setattr__(self, "indices", modes)
        object.__setattr__(self, "temperature_kelvin", temperature)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "delta_mu_ev", offset)
        object.__setattr__(self, "full_bandwidth", metadata["full_bandwidth"])
        object.__setattr__(self, "context_json", encoded)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class EPWAnisotropicData:
    """Owned real EPW fields, shape (frequency, native state ordinal).

    ``context`` is the source record, not a claim that a newly reconstructed
    field still has the source file's bytes. Only the native reader marks
    source bytes verified. The native writer returns the new actual manifest.
    """

    _layout: _Layout
    z: FloatArray
    delta_ev: FloatArray
    chi_ev: FloatArray
    _source_bytes: bytes | None

    def __init__(
        self,
        frequencies_ev: npt.ArrayLike,
        state_energies_ev: npt.ArrayLike,
        z: npt.ArrayLike,
        delta_ev: npt.ArrayLike,
        chi_ev: npt.ArrayLike | None,
        *,
        context: Mapping[str, Any],
    ) -> None:
        layout = _Layout(frequencies_ev, state_energies_ev, context)
        zz, gap = _array(z, "z", 2), _array(delta_ev, "delta_ev", 2)
        shape = (len(layout.indices), len(layout.state_energies_ev))
        if zz.shape != shape or gap.shape != shape or np.any(zz <= 0):
            raise ValueError(
                "positive Z and matching frequency-by-state fields required"
            )
        if chi_ev is None:
            if layout.full_bandwidth:
                raise ValueError("full-bandwidth native shift column is required")
            shift = np.zeros(shape)
        else:
            shift = _array(chi_ev, "chi_ev", 2)
        if shift.shape != shape or (not layout.full_bandwidth and np.any(shift != 0)):
            raise ValueError(
                "shift shape/convention disagrees with native bandwidth mode"
            )
        object.__setattr__(self, "_layout", layout)
        object.__setattr__(self, "z", _own(zz))
        object.__setattr__(self, "delta_ev", _own(gap))
        object.__setattr__(self, "chi_ev", _own(shift))
        object.__setattr__(self, "_source_bytes", None)

    @property
    def frequencies_ev(self) -> FloatArray:
        return self._layout.frequencies_ev

    @property
    def state_energies_ev(self) -> FloatArray:
        return self._layout.state_energies_ev

    @property
    def indices(self) -> IntArray:
        return self._layout.indices

    @property
    def beta(self) -> float:
        return self._layout.beta

    @property
    def temperature_kelvin(self) -> float:
        return self._layout.temperature_kelvin

    @property
    def context(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._layout.context_json))

    @property
    def source_bytes_verified(self) -> bool:
        return self._source_bytes is not None

    def replace_fields(
        self, z: npt.ArrayLike, delta_ev: npt.ArrayLike, chi_ev: npt.ArrayLike
    ) -> EPWAnisotropicData:
        """Own new fields retaining the source context without its source bytes."""
        return EPWAnisotropicData(
            self.frequencies_ev,
            self.state_energies_ev,
            z,
            delta_ev,
            chi_ev,
            context=self.context,
        )


def read_epw_anisotropic(
    path: str | Path, manifest: Mapping[str, Any] | str | Path
) -> EPWAnisotropicData:
    """Read strict fixed-width native rows and verify their mandatory source hash."""
    context = (
        _load_json(Path(manifest))
        if isinstance(manifest, (str, Path))
        else dict(manifest)
    )
    raw = Path(path).read_bytes()
    if context.get("native_sha256") != _digest(raw):
        raise ValueError("native file SHA256 is missing or mismatched")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("EPW native data must be ASCII") from error
    full = context.get("full_bandwidth")
    if type(full) is not bool:
        raise ValueError("full_bandwidth must be an explicit boolean")
    fields = 5 if full else 4
    if not lines or lines[0] != "".join(
        f"{label:>20}" for label in _HEADER_FIELDS[:fields]
    ):
        raise ValueError("native anisotropic header does not match the declared format")
    count = _integer(context.get("state_count"), "state_count", positive=True)
    modes = _indices(context.get("matsubara_indices", []))
    if len(lines) - 1 != len(modes) * count:
        raise ValueError("native row count is truncated or disagrees with context")
    parsed = np.empty((len(lines) - 1, fields), dtype=np.float64)
    for row, line in enumerate(lines[1:]):
        if len(line) != 20 * fields:
            raise ValueError("native row must contain exactly fixed ES20.10 fields")
        for column in range(fields):
            text = line[column * 20 : (column + 1) * 20]
            if _FIELD.fullmatch(text) is None:
                raise ValueError("malformed or unsupported native ES20.10 value")
            parsed[row, column] = float(text.replace("D", "E").replace("d", "e"))
    if not np.all(np.isfinite(parsed)):
        raise ValueError("native observations must be finite")
    array = parsed.reshape(len(modes), count, fields)
    frequencies, energies = array[:, 0, 0], array[0, :, 1]
    if not np.all(array[:, :, 0] == frequencies[:, None]) or not np.all(
        array[:, :, 1] == energies[None, :]
    ):
        raise ValueError("native frequency/state ordinal order changes between rows")
    result = EPWAnisotropicData(
        frequencies,
        energies,
        array[:, :, 2],
        array[:, :, 3],
        array[:, :, 4] if full else None,
        context=context,
    )
    object.__setattr__(result, "_source_bytes", bytes(raw))
    return result


def _native_bytes(data: EPWAnisotropicData, preserve_source: bool) -> bytes:
    if preserve_source and data._source_bytes is not None:
        return data._source_bytes
    count = 5 if data._layout.full_bandwidth else 4
    lines = ["".join(f"{label:>20}" for label in _HEADER_FIELDS[:count])]
    for i, frequency in enumerate(data.frequencies_ev):
        for j, energy in enumerate(data.state_energies_ev):
            row = [frequency, energy, data.z[i, j], data.delta_ev[i, j]]
            if count == 5:
                row.append(data.chi_ev[i, j])
            fields = [f"{value:20.10E}" for value in row]
            if any(
                len(field) != 20 or _FIELD.fullmatch(field) is None for field in fields
            ):
                raise ValueError(
                    "value is not representable by supported native ES20.10"
                )
            lines.append("".join(fields))
    return ("\n".join(lines) + "\n").encode("ascii")


def write_epw_anisotropic(
    data: EPWAnisotropicData,
    path: str | Path,
    *,
    preserve_source: bool = True,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write native bytes and return a context with their actual new SHA256."""
    if not isinstance(data, EPWAnisotropicData):
        raise TypeError("data must be EPWAnisotropicData")
    if type(preserve_source) is not bool or type(overwrite) is not bool:
        raise ValueError("writer controls must be booleans")
    target = Path(path)
    if target.exists() and not overwrite:
        raise FileExistsError(target)
    raw = _native_bytes(data, preserve_source)
    result = data.context
    if "native_sha256" in result:
        result["source_native_sha256"] = result["native_sha256"]
    result["native_sha256"], result["native_filename"] = _digest(raw), target.name
    target.write_bytes(raw)
    return result


def nambu_green(
    data: EPWAnisotropicData, *, rows: npt.ArrayLike | None = None
) -> ComplexArray:
    """Solve the declared real-Eliashberg system on explicit native rows.

    Supplying rows isolates observations from reserved validation values;
    omitted rows are never used to construct a Green matrix.
    """
    if not isinstance(data, EPWAnisotropicData):
        raise TypeError("data must be EPWAnisotropicData")
    selected = np.arange(len(data.indices)) if rows is None else _indices(rows, "rows")
    if np.any(selected < 0) or np.any(selected >= len(data.indices)):
        raise ValueError("rows must select native frequency rows")
    z, gap, shift = data.z[selected], data.delta_ev[selected], data.chi_ev[selected]
    matrix = np.zeros((*z.shape, 2, 2), dtype=np.complex128)
    with np.errstate(all="ignore"):
        iwz = 1j * data.frequencies_ev[selected, None] * z
        xi = data.state_energies_ev[None, :] - data._layout.delta_mu_ev + shift
        matrix[:, :, 0, 0], matrix[:, :, 1, 1] = iwz - xi, iwz + xi
        matrix[:, :, 0, 1] = matrix[:, :, 1, 0] = -z * gap
    if not np.all(np.isfinite(matrix)):
        raise NumericalBreakdownError("nonfinite Nambu system")
    try:
        result = np.linalg.solve(matrix, np.broadcast_to(np.eye(2), matrix.shape))
    except np.linalg.LinAlgError as error:
        raise NumericalBreakdownError("singular Nambu system") from error
    if not np.all(np.isfinite(result)):
        raise NumericalBreakdownError("nonfinite Nambu Green function")
    return cast(ComplexArray, _immutable_copy(np.asarray(result, dtype=np.complex128)))


@dataclass(frozen=True, slots=True)
class EPWRecoveryDiagnostics:
    condition_max: float
    inverse_derivative_gain_max: float
    z_derivative_gain_max: float
    solve_residual_max: float
    offdiagonal_relative_defect: float
    imaginary_z_max: float
    imaginary_delta_ev_max: float
    imaginary_chi_ev_max: float
    non_fbw_shift_zero_defect_ev: float = 0.0
    precision_bits: int = 53
    certification_status: str = "finite_float64_diagnostics_not_a_certificate"


@dataclass(frozen=True, slots=True)
class EPWRecoveryResult:
    data: EPWAnisotropicData
    diagnostics: EPWRecoveryDiagnostics


def recover_epw_nambu(
    data: EPWAnisotropicData,
    values: npt.ArrayLike,
    *,
    max_condition: float = 1e10,
    max_z_amplification: float = 1e10,
    real_tolerance: float = 1e-8,
) -> EPWRecoveryResult:
    """Recover real native fields with solves and explicit sensitivity rejection.

    The derivative gains are local binary64 diagnostics, not certified finite-
    perturbation bounds. Accepted imaginary/offdiagonal defects are reported.
    """
    if not isinstance(data, EPWAnisotropicData):
        raise TypeError("data must be EPWAnisotropicData")
    return _recover_layout(
        data._layout,
        values,
        max_condition=max_condition,
        max_z_amplification=max_z_amplification,
        real_tolerance=real_tolerance,
    )


def _recover_layout(
    layout: _Layout,
    values: npt.ArrayLike,
    *,
    max_condition: float = 1e10,
    max_z_amplification: float = 1e10,
    real_tolerance: float = 1e-8,
) -> EPWRecoveryResult:
    condition_limit = _real(max_condition, "max_condition", positive=True)
    gain_limit = _real(max_z_amplification, "max_z_amplification", positive=True)
    real_limit = _real(real_tolerance, "real_tolerance", positive=True)
    if condition_limit < 1 or real_limit >= 1:
        raise ValueError("invalid recovery condition/real tolerances")
    array = _numeric_array(values, "Nambu values")
    if array.shape != (len(layout.indices), len(layout.state_energies_ev), 2, 2):
        raise ValueError(
            "Nambu values must preserve native frequency/state/matrix axes"
        )
    singular = np.linalg.svd(array, compute_uv=False)
    with np.errstate(all="ignore"):
        condition = singular[..., 0] / singular[..., -1]
        gain = 1 / singular[..., -1] ** 2
        z_gain = gain / np.abs(layout.frequencies_ev[:, None])
    if (
        not np.all(np.isfinite(condition))
        or not np.all(np.isfinite(z_gain))
        or np.any(condition > condition_limit)
        or np.any(z_gain > gain_limit)
    ):
        raise NumericalBreakdownError(
            "unsafe Nambu recovery conditioning/amplification"
        )
    try:
        inverse = np.linalg.solve(array, np.broadcast_to(np.eye(2), array.shape))
    except np.linalg.LinAlgError as error:
        raise NumericalBreakdownError("singular Nambu recovery") from error
    with np.errstate(all="ignore"):
        residual = float(np.max(np.abs(array @ inverse - np.eye(2))))
        z = (inverse[..., 0, 0] + inverse[..., 1, 1]) / (
            2j * layout.frequencies_ev[:, None]
        )
        shift = (
            (inverse[..., 1, 1] - inverse[..., 0, 0]) / 2
            - layout.state_energies_ev[None, :]
            + layout.delta_mu_ev
        )
        gap = -(inverse[..., 0, 1] + inverse[..., 1, 0]) / (2 * z)
        relative = np.abs(inverse[..., 0, 1] - inverse[..., 1, 0]) / np.maximum(
            1, np.max(np.abs(inverse), axis=(-1, -2))
        )
    if (
        not all(np.all(np.isfinite(value)) for value in (inverse, z, shift, gap))
        or not math.isfinite(residual)
        or residual > 1e-10
    ):
        raise NumericalBreakdownError("Nambu recovery solve is inconclusive")
    discarded = [float(np.max(np.abs(value.imag))) for value in (z, gap, shift)]
    scales = [max(1.0, float(np.max(np.abs(value.real)))) for value in (z, gap, shift)]
    defect = float(np.max(relative))
    if defect > real_limit or any(
        value > real_limit * scale
        for value, scale in zip(discarded, scales, strict=True)
    ):
        raise NumericalBreakdownError(
            "Nambu recovery is materially non-real/asymmetric"
        )
    zero_shift_defect = 0.0
    real_shift = shift.real
    if not layout.full_bandwidth:
        zero_shift_defect = float(np.max(np.abs(real_shift)))
        if zero_shift_defect > real_limit:
            raise NumericalBreakdownError("recovery violates explicit non-FBW chi=0")
        real_shift = np.zeros_like(real_shift)
    output = EPWAnisotropicData(
        layout.frequencies_ev,
        layout.state_energies_ev,
        z.real,
        gap.real,
        real_shift,
        context=json.loads(layout.context_json),
    )
    diagnostics = EPWRecoveryDiagnostics(
        condition_max=float(np.max(condition)),
        inverse_derivative_gain_max=float(np.max(gain)),
        z_derivative_gain_max=float(np.max(z_gain)),
        solve_residual_max=residual,
        offdiagonal_relative_defect=defect,
        imaginary_z_max=discarded[0],
        imaginary_delta_ev_max=discarded[1],
        imaginary_chi_ev_max=discarded[2],
        non_fbw_shift_zero_defect_ev=zero_shift_defect,
    )
    return EPWRecoveryResult(output, diagnostics)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class EPWNambuFit:
    """A shared version-2 FermionicFit plus owning native layout/context.

    Archive source metadata is historical provenance, not a new execution of
    the upstream solver. The generic complex fitter never infers symmetry.
    """

    fit: FermionicFit
    _layout: _Layout
    fit_rows: IntArray
    completion_policy: str
    _analysis_json: str

    def __init__(
        self,
        fit: FermionicFit,
        *,
        frequencies_ev: npt.ArrayLike,
        state_energies_ev: npt.ArrayLike,
        context: Mapping[str, Any],
        fit_rows: npt.ArrayLike,
        completion_policy: str,
        analysis_provenance: Mapping[str, Any] | None = None,
    ) -> None:
        if not isinstance(fit, FermionicFit) or fit.observation_domain != "matsubara":
            raise ValueError(
                "EPW Nambu fits require the shared Matsubara fit interface"
            )
        layout = _Layout(frequencies_ev, state_energies_ev, context)
        if not math.isclose(
            fit.basis.beta, layout.beta, rel_tol=32 * np.finfo(float).eps
        ):
            raise ValueError("basis beta must equal the declared physical temperature")
        if fit.axis != 0 or fit.coefficients.shape != (
            fit.basis.rank,
            len(layout.state_energies_ev),
            2,
            2,
        ):
            raise ValueError("fit must preserve coefficient/state/2x2 axes")
        rows = _indices(fit_rows, "fit_rows")
        if np.any(rows < 0) or np.any(rows >= len(layout.indices)):
            raise ValueError("fit_rows must select native frequency rows")
        if completion_policy not in ("none", "explicit_adjoint_real_eliashberg"):
            raise ValueError("unknown explicit completion policy")
        if completion_policy != "none" and not np.all(layout.indices >= 0):
            raise ValueError("adjoint completion requires a positive-only source mesh")
        provenance = _json_object(analysis_provenance or {}, "analysis_provenance")
        object.__setattr__(self, "fit", fit)
        object.__setattr__(self, "_layout", layout)
        object.__setattr__(self, "fit_rows", rows)
        object.__setattr__(self, "completion_policy", completion_policy)
        object.__setattr__(self, "_analysis_json", provenance)

    @property
    def source_context(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._layout.context_json))

    @property
    def analysis_provenance(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._analysis_json))

    @property
    def indices(self) -> IntArray:
        return self._layout.indices

    def evaluate_matsubara(self, indices: npt.ArrayLike) -> ComplexArray:
        result = self.fit.evaluate_matsubara(indices)
        return cast(
            ComplexArray, _immutable_copy(np.asarray(result, dtype=np.complex128))
        )

    def reconstruct_native(
        self,
        *,
        max_condition: float = 1e10,
        max_z_amplification: float = 1e10,
        real_tolerance: float = 1e-8,
    ) -> EPWRecoveryResult:
        """Recover fields on the native mesh, with explicit diagnostics."""
        return _recover_layout(
            self._layout,
            self.evaluate_matsubara(self.indices),
            max_condition=max_condition,
            max_z_amplification=max_z_amplification,
            real_tolerance=real_tolerance,
        )

    def save(self, path: str | Path, *, overwrite: bool = False) -> None:
        """Save an owning v2 fit NPZ plus a hash-bound versioned application JSON."""
        if type(overwrite) is not bool:
            raise ValueError("overwrite must be a boolean")
        manifest = Path(path)
        archive = manifest.with_suffix(".npz")
        if manifest == archive:
            raise ValueError("manifest path must differ from its .npz component")
        if not overwrite and (manifest.exists() or archive.exists()):
            raise FileExistsError(manifest)
        self.fit.save(archive)
        payload = {
            "schema": _ARCHIVE_SCHEMA,
            "fit_archive": {
                "file": archive.name,
                "sha256": _digest(archive.read_bytes()),
            },
            "layout": {
                "context": self.source_context,
                "frequencies_ev": self._layout.frequencies_ev.tolist(),
                "state_energies_ev": self._layout.state_energies_ev.tolist(),
            },
            "fit_rows": self.fit_rows.tolist(),
            "completion_policy": self.completion_policy,
            "analysis_provenance": self.analysis_provenance,
        }
        payload["record_sha256"] = _digest(
            _json_object(payload, "archive record").encode()
        )
        manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )

    @classmethod
    def load(cls, path: str | Path) -> EPWNambuFit:
        """Check component integrity and conventions; never execute archive contents."""
        manifest = Path(path)
        record = _load_json(manifest)
        if (
            set(record)
            != {
                "schema",
                "fit_archive",
                "layout",
                "fit_rows",
                "completion_policy",
                "analysis_provenance",
                "record_sha256",
            }
            or record["schema"] != _ARCHIVE_SCHEMA
        ):
            raise ValueError("unsupported EPW Nambu archive manifest")
        expected_record_hash = record.pop("record_sha256")
        if expected_record_hash != _digest(
            _json_object(record, "archive record").encode()
        ):
            raise ValueError("EPW Nambu manifest metadata checksum mismatched")
        component = record["fit_archive"]
        if not isinstance(component, dict) or set(component) != {"file", "sha256"}:
            raise ValueError("invalid fit component record")
        filename = component["file"]
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise ValueError("fit component must be a local filename")
        archive = (manifest.parent / filename).resolve()
        if (
            archive.parent != manifest.parent.resolve()
            or _digest(archive.read_bytes()) != component["sha256"]
        ):
            raise ValueError(
                "fit component escaped manifest directory or hash mismatched"
            )
        layout = record["layout"]
        if not isinstance(layout, dict) or set(layout) != {
            "context",
            "frequencies_ev",
            "state_energies_ev",
        }:
            raise ValueError("invalid native layout record")
        fit = FermionicFit.load(archive)
        return cls(
            fit,
            frequencies_ev=layout["frequencies_ev"],
            state_energies_ev=layout["state_energies_ev"],
            context=layout["context"],
            fit_rows=record["fit_rows"],
            completion_policy=record["completion_policy"],
            analysis_provenance=record["analysis_provenance"],
        )


def fit_epw_nambu(
    data: EPWAnisotropicData,
    basis: FermionicBasis,
    *,
    fit_rows: npt.ArrayLike | None = None,
    rcond: float | None = None,
    complete_adjoint: bool = False,
    max_retained_condition: float = 1e12,
    analysis_provenance: Mapping[str, Any] | None = None,
) -> EPWNambuFit:
    """Fit source rows with optional explicitly requested adjoint completion."""
    if not isinstance(data, EPWAnisotropicData) or not isinstance(
        basis, FermionicBasis
    ):
        raise TypeError("EPWAnisotropicData and FermionicBasis are required")
    if type(complete_adjoint) is not bool:
        raise ValueError("complete_adjoint must be an explicit boolean")
    rows = _indices(
        np.arange(len(data.indices)) if fit_rows is None else fit_rows, "fit_rows"
    )
    if np.any(rows < 0) or np.any(rows >= len(data.indices)):
        raise ValueError("fit_rows must select source rows")
    if not math.isclose(basis.beta, data.beta, rel_tol=32 * np.finfo(float).eps):
        raise ValueError("basis beta must equal the declared physical temperature")
    modes = data.indices[rows]
    values = nambu_green(data)[rows]
    policy = "none"
    if complete_adjoint:
        if not np.all(data.indices >= 0):
            raise ValueError("adjoint completion requires a positive-only source mesh")
        modes = np.concatenate((modes, -modes - 1))
        values = np.concatenate((values, values.swapaxes(-1, -2).conj()))
        policy = "explicit_adjoint_real_eliashberg"
    condition_limit = _real(
        max_retained_condition, "max_retained_condition", positive=True
    )
    if condition_limit < 1:
        raise ValueError("max_retained_condition must be at least one")
    if len(modes) < basis.rank:
        raise NumericalBreakdownError("undersampled EPW Matsubara design")
    plan = basis.prepare_matsubara(modes, rcond=rcond)
    diagnostic = plan.diagnostics
    if (
        diagnostic.effective_rank == 0
        or diagnostic.retained_condition_number is None
        or diagnostic.retained_condition_number > condition_limit
    ):
        raise NumericalBreakdownError(
            "unsafe retained EPW Matsubara design conditioning"
        )
    fit = basis.fit_matsubara(values, modes, rcond=rcond)
    return EPWNambuFit(
        fit,
        frequencies_ev=data.frequencies_ev,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=rows,
        completion_policy=policy,
        analysis_provenance=analysis_provenance,
    )
