"""Explicit real-Nambu parity fitting with a frequency-native Cauchy cross.

This EPW-specific real-linear map is independent of the generic complex-linear
FermionicFit interface. Its safe numeric archive uses the shared transport with
an explicitly different representation schema; existing fit archives are intact.
"""

from __future__ import annotations

import json
import math
import time
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import numpy.typing as npt

from .._basis_io import read_archive, write_archive
from ..basis import _indices_array
from ..fitting import LeastSquaresDiagnostics, _immutable_copy
from ..types import FloatArray, NumericalBreakdownError
from ._frequency_cross import CauchyFrequencyPlan
from .epw import (
    ComplexArray,
    EPWAnisotropicData,
    EPWRecoveryResult,
    IntArray,
    _array,
    _digest,
    _indices,
    _json_object,
    _Layout,
    _load_json,
    _real,
    _recover_layout,
    nambu_green,
)

_FORMAT = "kernelgecp.epw_nambu_parity_cross"
_SCHEMA = "kernelgecp.epw-parity-cross-application.v1"
_POLICY = "explicit_real_eliashberg_nambu_parity"
_REPRESENTATION_UNITS = {
    "observation_x_nodes_poles": "eV^2",
    "cardinal_weights": "eV^2",
    "channels_and_coefficients": ["eV^-1", "eV^-2", "eV^-1"],
    "green_function": "eV^-1",
    "z": "dimensionless",
}


@dataclass(frozen=True, slots=True, init=False, eq=False)
class EPWFrequencyCrossFit:
    """Immutable Green-channel coefficients and owned native state layout.

    Channels are Re G11, Im G11/omega and Re G12. The declared real Nambu
    symmetry reconstructs G22=-Re G11+i Im G11 and G21=G12. This operation is
    explicit and cannot represent arbitrary complex Nambu data.
    """

    plan: CauchyFrequencyPlan
    coefficients: FloatArray
    fit_rows: IntArray
    _layout: _Layout
    _diagnostics_json: str
    _provenance_json: str

    def __init__(
        self,
        plan: CauchyFrequencyPlan,
        coefficients: npt.ArrayLike,
        *,
        frequencies_ev: npt.ArrayLike,
        state_energies_ev: npt.ArrayLike,
        context: Mapping[str, Any],
        fit_rows: npt.ArrayLike,
        diagnostics: Mapping[str, Any],
        analysis_provenance: Mapping[str, Any] | None = None,
    ):
        if not isinstance(plan, CauchyFrequencyPlan):
            raise TypeError("plan must be a CauchyFrequencyPlan")
        layout = _Layout(frequencies_ev, state_energies_ev, context)
        rows = _indices(fit_rows, "fit_rows")
        if (
            not json.loads(layout.context_json)["positive_only"]
            or np.any(rows < 0)
            or np.any(rows >= len(layout.indices))
        ):
            raise ValueError("real-parity cross fitting requires positive native rows")
        expected = ((2 * layout.indices[rows] + 1) * (np.pi / layout.beta)) ** 2
        if plan.observation_x.shape != expected.shape or not np.allclose(
            plan.observation_x, expected, rtol=32 * np.finfo(float).eps, atol=0
        ):
            raise ValueError(
                "cross observation coordinates must match native indices/beta"
            )
        coeff = _array(coefficients, "coefficients", 3)
        if coeff.shape != (plan.rank, len(layout.state_energies_ev), 3):
            raise ValueError(
                "coefficients must preserve cross/state/three-channel axes"
            )
        object.__setattr__(self, "plan", plan)
        object.__setattr__(
            self, "coefficients", cast(FloatArray, _immutable_copy(coeff))
        )
        object.__setattr__(self, "fit_rows", rows)
        object.__setattr__(self, "_layout", layout)
        object.__setattr__(
            self, "_diagnostics_json", _json_object(diagnostics, "diagnostics")
        )
        object.__setattr__(
            self,
            "_provenance_json",
            _json_object(analysis_provenance or {}, "analysis_provenance"),
        )

    @property
    def indices(self) -> IntArray:
        return self._layout.indices

    @property
    def source_context(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._layout.context_json))

    @property
    def diagnostics(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._diagnostics_json))

    @property
    def analysis_provenance(self) -> dict[str, Any]:
        return cast(dict[str, Any], json.loads(self._provenance_json))

    def evaluate_matsubara(self, indices: npt.ArrayLike) -> ComplexArray:
        modes = _indices_array(indices)
        if modes.size * self.coefficients.shape[1] * 4 * 16 > 64 * 1024**2:
            raise ValueError("Nambu evaluation exceeds 64 MiB output resource limit")
        with np.errstate(all="ignore"):
            omega = (2 * modes + 1).astype(float) * (np.pi / self._layout.beta)
            squared = omega * omega
        if not np.all(np.isfinite(squared)) or np.any(squared <= 0):
            raise ValueError("unrepresentable Matsubara coordinates")
        if not modes.size:
            result = np.empty((0, self.coefficients.shape[1], 2, 2), complex)
        else:
            design = self.plan.design(squared)
            channels = np.einsum("ij,jsc->isc", design, self.coefficients)
            result = np.empty((*channels.shape[:2], 2, 2), complex)
            real, imaginary, anomalous = (
                channels[..., 0],
                omega[:, None] * channels[..., 1],
                channels[..., 2],
            )
            result[..., 0, 0] = real + 1j * imaginary
            result[..., 1, 1] = -real + 1j * imaginary
            result[..., 0, 1] = result[..., 1, 0] = anomalous
        if not np.all(np.isfinite(result)):
            raise NumericalBreakdownError("nonfinite real-parity Nambu evaluation")
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
        return _recover_layout(
            self._layout,
            self.evaluate_matsubara(self.indices),
            max_condition=max_condition,
            max_z_amplification=max_z_amplification,
            real_tolerance=real_tolerance,
        )

    def save(self, path: str | Path, *, overwrite: bool = False) -> None:
        """Write shared bounded numeric transport and hash-bound EPW JSON.

        The version-2 parity representation is loaded by this facade; it is
        deliberately distinguished from the time-cardinal FermionicFit layout.
        """
        if type(overwrite) is not bool:
            raise ValueError("overwrite must be a boolean")
        manifest = Path(path)
        archive = manifest.with_suffix(".npz")
        if manifest == archive:
            raise ValueError("application manifest must differ from .npz component")
        if not overwrite and (manifest.exists() or archive.exists()):
            raise FileExistsError(manifest)
        plan_metadata, arrays = self.plan.to_payload()
        arrays.update(
            coefficients=self.coefficients,
            fit_rows=self.fit_rows,
            frequencies_ev=self._layout.frequencies_ev,
            state_energies_ev=self._layout.state_energies_ev,
        )
        metadata: dict[str, object] = {
            "format": _FORMAT,
            "version": 2,
            "observation_domain": "matsubara",
            "completion_policy": _POLICY,
            "representation_units": _REPRESENTATION_UNITS,
            "plan": plan_metadata,
            "diagnostics": self.diagnostics,
            "context": self.source_context,
            "analysis_provenance": self.analysis_provenance,
        }
        write_archive(archive, metadata, arrays)
        record = {
            "schema": _SCHEMA,
            "archive": {"file": archive.name, "sha256": _digest(archive.read_bytes())},
            "representation": _POLICY,
            "units": "eV",
            "representation_units": _REPRESENTATION_UNITS,
            "provider": self.source_context["provider"],
            "state_order": self.source_context["state_order"],
            "source_context": self.source_context,
            "analysis_provenance": self.analysis_provenance,
        }
        record["record_sha256"] = _digest(_json_object(record, "record").encode())
        manifest.write_text(
            json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )

    @classmethod
    def load(cls, path: str | Path) -> EPWFrequencyCrossFit:
        manifest = Path(path)
        record = _load_json(manifest)
        if set(record) != {
            "schema",
            "archive",
            "representation",
            "units",
            "representation_units",
            "provider",
            "state_order",
            "source_context",
            "analysis_provenance",
            "record_sha256",
        }:
            raise ValueError("invalid parity-cross application manifest keys")
        expected = record.pop("record_sha256")
        if (
            record["schema"] != _SCHEMA
            or record["representation"] != _POLICY
            or record["units"] != "eV"
            or record["representation_units"] != _REPRESENTATION_UNITS
            or expected != _digest(_json_object(record, "record").encode())
        ):
            raise ValueError("unsupported parity-cross manifest or checksum mismatch")
        component = record["archive"]
        if not isinstance(component, dict) or set(component) != {"file", "sha256"}:
            raise ValueError("invalid parity-cross archive component")
        filename = component["file"]
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise ValueError("parity-cross component must be a local filename")
        archive = (manifest.parent / filename).resolve()
        if archive.stat().st_size > 65 * 1024**2:
            raise ValueError("parity-cross archive exceeds 65 MiB resource limit")
        if (
            archive.parent != manifest.parent.resolve()
            or _digest(archive.read_bytes()) != component["sha256"]
        ):
            raise ValueError("parity-cross component escaped or hash mismatched")
        metadata, arrays = read_archive(archive)
        if (
            set(metadata)
            != {
                "format",
                "version",
                "observation_domain",
                "completion_policy",
                "representation_units",
                "plan",
                "diagnostics",
                "context",
                "analysis_provenance",
            }
            or metadata["format"] != _FORMAT
            or type(metadata["version"]) is not int
            or metadata["version"] != 2
            or metadata["observation_domain"] != "matsubara"
            or metadata["completion_policy"] != _POLICY
            or metadata["representation_units"] != _REPRESENTATION_UNITS
            or metadata["context"] != record["source_context"]
            or metadata["analysis_provenance"] != record["analysis_provenance"]
        ):
            raise ValueError("unsupported/inconsistent parity-cross archive metadata")
        required = {
            "observation_x",
            "nodes",
            "poles",
            "weights",
            "coefficients",
            "fit_rows",
            "frequencies_ev",
            "state_energies_ev",
        }
        if set(arrays) != required:
            raise ValueError("unsupported parity-cross archive arrays")
        plan = CauchyFrequencyPlan.from_payload(
            cast(dict[str, Any], metadata["plan"]),
            {k: arrays[k] for k in ("observation_x", "nodes", "poles", "weights")},
        )
        fit = cls(
            plan,
            arrays["coefficients"],
            frequencies_ev=arrays["frequencies_ev"],
            state_energies_ev=arrays["state_energies_ev"],
            context=cast(dict[str, Any], metadata["context"]),
            fit_rows=arrays["fit_rows"],
            diagnostics=cast(dict[str, Any], metadata["diagnostics"]),
            analysis_provenance=cast(dict[str, Any], metadata["analysis_provenance"]),
        )
        if (
            record["provider"] != fit.source_context["provider"]
            or record["state_order"] != fit.source_context["state_order"]
        ):
            raise ValueError(
                "parity-cross manifest provider/state ordering inconsistent"
            )
        return fit


def fit_epw_frequency_cross(
    data: EPWAnisotropicData,
    *,
    real_nambu_parity: bool,
    fit_rows: npt.ArrayLike | None = None,
    rank: int = 60,
    pole_cutoff_ev: float = 13,
    dictionary_size: int = 1201,
    precision_bits: int = 256,
    rcond: float | None = None,
    max_evaluation_condition: float = 1e4,
    max_retained_condition: float = 1e8,
    analysis_provenance: Mapping[str, Any] | None = None,
) -> EPWFrequencyCrossFit:
    """Fit actual Green channels with explicitly requested real Nambu parity.

    Only supplied frequency rows construct Green matrices. The dictionary
    cutoff is not a physical-support certificate; validate held-out native data.
    No self-energy subtraction, tail correction or another fitting method is
    applied. Ill-conditioned selection/evaluation fails without substitution.
    """
    if type(real_nambu_parity) is not bool or not real_nambu_parity:
        raise ValueError("explicit real_nambu_parity=True is required")
    if not isinstance(data, EPWAnisotropicData):
        raise TypeError("data must be EPWAnisotropicData")
    if not data.context["positive_only"]:
        raise ValueError("real-parity frequency cross supports positive native meshes")
    rows = (
        np.arange(len(data.indices))
        if fit_rows is None
        else _indices(fit_rows, "fit_rows")
    )
    if np.any(rows < 0) or np.any(rows >= len(data.indices)):
        raise ValueError("fit_rows must select native frequencies")
    limit = _real(max_retained_condition, "max_retained_condition", positive=True)
    if limit < 1:
        raise ValueError("max_retained_condition must be at least1")
    omega = (2 * data.indices[rows] + 1) * (np.pi / data.beta)
    start = time.perf_counter()
    plan = CauchyFrequencyPlan(
        omega**2,
        pole_cutoff_ev=pole_cutoff_ev,
        rank=rank,
        dictionary_size=dictionary_size,
        precision_bits=precision_bits,
        max_evaluation_condition=max_evaluation_condition,
    )
    constructed = time.perf_counter()
    green = nambu_green(data, rows=rows)
    channels = np.stack(
        (
            green[..., 0, 0].real,
            green[..., 0, 0].imag / omega[:, None],
            green[..., 0, 1].real,
        ),
        axis=-1,
    )
    observed = time.perf_counter()
    prepared = plan.prepare(rcond=rcond)
    prepared_at = time.perf_counter()
    diagnostic: LeastSquaresDiagnostics = prepared.diagnostics
    condition = diagnostic.retained_condition_number
    if (
        diagnostic.effective_rank == 0
        or condition is None
        or not math.isfinite(condition)
        or condition > limit
    ):
        raise NumericalBreakdownError(
            "unsafe frequency-cross retained SVD conditioning"
        )
    coefficients = plan.fit(channels, rcond=rcond)
    fitted_at = time.perf_counter()
    predicted = np.einsum("ij,jsc->isc", plan.design(omega**2), coefficients)
    diagnostics = {
        "least_squares": asdict(diagnostic),
        "training_channel_max_residual": float(np.max(np.abs(predicted - channels))),
        "parity_policy": _POLICY,
        "observation_domain": "matsubara",
        "channel_order": ["ReG11", "ImG11_over_omega", "ReG12"],
        "native_frequency_rounding": (
            "native_G_values_fitted_at_canonical_ordinary_indices"
        ),
        "selection_scope": "fixed_rank_finite_dictionary_log_schur_gecp",
        "construction_precision_bits": precision_bits,
        "evaluation_precision_bits": 53,
        "certification_status": "not_certified",
        "timing_seconds": {
            "basis_construction": constructed - start,
            "nambu_observation_construction": observed - constructed,
            "least_squares_preparation": prepared_at - observed,
            "coefficient_fit": fitted_at - prepared_at,
            "training_replay": time.perf_counter() - fitted_at,
        },
    }
    return EPWFrequencyCrossFit(
        plan,
        coefficients,
        frequencies_ev=data.frequencies_ev,
        state_energies_ev=data.state_energies_ev,
        context=data.context,
        fit_rows=rows,
        diagnostics=diagnostics,
        analysis_provenance=analysis_provenance,
    )
