"""Explicit model plus lossless native-field XOR restart representation.

The corrections encode all differences from a deterministic binary64 model;
they are not evidence that the uncorrected model passes a native restart.
SHA-256 checks bind contents for integrity; they do not authenticate a producer.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np

from .._basis_io import read_archive, write_archive
from ..fitting import _immutable_copy
from ..types import NumericalBreakdownError
from ._frequency_cross import CauchyFrequencyPlan
from .epw import (
    EPWAnisotropicData,
    _digest,
    _json_object,
    _load_json,
    write_epw_anisotropic,
)
from .epw_frequency import EPWFrequencyCrossFit

_POLICY = "lossless_native_xor_v1"
_SCHEMA = "kernelgecp.epw-restart-xor-application.v1"
_FORMAT = "kernelgecp.epw-restart-model-xor.v1"
_MAX_BYTES = 64 * 1024**2
_MAX_WORK = 50_000_000
# Local derivatives of gap in eV with respect to c (eV^-1) and b
# (eV^-2), respectively; these are diagnostic caps, not certificates.
_GAP_C_GAIN_LIMIT_EV2 = 1e10
_GAP_B_GAIN_LIMIT_EV3 = 1e10
_PLAN_KEYS = ("observation_x", "nodes", "poles", "weights")
_ARRAY_KEYS = set(_PLAN_KEYS) | {
    "coefficients",
    "fit_rows",
    "frequencies_ev",
    "state_energies_ev",
    "corrections",
}


def _bits(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype="<f8", order="C").view("<u8")


def _hash(values: np.ndarray) -> str:
    return _digest(np.asarray(values, dtype="<u8", order="C").tobytes())


def _owned_fit(fit: EPWFrequencyCrossFit) -> EPWFrequencyCrossFit:
    metadata, arrays = fit.plan.to_payload()
    plan = CauchyFrequencyPlan.from_payload(metadata, arrays)
    return EPWFrequencyCrossFit(
        plan,
        fit.coefficients,
        frequencies_ev=fit._layout.frequencies_ev,
        state_energies_ev=fit._layout.state_energies_ev,
        context=fit.source_context,
        fit_rows=fit.fit_rows,
        diagnostics=fit.diagnostics,
        analysis_provenance=fit.analysis_provenance,
    )


def _baseline(fit: EPWFrequencyCrossFit) -> np.ndarray:
    """Use scalar Python operations and fsum; no BLAS, FMA, or reselection."""
    n, states, rank = (
        len(fit.indices),
        len(fit._layout.state_energies_ev),
        fit.plan.rank,
    )
    if n * states * 3 * 8 > _MAX_BYTES or n * states * rank * 3 > _MAX_WORK:
        raise ValueError("restart evaluation exceeds resource limit")
    plan = fit.plan
    nodes = [float(x) for x in plan.nodes]
    weights = [float(x) for x in plan.weights]
    maximum = max(abs(x) for x in weights)
    normalized = [x / maximum for x in weights]
    if any(x == 0 for x in normalized):
        raise NumericalBreakdownError("scaled cardinal weights underflow")
    result = np.empty((n, states, 3), dtype="<f8")
    for i, mode in enumerate(fit.indices):
        omega = float(2 * int(mode) + 1) * (math.pi / fit._layout.beta)
        x = omega * omega
        if not math.isfinite(x) or x <= 0:
            raise NumericalBreakdownError("nonfinite frequency coordinate")
        delta = [x - node for node in nodes]
        if 0.0 in delta:
            design = [0.0] * rank
            design[delta.index(0.0)] = 1.0
        else:
            nearest = min(abs(d) for d in delta)
            terms = [w * (nearest / d) for w, d in zip(normalized, delta, strict=True)]
            constant = nearest / maximum
            den = math.fsum([constant, *terms])
            num = math.fsum([abs(constant), *(abs(t) for t in terms)])
            if not math.isfinite(den) or den == 0 or num == 0:
                raise NumericalBreakdownError("invalid cardinal denominator")
            if constant == 0:
                tiny = math.ulp(0.0)
                lower = math.nextafter(abs(den) - tiny, 0.0)
                upper = math.nextafter(num + tiny, math.inf)
                if lower <= 0 or upper / lower > plan.max_evaluation_condition:
                    raise NumericalBreakdownError(
                        "underflow prevents condition assessment"
                    )
            condition = num / abs(den)
            if (
                not math.isfinite(condition)
                or condition > plan.max_evaluation_condition
            ):
                raise NumericalBreakdownError(
                    "cardinal evaluation condition limit exceeded"
                )
            design = [t / den for t in terms]
        for j in range(states):
            channels = [
                math.fsum(
                    design[k] * float(fit.coefficients[k, j, c]) for k in range(rank)
                )
                for c in range(3)
            ]
            a, b, c = channels
            wb = omega * b
            magnitude2 = math.fsum([a * a, wb * wb, c * c])
            if not math.isfinite(magnitude2) or magnitude2 <= 0 or b == 0:
                raise NumericalBreakdownError("unsafe Nambu recovery")
            # The real-parity matrix is normal with equal singular values.
            native_omega = float(fit._layout.frequencies_ev[i])
            gain = (1.0 / magnitude2) / abs(native_omega)
            if not math.isfinite(gain) or gain > 1e10:
                raise NumericalBreakdownError("unsafe Nambu recovery amplification")
            ratio = native_omega / omega
            b_abs = abs(b)
            b_squared = b_abs * b_abs
            if b_squared == 0 or not math.isfinite(b_squared):
                raise NumericalBreakdownError("gap sensitivity underflow/overflow")
            c_gain = ratio / b_abs
            b_gain = (ratio * abs(c)) / b_squared
            if (
                not math.isfinite(c_gain)
                or not math.isfinite(b_gain)
                or c_gain > _GAP_C_GAIN_LIMIT_EV2
                or b_gain > _GAP_B_GAIN_LIMIT_EV3
                or (c != 0 and b_gain == 0)
            ):
                raise NumericalBreakdownError("unsafe gap sensitivity")
            det = -magnitude2
            z = (omega / native_omega) * (b / det)
            gap = ratio * (c / b)
            chi = (
                a / det - float(fit._layout.state_energies_ev[j])
            ) + fit._layout.delta_mu_ev
            if not all(math.isfinite(v) for v in (z, gap, chi)) or z <= 0:
                raise NumericalBreakdownError(
                    "nonfinite or nonpositive recovered native fields"
                )
            if not fit._layout.full_bandwidth:
                if abs(chi) > 1e-8:
                    raise NumericalBreakdownError(
                        "recovery violates non-FBW chi convention"
                    )
                chi = 0.0
            result[i, j] = (z, gap, chi)
    return result


@dataclass(frozen=True, slots=True, init=False, eq=False)
class EPWRestartArchive:
    """Immutable full model and XOR field corrections; not a pure fitted restart."""

    fit: EPWFrequencyCrossFit
    _corrections: np.ndarray
    _baseline_sha256: str
    _target_sha256: str

    def __init__(
        self,
        fit: EPWFrequencyCrossFit,
        corrections: np.ndarray,
        baseline_sha256: str,
        target_sha256: str,
    ):
        owned = _owned_fit(fit)
        shape = (len(owned.indices), len(owned._layout.state_energies_ev), 3)
        if (
            not isinstance(corrections, np.ndarray)
            or corrections.dtype != np.dtype("<u8")
            or corrections.shape != shape
            or corrections.nbytes > _MAX_BYTES
        ):
            raise ValueError("invalid uint64 correction shape or dtype")
        for digest in (baseline_sha256, target_sha256):
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)
            ):
                raise ValueError("invalid field SHA256")
        object.__setattr__(self, "fit", owned)
        object.__setattr__(self, "_corrections", _immutable_copy(corrections))
        object.__setattr__(self, "_baseline_sha256", baseline_sha256)
        object.__setattr__(self, "_target_sha256", target_sha256)

    @classmethod
    def from_data(
        cls, data: EPWAnisotropicData, fit: EPWFrequencyCrossFit, *, policy: str
    ) -> EPWRestartArchive:
        if policy != _POLICY:
            raise ValueError("explicit lossless_native_xor_v1 policy required")
        if not isinstance(data, EPWAnisotropicData) or not isinstance(
            fit, EPWFrequencyCrossFit
        ):
            raise TypeError("EPW data and frequency cross fit required")
        if (
            data.context != fit.source_context
            or not np.array_equal(data.frequencies_ev, fit._layout.frequencies_ev)
            or not np.array_equal(data.state_energies_ev, fit._layout.state_energies_ev)
        ):
            raise ValueError("data and fit native layouts differ")
        base = _bits(_baseline(fit))
        target = _bits(np.stack((data.z, data.delta_ev, data.chi_ev), axis=-1))
        return cls(fit, np.bitwise_xor(base, target), _hash(base), _hash(target))

    def reconstruct_native(self) -> EPWAnisotropicData:
        base = _bits(_baseline(self.fit))
        if _hash(base) != self._baseline_sha256:
            raise ValueError(
                "canonical baseline SHA256 mismatch; platform or model differs"
            )
        target = np.bitwise_xor(base, self._corrections)
        if _hash(target) != self._target_sha256:
            raise ValueError("target native field SHA256 mismatch")
        fields = target.view("<f8")
        return EPWAnisotropicData(
            self.fit._layout.frequencies_ev,
            self.fit._layout.state_energies_ev,
            fields[..., 0],
            fields[..., 1],
            fields[..., 2],
            context=self.fit.source_context,
        )

    def write_native_restart(
        self, path: str | Path, *, overwrite: bool = False
    ) -> dict[str, Any]:
        """Write native ES20.10 text; arbitrary binary64 values may round in text."""
        return write_epw_anisotropic(
            self.reconstruct_native(), path, preserve_source=False, overwrite=overwrite
        )

    def save(self, path: str | Path, *, overwrite: bool = False) -> None:
        if type(overwrite) is not bool or overwrite:
            raise ValueError("restart archive overwrite is unsupported")
        manifest = Path(path)
        archive = manifest.with_suffix(".npz")
        if manifest == archive:
            raise ValueError("manifest must differ from archive")
        if manifest.exists() or archive.exists():
            raise FileExistsError(manifest)
        plan, arrays = self.fit.plan.to_payload()
        arrays.update(
            coefficients=self.fit.coefficients,
            fit_rows=self.fit.fit_rows,
            frequencies_ev=self.fit._layout.frequencies_ev,
            state_energies_ev=self.fit._layout.state_energies_ev,
            corrections=self._corrections,
        )
        metadata: dict[str, object] = {
            "format": _FORMAT,
            "policy": _POLICY,
            "plan": plan,
            "context": self.fit.source_context,
            "diagnostics": self.fit.diagnostics,
            "analysis_provenance": self.fit.analysis_provenance,
            "baseline_sha256": self._baseline_sha256,
            "target_sha256": self._target_sha256,
        }
        # Prepare both files on the destination filesystem. Hard links publish
        # exclusively, so existing pairs cannot be overwritten by a racing save.
        with tempfile.TemporaryDirectory(
            prefix=".epw-restart-", dir=manifest.parent
        ) as directory:
            staged_archive = Path(directory) / archive.name
            staged_manifest = Path(directory) / manifest.name
            write_archive(staged_archive, metadata, arrays)
            record = {
                "schema": _SCHEMA,
                "policy": _POLICY,
                "units": "eV",
                "provider": self.fit.source_context["provider"],
                "state_order": self.fit.source_context["state_order"],
                "extra_encoded_information": (
                    "uint64 XOR corrections for all native Z, gap, and chi fields"
                ),
                "archive": {
                    "file": archive.name,
                    "sha256": _digest(staged_archive.read_bytes()),
                },
            }
            record["record_sha256"] = _digest(_json_object(record, "record").encode())
            staged_manifest.write_text(
                json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + "\n"
            )
            published_archive = False
            try:
                os.link(staged_archive, archive)
                published_archive = True
                os.link(staged_manifest, manifest)
            except BaseException:
                if published_archive:
                    archive.unlink(missing_ok=True)
                raise

    @classmethod
    def load(cls, path: str | Path) -> EPWRestartArchive:
        manifest = Path(path)
        record = _load_json(manifest)
        if set(record) != {
            "schema",
            "policy",
            "units",
            "provider",
            "state_order",
            "extra_encoded_information",
            "archive",
            "record_sha256",
        }:
            raise ValueError("invalid restart manifest keys")
        checksum = record.pop("record_sha256")
        if (
            record["schema"] != _SCHEMA
            or record["policy"] != _POLICY
            or record["units"] != "eV"
            or record["extra_encoded_information"]
            != "uint64 XOR corrections for all native Z, gap, and chi fields"
            or checksum != _digest(_json_object(record, "record").encode())
        ):
            raise ValueError("unsupported restart manifest or checksum mismatch")
        component = record["archive"]
        if not isinstance(component, dict) or set(component) != {"file", "sha256"}:
            raise ValueError("invalid archive component")
        filename = component["file"]
        if not isinstance(filename, str) or Path(filename).name != filename:
            raise ValueError("archive component must be local filename")
        archive = (manifest.parent / filename).resolve()
        if (
            archive.parent != manifest.parent.resolve()
            or archive.stat().st_size > 65 * 1024**2
            or _digest(archive.read_bytes()) != component["sha256"]
        ):
            raise ValueError(
                "archive escaped, exceeds resource limit, or hash mismatched"
            )
        metadata, arrays = read_archive(archive)
        if (
            set(metadata)
            != {
                "format",
                "policy",
                "plan",
                "context",
                "diagnostics",
                "analysis_provenance",
                "baseline_sha256",
                "target_sha256",
            }
            or metadata["format"] != _FORMAT
            or metadata["policy"] != _POLICY
            or set(arrays) != _ARRAY_KEYS
        ):
            raise ValueError("invalid restart archive schema")
        for key in ("coefficients", "frequencies_ev", "state_energies_ev"):
            if arrays[key].dtype != np.dtype("float64"):
                raise ValueError(f"{key} must be binary64")
        if arrays["fit_rows"].dtype != np.dtype("int64"):
            raise ValueError("fit_rows must be int64")
        plan = CauchyFrequencyPlan.from_payload(
            cast(dict[str, Any], metadata["plan"]), {k: arrays[k] for k in _PLAN_KEYS}
        )
        fit = EPWFrequencyCrossFit(
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
            raise ValueError("restart manifest context mismatch")
        result = cls(
            fit,
            arrays["corrections"],
            cast(str, metadata["baseline_sha256"]),
            cast(str, metadata["target_sha256"]),
        )
        result.reconstruct_native()
        return result
