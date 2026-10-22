"""Shared, explicitly scoped output for spectral-atom proposal algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .types import FloatArray, NumericalMetadata


@dataclass(frozen=True, slots=True)
class AtomProposal:
    """Candidate real spectral frequencies, not a fitted or certified measure.

    ``metadata.converged`` refers only to the algorithm's native data domain.
    A caller must fit weights in the original kernel and assess its returned
    time-domain representation separately. Diagnostics preserve rejected poles
    and native errors; no complex or out-of-band pole is silently projected.
    """

    frequencies: FloatArray
    metadata: NumericalMetadata
    diagnostics: dict[str, Any]
