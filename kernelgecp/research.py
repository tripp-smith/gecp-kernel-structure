"""Claim-disciplined diagnostics for the open fermionic GECP rate problem."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal, localcontext
from pathlib import Path

import mpmath as mp
import numpy as np

from .certified import interval_certified_pivot
from .gecp import _FermionicResidualEvaluator
from .types import PivotCertificate


@dataclass(frozen=True, slots=True)
class BlockContractionEvidence:
    """Finite-grid evidence for cutoff-scaled residual block contraction.

    ``continuous_certified`` is deliberately always false here.  The canonical
    census uses high-precision finite grids; it does not certify a continuous
    GECP trajectory.
    """

    cutoff: str
    tolerance: str
    scale: int
    block_size: int
    first_step_observed: str
    first_step_exact: str
    complete_block_ratios: tuple[str, ...]
    trailing_steps: int
    trailing_ratio: str | None
    maximum_complete_block_ratio: str | None
    continuous_certified: bool = False


@dataclass(frozen=True, slots=True)
class CertifiedBlockContractionEvidence:
    """Continuous interval evidence for one approximate-pivot trajectory."""

    cutoff: float
    requested_steps: int
    completed_steps: int
    initial_lower_bound: float
    final_upper_bound: float
    contraction_ratio_upper: float
    half_contraction_certified: bool
    trajectory_certified: bool
    stop_reason: str
    pivot_certificates: tuple[PivotCertificate, ...]


def dyadic_cutoff_scale(cutoff: Decimal) -> int:
    """Return the least ``s`` with ``cutoff <= 2**s`` for ``cutoff >= 1``."""

    if not cutoff.is_finite() or cutoff < 1:
        raise ValueError("cutoff must be finite and at least one")
    scale = 0
    bound = Decimal(1)
    while bound < cutoff:
        bound *= 2
        scale += 1
    return scale


def fermionic_first_step_ratio(cutoff: Decimal, *, precision: int = 80) -> Decimal:
    """Return the exact-form first-step ratio ``1 - exp(-2*cutoff)``.

    This is the reflected-corner residual divided by the initial complete
    pivot.  It tends to one with the cutoff, obstructing a fixed one-step
    contraction factor.
    """

    if not cutoff.is_finite() or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if precision < 16:
        raise ValueError("precision must be at least 16 decimal digits")
    with localcontext() as context:
        context.prec = precision
        return +(Decimal(1) - (-2 * cutoff).exp())


def certify_fermionic_block_contraction(
    cutoff: float,
    steps: int,
    *,
    precision_bits: int = 128,
    abs_tol: float = 1e-8,
    rel_tol: float = 1e-6,
    max_cells: int = 100_000,
) -> CertifiedBlockContractionEvidence:
    """Certify a continuous residual block along interval-certified pivots.

    The result concerns the returned approximate-pivot trajectory.  It does
    not assert that every exact GECP tie choice follows the same trajectory.
    Unlike the production GECP result builder, this focused diagnostic omits
    core SVD computations and certifies only the residual suprema needed by
    the contraction question.
    """

    if not math.isfinite(cutoff) or cutoff < 1:
        raise ValueError("cutoff must be finite and at least one")
    if steps < 1:
        raise ValueError("steps must be positive")
    if precision_bits < 53 or abs_tol < 0 or rel_tol < 0 or max_cells < 1:
        raise ValueError("invalid certificate precision, tolerance, or budget")

    evaluator = _FermionicResidualEvaluator(cutoff)
    certificates: list[PivotCertificate] = []
    initial_lower = 0.0
    final_upper = math.inf
    completed = 0
    stop_reason = "requested_steps"
    trajectory_certified = True

    with mp.workprec(precision_bits):
        for step in range(steps + 1):
            certificate = interval_certified_pivot(
                evaluator.point,
                evaluator.interval,
                omega_bounds=(-cutoff, cutoff),
                precision_bits=precision_bits,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
                max_cells=max_cells,
                gradient_interval_function=evaluator.interval_gradient,
            )
            certificates.append(certificate)
            if step == 0:
                initial_lower = certificate.lower_bound
            if step == steps:
                final_upper = certificate.upper_bound
                stop_reason = certificate.termination_reason
                break
            if not certificate.certified:
                trajectory_certified = False
                final_upper = certificate.upper_bound
                stop_reason = certificate.termination_reason
                break
            t = mp.mpf(repr(certificate.t))
            omega = mp.mpf(repr(certificate.omega))
            pivot = evaluator.point(t, omega)
            if not mp.isfinite(pivot) or pivot == 0:
                trajectory_certified = False
                final_upper = certificate.upper_bound
                stop_reason = "unusable_pivot"
                break
            evaluator.add_pivot(t, omega, pivot)
            completed += 1

    ratio = (
        float(np.nextafter(final_upper / initial_lower, math.inf))
        if initial_lower > 0 and math.isfinite(final_upper)
        else math.inf
    )
    return CertifiedBlockContractionEvidence(
        cutoff=cutoff,
        requested_steps=steps,
        completed_steps=completed,
        initial_lower_bound=initial_lower,
        final_upper_bound=final_upper,
        contraction_ratio_upper=ratio,
        half_contraction_certified=trajectory_certified
        and completed == steps
        and ratio <= 0.5,
        trajectory_certified=trajectory_certified and completed == steps,
        stop_reason=stop_reason,
        pivot_certificates=tuple(certificates),
    )


def _decimal_sequence(value: object, *, field: str) -> list[Decimal]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{field} must be a sequence")
    result: list[Decimal] = []
    for item in value:
        try:
            parsed = Decimal(str(item))
        except Exception as error:  # pragma: no cover - Decimal error taxonomy varies
            raise ValueError(f"{field} contains a non-decimal value") from error
        if not parsed.is_finite():
            raise ValueError(f"{field} contains a non-finite value")
        result.append(abs(parsed))
    return result


def analyze_census_record(
    record: Mapping[str, object], *, block_multiplier: int = 2, precision: int = 80
) -> BlockContractionEvidence:
    """Analyze one canonical census record without promoting it to a theorem."""

    if block_multiplier < 1:
        raise ValueError("block_multiplier must be positive")
    try:
        cutoff_text = str(record["cutoff"])
        tolerance_text = str(record["tolerance"])
        residual = abs(Decimal(str(record["residual"])))
    except (KeyError, ArithmeticError) as error:
        raise ValueError(
            "record is missing valid cutoff, tolerance, or residual"
        ) from error
    cutoff = Decimal(cutoff_text)
    if not residual.is_finite():
        raise ValueError("residual must be finite")
    pivots = _decimal_sequence(record.get("pivots"), field="pivots")
    if not pivots or pivots[0] == 0:
        raise ValueError("record must contain a nonzero first pivot")
    history = [*pivots, residual]
    if any(value == 0 for value in history[:-1]):
        raise ValueError("preterminal residual history must be nonzero")

    scale = dyadic_cutoff_scale(cutoff)
    block_size = block_multiplier * (scale + 1)
    first_observed = history[1] / history[0]
    exact_ratio = fermionic_first_step_ratio(cutoff, precision=precision)

    full_blocks = (len(history) - 1) // block_size
    complete_ratios = tuple(
        history[(block + 1) * block_size] / history[block * block_size]
        for block in range(full_blocks)
    )
    trailing_steps = (len(history) - 1) - full_blocks * block_size
    trailing_ratio = (
        history[-1] / history[full_blocks * block_size] if trailing_steps else None
    )
    maximum = max(complete_ratios) if complete_ratios else None
    return BlockContractionEvidence(
        cutoff=cutoff_text,
        tolerance=tolerance_text,
        scale=scale,
        block_size=block_size,
        first_step_observed=str(first_observed),
        first_step_exact=str(exact_ratio),
        complete_block_ratios=tuple(str(value) for value in complete_ratios),
        trailing_steps=trailing_steps,
        trailing_ratio=None if trailing_ratio is None else str(trailing_ratio),
        maximum_complete_block_ratio=None if maximum is None else str(maximum),
    )


def analyze_census_file(path: Path) -> list[BlockContractionEvidence]:
    """Analyze the strictest-tolerance record for each cutoff in a JSONL census."""

    records: list[Mapping[str, object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        decoded = json.loads(line)
        if not isinstance(decoded, dict):
            raise ValueError("each census line must be a JSON object")
        records.append(decoded)
    strictest: dict[Decimal, Mapping[str, object]] = {}
    for record in records:
        try:
            cutoff = Decimal(str(record["cutoff"]))
            tolerance = Decimal(str(record["tolerance"]))
        except (KeyError, ArithmeticError) as error:
            raise ValueError("census record has invalid cutoff or tolerance") from error
        current = strictest.get(cutoff)
        if current is None or tolerance < Decimal(str(current["tolerance"])):
            strictest[cutoff] = record
    return [analyze_census_record(strictest[key]) for key in sorted(strictest)]
