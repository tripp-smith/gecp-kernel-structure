"""Stable evaluations of kernels used by the project."""

from __future__ import annotations

from dataclasses import dataclass
from typing import overload

import mpmath as mp
import numpy as np
import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class FermionicKernel:
    r"""The dimensionless fermionic kernel on ``[0, 1] × [-cutoff, cutoff]``.

    The mathematical expression is

    ``exp(-t * omega) / (1 + exp(-omega))``.

    Evaluation uses sign-dependent but algebraically equal forms so neither
    numerator nor denominator overflows at large negative frequency.
    """

    cutoff: float
    strict: bool = True

    def __post_init__(self) -> None:
        if not np.isfinite(self.cutoff) or self.cutoff < 1:
            raise ValueError("cutoff must be finite and at least one")

    @overload
    def __call__(self, t: float, omega: float) -> float: ...

    @overload
    def __call__(
        self, t: npt.ArrayLike, omega: npt.ArrayLike
    ) -> npt.NDArray[np.float64]: ...

    def __call__(
        self, t: float | npt.ArrayLike, omega: float | npt.ArrayLike
    ) -> float | npt.NDArray[np.float64]:
        t_array = np.asarray(t, dtype=np.float64)
        omega_array = np.asarray(omega, dtype=np.float64)
        t_broadcast, omega_broadcast = np.broadcast_arrays(t_array, omega_array)
        if self.strict:
            if np.any(~np.isfinite(t_broadcast)) or np.any(
                ~np.isfinite(omega_broadcast)
            ):
                raise ValueError("kernel coordinates must be finite")
            if np.any((t_broadcast < 0) | (t_broadcast > 1)):
                raise ValueError("t must lie in [0, 1]")
            if np.any(np.abs(omega_broadcast) > self.cutoff):
                raise ValueError("omega lies outside the configured cutoff")

        positive = omega_broadcast >= 0
        values = np.empty_like(omega_broadcast, dtype=np.float64)
        values[positive] = np.exp(
            -t_broadcast[positive] * omega_broadcast[positive]
        ) / (1.0 + np.exp(-omega_broadcast[positive]))
        values[~positive] = np.exp(
            (1.0 - t_broadcast[~positive]) * omega_broadcast[~positive]
        ) / (1.0 + np.exp(omega_broadcast[~positive]))
        if values.ndim == 0:
            return float(values)
        return values

    def centered(self, t: float, omega: float) -> float:
        """Evaluate ``exp(-(t-1/2)ω) / (2 cosh(ω/2))`` stably."""

        return float(self(t, omega))

    def high_precision(
        self, t: float, omega: float, precision_bits: int = 128
    ) -> mp.mpf:
        """Evaluate using mpmath at an explicit binary precision."""

        if precision_bits < 53:
            raise ValueError("precision_bits must be at least 53")
        if self.strict and not (0 <= t <= 1 and abs(omega) <= self.cutoff):
            raise ValueError("coordinates lie outside the kernel domain")
        with mp.workprec(precision_bits):
            t_mp = mp.mpf(t)
            omega_mp = mp.mpf(omega)
            if omega >= 0:
                return mp.exp(-t_mp * omega_mp) / (1 + mp.exp(-omega_mp))
            return mp.exp((1 - t_mp) * omega_mp) / (1 + mp.exp(omega_mp))

    def reflected(self, t: float, omega: float) -> float:
        """Evaluate the reflected coordinate ``K(1-t, -omega)``."""

        return float(self(1.0 - t, -omega))
