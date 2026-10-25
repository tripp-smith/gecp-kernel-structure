"""Cached numerical operations for a frozen physical fermionic basis.

Working precision applies to kernel evaluation and the selected-core solve.
Coordinates and returned arrays are binary64. Conditioning below describes
the exported binary64 core, including when the factorization uses mpmath.
"""

from __future__ import annotations

import math
import warnings
from numbers import Real
from typing import Any

import mpmath as mp
import numpy as np
import numpy.typing as npt
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve

from ._reference_backend import ReferenceBackend, RustReferenceBackend
from .kernels import FermionicKernel
from .types import FloatArray, NumericalBreakdownError

ComplexArray = npt.NDArray[np.complex128]


def _real_vector(values: npt.ArrayLike, name: str) -> FloatArray:
    raw = np.asarray(values)
    if raw.ndim != 1 or raw.dtype.kind not in "fiu":
        raise ValueError(f"{name} must be a one-dimensional real numeric array")
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.asarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values")
    return result


def _immutable_copy(values: FloatArray) -> FloatArray:
    # Bytes provide an owned copy whose write flag cannot be re-enabled.
    return np.frombuffer(values.tobytes(), dtype=np.float64)


class BasisNumerics:
    """Own a selected core and one cached transpose LU factorization.

    ``omega_nodes`` are dimensionless frequencies. Physical time coordinates
    are converted by binary64 division by beta both when freezing the core
    and when evaluating. Reference arithmetic uses an independent mpmath
    context; no caller precision setting is changed.
    """

    def __init__(
        self,
        beta: float,
        cutoff: float,
        tau_nodes: npt.ArrayLike,
        omega_nodes: npt.ArrayLike,
        *,
        precision_bits: int,
        reference_backend: ReferenceBackend = "python",
    ) -> None:
        for name, value in (("beta", beta), ("cutoff", cutoff)):
            if (
                isinstance(value, (bool, np.bool_))
                or not isinstance(value, Real)
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ValueError(f"{name} must be finite and positive")
        if cutoff < 1:
            raise ValueError("cutoff must be at least one")
        if (
            isinstance(precision_bits, (bool, np.bool_))
            or not isinstance(precision_bits, (int, np.integer))
            or precision_bits < 53
        ):
            raise ValueError("precision_bits must be an integer at least 53")
        self.beta = float(beta)
        self.cutoff = float(cutoff)
        self.precision_bits = int(precision_bits)
        if reference_backend not in ("python", "rust"):
            raise ValueError("reference_backend must be 'python' or 'rust'")
        if reference_backend == "rust" and self.precision_bits not in (53, 128):
            raise ValueError("Rust reference backend requires exactly 128 working bits")
        self._rust_backend: RustReferenceBackend | None = None
        self.native_backend_info: tuple[tuple[str, str | int | bool], ...] = ()
        self.time_evaluation_backend = (
            "scipy_lapack_float64"
            if self.precision_bits == 53
            else "python_mpmath_reference"
        )
        self.frequency_evaluation_backend = self.time_evaluation_backend
        times = _real_vector(tau_nodes, "tau_nodes")
        frequencies = _real_vector(omega_nodes, "omega_nodes")
        if times.size == 0 or times.shape != frequencies.shape:
            raise ValueError("selected nodes must have matching nonempty shapes")
        if np.any(times < 0) or np.any(times > self.beta):
            raise ValueError("tau_nodes must lie in [0, beta]")
        if np.any(np.abs(frequencies) > self.cutoff):
            raise ValueError("omega_nodes must lie in [-cutoff, cutoff]")
        if (
            np.unique(times).size != times.size
            or np.unique(frequencies).size != frequencies.size
        ):
            raise NumericalBreakdownError(
                "selected basis core LU failed: repeated nodes"
            )
        self.tau_nodes = _immutable_copy(times)
        self.omega_nodes = _immutable_copy(frequencies)
        self.rank = int(times.size)
        self._kernel = FermionicKernel(self.cutoff)
        self._context: Any = None
        self._float_lu: Any = None
        self._mp_lu: Any = None
        self._mp_omega_nodes: list[Any] = []
        with np.errstate(under="ignore"):
            normalized_times = self.tau_nodes / self.beta
        if self.precision_bits == 53:
            with np.errstate(under="ignore"):
                core = self._kernel(
                    normalized_times[:, None], self.omega_nodes[None, :]
                )
            if (
                np.unique(core, axis=0).shape[0] != self.rank
                or np.unique(core, axis=1).shape[1] != self.rank
            ):
                raise NumericalBreakdownError(
                    "selected basis core LU failed: repeated rounded rows or columns"
                )
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("error", LinAlgWarning)
                    self._float_lu = lu_factor(core.T)
            except (LinAlgWarning, ValueError, np.linalg.LinAlgError) as exc:
                raise NumericalBreakdownError("selected basis core LU failed") from exc
        else:
            self._context = mp.mp.clone()
            self._context.prec = self.precision_bits
            self._mp_omega_nodes = [
                self._context.mpf(float(value)) for value in self.omega_nodes
            ]
            reference_core = self._context.matrix(
                self._reference_kernel_rows(normalized_times)
            )
            try:
                permutation, lower, upper = self._context.lu(reference_core.T)
            except (ZeroDivisionError, ValueError) as exc:
                raise NumericalBreakdownError("selected basis core LU failed") from exc
            ordering = [
                next(j for j in range(self.rank) if permutation[i, j])
                for i in range(self.rank)
            ]
            self._mp_lu = ordering, lower, upper
            core = np.asarray(reference_core.tolist(), dtype=np.float64)
            if reference_backend == "rust":
                self._rust_backend = RustReferenceBackend(
                    self.omega_nodes, ordering, lower, upper, self._context
                )
                self.time_evaluation_backend = "rust_mpfr_128"
                self.native_backend_info = self._rust_backend.info
        self.exported_core_condition = self._condition(core)
        # Kept as the requested scalar interface; its explicit scope is above.
        self.core_condition = self.exported_core_condition
        self.nodal_error = float(
            np.max(np.abs(self.design_tau(self.tau_nodes) - np.eye(self.rank)))
        )

    @staticmethod
    def _condition(core: FloatArray) -> float | None:
        try:
            singular = np.linalg.svd(core, compute_uv=False)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(singular)) or singular[-1] <= 0:
            return None
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            condition = float(singular[0] / singular[-1])
        return condition if math.isfinite(condition) else None

    def _reference_kernel_rows(self, times: FloatArray) -> list[list[Any]]:
        context = self._context
        rows = []
        for value in times:
            time = context.mpf(float(value))
            rows.append(
                [
                    context.exp(-time * omega) / (1 + context.exp(-omega))
                    if omega >= 0
                    else context.exp((1 - time) * omega) / (1 + context.exp(omega))
                    for omega in self._mp_omega_nodes
                ]
            )
        return rows

    def _solve_reference(self, rows: list[list[Any]]) -> list[list[Any]]:
        """Batched substitution using the previously computed transpose LU."""
        ordering, lower, upper = self._mp_lu
        context = self._context
        count = len(rows)
        forward = [[context.zero for _ in rows] for _ in range(self.rank)]
        solution = [[context.zero for _ in rows] for _ in range(self.rank)]
        for i in range(self.rank):
            for column in range(count):
                forward[i][column] = (
                    rows[column][ordering[i]]
                    - context.fsum(lower[i, j] * forward[j][column] for j in range(i))
                ) / lower[i, i]
        for i in range(self.rank - 1, -1, -1):
            for column in range(count):
                solution[i][column] = (
                    forward[i][column]
                    - context.fsum(
                        upper[i, j] * solution[j][column]
                        for j in range(i + 1, self.rank)
                    )
                ) / upper[i, i]
        return [
            [solution[i][column] for i in range(self.rank)] for column in range(count)
        ]

    @staticmethod
    def _check_output(values: FloatArray | ComplexArray) -> None:
        if not np.all(np.isfinite(values)):
            raise NumericalBreakdownError(
                "basis solve produced nonfinite binary64 values"
            )

    def design_tau(self, tau: npt.ArrayLike) -> FloatArray:
        """Evaluate cardinal columns at physical times, exporting binary64."""
        times = _real_vector(tau, "tau")
        if np.any(times < 0) or np.any(times > self.beta):
            raise ValueError("tau must lie in [0, beta]")
        if times.size == 0:
            return np.empty((0, self.rank), dtype=np.float64)
        with np.errstate(under="ignore"):
            normalized_times = times / self.beta
            if self.precision_bits == 53:
                columns = self._kernel(
                    normalized_times[:, None], self.omega_nodes[None, :]
                )
                result = np.asarray(
                    lu_solve(self._float_lu, columns.T).T, dtype=np.float64
                )
            elif self._rust_backend is not None:
                result = self._rust_backend.design_tau(normalized_times)
            else:
                result = np.asarray(
                    self._solve_reference(
                        self._reference_kernel_rows(normalized_times)
                    ),
                    dtype=np.float64,
                )
        self._check_output(result)
        return result

    def design_matsubara(self, indices: npt.ArrayLike) -> ComplexArray:
        """Integrate cardinal columns using ordinary fermionic indices n.

        For z = i*pi*(2*n+1), the positive kernel column integrates to
        -beta/(z-Omega). Coefficients fitted to physical negative G(tau)
        already contain their sign; no additional sign is applied here.
        """
        modes = np.asarray(indices)
        if modes.ndim != 1 or modes.dtype.kind not in "iu":
            raise ValueError("indices must be a one-dimensional integer array")
        if modes.size == 0:
            return np.empty((0, self.rank), dtype=np.complex128)
        with np.errstate(under="ignore"):
            if self.precision_bits == 53:
                odd = np.asarray([2 * int(n) + 1 for n in modes], dtype=np.float64)
                points = 1j * np.pi * odd
                columns = -self.beta / (points[:, None] - self.omega_nodes[None, :])
                result = np.asarray(
                    lu_solve(self._float_lu, columns.T).T, dtype=np.complex128
                )
            else:
                context = self._context
                beta = context.mpf(self.beta)
                rows = [
                    [
                        -beta / (context.j * context.pi * (2 * int(n) + 1) - omega)
                        for omega in self._mp_omega_nodes
                    ]
                    for n in modes
                ]
                result = np.asarray(self._solve_reference(rows), dtype=np.complex128)
        self._check_output(result)
        return result
