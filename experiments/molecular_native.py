"""Native Green IR reconstruction independent of kernelgecp compression.

The saved Green grid supplies its actual physical basis normalization. Native
interior time values are solved against that basis once. Arbitrary queries use
the independently shipped sparse-ir 2.1.6 functions, with an explicit consistency
check against the native grid tables. No Green inverse matrix is copied.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from kernelgecp.integrations.green import GreenTauData


class NativeIRSolution:
    def __init__(self, data: GreenTauData, grid_path: Path, native_lambda: float):
        import h5py
        import sparse_ir

        self.data = data
        self.native_lambda = native_lambda
        self.basis = sparse_ir.FiniteTempBasis("F", beta=1.0, wmax=native_lambda)
        with h5py.File(grid_path, "r") as source:
            if source.attrs.get("__grids_version__") != data.provider.grids:
                raise ValueError("native grid and output writer versions differ")
            if float(source["fermi/metadata/lambda"][()]) != native_lambda:
                raise ValueError("native grid cutoff metadata differs from the case")
            self.indices = np.asarray(source["fermi/ngrid"], dtype=np.int64)
            x = np.asarray(source["fermi/xgrid"], dtype=float)
            interior = np.asarray(source["fermi/uxl"], dtype=float)
            frequency = np.asarray(source["fermi/uwl"], dtype=complex)
            edges = [
                np.asarray(source["fermi/" + name], dtype=float)
                for name in ("u1l_neg", "u1l_pos")
            ]
        self.dimension = interior.shape[1]
        if (
            not all(
                np.all(np.isfinite(value)) for value in (x, interior, frequency, *edges)
            )
            or self.indices.ndim != 1
            or len(np.unique(self.indices)) != len(self.indices)
            or np.any(self.indices < -(2**62))
            or np.any(self.indices >= 2**62)
        ):
            raise ValueError(
                "native grid tables must be finite with unique safe integer indices"
            )
        self.time_table = np.vstack((edges[0], interior, edges[1])) * np.sqrt(
            2 / data.beta
        )
        self.frequency_table = frequency * np.sqrt(data.beta)
        tau = np.r_[0.0, (x + 1) * data.beta / 2, data.beta]
        if data.tau.shape != tau.shape or np.max(abs(data.tau - tau)) > 1e-12:
            raise ValueError(
                "native output time mesh does not match the retained IR grid"
            )
        if interior.shape[0] != interior.shape[1] or self.basis.size < self.dimension:
            raise ValueError(
                "native grid exceeds the independently installed sparse-ir basis"
            )
        # Native tables are defined at xgrid, before their physical-tau
        # conversion. At large Lambda, converting x->tau->tau/beta can move a
        # near-endpoint argument by one ulp and magnify a basis column's steep
        # derivative. Check normalization at the original defining arguments;
        # record the round-trip matrix effect separately from field errors.
        defining_times = np.r_[0.0, (x + 1) / 2, 1.0]
        independent_time = np.asarray(
            self.basis.u[: self.dimension](defining_times)
        ).T / np.sqrt(data.beta)
        independent_frequency = self._frequency_design(self.indices)
        self.table_errors = {
            "tau": float(np.max(abs(independent_time - self.time_table))),
            "frequency": float(
                np.max(abs(independent_frequency - self.frequency_table))
            ),
        }
        if max(self.table_errors.values()) > 2e-10:
            raise ValueError(
                "independent sparse-ir does not reproduce the native grid normalization"
            )
        self.condition_number = float(np.linalg.cond(self.time_table[1:-1]))
        self.physical_tau_roundtrip_matrix_error = float(
            np.max(abs(self._time_design(data.tau) - self.time_table))
        )
        factor = lu_factor(self.time_table[1:-1])
        self.coefficients = {
            "green": lu_solve(factor, data.g_tau[1:-1].reshape(self.dimension, -1)),
            "sigma": lu_solve(factor, data.sigma_tau[1:-1].reshape(self.dimension, -1)),
        }
        if not all(np.all(np.isfinite(value)) for value in self.coefficients.values()):
            raise ValueError("native IR reconstruction is nonfinite")
        self.endpoint_reconstruction_errors = {
            name: float(
                np.max(
                    abs(
                        (self.time_table @ coefficients).reshape(data.g_tau.shape)
                        - getattr(data, "g_tau" if name == "green" else "sigma_tau")
                    )
                )
            )
            for name, coefficients in self.coefficients.items()
        }
        self.independent_field_reconstruction_errors = {
            field: float(
                np.max(
                    abs(
                        self.evaluate_tau(data.tau, field=field)
                        - getattr(data, "g_tau" if field == "green" else "sigma_tau")
                    )
                )
            )
            for field in self.coefficients
        }
        if max(self.independent_field_reconstruction_errors.values()) > 2e-9:
            raise ValueError("independent native IR fields fail reconstruction")

    def _time_design(self, tau: Any) -> np.ndarray:
        points = np.asarray(tau, dtype=float)
        if (
            points.ndim != 1
            or not np.all(np.isfinite(points))
            or np.any(points < 0)
            or np.any(points > self.data.beta)
        ):
            raise ValueError("native time queries must be finite in [0,beta]")
        return np.asarray(
            self.basis.u[: self.dimension](points / self.data.beta)
        ).T / np.sqrt(self.data.beta)

    def _frequency_design(self, indices: Any) -> np.ndarray:
        points = np.asarray(indices)
        if points.ndim != 1 or points.dtype.kind not in "iu":
            raise ValueError(
                "native frequency queries must be ordinary integer indices"
            )
        return np.asarray(
            self.basis.uhat[: self.dimension](2 * points + 1), dtype=complex
        ).T * np.sqrt(self.data.beta)

    def evaluate_tau(self, tau: Any, *, field: str = "green") -> np.ndarray:
        if field not in self.coefficients:
            raise ValueError("field must be green or sigma")
        design = self._time_design(tau)
        return (design @ self.coefficients[field]).reshape(
            (len(design), *self.data.g_tau.shape[1:])
        )

    def evaluate_matsubara(self, indices: Any, *, field: str = "green") -> np.ndarray:
        if field not in self.coefficients:
            raise ValueError("field must be green or sigma")
        design = self._frequency_design(indices)
        return (design @ self.coefficients[field]).reshape(
            (len(design), *self.data.g_tau.shape[1:])
        )

    def dyson_residual(self, indices: Any) -> float:
        points = np.asarray(indices, dtype=np.int64)
        z = 1j * (2 * points + 1) * np.pi / self.data.beta
        green = self.evaluate_matsubara(points)
        sigma = self.evaluate_matsubara(points, field="sigma")
        operator = (
            (z + self.data.mu)[:, None, None, None, None] * self.data.overlap
            - self.data.h_core
            - self.data.sigma_static
            - sigma
        )
        residual = operator @ green - np.eye(green.shape[-1])
        return float(np.max(abs(residual)))
