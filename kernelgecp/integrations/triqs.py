"""Owning single-Gf interoperability with the pinned TRIQS 4.0.2 runtime.

TRIQS is imported only when an adapter is called. Green functions retain the
physical minus sign and exp(+i*wn*tau) Fourier convention of FermionicBasis.
No symmetry, conjugation, moment repair or static-term subtraction is applied.
BlockGf callers adapt each block in their original order, then construct a new
BlockGf from the ordered names and output blocks.
"""

from __future__ import annotations

import importlib
import operator
from typing import Any, cast

import numpy as np
import numpy.typing as npt

from ..basis import FermionicBasis, FermionicFit
from ..types import FloatArray

__all__ = ["fit_triqs_tau", "to_triqs_tau", "to_triqs_iw"]


def _provider() -> tuple[Any, Any]:
    try:
        version = importlib.import_module("triqs.version")
        gfs = importlib.import_module("triqs.gfs")
        meshes = importlib.import_module("triqs.mesh")
    except ImportError as error:
        raise ImportError(
            "TRIQS 4.0.2 is required for this optional integration; "
            "use the isolated environments/triqs.yml runtime"
        ) from error
    if getattr(version, "version", None) != "4.0.2":
        raise ImportError("this integration supports the pinned TRIQS 4.0.2 runtime")
    return gfs, meshes


def _mesh(mesh: object, beta: float, mesh_type: Any, name: str) -> Any:
    if not isinstance(mesh, mesh_type):
        raise TypeError(f"{name} requires an actual TRIQS {mesh_type.__name__} mesh")
    value = cast(Any, mesh)
    actual_beta = float(value.beta)
    if not np.isfinite(actual_beta) or actual_beta <= 0 or actual_beta != beta:
        raise ValueError("TRIQS mesh beta must equal the basis inverse temperature")
    if value.statistic != "Fermion":
        raise ValueError("TRIQS mesh statistics must be Fermion")
    if len(value) == 0:
        raise ValueError("TRIQS mesh must contain observations")
    return value


def _tau_coordinates(mesh: object, beta: float, meshes: Any) -> FloatArray:
    value = _mesh(mesh, beta, meshes.MeshImTime, "time conversion")
    coordinates = np.asarray([point.value for point in value], dtype=np.float64)
    if (
        coordinates.ndim != 1
        or coordinates.size != len(value)
        or coordinates.size < 2
        or not np.all(np.isfinite(coordinates))
        or not np.all(np.diff(coordinates) > 0)
        or coordinates[0] != 0.0
        or coordinates[-1] != beta
    ):
        raise ValueError("TRIQS time mesh must be ordered from 0+ through beta-")
    return coordinates


def _iw_indices(mesh: object, beta: float, meshes: Any) -> npt.NDArray[np.int64]:
    value = _mesh(mesh, beta, meshes.MeshImFreq, "frequency conversion")
    points = list(value)
    try:
        raw = [operator.index(point.index) for point in points]
        if any(isinstance(point.index, (bool, np.bool_)) for point in points):
            raise ValueError("boolean frequency indices are invalid")
        indices = np.asarray(raw, dtype=np.int64)
        frequencies = np.asarray([complex(point.value) for point in points])
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("TRIQS mesh must expose ordinary integer indices") from error
    expected = (
        1j
        * np.asarray([2 * index + 1 for index in raw], dtype=np.float64)
        * np.pi
        / beta
    )
    slack = 32 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(expected))))
    if (
        indices.shape != (len(value),)
        or any(right != left + 1 for left, right in zip(raw[:-1], raw[1:], strict=True))
        or indices[0] != operator.index(value.first_index)
        or indices[-1] != operator.index(value.last_index)
        or (bool(value.positive_only) and indices[0] != 0)
        or not np.all(np.isfinite(frequencies))
        or not np.all(np.isfinite(expected))
        or float(np.max(np.abs(frequencies - expected))) > slack
    ):
        raise ValueError("TRIQS frequency coordinates must match ordered fermionic n")
    return indices


def _output(fit: FermionicFit, mesh: object, values: npt.ArrayLike, gfs: Any) -> Any:
    array = np.moveaxis(np.asarray(values), fit.axis, 0)
    target_shape = array.shape[1:]
    if any(size == 0 for size in target_shape) or not np.all(np.isfinite(array)):
        raise ValueError("TRIQS output must have finite nonempty target axes")
    # MeshImTime and MeshImFreq are mutable through copy_from; Gf retains the
    # exact mesh object supplied to its constructor. Own the mesh as well as
    # the data so later caller/output mesh changes cannot reach the other.
    result = gfs.Gf(mesh=cast(Any, mesh).copy(), target_shape=list(target_shape))
    if result.data.shape != array.shape:
        raise ValueError("TRIQS output layout does not match the evaluated fit")
    # Allocate through Gf first, then assign. Never hand it a NumPy buffer that
    # could remain aliased to the caller's fit, input Gf or cached design.
    result.data[...] = array
    return result


def fit_triqs_tau(
    basis: FermionicBasis, gf: object, *, rcond: float | None = None
) -> FermionicFit:
    """Fit one actual fermionic imaginary-time Gf using its physical mesh.

    Scalar target_shape=() and arbitrary trailing target axes are preserved.
    For BlockGf, iterate blocks in caller order and invoke this function for
    each block. Input data are copied before fitting and never conjugated.
    """
    if not isinstance(basis, FermionicBasis):
        raise TypeError("basis must be a FermionicBasis")
    gfs, meshes = _provider()
    if not isinstance(gf, gfs.Gf):
        raise TypeError("fit_triqs_tau accepts one Gf; adapt BlockGf blocks in order")
    source = cast(Any, gf)
    tau = _tau_coordinates(source.mesh, basis.beta, meshes)
    values = np.array(source.data, copy=True)
    target_shape = tuple(source.target_shape)
    if (
        values.shape != (tau.size, *target_shape)
        or values.dtype.kind not in "fc"
        or any(size == 0 for size in target_shape)
        or not np.all(np.isfinite(values))
    ):
        raise ValueError("TRIQS Gf data must match its mesh and finite target shape")
    return basis.fit_tau(values, tau, rcond=rcond)


def to_triqs_tau(fit: FermionicFit, mesh: object) -> Any:
    """Evaluate into a new owning Gf on an actual 0+-through-beta- mesh."""
    if not isinstance(fit, FermionicFit):
        raise TypeError("fit must be a FermionicFit")
    gfs, meshes = _provider()
    tau = _tau_coordinates(mesh, fit.basis.beta, meshes)
    return _output(fit, mesh, fit.evaluate_tau(tau), gfs)


def to_triqs_iw(fit: FermionicFit, mesh: object) -> Any:
    """Evaluate into an owning Gf in the supplied mesh's actual storage order.

    MeshPoint.index is ordinary signed n, not the reduced odd integer 2*n+1.
    Both full signed and positive-only meshes are supported. No values are
    synthesized by conjugating data from the other frequency sign.
    """
    if not isinstance(fit, FermionicFit):
        raise TypeError("fit must be a FermionicFit")
    gfs, meshes = _provider()
    indices = _iw_indices(mesh, fit.basis.beta, meshes)
    return _output(fit, mesh, fit.evaluate_matsubara(indices), gfs)
