"""Optional TRIQS MaxEnt adapter for a mass-constrained hidden image.

The provider is imported only when the factory is called.  The adapter uses
TRIQS's general cost function and leaves the choice of minimizer to its caller.
"""

from __future__ import annotations

import copy
import importlib
from collections.abc import Callable
from typing import Any, TypeVar, cast

import numpy as np
from numpy.typing import NDArray

from kernelgecp.integrations._maxent_mass import NormalizedMassMap
from kernelgecp.types import NumericalBreakdownError

FloatArray = NDArray[np.float64]
_F = TypeVar("_F", bound=Callable[..., Any])


def make_native_mass_cost(mass: float) -> tuple[Any, Any]:
    """Create a native cost function and its normalized H(v) component.

    TRIQS reduces the kernel's singular space immediately before calling inv;
    refreshing there ensures the map uses that reduced space.  The state holder
    is shared by the shallow copies created by the provider's cached interface.
    """
    functions = importlib.import_module("triqs_maxent.functions")
    costs = importlib.import_module("triqs_maxent.cost_functions")
    generic_h: Any = functions.GenericH_of_v

    def cached(func: _F) -> _F:
        return cast(_F, functions.cached(func))

    class NativeNormalizedH(generic_h):  # type: ignore[misc]
        def __init__(self) -> None:
            self._mass_state: dict[str, Any] = {"map": None}
            super().__init__()

        def parameter_change(self) -> None:
            self._mass_state["map"] = None

        def _map(self, *, refresh: bool = False) -> NormalizedMassMap:
            current = self._mass_state["map"]
            if refresh or current is None:
                current = NormalizedMassMap(
                    np.asarray(self.D.D, dtype=np.float64),
                    np.asarray(self.K.V, dtype=np.float64),
                    mass,
                )
                self._mass_state["map"] = current
            return current  # type: ignore[no-any-return]

        @property
        def input_size(self) -> tuple[int]:
            return (self._map().rank,)

        @property
        def diagnostics(self) -> dict[str, Any]:
            return copy.deepcopy(dict(self._map().diagnostics))

        @cached
        def f(self, v: FloatArray) -> FloatArray:
            # Native NormalEntropy writes into its argument, including empty
            # floor selections. Its private value must be writable and owned;
            # the immutable mathematical map remains isolated from that state.
            mapping = self._map()
            image = mapping.f(v)
            # NormalEntropy floors H, D and H/D at 1e-100. Reject that
            # numerical regime before its in-place changes can invalidate the
            # chain-rule derivatives or alter the caller's default model.
            with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                ratio = image / mapping.default_weights
            if (
                np.any(image <= 1e-100)
                or np.any(mapping.default_weights <= 1e-100)
                or np.any(ratio <= 1e-100)
                or not np.isfinite(ratio).all()
            ):
                raise NumericalBreakdownError(
                    "constrained image reaches native NormalEntropy floor or overflow"
                )
            return image.copy()

        @cached
        def d(self, v: FloatArray) -> FloatArray:
            return self._map().d(v)

        @cached
        def dd(self, v: FloatArray) -> FloatArray:
            return self._map().dd(v)

        @cached
        def inv(self, hidden: FloatArray) -> FloatArray:
            # LevenbergMinimizer owns and updates this vector with v -= dv.
            return self._map(refresh=True).inv(hidden).copy()

    h = NativeNormalizedH()
    return costs.MaxEntCostFunction(d_dv=True, H_of_v=h), h
