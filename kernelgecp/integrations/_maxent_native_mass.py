"""Optional TRIQS MaxEnt adapter for a mass-constrained hidden image.

The provider is imported only when the factory is called.  The adapter uses
TRIQS's general cost function and leaves the choice of minimizer to its caller.
"""

from __future__ import annotations

import copy
import importlib
import math
from collections.abc import Callable
from numbers import Real
from typing import Any, TypeVar, cast

import numpy as np
from numpy.typing import NDArray

from kernelgecp.integrations._maxent_curvature import (
    shifted_projected_metric,
    validate_floor_relative,
)
from kernelgecp.integrations._maxent_mass import NormalizedMassMap
from kernelgecp.types import NumericalBreakdownError

FloatArray = NDArray[np.float64]
_F = TypeVar("_F", bound=Callable[..., Any])


def make_native_mass_cost(
    mass: float,
    *,
    curvature: str = "exact",
    interior_fraction: float = 0.0,
    curvature_floor_relative: float = 0.0,
) -> tuple[Any, Any]:
    """Create a native cost function and its normalized H(v) component.

    TRIQS reduces the kernel's singular space immediately before calling inv;
    refreshing there ensures the map uses that reduced space.  The state holder
    is shared by the shallow copies created by the provider's cached interface.
    """
    if curvature not in ("exact", "projected"):
        raise ValueError("curvature must be exact or projected")
    if (
        isinstance(interior_fraction, bool)
        or not isinstance(interior_fraction, Real)
        or not math.isfinite(interior_fraction)
        or not 0 <= interior_fraction < 1
    ):
        raise ValueError("interior_fraction must be finite and in [0, 1)")
    if interior_fraction > 0 and 1.0 - interior_fraction == 1.0:
        raise ValueError("interior_fraction is not representable in mixture scale")
    floor = validate_floor_relative(curvature_floor_relative)
    if floor and curvature != "projected":
        raise ValueError("mass_curvature_floor_relative requires projected curvature")
    functions = importlib.import_module("triqs_maxent.functions")
    costs = importlib.import_module("triqs_maxent.cost_functions")
    generic_h: Any = functions.GenericH_of_v

    def cached(func: _F) -> _F:
        return cast(_F, functions.cached(func))

    class NativeNormalizedH(generic_h):  # type: ignore[misc]
        def __init__(self) -> None:
            self._mass_state: dict[str, Any] = {
                "map": None,
                "curvature_metric": {
                    "evaluation_count": 0,
                    "last_evaluation": None,
                    "certified": False,
                },
                "last_evaluation": None,
                "max_latent_zero_count": 0,
                "evaluation_count": 0,
                "min_observed_H": None,
                "min_observed_mixture_floor": None,
                "min_observed_q": None,
            }
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
                    interior_fraction=interior_fraction,
                )
                self._mass_state["map"] = current
            return current

        @property
        def input_size(self) -> tuple[int]:
            return (self._map().rank,)

        @property
        def diagnostics(self) -> dict[str, Any]:
            result = copy.deepcopy(dict(self._map().diagnostics))
            result.update(
                curvature=curvature,
                d_dv=curvature == "exact",
                dA_projection=2,
                chi2_factor=1.0,
                floor_policy="reject_native_entropy_floor_or_overflow",
                last_evaluation=copy.deepcopy(self._mass_state["last_evaluation"]),
                max_latent_zero_count=self._mass_state["max_latent_zero_count"],
                evaluation_count=self._mass_state["evaluation_count"],
                min_observed_H=self._mass_state["min_observed_H"],
                min_observed_mixture_floor=self._mass_state[
                    "min_observed_mixture_floor"
                ],
                min_observed_q=self._mass_state["min_observed_q"],
            )
            if floor:
                result["curvature_metric"] = copy.deepcopy(
                    self._mass_state["curvature_metric"]
                )
            return result

        @cached
        def f(self, v: FloatArray) -> FloatArray:
            # Native NormalEntropy writes into its argument, including empty
            # floor selections. Its private value must be writable and owned;
            # the immutable mathematical map remains isolated from that state.
            mapping = self._map()
            image = mapping.f(v)
            observation = dict(mapping.evaluation_diagnostics(v))
            observation["min_q"] = float(
                np.min(mapping.probabilities(np.zeros(mapping.rank)))
            )
            for key, value in (
                ("min_observed_H", observation["min_H"]),
                ("min_observed_mixture_floor", observation["actual_mixture_floor"]),
                ("min_observed_q", observation["min_q"]),
            ):
                previous = self._mass_state[key]
                self._mass_state[key] = (
                    value if previous is None else min(previous, value)
                )
            self._mass_state["last_evaluation"] = observation
            self._mass_state["evaluation_count"] += 1
            self._mass_state["max_latent_zero_count"] = max(
                self._mass_state["max_latent_zero_count"],
                observation["latent_zero_count"],
            )
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
    base_cost: Any = costs.MaxEntCostFunction
    cost_class = base_cost
    if floor:

        class ShiftedProjectedCost(base_cost):  # type: ignore[misc]
            @cached
            def dd(self, v: FloatArray) -> FloatArray:
                metric, info = shifted_projected_metric(super().dd(v), floor)
                state = h._mass_state["curvature_metric"]
                state["evaluation_count"] += 1
                state["last_evaluation"] = info
                return metric

        cost_class = ShiftedProjectedCost
    return cost_class(d_dv=curvature == "exact", dA_projection=2, H_of_v=h), h
