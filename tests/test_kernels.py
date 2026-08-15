import math

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from kernelgecp import FermionicKernel


@given(
    t=st.floats(min_value=0, max_value=1, allow_nan=False, allow_infinity=False),
    omega=st.floats(
        min_value=-100, max_value=100, allow_nan=False, allow_infinity=False
    ),
)
def test_reflection_identity(t: float, omega: float) -> None:
    kernel = FermionicKernel(100)
    assert kernel(t, omega) == pytest.approx(kernel(1 - t, -omega), rel=2e-14)


def test_boundary_evaluation_is_finite() -> None:
    kernel = FermionicKernel(1_000_000)
    values = kernel(
        np.array([0.0, 1.0, 0.25, 0.75]),
        np.array([1_000_000.0, -1_000_000.0, 1_000_000.0, -1_000_000.0]),
    )
    assert np.all(np.isfinite(values))
    assert values[0] == pytest.approx(1.0)
    assert values[1] == pytest.approx(1.0)


def test_float64_matches_high_precision() -> None:
    kernel = FermionicKernel(1_000)
    for t, omega in ((0.2, -40.0), (0.7, 50.0), (0.5, 0.0)):
        expected = float(kernel.high_precision(t, omega, 256))
        assert kernel(t, omega) == pytest.approx(expected, rel=4e-15, abs=1e-300)


def test_validation() -> None:
    with pytest.raises(ValueError):
        FermionicKernel(math.inf)
    kernel = FermionicKernel(10)
    with pytest.raises(ValueError):
        kernel(1.1, 0.0)
    with pytest.raises(ValueError):
        kernel(0.5, 11.0)
