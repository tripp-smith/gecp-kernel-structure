import numpy as np
import pytest

from kernelgecp.integrations._maxent_mass import NormalizedMassMap
from kernelgecp.types import NumericalBreakdownError


def example():
    x = np.linspace(-1, 1, 9)
    return NormalizedMassMap(np.exp(x), np.column_stack((x, x**2, np.ones(9))), 2.3)


def test_derivatives_against_central_differences_and_mass():
    model = example()
    v = np.array([0.3, -0.4])
    eps = 1e-5
    jac = model.d(v)
    hess = model.dd(v)
    for j in range(model.rank):
        step = np.eye(model.rank)[j] * eps
        np.testing.assert_allclose(
            (model.f(v + step) - model.f(v - step)) / (2 * eps),
            jac[:, j],
            rtol=2e-8,
            atol=1e-10,
        )
        np.testing.assert_allclose(
            (model.d(v + step) - model.d(v - step)) / (2 * eps),
            hess[:, :, j],
            rtol=2e-8,
            atol=1e-10,
        )
    assert np.isclose(model.f(v).sum(), model.mass)
    np.testing.assert_allclose(jac.sum(axis=0), 0, atol=1e-15)
    np.testing.assert_allclose(hess.sum(axis=0), 0, atol=1e-15)


def test_roundtrip_default_gauge_and_ownership():
    d = np.array([1.0, 2.0, 3.0, 4.0])
    f = np.column_stack(([0.0, 1.0, 2.0, 3.0], [1.0, -1.0, -1.0, 1.0], np.ones(4)))
    model = NormalizedMassMap(d, f, 1.0)
    assert model.rank == model.input_size == 2
    np.testing.assert_allclose(model.f(np.zeros(2)), d / d.sum())
    np.testing.assert_allclose(model.inv(d / d.sum()), np.zeros(2), atol=1e-14)
    v = np.array([0.7, -0.9])
    np.testing.assert_allclose(model.inv(model.f(v)), v, atol=1e-12)
    d[:] = 100
    f[:] = 0
    np.testing.assert_allclose(model.f(np.zeros(2)), [0.1, 0.2, 0.3, 0.4])
    for array in (
        model.basis,
        model.default_weights,
        model.singular_values,
        model.f(v),
        model.d(v),
        model.dd(v),
        model.inv(model.f(v)),
    ):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(TypeError):
        model.diagnostics["effective_rank"] = 1
    with pytest.raises(ValueError, match="outside"):
        model.inv(np.array([0.1, 0.1, 0.1, 0.7]))


@pytest.mark.parametrize(
    "d,f,m",
    [
        ([1.0, 0.0], [[0.0], [1.0]], 1.0),
        ([1.0, np.nan], [[0.0], [1.0]], 1.0),
        ([1.0, 1.0], [[0.0], [np.inf]], 1.0),
        ([1.0, 1.0], [[1.0], [1.0]], 1.0),
        ([1.0, 1.0], [[0.0], [1.0]], 0.0),
        ([1.0, 1.0], [[0.0]], 1.0),
    ],
)
def test_invalid_input(d, f, m):
    with pytest.raises(ValueError):
        NormalizedMassMap(d, f, m)


def test_resource_cap_before_svd(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("SVD must not run")

    monkeypatch.setattr(np.linalg, "svd", forbidden)
    for rows, columns in [(2049, 1), (2, 129), (2048, 128)]:
        with pytest.raises(ValueError, match="cap"):
            NormalizedMassMap(np.ones(rows), np.ones((rows, columns)), 1)


def test_numerical_breakdown_and_bad_inverse():
    model = example()
    with pytest.raises(NumericalBreakdownError):
        model.f(np.full(model.rank, 1e10))
    with pytest.raises(ValueError):
        model.inv(np.zeros(9))
    with pytest.raises(ValueError):
        model.f([np.nan, 0])
    with pytest.raises(NumericalBreakdownError):
        NormalizedMassMap([1e-320, 1e300], [[0.0], [1.0]], 1)


def test_near_constant_feature_is_removed():
    x = np.linspace(-1, 1, 6)
    model = NormalizedMassMap(np.ones(6), np.column_stack((x, 1 + 1e-15 * x)), 1)
    assert model.rank == 1
    assert model.diagnostics["original_rank"] == 2
    assert NormalizedMassMap(np.ones(6), (1 + 1e-15 * x)[:, None], 1).rank == 1


def test_offset_independence_and_exact_constant():
    x = np.array([-1.0, 0.0, 1.0])
    plain = NormalizedMassMap(np.ones(3), x[:, None], 1.0)
    shifted = NormalizedMassMap(np.ones(3), (1e13 + x)[:, None], 1.0)
    h = plain.f([0.7])
    np.testing.assert_allclose(shifted.f(shifted.inv(h)), h, rtol=1e-13)
    with pytest.raises(ValueError, match="zero effective rank"):
        NormalizedMassMap(np.ones(3), np.full((3, 1), 1e13), 1.0)


def test_cap_precedes_coercion(monkeypatch):
    integer_features = np.ones((2049, 1), dtype=np.int64)
    integer_default = np.ones(2049, dtype=np.int64)

    def forbidden(*args, **kwargs):
        raise AssertionError("must reject before coercion")

    monkeypatch.setattr(np, "asarray", forbidden)
    with pytest.raises(ValueError, match="cap"):
        NormalizedMassMap(integer_default, integer_features, 1.0)
    with pytest.raises(ValueError, match="cap"):
        NormalizedMassMap([1] * 2049, [[1]] * 2049, 1.0)


def test_reject_unsupported_types_and_subnormal_mass():
    for weights, features in [
        ([1 + 1j, 1], [[0.0], [1.0]]),
        ([1, 1], np.array([[0j], [1j]])),
        (np.array([1, 1], dtype=object), [[0], [1]]),
    ]:
        with pytest.raises(ValueError, match="real numeric"):
            NormalizedMassMap(weights, features, 1.0)
    for mass in (1e-320, True, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            NormalizedMassMap([1, 1], [[0], [1]], mass)
    with pytest.raises(ValueError, match="real numeric"):
        example().f([1j, 0])
