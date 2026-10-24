"""API routing/failure isolation plus optional real-package smoke checks."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from experiments import dlr_adapters as adapters


def _fake_module(method: str, calls: dict) -> SimpleNamespace:
    class Basis:
        size = 2

        def __init__(self, *args: object, **kwargs: object) -> None:
            calls["construction"] = (args, kwargs)

        def __len__(self) -> int:
            return self.size

        def get_tau(self, beta: float) -> np.ndarray:
            return np.asarray([0.1, 0.9]) * beta

        def u(self, tau: np.ndarray) -> np.ndarray:
            return np.asarray([np.ones(tau.size), tau])

        def eval_dlr_tau(
            self, coeffs: np.ndarray, tau: np.ndarray, beta: float
        ) -> np.ndarray:
            assert coeffs.ndim in (1, 3)
            calls["evaluation_beta"] = beta
            return np.tensordot(self.u(tau).T, coeffs, axes=(1, 0))

        def lstsq_dlr_from_tau(
            self, tau: np.ndarray, values: np.ndarray, beta: float
        ) -> np.ndarray:
            calls["fit_tau"] = tau.copy()
            calls["fit_values"] = values.copy()
            calls["fit_beta"] = beta
            assert values.ndim == 3
            return np.linalg.lstsq(self.u(tau).T, values[:, :, 0], rcond=None)[0][
                :, :, None
            ]

    class Sampling:
        def __init__(
            self,
            basis: Basis,
            sampling_points: np.ndarray | None = None,
            use_positive_taus: bool = False,
        ) -> None:
            self.basis = basis
            self.sampling_points = (
                np.array([0.3, 2.7])
                if sampling_points is None
                else sampling_points.copy()
            )
            if sampling_points is None:
                calls["positive_taus"] = use_positive_taus

        def fit(self, values: np.ndarray, axis: int) -> np.ndarray:
            assert axis == 0
            calls["fit_tau"] = self.sampling_points.copy()
            calls["fit_values"] = values.copy()
            shape = (2, *values.shape[1:])
            return np.linalg.lstsq(
                self.basis.u(self.sampling_points).T,
                values.reshape(values.shape[0], -1),
                rcond=None,
            )[0].reshape(shape)

    if method == "pydlr":
        return SimpleNamespace(dlr=Basis, __version__="fake-pydlr")
    return SimpleNamespace(
        FiniteTempBasis=Basis, TauSampling=Sampling, __version__="fake-sparse-ir"
    )


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_external_adapter_uses_supplied_physical_times_and_rhs(
    method: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: dict = {}
    module = _fake_module(method, calls)
    monkeypatch.setattr(adapters.importlib, "import_module", lambda name: module)
    adapter = adapters.build_adapter(method, beta=3, wmax=2, tol=1e-7)
    args, kwargs = calls["construction"]
    if method == "pydlr":
        assert kwargs == {"lamb": 6.0, "eps": 1e-7, "xi": -1, "python_impl": True}
    else:
        assert args == ("F",)
        assert kwargs == {"beta": 3.0, "wmax": 2.0, "eps": 1e-7}
        assert calls["positive_taus"]
    assert adapter.rank == 2
    np.testing.assert_allclose(adapter.tau_nodes, [0.3, 2.7])
    assert adapter.metadata["lambda"] == 6
    assert adapter.metadata["certification_status"] == "not_certified"
    times = np.array([0.0, 0.4, 1.7, 3.0])
    expected = np.array([[-0.7, 0.2j], [0.1, -0.05j]])
    values = np.column_stack([np.ones(times.size), times]) @ expected
    fitted = adapter.fit(values, times)
    np.testing.assert_allclose(fitted, expected, atol=1e-14)
    np.testing.assert_array_equal(calls["fit_tau"], times)
    fit_values = calls["fit_values"]
    if method == "pydlr":
        fit_values = fit_values[:, :, 0]
        assert calls["fit_beta"] == 3
    np.testing.assert_array_equal(fit_values, values)
    heldout = np.array([0.13, 2.41])
    np.testing.assert_allclose(
        adapter.evaluate(fitted, heldout), adapter.design_matrix(heldout) @ expected
    )
    scalar = adapter.fit(values[:, 0].real, times)
    assert scalar.shape == (2,)
    assert adapter.evaluate(scalar, heldout).shape == (2,)
    tensor = adapter.fit(np.stack([values, values], axis=2), times)
    assert tensor.shape == (2, 2, 2)
    assert adapter.evaluate(tensor, heldout).shape == (2, 2, 2)


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_external_failures_distinguish_import_from_runtime(
    method: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unavailable(name: str) -> None:
        raise ModuleNotFoundError("dependency missing")

    monkeypatch.setattr(adapters.importlib, "import_module", unavailable)
    with pytest.raises(adapters.AdapterUnavailableError) as missing:
        adapters.build_adapter(method, 3, 2, 1e-7)
    assert missing.value.status == "import_unavailable"

    def broken(name: str) -> None:
        raise OSError("shared library could not load")

    monkeypatch.setattr(adapters.importlib, "import_module", broken)
    with pytest.raises(adapters.AdapterRuntimeError) as failed:
        adapters.build_adapter(method, 3, 2, 1e-7)
    assert failed.value.status == "runtime_failed"
    monkeypatch.setattr(
        adapters.importlib, "import_module", lambda name: SimpleNamespace()
    )
    with pytest.raises(adapters.AdapterRuntimeError, match="construction"):
        adapters.build_adapter(method, 3, 2, 1e-7)


def test_sparse_ir_legacy_positive_time_constructor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: dict = {}
    module = _fake_module("sparse_ir", calls)
    modern_sampling = module.TauSampling

    class LegacySampling:
        def __init__(self, basis: object, sampling_points: object = None) -> None:
            self._inner = modern_sampling(basis, sampling_points)
            self.sampling_points = self._inner.sampling_points

        def fit(self, values: np.ndarray, axis: int = 0) -> np.ndarray:
            return self._inner.fit(values, axis)

    module.TauSampling = LegacySampling
    monkeypatch.setattr(adapters.importlib, "import_module", lambda name: module)
    adapter = adapters.build_adapter("sparse_ir", 3, 2, 1e-7)
    np.testing.assert_allclose(adapter.tau_nodes, [0.3, 2.7])
    values = np.array([-0.7, -0.4])
    coeffs = adapter.fit(values, np.array([0.0, 3.0]))
    np.testing.assert_allclose(adapter.evaluate(coeffs, np.array([0.0, 3.0])), values)


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_adapter_validation_and_runtime_output_checks(
    method: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        adapters.importlib, "import_module", lambda name: _fake_module(method, {})
    )
    adapter = adapters.build_adapter(method, 3, 2, 1e-7)
    for tau in ([float("nan")], [-0.01], [3.01], [[0.5]], [0.5j]):
        with pytest.raises(ValueError):
            adapter.design_matrix(np.asarray(tau))
    for tau, values in (([0.2, 0.2], [1, 2]), ([0.2], [1]), ([0, 1], [1])):
        with pytest.raises(ValueError):
            adapter.fit(np.asarray(values), np.asarray(tau))
    with pytest.raises(ValueError):
        adapter.evaluate(np.array([1.0]), np.array([0.5]))
    with pytest.raises(ValueError):
        adapter.fit(np.array([1.0, np.inf]), np.array([0.0, 1.0]))
    # Positive-time endpoint semantics must not depend on the sign bit of zero.
    if method == "pydlr":
        original = adapter._basis.eval_dlr_tau

        def evaluate(coeffs: np.ndarray, tau: np.ndarray, beta: float) -> np.ndarray:
            assert not np.any(np.signbit(tau[tau == 0]))
            return original(coeffs, tau, beta)

        monkeypatch.setattr(adapter._basis, "eval_dlr_tau", evaluate)
    adapter.design_matrix(np.array([-0.0, 3.0]))
    monkeypatch.setattr(adapter._basis, "u", lambda tau: np.full((2, tau.size), np.nan))
    with pytest.raises(adapters.AdapterRuntimeError) as invalid:
        adapter.design_matrix(np.array([0.3]))
    assert invalid.value.status == "runtime_failed"


@pytest.mark.parametrize(
    "method,beta,wmax,tol",
    [
        ("unknown", 3, 2, 1e-7),
        ("pydlr", True, 2, 1e-7),
        ("pydlr", 0, 2, 1e-7),
        ("pydlr", 3, float("inf"), 1e-7),
        ("pydlr", 3, 2, 0),
        ("pydlr", 3, 2, 1),
        ("pydlr", 1e308, 1e308, 1e-7),
        ("pydlr", 1e-300, 1e-300, 1e-7),
    ],
)
def test_invalid_configuration_does_not_import(
    method: str, beta: float, wmax: float, tol: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unexpected(name: str) -> None:
        pytest.fail("invalid configuration attempted an optional import")

    monkeypatch.setattr(adapters.importlib, "import_module", unexpected)
    with pytest.raises(ValueError):
        adapters.build_adapter(method, beta, wmax, tol)


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_actual_optional_adapter_signed_finite_poles(method: str) -> None:
    pytest.importorskip(method)
    beta, wmax = 3.0, 2.0
    adapter = adapters.build_adapter(method, beta, wmax, 1e-9)
    poles = np.array([-1.5, 0.35, 1.2])
    weights = np.array([0.2, 0.5, 0.3])

    def truth(tau: np.ndarray) -> np.ndarray:
        return -np.exp(-tau[:, None] * poles - np.logaddexp(0, -beta * poles)) @ weights

    training = np.linspace(0, beta, 101)
    heldout = (np.arange(87) + 0.31) * beta / 87
    coeffs = adapter.fit(truth(training), training)
    assert np.max(np.abs(adapter.evaluate(coeffs, heldout) - truth(heldout))) < 5e-8
    common = np.linalg.lstsq(
        adapter.design_matrix(training), truth(training), rcond=None
    )[0]
    assert (
        np.max(np.abs(adapter.design_matrix(heldout) @ common - truth(heldout))) < 5e-8
    )
    multiplier = np.array([[1.0, 2.0j], [0.5j, -0.25]])
    tensor_values = truth(training)[:, None, None] * multiplier
    tensor_coeffs = adapter.fit(tensor_values, training)
    tensor_prediction = adapter.evaluate(tensor_coeffs, heldout)
    np.testing.assert_allclose(
        tensor_prediction, truth(heldout)[:, None, None] * multiplier, atol=1e-7
    )
