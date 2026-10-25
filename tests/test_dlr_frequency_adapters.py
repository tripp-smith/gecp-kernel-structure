"""Physical Matsubara conventions, complex linearity, and optional native checks."""

from __future__ import annotations

from types import SimpleNamespace

import mpmath as mp
import numpy as np
import pytest

from experiments import dlr_adapters as adapters


def _fake_frequency_module(method: str, calls: dict) -> SimpleNamespace:
    class Basis:
        size = 2

        def __init__(self, *args: object, **kwargs: object) -> None:
            self.beta = float(kwargs.get("beta", 3.7))
            self.poles = np.array([-0.8, 1.1])

        def __len__(self) -> int:
            return self.size

        def get_tau(self, beta: float) -> np.ndarray:
            return np.array([0.1, 0.9]) * beta

        def u(self, tau: np.ndarray) -> np.ndarray:
            return np.exp(
                -tau[None, :] * self.poles[:, None]
                - np.logaddexp(0, -self.beta * self.poles[:, None])
            )

        def eval_dlr_tau(
            self, coeffs: np.ndarray, tau: np.ndarray, beta: float
        ) -> np.ndarray:
            assert beta == self.beta
            return np.tensordot(self.u(tau).T, coeffs, axes=(1, 0))

        def eval_dlr_freq(
            self, coeffs: np.ndarray, z: np.ndarray, beta: float
        ) -> np.ndarray:
            calls["frequency_coefficients"] = coeffs.copy()
            calls["frequency_z"] = z.copy()
            calls["frequency_beta"] = beta
            result = np.tensordot(
                1 / (z[:, None] + self.poles[None, :]), coeffs, axes=(1, 0)
            )
            if result.ndim == 3:
                result = result.transpose(0, 2, 1)
            return result.conj()

        def uhat(self, reduced: np.ndarray) -> np.ndarray:
            calls["frequency_reduced"] = reduced.copy()
            return 1 / (self.poles[:, None] - 1j * reduced[None, :] * np.pi / self.beta)

    class Sampling:
        def __init__(self, basis: Basis, use_positive_taus: bool = False) -> None:
            self.sampling_points = np.array([0.1, 0.9]) * basis.beta

    if method == "pydlr":
        return SimpleNamespace(dlr=Basis)
    return SimpleNamespace(FiniteTempBasis=Basis, TauSampling=Sampling)


def _fake_adapter(
    method: str, monkeypatch: pytest.MonkeyPatch, calls: dict | None = None
) -> adapters.FermionicAdapter:
    module = _fake_frequency_module(method, {} if calls is None else calls)
    monkeypatch.setattr(adapters.importlib, "import_module", lambda name: module)
    return adapters.build_adapter(method, beta=3.7, wmax=2.4, tol=1e-9)


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_frequency_design_routes_ordinary_indices_and_preserves_complex_batches(
    method: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: dict = {}
    adapter = _fake_adapter(method, monkeypatch, calls)
    indices = np.array([-4, -1, 0, 2], dtype=np.int64)
    frequencies = (2 * indices + 1) * np.pi / adapter.beta
    expected = 1 / (np.array([-0.8, 1.1])[None, :] - 1j * frequencies[:, None])
    design = adapter.design_matsubara(indices)
    np.testing.assert_allclose(design, expected)
    if method == "pydlr":
        np.testing.assert_array_equal(calls["frequency_coefficients"], np.eye(2))
        assert not np.iscomplexobj(calls["frequency_coefficients"])
        np.testing.assert_allclose(calls["frequency_z"], 1j * frequencies)
        assert calls["frequency_beta"] == adapter.beta
    else:
        np.testing.assert_array_equal(calls["frequency_reduced"], 2 * indices + 1)
    coefficients = np.array(
        [[[1 + 2j, -0.3j], [0.1, 0.7j]], [[-0.2j, 3.0], [0.8 + 0.4j, -1.1]]]
    )
    prediction = np.tensordot(design, coefficients, axes=(1, 0))
    for row in range(2):
        for column in range(2):
            np.testing.assert_allclose(
                prediction[:, row, column], expected @ coefficients[:, row, column]
            )
    times = np.array([0.0, 0.7, adapter.beta])
    np.testing.assert_array_equal(
        adapter.design_tau(times), adapter.design_matrix(times)
    )
    assert adapter.metadata["matsubara_index_convention"] == "ordinary_integer_n"
    assert "+i*omega_n*tau" in adapter.metadata["fourier_convention"]


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
@pytest.mark.parametrize(
    "indices",
    [
        0,
        [[0, 1]],
        [True],
        [0, True],
        [np.bool_(False), np.int64(1)],
        [0.0],
        [0.5],
        [complex(0)],
        ["0"],
        [None],
        [np.nan],
        [np.inf],
        [-(1 << 62) - 1],
        [1 << 62],
        [np.uint64(1 << 63)],
    ],
)
def test_frequency_design_rejects_invalid_indices_before_native_evaluation(
    method: str, indices: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: dict = {}
    adapter = _fake_adapter(method, monkeypatch, calls)
    with pytest.raises(ValueError):
        adapter.design_matsubara(indices)
    assert not any(key.startswith("frequency_") for key in calls)


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_frequency_design_empty_extreme_and_overflow_inputs(
    method: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: dict = {}
    adapter = _fake_adapter(method, monkeypatch, calls)
    assert adapter.design_matsubara(np.array([], dtype=np.int64)).shape == (0, 2)
    assert not calls
    indices = np.array([-(1 << 62), (1 << 62) - 1], dtype=np.int64)
    assert np.all(np.isfinite(adapter.design_matsubara(indices)))
    if method == "sparse_ir":
        np.testing.assert_array_equal(
            calls["frequency_reduced"], [-(1 << 63) + 1, (1 << 63) - 1]
        )
    # Construction remains valid, but the smallest nonzero physical frequency
    # cannot be represented for this beta. Do not pass infinity to a package.
    adapter.beta = np.nextafter(0.0, 1.0)
    calls.clear()
    with pytest.raises(ValueError, match="frequencies must be finite"):
        adapter.design_matsubara(np.array([0]))
    assert not calls


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
@pytest.mark.parametrize("failure", ["shape", "nonfinite", "exception"])
def test_frequency_design_native_failures_are_classified(
    method: str, failure: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    adapter = _fake_adapter(method, monkeypatch)

    def invalid(*args: object) -> np.ndarray:
        if failure == "exception":
            raise RuntimeError("native evaluation failed")
        if failure == "shape":
            return np.ones((3, 1))
        return np.full((2, 2), complex(np.nan))

    monkeypatch.setattr(
        adapter._basis, "eval_dlr_freq" if method == "pydlr" else "uhat", invalid
    )
    with pytest.raises(adapters.AdapterRuntimeError, match="Matsubara design") as error:
        adapter.design_matsubara(np.array([-1, 0]))
    assert error.value.status == "runtime_failed"


@pytest.mark.parametrize("method", adapters.EXTERNAL_METHODS)
def test_actual_frequency_design_physical_poles_and_complex_matrix_batches(
    method: str,
) -> None:
    pytest.importorskip(method)
    beta, wmax = 3.7, 2.4
    adapter = adapters.build_adapter(method, beta, wmax, 1e-11)
    poles = np.array([-1.9, 0.35, 1.6])
    weights = np.array([0.2, 0.5, 0.3])
    times = np.linspace(0, beta, 137)
    indices = np.array([-19, -5, -1, 0, 2, 13])
    multiplier = np.array([[1.0 + 0.2j, 2.0j], [0.5j, -0.25 + 0.1j]])
    physical_time = (
        -np.exp(-times[:, None] * poles - np.logaddexp(0, -beta * poles)) @ weights
    )
    values = physical_time[:, None, None] * multiplier
    # Use the common solve to exercise exactly the design API consumed by the
    # comparison. The imaginary-time observations already carry the minus sign.
    coefficients = np.linalg.lstsq(
        adapter.design_tau(times), values.reshape(times.size, -1), rcond=None
    )[0].reshape(adapter.rank, 2, 2)
    prediction = np.tensordot(
        adapter.design_matsubara(indices), coefficients, axes=(1, 0)
    )
    with mp.workdps(70):
        exact = []
        for index in indices:
            omega = (2 * int(index) + 1) * mp.pi / mp.mpf(str(beta))
            exact.append(
                complex(
                    sum(
                        mp.mpf(str(weight)) / (mp.j * omega - mp.mpf(str(pole)))
                        for pole, weight in zip(poles, weights, strict=True)
                    )
                )
            )
        # Independent high-precision Fourier quadrature checks the sign,
        # physical beta scaling, and ordinary-index convention on one pole.
        pole = mp.mpf("0.35")
        for index in (-3, 0, 2):
            omega = (2 * index + 1) * mp.pi / mp.mpf(str(beta))
            integral = mp.quad(
                lambda tau, omega=omega: (
                    -mp.exp(mp.j * omega * tau - tau * pole)
                    / (1 + mp.exp(-mp.mpf(str(beta)) * pole))
                ),
                [0, mp.mpf(str(beta))],
            )
            assert abs(integral - 1 / (mp.j * omega - pole)) < mp.mpf("1e-60")
    expected = np.asarray(exact)[:, None, None] * multiplier
    np.testing.assert_allclose(prediction, expected, atol=2e-9, rtol=2e-9)
    # Native time fit/evaluation remains compatible with the shared frequency
    # design, including a complex nonsymmetric matrix RHS.
    native_coefficients = adapter.fit(values, times)
    native_prediction = np.tensordot(
        adapter.design_matsubara(indices), native_coefficients, axes=(1, 0)
    )
    np.testing.assert_allclose(native_prediction, expected, atol=2e-9, rtol=2e-9)
    np.testing.assert_allclose(
        adapter.evaluate(native_coefficients, times), values, atol=2e-9
    )


def test_actual_pydlr_raw_kernel_transform_orientation() -> None:
    pytest.importorskip("pydlr")
    beta = 3.7
    adapter = adapters.build_adapter("pydlr", beta, 2.4, 1e-10)
    indices = np.array([-4, -1, 0, 3])
    # A single native DLR basis column is a positive raw kernel. Its transform
    # has the opposite sign to the physical negative single-pole Green function.
    column = adapter.rank // 3
    pole = float(adapter._basis.dlrrf[column] / beta)
    transformed = adapter.design_matsubara(indices)[:, column]
    with mp.workdps(70):
        b, p = mp.mpf(str(beta)), mp.mpf(str(pole))
        for row, index in enumerate(indices):
            omega = (2 * int(index) + 1) * mp.pi / b
            integral = mp.quad(
                lambda tau, omega=omega: (
                    mp.exp(mp.j * omega * tau - tau * p) / (1 + mp.exp(-b * p))
                ),
                [0, b],
            )
            np.testing.assert_allclose(transformed[row], complex(integral), atol=2e-14)
    sample = np.array([0.0, 0.37, beta])
    raw_time = np.exp(-sample * pole - np.logaddexp(0, -beta * pole))
    np.testing.assert_allclose(
        adapter.design_tau(sample)[:, column], raw_time, atol=2e-14
    )
