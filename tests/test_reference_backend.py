"""Optional reference backend contracts, provider errors and numerical parity.

Fake native modules test error routing and bounded transfer, not numerical proof.
Actual-extension tests skip when only the universal Python package is installed.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import mpmath as mp
import numpy as np
import pytest

from kernelgecp import FermionicBasis, FermionicFit, NumericalBreakdownError
from kernelgecp._basis_numerics import BasisNumerics
from kernelgecp._reference_backend import RustReferenceBackend


def factors() -> tuple[Any, list[int], Any, Any]:
    context = mp.mp.clone()
    context.prec = 128
    return context, [0, 1], context.eye(2), context.eye(2)


class FakeHandle:
    calls: ClassVar[list[list[float]]] = []
    payload: ClassVar[bytes | None] = None

    def __init__(self, bits: int, omega: list[float], *args: object) -> None:
        assert bits == 128
        self.rank = len(omega)

    def design_tau(self, times: list[float]) -> bytes:
        self.calls.append(times)
        if self.payload is not None:
            return self.payload
        return np.asarray([[t, 1 - t] for t in times], dtype="<f8").tobytes()


def fake_module(**overrides: object) -> Any:
    values: dict[str, object] = {
        "NATIVE_PROTOCOL_VERSION": 1,
        "WORKING_PRECISION_BITS": 128,
        "ReferenceCardinal128": FakeHandle,
        "build_info": lambda: {
            "protocol_version": 1,
            "working_precision_bits": 128,
            "backend": "fake_for_behavior_tests",
        },
    }
    values.update(overrides)
    return SimpleNamespace(**values)


@pytest.fixture(autouse=True)
def clean_fake() -> None:
    FakeHandle.calls = []
    FakeHandle.payload = None


def make_backend(
    monkeypatch: pytest.MonkeyPatch, **overrides: object
) -> RustReferenceBackend:
    monkeypatch.setattr(
        importlib, "import_module", lambda name: fake_module(**overrides)
    )
    context, ordering, lower, upper = factors()
    return RustReferenceBackend(np.array([-1.0, 1.0]), ordering, lower, upper, context)


def test_adapter_batches_and_decodes_little_endian(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend = make_backend(monkeypatch)
    times = np.linspace(0, 1, 600)
    result = backend.design_tau(times)
    np.testing.assert_array_equal(result, np.column_stack((times, 1 - times)))
    assert [len(batch) for batch in FakeHandle.calls] == [256, 256, 88]
    assert backend.rank == 2
    assert backend.design_tau(np.array([], dtype=float)).shape == (0, 2)


@pytest.mark.parametrize(
    "key,value",
    [
        ("NATIVE_PROTOCOL_VERSION", 2),
        ("NATIVE_PROTOCOL_VERSION", True),
        ("WORKING_PRECISION_BITS", 53),
    ],
)
def test_adapter_rejects_wrong_protocol(
    monkeypatch: pytest.MonkeyPatch, key: str, value: object
) -> None:
    with pytest.raises(NumericalBreakdownError, match="incompatible"):
        make_backend(monkeypatch, **{key: value})


@pytest.mark.parametrize(
    "info",
    [
        None,
        {},
        {"protocol_version": True, "working_precision_bits": 128},
        {"protocol_version": 1, "working_precision_bits": 53},
        {"x": []},
    ],
)
def test_adapter_rejects_invalid_build_info(
    monkeypatch: pytest.MonkeyPatch, info: object
) -> None:
    with pytest.raises(NumericalBreakdownError, match="build info"):
        make_backend(monkeypatch, build_info=lambda: info)


@pytest.mark.parametrize(
    "error",
    [
        ImportError("missing"),
        OSError("bad loader"),
        RuntimeError("initialization failed"),
    ],
)
def test_missing_native_has_actionable_error(
    monkeypatch: pytest.MonkeyPatch, error: Exception
) -> None:
    def fail(name: str) -> Any:
        raise error

    monkeypatch.setattr(importlib, "import_module", fail)
    context, ordering, lower, upper = factors()
    with pytest.raises(NumericalBreakdownError, match="reference_backend='python'"):
        RustReferenceBackend(np.array([-1.0, 1.0]), ordering, lower, upper, context)


@pytest.mark.parametrize(
    "payload", [b"", b"x" * 17, np.array([np.nan, 1.0], dtype="<f8").tobytes()]
)
def test_adapter_rejects_bad_native_payload(
    monkeypatch: pytest.MonkeyPatch, payload: bytes
) -> None:
    backend = make_backend(monkeypatch)
    FakeHandle.payload = payload
    with pytest.raises(NumericalBreakdownError):
        backend.design_tau(np.array([0.5]))


def test_native_runtime_failure_never_falls_back(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend = make_backend(monkeypatch)

    def fail(times: list[float]) -> bytes:
        raise ArithmeticError("MPFR underflow")

    monkeypatch.setattr(backend._handle, "design_tau", fail)
    with pytest.raises(NumericalBreakdownError, match="no Python fallback"):
        backend.design_tau(np.array([0.5]))


@pytest.mark.parametrize(
    "times",
    [
        np.array([True]),
        np.array([0.2, np.inf]),
        np.array([-0.01]),
        np.array([1.01]),
        np.ones((2, 1)),
    ],
)
def test_adapter_validates_queries(monkeypatch: pytest.MonkeyPatch, times: Any) -> None:
    backend = make_backend(monkeypatch)
    with pytest.raises(ValueError, match="normalized times"):
        backend.design_tau(times)


def test_reference_factor_transfer_is_exact(monkeypatch: pytest.MonkeyPatch) -> None:
    observed: list[object] = []

    class Packed(FakeHandle):
        def __init__(self, bits: int, omega: list[float], *args: object) -> None:
            observed.extend(args)
            super().__init__(bits, omega, *args)

    make_backend(monkeypatch, ReferenceCardinal128=Packed)
    assert observed == [
        [0, 1],
        [["1.0", "0.0"], ["0.0", "1.0"]],
        [["1.0", "0.0"], ["0.0", "1.0"]],
    ]


@pytest.mark.parametrize("bits", [64, 256])
def test_incompatible_native_precision_rejected_before_selection(
    monkeypatch: pytest.MonkeyPatch, bits: int
) -> None:
    def fail(*args: object, **kwargs: object) -> Any:
        raise AssertionError("selection must not run")

    monkeypatch.setattr(FermionicBasis, "_select", fail)
    with pytest.raises(ValueError, match="precision_bits=128"):
        FermionicBasis(2, 1, reference_backend="rust", precision_bits=bits)


def test_python_default_and_float64_path_do_not_import_native(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(name: str) -> Any:
        raise AssertionError("native import must not run")

    monkeypatch.setattr(importlib, "import_module", fail)
    for backend in ("python", "rust"):
        numerics = BasisNumerics(
            2, 2, [0, 2], [-1, 1], precision_bits=53, reference_backend=backend
        )
        assert numerics.time_evaluation_backend == "scipy_lapack_float64"
    basis = FermionicBasis(2, 1, tol=1e-6, precision="reference", grid_order=4)
    assert basis.diagnostics.requested_reference_backend == "python"
    assert basis.diagnostics.native_backend_info == ()


def test_archive_backend_history_does_not_import_native(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    basis = FermionicBasis(2, 1, tol=1e-6, precision="reference", grid_order=4)
    fit = FermionicFit(basis, np.ones(basis.rank))
    path = tmp_path / "fit.npz"
    fit.save(path)

    def fail(name: str) -> Any:
        raise AssertionError("archive history must not import native")

    monkeypatch.setattr(importlib, "import_module", fail)
    restored = FermionicFit.load(path)
    assert restored.basis.diagnostics.requested_reference_backend == "python"
    assert restored.basis.diagnostics.history_status == "restored_unverified_history"
    np.testing.assert_array_equal(restored.evaluate_tau([0.3]), fit.evaluate_tau([0.3]))


@pytest.mark.parametrize("beta,wmax", [(2.0, 1.0), (3.7, 2.4)])
def test_actual_native_time_designs_and_matrix_axes_match_python(
    beta: float, wmax: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("kernelgecp._native")
    options = dict(tol=1e-9, precision="reference", grid_order=6)
    python = FermionicBasis(beta, wmax, **options)
    rust = FermionicBasis(beta, wmax, **options, reference_backend="rust")
    np.testing.assert_array_equal(python.tau_nodes, rust.tau_nodes)
    np.testing.assert_array_equal(python.omega_nodes, rust.omega_nodes)
    tau = beta * (np.arange(601) + 0.37) / 602
    tau = np.concatenate(([0.0], tau, [beta]))

    def fail_python_solve(*args: object) -> Any:
        raise AssertionError("native time path must not call the Python solve")

    monkeypatch.setattr(rust._numerics, "_solve_reference", fail_python_solve)
    np.testing.assert_array_equal(rust.design_tau(tau), python.design_tau(tau))
    coefficients = np.arange(4 * rust.rank, dtype=float).reshape(
        2, rust.rank, 2
    ) + 1j * np.linspace(0.1, 0.7, 4 * rust.rank).reshape(2, rust.rank, 2)
    np.testing.assert_array_equal(
        rust.evaluate_tau(coefficients, tau, axis=1),
        python.evaluate_tau(coefficients, tau, axis=1),
    )
    assert rust.diagnostics.time_evaluation_backend == "rust_mpfr_128"
    assert rust.diagnostics.frequency_evaluation_backend == "python_mpmath_reference"
    assert rust.diagnostics.fit_backend == "numpy_svd_float64_or_complex128"


def test_actual_native_fit_and_archive_are_portable(tmp_path: Path) -> None:
    pytest.importorskip("kernelgecp._native")
    basis = FermionicBasis(
        3.7,
        2.4,
        tol=1e-9,
        precision="reference",
        grid_order=6,
        reference_backend="rust",
    )
    tau = np.linspace(0, basis.beta, 73)
    values = -np.exp(-0.31 * tau) / (1 + np.exp(-0.31 * basis.beta))
    amplitude = np.array([[1, 0.2 + 0.5j], [0.2 - 0.5j, 0.75]])
    fit = basis.fit_tau(values[:, None, None] * amplitude, tau)
    path = tmp_path / "native-fit.npz"
    fit.save(path)
    python_restore = FermionicFit.load(path)
    native_restore = FermionicFit.load(path, reference_backend="rust")
    assert python_restore.basis.diagnostics.requested_reference_backend == "python"
    assert (
        python_restore.basis.diagnostics.time_evaluation_backend
        == "python_mpmath_reference"
    )
    assert native_restore.basis.diagnostics.time_evaluation_backend == "rust_mpfr_128"
    query = np.array([0, 0.67, 1.413, 3.7])
    np.testing.assert_array_equal(
        fit.evaluate_tau(query), python_restore.evaluate_tau(query)
    )
    np.testing.assert_array_equal(
        fit.evaluate_tau(query), native_restore.evaluate_tau(query)
    )
    np.testing.assert_array_equal(
        fit.evaluate_matsubara([-3, -1, 0, 5]),
        python_restore.evaluate_matsubara([-3, -1, 0, 5]),
    )
