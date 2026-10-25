"""Optional fermionic DLR/IR adapters using physical imaginary time.

Only ``build_adapter`` imports external packages.  Neither this module nor an
adapter installs dependencies or accesses a truth function.  Values are already
the physical Green's function, including its minus sign; no sign is inserted.
``fit`` uses only the explicitly supplied observations, with the first array
axis representing time and any trailing axes representing independent RHSs.
The design methods share these coefficients: ``design_tau`` uses physical
``tau`` in ``[0, beta]``, and ``design_matsubara`` uses ordinary integer indices
``n`` for ``omega_n = (2*n + 1)*pi/beta`` and the Fourier factor
``exp(+1j*omega_n*tau)``. Both return point-by-rank matrices.

API references (checked with pydlr 1.0.1 and sparse-ir 2.1.6):
https://libdlr.readthedocs.io/en/latest/pydlr.html
https://libdlr.readthedocs.io/en/latest/_modules/pydlr/pydlr.html
https://sparse-ir.readthedocs.io/en/latest/basis.html
https://sparse-ir.readthedocs.io/en/latest/sampling.html
https://github.com/SpM-lab/sparse-ir/blob/mainline/src/sparse_ir/sampling.py
"""

from __future__ import annotations

import importlib
import importlib.metadata
import inspect
import math
from typing import Any, Protocol

import numpy as np

EXTERNAL_METHODS = ("pydlr", "sparse_ir")


class AdapterUnavailableError(ImportError):
    """An optional package or one of its import dependencies is unavailable."""

    status = "import_unavailable"


class AdapterRuntimeError(RuntimeError):
    """An imported package could not construct, fit, or evaluate its basis."""

    status = "runtime_failed"


class FermionicAdapter(Protocol):
    rank: int
    beta: float
    tau_nodes: np.ndarray
    metadata: dict[str, Any]

    def design_matrix(self, tau: np.ndarray) -> np.ndarray: ...

    def design_tau(self, tau: np.ndarray) -> np.ndarray: ...

    def design_matsubara(self, indices: np.ndarray) -> np.ndarray: ...

    def fit(self, values: np.ndarray, tau: np.ndarray) -> np.ndarray: ...

    def evaluate(self, coeffs: np.ndarray, tau: np.ndarray) -> np.ndarray: ...


def _numeric_array(value: object, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype.kind not in "buifc":
        raise ValueError(f"{name} must contain numeric values")
    array = np.asarray(array, dtype=complex if np.iscomplexobj(array) else float)
    if array.ndim == 0 or not array.size or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite array")
    return array


def _tau_array(tau: object, beta: float, *, fitting: bool = False) -> np.ndarray:
    array = _numeric_array(tau, "tau")
    if array.ndim != 1 or np.iscomplexobj(array):
        raise ValueError("tau must be a one-dimensional real array")
    if np.any(array < 0) or np.any(array > beta):
        raise ValueError("tau must lie in [0, beta]")
    if fitting and np.unique(array).size != array.size:
        raise ValueError("fitting tau must be pairwise distinct")
    # sparse-ir distinguishes signed zero at the fermionic discontinuity.
    # This interface always uses the positive-time one-sided endpoint.
    array = array.copy()
    array[array == 0] = 0.0
    return array


def _matsubara_arrays(indices: object, beta: float) -> tuple[np.ndarray, np.ndarray]:
    # Inspect before numeric coercion: np.asarray([0, True]) otherwise silently
    # converts the boolean to an integer. Floats, including integral floats,
    # are not ordinary integer indices in this interface.
    array = np.asarray(indices, dtype=object)
    if array.ndim != 1:
        raise ValueError("Matsubara indices must be a one-dimensional integer array")
    if any(
        isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
        for value in array
    ):
        raise ValueError("Matsubara indices must be integers, excluding booleans")
    # sparse-ir's reduced index is a signed int64. Compute it only after
    # checking the ordinary index, so 2*n + 1 cannot wrap around.
    lower, upper = -(1 << 62), (1 << 62) - 1
    if any(value < lower or value > upper for value in array):
        raise ValueError("Matsubara reduced indices 2*n+1 must fit signed int64")
    ordinary = np.asarray(array, dtype=np.int64)
    reduced = 2 * ordinary + 1
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        frequencies = reduced.astype(float) * (np.pi / beta)
    if not np.all(np.isfinite(frequencies)):
        raise ValueError("physical Matsubara frequencies must be finite")
    return reduced, frequencies


class _Adapter:
    def __init__(
        self, method: str, module: Any, beta: float, wmax: float, tol: float
    ) -> None:
        self.method = method
        self.beta = beta
        self.wmax = wmax
        self.tol = tol
        distribution = "sparse-ir" if method == "sparse_ir" else method
        try:
            version = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            version = str(getattr(module, "__version__", "unknown"))
        self.metadata = {
            "package": distribution,
            "package_version": version,
            "statistics": "fermionic",
            "beta": beta,
            "wmax": wmax,
            "lambda": beta * wmax,
            "requested_tolerance": tol,
            "precision_bits": 53,
            "precision_scope": "exported_design_and_fit_not_internal_construction",
            "native_tolerance_meaning": (
                "native_interpolative_decomposition_eps"
                if method == "pydlr"
                else "relative_singular_value_truncation"
            ),
            "certification_status": "not_certified",
            "matsubara_index_convention": "ordinary_integer_n",
            "matsubara_frequency_convention": "omega_n=(2*n+1)*pi/beta",
            "fourier_convention": "integral_0^beta exp(+i*omega_n*tau) G(tau) dtau",
            "physical_green_function_sign": "observations_already_include_minus_sign",
        }

    def _check_basis(self) -> None:
        if self.rank < 1:
            raise ValueError("external basis must have positive rank")
        self.tau_nodes = _tau_array(self.tau_nodes, self.beta, fitting=True)
        if self.tau_nodes.size < self.rank:
            raise ValueError("external basis has fewer default nodes than its rank")
        self.metadata["basis_rank"] = self.rank

    def design_tau(self, tau: np.ndarray) -> np.ndarray:
        """Alias for the existing physical-time design matrix."""
        return self.design_matrix(tau)

    def _values(self, values: object, tau: np.ndarray) -> np.ndarray:
        array = _numeric_array(values, "values")
        if array.shape[0] != tau.size:
            raise ValueError("values first axis must match tau")
        if tau.size < self.rank:
            raise ValueError("fitting requires at least rank observations")
        return array

    def _coefficients(self, coeffs: object) -> np.ndarray:
        array = _numeric_array(coeffs, "coeffs")
        if array.shape[0] != self.rank:
            raise ValueError("coeffs first axis must match basis rank")
        return array

    def _output(self, value: object, shape: tuple[int, ...]) -> np.ndarray:
        array = _numeric_array(value, "external output")
        if array.shape != shape:
            raise ValueError(f"external output shape {array.shape} != {shape}")
        return array


class _PydlrAdapter(_Adapter):
    def __init__(self, module: Any, beta: float, wmax: float, tol: float) -> None:
        super().__init__("pydlr", module, beta, wmax, tol)
        self._basis = module.dlr(lamb=beta * wmax, eps=tol, xi=-1, python_impl=True)
        self.rank = int(len(self._basis))
        self.tau_nodes = np.asarray(self._basis.get_tau(beta))
        self._check_basis()

    def design_matrix(self, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta)
        try:
            matrix = self._basis.eval_dlr_tau(
                np.eye(self.rank)[:, :, None], times, self.beta
            )
            return self._output(matrix, (times.size, self.rank, 1))[:, :, 0]
        except Exception as exc:
            raise AdapterRuntimeError(f"pydlr design matrix failed: {exc}") from exc

    def design_matsubara(self, indices: np.ndarray) -> np.ndarray:
        """Fourier design at ordinary fermionic Matsubara indices."""
        _, frequencies = _matsubara_arrays(indices, self.beta)
        if not frequencies.size:
            return np.empty((0, self.rank), dtype=complex)
        try:
            # Native eval_dlr_freq conjugates its result and transposes a 3D
            # RHS. Real 2D identity coefficients extract only the correctly
            # oriented transformed basis; callers can then multiply arbitrary
            # complex coefficients without conjugating or transposing them.
            matrix = self._basis.eval_dlr_freq(
                np.eye(self.rank), 1j * frequencies, self.beta
            )
            return self._output(matrix, (frequencies.size, self.rank))
        except Exception as exc:
            raise AdapterRuntimeError(f"pydlr Matsubara design failed: {exc}") from exc

    def fit(self, values: np.ndarray, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta, fitting=True)
        array = self._values(values, times)
        try:
            # pydlr's arbitrary-time fitter requires exactly three dimensions.
            fitted = self._basis.lstsq_dlr_from_tau(
                times, array.reshape(times.size, -1, 1), self.beta
            )
            shape = (self.rank, *array.shape[1:])
            fitted = self._output(fitted, (self.rank, array.size // times.size, 1))
            return fitted.reshape(shape)
        except Exception as exc:
            raise AdapterRuntimeError(f"pydlr fit failed: {exc}") from exc

    def evaluate(self, coeffs: np.ndarray, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta)
        array = self._coefficients(coeffs)
        try:
            values = self._basis.eval_dlr_tau(
                array.reshape(self.rank, -1, 1), times, self.beta
            )
            values = self._output(values, (times.size, array.size // self.rank, 1))
            return values.reshape(times.size, *array.shape[1:])
        except Exception as exc:
            raise AdapterRuntimeError(f"pydlr evaluation failed: {exc}") from exc


class _SparseIRAdapter(_Adapter):
    def __init__(self, module: Any, beta: float, wmax: float, tol: float) -> None:
        super().__init__("sparse_ir", module, beta, wmax, tol)
        self._module = module
        self._basis = module.FiniteTempBasis("F", beta=beta, wmax=wmax, eps=tol)
        self.rank = int(self._basis.size)
        # The 1.x constructor always uses positive times; newer constructors
        # expose a choice, which we set explicitly instead of relying on default.
        options = (
            {"use_positive_taus": True}
            if "use_positive_taus" in inspect.signature(module.TauSampling).parameters
            else {}
        )
        self.tau_nodes = np.asarray(
            module.TauSampling(self._basis, **options).sampling_points
        )
        self._check_basis()

    def design_matrix(self, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta)
        try:
            return self._output(self._basis.u(times).T, (times.size, self.rank))
        except Exception as exc:
            raise AdapterRuntimeError(f"sparse_ir design matrix failed: {exc}") from exc

    def design_matsubara(self, indices: np.ndarray) -> np.ndarray:
        """Fourier design, converting ordinary n to sparse-ir's odd 2*n+1."""
        reduced, _ = _matsubara_arrays(indices, self.beta)
        if not reduced.size:
            return np.empty((0, self.rank), dtype=complex)
        try:
            matrix = self._basis.uhat(reduced).T
            return self._output(matrix, (reduced.size, self.rank))
        except Exception as exc:
            raise AdapterRuntimeError(
                f"sparse_ir Matsubara design failed: {exc}"
            ) from exc

    def fit(self, values: np.ndarray, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta, fitting=True)
        array = self._values(values, times)
        try:
            sampling = self._module.TauSampling(self._basis, sampling_points=times)
            fitted = sampling.fit(array, axis=0)
            return self._output(fitted, (self.rank, *array.shape[1:]))
        except Exception as exc:
            raise AdapterRuntimeError(f"sparse_ir fit failed: {exc}") from exc

    def evaluate(self, coeffs: np.ndarray, tau: np.ndarray) -> np.ndarray:
        times = _tau_array(tau, self.beta)
        array = self._coefficients(coeffs)
        try:
            values = np.tensordot(self.design_matrix(times), array, axes=(1, 0))
            return self._output(values, (times.size, *array.shape[1:]))
        except Exception as exc:
            raise AdapterRuntimeError(f"sparse_ir evaluation failed: {exc}") from exc


def build_adapter(
    method: str, beta: float, wmax: float, tol: float
) -> FermionicAdapter:
    """Build an optional adapter, distinguishing missing imports from failures.

    The external tolerance chooses a basis and is not an error certificate.
    The caller may use ``design_tau`` (also named ``design_matrix``) and
    ``design_matsubara`` with one common fitting algorithm for every method;
    ``fit`` separately exposes each package's native time-domain fitter.
    """
    if method not in EXTERNAL_METHODS:
        raise ValueError(f"unknown external method: {method!r}")
    for name, value in (("beta", beta), ("wmax", wmax), ("tol", tol)):
        if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if tol >= 1 or not math.isfinite(beta * wmax) or beta * wmax <= 0:
        raise ValueError("require tol < 1 and finite positive beta*wmax")
    try:
        module = importlib.import_module(method)
    except ImportError as exc:
        raise AdapterUnavailableError(f"{method} import unavailable: {exc}") from exc
    except Exception as exc:
        raise AdapterRuntimeError(f"{method} import failed: {exc}") from exc
    try:
        cls = _PydlrAdapter if method == "pydlr" else _SparseIRAdapter
        return cls(module, float(beta), float(wmax), float(tol))
    except Exception as exc:
        raise AdapterRuntimeError(f"{method} construction failed: {exc}") from exc
