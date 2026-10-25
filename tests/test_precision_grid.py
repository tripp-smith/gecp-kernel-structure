"""Independent exact-input QR checks of finite-grid pivot certificates."""

from dataclasses import FrozenInstanceError
from fractions import Fraction
from typing import Any

import mpmath as mp
import numpy as np
import pytest

import kernelgecp.precision_grid as grid
from kernelgecp import FermionicBasis
from kernelgecp.kernels import FermionicKernel


@pytest.fixture(scope="module")
def basis() -> FermionicBasis:
    return FermionicBasis(
        beta=3.7, wmax=2.4, tol=1e-7, grid_order=6, precision="reference"
    )


def _independent_residuals(basis: FermionicBasis, tau: Any, omega: Any) -> Any:
    context = mp.mp.clone()
    context.prec = 320
    beta = context.mpf(basis.beta)

    def kernel(time: float, frequency: float) -> Any:
        return context.exp(-context.mpf(float(time)) / beta * float(frequency)) / (
            1 + context.exp(-context.mpf(float(frequency)))
        )

    core = context.matrix(
        [[kernel(t, w) for w in basis.omega_nodes] for t in basis.tau_nodes]
    )
    residuals = []
    for time in tau:
        # QR is independent of the production interval Gaussian LU.
        coefficients, _ = context.qr_solve(
            core.T, context.matrix([kernel(time, w) for w in basis.omega_nodes])
        )
        residuals.append(
            [
                kernel(time, w)
                - context.fsum(
                    coefficients[k] * kernel(node, w)
                    for k, node in enumerate(basis.tau_nodes)
                )
                for w in omega
            ]
        )
    return context, residuals


def test_true_grid_quality_and_discrepancy_against_independent_qr(
    basis: FermionicBasis,
) -> None:
    tau = basis.beta * np.array([0.071, 0.23, 0.41, 0.63, 0.91])
    omega = basis.cutoff * np.array([-0.91, -0.41, 0.03, 0.37, 0.83])
    report = grid.certify_grid_pivot_precision(basis, tau, omega)
    assert report.certified, report.termination_reason
    assert report.eta is not None and 0.99 < report.eta <= 1
    assert report.discrepancy_upper is not None
    assert report.discrepancy_upper < 1e-12
    assert report.candidate_index is not None
    context, ideal = _independent_residuals(basis, tau, omega)
    kernel = FermionicKernel(basis.cutoff)
    actual = kernel(tau[:, None] / basis.beta, omega[None, :]) - basis.design_tau(
        tau
    ) @ kernel(basis.tau_nodes[:, None] / basis.beta, omega[None, :])
    p, q = report.candidate_index
    for i in range(len(tau)):
        for j in range(len(omega)):
            assert abs(ideal[i][j] - float(actual[i, j])) <= context.mpf(
                report.discrepancy_upper
            )
            assert context.mpf(report.eta) * abs(ideal[i][j]) <= abs(ideal[p][q])
    assert report.objective_scope.endswith("supplied_finite_grid")
    assert not report.continuous_certificate


def test_deterministic_ties_context_independence_and_input_binding(
    basis: FermionicBasis,
) -> None:
    tau = np.full(2, basis.beta * 0.23)
    omega = np.full(2, basis.cutoff * 0.41)
    old_mp, old_iv = mp.mp.prec, mp.iv.prec
    try:
        mp.mp.prec, mp.iv.prec = 71, 91
        first = grid.certify_grid_pivot_precision(basis, tau, omega)
        assert (mp.mp.prec, mp.iv.prec) == (71, 91)
        mp.mp.prec, mp.iv.prec = 140, 200
        second = grid.certify_grid_pivot_precision(basis, tau, omega)
        assert (mp.mp.prec, mp.iv.prec) == (140, 200)
    finally:
        mp.mp.prec, mp.iv.prec = old_mp, old_iv
    assert first == second
    assert first.certified
    assert first.candidate_index == (0, 0)
    assert first.exact_tie_count == first.near_tie_count == 4
    tau[0] = basis.beta * 0.24
    assert (
        first.content_sha256
        != grid.certify_grid_pivot_precision(basis, tau, omega).content_sha256
    )
    with pytest.raises(FrozenInstanceError):
        first.eta = 0.0  # type: ignore[misc]


@pytest.mark.parametrize(
    "budget,reason", [("rank", "rank_budget"), ("work", "work_budget")]
)
def test_budgets_are_explicit_inconclusive_reports(
    basis: FermionicBasis, budget: str, reason: str
) -> None:
    options = {"max_rank": 1} if budget == "rank" else {"max_work": 1}
    report = grid.certify_grid_pivot_precision(basis, [0.1], [0.2], **options)
    assert not report.certified and report.termination_reason == reason
    assert report.discrepancy_upper is report.eta is None


def test_interval_failure_is_not_singularity_or_zero_error(
    basis: FermionicBasis, monkeypatch: pytest.MonkeyPatch
) -> None:
    def failed(*args: Any) -> None:
        raise ArithmeticError("pivot interval contains zero")

    monkeypatch.setattr(grid, "_IntervalBasis", failed)
    report = grid.certify_grid_pivot_precision(basis, [0.1], [0.2])
    assert not report.certified
    assert report.termination_reason == "interval_inconclusive:ArithmeticError"
    assert report.discrepancy_upper is report.eta is None


def test_malformed_complex_design_cannot_silently_drop_imaginary_residual(
    basis: FermionicBasis, monkeypatch: pytest.MonkeyPatch
) -> None:
    original = FermionicBasis.design_tau

    def complex_design(self: FermionicBasis, tau: Any) -> Any:
        return original(self, tau) + 0.001j

    monkeypatch.setattr(FermionicBasis, "design_tau", complex_design)
    report = grid.certify_grid_pivot_precision(basis, [0.1, 0.8], [-0.3, 0.2])
    assert not report.certified
    assert report.termination_reason == "invalid_real_rounded_residual"
    assert report.discrepancy_upper is None


def test_exact_zero_grid_has_explicit_zero_branch(
    basis: FermionicBasis, monkeypatch: pytest.MonkeyPatch
) -> None:
    from mpmath.ctx_iv import MPIntervalContext

    class ExactConstant:
        def __init__(self, *args: Any) -> None:
            self.ctx = MPIntervalContext()
            self.ctx.prec = 256
            self.transpose_lu = self
            self.pivot_lower_bounds = [1.0] * basis.rank

        def kernel(self, *args: Any) -> Any:
            return self.ctx.mpf(1)

        def tau_rows(self, tau: Any) -> Any:
            return [[self.ctx.mpf(j == 0) for j in range(basis.rank)] for _ in tau]

    def constant_design(self: Any, tau: Any) -> Any:
        result = np.zeros((len(tau), basis.rank))
        result[:, 0] = 1
        return result

    class ConstantKernel:
        def __init__(self, *args: Any) -> None:
            pass

        def __call__(self, t: Any, w: Any) -> Any:
            return np.ones(np.broadcast_shapes(t.shape, w.shape))

    monkeypatch.setattr(grid, "_IntervalBasis", ExactConstant)
    monkeypatch.setattr(FermionicBasis, "design_tau", constant_design)
    monkeypatch.setattr(grid, "FermionicKernel", ConstantKernel)
    result = grid.certify_grid_pivot_precision(basis, [0.1], [0.2])
    assert result.certified and result.eta == 1
    assert result.discrepancy_upper == result.true_global_upper == 0
    assert result.termination_reason == "exact_zero_grid_residual"


@pytest.mark.parametrize(
    "tau,omega,options",
    [
        ([], [0], {}),
        ([0], [], {}),
        ([False, 0.1], [0], {}),
        ([0.1], [True, 0], {}),
        ([-0.1], [0], {}),
        ([0.1], [1e3], {}),
        ([0.1], [1j], {}),
        ([0.1], [0], {"precision_bits": True}),
        ([0.1], [0], {"max_rank": 0}),
        ([0.1], [0], {"max_work": 1.5}),
    ],
)
def test_invalid_inputs_fail_before_certification(
    basis: FermionicBasis, tau: Any, omega: Any, options: Any
) -> None:
    with pytest.raises(ValueError):
        grid.certify_grid_pivot_precision(basis, tau, omega, **options)


def test_directed_exports_preserve_positive_subnormal_and_ratio_bounds() -> None:
    half_subnormal = Fraction(1, 2**1075)
    assert grid._directed(half_subnormal, True) == np.nextafter(0.0, 1.0)
    assert grid._directed(half_subnormal, False) == 0
    exact = Fraction(765, 1027)
    assert Fraction(grid._directed(exact, False)) <= exact
    assert Fraction(grid._directed(exact, True)) >= exact
