import mpmath as mp
import numpy as np

from kernelgecp import (
    FermionicKernel,
    GECPConfig,
    certified_pivot,
    gecp,
    interval_certified_pivot,
)


def objective(t: float, omega: float) -> float:
    return 1.0 - (t - 0.25) ** 2 - (omega + 0.4) ** 2


def test_certificate_contains_dense_grid_maximum() -> None:
    certificate = certified_pivot(
        objective,
        omega_bounds=(-1.0, 1.0),
        lipschitz=(1.5, 2.8),
        abs_tol=2e-3,
        rel_tol=0,
        max_cells=20_000,
    )
    t = np.linspace(0, 1, 401)
    omega = np.linspace(-1, 1, 801)
    dense_maximum = float(np.max(np.abs(objective(t[:, None], omega[None, :]))))
    assert certificate.certified
    assert certificate.lower_bound <= dense_maximum <= certificate.upper_bound


def test_exhausted_budget_is_not_certified() -> None:
    certificate = certified_pivot(
        objective,
        omega_bounds=(-1.0, 1.0),
        lipschitz=(1.5, 2.8),
        abs_tol=0,
        rel_tol=0,
        max_cells=1,
    )
    assert not certificate.certified
    assert certificate.termination_reason == "cell_budget"


def interval_objective(t: object, omega: object) -> object:
    return 1 - (t - mp.mpf("0.25")) ** 2 - (omega + mp.mpf("0.4")) ** 2


def interval_objective_gradient(t: object, omega: object) -> tuple[object, object]:
    return -2 * (t - mp.mpf("0.25")), -2 * (omega + mp.mpf("0.4"))


def test_interval_certificate_contains_known_maximum() -> None:
    certificate = interval_certified_pivot(
        interval_objective,
        interval_objective,
        omega_bounds=(-1.0, 1.0),
        abs_tol=2e-4,
        rel_tol=0,
        max_cells=20_000,
        gradient_interval_function=interval_objective_gradient,
    )
    assert certificate.certified
    assert certificate.lower_bound <= 1.5225 <= certificate.upper_bound


def test_interval_certificate_reports_budget_exhaustion() -> None:
    certificate = interval_certified_pivot(
        interval_objective,
        interval_objective,
        omega_bounds=(-1.0, 1.0),
        abs_tol=0,
        rel_tol=0,
        max_cells=2,
    )
    assert not certificate.certified
    assert certificate.termination_reason == "cell_budget"


def test_adaptive_fermionic_gecp_retains_certificates() -> None:
    result = gecp(
        FermionicKernel(1.0),
        config=GECPConfig(
            tol=1e-2,
            max_rank=2,
            pivot="adaptive",
            precision_bits=96,
            certificate_abs_tol=2e-3,
            certificate_rel_tol=0,
            max_cells=20_000,
        ),
    )
    assert result.precision_bits == 96
    assert result.high_precision is not None
    assert result.rank == len(result.pivot_certificates)
    assert all(certificate.certified for certificate in result.pivot_certificates)
    assert all(
        certificate.lower_bound <= certificate.upper_bound
        for certificate in result.pivot_certificates
    )
