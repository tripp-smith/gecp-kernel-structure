import numpy as np

from kernelgecp import certified_pivot


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
