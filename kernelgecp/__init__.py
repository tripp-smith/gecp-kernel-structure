"""GECP algorithms and fermionic-kernel research utilities."""

from .applications import (
    GreenCompressionResult,
    PSDLandmarkResult,
    SparseRecoveryResult,
    SyntheticApplicationConfig,
    SyntheticApplicationSuite,
    fermionic_green_from_density,
    quadrature_weights,
    run_synthetic_application_suite,
    synthetic_sensor_points,
    synthetic_spectral_density,
)
from .approximation import (
    DyadicTaylorApproximation,
    SeparatedApproximation,
    dlr_rank_bound,
    dyadic_taylor_rank_bound,
    exp_family_rank_bound,
    fermionic_dyadic_taylor_approximation,
    fermionic_separated_approximation,
)
from .certified import certified_pivot, interval_certified_pivot
from .cholesky import pivoted_cholesky
from .gecp import cross_approximation, evaluate_residual, gecp, gecp_matrix
from .kernels import FermionicKernel
from .sparse import sparse_representation
from .types import (
    CholeskyResult,
    GECPConfig,
    GECPResult,
    HighPrecisionData,
    NumericalBreakdownError,
    PivotCertificate,
    SparseResult,
)

__all__ = [
    "CholeskyResult",
    "DyadicTaylorApproximation",
    "FermionicKernel",
    "GECPConfig",
    "GECPResult",
    "GreenCompressionResult",
    "HighPrecisionData",
    "NumericalBreakdownError",
    "PSDLandmarkResult",
    "PivotCertificate",
    "SeparatedApproximation",
    "SparseRecoveryResult",
    "SparseResult",
    "SyntheticApplicationConfig",
    "SyntheticApplicationSuite",
    "certified_pivot",
    "interval_certified_pivot",
    "cross_approximation",
    "dlr_rank_bound",
    "dyadic_taylor_rank_bound",
    "evaluate_residual",
    "exp_family_rank_bound",
    "fermionic_dyadic_taylor_approximation",
    "fermionic_green_from_density",
    "fermionic_separated_approximation",
    "gecp",
    "gecp_matrix",
    "pivoted_cholesky",
    "quadrature_weights",
    "run_synthetic_application_suite",
    "sparse_representation",
    "synthetic_sensor_points",
    "synthetic_spectral_density",
]

__version__ = "1.0.0"
