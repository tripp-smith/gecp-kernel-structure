"""GECP algorithms and fermionic-kernel research utilities."""

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
    "HighPrecisionData",
    "NumericalBreakdownError",
    "PivotCertificate",
    "SeparatedApproximation",
    "SparseResult",
    "certified_pivot",
    "interval_certified_pivot",
    "cross_approximation",
    "dlr_rank_bound",
    "dyadic_taylor_rank_bound",
    "evaluate_residual",
    "exp_family_rank_bound",
    "fermionic_dyadic_taylor_approximation",
    "fermionic_separated_approximation",
    "gecp",
    "gecp_matrix",
    "pivoted_cholesky",
    "sparse_representation",
]

__version__ = "1.0.0"
