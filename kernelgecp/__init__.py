"""GECP algorithms and fermionic-kernel research utilities."""

from .approximation import (
    SeparatedApproximation,
    dlr_rank_bound,
    exp_family_rank_bound,
    fermionic_separated_approximation,
)
from .certified import certified_pivot
from .cholesky import pivoted_cholesky
from .gecp import cross_approximation, evaluate_residual, gecp, gecp_matrix
from .kernels import FermionicKernel
from .sparse import sparse_representation
from .types import (
    CholeskyResult,
    GECPConfig,
    GECPResult,
    NumericalBreakdownError,
    PivotCertificate,
    SparseResult,
)

__all__ = [
    "CholeskyResult",
    "FermionicKernel",
    "GECPConfig",
    "GECPResult",
    "NumericalBreakdownError",
    "PivotCertificate",
    "SeparatedApproximation",
    "SparseResult",
    "certified_pivot",
    "cross_approximation",
    "dlr_rank_bound",
    "evaluate_residual",
    "exp_family_rank_bound",
    "fermionic_separated_approximation",
    "gecp",
    "gecp_matrix",
    "pivoted_cholesky",
    "sparse_representation",
]

__version__ = "0.1.0"
