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
from .research import (
    BlockContractionEvidence,
    CertifiedBlockContractionEvidence,
    analyze_census_file,
    analyze_census_record,
    certify_fermionic_block_contraction,
    dyadic_cutoff_scale,
    fermionic_first_residual,
    fermionic_first_step_ratio,
    fermionic_two_corner_residual,
)
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
    "BlockContractionEvidence",
    "CertifiedBlockContractionEvidence",
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
    "analyze_census_file",
    "analyze_census_record",
    "certified_pivot",
    "certify_fermionic_block_contraction",
    "interval_certified_pivot",
    "cross_approximation",
    "dlr_rank_bound",
    "dyadic_taylor_rank_bound",
    "dyadic_cutoff_scale",
    "evaluate_residual",
    "exp_family_rank_bound",
    "fermionic_dyadic_taylor_approximation",
    "fermionic_green_from_density",
    "fermionic_first_residual",
    "fermionic_first_step_ratio",
    "fermionic_two_corner_residual",
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
