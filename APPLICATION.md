# Application guide

The public package supplies stable fermionic-kernel evaluation, deterministic
finite-grid GECP, solve-based cross approximation, diagonal-pivoted Cholesky,
conditional Lipschitz pivot certificates, numerical separated interpolation,
and sparse spectral recovery.

```python
import numpy as np

from kernelgecp import FermionicKernel, GECPConfig, gecp, evaluate_residual

kernel = FermionicKernel(cutoff=100)
result = gecp(kernel, config=GECPConfig(pivot="grid", tol=1e-10))
residual = evaluate_residual(
    kernel,
    result,
    np.linspace(0, 1, 101),
    np.linspace(-100, 100, 401),
)
print(result.rank, result.converged, result.stop_reason, abs(residual).max())
```

`cross_approximation` solves with the selected core; it never forms an explicit
inverse. Results include residual history, determinants, smallest singular
values, near-tie counts, convergence, stop reason, and actual precision.
For PSD matrices, `pivoted_cholesky` and `gecp_matrix` use the same
diagonal-preferring lexicographic tie rule; the recursive equivalence of those
pivot policies is separately proved in Lean.

## Certified pivots

`certified_pivot` performs deterministic cell subdivision from a supplied
coordinatewise Lipschitz bound. The returned bounds are conditional on that
analytic bound. Cell-budget exhaustion returns `certified=False`; it never
silently becomes success. The implementation currently uses float64 arithmetic
without directed rounding, so these objects are numerical certificates, not
formal proofs.

## Separated interpolation

`fermionic_separated_approximation` uses piecewise barycentric interpolation on
dyadic Chebyshev frequency bands. `dlr_rank_bound` and
`exp_family_rank_bound` report explicit arithmetic counts for research
comparison. They do not certify that the planned Lean uniform-error theorem has
been proved.

## Sparse spectral recovery

`sparse_representation` performs orthogonal matching pursuit and optional
bounded nonlinear frequency refinement. It returns `converged=False` when the
requested tolerance is not met. The fixed two-delta regression reaches
`1e-8` with two atoms. AAA and vector fitting are deferred extension methods,
not placeholder APIs.

## Reproducible experiments

```bash
uv run python -m kernelgecp.census \
  --config experiments/census.toml \
  --output experiments/data/gecp_census.jsonl
uv run python experiments/exact_surrogates.py
```

Runtime algorithms require no network access. The canonical census is
finite-grid float64; arbitrary-precision continuous certification remains open.
