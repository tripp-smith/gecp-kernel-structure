# Migration from 1.1.0 to the 1.2.0 integration candidate

Version 1.2.0 remains unreleased until all domain and publication gates pass.
Published 1.1.0 wheels, source archives, release notices and fit archives remain
immutable. Keep the old environment if reproducing an old result.

## Fitting and archives

The existing time-domain fitting interface and coefficient meaning remain.
The new `basis.prepare_matsubara(indices, rcond=None)` and
`basis.fit_matsubara(values, indices, axis=0, rcond=None)` use ordinary integer
indices with omega_n=(2n+1)pi/beta. They preserve supplied ordering, complex
linearity and trailing matrix axes. They do not infer symmetry, subtract static
terms or correct tails.

Fit diagnostics now name their observation domain. New fits use archive format
2; the loader continues interpreting format-1 fits as time-domain observations.
Loading historical diagnostics does not rerun or certify the original fit.
Retain provenance and independently validate transferred outputs. Existing
size/rank/work limits and pickle-free validation still apply.

## Optional adapters and numerical workflows

New integrations live under `kernelgecp.integrations` and import optional
providers lazily. Use the [tested compatibility table](INSTALLATION_COMPATIBILITY.md)
and exact locks. Installing the base package alone does not promise a native
solver runtime. TRIQS support starts at the tested version-4 interface.

Green molecular import requires a versioned bound run manifest and preserves
spin/orbital axes, overlap, static self-energy and chemical potential. EPW
interchange requires explicit native ordering, units and chemical-potential
convention. Its supported restart archive retains explicit lossless corrections;
the failed fit-only route stays unsupported. Neither adapter replaces its
provider's internal solver or proves a performance advantage.

MaxEnt datasets carry complete covariance, projection and source-time
coordinates. Only scalar or diagonal positive spectral functions are supported.
Initial minimizer damping is explicit; it is not a floor or an uncertainty
model. Use the accepted campaign configuration and preserve failed native trials.
No unique-peak or full bootstrap-coverage guarantee is introduced.

The dense pivoted-Cholesky API remains. Column-oracle `KernelCholeskyFactor`
adds O(nr) retained state, out-of-sample extension, SciPy operators and shifted
solve/preconditioner objects. Optional GPyTorch conversion is CPU float64 and
does not differentiate through pivot selection. Numerical residuals and
resource limits remain acceptance obligations.

## Reproduction and publication

Read [PRODUCTIZATION.md](PRODUCTIZATION.md#production-integrations-12) for
phase-specific scientific and delivered-main receipts. Use the source bound
to each experiment, rather than relabeling old data with a new environment.
A later source or security update does not create fresh blind validation.
Release payload verification, scientific eligibility, installed-package tests
and public PyPI filename/hash verification are separate gates. No universal
speedup or extension of the original Lean theorem scope is implied.
