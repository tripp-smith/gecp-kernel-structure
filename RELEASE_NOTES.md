# kernelgecp 1.2.0

This release adds optional scientific integrations and reusable numerical
interfaces to the 1.1.0 package.

- Matsubara fitting and versioned archives preserve complex matrix axes and
  support the existing time-domain interface.
- Provider-isolated workflows cover TRIQS and CT-HYB, molecular Green, EPW,
  positive scalar/diagonal MaxEnt continuation, and CPU GPyTorch conversion.
- Column-oracle pivoted Cholesky provides bounded retained state, extension,
  operators, and shifted solves with explicit numerical diagnostics.

The [canonical EH run](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/42978469706)
completed successfully at frozen source `6caf85f`. The independent
[saved-data replay](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/43135341621)
passed all five audits, with its origin pinned at `3945d503a7b2cda261a8c47fc159587890d165a6`.
The [EH receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/146#issuecomment-6569177201)
and [EJ receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/144#issuecomment-6569479579)
track scientific and release gates. See the
[compatibility guide](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/INSTALLATION_COMPATIBILITY.md)
and [migration guide](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/MIGRATION.md).

Provider locks and actual execution receipts define supported combinations.
No universal speedup, unique spectral recovery, continuous inexact-pivot rate,
or extension of the original Lean theorem scope is claimed.

---

# kernelgecp 1.1.0

This release packages the practical Green-function interface developed after
the original mathematical research closure.

- Reusable `FermionicBasis`, `FermionicFit` and `LeastSquaresPlan` support cached
  factors, complex matrix batches, physical time/ordinary Matsubara indices and
  bounded, versioned NPZ/JSON archives.
- Explicit finite-query precision certificates separate arithmetic/solve
  discrepancy from physical/model error and retain inconclusive outcomes.
- Optional Rust/MPFR 128-bit reference time evaluation passes the frozen
  streaming-workload retention gate; Python remains the default. Small reuse
  and warm cache hits show much less benefit.
- Universal Python and Linux x86_64/macOS arm64 native wheels include typing.
  Native wheels include complete notices, dependency inventory and a matching
  source asset with offline rebuild instructions. Windows uses Python.
- Shared-accuracy DLR/IR comparisons, a small noninteracting physical solver,
  independent numerical checks, preserved failures, METHODS and CITATION make
  the evidence reproducible at its recorded scope.

Research GECP defaults and the original Lean theorem statements retain their
scope. No universal fastest method, native Windows support, continuous
inexact-pivot rate or upstream maintainer endorsement is claimed.
