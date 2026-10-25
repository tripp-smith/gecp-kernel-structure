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
