# Production application integration contract

This records the user's accepted EB–EJ implementation scope. Progress and
canonical phase receipts are tracked in PRODUCTIZATION.md; this document
defines requirements rather than a competing completion ledger.

## Goal and boundaries

Deliver tested public interfaces, actual provider workflows, reproducible
environments, realistic datasets, matched accuracy/performance evidence, and
kernelgecp 1.2.0 for correlated materials/DMFT, finite-temperature quantum
chemistry, superconductivity, spectroscopy, and Gaussian-process covariance.
Users are computational scientists and scientific-software developers.
Support adapters plus a closed-loop DMFT integration. Preserve Green/EPW
internal solvers; downstream compression/archive/evaluation/interchange does
not imply internal solver acceleration. Upstream maintainer acceptance is
external and not required for controlled delivery.

Use the 16 GiB workstation and existing GitHub CI, bounded overnight CPU runs,
no new paid compute or GPU requirement. Default dependencies remain
NumPy/SciPy/mpmath. Optional providers are isolated and lazily imported.
Support only actual tested versions/workflows. Full ab initio regeneration,
Green/EPW internal transform replacement, off-diagonal spectral continuation,
GPU training, universal speedup/priority claims and new continuous numerical
theorems are outside this release. Published 1.1.0 artifacts stay immutable.

## Shared interfaces

- `FermionicBasis.prepare_matsubara(indices, rcond=None)` returns the cached
  SVD `LeastSquaresPlan`; `fit_matsubara(values, indices, axis=0, rcond=None)`
  returns `FermionicFit`. Preserve complex linearity and trailing axes.
- Ordinary integer n means `(2*n+1)*pi/beta`. Add observation-domain diagnostics.
  Write archive version 2; read v1 with its original time-fit interpretation.
- Optional integrations under `kernelgecp.integrations`: TRIQS tau import and
  tau/iw export; molecular Green HDF5 import; EPW Nambu/native restart
  interchange; covariance-aware MaxEnt; SciPy and CPU GPyTorch operators.
- Preserve per-fit NPZ plus versioned JSON application manifests with hashes,
  units, provider identity, native block/state ordering and parameters.
  No pickle, implicit conjugation/symmetry, static-term subtraction, tail
  repair, dense fallback or silent provider substitution.
- Invalid/unsupported input fails clearly. Numerical/resource exhaustion is
  failed or inconclusive. Training residual, controller convergence,
  independent physical error and mathematical certification remain distinct.

## EB — runtime environments

Start with TRIQS 4.0.2, CT-HYB 4.0.1, MaxEnt 4.0.1, Green 1.0.0,
PySCF 2.14.0, qe-7.5/EPW 6.0 and CPU GPyTorch. Linux x86_64 is canonical for
native solvers; macOS arm64 support requires actual provider execution.
Build/cache pinned Linux environments in CI without assuming local Docker.
Freeze source, transitive/build locks, licenses, MPI/HDF5 and thread settings
after real imports and minimal calculations pass. Resolution/build success
alone does not establish runtime acceptance.

## EC — frequency fitting and TRIQS

Implement shared frequency fitting and `fit_triqs_tau`, `to_triqs_tau`,
`to_triqs_iw`, plus ordered BlockGf recipes. Validate actual mesh coordinates,
beta, fermionic statistics, finite data, shapes, order and ownership.
Keep 0+ and beta- distinct; initial support is TRIQS 4.
Acceptance: actual-provider scalar/complex matrix pole errors <=1e-8;
independent Fourier checks; positive-only/signed meshes; cache/ownership/block
order tests; v1/v2 archives; ordinary imports without optional providers.

## ED — deterministic interacting DMFT

Half-filled paramagnetic Bethe IPT: t=1, beta=10, U=0/1/2, mixing 0.5,
200 iterations maximum. Cancel Hartree analytically and use separate Green
and dynamic-self-energy cutoffs. Compare dense TRIQS, kernelgecp, pydlr and
sparse-ir through frozen-trajectory replay and independent loops.
Acceptance: exact atom/finite-bath Lehmann and finite-pole IPT-update oracles;
sampled G/Sigma errors <=1e-6; reference refinement change <1e-7;
unmixed equation defect <1e-8 twice consecutively. Report density, causality,
moments, projection and convergence separately. IPT is an approximate solver.

## EE — CT-HYB DMFT

Run actual CT-HYB dense/kernelgecp closed loops at beta=10,U=2 using measured
G_tau -> fit -> reconstructed G_iw in the bath update and native moment
handling. Replay identical measurements through GECP/DLR/IR. Calibrate with
cycle budgets 100000/400000/1600000, then freeze the choice. Retain independent
fixed-bath replicas for uncertainty and EH. Acceptance: density, double
occupancy and selected G values agree within 1e-3 plus three combined standard
errors; observable standard errors <=3e-4; thermalization/sampling recorded.
Separate representation, Monte Carlo and self-consistency error. Same seeds
do not establish paired paths once baths diverge; budget misses are inconclusive.

## EF — molecular GF2

Actual restricted Green GF2 H2/STO-3G at 0.74/1.5 Angstrom, beta=20 Hartree^-1,
declared symmetric Lowdin basis. Independent finite-temperature Lehmann
references use PySCF integrals and FCI sector-energy checks.
`read_green_tau(input_path, results_path, run_manifest, iteration=None)`
preserves spin/orbital axes, overlap, mu, static Sigma, iteration and energies.
Require a manifest and actual convergence; reject incomplete checkpoints.
Acceptance: real native GF2 converges; same-solution compression error <=1e-6;
native-grid/physical-cutoff refinement; independent Dyson defect. GF2-vs-FCI
disagreement is model error, not representation failure. Normal CI uses frozen
attributed actual outputs; release gates rerun the solver.

## EG — EPW superconductivity

Use official coarse MgB2 epw_super outputs from qe-7.5/EPW6 for native output,
Nambu construction, Matsubara fit, compact archive and anisotropic restart
export. Preserve state ordinal, frequency order, temperature, eV units and mu
convention; require explicit offset for moving mu. Do not impersonate EPW's
specialized fermionic/bosonic/parity/singular-value IR object format.
Acceptance: native parse/write order/rounding; held-out scaled Nambu error
<=1e-6; recovered Z error <=1e-5; Delta/chi errors <=1e-6 eV. Reject unsafe
conditioning. Unmodified EPW must read exported anisotropic files and pass a
bounded restart comparison. Coarse grids demonstrate interoperability, not
converged experimental predictions.

## EH — covariance-aware spectroscopy

Actual TauMaxEnt for scalar/diagonal positive spectra: off-grid atoms, smooth
densities, independent/correlated noise and EE CT-HYB observations. Compare
raw and GECP/DLR/IR-compressed continuation. Propagate the complete projection
covariance; retain declared eigenmodes/effective rank. Interpolated points
are not independent measurements. Acceptance: positivity/normalization to
1e-6; independent held-out forward predictions at the declared noise budget
(whitened RMS <=2 for noisy fixtures); spectral-window/CDF errors, bootstrap
variation and default-model sensitivity. Existing sparse methods remain
separately labeled comparisons. No unique peaks/full uncertainty-coverage claim.

## EI — scalable covariance and GP

Add immutable column-oracle KernelCholeskyFactor with out-of-sample extension,
SciPy LinearOperator, solve-based shifted preconditioner and optional CPU
float64 GPyTorch conversion; no autograd-through-selection claim. Preserve
the dense Cholesky API. Construction retains O(n*r) state, never an n*n matrix.
Use hash-pinned attributed UCI Beijing Aotizhongxin PM2.5/meteorology data,
chronological train/calibration/test splits, training-only normalization.
Sizes 512 through 16384; ranks 32/64/128/256/512. At n<=2048, means and
standard deviations match dense SciPy normalized RMS <=1e-3/max <=1e-2.
Approximate solve residual <=1e-10; original-kernel preconditioned CG <=1e-8.
At n=16384 retain <=128 MiB numerical state and fresh-process RSS <=4 GiB.
Test identical-factor SciPy/GPyTorch, duplicates, rank deficiency, callbacks,
ties and tiny noise. No optimal geographical sensor-placement claim.

## EJ — release and common delivery gates

Publish 1.2.0 only after all mandatory real-provider domain gates pass. Ship
tested installation/compatibility guides, runnable examples, troubleshooting,
reports, source/license inventories, migration notes and final handoff.
Run ./scripts/verify.sh each phase, with fail-closed dedicated provider lanes.
Small real fixtures belong in PR CI; full scientific runs gate domain closure.
Hosted jobs <=330 minutes; checkpoint longer campaigns across jobs and run one
memory-intensive local workload at a time with fixed threads. Freeze source,
dependencies/config/seeds before final generation and separate calibration
from untouched validation. Measure construction/fitting/evaluation/full
workflow time, RSS, storage, ranks, conditioning and errors; three fresh
process repeats for performance claims. Matched accuracy precedes any speedup
claim; retain slower and failed cases. No domain closes on synthetic tests only.
After each phase update current docs, verify exact head, merge under standing
authorization, push/synchronize main and verify main. Final installed-package,
public publication, exact-main CI and canonical receipt gates are mandatory.
