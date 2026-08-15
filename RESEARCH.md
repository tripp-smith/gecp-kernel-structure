# Research contracts

## Phase 0 — bootstrap

State: complete; merged in [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1).

Deliverables are the repository-scoped skills, pinned Lean and Python
toolchains, package roots, verification scripts, CI, and progress metadata.
This phase makes no mathematical claim.

## Milestone A — exact selected core closure

State: complete; merged in [PR #4](https://github.com/tripp-smith/gecp-kernel-structure/pull/4).

The phase closes the remaining gap between the dependent `GECP.Run` object and
the determinant/nonsingularity API. The exact public contract is:

- `Run.SelectedIndex` enumerates the successful pivots in their run order;
- `Run.selectedCore` evaluates the original kernel at those selected rows and
  columns;
- `Run.finSelectedCore` reindexes that same matrix by `Fin run.pivots.length`;
- `gecp_core_det_eq_prod_pivots` proves that the determinant of the finite
  selected core is the ordered product of the actual residual pivots carried
  by `run`;
- `gecp_core_nonsingular` proves the selected core determinant is nonzero.

The proof must derive its block elimination step from `residualUpdate` and the
nonzero pivot witness already stored in `Run.step`. It may not assume external
LDU data. This phase does not add a GECP convergence or pivot-selection claim.

## Milestone B — recursive PSD/Cholesky equivalence

State: complete; merged in [PR #5](https://github.com/tripp-smith/gecp-kernel-structure/pull/5).

The phase strengthens the existing one-step diagonal-domination result to the
following exact contract for finite real PSD matrices:

- `diagonalResidual` is the symmetric same-domain Schur update at a positive
  diagonal pivot;
- `posDef_gecp_residual_posSemidefinite` proves that update remains PSD;
- `IsDiagonalPivot` and `IsCompletePivot` state the pivoted-Cholesky and GECP
  maximizer predicates;
- `PivotedCholeskyTrace` and `GECPTrace` record the same ordered residual
  updates with their respective pivot predicates;
- `gecp_eq_pivotedCholesky_of_posDef` proves the two recursive trace relations
  equivalent for every pivot list.

The complete-pivot relation deliberately uses a diagonal pivot location and
the repository's diagonal-preferring canonical tie policy. The theorem does
not claim that every arbitrary off-diagonal tie choice produces the same run,
and it does not claim the stronger published Lipschitz convergence rate.

## Milestone C — high-precision census and certified pivot closure

State: complete; merged in [PR #6](https://github.com/tripp-smith/gecp-kernel-structure/pull/6).

This closure phase delivers:

- a genuine `mpmath` finite-grid GECP backend that executes at
  `GECPConfig.precision_bits` and preserves canonical quantities as decimal
  strings;
- an outward-rounded interval branch-and-bound pivot certificate for the
  fermionic residual on the continuous rectangle;
- adaptive GECP results that retain every pivot certificate and never report
  certification after a cell-budget exhaustion;
- a canonical endpoint-resolved finite-grid census whose recorded 128-bit
  precision and grid algorithm match the backend actually used, with
  repeat-run byte equality;
- the kernel derivative bounds and grid-to-continuous theorem/API constants
  needed to interpret approximate pivots.

The canonical census uses finite grids because a complete continuous adaptive
census is substantially more expensive. The supported adaptive mode retains
interval certificates and reports budget exhaustion honestly. Neither the
census nor those certificates are substituted for Lean proofs of the
Milestone D or E rate theorems.

## Milestone D — constructive separated approximation

State: complete; merged in [PR #7](https://github.com/tripp-smith/gecp-kernel-structure/pull/7).

Delivered:

- `expNegTaylor_error_le_next` proves the exact Lagrange remainder for the
  truncated negative exponential;
- `pow_div_factorial_eight_mul_le` proves that `8p` terms give error at most
  `2⁻ᵖ` whenever the local product is at most `2p`;
- `expFamily_separatedApprox` constructs `8p(s+1)` separated terms and proves
  uniform error `2⁻ᵖ` on `[0,1] × [0,2ˢ]`;
- `fermionicKernel_separatedApprox` uses the exact reflection identity and the
  positive denominator to construct `16p(s+1)` terms with the same uniform
  error on `[0,1] × [-2ˢ,2ˢ]`;
- `fermionic_dyadic_taylor_approximation` independently implements the same
  band/time-cutoff rules at arbitrary precision and tests band boundaries and
  cutoff `10⁶`.

The primary formal construction is a dyadic truncated-Taylor construction,
not the selected-exponential Chebyshev construction in Gimbutas–Marshall–
Rokhlin. It proves the same required `O(log Λ log(1/ε))` separated-rank scale
with conservative constants and has a substantially smaller formal
interpolation dependency surface. The published selected-exponential variant
remains a possible strengthening, not a prerequisite for the delivered
low-rank theorem. Nothing in this milestone proves that GECP selects these
terms or converges at this rate.

## Milestone E — structural GECP

State: complete; merged in [PR #8](https://github.com/tripp-smith/gecp-kernel-structure/pull/8).

Exact checks completed for all requested geometric surrogates. There were no
zero minors and no sign mismatches. However, the exact pivot sequences depend
on `q`; for size eight the three requested ratios produce different row and
column orders.

The rigorous closure contract is implemented in two parts:

- `signRegular_not_sufficient_for_parameterIndependent_pivots` gives a
  minimized exact `2 × 2` counterexample: two kernels have the same strict
  one- and two-minor signs but disjoint complete-pivot sets. Minor-sign data
  alone therefore cannot determine a universal GECP pivot path.
- `CrossRatioControl` bounds the post-update cross numerator relative to a
  uniform residual bound. `residualUpdate_le_of_crossRatioControl` proves
  one-step contraction, and
  `gecp_error_le_geometric_of_crossRatioControl` proves the resulting `θⁿ`
  residual bound for an exact update sequence.

This is the SPEC-approved “certified obstruction plus refined sufficient
condition” outcome. It does not establish that fermionic GECP satisfies a
uniform `CrossRatioControl` with `θ < 1`; therefore it is not the target rate,
a weaker fermionic rate, or a solution of Simons Problem 4.2.

The endpoint-resolved 128-bit finite-grid census exhibits rank growth across
all seven cutoffs, but remains numerical evidence rather than a continuous
rate theorem. The next optional strengthening is:

1. formalize sign regularity for geometric/exponential collocation minors;
2. derive `CrossRatioControl` from a more intrinsic determinant or
   near-max-volume condition;
3. test whether the fermionic residual satisfies the hypothesis using
   outward-rounded or exact bounds.

## Phase G — fermionic GECP rate research

State: verified obstruction-plus-base-case outcome on
`codex/phase-g-gecp-rate`; Conjecture G1 remains open.

The target is Conjecture G1: an exact continuous GECP residual bound of the
form

\[
e_k(\Lambda) \le C\exp\!\left(-c k/\log(1+\Lambda)\right)
\]

with universal positive constants. A cutoff-independent contraction at every
single pivot would be stronger than this target and is not a viable contract:
after selecting the complete pivot `(0, Λ)`, the reflected-corner residual is
exactly `1 - exp(-Λ)`, while the initial maximum is
`1 / (1 + exp(-Λ))`. Their ratio is `1 - exp(-2Λ)` and therefore approaches
one. Phase G consequently tests contraction over blocks of
`O(log(1 + Λ))` pivots, or an equivalent determinant/near-volume condition.

Frozen public Lean targets:

- `fermionicKernel_firstPivot_reflected_residual`: the exact first-update
  formula above, with no numerical assumptions;
- `fermionicKernel_firstPivot_abs_le_reflectedCorner`: the reflected corner is
  the actual second complete pivot for every positive cutoff;
- `fermionicKernel_no_uniform_firstStep_contraction`: a quantified obstruction
  to any cutoff-independent one-step factor below one;
- `gecp_error_le_product_of_crossRatioControl`: variable-factor exact residual
  control, so a block product rather than every individual factor may contract;
- `fermionicKernel_twoCornerResidual_le_half_initial`: the actual two-step
  continuous GECP residual contracts by one half on `0 < Λ ≤ 1`;
- `fermionicKernel_gecp_error_le_exp`: added only if a fermionic-specific block,
  determinant, or structural lemma genuinely proves Conjecture G1.

The research pass has proved the infrastructure targets and the following
fermionic-specific refinements:

- `fermionicKernel_le_cutoffCorner` proves that `(0, Λ)` and its reflected
  partner are genuine complete pivots on the cutoff rectangle;
- `fermionicKernel_firstPivot_abs_le_reflectedCorner` strengthens this to the
  recursive statement: after `(0,Λ)`, `(1,-Λ)` is the actual complete pivot
  for every `Λ>0`;
- `fermionicKernel_firstPivot_reflected_ratio` proves the exact ratio
  `1 - exp(-2Λ)`, and
  `fermionicKernel_no_uniform_firstStep_contraction` constructs a positive
  cutoff defeating every proposed fixed one-step factor below one;
- `gecp_error_le_dyadic_of_crossRatioProduct` shows that a cumulative factor
  bound at ranks `block * p` gives error `2⁻ᵖ` times the initial bound;
- `PivotCrossProductSignCoherent` captures the selected-cross sign pattern suggested by
  totally positive exponential kernels.
  `crossRatioControl_one_of_signCoherent` and
  `residualUpdate_le_of_signCoherent` prove that exact complete pivoting is
  nonexpansive whenever the current residual has this property.
- `centeredExpSecant_error_le_quarter` and
  `twoCornerApproximation_error_le_eighth` construct and bound an explicit
  rank-two comparison space on `0<Λ≤1`;
- `cornerWeights_nonnegative_sum_le_one` proves stability of its endpoint
  interpolation operator, while
  `fermionicKernel_twoCornerResidual_eq_sub_interpolate` identifies that error
  with the actual recursive residual;
- `fermionicKernel_twoCornerResidual_le_half_initial` closes the continuous
  cutoff-one block exactly. The interval certificate below is now an
  independent numerical check of a formal theorem, not the source of the
  claim.

The canonical finite-grid data suggests the sharper block contract

\[
  \|R_{n+2(s+1)}\|_\infty \le \tfrac12\|R_n\|_\infty,
  \qquad \Lambda\le 2^s.
\]

At the strictest stored tolerance, every complete block satisfies this test.
The maximum observed complete-block ratios for cutoffs
`1, 10, ..., 10^6` are respectively approximately
`7.74e-2`, `4.97e-6`, `2.04e-4`, `1.13e-3`, `1.00e-3`,
`2.26e-3`, and `2.48e-3`. The interval engine independently certifies the
continuous two-step ratio below `0.078` at cutoff one along its returned
certified approximate-pivot trajectory.

Two tempting stronger routes were rejected rather than promoted:

- per-step cutoff-uniform contraction fails by the exact first-step theorem;
- after reflecting and scaling frequency, the kernel is diagonally similar to
  a symmetric positive-definite exponential kernel, but stored pivots move far
  off the transformed diagonal for nontrivial cutoffs. The existing
  pivoted-Cholesky theorem therefore does not apply directly.

The current missing scaling lemma is now precise: prove that every exact fermionic
GECP residual is `PivotCrossProductSignCoherent` at its selected pivot (or an equally strong
nonexpansiveness invariant), and prove a cutoff-uniform half reduction within
`C(s+1)` subsequent pivots by transporting the proved `Λ≤1` base mechanism
through dyadic frequency layers. The former is suggested by strict sign regularity
of exponential collocation determinants; the latter still needs a dyadic
pivot-localization, determinant, or near-volume argument. Neither numerical
observation is being used as that proof.

Python research targets:

- record every residual ratio and its cutoff-normalized block product using
  arbitrary precision;
- independently verify the exact first-step formula;
- use outward-rounded interval residual bounds for any observation promoted to
  a fermionic structural lemma;
- minimize and retain any failure of the proposed block condition before the
  theorem contract is revised.

Non-claims for this verified partial result:

- finite-grid collapse of residual curves is not a continuous theorem;
- a factor fitted from the seven canonical cutoffs is not cutoff-uniform;
- the generic variable-factor theorem does not assert that fermionic GECP
  satisfies its hypotheses;
- The verified Phase G result is a rigorous obstruction plus a strict base
  block, not completion of Conjecture G1; Simons Problem 4.2 is not solved
  unless the dyadic scaling/dependence gap is closed rigorously.

## Numerical evidence identifiers

- `experiments/data/exact_surrogates.json`:
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`
- `experiments/data/gecp_census.jsonl`:
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`

## Run provenance

The model-identification boundary, collaboration mode, elapsed-time and token
snapshot, and explicitly non-billing cost estimate for the implementation run
are recorded in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). These
operational metadata do not change any proved, observed, conjectured, or
not-claimed research status above.

## Phase S — realistic synthetic applications

State: complete; implementation commit `6404311`, with canonical evidence and
final documentation delivered directly to `main` as requested by the user.

This post-v1 phase applies the delivered Python algorithms to stylized but
scientifically recognizable workloads. It depends on the completed fermionic
kernel, grid GECP, cross approximation, sparse recovery, pivoted-Cholesky, and
Green-error-transfer work. It adds no Lean theorem and does not reopen any
completed mathematical milestone.

The exact contract is:

- add a public, typed `SyntheticApplicationConfig` and
  `run_synthetic_application_suite` API;
- compress normalized Hubbard-like three-peak and gapped two-band spectral
  densities with one universal fermionic GECP basis, then verify the discrete
  `L∞`-to-`L¹` Green-error transfer on held-out time/frequency grids;
- recover a multi-line quasiparticle/satellite spectrum from a known
  transition library and represent a second spectrum through a dense blind
  scan of deterministic noisy imaginary-time data, reporting dictionary
  strategy, held-out errors, and stop reasons;
- use GECP and canonical pivoted Cholesky as matching landmark selectors for a
  clustered synthetic spatial covariance matrix;
- commit configuration-addressed JSON and a summary plot, and test the public
  API, deterministic serialization, transfer inequality, recovery behavior,
  and PSD pivot agreement.

Non-claims:

- the spectral fixtures are synthetic and are not fits to a named material,
  impurity calculation, experiment, or quantum Monte Carlo dataset;
- the noisy sparse example is not an analytic-continuation uniqueness,
  uncertainty-quantification, or minimax-stability theorem;
- finite quadrature verifies a discrete instance of the formal transfer
  theorem and is not a replacement for its continuous Lean proof;
- the covariance example is a control application of the PSD baseline, not a
  formal theorem about optimal sensor placement;
- none of the new evidence proves the target continuous fermionic GECP rate.

Delivered observations:

- the common 12-pivot fermionic basis has validation kernel error
  `5.404231354739705e-9`; the Hubbard-like and gapped spectra have Green errors
  `4.162832301091157e-10` and `2.909588125987739e-10`, respectively, below the
  discrete transfer bound;
- the known-transition dictionary recovers four atoms to machine precision,
  while the noisy dense scan returns eight effective atoms with
  `8.559113029438237e-9` noiseless held-out error and makes no identification
  claim;
- GECP and Cholesky agree on the 31 selected landmarks for the 72-site
  covariance, whose cross error is `7.782740553130552e-7` in maximum norm;
- `experiments/data/synthetic_applications.json` is byte-reproducible for
  implementation commit `6404311`, configuration hash
  `b41f4c17e76ccb108f2999b0e07af1003eb011bf8d96f2c7456e92b9d36f8d35`,
  and SHA-256
  `f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2`;
- the reviewed summary plot has SHA-256
  `0529263affd2fe10773f71784f06007650bff833374fbbebbb16f2507b2ab169`.

Scientific motivation comes from the DLR effective-delta representation of
imaginary-time Green functions, the quasiparticle/spectral-weight structure
studied in DMFT, and the noisy analytic-continuation setting. The experiments
remain stylized benchmarks rather than replications of those prior works.
