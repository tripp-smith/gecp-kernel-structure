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

State: complete as a bounded obstruction-plus-base-case research phase; merged
in [PR #11](https://github.com/tripp-smith/gecp-kernel-structure/pull/11).
Conjecture G1 remains open.

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

## Phase H — verification integrity refresh

State: verified locally with green CI; merged
[PR #12](https://github.com/tripp-smith/gecp-kernel-structure/pull/12).

This bounded maintenance phase restores the repository's root definition of
done before the next theorem phase. It depends on the merged Phase G library
and test suite but does not reopen or alter any mathematical result.

The exact contract is:

- make the local root test invocation robust to editable-install entry-point
  state by launching pytest through the synchronized project interpreter;
- keep the local root command and GitHub Actions invocation identical;
- refresh the development dependency lock so `pip-audit` has no known
  vulnerability finding in the resolved environment;
- verify all Lean builds, public-theorem axioms, proof-placeholder rejection,
  formatting, typing, Python tests, and dependency audit through
  `./scripts/verify.sh` from a fresh synchronized environment.

Independent checks:

- invoke the Python suite through the project interpreter and confirm the
  existing 45 tests pass;
- inspect the reverse dependency path for any audited package upgrade;
- run the complete root command after the lock refresh rather than treating
  focused checks as completion.

Non-claims:

- this phase proves no new Lean theorem and adds no numerical evidence;
- it does not strengthen Phase G, prove `PivotCrossProductSignCoherent` for
  fermionic residuals, establish dyadic block contraction, or resolve
  Conjecture G1;
- a clean dependency audit is repository-maintenance evidence, not evidence
  for any mathematical or numerical claim.

Delivered verification evidence:

- `scripts/verify.sh` and GitHub Actions launch pytest as
  `uv run python -m pytest`, binding collection to the synchronized project
  interpreter even when an older editable-install console script is stale;
- the development lock resolves `urllib3` 2.8.0 through
  `pip-audit -> requests -> urllib3`, replacing the audited 2.7.0 lock;
- `./scripts/verify.sh` passes both in the maintained environment and in a
  newly created `uv sync --all-extras --frozen` environment: the Lean build,
  axiom audit, placeholder scan, Ruff, formatting, mypy, all 45 tests, and
  `pip-audit` are green.

## Phase I — determinantal sign-coherence bridge

State: complete; merged PR #13 with green CI.

This phase attacks the selected-pivot nonexpansiveness route left open by
Phase G. Mathlib supplies determinant, reindexing, and Schur-complement
infrastructure but no existing total-positivity or strict-sign-regularity
abstraction. The phase therefore freezes the finite algebraic bridge needed
before any fermionic analytic sign theorem is attempted.

Frozen public Lean targets:

- `GECP.Run.borderedCore`: the original kernel evaluated on a successful
  run's selected rows and columns with one additional row and column;
- `GECP.Run.borderedCore_det_eq_selectedCore_det_mul_finalResidual`: the exact
  determinant identity
  `det(border(run,x,y)) = det(core(run)) * run.finalResidual x y`;
- `BorderedMinorSignCoherent`: the four-bordered-minor product condition at a
  proposed pivot;
- `pivotCrossProductSignCoherent_iff_borderedMinorSignCoherent`: equivalence
  between the residual selected-cross condition and the determinant condition
  for every successful run.

Independent checks:

- retain exact `Fraction` checking of every selected cross on the geometric
  surrogates, and return a minimized step/coordinate witness if the condition
  fails;
- use arbitrary-precision fermionic finite-grid residuals only to test the
  theorem direction and choose the next analytic claim, never as proof;
- add every public structural theorem to the axiom audit and run the full root
  verification command.

Observed before implementation:

- exact geometric surrogates pass selected-pivot sign coherence for the
  canonical sizes 2--8 and `q in {1/2,2/3,3/4}`;
- a 70-decimal-digit scan through 24 pivots at cutoffs `1,2,10,100` found no
  substantive negative cross product; normalized negative roundoff at exact
  zero rows or columns was between approximately `1e-89` and `1e-131`.

Implemented result:

- the recursive bordered index keeps pivot order and appends the requested
  border, allowing an explicit unit-lower block elimination proof at each run
  step;
- `borderedCore_det_eq_selectedCore_det_mul_finalResidual` identifies every
  bordered determinant with the selected-core determinant times the exact
  final residual;
- selected-core nonsingularity makes the fourth power of its determinant
  strictly positive, so
  `pivotCrossProductSignCoherent_iff_borderedMinorSignCoherent` follows by
  substituting the four determinant identities;
- exact regression checks verify the determinant identity along a complete
  geometric-surrogate run, all 21 canonical surrogate paths pass the sign
  criterion, and a generic `2 x 2` matrix supplies the exact negative witness
  `(step,pivot,point,product) = (0,(0,0),(1,1),-48)`.
- `./scripts/verify.sh` passes with the full Lean build, axiom audit,
  placeholder scan, Ruff, formatting, mypy, all 47 tests, and no known
  dependency vulnerabilities; the new public theorems use only the audited
  Lean defaults `propext`, `Classical.choice`, and `Quot.sound`.

Next decision:

- seek a determinant-sign factorization for ordered exponential collocation
  matrices strong enough to prove the four-minor product nonnegative under
  the row and column order induced by GECP; do not attempt Conjecture G1 until
  this condition, or a different nonexpansiveness invariant, is proved.

Non-claims:

- the bordered-minor criterion is not itself proof that fermionic residuals
  satisfy it;
- finite exact surrogates and high-precision grids do not establish a
  continuous-domain invariant;
- this phase does not prove strict contraction, dyadic localization,
  Conjecture G1, or the Simons Problem 4.2 rate.

## Phase J — strict-sign-regularity transfer

State: complete; green CI; merged
[PR #14](https://github.com/tripp-smith/gecp-kernel-structure/pull/14).

Phase I reduced selected-residual sign coherence to four original-kernel
bordered minors. The remaining algebraic gap is to connect that condition to
the ordered-minor language in the SPEC before attempting the all-orders
analytic theorem for the exponential kernel. A mathlib search found tuple
sorting, determinant permutation, and row/column scaling lemmas, but no
total-positivity, Chebyshev-system, or strict-sign-regularity abstraction.

Frozen public Lean targets:

- `StrictSignRegularAtOrder`: every square minor on strictly increasing row
  and column tuples has a prescribed strict determinant sign;
- `StrictSignRegular`: the corresponding all-orders signature sequence;
- `strictSignRegular_borderedMinorSignCoherent`: all-orders strict sign
  regularity implies the Phase I four-bordered-minor condition for every
  successful run and proposed pivot;
- `strictSignRegular_pivotCrossProductSignCoherent`: the resulting exact
  selected-residual sign coherence;
- `strictSignRegularAtOrder_columnScale`: multiplication by a positive column
  weight preserves the prescribed ordered-minor signs;
- `expKernel_strictSignRegularAtOrder_two` and
  `fermionicKernel_strictSignRegularAtOrder_two`: the first nontrivial
  exponential and fermionic minor-sign cases.

Implementation route:

- sort arbitrary injective finite row and column tuples, track both
  determinant permutation signs, and prove that the four border determinants
  cancel every orientation sign in pairs;
- prove selected and bordered tuples with a repeated node have zero
  determinant, so the transfer theorem covers duplicate proposed points
  without adding distinctness assumptions;
- express the fermionic kernel as the exponential interaction times its
  positive frequency-only weight and use determinant column scaling;
- independently enumerate exact geometric-surrogate minors under row and
  column permutations to verify the orientation convention.

Implemented result:

- the frozen definitions and transfer theorems above are implemented and
  re-exported by `GECPKernelStructure.Theorems`;
- repeated proposed rows or columns are handled by zero-determinant lemmas,
  rather than excluded by an extra run hypothesis;
- `fermionicKernel_strictSignRegularAtOrder_of_expKernel` transfers any fixed
  exponential order, while
  `fermionicKernel_strictSignRegular_of_expKernel` transfers the all-orders
  statement;
- orders one and two are proved directly for `expKernel` and transferred to
  `fermionicKernel`; the all-orders premise remains open;
- exhaustive exact permutation checks for the `q = 2/3` geometric surrogate
  at sizes two through four independently confirm the orientation convention.

Next analytic obligation:

```lean
StrictSignRegular expKernel expKernelSignature
```

Johnson and Richards identify `exp(x*y)` as strictly totally positive of every
order and give the relevant Vandermonde--Schur expansion in equation (3.9) of
[Hyperdeterminantal Total Positivity](https://arxiv.org/abs/2412.03000). A
Lean-oriented reconstruction can avoid importing that theory wholesale:

1. shift each finite row and column tuple into the nonnegative reals; the
   shift changes `exp(x*y)` only by positive row and column factors;
2. expand a truncated exponential kernel and apply finite Cauchy--Binet;
3. prove nonnegativity of the resulting generalized Vandermonde products and
   retain the strictly positive contribution from exponents `0, ..., n-1`;
4. pass to the exponential-series limit using continuity of the finite
   determinant;
5. reverse the frequency tuple to convert `exp(t*y)` positivity into the
   signature `(-1)^(n.choose 2)` for `exp(-t*omega)`.

The current mathlib search found ordinary Vandermonde determinant machinery
but no generalized-Vandermonde positivity or total-positivity library. That
missing lemma family, not the GECP transfer, is now the exact obstruction.

Verification:

- `lake build GECPKernelStructure.Theorems` passes;
- the public axiom audit reports only Lean defaults;
- the focused surrogate suite passes all seven tests;
- `./scripts/verify.sh` passes all repository checks with 48 Python tests and
  no known dependency vulnerabilities.

Non-claims:

- this phase does not prove strict sign regularity of `exp(-t*omega)` or the
  fermionic kernel at arbitrary order;
- order-two minor signs do not imply the all-orders hypothesis used by the
  transfer theorem;
- the phase does not establish strict contraction, dyadic localization,
  Conjecture G1, or the Simons Problem 4.2 rate.

## Phase K — exponential total-positivity scaffolding

State: complete; green CI; merged
[PR #15](https://github.com/tripp-smith/gecp-kernel-structure/pull/15).

Phase J reduced the remaining analytic work to
`StrictSignRegular expKernel expKernelSignature`. A fresh mathlib search found
the ordinary Vandermonde determinant formula and general matrix determinant
machinery, but no total-positivity, Chebyshev-system, Schur-positivity, or
generalized-Vandermonde positivity theorem. The next bounded step is therefore
to formalize the strictly positive principal term in the exponential
Cauchy--Binet expansion before attempting the full nonnegative-term sum.

Frozen public Lean targets:

- `strictSignRegularAtOrder_rowScale`: positive row scaling preserves an
  ordered-minor signature, complementing the Phase J column theorem;
- `vandermonde_det_pos_of_strictMono`: an increasing real tuple has strictly
  positive ordinary Vandermonde determinant;
- `expTaylorPrincipalMatrix`: the square truncation using exponents
  `0, ..., n-1` and weights `1/k!`;
- `expTaylorPrincipalMatrix_det_pos`: its determinant is strictly positive for
  increasing row and column tuples;
- an exact rational regression comparing the determinant with the product of
  both Vandermonde determinants and reciprocal factorial weights.

Decision rule:

- if these targets verify, close the bounded phase and make generalized
  Vandermonde nonnegativity the next exact theorem;
- do not claim the exponential determinant positive merely because one
  Cauchy--Binet term is positive: every remaining term must first be proved
  nonnegative and the finite-to-infinite limit justified.

Implemented result:

- `strictSignRegularAtOrder_rowScale` is implemented alongside the Phase J
  positive-column-scaling theorem;
- `vandermonde_det_pos_of_strictMono` specializes mathlib's determinant formula
  to strict positivity on increasing real tuples;
- `expTaylorPrincipalMatrix_apply` identifies the matrix product with the
  first `n` terms of the `exp(x*y)` series;
- `expTaylorPrincipalMatrix_det_pos` proves its determinant strictly positive
  by factoring it into two positive Vandermonde determinants and the positive
  product of reciprocal factorials;
- exact rational checks through size six reproduce that determinant product.

Next analytic obligation:

For nonnegative strictly increasing nodes `x` and a strictly increasing natural
exponent tuple `m`, prove

```text
0 ≤ det (x_i ^ m_j),
```

with strict positivity under the appropriate nondegeneracy assumptions. This
is the generalized Vandermonde/Schur-positivity lemma needed to show that every
finite Cauchy--Binet term is nonnegative. After that, the already-positive
principal term supplies strictness and only the exponential-series determinant
limit remains.

Verification:

- the public theorem build and axiom audit pass;
- all eight focused surrogate tests pass;
- `./scripts/verify.sh` passes all repository checks with 49 Python tests and
  no known dependency vulnerabilities.

Non-claims:

- this phase does not prove arbitrary-order strict total positivity of
  `exp(x*y)` or strict sign regularity of `exp(-t*omega)`;
- it does not prove the generalized Vandermonde sign theorem or the infinite
  exponential determinant limit;
- it does not establish fermionic GECP nonexpansiveness, dyadic localization,
  Conjecture G1, or the Simons Problem 4.2 rate.

## Phase L — generalized Vandermonde positivity

State: complete with green CI in merged PR #16.

Phase K proved that the principal `0, ..., n-1` exponential Taylor term is
strictly positive. To make every finite Cauchy--Binet term nonnegative, the
next exact target is the generalized Vandermonde determinant on positive
ordered nodes and increasing natural exponents. Mathlib has Descartes' rule of
signs, polynomial support/leading-term induction, ordinary Vandermonde
positivity, determinant continuity, and the intermediate value theorem, but no
assembled generalized result.

Frozen proof chain:

- `signVariations_lt_card_support`: a nonzero real polynomial has fewer sign
  variations than nonzero coefficients;
- `positiveRoots_lt_card_support`: Descartes' rule bounds distinct positive
  roots by support cardinality minus one;
- `generalizedVandermonde_det_ne_zero`: distinct positive nodes and distinct
  natural exponents give a nonsingular power matrix;
- `geometric_generalizedVandermonde_det_pos`: for nodes `2^i`, transpose turns
  the power matrix into an ordinary Vandermonde on `2^m_j`;
- `generalizedVandermonde_det_pos`: linearly homotope arbitrary positive
  increasing nodes to the geometric nodes, use nonsingularity along the path,
  and rule out a sign change by the intermediate value theorem.

Implementation outcome:

- the full frozen proof chain is implemented with the stated assumptions;
- the root-count theorem counts positive roots with multiplicity, so it is
  stronger than the distinct-root bound needed for nonsingularity;
- an exact rational regression checks all increasing exponent choices from
  `0, ..., n+2` at matrix sizes one through five;
- no failed hypothesis or weakened theorem was required.

Verification:

- the focused Lean build and public axiom audit pass, with only Lean's default
  axioms reported;
- the exact generalized-Vandermonde regression passes;
- `./scripts/verify.sh` passes all repository checks with 50 Python tests and
  no known dependency vulnerabilities.

Assumptions and intended strength:

- nodes are strictly positive and strictly increasing;
- exponents are strictly increasing natural numbers;
- the target is strict determinant positivity at arbitrary finite order;
- a later boundary argument is still required for merely nonnegative nodes,
  which can occur in an exponential Cauchy--Binet expansion after shifting.

Decision rule:

- prove the whole chain if the root-count and homotopy interfaces compose;
- after three failures at the same interface, retain the strongest verified
  prefix and name the exact missing lemma rather than weakening assumptions
  until they contain the conclusion.

Non-claims:

- this phase does not yet prove nonnegativity when a node equals zero;
- it does not prove the finite Cauchy--Binet formula or exponential-series
  determinant limit;
- it does not prove all-orders strict sign regularity of the exponential or
  fermionic kernel, dyadic localization, Conjecture G1, or Problem 4.2.

## Phase M — generalized Vandermonde boundary

State: complete with green CI in merged PR #17.

Phase L proves strict determinant positivity when every ordered node is
positive. The finite exponential Cauchy--Binet expansion also produces the
boundary case in which the first shifted node is zero. The exact target is:

- `generalizedVandermonde_det_nonneg`: nonnegative strictly increasing real
  nodes and strictly increasing natural exponents give a nonnegative power
  determinant.

Proof contract:

- shift every node by `1/(k+1)`, preserving strict order and making all nodes
  positive;
- apply `generalizedVandermonde_det_pos` at every positive shift;
- use coordinatewise convergence, determinant continuity, and closedness of
  the nonnegative ray to pass to the zero-shift limit;
- add the public theorem to the axiom audit and retain an exact zero-node
  regression, including a genuinely zero determinant when every exponent is
  positive.

Implementation outcome:

- `generalizedVandermonde_det_nonneg` is implemented with the frozen
  assumptions and proof route;
- exact rational zero-node matrices of sizes one through five are nonnegative
  over every increasing exponent choice from `0, ..., n+2`; they are strictly
  positive exactly when the least exponent is zero.

Verification:

- the focused Lean build and public axiom audit pass, with only Lean's default
  axioms reported;
- the exact zero-node regression passes;
- `./scripts/verify.sh` passes all repository checks with 51 Python tests and
  no known dependency vulnerabilities.

Non-claims:

- this phase does not yet assemble the finite Cauchy--Binet determinant
  identity or its infinite-series limit;
- it does not prove all-orders exponential or fermionic strict sign
  regularity, dyadic localization, Conjecture G1, or Problem 4.2.

## Phase N — finite exponential Cauchy--Binet

State: complete with green CI in merged
[PR #18](https://github.com/tripp-smith/gecp-kernel-structure/pull/18).

Phases K--M now provide exactly the sign information for every generalized
Vandermonde factor in a finite exponential Taylor expansion. Mathlib proves
the square `Matrix.det_mul` identity and exposes determinant multilinearity,
finite subset ordering, and submatrices, but the current dependency contains
no assembled rectangular Cauchy--Binet theorem.

Frozen targets:

- a reusable rectangular Cauchy--Binet identity, indexed by cardinality-`n`
  subsets of the intermediate finite ordered type;
- `expTaylorMatrix`: the square matrix whose entries truncate
  `exp(row_i * column_j)` to a prescribed finite term count;
- `expTaylorMatrix_apply`: the entrywise finite Taylor-sum formula;
- `expTaylorMatrix_det_pos`: strict determinant positivity when the term count
  is at least the matrix order and both node tuples are nonnegative and
  strictly increasing.

Proof contract:

- expand the rectangular product determinant and eliminate noninjective index
  maps using the determinant sign-reversing involution underlying mathlib's
  square `det_mul` proof;
- group injective maps by their range and the canonical increasing enumeration
  `Finset.orderEmbOfFin`, producing the two ordered minor determinants;
- apply `generalizedVandermonde_det_nonneg` to every Cauchy--Binet term and
  `expTaylorPrincipalMatrix_det_pos` to the principal subset for strictness;
- if the grouping interface fails three times, retain the verified injective-
  map expansion and record canonical-range grouping as the exact obstruction.

Implementation outcome:

- `Matrix.det_mul_rect_eq_sum_injective` expands a rectangular product
  determinant over injective intermediate tuples after a sign-reversing
  involution cancels every noninjective tuple;
- `Matrix.injectiveFunctionEquivStrictMonoPerm` canonically decomposes an
  injective finite tuple into its increasing rearrangement and a permutation;
- `Matrix.det_mul_rect` groups the injective expansion into the ordered
  Cauchy--Binet product of the two corresponding minors;
- `expTaylorLeft`, `expTaylorRight`, `expTaylorMatrix`, and
  `expTaylorMatrix_apply` expose the finite exponential feature product and
  its entrywise Taylor sum;
- `expTaylorMatrix_det_pos` applies generalized Vandermonde nonnegativity to
  every ordered Cauchy--Binet term and proves the principal term positive when
  the term count is at least the matrix order;
- exact rational matrices of sizes one through five, with zero boundary nodes
  and term counts from `n` through `n+3`, independently have positive
  determinants.

Verification:

- the focused Lean build and public axiom audit pass, with only Lean's default
  axioms reported;
- the focused exact surrogate suite passes all eleven tests;
- `./scripts/verify.sh` passes all repository checks with 52 Python tests and
  no known dependency vulnerabilities.

Next analytic obligation:

Prove strict positivity of the full matrix `exp(rows_i * columns_j)` by
entrywise Taylor convergence and determinant continuity. Merely taking the
limit of the positive finite determinants is insufficient: the proof must
retain the fixed positive principal Cauchy--Binet contribution as a uniform
lower bound. Positive row and column factors can then undo shifts to
nonnegative nodes, after which frequency reversal gives the prescribed
`expKernelSignature`.

Non-claims:

- this phase does not pass from finite Taylor matrices to the exponential
  kernel determinant;
- it does not yet prove all-orders exponential or fermionic strict sign
  regularity, dyadic localization, Conjecture G1, or Problem 4.2.

## Phase O — exponential determinant limit

State: complete with green CI in merged
[PR #19](https://github.com/tripp-smith/gecp-kernel-structure/pull/19).

Phase N proves every sufficiently long finite exponential Taylor determinant
positive, but strict positivity is not closed under limits. The next bounded
phase therefore retains the already-proved principal Cauchy--Binet summand as
a fixed positive lower bound while the remaining nonnegative summands grow.

Frozen public Lean targets:

- `expTaylorMatrix_det_ge_principal`: every truncation with at least `n` terms
  dominates the determinant of the principal `0, ..., n-1` feature block;
- `expTaylorMatrix_tendsto_expMatrix`: the finite Taylor matrices converge
  entrywise to the matrix `(i,j) ↦ exp(rows_i * columns_j)`;
- `expMatrix_det_pos`: that full exponential matrix has positive determinant
  for nonnegative strictly increasing row and column tuples.

Proof contract:

- identify the principal ordered tuple inside `Matrix.det_mul_rect` and bound
  its positive summand by the nonnegative finite sum;
- use `NormedSpace.expSeries_div_hasSum_exp` and finite-sum convergence for
  every matrix entry, then determinant continuity;
- pass the fixed lower bound to the limit with closed-order reasoning and use
  `expTaylorPrincipalMatrix_det_pos` for strictness;
- if the series-to-`Fin` indexing interface fails three times, isolate and
  prove the entrywise range-sum convergence lemma rather than weakening the
  determinant claim.

Implementation outcome:

- `expTaylorMatrix_det_ge_principal` identifies the canonical increasing
  intermediate tuple in the ordered Cauchy--Binet sum and proves that every
  truncation with at least `n` terms dominates the fixed principal determinant;
- `expMatrix` names the full positive exponential interaction matrix;
- `expTaylorMatrix_tendsto_expMatrix` derives entrywise convergence from
  mathlib's Banach-algebra exponential series, and
  `expTaylorMatrix_det_tendsto_expMatrix_det` composes it with determinant
  continuity;
- `expMatrix_det_pos` passes the fixed lower bound to the limit and combines it
  with `expTaylorPrincipalMatrix_det_pos`, avoiding the invalid inference that
  an arbitrary limit of positive numbers is positive;
- a 100-decimal regression through size five independently confirms positivity
  and the principal Vandermonde-factorial lower bound.

Verification:

- the focused theorem build and public axiom audit pass with only Lean's
  default axioms;
- all twelve focused surrogate tests pass;
- `./scripts/verify.sh` passes all repository checks with 53 Python tests and
  no known dependency vulnerabilities.

Next analytic obligation:

Translate arbitrary increasing row and column tuples to nonnegative tuples and
factor `exp(x*y)` into the translated kernel times positive row and column
weights. Then apply the positive determinant theorem to the negated reversed
frequency tuple and prove that the reversal permutation has sign
`(-1)^(n.choose 2)`. These finite algebraic steps are the remaining path to
`StrictSignRegular expKernel expKernelSignature`.

Non-claims:

- this phase does not yet remove the nonnegative-node restriction by
  exponential row/column scaling;
- it does not yet reverse frequency orientation or prove
  `StrictSignRegular expKernel expKernelSignature`;
- it does not establish dyadic localization, Conjecture G1, or Problem 4.2.

## Phase P — all-orders exponential sign regularity

State: complete with green CI in merged
[PR #20](https://github.com/tripp-smith/gecp-kernel-structure/pull/20).

Phase O proves strict positivity for `exp(x*y)` when both ordered node tuples
are nonnegative. Only finite scaling and orientation steps remain before the
all-orders premise already consumed by the Phase J fermionic transfer.

Frozen public Lean targets:

- `expMatrix_det_pos_of_strictMono`: arbitrary strictly increasing real row
  and column tuples give a positive `exp(x*y)` determinant;
- `Fin.sign_revPerm`: reversal on `Fin n` has sign
  `(-1)^(n.choose 2)`;
- `expKernel_strictSignRegular`: the kernel `exp(-t*omega)` is strictly sign
  regular at every order with signature `expKernelSignature`;
- `fermionicKernel_strictSignRegular`: the existing positive column-scaling
  transfer yields the same all-orders theorem for the fermionic kernel.

Proof contract:

- for nonempty tuples, subtract the first row and column nodes, apply
  `expMatrix_det_pos`, and recover the original matrix by positive exponential
  row and column scalings; discharge the empty order directly;
- compute `Fin.revPerm.sign` from mathlib's inversion-product formula, where
  every pair is inverted;
- apply arbitrary-node exponential positivity to the increasing tuple
  `j ↦ -columns (Fin.rev j)`, then permute columns back and account for the
  reversal sign;
- use the already-proved fermionic positive column scaling without adding a
  new analytic assumption.

Implementation outcome:

- `expMatrix_det_pos_of_strictMono` handles the empty order directly and, for
  nonempty tuples, subtracts the first nodes, applies `expMatrix_det_pos`, and
  reconstructs the original determinant through positive diagonal row and
  column factors;
- `Fin.sign_revPerm` computes the reversal sign from mathlib's product over
  inversions and the identity `sum(range n) = n.choose 2`;
- `expKernel_strictSignRegular` applies arbitrary-node positivity to
  `j ↦ -columns (Fin.rev j)` and permutes the columns back with exactly the
  signature `(-1)^(n.choose 2)`;
- `fermionicKernel_strictSignRegular` instantiates the existing positive
  frequency-weight transfer, and
  `fermionicKernel_pivotCrossProductSignCoherent` closes the selected-residual
  sign premise for every successful fermionic run;
- a 100-decimal regression through order six independently confirms the
  predicted signs on tuples containing both negative and positive nodes.

Verification:

- the focused theorem build and expanded public axiom audit pass with only
  Lean's default axioms;
- all thirteen focused surrogate tests pass;
- `./scripts/verify.sh` passes all repository checks with 54 Python tests and
  no known dependency vulnerabilities.

Next analytic obligation:

The Phase J--P chain now proves exact nonexpansiveness at every selected
fermionic pivot. To obtain Conjecture G1, combine this invariant with a strict
contraction mechanism on each dyadic frequency block, extending the proved
two-corner `0 < Λ ≤ 1` base case through cutoff rescaling or localization.

Non-claims:

- strict sign regularity supplies selected-residual nonexpansiveness, not the
  dyadic strict contraction needed for Conjecture G1;
- this phase does not prove cutoff localization, the full GECP rate, or
  Problem 4.2.

## Phase Q — complete-run nonexpansiveness

State: complete with green CI in merged
[PR #21](https://github.com/tripp-smith/gecp-kernel-structure/pull/21).

Phase P proves sign coherence for the final residual of every finite prefix of
a fermionic run. The next exact obligation is to compose prefixes with
continuations and turn that pointwise statement into a theorem for an entire
exact complete-pivot run on a restricted domain.

Frozen public Lean targets:

- `Run.append` and `Run.finalResidual_append`: compose dependent successful
  runs without forgetting the original kernel or changing the final residual;
- `GECP.CompletePivotOn` and `GECP.Run.CompleteOn`: record that selected pivots
  belong to the chosen row and column domains and maximize the current
  residual there;
- `residualUpdate_le_of_signCoherentOn`: derive the factor-one update bound
  from bounds on only the four points in the selected residual cross;
- `strictSignRegular_gecp_error_nonincreasing`: every realized exact
  complete-pivot residual sequence of a strictly sign-regular kernel preserves
  any initial uniform residual bound on the restricted domain;
- `fermionicKernel_gecp_error_le_cutoffCorner`: on
  `[0,1] x [-Lambda,Lambda]`, every finite exact complete-pivot run has final
  residual bounded by the initial cutoff-corner maximum.

Proof contract:

- retain every residual's successful original-kernel prefix and induct on the
  residual index;
- obtain sign coherence at the current residual from strict sign regularity of
  the original kernel, then apply a domain-local factor-one update theorem;
- use domain membership of the selected pivot to compare its magnitude with
  the inherited uniform bound;
- specialize to the fermionic kernel using positivity and the proved cutoff
  corner maximum, without adding compactness or supremum assumptions.

Independent checks:

- a finite exact rational run verifies the final residual after appending two
  one-step runs;
- focused Lean builds and the public axiom audit cover every new structural
  theorem before the full root verification command.

Decision rule:

- if dependent run composition cannot expose the required prefix to the Phase
  P theorem, stop after three distinct proof routes and record the exact type
  obstruction rather than weakening complete pivoting into a conclusion-shaped
  assumption.

Implementation outcome:

- `Run.append` composes a dependent successful prefix with a continuation, and
  `Run.finalResidual_append` proves that the composed run ends at exactly the
  continuation residual;
- `GECP.CompletePivotOn` and `GECP.Run.CompleteOn` express domain membership
  and exact residual maximization without asserting a global maximum outside
  the GECP rectangle;
- `residualUpdate_le_of_signCoherentOn` closes the domain mismatch discovered
  during the first focused build: only the four entries in the Schur cross
  need the common pivot bound;
- `strictSignRegular_gecp_error_nonincreasing` uses a successful original-
  kernel prefix realizing every residual, so Phase P supplies sign coherence
  at every step and the initial domain bound propagates to every finite rank;
- `fermionicKernel_gecp_error_le_cutoffCorner` instantiates the generic theorem
  on `[0,1] x [-Lambda,Lambda]` using fermionic positivity and the exact cutoff-
  corner maximum.

Verification:

- the focused module and public re-export build;
- the four new structural theorems in the axiom audit use only Lean's default
  axioms;
- an exact two-by-two rational regression composes two one-step dependent runs
  and checks that every entry of the final residual is zero;
- `./scripts/verify.sh` passes with 54 Python tests and no known dependency
  vulnerabilities.

Next analytic obligation:

Nonexpansiveness is now a theorem for the whole realized sequence, rather than
only a one-step consequence. The remaining rate problem is therefore isolated
to a strict statement: prove a half reduction within `O(s+1)` complete pivots
when `Lambda <= 2^s`, most plausibly by dyadic pivot localization or by a
determinant/near-volume estimate that transports the cutoff-one two-corner
mechanism.

Non-claims:

- factor-one nonexpansiveness is not strict contraction;
- this phase does not localize pivots to dyadic frequency bands, extend the
  cutoff-one half-contraction, prove Conjecture G1, or solve Problem 4.2.
