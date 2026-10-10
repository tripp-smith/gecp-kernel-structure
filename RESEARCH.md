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

## Phase T — geometric-mean determinant reduction

State: complete with green CI in merged
[PR #22](https://github.com/tripp-smith/gecp-kernel-structure/pull/22).

The finite-grid pivot locations jump across dyadic bands after the two endpoint
pivots, so a simple outer-to-inner localization order is not a credible next
contract. A sharper route is the geometric-mean determinant estimate for
pivoted LU developed in Marc Aurèle Gilles,
[Convergence rates for pivoted QR and LU](https://arxiv.org/abs/2607.26863)
(2026). This repository already proves the pivot-product identity, all-orders
fermionic sign regularity, and complete-sequence nonexpansiveness needed to
make the exact-pivot functional reduction especially direct.

Frozen public Lean targets:

- `strictSignRegular_gecp_pivotMagnitude_antitone`: selected complete-pivot
  magnitudes in a realized strictly sign-regular residual sequence are
  antitone;
- `GECP.Run.abs_finSelectedCore_det_eq_prod_abs_pivots`: the absolute selected-
  core determinant is the product of absolute pivot magnitudes;
- `strictSignRegular_gecp_error_pow_le_selectedCore_det`: if a successful run
  records the first `n` pivots of the realized sequence, every domain residual
  entry at rank `n`, raised to the `n`th power, is bounded by the absolute
  selected-core determinant.

Proof contract:

- derive consecutive pivot monotonicity from the Phase Q localized update
  theorem, using completeness at the current and next selected points;
- compare the constant finite product with the pivot product using mathlib's
  `Finset.prod_le_prod`, rather than introduce a local geometric-mean
  abstraction;
- combine the exact pivot list with `gecp_core_det_eq_prod_pivots`, distributing
  absolute value over the finite product;
- retain the power-form inequality and do not introduce real roots merely for
  presentation.

Independent checks:

- the exact rational two-pivot run has pivot magnitudes `2` and `1/2`, so the
  final-pivot square `1/4` is below the selected-core determinant `1`;
- focused builds and the public axiom audit cover all new structural results.

Decision rule:

- if the dependent run API cannot connect the sequence pivot list to the
  selected core without a conclusion-shaped assumption, stop at the antitone
  pivot-product theorem and record the missing prefix-coherence field exactly.

Implementation outcome:

- `GECP.Run.RealizesPivotPrefix` records exactly the alignment needed between
  a dependent run's stored pivots and the first `n` pivots of a residual
  sequence, without assuming the determinant conclusion;
- `strictSignRegular_gecp_pivotMagnitude_antitone` combines Phase Q's localized
  sign-coherent update bound with exact completeness at consecutive selected
  points, proving that complete-pivot magnitudes cannot increase;
- `strictSignRegular_gecp_pivotMagnitude_pow_le_product` compares the rank-`n`
  pivot raised to the `n`th power with the product of the preceding pivot
  magnitudes;
- `GECP.Run.abs_finSelectedCore_det_eq_prod_abs_pivots` converts the existing
  signed determinant identity into the required absolute product identity;
- `strictSignRegular_gecp_error_pow_le_selectedCore_det` uses rank-`n`
  completeness to extend the pivot inequality to every point in the chosen
  domain, then identifies the product with the selected-core determinant.

Verification:

- the focused module, public re-export, exact small-instance check, and public
  axiom audit pass;
- every new public theorem uses only Lean's permitted default axioms;
- the exact rational two-pivot run verifies the power-form determinant bound;
- `./scripts/verify.sh` passes with 54 Python tests and no known dependency
  vulnerabilities.

Next analytic obligation:

The GECP-specific part of the geometric-mean argument is now formalized. The
remaining strict-rate problem is a determinant estimate: combine the proved
dyadic separated approximation with a quantitative selected-core determinant
perturbation or near-volume bound strong enough to force decay in the right
side of the power inequality.

Non-claims:

- a determinant-power reduction is not a determinant decay bound;
- this phase does not yet combine the dyadic separated approximation with a
  determinant perturbation estimate, prove strict block contraction,
  Conjecture G1, or Problem 4.2.

## Phase U — low-rank perturbation determinant bound

State: complete with green CI in merged
[PR #23](https://github.com/tripp-smith/gecp-kernel-structure/pull/23).

Phase T reduces exact complete-pivot residuals to selected-core determinants.
The next bridge is suggested independently by several neighboring theories:

- Gilles's
  [geometric-mean analysis of pivoted LU](https://arxiv.org/abs/2607.26863)
  controls functional LU by arbitrary sampled determinants;
- reduced-basis greedy theory compares greedy errors with Kolmogorov widths,
  including geometric means of width sequences
  ([DeVore--Petrova--Wojtaszczyk](https://arxiv.org/abs/1204.2290),
  [Nguyen](https://arxiv.org/abs/1804.03935));
- max-volume cross-approximation theory likewise isolates determinant and
  growth-factor control
  ([Cortinovis--Kressner--Massei](https://arxiv.org/abs/1902.02283));
- in exterior-algebra language, the top wedge of a rank-`r` approximant plus
  an error has no term with more than `r` approximant columns. Therefore every
  surviving `n`-column term contains at least `n-r` error columns.

The last observation is the strongest fit to the repository: it uses the
already formalized separated factors directly, avoids singular-value
infrastructure, and preserves the crucial power of the uniform error.

Frozen public Lean targets:

- `Matrix.det_mul_rect_eq_zero_of_card_lt`: a rectangular product through a
  strictly smaller finite index type has zero square determinant;
- `Matrix.det_add_eq_sum_columnChoices`: expand the determinant of `A + E`
  over all choices of approximation and error columns;
- `Matrix.abs_det_le_factorial_mul_prod_columnBounds`: a columnwise version of
  mathlib's uniform `Matrix.det_le` bound;
- `Matrix.abs_det_le_of_factors_approx`: if `B = L * R` factors through `r`
  coordinates, `A-B` is entrywise at most `epsilon`, the entries of `B` are at
  most `C`, `epsilon ≤ 1`, `1 ≤ C`, and `r ≤ n`, then
  `|det A|` is bounded by
  `2^n * n! * C^r * epsilon^(n-r)`;
- `fermionicKernel_sample_det_le_separatedApprox`: specialize the generic
  theorem to arbitrary samples in `[0,1] x [-2^s,2^s]`, using the delivered
  `16p(s+1)`-term approximation and the kernel bound by one.

Proof contract:

- define a column-choice matrix and prove the determinant expansion from
  `Matrix.det_apply'` plus `Fintype.prod_add`;
- factor every mixed column-choice matrix through the sum of the separated
  coordinate type and the chosen error-column subtype;
- use the repository's rectangular Cauchy--Binet theorem and finite-cardinality
  pigeonhole principle to kill choices with fewer than `n-r` error columns;
- bound every remaining mixed determinant columnwise and sum over at most
  `2^n` choices;
- represent `evalSeparated` as a rectangular matrix product indexed by
  `Fin terms.length`, then instantiate the fermionic pointwise error theorem;
- keep the result in power form. Root extraction and asymptotic simplification
  belong to a later phase.

Independent checks:

- exactly enumerate rational matrices of sizes two through five, random
  low-rank factorizations, and rational perturbations; verify both the
  column-choice identity and determinant bound without floating point;
- include a Lean small-instance theorem where a rank-one `3 x 3` approximant
  forces at least two perturbation columns;
- focused builds and the public axiom audit cover every public structural
  theorem.

Decision rule:

- after three distinct failed routes for mixed-column rank vanishing, stop and
  deliver the exact column-choice expansion plus the minimal missing
  factorization lemma; do not replace rank structure with a conclusion-shaped
  assumption.

Non-claims:

- this phase does not yet combine the determinant bound with a realized GECP
  run or simplify its constants into Conjecture G1;
- the factorial and `2^n` constants are deliberately explicit and coarse;
- related-work analogies guide the proof architecture but do not substitute
  for the Lean theorem or exact checks.

Delivered identifiers:

- `Matrix.det_mul_rect_eq_zero_of_card_lt` proves singularity of a square
  product through a strictly smaller finite coordinate type;
- `Matrix.det_add_eq_sum_columnChoices` is the exact exterior-power expansion
  over approximation/error column choices;
- `Matrix.mixedLeft_mul_mixedRight` factors each mixed choice through the sum
  of the low-rank coordinate type and the selected error-column subtype;
- `Matrix.abs_det_le_factorial_mul_prod_columnBounds` and
  `Matrix.abs_det_columnChoice_le` retain separate bounds for error and
  approximant columns;
- `Matrix.abs_det_le_of_factors_approx` gives the exact surviving-subset sum,
  while `Matrix.abs_det_le_two_pow_mul_of_factors_approx` gives the frozen
  coarse bound `2^n n! epsilon^(n-r) C^r`;
- `fermionicKernel_sample_det_le_separatedApprox` instantiates the exact sum
  at rank `16p(s+1)`, and `fermionicKernel_sample_det_le_two_pow` gives the
  closed arbitrary-sample determinant decay bound.

The proof route succeeded without a fallback. Exact `Fraction` regressions
enumerate all column subsets for deterministic low-rank rational examples of
sizes two through five, check the expansion, rank-forced vanishing, individual
column bounds, and both summed bounds. A Lean `3 x 3` instance checks that a
rank-one background plus one error column has zero determinant.

Verification:

- focused Lean builds and the public axiom audit pass, with only Lean's
  permitted default axioms;
- the exact rational subset census passes for sizes two through five;
- `./scripts/verify.sh` passes with 55 Python tests and no known dependency
  vulnerabilities.

Next analytic obligation:

Combine the closed fermionic determinant bound with Phase T's residual-power
inequality, select `p` as a function of sample order and cutoff scale, and
determine whether the factorial and subset constants can be absorbed into a
cutoff-uniform strict rate. This arithmetic optimization is separate from the
now-complete determinant perturbation mechanism.

## Phase V — GECP/determinant composition

State: complete with green CI in merged
[PR #24](https://github.com/tripp-smith/gecp-kernel-structure/pull/24).

Phase T controls a complete-pivot residual power by the determinant of a
realized run's selected core. Phase U controls an arbitrary physical-domain
sample determinant by the explicit fermionic separated approximation. The
remaining structural composition must verify that the run's recursively
selected coordinates are physical-domain samples; equality of pivot values
alone does not encode that fact.

Frozen public Lean targets:

- `GECP.Run.selectedRow_mem_of_completeOn` and
  `GECP.Run.selectedColumn_mem_of_completeOn`: every recursive selected
  coordinate of a domain-complete run lies in its declared domain;
- `GECP.Run.selectedCore_eq_sampleMatrix`: identify the original-kernel
  selected core with the generic sampled-matrix representation;
- `fermionicKernel_gecp_error_pow_le_separatedApprox`: compose Phase T with
  Phase U's exact surviving-subset determinant sum;
- `fermionicKernel_gecp_error_pow_le_two_pow`: give the closed bound
  `|R_n(x,y)|^n <= 2^n n! 2^(-p(n-r)) 2^r`, with
  `r = 16p(s+1)`, for a realized physical-domain exact complete-pivot run.

Proof contract:

- prove selected-coordinate membership by induction on the dependent run and
  case analysis on its recursive sum index;
- unfold `Run.selectedCore` only in the small extensional bridge theorem;
- derive `run.pivots.length = n` from `Run.RealizesPivotPrefix` rather than
  adding it as an assumption;
- reuse the Phase T residual-power theorem and Phase U sample-determinant
  theorems without reproving either estimate;
- keep the free approximation order `p` explicit. Optimization and root
  extraction are admitted into this phase only if they reduce to stable
  natural-number and real-power lemmas without changing the contract.

Independent checks:

- instantiate the composition theorem on the exact two-pivot rational-style
  run pattern at a small fermionic sample, or add a structurally equivalent
  finite check if transcendental normalization prevents `norm_num` closure;
- focused builds, public re-export, axiom audit, and root verification.

Non-claims:

- equality of stored pivot values is not asserted to identify selected
  coordinates; `Run.CompleteOn` supplies the missing domain invariant;
- this phase does not claim an optimized nth-root rate, Conjecture G1, or
  Problem 4.2 unless the explicit constants are actually discharged in Lean.

Delivered identifiers:

- `GECP.Run.selectedRow_mem_of_completeOn` and
  `GECP.Run.selectedColumn_mem_of_completeOn` prove the recursive domain
  invariant, and `GECP.Run.selectedCore_eq_sampleMatrix` supplies the exact
  matrix bridge;
- `fermionicKernel_gecp_error_pow_le_separatedApprox` and
  `fermionicKernel_gecp_error_pow_le_two_pow` compose the Phase T residual
  power with the Phase U determinant estimates;
- `fermionicKernel_gecp_error_le_oddBlock` chooses approximation order
  `p = 2m+1` and block length `n = 32(2m+1)(s+1)`, then extracts the exact
  root to prove `|R_n(x,y)| <= 2n 2^-m`;
- `oddBlock_factor_le_half` reduces half contraction to
  `128(2m+1)(s+1) <= 2^m`, while `oddBlock_scale_quadratic` proves that
  `m = 16(s+1)` always satisfies this condition;
- `fermionicKernel_gecp_error_le_half_quadraticBlock` therefore proves error
  at most one half by
  `n = 32(32(s+1)+1)(s+1)` for every realized physical-domain exact run.

The odd-order choice is a divisibility device: after the separated rank is
set to half the sample order, the remaining dyadic determinant exponent is an
exact multiple of `n`. This converts the geometric-mean power estimate to an
ordinary residual estimate without logarithms or real roots. The quadratic
choice is deliberately elementary; sharper choices of `m` should recover an
`O((s+1)(log(s+1)+log(1/epsilon)))` rank shape, but that asymptotic inversion
is not claimed here.

Updated non-claims:

- Phase V proves a strict rate, but not the conjectured `2(s+1)` half-
  contraction block or a cutoff-uniform constant with the target linear scale;
- no claim is made that the quadratic block is sharp.

Verification:

- focused Lean builds and the expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- an independent integer regression checks the quadratic scale inequality for
  scales zero through 32;
- `./scripts/verify.sh` passes with 56 Python tests and no known dependency
  vulnerabilities.

## Phase W — logarithmic-scale contraction block

State: complete with green CI on 2026-10-09; delivered in merged PR #25.

Phase V deliberately used `m = 16(s+1)` to discharge
`128(2m+1)(s+1) <= 2^m` with elementary polynomial arithmetic. The exact
odd-block theorem only needs this inequality, so the natural next step is the
ceil-log choice
`m = 16 + 2 * Nat.clog 2 (s+1)`. Writing
`ell = Nat.clog 2 (s+1)`, mathlib supplies `s+1 <= 2^ell`; the remaining
inequality is a uniform comparison of the linear factor `33+4ell` with an
exponential in `ell`.

Frozen public Lean targets:

- `oddBlock_scale_logarithmic`: prove
  `128(2m+1)(s+1) <= 2^m` for
  `m = 16 + 2 * Nat.clog 2 (s+1)`;
- `fermionicKernel_gecp_error_le_half_logarithmicBlock`: instantiate Phase V's
  half-contraction theorem at that `m`, giving block length
  `32(2m+1)(s+1)`;
- retain the exact natural-number expression in the theorem statement and
  document its `O((s+1) log(s+1))` interpretation separately.

Proof contract:

- use `Nat.le_pow_clog` rather than introducing a local logarithm abstraction;
- separate the bound `s+1 <= 2^ell`, the elementary
  `33+4ell <= 512 * 2^ell`, and the final power normalization;
- use induction only for the one-variable auxiliary exponential inequality;
- reuse `fermionicKernel_gecp_error_le_half_of_oddBlock` unchanged.

Independent checks:

- exact integer evaluation over a broad deterministic scale range, including
  dyadic boundaries where `Nat.clog` changes;
- focused Lean build, public axiom audit, and root verification.

Decision rule:

- after three failed normalizations of the same ceil-log inequality, expose
  the smallest missing arithmetic helper as the phase result rather than
  weakening the GECP assumptions.

Non-claims:

- this phase does not prove the conjectured `2(s+1)` block, remove all
  logarithmic overhead, establish optimal constants, or close Conjecture G1.

Delivered identifiers:

- `oddBlock_scale_logarithmic` proves the exact odd-block scale condition for
  `m = 16 + 2 * Nat.clog 2 (s+1)` using `Nat.le_pow_clog` and the elementary
  auxiliary comparison `33+4ell <= 512 * 2^ell`;
- `fermionicKernel_gecp_error_le_half_logarithmicBlock` instantiates the Phase
  V contraction theorem at that scale and proves error at most one half by
  `n = 32(33 + 4 * Nat.clog 2 (s+1))(s+1)`.

Result:

- the certified general block improves from quadratic in `s+1` to
  `O((s+1) log(s+1))` without strengthening any GECP premise;
- the exact expression, rather than asymptotic notation, remains the public
  theorem interface;
- the conjectured linear block, optimal constants, and Conjecture G1 remain
  open.

Verification:

- focused Lean builds and the expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- an independent exact-integer regression covers scales through 1024 and the
  neighbors of every dyadic boundary through `2^64`;
- `./scripts/verify.sh` passes with 57 Python tests and no known dependency
  vulnerabilities.

## Phase X — arbitrary dyadic accuracy

State: complete with green CI on 2026-10-09; delivered in merged PR #26.

Phase W proves one half contraction using an accuracy order large enough to
absorb the dimension factor in Phase V's estimate
`|R_n(x,y)| <= 2n 2^-m`. The same arithmetic has more content: for a requested
dyadic accuracy order `q`, taking
`m = 16 + 2 * Nat.clog 2 (s+1) + 2q` should absorb both the scale startup and
the target factor `2^-q`. This exposes the strongest direct convergence
theorem currently implied by the determinant route before attempting to
remove its logarithmic startup overhead.

Frozen public Lean targets:

- `oddBlock_factor_le_accuracy`: reduce the residual factor bound
  `2n 2^-m <= 2^-q` to the exact natural-number condition
  `64(2m+1)(s+1)2^q <= 2^m`;
- `oddBlock_scale_accuracy`: prove that condition for
  `m = 16 + 2 * Nat.clog 2 (s+1) + 2q`;
- `fermionicKernel_gecp_error_le_dyadicAccuracyBlock`: instantiate the Phase V
  odd-block theorem and prove error at most `(1/2)^q` by
  `n = 32(33 + 4 * Nat.clog 2 (s+1) + 4q)(s+1)`;
- keep the exact natural-number rank in the theorem and state its
  `O((s+1)(log(s+1)+q))` interpretation only in documentation.

Proof contract:

- reuse the Phase W exponential comparison at `ell+q` rather than introduce a
  second logarithmic abstraction;
- separate the natural-number scale inequality from the real-power factor
  bound and the GECP corollary;
- preserve the exact complete-pivot and physical-domain assumptions unchanged;
- search mathlib first for power/division normalization lemmas and expose a
  small helper if three normalization attempts fail.

Independent checks:

- exact integer evaluation across ordinary scales, dyadic boundaries, and a
  broad range of requested accuracy orders;
- focused Lean build, public axiom audit, and root verification.

Decision rule:

- if the proposed linear-in-`q` parameter fails, minimize the exact integer
  witness and deliver the weakest explicit affine choice that passes; do not
  weaken the GECP hypotheses.

Non-claims:

- the additive `(s+1) log(s+1)` startup remains, so this phase does not prove
  Conjecture G1, the conjectured `2(s+1)` half-contraction block, optimal
  constants, or Simons Problem 4.2.

Delivered identifiers:

- `oddBlock_factor_le_accuracy` converts the exact natural-number condition
  `64(2m+1)(s+1)2^q <= 2^m` into the real residual-factor bound
  `2n 2^-m <= 2^-q`;
- `oddBlock_scale_accuracy` proves that condition for
  `m = 16 + 2 * Nat.clog 2 (s+1) + 2q` by reusing the Phase W exponential
  comparison at `Nat.clog 2 (s+1) + q`;
- `fermionicKernel_gecp_error_le_dyadicAccuracyBlock` composes those arithmetic
  results with the Phase V odd-block theorem to prove error at most `(1/2)^q`
  by `n = 32(33 + 4 * Nat.clog 2 (s+1) + 4q)(s+1)`.

Result:

- arbitrary accuracy is obtained directly from the original complete-pivot
  run, without assuming that a contraction theorem can be restarted on an
  arbitrary residual;
- the rank is `O((s+1)(log(s+1)+q))`, retaining linear dependence on the
  requested dyadic accuracy order;
- for `s` comparable to `log(1+Lambda)` and `q` comparable to
  `log(1/epsilon)`, this gives
  `O(log(1+Lambda)(log log(1+Lambda)+log(1/epsilon)))`;
- the remaining gap to Conjecture G1 is now isolated to the additive
  `(s+1) log(s+1)` startup and constants, rather than the accuracy dependence.

Verification:

- focused Lean builds and the expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- an independent exact-integer regression checks 65 accuracy orders across
  scales through 1024 and every dyadic transition through `2^64`;
- `./scripts/verify.sh` passes with 58 Python tests and no known dependency
  vulnerabilities.

## Phase Y — determinant-prefactor obstruction

State: complete with green CI on 2026-10-09; delivered in merged PR #27.

Phase X isolates the remaining gap to an additive `(s+1) log(s+1)` startup.
In the current determinant route that loss comes from taking an `n`th root of
a generic entrywise perturbation estimate with a dimension-growing
determinant prefactor. A tempting next step is to replace that prefactor by
`C^n` for a universal constant `C`. Before investing in a fermionic proof, this
phase tests whether such an improvement can possibly follow from entrywise
low-rank approximation alone.

Two neighboring theories sharpen the decision. Gilles's geometric-mean LU
analysis (arXiv:2607.26863) confirms sampled determinants as the correct pivot
quantity, while modern empirical-interpolation estimates
(arXiv:2401.13985) retain Lebesgue and finite-dimensional volume factors.
Neither removes the dimension dependence without additional kernel structure.
The exact finite obstruction is the Sylvester--Hadamard family: bounded
entries coexist with determinant root growing like the square root of the
dimension.

Frozen public Lean targets:

- `SylvesterIndex` and `sylvesterHadamard`: construct the power-of-two
  Kronecker family over real matrices;
- `sylvesterHadamard_isHadamard` and `sylvesterIndex_card`: certify
  orthogonality, unit-modulus entries, and order `2^(k+1)`;
- `sylvesterHadamard_det_square`: prove the exact identity
  `det(H_k)^2 = (2^(k+1))^(2^(k+1))`;
- `entrywise_determinant_constant_base_obstruction`: for every natural base
  `C`, exhibit a bounded-entry matrix whose absolute determinant exceeds
  `C^n`, ruling out a universal constant-base determinant bound under only
  entrywise hypotheses.

Proof contract:

- use mathlib's `Matrix.IsHadamard.kronecker` and determinant identity rather
  than reprove row orthogonality or expand permutations;
- keep the obstruction over exact real/integer arithmetic;
- separate family construction, cardinality, determinant identity, and the
  asymptotic witness;
- after three failures on the quantified witness arithmetic, retain the exact
  family identity as the certified obstruction and state the remaining
  elementary corollary explicitly rather than weaken the matrix assumptions.

Independent checks:

- recursively construct small Sylvester matrices with exact integers and
  compare their determinants with the closed form;
- focused Lean build, public axiom audit, and root verification.

Decision rule:

- if the obstruction is certified, stop pursuing generic entrywise
  determinant improvements and make a fermionic-specific total-positivity,
  divided-difference, or residual-relative estimate the next research target.

Non-claims:

- this phase does not show that the fermionic determinant has Hadamard growth;
  the obstruction matrices are not totally positive;
- it does not disprove a dimension-free bound under fermionic sign regularity
  or residual-selected sampling, and it does not prove Conjecture G1.

Delivered identifiers:

- `SylvesterIndex`, `hadamardTwo`, and `sylvesterHadamard` define exact real
  Sylvester--Hadamard matrices of order `2^(k+1)` by repeated Kronecker
  products;
- `hadamardTwo_isHadamard`, `sylvesterHadamard_isHadamard`, and
  `sylvesterIndex_card` certify unit-modulus entries, orthogonality, and the
  exact matrix order;
- `sylvesterHadamard_det_square` proves
  `det(H_k)^2 = (2^(k+1))^(2^(k+1))` using mathlib's Hadamard determinant
  identity;
- `entrywise_determinant_constant_base_obstruction` proves that for every
  natural `C` some unit-entry square matrix of order `n` satisfies
  `C^n < |det A|`.

Result:

- a determinant bound of the form `|det A| <= C^n` with fixed `C` cannot be
  derived from a uniform entrywise bound alone: the Hadamard family has
  determinant root `sqrt(n)`;
- therefore replacing Phase U's dimension factor by a universal exponential
  base requires additional fermionic or residual-selected hypotheses, not a
  sharper generic determinant inequality;
- the next promising proof target is a determinant or divided-difference
  estimate for totally positive exponential collocation matrices, ideally
  localized to the rows and columns selected by the GECP residual.

Verification:

- the focused Lean library build and expanded public axiom audit pass, with
  only Lean's permitted default axioms;
- an independent exact-integer regression checks the determinant identity for
  Sylvester orders two through 64 and the arithmetic witness for bases zero
  through 256;
- `./scripts/verify.sh` passes with 59 Python tests and no known dependency
  vulnerabilities.

## Phase Z — multiplicative-Monge corner localization

State: complete with green CI on 2026-10-09; delivered in merged PR #28.

Phase Y proves that bounded entries and low-rank approximation alone cannot
remove the determinant route's dimension loss. The next target must use an
identity special to the fermionic kernel. Its logarithm is affine in the
time-frequency interaction, so every four-point multiplicative cross ratio
collapses to one scalar exponential. Consequently, elimination at a corner of
an oppositely ordered rectangle has no cancellation ambiguity: the updated
residual is the original positive kernel times an explicit area factor.

Frozen public Lean targets:

- `fermionicKernel_cross_product`: prove the exact four-point exponential
  cross-product identity;
- `fermionicKernel_residualUpdate_factor`: factor an arbitrary one-step
  fermionic residual as
  `K(t,omega) * (1 - exp((t-t0)*(omega-omega0)))`;
- `fermionicKernel_cornerResidual_nonneg`: on `t0 <= t` and `omega <= omega0`,
  certify that the corner residual is nonnegative;
- `fermionicKernel_cornerResidual_le_area`: bound its absolute value by the
  dimensionless rectangle area `(t-t0)*(omega0-omega)` on `0 <= t <= 1`;
- `fermionicKernel_cornerResidual_le_half`: obtain a local half contraction
  whenever that area is at most `log 2`.

Proof contract:

- derive the factorization from the kernel definition and exponential laws,
  not by adding a cross-ratio hypothesis;
- keep the arbitrary pivot identity separate from the ordered-corner
  inequalities;
- use `1-u <= exp(-u)` for the linear area bound and exponential monotonicity
  for the `log 2` threshold;
- expose only exact real inequalities; asymptotic and covering language stays
  in documentation.

Independent checks:

- recompute the update and factorized forms independently at high precision
  on deterministic interior points and both reflected corner orientations;
- verify nonnegativity, the area bound, and the `log 2` half threshold on a
  dense deterministic sample of admissible local rectangles;
- focused Lean build, public axiom audit, and root verification.

Decision rule:

- if the exact factorization holds but the proposed area constants fail,
  preserve the identity and deliver the sharpest direct monotone-exponential
  bound; do not add assumptions about later GECP residuals.

Non-claims:

- a local corner theorem for the original kernel does not show that every
  later GECP residual has the same multiplicative cross ratio;
- this phase does not yet cover the cutoff rectangle by selected residual
  tiles, prove the conjectured `2(s+1)` block, prove Conjecture G1, or solve
  Simons Problem 4.2.

Delivered identifiers:

- `fermionicKernel_cross_product` proves that the four kernel factors differ
  by exactly `exp((t-t0)*(omega-omega0))`; the logistic column normalization
  cancels completely;
- `fermionicKernel_residualUpdate_factor` converts the arbitrary one-step
  Schur update into `K(t,omega) * (1-exp((t-t0)*(omega-omega0)))`;
- `fermionicKernel_residualUpdate_nonneg_of_cross_nonpos`,
  `fermionicKernel_residualUpdate_le_crossArea`, and
  `fermionicKernel_residualUpdate_le_half_of_crossArea` express the result
  without choosing an orientation, so both reflected corners are covered;
- `fermionicKernel_cornerResidual_nonneg` specializes the identity to an
  oppositely ordered corner and removes cancellation from the update;
- `fermionicKernel_cornerResidual_le_area` uses `K <= 1` and
  `1-exp(-u) <= u` to bound the residual by the time-frequency cross area;
- `fermionicKernel_cornerResidual_le_half` sharpens the monotone exponential
  bound to one half when that cross area is at most `log 2`.

Result:

- the original fermionic kernel has an exact local contraction law unavailable
  to the arbitrary bounded matrices from Phase Y;
- a corner pivot controls an entire oppositely ordered tile using only its
  dimensionless time-frequency area, giving a scale-invariant target for a
  dyadic tiling argument;
- the remaining global obligation is to show that GECP-selected residual
  pivots induce or dominate a collection of such tiles across only `O(s+1)`
  frequency scales.

Verification:

- the focused Lean build and expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- a 100-decimal regression independently compares the Schur update and
  factorized formulas at arbitrary and reflected corner points, then checks
  nonnegativity, the area bound, and the `log 2` threshold on 4,590
  deterministic oriented tile samples;
- `./scripts/verify.sh` passes with 60 Python tests and no known dependency
  vulnerabilities.

## Phase AA — two-corner power-secant residual

State: complete with green CI on 2026-10-09; delivered in merged PR #29.

Phase Z gives an exact area law for one update of the original kernel, but a
global GECP proof needs structure that survives beyond one pivot. The smallest
nontrivial case is the run selecting time endpoints `0,1` and two arbitrary
frequency endpoints `omegaLow < omegaHigh`. Under the positive coordinate
change `z = exp(-omega)`, the kernel becomes `z^t/(1+z)`. Direct elimination
then predicts that the two-step residual is the concave power `z^t` minus its
linear secant between the transformed frequency endpoints, divided by `1+z`.

Frozen public Lean targets:

- `fermionicKernel_eq_rpow`: prove the exact logistic-power coordinate
  identity `K(t,omega) = exp(-omega)^t / (1+exp(-omega))`;
- `powerSecant`: define linear endpoint interpolation of `z^t`;
- `powerSecant_le_rpow`: use mathlib's `Real.concaveOn_rpow` to prove that the
  secant lies below `z^t` for `0 <= t <= 1` and `a <= z <= b`;
- `fermionicKernel_firstCornerResidual_eq_rpow_sub`: express the first corner
  residual as `(z^t-a^t)/(1+z)`;
- `fermionicKernel_firstCorner_secondPivot_pos`: certify the arbitrary-band
  second pivot at `(1,omegaLow)` when `omegaLow < omegaHigh`;
- `fermionicKernel_twoCornerResidual_eq_powerSecantError`: identify the exact
  nested residual with `(z^t-powerSecant a b t z)/(1+z)`;
- `fermionicKernel_twoCornerResidual_nonneg`: certify nonnegativity throughout
  the frequency band for physical times.

Proof contract:

- prove the coordinate and residual identities algebraically before invoking
  concavity;
- use the existing `residualUpdate` definition and actual nonzero pivot proof,
  not a standalone interpolation surrogate;
- support asymmetric frequency bands rather than only the existing symmetric
  cutoff pair;
- keep quantitative secant-error bounds and dyadic covering as subsequent
  obligations unless they follow without adding derivative machinery.

Independent checks:

- compare the nested Schur update and secant-error formula at 100 decimal
  digits on asymmetric and reflected bands;
- densely sample admissible physical times and interior frequencies to check
  the proved sign and endpoint zeros;
- focused Lean build, public axiom audit, and root verification.

Decision rule:

- if exact reduction succeeds but formal concavity cannot be connected after
  three proof attempts, deliver the reduction and record the precise missing
  convex-combination lemma rather than replace it with a numerical claim.

Non-claims:

- nonnegativity of the two-corner residual is not yet a cutoff-uniform error
  bound or a statement about the actual third complete pivot;
- this phase does not prove inheritance for arbitrary pivot histories, the
  dyadic tile covering, Conjecture G1, or Simons Problem 4.2.

Delivered identifiers:

- `fermionicKernel_eq_rpow` moves the kernel exactly to the positive
  coordinate `z = exp(-omega)` as `z^t/(1+z)`;
- `powerSecant_le_rpow` proves the concave-power secant inequality on every
  positive interval, while `mul_rpow_sub_one_le_powerSecant` supplies the
  complementary lower line through the right endpoint;
- `fermionicKernel_firstCornerResidual_eq_rpow_sub` and
  `fermionicKernel_secondCornerResidual_eq_rpow_sub` identify the two
  reflected one-corner power errors;
- `fermionicKernel_firstCorner_secondPivot_pos` certifies the actual nonzero
  second Schur pivot for every strictly ordered frequency band;
- `fermionicKernel_twoCornerResidual_eq_powerSecantError` identifies the
  nested residual exactly with the normalized power-secant error;
- `fermionicKernel_twoCornerResidual_nonneg` proves its sign from concavity;
- `fermionicKernel_twoCornerResidual_le_upperCrossArea` and
  `fermionicKernel_twoCornerResidual_le_lowerCrossArea` compare the same
  residual to the two reflected Phase Z corner errors;
- `fermionicKernel_twoCornerResidual_mem_Icc_minCrossArea` packages the result
  as the tent bound
  `0 <= R2 <= min(t*(omegaHigh-omega), (1-t)*(omega-omegaLow))`.

Result:

- the two prescribed endpoint-time pivots retain exact structure after the
  first elimination, rather than merely inheriting a coarse factor-one bound;
- the residual vanishes on all four edges of its time-frequency band and is
  quantitatively localized by distance to both transformed corners;
- unlike the earlier symmetric small-cutoff base case, the theorem holds on
  every asymmetric finite band with no cutoff-size restriction;
- the next research obligation is no longer a secant estimate: it is to
  connect actual complete-pivot selection or a dyadic band decomposition to a
  controlled collection of these residual tents.

Verification:

- the focused Lean build and expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- a 100-decimal regression independently compares the nested Schur update and
  power-secant formula, then checks sign, all edge zeros, and both tent sides
  at 1,156 deterministic points over asymmetric and reflected bands;
- `./scripts/verify.sh` passes with 61 Python tests and no known dependency
  vulnerabilities.

## Phase AB — two-corner half-contraction obstruction

State: complete with green CI on 2026-10-09; delivered in merged PR #30.

Phase AA suggests the tempting shortcut that the first two symmetric corner
pivots might reduce the residual by one half for every cutoff, which would be
far stronger than Conjecture G1. Exact exploration of the power-secant formula
shows that this shortcut is false. The smallest clean witness found uses
transformed symmetric endpoints `a = 24^-3`, `b = 24^3`, physical time
`t = 2/3`, and interior coordinate `z = (5/4)^3`. All fractional powers are
rational, so the obstruction can be audited without floating-point bounds.

Frozen public Lean targets:

- `powerSecant_twoThirds_counterexample_value`: evaluate the normalized
  power-secant error at the rational witness as exactly
  `4495348 / 8973531`;
- `powerSecant_twoThirds_counterexample_gt_half`: prove that value exceeds
  half the corresponding initial symmetric-corner pivot
  `(1/2) * (13824/13825)` by the exact positive margin
  `222676 / 224338275`;
- `fermionicKernel_twoCornerResidual_counterexample_value`: transport the
  rational-coordinate identity to the actual nested `residualUpdate` at
  cutoff `log 13824` and frequency `-log(125/64)`;
- `fermionicKernel_twoCornerResidual_not_half_contraction`: state the resulting
  strict failure of two-pivot half contraction relative to the initial
  complete-pivot value.

Proof contract:

- derive the witness through Phase AA's exact power-secant identity;
- discharge every fractional power and final inequality exactly with rational
  arithmetic and mathlib's real-power normalization;
- retain the actual nonzero second-pivot proof in the nested residual;
- make no conclusion about longer complete-pivot blocks from this two-pivot
  obstruction.

Independent checks:

- recompute the exact secant, residual, initial pivot, and positive margin with
  Python `Fraction` arithmetic;
- scan a deterministic high-precision cutoff family to confirm the witness is
  part of the large-cutoff failure regime rather than an isolated transcription
  error;
- run the focused Lean build, public axiom audit, and root verification.

Decision rule:

- if the direct logarithmic transport is awkward but the positive-coordinate
  theorem compiles, isolate the missing exponential rewrite as a helper lemma;
  do not weaken the exact inequality to a decimal approximation.

Non-claims:

- failure after two prescribed pivots does not refute Conjecture G1, whose
  block length grows with cutoff scale;
- the witness does not rule out a uniform constant after more pivots, dyadic
  localization, residual dominance, or the existing proved
  `O((s+1) log(s+1))` contraction block;
- this phase does not solve Simons Problem 4.2.

Delivered identifiers:

- `powerSecant_twoThirds_counterexample_value` proves that the normalized
  power-secant error at the rational witness is exactly
  `4495348 / 8973531`;
- `powerSecant_twoThirds_counterexample_margin` identifies its exact excess
  over half the initial pivot as `222676 / 224338275`, and
  `powerSecant_twoThirds_counterexample_gt_half` records strict positivity;
- `fermionicKernel_twoCornerResidual_counterexample_value` transports the
  rational identity through exponentials and logarithms to the actual nested
  fermionic Schur update;
- `fermionicKernel_twoCornerResidual_counterexample_gt_half_initial` gives the
  pointwise strict failure relative to the initial pivot;
- `fermionicKernel_twoCornerResidual_not_half_contraction` proves that the
  universal half-contraction property on the physical rectangle at cutoff
  `log 13824` is false, including exact proofs that the witness lies inside
  the domain.

Result:

- the tempting all-cutoff extension of the cutoff-one two-corner theorem is
  formally closed as a false route;
- because earlier phases identify the two symmetric corners as the first two
  complete pivots, the failure concerns the actual canonical GECP trajectory,
  not an artificial pivot choice;
- the obstruction does not touch Conjecture G1's scale-dependent block length:
  it instead shows why intermediate, localized pivots are essential at large
  cutoff;
- the next promising target is the actual third-pivot geometry or a dyadic
  residual-dominance lemma that decomposes the large band into smaller secant
  tents.

Verification:

- the focused Lean build and expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- exact Python `Fraction` arithmetic reproduces the secant, residual, initial
  pivot, and margin, and confirms the same strict failure for every integer
  endpoint base from 24 through 64;
- a separate 100-decimal calculation forms the two nested Schur updates from
  kernel evaluations and agrees with the exact rational value and margin to
  `1e-90`;
- `./scripts/verify.sh` passes with 62 Python tests and no known dependency
  vulnerabilities.

## Phase AC — third-pivot frequency localization

State: complete with green CI on 2026-10-09; delivered in merged PR #31.

Phase AB shows that the two outer corners alone cannot halve every large-cutoff
residual. The Phase AA secant representation nevertheless implies a dyadic
localization mechanism. On the positive outer half-band, the two-corner
residual is bounded by the first-corner error
`exp(-t*omega)-exp(-t*Lambda)`. When `omega >= Lambda/2`, this is at most
`x*(1-x) <= 1/4` for `x=exp(-t*Lambda/2)`. Reflection gives the same bound on
the negative outer half-band. At the fixed center `(1/2,0)`, the residual is
strictly greater than `1/4` once `2*log 4 <= Lambda`. Therefore no complete
third pivot can remain in either outer half-band.

Frozen public Lean targets:

- `fermionicKernel_twoCornerResidual_le_firstCornerResidual` and
  `fermionicKernel_twoCornerResidual_le_secondCornerResidual`: expose the two
  exact comparison inequalities implicit in Phase AA;
- `fermionicKernel_secondCornerResidual_reflection`: identify the lower-corner
  comparison with the reflected upper-corner error;
- `fermionicKernel_firstCornerResidual_le_quarter_of_halfCutoff_le_frequency`:
  prove the `x*(1-x)` outer-half estimate;
- `fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_absFrequency`:
  combine both outer orientations for the actual symmetric two-corner
  residual;
- `fermionicKernel_symmetricTwoCornerResidual_center` and
  `fermionicKernel_symmetricTwoCornerResidual_center_gt_quarter`: compute the
  exact center value and prove the large-cutoff strict lower witness;
- `fermionicKernel_thirdCompletePivot_frequency_lt_halfCutoff`: prove that any
  third complete pivot dominating the center has `abs omega < Lambda/2`.

Proof contract:

- retain the actual nested `residualUpdate` and nonzero second-pivot witness;
- prove the positive and reflected outer estimates separately before using
  absolute frequency;
- derive the center comparison from an exact closed form, not a decimal bound;
- formulate complete-pivot localization through the defining domination of
  the center point, without assuming uniqueness or a tie-breaking rule.

Independent checks:

- evaluate the exact residual on deterministic high-precision grids bracketing
  the threshold `2*log 4` and verify that outer-half samples stay below one
  quarter while the center crosses above it;
- compare sampled maximizers with the proved central-half frequency region;
- run the focused Lean build, public axiom audit, and root verification.

Decision rule:

- if the center formula is the only blocked step after three attempts, retain
  the proved outer-half theorem and isolate the exact hyperbolic/exponential
  identity as the next missing lemma; do not replace strict localization with
  a sampled claim.

Non-claims:

- frequency localization does not identify the third pivot time, prove a
  unique maximizer, or prescribe tie-breaking between reflected maxima;
- one halving of the admissible frequency band does not yet prove iteration of
  the mechanism through later residuals or the full `O(s+1)` block;
- this phase does not prove Conjecture G1 or solve Simons Problem 4.2.

Delivered identifiers:

- `symmetricTwoCornerResidual` names the actual residual after the first two
  symmetric cutoff-corner pivots;
- `fermionicKernel_twoCornerResidual_le_firstCornerResidual` and
  `fermionicKernel_twoCornerResidual_le_secondCornerResidual` expose the two
  reflected comparison inequalities behind the Phase AA tent;
- `fermionicKernel_secondCornerResidual_reflection` identifies the lower
  comparison exactly with the reflected upper-corner residual;
- `fermionicKernel_firstCornerResidual_le_quarter_of_halfCutoff_le_frequency`
  proves the elementary `x*(1-x) <= 1/4` estimate on the positive outer
  half-band;
- the three `fermionicKernel_symmetricTwoCornerResidual_le_quarter...`
  theorems transfer that bound to each outer orientation and finally to
  `Lambda/2 <= abs omega`;
- `fermionicKernel_symmetricTwoCornerResidual_center` computes the fixed
  center as `1/2 - 1/(2*cosh(Lambda/2))`, while its `center_gt_quarter`
  corollary proves strict separation above the threshold `2*log 4`;
- `fermionicKernel_thirdCompletePivot_frequency_lt_halfCutoff` consumes the
  repository's actual `CompletePivotOn` predicate and proves
  `abs omega < Lambda/2` for every third complete pivot.

Result:

- the exact two-corner residual has a proved outer-versus-center gap, not just
  a sampled concentration pattern;
- the third pivot need not be unique and its time remains unspecified, but
  every valid complete-pivot choice is forced into the central half-frequency
  band;
- Phase AB's obstruction and Phase AC's localization fit together: the outer
  corners do not finish the contraction, but they eliminate both outer halves
  from the next complete-pivot search;
- the next formal obligation is an inheritance or residual-dominance theorem
  showing that the third update converts this positional localization into a
  reusable smaller-band problem.

Verification:

- the focused Lean build and expanded public axiom audit pass, with only
  Lean's permitted default axioms;
- a 100-decimal regression independently forms the nested Schur residual at
  five cutoffs beginning at `2*log 4`, verifies the exact center formula and
  outer-half quarter bound on 16,605 deterministic samples, and confirms each
  sampled maximizer lies in the central half-band;
- `./scripts/verify.sh` passes with 63 Python tests and no known dependency
  vulnerabilities.

## Phase AD — persistent outer-half exclusion

State: complete with green CI on 2026-10-09; delivered in merged PR #32.

Phase AC excludes both outer frequency half-bands from the third complete
pivot, but its quarter estimate is stated only for the two-corner residual.
The next reusable step is to show that a sign-coherent complete-pivot update
preserves any column-strip bound. The mechanism is asymmetric: the two
factors sampled inside the strip are bounded by the strip constant, while the
other two are bounded by the complete pivot. Same-sign cross products prevent
these mixed bounds from adding. Iteration then keeps the two-corner outer-half
quarter bound valid at every later residual.

Frozen public Lean targets:

- `residualUpdate_le_stripBound_of_signCoherent`: prove the one-step mixed
  strip/global estimate for a sign-coherent complete pivot;
- `stripBound_preserved_of_signCoherentCompletePivot`: iterate that estimate
  over an exact residual sequence;
- `fermionicKernel_outerHalfBound_preserved`: specialize the invariant to
  residuals realized by finite fermionic-kernel runs after the two symmetric
  cutoff corners;
- `fermionicKernel_laterCompletePivot_frequency_lt_halfCutoff_of_quarter_lt`:
  prove that every later complete pivot whose magnitude is strictly above one
  quarter has frequency strictly inside the central half-band.

Proof contract:

- use the existing `PivotCrossProductSignCoherent` and `CompletePivotOn`
  interfaces rather than introducing a fermionic-specific update identity;
- keep the strip estimate local in the frequency variable and use complete
  pivot maximality only for the off-strip column factor;
- require realization by a finite run from `fermionicKernel` so strict total
  positivity supplies sign coherence at every stage;
- obtain the initial absolute quarter bound from Phase AC together with the
  proved nonnegativity of the exact two-corner residual.

Independent checks:

- add a deterministic exact-arithmetic regression for the mixed-bound
  cancellation lemma, including both common-sign orientations;
- sample continued complete-pivot updates and verify that the outer-half
  maximum never exceeds the two-corner quarter bound;
- run the focused Lean build, public axiom audit, and root verification.

Decision rule:

- if direct sequence induction is blocked by index or run-realization
  bookkeeping after three attempts, ship the generic one-step invariant and
  isolate only the finite-run composition lemma as the next target; do not
  replace the realization hypothesis with an unproved preservation claim.

Non-claims:

- persistent outer-half exclusion does not prove contraction below one
  quarter, bound how many central pivots occur, or show self-similarity of the
  central restriction;
- this phase does not identify or pair later pivots, prove Conjecture G1, or
  solve Simons Problem 4.2.

Delivered identifiers:

- `residualUpdate_le_stripBound_of_signCoherent` proves the asymmetric
  one-step estimate with a protected-column bound and a pivot-sized global
  factor;
- `stripBound_preserved_of_signCoherentCompletePivot` iterates that estimate
  for an arbitrary sign-coherent complete-pivot residual sequence;
- `fermionicKernel_outerHalfBound_preserved` obtains coherence from realized
  fermionic runs and preserves the exact two-corner outer-half quarter bound;
- `fermionicKernel_laterCompletePivot_frequency_lt_halfCutoff_of_quarter_lt`
  converts the invariant into persistent frequency localization for every
  later selected pivot above one quarter.

Result:

- the third-pivot localization is no longer a one-step phenomenon: neither
  the third update nor any later exact complete-pivot update can recreate a
  residual larger than one quarter in the original outer half-bands;
- the proof does not assume a reflected pivot pair, uniqueness, or a
  tie-breaking rule;
- the next missing mechanism is relative rescaling on the central half-band:
  a new witness and protected-strip threshold must shrink with that band to
  turn the fixed invariant into a dyadic recurrence.

Verification:

- the focused Lean build and expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` cases cover both common-sign orientations and both pivot
  signs in the mixed-bound cancellation;
- an 80-decimal regression follows five further sampled complete-pivot
  updates at cutoffs 4, 8, and 16, checking the outer-half quarter bound and
  the location of every sampled pivot above it;
- `./scripts/verify.sh` passes with 64 Python tests and no known dependency
  vulnerabilities.

## Phase AE — exceptional-column determinant restart

State: complete with green CI on 2026-10-09; delivered in merged PR #33.

Phase AD shows that every continuation column selected while the pivot remains
above one quarter lies in the central half-band. Reusing the existing
smaller-cutoff separated approximation is blocked only by the two original
corner columns at `-Lambda` and `Lambda`. Exterior-power reasoning suggests
the correct repair: represent those exceptional columns exactly, augmenting
the approximating factorization rank by their cardinality while retaining an
error factor for every other column.

Frozen public Lean targets:

- `Matrix.abs_det_le_of_factors_approx_except`: extend the cardinality-sensitive
  low-rank determinant sum to a finite set of columns represented exactly;
- `Matrix.abs_det_le_two_pow_mul_of_factors_approx_except`: provide the closed
  coarse form with effective rank `card k + exceptional.card`;
- `fermionicKernel_sample_det_le_separatedApprox_except`: apply the explicit
  fermionic separated approximation only outside the exceptional set;
- `fermionicKernel_sample_det_le_two_pow_except`: expose the resulting closed
  determinant decay with the exact exceptional-column rank shift.

Proof contract:

- reuse `Matrix.mixedLeft`, `Matrix.mixedRight`, and `columnChoice` to encode
  exact exceptional columns rather than duplicating the factorization API;
- require cutoff membership only for nonexceptional frequencies; exceptional
  fermionic columns are bounded directly by the kernel's global unit bound;
- charge exactly `exceptional.card` additional factor coordinates and retain
  the same `2^-p` power for every column beyond the augmented rank;
- keep the theorem independent of how the exceptional columns were selected,
  so Phase AF can instantiate them with the two cutoff corners and Phase AD's
  localized continuation.

Independent checks:

- construct exact rational matrices consisting of a low-rank background plus
  two arbitrary exact columns and verify the predicted rank shift and
  determinant vanishing threshold;
- check the zero-exception specialization agrees with the existing bound;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if the augmented-factor cardinality cannot be normalized through the
  existing `Sum`/`Finset` APIs after three attempts, deliver the exact
  factorization identity and leave only its closed cardinal arithmetic as the
  next lemma; do not weaken exceptional columns to an unproved numerical
  approximation.

Non-claims:

- determinant decay with two exceptional columns does not itself prove that
  an actual GECP continuation is complete on the central band or bound the
  number of above-quarter pivots;
- this phase does not remove the logarithmic startup, prove Conjecture G1, or
  solve Simons Problem 4.2.

Delivered identifiers:

- `Matrix.abs_det_le_of_factors_approx_except` reuses the existing mixed
  factorization to make exact exceptional columns part of the low-rank
  background;
- `Matrix.abs_det_le_two_pow_mul_of_factors_approx_except` gives the closed
  determinant bound with effective rank `card k + exceptional.card`;
- `fermionicKernel_sample_det_le_separatedApprox_except` applies the explicit
  separated approximation only to nonexceptional frequency columns;
- `fermionicKernel_sample_det_le_two_pow_except` exposes the corresponding
  closed fermionic decay estimate.

Result:

- the two original cutoff-corner columns can be retained exactly while every
  later central-half column uses the cutoff-halved separated approximation;
- exceptional frequencies need no cutoff hypothesis and cost exactly two
  rank coordinates when the exceptional set contains the two corners;
- the next formal obligation is to identify those two columns inside the
  appended GECP selected core and combine the shifted determinant estimate
  with the continuation pivot product.

Verification:

- the focused Lean build and expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` matrices independently verify determinant expansion,
  effective-rank vanishing, and the coarse error exponent for both zero and
  two exceptional columns;
- `./scripts/verify.sh` passes with 65 Python tests and no known dependency
  vulnerabilities.

## Phase AF — GECP exceptional-column composition

State: complete with green CI on 2026-10-09; delivered in merged PR #34.

Phase AE proves determinant decay for a sample matrix with finitely many exact
out-of-band columns. The next bridge is to apply that theorem to the selected
core of an actual realized fermionic GECP prefix. The complete-pivot sequence
may run on a larger dyadic rectangle, while the determinant approximation uses
a smaller scale for every selected column outside an explicit exceptional
index set.

Frozen public Lean targets:

- `fermionicKernel_gecp_error_pow_le_separatedApprox_except`: compose the
  selected-core residual-power inequality with Phase AE's cardinality-
  sensitive determinant sum;
- `fermionicKernel_gecp_error_pow_le_two_pow_except`: expose the closed power
  bound with effective rank `16*p*(s+1) + exceptional.card`, while complete
  pivoting and the evaluated residual may live at a larger scale.

Proof contract:

- use the actual `Run.SelectedIndex`, `Run.selectedCore`, and
  `Run.RealizesPivotPrefix` interfaces;
- keep the complete-pivot domain scale independent from the smaller separated-
  approximation scale;
- require physical time membership for every selected row and smaller-band
  frequency membership only for nonexceptional selected columns;
- retain Phase T's exact residual-power-to-selected-determinant inequality and
  Phase AE's exact exceptional cardinality shift without replacing either by
  a sampled surrogate.

Independent checks:

- use exact rational GECP matrices to verify that two out-of-band selected
  columns are charged as two rank coordinates and that the pivot-product/
  selected-core identity remains exact;
- verify that an empty exceptional set reduces arithmetically to the existing
  GECP determinant exponent;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if dependent selected-index rewriting blocks the closed theorem after three
  attempts, deliver the cardinality-sensitive selected-core theorem and
  isolate only the normalization of `Fintype.card run.SelectedIndex` as the
  next lemma; do not replace the actual selected core by an unrelated matrix.

Non-claims:

- this phase does not yet prove that the exceptional set consists of exactly
  the first two cutoff-corner indices or derive smaller-band membership from
  the Phase AD pivot threshold;
- it does not prove a shrinking-band recurrence, remove the logarithmic
  startup, prove Conjecture G1, or solve Simons Problem 4.2.

Delivered identifiers:

- `fermionicKernel_gecp_error_pow_le_separatedApprox_except` composes the
  actual selected-core residual-power inequality with Phase AE's exact
  exceptional-column determinant sum;
- `fermionicKernel_gecp_error_pow_le_two_pow_except` gives the closed form at
  an approximation scale independent of the larger complete-pivot domain.

Result:

- an actual realized fermionic GECP prefix can now use a cutoff-halved
  separated approximation for all but explicitly named selected columns;
- the residual-power exponent retains the exact exceptional-set cardinality
  rather than charging the larger complete-pivot scale;
- the next obligation is dependent-index bookkeeping: construct the two-
  element exceptional set for the initial corner steps and prove every other
  selected column lies in the central half while the relevant pivots exceed
  one quarter.

Verification:

- the focused Lean build and expanded public axiom audit pass with only
  Lean's permitted default axioms;
- an exact full-rank rational GECP run independently verifies the selected-
  core determinant/pivot-product identity and the two-exception rank shift;
- `./scripts/verify.sh` passes with 66 Python tests and no known dependency
  vulnerabilities.

## Phase AG — two-corner exceptional indices

State: complete with green CI on 2026-10-09; delivered in merged PR #35.

Phase AF accepts an arbitrary exceptional subset of an actual GECP selected
core. The canonical application needs a dependent run whose first two steps
are the prescribed symmetric cutoff corners and a concrete exceptional set
containing exactly those recursive selected indices. Every other index then
belongs definitionally to the continuation run, so central-band completeness
of the continuation supplies the smaller-cutoff hypotheses required by Phase
AE without list-position reconstruction.

Frozen public Lean targets:

- `symmetricTwoCornerRun`: prepend the exact pivots `(0,Lambda)` and
  `(1,-Lambda)` to any successful continuation of
  `symmetricTwoCornerResidual`;
- `symmetricTwoCornerExceptional`: define the two first recursive selected
  indices of that run;
- `symmetricTwoCornerExceptional_card`: prove the exceptional set has exactly
  two elements;
- `symmetricTwoCornerRun_selectedColumn_mem_of_not_exceptional`: transfer
  continuation column-domain membership to every nonexceptional full-run
  selected index;
- `fermionicKernel_symmetricTwoCornerRun_sample_det_le_two_pow`: instantiate
  Phase AE's determinant estimate with the two exact corner columns and a
  smaller-band complete continuation.

Proof contract:

- retain the exact nonzero witnesses already used by
  `symmetricTwoCornerResidual` so the run's final residual is definitionally
  the continuation kernel;
- define the exceptional set in the recursive `Fin 1 + Fin 1 + tail` index
  type, not by an assumed list ordering;
- prove nonexceptional membership by dependent sum elimination and the
  existing `Run.selectedColumn_mem_of_completeOn` theorem;
- charge exactly two exceptional coordinates in the closed determinant
  exponent.

Independent checks:

- run an exact rational prescribed two-pivot prefix followed by a central
  continuation and verify selected-column order, two exceptions, determinant,
  and pivot product;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if `Finset` simplification over the dependent recursive index fails after
  three attempts, expose the equivalent two-constructor membership predicate
  and prove its finite cardinality separately; do not assume the exception
  count or selected-column order.

Non-claims:

- continuation completeness on the central band remains a hypothesis here;
  deriving it from Phase AD's above-quarter condition is the next phase;
- this phase does not prove a shrinking-band recurrence, remove the
  logarithmic startup, prove Conjecture G1, or solve Simons Problem 4.2.

Delivered identifiers:

- `symmetricTwoCornerRun` embeds the exact prescribed corner steps before an
  arbitrary successful continuation;
- `symmetricTwoCornerExceptional` and
  `symmetricTwoCornerExceptional_card` identify exactly two recursive selected
  indices without converting through lists;
- `symmetricTwoCornerRun_selectedRow_mem` and
  `symmetricTwoCornerRun_selectedColumn_mem_of_not_exceptional` transfer the
  physical row and continuation-column domains to the full run;
- `fermionicKernel_symmetricTwoCornerRun_sample_det_le_two_pow` specializes
  Phase AE with the exact `+2` effective-rank shift.

Result:

- the two cutoff corners are now formally isolated as the only selected
  columns exempt from a smaller-band approximation;
- no ordering assumption or cardinality hypothesis is used for the exception
  set: both follow from the dependent run constructors;
- the next obligation is to build a finite continuation from the Phase AD
  residual sequence and prove its `CompleteOn` property on the central band
  while its pivots remain above one quarter.

Verification:

- the focused Lean build and expanded public axiom audit pass with only
  Lean's permitted default axioms;
- an exact rational prescribed-corner prefix followed by a central complete-
  pivot continuation verifies selected-column order, the two exceptions,
  selected-core determinant, and pivot product;
- `./scripts/verify.sh` passes with 67 Python tests and no known dependency
  vulnerabilities.

## Phase AH — finite central-prefix runs

State: complete with green CI on 2026-10-09; delivered in merged PR #36.

Phase AD localizes each later pivot above one quarter to the central frequency
half-band, while Phase AG consumes a continuation whose complete pivots are
already known to lie in that band. The missing bridge is structural: assemble
a finite dependent `Run` from the indexed residual recurrence, and lift the
pointwise localization/maximality facts to `Run.CompleteOn` for every finite
above-quarter prefix.

Frozen public Lean targets:

- `Run.ofResidualSequenceFrom`: construct a successful run of any requested
  finite length from a residual recurrence and its nonzero pivot witnesses;
- `Run.finalResidual_ofResidualSequenceFrom`: identify its final residual with
  the indexed residual at `start + length`;
- `Run.pivots_ofResidualSequenceFrom` and
  `Run.card_selectedIndex_ofResidualSequenceFrom`: identify the consecutive
  pivot values and exact dependent selected-index cardinality;
- `Run.completeOn_ofResidualSequenceFrom`: assemble bounded pointwise
  `CompletePivotOn` witnesses into dependent `Run.CompleteOn` evidence;
- `fermionicKernel_laterCompletePivot_completeOn_centralHalf_of_quarter_lt`:
  restrict a later full-cutoff complete pivot above one quarter to the central
  half-band;
- `fermionicKernel_aboveQuarterPrefix_completeOn_centralHalf`: package every
  finite above-quarter continuation prefix as central-half complete.
- `symmetricTwoCornerResidualSequence_realized`: derive the formerly separate
  realization hypothesis from the recurrence itself by prepending the exact
  two-corner run.
- `symmetricTwoCornerPrefixRun_card`: expose the exact `length + 2` cardinality
  needed by Phase AG's determinant bound.
- `fermionicKernel_aboveQuarterPrefix_sample_det_le_two_pow`: discharge Phase
  AG's continuation-domain and cardinality hypotheses from the indexed
  above-quarter trajectory, completing the Phase AD-to-AG composition.

Proof contract:

- transport the recursively constructed tail only across the stated residual
  update equality; do not assume definitional equality of indexed residuals;
- prove the final-residual and completeness transport lemmas by equality
  elimination, so no proof-irrelevance or list-position assumption enters the
  public statements;
- obtain central frequency membership from Phase AD and obtain maximality on
  the smaller band solely by restricting the existing full-band complete
  pivot witness;
- keep the theorem finite and threshold-conditional: the first pivot at or
  below one quarter ends the prefix and is not silently discarded.
- do not retain an independent realization assumption in the final prefix
  theorem when the constructed dependent run already proves it.

Independent checks:

- exactly execute a finite residual recurrence and verify the reconstructed
  prefix length, pivot order, final residual, and central-band restriction;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if dependent equality transport becomes brittle after three proof attempts,
  expose a private cast helper proved by `cases` on the update equality; do not
  weaken the run constructor to a list or add an equality assumption containing
  the desired conclusion.

Non-claims:

- the phase does not prove that the above-quarter prefix has a prescribed
  length, nor that the first following pivot is at most one quarter;
- it does not yet iterate the central-half restart, remove the logarithmic
  startup, prove Conjecture G1, or solve Simons Problem 4.2.

Delivered identifiers:

- `Run.ofResidualSequenceFrom` constructs a dependent run from an indexed
  exact update recurrence, while its final-residual, pivot-list, and selected-
  index-cardinality theorems recover the exact requested prefix data;
- `Run.completeOn_ofResidualSequenceFrom` converts bounded pointwise complete-
  pivot evidence into recursive `Run.CompleteOn` evidence;
- `symmetricTwoCornerContinuationPrefix` transports that run to the exact
  symmetric two-corner residual, and
  `symmetricTwoCornerResidualSequence_realized` prepends the prescribed
  corners to realize every continuation residual from the original kernel;
- `fermionicKernel_aboveQuarterPrefix_completeOn_centralHalf` turns Phase AD's
  per-step localization into central-half completeness for the whole prefix;
- `fermionicKernel_aboveQuarterPrefix_sample_det_le_two_pow` discharges Phase
  AG's continuation-domain, realization, and cardinality bookkeeping and gives
  the closed smaller-band determinant bound with exact `length + 2` order.

Result:

- Phase AD and Phase AG now compose on the actual dependent GECP trajectory,
  rather than through a separately hypothesized continuation;
- the remaining quantitative task is a stopping argument: compare the lower
  determinant product forced by an all-above-quarter prefix with the new upper
  determinant estimate, choose the approximation order, and prove that such a
  prefix cannot exceed an explicit scale-dependent length;
- this is still not a shrinking-band recurrence below the first quarter-scale
  stopping point, Conjecture G1, or a solution of Problem 4.2.

Verification:

- the focused Lean build and expanded public axiom audit pass with only Lean's
  permitted default axioms;
- an exact rational recurrence reconstructs a three-step above-quarter prefix,
  its complete pivots, central columns, and final residual without rounding;
- `./scripts/verify.sh` passes with 68 Python tests and no known dependency
  vulnerabilities.

## Phase AI — above-quarter stopping bound

State: complete with green CI on 2026-10-09; delivered in merged PR #37.

Phase AH gives an upper bound for the selected-core determinant of any finite
continuation prefix whose pivots all remain above one quarter. The same
determinant is exactly the product of the two prescribed corner pivots and the
continuation pivots. A volume-sandwich argument can therefore rule out a long
all-above-quarter prefix once an explicit dyadic approximation order makes the
upper determinant smaller than the product lower bound.

Frozen public Lean targets:

- lower-bound the two prescribed corner pivots uniformly at dyadic cutoffs;
- identify the absolute selected-core determinant of the constructed full run
  with the product of its absolute pivots and derive a strict quarter-power
  lower bound from an all-above-quarter continuation;
- prove a reusable arithmetic criterion that bounds Phase AH's determinant
  upper estimate by the same quarter power using `n! <= n^n`, a binary size
  bound, and an exponent-budget inequality;
- instantiate the criterion with explicit polynomial parameters and conclude
  that some continuation pivot in the certified prefix is at most one quarter.
- transfer that stopping-pivot bound through complete-pivot maximality to a
  global quarter bound for the corresponding residual.

Proof contract:

- derive the determinant lower bound from actual stored pivots, not a fresh
  nonsingularity or volume assumption;
- keep the two strict corner inequalities explicit so the contradiction
  remains strict even for a zero-length continuation product;
- separate analytic, determinant, and natural-number exponent arithmetic into
  independently reusable lemmas;
- state the resulting stopping index explicitly and preserve complete-pivot
  tie semantics (`<= 1/4`, not an assumed strict drop).

Independent checks:

- evaluate the proposed parameter family with exact integers across ordinary
  scales and powers-of-two transition neighborhoods;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if the first explicit polynomial constants make the exponent budget awkward
  after three attempts, enlarge them and retain the clean proof; optimization
  belongs after the stopping mechanism is certified.

Non-claims:

- this phase forces the first quarter-scale stopping event but does not yet
  turn it into a repeated shrinking-band recurrence;
- it does not remove the logarithmic startup in the existing global rate,
  prove Conjecture G1, or solve Simons Problem 4.2.

Delivered identifiers:

- `Run.abs_selectedCore_det_eq_prod_abs_pivots` exposes the absolute pivot
  product on the recursively indexed selected core;
- `fermionicKernel_cutoffCorner_gt_half` and
  `fermionicKernel_firstPivot_reflected_gt_half` give strict uniform lower
  bounds for the two prescribed corner pivots at cutoff at least one;
- `fermionicKernel_aboveQuarterPrefix_sample_det_gt_quarter_pow` proves the
  strict determinant lower bound `(1/4)^(length+1)` from an actual all-above-
  quarter trajectory prefix;
- `aboveQuarter_determinant_upper_le_quarter_pow` reduces the determinant
  contradiction to a binary size bound and a natural-number exponent budget;
- `fermionicKernel_exists_pivot_le_quarter_of_exponent_budget` packages the
  reusable stopping criterion;
- `aboveQuarterStoppingOrder`, `aboveQuarterStoppingLength`, and
  `aboveQuarterStoppingBits` instantiate it with `24(s+1)`, `2048(s+1)^2`,
  and `2(s+1)+10`;
- `fermionicKernel_exists_residual_le_quarter_explicit` concludes that at some
  continuation index below `2048(s+1)^2`, complete-pivot maximality bounds the
  entire cutoff-domain residual by one quarter.

Result:

- the Phase AD localization is no longer merely conditional for arbitrarily
  long prefixes: the Phase AH volume upper bound and exact pivot-product lower
  bound force a finite quarter-scale stopping event;
- the constants are intentionally conservative and the quadratic scale cost
  is not competitive with the earlier global logarithmic block, but the new
  theorem tracks the localized dependent continuation and is suitable for a
  renormalized or nested-band argument;
- the next structural obligation is to generalize the fixed outer-quarter
  strip invariant to a threshold/band family that can restart after this
  stopping event, rather than proving only the first absolute threshold.

Verification:

- the focused Lean build and expanded public axiom audit pass with only Lean's
  permitted default axioms;
- exact integer regression checks the rank, binary-size, and exponent budgets
  through scale 256 and neighborhoods of powers of two through `2^16`;
- `./scripts/verify.sh` passes with 69 Python tests and no known dependency
  vulnerabilities.

## Phase AJ — multi-bordered residual minors

State: complete with green CI on 2026-10-10; delivered in merged PR #38.

The current bordered determinant API handles one appended row and column,
which is exactly enough for pointwise residual identities and order-two cross
sign coherence. A renormalization argument needs the stronger Schur-complement
closure fact: every finite minor of a run's final residual is the corresponding
multi-bordered original-kernel determinant divided by the selected-core
determinant. This is the all-orders analogue of the one-border identity.

Frozen public Lean targets:

- define a recursive arbitrary-border index, row map, column map, and augmented
  core that preserve the run's dependent pivot order;
- prove the one-step augmented determinant factorization by exact block
  elimination, without introducing an inverse;
- prove the full identity `det augmented = det selectedCore * det residualMinor`;
- transfer nonsingularity from an injective augmented original-kernel minor to
  the corresponding final-residual minor;
- specialize the transfer to strictly sign-regular kernels and the fermionic
  kernel at every finite order.

Proof contract:

- recurse on `Run`, retaining the selected coordinates before all appended
  border coordinates; do not flatten through lists;
- use factorization/elimination rather than an explicit matrix inverse;
- keep arbitrary finite border index types in the generic algebraic theorem;
- require injectivity of the combined augmented coordinates explicitly when
  invoking strict sign regularity, so repeated selected/border nodes correctly
  remain a zero-minor obstruction.

Independent checks:

- verify the multi-border identity exactly for rational runs with border sizes
  zero through three;
- run the focused Lean build, expanded public axiom audit, and root
  verification.

Decision rule:

- if transporting an arbitrary border through the recursive sum index becomes
  brittle after three attempts, introduce a dedicated recursive
  `AugmentedIndex` rather than reindexing `(SelectedIndex + border)` at every
  step.

Non-claims:

- this phase establishes the determinant mechanism for all residual orders;
  it does not yet prove a globally fixed residual signature after arbitrary
  pivot interleavings;
- it does not by itself create the threshold-parametric nested-band restart,
  prove Conjecture G1, or solve Simons Problem 4.2.

Delivered identifiers:

- `Run.selectedCore_det_ne`, `Run.selectedRow_injective`, and
  `Run.selectedColumn_injective` expose the nonsingularity and coordinate
  distinctness already forced by a valid dependent run;
- `Run.AugmentedIndex`, `Run.augmentedRow`, `Run.augmentedColumn`, and
  `Run.augmentedCore` retain every pivot before an arbitrary finite border;
- `Run.augmentedCore_step_det` proves the exact one-step block factorization,
  and `Run.augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor`
  iterates it through the full run;
- `GECP.strictSignRegular_finalResidual_minor_oriented_pos` transports the
  strict sign of an injective augmented minor through the selected-core and
  explicit row/column orientation gauges;
- `GECP.strictSignRegular_finalResidual_minor_ne_zero` and the fermionic
  specializations in `Fermionic.ResidualMinors` prove all corresponding
  residual minors nonzero.

Result:

- the former one-border/order-two mechanism now holds for every finite minor
  order without an inverse or a fixed-size reindexing;
- sign regularity is preserved exactly after accounting for the permutation
  gauge created by interleaving new border coordinates with the run's pivot
  order;
- the next useful target is to factor those tuple orientations into pointwise
  row and column gauges. That would package a gauge-scaled final residual as a
  strictly sign-regular kernel and make recursive band restarts substantially
  easier to state and compose.

Verification:

- focused generic and fermionic Lean builds and the expanded public axiom
  audit pass with only Lean's permitted default axioms;
- exact `Fraction` arithmetic checks the determinant identity and nonzero
  residual minors for border sizes zero through three;
- `./scripts/verify.sh` passes with 70 Python tests and no known dependency
  vulnerabilities.

## Phase AK — balanced-sign residual gauges

State: complete with green CI on 2026-10-10; delivered in merged PR #39.

Phase AJ identifies the exact orientation carried by every residual minor, but
the tuple-level sorting signs are awkward for recursive analysis. The order-one
part has a cleaner interpretation from signed bipartite graph theory: Phase P's
cross-product sign coherence says the nonzero entry signs form a balanced
bipartite signing, so one anchor row and column determine diagonal row/column
sign gauges that make every remaining entry positive.

Frozen public Lean targets:

- define anchor-based row and column weights and their magnitude-one sign
  gauges for an arbitrary real kernel;
- prove that cross-product sign coherence plus four nonzero entries makes the
  sign-gauged entry strictly positive;
- prove that sign gauging preserves every pointwise absolute value and expose
  the exact determinant scaling for arbitrary finite minors;
- define fresh row and column coordinates for a dependent run and derive the
  required nonvanishing from Phase AJ;
- specialize the positive gauge normalization to every fermionic final
  residual on fresh coordinates.

Proof contract:

- derive sign balance from the existing cross-product theorem rather than
  reproving permutation parity;
- use only diagonal sign scalings, so complete-pivot magnitudes are unchanged;
- retain explicit freshness assumptions because selected rows and columns are
  exact residual zeros;
- do not claim all-orders strict total positivity until the compound-minor
  gauges are also factored.

Independent checks:

- verify exact gauge positivity and absolute-value preservation on rational
  geometric-surrogate residuals for every unselected anchor and entry;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if direct `Real.sign` algebra becomes brittle after three attempts, prove the
  raw nonzero diagonal-scaling statement first and recover the magnitude-one
  version from `Real.sign_mul_pos_of_ne_zero`.

Non-claims:

- entrywise positive normalization is the order-one checkerboard closure, not
  yet an all-orders strictly totally positive residual theorem;
- this phase does not itself sharpen the determinant prefactor, construct a
  linear-size dyadic contraction block, prove Conjecture G1, or solve Simons
  Problem 4.2.

Delivered identifiers:

- `anchoredRowWeight`, `anchoredColumnWeight`, and `anchoredGaugeKernel`
  express the balanced bipartite signing as diagonal scaling from one anchor;
- `anchoredSignGaugeKernel` replaces those weights by magnitude-one signs;
- `anchoredGaugeKernel_pos_of_crossProductSignCoherent` and
  `anchoredSignGaugeKernel_pos_of_crossProductSignCoherent` prove strict
  positivity from cross-product coherence and four nonzero entries;
- `abs_anchoredSignGaugeKernel` proves pointwise magnitude preservation, while
  `anchoredSignGaugeKernel_minor_det` records the exact row/column gauge
  products multiplying every finite minor;
- `Run.FreshRow`, `Run.FreshColumn`, and
  `strictSignRegular_finalResidual_ne_zero_of_fresh` package the Phase AJ
  nonvanishing premise at order one;
- `fermionicKernel_finalResidual_anchoredSignGauge_pos` gives the positive
  normalization for every fermionic final residual and every fresh anchor and
  entry.

Result:

- the order-one checkerboard ambiguity is completely removed by diagonal
  signs of absolute value one, so exact complete-pivot magnitudes are
  unchanged;
- the construction avoids explicit sorting-permutation arithmetic and follows
  the balanced-sign characterization familiar from signed bipartite graphs;
- the next target is the exterior-power analogue: prove that each fixed-order
  compound kernel of residual minors is itself cross-product sign coherent,
  then apply the same anchor-gauge theorem at the compound level. This should
  expose all-orders positivity without forcing one pointwise gauge to encode
  every tuple permutation at once.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic checks every unselected anchor and entry of a
  three-pivot geometric-surrogate residual for positivity and magnitude
  preservation;
- `./scripts/verify.sh` passes with 71 Python tests and no known dependency
  vulnerabilities.

## Phase AL — compound residual gauges

State: complete with green CI on 2026-10-10; delivered in merged PR #40.

Phase AK solves the sign normalization problem for individual residual
entries. Exterior algebra suggests the correct all-orders lift: for each
minor order, regard residual determinants as the entries of a compound kernel
whose row and column vertices are fresh tuples. Phase AJ's orientation formula
then makes the compound signing balanced because four row/column orientation
factors occur in squares.

Frozen public Lean targets:

- define fresh row-tuple and column-tuple types for an arbitrary finite minor
  order and a dependent run;
- define the final-residual compound kernel whose entries are the corresponding
  minors;
- prove every compound entry nonzero under strict sign regularity;
- prove the compound kernel is cross-product sign coherent by multiplying four
  Phase AJ orientation-aware inequalities and cancelling squared signs;
- apply Phase AK's anchor gauge to obtain a positive, magnitude-preserving
  normalization of every fixed-order fermionic compound kernel.

Proof contract:

- keep the minor order as an arbitrary finite index type rather than reduce to
  lists or fixed natural sizes;
- use the existing all-orders residual identity and orientation theorem; do
  not reprove Schur complement algebra;
- expose strict positivity of the four-minor compound product, not merely
  nonnegativity, because freshness supplies nonvanishing;
- do not claim that the independent gauge at each compound order is induced by
  one base-kernel row/column gauge.

Independent checks:

- for exact rational geometric-surrogate residuals, build compound matrices at
  orders one through three and verify every anchor gauge, magnitude, and
  four-entry product;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if the four-factor sign algebra becomes brittle after three attempts, first
  package the common signed determinant factor and prove its square is one,
  then finish with a separate ordered-ring lemma.

Non-claims:

- compound positivity is a sign-structure theorem; it does not yet provide a
  quantitative determinant decay estimate or a linear-size contraction block;
- this phase does not prove Conjecture G1 or solve Simons Problem 4.2.

Delivered identifiers:

- `Run.FreshRows` and `Run.FreshColumns` package arbitrary finite tuples whose
  augmentation by the selected coordinates remains injective;
- `Run.finalResidualMinorKernel` treats fixed-order residual determinants as a
  kernel on those tuple types;
- `strictSignRegular_finalResidualMinorKernel_ne_zero` proves every compound
  entry nonzero;
- `strictSignRegular_finalResidualMinorKernel_crossProductSignCoherent`
  multiplies four Phase AJ signed-minor inequalities and cancels the squared
  orientation factors, proving balanced signs at every order;
- `strictSignRegular_finalResidualMinorKernel_anchoredSignGauge_pos` and its
  absolute-value companion apply Phase AK to expose a positive compound
  kernel without changing any minor magnitude;
- the fermionic cross-product and positive-gauge corollaries instantiate the
  result for every successful fermionic run.

Result:

- the order-one balanced-sign theorem now holds on every exterior power of a
  final residual, with the minor order represented by an arbitrary finite
  index type;
- all fresh fermionic minors lie in a single positive sign chamber after an
  independently chosen magnitude-one compound gauge at each order;
- this supplies the sign premise needed for subtraction-sensitive minor
  identities. The next quantitative target is a Desnanot--Jacobi/Dodgson
  condensation inequality for these positive compound coordinates, seeking a
  recursive determinant bound without the generic factorial prefactor.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic constructs compound matrices of orders one
  through three after a three-pivot geometric-surrogate run and checks every
  four-entry product, anchor gauge, and preserved magnitude;
- `./scripts/verify.sh` passes with 72 Python tests and no known dependency
  vulnerabilities.

## Phase AM — two-border condensation inequality

State: complete; delivered in merged PR #41 with green CI.

Phase AL supplies positive sign chambers for compound minors. The first
quantitative identity to exploit is Desnanot--Jacobi condensation: a selected
core, its four one-border extensions, and a two-border extension satisfy an
exact quadratic relation. Strict sign regularity makes the two products on the
right have the same sign, replacing the generic triangle factor two by a
factor-one maximum bound.

Frozen public Lean targets:

- prove a generic same-sign inequality
  `abs (a-b) <= max (abs a) (abs b)`;
- derive the exact two-border selected-core condensation identity directly
  from the Phase I and Phase AJ determinant factorizations;
- combine strict sign regularity with bordered-minor sign coherence to prove a
  factor-one absolute condensation inequality;
- specialize the result to every fermionic run;
- state explicitly which additional maximal-minor or approximation estimate is
  still required before condensation yields a rate improvement.

Proof contract:

- use `Fin 2` borders and `Matrix.det_fin_two`; do not introduce inverses or a
  general adjugate API merely for this identity;
- preserve arbitrary row and column coordinates, including duplicate cases;
- derive the sign premise from the existing bordered-minor theorem rather than
  from numerical positivity;
- do not claim that factor-one condensation alone removes the determinant
  prefactor.

Independent checks:

- verify the exact identity and factor-one bound over every two-border choice
  of an exact rational geometric-surrogate run;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if direct rewriting through the bordered identities fails three times,
  first prove the residual-level `2 x 2` determinant formula and compose it
  with the selected-core factors in a separate theorem.

Non-claims:

- this phase packages the first subtraction-sensitive consequence of compound
  sign balance; it does not yet prove a dimension-free approximation bound,
  Conjecture G1, or Simons Problem 4.2.

Delivered:

- `GECP.orderedPair` provides the displayed `Fin 2` border order used by the
  exact identity;
- `GECP.abs_sub_le_max_of_mul_nonneg` replaces the generic triangle factor two
  by a factor-one maximum whenever the two products have the same sign;
- `Run.augmentedCore_fin_two_mul_selectedCore_det` proves the selected-core
  Desnanot--Jacobi identity from the existing one-border and arbitrary-border
  Schur factorizations, without an inverse;
- `strictSignRegular_augmentedCore_fin_two_condensation` derives the same-sign
  premise from bordered-minor coherence and proves the factor-one absolute
  bound, with
  `fermionicKernel_augmentedCore_fin_two_condensation` specializing it to every
  fermionic run.

Result:

- compound sign balance has a genuine quantitative payoff: subtraction in the
  two-border identity costs the larger product rather than the sum of both
  products;
- the estimate is scale-exact but homogeneous. It propagates a bound already
  known for the four one-border minors and therefore does not itself create the
  missing decay or remove the determinant prefactor;
- the next useful target must control a maximal one-border-minor envelope or
  prove stable interpolation/approximation on a sign chamber. Repackaging
  condensation at higher order without such a magnitude input would be
  structurally correct but quantitatively circular.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic checks the identity, sign premise, and
  factor-one bound for all 256 ordered two-row/two-column fresh-border choices
  after a three-pivot geometric-surrogate run;
- `./scripts/verify.sh` passes with 73 Python tests and no known dependency
  vulnerabilities.

## Phase AN — factorial-free continuation determinants

State: complete; delivered in merged PR #42 with green CI.

Phase AM shows that compound sign balance removes a local triangle factor, but
also shows why condensation alone is homogeneous and cannot manufacture
decay. The next useful abstraction comes from total-positive elimination: if a
strictly sign-regular residual is already uniformly bounded by `B`, every
complete continuation pivot remains bounded by `B`, so its selected-core
determinant is bounded by `B^m` with no Leibniz factorial. This is the exact
determinant estimate needed if the dyadic approximation remainder can be shown
to retain compatible sign regularity.

Frozen public Lean targets:

- express a finite continuation whose current residual is realized by a prefix
  of the original strictly sign-regular kernel;
- prove that a domain-complete continuation preserves an initial uniform bound
  at every recursive step using selected-cross sign coherence;
- bound every continuation pivot by that same envelope;
- combine the pivot-product identity with the pointwise bounds to prove
  `abs(det selectedCore) <= B ^ card` without `card.factorial`;
- specialize the result to fermionic prefixes and cutoff-domain continuations.

Proof contract:

- use the existing dependent `Run`, `Run.append`, complete-pivot, and
  sign-coherent nonexpansiveness APIs; do not introduce an inverse;
- keep the theorem finite and continuation-local rather than encoding an
  artificial infinite residual sequence;
- make nonnegativity of `B` explicit and retain the exact selected-core
  exponent;
- do not claim strict decay: the result converts a small structured remainder
  envelope into a determinant power but does not prove that the current
  fermionic dyadic remainder satisfies the required compatibility.

Independent checks:

- use exact `Fraction` arithmetic on geometric sign-regular surrogates to check
  every suffix determinant against the appropriate residual sup power;
- include an exact witness that complete-pivot interpolation coefficients need
  not be bounded by one, ruling out that stronger shortcut;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if dependent prefix transport blocks direct recursion three times, isolate a
  reusable bounded-continuation predicate and prove determinant control from
  that predicate before reconnecting strict sign regularity.

Non-claims:

- this phase does not prove total positivity of the dyadic Taylor remainder,
  compatibility of mixed approximation/remainder columns, Conjecture G1, or
  Simons Problem 4.2.

Delivered:

- `Run.SignCoherent` and `Run.PivotsBounded` record the recursive hypotheses
  needed by a finite dependent run without inventing an infinite sequence;
- `Run.pivotsBounded_of_signCoherent_completeOn` proves by recursive Schur
  updates that complete sign-coherent pivoting preserves the initial uniform
  envelope at every step;
- `Run.abs_selectedCore_det_le_pow_of_pivotsBounded` combines those pointwise
  bounds with the exact pivot-product identity, producing `B^m` with no
  `m.factorial`;
- `strictSignRegular_continuation_signCoherent` transports the original
  all-orders sign theorem through an arbitrary successful prefix and finite
  continuation;
- `strictSignRegular_completeContinuation_pivotsBounded` and
  `strictSignRegular_completeContinuation_selectedCore_det_le_pow` close the
  generic result, while
  `fermionicKernel_completeContinuation_selectedCore_det_le_pow` supplies the
  fermionic specialization.

Result:

- total-positive elimination does remove the generic determinant factorial
  once the matrix being eliminated is itself a complete sign-coherent
  continuation with a small entry envelope;
- this is precisely the `epsilon^m` mechanism needed for a structured
  approximation remainder, but the current dyadic decomposition has not yet
  been proved compatible with that continuation structure;
- a stronger shortcut is false: an exact `8 x 8` geometric surrogate at rank
  six has a row-replacement interpolation coefficient
  `6243374306 / 4938550965 > 1`, so complete pivoting plus strict sign
  regularity does not make the selected core 1-dominant;
- high-precision diagnostics indicate that even Taylor remainders retain the
  expected minor signs in the approximation regime. This is consistent with
  neighboring beta--gamma total-positivity results for integer parameters
  ([Simon, 2012](https://arxiv.org/abs/1207.6464)), and makes remainder
  sign-regularity plus mixed-column compatibility the next disciplined target.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic checks every suffix length after every prefix of
  an eight-point geometric-surrogate GECP run against the corresponding
  residual-supremum power;
- the same exact regression certifies the explicit failure of factor-one core
  dominance;
- `./scripts/verify.sh` passes with 75 Python tests and no known dependency
  vulnerabilities.

## Phase AO — masked dyadic remainder obstruction

State: complete; delivered in merged PR #43 with green CI.

Phase AN removes the determinant factorial if the small approximation
remainder itself supplies compatible sign-coherent complete elimination. The
unmasked even Taylor tail has strong numerical and neighboring theoretical
evidence for the expected total-positive signs, but the approximation actually
used by this repository contains discontinuous dyadic masks. Before attempting
an all-orders tail theorem, this phase tests the first scale transition exactly.

Frozen public Lean targets:

- define the positive-frequency error of the existing masked dyadic Taylor
  approximation;
- reduce its values at `t = (1/2, 3/5)` and `omega = (2, 5/2)` to three
  eighth-order Taylor remainders and one masked full exponential;
- prove by explicit exponential-series bounds that the resulting `2 x 2`
  determinant is positive;
- conclude that the masked remainder is not strictly sign regular at order two
  with `expKernelSignature 2`.

Proof contract:

- prove the sign with rational inequalities and mathlib exponential bounds;
  floating-point output may guide constants but is not the theorem source;
- reuse `eval_dyadicSeparatedTerms_succ` and `expNegTaylor`; do not restate a
  simplified surrogate mask;
- keep the obstruction on the positive-frequency exponential numerator, since
  the fermionic logistic denominator is a positive column scaling and cannot
  repair the minor sign;
- do not infer failure of a different smooth dyadic approximation or of the
  unmasked local Taylor tail.

Independent checks:

- evaluate the four entries and determinant at high precision using a stable
  Taylor-tail expression;
- sweep nearby rational/grid points across the same mask boundary;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if direct `Real.exp_bound` arithmetic becomes unwieldy after three attempts,
  prove a reusable alternating-tail interval lemma first and instantiate it at
  the four rational products.

Non-claims:

- this obstruction does not rule out a smooth partition, an unmasked
  remainder theorem, another compatible low-rank decomposition, Conjecture G1,
  or Simons Problem 4.2.

Delivered:

- `dyadicTaylorError` is the positive-frequency error of the repository's
  actual `dyadicSeparatedTerms (8 * p) p s` construction;
- at `p = 1`, `s = 2`, rows `(1/2, 3/5)`, and columns `(2, 5/2)`, four
  reduction lemmas expose three eighth-order Taylor tails and the masked value
  `exp(-3/2)` directly through `eval_dyadicSeparatedTerms_succ`;
- `eighthOrder_maskTransition_det_pos` proves the resulting determinant is
  positive using rational Taylor-remainder bounds and elementary exponential
  inequalities, with no floating-point premise;
- `dyadicTaylorError_not_strictSignRegularAtOrder_two` proves that this
  positive ordered minor contradicts the negative signature required by
  `expKernelSignature 2`.

Result:

- the Phase AN factorial-free continuation mechanism cannot be applied
  directly to the current hard-masked dyadic approximation remainder;
- the failure occurs at the first transition and at the smallest nontrivial
  minor order, so an all-orders proof search for this particular remainder is
  closed rather than merely deferred;
- the next approximation-level route should remove the discontinuous mask,
  for example by a smooth overlap or a globally sign-compatible remainder,
  and must retest order two before attempting an all-orders theorem.

Verification:

- the focused module and expanded public axiom audit pass with only Lean's
  permitted default axioms;
- a 100-decimal evaluator built from the convergent exponential tail agrees
  with the repository implementation and places the witness determinant in
  `(4.96e-6, 4.97e-6)`;
- all 36 nearby transition-crossing samples in the independent sweep have the
  same positive determinant sign;
- `./scripts/verify.sh` passes with 76 Python tests and no known dependency
  vulnerabilities.

## Phase AP — unmasked Taylor-tail beta structure

State: complete; delivered in merged PR #44 with green CI.

Phase AO rules out the discontinuously masked remainder, not the local Taylor
tail itself. Taylor's integral remainder rewrites the positive eighth-order
tail as a monomial row/column scaling of a beta-weighted Laplace transform.
This is the natural bridge to the basic composition and variation-diminishing
ideas of total-positivity theory, while retaining an exact Lean target.

Frozen public Lean targets:

- define the unmasked eighth-order negative-exponential tail and its normalized
  beta moment;
- prove the exact unit-interval integral representation from mathlib's Taylor
  theorem with integral remainder;
- prove strict positivity for positive arguments;
- factor every ordered `2 x 2` tail minor into positive monomial row/column
  factors and the corresponding normalized beta-moment minor, isolating the
  one remaining sign inequality without assuming it.

Proof contract:

- derive the representation from `map_add_eq_sum_add_integral_iteratedFDeriv`
  or `taylor_integral_remainder`, not from a new axiom or an equality fitted
  numerically;
- keep the beta weight and factorial normalization explicit;
- treat total-positivity composition as motivation only until its hypotheses
  are formalized for this product kernel;
- do not reintroduce hard masks or claim that local order-two evidence closes
  the all-orders, mixed-column, or global dyadic problem.

Independent checks:

- compare the integral, convergent tail, and direct implementation at high
  precision over the approximation regime `0 < x <= 2`;
- sweep ordered positive row/column pairs over logarithmic and boundary-focused
  grids and record the smallest signed order-two margin;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if the exact integral representation does not reduce cleanly after three
  proof attempts, retain the interval `[0,x]` representation and defer the
  unit-interval change of variables explicitly;
- if order-two signs fail numerically, minimize and formalize the witness
  instead of pursuing the beta-moment inequality.

Non-claims:

- this phase does not prove all-orders strict sign regularity of an unmasked
  Taylor tail, construct a global smooth dyadic approximation, prove
  Conjecture G1, or solve Simons Problem 4.2.

Delivered:

- `eighthOrderTail_eq_betaMoment` derives directly from mathlib's Taylor
  theorem with integral remainder the identity
  `tail(x) = x^8 / 7! * integral_0^1 (1-u)^7 exp(-ux) du`;
- `eighthOrderBetaMoment_pos` and `eighthOrderTail_pos` prove the normalized
  moment and the positive-argument tail are strictly positive;
- `eighthOrderTail_fin_two_det_factor` removes all monomial row/column factors
  from an arbitrary order-two tail minor;
- `eighthOrderTail_fin_two_det_neg_iff_betaMoment` proves that on positive row
  and column arguments the tail minor is negative exactly when its normalized
  beta-moment minor is negative.

Result:

- the unmasked tail retains a smooth positive integral structure absent from
  the hard-masked approximation;
- the next exact analytic target is no longer a four-entry transcendental
  determinant: it is the order-two submodularity inequality for the single
  scalar beta moment;
- the representation is compatible with classical total-positivity
  composition ideas, but the product argument `t*omega` means the basic
  composition formula cannot be invoked without an additional kernel theorem.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- at 100 decimal digits, the beta integral agrees with both the convergent
  tail and direct Taylor subtraction at six points from `0.001` through `2`;
- the beta representation agrees with the repository evaluator at every node
  of a positive `6 x 7` grid, and all 315 ordered order-two minors have the
  expected negative sign with normalized margin above `4.3e-6`;
- `./scripts/verify.sh` passes with 77 Python tests and no known dependency
  vulnerabilities.

## Phase AQ — tilted-beta moment gap

State: verified locally; delivery pending.

Phase AP reduces the local tail's order-two sign to the normalized moment
`M(x) = integral_0^1 (1-u)^7 exp(-ux) du`. Differentiating its logarithmic
elasticity produces the tilted-beta variance numerator. On the local Taylor
regime `0 <= x <= 2`, it is enough to prove the sharper elementary gap
`2 J_2(x) < J_1(x)`, where `J_k` inserts `u^k` in the same integral.

Frozen public Lean targets:

- define the first and second tilted-beta moments consistently with Phase AP's
  normalized moment;
- prove the exact unweighted identity
  `integral_0^1 u(1-2u)(1-u)^7 du = 1/120`;
- prove `2 J_2(x) < J_1(x)` for every `x >= 0` by comparing the decreasing
  exponential weight on the two sides of `u = 1/2`;
- derive the strict variance-numerator inequality
  `x (J_2 M - J_1^2) < J_1 M` for `0 <= x <= 2`.

Proof contract:

- use interval-integral monotonicity and exact polynomial calculus; do not
  replace the continuous moment inequality by a sampled certificate;
- keep all denominators and strictness explicit at `x = 0` and `x = 2`;
- expose the variance numerator without yet assuming differentiation under the
  integral or the final monotonic-ratio theorem;
- do not broaden the result beyond the local eighth-order regime in this phase.

Independent checks:

- recompute the polynomial integral exactly with `Fraction` arithmetic;
- verify the moment gap and variance inequality at 100 digits on endpoint,
  logarithmic, and dense local grids;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if direct polynomial integration is brittle after three attempts, introduce
  the explicit antiderivative and prove its derivative by `ring`;
- if the midpoint comparison is insufficient, formalize the single-crossing
  cumulative-integral argument rather than weakening strictness.

Non-claims:

- this phase does not yet prove the beta-moment kernel's order-two sign,
  all-orders tail sign regularity, a global smooth dyadic approximation,
  Conjecture G1, or Simons Problem 4.2.

Delivered:

- `eighthOrder_unweighted_moment_gap` evaluates the exact polynomial integral
  to `1/120` through an explicit degree-ten antiderivative;
- `eighthOrderBetaMoment_two_mul_two_lt_one` uses the one sign change at
  `u = 1/2` and the decreasing exponential tilt to prove
  `2 J_2(x) < J_1(x)` for every `x >= 0`;
- `eighthOrderBetaMoment_variance_numerator_lt` combines that gap with Phase
  AP's strict positivity to prove the local variance-numerator inequality for
  every `0 <= x <= 2`.

Result:

- the numerical order-two pattern from Phase AP now has its decisive local
  moment inequality proved exactly;
- after differentiating `M`, the remaining route is to show the logarithmic
  elasticity `-x J_1(x)/M(x)` is strictly decreasing and transfer that
  monotonicity to ordered products `t*omega`;
- the proof's midpoint comparison is a one-dimensional monotone-likelihood-
  ratio argument, a concrete variation-diminishing mechanism rather than an
  appeal to total-positivity terminology alone.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic independently recovers the unweighted value
  `1/120`;
- a 100-decimal grid containing both endpoints, logarithmic near-zero values,
  and 101 equally spaced local points verifies the moment gap and variance
  bound, with observed margins above `0.0062` and `0.0008` respectively;
- `./scripts/verify.sh` passes with 78 Python tests and no known dependency
  vulnerabilities.
