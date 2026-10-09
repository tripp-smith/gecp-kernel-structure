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
