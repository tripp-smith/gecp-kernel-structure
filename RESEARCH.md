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

State: complete; delivered in merged PR #45 with green CI.

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

## Phase AR — beta-moment logarithmic elasticity

State: complete; delivered in merged PR #46 with green CI.

Phase AQ proves the strict numerator inequality obtained by formally
differentiating the logarithmic elasticity of the normalized beta moment. This
phase supplies that missing differentiation step and converts the pointwise
inequality into a strict monotonicity theorem on the full local interval.

Frozen public Lean targets:

- prove that `eighthOrderBetaMoment` has derivative
  `-eighthOrderBetaMomentOne` at every real argument;
- prove that `eighthOrderBetaMomentOne` has derivative
  `-eighthOrderBetaMomentTwo` at every real argument;
- define
  `eighthOrderBetaElasticity x = -x * eighthOrderBetaMomentOne x /
  eighthOrderBetaMoment x` and prove its exact quotient-rule derivative;
- prove `eighthOrderBetaElasticity_strictAntiOn` on `Set.Icc 0 2` from the
  Phase AQ variance-numerator inequality.

Proof contract:

- justify differentiation under the interval integral with mathlib's
  dominated parametric-integral API; do not introduce the derivative identity
  as an assumption or derive it from numerical evidence;
- use the exact quotient-rule numerator and the previously proved strict
  moment inequality, keeping positivity of the denominator explicit;
- prove continuity at both endpoints and strict decrease on the closed
  interval through the derivative sign on its interior;
- keep the result one-dimensional: the ordered-product transfer belongs to a
  separate phase.

Independent checks:

- compare analytic derivatives against 100-digit differentiated integrals at
  endpoint, logarithmic, and dense local points;
- verify strict monotonicity and retain the smallest observed derivative
  margin as numerical evidence only;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if the specialized dominated-integral proof repeats substantially, factor a
  private generic beta-power-moment derivative lemma and instantiate powers
  zero and one;
- if the compact domination bound becomes the only obstruction after three
  proof attempts, freeze that bound as the next exact helper instead of
  weakening the derivative statement.

Non-claims:

- this phase does not yet prove the beta-moment kernel's order-two sign, the
  unmasked tail's order-two sign, all-orders tail sign regularity, a global
  smooth dyadic approximation, Conjecture G1, or Simons Problem 4.2.

Delivered:

- a private generic beta-power-moment lemma applies dominated differentiation
  under the interval integral for every natural moment order;
- `eighthOrderBetaMoment_hasDerivAt` and
  `eighthOrderBetaMomentOne_hasDerivAt` instantiate the zeroth and first
  moments;
- `eighthOrderBetaElasticity_hasDerivAt` gives the exact quotient-rule
  derivative whose numerator is Phase AQ's strict inequality;
- `eighthOrderBetaElasticity_strictAntiOn` proves strict decrease on the closed
  local interval `Set.Icc 0 2`.

Result:

- the one-dimensional analytic chain from the normalized beta integral to a
  strict elasticity order is now formalized without numerical assumptions;
- the remaining order-two task is algebraic/order-theoretic: transfer strict
  elasticity decrease across the four ordered products and recover the beta-
  moment determinant sign.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- 100-decimal numerical differentiation agrees with both moment identities at
  six points from `0` through `2`;
- a 104-point local sweep verifies strict decrease and keeps the observed
  derivative margin `-E'(x)` above `0.078`;
- `./scripts/verify.sh` passes with 79 Python tests and no known dependency
  vulnerabilities.

## Phase AS — local beta-moment order-two sign

State: complete; delivered in merged PR #47 with green CI.

Phase AR proves that the logarithmic elasticity
`E(x) = -x J_1(x) / M(x)` is strictly decreasing on `0 <= x <= 2`. For
positive ordered scales `t₀ < t₁`, this makes the ratio
`M(t₀ * ω) / M(t₁ * ω)` strictly increasing in positive `ω`, which is exactly
the cross-multiplied order-two determinant inequality.

Frozen public Lean targets:

- define the beta-moment scale ratio and prove its exact derivative from the
  Phase AR moment derivative;
- prove strict monotonicity of the ratio on every positive frequency interval
  whose largest product `t₁ * ω₁` is at most `2`;
- prove the corresponding strict negative order-two beta-moment determinant;
- transfer the sign through Phase AP's positive factorization to the unmasked
  eighth-order Taylor tail.

Proof contract:

- obtain the derivative numerator from the quotient rule and derive its strict
  positivity by cross-multiplying the Phase AR elasticity inequality through
  positive moments;
- keep every positivity and local-product condition explicit; do not use a
  logarithm whose domain obligations obscure the argument;
- use the monotone scale ratio to prove the determinant sign, not an sampled
  four-point certificate;
- keep the theorem local to the eighth-order regime `t₁ * ω₁ <= 2`.

Independent checks:

- verify the scale-ratio derivative against 100-digit numerical
  differentiation on asymmetric row pairs and endpoint/logarithmic frequency
  grids;
- verify all ordered beta-moment and tail minors on an independent product-
  bounded grid and record the smallest signed margins;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if direct derivative positivity is algebraically brittle after three
  attempts, isolate the exact cross-multiplication lemma between elasticity
  order and the scale-ratio numerator;
- if closed-interval strict monotonicity creates only endpoint obligations,
  prove continuity globally and derivative positivity on the open interior
  rather than weakening the determinant result.

Non-claims:

- this phase does not prove all-orders tail sign regularity, sign regularity of
  a global smooth dyadic remainder, a factorial-free approximation theorem,
  Conjecture G1, or Simons Problem 4.2.

Delivered:

- `eighthOrderBetaMomentScaleRatio_hasDerivAt` gives the exact quotient-rule
  derivative of `M(t₀ * ω) / M(t₁ * ω)`;
- `eighthOrderBetaMomentScaleRatio_strictMonoOn` derives strict increase on
  every positive product-bounded frequency interval from Phase AR elasticity;
- `eighthOrderBetaMoment_fin_two_det_neg` proves the normalized local
  order-two minor is strictly negative;
- `eighthOrderTail_fin_two_det_neg` transfers that sign through Phase AP's
  exact positive tail factorization.

Result:

- the entire local order-two sign argument for the unmasked eighth-order tail
  is now formal: beta integral, moment inequality, differentiation,
  elasticity order, scale-ratio order, and determinant sign;
- the next structural question is whether this order-two mechanism extends to
  all compound orders or can support a smooth overlapping approximation whose
  remainder is compatible with Phase AN's factorial-free continuation bound.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- 100-decimal numerical differentiation validates the scale-ratio derivative
  on three asymmetric row pairs, with observed derivative above `0.0015`;
- an independent product-bounded grid verifies all beta-moment and tail minors,
  with signed margins above `2.4e-5` and `1e-46` respectively;
- `./scripts/verify.sh` passes with 80 Python tests and no known dependency
  vulnerabilities.

## Phase AT — truncated-power composition bridge

State: complete; delivered in merged PR #48 with green CI.

The Phase AP beta integral admits the scaling substitution `v = ω u`:
`ω^8 M(tω) = integral_0^ω (ω-v)^7 exp(-tv) dv`. Consequently the unmasked
tail is, up to the positive row factor `t^8 / 7!`, the composition of the
sign-regular exponential kernel with the truncated-power kernel. This is the
classical spline/variation-diminishing structure suggested by the successful
order-three through order-five numerical minors.

Frozen public Lean targets:

- define the nonnegative seventh-degree truncated-power kernel
  `max (ω-v) 0 ^ 7`;
- prove the exact scaled-beta identity on positive frequencies;
- prove a common-cutoff composition identity that replaces the column-varying
  upper limit by any shared cutoff above `ω`;
- derive the corresponding composition representation of `eighthOrderTail`.

Proof contract:

- obtain the scaling identity from mathlib's interval-integral change of
  variables, with exact powers and no asymptotic argument;
- prove the common-cutoff extension by interval additivity and exact vanishing
  of the truncated power above its knot;
- expose positive row scaling separately so future determinant composition can
  reuse the kernel identity directly;
- do not claim total positivity from the representation alone.

Independent checks:

- recompute the beta, scaled truncated-power, and common-cutoff integrals at
  100 digits on asymmetric points;
- sweep ordered minors through orders two to five on a product-bounded grid to
  test the all-orders direction without promoting it to proof;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if a common-cutoff identity is blocked only by max-function integration,
  first prove exact interval vanishing above the knot and then compose with
  interval additivity;
- if strict all-orders sign needs a missing continuous Cauchy--Binet theorem,
  record that theorem as the next contract rather than asserting sign
  regularity from numerical evidence.

Non-claims:

- this phase does not prove total positivity of the truncated-power kernel,
  an integral Cauchy--Binet composition theorem, all-orders tail sign
  regularity, a globally sign-regular dyadic remainder, Conjecture G1, or
  Simons Problem 4.2.

Delivered:

- `eighthOrderTruncatedPowerKernel` defines the nonnegative kernel
  `max (ω-v) 0 ^ 7`;
- `eighthOrderBetaMoment_scaled_eq_truncatedPower` proves the exact scaling
  substitution `ω^8 M(tω) = integral_0^ω (ω-v)^7 exp(-tv) dv`;
- `eighthOrderTruncatedPower_integral_commonCutoff` proves that any shared
  cutoff above the knot gives the same composition integral;
- `eighthOrderTail_eq_truncatedPowerComposition_commonCutoff` exposes the
  unmasked tail as a positive row scaling of that shared-domain composition.

Result:

- the product argument that blocked a direct application of the basic
  composition formula in Phase AP has been removed by the scaled variable
  `v = ωu`;
- the next exact target is now structural and reusable: total nonnegativity of
  the truncated-power kernel plus a finite-order continuous Cauchy--Binet
  formula with strictness conditions.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- the beta integral, scaled truncated-power integral, shared-cutoff integral,
  and tail composition agree at 100 digits on four asymmetric configurations;
- an independent product-bounded grid has the expected strict signature at
  orders two through five, with observed margins above `6e-6`, `4e-12`,
  `1.7e-18`, and `2e-24` respectively;
- `./scripts/verify.sh` passes with 81 Python tests and no known dependency
  vulnerabilities.

## Phase AU — finite composition strictness

State: complete; delivered in merged PR #49 with green CI.

Phase AT reduces the eighth-order unmasked tail to a shared-domain
composition of the exponential kernel and the seventh-degree truncated-power
kernel. Before introducing product-measure integration, the exact finite
algebraic strictness mechanism should be isolated: an ordered Cauchy--Binet
sum is bounded below by any selected summand when all signed summands are
nonnegative, and is therefore strictly signed when one selected summand is
positive.

Frozen public Lean targets:

- define a readable signed ordered Cauchy--Binet summand;
- prove nonnegativity of the signed determinant from nonnegativity of every
  ordered summand;
- prove a selected-summand lower bound and strict signed determinant theorem;
- keep the result generic over rectangular real matrices so finite
  quadratures and later simple-function approximations can reuse it.

Proof contract:

- use the existing ordered expansion `Matrix.det_mul_rect` rather than a new
  determinant expansion;
- distribute the external signature through the finite sum exactly;
- retain an explicit strictly positive summand instead of relying only on
  positivity of a sequence of determinants;
- inspect the assumptions to ensure they constrain individual composition
  minors and do not restate the determinant conclusion.

Independent checks:

- use exact rational arithmetic on a signed `2 x 3` by `3 x 2` composition;
- verify every ordered signed summand, the Cauchy--Binet sum, the selected
  lower bound, and strictness independently;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if the ordered-sum proof needs only a missing finite-sum lemma, prove the
  smallest local helper;
- if the actual truncated-power factor has a negative ordered minor, record a
  minimized exact obstruction rather than weakening the generic algebraic
  theorem;
- otherwise make total nonnegativity of the truncated-power factor the next
  theorem contract.

Non-claims:

- this phase does not prove total nonnegativity of the truncated-power kernel,
  continuous or integral Cauchy--Binet, all-orders tail sign regularity, a
  globally sign-regular dyadic remainder, Conjecture G1, or Simons Problem
  4.2.

Delivered:

- `Matrix.signedOrderedCauchyBinetTerm` names the signed product of the two
  ordered minors at one increasing intermediate tuple;
- `Matrix.signature_mul_det_mul_rect_eq_sum` distributes the orientation
  through the existing ordered Cauchy--Binet expansion;
- `Matrix.signature_mul_det_mul_rect_nonneg` proves signed determinant
  nonnegativity from termwise nonnegativity;
- `Matrix.signedOrderedCauchyBinetTerm_le_signature_mul_det_mul_rect` retains
  any selected summand as an exact lower bound, and
  `Matrix.signature_mul_det_mul_rect_pos` converts one positive witness into
  strict signed determinant positivity.

Result:

- finite positive quadratures of a sign-compatible kernel composition now
  inherit both the desired determinant orientation and an explicit strict
  lower bound;
- the hypotheses are local to the factor minors and do not assume the product
  determinant conclusion;
- the next structural target is total nonnegativity of
  `max (ω-v) 0 ^ 7`, including the support-induced zero cases and a strict
  witness region suitable for the retained-term argument.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact rational arithmetic gives the three signed ordered summands
  `1, 4, 1`, their sum and the signed product determinant `6`, and the retained
  strict lower bound `4 <= 6`;
- `./scripts/verify.sh` passes with 82 Python tests and no known dependency
  vulnerabilities.

## Phase AV — truncated-power total-nonnegativity base

State: complete; delivered in merged PR #50 with green CI.

Phase AU reduces strictness of a finite composition determinant to common-sign
factor minors plus one positive summand. The factor still missing is the
seventh-degree truncated-power kernel `max (ω-v) 0 ^ 7`. Classical Pólya-
frequency theory predicts total nonnegativity at every order for positive
integer powers, with strictness controlled by the Schoenberg--Whitney support
geometry. The Lean research loop begins with the first nontrivial order so
that support zeros and the orientation convention are fixed before an
all-orders proof is attempted.

Frozen public Lean targets:

- add reusable `TotallyNonnegativeAtOrder` and `TotallyNonnegative` kernel
  predicates alongside the existing strict-sign-regular predicates;
- prove total nonnegativity of `eighthOrderTruncatedPowerKernel` at orders one
  and two on arbitrary strictly increasing real tuples;
- prove the exact order-two strict-support implication: the determinant is
  positive when each ordered row knot lies strictly below its paired column
  knot;
- expose the scalar cross-product inequalities used by the matrix theorems.

Proof contract:

- split the support geometry exactly rather than replacing `max` by a
  polynomial globally;
- in the full-support case reduce the seventh powers to the base cross-product
  identity `(v₁-v₀)(ω₁-ω₀) > 0`;
- distinguish nonnegative minors from strict minors in the public API;
- do not import the classical all-orders theorem as an axiom.

Independent checks:

- enumerate exact rational separated, interlacing, boundary-zero, and
  full-support configurations;
- verify the determinant sign and strict-support criterion directly;
- run focused Lean builds, the expanded public axiom audit, and root
  verification.

Decision rule:

- if the two-by-two claim fails, minimize and record the exact support
  configuration;
- if it succeeds, use the support classification to choose between a
  Schoenberg--Whitney determinant induction and a repeated order-kernel
  composition proof for all orders.

Non-claims:

- this phase does not prove all-orders total nonnegativity of the
  truncated-power kernel, continuous Cauchy--Binet, all-orders tail sign
  regularity, a globally sign-regular dyadic remainder, Conjecture G1, or
  Simons Problem 4.2.

Delivered:

- `GECP.TotallyNonnegativeAtOrder` and `GECP.TotallyNonnegative` distinguish
  nonnegative ordered minors from the repository's existing strict-sign-
  regular predicates;
- `eighthOrderTruncatedPowerKernel_cross_nonneg` proves the global ordered
  cross-product inequality with all support-zero cases included;
- `eighthOrderTruncatedPowerKernel_cross_pos_iff` proves that, for strict row
  and column pairs, positivity is equivalent to the two paired support
  inequalities;
- `eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_one` and
  `eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_two` establish
  the first two total-nonnegativity orders;
- `eighthOrderTruncatedPowerKernel_fin_two_det_pos_iff` packages the exact
  order-two Schoenberg--Whitney criterion at the matrix-minor level.

Result:

- support zeros are now part of the formal theorem rather than excluded by a
  full-support assumption;
- in the full-support case the proof reduces seventh powers to the strict base
  identity controlled by `(v₁-v₀)(ω₁-ω₀)`;
- a [primary-source check of Karlin's Pólya-frequency theorem](https://arxiv.org/html/2008.05121v2#S1.SS2)
  confirms total nonnegativity for positive integer truncated powers after
  removing the theorem's harmless positive exponential row and column factors,
  so the most promising all-orders route is a normalized truncated-power
  convolution ladder rather than separate determinant expansions at each
  order;
- the next exact target is the shared-domain recurrence for normalized
  truncated powers and the all-orders total nonnegativity of its order-kernel
  base, followed by continuous composition closure.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact rational arithmetic checks separated, boundary, interlacing, and
  full-support witnesses, including the full-support determinant `14197`;
- all 784 ordered pair configurations on the integer grid `[-3,4]` satisfy
  nonnegativity and the exact paired-support strictness equivalence;
- `./scripts/verify.sh` passes with 83 Python tests and no known dependency
  vulnerabilities.

## Phase AW contract: normalized truncated-power convolution ladder

Phase AW freezes the analytic recurrence suggested by Phase AV's external
total-nonnegativity check.  The purpose is to express every normalized
positive-integer truncated power as one more composition with the elementary
order kernel, on a domain that is shared by all matrix columns.

Frozen public Lean targets:

- define `normalizedTruncatedPowerKernel degree v omega` with exact factorial
  normalization and the same support convention as the existing kernel;
- identify degree zero with the elementary order kernel away from its
  measure-zero boundary;
- prove the successor integral on its natural ordered interval;
- prove the exact shared-domain successor recurrence when the common bounds
  contain both endpoints;
- identify `eighthOrderTruncatedPowerKernel` with `7!` times normalized degree
  seven.

Proof contract:

- derive the recurrence from the exact power integral and the factorial
  successor identity, not from a classical total-positivity theorem;
- handle the order-kernel discontinuity by an explicit almost-everywhere
  argument and keep endpoint conventions visible;
- keep all domain-order assumptions explicit so the result can feed a later
  continuous Cauchy--Binet theorem.

Independent checks:

- verify the normalized recurrence with exact rational arithmetic for several
  degrees, ordered endpoints, and expanded common domains;
- run focused Lean builds, extend the public axiom audit, and run root
  verification.

Decision rule:

- if the shared-domain identity is formalized, the next phase targets total
  nonnegativity of the degree-zero order kernel and continuous composition
  closure;
- if interval-endpoint measurability blocks the shared-domain theorem, retain
  only the ordered-interval recurrence and record the precise missing mathlib
  lemma rather than weakening assumptions silently.

Non-claims:

- this phase does not prove all-orders total nonnegativity, continuous
  Cauchy--Binet, all-orders tail sign regularity, a globally sign-regular
  dyadic remainder, Conjecture G1, or Simons Problem 4.2.

Delivered:

- `normalizedTruncatedPowerKernel` packages every nonnegative integer degree
  with factorial normalization and strict support;
- `normalizedTruncatedPowerKernel_zero` identifies degree zero with the
  elementary strict order kernel;
- `normalizedTruncatedPowerKernel_succ_integral` proves the exact successor
  recurrence on the natural interval;
- `normalizedTruncatedPowerKernel_succ_integral_commonDomain` proves the same
  recurrence on any shared interval containing the row and column endpoints;
- `eighthOrderTruncatedPowerKernel_eq_factorial_mul_normalized` identifies the
  Phase AT kernel with `7!` times normalized degree seven.

Result:

- the degree-seven factor is now connected by an exact seven-step convolution
  ladder to the degree-zero order kernel;
- the shared-domain theorem explicitly discharges the two zero-support outer
  intervals and the measure-zero endpoint convention, so no column-dependent
  integration domain remains;
- the next obstruction is algebraic/measure-theoretic rather than analytic:
  prove all-orders total nonnegativity of the degree-zero order kernel and a
  continuous determinant-composition closure theorem.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic independently expands and integrates
  `(omega-z)^degree` term by term for degrees zero through seven, including
  rational and coincident endpoints and enlarged common domains;
- `./scripts/verify.sh` passes with 84 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #51.

## Phase AX contract: all-orders order-kernel total nonnegativity

Phase AX attacks the algebraic base of the convolution ladder.  On increasing
row and column tuples, the degree-zero matrix is a staircase of zeros and
ones.  Expanding its first column yields a three-way induction: a zero column,
duplicate all-one rows, or one leading pivot followed by the smaller staircase
minor.

Frozen public Lean targets:

- name the strict order kernel underlying normalized degree zero;
- prove that every ordered square minor of this kernel is exactly zero or one;
- derive all-orders `TotallyNonnegative` for both the named order kernel and
  `normalizedTruncatedPowerKernel 0`.

Proof contract:

- induct on the matrix order using the first-column staircase geometry;
- in the two-support-row case prove determinant zero from duplicate rows;
- in the one-support-row case reduce the Laplace expansion exactly to the
  tail minor with positive sign;
- do not appeal to an external total-unimodularity or total-positivity axiom.

Independent checks:

- exhaust strictly increasing integer row and column tuples through order
  five and verify every determinant is zero or one;
- run focused Lean builds, extend the public axiom audit, and run root
  verification.

Decision rule:

- if the exact induction closes, the next phase targets a bounded continuous
  composition theorem capable of propagating this result up the normalized
  truncated-power ladder;
- if the Laplace reduction fails, record the exact missing finite-index lemma
  rather than weakening to sampled orders.

Non-claims:

- this phase does not prove higher-degree total nonnegativity, continuous
  Cauchy--Binet, all-orders tail sign regularity, Conjecture G1, or Simons
  Problem 4.2.

Delivered:

- `strictOrderKernel` names the strict zero/one kernel at normalized degree
  zero;
- `strictOrderKernel_minor_det_eq_zero_or_one` proves by induction that every
  ordered square minor is exactly zero or one;
- `strictOrderKernel_totallyNonnegative` and
  `normalizedTruncatedPowerKernel_zero_totallyNonnegative` establish the
  all-orders base theorem in both named forms.

Result:

- the first-column induction exposes the exact finite mechanism: no supported
  row gives a zero column, at least two supported rows give duplicate all-one
  rows, and exactly one supported row leaves the tail minor with positive
  cofactor sign;
- this is a stronger finite statement than nonnegativity alone and confirms
  that the degree-zero kernel is a staircase/network matrix with no hidden
  negative orientation;
- together with Phase AW, the remaining gap to all normalized degrees is now
  a continuous composition closure theorem with explicit measurability and
  interval-integrability hypotheses.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic exhausts 3,381 pairs of strictly increasing
  integer tuples on `[-3,3]` through order five and finds every determinant in
  `{0,1}`;
- `./scripts/verify.sh` passes with 85 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #52.

## Phase AY contract: determinant/interval-integral commutation

Phase AY isolates the analytic atom needed by continuous Cauchy--Binet.  With
all but one matrix row fixed, the determinant is a linear functional of the
remaining row; finite dimensionality makes that functional continuous, so it
commutes with the Bochner interval integral.

Frozen public Lean targets:

- define the continuous linear functional obtained by updating one matrix row
  and taking the determinant;
- prove `Matrix.det_updateRow_intervalIntegral`, commuting that functional
  with an interval integral of row vectors;
- derive nonnegative transfer from pointwise nonnegative updated determinants;
- derive strict transfer for a continuous row path that is nonnegative on the
  interval and positive somewhere.

Proof contract:

- use `Matrix.det_updateRow_add`, `Matrix.det_updateRow_smul`, finite-
  dimensional continuity, and `ContinuousLinearMap.intervalIntegral_comp_comm`;
- state interval order, integrability, continuity, and pointwise sign
  hypotheses explicitly;
- do not encode the desired determinant sign inside a certificate wrapper.

Independent checks:

- integrate an affine rational row in an exact `3 x 3` example and recompute
  both determinant sides independently;
- run focused Lean builds, extend the public axiom audit, and run root
  verification.

Decision rule:

- if the row-integral bridge closes, the next phase iterates it over ordered
  disjoint row intervals and combines it with determinant-preserving adjacent
  row differences;
- if finite-dimensional continuity is insufficient, record the exact missing
  normed-space instance rather than specializing the theorem to one order.

Non-claims:

- this phase does not yet prove a multi-row continuous Cauchy--Binet theorem,
  higher-degree total nonnegativity, all-orders tail sign regularity,
  Conjecture G1, or Simons Problem 4.2.

Delivered:

- `Matrix.detUpdateRowCLM` packages the determinant after one row update as a
  continuous linear functional;
- `Matrix.det_updateRow_intervalIntegral` proves exact commutation with the
  Bochner interval integral of a row-vector path;
- `Matrix.det_updateRow_intervalIntegral_nonneg` transfers pointwise
  nonnegativity on an ordered interval;
- `Matrix.det_updateRow_intervalIntegral_pos` transfers strictness from a
  continuous nonnegative path that is positive somewhere.

Result:

- determinant multilinearity is now connected to the interval-integral API
  without expanding permutations or constructing a product measure;
- the strict corollary preserves the witness mechanism needed later for
  Schoenberg--Whitney support regions;
- the next exact target is to iterate this one-row theorem over pairwise
  ordered row intervals after determinant-preserving adjacent row
  differences.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic integrates the affine row
  `[z, 5+z, 1]` on `[0,2]` inside a `3 x 3` determinant and independently
  obtains `6` from both sides, with pointwise determinant `4-z > 0`;
- `./scripts/verify.sh` passes with 86 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #53.

## Phase AZ contract: ordered multi-row determinant integration

Phase AZ iterates the one-row analytic transfer without introducing a product
measure.  The determinant is integrated one row at a time, while all remaining
rows are frozen at admissible samples.  This is the finite-dimensional analogue
of iterated expectation and keeps the proof aligned with determinant
multilinearity.

Frozen public Lean targets:

- extend total nonnegativity from strictly increasing tuples to monotone tuples,
  with repeated rows or columns handled by a zero determinant;
- define a partially integrated row matrix indexed by a finite set of rows;
- prove nonnegativity after integrating every row whenever every simultaneous
  pointwise sample has nonnegative determinant;
- specialize the result to a totally nonnegative kernel on ordered disjoint row
  intervals.

Proof contract:

- induct over the finite set of integrated rows and use
  `Matrix.det_updateRow_intervalIntegral_nonneg` at each step;
- retain explicit interval-order, bound, and Bochner-integrability hypotheses;
- derive weakly ordered pointwise minors from total nonnegativity, using the
  existing repeated-row and repeated-column determinant lemmas at shared
  endpoints;
- do not claim a product-measure Andreief identity or higher-degree closure yet.

Independent checks:

- integrate the strict order kernel exactly over three adjacent rational
  intervals, obtaining the unit upper-triangular matrix and determinant `1`;
- collapse one interval and confirm the boundary determinant is exactly zero;
- run focused Lean builds, the public axiom audit, and root verification.

Decision rule:

- if finite row iteration closes, the next phase proves determinant-preserving
  adjacent row differences for common-upper integrals and combines the two
  results into total-nonnegativity propagation along the normalized
  truncated-power ladder;
- if row-update bookkeeping blocks the induction, retain the weak-order lemma
  and record the smallest missing matrix-update identity.

Non-claims:

- this phase does not yet prove that every normalized truncated-power degree is
  totally nonnegative, all-orders tail sign regularity, Conjecture G1, or
  Simons Problem 4.2.

Delivered:

- `TotallyNonnegativeAtOrder.det_minorMatrix_nonneg_of_monotone` extends the
  strict-tuple definition to weakly ordered samples;
- `partialRowIntervalIntegral` records exactly which rows have already been
  integrated;
- `det_rowIntervalIntegral_nonneg` performs finite iterated row integration;
- `det_rowIntervalIntegral_kernel_nonneg` transfers a totally nonnegative minor
  through ordered disjoint row intervals.

Result:

- multi-row determinant integration now follows from the one-row theorem with
  no Fubini or product-measure dependency;
- shared interval endpoints are covered exactly by zero-minor degeneracy;
- the remaining algebraic bridge is the adjacent-row difference reduction from
  common-upper integrals to ordered disjoint intervals.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic integrates the strict order kernel over three
  adjacent intervals to a unit upper-triangular matrix with determinant `1`,
  while a collapsed middle interval gives determinant zero;
- `./scripts/verify.sh` passes with 87 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #54.

## Phase BA contract: all-orders truncated-power total nonnegativity

Phase BA closes the algebraic bridge from common-upper row integrals to the
ordered disjoint intervals handled by Phase AZ. Successive subtraction of the
next row preserves the determinant and converts nested tails into adjacent
intervals. Applied to the normalized truncated-power recurrence, this raises
total nonnegativity by one degree.

Frozen public Lean targets:

- define successive adjacent row differences and prove determinant invariance;
- prove that a matrix of common-upper row integrals has the same determinant as
  its adjacent-interval row-integral matrix;
- transfer total nonnegativity through common-upper row integration;
- prove integrability and the common-upper recurrence for every normalized
  truncated-power degree;
- prove degree-successor total-nonnegativity closure, all-degree closure, and
  the all-orders result for the seventh-degree factor.

Proof contract:

- reverse the row order and reuse mathlib's predecessor-row determinant
  invariance rather than developing a second elimination recursion;
- use interval additivity to identify row differences exactly;
- choose the maximum of the final row and column knots as a shared upper bound;
- handle order zero separately and derive all positive orders from common-upper
  integration;
- scale the normalized degree-seven result by `7!` only through an explicit
  nonnegative-scaling lemma.

Independent checks:

- enumerate exact rational minors for degrees zero through seven and orders one
  through four on a six-point grid;
- require all 7,088 determinants to be nonnegative;
- run focused Lean builds, the public axiom audit, and root verification.

Decision rule:

- if the degree-successor theorem closes, replace the former order-one/two
  factor boundary with the all-orders theorem and next propagate it through the
  exponential/truncated-power composition;
- if common-upper integrability blocks the result, retain the generic
  determinant reduction and record the missing scalar-to-vector integration
  bridge.

Non-claims:

- total nonnegativity alone does not identify exact strict-support regions;
- this phase does not yet prove all-orders tail sign regularity, Conjecture G1,
  or Simons Problem 4.2.

Delivered:

- `Matrix.successiveRowDifference` and
  `Matrix.det_successiveRowDifference` formalize determinant-preserving
  subtraction of each next row;
- `Matrix.det_commonUpperRowIntegral_eq_adjacent` identifies nested-tail and
  adjacent-interval determinants;
- `det_commonUpperRowIntegral_kernel_nonneg` transfers total nonnegativity
  through common-upper row integration;
- `normalizedTruncatedPowerKernel_intervalIntegrable_to_commonUpper` and
  `normalizedTruncatedPowerKernel_succ_integral_to_commonUpper` supply the
  analytic hypotheses and recurrence;
- `normalizedTruncatedPowerKernel_succ_totallyNonnegativeAtOrder` propagates
  one degree, `normalizedTruncatedPowerKernel_totallyNonnegative` closes every
  degree, and `eighthOrderTruncatedPowerKernel_totallyNonnegative` delivers the
  all-orders seventh-degree factor theorem.

Result:

- the truncated-power factor is no longer limited to minors of orders one and
  two: every ordered square minor is nonnegative;
- the proof is constructive along the convolution ladder and isolates reusable
  common-upper integration machinery;
- exact strictness/support geometry and the continuous tail-composition sign
  theorem remain the next boundary.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic checks 7,088 ordered minors across normalized
  degrees zero through seven and orders one through four on a six-point grid;
- `./scripts/verify.sh` passes with 88 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #55.

## Phase BB contract: Andréief integrand symmetrization

Phase BB isolates the algebraic heart of continuous Cauchy--Binet. Before
integration, sum the row-wise composition integrand over every permutation of
the intermediate sample. Alternation of the right determinant supplies the
permutation sign, and the sum becomes exactly the product of the sampled left
and right determinants.

Frozen public Lean targets:

- define the row-wise composition sample integrand;
- prove its sum over all coordinate permutations equals the product of the two
  sampled determinants;
- specialize to the exponential and seventh-degree truncated-power kernels and
  prove the prescribed weak sign on every strictly ordered intermediate sample.

Proof contract:

- use `Matrix.det_permute` for the right factor and `Matrix.det_apply'` on the
  transpose of the left factor;
- consume the all-orders exponential strict-sign theorem and Phase BA
  truncated-power total nonnegativity directly;
- keep integration and strictness outside the theorem statement so the exact
  remaining analytic obligation is visible.

Independent checks:

- evaluate the six permutation terms of an exact rational `3 x 3` example and
  independently multiply the two determinants;
- require both calculations to equal `-2205`;
- run focused Lean builds, the public axiom audit, and root verification.

Decision rule:

- if symmetrization closes, the next phase proves that the finite product
  measure on a shared interval is invariant under coordinate permutations and
  integrates the identity into continuous Cauchy--Binet;
- if determinant orientations disagree, record the exact permutation
  convention rather than altering the kernel signature.

Non-claims:

- this phase does not yet integrate the symmetrized identity, establish strict
  support for the truncated-power determinant, prove all-orders tail sign
  regularity, Conjecture G1, or Simons Problem 4.2.

Delivered:

- `Matrix.compositionSampleIntegrand` records the unsymmetrized row-wise
  composition term;
- `Matrix.sum_compositionSampleIntegrand_comp_perm` proves exact finite
  permutation symmetrization;
- `eighthOrderComposition_symmetrized_nonneg` combines exponential strict sign
  regularity with all-orders truncated-power total nonnegativity on ordered
  intermediate samples.

Result:

- the determinant algebra needed by Andréief is complete at every finite order;
- the remaining bridge is measure-theoretic: permutation invariance of the
  shared cube and finite-sum/integral exchange;
- strict tail signs additionally require a positive-measure region where the
  truncated-power determinant is positive.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- exact `Fraction` arithmetic sums all six permutation terms of a rational
  `3 x 3` example and independently multiplies its two determinants, obtaining
  `-2205` by both routes;
- `./scripts/verify.sh` passes with 89 Python tests and no known dependency
  vulnerabilities.
- The change passed CI and merged as PR #56.

## Phase BF contract: factorial-free Taylor-tail GECP envelope

State: complete; green CI; merged in PR #60.

Phase BF converts the Phase BE ordered-minor theorem into the exact interface
needed by the earlier continuation-determinant machinery. It is intentionally
restricted to the positive unit band and the eighth-order local Taylor split;
arbitrary dyadic scales and order `8p` are a later phase.

Frozen public Lean targets:

- package `eighthOrderTail (t * omega)` on positive-real subtypes as a global
  `StrictSignRegular` kernel;
- transfer strict sign regularity through the positive fermionic denominator;
- expose the fermionic kernel as an exact eight-separated-term Taylor
  component plus that compatible tail;
- prove the unit-band tail and approximation error bound `1 / 8!`;
- combine strict sign regularity, complete pivoting, and the uniform tail bound
  into the factorial-free selected-core estimate `(1 / 8!)^rank`.

Independent check:

- use 120-digit arithmetic on an eight-by-eight positive unit grid, verify the
  exact decomposition numerically, execute deterministic complete pivoting on
  the fermionic tail, and check every pivot and cumulative selected-core
  determinant against the proved envelope.

Decision rule:

- if subtype packaging or column scaling fails, isolate the missing transport
  lemma without weakening the Phase BE strict-minor theorem;
- if the complete-run composition fails, retain the exact decomposition and
  state the precise mismatch with `ContinuationDeterminant`;
- after success, generalize the same interface from order eight to order
  `8p`, then attack the mixed low-rank-plus-tail determinant estimate.

Non-claims:

- the phase does not yet give a factorial-free bound for determinants of the
  sum of the Taylor component and its tail;
- it does not remove dyadic masks, cover negative frequencies, prove a
  cutoff-uniform block contraction, Conjecture G1, or Simons Problem 4.2.

Delivered:

- `positiveEighthOrderTailKernel_strictSignRegular` packages the Phase BE
  theorem on positive-real row and column subtypes;
- `positiveEighthOrderFermionicTailKernel_strictSignRegular` transfers it
  through the positive denominator scaling;
- `positiveEighthOrderFermionicTaylorApprox_eq_evalSeparated` and
  `fermionicKernel_eq_positiveEighthOrderTaylorApprox_add_tail` give the exact
  eight-term-plus-compatible-tail split;
- the unit-band error theorems prove the explicit `1 / 8!` envelope;
- the two `_completeRun_selectedCore_det_le` theorems compose the envelope
  with complete-pivot sign coherence to obtain `(1 / 8!)^rank` without a
  determinant factorial.

Result:

- the earlier Phase AN continuation theorem now applies to a concrete
  fermionic approximation remainder rather than an abstract compatible tail;
- the remaining determinant problem is specifically mixed: control a sampled
  determinant of the sum from the rank-eight part and the compatible tail;
- the next analytic generalization is from the local eighth-order split to
  the `8p` orders used by the dyadic approximation.

Verification:

- focused Lean builds and the expanded public axiom audit pass with only
  Lean's permitted default axioms;
- a deterministic 120-digit eight-by-eight complete-pivot computation checks
  the exact split, every tail pivot, and every cumulative selected-core bound;
- `./scripts/verify.sh` passes with 93 Python tests and no known dependency
  vulnerabilities.

## Phase BG contract: positive-series compatible fermionic split

State: complete; green CI; merged in PR #61.

The direct negative-exponential split from Phase BF is compatible at the pure
tail endpoint but not coefficientwise compatible with its low-rank Taylor
part. Exact high-precision diagnostics show that the one-tail-column
coefficient of a `2 x 2` determinant has the opposite sign. Phase BG changes
to the equivalent positive coordinate `u = 1 - t`, where
`fermionicKernel t omega = exp (u * omega) / (1 + exp omega)` for positive
frequency and every Taylor feature has positive coefficient.

Frozen public Lean targets:

- define the positive eighth-order exponential tail and its rank-eight Taylor
  component;
- prove the exact positive-series fermionic identity and decomposition;
- prove strict total positivity of the positive tail on positive arguments,
  using the Phase BE composition theorem with the necessary row-orientation
  transport;
- prove the unit-square `1 / 8!` error envelope after positive denominator
  scaling.

Independent checks:

- record the sign failure of the negative-series `2 x 2` mixed coefficient;
- verify on deterministic high-precision grids that every coefficient of
  `det(P + lambda E)` for the positive-series split has the predicted common
  sign through the tested orders.

Decision rule:

- if positive-tail strictness fails under row reversal, isolate the exact
  orientation factor rather than weakening the determinant statement;
- if positive-series mixed coefficients still disagree, minimize and retain a
  counterexample instead of claiming compatibility;
- after success, formalize a coefficientwise mixed-determinant bound and then
  generalize from order eight to `8p`.

Non-claims:

- this phase does not itself bound every mixed determinant or prove a global
  GECP contraction;
- the endpoint `u = 0`, negative frequencies, dyadic assembly, Conjecture G1,
  and Simons Problem 4.2 remain outside this contract.

Failed hypothesis retained:

- for rows and columns `(1/2, 1)`, the one-tail-column coefficient in the
  `2 x 2` negative-exponential split is positive even though the full ordered
  determinant orientation is negative. Thus Phase BF's direct split is not
  coefficientwise sign compatible, despite both endpoint kernels having the
  prescribed pure-minor sign.

Delivered:

- `positiveTaylorSeparatedTerms` gives the positive monomial feature family,
  and `dividePositiveSeparatedRight` absorbs the reflected fermionic
  denominator;
- `positiveExpTaylor_error_le_next_mul_exp` supplies the positive-exponential
  Lagrange envelope;
- `eighthOrderTail_minor_signed_pos_of_rows_ne_zero` exposes the even-power
  strengthening needed for row reversal;
- `positiveSeriesEighthOrderTailKernel_strictSignRegular` proves strict total
  positivity of the positive tail;
- `fermionicKernel_one_sub_eq_positiveSeriesTaylorApprox_add_tail` gives the
  exact reflected rank-eight-plus-tail decomposition;
- the positive-series tail error and complete-run theorems retain the sharp
  `1 / 8!` and `(1 / 8!)^rank` bounds after denominator scaling.

Result:

- the low and high feature families now use the same positive monomial
  orientation, removing the exact order-two obstruction in the original
  coordinate;
- a 140-digit sweep finds every coefficient of `det(P + lambda E)` positive
  through order six for deterministic positive nodes;
- coefficientwise positivity is numerical evidence here, not yet a Lean
  theorem. Its finite Cauchy--Binet formulation is the next contract.

Verification:

- focused Lean builds and the public axiom audit pass with only Lean's
  permitted default axioms;
- the independent mixed-coefficient regression checks both the negative-route
  obstruction and the positive-series orientation through order six;
- `./scripts/verify.sh` passes with 94 Python tests and no known dependency
  vulnerabilities.

## Phase BH contract: block-feature coefficientwise Cauchy--Binet

State: complete; merged as PR #62 with green CI.

Phase BH freezes the finite algebra suggested by Phase BG. A rectangular
feature product is split into low and high feature blocks; scaling every high
feature by `lambda` makes each ordered Cauchy--Binet summand acquire exactly
`lambda^m`, where `m` is the number of selected high features.

Frozen public Lean targets:

- define the high-feature count of an ordered feature selection;
- prove the exact weighted rectangular Cauchy--Binet expansion with exponent
  equal to that count;
- define the coefficient obtained by grouping selections of a fixed count and
  prove it nonnegative whenever every unweighted ordered summand is
  nonnegative;
- specialize the hypotheses to finite positive exponential Taylor features
  on nonnegative strictly ordered nodes.

Independent check:

- use exact rational arithmetic to expand a finite positive Taylor feature
  determinant as a polynomial in the high-block weight and independently
  group Cauchy--Binet terms by high-feature count.

Decision rule:

- if the general sum regrouping is disproportionately costly, retain the
  exact selection-indexed monomial expansion plus fixed-count coefficient
  nonnegativity; that is sufficient for the later limit theorem;
- if a finite positive-Taylor summand changes sign, minimize the exact
  rational witness instead of weakening the generic algebra.

Non-claims:

- this phase does not pass to the infinite positive-exponential tail;
- it does not yet provide a closed determinant magnitude bound, dyadic
  assembly, Conjecture G1, or Simons Problem 4.2.

Delivered:

- `Matrix.highFeatureCount` counts marked features in an ordered rectangular
  selection;
- `Matrix.signature_mul_det_weightedFeatureProduct_eq_sum` gives the exact
  determinant expansion after scaling the marked feature block by `lambda`;
- `Matrix.blockFeatureCoefficient_nonneg` groups any fixed high-feature count
  and proves its coefficient nonnegative from termwise nonnegativity;
- `expTaylor_signedOrderedCauchyBinetTerm_nonneg` and
  `expTaylor_blockFeatureCoefficient_nonneg` discharge those hypotheses for
  positive exponential Taylor features on nonnegative ordered nodes.

Result:

The positive-series observation from Phase BG is now a theorem at every finite
Taylor cutoff and every block degree: its mixed determinant coefficients are
nonnegative sums of ordered Cauchy--Binet terms. This is the finite algebraic
input for a limit argument; it does not yet identify the coefficients of the
exact infinite tail or bound their sum quantitatively.

Verification:

- every new public declaration reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- an exact rational `3 x 3`, six-feature regression independently agrees
  between determinant multilinearity and feature-selection grouping and finds
  all four coefficients strictly positive;
- `./scripts/verify.sh` passes with 95 Python tests and no known dependency
  vulnerabilities.

## Phase BJ contract: infinite-tail mixed coefficients

State: complete; merged as PR #64 with green CI.

Phase BJ passes the finite physical-column coefficient sign from Phase BI to
the exact positive exponential tail. Because a mixed coefficient is a finite
sum of finite determinants, entrywise matrix convergence is sufficient; no
infinite-dimensional Cauchy--Binet theorem is needed.

Frozen public Lean targets:

- prove continuity of `mixedColumnCoefficient` in its high and low matrix
  arguments;
- define the fixed eighth-order Taylor matrix, finite tail truncations, and
  exact exponential tail matrix;
- prove finite tail truncations converge entrywise to the exact tail matrix;
- identify the finite truncations with the Phase BI low/high feature products;
- prove every fixed-cardinality exact-tail mixed-column coefficient
  nonnegative on nonnegative strictly ordered nodes.

Independent check:

- use high precision to show the finite coefficient vector converges to the
  direct exact-tail coefficient vector on a deterministic ordered sample.

Decision rule:

- factor the proof through generic finite-matrix continuity before specializing
  exponential Taylor sums;
- if the feature-product identification is the only obstruction, expose it as
  a separate exact finite-sum lemma rather than weakening the limit theorem.

Non-claims:

- coefficientwise nonnegativity does not by itself give a closed determinant
  magnitude bound;
- dyadic assembly, Conjecture G1, and Simons Problem 4.2 remain open.

Delivered:

- `Matrix.mixedColumnCoefficient_tendsto` proves joint continuity of every
  fixed-cardinality mixed-column coefficient;
- `positiveEighthOrderFiniteTailMatrix_tendsto` passes the finite Taylor tail
  matrices to the exact positive exponential tail;
- `lowFeatureProduct_eight` and `highFeatureProduct_eight` identify those
  matrices with the Phase BI feature split exactly;
- `positiveEighthOrderTailMatrix_apply` identifies the limit entrywise with
  `positiveEighthOrderTail`;
- `positiveEighthOrderTail_mixedColumnCoefficient_nonneg` proves every exact
  mixed coefficient sum nonnegative on nonnegative strictly ordered nodes.

Result:

The coefficientwise sign pattern observed numerically in Phase BG is now
proved for the exact infinite positive exponential tail at every matrix order
and tail-column count. The proof uses only finite determinant continuity after
the polynomial bridge, avoiding an infinite feature-index formalization.

Verification:

- all new public declarations report only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- a 120-digit ordered `3 x 3` regression checks strict positivity of the exact
  coefficient vector and monotone convergence from cutoffs 16, 24, 40, and 64;
- `./scripts/verify.sh` passes with 96 Python tests and no known dependency
  vulnerabilities.

## Phase BK contract: low-rank tail amplification

State: complete; merged as PR #65 with green CI.

Phase BK converts exact coefficientwise nonnegativity into the first
quantitative determinant inequality. If the low block factors through `k`
features, every mixed coefficient below degree `n-k` vanishes. For a tail
weight `lambda >= 1`, every remaining monomial gains at least
`lambda^(n-k)`.

Frozen public Lean targets:

- prove `mixedColumnCoefficient` vanishes when the chosen tail-column count
  plus the low feature rank is smaller than the matrix order;
- express determinant-pencil evaluation through its coefficient support;
- prove the generic amplification inequality
  `lambda^(n-k) * det (E+B) <= det (lambda • E+B)` under low-rank
  factorization and coefficient nonnegativity;
- specialize to the exact eighth-order positive exponential Taylor/tail split
  on nonnegative strictly ordered nodes.

Independent check:

- use high precision to compare both sides across orders 9--16 and several
  `lambda >= 1`, recording the smallest observed slack.

Decision rule:

- keep the comparison division-free so it remains valid without proving the
  amplified determinant positive separately;
- if polynomial support bookkeeping dominates, isolate it as a generic lemma
  instead of weakening the rank threshold.

Non-claims:

- this phase does not yet bound the amplified determinant by a fixed-base
  exponential;
- the GECP contraction, dyadic assembly, Conjecture G1, and Problem 4.2 remain
  open until that envelope is supplied.

Delivered:

- `Matrix.mixedColumnCoefficient_eq_zero_of_rank` proves exact vanishing below
  the codimension of a factored low block;
- `Matrix.mixedColumnPolynomial_rank_amplification` proves the division-free
  generic inequality
  `lambda^(n-k) * det(E + LR) <= det(lambda E + LR)` for `lambda >= 1` and
  nonnegative mixed coefficients;
- `positiveEighthOrderTail_rank_amplification` specializes the comparison to
  the exact positive exponential Taylor/tail split;
- `positiveSeriesFermionic_rank_amplification` preserves it under the positive
  fermionic denominator scaling.

Result:

Every determinant term that survives the rank-eight low block contains at
least `n - 8` exact tail columns. Increasing the tail weight therefore
amplifies the full determinant by at least `lambda^(n-8)`, without division,
factorials, or a finite Taylor cutoff. This isolates the remaining quantitative
problem: bound the amplified totally-positive determinant by `C^n` for a
fixed `C`.

Verification:

- each public declaration reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- a 350-digit sweep checks strict inequality at orders 9--16 for weights 2,
  8, and `8!`;
- `./scripts/verify.sh` passes with 97 Python tests and no known dependency
  vulnerabilities.

## Phase BL contract: amplified-kernel envelope

State: complete; merged as PR #66 with green CI.

Phase BK reduces the determinant-decay problem to an upper envelope for the
tail-amplified positive-series fermionic kernel. Phase BL freezes the first
quantitative half of that envelope. It treats the amplification parameter as
a genuine kernel pencil, proves that increasing a nonnegative mixed-column
pencil preserves strict ordered-minor positivity, and fixes the natural
weight `8!`, where the existing unit-band tail estimate balances exactly.

Frozen public Lean targets:

- prove monotonicity of a determinant pencil between nonnegative scalar
  weights whenever every mixed-column coefficient is nonnegative;
- define the positive-series amplified fermionic kernel and identify each
  ordered sample matrix with `positiveSeriesAmplifiedFermionicMatrix`;
- prove strict total positivity of that kernel for every amplification weight
  at least one;
- prove the pointwise unit-band envelope `abs K <= 2` at weight `8!`;
- derive the existing factorial-free `2^n` selected-core determinant bound
  for complete runs of the amplified kernel.

Independent check:

- evaluate the `8!`-amplified kernel at high precision on endpoint-resolved
  unit-band grids and check ordered determinants through several orders,
  recording the largest entry and smallest positive determinant.

Decision rule:

- reuse mixed-column coefficient positivity rather than attempting a second
  total-positivity representation of the weighted tail;
- keep the envelope at `2`, avoiding numerical optimization of an irrelevant
  constant;
- if complete-run machinery cannot apply directly, retain the strict
  total-positivity and pointwise theorem and record the exact type mismatch,
  without weakening either analytic claim.

Non-claims:

- a complete run of the amplified kernel is not assumed to share pivot nodes
  with a GECP run of the original kernel;
- this phase does not yet prove the arbitrary ordered-sample inequality
  `det K <= 2^n`, compose determinant decay into a GECP residual contraction,
  or close Conjecture G1 or Simons Problem 4.2.

Delivered:

- `Matrix.mixedColumnPolynomial_det_mono` proves monotonicity of every
  nonnegative-coefficient determinant pencil on nonnegative weights;
- `positiveSeriesAmplifiedFermionicKernel_minorMatrix` identifies the kernel
  and matrix formulations exactly;
- `positiveSeriesAmplifiedFermionicKernel_strictSignRegular` proves strict
  total positivity for every weight at least one;
- `positiveSeriesAmplifiedFermionicKernel_factorial_abs_le_two` proves the
  fixed entry envelope at weight `8!`;
- `positiveSeriesAmplifiedFermionicKernel_factorial_completeRun_det_le`
  supplies the corresponding `2^n` envelope for complete runs of that kernel.

Result:

The amplification comparison from Phase BK no longer points at an uncontrolled
kernel: at the balancing weight `8!`, the target kernel is strictly totally
positive and uniformly bounded by two on the unit band. Its own complete GECP
runs therefore have factorial-free selected-core determinants. The remaining
bridge is now sharply isolated: prove the same `2^n` determinant envelope for
an arbitrary ordered sample minor, because the original GECP nodes need not be
complete pivots for the amplified kernel.

Verification:

- every new public declaration reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- a 250-digit endpoint-resolved regression verifies the entry envelope and
  positive ordered determinants through order eight; its largest sampled entry
  is approximately `1.0331592962`, and the order-eight determinant is
  approximately `2.2972915790e-38`;
- `./scripts/verify.sh` passes with 98 Python tests and no known dependency
  vulnerabilities.

## Phase BM contract: total-positive Hadamard envelope

State: complete; merged as PR #67 with green CI.

Phase BL proves that the `8!`-amplified positive-series kernel is strictly
totally positive and bounded by two, but its complete-run envelope does not
apply to arbitrary nodes selected by a different kernel. Phase BM removes
that path-dependence. The proposed finite theorem is the total-positive
analogue of Hadamard's inequality: a strictly totally positive matrix has
determinant at most the product of its diagonal entries.

Frozen public Lean targets:

- define finite strict total positivity in terms of every ordered square
  subminor;
- prove that the northwest scalar Schur complement of a strictly totally
  positive matrix is again strictly totally positive;
- prove by induction that a strictly totally positive matrix satisfies
  `det A <= product i, A i i`;
- transfer this to arbitrary ordered kernel samples with positive minor
  signature and a pointwise bound `B`, obtaining `det <= B^n`;
- specialize to the `8!`-amplified positive-series fermionic kernel on the
  unit band, obtaining the arbitrary-sample bound `det <= 2^n`;
- compose that bound with Phase BK's amplification inequality to obtain
  `(8!)^(n-8) det K_n <= 2^n` for every order at least eight.

Independent check:

- use exact rational arithmetic on finite positive-exponential feature
  matrices through order eight, checking every ordered subminor, the diagonal-
  product inequality, and each northwest Schur-complement stage.

Decision rule:

- prove the generic finite theorem before the fermionic specialization;
- use the existing one-step augmented-determinant factorization for the Schur
  identity instead of introducing an inverse or an LU axiom;
- after three failures around one indexing fact, isolate that fact as an
  explicit helper rather than weakening the determinant claim.

Non-claims:

- this phase does not assert that arbitrary entrywise-bounded matrices have a
  fixed-base determinant envelope;
- it does not identify pivot paths between the original and amplified kernels;
- the final GECP residual contraction, dyadic assembly, Conjecture G1, and
  Simons Problem 4.2 remain open until the determinant decay is composed with
  the existing geometric-mean machinery.

Delivered:

- `Matrix.StrictlyTotallyPositive` packages positivity of every ordered square
  subminor, and `Matrix.StrictlyTotallyPositive.northwestSchur` proves that a
  scalar northwest Schur step preserves it using the existing augmented-core
  determinant factorization;
- `Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive` iterates those
  Schur steps to prove the finite diagonal-product inequality;
- `minorMatrix_det_le_pow_of_strictSignRegular_one` transfers the result to
  every ordered sample of a positive-minor kernel, independently of how the
  sample nodes were selected;
- `positiveSeriesAmplifiedFermionicKernel_factorial_minor_det_le_two_pow`
  gives the arbitrary unit-band amplified-minor envelope `det <= 2^n`;
- `positiveSeriesFermionic_factorial_amplified_det_le_two_pow` and
  `positiveSeriesFermionic_det_le_two_pow_div_factorial_pow` compose the
  envelope with Phase BK, giving
  `(8!)^(n-8) det K_n <= 2^n` and its quotient form for every `n >= 8`.

Verification:

- every new public theorem reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the public axiom audit;
- an exact `Fraction` regression checks all 12,869 nonempty ordered square
  minors of an `8 x 8` positive-feature matrix, plus the determinant/product
  inequality and every northwest Schur stage;
- `./scripts/verify.sh` passes with 99 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #67.

## Phase BN contract: positive-series GECP composition

State: complete; merged as PR #68 with green CI.

Phase BM supplies exponential determinant decay for ordered samples of the
reflected positive-series kernel, while Phase T bounds a complete-pivot
residual power by the selected-core determinant. The selected pivots occur in
algorithmic rather than sorted order, so the exact remaining bridge is an
absolute determinant bound invariant under independent row and column
permutations.

Frozen public Lean targets:

- prove a generic sorting lemma that transfers an ordered-sample determinant
  upper bound to the absolute determinant of every injective sample of a
  positive-minor kernel;
- define the unamplified reflected positive-series fermionic kernel and identify
  its sample matrices with `positiveSeriesFermionicMatrix`;
- transfer Phase BM to arbitrary injective unit-band samples, obtaining
  `|det K_n| <= 2^n / (8!)^(n-8)` for `n >= 8`;
- compose that estimate with
  `strictSignRegular_gecp_error_pow_le_selectedCore_det` for exact complete-
  pivot runs on the positive unit band;
- discharge the explicit arithmetic at `n = 10` and prove that every residual
  entry after ten pivots has magnitude at most `1/2`.

Independent check:

- run high-precision complete pivoting on deterministic positive-unit-band
  grids, verify the selected-core determinant inequality at rank ten, and
  confirm the residual maximum is below the formal half bound.

Decision rule:

- sort the two selected tuples independently and track determinant signs only
  through absolute values; do not assume the GECP pivot order is monotone;
- reuse `Run.selectedRow_injective`, `Run.selectedColumn_injective`, and the
  Phase T residual-power theorem instead of restating pivot-product algebra;
- if endpoint closure is needed, isolate it as a later continuity theorem
  rather than silently replacing strict positive arguments by nonnegative
  ones.

Non-claims:

- this phase treats complete pivoting on the reflected positive-series kernel
  over positive unit-band arguments; it does not yet combine positive and
  negative frequency columns from a full fermionic run;
- it does not yet prove endpoint closure, dyadic rescaling, Conjecture G1, or
  Simons Problem 4.2.

Delivered:

- `abs_minorMatrix_det_le_pow_of_strictSignRegular_one` makes the Phase BM
  diagonal-product envelope invariant under independent permutations of every
  injective row and column sample;
- `positiveSeriesFermionicKernel` names the reflected unamplified kernel,
  `positiveSeriesFermionicKernel_eq_amplified_one` identifies it with unit
  tail weight, and its sample-matrix and strict-positivity theorems expose the
  exact finite representation;
- `positiveSeriesFermionicKernel_abs_minor_det_le` transports the explicit
  `2^n / (8!)^(n-8)` determinant decay from sorted samples to the algorithmic
  order of any injective sample;
- `positiveSeriesFermionicKernel_gecp_error_pow_le` composes that bound with
  the exact complete-pivot geometric-mean inequality;
- `positiveSeriesFermionicKernel_gecp_error_le_half_ten` discharges the
  factorial arithmetic at rank ten and proves a uniform `1/2` residual bound
  throughout the positive unit band.

Verification:

- every new public result reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the public axiom audit;
- a deterministic 180-digit complete-pivot regression on `16 x 16`, `24 x
  24`, and `32 x 32` positive-unit-band grids checks determinant magnitude
  invariance under sorting, the rank-ten determinant envelope, and the final
  residual maximum;
- `./scripts/verify.sh` passes with 100 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #68.

## Phase BO contract: signed unit-band amplification

State: complete; merged as PR #69 with green CI.

The full signed unit band admits a single positive-coordinate representation,
without a sign partition.  For `0 <= t <= 1` and `-1 <= omega <= 1`, setting
`z = 1 - omega` gives
`K(t, omega) = exp(-t) * exp(t*z) / (1 + exp(z-1))`, with
`0 <= z <= 2`.  This affine shift turns the remaining mixed-sign problem into
the same positive exponential tail-amplification problem as Phases BK--BN,
but on a product interval of width two.

Frozen public Lean targets:

- define the shifted signed-band kernel and prove its exact identity with the
  original fermionic kernel;
- prove an order-eight, weight-`128` tail amplification inequality for every
  ordered signed-unit-band sample;
- prove the amplified finite matrix is strictly totally positive and uniformly
  bounded by `2`, then obtain the arbitrary-sample determinant estimate
  `|det K_n| <= 2^n / 128^(n-8)` for `n >= 8`;
- compose the determinant estimate with the existing exact complete-pivot
  residual-power theorem;
- discharge the arithmetic at `n = 12`, proving that twelve exact complete
  pivots contract every residual entry on `[0,1] x [-1,1]` by one half.

Independent check:

- run deterministic high-precision complete pivoting on full signed-unit-band
  grids, check the selected-core determinant envelope at rank twelve, and
  confirm the residual maximum is below the formal half bound.

Decision rule:

- use the affine shift `z = 1 - omega`, not a mixed-column Laplace expansion;
- keep the explicit amplification weight `128`, whose eighth-order tail cost
  is absorbed by `2^8 / 8! < 1/128`;
- sort arbitrary selected tuples only when passing from ordered determinants to
  absolute determinants; do not assume GECP chooses monotone coordinates.

Non-claims:

- this phase treats only the full signed unit band;
- it does not yet prove the shrinking-band recurrence, dyadic rescaling,
  Conjecture G1, or Simons Problem 4.2.

Delivered:

- `fermionicKernel_eq_signedUnitShiftedKernel` replaces the proposed
  positive/negative column partition by the single affine coordinate
  `z = 1 - omega`, valid at both signed-band endpoints;
- `signedUnitShiftedAmplifiedMatrix_det_pos` proves strict total positivity of
  every closed nonnegative shifted sample after arbitrary tail amplification
  `lambda >= 1`;
- `signedUnitShiftedAmplifiedMatrix_one_two_abs_le_two` uses the weighted
  Taylor remainder and the exact choice `lambda = 128` to keep the amplified
  kernel inside the entry envelope `2` on `[0,1] x [0,2]`;
- `signedUnitShiftedMatrix_abs_det_le_two_pow_div_128_pow` and
  `fermionicKernel_signedUnit_abs_sample_det_le` give
  `|det K_n| <= 2^n / 128^(n-8)` for every injective sample on the full closed
  signed unit band, independently of sample order;
- `fermionicKernel_signedUnit_gecp_error_pow_le` composes that bound with the
  exact complete-pivot geometric-mean inequality, and
  `fermionicKernel_signedUnit_gecp_error_le_half_twelve` proves the uniform
  twelve-pivot half contraction.

Verification:

- every new public theorem reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the public axiom audit;
- a deterministic 180-digit complete-pivot regression on full closed
  `16 x 16`, `24 x 24`, and `32 x 32` signed-unit-band grids checks determinant
  magnitude invariance, the rank-twelve envelope, and the final residual;
- `./scripts/verify.sh` passes with 101 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #69.

## Phase BP contract: translated unit-band amplification

State: complete; merged as PR #70 with green CI.

Phase BO centers its affine positive coordinate at frequency zero, but none of
the tail-amplification constants should depend on that center.  For an
arbitrary `center` and `omega` in `[center - 1, center + 1]`, setting
`z = center + 1 - omega` again gives `0 <= z <= 2` and expresses the original
fermionic kernel as a positive row scaling, a positive exponential product,
and a positive column scaling.  The key quantitative observation is uniform:
the weighted exponential term is the original fermionic kernel itself, hence
is at most one on the full time strip independently of `center`.

Frozen public Lean targets:

- define the center-parametric shifted kernel and prove its exact identity
  with the original fermionic kernel;
- prove strict positivity of every ordered tail-amplified sample for arbitrary
  center and amplification weight at least one;
- prove that weight `128` stays inside the entry envelope `2` uniformly in the
  center on `0 <= t <= 1` and `0 <= z <= 2`;
- derive the arbitrary-sample determinant estimate
  `|det K_n| <= 2^n / 128^(n-8)` on every translated width-two band;
- compose the estimate with the exact complete-pivot residual-power theorem
  and prove the same twelve-pivot half contraction on
  `[0,1] x [center-1,center+1]`.

Independent check:

- run deterministic high-precision complete pivoting at several widely
  separated positive and negative centers, checking the amplified entry
  envelope, determinant estimate, and twelve-step residual contraction.

Decision rule:

- retain the center in row and column diagonal factors and reuse the existing
  positive-exponential mixed-coefficient theorem; do not reduce to Phase BO by
  a diagonal comparison whose products would obscure the uniform constants;
- use `fermionicKernel_le_one` for the center-independent tail envelope;
- if direct duplicated matrix plumbing dominates, factor only the shared
  diagonal determinant identity rather than refactoring the verified Phase BO
  API during this phase.

Non-claims:

- a translated original-kernel block is not yet a restart theorem for a Schur
  residual conditioned by earlier pivots;
- this phase does not prove the exceptional-column coefficient sign law, the
  shrinking-band recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered identifiers:

- `translatedUnitShiftedKernel` and
  `fermionicKernel_eq_translatedUnitShiftedKernel` give the exact affine
  coordinate at an arbitrary band center;
- `translatedUnitShiftedAmplifiedMatrix_det_pos` and
  `translatedUnitShiftedAmplifiedMatrix_strictlyTotallyPositive` preserve
  strict total positivity at every tail weight at least one;
- `translatedUnitShiftedAmplifiedMatrix_abs_le_two` proves the center-uniform
  weight-`128` entry envelope;
- `fermionicKernel_translatedUnit_abs_sample_det_le` proves
  `|det K_n| <= 2^n / 128^(n-8)` for every injective sample in any translated
  width-two band;
- `fermionicKernel_translatedUnit_gecp_error_pow_le` and
  `fermionicKernel_translatedUnit_gecp_error_le_half_twelve` compose the
  determinant decay with exact complete pivoting.

Result:

- absolute frequency location has been removed from the local-block theorem:
  the determinant and twelve-pivot contraction constants are identical on
  `[center-1,center+1]` for every real center;
- the remaining transition problem is therefore purely about conditioning by
  the prior selected core, not about recentering the fermionic denominator;
- the next target is a two-flank exceptional-column or equivalent Schur-
  residual amplification theorem that transports this local block through an
  existing outer-scale prefix.

Verification:

- every new public theorem reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the public axiom audit;
- a deterministic 180-digit regression checks the amplified entry envelope,
  selected determinant bound, and twelve-step residual at centers `-8`, `-2`,
  `3`, and `9`;
- `./scripts/verify.sh` passes with 102 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #70.

## Phase CW contract: aggregate mean-shift Lorentzian bridge

State: complete; merged as PR #103 with green CI.

Phase CV refutes every route that requires the left and right augmented
coordinates to have a common sign one pair at a time.  The fresh-order-two
target is only the aggregate quadratic discriminant, however.  In the affine
rank-two chart, the two within-block determinants are weighted covariance
numerators and the mixed coefficient is a cross increment.  The exact law of
total covariance already splits that cross increment into the two rescaled
within-block covariances plus one product of block-mean shifts.  Thus
pointwise cross-block ordering is stronger than the final inequality needs.

Frozen public Lean targets:

- prove the sharp reverse-covariance inequality from positive block masses,
  nonnegative within-block covariance numerators, and nonnegativity of the
  single product of the two aggregate block-mean shifts;
- lift that criterion to affine rank-two coefficients and then to the
  `Fin 2` determinant-pencil discriminant, without division or an inverse;
- keep the existing pointwise-order theorem as a corollary-compatible route,
  while exposing the new aggregate premise as the preferred interface after
  Phase CV.

Independent check:

- condition monomial feature factors through the two exact complete pivots of
  the Phase CV rational exponential kernel and check every cutoff for feature
  truncations four through ten;
- rerun a broader deterministic exact-arithmetic corpus with deeper prefixes,
  retaining coefficient signs, discriminants, and the smallest normalized
  discriminant margin.

Decision rule:

- if the Phase CV aggregate pencil has a negative discriminant, minimize and
  formalize that obstruction instead of adding a false bridge;
- if it survives, prove the aggregate theorem first and make the remaining
  fermionic obligation exactly the block-mean-shift product, not a collection
  of coordinatewise signs;
- do not infer the infinite-tail sign or the order-eighteen companion envelope
  from finite truncations alone.

Non-claims:

- this phase does not prove the aggregate mean-shift sign for every fermionic
  prefix or cutoff;
- it does not prove the nine amplified companion envelopes, order-eighteen
  contraction, dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq_of_meanShiftProduct_nonneg`
  proves the sharp reverse-covariance inequality from positive block masses,
  nonnegative within-block covariance numerators, and one nonnegative product
  of aggregate block-mean shifts. It does not assume nonnegative singleton
  weights or any pointwise cross-block order;
- `affineRankTwoCoefficients_lorentzian_of_meanShiftProduct_nonneg` transports
  the criterion to affine rank-two coefficients, and
  `rankTwoBlockPencil_lorentzian_of_affine_meanShift` exposes the exact
  determinant-pencil interface without division or an inverse;
- the older pointwise-order theorem now factors through the aggregate result,
  making the latter a strict weakening of the previous proof contract rather
  than a parallel lemma;
- exact conditioning through both Phase CV complete pivots checks all 42
  cutoffs from feature counts four through ten. Forty use the new mean-shift
  certificate directly; the two exceptional four-feature endpoint splits
  have a rank-one endpoint and therefore satisfy the discriminant trivially;
- a deterministic extended corpus checks 16,424 exact conditioned pencils
  with prefix depths through six. Every oriented constant, mixed, and
  quadratic coefficient is nonnegative, every discriminant is nonnegative,
  and the smallest exact normalized margin lies strictly between `0.37` and
  `0.38`;
- all three new public theorems use only `propext`, `Classical.choice`, and
  `Quot.sound`, and `./scripts/verify.sh` passes a 3,600-job Lean build and
  122 Python tests with no known dependency vulnerabilities.

Next decision:

The remaining fresh-order-two obligation is now one scalar aggregate sign:
construct, for every conditioned fermionic prefix and cutoff, an affine or
projective chart with positive block masses, nonnegative within-block
covariances, and a nonnegative product of the two block-mean shifts. The most
promising route is a variation-diminishing or continued-fraction statement
for the two conditioned ratio sequences; unlike common chirotopes, it may
permit finitely many reversed coordinate pairs so long as their weighted mean
shifts agree. The Phase CV witness and the 16,424-pencil corpus are mandatory
regressions for that theorem.

## Phase CV contract: complete-pivot chirotope obstruction

State: complete; merged as PR #102 with green CI.

Phase CU refutes nonvanishing on the full ordered chamber but leaves open the
paired-subtuple endpoint geometry observed in the conditioned GECP corpus.
That containment condition is also insufficient, even with the endpoint
selection pattern produced by two complete pivots. A sharper construction
uses the exact rational exponential kernel

`K(i,j) = (101/100)^((i+1)(j+21))`, for `i,j : Fin 4`.

It is sampled by open-unit nodes whose products are integer multiples of
`log(101/100)`. The first complete pivot is `(3,3)`, the next complete pivot of
its Schur residual is `(0,0)`, and the selected subtuples are therefore the two
endpoints on both axes. Nevertheless, the crossed exponential-cubic
determinants for appended degrees `(0,3)` have opposite signs.

Frozen public Lean targets:

- define the real row and column nodes, prove that they are strictly increasing
  in `(0,1)`, and prove that their exponential kernel equals the exact rational
  power matrix above;
- certify `(3,3)` as a complete pivot on the full `Fin 4` domain and `(0,0)` as
  a complete pivot of the exact first residual;
- define the two endpoint-selected crossed augmented matrices with appended
  monomials `1` and `x^3`;
- reduce both augmented determinants by positive column scaling to rational
  `4 x 4` matrices and prove that their signs are opposite.

Independent check:

- recompute all sixteen kernel entries, the first Schur residual, and both
  complete-pivot maxima in exact `Fraction` arithmetic;
- recompute both crossed augmented determinant signs exactly and verify the
  negative product;
- retain the Phase CT randomized corpus to document why its 11,024 successes
  did not cover this narrow boundary-scale witness.

Decision rule:

- if Lean certifies both complete pivots and the opposite augmented signs,
  retire arbitrary-run chirotope synchronization, including its
  paired-subtuple strengthening, as a universal route;
- if the pivot certificate fails, minimize the discrepancy between the exact
  residual and the proposed complete-pivot order before changing the witness;
- do not infer failure of the ultimate determinant envelope: the obstruction
  targets this sign-based sufficient condition, not every cancellation or
  absolute-bound proof.

Non-claims:

- this phase does not refute a quantitative cancellation theorem, an
  absolute-determinant bound, or a bound specialized to later amplification
  coefficients;
- it does not prove higher fresh order, the normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `phaseCVKernelQ` is the exact rational power kernel, while `phaseCVRows`,
  `phaseCVColumns`, and `phaseCV_expKernel_eq` realize it as a genuine
  exponential kernel on strictly increasing open-unit nodes;
- `phaseCV_firstPivot_complete` certifies the bottom-right pivot `(3,3)`, and
  `phaseCV_secondPivot_complete` certifies `(0,0)` after its exact Schur
  update;
- `phaseCVLeftReduced_det_pos` and `phaseCVRightReduced_det_neg` prove the
  opposite signs of the rationally reduced endpoint-selected degree-`(0,3)`
  augmented determinants;
- `phaseCVLeftAugmented_det_pos`, `phaseCVRightAugmented_det_neg`, and
  `phaseCV_crossedAugmented_det_mul_neg` lift the obstruction to the actual
  exponential-monomial determinants through positive column scalings;
- all public results use only `propext`, `Classical.choice`, and `Quot.sound`,
  and `./scripts/verify.sh` passes a 3,600-job Lean build and 121 Python tests
  with no known dependency vulnerabilities.

Next decision:

Retire arbitrary-run common-chirotope synchronization as a universal route.
The remaining rank-two target must keep cancellation between degree pairs
inside the aggregate mixed coefficient or discriminant, or prove an absolute
determinant envelope without coordinatewise signs. The Phase CV exact kernel
is now the mandatory regression for any such replacement: a candidate must
survive its negative degree-`(0,3)` coordinate rather than assuming it away.

## Phase CU contract: ordered-chamber obstruction

State: complete; merged as PR #101 with green CI.

Phase CT isolates nonvanishing of a canonical confluent exponential
determinant along a straight ordered-node homotopy. The tempting strengthening
to the entire ordered positive unit chamber is false. For four increasing row
nodes, two increasing exponential column nodes, and appended monomials `1`
and `x^3`, exact reconnaissance finds opposite determinant signs inside that
chamber. Choosing the column nodes as `9*log(12/11)` and `9*log(11/10)` and
equally spaced rows makes successive exponential ratios rational, so the
obstruction can be proved in Lean without floating-point assumptions.

Frozen public Lean targets:

- define the exponential-cubic augmented matrix on four row nodes and two
  exponential column nodes;
- prove that the row tuples `(1,2,3,4)/9` and `(5,6,7,8)/9` and the two
  logarithmic column nodes are strictly increasing and lie in `(0,1)`;
- reduce both exponential determinants by positive column scaling to exact
  rational matrices, proving positive sign at the first tuple and negative
  sign at the second;
- prove that a continuous scalar path with opposite nonzero endpoint signs
  has a zero, then instantiate it on the straight row-node interpolation.

Independent check:

- recompute the two reduced `4 x 4` determinants in exact `Fraction`
  arithmetic and retain their values `83/194059800` and `-1/194059800`;
- verify symbolically that all straight-path row tuples remain strictly
  increasing in the open unit interval and that the columns remain fixed;
- retain the Phase CT 11,024-case midpoint audit as a guard that the physical
  paired-endpoint corpus is narrower than this counterexample.

Decision rule:

- if Lean certifies the opposite signs and forced zero, retire universal
  ordered-chamber nonvanishing and make the next contract the paired-subtuple
  theorem reflecting `X_sel ⊂ X_aug` and `Y_sel ⊂ Y_aug` at opposite endpoints;
- if the logarithmic normalization fails, keep the exact rational reduced
  matrices and isolate the missing exponential identity before changing the
  witness;
- do not transfer the obstruction to the GECP-correlated path without a
  witness satisfying its two endpoint subtuple constraints.

Non-claims:

- this phase does not refute nonvanishing for the paired-subtuple paths
  generated by successful GECP runs;
- it does not prove arbitrary-run chirotope synchronization, separator
  existence, the infinite-tail passage, higher fresh order, the normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `expCubicAugmented` defines the four-row determinant with two exponential
  columns and appended monomials `1` and `x^3`;
- `phaseCURows_strictMono`, `phaseCUColumns_strictMono`, and the open-unit
  membership theorems certify that both endpoint configurations and the whole
  row path stay inside the proposed ordered chamber;
- `phaseCUReduced_one_det` and `phaseCUReduced_five_det` prove the exact
  rational reduced determinants `83/194059800` and `-1/194059800`;
- `expCubicAugmented_phaseCURows_one_det_pos` and its negative endpoint
  counterpart lift those signs through positive exponential column scaling;
- `Matrix.exists_eq_zero_of_continuousOn_of_mul_neg` supplies the generic
  opposite-sign intermediate-value interface, and
  `phaseCUMatrixPath_exists_det_eq_zero` proves that the explicit straight
  ordered-chamber path contains a singular determinant;
- all public results use only `propext`, `Classical.choice`, and `Quot.sound`,
  and `./scripts/verify.sh` passes a 3,598-job Lean build and 120 Python tests
  with no known dependency vulnerabilities.

Next decision:

Replace the false full-chamber target with the endpoint-coupled statement
actually produced by GECP: `X_sel` is a subtuple of `X_aug`, `Y_sel` is a
subtuple of `Y_aug`, and the two determinants cross these pairs as
`D(X_aug,Y_sel)` and `D(Y_aug,X_sel)`. The Phase CU witness does not satisfy
those constraints, while all 11,024 Phase CT pairs do.

## Phase CT contract: homotopy sign lock

State: complete; merged as PR #100 with green CI.

Phase CS shows that an arbitrary common chirotope schedule is sufficient but
does not construct it.  After tuple reorientation and positive denominator
scaling, the left augmented determinant is a canonical confluent exponential
determinant `D(X_aug,Y_sel;i,j)`, while the transpose of the right determinant
is `D(Y_aug,X_sel;i,j)`.  Componentwise convex interpolation simultaneously
moves `X_aug` to `Y_aug` and `Y_sel` to `X_sel`; strict tuple order and the
positive unit chamber are convex.  Thus continuity plus nonvanishing along
this single path locks the two endpoint signs without an explicit formula.

Frozen public Lean targets:

- prove that a real continuous function nonzero on `[0,1]` has endpoint
  product strictly positive;
- specialize the result to entrywise-continuous square-matrix paths whose
  determinants never vanish;
- keep nonvanishing as an explicit analytic hypothesis, not a disguised sign
  assumption, so a later confluent-Chebyshev theorem can discharge it.

Independent check:

- put both exact augmented determinants into the canonical sorted/scaled form;
- evaluate the componentwise rational midpoint of the row/column node tuples
  for every one of the 11,024 feature pairs;
- verify endpoint normalization, midpoint nonvanishing, and endpoint sign
  agreement while retaining the fixed-schedule obstruction.

Decision rule:

- if the homotopy interface and midpoint audit close, make the next phase only
  the ordered-unit-chamber nonvanishing theorem for the canonical confluent
  exponential determinant;
- if midpoint zeros occur, abandon the straight-line path and search within
  the same ordered chamber before weakening the invariant;
- do not infer all-path nonvanishing from finitely many samples.

Non-claims:

- this phase does not prove homotopy nonvanishing, the arbitrary-run
  chirotope theorem, separator existence, or the infinite-tail passage;
- it does not prove higher fresh order, the normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `Matrix.endpoint_mul_pos_of_continuousOn_of_ne_zero` proves that a
  continuous real path avoiding zero on the unit interval has endpoint
  product strictly positive;
- `Matrix.det_endpoint_mul_pos_of_entrywiseContinuous_of_ne_zero` specializes
  the sign lock to entrywise-continuous square-matrix paths with nonvanishing
  determinant;
- exact canonical normalization and midpoint evaluation succeed for all
  11,024 left/right augmented determinant pairs, preserving the Phase CS
  fixed-schedule obstruction;
- `./scripts/verify.sh` passes a 3,597-job Lean build and 119 Python tests with
  no known dependency vulnerabilities.

The next target is precisely the missing analytic premise: nonvanishing of the
canonical confluent exponential determinant throughout the ordered positive
unit chamber, or a rigorously delimited obstruction if that statement is
false.

## Phase CS contract: common chirotope schedule

State: complete; merged as PR #99 with green CI.

Phase CR turns the conditioned matched coordinate into a product of original
augmented determinants. The initial contract conjectured that their common raw
sign schedule depended only on pivot depth and the two degrees. The exact audit
refutes that stronger claim: already at depth one, pairs `(0,1)` and `(0,2)`
need opposite schedules from the proposed rule. The robust positive-
Grassmannian invariant is that, after removing the independent row/column tuple
orientations, the left and right augmented determinants have the same nonzero
sign. That common sign may depend on the run and pair. Its square still cancels,
while the two tuple orientations multiply to the global Phase CQ orientation.

Frozen public Lean targets:

- prove an abstract common-schedule cancellation lemma for two scalar Plucker
  coordinates with independent reorientations;
- compose it with the Phase CR selected-core-square equivalence to derive the
  conditioned oriented matched-coordinate sign from separate augmented-minor
  sign certificates;
- keep the schedule scalar explicit, nonzero, and independent of the two
  factors so the next confluent-Vandermonde theorem can instantiate it.

Independent check:

- compute tuple orientations of the selected-plus-fresh row and column lists;
- verify that both separately reoriented augmented determinants are nonzero and
  have the same sign for all 11,024 feature pairs in the 500 deterministic
  arbitrary-pivot cases, with both schedule signs occurring;
- retain the failed fixed pivot-depth rule as an explicit obstruction;
- retain the existing core-square, cyclic-order, separator, and Lorentzian
  checks.

Decision rule:

- if a common schedule exists throughout the exact audit, make the next phase
  the analytic left/right chirotope-synchronization theorem via compound total
  positivity or positive-Grassmannian duality;
- if either side needs a different schedule, record the precise discrepancy
  and do not hide it inside a product-level orientation;
- do not infer the analytic sign theorem from the finite audit.

Non-claims:

- this phase does not prove the augmented Taylor determinant sign schedule,
  cyclic-prefix separator existence, or infinite-tail passage;
- it does not prove higher fresh order, the normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered Lean identifiers:

- `Run.mul_nonneg_of_common_sign_schedule`;
- `Run.orientedConditionedRankTwoMatchedCoordinate_nonneg_of_common_sign_schedule`.

Verification:

- the public axiom audit reports only `propext`, `Classical.choice`, and
  `Quot.sound` for both new theorems;
- exact rational reconstruction confirms a shared nonzero left/right schedule
  for all 11,024 augmented pairs, observes both schedule signs, and finds
  counterexamples to the proposed fixed depth-only schedule;
- `./scripts/verify.sh` passes a 3,596-job Lean build and 119 Python tests with
  no known dependency vulnerabilities.

## Phase CR contract: conditioned rank-two Plucker bridge

State: complete; merged as PR #98 with green CI.

Phase CQ isolates oriented signs of conditioned rank-two matched coordinates
as the remaining within-block input.  Those coordinates are determinants of
the two conditioned `2 x 2` feature factors, so the companion Plucker
identities should clear every conditioning denominator at once: multiplication
by the square of the selected-core determinant turns one matched coordinate
into the product of the corresponding original left and right feature-
augmented core determinants.  The square is strictly positive for every
successful run, making oriented sign transport an equivalence.

Frozen public Lean targets:

- prove the exact selected-core-square identity for arbitrary kernels, runs,
  feature factors, two fresh rows/columns, and two feature indices;
- derive unoriented and scalar-oriented nonnegativity equivalences without
  assuming either target sign;
- expose the bridge from a dedicated GECP module and the public theorem
  surface so the conditioned-Taylor specialization can use it directly.

Independent check:

- reconstruct both original augmented feature cores in exact rational
  arithmetic for every feature pair in the 500 deterministic arbitrary-pivot
  conditioned-Taylor cases;
- verify the identity against the already audited conditioned minors, including
  all 11,024 matched coordinates and both pivot-depth orientations.

Decision rule:

- if the generic identity closes, make the next phase a pure augmented-
  generalized-Vandermonde sign/cyclic-interval theorem;
- if feature-index ordering introduces an unexpected determinant sign, record
  it explicitly rather than absorbing it into a sample-dependent orientation;
- do not infer the augmented-minor sign theorem from the finite audit.

Non-claims:

- this phase does not prove the augmented Taylor minor signs, cyclic-prefix
  order, or separator existence;
- it does not prove the infinite-tail passage, higher fresh order, normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered Lean identifiers:

- `Run.selectedCore_sq_mul_conditionedRankTwoMatchedCoordinate`;
- `Run.selectedCore_sq_mul_orientedConditionedRankTwoMatchedCoordinate`;
- `Run.orientedConditionedRankTwoMatchedCoordinate_nonneg_iff`;
- `Run.conditionedRankTwoMatchedCoordinate_nonneg_iff`.

Verification:

- the public axiom audit reports only `propext`, `Classical.choice`, and
  `Quot.sound` for every new theorem;
- exact rational reconstruction verifies both original augmented cores, the
  core-square identity, and oriented sign equivalence for all 11,024 feature
  pairs in the 500 arbitrary-pivot conditioned-Taylor cases;
- `./scripts/verify.sh` passes a 3,596-job Lean build and 119 Python tests with
  no known dependency vulnerabilities.

## Phase CQ contract: matched-minor separator reduction

State: complete; merged as PR #97 with green CI.

Phase CP still asks an analytic caller to restate within-block comonotonicity,
positive masses, and reciprocal cross order in a chosen pole chart. The exact
Phase CN data indicate that the first two are not independent hypotheses.
Matched rank-two coordinates are projective invariants: after the two pole
bases, each coordinate is the transformed singleton-weight product times the
two reciprocal increments. The audit also reveals two independent signs. An
orientation reverses one reciprocal projective circle, while a common chart
scale flips all singleton weights without changing ratios. Thus strict
positivity of the scaled transformed weights turns oriented nonnegative
matched minors directly into `MonovaryOn`; nonempty blocks then supply
positive masses. This is the rank-two positive-Grassmannian/Gale-projection
invariant behind the observed cyclic order.

Frozen public Lean targets:

- prove exact determinant scaling of `affineRankTwoMatchedCoordinate` under
  independently oriented and scaled pole charts;
- derive `MonovaryOn` for the reciprocal ratios from positive transformed
  singleton weights and nonnegative original matched coordinates;
- package a pole-to-Lorentzian theorem whose remaining geometric hypotheses
  are only noncollision, nonzero orientation/scale, positive transformed
  weights, nonempty blocks, oriented within-block matched-coordinate signs,
  and reciprocal cross order;
- keep all hypotheses pointwise and finite so the conditioned-Taylor proof can
  discharge them with generalized-Vandermonde signs.

Independent check:

- for every one of the 2,868 nontrivial Phase CN degree cutoffs, construct
  exact rational left/right poles from gaps in the observed projective orders;
- verify noncollision, strictly positive transformed weights, reciprocal cross
  order, within-block matched-coordinate signs, and the original Lorentzian
  discriminant;
- record the number of candidate gap pairs tried and retain a deliberate pole
  collision as a rejected boundary.

Decision rule:

- if the reduction and constructive audit close, make the next phase solely
  the arbitrary-pivot conditioned-Taylor matched-minor/separator theorem;
- if common cyclic order does not suffice to align positive transformed
  weights, strengthen the analytic target to an oriented cyclic order rather
  than adding sample-specific sign choices;
- do not infer the arbitrary-pivot theorem from the finite construction.

Non-claims:

- this phase does not prove projective matched-minor signs or separator
  existence for arbitrary conditioned Taylor factors;
- it does not prove the infinite-tail passage, higher fresh order, normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered Lean identifiers:

- `Matrix.finTwoOrientedPoleBasis`,
  `Matrix.finTwoScaledOrientedPoleBasis`, and their determinant/block-action
  theorems;
- `Matrix.affineRankTwoMatchedCoordinate_scaledOrientedPoleChart`;
- `Matrix.scaledOrientedPoleChartRatio_monovaryOn_of_matched_nonneg`;
- `Matrix.rankTwoBlockPencil_lorentzian_of_scaled_oriented_poles`.

Verification:

- the public axiom audit reports only `propext`, `Classical.choice`, and
  `Quot.sound` for every new theorem;
- exact rational gap search constructs an oriented and scaled pole pair on
  the first eligible candidate for all 2,868 nontrivial Phase CN cutoffs,
  verifies positive weights, both reciprocal cross orders, matched-coordinate
  invariance, and the Lorentzian discriminant, and retains the collision
  boundary; negative orientations and negative raw chart scales both occur;
- `./scripts/verify.sh` passes a 3,595-job Lean build and 119 Python tests with
  no known dependency vulnerabilities.

## Phase CP contract: separating projective pole chart

State: complete; merged as PR #96 with green CI.

Phase CO permits an arbitrary invertible basis change but leaves a later proof
to construct the changed affine coordinates. For a homogeneous feature column
`a_k * (1,f_k)`, a pole `c` uses the basis
`[[-c,1],[1,0]]`, producing base `a_k*(f_k-c)` and ratio
`1/(f_k-c)`. The analogous right change at pole `d` gives weight
`a_k*b_k*(f_k-c)*(g_k-d)`. This is exactly the chart obtained by placing the
projective point at infinity between the two cyclic degree blocks.

Frozen public Lean targets:

- define the explicit invertible `2 x 2` pole basis and compute its
  determinant;
- prove the exact two-sided pole transformation of a rank-two block matrix;
- prove the transformed affine factorization using only nonvanishing of the
  two pole differences;
- combine the construction with Phases CM--CO to prove the original pencil
  Lorentzian from positivity and block order of the reciprocal pole ratios.

Independent check:

- use exact rational feature factors, poles on both sides and between ratio
  clusters, and all nontrivial degree cutoffs;
- check the block-matrix transformation, transformed weights, reciprocal
  ratios, and discriminant equality directly;
- retain a pole colliding with one projective point as the rejected boundary.

Decision rule:

- if the pole construction closes, the next contract is solely the
  conditioned Taylor cyclic-prefix/separating-pole theorem via confluent
  generalized Vandermonde signs;
- if the explicit chart changes orientation, retain the exact determinant
  sign and use Phase CO rather than adding ad hoc reorientation assumptions;
- do not claim a pole exists merely from the finite audit.

Non-claims:

- this phase does not prove cyclic-prefix order or separating-pole existence
  for arbitrary conditioned Taylor factors;
- it does not prove the infinite-tail passage, higher fresh order, normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered Lean identifiers:

- `Matrix.finTwoPoleBasis` and `Matrix.finTwoPoleBasis_det`;
- `Matrix.finTwoPoleBasis_mul_rankTwoBlockMatrix_mul`;
- `Matrix.poleChartBase`, `Matrix.poleChartRatio`, and
  `Matrix.poleChartBase_mul_poleChartRatio`;
- `Matrix.rankTwoBlockPencil_lorentzian_of_poles`.

Verification:

- the public axiom audit reports only `propext`, `Classical.choice`, and
  `Quot.sound` for every new theorem;
- exact rational checks cover pole charts exterior to either side and between
  ratio values, every one of the five nontrivial cutoffs, exact discriminant
  preservation, positive transformed weights, reciprocal block order, and a
  deliberately rejected pole collision;
- `./scripts/verify.sh` passes a 3,595-job Lean build and 119 Python tests with
  no known dependency vulnerabilities.

## Phase CO contract: projective chart invariance

State: complete; merged as PR #95 with green CI.

Phase CN proves the Lorentzian inequality in any positive ordered affine
chart, while its exact audit shows that the useful chart must depend on the
degree cutoff. Changing the left and right homogeneous bases sends each block
matrix `M` to `P*M*Q`. Both endpoint determinants and the mixed determinant
coefficient should scale by `det(P)*det(Q)`, so the discriminant inequality
scales by its square and is invariant for invertible chart changes.

Frozen public Lean targets:

- recover the mixed `2 x 2` pencil coefficient from determinant evaluation at
  unit weight;
- prove exact mixed-coefficient scaling under arbitrary two-sided basis
  changes;
- prove forward and reverse Lorentzian transfer, with invertibility required
  only for the reverse direction;
- compose the transfer with Phase CN so a positive ordered transformed chart
  certifies the original block pencil.

Independent check:

- use deterministic exact rational invertible basis changes on the Phase CL
  rank-two pencils;
- check endpoint, mixed, and discriminant scaling separately, including
  negative determinant products;
- retain a singular transformation as a one-way boundary case.

Decision rule:

- if chart invariance closes, attack the conditioned Taylor cyclic-prefix
  theorem directly through bordered generalized-Vandermonde signs;
- if mixed-coefficient scaling requires stronger hypotheses, isolate the
  minimal matrix identity rather than constraining the projective chart;
- do not treat basis invariance as existence of a separating chart.

Non-claims:

- this phase does not prove the cyclic-prefix invariant or construct the
  cutoff-specific chart for conditioned Taylor factors;
- it does not prove the infinite-tail passage, higher fresh order, normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered: `finTwoPencilMixed_eq_det_add_sub` recovers the mixed coefficient
from the determinant at unit weight. `finTwoPencilMixed_mul_mul` and
`finTwoPencilDiscriminant_mul_mul` prove exact two-sided basis-change scaling,
including negative and singular determinant products.
`finTwoPencil_lorentzian_mul_mul` gives one-way transfer through arbitrary
changes, while `finTwoPencil_lorentzian_mul_mul_iff` gives equivalence for
invertible changes. Finally,
`rankTwoBlockPencil_lorentzian_of_projective_affine` composes this invariance
with Phase CN: any invertible cutoff-specific chart satisfying the positive
ordered affine certificate proves the original pencil inequality. Exact
rational checks cover 512 invertible changes, including 271 negative scale
products, plus one singular boundary degeneration. The root command passes a
3,595-job Lean build and 118 tests.

## Phase CN contract: projective rank-two ratio bridge

State: complete; merged as PR #94 with green CI.

Phase CM proves the sharp Lorentzian inequality once a rank-two feature pair
is written in a positive affine chart with two comonotone ratio coordinates.
The naive choice of the first fresh row and column as denominators is not
degree-monotone: only 42 of 500 exact conditioned-Taylor cases are monotone in
that chart, and the two aggregate block shifts have opposite signs in 276 of
2,868 nontrivial threshold splits. That global affine-chart target is false.

The projective invariant is substantially cleaner. In all 500 cases the left
and right feature ratios have the same cyclic order, and every one of the
2,868 degree prefixes is a circular interval in both orders. The dominant
orders are the pivot-depth patterns
`p,p-1,...,0,p+1,...,n-1`, up to cyclic rotation and reversal. Thus each
cutoff can plausibly choose its own projective point at infinity, unwrap the
prefix/suffix interval into two ordered affine blocks, and invoke Phase CM.

Frozen public Lean targets:

- define the matched rank-two pair coordinate associated with left/right
  affine ratio factorizations;
- prove its exact product-of-weights/product-of-ratio-increments identity;
- identify the within-block sum with `weightedCovarianceNumerator` and the
  cross-block sum with `weightedCrossIncrement`;
- combine the identities with Phase CM to expose a reusable
  ratio-factorized threshold Lorentzian theorem;
- keep the hypotheses certificate-oriented so a later projective chart can
  instantiate them without division inside the final theorem.

Independent check:

- count the direct-chart monotonicity and mean-shift failures rather than
  silently choosing a favorable sample;
- verify common cyclic order and prefix circular-interval separation for all
  deterministic Phase CL cases and all nontrivial cutoffs;
- retain exact rational arithmetic and boundary feature counts.

Decision rule:

- if the algebraic ratio bridge closes, next prove the pivot-depth cyclic
  order for conditioned Taylor feature factors using generalized
  Vandermonde signs or variation diminution;
- if a prefix fails to be a circular interval, replace the two-block chart by
  a multi-interval Lorentzian decomposition before attempting specialization;
- do not claim the finite audit as a proof of the projective-order theorem.

Non-claims:

- this phase does not prove existence of the chart for arbitrary conditioned
  Taylor or infinite positive-series factors;
- it does not prove higher fresh order, the normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered: `rankTwoMatchedCoordinate_eq_affine` proves the exact chart
factorization of each matched Plucker product.
`affineRankTwoWithinCoefficient_eq_weightedCovarianceNumerator` and
`affineRankTwoCrossCoefficient_eq_weightedCrossIncrement` identify the three
threshold coefficients with the Phase CM objects. `rankTwoBlockMatrix_det`
and `finTwoPencilMixed_rankTwoBlockMatrix` connect those scalar coefficients
to the original `2 x 2` block matrices, and
`rankTwoBlockPencil_lorentzian_of_affine` exposes the final division-free
certificate consumed by Phase CL. The exact audit confirms sign-coherent
direct-chart weights in all 500 cases but degree monotonicity in only 42 and
opposite mean shifts in 276 of 2,868 splits. In contrast, all 500 left/right
configurations have the same unoriented cyclic order and every one of the
2,868 prefixes is a circular interval. The root command passes a 3,595-job
Lean build and 117 tests.

## Phase CM contract: two-block Lorentzian covariance

State: complete; merged as PR #93 with green CI.

Phase CL reduces fresh-order-two real-rootedness to
`4*c0*c2 ≤ c1^2`. The exact audit reveals an additional structural fact:
every oriented two-feature coordinate product is nonnegative. A rank-two
factor with a nonzero reference coordinate writes each Plucker coordinate as
a positive pair weight times a difference of ratios. Thus the three pencil
coefficients become the within-lower, cross-block, and within-upper weighted
increment sums of two comonotone sequences.

The relevant reverse covariance inequality has a division-free proof. If
`mL,mU` are the two block masses, `c0,c2` their covariance numerators, and
`c1` the cross-block increment sum, then

`mL*mU*c1 = mU^2*c0 + mL^2*c2 + rowShift*columnShift`.

Ordered blocks make the final product nonnegative. The first two terms then
imply `4*c0*c2 ≤ c1^2` because their AM--GM remainder is one square.

Frozen public Lean targets:

- define the weighted cross-block increment sum for two finite index blocks;
- prove its exact division-free two-block covariance decomposition;
- prove nonnegativity of the two within-block covariance numerators and both
  mean-shift factors from nonnegative weights and internal/cross-block order;
- derive the sharp reverse covariance inequality
  `4*covLower*covUpper ≤ crossIncrement^2` for positive block masses;
- expose a ratio-coordinate corollary ready to consume rank-two conditioned
  feature factors.

Independent check:

- exhaust small positive integer weights and two independently comonotone
  rational sequences, checking the decomposition and sharp inequality;
- extend the Phase CL audit to count oriented individual two-feature products,
  distinguishing the true fresh-order-two coordinate sign from the false
  higher-order coordinatewise claim;
- retain boundary cases with singleton blocks, repeated sequence values, and
  zero within-block covariance.

Decision rule:

- if the reverse covariance theorem closes, next prove that arbitrary
  conditioned Taylor rank-two feature factors admit its positive-weight,
  comonotone ratio representation;
- if individual fresh-order-two coordinate products can be negative, weaken
  the bridge to grouped ratio intervals before attempting the fermionic
  specialization;
- do not extrapolate the rank-two result to higher compound order.

Non-claims:

- this phase does not prove the conditioned-factor ratio representation or
  the full fermionic fresh-order-two theorem;
- it does not extend the Lorentzian inequality to fresh order three or pass
  from finite Taylor features to the infinite positive-series tail;
- it does not prove the normalized tangent envelope, order-eighteen
  contraction, dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered: `weightedCrossIncrement_eq` expands the mixed block coefficient,
`weightedMass_mul_weightedMass_mul_crossIncrement` proves the division-free
two-block law of total covariance, and
`four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq` derives the
sharp reverse covariance inequality from ordered nonnegative weights.
`geometricCrossIncrement_lorentzian` specializes the result to geometric
ratio coordinates on a degree prefix/suffix split. An exhaustive exact audit
checks the identity and inequality on 1,764 weighted splits, including
singleton and zero-covariance boundaries. The conditioned-Taylor audit checks
11,024 oriented two-feature coordinate products without a negative case.
The root command passes a 3,594-job Lean build and 117 tests.

## Phase CL contract: fresh-order-two Lorentzian bridge

State: complete; merged as PR #92 with green CI.

Phase CK isolates nonnegative generalized spectrum as the promising compound
invariant. At fresh order two the threshold determinant is quadratic, so this
claim has a sharp scalar form. If its oriented coefficients are
`c0,c1,c2`, real nonpositive roots are equivalent (away from degenerate
endpoints) to the Lorentzian/Alexandrov--Fenchel inequality
`4*c0*c2 ≤ c1^2`. A first attempt to force this from one row or column
cross-minor sign is false on most exact samples; the proof must preserve the
full mixed-minor sum, as in a Hodge--Riemann or mixed-discriminant argument.

Frozen public Lean targets:

- define the constant, mixed, and quadratic coefficients of a `2 x 2`
  determinant pencil and prove its exact polynomial expansion;
- identify the polynomial discriminant with `c1^2 - 4*c0*c2` and expose the
  exact Lorentzian inequality required for real roots;
- construct nonnegative spectral factors from nonnegative coefficients,
  positive quadratic coefficient, and the Lorentzian inequality, proving the
  factorization for every scalar weight;
- compose that scalar factorization with the matrix-pencil expansion, without
  introducing an inverse or an eigenvalue axiom.

Independent check:

- audit deterministic exact finite-Taylor complete-pivot prefixes through
  fresh order two, checking every threshold coefficient and discriminant;
- evaluate both tempting single row- and column-cross-minor sufficient signs,
  and retain their failure count so they cannot be mistaken for the missing
  theorem;
- independently substitute the constructed quadratic factors back into each
  audited polynomial using exact arithmetic plus certified real-square-root
  identities in Lean.

Decision rule:

- if the exact Lorentzian inequality survives while the single-cross-minor
  criteria fail, next seek a compound mixed-discriminant/Hodge--Riemann proof;
- if the inequality fails, minimize the exact fresh-order-two witness and
  retreat to the nine finite Phase BS weights;
- do not infer arbitrary fresh order or the infinite-tail result from a
  quadratic theorem.

Non-claims:

- this phase does not prove the fermionic Lorentzian inequality itself or the
  compound grouped-sign theorem beyond the stated discriminant premise;
- it does not prove the infinite-tail passage, normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `finTwoPencilMixed` names the mixed coefficient and
  `finTwoPencil_det_expansion` proves the exact `2 x 2` determinant quadratic;
- `finTwoPencilDiscriminant` and
  `finTwoPencilDiscriminant_nonneg_iff` identify real-rootedness with the sharp
  Lorentzian coefficient inequality;
- `nonnegative_quadratic_factorization` constructs two nonnegative affine
  spectral offsets using the exact square root of the discriminant, and
  `finTwoPencil_nonnegative_factorization` composes that construction with the
  matrix pencil without an inverse;
- `finTwoPencil_det_nonneg_of_lorentzian` gives the resulting sign certificate
  at every nonnegative weight;
- 500 deterministic exact finite-Taylor runs check 3,868 fresh-order-two
  threshold pencils. Every oriented coefficient is nonnegative and every
  discriminant satisfies `4*c0*c2 ≤ c1^2`;
- only 1,030 pencils satisfy at least one naive row/column cross-minor sign, so
  2,838 exact cases rule out that simpler proof route while retaining the
  full Lorentzian inequality.
- `./scripts/verify.sh` passes a 3,594-job Lean build, the expanded public
  axiom audit, 117 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #92.

Next decision:

Express the three coefficients as mixed compound pairings and prove their
Lorentzian inequality by a reverse Cauchy--Schwarz/Hodge--Riemann argument.
The exact target is now scalar and fresh-order-two; no eigenvalue API is
needed until that mixed-minor inequality is understood.

## Phase CK contract: compound generalized-spectrum reduction

State: complete; merged as PR #91 with green CI.

Phase CJ proves the complete one-pivot grouped theorem, but its scalar
conditional-covariance argument does not extend coordinatewise to feature
subsets.  The natural subset-lattice route is also too strong: exact Phase BV
data violate absolute-value `MTP2`.  A different invariant survives.  For a
fixed threshold, the grouped coefficients form the determinant polynomial of
the conditioned low/high matrix pair; exact compound samples indicate that
this polynomial has only nonpositive real zeros.  This is the
Pólya-frequency/generalized-eigenvalue analogue of scalar association.

Frozen public Lean targets:

- factor a determinant pencil when its low and high matrices admit common
  left and right factors with diagonal middle coordinates;
- prove nonnegativity, and a sharp strict-positivity variant, for every
  nonnegative pencil weight from nonnegative diagonal generalized-spectrum
  data and the orientation of the common factors;
- expose the result without an explicit matrix inverse, so a future
  oscillatory-matrix or generalized-eigenvalue theorem can discharge exactly
  the remaining simultaneous-factorization premise.

Independent check:

- on the exact Phase BV 56-coordinate witness, evaluate every threshold
  polynomial, its Newton log-concavity inequalities, and its quadratic or
  cubic discriminant using `Fraction` arithmetic;
- count failures of absolute-value `MTP2` on the ordered three-subset lattice;
- repeat the oriented coefficient, log-concavity, discriminant, and `MTP2`
  audit on deterministic exact compound Taylor samples through fresh order
  three, using factorization and determinant solves only, never an explicit
  inverse.

Decision rule:

- if the oriented threshold polynomials remain nonnegative and real-rooted
  while `MTP2` fails, next target a sign-regular or oscillatory generalized-
  eigenvalue theorem for the conditioned low/high pair;
- if real-rootedness fails, retain only the common-factor certificate and use
  the minimized exact witness to choose a weaker sector or finite-weight
  invariant;
- do not infer all-order compound positivity from the finite audit.

Non-claims:

- this phase does not construct the simultaneous diagonal factorization for
  arbitrary fermionic prefixes or prove the compound grouped-sign theorem;
- it does not prove the infinite-tail passage, normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `det_add_smul_eq_prod_of_common_diagonal` gives a division-free determinant
  factorization for a matrix pair with common left/right factors and diagonal
  middle coordinates;
- `det_add_smul_nonneg_of_common_diagonal` and the two strict variants turn
  nonnegative generalized-spectrum coordinates and factor orientation into an
  all-nonnegative-weight determinant sign;
- the exact Phase BV witness has 372 failures among 3,136 absolute-value
  `MTP2` comparisons, while all nine threshold coefficient sequences are
  nonnegative, log-concave, and have nonnegative exact discriminant;
- 200 deterministic exact compound Taylor cases produce 1,481 oriented
  threshold polynomials through fresh order three. Every polynomial retains
  nonnegative coefficients, Newton log-concavity, and nonnegative exact
  discriminant, although 147 cases violate absolute-value `MTP2`;
- the audit therefore rejects the subset-lattice FKG route and isolates the
  missing theorem as nonnegative real generalized spectrum, or an equivalent
  common diagonal factorization, for the conditioned low/high pair.
- `./scripts/verify.sh` passes a 3,593-job Lean build, the expanded public
  axiom audit, 116 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #91.

Next decision:

Prove the generalized-spectrum statement first at fresh order two. The exact
criterion is a nonnegative quadratic discriminant for every cutoff; a
sign-regular matrix-pencil or oscillatory generalized-eigenvalue argument can
then be tested before attempting arbitrary compound order.

## Phase CJ contract: weighted centered-group association

State: complete; merged as PR #90 with green CI.

Phase CI's positive-coefficient polynomials are the first instance of a more
general association principle.  At one pivot, normalize the positive pivot
feature products into weights `w_k`; the row and column feature ratios become
two increasing sequences `f_k,g_k`.  Every matched bordered product is a
positive common scale times
`w_k (f_k-E f)(g_k-E g)`.  A prefix or suffix sum is nonnegative because its
conditional covariance is nonnegative and both conditional means shift from
the global means in the same direction.

Frozen public Lean targets:

- define a finite weighted covariance numerator and prove its nonnegativity
  for monovary sequences by an exact pairwise-sum identity;
- prove the algebraic two-block decomposition of a globally centered product
  sum into a local covariance numerator and the product of two block-mean
  shifts;
- prove nonnegativity for a lower block and an upper block when both sequences
  have the corresponding cross-block order;
- specialize the result to positive weights and two increasing geometric
  feature-ratio sequences on arbitrary finite degree sets;
- expose the resulting arbitrary-feature one-pivot prefix and suffix grouped
  signs in a form composable with the Phase CG coefficient bridge.

Independent check:

- verify the pairwise covariance identity and two-block decomposition in exact
  rational arithmetic on exhaustive small weighted sequences;
- test arbitrary feature counts, cutoffs, positive rational weights, and
  increasing geometric ratios, including groups containing negative
  individual centered products;
- compare the generic centered products with direct one-pivot Taylor bordered
  determinants for deterministic physical and super-physical node grids.

Decision rule:

- if the weighted association theorem closes arbitrary finite one-pivot
  groups, next seek its compound/exterior-power analogue for multiple fresh
  rows and columns;
- if Taylor bordered products fail to match the centered form, isolate the
  missing factorial or pivot scaling before weakening the theorem;
- do not infer the higher compound-order sign or infinite-tail passage from a
  one-pivot theorem.

Non-claims:

- this phase does not prove grouped signs for more than one fresh row/column;
- it does not prove the infinite-tail passage, normalized tangent envelope,
  residual half contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `two_mul_weightedCovarianceNumerator` proves the exact pairwise-increment
  formula for weighted covariance, and
  `weightedCovarianceNumerator_nonneg` gives a weighted Chebyshev theorem for
  arbitrary finite monovary data;
- `weightedMass_mul_weightedCenteredGroup` is the division-free conditional-
  covariance decomposition; `weightedCenteredLowerGroup_nonneg` and
  `weightedCenteredUpperGroup_nonneg` prove both ordered block signs from
  local covariance and same-direction mean shifts;
- `geometricCenteredProduct_prefix_nonneg` and
  `geometricCenteredProduct_suffix_nonneg` specialize the abstraction to
  arbitrary finite positive weights and two geometric sequences with ratios
  at least one;
- `onePivotTaylorProduct_eq_geometricCenteredProduct` identifies every
  finite-Taylor bordered product with that centered model, including all
  factorial and pivot scalings;
- `onePivotTaylorProduct_prefix_nonneg` and
  `onePivotTaylorProduct_suffix_nonneg` prove every nonempty finite one-pivot
  Taylor prefix and suffix group nonnegative at arbitrary ordered positive
  nodes, for every feature count;
- `onePivotTaylorProduct_step_nonneg` combines those groups to prove every
  nontrivial two-level contiguous step nonnegative for arbitrary nonnegative
  scalar weight;
- exact arithmetic exhausts 4,644 small weighted groups and checks 6,782
  Taylor prefix/suffix groups across 512 randomized cases through 12 features;
  all group signs survive despite 218 negative individual coordinates;
- `./scripts/verify.sh` passes a 3,592-job Lean build, the expanded public
  axiom audit, 115 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #90.

Next decision:

Lift the centered-association model from single degrees to feature subsets.
For multiple fresh rows and columns, the natural index space is the ordered
subset lattice; oriented generalized-Vandermonde minors supply positive base
weights and Schur-like ratio functions.  The promising analogue of scalar
Chebyshev is multivariate total positivity (`MTP2`) or an FKG association
theorem for the resulting biorthogonal ensemble, conditioned on tail count.
First identify the exact subset weight and test its log-supermodular minors on
the Phase BV witness before committing to a formal compound theorem.

## Phase CI contract: three-feature grouped sign theorem

State: complete; merged as PR #89 with green CI.

Phase CH exhibits the smallest observed coordinate cancellation but also shows
that both neighboring groups stay positive.  Symbolic elimination for general
ordered nodes reveals the reason: each adjacent two-coordinate sum factors as
the product of the row gap, the column gap, and a polynomial with positive
coefficients.  This is a concrete variation-diminishing theorem rather than a
finite witness.

Frozen public Lean targets:

- define the generic three-feature positive-Taylor kernel and its one-pivot
  left, right, and matched bordered coordinates;
- prove closed factorizations of the three individual matched coordinates;
- prove the adjacent prefix and suffix sums equal the ordered-gap product
  times their explicit positive-coefficient polynomials;
- prove both grouped sums strictly positive for arbitrary strictly ordered
  nonnegative row and column nodes;
- deduce positivity of both nontrivial contiguous-step determinants for every
  nonnegative scalar weight.

Independent check:

- compare the symbolic closed forms with direct exact bordered determinants;
- exhaust a deterministic rational boundary/interior grid, including zero
  endpoints and nodes beyond the physical row interval;
- use deterministic randomized exact cases to verify both factorizations and
  all-weight consequences independently of Lean.

Decision rule:

- if the positive-coefficient factorizations formalize, seek a recurrence or
  sign-reversing involution that extends consecutive grouped sums to arbitrary
  feature count and cutoff;
- if either factorization needs an unstated bound, record the exact domain and
  test whether the physical `0 < x <= 1`, `0 < y` hypotheses suffice;
- do not infer arbitrary-prefix, arbitrary-fresh-order, or infinite-tail signs
  from the one-pivot order-three theorem.

Non-claims:

- this phase does not prove the grouped sign theorem for more than one fresh
  row/column or more than three features;
- it does not prove the infinite-tail passage, normalized tangent envelope,
  residual half contraction, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2.

Delivered result:

- `positiveTaylorThreeProductZero_factor`,
  `positiveTaylorThreeProductOne_factor`, and
  `positiveTaylorThreeProductTwo_factor` compute all three generic one-pivot
  matched bordered coordinates;
- `positiveTaylorThree_prefix_factor` and
  `positiveTaylorThree_suffix_factor` factor the two adjacent groups into the
  ordered row gap, ordered column gap, and explicit polynomials whose
  coefficients are all positive;
- `positiveTaylorThree_prefix_pos` and
  `positiveTaylorThree_suffix_pos` prove both groups strictly positive for
  arbitrary strictly ordered nonnegative nodes, without the physical
  restriction `x <= 1`;
- `positiveTaylorThree_firstStep_pos` and
  `positiveTaylorThree_secondStep_pos` deduce strict positivity of both
  nontrivial contiguous-step determinants for every nonnegative scalar
  weight, with the necessary positive selected-node condition at the first
  cut when the weight can vanish;
- exact arithmetic verifies every coordinate and group factorization on 150
  boundary/interior grid cases and 2,048 deterministic randomized rational
  cases, including zero endpoints and row nodes larger than one;
- `./scripts/verify.sh` passes a 3,588-job Lean build, the expanded public
  axiom audit, 114 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #89.

Next decision:

Lift from three features to arbitrary finite one-pivot Taylor truncations by
normalizing the pivot feature weights into a finite probability measure.  Each
matched coordinate becomes a positive common factor times
`(r^k-E[r^J])(s^k-E[s^J])` for `r,s>1`.  For a prefix or suffix, conditional-
mean decomposition writes its sum as conditional covariance plus a product of
two same-sign mean shifts.  Monotonicity makes both terms nonnegative.  This
probabilistic association identity is the most promising route to the full
one-pivot grouped theorem before attempting compound orders.

## Phase CH contract: minimal grouped cancellation witness

State: complete; merged as PR #88 with green CI.

Phase CG reduces the analytic problem to signs of tail-count groups.  A search
for the smallest cancellation mechanism finds it already with one pivot, one
fresh row and column, and three consecutive positive-Taylor features at
physical-domain nodes `x=(1/8,1/4)`, `y=(6,22)`.  The individual bordered
products by degree are `561/64`, `-7/256`, and `801/32`; every contiguous
prefix and suffix sum is nevertheless positive.  This isolates adjacent-
degree compensation without the combinatorics of a long prefix.

Frozen public Lean targets:

- define the exact two-node, three-feature Taylor factors and their finite
  base kernel;
- compute the three one-pivot bordered Plucker products and formally certify
  the negative middle coordinate;
- prove positivity of every contiguous prefix and suffix group at the two
  nontrivial cutoffs;
- prove the resulting scalar contiguous-step determinant positive for every
  nonnegative weight, despite the negative coordinate.

Independent check:

- recompute all factor entries, base-kernel entries, conditioned factors,
  bordered products, grouped sums, and step evaluations in exact `Fraction`
  arithmetic;
- exhaust a deterministic small rational grid for two- and three-feature
  one-pivot witnesses, confirming that three features are the first observed
  coordinate-level cancellation in the physical row interval `(0,1]`;
- retain the Phase CF sparse-projection failures and Phase CG grouped
  identities as regression guards.

Decision rule:

- if the exact witness formalizes, target a prefix/suffix-sum or adjacent-
  degree sign-reversing theorem for one-pivot Taylor features before lifting
  to arbitrary prefixes and fresh order;
- if either grouped sum is nonpositive, abandon coefficientwise grouped
  positivity and return to the nine finite evaluations;
- do not infer the all-orders grouped theorem from this witness alone.

Non-claims:

- this phase does not prove grouped signs for arbitrary nodes, prefixes, or
  fresh order, nor the infinite-tail passage or normalized tangent envelopes;
- it does not prove the residual half contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `phaseCHKernel_eq` computes the exact three-feature Taylor kernel at
  `x=(1/8,1/4)` and `y=(6,22)`;
- `phaseCHBorderedProduct_zero`, `phaseCHBorderedProduct_one`, and
  `phaseCHBorderedProduct_two` compute the three bordered products as
  `561/64`, `-7/256`, and `801/32`, while
  `phaseCH_middle_borderedProduct_neg` formally certifies the coordinatewise
  obstruction;
- `phaseCH_proper_contiguous_groups_pos` proves all proper nonempty prefix and
  suffix groups positive, and `phaseCH_firstStep_pos` and
  `phaseCH_secondStep_pos` prove both nontrivial contiguous-step determinants
  positive for every nonnegative weight;
- exact `Fraction` arithmetic independently recomputes the factors, kernel,
  conditioned factors, bordered products, groups, and step evaluations;
- an exhaustive ordered grid of 6,468 physical-row/positive-column cases has
  no two-feature negative coordinate, while 800 three-feature cases do; the
  formal witness is the first in lexicographic enumeration;
- `./scripts/verify.sh` passes a 3,587-job Lean build, the expanded public
  axiom audit, 113 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #88.

Next decision:

Derive the one-pivot three-feature products symbolically and prove the two
adjacent grouped inequalities on the full ordered physical domain.  The
promising mechanism is variation diminution in the feature degree: the middle
coordinate may reverse sign, but its magnitude is controlled by a neighboring
endpoint coordinate.  A successful order-three theorem should then be stated
as a prefix/suffix-sum invariant that can plausibly lift by compound-matrix or
sign-reversing-involution arguments, rather than as the already false
coordinatewise positivity claim.

## Phase CG contract: grouped conditioned Plucker bridge

State: complete; merged as PR #87 with green CI.

Phase CF narrows the viable diagonal actions to the contiguous two-level
weights used by the physical tail amplification.  The existing block-feature
polynomial already groups Cauchy--Binet terms by the number of tail features,
while the companion Plucker API represents each conditioned term by original
bordered determinants.  Their composition is the exact algebraic statement
needed before attacking the remaining grouped signs.

Frozen public Lean targets:

- prove the generic squared-core identity for one conditioned signed ordered
  Cauchy--Binet term and its two original feature-bordered determinants;
- sum that identity at one fixed tail count to express every conditioned
  block-feature coefficient as a grouped bordered Plucker pairing;
- sum with the powers of a scalar step weight to express the entire
  conditioned contiguous-step determinant as the corresponding weighted
  bordered pairing;
- expose the exact nonnegativity premise that would imply every required
  finite step-weight sign, without assuming individual coordinate signs.

Independent check:

- on the Phase BV witness, compare every grouped coefficient computed from
  conditioned feature factors with the same grouping of bordered products;
- verify that the two negative sparse coordinates are absorbed into positive
  grouped coefficients for each of the nine contiguous cutoffs;
- compare direct step-weight determinants with evaluation of the grouped
  coefficient polynomial at all Phase BS nodes.

Decision rule:

- if all grouped coefficients remain nonnegative, make their analytic sign
  the next formal theorem and reuse the existing rank-amplification API;
- if a grouped coefficient is negative but the nine required evaluations stay
  positive, target the finite-node envelope directly and retain the
  coefficient obstruction;
- do not strengthen the claim back to coordinatewise or arbitrary-diagonal
  positivity.

Non-claims:

- this phase does not prove the grouped coefficient signs or the nine
  normalized tangent envelopes;
- it does not prove the residual half contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `Run.selectedCore_sq_mul_conditionedSignedTerm` gives the generic one-term
  squared-core identity for arbitrary real feature factors;
- `Run.selectedCore_sq_mul_conditionedBlockFeatureCoefficient` sums that
  identity at one fixed tail count and proves that the conditioned polynomial
  coefficient is exactly the grouped original bordered Plucker pairing;
- `Run.selectedCore_sq_mul_conditionedStep_det_eq_sum` evaluates the whole
  contiguous-step determinant as the same bordered products weighted by
  `lambda` to their tail counts;
- `Run.conditionedBlockFeatureCoefficient_nonneg_of_grouped` proves that a
  nonnegative grouped bordered pairing is sufficient for the corresponding
  conditioned coefficient, using the already certified nonzero selected-core
  determinant rather than adding an assumption;
- in the Phase BV witness, every grouped coefficient is nonnegative for all
  nine cutoffs `0,...,8`, despite the two negative individual coordinates;
- all 90 cutoff/node determinant evaluations exactly equal the associated
  grouped coefficient polynomial and remain strictly positive;
- `./scripts/verify.sh` passes a 3,586-job Lean build, the expanded public
  axiom audit, 112 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #87.

Next decision:

Prove nonnegativity of the grouped bordered pairing for finite positive-Taylor
features with a contiguous degree cutoff, then pass it to the infinite tail.
The most promising route is a discrete variation-diminishing or sign-
reversing-involution theorem indexed by tail count: it must pair the sparse
negative coordinates inside each group while preserving the selected/fresh
cross-inversion orientation.

## Phase CF contract: sparse diagonal obstruction

State: complete; merged as PR #86 with green CI.

Phase CE's dense exact census suggested that arbitrary nonnegative diagonal
weights might preserve the distinguished orientation.  Phase BV already
contains the decisive stress test: two matched conditioned Plucker
coordinates are negative.  A diagonal projection supported exactly on either
feature subset isolates that coordinate, so the broad nonnegative diagonal-
action conjecture is false.  The physical amplification is much narrower: its
weights are the contiguous two-level step `1,...,1,lambda,...,lambda`.

Frozen public Lean targets:

- define the rectangular inclusion matrix attached to an injective ordered
  feature tuple and prove its left and right multiplication identities;
- define the induced sparse feature projection and prove it is entrywise
  nonnegative with zero off-diagonal entries;
- prove that inserting this projection between arbitrary feature factors
  isolates exactly the selected matched Cauchy--Binet coordinate in the final
  determinant;
- conclude abstractly that any negative matched coordinate produces a
  nonnegative diagonal-action counterexample.

Independent check:

- revisit the exact Phase BV witness and identify the negative feature subsets
  `(0,1,5)` and `(0,2,5)`;
- verify directly that the two corresponding sparse projections make the
  conditioned determinant negative;
- enumerate every contiguous two-level cutoff for representative physical
  weights, including all nine Phase BS nodes, and retain their positive signs
  as the narrower next target;
- record why the Phase CE dense randomized census missed the sparse boundary.

Decision rule:

- retire arbitrary nonnegative diagonal action once the sparse isolation
  theorem and witness agree;
- pursue only the contiguous step-weight cone used by positive-series tail
  amplification, preferably through grouped Plucker sums or a variation-
  diminishing theorem rather than coordinatewise signs;
- if a contiguous step-weight counterexample appears, minimize it before
  proposing any new positivity theorem.

Non-claims:

- this phase does not prove positivity for contiguous step weights or any of
  the nine normalized tangent envelopes;
- it does not prove the residual half contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `orderedFeatureInclusion` is the one-hot rectangular inclusion attached to
  a feature tuple; `mul_orderedFeatureInclusion` and its transposed companion
  identify multiplication with the corresponding matrix submatrices;
- `sparseFeatureProjection` is the resulting feature-space projection;
  `sparseFeatureProjection_nonneg` proves every entry nonnegative and
  `sparseFeatureProjection_offDiagonal_eq_zero` proves all off-diagonal
  entries vanish;
- `det_mul_sparseFeatureProjection_mul` proves that inserting the projection
  isolates exactly the product of the selected left and right minors;
- `mul_sparseFeatureProjection_mul_det_neg` turns any negative matched
  Plucker coordinate into a formal sign counterexample for the projected
  determinant;
- in the exact Phase BV witness, the feature subsets `(0,1,5)` and `(0,2,5)`
  are the only negative coordinates among 56, and both corresponding sparse
  projected determinants are negative;
- all 90 combinations of the nine Phase BS amplification nodes plus weight
  one with the nine contiguous feature cutoffs `0,...,8` remain strictly
  positive in the same witness;
- the Phase CE randomized census used dense independent weights and therefore
  almost never sampled an exact sparse support; it was evidence, not a proof,
  and the sparse projection is the decisive adversarial test;
- `./scripts/verify.sh` passes a 3,585-job Lean build, the expanded public
  axiom audit, 112 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #86.

Next decision:

Target the contiguous step-weight determinant directly. Its Cauchy--Binet
expansion groups Plucker coordinates by the number of indices above the
cutoff, so the plausible invariant is positivity of every grouped tail-count
coefficient, not positivity of each coordinate. Phase BV's graded pairing
and the mixed-column coefficient API already provide the algebraic container;
the missing work is an analytic sign theorem for those grouped conditioned
coefficients.

## Phase CE contract: conditioned tangent cancellation collapse

State: complete; merged as PR #85 with green CI.

Phase CD proves that separate left, middle, and right flag signs cannot explain
the positive distinguished tangent determinant.  The fixed-run companion API
already contains the correct cancellation-preserving object: when
`K = A*B` and `H = A*W*B`, two-sided conditioning sends `H` to the single
product `conditionRows(A) * W * conditionColumns(B)`.  Combining this with
the squared-core tangent identity regroups the complete double-flag sum before
any sign estimate is attempted.

Frozen public Lean targets:

- prove that every distinguished tangent minor of a shared factorization
  `K=A*B`, `H=A*W*B` is `det(selectedCore)^2` times the determinant of the
  sampled conditioned feature product;
- combine this identity with ordered triple Cauchy--Binet to identify the
  entire signed double-flag sum with that single conditioned determinant;
- specialize the collapse to the positive-series amplified kernel so the nine
  Phase BS nodes are expressed directly through conditioned low/tail feature
  weights, without separate flag-orientation hypotheses;
- publicly export and axiom-audit the new identities.

Independent check:

- exhaustively or deterministically test generalized positive monomial feature
  systems with independently ordered row and column grids, shuffled selected
  prefixes, sorted fresh samples, and arbitrary nonnegative diagonal weights
  in exact rational arithmetic;
- minimize and retain an exact counterexample once a diagonal weight is
  allowed to be negative, while guarding against the spurious sign change
  caused by presenting a fresh tuple out of canonical sorted order;
- independently compare the double-flag sum, direct tangent determinant, and
  conditioned-product determinant on the Phase CD witness.

Decision rule:

- if the nonnegative-weight census stays sign coherent, formulate the next
  analytic target as preservation of the distinguished orientation under a
  nonnegative diagonal action on two totally positive feature flags;
- if a nonnegative counterexample appears, minimize it and replace the target
  by the narrower step weight `1,...,1,lambda,...,lambda` used by the
  fermionic amplification;
- do not return to termwise double-flag positivity or assume the desired
  determinant sign.

Non-claims:

- this phase does not prove the nonnegative diagonal-action positivity theorem or
  any of the nine normalized tangent envelopes;
- it does not prove the residual half contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `Run.upperTangentMatrix_factorized_distinguishedMinor_det` proves that the
  distinguished tangent minor for `K=A*B`, `H=A*W*B` is the squared selected-
  core determinant times the sampled determinant of
  `conditionRows(A) * W * conditionColumns(B)`;
- `Run.signature_mul_conditionedFeature_det_eq_doubleFlagSum` reindexes the
  complete double border by one common equivalence and proves that the whole
  ordered double-flag sum equals that single conditioned determinant, with no
  permutation sign because both determinant axes use the same reindexing;
- `expTaylor_upperTangent_diagonal_distinguishedMinor_det` specializes the
  cancellation-preserving identity to finite positive-Taylor monomial flags
  and arbitrary diagonal degree weights;
- all three public theorems use only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- deterministic exact rational arithmetic checks 5,000 generalized monomial
  cases with arbitrary nonnegative weights, independently ordered grids,
  shuffled prefixes, and fresh orders through three: 3,964 are strictly
  oriented-positive and 1,036 vanish for rank or zero-weight reasons;
- allowing a negative degree weight immediately fails: weights `(-2,0)` on
  the two-feature sample give the conditioned scalar `-4/169`;
- the Phase CD witness independently satisfies the new collapse identity, so
  its two negative double flags cancel inside the positive conditioned
  determinant rather than disappearing under a termwise argument;
- `./scripts/verify.sh` passes a 3,584-job Lean build, the expanded public
  axiom audit, 112 Python tests, formatting, typing, and the dependency audit.
  The change passed CI and merged as PR #85.

Next decision:

Prove the oriented nonnegativity of the conditioned product determinant for
two totally positive monomial flags and a nonnegative diagonal action.  The
right abstraction is a cancellation-preserving compound/exterior statement,
or equivalently a planar-network theorem for an oblique Schur projection; it
must conclude the sign of the whole determinant and must not decompose it into
the false Phase CA termwise sign premises.

Phase CF subsequently refutes this target for arbitrary nonnegative diagonal
actions by activating the sparse supports omitted by the dense census. Only
the contiguous two-level step weights used by the physical amplification
survive as the next sign family.

## Phase CD contract: factorwise orientation obstruction

State: complete; merged as PR #84 with green CI.

Phase CC closes the middle matrix's all-orders total nonnegativity, but the
Phase CA criterion additionally requires the three factor orientations to
multiply to one external sign for every pair of intermediate flags.  Exact
enumeration of the smallest distinguished positive-monomial tangent example
shows that this stronger termwise condition fails: with two physical nodes,
three features, positive diagonal weights, and one selected plus one fresh
row and column, 12 double-flag terms are nonzero, 10 are positive, and 2 are
negative, while their complete sum is positive.  The two-feature case has no
mixed sign, making three features the first obstruction in this family.

Frozen public Lean targets:

- define the exact two-node, three-feature positive monomial factors,
  positive diagonal tangent weights, distinguished row/column sample, and
  explicit negative double flag;
- prove the left, middle, and right flag minors and their negative product by
  exact rational normalization;
- prove that no left/middle/right orientation functions can satisfy Phase
  CA's factorwise nonnegativity interface with external signature `1` on this
  witness;
- compute the complete distinguished tangent determinant and prove it is
  positive, separating failure of the proof interface from failure of the
  desired global sign.

Independent check:

- enumerate every double flag in exact `Fraction` arithmetic for feature
  counts two and three;
- verify that the two-feature witness has two positive nonzero terms, while
  the three-feature witness has exactly 12 nonzero terms, 10 positive and 2
  negative, summing to `24929 / 27648`;
- retain the explicit negative term `-1 / 192` and independently compare the
  full double sum with the directly sampled tangent determinant.

Decision rule:

- if Lean confirms the obstruction, retire termwise factor-orientation
  composition as the route to the positive-series envelope and target a
  cancellation-preserving confluent/jet determinant theorem next;
- if the exact witness fails to formalize, identify the discrepancy between
  the Python flags and the canonical Lean reindexing before changing claims;
- do not weaken the final determinant sign objective merely because the
  sufficient Phase CA interface is too strong.

Non-claims:

- this phase does not refute positivity of the complete distinguished tangent
  determinant or its normalized envelope;
- it does not prove the order-eighteen contraction, dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `phaseCD_negative_triple_term` proves that the flags `(0,1,5)` and
  `(0,4,5)` contribute exactly `-1/192` to the ordered triple
  Cauchy--Binet expansion;
- `phaseCD_no_factorwise_positive_orientation` turns that witness into a
  formal impossibility theorem for Phase CA's separate factor-orientation
  interface with external signature one;
- `phaseCD_complete_tangent_det` and its positivity corollary prove that the
  complete distinguished determinant is nevertheless
  `24929/27648 > 0`;
- exact `Fraction` enumeration confirms that two features have only the
  positive terms `1/2` and `1/36`, whereas three features have 12 nonzero
  terms, 10 positive and 2 negative, whose sum matches the direct
  determinant;
- `./scripts/verify.sh` passes a 3,582-job Lean build, the expanded public
  axiom audit, 111 Python tests, formatting, typing, and the dependency audit.

Next decision:

Retire separate factorwise orientation as the main route. The next target is
a cancellation-preserving theorem for the complete confluent/jet determinant:
either a planar-network/Lindstrom--Gessel--Viennot model whose path families
encode the whole double-flag sum, or an exterior-algebra/sign-reversing
involution that pairs the negative flags before estimating the determinant.

## Phase CC contract: upper-bidiagonal total nonnegativity

State: complete; merged as PR #83 with green CI.

Phase CB exposes the feature-major tangent middle matrix as a nonnegative
upper-bidiagonal network: every row has its diagonal unit edge and at most one
forward edge, and the forward edges are disjoint.  The exhaustive exact
checks through five cells find no negative ordered minor.  The strongest
remaining middle-factor statement is therefore not a collection of bounded
checks but all-orders total nonnegativity.

Frozen public Lean targets:

- isolate a reusable total-nonnegativity theorem for the relevant sparse
  upper-bidiagonal or ordered direct-sum structure;
- prove `interleavedTangentMiddle` totally nonnegative for arbitrary finite
  feature count and arbitrary nonnegative weights;
- remove the conditional interleaved-TN premise from the diagonal tangent
  middle-factor orientation certificate;
- keep the theorem stated with explicit tuple orientations, since ordinary
  total nonnegativity remains false in layer-major order.

Independent check:

- retain the exhaustive exact minor census through five cells;
- add a deterministic theorem-facing regression that exercises zero and
  positive weights and compares the unconditional oriented certificate with
  the feature-major minors;
- retain the negative ordinary layer-major minor as a guard against erasing
  the shuffle sign.

Decision rule:

- prefer an order-theoretic determinant proof showing that every nonzero
  matching is order-preserving and hence unique;
- if direct determinant expansion creates disproportionate proof overhead,
  factor the matrix into adjacent nonnegative Jacobi cells and prove closure
  under ordered Cauchy--Binet;
- if either route reveals a false generic abstraction, minimize the exact
  obstruction and specialize only as far as the tangent-cell matrix requires.

Non-claims:

- this phase does not prove either outer confluent/generalized Vandermonde
  flag-minor family;
- it does not prove the normalized tangent envelope, the Phase BS contraction,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `HasNoncrossingSupport` isolates the ordered planar-support condition that
  forbids two crossed nonzero matrix entries;
- `totallyNonnegative_of_nonnegative_of_noncrossing` proves that a
  nonnegative noncrossing kernel is totally nonnegative at every order by
  showing every nonidentity Leibniz term vanishes;
- `interleavedTangentMiddle_hasNoncrossingSupport` proves the feature-major
  tangent cells have that support, including at zero weights;
- `interleavedTangentMiddle_totallyNonnegative` closes the Phase CB
  all-orders hypothesis for arbitrary finite feature counts;
- `upperTangentFeatureMiddle_diagonal_oriented_minor_nonneg` supplies the
  exact Phase CA middle-factor shuffle orientation directly from
  nonnegative diagonal weights;
- exact arithmetic retains all 397,252 Phase CB minor checks with mixed zero
  and positive weights and separately verifies 225 crossing pairs and 12,403
  nonidentity matching terms in a three-cell matrix;
- `./scripts/verify.sh` passes a 3,581-job Lean build, the expanded public
  axiom audit, 110 Python tests, formatting, typing, and the dependency audit.

Next decision:

Prove the sampled left and right block-diagonal feature-factor signs as
shuffle-oriented confluent/generalized Vandermonde flag minors. With the
middle factor now unconditional, those are the only missing sign families in
the Phase CA three-factor expansion.

## Phase CB contract: interleaved tangent cells

State: complete; merged as PR #82 with green CI.

Phase CA separates the tangent determinant into left, middle, and right flag
minors. For the diagonal feature weights used by the exponential family, the
middle matrix is a disjoint collection of `2 x 2` upper-unitriangular cells;
its negative layer-major minors arise only because those cells are not
contiguous in the layer ordering. Interleaving the two copies of each feature
should turn the middle problem into ordinary total nonnegativity of a sparse
upper-bidiagonal matrix, with the lost sign recovered by tuple orientation.

Frozen public Lean targets:

- define the layer/feature equivalence and the local weighted tangent cell;
- prove that the diagonal-weight upper-tangent middle matrix reindexes to a
  block diagonal family of those cells in feature-major order;
- prove every local cell minor is nonnegative and compute the determinant of
  the complete interleaved middle matrix;
- add the nonstrict analogue of arbitrary-tuple orientation transfer for a
  totally nonnegative kernel;
- derive a middle-factor sign certificate from total nonnegativity of the
  interleaved cell matrix, with the exact shuffle orientation exposed in the
  theorem statement.

Independent check:

- exhaustively enumerate all square minors for deterministic diagonal weights
  and feature counts through five;
- verify that feature-major minors are nonnegative and layer-major signs agree
  with the row/column sorting orientations;
- retain a negative layer-major order-two witness so no ordinary-TN claim can
  enter through the original ordering.

Decision rule:

- if the interleaved matrix is uniformly TN, make its sparse upper-bidiagonal
  proof the next exact theorem (or close it in this phase if the Lean proof is
  short);
- if feature-major TN fails, minimize the first exact cell-level obstruction
  before changing representations;
- do not assume a layer-major minor sign without its explicit orientation.

Non-claims:

- this phase need not prove all-orders TN of the interleaved middle matrix;
- it does not prove the outer confluent Vandermonde flag signs, the normalized
  tangent envelope, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `oriented_minor_nonneg_of_totallyNonnegativeAtOrder` is the nonstrict
  arbitrary-tuple counterpart of the existing strict-sign transfer theorem;
- `layerFeatureEquiv` and `featureInterleaveEquiv` make the layer-major to
  feature-major shuffle explicit;
- `tangentCell` packages `[[1,w],[0,1]]`, and
  `tangentCell_totallyNonnegative` proves every ordered minor of one cell is
  nonnegative for `w >= 0`;
- `reindex_upperTangent_diagonal_eq_blockDiagonal` proves that a diagonal
  tangent middle factor is exactly a block diagonal family of those cells;
- `interleavedTangentMiddle_det` proves the complete feature-major determinant
  is one;
- `upperTangentFeatureMiddle_oriented_minor_nonneg` transfers any
  interleaved-TN proof into the precise tuple-orientation certificate required
  by Phase CA's middle factor;
- exact rational arithmetic exhaustively verifies 198,626 feature-major
  minors and 198,626 corresponding oriented layer-major minors for one through
  five cells, and retains a negative unoriented layer-major witness;
- `./scripts/verify.sh` passes a 3,580-job Lean build, the expanded public
  axiom audit, 109 Python tests, formatting, typing, and the dependency audit.

The experiment supports, but does not replace, the remaining all-orders TN
theorem for `interleavedTangentMiddle`. Its entry formula is a nonnegative
upper-bidiagonal matrix with disjoint superdiagonal weights. A proof by the
Leibniz expansion should show that every nonzero matching is order-preserving,
so every ordered minor has at most one nonzero permutation term.

Next decision:

Prove all-orders total nonnegativity of a generic nonnegative upper-bidiagonal
matrix, then specialize it to `interleavedTangentMiddle`. This closes the
middle orientation hypothesis of Phase CA outright and leaves only the two
outer confluent/generalized Vandermonde flag-minor families.

## Phase CA contract: shuffle-aware tangent Cauchy--Binet

State: complete; merged as PR #81 with green CI.

Phase BZ identifies the required companion determinant with a distinguished
minor of a three-factor upper-tangent lift.  The remaining sign problem is
therefore not an ordinary total-nonnegativity claim: it is a flag-minor claim
whose left, middle, and right factors carry separate shuffle orientations.
The next reusable step is to expose all three minors in one ordered
Cauchy--Binet expansion and make their orientation certificates composable.

Frozen public Lean targets:

- define the ordered three-factor Cauchy--Binet summand for rectangular real
  matrices `A * M * B`;
- prove the exact double-sum expansion of an externally oriented determinant;
- prove a nonnegative-summand criterion and a modular criterion that combines
  separate left, middle, and right orientation certificates;
- specialize the expansion to finite upper-tangent feature factorizations,
  including arbitrary distinguished row and column samples;
- expose a single theorem whose remaining hypotheses are precisely the three
  families of oriented flag-minor inequalities needed for the Phase BZ
  envelope.

Independent check:

- add deterministic exact rational checks of the three-factor expansion under
  nontrivial sample shuffles and rectangular feature dimensions;
- retain the Phase BZ negative ordinary minor as a regression guard;
- compare the summand signs with the recorded selected/fresh cross-inversion
  parity rather than inferring positivity from determinant values alone.

Decision rule:

- if the upper-tangent specialization isolates uniform factorwise signs,
  pursue those finite exponential flag-minor inequalities next;
- if a factor family still changes sign after the explicit shuffle
  orientation, minimize that factor-level obstruction and move to a
  cancellation-preserving network or compound-matrix formulation;
- do not assume the desired determinant sign as a hypothesis.

Non-claims:

- this phase does not prove the finite exponential flag-minor inequalities;
- it does not prove the nine normalized tangent-minor bounds, residual half
  contraction, the dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `signedOrderedTripleCauchyBinetTerm` exposes the left, middle, and right
  minors of one ordered three-factor summand;
- `signature_mul_det_mul_mul_rect_eq_sum` proves the exact double-flag
  expansion for arbitrary rectangular real factors;
- `signature_mul_det_mul_mul_rect_nonneg_of_factor_orientations` composes
  independent factor orientations without assuming the final determinant
  sign;
- `upperTangentSampleLeft`, `upperTangentFeatureMiddle`, and
  `upperTangentSampleRight` canonically reindex the two feature layers;
- `signature_mul_det_submatrix_upperTangentMatrix_eq_sum` specializes the
  expansion to arbitrary sampled upper-tangent feature factorizations;
- `signature_mul_det_submatrix_upperTangentMatrix_nonneg_of_flag_orientations`
  states the missing sign theorem using exactly three factor-minor families;
- exact rational arithmetic verifies 24 randomized rectangular
  factorizations and 24 shuffled three-factor expansions, while the Phase BZ
  negative ordinary-minor and 360 oriented-coefficient checks remain active;
- `./scripts/verify.sh` passes a 3,579-job Lean build, the expanded public
  axiom audit, 108 Python tests, formatting, typing, and the dependency audit.

Next decision:

Classify minors of the canonically reindexed middle matrix
`[[1,W],[0,1]]` for nonnegative diagonal `W`.  Interleaving the two copies of
each feature turns it into contiguous `2 x 2` unitriangular blocks; the
conversion back to layer-major order should supply the middle shuffle sign.
The remaining left and right minors are then confluent/generalized
Vandermonde flag minors, for which the existing exponential total-positivity
library is the natural source of the two outer orientation certificates.

## Phase BZ contract: upper-tangent companion lift

State: complete; merged as PR #80 with green CI.

Phase BY makes every recursive companion determinant a normalized Phase BX
double border. The most useful next representation comes from tangent and
dual-number linear algebra: the double-bordered matrix is a distinguished
minor of the two-layer kernel

`upperTangent(K,H) = [[K,H],[0,K]]`.

When `K = A*B` and `H = A*W*B`, this lift factors through the same left and
right feature systems with the upper-unitriangular middle block
`[[1,W],[0,1]]`. Thus the missing envelope becomes a structured layered-minor
problem closely related to confluent total positivity, planar networks, and
compound matrices. This keeps all cancellations inside one determinant while
exposing the positive diagonal feature weights of the amplified exponential
family.

Frozen public Lean targets:

- define the generic rectangular upper-tangent matrix `[[K,H],[0,K]]`;
- define the distinguished row and column embeddings selected by a GECP run
  and arbitrary fresh samples;
- prove that the induced submatrix is exactly Phase BX's double-bordered
  companion matrix;
- combine that equality with Phase BY to identify its determinant with
  `det(selectedCore)^2` times the sampled recursive `conditionBoth`
  determinant;
- prove the generic feature-space factorization of the upper-tangent matrix
  when `K=A*B` and `H=A*W*B`;
- expose the positive-series amplified companion as
  `lambda • positiveSeriesTailMatrix + positiveSeriesLowLeft *
  positiveSeriesLowRight`.

Independent check:

- verify the layered-minor and feature-factorization identities in exact
  rational arithmetic across deterministic randomized dimensions and weights;
- enumerate the distinguished minors and their polynomial coefficients along
  complete-pivot prefixes, testing whether a fixed shuffle orientation makes
  them nonnegative;
- if that sign claim fails, minimize and retain the first exact obstruction,
  separating failure of ordinary total positivity from the narrower
  distinguished-minor class actually needed by Phase BS.

Decision rule:

- if distinguished coefficients retain one orientation, make matrix-valued or
  confluent total positivity the next theorem contract;
- if they do not, use the minimized obstruction to select an absolute-norm,
  complementary-minor, or mixed-discriminant bound instead;
- do not infer a determinant envelope from finite sign evidence.

Non-claims:

- this phase does not itself prove the nine normalized double-border bounds;
- it does not prove residual half contraction, the dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `Matrix.upperTangentMatrix` defines the rectangular layered matrix
  `[[K,H],[0,K]]`, with explicit distinguished double-border row and column
  embeddings;
- `Matrix.submatrix_upperTangentMatrix_eq_doubleBorderedCompanionMatrix` and
  `Run.upperTangentMatrix_submatrix_eq_doubleBorderedCompanionMatrix` prove
  that the selected/fresh layered minor is definitionally the Phase BX double
  border;
- `Matrix.upperTangentMatrix_mul_factorization` proves
  `upperTangent(A*B,A*W*B) = diag(A,A) * [[1,W],[0,1]] * diag(B,B)`, exposing
  the common feature system and the sole companion weight block;
- `Run.upperTangentMatrix_distinguishedMinor_det` identifies the distinguished
  determinant with `det(selectedCore)^2` times the recursive companion sample;
- `positiveSeriesAmplifiedFermionicKernel_eq_smul_tail_add_low` and
  `conditionBoth_positiveSeriesAmplifiedFermionicKernel` show that the
  positive-series amplified kernel is exactly the two-sided Phase BS pencil;
- `positiveSeriesAmplified_upperTangentMinor_det` is its specialized tangent
  identity, and
  `positiveSeriesFermionicKernel_det_le_halfPow_of_upperTangent` proves that
  normalized `2^18` bounds for the nine tangent minors already imply the
  required order-eighteen half-power determinant estimate;
- exact rational arithmetic verifies 24 randomized feature factorizations and
  360 complete-prefix determinant polynomials. Every polynomial coefficient
  becomes nonnegative after multiplying by the explicit shuffle sign
  `(-1)^(row cross-inversions + column cross-inversions)`; at the same time,
  ordinary order-two minors of the full lift remain negative.
- `./scripts/verify.sh` passes a 3,577-job Lean build, the expanded public
  axiom audit, 108 Python tests, formatting, typing, and the dependency audit.

The coefficient observation is evidence, not a theorem. It is materially
narrower than total nonnegativity and survives precisely where the ordinary
minor property fails.

Next decision:

Formalize the cross-inversion-oriented coefficient statement as a flag-minor
or matrix-valued total-positivity theorem for the exponential feature lift.
The feature factorization suggests two disciplined routes: a shuffle-aware
Cauchy--Binet expansion of `diag(A,A) * [[1,W],[0,1]] * diag(B,B)`, or a
two-layer planar network whose path-family crossing sign is exactly the
recorded parity. Either theorem should target only the distinguished
selected/fresh minors, since ordinary total nonnegativity is exactly false.

## Phase BY contract: selected-core companion closed form

State: complete; merged as PR #79 with green CI.

Phase BX proves the correct generic double-border determinant formula, and
exact arithmetic confirms that it matches a three-pivot companion run. The
remaining algebraic bridge is to identify the recursively defined
`Run.conditionBoth` with that selected-core Schur expression for every run.
This can be proved without expanding Plücker coordinates: row conditioning is
the unique correction by selected kernel columns that vanishes on selected
rows, the column statement is dual, and left/right elementary projectors
commute.

Frozen public Lean targets:

- define selected kernel and companion row/column block matrices for a run;
- prove closed selected-core inverse formulas for `conditionRows` and
  `conditionColumns`;
- prove fixed left and right companion updates commute and deduce
  `conditionBoth = conditionColumns ∘ conditionRows`;
- derive the four-term selected-core formula for `conditionBoth`;
- specialize Phase BX's squared-core determinant identity to arbitrary finite
  samples of a successful run.

Independent check:

- retain the exact Phase BX three-pivot entrywise and determinant checks;
- add deterministic exact randomized runs comparing the recursive and closed
  formulas for multiple prefix lengths, companion weights, and fresh blocks;
- use solve-based rational arithmetic and never form an explicit inverse.

Decision rule:

- construct the selected-core inverse instance only from the already proved
  `selectedCore_det_ne` theorem;
- use annihilation and interpolation uniqueness instead of recursively
  expanding a block inverse when possible;
- stop at a precise missing lemma if dependent selected-index bookkeeping
  prevents the full sampled specialization.

Non-claims:

- this phase does not prove the normalized double-border envelope;
- it does not prove residual half contraction, the dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `Run.conditionRows_eq_selectedCoreSchur` and
  `Run.conditionColumns_eq_selectedCoreSchur` identify both iterated
  one-sided operators with their selected-core oblique-projector formulas;
- `companionRowUpdate_companionColumnUpdate` and
  `Run.conditionBoth_eq_conditionColumns_conditionRows` prove that the
  recursively interleaved left/right operations may be globally composed;
- `Run.conditionBoth_eq_twoSidedCompanionSchur` proves the exact four-term
  formula for every successful run over an arbitrary field;
- `Run.conditionBoth_sample_eq_twoSidedCompanionSchur` transports the formula
  to arbitrary finite row and column samples, while
  `Run.doubleBorderedCompanionMatrix_det_eq_conditionBoth_sample` specializes
  Phase BX to show that the double-border determinant is exactly
  `det(selectedCore)^2` times the sampled recursive companion determinant;
- the selected-core inverse is constructed solely from
  `selectedCore_det_ne`; no extra invertibility assumption is exposed;
- exact solve-based rational arithmetic checks 144 deterministic randomized
  combinations of Cauchy-kernel runs, prefix lengths one through three, four
  companion weights, and varying fresh blocks, in addition to retaining the
  Phase BX witness checks.
- `./scripts/verify.sh` passes a 3,574-job Lean build, the expanded public
  axiom audit, 107 Python tests, formatting, typing, and the dependency audit.

Next decision:

Control the normalized double-bordered determinant for the positive-series
amplified companions at the nine Phase BS weights. The recursive/closed-form
ambiguity is gone: this is now the exact remaining analytic premise. Ordinary
total nonnegativity is unavailable because Phase BX's negative order-two
minors persist, so the next phase should test cancellation-preserving routes
such as complementary-minor condensation, planar-network realizations with
internal signed paths, or a mixed-discriminant inequality stable under the
selected-core Schur prefix.

## Phase BX contract: double-bordered companion determinant

State: complete; merged as PR #78 with green CI.

Phase BW linearizes a rank-update pencil, but Phase BS requires bounds on the
amplified conditioned tail plus the conditioned rank-eight block. For a fixed
GECP prefix, two-sided companion conditioning has a more direct global form.
If `C` is the selected core, the fresh companion block is the two-sided Schur
expression

`H - K_L C⁻¹ H_R - H_L C⁻¹ K_R + K_L C⁻¹ H_C C⁻¹ K_R`.

This is the Schur complement of one double-bordered matrix whose leading block
is upper triangular with two copies of `C`. Its determinant therefore carries
the factor `det(C)^2` seen termwise in Phase BU, but without expanding any
Plücker coordinates. This is the object that directly matches the two-sided
finite-weight companion envelope.

Frozen public Lean targets:

- define the generic two-sided companion Schur expression for an invertible
  finite core `C`;
- define its double-bordered matrix with leading block
  `[[C, H_C], [0, C]]` and the corresponding `K/H` cross blocks;
- prove that the double-bordered determinant equals
  `det(C)^2 * det(twoSidedCompanionSchur ...)`;
- keep the identity generic over finite real matrices and make every
  invertibility assumption explicit.

Independent check:

- for the exact Phase BV three-pivot witness, compare direct iterative
  `conditionBoth` against the closed two-sided Schur formula entrywise;
- verify the squared-core determinant identity at weight one and all nine
  Phase BS weights;
- inspect whether the double border removes or merely relocates the Phase BW
  signed-cycle obstruction.

Decision rule:

- use Mathlib's upper-block-triangular inverse and Schur determinant APIs;
- retain the complete determinant rather than expanding feature subsets;
- if the fixed-run equality needs a separate selected-core inverse theorem,
  record that theorem precisely as the next contract instead of hiding it in
  a numerical premise.

Non-claims:

- this phase does not yet identify every recursive `Run.conditionBoth` entry
  with the generic Schur expression in Lean;
- it does not prove the double-bordered determinant envelope, residual half
  contraction, the dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `twoSidedCompanionSchur` is the exact four-term projection formula for an
  arbitrary invertible finite core and separate kernel/companion cross blocks;
- `doubleBorderedCompanionMatrix` packages the core, companion core, four
  cross blocks, and fresh companion block into one matrix;
- `doubleBorderedCompanionMatrix_det` and
  `abs_doubleBorderedCompanionMatrix_det` prove the exact and absolute
  factorizations by `det(C)^2`;
- exact solve-based rational arithmetic shows that three iterative companion
  pivots equal the closed Schur expression entrywise at weight one and every
  Phase BS norming weight, and independently verifies the determinant
  factorization;
- negative order-two minors remain in every checked double border, so the new
  representation preserves the required cancellation but is not totally
  nonnegative;
- `./scripts/verify.sh` passes a 3,573-job Lean build, the public axiom audit,
  107 Python tests, formatting, typing, and the dependency audit.

Next decision:

Prove the selected-core closed form for `Run.conditionBoth` in Lean and
specialize the double-border factorization to the positive-series amplified
companion. That makes the nine Phase BS premises equivalent to normalized
double-bordered determinant bounds, with no ambiguity between the coefficient
surrogates of Phases BV/BW and the actual conditioned-tail envelope.

## Phase BW contract: bordered-pencil linearization

State: complete; merged as PR #77 with green CI.

Phase BV identifies the surviving signed Plücker sum with a coefficient of
`det(mu • (U*V) + A*B)`. Expanding that coefficient returns to the cancelling
coordinates. A standard linearization from mixed-discriminant and systems
theory instead embeds the entire pencil into one bordered determinant with an
identity corner. This keeps all cancellations internal and creates a single
matrix on which total-positivity, planar-network, or Schur-complement methods
can operate.

Frozen public Lean targets:

- define the bordered rank-update pencil with blocks `1`, `V`, `-mu • U`,
  and `B`;
- prove its determinant is exactly `det(mu • (U*V) + B)` without any
  invertibility hypothesis;
- compose the identity with Phase BV's combined fixed-plus-feature matrices
  and the existing block-feature polynomial;
- expose the result as a reusable linearization of every finite weight needed
  by Phase BS.

Independent check:

- verify the bordered determinant identity in exact rational arithmetic on
  the Phase BV three-pivot witness at weights `1, 56, 65, ..., 128`;
- preserve the exact 2-of-56 coordinate-sign obstruction while confirming
  that the unexpanded bordered determinant agrees at every weight;
- test whether a fixed row/column sign gauge makes the bordered matrices
  totally nonnegative; if not, record the obstruction and keep the next target
  at the aggregate determinant level.

Decision rule:

- use Mathlib's Schur-complement determinant identity with the identity block;
- do not assume `B`, `U`, `V`, or any selected core is invertible;
- prefer one determinant identity over another coordinate expansion.

Non-claims:

- this phase does not prove positivity or the required envelope for the
  bordered determinant;
- it does not prove the nine Phase BS companion envelopes, residual half
  contraction, the dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `borderedRankUpdatePencil` places `1`, `V`, `-mu • U`, and `B` in one
  square block matrix;
- `borderedRankUpdatePencil_det` proves its determinant is exactly
  `det(mu • (U*V) + B)` without invertibility assumptions;
- `mixedColumnPolynomial_eval_eq_borderedRankUpdatePencil` and
  `blockFeaturePolynomial_eval_eq_borderedRankUpdatePencil` connect the
  generic and Phase BV graded polynomials to that single determinant;
- exact rational arithmetic checks the identity at weight one and all nine
  Phase BS weights `56, 65, ..., 128` on the Phase BV witness;
- a negative four-cycle involving the identity, positive `V`, negative `U`,
  and positive low block is invariant under row/column sign changes, proving
  that no such gauge can make the bordered matrix entrywise nonnegative;
- `./scripts/verify.sh` passes a 3,572-job Lean build, the public axiom audit,
  107 Python tests, formatting, typing, and the dependency audit.

Next decision:

Exploit the bordered determinant without demanding total nonnegativity of its
entries. The signed four-cycle points toward a cancellation-preserving planar
network, a complementary-minor identity, or a mixed-discriminant inequality.
Any next abstraction should target the determinant or its Schur complement as
a whole and retain the exact finite-weight envelope required by Phase BS.

## Phase BV contract: graded Plücker pairing

State: complete; merged as PR #76 with green CI.

Phase BU converts conditioned matched terms into products of original
feature-bordered determinants. A stronger exact search then falsifies the
tempting termwise continuation: after three positive-unit finite-Taylor GECP
pivots, two of 56 order-three feature subsets have the opposite bordered-
product sign. The fixed-degree sums nevertheless retain a common orientation
in every search so far. The correct object is therefore the graded exterior-
algebra pairing of the two Plücker vectors, not each coordinate product.

Frozen public Lean targets:

- define the finite combined feature matrices obtained by adjoining a fixed
  coordinate block `U,V` to a feature factorization `A,B`;
- prove that marking the fixed block by `mu` produces exactly
  `mu • (U*V) + A*B`;
- identify the block-feature coefficient at the full fixed-coordinate degree
  with the corresponding physical mixed-column coefficient;
- expose the ordered Cauchy--Binet sum at that degree as the graded Plücker
  pairing which survives the termwise obstruction;
- record the exact finite-Taylor counterexample to coordinatewise sign
  coherence and verify that its full graded sum has the predicted sign.

Independent check:

- reproduce the three-pivot rational witness with selected rows
  `[6,0,3]`, selected columns `[6,0,2]`, and all 56 order-three subsets of
  eight low features;
- verify that exactly two individual bordered products reverse sign while the
  total fixed-degree pairing agrees with the determinant coefficient;
- continue deterministic randomized checks of the graded coefficient rather
  than treating the failed coordinatewise property as evidence.

Decision rule:

- preserve the full signed sum: do not take absolute values term by term;
- use existing block-feature Cauchy--Binet and mixed-column polynomial APIs
  instead of creating a parallel determinant expansion;
- keep the algebra generic over finite real matrices and isolate all
  fermionic sign claims as later analytic premises.

Non-claims:

- this phase does not prove positivity of the graded pairing for arbitrary
  fermionic prefixes;
- it does not prove the nine Phase BS companion envelopes, residual half
  contraction, the dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered result:

- `combinedFeatureLeft` and `combinedFeatureRight` adjoin a fixed block to an
  arbitrary rectangular feature factorization;
- `weightedCombinedFeatureProduct` identifies fixed-coordinate weighting with
  the physical pencil `mu • (U*V) + A*B`;
- `fixedBlockPluckerPairing_eq_mixedColumnCoefficient` and
  `fixedBlockPluckerPairing_eq_sum` identify its full-fixed-degree coefficient
  with the surviving signed ordered Cauchy--Binet sum;
- the exact three-pivot witness recovers selected rows `[6,0,3]`, selected
  columns `[6,0,2]`, and 2 negative versus 54 positive feature-coordinate
  products, while their total is positive and exactly equals the degree-three
  mixed-column coefficient;
- `./scripts/verify.sh` passes a 3,571-job Lean build, the public axiom audit,
  107 Python tests, formatting, typing, and the dependency audit.

Next decision:

Treat coefficient positivity as a single compound-matrix pairing problem.
The most promising mechanisms are a sign-reversing cancellation/involution on
the exceptional Plücker coordinates, or a stability/Lorentzian property of
the determinant pencil strong enough to control its full-fixed-degree
coefficient. Coordinatewise total positivity is ruled out by the exact
witness and should not be revisited without a genuinely stronger invariant.

## Phase BU contract: companion Plücker factorization

State: complete; merged as PR #75 with green CI.

Phase BT shows that the conditioned low block factors through separately
projected left and right features, and exact evidence suggests their matched
minors acquire opposite sign changes which cancel in each Cauchy--Binet term.
The natural exact invariant comes from Schur-complement and Grassmannian
geometry: a conditioned feature minor should be the corresponding original
feature-bordered determinant divided by the selected-core determinant.

Frozen public Lean targets:

- define a left feature-augmented core whose selected columns come from the
  pivoting kernel and whose finite border columns come from an arbitrary left
  companion factor;
- define the transposed right analogue with selected kernel rows and arbitrary
  right companion feature rows;
- prove both determinant factorizations through `run.selectedCore.det` and
  the sampled minors of `run.conditionRows` and `run.conditionColumns`;
- combine the two formulas into an exact matched-product identity with the
  square of the selected-core determinant;
- specialize the interface to the positive-series fermionic left and right
  rank-eight features.

Independent check:

- use exact `Fraction` arithmetic across multi-pivot finite positive-Taylor
  runs to compare both augmented determinants with the conditioned minors;
- verify the matched-product identity for every feature subset of each tested
  order, including the Phase BT mixed-sign witness.

Decision rule:

- define the feature-augmented cores recursively in pivot order so the proof
  uses elementary elimination and never an explicit inverse;
- keep row and column identities separate until their determinant signs and
  index orientations are explicit;
- treat the formulas as Plücker-coordinate identities only after exact
  determinant equalities are proved.

Non-claims:

- this phase does not prove the signs of the original bordered feature
  determinants or all-orders matched-minor coherence;
- it does not prove the nine Phase BS companion envelopes, residual half
  contraction, the dyadic recurrence, Conjecture G1, or Simons Problem 4.2.

Delivered identifiers:

- `GECP.Run.leftFeatureColumn`, `leftFeatureCore`, `rightFeatureRow`, and
  `rightFeatureCore`;
- `leftFeatureCore_step_det` and `rightFeatureCore_step_det`;
- `leftFeatureCore_det_eq_selectedCore_det_mul_conditionRows_minor`;
- `rightFeatureCore_det_eq_selectedCore_det_mul_conditionColumns_minor`;
- `leftFeatureCore_det_mul_rightFeatureCore_det`;
- `positiveSeriesLowLeft_featureCore_det` and
  `positiveSeriesLowRight_featureCore_det`;
- `positiveSeriesConditionedLow_core_sq_mul_signedTerm`;
- `positiveSeriesConditionedLow_signedTerm_nonneg_iff`.

Result:

- every left or right conditioned feature minor is exactly its original
  feature-augmented determinant divided algebraically by the nonzero selected
  core, with the proof performed by elementary elimination rather than an
  inverse;
- multiplying matched left/right identities makes the selected-core factor a
  square, so its pivot-order sign disappears;
- for the positive-series rank-eight features, nonnegativity of every
  conditioned signed Cauchy--Binet term is equivalent to nonnegativity of the
  corresponding product of two original bordered determinants;
- the remaining sign problem is therefore unconditioned. Each bordered
  family mixes exponential columns or rows with monomials, and the monomials
  are confluent derivatives of the exponential family at parameter zero.

Independent checks:

- exact `Fraction` arithmetic verifies 2,382 left, right, and matched-product
  determinant identities over the Phase BT witness and deterministic
  randomized multi-pivot finite positive-Taylor runs;
- every tested feature subset is enumerated at each tested minor order, and
  the earlier 2,760-sample two-sided sign regression remains intact.

Verification:

- every new public theorem reports only the accepted Lean foundations in the
  axiom audit, with no `sorry`, `admit`, or custom axiom;
- `./scripts/verify.sh` passes a 3,570-job Lean build, 106 Python tests, and
  the dependency audit with no known vulnerabilities; the change passed CI
  and merged as PR #75.

## Phase BT contract: two-sided companion projection

State: complete; merged as PR #74 with green CI.

Phase BS leaves nine determinant envelopes for the one-sided conditioned
companion. Exact searches show that those envelopes are plausible, but the
one-sided companion loses coefficient signs and strict sign regularity. The
natural construction from oblique projection theory is two-sided: apply the
fixed pivot's left row projector and right column projector to the companion
matrix. This map is linear in the companion, recovers the ordinary Schur
update on the pivoting matrix, annihilates the selected row and column, and
maps `L*R` to a product of separately conditioned feature factors.

Frozen public Lean targets:

- define fixed-pivot left and right companion feature updates and the
  two-sided companion update;
- prove the expanded four-term formula and show that self-application is the
  ordinary residual update;
- prove linearity, scalar compatibility, selected-row/column annihilation,
  and exact factor transport
  `biUpdate(K, L*R) = leftUpdate(K,L) * rightUpdate(K,R)`;
- iterate the construction through an arbitrary fixed `GECP.Run`, proving
  self-conditioning and transport of an initial `K = E + L*R` split with the
  same feature dimension on both sides;
- expose the finite sampled factorization needed by the mixed-column and
  block-feature determinant APIs.

Independent check:

- use exact `Fraction` arithmetic on multi-pivot finite positive-Taylor
  matrices to verify the two-sided split, annihilation, and feature product;
- exhaustively check that all nonzero mixed coefficients share one sign over
  deterministic randomized prefixes and fresh samples, while retaining a
  fixed witness where the one-sided coefficients have mixed signs.

Decision rule:

- keep pivot data fixed by the original run and preserve linearity in the
  companion;
- represent the two-sided map through separately conditioned feature factors,
  not an explicit inverse;
- expose matched left/right feature-minor sign coherence as the next analytic
  theorem; do not promote the exact search to a proof.

Non-claims:

- this phase does not prove all-orders matched-minor sign coherence or the
  nine Phase BS determinant envelopes;
- it does not yet prove residual half contraction, the dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered identifiers:

- `GECP.companionColumnUpdate` and `GECP.companionBiUpdate`;
- `companionBiUpdate_apply`, `companionBiUpdate_self`,
  `companionBiUpdate_selectedRow`, and `companionBiUpdate_selectedColumn`;
- `companionBiUpdate_add`, `companionBiUpdate_smul`, and
  `companionBiUpdate_mul`;
- `GECP.Run.conditionColumns` and `GECP.Run.conditionBoth` with their
  self-conditioning, linearity, and factor-transport theorems;
- `GECP.Run.finalResidual_eq_conditionBoth_add_mul` and its sampled form;
- `positiveSeriesFermionicKernel_finalResidual_biFeatureSplit` and its
  sampled form;
- `positiveSeriesConditionedLow_signature_mul_det_eq_sum`.

Result:

- one fixed pivot acts on a companion as the explicit four-term two-sided
  oblique projection and annihilates both its selected row and selected
  column;
- self-application is exactly the ordinary Schur residual update, while a
  rank-`k` product is transported as separately conditioned left and right
  factors, so an arbitrary prefix preserves the feature dimension on both
  sides;
- the positive-series fermionic residual therefore has a canonical
  two-sided rank-eight split, and its conditioned low determinant is exposed
  as a sum of products of matched left/right feature minors;
- exact arithmetic retains the one-sided witness with coefficient signs
  `0, 0, +, -, +`, whereas its two-sided counterpart has signs
  `0, 0, +, +, +`. Across 2,760 exhaustively enumerated fresh samples from
  deterministic randomized prefixes, every set of nonzero two-sided mixed
  coefficients has one common orientation. This is evidence, not the
  all-orders sign theorem.

Verification:

- all new public results report only the accepted Lean foundations in the
  axiom audit, with no `sorry`, `admit`, or custom axiom;
- the exact regression verifies the residual split, selected row and column
  zeros, and separately transported feature product after every pivot;
- `./scripts/verify.sh` passes a 3,568-job Lean build, 106 Python tests, and
  the dependency audit with no known vulnerabilities; the change passed CI
  and merged as PR #74.

## Phase BS contract: finite-weight conditioned-pencil amplification

State: complete; merged as PR #73 with green CI.

Phase BR reduces the restart problem to a rank-eight determinant pencil, but
its exact regression shows that conditioning can make individual mixed-column
coefficients negative. Further exact searches rule out two tempting repairs:
the conditioned pencil need not be real-rooted, and amplification at the
single weight `128` need not dominate its value at one. The degree after
factoring the forced rank-eight zero is nevertheless at most eight. Polynomial
norming sets therefore offer a sign-free replacement: a degree-eight
polynomial is controlled by its values at nine fixed weights.

Frozen public Lean targets:

- define the nine integer weights `56, 65, ..., 128` and their exact Lagrange
  coefficients for evaluation at one;
- prove the degree-eight interpolation identity and the explicit coefficient-
  mass bound by `518173`;
- combine low-rank coefficient vanishing with that identity to prove a
  sign-free determinant amplification theorem from uniform bounds at those
  nine weights;
- specialize to rank eight and fresh order eighteen, where exact arithmetic
  turns a `2^18` envelope at every norming weight into a `2^-18` bound for the
  physical residual determinant;
- compose the result with Phase BR's canonical conditioned fermionic split,
  leaving only the finite companion-envelope premise explicit.

Independent check:

- use exact `Fraction` arithmetic to verify every Lagrange moment through
  degree eight, the constant bound, and the order-eighteen half-power
  inequality;
- search exact finite positive-Taylor companion pencils for failures of
  coefficient positivity, real-rootedness, and single-weight amplification,
  while checking the nine-weight certificate on every case.

Decision rule:

- use a finite positive norming set, not an unproved sign, root-location, or
  stability hypothesis;
- keep all weights at most `128`, so every future analytic envelope can reuse
  the translated-band entry bound `2`;
- expose the nine companion determinant envelopes as the next analytic
  obligation rather than hiding them in a new certificate predicate.

Non-claims:

- this phase does not prove the companion determinant envelopes themselves;
- it does not yet prove residual half contraction, the dyadic recurrence,
  Conjecture G1, or Simons Problem 4.2.

Delivered identifiers:

- `rankEightAmplificationNode` and `rankEightLagrangeWeight` encode the exact
  nine-node norming set `56, 65, ..., 128` and its evaluation-at-one weights;
- `rankEightLagrange_moment`, `rankEightLagrange_eval_one`, and
  `rankEightLagrange_abs_sum_le` prove the degree-eight interpolation identity
  and the explicit coefficient-mass bound `518173`;
- `rankEightLagrange_norming_bound` removes every coefficient-sign and
  root-location assumption from the polynomial comparison;
- `mixedColumnPolynomial_rankEight_finiteWeight_bound` combines that norming
  theorem with rank-eight coefficient vanishing;
- `rankEight_finiteWeight_det_le_halfPow` proves that nine `2^18` envelopes at
  fresh order eighteen imply the exact physical determinant bound `2^-18`;
- `positiveSeriesFermionicKernel_conditioned_det_le_halfPow` applies the
  result directly to Phase BR's canonical conditioned fermionic split.

Result:

- conditioned coefficient positivity, alternating signs, real-rootedness,
  and domination at one preferred weight are no longer proof obligations;
- the entire determinant-pencil question is reduced to nine fixed companion
  matrices, all with weights at most `128` and therefore compatible with the
  existing translated-band entry envelope;
- the order-eighteen constants have ample exact slack:
  `518173 * 4^18 < 56^10`;
- the remaining analytic target is to prove a `2^18` determinant envelope for
  each of those nine companion residual samples. This is strictly weaker and
  more stable than coefficientwise positivity and is supported by exhaustive
  exact finite-Taylor minor checks performed during the research loop.

Verification:

- every new public theorem reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the public axiom audit;
- exact `Fraction` arithmetic verifies all nine interpolation moments, the
  constant and half-power inequalities, a genuine failure of fixed-weight
  `+128` amplification, a negative-discriminant conditioned pencil, and the
  surviving nine-weight certificate;
- `./scripts/verify.sh` passes with 105 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #73.

## Phase BR contract: Schur-aware feature elimination

State: complete; merged as PR #72 with green CI.

Phase BQ leaves one exact premise: the final residual sample must be written as
a conditioned tail plus a uniformly low-rank block.  The Schur update is
nonlinear in the kernel that selects the pivot, but it is linear in a companion
kernel once the pivot run is fixed.  Applying those fixed row operations to an
initial identity `K = E + L*R` therefore transports `E` and `L` while leaving
the right feature factor `R` unchanged.  This is the feature-level analogue of
Gaussian elimination with multiple right-hand sides.

Frozen public Lean targets:

- define the one-step and full-run row-conditioning operators attached to an
  exact `GECP.Run`;
- prove that self-conditioning recovers `run.finalResidual`;
- prove additivity and scalar compatibility of the conditioning operator;
- prove that conditioning commutes with multiplication by a fixed right
  feature matrix;
- conclude that every initial split `K = E + L*R` yields the canonical final
  split `run.finalResidual = condition(E) + condition(L)*R`, so the low-rank
  feature dimension does not grow with prefix length;
- sample that identity on arbitrary fresh rows and columns and compose it with
  Phase BQ's residual amplification theorem under the remaining mixed-
  coefficient sign premise.

Independent check:

- use exact `Fraction` arithmetic on a nontrivial multi-pivot matrix to compare
  recursive GECP residuals with companion row conditioning, verify the
  transported feature split entrywise, and check that the conditioned low
  block retains its original feature rank bound.

Decision rule:

- keep the elimination coefficients fixed by the original run; do not define
  a new run for the tail or amplified companion kernel;
- preserve the right feature factor exactly, as in elimination with multiple
  right-hand sides;
- expose mixed-coefficient nonnegativity as the next analytic obligation
  rather than assuming it through a renamed certificate.

Non-claims:

- this phase does not prove that the canonical conditioned fermionic tail has
  nonnegative mixed coefficients;
- it does not prove a uniform amplified determinant envelope for the
  conditioned tail, the dyadic recurrence, Conjecture G1, or Simons Problem
  4.2.

Delivered identifiers:

- `GECP.companionRowUpdate` and `GECP.Run.conditionRows` apply the pivot
  coefficients selected by one fixed run to an arbitrary companion matrix;
- `companionRowUpdate_self` and `conditionRows_self` recover the ordinary
  residual update and final residual;
- `conditionRows_add`, `conditionRows_smul`, and `conditionRows_mul` establish
  linearity and exact commutation with a fixed right feature factor;
- `finalResidual_eq_conditionRows_add_mul` and its sampled form transport any
  initial identity `K = E + L*R` through an arbitrary selected prefix without
  enlarging the feature dimension;
- `positiveSeriesFermionicKernel_finalResidual_featureSplit` gives the
  canonical rank-eight split for the positive-series fermionic kernel;
- `conditionRows_rank_amplification` and
  `positiveSeriesFermionicKernel_conditioned_rank_amplification` compose that
  split with Phase BQ once the conditioned-tail coefficient signs are known.

Result:

- the algebraic compatibility premise of Phase BQ is now constructed rather
  than assumed: every fermionic prefix has an exact conditioned-tail plus
  rank-eight feature split;
- the right feature factor is literally unchanged through the full prefix,
  so earlier pivots do not consume the fresh determinant exponent;
- the remaining sign premise is not automatic. In the exact finite
  positive-Taylor regression, the first conditioned mixed coefficients have
  signs `0, 0, +, -, +`; in particular the degree-three coefficient is
  `-846159 / 4172800000000000`;
- this rules out a direct reuse of the unconditioned coefficientwise-
  nonnegative pencil. The promising replacement is a signed, alternating, or
  real-rooted determinant-pencil theorem strong enough to compare the
  physical and amplified values at the two positive weights actually used.

Verification:

- every new public theorem reports only the accepted Lean foundations in the
  public axiom audit;
- an exact `Fraction` regression follows three GECP pivots, checks the
  companion split entrywise after every step, and verifies every `3 x 3`
  minor of the conditioned rank-two low block vanishes;
- `./scripts/verify.sh` passes with 104 Python tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #72.

## Phase BQ contract: residual amplification bridge

State: complete; merged as PR #71 with green CI.

Phase BP supplies a center-uniform local determinant contraction for the
original kernel.  Restarting it after earlier pivots is not an exceptional-
column problem: changing a fresh column also changes its entries on previously
selected rows, so the selected Schur data are no longer fixed.  Exact finite
positive-Taylor experiments exhibit the resulting sign failure.  The correct
object is the final residual itself.  Its finite sample is already the Schur
complement encoded by the run, and Phase AJ identifies every augmented
original-kernel determinant with the selected-core determinant times that
residual minor.

Frozen public Lean targets:

- give an exact rational finite positive-Taylor witness showing that keeping
  one old column exact while tail-amplifying the remaining columns can reverse
  the amplified determinant sign and violate the proposed rank amplification;
- prove that a sampled final residual admitting an exact decomposition
  `E + L*R`, nonnegative mixed coefficients, and rank `k` satisfies the same
  `lambda^(n-k)` determinant amplification as an unconditioned kernel;
- prove nonnegativity of the original and amplified residual determinants from
  the coefficient and rank hypotheses, so the result is available in absolute
  value without an orientation assumption;
- compose the residual theorem with the arbitrary-border determinant identity,
  exposing the selected-core determinant only as a multiplicative factor and
  keeping the decay exponent independent of the prefix length;
- specialize the bridge to fermionic runs while leaving the analytic residual
  decomposition as an explicit premise.

Delivered identifiers:

- `Matrix.mixedColumnPolynomial_rank_amplification_abs`;
- `GECP.Run.finalResidual_rank_amplification`;
- `GECP.Run.augmentedCore_rank_amplification`;
- `frozenColumnLow_eq_mul`;
- `frozenColumnTaylor_original_det`;
- `frozenColumnTaylor_amplified_det`;
- `frozenColumnTaylor_rankAmplification_fails`;
- `fermionicKernel_finalResidual_rank_amplification`;
- `fermionicKernel_augmentedCore_rank_amplification`.

Result:

- an exact six-term positive-Taylor sample has original determinant
  `91063777 / 85899345920 > 0`, but amplifying the three fresh-column tails by
  `128` gives determinant `-32176499 / 983040`; its frozen background factors
  through two features, so one-sided frozen-column rank amplification is
  false even in this finite positive-feature setting;
- for an exact final-residual split `E + L*R` with nonnegative mixed
  coefficients and residual low rank `k < n`, weight `lambda >= 1` gives
  `lambda^(n-k) * |det(E+L*R)| <= |det(lambda E+L*R)|`;
- composing with the Phase AJ augmented-core identity leaves the selected
  prefix only in the multiplicative factor `|det selectedCore|`; it does not
  reduce the fresh-sample amplification exponent.

Independent checks:

- recompute the obstruction entirely with `Fraction` arithmetic, including
  the original positive determinant, the negative weight-`128` determinant,
  and the failed comparison;
- verify the augmented/residual determinant comparison on exact rational
  geometric-surrogate runs with nontrivial prefixes and low-rank splits.

Verification:

- the rational obstruction is kernel-checked from explicit LU
  factorizations; the public axiom audit contains only the accepted Lean
  foundations and no evaluator or sorry axiom;
- the exact `Fraction` regression independently checks the obstruction, a
  compatible one-pivot residual split, coefficient signs, rank amplification,
  and the augmented/residual determinant identity;
- `./scripts/verify.sh` passes all 3,563 Lean jobs and 103 Python tests with no
  known dependency vulnerabilities; the change passed CI and merged as PR
  #71.

Decision rule:

- formulate the bridge at the residual sample level rather than masking rows
  or columns of the original kernel;
- reuse `mixedColumnPolynomial_rank_amplification` and Phase AJ's exact
  augmented determinant factorization; do not introduce an inverse;
- if absolute-value transport requires extra orientation plumbing, first prove
  determinant nonnegativity directly from nonnegative coefficients and strict
  rank deficiency.

Non-claims:

- the algebraic bridge does not assert that a fermionic final residual already
  has the required low-rank-plus-compatible-tail decomposition;
- this phase does not prove the shrinking-band recurrence, Conjecture G1, or
  Simons Problem 4.2.

## Phase BI contract: mixed-column coefficient bridge

State: complete; merged as PR #63 with green CI.

Phase BI connects the feature-index coefficient from Phase BH to the
coordinate-space mixed determinants needed for the exact tail limit. The
low/high feature split produces a matrix pencil `B + lambda E`; determinant
multilinearity produces its coefficients as sums over physical column choices,
while weighted Cauchy--Binet produces the same polynomial by feature count.

Frozen public Lean targets:

- define low- and high-block feature matrices and prove their weighted
  recombination identity;
- define the fixed-cardinality sum of physical mixed-column determinants;
- prove, by equality of determinant polynomials, that this physical
  coefficient equals the grouped feature coefficient;
- conclude nonnegativity of every physical mixed-column coefficient from
  termwise nonnegative ordered Cauchy--Binet summands;
- specialize the result to finite positive exponential Taylor features.

Independent check:

- extend the exact rational regression to compare each feature-count
  coefficient with the corresponding fixed-cardinality physical-column sum.

Decision rule:

- use polynomial extensionality over `Real` rather than differentiating the
  determinant pencil; this exposes the algebraic coefficient identity directly;
- if the generic polynomial bridge becomes dominated by coercion plumbing,
  retain the theorem over real matrices rather than prematurely generalizing
  the coefficient ring.

Non-claims:

- this phase remains finite-dimensional in the feature cutoff;
- it does not yet pass the positive Taylor tail to infinity or prove the
  quantitative determinant, dyadic, Conjecture G1, or Problem 4.2 bounds.

Delivered:

- `Matrix.mixedColumnPolynomial` packages determinant multilinearity by the
  number of selected high-block columns;
- `Matrix.blockFeaturePolynomial_eq_mixedColumnPolynomial` proves equality of
  the feature-selection and physical-column generating polynomials by
  evaluation at every real weight;
- `Matrix.blockFeatureCoefficient_eq_mixedColumnCoefficient` extracts the
  coefficient identity by polynomial uniqueness;
- `expTaylor_mixedColumnCoefficient_nonneg` proves every fixed-cardinality sum
  of finite positive-Taylor mixed-column determinants nonnegative.

Result:

Phase BH's abstract feature-count coefficient is exactly the coordinate-space
quantity needed for a tail limit. This matters because the latter is a finite
sum of ordinary determinants, so entrywise convergence of finite Taylor tails
will preserve its sign without an infinite-dimensional Cauchy--Binet theorem.

Verification:

- every new public result reports only `propext`, `Classical.choice`, and
  `Quot.sound` in the axiom audit;
- the exact rational regression compares determinant multilinearity,
  fixed-cardinality physical-column sums, and ordered feature grouping;
- `./scripts/verify.sh` passes with 95 Python tests and no known dependency
  vulnerabilities.

## Phase CY contract: adjacent-transfer interlacing bridge

State: complete; merged as [PR #105](https://github.com/tripp-smith/gecp-kernel-structure/pull/105) with green CI.

Phase CY targets the one-sign-change prefix cone isolated by Phase CX. Write a
prefix/suffix determinant pencil as `A + B*t + C*t^2`. Moving one rank-one
Taylor feature from the suffix to the prefix changes its coefficients to
`A+l`, `B-l+r`, and `C-r`, where `l` and `r` are the feature's matched sums
against the lower and remaining upper blocks.

Frozen public Lean targets:

- `quadraticTransfer_eq` proves the exact polynomial update
  `P_next(t) = P(t) + (1-t)*(l+r*t)`;
- `quadraticTransferNodeNumerator_invariant` proves that the division-free
  scaled value `A*r^2 - B*l*r + C*l^2` is shared by both adjacent cutoffs;
- `quadratic_discriminant_nonneg_of_transferNode` proves one Lorentzian
  discriminant from positive transfer sums, nonnegative endpoints, and a
  nonpositive canonical transfer-node numerator;
- `quadraticTransfer_discriminants_nonneg_of_transferNode` proves both adjacent
  cutoff discriminants from the same certificate.

Structural refinement after the scalar proof audit:

- `quadraticTransferNodeNumerator_outerFeature_eq` factors the canonical
  numerator as the negative product of projected row and column cross
  determinants of the lower and remaining upper blocks;
- `finTwoPencil_transfer_lorentzian_of_projectedCross` derives both adjacent
  inequalities from aggregate projected-cross sign agreement and explicit
  endpoint/transfer hypotheses. This reduces the next analytic problem to
  cumulative projections, without assumptions about individual coordinates.

The phase name refers to an interlacing route. The delivered theorems prove
the discriminant inequalities; an ordering theorem for adjacent roots is not
asserted.

Independent exact check:

- extend the deterministic conditioned-Taylor corpus to reconstruct `l` and
  `r` from consecutive cutoff coefficients, verify the transfer identities,
  and audit the canonical numerator on every interior transition;
- retain the Phase CV rank-one boundary as an explicit exception rather than
  weakening or hiding the certificate.

Dependencies:

- Phase CL's quadratic discriminant interface;
- Phase CW's aggregate fresh-order-two Lorentzian bridge;
- Phase CX's proof that only the prefix cone, not the full feature space, can
  support the target reverse Cauchy--Schwarz law.

Non-claims:

- the phase does not prove the canonical transfer-node sign for arbitrary
  complete-pivot prefixes or the infinite Taylor tail;
- the exact corpus remains evidence, not a theorem;
- no higher-fresh-order envelope, dyadic recurrence, Conjecture G1, or Simons
  Problem 4.2 conclusion is claimed.

Focused evidence and revised decision:

- The 2,500-case deterministic finite-Taylor corpus checks 13,924 adjacent
  transfer identities and positive transfer sums. All 9,170 transitions with
  `2 ≤ moved_degree < feature_count - 2` have negative canonical numerator.
  The projected-cross factorization is checked exactly on every transition,
  with the same cross-inversion orientation as the determinant coefficients.
  This corpus uses sizes four through eight, three through twelve retained
  features, five additional Taylor terms in the pivot kernel, and the first
  two available rows and columns after each selected complete-pivot prefix.
  It does not exhaust all remaining row/column pairs or infinite expansions.
- A universal transfer-node sign already fails next to rank-one boundary
  blocks. The Phase CV witness has 35 transitions; its sole positive numerator
  occurs at feature count four, moved degree one, where the old lower
  determinant vanishes. Its 21 strict interior transitions pass.
- The next decision is to prove projected-cross sign agreement for arbitrary
  interior conditioned Taylor blocks, then pass to the infinite tail. Neither
  the finite corpus nor the generic bridge establishes that analytic sign.

Local verification on 2026-10-19: `./scripts/verify.sh` passes a 3,602-job Lean
build, the public axiom audit (only `propext`, `Classical.choice`, and
`Quot.sound`), Ruff, formatting, mypy, 123 tests in 101.16 seconds, and the
dependency audit with no known vulnerabilities. The local package is skipped
by pip-audit because it is not a PyPI dependency. Both push and PR CI passed;
the change merged as PR #105.

## Phase CZ contract: interior transfer-node obstruction

State: complete; merged as [PR #106](https://github.com/tripp-smith/gecp-kernel-structure/pull/106) with green CI.

The proposed universal interior sign after Phase CY is false even for a
complete-pivot prefix of a physical-unit-domain finite Taylor kernel. This
phase freezes the minimized exact obstruction before further analytic work.

Named Lean targets in `Fermionic.InteriorTransferNodeObstruction`:

- `phaseCZ_firstPivot_complete` and `phaseCZ_secondPivot_complete` certify
  the two endpoint pivots of the nine-term Taylor kernel on row numerators
  `(17,18,19,20)/21` and column numerators `(15,16,17,19)/20`;
- `phaseCZ_interior_transferNode_pos` proves that the canonical numerator
  for moving degree three is strictly positive, with both adjacent block
  determinants positive;
- `phaseCZ_alternate_transferNodes_neg` proves that degrees two and four
  supply negative canonical numerators, so alternate adjacent certificates
  remain available;
- `phaseCZ_neighbor_discriminants_pos` verifies that the failed certificate
  does not imply failure of the two target discriminants.

Independent evidence:

- reconstruct the kernel, selected pivots, two conditioned feature factors,
  coefficient transfers, failed interior node, and alternate certificates
  using exact rational arithmetic;
- extend sampling beyond the earlier first-fresh-pair, size-at-most-eight
  corpus, retaining failed signs explicitly;
- a 500-case exploratory audit over every remaining fresh pair observed
  339,950 transfers, six positive canonical numerators, and neighboring-node
  coverage for all 292,076 nontrivial interior cutoffs. A reproducible bounded
  regression will accompany the minimized witness.

Decision: stop pursuing a universal fixed-node projected sign. Test and then
prove a weaker adaptive coverage principle, or a biorthogonal/oscillation
representation of the determinant pencil. Related primary literature:
[Bertola, Gekhtman, Szmigielski, Cauchy Biorthogonal Polynomials](https://arxiv.org/abs/0904.2602).
Its general total-positive-kernel zero theorem and Cauchy-specific interlacing
are useful context; neither is applied as a theorem about these conditioned
Taylor pencils without a verified representation.

Non-claims: this finite obstruction does not refute the full exponential
infinite-tail inequality, the adaptive coverage conjecture, or Simons Problem
4.2. Numerical tail stabilization is evidence only. Adjacent-node coverage on
a finite corpus does not prove it for arbitrary complete-pivot runs.
Complete pivoting here means maximization on the specified finite node sets;
no continuous-domain pivot-maximization claim is made.

Delivered evidence:

- `phaseCZLeftQ_mul_rightQ` proves that all retained features reconstruct the
  pivot kernel, and `phaseCZTransferQ_cast` connects the rational normalization
  to the real-valued Phase CY certificate. The exact obstruction, positive
  endpoints, alternate signs, and surviving discriminants all compile with
  only the accepted Lean axioms.
- The committed default audit covers 54,845 transfers and 47,496 interior
  cutoffs, retaining four failed fixed nodes. Its seed is `0xC212C0DE` and it
  tests every remaining row/column pair in 100 configurations.
- A separate 80/160-decimal regression uses the actual exponential kernel and
  obtains the tail block by subtracting the low features from its two-pivot
  residual. The normalized degree-three numerator agrees to 35 significant
  digits at `0.00010354916980416281785667342918876537`, while the neighboring
  nodes remain negative. This is numerical evidence about the infinite tail,
  not a formal extension of the finite obstruction.

Local verification on 2026-10-19: `./scripts/verify.sh` passes a 3,603-job Lean
build, the public axiom audit with only the accepted foundations, Ruff,
formatting, mypy, 126 tests in 132.87 seconds, and the dependency audit with no
known vulnerabilities. The local package is skipped by pip-audit because it
is not a PyPI dependency. The committed 500-case audit command independently
reproduces 339,950 transfers, six failed fixed nodes, and coverage of 292,076
interior cutoffs. Both push and PR CI passed; the change merged as PR #106.

## Phase DA contract: one-pivot exponential Chebyshev sign

State: complete; merged as PR #107 with green push and PR CI.

The proposed localization of bad transfer nodes at `pivot_depth + 1` fails
in deeper exact tests. At depth ten, a twelve-node Taylor case has a positive
degree-ten numerator; neighboring nodes still cover its cutoffs. The expanded
400-case exploratory audit through depth eighteen checks 6,364 transfers and
covers 5,964 interior cutoffs, with sixteen positive nodes, ten outside the
proposed exceptional degree. These observations do not prove adaptive
coverage. The explicit depth-ten configuration will be retained as a bounded
exact regression.

The analytic base case uses the actual exponential rather than a Taylor
truncation. For `0 < a ≤ 1`, `0 ≤ x < y < z ≤ 1`, and natural degrees `i < j`,
the determinant with columns `exp(a*t)`, `t^i`, `t^j` has a sign depending only
on the degree pair: negative for `i=0` and `j≥2`, positive otherwise. Thus two
independent ordered augmented triples have a positive matched product at one
pivot, uniformly over every pair of degrees. The claim is false as a generic
multi-pivot sign rule, as Phase CV already demonstrates.

Named Lean targets:

- `Matrix.threePointCross_pos_of_derivativeRatio`: Cauchy mean-value transfer
  from strictly increasing derivative ratio to a three-point cross sign;
- `expMonomialMinor_pos_of_pos_degrees`: the positive-degree case;
- `expMonomialMinor_neg_of_zero_left`: the zero-degree/high-degree case;
- `expMonomialMinor_zero_one_pos`: the remaining degree `(0,1)` case;
- `expMonomialMinor_signed_pos` and `expMonomialMinor_matched_mul_pos`: unified
  degree schedule and matched-product sign for two independent unit triples.

Method: divide each sampled row by its positive exponential factor and apply
Cauchy's mean value theorem on the two adjacent intervals. For positive
degrees the derivative ratio is
`t^(j-i)*(j-a*t)/(i-a*t)`. Its derivative has the positive factor
`i*j-(i+j-1)*a*t+(a*t)^2` on the open unit interval. The zero-degree cases have
explicit derivative ratios of opposite monotonicity.

Independent checks: high-precision actual-exponential determinants over
boundary and interior unit triples, all degree-pair sign branches, and exact
rational reconstruction of the failed localization example.

Scope: this closes an analytic one-pivot augmented-minor base case. It does
not establish arbitrary-prefix adaptive coverage, the higher-order companion
envelope, dyadic recurrence, G1, or Simons Problem 4.2. No percentage for the
remaining research effort is inferred from completing this base case.

Delivered proof: `OnePivotChebyshevSign.lean` establishes all six named targets
using the actual exponential, with no Taylor-tail limit assumption. The two
matched minors may even use independent rates in `(0,1]`. The normalized
derivative bracket is strictly positive because, with `i=m+1`, `j=m+n+2`,
and `s=a*t`, it equals
`1 + m*(m+n+1) + (2*m+n)*(1-s) + (1-s)^2`.
The zero-left/high-right branch is strictly negative; omitting it would give
an incorrect all-positive Chebyshev assertion. The boundary nodes `0` and `1`
are allowed: the derivatives need their strict sign only on the open interval.

Independent regressions in `tests/test_one_pivot_chebyshev.py` compare 324
raw exponential determinants and their row-normalized form at both 80 and
160 decimal digits, agreeing to 35 significant digits. Another 66 pairs use
independent rates and independent node triples. These are numerical checks,
not formal certificates or an infinite-family proof.

The exact localization witness uses rows
`(3,4,7,12,13,14,18,20,22,25,27,31)/34`, columns
`(1,2,3,6,11,14,16,19,21,26,28,31)/33`, 23 pivot-kernel Taylor terms,
and 18 retained features. Its depth-ten complete-pivot rows are
`(11,0,6,9,2,3,10,1,8,5)` and columns `(11,0,7,4,9,3,10,2,6,8)`.
The remaining pairs are `(4,7)` and `(1,5)`. Its positive degree-ten node
has normalized numerator between `0.0024` and `0.0025`; degrees nine and
eleven are negative and all nontrivial cutoff discriminants remain positive.
This exact-rational regression is not a new Lean-certified obstruction.
The 18 features do not reconstruct the 23-term pivot kernel, so this witness
was independently retested with all 23 features retained. The same selected
pivots and positive degree-ten numerator persist, now normalized to
approximately `0.0024093902186396532`. The committed regression uses this
coherent 23-term variant and checks exact fresh-block reconstruction, positive
endpoint determinants, positive discriminants, and adjacent-node coverage.
It therefore rejects degree localization even for a coherent finite Taylor
kernel. It is still an exact-rational regression, not a Lean theorem or an
infinite-tail claim.

The exploratory 400-case sample uses seed `0xDA10CA7E`, depth `2+case%17`,
size `depth+2`, `depth+8` retained features, and five additional pivot-kernel
terms. Every third case clusters nodes near one: each denominator is chosen
by `randrange(4*size,8*size)`, and nodes are a sorted `sample` of size `size`
from the final `2*size` positive numerators below that denominator. Other
cases use `randrange(size+1,4*size)` and sample the full positive numerator
range. Row then column denominators are drawn before row then column samples.
Only the explicit coherent witness, not this broader exploratory generator,
is installed as a bounded test.

Decision: retain the analytic one-pivot sign theorem as a base case, but do
not induct on a universal degree sign or a fixed exceptional degree. The next
critical lemma is still arbitrary-prefix cumulative neighboring-node coverage
(or a valid oscillation representation implying it). Neither the numerical
coverage observations nor this three-node theorem proves that lemma.

Local verification on 2026-10-19: `./scripts/verify.sh` passes a 3,605-job
Lean build, the public axiom audit using only `propext`, `Classical.choice`,
and `Quot.sound`, Ruff, formatting, mypy, 129 tests in 132.74 seconds, and the
dependency audit with no known vulnerabilities. The local package is skipped
by pip-audit because it is not a PyPI dependency. Push and PR CI both pass on
the final branch head; PR #107 merged on 2026-10-19. Phase DA is complete.

## Phase DB contract: adaptive transfer coverage decision

State: complete; merged as PR #108 with green push and PR CI.

Attack the remaining fresh-order-two invariant directly. Let a finite positive
Taylor kernel on independently ordered open-unit row and column nodes use all
of its `m` features for pivoting and reconstruction. After an exact complete-
pivot prefix, select any two fresh rows and two fresh columns, orient their
matched coordinates consistently, and form every degree-prefix pencil.
For every strict interior cutoff `2≤k≤m-2`, the proposed coverage principle is
`T_(k-1)≤0 or T_k≤0`, where `T` is Phase CY's canonical transfer numerator and
the corresponding left and right transfer sums must be positive. This would
provide the two-sided Lorentzian certificate needed before attempting higher
fresh order. It must not be replaced by the false fixed-node sign or the false
`depth+1` localization assertion.

First test the complete-pivot/coherent finite statement at larger depths and
under asymmetric and clustered samples, preserving exact rational arithmetic.
If it fails, minimize and independently reconstruct a witness, then certify
the exact failed coverage and surviving discriminant (if positive) in Lean.
If it survives, identify a structural recurrence or oscillation invariant
which implies coverage and prove it without assuming coverage itself.

Provisional named outputs: `conditionedTaylor_transferCoverage` if the general
finite theorem is proved; otherwise explicit witness-specific results such as
`phaseDB_adjacent_transferNodes_pos` and `phaseDB_cutoff_discriminant_pos`.
The committed Python regression must reconstruct the actual pivot kernel,
both conditioned factors, every claimed pivot, and the canonical numerator
from independent block determinants. No tautological certificate wrapper
counts as this phase's theorem.

Non-claims: a finite-Taylor obstruction need not refute the actual exponential
infinite tail; decimal continuation is evidence only. A surviving finite
corpus cannot prove the arbitrary-prefix theorem. Any one-pivot-only result
cannot close this contract. The higher-fresh-order envelope, dyadic recurrence,
G1, and Simons Problem 4.2 remain outside what a rank-two base decision alone
would establish.

Scope audit before drawing a conclusion: the actual restart interface
`positiveSeriesFermionicKernel_finalResidual_biFeatureSplit` accepts a
`Run` on `PositiveArgument`, not a run confined to the fresh unit band.
Consequently even a proof of all-unit coverage would leave outer-prefix
transport open. The same phase must stress-test outer-scale selected nodes
while requiring the fresh sample to remain local, and must distinguish
unweighted exponential pivots from the column-weighted fermionic pivots.
The all-unit coherent search is a base test, not a substitute for that
stronger application requirement.

Delivered decision: the proposed neighboring-node coverage is false even
before outer-prefix transport. The coherent 13-term Taylor witness has rows
`(7,11,12,18,19,21,22,23,27)/29` and columns
`(8,9,11,13,24,33,39,57,60)/63`. Its seven complete pivots select rows
`(8,0,3,2,7,5,1)` and columns `(8,0,6,4,7,3,5)`. The fresh sample is rows
`(6,4)` and columns `(1,2)`; the reversed row order fixes determinant
orientation. At cutoff eight, **both** adjacent canonical numerators `T7`
and `T8` are strictly positive, with all four transfer sums strictly positive.
The three neighboring pencils have positive endpoint and mixed coefficients,
and the middle discriminant is strictly positive. The failure therefore
refutes the sufficient coverage mechanism, not the underlying Lorentzian
inequality.

`AdaptiveCoverageObstructionData.lean` stores an integer fraction-free
certificate, independently regenerated by
`uv run python scripts/generate_coverage_certificate.py --check`.
`AdaptiveCoverageObstructionCore.lean` checks its exact Schur recurrence,
complete-pivot maxima, open-unit nodes, and coherent Taylor reconstruction.
The kernel is uniformly scaled by the positive integer `(29*63)^12 * 12!`;
this changes neither pivot locations nor the signs in question.
`AdaptiveCoverageObstruction.lean` checks both actual feature projectors and
the conditioned reconstruction. Its numerator matrices have a common nonzero
305-digit normalization factor; homogeneity transports the computed signs
back to the actual conditioned blocks. `phaseDBTransferQ_cast` identifies
the rational expression with Phase CY's real canonical numerator.
The closed arithmetic uses `decide +kernel`, not `native_decide` or a new
axiom. The public result is `phaseDB_not_adaptive_coverage` with
`phaseDB_actual_endpoints_and_transfers_pos` and
`phaseDB_actual_cutoff_discriminant_pos` certifying its nondegeneracy.

The exploratory coherent-unit search found the first failure at case 214
with seed `0xDB202610`, depth seven and 18 features. Removing unused features
retains failure down to 13 terms. Deleting one selected row and column did
not produce a smaller witness in that bounded deletion search; no global
minimality claim is made. The prior 100/500-case Phase CZ coverage counts
remain correct for those fixed corpora, but do not extrapolate to all runs.

Independent actual-exponential calculations at 160 and 320 decimal digits
agree to 45 significant digits on normalized numerators
`0.00733642663708844273392090806344167651521625274` and
`0.0154696286973923495855634025862919067083240956`.
The tail block is the full exponential residual minus its conditioned low
features, not a finite-tail approximation. Pivot choices have no near tie at
the tested precision. This is evidence, not a formal infinite-tail theorem.

The outer-prefix stress test also retains a separate exact weighted fixture:
rows `(4,6,8,12)/13`, columns `(1/5,2/5,13/15,22/3)`, and twelve features.
Column weights are `1/(1+P12(y))`, so all features reconstruct the same scaled
Taylor kernel used for complete pivoting. Pivots `(3,2)` and `(2,3)` leave a
local fresh sample in the unit square. Both degree-five and degree-six
numerators are positive and both left transfer sums are negative, while
the cutoff-six discriminant remains positive. The actual fermionic positive-
coordinate kernel reproduces these signs at 160/320 digits. This weighted
fixture is an exact regression plus numerical evidence, not a Lean theorem.

Related work considered: discrete multiple-orthogonal-polynomial interlacing
and Poisson/Charlier zero separation, especially
[Haneczok and Van Assche](https://arxiv.org/abs/1108.3917) and
[Ndayiragije and Van Assche](https://arxiv.org/abs/1108.3918).
Their recurrence and orthogonality hypotheses have not been established for
these conditioned factors. No theorem from those papers is applied here.

Next decision: stop pursuing neighboring canonical-node coverage as a general
finite-prefix invariant. Either prove the cumulative cutoff Lorentzian
inequality directly through a verified oscillation/biorthogonal
representation, or return to the sign-free nine-weight normalized companion
envelope. A proof may need continuous-domain complete-pivot inequalities:
these formal counterexamples use complete maxima on their stated finite
grids, not on the entire continuous domain. The full exponential inequality,
higher-fresh-order bound, dyadic recurrence, G1, and Problem 4.2 remain open.

Local verification on 2026-10-19: `./scripts/verify.sh` passes a 3,608-job
Lean build, the expanded public axiom audit with only `propext`,
`Classical.choice`, and `Quot.sound`, Ruff, formatting, mypy, 133 tests in
133.09 seconds, and the dependency audit with no known vulnerabilities.
The local package is skipped by pip-audit because it is not a PyPI dependency.
The certificate generator's `--check` command independently reproduces every
committed integer table. Push and PR CI both pass on the final branch head;
PR #108 merged on 2026-10-19.

## Phase DC contract: direct sign-free companion envelope

State: complete as the explicitly weaker all-orders one-step bridge below;
merged as PR #109 with green CI. The original companion envelope remains open.

Return to the quantitative premise consumed by the Phase BS/BZ norming
theorems. For the fixed original positive-series fermionic run and each of
the nine weights `56,65,...,128`, seek the bound
`abs(det(sample(conditionBoth(H_lambda)))) ≤ 2^18` on an eighteen-row,
eighteen-column fresh local unit-band sample. Previously selected nodes may
lie outside that band. Pivot choices must come from the original kernel,
never from the amplified companion. The one-sided Phase BS variant is a
separate interface and must not be silently identified with the two-sided
operator.

First stress-test the actual exponential and denominator weights, not just
short Taylor surrogates, at deeper cross-scale prefixes. Use independent
precision runs and factorization-based determinants; singular or unstable
numerical output is uncertified. Test interpolation-projector/entry bounds
as possible sufficient ingredients, but a scalar or rank-two result alone
cannot close this contract. A different explicit constant fresh order is
acceptable only if it still yields the SPEC's cutoff-uniform contraction
mechanism, not as a substitute for that mechanism.

Provisional output: a genuine quantitative theorem discharging a companion
envelope, or a certified counterexample to its stated scope. No wrapper
whose assumptions merely assert the desired envelope will count as progress
on the phase's theorem. If a promising sufficient estimate fails, retain
its witness and distinguish it from failure of the determinant target.

Non-claims: finite-grid complete maxima do not certify continuous-domain
maxima; numerical checks do not prove an arbitrary-prefix or infinite-tail
bound. Do not claim the dyadic recurrence, G1, or Problem 4.2 until the actual
quantitative hypotheses and their domains have been proved.

Promising alternative identified before implementation: a principal singleton
Hadamard--Fischer inequality for a strictly totally positive matrix gives
`det A ≤ A[i,i] * det A[delete i,delete i]`. For a scalar Schur step and a
fixed local column tuple, choose the row having the same sorted rank as the
pivot column in the augmented sample. Its pivot-column entry is at most the
chosen complete pivot. This suggests contraction of the maximum local
minor of every fixed order, rather than a companion estimate. If valid and
transported through the residual's sign gauge, it would preserve an initial
local analytic determinant budget without charging the prefix length.
This is a route toward the original convergence goal, not a change of that
goal. The precise sampling/order/gauge obligations must all be proved.

The first exact tests cover 1,000 diagonally scaled geometric TP matrices
and 400 positive Jacobi-product TP matrices, fresh orders two and three,
and all local column tuples avoiding the pivot. None violates contraction.
These finite tests select the next proof, but are not its proof.
The classical matrix inequality is discussed in
[Skandera and Soskin, Barrett--Johnson inequalities](https://arxiv.org/abs/2209.06466).
No external theorem will be imported as an axiom. The next named target is
`Matrix.det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive`,
followed by the local-minor contraction theorem if its full hypotheses hold.

DC proof decision: the direct companion envelope is not discharged by the
stress tests. Deliver an explicitly weaker, all-orders one-step bridge rather
than claim that envelope. The frozen finite statement is: for an ordered
`(n+1)`-sample of a strictly totally positive kernel and a column-maximal
pivot at sorted ranks `(p,q)`, the absolute determinant of the complementary
`n`-sample of its scalar Schur residual is at most the original minor deleting
row `q` and column `q`. The local columns therefore stay fixed; the old row
tuple can include the pivot row. Public outputs are
`Matrix.abs_det_pivotSchur_le_matchedComplement_of_columnMax`,
`GECP.localMinor_schur_le_of_strictTotalPositive`, and the actual positive-
series fermionic specialization. This is all orders, not a scalar/rank-two
surrogate, but it still loses arbitrary-prefix scope. The next precise lemma
is preservation of all-orders strict total positivity by the actual residual's
single magnitude-preserving row/column sign gauge. Independent compound-order
gauges from Phase AL do not establish that premise. Until that transport and
the initial dyadic determinant budget are proved, no cutoff-uniform residual
contraction or completion of DC's original quantitative target is claimed.

Delivered proof: singleton Fischer follows by northwest-Schur induction at
positive indices, then simultaneous reversal for index zero. The arbitrary
Schur determinant identity uses `cycleRange.symm` to move each pivot index to
zero; independent row/column permutations preserve its absolute determinant.
Fischer at index `q` and `A[q,q] ≤ A[p,q]` yield the matched complement bound.
The kernel theorem and actual positive-series fermionic specialization build
strict total positivity from ordered kernel samples, not from the desired
determinant conclusion. The separate budget corollary transports an **old**
same-column envelope into the new minor and is not a wrapper assuming a new
residual envelope.

Independent checks: 2,520 exact one-step cases through fresh order four, plus
1,376 exact complete-prefix local-minor/gauge cases. The latter support a
precise next contract: on coordinates fresh for a positive-series run, multiply
its residual by `prod_j sign(u-u_j) * prod_j sign(y-y_j)` and prove strict
total positivity at every order on the inherited coordinate order. Each
factor has magnitude one on fresh nodes. Arbitrary selected nodes may lie
outside the local band. A possible proof factors the orientation of a bordered
tuple using its Vandermonde cross differences; this is a research direction,
not an established lemma. Once available, apply the one-step theorem backwards
through the prefix to preserve every fixed-column maximum-minor budget.
Then prove an initial scaled dyadic-band budget and compose the existing
pigeonhole/continuation inequalities. Those are distinct remaining obligations.

Direct companion evidence: the nine configurations `(depth,cutoff)` are
`(0,1)`, `(2,16)`, `(8,16)`, `(16,16)`, `(32,16)`, `(8,256)`, `(16,256)`,
`(32,256)`, `(64,256)`. Both 320- and 640-digit runs select identical pivots
and reproduce all 45-significant-digit normalized order-eighteen determinant
ratios for all nine weights. No sampled ratio exceeds one. At depth 64 the
largest ratio is `1.0326550577e-1234`, the largest companion entry is about
`5.987823e-41`, and the minimum relative pivot gap is about `1.79556e-108`.
Pivots are never recomputed from an amplified kernel. Determinants use
factorizations, not explicit inverses. Zeros are rejected as uncertified;
conditioning is not estimated and decimal agreement is not interval proof.

Candidate initial dyadic budget (not yet formalized): on `omega in [A,2*A]`,
`A>0`, write the original kernel as
`exp(-2*A*t) * exp(t*(2*A-omega)) / (1+exp(-omega))` and reverse the frequency
order. For every Taylor degree `k>=r`,
`(t*(2*A-omega))^k <= 2^(-r) * (2*A*t)^k`. Positive term comparison therefore
suggests a damped degree-`r` tail bounded by `2^(-r)`, independently of `A`.
Amplifying that tail by `2^r` would keep the initial kernel inside envelope
two. A generalized rank-`r` amplification theorem would then give the initial
local determinant budget `2^n / (2^r)^(n-r)` for `n>=r`. This avoids a
Poisson-tail asymptotic estimate or subdivision into smaller subbands. A
fixed rank-eight budget alone does not give the desired arbitrary-accuracy
rate: the degree must vary with accuracy, and the mixed-coefficient and
rank-amplification interfaces must be generalized or located before claiming
that bound. The global gauge theorem must cover signed frequencies and the
closed time interval, not only `PositiveArgument` samples used by DC's
convenient concrete specialization. Finally an exact count of band-local
pivots above a threshold is needed; a same-column bound is not itself G1.

Verification: `./scripts/verify.sh` passes a 3,610-job Lean build, accepted
public theorem axioms, Ruff, formatting, mypy, 140 tests in 135.08 seconds,
and no known dependency vulnerabilities. Push and PR CI both pass on the
final branch head. PR #109 merged on 2026-10-19. This verifies the delivered
one-step bridge and evidence, not the unproved original envelope or G1.

## Phase DD contract: a common all-orders residual sign gauge

State: complete; merged as PR #110 with green CI. The common all-orders gauge
and actual arbitrary-prefix next-step bound below are proved; G1 remains open.

For a real-coordinate kernel `K` with `StrictSignRegular K (fun _ => 1)`
and any successful exact `Run K`, define
`P_run(x) = prod_i (x - selectedRow_i)` and the analogous column polynomial
`Q_run(y)`. On fresh coordinates, prove that the single kernel
`sign(P_run(x)) * finalResidual(x,y) * sign(Q_run(y))` is strictly totally
positive at **every** order, and preserves the absolute determinant of every
finite minor. Do not assume minor positivity of the gauged residual. Successful
pivots need not be complete, and selected nodes need not lie in a local band.
The concrete specialization must cover the actual reflected fermionic kernel
on all real coordinates, including time endpoints and signed frequencies.

Named targets: `GECP.Run.gaugedFinalResidual_strictSignRegular`, its absolute-
minor preservation theorem, and a real-coordinate reflected-fermionic
specialization. Planned sources are a reusable Vandermonde block factorization
and `GECP/ResidualTotalPositivity.lean`; independent rational regressions will
exercise complete and deliberately non-complete prefixes, all minor orders,
empty prefixes/minors, boundary time nodes, and signed frequency labels.

Proof route selected before implementation: the sign of any injective kernel
minor is the product of the signs of the two Vandermonde determinants.
Factor each augmented Vandermonde into the selected determinant, the fresh
determinant, and the selected/fresh cross differences. The exact multi-bordered
Schur identity then cancels the selected determinants and leaves the two
polynomial weights above. Normalize by their nonzero absolute values on fresh
coordinates. This uses one gauge for all orders, unlike Phase AL's independent
compound-order gauges. Search mathlib before adding product/indexing helpers.

Non-claims: this phase does not prove the initial variable-degree dyadic
determinant budget, the count of band-local pivots above a threshold, the
companion envelope, G1, or Problem 4.2. It discharges a structural premise of
the Phase DC route without replacing the full convergence objective.

Planned immediate application, before its implementation: combine the proved
common gauge with Phase DC and diagonal-scaling covariance of a scalar Schur
update. `GECP.Run.finalResidual_localMinor_schur_le` should give the actual
residual's matched-rank same-column bound after an arbitrary successful
prefix, assuming absolute column maximality for the next pivot. This must be
a bound for the actual next residual, not merely for an auxiliary signed
kernel. Iteration to an initial band budget and pivot counting remain separate.

Delivered proof path: `Matrix.det_vandermonde_append` factors the block
determinant over any commutative ring, without injectivity or sign premises.
Injective real tuples have Vandermonde determinant signs equal to their sorting
orientations. Multiplying an original-kernel minor by both Vandermonde
determinants therefore gives a positive number. Apply this to the selected
core and augmented core; the exact selected-core/residual identity cancels the
selected factors. The remaining raw weights are exactly `P_run` and `Q_run`.
Freshness makes them nonzero, so their absolute values provide a positive
normalization to one fixed sign gauge. Absolute minors are preserved even for
unordered or repeated fresh samples. This proves the frozen all-orders target,
not an additional assumption asserting it.

Diagonal-scaling covariance of a scalar Schur update then transfers Phase DC
back to the actual residual. `Run.finalResidual_localMinor_schur_le` requires
only the next absolute-column maximum on the sampled augmented rows. No
complete-pivot, localization, or positivity assumption is added for earlier
pivots. `Run.negateColumns` transports every original fermionic run exactly,
and `finalResidual_negateColumns` verifies the inverse reflection; this is not
a new run selected from a companion kernel or a new canonical tie policy.

Independent checks: 25 exact block identities; 3,825 exact gauge minors over
three pivot policies and all prefix depths on the six-by-six signed/boundary
fixture; 350 exact Schur-covariance/matched-rank bounds; repeated-coordinate
zero minors; exact original-run reflection; and 19 actual-exponential residual
minors at both 160 and 320 decimal digits with 45-significant-digit agreement.
Numerical values use factorization-based determinants and do not certify an
arbitrary-input or continuous-domain runtime result.

Next quantitative decision: first try the **existing generic rank-approximation
determinant estimate locally**. Its factorial depends on the fixed local minor
order, not on the dyadic cutoff scale. A variable-degree local Taylor error
`2^(-r)` with rank `r`, followed by the same-column contraction and a band-local
pivot count, may give a linear-in-accuracy-order count per band. This avoids
generalizing positive mixed-coefficient amplification unless the generic bound
actually fails. The existing global bound already has an extra `log(s+1)`
overhead; rederiving that bound is not the missing cutoff-uniform G1 result.
For the local approximation, either prove the positive-series term comparison
recorded under DC, or split each dyadic band into a constant number of narrower
bands and reuse the existing Lagrange remainder/`pow / factorial <= exp`
estimates. Both retain the original logarithmic cutoff target. Neither the
local budget nor its iteration/counting is claimed proved in DD.

Local verification: `./scripts/verify.sh` passes a 3,614-job Lean build,
accepted public theorem axioms, Ruff, formatting, mypy, 146 tests in 133.88
seconds, and no known dependency vulnerabilities. The final branch-head
repeat passes all 146 tests in 135.55 seconds. Push and PR CI both pass;
PR #110 merged on 2026-10-19.

Sharper next-target discovery at DD closure: the existing
`eval_dyadicSeparatedTerms_succ` already isolates **one rank-`8p` panel** on
`2^j < omega <= 2^(j+1)`. It is `taylorSeparatedTerms (8*p)` masked only by
`t*2^j <= p`; dividing the right factors by the fermionic denominator preserves
that rank. The existing `expFamily_separatedApprox_error` supplies error
`2^(-p)` on that panel, while the lower-frequency branch vanishes there.
Thus neither a new scaled Taylor remainder nor a variable-rank positive
companion theorem appears necessary for the initial local budget. Central
and negative-frequency panels can reuse the existing base and reflection.

Proposed explicit quantitative assembly, **not yet a Lean theorem**: put
`p=2*m+1`, local minor order `n=16*p`, and reuse the generic determinant bound
with rank `8*p`. The same odd-exponent cancellation as
`fermionicKernel_gecp_error_le_oddBlock` would give a local determinant budget
`(2*n*2^(-m))^n`, now without a factor `s+1` in `n`. Choosing `m=16+2*q`
suggests a budget at most `2^(-q*n)` with `n=16*(33+4*q)`.
There are `2*(s+1)` disjoint signed dyadic/central panels. If the new
same-column bound can be iterated backwards through intervening global
pivots, at most `n-1` pivots per panel can exceed `2^(-q)`. This suggests
the explicit target `32*(s+1)*(33+4*q)` global pivots for error `2^(-q)`,
removing the earlier `Nat.clog 2 (s+1)` overhead. The panel approximation,
backwards minor construction, partition bookkeeping, and final arithmetic
must all be proved before making that claim. The earlier global odd-block
theorem alone does not establish it.

## Phase DE contract: a cutoff-uniform local panel determinant budget

State: complete; merged as PR #111 with green CI. The initial local budget is
proved; its iteration/counting and G1 remain open.

Prove an initial determinant envelope on each signed dyadic or central
frequency panel, independent of its scale. The positive central panel is
`[0,1]`; positive panel `j+1` is `(2^j,2^(j+1)]`. Negative panels are their
frequency reflections, with zero assigned only to the positive central panel.
For every `p>0`, build exactly `8*p` separated terms with uniform error
`2^(-p)` on each panel and all times in `[0,1]`, reusing the existing masked
Taylor construction and denominator/reflection lemmas.

Named target `fermionicKernel_panel_sample_det_le_two_pow`: for a sample of
order `n>=8*p` whose columns lie in one panel and whose rows lie in `[0,1]`,
prove `abs(det K_n) <= 2^n * n! * 2^(-p*(n-8*p)) * 2^(8*p)`.
Then `fermionicKernel_panel_minor_det_le_dyadicBudget` must prove the explicit
scale-free budget `abs(det K_n) <= (2^(-q))^n` at
`n=16*(33+4*q)`, for every natural accuracy order `q` and every signed panel.
The determinant bound must concern the actual fermionic exponential, not a
Taylor surrogate or a hypothesis asserting the desired determinant estimate.
All sample tuples, including repeated or unsorted tuples, are permitted.

Planned source: `Fermionic/LocalPanelBudget.lean`. Independent checks will
verify the exact masked-feature factorization and integer exponent arithmetic,
exercise panel/time-cutoff boundaries and both frequency signs, and compare
actual-exponential approximation/determinant evaluations at independent
precisions. Large target-order floating-point determinants are not required
and must not be substituted for the all-orders proof. Search existing mathlib
and repository bounds before introducing a new analytic abstraction.

Non-claims: this is an **initial** budget. Its backward transport through
arbitrary prefixes, the count of band-local pivots, the full signed-panel
partition count, G1, and Problem 4.2 remain separate obligations. Do not
claim the proposed global `32*(s+1)*(33+4*q)` count until they are proved.

Delivered proof: `positivePanelTerms` selects the central Taylor list or the
single row-masked dyadic branch, then absorbs the fermionic denominator.
`eval_positivePanelTerms_eq_dyadic` identifies its value with the already
certified dyadic construction on the stated panel. Reflection handles negative
frequencies without changing the `8p` feature count. The approximant entries
are bounded by two from kernel bound one plus approximation error at most one.
The existing generic determinant estimate therefore applies with effective
rank `8p`. Its factorial is retained, bounded by `n^n` for the odd-block
specialization, and absorbed by the explicit accuracy choice. No companion
positivity, ordering, or injectivity premise is introduced.

The frozen determinant target is proved as
`fermionicKernel_panel_minor_det_le_dyadicBudget`, at the exact local order
`16*(33+4*q)` for every signed panel and every natural `q`. The scalar
accuracy step reuses `oddBlock_scale_accuracy 0 q`; its scale is zero only
in an arithmetic lemma, **not** a restriction on the physical panel index.

Independent checks: 345 exact feature-product identities; 435 actual-kernel
approximation errors compared at 160/320 digits; twelve small-order actual
determinants at 320/640 digits; 65 odd-exponent identities; 1,001 scalar
accuracy conditions; seventeen factorial-budget inequalities; and finite
signed-panel boundary membership checks. Precision comparisons use decimal
strings with 45 significant digits. The exact feature test uses rational
column weights only to isolate finite algebra; the actual denominator and
exponential are separately checked. These observations do not replace the
all-orders theorem or certify a continuous-domain runtime computation.

Next precise proof obligation: backwards construction of nonzero local minor
witnesses along an actual complete-pivot sequence. At a pivot inside a chosen
panel, prepend its row and column and use the exact Schur determinant identity
to gain its pivot factor. At an intervening pivot outside the panel, use DD's
matched-rank bound without changing the number of local columns. A nonzero
residual minor gives a nonzero augmented original-kernel minor by the selected-
core factorization; this can certify the needed freshness/injectivity instead
of transporting an entire run through successive signed subtypes. All
replacement rows must remain in the physical time domain. After that proof,
partition `[-2^s,2^s]` into the `2*(s+1)` signed central/dyadic panels and count
pivots above `2^(-q)`. The proposed `32*(s+1)*(33+4*q)` global count remains
unproved until this construction and partition/count bookkeeping are complete.

Local verification: `./scripts/verify.sh` passes a 3,615-job Lean build,
accepted public theorem axioms, Ruff, formatting, mypy, 151 tests in 139.28
seconds, and no known dependency vulnerabilities. The final branch-head repeat
passes all 151 tests in 139.20 seconds. Push and PR CI both pass;
PR #111 merged on 2026-10-19.

Refinement for the next witness proof: if a panel has at least the target
number `n` of above-threshold pivots, select **exactly `n` marked positions**.
Gain a pivot factor only at those positions; use the nonincreasing matched-
minor bound at every other step, including unmarked pivots inside the panel.
This keeps the initial witness order exactly equal to DE's budget order,
rather than incorrectly applying that fixed-order theorem to a larger count.
For freshness, a nonzero matrix minor has no identically zero sampled row or
column. Nonzero residual entries cannot lie on previously selected coordinates
by interpolation, so they can certify the fresh-subtype inputs needed by DD.
These are proposed proof steps, not a substitute for the missing iteration.

## Phase DF contract: backward minor witnesses and the global local-panel rate

State: verified; merged as PR #112 from `codex/phase-df-local-pivot-count-rate`
on 2026-10-19 with green push and PR CI.

Prove backwards transport of a nonzero local minor through an actual scalar
Schur step. For a marked pivot whose column lies in the chosen panel, prepend
the pivot row/column and gain its absolute pivot factor. For every unmarked
pivot, preserve the local minor order using the matched-rank bound, even when
that pivot also lies in the panel. Replacement rows must remain in the original
physical row domain; replacement columns must remain in the chosen panel.
Freshness must be derived from nonzero minors and exact interpolation rather
than assumed as an independent invariant for the entire sequence.

Named generic target: `GECP.markedPivotProduct_le_initialMinor`, constructing
an actual original-kernel sample on exactly the marked-column count whose
absolute determinant dominates the marked pivot product. All previous pivots
are actual original-kernel complete pivots, not pivots of an amplified kernel.
Combine this with DE's initial local budget to bound the number of pivots above
`2^(-q)` in any panel. Prove the signed panel cover/finite count rather than
relying on the finite boundary regression.

Main finite target: `fermionicKernel_gecp_error_le_localPanelBlock`, with
`N=32*(s+1)*(33+4*q)`, proving `abs(R_N(t,omega)) <= 2^(-q)` on
`[0,1] x [-2^s,2^s]` for every realized exact complete-pivot sequence starting
from the actual fermionic kernel. The only algorithm hypotheses are its
initial kernel, nonzero successful pivots, exact updates, realized finite
prefixes, and domain-complete maxima. No residual determinant budget,
companion bound, or convergence conclusion may be added as an assumption.
Frequency-reflection transport must link the positive-minor presentation
to the original sequence without reselecting pivots.

Planned sources: `GECP/LocalMinorWitness.lean` and
`Fermionic/LocalPanelRate.lean`. Independently check backwards witness
construction and exact marked products on rational samples, unmarked pivots
inside and outside a panel, sorted-rank changes, empty witness/mark sets,
signed partition boundaries, and any finite closed-form count arithmetic.
Numerical checks can supplement, but cannot replace, the universal theorem.

Non-claims until proved: the proposed global count is not already implied by
DD or DE alone. A real-logarithmic/asymptotic corollary and a final full-SPEC
completion audit remain separate from the natural-number finite theorem.
If a genuine mathematical obstruction appears, preserve it explicitly rather
than weakening the intended rate silently or marking the full goal complete.

### DF outcome

The named finite target is proved, together with
`fermionicKernel_gecp_error_le_of_localPanelBlock` for every later rank.
The bound is `N=32*(s+1)*(33+4*q)`, on the actual original kernel and closed
physical rectangle. It removes the earlier `Nat.clog` startup overhead.
No determinant envelope, companion certificate, or convergence conclusion is
assumed in these fermionic rate theorems. They assume an exact successful
complete-pivot sequence, realized finite prefixes, and initialization.

`exists_minor_witness_before_step` sorts the augmented samples in DD's common
gauge and deletes the row at the pivot column's sorted rank, not necessarily
the pivot row. `exists_marked_minor_witness` iterates this construction
backwards. Marked pivots prepend their actual row/column; unmarked pivots
preserve order, even if in the same panel. The empty witness has determinant
one. Nonzero minors prove their own freshness using interpolation.
`panel_pivot_count_lt_of_initialBudget` selects exactly the budget order of
marked positions, avoiding any unsupported extrapolation to larger minors.
Its generic initial-budget premise is discharged by DE in the actual rate.

`fermionicPanelDomain_cover`, `fermionicPanelDomain_unique`, and
`fermionicPanel_count` prove the signed finite partition, assigning zero only
to the positive central panel. Pigeonhole counting over `2*(s+1)` panels and
pivot antitonicity give the global bound. Frequency reflection transports
the original chosen pivots; it does not choose a new sequence or tie policy.

Independent evidence: all 64 mark sets on a six-pivot exact rational sample
of the actual positive presentation; both inside/outside unmarked steps and
non-pivot matched-rank row deletion; 135 separately enumerated initial-minor
budgets; three actual-exponential witness determinants at 160/320 digits,
compared as 45-significant-digit decimal strings; signed boundary partition
checks through scale 11 and 780 exact count identities. Zero/nonfinite
determinants are rejected. Finite-grid evidence is not a continuous-domain
runtime certificate, and no enormous target-order numerical determinant is
used to justify the theorem. Root verification passes a 3,619-job Lean build,
the accepted public axioms, Ruff/formatting/mypy, all 155 tests in 139.22
seconds, and no known dependency vulnerabilities. The final branch-head repeat
passes all 155 tests in 139.98 seconds. Push and PR CI pass; PR #112 merged
on 2026-10-19.

Next: generalize the rectangle hypothesis to an arbitrary real cutoff bounded
by `2^s` (a sequence complete on a smaller rectangle is not necessarily
complete on the larger dyadic rectangle). Then derive explicit real-cutoff,
real-accuracy/logarithmic corollaries and audit their quantifiers against G1
and the full SPEC. The dyadic theorem alone must not be advertised as already
proving a theorem for all real cutoffs by padding the pivot domain.

## Phase DG contract: actual real cutoffs and explicit logarithmic rates

State: verified; merged as PR #113 from `codex/phase-dg-real-cutoff-rate`
on 2026-10-20 with green push and PR CI.

Generalize DF's count argument to complete pivots on the **actual** rectangle
`[0,1] x [-Lambda,Lambda]`, using `Lambda <= 2^s` only to label its columns.
Do not assume completeness on a larger rectangle. Derive realized prefixes
from initialization and successful exact updates using the existing
`Run.ofResidualSequenceFrom` API instead of requiring duplicate run evidence.
Named output: `fermionicKernel_gecp_error_le_realCutoffBlock`, with exactly
the same explicit count `32*(s+1)*(33+4*q)` and all-later-rank guarantee.

Source: `Fermionic/RealCutoffRate.lean`. Derive an explicit exponential
corollary at every natural rank, first with denominator `s+1`, then with
`log(1+Lambda)` for `Lambda >= 1`. Proposed safe constants are `C=1024` and
`c=(log 2)^2/256`. For all positive cutoffs, prove a `log(2+Lambda)` version
rather than silently treating `log(1+Lambda)` as bounded below. Give a
real-tolerance rank choice with explicit ceilings/additive constants.

Audit the literal small-cutoff wording of G1: the already proved first
complete pivot `(0,Lambda)` leaves reflected-corner residual
`1-exp(-Lambda)`. Determine rigorously whether universal constants in the
literal `log(1+Lambda)` denominator can cover arbitrarily small positive
cutoffs. Preserve SPEC's mathematical goals; record any obstruction as a
status annotation, not a rewritten conjecture.

Independent checks: exact rank/floor arithmetic, cutoff and accuracy boundary
choices, direct actual-kernel first-pivot values and logarithmic envelope
comparisons at independent precisions. Floating-point checks are evidence,
not a proof or a runtime complete-pivot certificate. Non-claims: sharp
constants, approximate-pivot robustness, automatic numerical certification,
literature priority, or full SPEC completion before the separate audit.

### DG outcome

All frozen mathematical targets compile. The real-cutoff block theorem keeps
maxima on `[-Lambda,Lambda]` and uses the dyadic cover only for column labels.
`GECP.realized_of_initial_updates` constructs each finite prefix from the
existing exact recurrence, eliminating the duplicate realization assumption.
The public rate assumptions are initialization, nonzero successful pivots,
exact updates, and complete maxima on the actual physical rectangle. They
still describe a sequence at all natural indices; finite-only run and exact-
termination interfaces belong in the final assumption audit.

`localPanelBlocks_le_exponential` sets `j=k/(128*(s+1))` and
`q=max(0,j-9)`. For `j>=9` the named block is at most `k`; for smaller ranks
the initial bound one suffices. In both cases `k/(128*(s+1)) < q+10`, yielding
the universal prefactor `2^10=1024`. The ceil-logarithmic cover satisfies
`(s+1)*log 2 <= 2*log(1+Lambda)` when `Lambda>=1`, by
`4*Lambda <= (1+Lambda)^2`. Thus
`fermionicKernel_gecp_error_le_logarithmic` proves
`abs(R_k) <= 1024*exp(-((log 2)^2/256)*k/log(1+Lambda))`.
For every positive cutoff, `fermionicKernel_gecp_error_le_allCutoff_logarithmic`
has denominator `log(2+Lambda)` and the same constants.
`fermionicKernel_gecp_error_le_realTolerance` uses the explicit natural
ceiling of `256/(log 2)^2 * log(2+Lambda) * log(1024/epsilon)`; every
positive tolerance is admitted, including a zero ceiling for large ones.

The literal all-positive-cutoff interpretation of G1 is false. For any
proposed `C,c>0`, choose `Lambda=min(1,c^2/(16*C))`. The actual first complete
pivot `(0,Lambda)` has reflected-corner residual `1-exp(-Lambda)`. Elementary
exponential/logarithm inequalities prove
`C*exp(-c/log(1+Lambda)) <= Lambda/4 < Lambda/2 <= 1-exp(-Lambda)`.
`fermionicKernel_smallCutoff_logarithmic_obstruction` includes the continuous
complete-pivot certificate, not just the scalar comparison. This is a named
obstruction to a parameter range, not a rewrite of SPEC or a refutation of
the large-cutoff near-optimal rate. No external asymptotic theorem is needed.

Independent checks: 6,400 exact integer rank/floor cases, 53 rational cutoff
covers (including both sides of dyadic boundaries), 212 logarithmic-envelope
values, 28 real-tolerance/ceiling cases including tolerances 1024 and 2048,
and nine actual-kernel small-cutoff witnesses. Analytic values are independently
evaluated at 160/320 digits and compared as 45-significant-digit strings.
Zero/nonfinite witness residuals are rejected. The observed finite first-pivot
grid checks are not substituted for Lean's continuous complete-pivot proof.
Root verification passes a 3,622-job Lean build, the permitted public axioms,
Ruff/formatting/mypy, all 159 tests in 139.51 seconds, and no known dependency
vulnerabilities. The final branch-head repeat passes all 159 tests in 136.89
seconds. Push and PR CI pass; PR #113 merged on 2026-10-20. One PR Lean
job failed while downloading the toolchain with a TLS error before compilation;
retrying the failed check passed, without source or gate changes.

Next decision: audit the full SPEC and the theorem's finite-prefix,
termination, and cross-approximation interfaces before declaring completion.
The main large-cutoff rate and small-cutoff obstruction are no longer open
research hypotheses. Earlier companion sign obstructions are not additional
requirements for this successful, different local-minor route.

## Phase DH contract: finite runs, exact termination, and physical error

State: verified; merged as PR #114 from `codex/phase-dh-finite-run-closure`
on 2026-10-20 with green push and PR CI.

Close the future-success assumption in DG. State the rate directly for an
arbitrary finite `Run fermionicKernel` with the existing recursive
`Run.CompleteOn` predicate on its actual physical rectangle. Its rank is its
stored successful pivot count. No independent residual sequence, realized
prefix evidence, extra complete pivot after the final residual, or successful
continuation beyond that rank may be assumed.

Build exact finite prefix extraction/indexed pivot interfaces, reusing the
existing `Run` and `CompleteOn` definitions. Transport DF's backwards marked
minor construction only over indices inside the stored run. Derive the final
residual bound from the last successful update and finite pivot inequalities,
not from a hypothetical next nonzero complete pivot. Use the same initial
panel budget and counts. Named targets in `Fermionic/FiniteRunRate.lean`:
`fermionicKernel_finiteRun_error_le_localPanelBlock` and
`fermionicKernel_finiteRun_error_le_exp`, with DG's explicit constants and
the same real-cutoff parameter range; also an all-positive-cutoff version.

Add the selected-core cross approximation as a named finite-run object and
prove `finalResidual = kernel - crossApproximation` from the existing closed
selected-core Schur identities, not a new inverse axiom. Formal matrix
inverse notation is algebraic only; numerical checks remain solve-based.
An actual complete pivot with zero value on the final residual certifies
exact termination on the physical domain, without attempting a zero-pivot
update. Export the corresponding exact cross reconstruction and the
logarithmic Green-error consequence through the existing integral transfer.
Planned sources: `GECP/FiniteLocalMinorWitness.lean`,
`GECP/RunCrossApproximation.lean`, and `Fermionic/FiniteRunRate.lean`.

Independent checks: every prefix and mark set of small exact rational traces,
empty/one-step/final successful prefixes, a rank-deficient exact-zero stop,
and solve-based selected-cross versus iterative residual/Green sums at
independent precision. Tests must not smuggle in pivots after the target rank.
No approximate-pivot robustness, interval-certified runtime convergence,
sharp constants, literature priority, or full-SPEC completion is claimed by
this phase. A requirement-by-requirement final audit remains separate.

### DH outcome

The finite interface is proved in the three named sources. `Run.take` stores
only a prefix, `length_rows_take` gives its truncated successful count, and
`completeOn_take` inherits completeness. `residualAt_succ`, `pivotAt_ne`, and
`completePivotAt` apply only below the stored length. At and beyond that
length, prefix extraction returns the same final run; no additional successful
pivot is asserted. The finite backwards witness and panel count retain all
domain, fixed-order, freshness, and marked-column invariants of DF.

`finalResidual_abs_le_pivot` derives the bound from the last successful update
and finite pivot comparisons, rather than an extra complete maximum after the
run. `fermionicKernel_finiteRun_error_le_localPanelBlock` consequently uses
only an actual `Run fermionicKernel`, its recursive `CompleteOn` predicate,
the dyadic cover and stored count. `fermionicKernel_finiteRun_error_le_exp`
has the same `C=1024`, `c=(log 2)^2/256`, `Lambda>=1` range as DG, at every
finite successful length including zero. An all-positive-cutoff version uses
`log(2+Lambda)`; the explicit real-tolerance ceiling is unchanged.

`Run.crossApproximation` is the selected-core algebraic cross and
`finalResidual_eq_sub_crossApproximation` proves exact equality to iterative
error using the existing closed Schur formula. No new inverse identity is
assumed. A zero-valued complete pivot on a final residual proves exact domain
reconstruction without attempting a zero-pivot update; the stopped-or-budget
corollary includes this early-stop case explicitly. A floating-point zero or
an exhausted certification budget is not such a proof.

`continuous_finalResidual` and `continuous_crossApproximation` derive finite
continuity. `greenError_le_uniformCrossError` then proves integrability of both
weighted kernels from `rho` integrability, the kernel bound and uniform error,
and bounds the actual **difference of their integrals**. The fermionic
`fermionicKernel_finiteRun_greenError_le_exp` specializes the explicit rate
for any measure supported almost everywhere on the cutoff interval.
The `_interval` corollary specializes to the restricted physical interval
measure and discharges its support property automatically.

Independent evidence: 127 witnesses across every prefix/mark set of a six-step
exact positive-presentation trace, with 642 backwards steps; 504 exact cross
entries and 84 signed Green-sum identities across both kernel presentations;
rank-zero and rank-one exact-stop cases; and 258 determinant/cross values
compared as 45-significant-digit strings at independent 160/320-digit precision.
At each precision, 42 actual Green-integral differences are also recomputed
directly and compared with both the iterative residual sum and exact reference.
The selected cores are solved, never explicitly inverted. Numerical zero or
nonfinite core determinants are rejected. Finite-grid observations do not
certify continuous complete pivots. Root verification passes a 3,625-job Lean
build, the permitted public axioms, Ruff/formatting/mypy, 163 tests in 138.64
seconds, and no known dependency vulnerabilities. The final-source repeat,
including the interval Green corollary, passes all 163 tests in 137.40 seconds.
The final branch-head gate including direct independent-precision Green sums
passes all 163 tests in 139.12 seconds. Push and PR CI pass; PR #114 merged
on 2026-10-20.

Next: requirement-by-requirement completion audit against the original SPEC,
including the optional/conditional approximate-pivot language and the precise
large-cutoff versus literal-small-cutoff claim boundary. The finite-run and
physical-error bridges are no longer missing mathematical premises.

## Phase DI contract: original-SPEC completion audit

State: verified; merged as PR #115 from `codex/phase-di-spec-completion-audit`
on 2026-10-20 with green push and PR CI.

Audit the original mathematical specification (sections 1–30, before the
historical status annotations), not just the recently delivered rate. Write
`COMPLETION_AUDIT.md` mapping every requirement, numbered success criterion,
named public API/artifact, verification gate and non-claim to authoritative
current evidence. Distinguish required outcomes from explicitly suggested
research routes or optional comparisons. Inspect theorem assumptions, not
just identifiers or successful builds. A conditional restatement of the
desired conclusion does not close a requirement.

Recheck all public Python exports and result metadata, exact/high-precision
and failure/certification paths, committed census/application artifacts,
documentation and implementation-run provenance, full root verification,
and branch/PR/main delivery. Add small interface or explicit-quantifier
corollaries only when they genuinely close an identified requirement and
their proof obligation is frozen before implementation. Preserve SPEC's
goals and the proved small-cutoff G1 obstruction; append status only.

If any required item is missing, circular, weaker than the required scope,
or supported only by indirect evidence, record it as not achieved and freeze
the next exact implementation/proof target. Do not label the full goal
complete simply because its principal rate is now proved. If the audit
proves every required item, synchronize final handoff and claim levels and
complete the goal only after merged/pushed exact-head verification.
Non-claims: literature priority, arbitrary approximate-pivot near-optimality
without a theorem, and billing/model details not exposed by the runtime.

### DI bounded closure contracts

The initial audit finds that D1 is currently stated only with natural dyadic
parameters. Add `fermionicKernel_separatedApprox_real` in
`Fermionic/RealSeparatedApprox.lean`: for every `Lambda>=1` and `0<epsilon<1`,
construct a separated list with uniform error at most `epsilon` and length
at most `128*log(1+Lambda)*log(128/epsilon)`. Reuse the actual dyadic Taylor
construction, DG's proved logarithmic cover, and natural ceiling accuracy.
No GECP or approximation conclusion may be added as a premise.

The audit also finds that `scripts/check_axioms.sh` prints dependencies but
does not itself fail on unexpected axioms. Add an output parser with an exact
allow-list (`propext`, `Classical.choice`, `Quot.sound`), expected-declaration
coverage checks, and fail-closed malformed/empty-output handling. Test positive,
custom-axiom, native-reduction-axiom, missing-declaration, and malformed cases.
Reject custom axiom declarations/proof escapes in root/CI source guards too.
This strengthens a gate; it must not hide any currently disallowed dependency.

The baseline currently covers finite PSD diagonal-preferring traces and a
conditional residual-Lipschitz bound, not the original-kernel continuous
baseline of sections 9–10. The census has no `1e-12` cases or normalized
residual-curve plot. Record these as required incomplete targets for subsequent
bounded phases; do not rewrite their requirements or call the goal complete.

Section 15's existing exact two-pivot run certifies its pivot product but not
the complete-pivot inequalities for that concrete grid. Add an explicit
real-valued cast of the same rational `2x2` kernel, a two-step `CompleteOn`
certificate, its entry/determinant positivity, and the exact nonincreasing
pivot list `[2,1/2]` in `SmallInstanceChecks.lean`. These close a finite
certificate interface, not a new general pivot-order conjecture.

The clean-build audit finds local duplicate working notes included by the
default source distribution. Exclude the exact untracked note filenames from
Hatch builds without reading, editing, staging or deleting their originals;
inspect rebuilt source/wheel member lists and smoke-test the wheel. The first
temporary source archive is removed, not published. Add a regression on these
known local-note exclusions.

Close section 7's compact-maximizer interface with
`fermionicKernel_finiteRun_exists_completePivot` in
`Fermionic/CompletePivotExistence.lean`: every finite successful residual has
a complete maximum on the closed physical rectangle for `Lambda>=0`.
Derive continuity from the finite run and apply mathlib's compact extreme-value
theorem; the chosen maximum is permitted to be zero (exact termination).

### DI outcome

`COMPLETION_AUDIT.md` inspects all original sections 1–30 and 40 concrete
requirement groups, including named interfaces/artifacts, gates and parameter
boundaries. After the bounded DI closures, 33 groups are supported and seven
are incomplete: R08/R09 (continuous SPD / original-kernel Lipschitz baseline),
R24 (metadata outside GECP), R25–R27 (tolerance/diagnostic/normalized census
scope), and R31 (listed missing sparse comparisons). About 80% is a checklist-
coverage estimate, not a weighted mathematical completion percentage. The
goal must remain active.

`fermionicKernel_separatedApprox_real` proves D1's exact real quantifiers and
rank bound `128*log(1+Lambda)*log(128/epsilon)` with an actual existing Taylor
list. Compact maximum existence is derived for every finite residual from
its proved continuity and mathlib's extreme-value theorem, permitting zero
termination. Exact `2x2` small-grid completeness, entry/determinant positivity
and pivot list `[2,1/2]` close the concrete Lean certificate surface.

The axiom gate now validates every declared report and the exact allowed
axiom set, rejecting custom/native-reduction axioms and empty/malformed/missing
reports. A shell integration test confirms that Lean's failing exit status
is not hidden even when its reports parse successfully. Root/CI source guards
also reject local custom axiom declarations and native proof escapes.

The first local source archive included untracked duplicate notes, so that
temporary archive was removed and exact filename exclusions were added.
Rebuilt source/wheel member lists contain neither notes nor local caches;
the isolated Python 3.12 wheel smoke imports all 42 exports and completes a
two-pivot GECP example outside the checkout. Original notes are unchanged;
no archive was published. New tests cover parser failure modes, D1 parameter
arithmetic at independent precision, note exclusions and audit accounting.
Final-source root verification passes 3,627 Lean jobs, all 838 declared axiom
reports under the enforced allow-list, Ruff/format/mypy, all 177 tests in
135.57 seconds, and no known dependency vulnerabilities. The final branch-head
repeat passes all 177 tests in 135.75 seconds. Push and PR CI pass; PR #115
merged on 2026-10-20. The seven original-SPEC gaps remain incomplete.

Next highest-priority exact proof: an arbitrary-domain PSD/SPD kernel baseline
and fill-distance control from **original** kernel regularity, not a residual-
Lipschitz premise. Preserve v1's valid finite canonical/conditional results
and add the missing stronger scope as a new phase. Numerical census/metadata
and Phase-II comparisons then need separate implementations and re-audit.

## Phase DJ — original-kernel fill-distance control

Frozen contract (SPEC sections 9–10, audit R09 and the PSD portion of R08):
add `PositiveDefinite/KernelGeometry.lean` and
`PositiveDefinite/KernelFillDistance.lean`. Define kernel positivity by all
finite sampled Gram matrices using mathlib's `Matrix.PosSemidef`; samples may
repeat. For any arbitrary-domain successful `GECP.Run` whose pivots are
diagonal, derive positivity of both the final residual and the removed
kernel `K-R`. No completeness or residual-regularity premise is required.
Derive the squared power-function bound
`R(x,x) <= K(x,x)-2*K(x,s)+K(s,s)` at every selected node `s`.
For a metric domain, original-kernel coordinatewise `LipschitzWith L` and a
selected-node cover of radius `h`, prove all residual entries are at most
`2*L*h` in absolute value. This must not assume a Lipschitz residual diagonal.

Also lift PSD diagonal/complete bound equivalence to arbitrary domains and
derive existence of a diagonal complete maximum for continuous residuals on
nonempty compact spaces. Generalize finite-run continuity without restricting
the ambient type to the real line. Preserve the old finite and conditional
interfaces. Public results are re-exported and axiom-audited.

Independent checks: exact rational feature-Gram kernels on repeated samples,
every prefix of positive diagonal updates, PSD/order quadratic forms,
power-function energy inequalities, and the original-kernel fill bound.
Include a rank-one weak-PSD off-diagonal tie counterexample so diagonal
availability is not confused with necessity of every maximizing pivot.

Mathlib search found `Matrix.PosSemidef.submatrix`,
`conjTranspose_mul_mul_same`, `posSemidef_vecMulVec_self_star`, scalar
positivity and compact extreme-value theorems. No existing general kernel
Schur-run/fill-distance theorem was found. Reuse matrix positivity, not a new
quadratic-form theory. Non-claims: arbitrary off-diagonal complete-pivot
necessity for weak PSD, strict-SPD fresh-domain reduction (R08 remains open
until that is proved), asymptotic P-greedy rates, practical approximate-pivot
rates, or completion of the seven-gap audit by annotation alone.

### DJ outcome and next decision

The frozen theorem is proved: `diagonalRun_residual_and_removed_posSemidef`,
`diagonalRun_power_le_original_energy`,
`diagonalRun_error_le_original_fillDistance`,
`continuous_finalResidual_kernel`, and
`diagonalRun_exists_diagonal_completePivot` cover arbitrary ambient domains
at their stated metric/topological scopes. Finite matrix positivity is lifted
through a sample enlarged by the pivot; the removed update is a positive
rank-one Gram matrix. Summing those updates derives `K-R` PSD without a
provided Loewner-order premise. A two-point quadratic form yields the energy
comparison, then two original-kernel coordinate Lipschitz bounds give `2*L*h`.
No residual regularity, asymptotic rate, or completeness premise is smuggled in.

`KernelGeometryChecks.lean` realizes the hypotheses for `1+x*y` on the real
line and computes the actual one-step remainder `x*y`. Exact independent
tests check two arbitrary diagonal orders, all six prefixes, repeated samples,
all principal minors of the sampled residual and removed matrices, 510 energy
comparisons, fill bounds, interpolation and exact rank-five stopping. The weak
PSD off-diagonal tie and factor-one failure examples prevent stronger claims.

The full root command passes 3,630 Lean jobs, all 857 axiom reports with only
permitted dependencies, Ruff/format/mypy, 180 tests in 138.60 seconds, and no
known dependency vulnerabilities. The exact branch-head root repeat passes
180 tests in 138.11 seconds; both push and PR CI pass at
`3d7313f5884c803788a80995b0035c33a61d8b48`. PR #116 merged on 2026-10-20
at 09:09 UTC as `e21510cf677277b634d8ca90eb62083c0097ebce`. Closure is pushed
to main and exact merged/pushed-head gates are repeated before final handoff.
R09 is supported; R08 remains incomplete,
so the audit is 34/40 (85% checklist coverage), not full completion.

Next exact claim: strict positivity on arbitrary injective fresh samples must
survive a diagonal update. Mathlib's
`Matrix.PosDef.conjTranspose_mul_mul_same` accepts an injective `mulVec` map;
the rectangular elimination map with an identity fresh-node block is a
credible route. Combined with strict two-point determinant positivity, it
should rule out off-diagonal successful complete pivots. This is a next proof
contract, not a delivered theorem or an assumed invariant in DJ.

## Phase DK — strict-SPD complete-pivot reduction

Frozen contract: close audit R08 / SPEC sections 3A and 9 on arbitrary
ambient kernel types, using the actual dependent `GECP.Run`, not just a
diagonal-preferring finite trace definition. Add
`PositiveDefinite/StrictKernel.lean` and
`PositiveDefinite/ContinuousCholesky.lean`. `StrictKernelOn` requires strict
positivity of every injective finite sample in a domain. Begin from mathlib's
`Matrix.PosDef K` on the arbitrary ambient type; derive finite-sample PSD and
strictness, never assume diagonal pivots or residual strictness as inputs to
the public complete-run theorem.

Prove the finite rectangular elimination-map congruence with an identity
fresh-node block and injective `mulVec`. Use it to preserve strict sample
positivity away from a successful diagonal pivot. For a PSD kernel, strict
two-point determinants and a nonzero complete maximum force equality of the
row and column pivot. Previously selected nodes have zero residual rows;
derive freshness, propagate the strict invariant, and prove
`completeRun_diagonal_of_posDef` for every finite exact complete run on the
whole domain (a compact physical domain may be the ambient subtype).

Define `DiagonalMaxRun` independently by equality of pivot coordinates and
diagonal maximality at every actual residual update. Prove
`posDef_completeRun_iff_pivotedCholesky`, plus residual PSD/diagonal domination,
compact diagonal-maximizer availability and the original-kernel fill bound
for actual SPD complete runs without a supplied `DiagonalRun` premise.
Continuity/compactness is needed for maximum existence, not finite algebra.

Independent exact rational checks: feature-Gram SPD kernels, every actual
complete-pivot branch including diagonal ties, every prefix/fresh principal
minor, strict two-point inequalities, elimination-map congruence/injectivity,
interpolation and exact stopping. Retain the weak-PSD off-diagonal tie
counterexample. Export and axiom-audit public results, run the full root gate,
update the original-SPEC audit only for the scope actually proved, and ship
one phase via draft PR / green CI / merge / pushed-main verification.

Mathlib search identified `Matrix.PosDef.submatrix`, `.det_pos`,
`.conjTranspose_mul_mul_same` with injective `mulVec`, and the existing local
finite PSD congruence. Its Schur block interface directly exposes semidefinite
preservation, so strictness needs the injective rectangular map. Non-claims:
diagonal necessity for weak PSD, a unique diagonal maximizer, numerical
approximate-pivot guarantees, or published asymptotic P-greedy rates. The five
non-baseline audit gaps remain required work.

### DK outcome and next decision

The frozen statement is proved. `posDef_diagonalSchur_fresh` constructs the
rectangular map with a lower identity block, proves its `mulVec` injective,
and transports strict positivity by mathlib congruence. `StrictKernelOn`
propagates this over every injective fresh sample. Strict two-point
determinants and zero previously selected rows then give
`completeRun_diagonal_of_posDef` on actual finite exact complete runs, with
no diagonal preference or supplied residual strictness. The independent
`DiagonalMaxRun` recursion gives
`posDef_completeRun_iff_pivotedCholesky`; reverse conversion needs only PSD.

Public corollaries supply residual PSD, strict fresh samples, squared diagonal
domination, compact continuous complete-maximum existence and original-kernel
fill control without a provided diagonal-run proof. The concrete existing
two-pivot real kernel is proved strictly SPD by an injective Gram factor and
the general equivalence applies to that actual run. Exact rational tests
independently traverse every complete branch of feature-Gram and identity
fixtures, every prefix/fresh principal minor and strict two-point inequality,
the rectangular congruence, interpolation and full stopping. All six identity
tie orders survive; a signed rank-one weak-PSD tie refutes stronger necessity.

The full root gate passes 3,633 Lean jobs, all 876 axiom reports with only
permitted dependencies, Ruff/format/mypy, 183 tests in 136.77 seconds and no
known dependency vulnerabilities. The exact branch-head repeat passes 183
tests in 136.15 seconds; push and PR CI pass at
`c5c367400b02555a3c9ffdf33028770823e7d707`. PR #117 merged on 2026-10-20
at 09:29 UTC as `39226a2dbe2ee0760afd287a03adc745715c655c`. Closure is pushed
to main, and exact merged/pushed-head gates are repeated before final handoff.
R08's staged baseline is supported, bringing
the audit to 35/40 (87.5% checklist coverage, not weighted effort).

The remaining original-SPEC obligations are computational: R24 metadata
outside GECP, R25 census through `1e-12`, R26 complete serialized diagnostics,
R27 normalized residual curves, and R31 listed sparse comparisons. Next
freeze a complete result-metadata contract and then the canonical census
record/plot contract. Do not relabel uncertified numerical quantities as
proofs, or infer exact complete-pivot theorems for approximate maxima.

## Phase DL — complete numerical result metadata

Frozen contract: close audit R24 / numerical agent rules and SPEC sections
23–25. Add typed, JSON-safe `NumericalMetadata`, `ConditioningDiagnostic` and
`RefinementDiagnostic` records in `types.py`, shared diagnostics helpers in
`diagnostics.py`, complete metadata for every `CholeskyResult` and
`SparseResult`, and a metadata view over the already recorded GECP evidence.
Expose actual arithmetic/precision separately from requested precision,
conditioning method/status/history, deterministic tie policy/counts/tolerance,
certification status/scope, convergence/error scope and effective settings.
Use null plus an explicit status for unavailable/non-applicable quantities,
never fake finite condition numbers, certification, or high precision.

Cholesky measures original selected-core conditioning and reports float64
execution even when its shared config requests more. Sparse recovery measures
selected/final-design conditioning and separately reports nonlinear optimizer
termination, budget, acceptance and Jacobian conditioning. Keep the legacy
`validation_error` field but label it training-only and expose `training_error`.
Propagate metadata into synthetic application summaries under schema version 2;
preserve historical schema-1 committed artifacts and their hashes/provenance.

Audit failure paths while adding metadata: empty/rank-zero and budget stops,
finite input validation, exact-zero tie tolerance, exhausted sparse libraries,
rank deficiency, materially invalid PSD residuals, nonfinite arithmetic,
nonlinear optimizer budget/failure and rejected worsening refinement. Final
`converged` must describe the returned approximation's actual sampled error,
not a stale pre-refinement flag. Do not select an already exhausted atom or
selected zero Cholesky pivot. Convergence is sampled numerical evidence, not
an analytic or rounding certificate. No explicit numerical inverse is allowed.

Test metadata against independently recomputed residuals/SVDs, serialize with
`allow_nan=False`, repeat deterministic runs, exercise all relevant failure/
not-applicable paths and retain existing successful fixtures. Root verification,
green CI, draft PR, merge and pushed-main gates are required. Non-claims:
arbitrary precision Cholesky/sparse execution, continuous certification for
sampled results, new Lean results, minimum sparse support or the remaining
census/comparison deliverables. No Lean research-loop skill is needed because
this phase adds no open formal theorem.

DL type-audit clarification: fixed-order interpolation and dyadic Taylor
construction objects are also public approximation outputs. Give them explicit
non-iterative/not-applicable metadata, distinguishing internal mpmath precision
from float64 returned evaluation arrays. Do not invent convergence tests or
uniform numerical certificates for these fixed constructions.

DL verification found newly disclosed dependency advisories on 2026-10-22:
the first full gate passed 232 tests but pip-audit reported 25 findings in
locked certifi 2026.7.22 and pillow 12.3.0, with fixed versions 2026.10.19 and
13.0.0 respectively. Update those transitive lock entries and repeat root/CI
gates; do not waive the audit or misreport the first full gate as successful.

### DL outcome and next decision

The metadata contract is implemented for GECP, Cholesky, sparse recovery,
fixed constructions and schema-2 application summaries. Conditioning estimates
are scoped by stage, unavailable values are null with a status, and optimizer
termination is not confused with returned-fit convergence. Independent review
and 49 new test cases cover singular values, actual returned errors, tie
thresholds, rank-zero/non-applicable cases, exhaustion, invalid inputs,
optimizer budget/failure/worsening/nonfinite rejection and deterministic JSON.
Cholesky's actual error and accepted near-symmetry floor are exposed; stable
symmetrization also handles finite SPD input at scale `1e308` without an
avoidable overflow. Fixed barycentric construction now uses a fixed random
state. No numerical inverse or new formal theorem is introduced.

The transitive lock updates are exactly certifi 2026.10.19 and pillow 13.0.0;
pip-audit reports no known vulnerabilities. A PNG render/read smoke passes.
The repeated full root command passes 3,633 Lean jobs, all 876 axiom reports
with permitted dependencies only, Ruff/format/mypy (14 source files), 232
tests in 137.28 seconds, and the clean dependency audit. Rebuilt local source
and wheel archives include the new diagnostics module and exclude the local
notes/caches. An isolated Python 3.12 wheel smoke outside the checkout imports
all 45 exports and executes Cholesky with metadata; nothing is published.

The final branch-head root repeat passes 232 tests in 135.97 seconds; push
and PR CI pass at `7e3cffede223a0016dc24b14b1bdbeaf298c6b29`. The final wheel
smoke outside the checkout also checks convergence/stop snapshots. PR #118
merged on 2026-10-22 at 09:53 UTC as
`1e95e5fad62cd7e13154e13ed82ee4b32dac7010`. Closure is pushed to main and
exact merged/pushed-head gates are repeated before final handoff.
R24 is supported; R25–R27 and R31 remain incomplete, so the audit is
36/40 (90% checklist coverage, not weighted effort). Next freeze a complete
census record contract: seven cutoffs through `1e-12`, all residual/core/
geometry histories and normalized `k/log(1+Lambda)` curves at deterministic
revision/config/precision. Sparse-method comparisons then remain independent
deliverables before a final original-SPEC audit can claim completion.

## Phase DM — complete census diagnostics and normalized curves

Frozen contract (SPEC sections 14/26, audit R25–R27): extend the seven-cutoff
canonical 128-bit grid census to `1e-6,1e-8,1e-10,1e-12` (28 records).
Add schema-2 `CensusRecord` output with full residual history (including the
terminal value), signed original-core determinant history, smallest singular
values, signed pivots and pivot magnitudes, per-pivot distances to all four
physical boundaries, explicit signed dyadic panels, reflected coordinates,
nearest selected reflection partner/distance and a stated matching tolerance.
Derive prefix metadata/settings from the actual strict trajectory without
claiming a prefix's convergence merely from a failed/uncertified certificate.
All high-precision quantities remain canonical decimal strings.

Record code revision, dirty state, effective configuration hash, grid sizes,
runtime/dependency versions and lock digest. Canonical config requires clean
tracked code before generation; exploratory dirty runs are labeled, not
silently attributed to a commit. Preserve legacy schema-1 JSONL/PNG artifacts
and hashes; write newly named schema-2 artifacts. Validate finite controls,
history alignment, exact-zero/budget/failure paths, deterministic ordering,
atomic output and byte-identical repeated small configurations. A second
full canonical run must match bytes at the same revision/configuration.

Extend `plot_census.py` with normalized residual curves using
`k/log(1+Lambda)`, one deepest trajectory per cutoff, accurately labeled as
sampled-grid numerical evidence. Keep rank-versus-cutoff plotting, disclose
numerical zero treatment and inspect rendered output. Independent exact
geometry/prefix tests and higher-precision sampled-core determinant/singular-
value checks supplement existing Lean identities. Add a denser-grid check at
representative cutoffs spanning `1…1e6`; record its limited scope rather than
claiming continuous certification or grid independence.

Numerical runtime may be substantial: resume/wait on confirmed live handles,
never restart merely after an observation timeout. No new Lean result or
near-optimal theorem for approximate numerical pivots is claimed. The remaining
Phase-II sparse comparisons (R31) and final original-SPEC audit stay required.
Use the phase skill for verification/delivery; the Lean research loop is not
needed because this phase introduces no open formal claim.

### DM implementation checks

The immutable numerical-generator/configuration snapshot is
`c0e8c5093134f369b93781c413750b7f1de46036`. Two complete canonical runs were
started from its clean tracked tree. Subsequent review repaired non-integer
outer frequency endpoints without changing any canonical integer-cutoff grid
byte, rejected adaptive/unsuccessful inputs in finite-grid plots, and rejected
negative schema-2 residual maxima in block analysis. These scope repairs do
not relabel the dataset's actual generation revision as a later commit.

The complementary held-out check is a fixed-selected-cross sample:
up to three source-disjoint order-32 interior coordinates in every endpoint/
dyadic panel per axis, Cartesian product, evaluated at 192 bits using one LU
factorization and batched triangular solves. It is not a denser GECP rerun,
full-denser-grid maximum, continuous certificate, or proof of grid independence.
Independent 192-bit core determinant, signed pivot-product and singular-value
checks use fixed acceptance thresholds; disagreements must be preserved and
reported rather than fixed by adapting thresholds after seeing the data.

Independent review requires the resolution component to include actual GECP
trajectories, not just held-out selected crosses. Add clean-revision
`experiments/census_resolution.toml`: order 32, 128 bits, tolerance `1e-12`,
cutoffs `1,1000,1000000`. Compare these three actual denser runs with the
canonical order-24 trajectories, including ranks, residuals, pivot-path/tie
changes and unchanged controls. This is a representative two-resolution
ladder, not an exhaustive all-cutoff sweep, residual-adaptive grid refinement,
grid independence or a continuous certificate. The existing interval-adaptive
mode remains separately tested. Retain held-out samples as complementary
evidence, not a replacement for the actual denser runs.

Next-target review also found that the historical synthetic fixture named
`dense_blind_scan` unions its dense candidate dictionary with the generating
frequencies. Preserve that historical artifact and measured function errors,
but correct the prose: it is support-augmented, not support-blind. R31's new
comparisons must separate known libraries from independent fixed dictionaries
and off-grid truth; no blind-recovery claim follows from the old fixture.

Promising next-phase bridges (to test independently, not new formal claims):
uniform samples have exponential nodes `z_j=exp(-h*omega_j)` and amplitudes
`a_j=g_j/(1+exp(-omega_j))`, so a Hankel/pencil Prony fit can yield candidate
frequencies followed by a solve in the original kernel design. At fermionic
Matsubara `z=i*pi*(2n+1)`, direct integration gives
`-integral_0^1 exp(z*t)*G(t) dt = sum_j g_j/(z-omega_j)` for finite atoms.
AAA or true vector-fitting pole relocation should use that rational bridge,
not mistake poles of a rational interpolant in time for spectral frequencies.
All six SPEC methods should share training observations/noise and held-out
time-domain errors; transform error, invalid poles, optimizer budgets and raw
versus refined outcomes remain visible. Native-frequency-data cohorts must be
separate from same-time-input comparisons. No global sparse optimality or
unique support recovery is inferred from the smallest tested passing fit.

### DM numerical results

The 28-record canonical schema-2 JSONL is generated from clean
`c0e8c5093134f369b93781c413750b7f1de46036`, configuration hash
`c0ac1bbaa7ddb39a2bb4a1c306985386f89f4103ff75cb4e21f99257e1e9f31b`,
with SHA-256
`c6854330fcb74007f8f9814fba48c314d34fe26c1cef5461d908ea5182b66bda`.
Both full runs match byte for byte; observed wall times were 1496.40 and
1508.07 seconds while running concurrently. All cases converge on their
sampled grids. For cutoffs `1,10,100,1000,10000,100000,1000000`, strictest
`1e-12` ranks are `9,17,37,63,91,119,146`. All 21 shared legacy tolerance
cases retain the exact previous signed pivot strings and ranks; the legacy
JSONL/PNG hashes are unchanged.

The rank and normalized plots render correctly and repeat byte-identically.
No numerical zero samples occur. Large-cutoff normalized trajectories lie on
a similar scale, but do not coincide; the cutoff-one curve is visibly
different. No common fitted law or new rate is inferred from that picture.

The preliminary 192-bit primary-core analysis passes all 21 selected-core
checks: maximum relative discrepancies are below `4.45e-27` for determinants,
`6.33e-27` for signed pivot products, and `8.70e-27` for smallest singular
values. Fixed acceptance thresholds remain `1e-18`, `1e-18`, and `1e-15`.
At representative cutoffs `1,1000,1000000`, the held-out Cartesian samples
contain `36,3960,15120` points; maximum observed errors are approximately
`6.31e-14,6.33e-13,8.61e-13`, respectively. These are point-arithmetic
observations, not full denser-grid maxima or continuous bounds.

The actual order-32 GECP ladder also converges at all three representative
cutoffs, with ranks `9,63,146`, matching order 24 in these cases. Grid shapes
are `63*63`, `621*683`, `1279*1303`; residual maxima are approximately
`1.784e-13,9.597e-13,9.717e-13`. Equal ranks are an observation, not a
grid-independent-path theorem. The denser dataset is generated from clean
`ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61`, configuration hash
`bf49e21ce208a1336bfc83d3e2ec955adf627e64e8056c8be9986d855ee2d3be`,
SHA-256
`7251cbb528179cc5cbfa6eb5635b8fa280e087d1e04db1348b56a9065026d588`,
in 2004.01 wall seconds. Core GECP, kernel, census and type source files are
identical between the two generation revisions; the intervening grid endpoint
repair preserves every canonical integer-cutoff grid byte. The combined
`gecp_census_checks.json` passes all 30 selected-core comparisons at 192 bits,
with the same discrepancy maxima as the primary checks above. Primary-core
and held-out results exactly repeat the preliminary computation. The combined
report SHA-256 is
`b7192e8bc57a72bf64802dcb9334d4c1a4761be5a2a2876368fd72aeda40247d`.
Denser-grid residual maxima are larger in all three cases and the common
selected-node prefix lengths are `3,0,0`: equal ranks do not mean equal paths.
All five committed-artifact audit tests pass. Full root and exact-head push/PR
CI pass; PR #119 is merged. Exact pushed-main closure gates repeat before
handoff. R31 and the final full-SPEC audit remain next.

### DM full local gate and CI finalization recovery

The full root command at `91f9b4a62c9440546d8b363da200f0c88f9bab39`
passes 3,633 Lean jobs, all 876 permitted-axiom reports, Ruff/format (69 files),
mypy (14 sources), 430 tests in 139.13 seconds, and the clean dependency audit.
Final source/wheel archives include every new numerical artifact and exclude
local notes/caches. An isolated installed-wheel smoke outside the checkout
loads all 45 exports and exercises MP128 census diagnostics and non-integer
physical endpoints. No package is published.

PR CI run `39814069068` passed on that exact head. Push CI run `39814062885`
had all successful job steps, including Python's `Complete job` at 16:29:58 UTC,
but its job/run summary remained in progress. GitHub rejected cancellation as
already completed and rejected a targeted rerun. This occurred while the
[official status page](https://www.githubstatus.com/) reported degraded
Actions/Pull Requests service. No test or audit was waived, and no numerical
process was restarted. The actual verification/provenance record was committed
as `4f99db923e89773e1d248eecde872f2350e0d853` to obtain fresh exact-head CI.

### DM outcome and next decision

At that final branch head, the full root repeat passes 430 tests in 139.64
seconds, the same 3,633 Lean jobs/876 axiom reports and clean dependency audit.
Push CI `39815369063` and PR CI `39815374660` both pass. PR #119 merged on
2026-10-22 at 16:44:20 UTC as
`b068e1070745d1308ee23101096e1197c345f057`. Main's pushed closure head repeats
the root and exact-head CI gates before final handoff; no artifact is published.

R25–R27 are supported. The audit is now 39/40 (97.5% checklist coverage, not
weighted effort). Next freeze R31's fair six-method comparison contract using
shared observations and held-out time-domain errors, actual spectral-atom
bridges and explicit failure/conditioning/precision metadata. Keep known
libraries separate from support-blind cases and retain failed fits. Only after
that implementation and a final original-SPEC audit can the full goal be
considered complete.
