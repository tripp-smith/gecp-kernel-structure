# Final Handoff Report

## Project summary

`gecp-kernel-structure` v1.0.0 delivers a `sorry`-free Lean 4 library and the
`kernelgecp` Python research package. The strongest fermionic result is a
machine-checked separated approximation with rank
`16p(s+1)` and uniform error `2⁻ᵖ` on
`[0,1] × [-2ˢ,2ˢ]`, which is the required
`O(log Λ log(1/ε))` low-rank scale.

Milestone E closes through the SPEC-approved obstruction outcome: exact strict
minor signs do not determine complete-pivot locations, while the refined
`CrossRatioControl` hypothesis is proved sufficient for geometric residual
decay. The repository does not claim that the fermionic residuals satisfy that
hypothesis with a cutoff-uniform contraction factor.

The post-v1 Python application phase demonstrates how those components work
together on realistic synthetic data: one universal fermionic basis compresses
two qualitatively different continuous spectra, sparse routines represent
clean and noisy Green functions under explicit dictionary assumptions, and
the PSD baseline selects spatial covariance landmarks. These are application
observations with held-out errors, not new formal convergence or inverse-
problem claims.

Post-v1 Phase G adds a rigorous rate obstruction and the first strict
fermionic block theorem. Lean proves that the first two complete pivots are
the two reflected cutoff corners for every positive cutoff, that no universal
factor below one can contract the first step, and that for `0<Λ≤1` the actual
two-step residual is at most one half of the initial complete-pivot magnitude.
This verifies the base block but does not supply the dyadic scaling lemma
needed for Conjecture G1.

Post-v1 Phase BL closes the qualitative and pointwise parts of the amplified-
determinant strategy. Lean proves that nonnegative mixed coefficients make the
tail-weight pencil determinant-monotone, that every amplification weight at
least one preserves strict total positivity, and that the `8!`-amplified
positive-series kernel is bounded by two on the unit band. Complete runs of
that amplified kernel therefore have a `2^n` selected-core determinant bound.
The remaining rate gap is the arbitrary ordered-sample determinant envelope:
the original GECP nodes are not assumed to be complete pivots of the amplified
kernel.

Post-v1 Phase BM removes that pivot-path restriction. A generic northwest-
Schur induction proves that every strictly totally positive finite matrix has
determinant at most its diagonal product, so an entry bound applies to every
ordered sample, not only a complete run of the bounded kernel. Combining the
arbitrary-sample `2^n` amplified envelope with the Phase BK rank-eight
amplification gives
`det K_n <= 2^n / (8!)^(n-8)` for positive unit-band samples and `n >= 8`.
The remaining rate gap is now composition with the original GECP selected-
core/geometric-mean theorem and assembly across reflections and dyadic bands.

Post-v1 Phase BN completes the first exact GECP composition of that decay.
Independent sorting removes the mismatch between monotone determinant
theorems and algorithmic pivot order, after which the selected-core
geometric-mean inequality gives a residual-power bound. Exact arithmetic at
rank ten proves a uniform one-half contraction for complete pivoting on the
reflected positive-series kernel over the positive unit band. The remaining
global step is to assemble positive and negative frequency columns, close the
zero endpoints, and transport the unit-band block across dyadic scales.

Post-v1 Phase BO closes the full signed-unit-band and endpoint obligations in
one step. The affine coordinate `z = 1 - omega` rewrites the original kernel
as positive diagonal weights times `exp(t*z)` on `[0,1] x [0,2]`, avoiding a
mixed-column Laplace expansion. A weight-`128` eighth-order tail amplification,
strict total positivity, and the diagonal-product envelope prove
`|det K_n| <= 2^n / 128^(n-8)` for arbitrary injective samples. Composing with
the exact selected-core residual-power theorem proves a uniform half
contraction after twelve complete pivots on `[0,1] x [-1,1]`. The remaining
global step is dyadic rescaling/continuation across larger cutoffs.

Post-v1 Phase BP removes absolute band location from that local theorem. For
every real center, the coordinate `z = center + 1 - omega` maps the closed
width-two band to `[0,2]`. The weight-`128` amplified kernel remains strictly
totally positive and bounded by two uniformly in the center, so arbitrary
injective samples satisfy the same `2^n / 128^(n-8)` determinant envelope and
twelve exact complete pivots give the same one-half residual contraction. The
remaining global step is now specifically a conditioned-residual restart
across the selected outer-scale prefix, not an affine recentering problem.

Post-v1 Phase BQ identifies the exact restart interface. It disproves
one-sided frozen-column amplification by a kernel-checked rational
positive-Taylor witness, then proves that an exact low-rank-plus-tail split of
the actual final Schur residual retains the full fresh-order amplification
exponent. Phase AJ's bordered factorization composes this with the selected
core without charging prefix length against that exponent. What remains is no
longer an algebraic composition ambiguity: it is the analytic construction of
a compatible split, with nonnegative mixed coefficients, for arbitrary
fermionic final residuals.

Post-v1 Phase BR constructs that compatible split exactly. The pivot
coefficients of a fixed run define a linear companion row-conditioning
operator which commutes with a fixed right feature factor. Consequently the
positive-series fermionic decomposition remains a conditioned tail plus a
rank-eight block after every prefix, with no prefix-length penalty in the
Phase BQ exponent. Exact arithmetic also shows that conditioning does not
preserve coefficientwise nonnegativity: a mixed coefficient becomes
negative. The remaining analytic target is therefore a signed or real-rooted
conditioned determinant-pencil comparison, not the previously anticipated
all-nonnegative coefficient theorem.

Post-v1 Phase BS removes that failed sign theorem from the critical path.
Rank-eight vanishing leaves a polynomial of degree at most eight, so exact
Lagrange interpolation at nine weights between `56` and `128` controls its
value at one without assumptions on coefficients or roots. The formal
order-eighteen corollary converts nine `2^18` companion determinant envelopes
into the `2^-18` physical determinant bound required for a half-contraction
block. The remaining restart problem is thus a finite constant-base envelope
for those companion residuals, not a coefficientwise positivity theorem.

Post-v1 Phase BT gives that finite-envelope problem its natural two-sided
coordinates. Each fixed pivot applies a left and right oblique projector to a
companion, annihilating both selected directions while preserving linearity.
The projector transports `L*R` by conditioning the two factors separately, so
every fermionic prefix retains an exact two-sided rank-eight split. The
conditioned low determinant is now an ordered Cauchy--Binet sum of matched
left/right feature minors. Exact finite searches consistently recover a
common coefficient orientation, but the corresponding all-orders
matched-minor sign theorem remains open.

Post-v1 Phase BU removes the conditioning from that sign theorem. Exact left
and right feature-augmented determinant identities express each conditioned
feature minor as a Schur quotient by the selected core. In a matched product
the core is squared, so pivot order and core sign cannot affect the result.
For the fermionic rank-eight block, every conditioned Cauchy--Binet term now
has the sign of two original bordered feature determinants. Those borders are
confluent exponential systems: the polynomial features are derivatives of
the exponential family at parameter zero. Proving their oriented signs is the
next analytic step.

Post-v1 Phase BV corrects that proposed next step. An exact three-pivot
finite-Taylor witness has two negative products among its 56 order-three
bordered-feature coordinates, so coordinatewise confluent positivity is
false. The complete signed sum nevertheless has a canonical invariant form:
adjoin the selected kernel blocks `U,V` to the two feature factors `A,B`, mark
the fixed coordinates, and take their full-degree block coefficient. Lean
proves this is exactly the mixed-column coefficient of
`det(mu • (U*V) + A*B)`. In the witness the coefficient equals the signed sum
and remains positive. The open analytic problem is positivity of this
aggregate graded pairing for arbitrary fermionic prefixes, not positivity of
each Plücker coordinate. The root command passes a 3,571-job Lean build with
107 tests and no known dependency vulnerabilities; the change passed CI and
merged as PR #76.

Post-v1 Phase H restores verification integrity after the local editable
pytest entry point became stale and the development lock acquired dependency
audit findings. The root and CI workflows now launch pytest through the
synchronized interpreter, the audited transitive lock resolves `urllib3`
2.8.0, and the complete root command passes in maintained and fresh frozen
environments. This is maintenance evidence and changes no mathematical claim.

Post-v1 Phase I isolates the next exact analytic obligation. Lean proves a
bordered selected-core determinant identity and uses core nonsingularity to
show that selected-residual cross-product sign coherence is equivalent to a
four-bordered-minor sign condition on the original kernel. Exact surrogate
checks retain first-failure witnesses. This is an algebraic reduction, not a
proof that the fermionic kernel satisfies the determinant-sign condition.

Post-v1 Phase J closes the algebraic part of that obligation. Lean formalizes
ordered strict sign regularity, transports it through arbitrary tuple
orientations and positive column scaling, and proves that an all-orders
exponential sign theorem would force selected-residual sign coherence for
every fermionic GECP run. Orders one and two are proved unconditionally. The
all-orders exponential theorem itself remains the precise analytic gap.

Post-v1 Phase K begins the constructive all-orders proof. Lean proves positive
row scaling, strict positivity of ordered ordinary Vandermonde determinants,
and strict positivity of the principal square exponential Taylor feature
block at every order. This certifies the principal Cauchy--Binet term while
leaving generalized Vandermonde positivity and the infinite-series limit open.

Post-v1 Phase L proves the positive-node generalized Vandermonde theorem at
every finite order. Lean derives a sparse-polynomial positive-root bound from
Descartes' rule, uses it to prove power-matrix nonsingularity, anchors the sign
at geometric nodes, and transports positivity through a zero-avoiding
homotopy. The zero-node boundary case and exponential-series determinant limit
remain open.

Post-v1 Phase M closes the zero-node boundary: nonnegative strictly increasing
nodes and strictly increasing natural exponents give a nonnegative generalized
Vandermonde determinant. This supplies the missing termwise sign input for a
finite exponential Cauchy--Binet expansion.

Post-v1 Phase N assembles that finite argument. A reusable rectangular
Cauchy--Binet identity groups injective intermediate tuples into increasing
tuples and permutations. Applied to finite exponential Taylor features, every
minor product is nonnegative and the principal product is strictly positive,
so every truncation with at least as many terms as the matrix order has
positive determinant. Passing to the infinite-series determinant while
retaining strictness remains open.

Post-v1 Phase O passes to that limit without losing strictness. Each finite
Taylor determinant dominates the same positive principal determinant; the
feature matrices converge entrywise to `exp(row_i * column_j)`, and
determinant continuity preserves the lower bound. Thus the full exponential
matrix has positive determinant on nonnegative strictly increasing node
tuples. Arbitrary-node scaling and the frequency-reversal orientation remain
before the all-orders sign-regularity theorem.

Post-v1 Phase P completes the all-orders sign theorem. Positive exponential
row and column scalings remove arbitrary node offsets, finite reversal has
sign `(-1)^(n.choose 2)`, and negated reversed frequencies convert positive
`exp(x*y)` determinants into the prescribed `exp(-t*omega)` signature. The
fermionic kernel inherits strict sign regularity at every order, so every
selected residual cross is sign coherent and exact complete-pivot updates are
nonexpansive. This still does not supply strict dyadic block contraction.

## Implementation run metadata

This is a provenance record for the Codex implementation goal, captured at
handoff on 2026-08-15 at 17:39 EDT (UTC-04:00).

| Field | Recorded value |
| --- | --- |
| Codex surface | Codex desktop app |
| Declared model identity | OpenAI Codex based on the GPT-5 family |
| Exact serving model/version | Not exposed to the task runtime; no model slug or snapshot can be verified |
| Collaboration mode | Default mode; primary/root agent |
| Reasoning effort | Not exposed to the task runtime |
| Workspace execution mode | Unrestricted local filesystem/network access; approval policy `never` |
| Goal state | Complete |
| Goal-accounting usage at completion | 1,801,132 aggregate tokens and 10,747 seconds (2:59:07) |
| Usage scope | The implementation goal through v1 completion; this post-completion metadata update is not included because its goal counter is no longer available |

The aggregate token counter does not split input, cached input, reasoning, or
output tokens and is not an OpenAI invoice. The runtime also does not expose
the exact serving model. For an order-of-magnitude comparison only, the
current [OpenAI API pricing](https://openai.com/api/pricing/) lists standard
GPT-5.6 Sol text rates of $5 per million input tokens and $30 per million
output tokens. Applying those rates to the 1.801132-million-token aggregate
gives:

- an illustrative 80% input / 20% output estimate of **$18.01**;
- an uncached standard-rate sensitivity envelope of **$9.01–$54.03** if the
  whole aggregate were priced at the input or output rate, respectively.

These figures are API-equivalent estimates, not the actual Codex subscription,
credit, or internal compute cost. They exclude cache discounts, long-context
premiums, service-tier adjustments, tool-call charges, and any usage that was
not represented by the completed goal counter.

## Phase status

| Phase | State | Delivery |
| --- | --- | --- |
| 0 — bootstrap | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| A — exact GECP core | complete | [PR #4](https://github.com/tripp-smith/gecp-kernel-structure/pull/4) |
| B — PSD baseline | complete | [PR #5](https://github.com/tripp-smith/gecp-kernel-structure/pull/5) |
| C — fermionic structure/census | complete | [PR #6](https://github.com/tripp-smith/gecp-kernel-structure/pull/6) |
| D — separated approximation | complete | [PR #7](https://github.com/tripp-smith/gecp-kernel-structure/pull/7) |
| E — structural GECP | complete | [PR #8](https://github.com/tripp-smith/gecp-kernel-structure/pull/8) |
| F — sparse spectral application | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| R — release readiness | complete | [PR #9](https://github.com/tripp-smith/gecp-kernel-structure/pull/9), [closure PR #10](https://github.com/tripp-smith/gecp-kernel-structure/pull/10) |
| S — realistic synthetic applications | complete | implementation `6404311`; canonical evidence delivered directly to `main` |
| G — fermionic GECP rate | complete bounded research phase; G1 open | merged [PR #11](https://github.com/tripp-smith/gecp-kernel-structure/pull/11) |
| BC — product-measure Andréief identity | complete; green CI | merged [PR #57](https://github.com/tripp-smith/gecp-kernel-structure/pull/57) |
| BD — arbitrary-sample and Taylor-tail weak signs | complete; green CI | merged [PR #58](https://github.com/tripp-smith/gecp-kernel-structure/pull/58) |
| BE — strict all-orders Taylor-tail signs | complete; green CI | merged [PR #59](https://github.com/tripp-smith/gecp-kernel-structure/pull/59) |
| BF — factorial-free Taylor-tail GECP envelope | complete; green CI | merged [PR #60](https://github.com/tripp-smith/gecp-kernel-structure/pull/60) |
| BG — positive-series compatible fermionic split | complete; green CI | merged [PR #61](https://github.com/tripp-smith/gecp-kernel-structure/pull/61) |
| BH — block-feature coefficientwise Cauchy--Binet | complete; green CI | merged [PR #62](https://github.com/tripp-smith/gecp-kernel-structure/pull/62) |
| BI — mixed-column coefficient bridge | complete; green CI | merged [PR #63](https://github.com/tripp-smith/gecp-kernel-structure/pull/63) |
| BJ — infinite-tail mixed coefficients | complete; green CI | merged [PR #64](https://github.com/tripp-smith/gecp-kernel-structure/pull/64) |
| BK — low-rank tail amplification | complete; green CI | merged [PR #65](https://github.com/tripp-smith/gecp-kernel-structure/pull/65) |
| BL — amplified-kernel envelope | complete; green CI | merged [PR #66](https://github.com/tripp-smith/gecp-kernel-structure/pull/66) |
| BM — total-positive Hadamard envelope | complete; green CI | merged [PR #67](https://github.com/tripp-smith/gecp-kernel-structure/pull/67) |
| BN — positive-series GECP composition | complete; green CI | merged [PR #68](https://github.com/tripp-smith/gecp-kernel-structure/pull/68) |
| BO — signed unit-band amplification | complete; green CI | merged [PR #69](https://github.com/tripp-smith/gecp-kernel-structure/pull/69) |
| BP — translated unit-band amplification | complete; green CI | merged [PR #70](https://github.com/tripp-smith/gecp-kernel-structure/pull/70) |
| H — verification integrity refresh | verified locally; green CI | merged [PR #12](https://github.com/tripp-smith/gecp-kernel-structure/pull/12) |
| I — determinantal sign-coherence bridge | complete; green CI | merged [PR #13](https://github.com/tripp-smith/gecp-kernel-structure/pull/13) |
| J — strict-sign-regularity transfer | complete; green CI | merged [PR #14](https://github.com/tripp-smith/gecp-kernel-structure/pull/14) |
| K — exponential total-positivity scaffolding | complete; green CI | merged [PR #15](https://github.com/tripp-smith/gecp-kernel-structure/pull/15) |

## Completed tasks

- T-001–T-003: Codex-native skills, pinned toolchains, repository skeleton,
  CI, axiom audit, and root verification command.
- T-004–T-005: exact GECP residual/run core plus finite-grid and
  arbitrary-precision numerical GECP.
- T-006: positive-semidefinite residual preservation and canonical recursive
  GECP/pivoted-Cholesky equivalence.
- T-007–T-008: fermionic identities, derivatives, endpoint-resolved grids,
  interval pivot certification, and deterministic census.
- T-009–T-010: explicit exponential-family and fermionic separated
  approximation theorems with an independent high-precision evaluator.
- T-011–T-012: exact structural census, minimized sign-regular obstruction,
  and cross-ratio sufficient condition.
- T-013–T-015: discrete-to-continuous certificate bridge, Green-function
  error transfer, and sparse two-atom application.
- T-016: v1.0.0 metadata, release smoke configuration, artifacts, status
  audit, and final handoff.
- T-017: typed synthetic-application API, Hubbard-like and gapped Green-
  function compression, clean/noisy sparse recovery, PSD covariance landmark
  selection, deterministic JSON, summary plot, and application documentation.
- T-018: variable-factor/block contraction infrastructure, exact one-step
  obstruction, exact first-two-pivot localization for every positive cutoff,
  rank-two secant/interpolation analysis, and the proved two-step half
  contraction on `0<Λ≤1`.
- T-019: interpreter-bound local/CI pytest entry point, refreshed audited
  development lock, and maintained-plus-fresh-environment root verification.
- T-020: bordered selected-core determinant recurrence and residual identity,
  equivalent four-minor sign criterion, exact selected-cross diagnostics, and
  a generic minimized failure witness.
- T-021: ordered-minor strict-sign-regularity definitions, orientation-aware
  transfer to bordered-minor and selected-residual coherence, positive column
  scaling, fermionic reduction, order-one/order-two proofs, and exact
  permutation regressions.
- T-022: positive row scaling, ordered Vandermonde determinant positivity,
  positive determinant of the principal exponential Taylor block, and an exact
  rational factorization regression.
- T-023: sparse-polynomial positive-root bounds, generalized Vandermonde
  nonsingularity, and positive-node determinant positivity.
- T-024: generalized Vandermonde nonnegativity at a zero boundary node by
  uniform positive shifts and determinant continuity.
- T-025: rectangular Cauchy--Binet expansion, finite exponential Taylor
  determinant positivity, and exact rational multi-term regressions.
- T-026: fixed principal determinant lower bound, entrywise exponential-series
  convergence, determinant-limit positivity, and a 100-decimal regression.
- T-027: arbitrary-node exponential determinant positivity, finite-reversal
  sign, all-orders exponential and fermionic strict sign regularity, and
  unconditional selected-residual sign coherence.

## Deferred research

No v1 implementation task is blocked. Open research and optional extensions
are:

- transport the proved `Λ≤1` two-corner half contraction through dyadic
  frequency layers, using the now-unconditional selected-pivot
  nonexpansiveness invariant between strict contractions, to obtain the
  `O(log Λ)` block length required by Conjecture G1;
- formalize the Gimbutas–Marshall–Rokhlin selected-exponential Chebyshev
  construction as an alternate low-rank theorem;
- formalize stronger published pivoted-Cholesky Lipschitz rates;
- add AAA, vector fitting, and multidimensional sparse recovery in a later
  version.

## Public Lean results

The public theorem surface includes:

- exact GECP: `gecp_interpolates_selected_rows`,
  `gecp_interpolates_selected_cols`, `gecp_core_nonsingular`, and
  `gecp_core_det_eq_prod_pivots`;
- PSD baseline: `posDef_gecp_residual_posSemidefinite`,
  `posDef_gecp_residual_abs_le_diag`,
  `posDef_gecp_residual_sup_eq_diag_sup`,
  `gecp_eq_pivotedCholesky_of_posDef`, and
  `posDef_gecp_error_le_fillDistance`;
- fermionic structure: `fermionicKernel_reflection`,
  `fermionicKernel_centered`, `fermionicKernel_continuous`, and the exact
  derivative/bound theorems;
- low rank: `expFamily_separatedApprox` and
  `fermionicKernel_separatedApprox`;
- structural outcome:
  `signRegular_not_sufficient_for_parameterIndependent_pivots`,
  `residualUpdate_le_of_crossRatioControl`, and
  `gecp_error_le_geometric_of_crossRatioControl`;
- Phase G structure:
  `gecp_error_le_product_of_crossRatioControl`,
  `gecp_error_le_dyadic_of_crossRatioProduct`,
  `fermionicKernel_no_uniform_firstStep_contraction`,
  `fermionicKernel_firstPivot_abs_le_reflectedCorner`, and
  `fermionicKernel_twoCornerResidual_le_half_initial`;
- Phase I structure: `Run.borderedCore`,
  `Run.borderedCore_det_eq_selectedCore_det_mul_finalResidual`,
  `BorderedMinorSignCoherent`, and
  `pivotCrossProductSignCoherent_iff_borderedMinorSignCoherent`;
- Phase J structure: `StrictSignRegularAtOrder`, `StrictSignRegular`,
  `oriented_minor_pos_of_strictSignRegularAtOrder`,
  `strictSignRegular_borderedMinorSignCoherent`,
  `strictSignRegular_pivotCrossProductSignCoherent`,
  `strictSignRegularAtOrder_columnScale`,
  `expKernel_strictSignRegularAtOrder_two`, and
  `fermionicKernel_strictSignRegularAtOrder_two`;
- Phase K structure: `strictSignRegularAtOrder_rowScale`,
  `vandermonde_det_pos_of_strictMono`, `expTaylorPrincipalMatrix_apply`, and
  `expTaylorPrincipalMatrix_det_pos`;
- Phase L structure: `signVariations_lt_card_support`,
  `positiveRoots_lt_card_support`, `generalizedVandermonde_det_ne_zero`,
  `geometric_generalizedVandermonde_det_pos`, and
  `generalizedVandermonde_det_pos`;
- Phase M structure: `generalizedVandermonde_det_nonneg`;
- Phase N structure: `Matrix.det_mul_rect_eq_sum_injective`,
  `Matrix.det_mul_rect`, `expTaylorMatrix_apply`, and
  `expTaylorMatrix_det_pos`;
- Phase O structure: `expTaylorMatrix_det_ge_principal`,
  `expTaylorMatrix_tendsto_expMatrix`,
  `expTaylorMatrix_det_tendsto_expMatrix_det`, and `expMatrix_det_pos`;
- Phase P structure: `Fin.sign_revPerm`,
  `expMatrix_det_pos_of_strictMono`, `expKernel_strictSignRegular`,
  `fermionicKernel_strictSignRegular`, and
  `fermionicKernel_pivotCrossProductSignCoherent`;
- application: `grid_sup_le_max_add_lipschitz`,
  `approxPivot_of_grid_certificate`, and
  `greenError_le_kernelError_mul_l1`.

The axiom audit reports only Lean defaults (`propext`, `Quot.sound`, and where
needed `Classical.choice`).

## Python and API results

The package exports stable `FermionicKernel` evaluation, `gecp`,
`gecp_matrix`, pivoted Cholesky, solve-based cross approximation, interval
pivot certificates, deterministic census tools, sparse spectral recovery,
piecewise Chebyshev interpolation, and the Lean-matched
`fermionic_dyadic_taylor_approximation`.

The research API additionally exports exact-form decimal evaluators
`fermionic_first_residual` and `fermionic_two_corner_residual`, census block
analysis, and focused interval-certified trajectories. The two-corner
evaluator independently checks the Lean base theorem; it is not used as the
source of that theorem.

It also exports `SyntheticApplicationConfig`,
`run_synthetic_application_suite`, typed result records,
`synthetic_spectral_density`, `fermionic_green_from_density`, quadrature
weights, and deterministic sensor-point generation. The application runner
combines existing public algorithms; it does not bypass their convergence,
conditioning, or stop metadata.

Normal non-convergence is represented by `converged` and `stop_reason`.
Certified-pivot budget exhaustion is never reported as certification, and
high-precision runs retain canonical decimal strings.

## Verification commands run

```bash
./scripts/verify.sh
uv build
uv run python -m kernelgecp.census \
  --config experiments/release-smoke.toml \
  --output /tmp/gecp_release_smoke.jsonl
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/data/synthetic_applications.json \
  --plot experiments/data/synthetic_applications.png
git status --short
```

## Verification results

- Lean build and public axiom audit pass; no `sorry`, `admit`, or custom axiom
  is present.
- Ruff, formatting, strict mypy, and all 54 pytest cases pass.
- `pip-audit` reports no known vulnerabilities; the local editable project is
  correctly skipped because it is not a PyPI dependency.
- The canonical 21-case census is byte-reproducible with SHA-256
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
- The exact-surrogate dataset has SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The 128-bit release smoke run converges at tolerance `1e-6` with rank 6 for
  cutoff 1 and rank 82 for cutoff `10⁶`. Its two-record JSONL has SHA-256
  `d7b5f53d61eb2ea49cd1ebde944b47a06f931b9d82687c8fa4458ac4c09907d5`.
- The synthetic-application JSON is byte-identical across repeated runs and
  has SHA-256
  `f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2`;
  the reviewed plot has SHA-256
  `0529263affd2fe10773f71784f06007650bff833374fbbebbb16f2507b2ab169`.
- At cutoff 8 and tolerance `1e-8`, a 12-pivot universal GECP basis gives
  Green errors `4.16e-10` and `2.91e-10` for the Hubbard-like and gapped
  densities, below the common `5.40e-9` discrete transfer bound.
- The known-transition fixture is recovered with four atoms to `3.33e-16`
  held-out error. The noisy blind scan uses eight effective atoms and reaches
  `8.56e-9` held-out function error without claiming atom identification.
- GECP and pivoted Cholesky agree on 31 of 72 covariance landmarks; the
  selected cross has maximum error `7.78e-7` and relative Frobenius error
  `1.14e-7`.

## Known limitations and non-claims

- No solution of Simons Problem 4.2 is claimed.
- The low-rank theorem is not a GECP theorem.
- Census curves and pivot certificates are evidence, not a continuous GECP
  convergence proof.
- The proved `0<Λ≤1` two-step half contraction is not extrapolated to larger
  cutoffs; the dyadic block scaling required by Conjecture G1 remains open.
- Phase BR constructs the exact compatible low-rank-plus-tail split for every
  fermionic prefix and preserves the rank-eight feature factor. It does not
  prove the conditioned-tail amplification inequality: exact arithmetic shows
  that the unconditioned coefficientwise-nonnegative premise can fail after
  conditioning.
- Phase BS proves the sign-free finite-weight amplification reduction, but it
  does not yet prove the nine amplified companion determinant envelopes. Those
  envelopes are the remaining premise for the order-eighteen half contraction.
- Phase BT proves the exact two-sided companion factorization and exposes the
  matched-minor interface, but the observed common orientation of those
  conditioned minor products is still numerical evidence rather than a
  theorem.
- Phase BU proves the exact bordered-determinant sign equivalence, but it does
  not make the resulting products sign coherent; Phase BV gives an exact
  counterexample with two negative products among 56 coordinates.
- Phase BV proves the aggregate graded-pairing identity, but the observed
  positivity of that mixed-column coefficient is not yet a theorem for
  arbitrary fermionic prefixes.
- All-orders exponential and fermionic minor signs, and selected-residual
  nonexpansiveness, do not imply a strict contraction factor below one.
- Full exponential determinant positivity is currently stated for
  nonnegative strictly increasing node tuples. The translation/scaling and
  frequency-reversal orientation needed for arbitrary ordered tuples and
  `expKernelSignature` are not yet assembled.
- Arbitrary PSD tie-breaking is not claimed equivalent to diagonal pivoting;
  the formal theorem uses the documented diagonal-preferring canonical rule.
- The formal separated construction is dyadic truncated Taylor, not the
  selected-exponential Chebyshev variant.
- The post-v1 spectral fixtures are synthetic analogues, not material-
  specific calculations, experimental fits, or uncertainty-quantified
  analytic continuation.
- Accurate held-out Green-function recovery in the noisy dense scan does not
  identify the four generating atoms uniquely; the eight recovered atoms are
  an effective representation.
- Matching covariance pivots do not prove optimal sensor placement.

## Numerical stability and performance

- Kernel evaluation uses sign-dependent overflow-safe formulas.
- GECP supports float64, exact rational surrogates, and arbitrary precision;
  the canonical census uses 128-bit arithmetic.
- Numerical cross approximation uses solves/factorizations rather than an
  explicit inverse.
- Every GECP result records pivots, residual history, conditioning, ties,
  convergence, and stopping metadata.
- The 21-case census converges in every tracked case; sampled rank at
  tolerance `1e-10` grows from 8 at cutoff 1 to 124 at cutoff `10⁶`.

## Security and dependency notes

The Lean revision and Python dependency graph are pinned. Runtime algorithms
need no network access, and the repository stores no credentials, personal
data, or telemetry. Generated environments, caches, and package builds remain
ignored.

## Files and artifacts delivered

- Lean library: `GECPKernelStructure`.
- Python package: `kernelgecp` v1.0.0.
- Canonical data and plot under `experiments/data/`.
- Synthetic application configuration, canonical JSON, and reviewed summary
  plot under `experiments/` and `experiments/data/`.
- Release artifacts: `dist/kernelgecp-1.0.0-py3-none-any.whl` and
  `dist/kernelgecp-1.0.0.tar.gz`.
- Claim and research documentation: `README.md`, `FINDINGS.md`, `RESEARCH.md`,
  `APPLICATION.md`, `MATHLIB.md`, and `SPEC.md`.
- Preserved workflow input: `autonomous-implementation.md`, explicitly labeled
  as a reference rather than current status and excluded from package builds.

## How to run locally

```bash
uv sync --all-extras
lake exe cache get
uv run python -c "import kernelgecp; print(kernelgecp.__version__)"
uv run python experiments/synthetic_applications.py
```

See `APPLICATION.md` for GECP, certified-pivot, low-rank, and sparse-recovery
examples.

## How to test

```bash
./scripts/verify.sh
```

Focused commands are documented in `AGENTS.md`; a phase is considered verified
only after the root command passes.

## Recommended next research step

Prove positivity of the single full-fixed-degree coefficient isolated by
Phase BV. Work at the compound-matrix level: seek either a sign-reversing
involution that pairs the exceptional Plücker coordinates with dominating
positive coordinates, or a stable/Lorentzian determinant-pencil theorem whose
coefficient inequalities survive the Schur prefix. The exact 2-of-56 witness
is the regression boundary for any candidate argument. Once this aggregate
coefficient is nonnegative at the nine Phase BS weights
`56, 65, ..., 128`, use the finite-weight norming theorem to obtain the
order-eighteen half contraction, then iterate across the
`O(log Lambda)` dyadic levels toward Conjecture G1 and Simons Problem 4.2.
