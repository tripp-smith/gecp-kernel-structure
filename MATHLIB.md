# Mathlib candidates

No upstream pull request is claimed.

Potential candidates after further use and generalization:

- the dependent successful-elimination run and preservation of previously zero
  rows/columns;
- the determinant recurrence for a dependent successful elimination run via
  an explicit unit lower block-elimination matrix;
- the bound-based formulation equating PSD diagonal and complete-pivot
  suprema;
- same-domain positive-semidefinite preservation of a positive diagonal Schur
  update via an explicit congruence factorization;
- the metric grid-cover plus Lipschitz supremum lemma.
- the generic cross-ratio-control lemma for one-step Schur residual contraction
  and its geometric-sequence induction theorem.
- the sign-coherent cross-product criterion for nonexpansive exact Schur
  updates, and the variable-factor/cumulative-product contraction induction.
- the determinant identity for adjoining one row and column to a successful
  elimination run's selected core, and the resulting equivalence between a
  residual selected-cross sign condition and a four-bordered-minor product
  condition on the original matrix.
- the ordered-minor strict-sign-regularity abstraction, its extension to
  arbitrary injective tuples by sorting and determinant orientation, and its
  preservation under positive row or column scaling;
- the generic transfer from all-orders strict sign regularity to
  four-bordered-minor and selected-residual sign coherence for successful
  elimination runs;
- dependent composition of successful elimination runs, domain-restricted
  complete-pivot predicates, and propagation of a domain bound through a
  realized strictly sign-regular complete-pivot residual sequence;
- recursive selected-coordinate membership for a domain-complete dependent
  elimination run and the generic selected-core/sample-matrix bridge;
- the generic finite-product geometric-mean pattern for an antitone
  nonnegative sequence, together with the absolute determinant/pivot-product
  identity for a successful elimination run;
- the determinant expansion over column choices for a sum, the vanishing of a
  square rectangular product through a smaller finite type, and the resulting
  entrywise low-rank perturbation bound retaining the full
  `epsilon^(n-r)` power;
- the reusable strict-positivity wrapper around `Matrix.det_vandermonde` for
  strictly increasing real tuples and the positive-determinant factorization
  of a square monomial-feature matrix with positive diagonal weights;
- the fewnomial bound `signVariations_lt_card_support` and its Descartes-rule
  corollary bounding positive roots by polynomial support cardinality;
- generalized Vandermonde nonsingularity for distinct positive nodes and
  exponents, together with strict determinant positivity for increasing nodes
  and exponents via a geometric-node homotopy, and nonnegativity at the
  nonnegative-node boundary by uniform positive shifts;
- a rectangular Cauchy--Binet determinant expansion for finite ordered index
  types, indexed by strictly increasing intermediate tuples, together with its
  injective-tuple precursor;
- the reversal-permutation identity
  `Fin.revPerm.sign = (-1)^(n.choose 2)`, derived from the finite inversion
  product;
- the generic endpoint-interpolation stability pattern: nonnegative cardinal
  weights of total mass at most one turn an error in a fixed comparison space
  into a factor-two residual bound.

The fermionic identities, DLR counts, and research-outcome enumeration are
project-specific and are not current upstream candidates.

Phase Y reuses mathlib's existing `Matrix.IsHadamard.kronecker` and
`Matrix.IsHadamard.det_mul_star_det` API. Its recursive Sylvester wrapper and
constant-base obstruction are research-routing results for this project, not
new upstream candidates at their present level of generality.

Phase Z's multiplicative cross-ratio identity and corner-area estimates are
specific to the normalized fermionic exponential kernel. No upstream mathlib
abstraction is proposed until the same proof pattern is needed for a second
kernel family.

Phase AA uses mathlib's existing `Real.concaveOn_rpow` and antitonicity of
nonpositive real powers to sandwich a power secant between two elementary
functions. The two short secant wrappers are useful locally, but the reusable
content is already present upstream; no new mathlib candidate is proposed.

Phase AB relies on mathlib's existing exact `norm_num` support for rational
real powers whose roots are rational. The counterexample and its logarithmic
transport are project-specific, so this phase adds no upstream candidate.

Phase AC combines existing exponential monotonicity, `x*(1-x) <= 1/4`, the
hyperbolic double-angle identity, and the project-specific fermionic Schur
formulas. No missing general-purpose mathlib abstraction was encountered.

Phase AD packages an asymmetric local/global Schur-update estimate from
elementary absolute-value, division, and same-sign inequalities. Its generic
statement is useful inside this repository's `Kernel`, `CompletePivotOn`, and
`PivotCrossProductSignCoherent` APIs, but no missing lower-level mathlib lemma
was encountered and no upstream candidate is proposed.

Phase AE's `Matrix.abs_det_le_of_factors_approx_except` and closed coarse
corollary are credible upstream candidates after API review. They extend a
low-rank determinant perturbation estimate by representing a finite set of
columns exactly through `k ⊕ exceptional`, charging exactly the exceptional
cardinality in the effective rank. The proof is generic matrix algebra and
reuses `mixedLeft`, `mixedRight`, and `columnChoice`; only the fermionic
specializations remain project-specific.

Phase AF composes the generic exceptional-column estimate with this project's
dependent GECP run and fermionic strict-sign-regularity APIs. The composition
introduces no new general-purpose matrix lemma beyond the Phase AE candidates,
so no additional upstream proposal is recorded.

Phase AG's dependent two-corner run and recursive selected-index exception set
are specific to the fermionic pivot prescription. They reuse existing `Sum`,
`Finset`, and `Run.SelectedIndex` APIs without exposing a missing general
abstraction, so no additional mathlib candidate is proposed.

Phase AH adds a reusable constructor from an indexed residual recurrence to
the project's dependent `Run`, together with exact final-residual, pivot-list,
cardinality, and `CompleteOn` transport theorems. These abstractions are generic
within the local GECP API, but that API is not part of mathlib; no standalone
upstream mathlib candidate is proposed. The central-half specialization and
determinant composition remain fermionic-project results.

Phase AI reuses mathlib's `Nat.factorial_le_pow`,
`Nat.two_mul_sq_add_one_le_two_pow_two_mul`, ordered-power lemmas, and existing
real exponential bounds. The determinant sandwich is formulated around the
project's dependent GECP run and fermionic corner pivots; it exposes no missing
standalone mathlib abstraction. The conservative exponent-budget arithmetic is
also project-specific, so no new upstream candidate is recorded.

Phase AJ's recursive `Run.AugmentedIndex` and its determinant factorization
reuse mathlib's finite `Sum` indices, block-matrix operations, determinant
multiplicativity, and permutation signs. The construction is generic inside
the project's dependent GECP API but is not a standalone mathlib abstraction;
the fermionic specialization is project-specific. No new upstream candidate
is recorded.

Phase AK reuses `Real.sign`, diagonal determinant scaling, and the project's
cross-product sign-coherence API to formalize the balanced signing of a
bipartite kernel. The anchor-gauge construction is short and generic, but its
current formulation is specialized to the local `Kernel` and GECP residual
interfaces; no standalone mathlib candidate is proposed.

Phase AL packages exterior-power entries as a local compound kernel and proves
their balanced signing by elementary ordered-field algebra. Mathlib supplies
the determinant and finite-type foundations but currently exposes no directly
applicable Desnanot--Jacobi or Dodgson-condensation theorem; that identity is a
credible upstream candidate if the next quantitative phase requires a generic
matrix formulation.

Phase AM proves the needed two-border Desnanot--Jacobi identity by composing
the project's selected-core Schur factorizations with `Matrix.det_fin_two`.
That wrapper is tied to the local dependent `Run` API, but a generic finite
matrix Desnanot--Jacobi theorem remains a credible upstream candidate. The
elementary real lemma `abs_sub_le_max_of_mul_nonneg` may also be useful more
broadly, although its proof is short and no upstream proposal is warranted
without a second independent consumer.

Phase AN adds two generic recursive run predicates and proves that a uniformly
bounded, sign-coherent complete continuation has determinant at most `B^m`
through its exact pivot product. The algebra is reusable, but the statements
are expressed in the project's dependent `Run`, `CompleteOn`, and kernel APIs;
they are therefore local candidates rather than standalone mathlib proposals.
The proof also reuses dependent run transport, whose pivot-list preservation
is already listed above.

Phase AO uses mathlib's `Real.exp_bound`, Lagrange Taylor remainder support,
and explicit `Matrix.det_fin_two` formula to certify a repository-specific
mask-transition obstruction. The witness depends on the local dyadic masking
policy and does not suggest a standalone mathlib theorem. A reusable
alternating exponential-tail interval lemma was considered but was not needed:
the existing exponential estimates close the exact rational bounds directly.

Phase AP uses `map_add_eq_sum_add_integral_iteratedFDeriv` to obtain the exact
unit-interval beta representation of the eighth-order exponential tail.
Mathlib's existing Taylor-integral, iterated-derivative, interval-integral
positivity, and finite determinant tools are sufficient; no missing upstream
lemma was encountered. The resulting beta-moment minor reduction is specialized
to the local approximation target and is not proposed for upstreaming.

Phase AQ proves a local tilted-beta moment inequality with existing interval-
integral monotonicity and the fundamental theorem of calculus. The explicit
degree-ten primitive and the split at `u = 1/2` are specialized to the eighth-
order tail. The single-crossing comparison pattern may be reusable, but no
upstream abstraction is proposed until a second independent consumer appears.

Phase AR uses
`intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le` together
with the derivative mean-value criterion `strictAntiOn_of_deriv_neg`. The only
new abstraction is a private beta-power-moment differentiation helper; its
weight and local domination bound are specific to the eighth-order tail, so no
upstream candidate is recorded.

Phase AS combines the standard quotient derivative, derivative-based strict
monotonicity on a convex interval, and ordered-field cross multiplication. The
scale-ratio argument is a project-specific specialization of the classical
monotone-ratio criterion for multiplicative kernels; mathlib's existing API is
sufficient, so no upstream candidate is recorded.

Phase AT uses `intervalIntegral.smul_integral_comp_mul_left`, interval
additivity, and exact support vanishing to put a truncated-power composition on
a shared domain. The next step exposes a credible upstream gap: a finite-order
continuous Cauchy--Binet/Andreief identity for determinants of integral-kernel
compositions, including measurable/integrable hypotheses and a strictness
corollary. That abstraction would apply beyond this repository and is recorded
as a candidate rather than rederived prematurely inside the fermionic API.

Phase AU builds only short order-theoretic corollaries of the repository's
ordered rectangular Cauchy--Binet formula: termwise nonnegativity, retention
of a selected summand, and strictness from one positive term. These corollaries
would be natural companions if the finite rectangular formula were upstreamed,
but they do not constitute a separate mathlib proposal. The continuous
Cauchy--Binet/Andreief candidate recorded in Phase AT remains the relevant
upstream gap.

Phase AV finds no existing mathlib vocabulary for total nonnegativity of a
two-variable kernel. The local `TotallyNonnegativeAtOrder` predicate and the
order-two truncated-power proof are small, but a reusable kernel-level API
with finite and integral composition closure would be a credible upstream
abstraction once the all-orders consumer is complete. The exact support
criterion also indicates that such an API must distinguish TN from TP rather
than encoding only strict minors.

Phase AW uses existing interval-integral primitives (`integral_pow`, affine
substitution, adjacent-interval additivity, and equality on open intervals) to
prove the normalized truncated-power recurrence. No new special-function
lemma is needed. The remaining reusable upstream gap is still a kernel-level
continuous composition/Cauchy--Binet theorem that transfers nonnegative
ordered minors through integration; that abstraction should be proposed only
after the order-kernel consumer fixes the necessary measurability and
integrability hypotheses.

Phase AX's first-column induction is a specialized consecutive-ones/Ferrers
matrix argument. Mathlib's `Matrix.IsTotallyUnimodular` API is adjacent, but
using it would still require proving the consecutive-ones criterion and would
lose the stronger nonnegative `0`/`1` orientation obtained directly here. The
local proof is therefore appropriate; a general consecutive-ones total-
unimodularity theorem is a possible upstream result only if another consumer
appears.

Phase AY composes existing determinant row-linearity, automatic continuity of
finite-dimensional linear maps, and Bochner integral commutation. The local
`Matrix.det_updateRow_intervalIntegral` lemma is short and broadly reusable;
it is a credible mathlib companion to `det_updateRow_add` and
`det_updateRow_smul`. A larger Andreief theorem should wait until the
repository's ordered-interval iteration determines the cleanest hypotheses.

Phase AZ supplies that ordered-interval iteration using only `Finset` induction
and the Phase AY one-row theorem. The generic
`det_rowIntervalIntegral_nonneg` is a useful intermediate abstraction: it
requires pointwise determinant nonnegativity for every simultaneous interval
sample and concludes nonnegativity of the determinant of row integrals. Its
ordered-kernel specialization is repository-local vocabulary, while the
generic finite-row theorem and the one-row commutation lemma together form a
credible small upstream API independent of a full Andreief identity.

Phase BA uses `Matrix.det_eq_of_forall_row_eq_smul_add_pred` after reversing
the row order to prove determinant invariance under successive next-row
subtraction. The resulting `Matrix.det_successiveRowDifference` and
`Matrix.det_commonUpperRowIntegral_eq_adjacent` are plausible general-purpose
companions to the existing row-operation and interval-additivity APIs. The
kernel-level total-nonnegativity consequences remain project vocabulary, while
the matrix identities are credible upstream candidates with consumers beyond
truncated powers.

Phase BB proves the finite algebraic symmetrization underlying Andréief's
identity using existing determinant expansion and row-permutation lemmas. The
generic `Matrix.sum_compositionSampleIntegrand_comp_perm` is independent of
analysis and is a plausible supporting lemma for an eventual mathlib
continuous Cauchy--Binet theorem. The next upstream-relevant question is which
finite-product-measure permutation API gives the cleanest integration proof.

Phase BC answers that question with `measurePreserving_piCongrLeft`, its
integrability and integral transport lemmas, finite-sum integration, and
`integral_fintype_prod_eq_prod`. The resulting generic division-free
finite-order Andreief identity is independent of the fermionic application
and is now a concrete upstream candidate, together with its coordinate-
permutation invariance helpers. Before proposing it upstream, the hypotheses
should be compared with a Bochner-valued or dependent-measure formulation;
the present real-valued common-measure statement is the minimal theorem needed
by this repository.

Phase BD uses `Tuple.sort` to isolate a generic sign-transfer lemma: a product
of two sampled determinants is invariant when the same permutation is applied
to the left columns and right rows, so an ordered-sample weak sign extends to
all samples after separating injective and repeated cases. Together with the
Phase BC Andreief theorem this yields a reusable weak-sign composition-closure
result. Both lemmas are plausible upstream companions to a continuous
Cauchy--Binet API; the Taylor-tail specialization remains project-local.

Phase BE adds strictness with a classical Schoenberg--Whitney interlacing
witness. The generic analytic pattern is broader than this kernel: a
continuous nonnegative Andreief integrand with one strict interior point is
localized by a compactly supported bump, giving positive integral and hence a
strict composition determinant. A reusable upstream theorem should abstract
that positive-support argument from the project-specific midpoint witness;
the current implementation keeps the abstraction local until its minimal
topological and measure hypotheses are clearer.

Phase BF is primarily an application-level composition of existing generic
results: subtype restriction and positive column scaling turn a strict-minor
family into a `StrictSignRegular` kernel, after which complete-pivot
nonexpansiveness yields a factorial-free determinant envelope. The reusable
upstream candidate remains the generic continuation theorem from Phase AN;
the new tail and rank-eight decomposition are fermionic-specific.

Phase BG uses a project-specific reflection to replace an alternating Taylor
feature family by positive monomial features. The reusable ingredient is a
potential coefficientwise Cauchy--Binet theorem for a feature matrix split
into two ordered index blocks with nonnegative weights. That abstraction is
not yet formalized; the current strict-tail and mixed-coefficient results
should remain local until the precise finite statement is proved.

Phase BH formalizes that reusable finite statement locally. The generic
upstream candidate is a block-weighted rectangular Cauchy--Binet theorem:
mark columns of one factor, scale them by a scalar, and group ordered minor
terms by the number of marked columns selected. The core theorem is independent
of the fermionic kernel; only the nonnegative generalized-Vandermonde
specialization is project-facing. Before upstreaming, it would be worth
generalizing the coefficient semiring and checking whether mathlib's
polynomial matrix API offers a more canonical coefficient formulation.

Phase BI finds that canonical formulation: the weighted feature expansion and
the physical mixed-column expansion are two presentations of the same
determinant polynomial. The generic upstream candidate now includes
`mixedColumnPolynomial` and its coefficient theorem, plus the equality with a
block-feature Cauchy--Binet polynomial. The proof uses existing polynomial
extensionality over `Real`; an upstream version should be stated over the
weakest infinite integral-domain or sufficiently-large-domain hypotheses.

Phase BJ adds a generic topological companion: fixed-cardinality sums of
mixed-column determinants are jointly continuous in both matrices. This is a
small reusable lemma over finite matrix types; an upstream version should be
generalized from real matrices to a topological commutative ring. The
finite-to-infinite Taylor-tail specialization remains project-specific.

Phase BK adds a generic coefficient-support consequence: if the background
matrix factors through `k`, multilinearity kills coefficients below degree
`n-k`; nonnegative coefficients then give a monotone determinant-pencil
amplification inequality. This is a plausible polynomial/matrix upstream
lemma once generalized beyond real matrices and phrased using an ordered
coefficient semiring. The eighth-order and fermionic instances remain local.

Phase BL adds the simpler companion
`Matrix.mixedColumnPolynomial_det_mono`: a determinant pencil is monotone on
nonnegative weights when all of its mixed-column coefficients are
nonnegative. This is a plausible upstream lemma beside the coefficient and
evaluation API, generalized from real matrices to a suitable ordered
commutative ring. The amplified fermionic kernel, its strict positivity, and
the numerical choice `8!` remain project-specific.

Phase BM adds a reusable finite-matrix abstraction around strict total
positivity. The main upstream candidates are preservation under a scalar
northwest Schur complement and
`Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive`. The current proof is
specialized to real square matrices and reuses this project's augmented GECP
determinant identity; an upstream version should be generalized to an ordered
field and connected to mathlib's existing Schur-complement API if that avoids
extra index equivalences. The positive-series kernel specialization and the
constants `8!` and `2` remain project-specific.

Phase BN adds the generic order-invariance companion
`abs_minorMatrix_det_le_pow_of_strictSignRegular_one`: independently sort two
injective sample tuples, apply the ordered positive-minor envelope, and remove
both permutation signs through absolute value. A more general upstream theorem
could accept an arbitrary ordered-sample determinant bound instead of the
current diagonal-product specialization. The positive-series GECP composition
and rank-ten arithmetic remain project-specific.

Phase BO needs no new general-purpose matrix primitive. Its key simplification
is kernel-specific: the affine coordinate `z = 1 - omega` converts the whole
signed unit band into a positive exponential product band, so the existing
mixed-column polynomial, total-positivity, and determinant-permutation APIs
compose directly. The reusable upstream candidates therefore remain the
generic amplification and strictly-totally-positive Hadamard lemmas identified
in Phases BK--BN; the constants `128`, `2`, and the affine fermionic identity
are project-specific.

Phase BP likewise needs no new general-purpose matrix primitive. It observes
that the Phase BO diagonal-factor and tail-amplification proof is uniform under
an arbitrary translation of the frequency window: the positive product
coordinate still has width two, while the weighted exponential is the
original fermionic kernel and remains bounded by one. The center-parametric
kernel identities and constants are project-specific; the reusable upstream
candidates remain the mixed-column rank amplification and strictly-total-
positive diagonal-product lemmas already identified in Phases BK--BN.

## Run provenance

Implementation-agent metadata is intentionally kept out of the upstream
candidate assessment. The model/version boundary, mode, elapsed time, token
snapshot, and API-equivalent cost assumptions are recorded in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata).
