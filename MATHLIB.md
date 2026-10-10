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

## Run provenance

Implementation-agent metadata is intentionally kept out of the upstream
candidate assessment. The model/version boundary, mode, elapsed time, token
snapshot, and API-equivalent cost assumptions are recorded in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata).
