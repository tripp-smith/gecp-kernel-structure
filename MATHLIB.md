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
- the reusable strict-positivity wrapper around `Matrix.det_vandermonde` for
  strictly increasing real tuples and the positive-determinant factorization
  of a square monomial-feature matrix with positive diagonal weights;
- the fewnomial bound `signVariations_lt_card_support` and its Descartes-rule
  corollary bounding positive roots by polynomial support cardinality;
- generalized Vandermonde nonsingularity for distinct positive nodes and
  exponents, together with strict determinant positivity for increasing nodes
  and exponents via a geometric-node homotopy, and nonnegativity at the
  nonnegative-node boundary by uniform positive shifts;
- the generic endpoint-interpolation stability pattern: nonnegative cardinal
  weights of total mass at most one turn an error in a fixed comparison space
  into a factor-two residual bound.

The fermionic identities, DLR counts, and research-outcome enumeration are
project-specific and are not current upstream candidates.

## Run provenance

Implementation-agent metadata is intentionally kept out of the upstream
candidate assessment. The model/version boundary, mode, elapsed time, token
snapshot, and API-equivalent cost assumptions are recorded in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata).
