import GECPKernelStructure.PositiveDefinite.ResidualPSD

namespace GECPKernelStructure
namespace PositiveDefinite

open Matrix

universe u

/-- A pivoted-Cholesky pivot maximizes the residual diagonal. -/
def IsDiagonalPivot {n : Type u} (M : Matrix n n ℝ) (pivot : n) : Prop :=
  ∀ i, M i i ≤ M pivot pivot

/-- A diagonal-preferring GECP pivot dominates every residual entry in absolute value. -/
def IsCompletePivot {n : Type u} (M : Matrix n n ℝ) (pivot : n) : Prop :=
  ∀ i j, |M i j| ≤ M pivot pivot

/-- The least diagonal maximizer under the deterministic node order. -/
def IsCanonicalDiagonalPivot {n : Type u} [LinearOrder n]
    (M : Matrix n n ℝ) (pivot : n) : Prop :=
  IsDiagonalPivot M pivot ∧ ∀ candidate, IsDiagonalPivot M candidate → pivot ≤ candidate

/-- The diagonal-preferring complete pivot with the same deterministic tie rule. -/
def IsCanonicalCompletePivot {n : Type u} [LinearOrder n]
    (M : Matrix n n ℝ) (pivot : n) : Prop :=
  IsCompletePivot M pivot ∧ ∀ candidate, IsDiagonalPivot M candidate → pivot ≤ candidate

/-- On a finite PSD residual, diagonal and diagonal-preferring complete pivots coincide. -/
theorem diagonalPivot_iff_completePivot {n : Type u} [Fintype n]
    (M : Matrix n n ℝ) (hM : M.PosSemidef) (pivot : n) :
    IsDiagonalPivot M pivot ↔ IsCompletePivot M pivot :=
  posDef_gecp_residual_sup_eq_diag_sup M hM (M pivot pivot)

/-- Canonical lexicographic tie-breaking is unchanged by the PSD reduction. -/
theorem canonicalDiagonalPivot_iff_canonicalCompletePivot {n : Type u}
    [Fintype n] [LinearOrder n] (M : Matrix n n ℝ) (hM : M.PosSemidef)
    (pivot : n) :
    IsCanonicalDiagonalPivot M pivot ↔ IsCanonicalCompletePivot M pivot := by
  rw [IsCanonicalDiagonalPivot, IsCanonicalCompletePivot,
    diagonalPivot_iff_completePivot M hM pivot]

/-- A matrix state carrying the invariant needed by recursive PSD pivoting. -/
structure PSDState (n : Type u) [Fintype n] [DecidableEq n] where
  matrix : Matrix n n ℝ
  posSemidefinite : matrix.PosSemidef

namespace PSDState

/-- The next PSD state after a positive diagonal Schur update. -/
noncomputable def next {n : Type u} [Fintype n] [DecidableEq n]
    (state : PSDState n) (pivot : n) (pivot_pos : 0 < state.matrix pivot pivot) :
    PSDState n where
  matrix := diagonalResidual state.matrix pivot
  posSemidefinite :=
    posDef_gecp_residual_posSemidefinite state.matrix state.posSemidefinite pivot pivot_pos

end PSDState

/-- Recursive canonical diagonal-pivoted Cholesky trace. -/
inductive PivotedCholeskyTrace {n : Type u} [Fintype n] [DecidableEq n]
    [LinearOrder n] : PSDState n → List n → Prop where
  | nil (state : PSDState n) : PivotedCholeskyTrace state []
  | step {state : PSDState n} (pivot : n)
      (pivot_pos : 0 < state.matrix pivot pivot)
      (pivot_rule : IsCanonicalDiagonalPivot state.matrix pivot)
      {pivots : List n}
      (tail : PivotedCholeskyTrace (state.next pivot pivot_pos) pivots) :
      PivotedCholeskyTrace state (pivot :: pivots)

/-- Recursive canonical diagonal-preferring complete-pivot GECP trace. -/
inductive GECPTrace {n : Type u} [Fintype n] [DecidableEq n]
    [LinearOrder n] : PSDState n → List n → Prop where
  | nil (state : PSDState n) : GECPTrace state []
  | step {state : PSDState n} (pivot : n)
      (pivot_pos : 0 < state.matrix pivot pivot)
      (pivot_rule : IsCanonicalCompletePivot state.matrix pivot)
      {pivots : List n}
      (tail : GECPTrace (state.next pivot pivot_pos) pivots) :
      GECPTrace state (pivot :: pivots)

theorem GECPTrace.toPivotedCholesky {n : Type u} [Fintype n]
    [DecidableEq n] [LinearOrder n] {state : PSDState n} {pivots : List n}
    (trace : GECPTrace state pivots) : PivotedCholeskyTrace state pivots := by
  induction trace with
  | nil state => exact PivotedCholeskyTrace.nil state
  | @step state pivot pivot_pos pivot_rule pivots tail ih =>
      exact PivotedCholeskyTrace.step pivot pivot_pos
        ((canonicalDiagonalPivot_iff_canonicalCompletePivot state.matrix
          state.posSemidefinite pivot).mpr pivot_rule) ih

theorem PivotedCholeskyTrace.toGECP {n : Type u} [Fintype n]
    [DecidableEq n] [LinearOrder n] {state : PSDState n} {pivots : List n}
    (trace : PivotedCholeskyTrace state pivots) : GECPTrace state pivots := by
  induction trace with
  | nil state => exact GECPTrace.nil state
  | @step state pivot pivot_pos pivot_rule pivots tail ih =>
      exact GECPTrace.step pivot pivot_pos
        ((canonicalDiagonalPivot_iff_canonicalCompletePivot state.matrix
          state.posSemidefinite pivot).mp pivot_rule) ih

/-- Canonical complete-pivot GECP and diagonal-pivoted Cholesky have exactly the same
recursive pivot traces on every finite PSD matrix. -/
theorem gecp_eq_pivotedCholesky_of_posDef {n : Type u} [Fintype n]
    [DecidableEq n] [LinearOrder n] (M : Matrix n n ℝ) (hM : M.PosSemidef)
    (pivots : List n) :
    GECPTrace ⟨M, hM⟩ pivots ↔ PivotedCholeskyTrace ⟨M, hM⟩ pivots := by
  exact ⟨GECPTrace.toPivotedCholesky, PivotedCholeskyTrace.toGECP⟩

end PositiveDefinite
end GECPKernelStructure
