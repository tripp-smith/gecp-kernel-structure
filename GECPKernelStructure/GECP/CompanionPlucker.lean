import GECPKernelStructure.GECP.CompanionBiConditioning
import GECPKernelStructure.GECP.AugmentedDeterminant

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w z t

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

/-- Selected kernel columns followed by arbitrary columns of a left feature
factor. -/
def leftFeatureColumn {γ : Type z} {δ : Type t} {K : Kernel α β 𝕜} :
    (run : Run K) → (γ → δ) → run.AugmentedIndex γ → β ⊕ δ
  | .nil _, features, i => Sum.inr (features i)
  | .step _ column _ tail, features, i =>
      Sum.elim (fun _ => Sum.inl column)
        (tail.leftFeatureColumn features) i

/-- The original selected kernel core bordered on the right by columns of an
arbitrary left feature factor. -/
def leftFeatureCore {γ : Type z} {δ : Type t} {K : Kernel α β 𝕜}
    (run : Run K) (rows : γ → α) (features : γ → δ)
    (L : Matrix α δ 𝕜) :
    Matrix (run.AugmentedIndex γ) (run.AugmentedIndex γ) 𝕜 :=
  fun i j =>
    Sum.elim
      (fun column => K (run.augmentedRow γ rows i) column)
      (fun feature => L (run.augmentedRow γ rows i) feature)
      (run.leftFeatureColumn features j)

/-- Selected kernel rows followed by arbitrary rows of a right feature
factor. -/
def rightFeatureRow {γ : Type z} {δ : Type t} {K : Kernel α β 𝕜} :
    (run : Run K) → (γ → δ) → run.AugmentedIndex γ → α ⊕ δ
  | .nil _, features, i => Sum.inr (features i)
  | .step row _ _ tail, features, i =>
      Sum.elim (fun _ => Sum.inl row)
        (tail.rightFeatureRow features) i

/-- The original selected kernel core bordered below by rows of an arbitrary
right feature factor. -/
def rightFeatureCore {γ : Type z} {δ : Type t} {K : Kernel α β 𝕜}
    (run : Run K) (features : γ → δ) (columns : γ → β)
    (R : Matrix δ β 𝕜) :
    Matrix (run.AugmentedIndex γ) (run.AugmentedIndex γ) 𝕜 :=
  fun i j =>
    Sum.elim
      (fun row => K row (run.augmentedColumn γ columns j))
      (fun feature => R feature (run.augmentedColumn γ columns j))
      (run.rightFeatureRow features i)

@[simp]
theorem leftFeatureCore_nil {γ : Type z} {δ : Type t}
    (K : Kernel α β 𝕜) (rows : γ → α) (features : γ → δ)
    (L : Matrix α δ 𝕜) :
    (Run.nil K).leftFeatureCore rows features L =
      fun i j => L (rows i) (features j) := by
  rfl

@[simp]
theorem rightFeatureCore_nil {γ : Type z} {δ : Type t}
    (K : Kernel α β 𝕜) (features : γ → δ) (columns : γ → β)
    (R : Matrix δ β 𝕜) :
    (Run.nil K).rightFeatureCore features columns R =
      fun i j => R (features i) (columns j) := by
  rfl

/-- One selected pivot factors a left-feature bordered determinant and applies
the corresponding fixed row operation to its feature border. -/
theorem leftFeatureCore_step_det
    {γ : Type z} [Fintype γ] [DecidableEq γ] {δ : Type t}
    {K : Kernel α β 𝕜} (row : α) (column : β)
    (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (rows : γ → α) (features : γ → δ) (L : Matrix α δ 𝕜) :
    ((Run.step row column pivot_ne tail).leftFeatureCore
        rows features L).det =
      K row column *
        (tail.leftFeatureCore rows features
          (companionRowUpdate K L row column pivot_ne)).det := by
  classical
  let ε := tail.AugmentedIndex γ
  let core : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    (Run.step row column pivot_ne tail).leftFeatureCore rows features L
  let pivotBlock : Matrix (Fin 1) (Fin 1) 𝕜 := fun _ _ => K row column
  let topRows : Matrix (Fin 1) ε 𝕜 := fun _ j =>
    Sum.elim (fun y => K row y) (fun d => L row d)
      (tail.leftFeatureColumn features j)
  let leftColumns : Matrix ε (Fin 1) 𝕜 := fun i _ =>
    K (tail.augmentedRow γ rows i) column
  let originalTail : Matrix ε ε 𝕜 := fun i j =>
    Sum.elim
      (fun y => K (tail.augmentedRow γ rows i) y)
      (fun d => L (tail.augmentedRow γ rows i) d)
      (tail.leftFeatureColumn features j)
  let multiplier : Matrix ε (Fin 1) 𝕜 := fun i _ =>
    -(K (tail.augmentedRow γ rows i) column / K row column)
  let lower : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    Matrix.fromBlocks (1 : Matrix (Fin 1) (Fin 1) 𝕜)
      (0 : Matrix (Fin 1) ε 𝕜) multiplier (1 : Matrix ε ε 𝕜)
  let eliminated : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix ε (Fin 1) 𝕜)
      (tail.leftFeatureCore rows features
        (companionRowUpdate K L row column pivot_ne))
  have hcoreBlocks :
      core = Matrix.fromBlocks pivotBlock topRows leftColumns originalTail := by
    ext (i | i) (j | j) <;> rfl
  have hleft : multiplier * pivotBlock + leftColumns = 0 := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, pivotBlock, leftColumns, pivot_ne]
  have htail : multiplier * topRows + originalTail =
      tail.leftFeatureCore rows features
        (companionRowUpdate K L row column pivot_ne) := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    cases hcolumn : tail.leftFeatureColumn features j with
    | inl y =>
        simp [multiplier, topRows, originalTail, leftFeatureCore,
          hcolumn, residualUpdate]
        ring
    | inr d =>
        simp [multiplier, topRows, originalTail, leftFeatureCore,
          hcolumn, companionRowUpdate]
        ring
  have hmul : lower * core = eliminated := by
    rw [hcoreBlocks]
    simp only [lower, Matrix.fromBlocks_multiply, Matrix.one_mul,
      Matrix.zero_mul, add_zero, hleft, htail, eliminated]
  have hlower : lower.det = 1 := by
    simp [lower]
  have heliminated : eliminated.det =
      K row column *
        (tail.leftFeatureCore rows features
          (companionRowUpdate K L row column pivot_ne)).det := by
    rw [show eliminated = Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix ε (Fin 1) 𝕜)
      (tail.leftFeatureCore rows features
        (companionRowUpdate K L row column pivot_ne)) by rfl]
    rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fin_one]
  calc
    ((Run.step row column pivot_ne tail).leftFeatureCore
        rows features L).det = core.det := rfl
    _ = lower.det * core.det := by rw [hlower, one_mul]
    _ = (lower * core).det := by rw [Matrix.det_mul]
    _ = eliminated.det := congrArg Matrix.det hmul
    _ = _ := heliminated

/-- One selected pivot factors a right-feature bordered determinant and
applies the corresponding fixed column operation to its feature border. -/
theorem rightFeatureCore_step_det
    {γ : Type z} [Fintype γ] [DecidableEq γ] {δ : Type t}
    {K : Kernel α β 𝕜} (row : α) (column : β)
    (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (features : γ → δ) (columns : γ → β) (R : Matrix δ β 𝕜) :
    ((Run.step row column pivot_ne tail).rightFeatureCore
        features columns R).det =
      K row column *
        (tail.rightFeatureCore features columns
          (companionColumnUpdate K R row column pivot_ne)).det := by
  classical
  let ε := tail.AugmentedIndex γ
  let core : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    (Run.step row column pivot_ne tail).rightFeatureCore features columns R
  let pivotBlock : Matrix (Fin 1) (Fin 1) 𝕜 := fun _ _ => K row column
  let topRows : Matrix (Fin 1) ε 𝕜 := fun _ j =>
    K row (tail.augmentedColumn γ columns j)
  let leftColumns : Matrix ε (Fin 1) 𝕜 := fun i _ =>
    Sum.elim (fun x => K x column) (fun d => R d column)
      (tail.rightFeatureRow features i)
  let originalTail : Matrix ε ε 𝕜 := fun i j =>
    Sum.elim
      (fun x => K x (tail.augmentedColumn γ columns j))
      (fun d => R d (tail.augmentedColumn γ columns j))
      (tail.rightFeatureRow features i)
  let multiplier : Matrix ε (Fin 1) 𝕜 := fun i _ =>
    -(Sum.elim (fun x => K x column) (fun d => R d column)
        (tail.rightFeatureRow features i) / K row column)
  let lower : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    Matrix.fromBlocks (1 : Matrix (Fin 1) (Fin 1) 𝕜)
      (0 : Matrix (Fin 1) ε 𝕜) multiplier (1 : Matrix ε ε 𝕜)
  let eliminated : Matrix (Fin 1 ⊕ ε) (Fin 1 ⊕ ε) 𝕜 :=
    Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix ε (Fin 1) 𝕜)
      (tail.rightFeatureCore features columns
        (companionColumnUpdate K R row column pivot_ne))
  have hcoreBlocks :
      core = Matrix.fromBlocks pivotBlock topRows leftColumns originalTail := by
    ext (i | i) (j | j) <;> rfl
  have hleft : multiplier * pivotBlock + leftColumns = 0 := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, pivotBlock, leftColumns, pivot_ne]
  have htail : multiplier * topRows + originalTail =
      tail.rightFeatureCore features columns
        (companionColumnUpdate K R row column pivot_ne) := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    cases hrow : tail.rightFeatureRow features i with
    | inl x =>
        simp [multiplier, topRows, originalTail, rightFeatureCore,
          hrow, residualUpdate]
        ring
    | inr d =>
        simp [multiplier, topRows, originalTail, rightFeatureCore,
          hrow, companionColumnUpdate]
        ring
  have hmul : lower * core = eliminated := by
    rw [hcoreBlocks]
    simp only [lower, Matrix.fromBlocks_multiply, Matrix.one_mul,
      Matrix.zero_mul, add_zero, hleft, htail, eliminated]
  have hlower : lower.det = 1 := by
    simp [lower]
  have heliminated : eliminated.det =
      K row column *
        (tail.rightFeatureCore features columns
          (companionColumnUpdate K R row column pivot_ne)).det := by
    rw [show eliminated = Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix ε (Fin 1) 𝕜)
      (tail.rightFeatureCore features columns
        (companionColumnUpdate K R row column pivot_ne)) by rfl]
    rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fin_one]
  calc
    ((Run.step row column pivot_ne tail).rightFeatureCore
        features columns R).det = core.det := rfl
    _ = lower.det * core.det := by rw [hlower, one_mul]
    _ = (lower * core).det := by rw [Matrix.det_mul]
    _ = eliminated.det := congrArg Matrix.det hmul
    _ = _ := heliminated

/-- A left conditioned feature minor is the exact Schur quotient of its
original feature-bordered determinant. -/
theorem leftFeatureCore_det_eq_selectedCore_det_mul_conditionRows_minor
    {γ : Type z} [Fintype γ] [DecidableEq γ] {δ : Type t}
    {K : Kernel α β 𝕜} (run : Run K)
    (rows : γ → α) (features : γ → δ) (L : Matrix α δ 𝕜) :
    (run.leftFeatureCore rows features L).det =
      run.selectedCore.det *
        Matrix.det (fun i j => run.conditionRows L (rows i) (features j)) := by
  induction run generalizing L with
  | nil initial =>
      have selected_det : (Run.nil initial).selectedCore.det = 1 := by
        rw [selectedCore_nil]
        exact Matrix.det_fin_zero
      change Matrix.det (fun i j => L (rows i) (features j)) =
        (Run.nil initial).selectedCore.det *
          Matrix.det (fun i j => L (rows i) (features j))
      rw [selected_det, one_mul]
  | step row column pivot_ne tail ih =>
      rw [leftFeatureCore_step_det, selectedCore_step_det, ih]
      simp only [conditionRows_step]
      ring

/-- A right conditioned feature minor is the exact Schur quotient of its
original feature-bordered determinant. -/
theorem rightFeatureCore_det_eq_selectedCore_det_mul_conditionColumns_minor
    {γ : Type z} [Fintype γ] [DecidableEq γ] {δ : Type t}
    {K : Kernel α β 𝕜} (run : Run K)
    (features : γ → δ) (columns : γ → β) (R : Matrix δ β 𝕜) :
    (run.rightFeatureCore features columns R).det =
      run.selectedCore.det *
        Matrix.det (fun i j => run.conditionColumns R (features i) (columns j)) := by
  induction run generalizing R with
  | nil initial =>
      have selected_det : (Run.nil initial).selectedCore.det = 1 := by
        rw [selectedCore_nil]
        exact Matrix.det_fin_zero
      change Matrix.det (fun i j => R (features i) (columns j)) =
        (Run.nil initial).selectedCore.det *
          Matrix.det (fun i j => R (features i) (columns j))
      rw [selected_det, one_mul]
  | step row column pivot_ne tail ih =>
      rw [rightFeatureCore_step_det, selectedCore_step_det, ih]
      simp only [conditionColumns_step]
      ring

/-- The product of matched conditioned feature minors is represented entirely
by original bordered determinants, with the selected core entering as a
square factor. -/
theorem leftFeatureCore_det_mul_rightFeatureCore_det
    {γ : Type z} [Fintype γ] [DecidableEq γ] {δ : Type t}
    {K : Kernel α β 𝕜} (run : Run K)
    (rows : γ → α) (columns : γ → β) (features : γ → δ)
    (L : Matrix α δ 𝕜) (R : Matrix δ β 𝕜) :
    (run.leftFeatureCore rows features L).det *
        (run.rightFeatureCore features columns R).det =
      run.selectedCore.det ^ 2 *
        (Matrix.det (fun i j => run.conditionRows L (rows i) (features j)) *
          Matrix.det (fun i j =>
            run.conditionColumns R (features i) (columns j))) := by
  rw [run.leftFeatureCore_det_eq_selectedCore_det_mul_conditionRows_minor,
    run.rightFeatureCore_det_eq_selectedCore_det_mul_conditionColumns_minor]
  ring

end Run
end GECP
end GECPKernelStructure
