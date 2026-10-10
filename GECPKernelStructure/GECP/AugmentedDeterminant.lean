import GECPKernelStructure.GECP.SignRegularity

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w z

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

/-- The recursively indexed selected core is nonsingular without reindexing. -/
theorem selectedCore_det_ne {K : Kernel α β 𝕜} (run : Run K) :
    run.selectedCore.det ≠ 0 := by
  intro selected_zero
  apply gecp_core_nonsingular run
  rw [finSelectedCore, Matrix.det_reindex_self, selected_zero]

/-- Successful elimination never selects the same row coordinate twice. -/
theorem selectedRow_injective {K : Kernel α β 𝕜} (run : Run K) :
    Function.Injective run.selectedRow := by
  intro i j rows_equal
  by_contra indices_ne
  apply run.selectedCore_det_ne
  apply Matrix.det_zero_of_row_eq indices_ne
  funext k
  exact congrArg (fun row => K row (run.selectedColumn k)) rows_equal

/-- Successful elimination never selects the same column coordinate twice. -/
theorem selectedColumn_injective {K : Kernel α β 𝕜} (run : Run K) :
    Function.Injective run.selectedColumn := by
  intro i j columns_equal
  by_contra indices_ne
  apply run.selectedCore_det_ne
  apply Matrix.det_zero_of_column_eq indices_ne
  intro k
  exact congrArg (fun column => K (run.selectedRow k) column) columns_equal

/--
The selected positions of a successful run followed by an arbitrary finite
border, retaining the recursive pivot order without repeated reindexing.
-/
def AugmentedIndex (γ : Type z) {K : Kernel α β 𝕜} : Run K → Type z
  | .nil _ => γ
  | .step _ _ _ tail => Fin 1 ⊕ tail.AugmentedIndex γ

noncomputable instance augmentedIndexFintype (γ : Type z) [Fintype γ]
    {K : Kernel α β 𝕜} (run : Run K) : Fintype (run.AugmentedIndex γ) := by
  induction run with
  | nil =>
      dsimp [AugmentedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [AugmentedIndex]
      letI := ih
      infer_instance

noncomputable instance augmentedIndexDecidableEq (γ : Type z) [Fintype γ]
    [DecidableEq γ] {K : Kernel α β 𝕜} (run : Run K) :
    DecidableEq (run.AugmentedIndex γ) := by
  induction run with
  | nil =>
      dsimp [AugmentedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [AugmentedIndex]
      letI := ih
      infer_instance

/-- The selected rows followed recursively by all rows in an arbitrary border. -/
def augmentedRow (γ : Type z) {K : Kernel α β 𝕜} :
    (run : Run K) → (γ → α) → run.AugmentedIndex γ → α
  | .nil _, rows, i => rows i
  | .step selected _ _ tail, rows, i =>
      Sum.elim (fun _ => selected) (tail.augmentedRow γ rows) i

/-- The selected columns followed recursively by all columns in an arbitrary border. -/
def augmentedColumn (γ : Type z) {K : Kernel α β 𝕜} :
    (run : Run K) → (γ → β) → run.AugmentedIndex γ → β
  | .nil _, columns, i => columns i
  | .step _ selected _ tail, columns, i =>
      Sum.elim (fun _ => selected) (tail.augmentedColumn γ columns) i

/-- The original kernel sampled on the selected core plus an arbitrary border. -/
def augmentedCore (γ : Type z) {K : Kernel α β 𝕜} (run : Run K)
    (rows : γ → α) (columns : γ → β) :
    Matrix (run.AugmentedIndex γ) (run.AugmentedIndex γ) 𝕜 :=
  fun i j => K (run.augmentedRow γ rows i) (run.augmentedColumn γ columns j)

/-- Orientation of the combined selected and border row tuple. -/
noncomputable def augmentedRowOrientation (γ : Type z) [Fintype γ]
    [LinearOrder α] {K : Kernel α β 𝕜} (run : Run K) (rows : γ → α) : ℝ :=
  tupleOrientation
    (run.augmentedRow γ rows ∘ (Fintype.equivFin (run.AugmentedIndex γ)).symm)

/-- Orientation of the combined selected and border column tuple. -/
noncomputable def augmentedColumnOrientation (γ : Type z) [Fintype γ]
    [LinearOrder β] {K : Kernel α β 𝕜} (run : Run K) (columns : γ → β) : ℝ :=
  tupleOrientation
    (run.augmentedColumn γ columns ∘
      (Fintype.equivFin (run.AugmentedIndex γ)).symm)

@[simp]
theorem augmentedCore_nil (γ : Type z) (K : Kernel α β 𝕜)
    (rows : γ → α) (columns : γ → β) :
    (Run.nil K).augmentedCore γ rows columns =
      (fun i j => K (rows i) (columns j)) := by
  rfl

/-- One elimination step factors every arbitrary-border determinant. -/
theorem augmentedCore_step_det (γ : Type z) [Fintype γ] [DecidableEq γ]
    {K : Kernel α β 𝕜} (row : α) (column : β)
    (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (rows : γ → α) (columns : γ → β) :
    ((Run.step row column pivot_ne tail).augmentedCore γ rows columns).det =
      K row column * (tail.augmentedCore γ rows columns).det := by
  classical
  let δ := tail.AugmentedIndex γ
  let core : Matrix (Fin 1 ⊕ δ) (Fin 1 ⊕ δ) 𝕜 :=
    fun i j => K (Sum.elim (fun _ => row) (tail.augmentedRow γ rows) i)
      (Sum.elim (fun _ => column) (tail.augmentedColumn γ columns) j)
  have hcore :
      (Run.step row column pivot_ne tail).augmentedCore γ rows columns = core := rfl
  let pivotBlock : Matrix (Fin 1) (Fin 1) 𝕜 := fun _ _ => K row column
  let topRows : Matrix (Fin 1) δ 𝕜 :=
    fun _ j => K row (tail.augmentedColumn γ columns j)
  let leftColumns : Matrix δ (Fin 1) 𝕜 :=
    fun i _ => K (tail.augmentedRow γ rows i) column
  let originalTail : Matrix δ δ 𝕜 :=
    fun i j => K (tail.augmentedRow γ rows i)
      (tail.augmentedColumn γ columns j)
  let multiplier : Matrix δ (Fin 1) 𝕜 :=
    fun i _ => -(K (tail.augmentedRow γ rows i) column / K row column)
  let lower : Matrix (Fin 1 ⊕ δ) (Fin 1 ⊕ δ) 𝕜 :=
    Matrix.fromBlocks (1 : Matrix (Fin 1) (Fin 1) 𝕜)
      (0 : Matrix (Fin 1) δ 𝕜) multiplier (1 : Matrix δ δ 𝕜)
  let eliminated : Matrix (Fin 1 ⊕ δ) (Fin 1 ⊕ δ) 𝕜 :=
    Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix δ (Fin 1) 𝕜) (tail.augmentedCore γ rows columns)
  have hcoreBlocks :
      core = Matrix.fromBlocks pivotBlock topRows leftColumns originalTail := by
    ext (i | i) (j | j) <;> rfl
  have hleft : multiplier * pivotBlock + leftColumns = 0 := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, pivotBlock, leftColumns, pivot_ne]
  have htail :
      multiplier * topRows + originalTail = tail.augmentedCore γ rows columns := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, topRows, originalTail, augmentedCore, residualUpdate]
    ring
  have hmul : lower * core = eliminated := by
    rw [hcoreBlocks]
    simp only [lower, Matrix.fromBlocks_multiply, Matrix.one_mul, Matrix.zero_mul,
      add_zero, hleft, htail, eliminated]
  have hlower : lower.det = 1 := by
    simp [lower]
  have heliminated : eliminated.det =
      K row column * (tail.augmentedCore γ rows columns).det := by
    rw [show eliminated = Matrix.fromBlocks pivotBlock topRows
      (0 : Matrix δ (Fin 1) 𝕜) (tail.augmentedCore γ rows columns) by rfl]
    rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fin_one]
  calc
    ((Run.step row column pivot_ne tail).augmentedCore γ rows columns).det =
        core.det := congrArg Matrix.det hcore
    _ = lower.det * core.det := by rw [hlower, one_mul]
    _ = (lower * core).det := by rw [Matrix.det_mul]
    _ = eliminated.det := congrArg Matrix.det hmul
    _ = K row column * (tail.augmentedCore γ rows columns).det := heliminated

/--
Every multi-bordered determinant equals the selected-core determinant times
the corresponding finite minor of the exact final residual.
-/
theorem augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    {K : Kernel α β 𝕜} (run : Run K)
    (rows : γ → α) (columns : γ → β) :
    (run.augmentedCore γ rows columns).det =
      run.selectedCore.det *
        Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) := by
  induction run with
  | nil initial =>
      have selected_det : (Run.nil initial).selectedCore.det = 1 := by
        rw [selectedCore_nil]
        exact Matrix.det_fin_zero
      change Matrix.det (fun i j : γ => initial (rows i) (columns j)) =
        (Run.nil initial).selectedCore.det *
          Matrix.det (fun i j : γ => initial (rows i) (columns j))
      rw [selected_det, one_mul]
  | step selectedRow selectedColumn pivot_ne tail ih =>
      rw [augmentedCore_step_det, selectedCore_step_det, ih]
      simp only [finalResidual]
      ring

/-- A nonzero augmented original-kernel minor gives a nonzero residual minor. -/
theorem finalResidual_minor_ne_zero_of_augmentedCore_det_ne_zero
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    {K : Kernel α β 𝕜} (run : Run K)
    (rows : γ → α) (columns : γ → β)
    (augmented_ne : (run.augmentedCore γ rows columns).det ≠ 0) :
    Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) ≠ 0 := by
  intro residual_zero
  apply augmented_ne
  rw [augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor,
    residual_zero, mul_zero]

end Run

/--
Strict sign regularity transfers an exact orientation-aware sign to every
minor of a final residual.
-/
theorem strictSignRegular_finalResidual_minor_oriented_pos
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows : γ → α) (columns : γ → β)
    (rows_injective : Function.Injective (run.augmentedRow γ rows))
    (columns_injective : Function.Injective (run.augmentedColumn γ columns)) :
    0 < signature (Fintype.card (run.AugmentedIndex γ)) *
        run.augmentedRowOrientation γ rows *
        run.augmentedColumnOrientation γ columns *
        run.selectedCore.det *
        Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) := by
  let e := Fintype.equivFin (run.AugmentedIndex γ)
  let augmentedRows : Fin (Fintype.card (run.AugmentedIndex γ)) → α :=
    run.augmentedRow γ rows ∘ e.symm
  let augmentedColumns : Fin (Fintype.card (run.AugmentedIndex γ)) → β :=
    run.augmentedColumn γ columns ∘ e.symm
  have augmentedRows_injective : Function.Injective augmentedRows :=
    rows_injective.comp e.symm.injective
  have augmentedColumns_injective : Function.Injective augmentedColumns :=
    columns_injective.comp e.symm.injective
  have oriented_pos := oriented_minor_pos_of_strictSignRegularAtOrder
    (regular (Fintype.card (run.AugmentedIndex γ)))
    augmentedRows augmentedColumns augmentedRows_injective augmentedColumns_injective
  have reindexed_det :
      (minorMatrix K augmentedRows augmentedColumns).det =
        (run.augmentedCore γ rows columns).det := by
    exact Matrix.det_reindex_self e (run.augmentedCore γ rows columns)
  rw [reindexed_det,
    Run.augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor
      (𝕜 := ℝ) γ run rows columns] at oriented_pos
  simpa [Run.augmentedRowOrientation, Run.augmentedColumnOrientation,
    augmentedRows, augmentedColumns, e, mul_assoc] using oriented_pos

/--
Strict sign regularity makes every residual minor nonsingular whenever the
combined selected and border coordinates are injective.
-/
theorem strictSignRegular_finalResidual_minor_ne_zero
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows : γ → α) (columns : γ → β)
    (rows_injective : Function.Injective (run.augmentedRow γ rows))
    (columns_injective : Function.Injective (run.augmentedColumn γ columns)) :
    Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) ≠ 0 := by
  let e := Fintype.equivFin (run.AugmentedIndex γ)
  let augmentedRows : Fin (Fintype.card (run.AugmentedIndex γ)) → α :=
    run.augmentedRow γ rows ∘ e.symm
  let augmentedColumns : Fin (Fintype.card (run.AugmentedIndex γ)) → β :=
    run.augmentedColumn γ columns ∘ e.symm
  have augmentedRows_injective : Function.Injective augmentedRows :=
    rows_injective.comp e.symm.injective
  have augmentedColumns_injective : Function.Injective augmentedColumns :=
    columns_injective.comp e.symm.injective
  have oriented_pos := oriented_minor_pos_of_strictSignRegularAtOrder
    (regular (Fintype.card (run.AugmentedIndex γ)))
    augmentedRows augmentedColumns augmentedRows_injective augmentedColumns_injective
  have reindexed_det :
      (minorMatrix K augmentedRows augmentedColumns).det =
        (run.augmentedCore γ rows columns).det := by
    exact Matrix.det_reindex_self e (run.augmentedCore γ rows columns)
  apply run.finalResidual_minor_ne_zero_of_augmentedCore_det_ne_zero γ rows columns
  intro augmented_zero
  rw [reindexed_det, augmented_zero, mul_zero] at oriented_pos
  exact lt_irrefl 0 oriented_pos

end GECP
end GECPKernelStructure
