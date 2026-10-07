import GECPKernelStructure.GECP.Determinant

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

/--
The selected positions of a successful run followed by one additional border
position. The recursive definition keeps the original pivot order and places
the border last without requiring a reindexing equivalence in the determinant
recurrence.
-/
def BorderedIndex {K : Kernel α β 𝕜} : Run K → Type
  | .nil _ => Fin 1
  | .step _ _ _ tail => Fin 1 ⊕ tail.BorderedIndex

noncomputable instance borderedIndexFintype {K : Kernel α β 𝕜} (run : Run K) :
    Fintype run.BorderedIndex := by
  induction run with
  | nil =>
      dsimp [BorderedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [BorderedIndex]
      letI := ih
      infer_instance

noncomputable instance borderedIndexDecidableEq {K : Kernel α β 𝕜} (run : Run K) :
    DecidableEq run.BorderedIndex := by
  induction run with
  | nil =>
      dsimp [BorderedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [BorderedIndex]
      letI := ih
      infer_instance

/-- The selected rows followed by one additional row. -/
def borderedRow {K : Kernel α β 𝕜} :
    (run : Run K) → α → run.BorderedIndex → α
  | .nil _, row, _ => row
  | .step selected _ _ tail, row, i =>
      Sum.elim (fun _ => selected) (tail.borderedRow row) i

/-- The selected columns followed by one additional column. -/
def borderedColumn {K : Kernel α β 𝕜} :
    (run : Run K) → β → run.BorderedIndex → β
  | .nil _, column, _ => column
  | .step _ selected _ tail, column, i =>
      Sum.elim (fun _ => selected) (tail.borderedColumn column) i

/--
The original kernel on the run's selected core with one additional row and
column appended after the recursive pivot positions.
-/
def borderedCore {K : Kernel α β 𝕜} (run : Run K) (row : α) (column : β) :
    Matrix run.BorderedIndex run.BorderedIndex 𝕜 :=
  fun i j => K (run.borderedRow row i) (run.borderedColumn column j)

@[simp]
theorem borderedCore_nil (K : Kernel α β 𝕜) (row : α) (column : β) :
    (Run.nil K).borderedCore row column =
      (fun _ _ : Fin 1 => K row column) := by
  rfl

@[simp]
theorem borderedCore_nil_det (K : Kernel α β 𝕜) (row : α) (column : β) :
    ((Run.nil K).borderedCore row column).det = K row column := by
  rw [borderedCore_nil]
  exact Matrix.det_fin_one _

/-- One elimination step factors the bordered determinant through the pivot. -/
theorem borderedCore_step_det {K : Kernel α β 𝕜} (row x : α) (column y : β)
    (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne)) :
    ((Run.step row column pivot_ne tail).borderedCore x y).det =
      K row column * (tail.borderedCore x y).det := by
  classical
  let core : Matrix (Fin 1 ⊕ tail.BorderedIndex) (Fin 1 ⊕ tail.BorderedIndex) 𝕜 :=
    fun i j => K (Sum.elim (fun _ => row) (tail.borderedRow x) i)
      (Sum.elim (fun _ => column) (tail.borderedColumn y) j)
  have hcore : (Run.step row column pivot_ne tail).borderedCore x y = core := rfl
  let pivotBlock : Matrix (Fin 1) (Fin 1) 𝕜 := fun _ _ => K row column
  let topRow : Matrix (Fin 1) tail.BorderedIndex 𝕜 :=
    fun _ j => K row (tail.borderedColumn y j)
  let leftColumn : Matrix tail.BorderedIndex (Fin 1) 𝕜 :=
    fun i _ => K (tail.borderedRow x i) column
  let originalTail : Matrix tail.BorderedIndex tail.BorderedIndex 𝕜 :=
    fun i j => K (tail.borderedRow x i) (tail.borderedColumn y j)
  let multiplier : Matrix tail.BorderedIndex (Fin 1) 𝕜 :=
    fun i _ => -(K (tail.borderedRow x i) column / K row column)
  let lower : Matrix (Fin 1 ⊕ tail.BorderedIndex) (Fin 1 ⊕ tail.BorderedIndex) 𝕜 :=
    Matrix.fromBlocks (1 : Matrix (Fin 1) (Fin 1) 𝕜)
      (0 : Matrix (Fin 1) tail.BorderedIndex 𝕜) multiplier
      (1 : Matrix tail.BorderedIndex tail.BorderedIndex 𝕜)
  let eliminated : Matrix (Fin 1 ⊕ tail.BorderedIndex)
      (Fin 1 ⊕ tail.BorderedIndex) 𝕜 :=
    Matrix.fromBlocks pivotBlock topRow
      (0 : Matrix tail.BorderedIndex (Fin 1) 𝕜) (tail.borderedCore x y)
  have hcoreBlocks : core = Matrix.fromBlocks pivotBlock topRow leftColumn originalTail := by
    ext (i | i) (j | j) <;> rfl
  have hleft : multiplier * pivotBlock + leftColumn = 0 := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, pivotBlock, leftColumn, pivot_ne]
  have htail : multiplier * topRow + originalTail = tail.borderedCore x y := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, topRow, originalTail, borderedCore, residualUpdate]
    ring
  have hmul : lower * core = eliminated := by
    rw [hcoreBlocks]
    simp only [lower, Matrix.fromBlocks_multiply, Matrix.one_mul, Matrix.zero_mul,
      add_zero, hleft, htail, eliminated]
  have hlower : lower.det = 1 := by
    simp [lower, Matrix.det_fromBlocks_zero₁₂]
  have heliminated : eliminated.det =
      K row column * (tail.borderedCore x y).det := by
    rw [show eliminated = Matrix.fromBlocks pivotBlock topRow
      (0 : Matrix tail.BorderedIndex (Fin 1) 𝕜) (tail.borderedCore x y) by rfl]
    rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fin_one]
  calc
    ((Run.step row column pivot_ne tail).borderedCore x y).det = core.det :=
      congrArg Matrix.det hcore
    _ = lower.det * core.det := by rw [hlower, one_mul]
    _ = (lower * core).det := by rw [Matrix.det_mul]
    _ = eliminated.det := congrArg Matrix.det hmul
    _ = K row column * (tail.borderedCore x y).det := heliminated

/--
The determinant of a selected core with one appended row and column equals the
selected-core determinant times the exact final residual at that border.
-/
theorem borderedCore_det_eq_selectedCore_det_mul_finalResidual
    {K : Kernel α β 𝕜} (run : Run K) (row : α) (column : β) :
    (run.borderedCore row column).det =
      run.selectedCore.det * run.finalResidual row column := by
  induction run with
  | nil initial =>
      rw [borderedCore_nil_det]
      have selected_det : (Run.nil initial).selectedCore.det = 1 := by
        rw [selectedCore_nil]
        exact Matrix.det_fin_zero
      rw [selected_det]
      simp [finalResidual]
  | step selectedRow selectedColumn pivot_ne tail ih =>
      rw [borderedCore_step_det, selectedCore_step_det, ih]
      simp only [finalResidual]
      ring

end Run
end GECP
end GECPKernelStructure
