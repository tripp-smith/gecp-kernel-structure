import GECPKernelStructure.GECP.Definitions
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

/-- Finite pivot positions, retaining the recursive order of a successful run. -/
def SelectedIndex {K : Kernel α β 𝕜} : Run K → Type
  | .nil _ => Fin 0
  | .step _ _ _ tail => Fin 1 ⊕ tail.SelectedIndex

noncomputable instance selectedIndexFintype {K : Kernel α β 𝕜} (run : Run K) :
    Fintype run.SelectedIndex := by
  induction run with
  | nil =>
      dsimp [SelectedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [SelectedIndex]
      letI := ih
      infer_instance

noncomputable instance selectedIndexDecidableEq {K : Kernel α β 𝕜} (run : Run K) :
    DecidableEq run.SelectedIndex := by
  induction run with
  | nil =>
      dsimp [SelectedIndex]
      infer_instance
  | step row column pivot_ne tail ih =>
      dsimp [SelectedIndex]
      letI := ih
      infer_instance

/-- The row selected at a recursive pivot position. -/
def selectedRow {K : Kernel α β 𝕜} : (run : Run K) → run.SelectedIndex → α
  | .nil _, i => Fin.elim0 i
  | .step row _ _ tail, i => Sum.elim (fun _ => row) tail.selectedRow i

/-- The column selected at a recursive pivot position. -/
def selectedColumn {K : Kernel α β 𝕜} : (run : Run K) → run.SelectedIndex → β
  | .nil _, i => Fin.elim0 i
  | .step _ column _ tail, i => Sum.elim (fun _ => column) tail.selectedColumn i

/-- The original kernel evaluated on the ordered rows and columns selected by `run`. -/
def selectedCore {K : Kernel α β 𝕜} (run : Run K) :
    Matrix run.SelectedIndex run.SelectedIndex 𝕜 :=
  fun i j => K (run.selectedRow i) (run.selectedColumn j)

/-- The order-preserving enumeration of recursive pivot positions by `Fin`. -/
noncomputable def selectedIndexEquivFin {K : Kernel α β 𝕜} :
    (run : Run K) → run.SelectedIndex ≃ Fin run.pivots.length
  | .nil _ => Equiv.refl _
  | .step _ _ _ tail =>
      (Equiv.sumCongr (Equiv.refl (Fin 1)) tail.selectedIndexEquivFin).trans
        (finSumFinEquiv.trans (finCongr (by simp [pivots, Nat.add_comm])))

/-- The selected core with the public `Fin k` indexing promised by the SPEC. -/
noncomputable def finSelectedCore {K : Kernel α β 𝕜} (run : Run K) :
    Matrix (Fin run.pivots.length) (Fin run.pivots.length) 𝕜 :=
  Matrix.reindex run.selectedIndexEquivFin run.selectedIndexEquivFin run.selectedCore

@[simp]
theorem selectedCore_nil (K : Kernel α β 𝕜) :
    (Run.nil K).selectedCore = (0 : Matrix (Fin 0) (Fin 0) 𝕜) := by
  ext i
  exact Fin.elim0 i

/-- One exact GECP step factors the selected-core determinant into its pivot and tail core. -/
theorem selectedCore_step_det {K : Kernel α β 𝕜} (row : α) (column : β)
    (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne)) :
    (Run.step row column pivot_ne tail).selectedCore.det =
      K row column * tail.selectedCore.det := by
  classical
  let core : Matrix (Fin 1 ⊕ tail.SelectedIndex) (Fin 1 ⊕ tail.SelectedIndex) 𝕜 :=
    fun i j => K (Sum.elim (fun _ => row) tail.selectedRow i)
      (Sum.elim (fun _ => column) tail.selectedColumn j)
  have hcore : (Run.step row column pivot_ne tail).selectedCore = core := rfl
  let pivotBlock : Matrix (Fin 1) (Fin 1) 𝕜 := fun _ _ => K row column
  let topRow : Matrix (Fin 1) tail.SelectedIndex 𝕜 :=
    fun _ j => K row (tail.selectedColumn j)
  let leftColumn : Matrix tail.SelectedIndex (Fin 1) 𝕜 :=
    fun i _ => K (tail.selectedRow i) column
  let originalTail : Matrix tail.SelectedIndex tail.SelectedIndex 𝕜 :=
    fun i j => K (tail.selectedRow i) (tail.selectedColumn j)
  let multiplier : Matrix tail.SelectedIndex (Fin 1) 𝕜 :=
    fun i _ => -(K (tail.selectedRow i) column / K row column)
  let lower : Matrix (Fin 1 ⊕ tail.SelectedIndex) (Fin 1 ⊕ tail.SelectedIndex) 𝕜 :=
    Matrix.fromBlocks (1 : Matrix (Fin 1) (Fin 1) 𝕜)
      (0 : Matrix (Fin 1) tail.SelectedIndex 𝕜) multiplier
      (1 : Matrix tail.SelectedIndex tail.SelectedIndex 𝕜)
  let eliminated : Matrix (Fin 1 ⊕ tail.SelectedIndex) (Fin 1 ⊕ tail.SelectedIndex) 𝕜 :=
    Matrix.fromBlocks pivotBlock topRow
      (0 : Matrix tail.SelectedIndex (Fin 1) 𝕜) tail.selectedCore
  have hcoreBlocks : core = Matrix.fromBlocks pivotBlock topRow leftColumn originalTail := by
    ext (i | i) (j | j) <;> rfl
  have hleft : multiplier * pivotBlock + leftColumn = 0 := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, pivotBlock, leftColumn, pivot_ne]
  have htail : multiplier * topRow + originalTail = tail.selectedCore := by
    ext i j
    rw [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_succ]
    simp [multiplier, topRow, originalTail, selectedCore, residualUpdate]
    ring
  have hmul : lower * core = eliminated := by
    rw [hcoreBlocks]
    simp only [lower, Matrix.fromBlocks_multiply, Matrix.one_mul, Matrix.zero_mul,
      add_zero, hleft, htail, eliminated]
  have hlower : lower.det = 1 := by
    simp [lower, Matrix.det_fromBlocks_zero₁₂]
  have heliminated : eliminated.det = K row column * tail.selectedCore.det := by
    rw [show eliminated = Matrix.fromBlocks pivotBlock topRow
      (0 : Matrix tail.SelectedIndex (Fin 1) 𝕜) tail.selectedCore by rfl]
    rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fin_one]
  calc
    (Run.step row column pivot_ne tail).selectedCore.det = core.det := congrArg Matrix.det hcore
    _ = lower.det * core.det := by rw [hlower, one_mul]
    _ = (lower * core).det := by rw [Matrix.det_mul]
    _ = eliminated.det := congrArg Matrix.det hmul
    _ = K row column * tail.selectedCore.det := heliminated

end Run

/-- The ordered selected-core determinant is the product of the actual GECP pivots. -/
theorem gecp_core_det_eq_prod_pivots {K : Kernel α β 𝕜} (run : Run K) :
    run.finSelectedCore.det = run.pivots.prod := by
  rw [Run.finSelectedCore, Matrix.det_reindex_self]
  induction run with
  | nil => simp [Run.SelectedIndex, Run.pivots]
  | step row column pivot_ne tail ih =>
      rw [Run.selectedCore_step_det, ih]
      rfl

/-- The selected core of every successful GECP run is nonsingular. -/
theorem gecp_core_nonsingular {K : Kernel α β 𝕜} (run : Run K) :
    run.finSelectedCore.det ≠ 0 := by
  rw [gecp_core_det_eq_prod_pivots]
  induction run with
  | nil => simp [Run.pivots]
  | step row column pivot_ne tail ih =>
      simp only [Run.pivots, List.prod_cons]
      exact mul_ne_zero pivot_ne ih

end GECP
end GECPKernelStructure
