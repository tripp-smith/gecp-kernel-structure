import GECPKernelStructure.GECP.Definitions

namespace GECPKernelStructure
namespace GECP

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

@[simp]
theorem residualUpdate_selected_row (R : Kernel α β 𝕜) (row : α) (column : β)
    (pivot_ne : R row column ≠ 0) (y : β) :
    residualUpdate R row column pivot_ne row y = 0 := by
  simp [residualUpdate, pivot_ne]

@[simp]
theorem residualUpdate_selected_column (R : Kernel α β 𝕜) (row : α) (column : β)
    (pivot_ne : R row column ≠ 0) (x : α) :
    residualUpdate R row column pivot_ne x column = 0 := by
  simp [residualUpdate, pivot_ne]

theorem residualUpdate_preserves_zero_row (R : Kernel α β 𝕜)
    (row : α) (column : β) (pivot_ne : R row column ≠ 0) {x : α}
    (hzero : ∀ y, R x y = 0) :
    ∀ y, residualUpdate R row column pivot_ne x y = 0 := by
  intro y
  simp [residualUpdate, hzero]

theorem residualUpdate_preserves_zero_column (R : Kernel α β 𝕜)
    (row : α) (column : β) (pivot_ne : R row column ≠ 0) {y : β}
    (hzero : ∀ x, R x y = 0) :
    ∀ x, residualUpdate R row column pivot_ne x y = 0 := by
  intro x
  simp [residualUpdate, hzero]

namespace Run

theorem finalResidual_preserves_zero_row {K : Kernel α β 𝕜} (run : Run K) {x : α}
    (hzero : ∀ y, K x y = 0) : ∀ y, run.finalResidual x y = 0 := by
  induction run with
  | nil => exact hzero
  | step row column pivot_ne tail ih =>
      exact ih (residualUpdate_preserves_zero_row _ row column pivot_ne hzero)

theorem finalResidual_preserves_zero_column {K : Kernel α β 𝕜} (run : Run K) {y : β}
    (hzero : ∀ x, K x y = 0) : ∀ x, run.finalResidual x y = 0 := by
  induction run with
  | nil => exact hzero
  | step row column pivot_ne tail ih =>
      exact ih (residualUpdate_preserves_zero_column _ row column pivot_ne hzero)

end Run
end GECP
end GECPKernelStructure
