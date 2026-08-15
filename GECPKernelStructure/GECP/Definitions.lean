import GECPKernelStructure.Definitions
import Mathlib.Algebra.Field.Basic

namespace GECPKernelStructure
namespace GECP

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w}

/-- One exact Gaussian-elimination residual update. -/
def residualUpdate [Field 𝕜] (R : Kernel α β 𝕜) (row : α) (column : β)
    (_pivot_ne : R row column ≠ 0) : Kernel α β 𝕜 :=
  fun x y ↦ R x y - R x column * R row y / R row column

/-- A successful ordered run carries the nonzero proof required at every pivot. -/
inductive Run [Field 𝕜] : Kernel α β 𝕜 → Type (max u v w) where
  | nil (K : Kernel α β 𝕜) : Run K
  | step {K : Kernel α β 𝕜} (row : α) (column : β)
      (pivot_ne : K row column ≠ 0)
      (tail : Run (residualUpdate K row column pivot_ne)) : Run K

namespace Run

variable [Field 𝕜] {K : Kernel α β 𝕜}

/-- The residual after all successful steps in a run. -/
noncomputable def finalResidual (run : Run K) : Kernel α β 𝕜 :=
  Run.rec (motive := fun _ _ ↦ Kernel α β 𝕜)
    (fun initial ↦ initial)
    (fun _ _ _ _ final ↦ final)
    run

/-- Ordered selected row nodes. -/
noncomputable def rows (run : Run K) : List α :=
  Run.rec (motive := fun _ _ ↦ List α)
    (fun _ ↦ [])
    (fun row _ _ _ tailRows ↦ row :: tailRows)
    run

/-- Ordered selected column nodes. -/
noncomputable def columns (run : Run K) : List β :=
  Run.rec (motive := fun _ _ ↦ List β)
    (fun _ ↦ [])
    (fun _ column _ _ tailColumns ↦ column :: tailColumns)
    run

/-- Ordered nonzero pivots, embedded in the original run. -/
noncomputable def pivots (run : Run K) : List 𝕜 :=
  Run.rec (motive := fun _ _ ↦ List 𝕜)
    (fun _ ↦ [])
    (fun {K} row column _ _ tailPivots ↦ K row column :: tailPivots)
    run

@[simp]
theorem length_rows_eq_length_columns (run : Run K) :
    run.rows.length = run.columns.length := by
  induction run with
  | nil => rfl
  | step _ _ _ _ ih =>
      change Nat.succ _ = Nat.succ _
      exact congrArg Nat.succ ih

@[simp]
theorem length_pivots_eq_length_rows (run : Run K) :
    run.pivots.length = run.rows.length := by
  induction run with
  | nil => rfl
  | step _ _ _ _ ih =>
      change Nat.succ _ = Nat.succ _
      exact congrArg Nat.succ ih

end Run
end GECP
end GECPKernelStructure
