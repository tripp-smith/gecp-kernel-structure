import GECPKernelStructure.GECP.Definitions
import Mathlib.Data.Real.Basic

namespace GECPKernelStructure

universe u v

namespace GECP

/-- A pivot belongs to the chosen domain and maximizes the residual there. -/
def CompletePivotOn {α : Type u} {β : Type v} (rowDomain : α → Prop)
    (columnDomain : β → Prop) (R : Kernel α β ℝ) (row : α) (column : β) : Prop :=
  rowDomain row ∧ columnDomain column ∧
    ∀ x y, rowDomain x → columnDomain y → abs (R x y) ≤ abs (R row column)

namespace Run

/-- Every pivot in a successful finite run is complete on the chosen domain. -/
def CompleteOn {α : Type u} {β : Type v} (rowDomain : α → Prop)
    (columnDomain : β → Prop) : {K : Kernel α β ℝ} → Run K → Prop
  | _, .nil _ => True
  | K, .step row column _ tail =>
      CompletePivotOn rowDomain columnDomain K row column ∧
        CompleteOn rowDomain columnDomain tail

end Run

end GECP
end GECPKernelStructure
