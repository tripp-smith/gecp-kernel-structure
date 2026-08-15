import GECPKernelStructure.GECP.Residual

namespace GECPKernelStructure
namespace GECP

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

/-- Every selected row vanishes in the residual after a successful run. -/
theorem gecp_interpolates_selected_rows {K : Kernel α β 𝕜} (run : Run K)
    {x : α} (hx : x ∈ run.rows) : ∀ y, run.finalResidual x y = 0 := by
  induction run with
  | nil => simp [rows] at hx
  | step row column pivot_ne tail ih =>
      simp only [rows, List.mem_cons] at hx
      rcases hx with rfl | hx
      · exact tail.finalResidual_preserves_zero_row
          (residualUpdate_selected_row _ _ _ pivot_ne)
      · exact ih hx

/-- Every selected column vanishes in the residual after a successful run. -/
theorem gecp_interpolates_selected_cols {K : Kernel α β 𝕜} (run : Run K)
    {y : β} (hy : y ∈ run.columns) : ∀ x, run.finalResidual x y = 0 := by
  induction run with
  | nil => simp [columns] at hy
  | step row column pivot_ne tail ih =>
      simp only [columns, List.mem_cons] at hy
      rcases hy with rfl | hy
      · exact tail.finalResidual_preserves_zero_column
          (residualUpdate_selected_column _ _ _ pivot_ne)
      · exact ih hy

end Run
end GECP
end GECPKernelStructure
