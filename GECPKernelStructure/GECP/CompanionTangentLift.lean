import GECPKernelStructure.GECP.CompanionClosedForm
import GECPKernelStructure.Matrix.UpperTangent

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

namespace Run

variable {K : Kernel α β 𝕜}

/-- The Phase BX double border of a run and a finite sample is exactly the
distinguished upper-tangent minor selected by the run. -/
theorem upperTangentMatrix_submatrix_eq_doubleBorderedCompanionMatrix
    {γ : Type*} (run : Run K) (H : Matrix α β 𝕜)
    (rows : γ → α) (columns : γ → β) :
    Matrix.submatrix
        (Matrix.upperTangentMatrix
          (fun x y => K x y : Matrix α β 𝕜) H)
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns) =
      Matrix.doubleBorderedCompanionMatrix run.selectedCore
        (Matrix.submatrix run.selectedKernelColumns rows id)
        (Matrix.submatrix run.selectedKernelRows id columns)
        (run.selectedCompanionCore H)
        (Matrix.submatrix (run.selectedCompanionColumns H) rows id)
        (Matrix.submatrix (run.selectedCompanionRows H) id columns)
        (Matrix.submatrix H rows columns) := by
  ext (i | i) (j | j)
  · rcases i with i | i <;> rcases j with j | j <;> rfl
  · rcases i with i | i <;> rfl
  · rcases j with j | j <;> rfl
  · rfl

set_option linter.style.haveILetI false in
/-- Every distinguished real upper-tangent minor selected by a successful run
factors into the squared selected-core determinant and the actual recursively
conditioned companion sample. -/
theorem upperTangentMatrix_distinguishedMinor_det
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    {Kℝ : Kernel α β ℝ} (run : Run Kℝ) (H : Matrix α β ℝ)
    (rows : γ → α) (columns : γ → β) :
    (Matrix.submatrix
        (Matrix.upperTangentMatrix
          (fun x y => Kℝ x y : Matrix α β ℝ) H)
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det =
      run.selectedCore.det ^ 2 *
        (Matrix.submatrix (run.conditionBoth H) rows columns).det := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  rw [run.upperTangentMatrix_submatrix_eq_doubleBorderedCompanionMatrix]
  exact doubleBorderedCompanionMatrix_det_eq_conditionBoth_sample
    run H rows columns

end Run
end GECP
end GECPKernelStructure
