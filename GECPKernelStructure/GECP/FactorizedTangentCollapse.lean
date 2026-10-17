import GECPKernelStructure.GECP.CompanionTangentLift
import GECPKernelStructure.Matrix.UpperTangentCauchyBinet

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

namespace Run

/-- A distinguished upper-tangent minor of a shared feature factorization
collapses before any Cauchy--Binet expansion: all double-flag cancellation is
contained in one determinant of the separately conditioned feature factors. -/
theorem upperTangentMatrix_factorized_distinguishedMinor_det
    {α : Type u} {β : Type v} {feature : Type w} {fresh : Type*}
    [Fintype feature] [DecidableEq feature]
    [Fintype fresh] [DecidableEq fresh]
    (A : Matrix α feature ℝ) (W : Matrix feature feature ℝ)
    (B : Matrix feature β ℝ)
    (run : Run (fun x y => (A * B) x y))
    (rows : fresh → α) (columns : fresh → β) :
    (Matrix.submatrix
        (Matrix.upperTangentMatrix (A * B) (A * W * B))
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det =
      run.selectedCore.det ^ 2 *
        (Matrix.submatrix
          (run.conditionRows A * W * run.conditionColumns B)
          rows columns).det := by
  have conditioned :
      run.conditionBoth (A * W * B) =
        run.conditionRows A * W * run.conditionColumns B := by
    rw [show A * W * B = A * (W * B) by rw [Matrix.mul_assoc],
      run.conditionBoth_mul, run.conditionColumns_mul]
    rw [Matrix.mul_assoc]
  rw [run.upperTangentMatrix_distinguishedMinor_det, conditioned]

/-- The ordered double-flag expansion of a distinguished factorized tangent
minor is exactly the same single conditioned-product determinant.  An
equivalence chooses a finite order for the double-border indices; using it on
both axes introduces no determinant sign. -/
theorem signature_mul_conditionedFeature_det_eq_doubleFlagSum
    {α : Type u} {β : Type v} {fresh : Type*}
    {a b : ℕ}
    [Fintype fresh] [DecidableEq fresh]
    (signature : ℝ)
    (A : Matrix α (Fin b) ℝ) (W : Matrix (Fin b) (Fin b) ℝ)
    (B : Matrix (Fin b) β ℝ)
    (run : Run (fun x y => (A * B) x y))
    (rows : fresh → α) (columns : fresh → β)
    (e : ((run.SelectedIndex ⊕ run.SelectedIndex) ⊕ fresh) ≃ Fin a) :
    signature * run.selectedCore.det ^ 2 *
        (Matrix.submatrix
          (run.conditionRows A * W * run.conditionColumns B)
          rows columns).det =
      ∑ q : {q : Fin a → Fin (b + b) // StrictMono q},
        ∑ r : {r : Fin a → Fin (b + b) // StrictMono r},
          Matrix.signedOrderedTripleCauchyBinetTerm signature
            (Matrix.upperTangentSampleLeft A
              (fun i => Matrix.upperTangentDoubleBorderRow
                run.selectedRow rows (e.symm i)))
            (Matrix.upperTangentFeatureMiddle W)
            (Matrix.upperTangentSampleRight B
              (fun i => Matrix.upperTangentDoubleBorderColumn
                run.selectedColumn columns (e.symm i))) q r := by
  let tangentRows : Fin a → α ⊕ α := fun i =>
    Matrix.upperTangentDoubleBorderRow run.selectedRow rows (e.symm i)
  let tangentColumns : Fin a → β ⊕ β := fun i =>
    Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns (e.symm i)
  have expansion :=
    Matrix.signature_mul_det_submatrix_upperTangentMatrix_eq_sum
      signature A W B tangentRows tangentColumns
  have reindexed :
      Matrix.submatrix
          (Matrix.upperTangentMatrix (A * B) (A * W * B))
          tangentRows tangentColumns =
        Matrix.reindex e e
          (Matrix.submatrix
            (Matrix.upperTangentMatrix (A * B) (A * W * B))
            (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
            (Matrix.upperTangentDoubleBorderColumn
              run.selectedColumn columns)) := by
    rfl
  rw [reindexed, Matrix.det_reindex_self,
    run.upperTangentMatrix_factorized_distinguishedMinor_det] at expansion
  simpa only [tangentRows, tangentColumns, mul_assoc] using expansion

end Run
end GECP
end GECPKernelStructure
