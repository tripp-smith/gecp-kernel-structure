import GECPKernelStructure.GECP.CompanionPlucker
import GECPKernelStructure.Matrix.BlockFeatureCauchyBinet

namespace GECPKernelStructure
namespace GECP

open Matrix Finset

universe u v

namespace Run

/-- After multiplying by the squared selected-core determinant, one signed
conditioned feature term is exactly the product of its two original bordered
Plucker coordinates. -/
theorem selectedCore_sq_mul_conditionedSignedTerm
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    {a b : ℕ} (run : Run K)
    (rows : Fin a → α) (columns : Fin a → β)
    (signature : ℝ) (A : Matrix α (Fin b) ℝ)
    (B : Matrix (Fin b) β ℝ)
    (q : {q : Fin a → Fin b // StrictMono q}) :
    run.selectedCore.det ^ 2 *
        signedOrderedCauchyBinetTerm signature
          (Matrix.submatrix (run.conditionRows A) rows id)
          (Matrix.submatrix (run.conditionColumns B) id columns) q =
      signature * (run.leftFeatureCore rows q.1 A).det *
        (run.rightFeatureCore q.1 columns B).det := by
  have matched := run.leftFeatureCore_det_mul_rightFeatureCore_det
    rows columns q.1 A B
  unfold signedOrderedCauchyBinetTerm
  change run.selectedCore.det ^ 2 *
      (signature *
        (Matrix.det (fun i j => run.conditionRows A (rows i) (q.1 j)) *
          Matrix.det (fun i j =>
            run.conditionColumns B (q.1 i) (columns j)))) = _
  calc
    _ = signature *
        (run.selectedCore.det ^ 2 *
          (Matrix.det (fun i j => run.conditionRows A (rows i) (q.1 j)) *
            Matrix.det (fun i j =>
              run.conditionColumns B (q.1 i) (columns j)))) := by ring
    _ = signature *
        ((run.leftFeatureCore rows q.1 A).det *
          (run.rightFeatureCore q.1 columns B).det) := by rw [matched]
    _ = _ := by ring

/-- A fixed tail-count coefficient of the conditioned feature determinant is
the corresponding grouped sum of original bordered Plucker products.  The
identity preserves cancellations inside each group. -/
theorem selectedCore_sq_mul_conditionedBlockFeatureCoefficient
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    {a b : ℕ} (run : Run K)
    (rows : Fin a → α) (columns : Fin a → β)
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix α (Fin b) ℝ) (B : Matrix (Fin b) β ℝ)
    (degree : ℕ) :
    run.selectedCore.det ^ 2 *
        blockFeatureCoefficient signature high
          (Matrix.submatrix (run.conditionRows A) rows id)
          (Matrix.submatrix (run.conditionColumns B) id columns) degree =
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        if highFeatureCount high q.1 = degree then
          signature * (run.leftFeatureCore rows q.1 A).det *
            (run.rightFeatureCore q.1 columns B).det
        else 0 := by
  unfold blockFeatureCoefficient
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro q _
  by_cases hdegree : highFeatureCount high q.1 = degree
  · rw [if_pos hdegree, if_pos hdegree,
      run.selectedCore_sq_mul_conditionedSignedTerm]
  · rw [if_neg hdegree, if_neg hdegree, mul_zero]

/-- The whole conditioned contiguous-step determinant is the weighted sum of
the original bordered Plucker products, with the exponent equal to the number
of selected tail coordinates. -/
theorem selectedCore_sq_mul_conditionedStep_det_eq_sum
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    {a b : ℕ} (run : Run K)
    (rows : Fin a → α) (columns : Fin a → β)
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (lambda : ℝ) (A : Matrix α (Fin b) ℝ)
    (B : Matrix (Fin b) β ℝ) :
    run.selectedCore.det ^ 2 *
        (signature *
          (weightedFeatureLeft high lambda
              (Matrix.submatrix (run.conditionRows A) rows id) *
            Matrix.submatrix (run.conditionColumns B) id columns).det) =
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        lambda ^ highFeatureCount high q.1 *
          (signature * (run.leftFeatureCore rows q.1 A).det *
            (run.rightFeatureCore q.1 columns B).det) := by
  rw [signature_mul_det_weightedFeatureProduct_eq_sum, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro q _
  calc
    run.selectedCore.det ^ 2 *
        (lambda ^ highFeatureCount high q.1 *
          signedOrderedCauchyBinetTerm signature
            (Matrix.submatrix (run.conditionRows A) rows id)
            (Matrix.submatrix (run.conditionColumns B) id columns) q) =
      lambda ^ highFeatureCount high q.1 *
        (run.selectedCore.det ^ 2 *
          signedOrderedCauchyBinetTerm signature
            (Matrix.submatrix (run.conditionRows A) rows id)
            (Matrix.submatrix (run.conditionColumns B) id columns) q) := by ring
    _ = _ := by
      rw [run.selectedCore_sq_mul_conditionedSignedTerm]

/-- Nonnegativity of one grouped bordered pairing is exactly sufficient for
nonnegativity of the corresponding conditioned block-feature coefficient. -/
theorem conditionedBlockFeatureCoefficient_nonneg_of_grouped
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    {a b : ℕ} (run : Run K)
    (rows : Fin a → α) (columns : Fin a → β)
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix α (Fin b) ℝ) (B : Matrix (Fin b) β ℝ)
    (degree : ℕ)
    (grouped_nonneg : 0 ≤
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        if highFeatureCount high q.1 = degree then
          signature * (run.leftFeatureCore rows q.1 A).det *
            (run.rightFeatureCore q.1 columns B).det
        else 0) :
    0 ≤ blockFeatureCoefficient signature high
      (Matrix.submatrix (run.conditionRows A) rows id)
      (Matrix.submatrix (run.conditionColumns B) id columns) degree := by
  have identity := run.selectedCore_sq_mul_conditionedBlockFeatureCoefficient
    rows columns signature high A B degree
  have core_sq_pos : 0 < run.selectedCore.det ^ 2 :=
    sq_pos_of_ne_zero run.selectedCore_det_ne
  rw [← identity] at grouped_nonneg
  apply nonneg_of_mul_nonneg_left _ core_sq_pos
  simpa only [mul_comm] using grouped_nonneg

end Run
end GECP
end GECPKernelStructure
