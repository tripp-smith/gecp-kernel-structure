import GECPKernelStructure.GECP.CompanionPlucker
import GECPKernelStructure.Fermionic.TwoSidedFeatureElimination

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The left conditioned rank-eight feature minor is represented by its
original kernel/feature bordered determinant. -/
theorem positiveSeriesLowLeft_featureCore_det
    {order : ℕ} (run : Run positiveSeriesFermionicKernel)
    (rows : Fin order → PositiveArgument) (features : Fin order → Fin 8) :
    (run.leftFeatureCore rows features positiveSeriesLowLeft).det =
      run.selectedCore.det * Matrix.det (fun i j =>
        run.conditionRows positiveSeriesLowLeft (rows i) (features j)) :=
  run.leftFeatureCore_det_eq_selectedCore_det_mul_conditionRows_minor
    rows features positiveSeriesLowLeft

/-- The right conditioned rank-eight feature minor is represented by its
original feature/kernel bordered determinant. -/
theorem positiveSeriesLowRight_featureCore_det
    {order : ℕ} (run : Run positiveSeriesFermionicKernel)
    (features : Fin order → Fin 8)
    (columns : Fin order → PositiveArgument) :
    (run.rightFeatureCore features columns positiveSeriesLowRight).det =
      run.selectedCore.det * Matrix.det (fun i j =>
        run.conditionColumns positiveSeriesLowRight
          (features i) (columns j)) :=
  run.rightFeatureCore_det_eq_selectedCore_det_mul_conditionColumns_minor
    features columns positiveSeriesLowRight

/-- After multiplying by the square of the selected-core determinant, each
conditioned ordered Cauchy--Binet term is exactly a product of two original
feature-bordered determinants.  Thus the selected-core orientation cancels. -/
theorem positiveSeriesConditionedLow_core_sq_mul_signedTerm
    {order : ℕ} (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin order → PositiveArgument) (signature : ℝ)
    (q : {q : Fin order → Fin 8 // StrictMono q}) :
    run.selectedCore.det ^ 2 *
        signedOrderedCauchyBinetTerm signature
          (Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id)
          (Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
            id columns) q =
      signature *
        (run.leftFeatureCore rows q.1 positiveSeriesLowLeft).det *
          (run.rightFeatureCore q.1 columns positiveSeriesLowRight).det := by
  have matched := run.leftFeatureCore_det_mul_rightFeatureCore_det
    rows columns q.1 positiveSeriesLowLeft positiveSeriesLowRight
  unfold signedOrderedCauchyBinetTerm
  change run.selectedCore.det ^ 2 *
      (signature *
        (Matrix.det (fun i j =>
          run.conditionRows positiveSeriesLowLeft (rows i) (q.1 j)) *
        Matrix.det (fun i j =>
          run.conditionColumns positiveSeriesLowRight (q.1 i) (columns j)))) = _
  calc
    _ = signature *
        (run.selectedCore.det ^ 2 *
          (Matrix.det (fun i j =>
            run.conditionRows positiveSeriesLowLeft (rows i) (q.1 j)) *
          Matrix.det (fun i j =>
            run.conditionColumns positiveSeriesLowRight
              (q.1 i) (columns j)))) := by ring
    _ = signature *
        ((run.leftFeatureCore rows q.1 positiveSeriesLowLeft).det *
          (run.rightFeatureCore q.1 columns positiveSeriesLowRight).det) := by
      rw [matched]
    _ = _ := by ring

/-- The sign of a conditioned matched-minor term is equivalent to the sign of
the corresponding original bordered-determinant product. -/
theorem positiveSeriesConditionedLow_signedTerm_nonneg_iff
    {order : ℕ} (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin order → PositiveArgument) (signature : ℝ)
    (q : {q : Fin order → Fin 8 // StrictMono q}) :
    0 ≤ signedOrderedCauchyBinetTerm signature
        (Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id)
        (Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
          id columns) q ↔
      0 ≤ signature *
        (run.leftFeatureCore rows q.1 positiveSeriesLowLeft).det *
          (run.rightFeatureCore q.1 columns positiveSeriesLowRight).det := by
  have identity := positiveSeriesConditionedLow_core_sq_mul_signedTerm
    run rows columns signature q
  have core_sq_pos : 0 < run.selectedCore.det ^ 2 :=
    sq_pos_of_ne_zero run.selectedCore_det_ne
  constructor <;> intro h <;> nlinarith

end Fermionic
end GECPKernelStructure
