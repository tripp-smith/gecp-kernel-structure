import GECPKernelStructure.GECP.CompanionTangentLift
import GECPKernelStructure.Fermionic.AmplifiedKernelEnvelope
import GECPKernelStructure.Fermionic.FiniteWeightAmplification
import GECPKernelStructure.Fermionic.SchurFeatureElimination
import GECPKernelStructure.Fermionic.TwoSidedFeatureElimination

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The positive-series amplified kernel is exactly the companion pencil whose
fixed-run conditioning appears in Phase BS. -/
theorem positiveSeriesAmplifiedFermionicKernel_eq_smul_tail_add_low
    (lambda : ℝ) :
    (positiveSeriesAmplifiedFermionicKernel lambda :
        Matrix PositiveArgument PositiveArgument ℝ) =
      lambda • positiveSeriesTailMatrix +
        positiveSeriesLowLeft * positiveSeriesLowRight := by
  ext u omega
  rw [Matrix.add_apply, Matrix.smul_apply, positiveSeriesLow_mul]
  simp only [positiveSeriesEighthOrderFermionicTailKernel,
    positiveSeriesEighthOrderTailKernel,
    positiveSeriesEighthOrderTaylorApprox, positiveSeriesAmplifiedFermionicKernel,
    positiveSeriesTailMatrix, smul_eq_mul]
  ring

/-- Conditioning the amplified kernel is precisely the two-sided Phase BS
companion pencil. -/
theorem conditionBoth_positiveSeriesAmplifiedFermionicKernel
    (run : Run positiveSeriesFermionicKernel) (lambda : ℝ) :
    run.conditionBoth (positiveSeriesAmplifiedFermionicKernel lambda) =
      lambda • run.conditionBoth positiveSeriesTailMatrix +
        run.conditionRows positiveSeriesLowLeft *
          run.conditionColumns positiveSeriesLowRight := by
  rw [positiveSeriesAmplifiedFermionicKernel_eq_smul_tail_add_low,
    run.conditionBoth_add, run.conditionBoth_smul, run.conditionBoth_mul]

/-- The distinguished tangent minor of the amplified positive-series kernel
is the squared selected-core determinant times the exact two-sided Phase BS
sample pencil. -/
theorem positiveSeriesAmplified_upperTangentMinor_det
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    (run : Run positiveSeriesFermionicKernel) (lambda : ℝ)
    (rows columns : γ → PositiveArgument) :
    (Matrix.submatrix
        (Matrix.upperTangentMatrix
          (positiveSeriesFermionicKernel :
            Matrix PositiveArgument PositiveArgument ℝ)
          (positiveSeriesAmplifiedFermionicKernel lambda))
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det =
      run.selectedCore.det ^ 2 *
        (lambda •
            Matrix.submatrix (run.conditionBoth positiveSeriesTailMatrix)
              rows columns +
          Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
            Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
              id columns).det := by
  have conditioned :=
    conditionBoth_positiveSeriesAmplifiedFermionicKernel run lambda
  have sampled := congrArg
    (fun M => Matrix.submatrix M rows columns) conditioned
  have sampled_eq :
      Matrix.submatrix
          (run.conditionBoth (positiveSeriesAmplifiedFermionicKernel lambda))
          rows columns =
        lambda •
            Matrix.submatrix (run.conditionBoth positiveSeriesTailMatrix)
              rows columns +
          Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
            Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
              id columns := by
    ext i j
    have entry := congrFun (congrFun sampled i) j
    simpa only [Matrix.submatrix_apply, Matrix.add_apply, Matrix.smul_apply,
      Matrix.mul_apply, id_eq] using entry
  calc
    (Matrix.submatrix
        (Matrix.upperTangentMatrix
          (positiveSeriesFermionicKernel :
            Matrix PositiveArgument PositiveArgument ℝ)
          (positiveSeriesAmplifiedFermionicKernel lambda))
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det =
        run.selectedCore.det ^ 2 *
          (Matrix.submatrix
            (run.conditionBoth (positiveSeriesAmplifiedFermionicKernel lambda))
            rows columns).det := by
      exact run.upperTangentMatrix_distinguishedMinor_det
        (fun u omega => positiveSeriesAmplifiedFermionicKernel lambda u omega)
        rows columns
    _ = _ := by rw [sampled_eq]

/-- Nine two-sided companion bounds imply the Phase BS order-eighteen
half-power determinant bound. -/
theorem positiveSeriesFermionicKernel_biConditioned_det_le_halfPow
    (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin 18 → PositiveArgument)
    (node_det_bound : ∀ j,
      abs (rankEightAmplificationNode j •
          Matrix.submatrix (run.conditionBoth positiveSeriesTailMatrix)
            rows columns +
        Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
            id columns).det ≤ (2 : ℝ) ^ 18) :
    abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      ((1 / 2 : ℝ) ^ 18) := by
  let conditionedTail : Matrix (Fin 18) (Fin 18) ℝ :=
    Matrix.submatrix (run.conditionBoth positiveSeriesTailMatrix) rows columns
  let conditionedLeft : Matrix (Fin 18) (Fin 8) ℝ :=
    Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id
  let conditionedRight : Matrix (Fin 8) (Fin 18) ℝ :=
    Matrix.submatrix (run.conditionColumns positiveSeriesLowRight) id columns
  have determinant_bound := Matrix.rankEight_finiteWeight_det_le_halfPow
    conditionedTail conditionedLeft conditionedRight node_det_bound
  have residual_split :=
    positiveSeriesFermionicKernel_finalResidual_sample_biFeatureSplit
      run rows columns
  change abs (Matrix.submatrix run.finalResidual rows columns).det ≤ _
  rw [residual_split]
  exact determinant_bound

/-- A normalized upper-tangent envelope at the nine Phase BS weights is
already sufficient for the order-eighteen half-power determinant bound. -/
theorem positiveSeriesFermionicKernel_det_le_halfPow_of_upperTangent
    (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin 18 → PositiveArgument)
    (node_tangent_bound : ∀ j,
      abs (Matrix.submatrix
          (Matrix.upperTangentMatrix
            (positiveSeriesFermionicKernel :
              Matrix PositiveArgument PositiveArgument ℝ)
            (positiveSeriesAmplifiedFermionicKernel
              (rankEightAmplificationNode j)))
          (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
          (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det ≤
        abs run.selectedCore.det ^ 2 * (2 : ℝ) ^ 18) :
    abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      ((1 / 2 : ℝ) ^ 18) := by
  apply positiveSeriesFermionicKernel_biConditioned_det_le_halfPow
  intro j
  have bound := node_tangent_bound j
  rw [positiveSeriesAmplified_upperTangentMinor_det, abs_mul, abs_pow] at bound
  have core_sq_pos : 0 < abs run.selectedCore.det ^ 2 :=
    pow_pos (abs_pos.mpr run.selectedCore_det_ne) 2
  exact le_of_mul_le_mul_left bound core_sq_pos

end Fermionic
end GECPKernelStructure
