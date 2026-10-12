import GECPKernelStructure.GECP.CompanionBiConditioning
import GECPKernelStructure.Fermionic.SchurFeatureElimination
import GECPKernelStructure.Matrix.CauchyBinet

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix Finset

/-- The global positive-series rank-eight split survives a fixed GECP prefix
with both feature factors transported by the corresponding oblique
projectors. -/
theorem positiveSeriesFermionicKernel_finalResidual_biFeatureSplit
    (run : Run positiveSeriesFermionicKernel) :
    (run.finalResidual : Matrix PositiveArgument PositiveArgument ℝ) =
      run.conditionBoth positiveSeriesTailMatrix +
        run.conditionRows positiveSeriesLowLeft *
          run.conditionColumns positiveSeriesLowRight :=
  run.finalResidual_eq_conditionBoth_add_mul
    positiveSeriesTailMatrix positiveSeriesLowLeft positiveSeriesLowRight
    positiveSeriesFermionicKernel_eq_tail_add_low

/-- Sampling the two-sided split produces the exact rank-eight factorization
used by the determinant APIs; no selected-prefix dimension is added. -/
theorem positiveSeriesFermionicKernel_finalResidual_sample_biFeatureSplit
    {n : Type*} [Fintype n] [DecidableEq n]
    (run : Run positiveSeriesFermionicKernel)
    (rows columns : n → PositiveArgument) :
    Matrix.submatrix run.finalResidual rows columns =
      Matrix.submatrix (run.conditionBoth positiveSeriesTailMatrix)
          rows columns +
        Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
            id columns :=
  run.finalResidual_sample_eq_conditionBoth_add_mul
    rows columns positiveSeriesTailMatrix positiveSeriesLowLeft
    positiveSeriesLowRight positiveSeriesFermionicKernel_eq_tail_add_low

/-- The conditioned low block is an ordered sum of products of matched left
and right feature minors.  Establishing a common sign for these terms is the
remaining analytic interface suggested by the exact Phase BT search. -/
theorem positiveSeriesConditionedLow_signature_mul_det_eq_sum
    {order : ℕ} (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin order → PositiveArgument) (signature : ℝ) :
    signature *
        (Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
            id columns).det =
      ∑ q : {q : Fin order → Fin 8 // StrictMono q},
        signedOrderedCauchyBinetTerm signature
          (Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id)
          (Matrix.submatrix (run.conditionColumns positiveSeriesLowRight)
            id columns) q :=
  signature_mul_det_mul_rect_eq_sum _ _ _

end Fermionic
end GECPKernelStructure
