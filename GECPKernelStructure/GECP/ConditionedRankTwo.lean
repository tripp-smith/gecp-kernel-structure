import GECPKernelStructure.GECP.CompanionPlucker
import GECPKernelStructure.Matrix.RankTwoRatio

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

namespace Run

/-- A conditioned rank-two matched coordinate is exactly the product of the
two corresponding original feature-augmented Plucker coordinates after
multiplication by the square of the selected-core determinant. -/
theorem selectedCore_sq_mul_conditionedRankTwoMatchedCoordinate
    {α : Type u} {β : Type v} {δ : Type w} {K : Kernel α β ℝ}
    (run : Run K) (rows : Fin 2 → α) (columns : Fin 2 → β)
    (features : Fin 2 → δ) (L : Matrix α δ ℝ) (R : Matrix δ β ℝ) :
    run.selectedCore.det ^ 2 *
        Matrix.rankTwoMatchedCoordinate
          (fun feature ↦ run.conditionRows L (rows 0) feature)
          (fun feature ↦ run.conditionRows L (rows 1) feature)
          (fun feature ↦ run.conditionColumns R feature (columns 0))
          (fun feature ↦ run.conditionColumns R feature (columns 1))
          (features 0) (features 1) =
      (run.leftFeatureCore rows features L).det *
        (run.rightFeatureCore features columns R).det := by
  have matched := run.leftFeatureCore_det_mul_rightFeatureCore_det
    rows columns features L R
  rw [matched]
  rw [Matrix.det_fin_two
    (fun i j : Fin 2 ↦ run.conditionRows L (rows i) (features j))]
  rw [Matrix.det_fin_two
    (fun i j : Fin 2 ↦ run.conditionColumns R (features i) (columns j))]
  unfold Matrix.rankTwoMatchedCoordinate
  ring

/-- The same Plucker bridge with an arbitrary common orientation scalar. -/
theorem selectedCore_sq_mul_orientedConditionedRankTwoMatchedCoordinate
    {α : Type u} {β : Type v} {δ : Type w} {K : Kernel α β ℝ}
    (run : Run K) (rows : Fin 2 → α) (columns : Fin 2 → β)
    (features : Fin 2 → δ) (L : Matrix α δ ℝ) (R : Matrix δ β ℝ)
    (orientation : ℝ) :
    run.selectedCore.det ^ 2 *
        (orientation * Matrix.rankTwoMatchedCoordinate
          (fun feature ↦ run.conditionRows L (rows 0) feature)
          (fun feature ↦ run.conditionRows L (rows 1) feature)
          (fun feature ↦ run.conditionColumns R feature (columns 0))
          (fun feature ↦ run.conditionColumns R feature (columns 1))
          (features 0) (features 1)) =
      orientation * ((run.leftFeatureCore rows features L).det *
        (run.rightFeatureCore features columns R).det) := by
  calc
    _ = orientation * (run.selectedCore.det ^ 2 *
        Matrix.rankTwoMatchedCoordinate
          (fun feature ↦ run.conditionRows L (rows 0) feature)
          (fun feature ↦ run.conditionRows L (rows 1) feature)
          (fun feature ↦ run.conditionColumns R feature (columns 0))
          (fun feature ↦ run.conditionColumns R feature (columns 1))
          (features 0) (features 1)) := by ring
    _ = _ := by
      rw [run.selectedCore_sq_mul_conditionedRankTwoMatchedCoordinate]

/-- Because a successful run has nonzero selected core, oriented
nonnegativity of a conditioned matched coordinate is equivalent to oriented
nonnegativity of the original augmented Plucker product. -/
theorem orientedConditionedRankTwoMatchedCoordinate_nonneg_iff
    {α : Type u} {β : Type v} {δ : Type w} {K : Kernel α β ℝ}
    (run : Run K) (rows : Fin 2 → α) (columns : Fin 2 → β)
    (features : Fin 2 → δ) (L : Matrix α δ ℝ) (R : Matrix δ β ℝ)
    (orientation : ℝ) :
    0 ≤ orientation * Matrix.rankTwoMatchedCoordinate
        (fun feature ↦ run.conditionRows L (rows 0) feature)
        (fun feature ↦ run.conditionRows L (rows 1) feature)
        (fun feature ↦ run.conditionColumns R feature (columns 0))
        (fun feature ↦ run.conditionColumns R feature (columns 1))
        (features 0) (features 1) ↔
      0 ≤ orientation * ((run.leftFeatureCore rows features L).det *
        (run.rightFeatureCore features columns R).det) := by
  have identity := run.selectedCore_sq_mul_orientedConditionedRankTwoMatchedCoordinate
    rows columns features L R orientation
  have core_sq_pos : 0 < run.selectedCore.det ^ 2 :=
    sq_pos_of_ne_zero run.selectedCore_det_ne
  constructor
  · intro conditioned_nonneg
    rw [← identity]
    exact mul_nonneg core_sq_pos.le conditioned_nonneg
  · intro augmented_nonneg
    rw [← identity] at augmented_nonneg
    apply nonneg_of_mul_nonneg_left _ core_sq_pos
    simpa only [mul_comm] using augmented_nonneg

/-- The unoriented sign bridge is the unit-orientation specialization. -/
theorem conditionedRankTwoMatchedCoordinate_nonneg_iff
    {α : Type u} {β : Type v} {δ : Type w} {K : Kernel α β ℝ}
    (run : Run K) (rows : Fin 2 → α) (columns : Fin 2 → β)
    (features : Fin 2 → δ) (L : Matrix α δ ℝ) (R : Matrix δ β ℝ) :
    0 ≤ Matrix.rankTwoMatchedCoordinate
        (fun feature ↦ run.conditionRows L (rows 0) feature)
        (fun feature ↦ run.conditionRows L (rows 1) feature)
        (fun feature ↦ run.conditionColumns R feature (columns 0))
        (fun feature ↦ run.conditionColumns R feature (columns 1))
        (features 0) (features 1) ↔
      0 ≤ (run.leftFeatureCore rows features L).det *
        (run.rightFeatureCore features columns R).det := by
  simpa using run.orientedConditionedRankTwoMatchedCoordinate_nonneg_iff
    rows columns features L R 1

end Run
end GECP
end GECPKernelStructure
