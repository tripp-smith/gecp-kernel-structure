import GECPKernelStructure.GECP.CompanionPlucker
import GECPKernelStructure.Matrix.RankTwoRatio

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

namespace Run

/-- Two separately reoriented scalar coordinates with the same nonzero sign
schedule have a nonnegative product after the two reorientations are
multiplied. The common schedule disappears through its positive square. -/
theorem mul_nonneg_of_common_sign_schedule
    {left right leftOrientation rightOrientation schedule : ℝ}
    (schedule_ne : schedule ≠ 0)
    (left_nonneg : 0 ≤ leftOrientation * schedule * left)
    (right_nonneg : 0 ≤ rightOrientation * schedule * right) :
    0 ≤ (leftOrientation * rightOrientation) * (left * right) := by
  have product_nonneg := mul_nonneg left_nonneg right_nonneg
  have schedule_sq_pos : 0 < schedule ^ 2 := sq_pos_of_ne_zero schedule_ne
  have scaled :
      0 ≤ schedule ^ 2 *
        ((leftOrientation * rightOrientation) * (left * right)) := by
    calc
      0 ≤ (leftOrientation * schedule * left) *
          (rightOrientation * schedule * right) := product_nonneg
      _ = schedule ^ 2 *
          ((leftOrientation * rightOrientation) * (left * right)) := by ring
  apply nonneg_of_mul_nonneg_left _ schedule_sq_pos
  simpa only [mul_comm] using scaled

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

/-- Separate augmented-minor sign certificates with a common nonzero
chirotope schedule imply the oriented conditioned matched-coordinate sign.
This is the interface consumed by a confluent generalized-Vandermonde sign
theorem. -/
theorem orientedConditionedRankTwoMatchedCoordinate_nonneg_of_common_sign_schedule
    {α : Type u} {β : Type v} {δ : Type w} {K : Kernel α β ℝ}
    (run : Run K) (rows : Fin 2 → α) (columns : Fin 2 → β)
    (features : Fin 2 → δ) (L : Matrix α δ ℝ) (R : Matrix δ β ℝ)
    (leftOrientation rightOrientation schedule : ℝ)
    (schedule_ne : schedule ≠ 0)
    (left_nonneg : 0 ≤ leftOrientation * schedule *
      (run.leftFeatureCore rows features L).det)
    (right_nonneg : 0 ≤ rightOrientation * schedule *
      (run.rightFeatureCore features columns R).det) :
    0 ≤ (leftOrientation * rightOrientation) *
      Matrix.rankTwoMatchedCoordinate
        (fun feature ↦ run.conditionRows L (rows 0) feature)
        (fun feature ↦ run.conditionRows L (rows 1) feature)
        (fun feature ↦ run.conditionColumns R feature (columns 0))
        (fun feature ↦ run.conditionColumns R feature (columns 1))
        (features 0) (features 1) := by
  apply (run.orientedConditionedRankTwoMatchedCoordinate_nonneg_iff
    rows columns features L R (leftOrientation * rightOrientation)).2
  exact mul_nonneg_of_common_sign_schedule schedule_ne left_nonneg right_nonneg

end Run
end GECP
end GECPKernelStructure
