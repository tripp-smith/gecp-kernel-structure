import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

namespace Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- A matrix pencil with common left and right factors and diagonal middle
coordinates has the expected product determinant.  This division-free form is
the generalized-spectrum interface used for conditioned feature pencils. -/
theorem det_add_smul_eq_prod_of_common_diagonal
    (low high left right : Matrix n n ℝ) (a b : n → ℝ) (weight : ℝ)
    (low_eq : low = left * diagonal a * right)
    (high_eq : high = left * diagonal b * right) :
    (low + weight • high).det =
      (left.det * right.det) * ∏ i, (a i + weight * b i) := by
  rw [low_eq, high_eq]
  have factorization :
      left * diagonal a * right + weight • (left * diagonal b * right) =
        left * (diagonal a + weight • diagonal b) * right := by
    rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul]
  rw [factorization, Matrix.det_mul, Matrix.det_mul]
  have middle_eq :
      diagonal a + weight • diagonal b =
        diagonal (fun i => a i + weight * b i) := by
    rw [← Matrix.diagonal_smul, Matrix.diagonal_add]
    rfl
  rw [middle_eq, Matrix.det_diagonal]
  ring

/-- Nonnegative common-factor orientation and nonnegative diagonal
generalized-spectrum coordinates certify a nonnegative determinant at every
nonnegative pencil weight. -/
theorem det_add_smul_nonneg_of_common_diagonal
    (low high left right : Matrix n n ℝ) (a b : n → ℝ) (weight : ℝ)
    (low_eq : low = left * diagonal a * right)
    (high_eq : high = left * diagonal b * right)
    (orientation_nonneg : 0 ≤ left.det * right.det)
    (a_nonneg : ∀ i, 0 ≤ a i) (b_nonneg : ∀ i, 0 ≤ b i)
    (weight_nonneg : 0 ≤ weight) :
    0 ≤ (low + weight • high).det := by
  rw [det_add_smul_eq_prod_of_common_diagonal low high left right a b weight
    low_eq high_eq]
  apply mul_nonneg orientation_nonneg
  apply Finset.prod_nonneg
  intro i _
  exact add_nonneg (a_nonneg i) (mul_nonneg weight_nonneg (b_nonneg i))

/-- If the common-factor orientation and every low generalized-spectrum
coordinate are positive, the same pencil determinant is strictly positive at
every nonnegative weight. -/
theorem det_add_smul_pos_of_common_diagonal_low_pos
    (low high left right : Matrix n n ℝ) (a b : n → ℝ) (weight : ℝ)
    (low_eq : low = left * diagonal a * right)
    (high_eq : high = left * diagonal b * right)
    (orientation_pos : 0 < left.det * right.det)
    (a_pos : ∀ i, 0 < a i) (b_nonneg : ∀ i, 0 ≤ b i)
    (weight_nonneg : 0 ≤ weight) :
    0 < (low + weight • high).det := by
  rw [det_add_smul_eq_prod_of_common_diagonal low high left right a b weight
    low_eq high_eq]
  apply mul_pos orientation_pos
  apply Finset.prod_pos
  intro i _
  exact add_pos_of_pos_of_nonneg (a_pos i)
    (mul_nonneg weight_nonneg (b_nonneg i))

/-- At a positive weight, strict positivity also follows when each diagonal
coordinate is positive in at least one of the two parts. -/
theorem det_add_smul_pos_of_common_diagonal
    (low high left right : Matrix n n ℝ) (a b : n → ℝ) (weight : ℝ)
    (low_eq : low = left * diagonal a * right)
    (high_eq : high = left * diagonal b * right)
    (orientation_pos : 0 < left.det * right.det)
    (a_nonneg : ∀ i, 0 ≤ a i) (b_nonneg : ∀ i, 0 ≤ b i)
    (coordinate_pos : ∀ i, 0 < a i ∨ 0 < b i)
    (weight_pos : 0 < weight) :
    0 < (low + weight • high).det := by
  rw [det_add_smul_eq_prod_of_common_diagonal low high left right a b weight
    low_eq high_eq]
  apply mul_pos orientation_pos
  apply Finset.prod_pos
  intro i _
  rcases coordinate_pos i with hai | hbi
  · exact add_pos_of_pos_of_nonneg hai (mul_nonneg weight_pos.le (b_nonneg i))
  · exact add_pos_of_nonneg_of_pos (a_nonneg i) (mul_pos weight_pos hbi)

end Matrix
