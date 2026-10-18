import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Topology.Instances.Matrix
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Tactic

namespace Matrix

open Set

/-- A continuous real path that never vanishes on the unit interval has
endpoints with the same strict sign. This sign-locking form avoids choosing
which of the two signs occurs. -/
theorem endpoint_mul_pos_of_continuousOn_of_ne_zero
    (f : ℝ → ℝ) (continuous : ContinuousOn f (Icc 0 1))
    (nonzero : ∀ t ∈ Icc (0 : ℝ) 1, f t ≠ 0) :
    0 < f 0 * f 1 := by
  have zero_mem : (0 : ℝ) ∈ Icc 0 1 := by constructor <;> norm_num
  have one_mem : (1 : ℝ) ∈ Icc 0 1 := by constructor <;> norm_num
  have f_zero_ne : f 0 ≠ 0 := nonzero 0 zero_mem
  have f_one_ne : f 1 ≠ 0 := nonzero 1 one_mem
  by_contra product_not_pos
  have product_le : f 0 * f 1 ≤ 0 := le_of_not_gt product_not_pos
  have product_ne : f 0 * f 1 ≠ 0 := mul_ne_zero f_zero_ne f_one_ne
  have product_neg : f 0 * f 1 < 0 := lt_of_le_of_ne product_le product_ne
  have between : (0 : ℝ) ∈ uIcc (f 0) (f 1) := by
    rw [Set.mem_uIcc]
    rcases mul_neg_iff.mp product_neg with signs | signs
    · exact Or.inr ⟨signs.2.le, signs.1.le⟩
    · exact Or.inl ⟨signs.1.le, signs.2.le⟩
  have continuous_uIcc : ContinuousOn f (uIcc (0 : ℝ) 1) := by
    simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using continuous
  obtain ⟨t, t_mem, path_zero⟩ :=
    (intermediate_value_uIcc continuous_uIcc) between
  have t_mem_unit : t ∈ Icc (0 : ℝ) 1 := by
    simpa [uIcc_of_le (by norm_num : (0 : ℝ) ≤ 1)] using t_mem
  exact nonzero t t_mem_unit path_zero

/-- Entrywise continuity and determinant nonvanishing lock the endpoint signs
of a square-matrix path. -/
theorem det_endpoint_mul_pos_of_entrywiseContinuous_of_ne_zero
    {n : Type*} [Fintype n] [DecidableEq n]
    (A : ℝ → Matrix n n ℝ)
    (continuous : ∀ i j, Continuous fun t ↦ A t i j)
    (nonzero : ∀ t ∈ Set.Icc (0 : ℝ) 1, (A t).det ≠ 0) :
    0 < (A 0).det * (A 1).det := by
  apply endpoint_mul_pos_of_continuousOn_of_ne_zero
    (fun t ↦ (A t).det) _ nonzero
  apply Continuous.continuousOn
  apply Continuous.matrix_det
  apply continuous_pi
  intro i
  apply continuous_pi
  intro j
  exact continuous i j

end Matrix
