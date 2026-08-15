import Mathlib.Data.Rat.Cast.Order

namespace GECPKernelStructure

/-- Exact positivity of the first nontrivial geometric surrogate minor. -/
theorem geometricTwoByTwo_det_pos {q : ℚ} (q_pos : 0 < q) (q_lt_one : q < 1) :
    0 < 1 - q ^ 2 := by
  apply sub_pos.mpr
  calc
    q ^ 2 = q * q := pow_two q
    _ < 1 * q := mul_lt_mul_of_pos_right q_lt_one q_pos
    _ < 1 * 1 := mul_lt_mul_of_pos_left q_lt_one zero_lt_one
    _ = 1 := one_mul 1

end GECPKernelStructure
