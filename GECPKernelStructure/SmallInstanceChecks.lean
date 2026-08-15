import GECPKernelStructure.GECP.Determinant
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

namespace GECP

/-- An exact two-by-two kernel used to regression-test the dependent-run determinant theorem. -/
def exactTwoByTwoKernel : Kernel (Fin 2) (Fin 2) ℚ :=
  fun i j => !![(2 : ℚ), 1; 1, 1] i j

/-- The exact run has residual pivots `2` and `1 / 2`. -/
noncomputable def exactTwoPivotRun : Run exactTwoByTwoKernel :=
  Run.step 0 0 (by norm_num [exactTwoByTwoKernel]) <|
    Run.step 1 1 (by norm_num [residualUpdate, exactTwoByTwoKernel]) <|
      Run.nil _

/-- Independent exact instance: the selected core determinant is `2 * (1 / 2) = 1`. -/
example : exactTwoPivotRun.finSelectedCore.det = 1 := by
  rw [gecp_core_det_eq_prod_pivots]
  norm_num [exactTwoPivotRun, Run.pivots, residualUpdate, exactTwoByTwoKernel]

end GECP

end GECPKernelStructure
