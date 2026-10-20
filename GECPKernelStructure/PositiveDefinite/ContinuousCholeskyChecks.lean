import GECPKernelStructure.PositiveDefinite.ContinuousCholesky
import GECPKernelStructure.SmallInstanceChecks

namespace GECPKernelStructure.PositiveDefinite

open Matrix GECP

/-- Strict SPD is certified by an actual injective Gram factor, not just positive entries. -/
theorem exactTwoByTwoRealKernel_posDef : Matrix.PosDef exactTwoByTwoRealKernel := by
  let A : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 1, 1]
  have hA : Function.Injective A.mulVec := by
    intro v w hvw
    have h0 := congrFun hvw 0
    have h1 := congrFun hvw 1
    simp [A, Matrix.mulVec, dotProduct, Fin.sum_univ_two] at h0 h1
    funext i
    fin_cases i
    · change v 0 = w 0
      exact h0
    · change v 1 = w 1
      linarith
  have h := Matrix.PosDef.conjTranspose_mul_self A hA
  have heq : Aᴴ * A = (exactTwoByTwoRealKernel : Matrix (Fin 2) (Fin 2) ℝ) := by
    ext i j
    change (∑ k : Fin 2, star (A k i) * A k j) = exactTwoByTwoRealKernel i j
    fin_cases i <;> fin_cases j <;>
      norm_num [A, Fin.sum_univ_two, exactTwoByTwoRealKernel, exactTwoByTwoKernel]
  rw [heq] at h
  exact h

/-- The existing two actual complete pivots satisfy the independently defined Cholesky rule. -/
theorem exactTwoPivotRealRun_pivotedCholesky : DiagonalMaxRun exactTwoPivotRealRun :=
  (posDef_completeRun_iff_pivotedCholesky _ exactTwoByTwoRealKernel_posDef).mp
    exactTwoPivotRealRun_complete

end GECPKernelStructure.PositiveDefinite
