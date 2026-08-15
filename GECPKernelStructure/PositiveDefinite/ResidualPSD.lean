import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Algebra.Order.Star.Real

namespace GECPKernelStructure
namespace PositiveDefinite

open Matrix

universe u v

variable {m : Type u} {n : Type v}

/-- The Schur residual of a positive pivot block remains positive semidefinite. -/
theorem posDef_gecp_residual_posSemidefinite [Fintype m] [DecidableEq m]
    [Finite n] {A : Matrix m m ℝ} (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    (pivotBlock_posDef : A.PosDef) [Invertible A]
    (block_posSemidefinite : (fromBlocks A B Bᴴ D).PosSemidef) :
    (D - Bᴴ * A⁻¹ * B).PosSemidef :=
  (Matrix.PosDef.fromBlocks₁₁ B D pivotBlock_posDef).mp block_posSemidefinite

/-- The off-diagonal square of a real PSD `2 × 2` matrix is diagonally bounded. -/
private theorem finTwo_offDiagonal_sq_le (M : Matrix (Fin 2) (Fin 2) ℝ)
    (hM : M.PosSemidef) : |M 0 1| ^ 2 ≤ M 0 0 * M 1 1 := by
  have hdet : 0 ≤ M.det := hM.det_nonneg
  have hsym : M 1 0 = M 0 1 := by
    simpa using hM.isHermitian.apply 0 1
  rw [Matrix.det_fin_two, hsym] at hdet
  rw [sq_abs]
  simpa only [pow_two] using sub_nonneg.mp hdet

/-- PSD residual entries obey the sharp two-point diagonal domination inequality. -/
theorem posDef_gecp_residual_abs_le_diag {n : Type u} [Fintype n]
    (M : Matrix n n ℝ) (hM : M.PosSemidef) (i j : n) :
    |M i j| ^ 2 ≤ M i i * M j j := by
  let select : Fin 2 → n := ![i, j]
  simpa [select] using finTwo_offDiagonal_sq_le (M.submatrix select select) (hM.submatrix select)

private theorem abs_entry_le_max_diag {n : Type u} [Fintype n]
    (M : Matrix n n ℝ) (hM : M.PosSemidef) (i j : n) :
    |M i j| ≤ max (M i i) (M j j) := by
  have hii : 0 ≤ M i i := hM.diag_nonneg
  have hjj : 0 ≤ M j j := hM.diag_nonneg
  have hiMax : M i i ≤ max (M i i) (M j j) := le_max_left _ _
  have hjMax : M j j ≤ max (M i i) (M j j) := le_max_right _ _
  have hMax : 0 ≤ max (M i i) (M j j) := hii.trans hiMax
  have hproduct : M i i * M j j ≤ max (M i i) (M j j) ^ 2 := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hiMax) (sub_nonneg.mpr hjMax)]
  have hsquare := (posDef_gecp_residual_abs_le_diag M hM i j).trans hproduct
  nlinarith [sq_nonneg (|M i j| - max (M i i) (M j j)), abs_nonneg (M i j)]

/-- Bounding every diagonal entry is equivalent to bounding every absolute entry.

This bound-based formulation states equality of the diagonal and complete-pivot
suprema without choosing an arbitrary maximizer.
-/
theorem posDef_gecp_residual_sup_eq_diag_sup {n : Type u} [Fintype n]
    (M : Matrix n n ℝ) (hM : M.PosSemidef) (bound : ℝ) :
    (∀ i, M i i ≤ bound) ↔ ∀ i j, |M i j| ≤ bound := by
  constructor
  · intro hdiag i j
    exact (abs_entry_le_max_diag M hM i j).trans (max_le (hdiag i) (hdiag j))
  · intro hall i
    simpa [abs_of_nonneg hM.diag_nonneg] using hall i i

end PositiveDefinite
end GECPKernelStructure
