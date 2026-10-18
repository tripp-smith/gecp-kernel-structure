import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

namespace Matrix

/-- The explicit Laplace expansion of a `4 x 4` determinant. -/
theorem det_fin_four_expansion {R : Type*} [CommRing R]
    (A : Matrix (Fin 4) (Fin 4) R) :
    A.det =
      A 0 0 * (A 1 1 * A 2 2 * A 3 3 - A 1 1 * A 2 3 * A 3 2 -
        A 1 2 * A 2 1 * A 3 3 + A 1 2 * A 2 3 * A 3 1 +
        A 1 3 * A 2 1 * A 3 2 - A 1 3 * A 2 2 * A 3 1) -
      A 0 1 * (A 1 0 * A 2 2 * A 3 3 - A 1 0 * A 2 3 * A 3 2 -
        A 1 2 * A 2 0 * A 3 3 + A 1 2 * A 2 3 * A 3 0 +
        A 1 3 * A 2 0 * A 3 2 - A 1 3 * A 2 2 * A 3 0) +
      A 0 2 * (A 1 0 * A 2 1 * A 3 3 - A 1 0 * A 2 3 * A 3 1 -
        A 1 1 * A 2 0 * A 3 3 + A 1 1 * A 2 3 * A 3 0 +
        A 1 3 * A 2 0 * A 3 1 - A 1 3 * A 2 1 * A 3 0) -
      A 0 3 * (A 1 0 * A 2 1 * A 3 2 - A 1 0 * A 2 2 * A 3 1 -
        A 1 1 * A 2 0 * A 3 2 + A 1 1 * A 2 2 * A 3 0 +
        A 1 2 * A 2 0 * A 3 1 - A 1 2 * A 2 1 * A 3 0) := by
  have h12 : (1 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h20 : (2 : Fin 4).succAbove (0 : Fin 3) = 0 := by decide
  have h21 : (2 : Fin 4).succAbove (1 : Fin 3) = 1 := by decide
  have h22 : (2 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h30 : (3 : Fin 4).succAbove (0 : Fin 3) = 0 := by decide
  have h31 : (3 : Fin 4).succAbove (1 : Fin 3) = 1 := by decide
  have h32 : (3 : Fin 4).succAbove (2 : Fin 3) = 2 := by decide
  have hs2 : Fin.succ (2 : Fin 3) = (3 : Fin 4) := by decide
  have hc2 : Fin.castSucc (2 : Fin 3) = (2 : Fin 4) := by decide
  rw [Matrix.det_succ_row_zero, Fin.sum_univ_succ, Fin.sum_univ_succ,
    Fin.sum_univ_succ, Fin.sum_univ_succ]
  simp only [Matrix.det_fin_three]
  norm_num [h12, h20, h21, h22, h30, h31, h32, hs2, hc2]
  ring

end Matrix
