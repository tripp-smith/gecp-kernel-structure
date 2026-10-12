import GECPKernelStructure.Matrix.GradedPluckerPairing
import Mathlib.LinearAlgebra.Matrix.SchurComplement

namespace Matrix

variable {p n : Type*} [Fintype p] [DecidableEq p]
  [Fintype n] [DecidableEq n]

/-- A bordered linearization of the rank-update determinant pencil
`mu • (U * V) + B`.  The identity corner makes the construction valid without
an invertibility assumption on any input block. -/
def borderedRankUpdatePencil (mu : ℝ) (U : Matrix n p ℝ)
    (V : Matrix p n ℝ) (B : Matrix n n ℝ) :
    Matrix (p ⊕ n) (p ⊕ n) ℝ :=
  Matrix.fromBlocks 1 V (-mu • U) B

/-- The bordered pencil has exactly the determinant of the corresponding
rank update. -/
theorem borderedRankUpdatePencil_det (mu : ℝ) (U : Matrix n p ℝ)
    (V : Matrix p n ℝ) (B : Matrix n n ℝ) :
    (borderedRankUpdatePencil mu U V B).det =
      (mu • (U * V) + B).det := by
  rw [borderedRankUpdatePencil, Matrix.det_fromBlocks_one₁₁]
  congr 1
  rw [Matrix.smul_mul]
  ext i j
  simp [sub_eq_add_neg, add_comm]

/-- Every evaluation of a mixed-column determinant polynomial is one
bordered determinant. -/
theorem mixedColumnPolynomial_eval_eq_borderedRankUpdatePencil
    (mu : ℝ) (U : Matrix n p ℝ) (V : Matrix p n ℝ) (B : Matrix n n ℝ) :
    (mixedColumnPolynomial (U * V) B).eval mu =
      (borderedRankUpdatePencil mu U V B).det := by
  rw [mixedColumnPolynomial_eval, borderedRankUpdatePencil_det]

/-- Phase BV's combined fixed-plus-feature determinant polynomial admits the
same bordered linearization at every weight. -/
theorem blockFeaturePolynomial_eval_eq_borderedRankUpdatePencil
    {fixed fresh features : ℕ} (mu : ℝ)
    (U : Matrix (Fin (fixed + fresh)) (Fin fixed) ℝ)
    (V : Matrix (Fin fixed) (Fin (fixed + fresh)) ℝ)
    (A : Matrix (Fin (fixed + fresh)) (Fin features) ℝ)
    (B : Matrix (Fin features) (Fin (fixed + fresh)) ℝ) :
    (blockFeaturePolynomial 1
        (isFixedCoordinate (p := fixed) (k := features))
        (combinedFeatureLeft U A) (combinedFeatureRight V B)).eval mu =
      (borderedRankUpdatePencil mu U V (A * B)).det := by
  rw [blockFeaturePolynomial_eval, weightedCombinedFeatureProduct,
    borderedRankUpdatePencil_det]

end Matrix
