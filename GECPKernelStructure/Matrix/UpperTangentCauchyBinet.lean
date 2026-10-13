import GECPKernelStructure.Matrix.TripleCauchyBinet
import GECPKernelStructure.Matrix.UpperTangent

namespace Matrix

open Finset

/-- The sampled left block-diagonal feature factor of an upper-tangent lift,
with the two feature layers put in the canonical finite order. -/
def upperTangentSampleLeft {a b : ℕ} {row : Type*}
    (A : Matrix row (Fin b) ℝ) (rows : Fin a → row ⊕ row) :
    Matrix (Fin a) (Fin (b + b)) ℝ :=
  (Matrix.fromBlocks A 0 0 A).submatrix rows finSumFinEquiv.symm

/-- The upper-unitriangular middle factor of an upper-tangent lift, reindexed
from two sum layers to one canonically ordered finite feature space. -/
def upperTangentFeatureMiddle {b : ℕ} (W : Matrix (Fin b) (Fin b) ℝ) :
    Matrix (Fin (b + b)) (Fin (b + b)) ℝ :=
  Matrix.reindex finSumFinEquiv finSumFinEquiv
    (upperTangentMatrix (1 : Matrix (Fin b) (Fin b) ℝ) W)

/-- The sampled right block-diagonal feature factor of an upper-tangent lift,
with the two feature layers put in the canonical finite order. -/
def upperTangentSampleRight {a b : ℕ} {column : Type*}
    (B : Matrix (Fin b) column ℝ) (columns : Fin a → column ⊕ column) :
    Matrix (Fin (b + b)) (Fin a) ℝ :=
  (Matrix.fromBlocks B 0 0 B).submatrix finSumFinEquiv.symm columns

/-- Sampling commutes with the three-factor upper-tangent feature
factorization after the two feature layers are canonically reindexed. -/
theorem submatrix_upperTangentMatrix_mul_factorization
    {a b : ℕ} {row column : Type*}
    (A : Matrix row (Fin b) ℝ) (W : Matrix (Fin b) (Fin b) ℝ)
    (B : Matrix (Fin b) column ℝ)
    (rows : Fin a → row ⊕ row) (columns : Fin a → column ⊕ column) :
    Matrix.submatrix (upperTangentMatrix (A * B) (A * W * B)) rows columns =
      upperTangentSampleLeft A rows * upperTangentFeatureMiddle W *
        upperTangentSampleRight B columns := by
  let L : Matrix (row ⊕ row) (Fin b ⊕ Fin b) ℝ :=
    Matrix.fromBlocks A 0 0 A
  let C : Matrix (Fin b ⊕ Fin b) (Fin b ⊕ Fin b) ℝ :=
    upperTangentMatrix (1 : Matrix (Fin b) (Fin b) ℝ) W
  let R : Matrix (Fin b ⊕ Fin b) (column ⊕ column) ℝ :=
    Matrix.fromBlocks B 0 0 B
  have lift_factorization :
      upperTangentMatrix (A * B) (A * W * B) = L * C * R := by
    simpa only [L, C, R] using upperTangentMatrix_mul_factorization A W B
  have left_middle :
      L.submatrix rows finSumFinEquiv.symm *
          C.submatrix finSumFinEquiv.symm finSumFinEquiv.symm =
        (L * C).submatrix rows finSumFinEquiv.symm :=
    Matrix.submatrix_mul_equiv L C rows finSumFinEquiv.symm
      finSumFinEquiv.symm
  have full_product :
      (L * C).submatrix rows finSumFinEquiv.symm *
          R.submatrix finSumFinEquiv.symm columns =
        (L * C * R).submatrix rows columns :=
    Matrix.submatrix_mul_equiv (L * C) R rows finSumFinEquiv.symm columns
  rw [lift_factorization]
  simp only [upperTangentSampleLeft, upperTangentFeatureMiddle,
    upperTangentSampleRight, Matrix.reindex_apply, L, C, R]
  rw [left_middle, full_product]

/-- The distinguished finite upper-tangent sample has an exact ordered
double-flag expansion.  Its summands expose one minor from each of the left,
middle, and right tangent factors. -/
theorem signature_mul_det_submatrix_upperTangentMatrix_eq_sum
    {a b : ℕ} {row column : Type*}
    (signature : ℝ) (A : Matrix row (Fin b) ℝ)
    (W : Matrix (Fin b) (Fin b) ℝ) (B : Matrix (Fin b) column ℝ)
    (rows : Fin a → row ⊕ row) (columns : Fin a → column ⊕ column) :
    signature *
        (Matrix.submatrix (upperTangentMatrix (A * B) (A * W * B))
          rows columns).det =
      ∑ q : {q : Fin a → Fin (b + b) // StrictMono q},
        ∑ r : {r : Fin a → Fin (b + b) // StrictMono r},
          signedOrderedTripleCauchyBinetTerm signature
            (upperTangentSampleLeft A rows)
            (upperTangentFeatureMiddle W)
            (upperTangentSampleRight B columns) q r := by
  rw [submatrix_upperTangentMatrix_mul_factorization]
  exact signature_mul_det_mul_mul_rect_eq_sum _ _ _ _

/-- A modular flag-minor interface for the upper-tangent sign problem.  It
requires only separate orientation certificates for the sampled left factor,
the unitriangular middle factor, and the sampled right factor. -/
theorem signature_mul_det_submatrix_upperTangentMatrix_nonneg_of_flag_orientations
    {a b : ℕ} {row column : Type*}
    (signature : ℝ) (A : Matrix row (Fin b) ℝ)
    (W : Matrix (Fin b) (Fin b) ℝ) (B : Matrix (Fin b) column ℝ)
    (rows : Fin a → row ⊕ row) (columns : Fin a → column ⊕ column)
    (leftOrientation :
      {q : Fin a → Fin (b + b) // StrictMono q} → ℝ)
    (middleOrientation :
      {q : Fin a → Fin (b + b) // StrictMono q} →
        {r : Fin a → Fin (b + b) // StrictMono r} → ℝ)
    (rightOrientation :
      {r : Fin a → Fin (b + b) // StrictMono r} → ℝ)
    (orientation_mul : ∀ q r, signature =
      leftOrientation q * middleOrientation q r * rightOrientation r)
    (left_nonneg : ∀ q, 0 ≤ leftOrientation q *
      ((upperTangentSampleLeft A rows).submatrix id q.1).det)
    (middle_nonneg : ∀ q r, 0 ≤ middleOrientation q r *
      ((upperTangentFeatureMiddle W).submatrix q.1 r.1).det)
    (right_nonneg : ∀ r, 0 ≤ rightOrientation r *
      ((upperTangentSampleRight B columns).submatrix r.1 id).det) :
    0 ≤ signature *
      (Matrix.submatrix (upperTangentMatrix (A * B) (A * W * B))
        rows columns).det := by
  rw [submatrix_upperTangentMatrix_mul_factorization]
  exact signature_mul_det_mul_mul_rect_nonneg_of_factor_orientations
    signature (upperTangentSampleLeft A rows)
      (upperTangentFeatureMiddle W) (upperTangentSampleRight B columns)
      leftOrientation middleOrientation rightOrientation orientation_mul
      left_nonneg middle_nonneg right_nonneg

end Matrix
