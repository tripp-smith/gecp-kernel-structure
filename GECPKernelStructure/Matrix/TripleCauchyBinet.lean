import GECPKernelStructure.Matrix.CauchyBinet

namespace Matrix

open Finset

/-- One summand in the ordered Cauchy--Binet expansion of a three-factor
rectangular product.  The external signature is deliberately kept outside
the three minors so that each factor may be given its own shuffle
orientation. -/
def signedOrderedTripleCauchyBinetTerm {a b c : ℕ} (signature : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ)
    (M : Matrix (Fin b) (Fin c) ℝ)
    (B : Matrix (Fin c) (Fin a) ℝ)
    (q : {q : Fin a → Fin b // StrictMono q})
    (r : {r : Fin a → Fin c // StrictMono r}) : ℝ :=
  signature * (A.submatrix id q.1).det *
    (M.submatrix q.1 r.1).det * (B.submatrix r.1 id).det

/-- Ordered Cauchy--Binet for three rectangular factors.  Unlike applying the
two-factor formula as a black box, this exposes both intermediate flags and
all three minors in every summand. -/
theorem signature_mul_det_mul_mul_rect_eq_sum {a b c : ℕ}
    (signature : ℝ) (A : Matrix (Fin a) (Fin b) ℝ)
    (M : Matrix (Fin b) (Fin c) ℝ)
    (B : Matrix (Fin c) (Fin a) ℝ) :
    signature * (A * M * B).det =
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        ∑ r : {r : Fin a → Fin c // StrictMono r},
          signedOrderedTripleCauchyBinetTerm signature A M B q r := by
  rw [Matrix.mul_assoc]
  rw [det_mul_rect, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro q _
  have sampled_product :
      (M * B).submatrix q.1 id = M.submatrix q.1 id * B := by
    ext i j
    simp only [Matrix.submatrix_apply, Matrix.mul_apply, id_eq]
  rw [sampled_product, det_mul_rect, Finset.mul_sum, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro r _
  have middle_sample :
      (M.submatrix q.1 id).submatrix id r.1 =
        M.submatrix q.1 r.1 := by rfl
  rw [middle_sample]
  simp only [signedOrderedTripleCauchyBinetTerm]
  ring

/-- Factorwise orientations certify the sign of a triple Cauchy--Binet term
when their product is the required external orientation. -/
theorem signedOrderedTripleCauchyBinetTerm_nonneg_of_factor_orientations
    {a b c : ℕ} (signature : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ)
    (M : Matrix (Fin b) (Fin c) ℝ)
    (B : Matrix (Fin c) (Fin a) ℝ)
    (leftOrientation : {q : Fin a → Fin b // StrictMono q} → ℝ)
    (middleOrientation :
      {q : Fin a → Fin b // StrictMono q} →
        {r : Fin a → Fin c // StrictMono r} → ℝ)
    (rightOrientation : {r : Fin a → Fin c // StrictMono r} → ℝ)
    (q : {q : Fin a → Fin b // StrictMono q})
    (r : {r : Fin a → Fin c // StrictMono r})
    (orientation_mul : signature =
      leftOrientation q * middleOrientation q r * rightOrientation r)
    (left_nonneg :
      0 ≤ leftOrientation q * (A.submatrix id q.1).det)
    (middle_nonneg :
      0 ≤ middleOrientation q r * (M.submatrix q.1 r.1).det)
    (right_nonneg :
      0 ≤ rightOrientation r * (B.submatrix r.1 id).det) :
    0 ≤ signedOrderedTripleCauchyBinetTerm signature A M B q r := by
  rw [signedOrderedTripleCauchyBinetTerm, orientation_mul]
  have product_nonneg := mul_nonneg (mul_nonneg left_nonneg middle_nonneg)
    right_nonneg
  simpa only [mul_assoc, mul_left_comm, mul_comm] using product_nonneg

/-- If all exposed three-factor summands have one orientation, the determinant
of the complete rectangular product has that orientation. -/
theorem signature_mul_det_mul_mul_rect_nonneg {a b c : ℕ}
    (signature : ℝ) (A : Matrix (Fin a) (Fin b) ℝ)
    (M : Matrix (Fin b) (Fin c) ℝ)
    (B : Matrix (Fin c) (Fin a) ℝ)
    (term_nonneg :
      ∀ (q : {q : Fin a → Fin b // StrictMono q})
        (r : {r : Fin a → Fin c // StrictMono r}),
        0 ≤ signedOrderedTripleCauchyBinetTerm signature A M B q r) :
    0 ≤ signature * (A * M * B).det := by
  rw [signature_mul_det_mul_mul_rect_eq_sum]
  exact Finset.sum_nonneg fun q _ => Finset.sum_nonneg fun r _ =>
    term_nonneg q r

/-- Separate left, middle, and right flag-minor sign certificates imply the
oriented sign of the complete three-factor determinant. -/
theorem signature_mul_det_mul_mul_rect_nonneg_of_factor_orientations
    {a b c : ℕ} (signature : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ)
    (M : Matrix (Fin b) (Fin c) ℝ)
    (B : Matrix (Fin c) (Fin a) ℝ)
    (leftOrientation : {q : Fin a → Fin b // StrictMono q} → ℝ)
    (middleOrientation :
      {q : Fin a → Fin b // StrictMono q} →
        {r : Fin a → Fin c // StrictMono r} → ℝ)
    (rightOrientation : {r : Fin a → Fin c // StrictMono r} → ℝ)
    (orientation_mul : ∀ q r, signature =
      leftOrientation q * middleOrientation q r * rightOrientation r)
    (left_nonneg : ∀ q,
      0 ≤ leftOrientation q * (A.submatrix id q.1).det)
    (middle_nonneg : ∀ q r,
      0 ≤ middleOrientation q r * (M.submatrix q.1 r.1).det)
    (right_nonneg : ∀ r,
      0 ≤ rightOrientation r * (B.submatrix r.1 id).det) :
    0 ≤ signature * (A * M * B).det := by
  apply signature_mul_det_mul_mul_rect_nonneg
  intro q r
  exact signedOrderedTripleCauchyBinetTerm_nonneg_of_factor_orientations
    signature A M B leftOrientation middleOrientation rightOrientation q r
    (orientation_mul q r) (left_nonneg q) (middle_nonneg q r)
    (right_nonneg r)

end Matrix
