import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

namespace Matrix

/-- The mixed coefficient in the determinant of a `2 x 2` matrix pencil. -/
noncomputable def finTwoPencilMixed (low high : Matrix (Fin 2) (Fin 2) ℝ) : ℝ :=
  low 0 0 * high 1 1 + high 0 0 * low 1 1 -
    low 0 1 * high 1 0 - high 0 1 * low 1 0

/-- The determinant of a `2 x 2` pencil is an explicit quadratic whose
constant and leading coefficients are the endpoint determinants. -/
theorem finTwoPencil_det_expansion
    (low high : Matrix (Fin 2) (Fin 2) ℝ) (weight : ℝ) :
    (low + weight • high).det =
      low.det + finTwoPencilMixed low high * weight + high.det * weight ^ 2 := by
  rw [Matrix.det_fin_two]
  simp only [Matrix.add_apply, smul_apply, smul_eq_mul]
  rw [Matrix.det_fin_two, Matrix.det_fin_two]
  unfold finTwoPencilMixed
  ring

/-- The discriminant of the fresh-order-two determinant pencil. -/
noncomputable def finTwoPencilDiscriminant
    (low high : Matrix (Fin 2) (Fin 2) ℝ) : ℝ :=
  finTwoPencilMixed low high ^ 2 - 4 * low.det * high.det

/-- The Lorentzian coefficient inequality is exactly nonnegativity of the
fresh-order-two pencil discriminant. -/
theorem finTwoPencilDiscriminant_nonneg_iff
    (low high : Matrix (Fin 2) (Fin 2) ℝ) :
    0 ≤ finTwoPencilDiscriminant low high ↔
      4 * low.det * high.det ≤ finTwoPencilMixed low high ^ 2 := by
  unfold finTwoPencilDiscriminant
  exact sub_nonneg

/-- The mixed coefficient can be recovered from the determinant at unit
pencil weight. -/
theorem finTwoPencilMixed_eq_det_add_sub
    (low high : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilMixed low high = (low + high).det - low.det - high.det := by
  have expansion := finTwoPencil_det_expansion low high 1
  norm_num at expansion
  linarith

/-- Two-sided changes of homogeneous basis scale the mixed determinant
coefficient by the product of the two basis determinants. -/
theorem finTwoPencilMixed_mul_mul
    (left low high right : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilMixed (left * low * right) (left * high * right) =
      left.det * right.det * finTwoPencilMixed low high := by
  rw [finTwoPencilMixed_eq_det_add_sub,
    finTwoPencilMixed_eq_det_add_sub]
  have add_eq :
      left * low * right + left * high * right =
        left * (low + high) * right := by
    rw [Matrix.mul_add, Matrix.add_mul]
  rw [add_eq]
  simp only [Matrix.det_mul]
  ring

/-- The discriminant scales by the square of the two-sided projective basis
factor. -/
theorem finTwoPencilDiscriminant_mul_mul
    (left low high right : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilDiscriminant (left * low * right) (left * high * right) =
      (left.det * right.det) ^ 2 * finTwoPencilDiscriminant low high := by
  unfold finTwoPencilDiscriminant
  rw [finTwoPencilMixed_mul_mul]
  simp only [Matrix.det_mul]
  ring

/-- Lorentzianity is preserved by every two-sided basis change, including a
singular degeneration. -/
theorem finTwoPencil_lorentzian_mul_mul
    (left low high right : Matrix (Fin 2) (Fin 2) ℝ)
    (lorentzian :
      4 * low.det * high.det ≤ finTwoPencilMixed low high ^ 2) :
    4 * (left * low * right).det * (left * high * right).det ≤
      finTwoPencilMixed (left * low * right) (left * high * right) ^ 2 := by
  rw [← finTwoPencilDiscriminant_nonneg_iff] at lorentzian ⊢
  rw [finTwoPencilDiscriminant_mul_mul]
  exact mul_nonneg (sq_nonneg _) lorentzian

/-- Invertible two-sided basis changes preserve and reflect the Lorentzian
discriminant inequality. -/
theorem finTwoPencil_lorentzian_mul_mul_iff
    (left low high right : Matrix (Fin 2) (Fin 2) ℝ)
    (left_det_ne : left.det ≠ 0) (right_det_ne : right.det ≠ 0) :
    4 * (left * low * right).det * (left * high * right).det ≤
        finTwoPencilMixed (left * low * right) (left * high * right) ^ 2 ↔
      4 * low.det * high.det ≤ finTwoPencilMixed low high ^ 2 := by
  rw [← finTwoPencilDiscriminant_nonneg_iff,
    ← finTwoPencilDiscriminant_nonneg_iff,
    finTwoPencilDiscriminant_mul_mul]
  have scale_sq_pos : 0 < (left.det * right.det) ^ 2 :=
    sq_pos_of_ne_zero (mul_ne_zero left_det_ne right_det_ne)
  constructor
  · intro scaled_nonneg
    exact nonneg_of_mul_nonneg_left
      (by simpa [mul_comm] using scaled_nonneg) scale_sq_pos
  · intro original_nonneg
    exact mul_nonneg scale_sq_pos.le original_nonneg

/-- A nonnegative quadratic coefficient sequence satisfying its sharp
Lorentzian discriminant inequality factors into two affine terms with
nonnegative offsets.  This is the scalar generalized-spectrum certificate. -/
theorem nonnegative_quadratic_factorization
    (constant mixed quadratic : ℝ)
    (constant_nonneg : 0 ≤ constant) (mixed_nonneg : 0 ≤ mixed)
    (quadratic_pos : 0 < quadratic)
    (discriminant_nonneg : 0 ≤ mixed ^ 2 - 4 * constant * quadratic) :
    ∃ lower upper : ℝ,
      0 ≤ lower ∧ 0 ≤ upper ∧
        ∀ weight : ℝ,
          constant + mixed * weight + quadratic * weight ^ 2 =
            quadratic * (weight + lower) * (weight + upper) := by
  let discriminant := mixed ^ 2 - 4 * constant * quadratic
  let root := Real.sqrt discriminant
  have root_nonneg : 0 ≤ root := Real.sqrt_nonneg _
  have root_sq : root ^ 2 = discriminant := by
    exact Real.sq_sqrt discriminant_nonneg
  have root_le_mixed : root ≤ mixed := by
    have product_nonneg : 0 ≤ constant * quadratic :=
      mul_nonneg constant_nonneg quadratic_pos.le
    dsimp [discriminant] at root_sq
    nlinarith
  refine ⟨(mixed - root) / (2 * quadratic),
    (mixed + root) / (2 * quadratic), ?_, ?_, ?_⟩
  · exact div_nonneg (sub_nonneg.mpr root_le_mixed)
      (mul_nonneg (by norm_num) quadratic_pos.le)
  · exact div_nonneg (add_nonneg mixed_nonneg root_nonneg)
      (mul_nonneg (by norm_num) quadratic_pos.le)
  · intro weight
    have quadratic_ne : quadratic ≠ 0 := ne_of_gt quadratic_pos
    dsimp [root, discriminant] at root_sq ⊢
    field_simp [quadratic_ne]
    ring_nf at root_sq ⊢
    nlinarith

/-- For a `2 x 2` pencil, nonnegative endpoint and mixed coefficients plus the
Lorentzian discriminant inequality construct the two nonnegative generalized
spectral factors directly. -/
theorem finTwoPencil_nonnegative_factorization
    (low high : Matrix (Fin 2) (Fin 2) ℝ)
    (low_det_nonneg : 0 ≤ low.det)
    (mixed_nonneg : 0 ≤ finTwoPencilMixed low high)
    (high_det_pos : 0 < high.det)
    (lorentzian :
      4 * low.det * high.det ≤ finTwoPencilMixed low high ^ 2) :
    ∃ lower upper : ℝ,
      0 ≤ lower ∧ 0 ≤ upper ∧
        ∀ weight : ℝ,
          (low + weight • high).det =
            high.det * (weight + lower) * (weight + upper) := by
  have discriminant_nonneg :
      0 ≤ finTwoPencilMixed low high ^ 2 - 4 * low.det * high.det := by
    linarith
  obtain ⟨lower, upper, lower_nonneg, upper_nonneg, factorization⟩ :=
    nonnegative_quadratic_factorization low.det
      (finTwoPencilMixed low high) high.det low_det_nonneg mixed_nonneg
      high_det_pos discriminant_nonneg
  refine ⟨lower, upper, lower_nonneg, upper_nonneg, ?_⟩
  intro weight
  rw [finTwoPencil_det_expansion]
  exact factorization weight

/-- The same Lorentzian hypotheses certify the determinant sign at every
nonnegative pencil weight. -/
theorem finTwoPencil_det_nonneg_of_lorentzian
    (low high : Matrix (Fin 2) (Fin 2) ℝ)
    (low_det_nonneg : 0 ≤ low.det)
    (mixed_nonneg : 0 ≤ finTwoPencilMixed low high)
    (high_det_pos : 0 < high.det)
    (lorentzian :
      4 * low.det * high.det ≤ finTwoPencilMixed low high ^ 2)
    (weight : ℝ) (weight_nonneg : 0 ≤ weight) :
    0 ≤ (low + weight • high).det := by
  obtain ⟨lower, upper, lower_nonneg, upper_nonneg, factorization⟩ :=
    finTwoPencil_nonnegative_factorization low high low_det_nonneg mixed_nonneg
      high_det_pos lorentzian
  rw [factorization weight]
  positivity

end Matrix
