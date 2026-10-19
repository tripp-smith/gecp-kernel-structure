import GECPKernelStructure.Matrix.FinTwoPencil

namespace Matrix

/-- A quadratic determinant pencil evaluated at one scalar weight. -/
def quadraticPencilValue (constant mixed quadratic weight : ℝ) : ℝ :=
  constant + mixed * weight + quadratic * weight ^ 2

/-- Moving one rank-one feature from the upper block to the lower block adds
`left` to the constant coefficient, removes `right` from the quadratic
coefficient, and changes the mixed coefficient by `right - left`. -/
theorem quadraticTransfer_eq
    (constant mixed quadratic left right weight : ℝ) :
    quadraticPencilValue
        (constant + left) (mixed - left + right) (quadratic - right) weight =
      quadraticPencilValue constant mixed quadratic weight +
        (1 - weight) * (left + right * weight) := by
  unfold quadraticPencilValue
  ring

/-- The numerator obtained by evaluating a quadratic pencil at the canonical
transfer node `-left / right`, with the denominator `right^2` cleared. -/
def quadraticTransferNodeNumerator
    (constant mixed quadratic left right : ℝ) : ℝ :=
  constant * right ^ 2 - mixed * left * right + quadratic * left ^ 2

theorem right_sq_mul_quadraticPencilValue_neg_div
    (constant mixed quadratic left right : ℝ) (right_ne : right ≠ 0) :
    right ^ 2 * quadraticPencilValue constant mixed quadratic (-left / right) =
      quadraticTransferNodeNumerator constant mixed quadratic left right := by
  unfold quadraticPencilValue quadraticTransferNodeNumerator
  field_simp [right_ne]
  ring

/-- The canonical transfer-node numerator is shared by both adjacent cutoff
pencils.  This is the division-free form of the common-value observation at
the second root of their difference. -/
theorem quadraticTransferNodeNumerator_invariant
    (constant mixed quadratic left right : ℝ) :
    quadraticTransferNodeNumerator
        (constant + left) (mixed - left + right) (quadratic - right)
        left right =
      quadraticTransferNodeNumerator constant mixed quadratic left right := by
  unfold quadraticTransferNodeNumerator
  ring

/-- A nonpositive value at the canonical negative transfer node forces the
quadratic determinant pencil to have nonnegative discriminant.  The statement
uses only a cleared numerator and therefore introduces no division. -/
theorem quadratic_discriminant_nonneg_of_transferNode
    (constant mixed quadratic left right : ℝ)
    (constant_nonneg : 0 ≤ constant)
    (quadratic_nonneg : 0 ≤ quadratic)
    (left_pos : 0 < left)
    (right_pos : 0 < right)
    (node_nonpos :
      quadraticTransferNodeNumerator
        constant mixed quadratic left right ≤ 0) :
    4 * constant * quadratic ≤ mixed ^ 2 := by
  have sum_nonneg :
      0 ≤ constant * right ^ 2 + quadratic * left ^ 2 :=
    add_nonneg
      (mul_nonneg constant_nonneg (sq_nonneg right))
      (mul_nonneg quadratic_nonneg (sq_nonneg left))
  have difference_nonneg :
      0 ≤ mixed * left * right -
        (constant * right ^ 2 + quadratic * left ^ 2) := by
    unfold quadraticTransferNodeNumerator at node_nonpos
    linarith
  have mixed_transfer_nonneg : 0 ≤ mixed * left * right := by
    linarith
  have sum_with_mixed_nonneg :
      0 ≤ mixed * left * right +
        (constant * right ^ 2 + quadratic * left ^ 2) :=
    add_nonneg mixed_transfer_nonneg sum_nonneg
  have scaled_nonneg :
      0 ≤ (mixed ^ 2 - 4 * constant * quadratic) * (left * right) ^ 2 := by
    calc
      0 ≤
          (mixed * left * right -
              (constant * right ^ 2 + quadratic * left ^ 2)) *
              (mixed * left * right +
                (constant * right ^ 2 + quadratic * left ^ 2)) +
            (constant * right ^ 2 - quadratic * left ^ 2) ^ 2 :=
        add_nonneg
          (mul_nonneg difference_nonneg sum_with_mixed_nonneg)
          (sq_nonneg _)
      _ = (mixed ^ 2 - 4 * constant * quadratic) *
          (left * right) ^ 2 := by ring
  have transfer_sq_pos : 0 < (left * right) ^ 2 :=
    sq_pos_of_pos (mul_pos left_pos right_pos)
  have discriminant_nonneg : 0 ≤ mixed ^ 2 - 4 * constant * quadratic :=
    nonneg_of_mul_nonneg_right
      (by simpa [mul_comm] using scaled_nonneg) transfer_sq_pos
  linarith

/-- One canonical transfer-node certificate proves the Lorentzian
discriminant inequalities for both adjacent cutoff pencils. -/
theorem quadraticTransfer_discriminants_nonneg_of_transferNode
    (constant mixed quadratic left right : ℝ)
    (constant_nonneg : 0 ≤ constant)
    (quadratic_nonneg : 0 ≤ quadratic)
    (transferred_constant_nonneg : 0 ≤ constant + left)
    (transferred_quadratic_nonneg : 0 ≤ quadratic - right)
    (left_pos : 0 < left)
    (right_pos : 0 < right)
    (node_nonpos :
      quadraticTransferNodeNumerator
        constant mixed quadratic left right ≤ 0) :
    4 * constant * quadratic ≤ mixed ^ 2 ∧
      4 * (constant + left) * (quadratic - right) ≤
        (mixed - left + right) ^ 2 := by
  constructor
  · exact quadratic_discriminant_nonneg_of_transferNode
      constant mixed quadratic left right constant_nonneg quadratic_nonneg
      left_pos right_pos node_nonpos
  · apply quadratic_discriminant_nonneg_of_transferNode
      (constant + left) (mixed - left + right) (quadratic - right)
      left right transferred_constant_nonneg transferred_quadratic_nonneg
      left_pos right_pos
    rwa [quadraticTransferNodeNumerator_invariant]

theorem finTwoPencilMixed_add_left
    (first second third : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilMixed (first + second) third =
      finTwoPencilMixed first third + finTwoPencilMixed second third := by
  unfold finTwoPencilMixed
  simp only [Matrix.add_apply]
  ring

theorem finTwoPencilMixed_add_right
    (first second third : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilMixed first (second + third) =
      finTwoPencilMixed first second + finTwoPencilMixed first third := by
  unfold finTwoPencilMixed
  simp only [Matrix.add_apply]
  ring

theorem finTwoPencilMixed_comm
    (first second : Matrix (Fin 2) (Fin 2) ℝ) :
    finTwoPencilMixed first second = finTwoPencilMixed second first := by
  unfold finTwoPencilMixed
  ring

theorem finTwo_det_add
    (first second : Matrix (Fin 2) (Fin 2) ℝ) :
    (first + second).det =
      first.det + finTwoPencilMixed first second + second.det := by
  rw [Matrix.det_fin_two, Matrix.det_fin_two, Matrix.det_fin_two]
  unfold finTwoPencilMixed
  simp only [Matrix.add_apply]
  ring

/-- Matrix form of the adjacent-cutoff bridge.  A rank-one feature matrix is
transferred from `upper` to `lower`; one nonpositive canonical-node numerator
then proves both neighboring determinant pencils Lorentzian. -/
theorem finTwoPencil_transfer_lorentzian_of_transferNode
    (lower moved upper : Matrix (Fin 2) (Fin 2) ℝ)
    (moved_rank_one : moved.det = 0)
    (lower_det_nonneg : 0 ≤ lower.det)
    (old_upper_det_nonneg : 0 ≤ (moved + upper).det)
    (new_lower_det_nonneg : 0 ≤ (lower + moved).det)
    (upper_det_nonneg : 0 ≤ upper.det)
    (left_pos : 0 < finTwoPencilMixed lower moved)
    (right_pos : 0 < finTwoPencilMixed moved upper)
    (node_nonpos :
      quadraticTransferNodeNumerator
        lower.det
        (finTwoPencilMixed lower (moved + upper))
        (moved + upper).det
        (finTwoPencilMixed lower moved)
        (finTwoPencilMixed moved upper) ≤ 0) :
    4 * lower.det * (moved + upper).det ≤
        finTwoPencilMixed lower (moved + upper) ^ 2 ∧
      4 * (lower + moved).det * upper.det ≤
        finTwoPencilMixed (lower + moved) upper ^ 2 := by
  have constant_transfer :
      lower.det + finTwoPencilMixed lower moved = (lower + moved).det := by
    rw [finTwo_det_add, moved_rank_one]
    ring
  have quadratic_transfer :
      (moved + upper).det - finTwoPencilMixed moved upper = upper.det := by
    rw [finTwo_det_add, moved_rank_one]
    ring
  have mixed_transfer :
      finTwoPencilMixed lower (moved + upper) -
          finTwoPencilMixed lower moved + finTwoPencilMixed moved upper =
        finTwoPencilMixed (lower + moved) upper := by
    rw [finTwoPencilMixed_add_right, finTwoPencilMixed_add_left]
    ring
  obtain ⟨old_lorentzian, new_lorentzian⟩ :=
    quadraticTransfer_discriminants_nonneg_of_transferNode
      lower.det
      (finTwoPencilMixed lower (moved + upper))
      (moved + upper).det
      (finTwoPencilMixed lower moved)
      (finTwoPencilMixed moved upper)
      lower_det_nonneg old_upper_det_nonneg
      (by rwa [constant_transfer])
      (by rwa [quadratic_transfer])
      left_pos right_pos node_nonpos
  exact ⟨old_lorentzian, by
    rwa [constant_transfer, quadratic_transfer, mixed_transfer] at new_lorentzian⟩

/-- A rank-one feature with its left and right scalar coordinates exposed. -/
def finTwoOuterFeature (left right : Fin 2 → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  fun i j ↦ left i * right j

/-- Cross determinant of the two block columns projected against the
perpendicular to the moved right feature. -/
def finTwoProjectedColumnCross
    (lower upper : Matrix (Fin 2) (Fin 2) ℝ) (right : Fin 2 → ℝ) : ℝ :=
  (lower 0 1 * right 0 - lower 0 0 * right 1) *
      (upper 1 1 * right 0 - upper 1 0 * right 1) -
    (upper 0 1 * right 0 - upper 0 0 * right 1) *
      (lower 1 1 * right 0 - lower 1 0 * right 1)

/-- The corresponding cross determinant of the two projected block rows. -/
def finTwoProjectedRowCross
    (lower upper : Matrix (Fin 2) (Fin 2) ℝ) (left : Fin 2 → ℝ) : ℝ :=
  (lower 1 0 * left 0 - lower 0 0 * left 1) *
      (upper 1 1 * left 0 - upper 0 1 * left 1) -
    (upper 1 0 * left 0 - upper 0 0 * left 1) *
      (lower 1 1 * left 0 - lower 0 1 * left 1)

/-- The transfer-node numerator factors into two cumulative projected cross
determinants. Thus the required sign is aggregate left/right sign agreement;
individual feature-coordinate signs are not required. -/
theorem quadraticTransferNodeNumerator_outerFeature_eq
    (lower upper : Matrix (Fin 2) (Fin 2) ℝ) (left right : Fin 2 → ℝ) :
    quadraticTransferNodeNumerator
        lower.det
        (finTwoPencilMixed lower (finTwoOuterFeature left right + upper))
        (finTwoOuterFeature left right + upper).det
        (finTwoPencilMixed lower (finTwoOuterFeature left right))
        (finTwoPencilMixed (finTwoOuterFeature left right) upper) =
      -(finTwoProjectedColumnCross lower upper right *
        finTwoProjectedRowCross lower upper left) := by
  simp only [quadraticTransferNodeNumerator, finTwoPencilMixed,
    Matrix.det_fin_two, finTwoOuterFeature, Matrix.add_apply,
    finTwoProjectedColumnCross, finTwoProjectedRowCross]
  ring

/-- Aggregate projected cross-sign agreement is sufficient for both adjacent
rank-two cutoff pencils to satisfy the sharp Lorentzian inequality. -/
theorem finTwoPencil_transfer_lorentzian_of_projectedCross
    (lower upper : Matrix (Fin 2) (Fin 2) ℝ) (left right : Fin 2 → ℝ)
    (lower_det_nonneg : 0 ≤ lower.det)
    (old_upper_det_nonneg : 0 ≤ (finTwoOuterFeature left right + upper).det)
    (new_lower_det_nonneg : 0 ≤ (lower + finTwoOuterFeature left right).det)
    (upper_det_nonneg : 0 ≤ upper.det)
    (left_pos : 0 < finTwoPencilMixed lower (finTwoOuterFeature left right))
    (right_pos : 0 < finTwoPencilMixed (finTwoOuterFeature left right) upper)
    (projected_cross_nonneg : 0 ≤ finTwoProjectedColumnCross lower upper right *
      finTwoProjectedRowCross lower upper left) :
    4 * lower.det * (finTwoOuterFeature left right + upper).det ≤
        finTwoPencilMixed lower (finTwoOuterFeature left right + upper) ^ 2 ∧
      4 * (lower + finTwoOuterFeature left right).det * upper.det ≤
        finTwoPencilMixed (lower + finTwoOuterFeature left right) upper ^ 2 := by
  apply finTwoPencil_transfer_lorentzian_of_transferNode lower
    (finTwoOuterFeature left right) upper
    (by simp only [Matrix.det_fin_two, finTwoOuterFeature]; ring)
    lower_det_nonneg old_upper_det_nonneg new_lower_det_nonneg upper_det_nonneg
    left_pos right_pos
  rw [quadraticTransferNodeNumerator_outerFeature_eq]
  exact neg_nonpos.mpr projected_cross_nonneg

end Matrix
