import GECPKernelStructure.Matrix.CauchyBinet
import Mathlib.Tactic

namespace Matrix

/-- The first three positive Taylor features at the two physical-domain row
nodes `1/8` and `1/4`. -/
noncomputable def phaseCHLeft : Matrix (Fin 2) (Fin 3) ℝ :=
  !![1, 1 / 8, 1 / 128;
     1, 1 / 4, 1 / 32]

/-- The matching monomial features at the two positive column nodes `6` and
`22`. -/
noncomputable def phaseCHRight : Matrix (Fin 3) (Fin 2) ℝ :=
  !![1, 1;
     6, 22;
     36, 484]

/-- The finite three-feature Taylor kernel of the Phase CH witness. -/
noncomputable def phaseCHKernel : Matrix (Fin 2) (Fin 2) ℝ :=
  phaseCHLeft * phaseCHRight

theorem phaseCHKernel_eq :
    phaseCHKernel =
      !![65 / 32, 241 / 32;
         29 / 8, 173 / 8] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCHKernel, phaseCHLeft, phaseCHRight, Matrix.mul_apply,
      Fin.sum_univ_succ]

/-- The left one-pivot bordered determinant for one Taylor feature. -/
noncomputable def phaseCHLeftBordered (degree : Fin 3) :
    Matrix (Fin 2) (Fin 2) ℝ :=
  !![phaseCHKernel 0 0, phaseCHLeft 0 degree;
     phaseCHKernel 1 0, phaseCHLeft 1 degree]

/-- The right one-pivot bordered determinant for the same Taylor feature. -/
noncomputable def phaseCHRightBordered (degree : Fin 3) :
    Matrix (Fin 2) (Fin 2) ℝ :=
  !![phaseCHKernel 0 0, phaseCHKernel 0 1;
     phaseCHRight degree 0, phaseCHRight degree 1]

/-- The matched bordered Plücker product at one feature coordinate. -/
noncomputable def phaseCHBorderedProduct (degree : Fin 3) : ℝ :=
  (phaseCHLeftBordered degree).det * (phaseCHRightBordered degree).det

theorem phaseCHBorderedProduct_zero :
    phaseCHBorderedProduct 0 = 561 / 64 := by
  rw [phaseCHBorderedProduct, phaseCHLeftBordered, phaseCHRightBordered,
    phaseCHKernel_eq, Matrix.det_fin_two]
  norm_num [phaseCHLeft, phaseCHRight, Matrix.det_fin_two]

theorem phaseCHBorderedProduct_one :
    phaseCHBorderedProduct 1 = -(7 / 256) := by
  rw [phaseCHBorderedProduct, phaseCHLeftBordered, phaseCHRightBordered,
    phaseCHKernel_eq]
  norm_num [phaseCHLeft, phaseCHRight, Matrix.det_fin_two]

theorem phaseCHBorderedProduct_two :
    phaseCHBorderedProduct 2 = 801 / 32 := by
  rw [phaseCHBorderedProduct, phaseCHLeftBordered, phaseCHRightBordered,
    phaseCHKernel_eq, Matrix.det_fin_two]
  norm_num [phaseCHLeft, phaseCHRight, Matrix.det_fin_two, cons_val_two,
    vecHead, vecTail]

/-- The middle individual coordinate is genuinely negative, so a proof by
coordinatewise Plücker positivity cannot cover even this minimal witness. -/
theorem phaseCH_middle_borderedProduct_neg :
    phaseCHBorderedProduct 1 < 0 := by
  rw [phaseCHBorderedProduct_one]
  norm_num

/-- At the first nontrivial cutoff, the two-coordinate suffix cancels the
negative middle coordinate and is strictly positive. -/
theorem phaseCH_first_suffix_pos :
    0 < phaseCHBorderedProduct 1 + phaseCHBorderedProduct 2 := by
  rw [phaseCHBorderedProduct_one, phaseCHBorderedProduct_two]
  norm_num

/-- At the second nontrivial cutoff, the two-coordinate prefix cancels the
negative middle coordinate and is strictly positive. -/
theorem phaseCH_second_prefix_pos :
    0 < phaseCHBorderedProduct 0 + phaseCHBorderedProduct 1 := by
  rw [phaseCHBorderedProduct_zero, phaseCHBorderedProduct_one]
  norm_num

/-- All proper nonempty prefix and suffix groups of the three coordinate
products are positive. -/
theorem phaseCH_proper_contiguous_groups_pos :
    0 < phaseCHBorderedProduct 0 ∧
      0 < phaseCHBorderedProduct 0 + phaseCHBorderedProduct 1 ∧
      0 < phaseCHBorderedProduct 1 + phaseCHBorderedProduct 2 ∧
      0 < phaseCHBorderedProduct 2 := by
  rw [phaseCHBorderedProduct_zero, phaseCHBorderedProduct_one,
    phaseCHBorderedProduct_two]
  norm_num

/-- The scalar two-level step polynomial at the first nontrivial cutoff. -/
noncomputable def phaseCHFirstStep (weight : ℝ) : ℝ :=
  phaseCHBorderedProduct 0 +
    weight * (phaseCHBorderedProduct 1 + phaseCHBorderedProduct 2)

/-- The scalar two-level step polynomial at the second nontrivial cutoff. -/
noncomputable def phaseCHSecondStep (weight : ℝ) : ℝ :=
  phaseCHBorderedProduct 0 + phaseCHBorderedProduct 1 +
    weight * phaseCHBorderedProduct 2

/-- Grouping makes the first-cut step determinant strictly positive for every
nonnegative amplification weight, despite its negative middle coordinate. -/
theorem phaseCH_firstStep_pos (weight : ℝ) (weight_nonneg : 0 ≤ weight) :
    0 < phaseCHFirstStep weight := by
  rw [phaseCHFirstStep, phaseCHBorderedProduct_zero,
    phaseCHBorderedProduct_one, phaseCHBorderedProduct_two]
  positivity

/-- Grouping likewise makes the second-cut step determinant strictly positive
for every nonnegative amplification weight. -/
theorem phaseCH_secondStep_pos (weight : ℝ) (weight_nonneg : 0 ≤ weight) :
    0 < phaseCHSecondStep weight := by
  rw [phaseCHSecondStep, phaseCHBorderedProduct_zero,
    phaseCHBorderedProduct_one, phaseCHBorderedProduct_two]
  positivity

end Matrix
