import GECPKernelStructure.Matrix.TripleCauchyBinet
import Mathlib.Tactic

namespace Matrix

/-- The sampled duplicated left monomial factor in the minimal Phase CD
obstruction. Rows are one selected upper-layer row, its lower-layer copy, and
one fresh upper-layer row. -/
noncomputable def phaseCDLeft : Matrix (Fin 3) (Fin 6) ℝ :=
  !![1, 1 / 3, 1 / 18, 0, 0, 0;
     0, 0, 0, 1, 1 / 3, 1 / 18;
     1, 2 / 3, 2 / 9, 0, 0, 0]

/-- The three-feature tangent middle factor with positive diagonal weights
`2, 3, 4`. -/
def phaseCDMiddle : Matrix (Fin 6) (Fin 6) ℝ :=
  !![1, 0, 0, 2, 0, 0;
     0, 1, 0, 0, 3, 0;
     0, 0, 1, 0, 0, 4;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1]

/-- The sampled duplicated right monomial factor in the same obstruction.
Columns are one selected upper-layer column, its lower-layer copy, and one
fresh lower-layer column. -/
noncomputable def phaseCDRight : Matrix (Fin 6) (Fin 3) ℝ :=
  !![1, 0, 0;
     1 / 4, 0, 0;
     1 / 16, 0, 0;
     0, 1, 1;
     0, 1 / 4, 3 / 4;
     0, 1 / 16, 9 / 16]

def phaseCDNegativeLeftFlagFn : Fin 3 → Fin 6 := ![0, 1, 5]

def phaseCDNegativeRightFlagFn : Fin 3 → Fin 6 := ![0, 4, 5]

theorem phaseCDNegativeLeftFlag_strictMono :
    StrictMono phaseCDNegativeLeftFlagFn := by
  rw [Fin.strictMono_iff_lt_succ]
  decide

theorem phaseCDNegativeRightFlag_strictMono :
    StrictMono phaseCDNegativeRightFlagFn := by
  rw [Fin.strictMono_iff_lt_succ]
  decide

def phaseCDNegativeLeftFlag :
    {q : Fin 3 → Fin 6 // StrictMono q} :=
  ⟨phaseCDNegativeLeftFlagFn, phaseCDNegativeLeftFlag_strictMono⟩

def phaseCDNegativeRightFlag :
    {r : Fin 3 → Fin 6 // StrictMono r} :=
  ⟨phaseCDNegativeRightFlagFn, phaseCDNegativeRightFlag_strictMono⟩

theorem phaseCDNegativeLeftMinor :
    (phaseCDLeft.submatrix id phaseCDNegativeLeftFlag.1).det =
      -(1 / 54 : ℝ) := by
  have sampled :
      phaseCDLeft.submatrix id phaseCDNegativeLeftFlag.1 =
        !![1, 1 / 3, 0; 0, 0, 1 / 18; 1, 2 / 3, 0] := by
    ext i j
    fin_cases i <;> fin_cases j <;> rfl
  rw [sampled]
  rw [Matrix.det_fin_three]
  change
    (1 : ℝ) * 0 * 0 - 1 * (1 / 18) * (2 / 3) -
      (1 / 3) * 0 * 0 + (1 / 3) * (1 / 18) * 1 +
      0 * 0 * (2 / 3) - 0 * 0 * 1 = -(1 / 54)
  norm_num

theorem phaseCDNegativeMiddleMinor :
    (phaseCDMiddle.submatrix phaseCDNegativeLeftFlag.1
      phaseCDNegativeRightFlag.1).det = (3 : ℝ) := by
  have sampled :
      phaseCDMiddle.submatrix phaseCDNegativeLeftFlag.1
          phaseCDNegativeRightFlag.1 =
        !![1, 0, 0; 0, 3, 0; 0, 0, 1] := by
    ext i j
    fin_cases i <;> fin_cases j <;> rfl
  rw [sampled]
  rw [Matrix.det_fin_three]
  change
    (1 : ℝ) * 3 * 1 - 1 * 0 * 0 - 0 * 0 * 1 +
      0 * 0 * 0 + 0 * 0 * 0 - 0 * 3 * 0 = 3
  norm_num

theorem phaseCDNegativeRightMinor :
    (phaseCDRight.submatrix phaseCDNegativeRightFlag.1 id).det =
      (3 / 32 : ℝ) := by
  have sampled :
      phaseCDRight.submatrix phaseCDNegativeRightFlag.1 id =
        !![1, 0, 0; 0, 1 / 4, 3 / 4; 0, 1 / 16, 9 / 16] := by
    ext i j
    fin_cases i <;> fin_cases j <;> rfl
  rw [sampled]
  rw [Matrix.det_fin_three]
  change
    (1 : ℝ) * (1 / 4) * (9 / 16) - 1 * (3 / 4) * (1 / 16) -
      0 * 0 * (9 / 16) + 0 * (3 / 4) * 0 +
      0 * 0 * (1 / 16) - 0 * (1 / 4) * 0 = 3 / 32
  norm_num

/-- A positive-monomial tangent factorization already has a negative ordered
double-flag summand at three features. -/
theorem phaseCD_negative_triple_term :
    signedOrderedTripleCauchyBinetTerm 1 phaseCDLeft phaseCDMiddle
      phaseCDRight phaseCDNegativeLeftFlag phaseCDNegativeRightFlag =
        -(1 / 192 : ℝ) := by
  rw [signedOrderedTripleCauchyBinetTerm, phaseCDNegativeLeftMinor,
    phaseCDNegativeMiddleMinor, phaseCDNegativeRightMinor]
  norm_num

/-- Consequently, no separate left/middle/right orientation functions can
satisfy Phase CA's termwise interface with the positive external signature
on this exact witness. -/
theorem phaseCD_no_factorwise_positive_orientation :
    ¬ ∃
      (leftOrientation : {q : Fin 3 → Fin 6 // StrictMono q} → ℝ)
      (middleOrientation :
        {q : Fin 3 → Fin 6 // StrictMono q} →
          {r : Fin 3 → Fin 6 // StrictMono r} → ℝ)
      (rightOrientation : {r : Fin 3 → Fin 6 // StrictMono r} → ℝ),
      (∀ q r, (1 : ℝ) = leftOrientation q *
        middleOrientation q r * rightOrientation r) ∧
      (∀ q, 0 ≤ leftOrientation q *
        (phaseCDLeft.submatrix id q.1).det) ∧
      (∀ q r, 0 ≤ middleOrientation q r *
        (phaseCDMiddle.submatrix q.1 r.1).det) ∧
      (∀ r, 0 ≤ rightOrientation r *
        (phaseCDRight.submatrix r.1 id).det) := by
  rintro ⟨leftOrientation, middleOrientation, rightOrientation,
    orientation_mul, left_nonneg, middle_nonneg, right_nonneg⟩
  have term_nonneg := signedOrderedTripleCauchyBinetTerm_nonneg_of_factor_orientations
    1 phaseCDLeft phaseCDMiddle phaseCDRight leftOrientation
      middleOrientation rightOrientation phaseCDNegativeLeftFlag
      phaseCDNegativeRightFlag
      (orientation_mul phaseCDNegativeLeftFlag phaseCDNegativeRightFlag)
      (left_nonneg phaseCDNegativeLeftFlag)
      (middle_nonneg phaseCDNegativeLeftFlag phaseCDNegativeRightFlag)
      (right_nonneg phaseCDNegativeRightFlag)
  rw [phaseCD_negative_triple_term] at term_nonneg
  norm_num at term_nonneg

/-- Failure of termwise sign composition is only a proof obstruction: the
complete distinguished tangent determinant in the witness remains positive. -/
theorem phaseCD_complete_tangent_det :
    (phaseCDLeft * phaseCDMiddle * phaseCDRight).det =
      (24929 / 27648 : ℝ) := by
  have product : phaseCDLeft * phaseCDMiddle * phaseCDRight =
      !![313 / 288, 163 / 72, 23 / 8;
         0, 313 / 288, 41 / 32;
         85 / 72, 23 / 9, 4] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [phaseCDLeft, phaseCDMiddle, phaseCDRight, Matrix.mul_apply,
        Fin.sum_univ_succ]
  rw [product]
  rw [Matrix.det_fin_three]
  change
    (313 / 288 : ℝ) * (313 / 288) * 4 -
      (313 / 288) * (41 / 32) * (23 / 9) -
      (163 / 72) * 0 * 4 + (163 / 72) * (41 / 32) * (85 / 72) +
      (23 / 8) * 0 * (23 / 9) - (23 / 8) * (313 / 288) * (85 / 72) =
        24929 / 27648
  norm_num

theorem phaseCD_complete_tangent_det_pos :
    0 < (phaseCDLeft * phaseCDMiddle * phaseCDRight).det := by
  rw [phaseCD_complete_tangent_det]
  norm_num

end Matrix
