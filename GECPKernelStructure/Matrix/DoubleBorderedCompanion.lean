import GECPKernelStructure.Matrix.BorderedPencil

namespace Matrix

variable {core fresh : Type*}
  [Fintype core] [DecidableEq core]
  [Fintype fresh] [DecidableEq fresh]

/-- The closed two-sided companion projection associated with an invertible
core and separate kernel/companion cross blocks. -/
noncomputable def twoSidedCompanionSchur
    (C : Matrix core core ℝ) [Invertible C]
    (KLeft : Matrix fresh core ℝ) (KRight : Matrix core fresh ℝ)
    (HCore : Matrix core core ℝ)
    (HLeft : Matrix fresh core ℝ) (HRight : Matrix core fresh ℝ)
    (H : Matrix fresh fresh ℝ) : Matrix fresh fresh ℝ :=
  H - KLeft * ⅟C * HRight - HLeft * ⅟C * KRight +
    KLeft * ⅟C * HCore * ⅟C * KRight

/-- The double-bordered matrix whose Schur complement is
`twoSidedCompanionSchur`.  Its leading block is upper triangular with two
copies of the selected core. -/
def doubleBorderedCompanionMatrix
    (C : Matrix core core ℝ)
    (KLeft : Matrix fresh core ℝ) (KRight : Matrix core fresh ℝ)
    (HCore : Matrix core core ℝ)
    (HLeft : Matrix fresh core ℝ) (HRight : Matrix core fresh ℝ)
    (H : Matrix fresh fresh ℝ) :
    Matrix ((core ⊕ core) ⊕ fresh) ((core ⊕ core) ⊕ fresh) ℝ :=
  Matrix.fromBlocks
    (Matrix.fromBlocks C HCore 0 C)
    (fun index j => Sum.elim (fun i => HRight i j) (fun i => KRight i j) index)
    (fun i index => Sum.elim (fun j => KLeft i j) (fun j => HLeft i j) index)
    H

set_option linter.style.haveILetI false in
/-- The complete double-bordered determinant is the squared core determinant
times the determinant of the two-sided companion Schur expression. -/
theorem doubleBorderedCompanionMatrix_det
    (C : Matrix core core ℝ) [Invertible C]
    (KLeft : Matrix fresh core ℝ) (KRight : Matrix core fresh ℝ)
    (HCore : Matrix core core ℝ)
    (HLeft : Matrix fresh core ℝ) (HRight : Matrix core fresh ℝ)
    (H : Matrix fresh fresh ℝ) :
    (doubleBorderedCompanionMatrix C KLeft KRight HCore HLeft HRight H).det =
      C.det ^ 2 *
        (twoSidedCompanionSchur C KLeft KRight HCore HLeft HRight H).det := by
  classical
  let Q : Matrix (core ⊕ core) (core ⊕ core) ℝ :=
    Matrix.fromBlocks C HCore 0 C
  let T : Matrix (core ⊕ core) fresh ℝ :=
    fun index j => Sum.elim (fun i => HRight i j) (fun i => KRight i j) index
  let L : Matrix fresh (core ⊕ core) ℝ :=
    fun i index => Sum.elim (fun j => KLeft i j) (fun j => HLeft i j) index
  letI : Invertible Q := Matrix.fromBlocksZero₂₁Invertible C HCore C
  have hQinv : ⅟Q = Matrix.fromBlocks (⅟C) (-(⅟C * HCore * ⅟C)) 0 (⅟C) := by
    exact Matrix.invOf_fromBlocks_zero₂₁_eq C HCore C
  let P : Matrix fresh (core ⊕ core) ℝ := fun i index =>
    Sum.elim
      (fun j => (KLeft * ⅟C) i j)
      (fun j =>
        (KLeft * (-(⅟C * HCore * ⅟C)) + HLeft * ⅟C) i j)
      index
  have hLQ : L * Matrix.fromBlocks (⅟C) (-(⅟C * HCore * ⅟C)) 0 (⅟C) = P := by
    ext i (j | j)
    · rw [Matrix.mul_apply, Fintype.sum_sum_type]
      simp [L, P, Matrix.mul_apply]
    · rw [Matrix.mul_apply, Fintype.sum_sum_type]
      simp [L, P, Matrix.mul_apply]
  let Sprod : Matrix fresh fresh ℝ :=
    (KLeft * ⅟C) * HRight +
      (KLeft * (-(⅟C * HCore * ⅟C)) + HLeft * ⅟C) * KRight
  have hPT : P * T = Sprod := by
    ext i j
    rw [Matrix.mul_apply, Fintype.sum_sum_type]
    simp [P, T, Sprod, Matrix.mul_apply]
  have hschur : H - L * ⅟Q * T =
      twoSidedCompanionSchur C KLeft KRight HCore HLeft HRight H := by
    rw [hQinv, hLQ, hPT]
    simp only [Sprod, twoSidedCompanionSchur, Matrix.neg_mul,
      Matrix.mul_neg, Matrix.add_mul, Matrix.mul_assoc]
    abel
  rw [show doubleBorderedCompanionMatrix C KLeft KRight HCore HLeft HRight H =
      Matrix.fromBlocks Q T L H by rfl]
  rw [Matrix.det_fromBlocks₁₁, hschur]
  rw [show Q.det = C.det * C.det by
    exact Matrix.det_fromBlocks_zero₂₁ C HCore C]
  ring

/-- Absolute-value form of the squared-core factorization, ready for a
companion determinant envelope. -/
theorem abs_doubleBorderedCompanionMatrix_det
    (C : Matrix core core ℝ) [Invertible C]
    (KLeft : Matrix fresh core ℝ) (KRight : Matrix core fresh ℝ)
    (HCore : Matrix core core ℝ)
    (HLeft : Matrix fresh core ℝ) (HRight : Matrix core fresh ℝ)
    (H : Matrix fresh fresh ℝ) :
    abs (doubleBorderedCompanionMatrix C KLeft KRight
        HCore HLeft HRight H).det =
      abs C.det ^ 2 *
        abs (twoSidedCompanionSchur C KLeft KRight
          HCore HLeft HRight H).det := by
  rw [doubleBorderedCompanionMatrix_det, abs_mul, abs_pow]

end Matrix
