import GECPKernelStructure.Matrix.SignLockedHomotopy
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Rat.Cast.Order
import Mathlib.Tactic

namespace GECPKernelStructure.Fermionic

open Set

/-- Four equally spaced row nodes, with the first node indexed by `c`. -/
noncomputable def phaseCURows (c : ℕ) : Fin 4 → ℝ := fun i => (c + i.1) / 9

/-- Two ordered unit-interval column nodes chosen so that their exponential
increments on the Phase CU row grid are rational. -/
noncomputable def phaseCUColumns : Fin 2 → ℝ := fun j =>
  if j.1 = 0 then 9 * Real.log (12 / 11) else 9 * Real.log (11 / 10)

/-- The two-exponential augmented matrix with appended monomials `1` and
`x^3`. -/
noncomputable def expCubicAugmented
    (rows : Fin 4 → ℝ) (columns : Fin 2 → ℝ) :
    Matrix (Fin 4) (Fin 4) ℝ := fun i j =>
  if j.1 = 0 then Real.exp (rows i * columns 0)
  else if j.1 = 1 then Real.exp (rows i * columns 1)
  else if j.1 = 2 then 1
  else (rows i) ^ 3

/-- The exact rational matrix remaining after the two positive exponential
column factors are removed. -/
def phaseCUReduced (c : ℕ) : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  if j.1 = 0 then (12 / 11 : ℚ) ^ i.1
  else if j.1 = 1 then (11 / 10 : ℚ) ^ i.1
  else if j.1 = 2 then 1
  else ((c + i.1 : ℚ) / 9) ^ 3

noncomputable def phaseCUWeights (c : ℕ) : Fin 4 → ℝ := fun j =>
  if j.1 = 0 then (12 / 11 : ℝ) ^ c
  else if j.1 = 1 then (11 / 10 : ℝ) ^ c
  else 1

noncomputable def phaseCUReducedReal (c : ℕ) : Matrix (Fin 4) (Fin 4) ℝ :=
  (Rat.castHom ℝ).mapMatrix (phaseCUReduced c)

theorem phaseCURows_strictMono (c : ℕ) : StrictMono (phaseCURows c) := by
  intro i j hij
  simp only [phaseCURows]
  apply (div_lt_div_iff_of_pos_right (by norm_num : (0 : ℝ) < 9)).mpr
  exact_mod_cast Nat.add_lt_add_left hij c

theorem phaseCURows_one_mem_openUnit (i : Fin 4) :
    phaseCURows 1 i ∈ Ioo (0 : ℝ) 1 := by
  fin_cases i <;> norm_num [phaseCURows]

theorem phaseCURows_five_mem_openUnit (i : Fin 4) :
    phaseCURows 5 i ∈ Ioo (0 : ℝ) 1 := by
  fin_cases i <;> norm_num [phaseCURows]

theorem phaseCUColumns_strictMono : StrictMono phaseCUColumns := by
  intro i j hij
  have log_order : Real.log (12 / 11) < Real.log (11 / 10) :=
    Real.strictMonoOn_log (by norm_num) (by norm_num) (by norm_num)
  change i.1 < j.1 at hij
  fin_cases i <;> fin_cases j
  all_goals norm_num at hij
  simpa [phaseCUColumns] using mul_lt_mul_of_pos_left log_order (by norm_num :
    (0 : ℝ) < 9)

theorem phaseCUColumns_mem_openUnit (j : Fin 2) :
    phaseCUColumns j ∈ Ioo (0 : ℝ) 1 := by
  have log_twelve_pos : 0 < Real.log (12 / 11) :=
    Real.log_pos (by norm_num)
  have log_eleven_pos : 0 < Real.log (11 / 10) :=
    Real.log_pos (by norm_num)
  have log_eleven_lt : Real.log (11 / 10) < (11 / 10 : ℝ) - 1 :=
    Real.log_lt_sub_one_of_pos (by norm_num) (by norm_num)
  have log_order : Real.log (12 / 11) < Real.log (11 / 10) :=
    Real.strictMonoOn_log (by norm_num) (by norm_num) (by norm_num)
  fin_cases j
  · constructor
    · simpa [phaseCUColumns] using mul_pos (by norm_num : (0 : ℝ) < 9)
        log_twelve_pos
    · simp only [phaseCUColumns]
      norm_num only [if_pos]
      nlinarith
  · constructor
    · simpa [phaseCUColumns] using mul_pos (by norm_num : (0 : ℝ) < 9)
        log_eleven_pos
    · simp only [phaseCUColumns]
      norm_num only [if_neg]
      nlinarith

/-- On the selected rational exponential increments, the augmented matrix is
a positive column scaling of the exact reduced matrix. -/
theorem expCubicAugmented_phaseCURows_eq (c : ℕ) :
    expCubicAugmented (phaseCURows c) phaseCUColumns =
      fun i j => phaseCUWeights c j * phaseCUReducedReal c i j := by
  ext i j
  fin_cases j
  · simp only [expCubicAugmented, phaseCURows, phaseCUWeights,
      phaseCUReducedReal, RingHom.mapMatrix_apply]
    norm_num [phaseCUColumns]
    rw [show ((c : ℝ) + i.1) / 9 * (9 * Real.log (12 / 11)) =
          (c + i.1 : ℕ) * Real.log (12 / 11) by push_cast; ring,
      Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ) < 12 / 11),
      pow_add]
    norm_num
    simp [phaseCUReduced]
  · simp only [expCubicAugmented, phaseCURows, phaseCUWeights,
      phaseCUReducedReal, RingHom.mapMatrix_apply]
    norm_num [phaseCUColumns]
    rw [show ((c : ℝ) + i.1) / 9 * (9 * Real.log (11 / 10)) =
          (c + i.1 : ℕ) * Real.log (11 / 10) by push_cast; ring,
      Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ) < 11 / 10),
      pow_add]
    norm_num
    simp [phaseCUReduced]
  · simp [expCubicAugmented, phaseCUWeights, phaseCUReducedReal,
      phaseCUReduced]
  · simp [expCubicAugmented, phaseCURows, phaseCUWeights, phaseCUReducedReal,
      phaseCUReduced]

private theorem det_fin_four_expansion {R : Type*} [CommRing R]
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

theorem phaseCUReduced_one_det :
    (phaseCUReduced 1).det = (83 / 194059800 : ℚ) := by
  rw [det_fin_four_expansion]
  norm_num [phaseCUReduced]

theorem phaseCUReduced_five_det :
    (phaseCUReduced 5).det = (-1 / 194059800 : ℚ) := by
  rw [det_fin_four_expansion]
  norm_num [phaseCUReduced]

theorem expCubicAugmented_phaseCURows_det (c : ℕ) :
    (expCubicAugmented (phaseCURows c) phaseCUColumns).det =
      (∏ j, phaseCUWeights c j) * (((phaseCUReduced c).det : ℚ) : ℝ) := by
  rw [expCubicAugmented_phaseCURows_eq]
  calc
    Matrix.det ((fun i j => phaseCUWeights c j * phaseCUReducedReal c i j) :
        Matrix (Fin 4) (Fin 4) ℝ) =
        (∏ j, phaseCUWeights c j) * (phaseCUReducedReal c).det := by
      exact Matrix.det_mul_row (phaseCUWeights c) (phaseCUReducedReal c)
    _ = (∏ j, phaseCUWeights c j) * (((phaseCUReduced c).det : ℚ) : ℝ) := by
      congr 1
      simpa [phaseCUReducedReal] using (Rat.cast_det (phaseCUReduced c)).symm

theorem expCubicAugmented_phaseCURows_one_det_pos :
    0 < (expCubicAugmented (phaseCURows 1) phaseCUColumns).det := by
  rw [expCubicAugmented_phaseCURows_det, phaseCUReduced_one_det]
  apply mul_pos
  · exact Finset.prod_pos fun j _ => by
      fin_cases j <;> norm_num [phaseCUWeights]
  · norm_num

theorem expCubicAugmented_phaseCURows_five_det_neg :
    (expCubicAugmented (phaseCURows 5) phaseCUColumns).det < 0 := by
  rw [expCubicAugmented_phaseCURows_det, phaseCUReduced_five_det]
  apply mul_neg_of_pos_of_neg
  · exact Finset.prod_pos fun j _ => by
      fin_cases j <;> norm_num [phaseCUWeights]
  · norm_num

/-- Straight interpolation between the two row tuples in the obstruction. -/
noncomputable def phaseCURowPath (t : ℝ) : Fin 4 → ℝ :=
  fun i => (1 + 4 * t + i.1) / 9

noncomputable def phaseCUMatrixPath (t : ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  expCubicAugmented (phaseCURowPath t) phaseCUColumns

theorem phaseCURowPath_strictMono (t : ℝ) : StrictMono (phaseCURowPath t) := by
  intro i j hij
  simp only [phaseCURowPath]
  apply (div_lt_div_iff_of_pos_right (by norm_num : (0 : ℝ) < 9)).mpr
  have hij_real : (i.1 : ℝ) < j.1 := by exact_mod_cast hij
  linarith

theorem phaseCURowPath_mem_openUnit (t : ℝ) (t_mem : t ∈ Icc (0 : ℝ) 1)
    (i : Fin 4) : phaseCURowPath t i ∈ Ioo (0 : ℝ) 1 := by
  fin_cases i <;> constructor <;> simp only [phaseCURowPath] <;> norm_num at t_mem ⊢ <;>
    nlinarith [t_mem.1, t_mem.2]

theorem phaseCUMatrixPath_entry_continuous (i j : Fin 4) :
    Continuous fun t => phaseCUMatrixPath t i j := by
  fin_cases j <;> simp [phaseCUMatrixPath, expCubicAugmented, phaseCURowPath] <;>
    fun_prop

/-- Universal ordered-unit-chamber nonvanishing is false: the straight path
between the two exact ordered configurations has a singular member. -/
theorem phaseCUMatrixPath_exists_det_eq_zero :
    ∃ t ∈ Icc (0 : ℝ) 1, (phaseCUMatrixPath t).det = 0 := by
  apply Matrix.exists_eq_zero_of_continuousOn_of_mul_neg
  · apply Continuous.continuousOn
    apply Continuous.matrix_det
    apply continuous_pi
    intro i
    apply continuous_pi
    intro j
    exact phaseCUMatrixPath_entry_continuous i j
  · have first := expCubicAugmented_phaseCURows_one_det_pos
    have last := expCubicAugmented_phaseCURows_five_det_neg
    have path_zero : phaseCUMatrixPath 0 =
        expCubicAugmented (phaseCURows 1) phaseCUColumns := by
      rw [phaseCUMatrixPath]
      congr 1
      ext i
      norm_num [phaseCURowPath, phaseCURows]
    have path_one : phaseCUMatrixPath 1 =
        expCubicAugmented (phaseCURows 5) phaseCUColumns := by
      rw [phaseCUMatrixPath]
      congr 1
      ext i
      norm_num [phaseCURowPath, phaseCURows]
    rw [path_zero, path_one]
    exact mul_neg_of_pos_of_neg first last

end GECPKernelStructure.Fermionic
