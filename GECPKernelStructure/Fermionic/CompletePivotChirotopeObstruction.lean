import GECPKernelStructure.Fermionic.CompleteRun
import GECPKernelStructure.Fermionic.OrderedChamberObstruction
import GECPKernelStructure.Matrix.FinFourDeterminant

namespace GECPKernelStructure.Fermionic

open GECP Set

def phaseCVBaseQ : ℚ := 101 / 100

/-- The exact rational exponential kernel used by the Phase CV obstruction. -/
def phaseCVKernelQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  phaseCVBaseQ ^ ((i.1 + 1) * (j.1 + 21))

noncomputable def phaseCVKernel : Kernel (Fin 4) (Fin 4) ℝ := fun i j =>
  (phaseCVKernelQ i j : ℚ)

/-- Open-unit row nodes whose products with `phaseCVColumns` turn the
exponential kernel into `phaseCVKernelQ`. -/
noncomputable def phaseCVRows : Fin 4 → ℝ := fun i =>
  (6 / 25) * (i.1 + 1)

noncomputable def phaseCVColumns : Fin 4 → ℝ := fun j =>
  (25 / 6) * Real.log (101 / 100) * (j.1 + 21)

theorem phaseCVRows_strictMono : StrictMono phaseCVRows := by
  intro i j hij
  simp only [phaseCVRows]
  have hij_real : (i.1 : ℝ) < j.1 := by exact_mod_cast hij
  nlinarith

theorem phaseCVRows_mem_openUnit (i : Fin 4) :
    phaseCVRows i ∈ Ioo (0 : ℝ) 1 := by
  fin_cases i <;> norm_num [phaseCVRows]

theorem phaseCVColumns_strictMono : StrictMono phaseCVColumns := by
  intro i j hij
  simp only [phaseCVColumns]
  have log_pos : 0 < Real.log (101 / 100) := Real.log_pos (by norm_num)
  have hij_real : (i.1 : ℝ) < j.1 := by exact_mod_cast hij
  nlinarith

theorem phaseCVColumns_mem_openUnit (j : Fin 4) :
    phaseCVColumns j ∈ Ioo (0 : ℝ) 1 := by
  have log_pos : 0 < Real.log (101 / 100) := Real.log_pos (by norm_num)
  have log_lt : Real.log (101 / 100) < (101 / 100 : ℝ) - 1 :=
    Real.log_lt_sub_one_of_pos (by norm_num) (by norm_num)
  fin_cases j <;> constructor <;> simp only [phaseCVColumns] <;> norm_num <;>
    nlinarith

theorem phaseCV_expKernel_eq (i j : Fin 4) :
    Real.exp (phaseCVRows i * phaseCVColumns j) = phaseCVKernel i j := by
  rw [show phaseCVRows i * phaseCVColumns j =
      ((i.1 + 1) * (j.1 + 21) : ℕ) * Real.log (101 / 100) by
    simp only [phaseCVRows, phaseCVColumns]
    push_cast
    ring]
  rw [Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ) < 101 / 100)]
  norm_num [phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]

theorem phaseCVKernel_pivot_ne : phaseCVKernel 3 3 ≠ 0 := by
  norm_num [phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]

theorem phaseCV_firstPivot_complete :
    CompletePivotOn (fun _ : Fin 4 => True) (fun _ : Fin 4 => True)
      phaseCVKernel 3 3 := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ, abs_of_nonneg]

noncomputable def phaseCVFirstResidual : Kernel (Fin 4) (Fin 4) ℝ :=
  residualUpdate phaseCVKernel 3 3 phaseCVKernel_pivot_ne

theorem phaseCVFirstResidual_pivot_ne : phaseCVFirstResidual 0 0 ≠ 0 := by
  norm_num [phaseCVFirstResidual, residualUpdate, phaseCVKernel,
    phaseCVKernelQ, phaseCVBaseQ]

theorem phaseCV_secondPivot_complete :
    CompletePivotOn (fun _ : Fin 4 => True) (fun _ : Fin 4 => True)
      phaseCVFirstResidual 0 0 := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCVFirstResidual, residualUpdate, phaseCVKernel,
      phaseCVKernelQ, phaseCVBaseQ, abs_of_nonneg, abs_of_nonpos]

noncomputable def phaseCVSelectedRows : Fin 2 → ℝ := fun j =>
  if j.1 = 0 then phaseCVRows 0 else phaseCVRows 3

noncomputable def phaseCVSelectedColumns : Fin 2 → ℝ := fun j =>
  if j.1 = 0 then phaseCVColumns 0 else phaseCVColumns 3

def phaseCVLeftReducedQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  if j.1 = 0 then phaseCVBaseQ ^ ((i.1 + 1) * 21)
  else if j.1 = 1 then phaseCVBaseQ ^ ((i.1 + 1) * 24)
  else if j.1 = 2 then 1
  else (i.1 + 1 : ℚ) ^ 3

def phaseCVRightReducedQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  if j.1 = 0 then phaseCVBaseQ ^ (i.1 + 21)
  else if j.1 = 1 then phaseCVBaseQ ^ (4 * (i.1 + 21))
  else if j.1 = 2 then 1
  else (i.1 + 21 : ℚ) ^ 3

noncomputable def phaseCVLeftWeights : Fin 4 → ℝ := fun j =>
  if j.1 = 3 then (6 / 25 : ℝ) ^ 3 else 1

noncomputable def phaseCVRightWeights : Fin 4 → ℝ := fun j =>
  if j.1 = 3 then ((25 / 6) * Real.log (101 / 100)) ^ 3 else 1

noncomputable def phaseCVLeftReduced : Matrix (Fin 4) (Fin 4) ℝ :=
  (Rat.castHom ℝ).mapMatrix phaseCVLeftReducedQ

noncomputable def phaseCVRightReduced : Matrix (Fin 4) (Fin 4) ℝ :=
  (Rat.castHom ℝ).mapMatrix phaseCVRightReducedQ

theorem phaseCVLeftAugmented_eq :
    expCubicAugmented phaseCVRows phaseCVSelectedColumns =
      fun i j => phaseCVLeftWeights j * phaseCVLeftReduced i j := by
  ext i j
  fin_cases j
  · simp [expCubicAugmented, phaseCVSelectedColumns, phaseCVLeftWeights,
      phaseCVLeftReduced, phaseCVLeftReducedQ, phaseCV_expKernel_eq,
      phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]
  · simp [expCubicAugmented, phaseCVSelectedColumns, phaseCVLeftWeights,
      phaseCVLeftReduced, phaseCVLeftReducedQ, phaseCV_expKernel_eq,
      phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]
  · simp [expCubicAugmented, phaseCVLeftWeights, phaseCVLeftReduced,
      phaseCVLeftReducedQ]
  · simp [expCubicAugmented, phaseCVRows, phaseCVLeftWeights,
      phaseCVLeftReduced, phaseCVLeftReducedQ]
    ring

theorem phaseCVRightAugmented_eq :
    expCubicAugmented phaseCVColumns phaseCVSelectedRows =
      fun i j => phaseCVRightWeights j * phaseCVRightReduced i j := by
  ext i j
  fin_cases j
  · simp [expCubicAugmented, phaseCVSelectedRows, phaseCVRightWeights,
      phaseCVRightReduced, phaseCVRightReducedQ, phaseCVBaseQ]
    rw [mul_comm, phaseCV_expKernel_eq 0 i]
    norm_num [phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]
  · simp [expCubicAugmented, phaseCVSelectedRows, phaseCVRightWeights,
      phaseCVRightReduced, phaseCVRightReducedQ, phaseCVBaseQ]
    rw [mul_comm, phaseCV_expKernel_eq 3 i]
    norm_num [phaseCVKernel, phaseCVKernelQ, phaseCVBaseQ]
  · simp [expCubicAugmented, phaseCVRightWeights, phaseCVRightReduced,
      phaseCVRightReducedQ]
  · simp [expCubicAugmented, phaseCVColumns, phaseCVRightWeights,
      phaseCVRightReduced, phaseCVRightReducedQ]
    ring

theorem phaseCVLeftReduced_det_pos : 0 < phaseCVLeftReducedQ.det := by
  rw [Matrix.det_fin_four_expansion]
  norm_num [phaseCVLeftReducedQ, phaseCVBaseQ]

theorem phaseCVRightReduced_det_neg : phaseCVRightReducedQ.det < 0 := by
  rw [Matrix.det_fin_four_expansion]
  norm_num [phaseCVRightReducedQ, phaseCVBaseQ]

theorem phaseCVLeftAugmented_det_pos :
    0 < (expCubicAugmented phaseCVRows phaseCVSelectedColumns).det := by
  rw [phaseCVLeftAugmented_eq]
  have determinant_eq :
      Matrix.det ((fun i j => phaseCVLeftWeights j * phaseCVLeftReduced i j) :
        Matrix (Fin 4) (Fin 4) ℝ) =
        (∏ j, phaseCVLeftWeights j) * phaseCVLeftReduced.det :=
    Matrix.det_mul_row phaseCVLeftWeights phaseCVLeftReduced
  rw [determinant_eq]
  apply mul_pos
  · exact Finset.prod_pos fun j _ => by
      fin_cases j <;> norm_num [phaseCVLeftWeights]
  · change 0 < (phaseCVLeftReducedQ.map fun x => (x : ℝ)).det
    rw [← Rat.cast_det]
    exact_mod_cast phaseCVLeftReduced_det_pos

theorem phaseCVRightAugmented_det_neg :
    (expCubicAugmented phaseCVColumns phaseCVSelectedRows).det < 0 := by
  rw [phaseCVRightAugmented_eq]
  have determinant_eq :
      Matrix.det ((fun i j => phaseCVRightWeights j * phaseCVRightReduced i j) :
        Matrix (Fin 4) (Fin 4) ℝ) =
        (∏ j, phaseCVRightWeights j) * phaseCVRightReduced.det :=
    Matrix.det_mul_row phaseCVRightWeights phaseCVRightReduced
  rw [determinant_eq]
  apply mul_neg_of_pos_of_neg
  · exact Finset.prod_pos fun j _ => by
      fin_cases j <;> simp [phaseCVRightWeights]
      positivity
  · change (phaseCVRightReducedQ.map fun x => (x : ℝ)).det < 0
    rw [← Rat.cast_det]
    exact_mod_cast phaseCVRightReduced_det_neg

/-- Even two exact complete pivots do not force the crossed augmented
chirotopes to agree. -/
theorem phaseCV_crossedAugmented_det_mul_neg :
    (expCubicAugmented phaseCVRows phaseCVSelectedColumns).det *
      (expCubicAugmented phaseCVColumns phaseCVSelectedRows).det < 0 :=
  mul_neg_of_pos_of_neg phaseCVLeftAugmented_det_pos
    phaseCVRightAugmented_det_neg

end GECPKernelStructure.Fermionic
