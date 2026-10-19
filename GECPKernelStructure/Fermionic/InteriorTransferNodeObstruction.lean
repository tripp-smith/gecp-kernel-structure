import GECPKernelStructure.Fermionic.CompleteRun
import GECPKernelStructure.GECP.CompanionBiConditioning
import GECPKernelStructure.Matrix.QuadraticTransfer

namespace GECPKernelStructure.Fermionic

open GECP Matrix

def phaseCZRowsQ : Fin 4 → ℚ := ![17 / 21, 18 / 21, 19 / 21, 20 / 21]
def phaseCZColumnsQ : Fin 4 → ℚ := ![15 / 20, 16 / 20, 17 / 20, 19 / 20]

def phaseCZTaylorQ (z : ℚ) : ℚ :=
  1 + z + z ^ 2 / 2 + z ^ 3 / 6 + z ^ 4 / 24 + z ^ 5 / 120 +
    z ^ 6 / 720 + z ^ 7 / 5040 + z ^ 8 / 40320

theorem phaseCZTaylorQ_eq_sum (z : ℚ) :
    phaseCZTaylorQ z = ∑ k : Fin 9, z ^ k.1 / (k.1.factorial : ℚ) := by
  norm_num [phaseCZTaylorQ, Fin.sum_univ_succ]
  ring

def phaseCZKernelQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j ↦
  phaseCZTaylorQ (phaseCZRowsQ i * phaseCZColumnsQ j)

noncomputable def phaseCZKernel : Kernel (Fin 4) (Fin 4) ℝ := fun i j ↦
  (phaseCZKernelQ i j : ℚ)

theorem phaseCZ_nodes_openUnit :
    (∀ i, 0 < phaseCZRowsQ i ∧ phaseCZRowsQ i < 1) ∧
      (∀ j, 0 < phaseCZColumnsQ j ∧ phaseCZColumnsQ j < 1) := by
  constructor <;> intro i <;> fin_cases i <;>
    norm_num [phaseCZRowsQ, phaseCZColumnsQ, cons_val_two, cons_val_three,
      vecHead, vecTail]

theorem phaseCZKernelQ_pivot_ne : phaseCZKernelQ 3 3 ≠ 0 := by
  norm_num [phaseCZKernelQ, phaseCZTaylorQ, phaseCZRowsQ, phaseCZColumnsQ,
    cons_val_two, cons_val_three, vecHead, vecTail]

theorem phaseCZKernel_pivot_ne : phaseCZKernel 3 3 ≠ 0 := by
  change ((phaseCZKernelQ 3 3 : ℚ) : ℝ) ≠ 0
  exact_mod_cast phaseCZKernelQ_pivot_ne

theorem phaseCZ_firstPivot_complete :
    CompletePivotOn (fun _ : Fin 4 ↦ True) (fun _ : Fin 4 ↦ True)
      phaseCZKernel 3 3 := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCZKernel, phaseCZKernelQ, phaseCZTaylorQ,
      phaseCZRowsQ, phaseCZColumnsQ, cons_val_two, cons_val_three,
      vecHead, vecTail, abs_of_nonneg]

def phaseCZFirstResidualQ : Matrix (Fin 4) (Fin 4) ℚ :=
  residualUpdate phaseCZKernelQ 3 3 phaseCZKernelQ_pivot_ne

theorem phaseCZFirstResidualQ_pivot_ne : phaseCZFirstResidualQ 0 0 ≠ 0 := by
  norm_num [phaseCZFirstResidualQ, residualUpdate, phaseCZKernelQ,
    phaseCZTaylorQ, phaseCZRowsQ, phaseCZColumnsQ, cons_val_two,
    cons_val_three, vecHead, vecTail]

noncomputable def phaseCZFirstResidual : Kernel (Fin 4) (Fin 4) ℝ :=
  residualUpdate phaseCZKernel 3 3 phaseCZKernel_pivot_ne

theorem phaseCZ_secondPivot_complete :
    CompletePivotOn (fun _ : Fin 4 ↦ True) (fun _ : Fin 4 ↦ True)
      phaseCZFirstResidual 0 0 := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCZFirstResidual, residualUpdate, phaseCZKernel,
      phaseCZKernelQ, phaseCZTaylorQ, phaseCZRowsQ, phaseCZColumnsQ,
      cons_val_two, cons_val_three, vecHead, vecTail,
      abs_of_nonneg, abs_of_nonpos]

def phaseCZLeftQ : Matrix (Fin 4) (Fin 9) ℚ := fun i k ↦
  phaseCZRowsQ i ^ k.1 / (k.1.factorial : ℚ)

def phaseCZRightQ : Matrix (Fin 9) (Fin 4) ℚ := fun k j ↦
  phaseCZColumnsQ j ^ k.1

/-- All nine retained features reconstruct the pivot kernel exactly. -/
theorem phaseCZLeftQ_mul_rightQ : phaseCZLeftQ * phaseCZRightQ = phaseCZKernelQ := by
  ext i j
  simp only [Matrix.mul_apply, phaseCZKernelQ, phaseCZTaylorQ_eq_sum,
    phaseCZLeftQ, phaseCZRightQ, mul_pow]
  apply Finset.sum_congr rfl
  intro k _
  ring

def phaseCZConditionedLeftQ : Matrix (Fin 4) (Fin 9) ℚ :=
  companionRowUpdate phaseCZFirstResidualQ
    (companionRowUpdate phaseCZKernelQ phaseCZLeftQ 3 3
      phaseCZKernelQ_pivot_ne) 0 0 phaseCZFirstResidualQ_pivot_ne

def phaseCZConditionedRightQ : Matrix (Fin 9) (Fin 4) ℚ :=
  companionColumnUpdate phaseCZFirstResidualQ
    (companionColumnUpdate phaseCZKernelQ phaseCZRightQ 3 3
      phaseCZKernelQ_pivot_ne) 0 0 phaseCZFirstResidualQ_pivot_ne

/-- The actual two-sided conditioned feature block on the two fresh rows
and columns. Both factor updates use the original selected pivot residuals. -/
def phaseCZBlockQ (degrees : Finset (Fin 9)) : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j ↦ ∑ k ∈ degrees,
    phaseCZConditionedLeftQ ⟨i.1 + 1, by omega⟩ k *
      phaseCZConditionedRightQ k ⟨j.1 + 1, by omega⟩

def phaseCZLowQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {0, 1, 2}
def phaseCZMovedQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {3}
def phaseCZUpperQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {4, 5, 6, 7, 8}

def phaseCZLowerTwoQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {0, 1}
def phaseCZMovedTwoQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {2}
def phaseCZUpperTwoQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {3, 4, 5, 6, 7, 8}
def phaseCZLowerFourQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {0, 1, 2, 3}
def phaseCZMovedFourQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {4}
def phaseCZUpperFourQ : Matrix (Fin 2) (Fin 2) ℚ := phaseCZBlockQ {5, 6, 7, 8}

def phaseCZMixedQ (low high : Matrix (Fin 2) (Fin 2) ℚ) : ℚ :=
  low 0 0 * high 1 1 + high 0 0 * low 1 1 -
    low 0 1 * high 1 0 - high 0 1 * low 1 0

def phaseCZTransferQ
    (low moved upper : Matrix (Fin 2) (Fin 2) ℚ) : ℚ :=
  low.det * phaseCZMixedQ moved upper ^ 2 -
    phaseCZMixedQ low (moved + upper) *
      phaseCZMixedQ low moved * phaseCZMixedQ moved upper +
    (moved + upper).det * phaseCZMixedQ low moved ^ 2

/-- Rational normalization computes exactly the real transfer-node quantity
used by the generic Phase CY theorem. -/
theorem phaseCZTransferQ_cast
    (low moved upper : Matrix (Fin 2) (Fin 2) ℚ) :
    (phaseCZTransferQ low moved upper : ℝ) =
      quadraticTransferNodeNumerator
        (low.map (fun x ↦ (x : ℝ))).det
        (finTwoPencilMixed (low.map (fun x ↦ (x : ℝ)))
          ((moved + upper).map (fun x ↦ (x : ℝ))))
        ((moved + upper).map (fun x ↦ (x : ℝ))).det
        (finTwoPencilMixed (low.map (fun x ↦ (x : ℝ)))
          (moved.map (fun x ↦ (x : ℝ))))
        (finTwoPencilMixed (moved.map (fun x ↦ (x : ℝ)))
          (upper.map (fun x ↦ (x : ℝ)))) := by
  simp only [phaseCZTransferQ, phaseCZMixedQ, quadraticTransferNodeNumerator,
    finTwoPencilMixed, Matrix.det_fin_two, Matrix.map_apply]
  push_cast
  rfl

local macro "phase_cz_compute" : tactic => `(tactic|
  (simp [phaseCZTransferQ, phaseCZMixedQ, phaseCZLowQ, phaseCZMovedQ,
    phaseCZUpperQ, phaseCZLowerTwoQ, phaseCZMovedTwoQ, phaseCZUpperTwoQ,
    phaseCZLowerFourQ, phaseCZMovedFourQ, phaseCZUpperFourQ,
    phaseCZBlockQ, phaseCZConditionedLeftQ,
    phaseCZConditionedRightQ, companionRowUpdate, companionColumnUpdate,
    phaseCZFirstResidualQ, residualUpdate, phaseCZLeftQ, phaseCZRightQ,
    phaseCZKernelQ, phaseCZTaylorQ, phaseCZRowsQ, phaseCZColumnsQ,
    Matrix.det_fin_two, cons_val_two, cons_val_three, vecHead, vecTail]
   norm_num))

set_option maxHeartbeats 4000000 in
/-- The shared-node certificate can fail at a strict interior degree. -/
theorem phaseCZ_interior_transferNode_pos :
    0 < phaseCZTransferQ phaseCZLowQ phaseCZMovedQ phaseCZUpperQ := by
  phase_cz_compute

set_option maxHeartbeats 4000000 in
/-- All four endpoints and both transfer sums are strictly positive; the
failed canonical certificate is not a rank-one endpoint phenomenon. -/
theorem phaseCZ_endpoints_and_transfers_pos :
    0 < phaseCZLowQ.det ∧
    0 < (phaseCZMovedQ + phaseCZUpperQ).det ∧
    0 < (phaseCZLowQ + phaseCZMovedQ).det ∧
    0 < phaseCZUpperQ.det ∧
    0 < phaseCZMixedQ phaseCZLowQ phaseCZMovedQ ∧
    0 < phaseCZMixedQ phaseCZMovedQ phaseCZUpperQ := by
  phase_cz_compute

set_option maxHeartbeats 4000000 in
/-- The nodes from moving degrees two and four supply valid alternatives on
either side of the failed degree-three shared-node certificate. -/
theorem phaseCZ_alternate_transferNodes_neg :
    phaseCZTransferQ phaseCZLowerTwoQ phaseCZMovedTwoQ phaseCZUpperTwoQ < 0 ∧
      phaseCZTransferQ phaseCZLowerFourQ phaseCZMovedFourQ phaseCZUpperFourQ < 0 := by
  phase_cz_compute

set_option maxHeartbeats 4000000 in
/-- Both neighboring cutoff discriminants remain strictly positive despite
failure of the degree-three canonical-node sign. -/
theorem phaseCZ_neighbor_discriminants_pos :
    0 < phaseCZMixedQ phaseCZLowQ (phaseCZMovedQ + phaseCZUpperQ) ^ 2 -
        4 * phaseCZLowQ.det * (phaseCZMovedQ + phaseCZUpperQ).det ∧
      0 < phaseCZMixedQ (phaseCZLowQ + phaseCZMovedQ) phaseCZUpperQ ^ 2 -
        4 * (phaseCZLowQ + phaseCZMovedQ).det * phaseCZUpperQ.det := by
  phase_cz_compute

end GECPKernelStructure.Fermionic
