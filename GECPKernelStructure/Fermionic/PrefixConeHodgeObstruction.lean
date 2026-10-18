import GECPKernelStructure.Fermionic.CompleteRun
import GECPKernelStructure.GECP.CompanionBiConditioning
import GECPKernelStructure.Matrix.RankTwoRatio

namespace GECPKernelStructure.Fermionic

open GECP Matrix

def phaseCXBaseQ : ℚ := 21 / 20

/-- A small exact exponential kernel whose first complete pivot is its
upper-right corner. -/
def phaseCXKernelQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  phaseCXBaseQ ^ ((i.1 + 1) * (j.1 + 1))

noncomputable def phaseCXKernel : Kernel (Fin 4) (Fin 4) ℝ := fun i j =>
  (phaseCXKernelQ i j : ℚ)

/-- Physical row nodes realizing `phaseCXKernelQ` as an exponential kernel. -/
noncomputable def phaseCXRows : Fin 4 → ℝ := fun i => (i.1 + 1) / 5

/-- Physical column nodes realizing `phaseCXKernelQ` as an exponential kernel. -/
noncomputable def phaseCXColumns : Fin 4 → ℝ := fun j =>
  5 * Real.log (21 / 20) * (j.1 + 1)

theorem phaseCXRows_mem_openUnit (i : Fin 4) :
    phaseCXRows i ∈ Set.Ioo (0 : ℝ) 1 := by
  fin_cases i <;> norm_num [phaseCXRows]

theorem phaseCXColumns_mem_openUnit (j : Fin 4) :
    phaseCXColumns j ∈ Set.Ioo (0 : ℝ) 1 := by
  have log_pos : 0 < Real.log (21 / 20) := Real.log_pos (by norm_num)
  have log_lt : Real.log (21 / 20) < (21 / 20 : ℝ) - 1 :=
    Real.log_lt_sub_one_of_pos (by norm_num) (by norm_num)
  fin_cases j <;> constructor <;> simp only [phaseCXColumns] <;> norm_num <;>
    nlinarith

theorem phaseCX_expKernel_eq (i j : Fin 4) :
    Real.exp (phaseCXRows i * phaseCXColumns j) =
      phaseCXKernel i j := by
  rw [show phaseCXRows i * phaseCXColumns j =
      ((i.1 + 1) * (j.1 + 1) : ℕ) * Real.log (21 / 20) by
    simp only [phaseCXRows, phaseCXColumns]
    push_cast
    ring]
  rw [Real.exp_nat_mul, Real.exp_log (by norm_num : (0 : ℝ) < 21 / 20)]
  norm_num [phaseCXKernel, phaseCXKernelQ, phaseCXBaseQ]

theorem phaseCXKernelQ_pivot_ne : phaseCXKernelQ 3 3 ≠ 0 := by
  norm_num [phaseCXKernelQ, phaseCXBaseQ]

theorem phaseCX_firstPivot_complete :
    CompletePivotOn (fun _ : Fin 4 => True) (fun _ : Fin 4 => True)
      phaseCXKernel 3 3 := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  fin_cases i <;> fin_cases j <;>
    norm_num [phaseCXKernel, phaseCXKernelQ, phaseCXBaseQ, abs_of_nonneg]

/-- The factorial-normalized monomial factor before the exact pivot update. -/
def phaseCXLeftQ : Matrix (Fin 4) (Fin 4) ℚ := fun i degree =>
  (i.1 + 1 : ℚ) ^ degree.1 / degree.1.factorial

/-- The matching monomial factor before the exact pivot update. -/
def phaseCXRightQ : Matrix (Fin 4) (Fin 4) ℚ := fun degree j =>
  (j.1 + 1 : ℚ) ^ degree.1

def phaseCXConditionedLeftQ : Matrix (Fin 4) (Fin 4) ℚ :=
  companionRowUpdate phaseCXKernelQ phaseCXLeftQ 3 3 phaseCXKernelQ_pivot_ne

def phaseCXConditionedRightQ : Matrix (Fin 4) (Fin 4) ℚ :=
  companionColumnUpdate phaseCXKernelQ phaseCXRightQ 3 3
    phaseCXKernelQ_pivot_ne

/-- The matched-coordinate Hessian for fresh rows `(0,1)`, fresh columns
`(1,2)`, and the first four normalized Taylor degrees after the complete
pivot `(3,3)`. -/
def phaseCXMatchedHessianQ : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  if i = j then 0 else
    (phaseCXConditionedLeftQ 0 i * phaseCXConditionedLeftQ 1 j -
        phaseCXConditionedLeftQ 0 j * phaseCXConditionedLeftQ 1 i) *
      (phaseCXConditionedRightQ i 1 * phaseCXConditionedRightQ j 2 -
        phaseCXConditionedRightQ i 2 * phaseCXConditionedRightQ j 1)

/-- Bilinear form of the matched-coordinate Hessian. -/
def phaseCXBilinearQ (x y : Fin 4 → ℚ) : ℚ :=
  ∑ i, ∑ j, x i * phaseCXMatchedHessianQ i j * y j

def phaseCXPositiveVectorQ : Fin 4 → ℚ := ![10, 177, 686, 706]

def phaseCXSecondPositiveVectorQ : Fin 4 → ℚ := ![999, -38, -6, 0]

/-- The normalized conditioned Hessian has a positive-definite restriction to
an explicit two-dimensional plane.  This is an exact obstruction to any
global one-positive-direction (Hodge--Riemann/Lorentzian) strengthening. -/
theorem phaseCX_positivePlane :
    0 < phaseCXBilinearQ phaseCXPositiveVectorQ phaseCXPositiveVectorQ ∧
    0 < phaseCXBilinearQ phaseCXSecondPositiveVectorQ
      phaseCXSecondPositiveVectorQ ∧
    phaseCXBilinearQ phaseCXPositiveVectorQ phaseCXSecondPositiveVectorQ ^ 2 <
      phaseCXBilinearQ phaseCXPositiveVectorQ phaseCXPositiveVectorQ *
        phaseCXBilinearQ phaseCXSecondPositiveVectorQ
          phaseCXSecondPositiveVectorQ := by
  simp [phaseCXBilinearQ, phaseCXMatchedHessianQ,
    phaseCXConditionedLeftQ, phaseCXConditionedRightQ,
    companionRowUpdate, companionColumnUpdate, phaseCXLeftQ, phaseCXRightQ,
    phaseCXKernelQ, phaseCXBaseQ, phaseCXPositiveVectorQ,
    phaseCXSecondPositiveVectorQ, Fin.sum_univ_four, cons_val_two,
    cons_val_three, vecHead, vecTail]
  norm_num

/-- Reverse Cauchy--Schwarz on every pair of positive vectors is a convenient
coordinate-free necessary condition for a one-positive-direction quadratic
form. -/
def ReverseCauchySchwarzOnPositive
    (B : (Fin 4 → ℚ) → (Fin 4 → ℚ) → ℚ) : Prop :=
  ∀ x y, 0 < B x x → 0 < B y y → B x x * B y y ≤ B x y ^ 2

theorem phaseCX_not_reverseCauchySchwarzOnPositive :
    ¬ ReverseCauchySchwarzOnPositive phaseCXBilinearQ := by
  intro h
  rcases phaseCX_positivePlane with ⟨hx, hy, hxy⟩
  exact (not_le_of_gt hxy) (h _ _ hx hy)

def phaseCXPrefixOneQ : Fin 4 → ℚ := ![1, 0, 0, 0]
def phaseCXSuffixOneQ : Fin 4 → ℚ := ![0, 1, 1, 1]
def phaseCXPrefixTwoQ : Fin 4 → ℚ := ![1, 1, 0, 0]
def phaseCXSuffixTwoQ : Fin 4 → ℚ := ![0, 0, 1, 1]
def phaseCXPrefixThreeQ : Fin 4 → ℚ := ![1, 1, 1, 0]
def phaseCXSuffixThreeQ : Fin 4 → ℚ := ![0, 0, 0, 1]

/-- In the same exact witness, every nontrivial contiguous prefix/suffix plane
still obeys the reverse Cauchy--Schwarz inequality needed by the determinant
pencil.  The surviving structure is therefore prefix-cone, not global. -/
theorem phaseCX_all_prefixPlanes_lorentzian :
    phaseCXBilinearQ phaseCXPrefixOneQ phaseCXSuffixOneQ ^ 2 ≥
        phaseCXBilinearQ phaseCXPrefixOneQ phaseCXPrefixOneQ *
          phaseCXBilinearQ phaseCXSuffixOneQ phaseCXSuffixOneQ ∧
    phaseCXBilinearQ phaseCXPrefixTwoQ phaseCXSuffixTwoQ ^ 2 ≥
        phaseCXBilinearQ phaseCXPrefixTwoQ phaseCXPrefixTwoQ *
          phaseCXBilinearQ phaseCXSuffixTwoQ phaseCXSuffixTwoQ ∧
    phaseCXBilinearQ phaseCXPrefixThreeQ phaseCXSuffixThreeQ ^ 2 ≥
        phaseCXBilinearQ phaseCXPrefixThreeQ phaseCXPrefixThreeQ *
          phaseCXBilinearQ phaseCXSuffixThreeQ phaseCXSuffixThreeQ := by
  simp [phaseCXBilinearQ, phaseCXMatchedHessianQ,
    phaseCXConditionedLeftQ, phaseCXConditionedRightQ,
    companionRowUpdate, companionColumnUpdate, phaseCXLeftQ, phaseCXRightQ,
    phaseCXKernelQ, phaseCXBaseQ, phaseCXPrefixOneQ, phaseCXSuffixOneQ,
    phaseCXPrefixTwoQ, phaseCXSuffixTwoQ, phaseCXPrefixThreeQ,
    phaseCXSuffixThreeQ, Fin.sum_univ_four, cons_val_two, cons_val_three,
    vecHead, vecTail]
  norm_num

end GECPKernelStructure.Fermionic
