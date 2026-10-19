import GECPKernelStructure.Fermionic.AdaptiveCoverageObstructionData
import GECPKernelStructure.Fermionic.CompleteRun
import GECPKernelStructure.GECP.CompanionBiConditioning
import GECPKernelStructure.Matrix.QuadraticTransfer

/-! A coherent finite-Taylor obstruction to neighboring canonical-node coverage.
Integer fraction-free tables certify the seven actual complete-pivot updates;
they are used to keep kernel checking independent of expression swell. -/

namespace GECPKernelStructure.Fermionic

open GECP Matrix

-- Closed certificate entries have up to 305 decimal digits before multiplication.
set_option maxRecDepth 16384

def phaseDBRowsQ (i : Fin 9) : ℚ := phaseDBRowNumerators i / 29
def phaseDBColumnsQ (j : Fin 9) : ℚ := phaseDBColumnNumerators j / 63

def phaseDBRawLeft (i : Fin 9) (k : Fin 13) : ℤ :=
  phaseDBRowNumerators i ^ k.1 * (29 * 63) ^ (12 - k.1) *
    (Nat.factorial 12 / Nat.factorial k.1 : ℕ)

def phaseDBRawRight (k : Fin 13) (j : Fin 9) : ℤ :=
  phaseDBColumnNumerators j ^ k.1

def phaseDBResidualQ (k : Fin 8) : Kernel (Fin 9) (Fin 9) ℚ := fun i j ↦
  (phaseDBState k i j : ℚ) / (phaseDBDivisors k : ℚ)

theorem phaseDB_nodes_openUnit :
    (∀ i, 0 < phaseDBRowsQ i ∧ phaseDBRowsQ i < 1) ∧
      (∀ j, 0 < phaseDBColumnsQ j ∧ phaseDBColumnsQ j < 1) := by
  decide +kernel

theorem phaseDB_nodes_strictMono : StrictMono phaseDBRowsQ ∧ StrictMono phaseDBColumnsQ := by
  decide +kernel

theorem phaseDB_kernelScale_pos : 0 < phaseDBKernelScale := by
  decide +kernel

theorem phaseDB_divisors_pos (k : Fin 8) : 0 < phaseDBDivisors k := by
  revert k
  decide +kernel

theorem phaseDB_pivot_eq_next_divisor (k : Fin 7) :
    phaseDBState k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) =
      phaseDBDivisors k.succ := by
  revert k
  decide +kernel

set_option maxHeartbeats 20000000 in
theorem phaseDB_fractionFree_recurrence (k : Fin 7) (i j : Fin 9) :
    phaseDBDivisors k.castSucc * phaseDBState k.succ i j =
      phaseDBState k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) *
          phaseDBState k.castSucc i j -
        phaseDBState k.castSucc i (phaseDBPivotColumns k) *
          phaseDBState k.castSucc (phaseDBPivotRows k) j := by
  revert k i j
  decide +kernel

set_option maxHeartbeats 20000000 in
theorem phaseDB_fractionFree_complete (k : Fin 7) (i j : Fin 9) :
    |phaseDBState k.castSucc i j| ≤
      phaseDBState k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) := by
  revert k i j
  decide +kernel

theorem phaseDB_pivot_ne (k : Fin 7) :
    phaseDBResidualQ k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) ≠ 0 := by
  unfold phaseDBResidualQ
  rw [phaseDB_pivot_eq_next_divisor]
  exact div_ne_zero (by exact_mod_cast (phaseDB_divisors_pos k.succ).ne')
    (by exact_mod_cast (phaseDB_divisors_pos k.castSucc).ne')

private theorem fractionFree_residualUpdate
    (current next : Kernel (Fin 9) (Fin 9) ℚ) (previous : ℚ)
    (row column : Fin 9) (previous_ne : previous ≠ 0)
    (pivot_ne : current row column ≠ 0)
    (identity : ∀ i j, previous * next i j =
      current row column * current i j - current i column * current row j) :
    residualUpdate (fun i j ↦ current i j / previous) row column
        (div_ne_zero pivot_ne previous_ne) =
      fun i j ↦ next i j / current row column := by
  funext i j
  have hn : next i j =
      (current row column * current i j - current i column * current row j) /
        previous := (eq_div_iff previous_ne).mpr (by simpa [mul_comm] using identity i j)
  rw [hn]
  dsimp [residualUpdate]
  field_simp

theorem phaseDB_residual_recurrence (k : Fin 7) :
    phaseDBResidualQ k.succ = residualUpdate (phaseDBResidualQ k.castSucc)
      (phaseDBPivotRows k) (phaseDBPivotColumns k) (phaseDB_pivot_ne k) := by
  have identity : ∀ i j,
      (phaseDBDivisors k.castSucc : ℚ) * (phaseDBState k.succ i j : ℚ) =
        (phaseDBState k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) : ℚ) *
          (phaseDBState k.castSucc i j : ℚ) -
        (phaseDBState k.castSucc i (phaseDBPivotColumns k) : ℚ) *
          (phaseDBState k.castSucc (phaseDBPivotRows k) j : ℚ) := by
    intro i j
    exact_mod_cast phaseDB_fractionFree_recurrence k i j
  have h := fractionFree_residualUpdate
    (fun i j ↦ (phaseDBState k.castSucc i j : ℚ))
    (fun i j ↦ (phaseDBState k.succ i j : ℚ)) (phaseDBDivisors k.castSucc : ℚ)
    (phaseDBPivotRows k) (phaseDBPivotColumns k)
    (by exact_mod_cast (phaseDB_divisors_pos k.castSucc).ne')
    (by rw [phaseDB_pivot_eq_next_divisor]; exact_mod_cast (phaseDB_divisors_pos k.succ).ne')
    identity
  unfold phaseDBResidualQ
  simpa only [phaseDB_pivot_eq_next_divisor] using h.symm

theorem phaseDB_complete_pivots (k : Fin 7) :
    CompletePivotOn (fun _ : Fin 9 ↦ True) (fun _ : Fin 9 ↦ True)
      (fun i j ↦ (phaseDBResidualQ k.castSucc i j : ℝ))
      (phaseDBPivotRows k) (phaseDBPivotColumns k) := by
  refine ⟨trivial, trivial, ?_⟩
  intro i j _ _
  have denom_pos : (0 : ℝ) < (phaseDBDivisors k.castSucc : ℝ) := by
    exact_mod_cast phaseDB_divisors_pos k.castSucc
  have pivot_pos : (0 : ℝ) <
      (phaseDBState k.castSucc (phaseDBPivotRows k) (phaseDBPivotColumns k) : ℝ) := by
    rw [phaseDB_pivot_eq_next_divisor]
    exact_mod_cast phaseDB_divisors_pos k.succ
  simp only [phaseDBResidualQ, Rat.cast_div, Rat.cast_intCast, abs_div,
    abs_of_pos denom_pos, abs_of_pos pivot_pos]
  exact div_le_div_of_nonneg_right (by exact_mod_cast phaseDB_fractionFree_complete k i j)
    denom_pos.le

set_option maxHeartbeats 20000000 in
theorem phaseDB_raw_feature_product (i j : Fin 9) :
    phaseDBState 0 i j = ∑ k : Fin 13, phaseDBRawLeft i k * phaseDBRawRight k j := by
  fin_cases i <;> fin_cases j <;> decide +kernel

theorem phaseDB_raw_feature_is_Taylor (i j : Fin 9) (k : Fin 13) :
    ((phaseDBRawLeft i k : ℤ) : ℚ) * (phaseDBRawRight k j : ℚ) =
      (phaseDBKernelScale : ℚ) * (phaseDBRowsQ i * phaseDBColumnsQ j) ^ k.1 /
        (k.1.factorial : ℚ) := by
  fin_cases k <;>
    norm_num [phaseDBRawLeft, phaseDBRawRight, phaseDBKernelScale,
      phaseDBRowsQ, phaseDBColumnsQ] <;> push_cast <;> ring

theorem phaseDB_coherent_Taylor_kernel (i j : Fin 9) :
    phaseDBResidualQ 0 i j = (phaseDBKernelScale : ℚ) *
      ∑ k : Fin 13, (phaseDBRowsQ i * phaseDBColumnsQ j) ^ k.1 / (k.1.factorial : ℚ) := by
  simp only [phaseDBResidualQ, show phaseDBDivisors 0 = 1 by rfl, Int.cast_one, div_one]
  rw [phaseDB_raw_feature_product]
  push_cast
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ ↦ by simpa [mul_div_assoc] using phaseDB_raw_feature_is_Taylor i j k

end GECPKernelStructure.Fermionic
