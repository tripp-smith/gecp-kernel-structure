import GECPKernelStructure.Fermionic.TruncatedPowerConvolution
import GECPKernelStructure.GECP.SignRegularity

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The strict order kernel at the base of the normalized truncated-power ladder. -/
noncomputable def strictOrderKernel (v ω : ℝ) : ℝ :=
  if v < ω then 1 else 0

/-- Normalized degree zero is exactly the strict order kernel. -/
theorem normalizedTruncatedPowerKernel_zero_eq_strictOrderKernel (v ω : ℝ) :
    normalizedTruncatedPowerKernel 0 v ω = strictOrderKernel v ω := by
  simp [normalizedTruncatedPowerKernel_zero, strictOrderKernel]

/-- Every ordered square minor of the strict order kernel is exactly zero or one. -/
theorem strictOrderKernel_minor_det_eq_zero_or_one
    (n : ℕ) (rows columns : Fin n → ℝ) (rows_mono : StrictMono rows)
    (columns_mono : StrictMono columns) :
    (minorMatrix strictOrderKernel rows columns).det = 0 ∨
      (minorMatrix strictOrderKernel rows columns).det = 1 := by
  induction n with
  | zero =>
      right
      simp
  | succ n ih =>
      cases n with
      | zero =>
          by_cases h : rows 0 < columns 0
          · right
            simp [minorMatrix, strictOrderKernel, h]
          · left
            simp [minorMatrix, strictOrderKernel, h]
      | succ n =>
          let A : Matrix (Fin (n + 2)) (Fin (n + 2)) ℝ :=
            minorMatrix strictOrderKernel rows columns
          by_cases h00 : rows 0 < columns 0
          · by_cases h10 : rows 1 < columns 0
            · left
              apply Matrix.det_zero_of_row_eq (show (0 : Fin (n + 2)) ≠ 1 by simp)
              funext j
              have hcj : columns 0 ≤ columns j :=
                columns_mono.monotone (Fin.zero_le j)
              have h0j : rows 0 < columns j := h00.trans_le hcj
              have h1j : rows 1 < columns j := h10.trans_le hcj
              simp [minorMatrix, strictOrderKernel, h0j, h1j]
            · have hcol : ∀ i : Fin (n + 1), A i.succ 0 = 0 := by
                intro i
                have hr1i : rows 1 ≤ rows i.succ := by
                  apply rows_mono.monotone
                  change (1 : ℕ) ≤ i.val + 1
                  omega
                have hle : columns 0 ≤ rows i.succ :=
                  (le_of_not_gt h10).trans hr1i
                simp [A, minorMatrix, strictOrderKernel, not_lt_of_ge hle]
              have hdet : A.det =
                  (A.submatrix Fin.succ Fin.succ).det := by
                rw [Matrix.det_succ_column_zero, Fin.sum_univ_succ]
                simp only [Fin.val_zero, pow_zero, one_mul]
                rw [Finset.sum_eq_zero]
                · simp [A, minorMatrix, strictOrderKernel, h00]
                · intro i _
                  rw [hcol]
                  ring
              have htail : A.submatrix Fin.succ Fin.succ =
                  minorMatrix strictOrderKernel (rows ∘ Fin.succ)
                    (columns ∘ Fin.succ) := by
                rfl
              rw [show minorMatrix strictOrderKernel rows columns = A by rfl, hdet]
              rw [htail]
              exact ih (rows ∘ Fin.succ) (columns ∘ Fin.succ)
                (rows_mono.comp Fin.strictMono_succ)
                (columns_mono.comp Fin.strictMono_succ)
          · left
            apply Matrix.det_eq_zero_of_column_eq_zero 0
            intro i
            have hle : columns 0 ≤ rows i :=
              (le_of_not_gt h00).trans (rows_mono.monotone (Fin.zero_le i))
            simp [minorMatrix, strictOrderKernel, not_lt_of_ge hle]

/-- The strict order kernel is totally nonnegative at every order. -/
theorem strictOrderKernel_totallyNonnegative : TotallyNonnegative strictOrderKernel := by
  intro n rows columns rows_mono columns_mono
  rcases strictOrderKernel_minor_det_eq_zero_or_one n rows columns
    rows_mono columns_mono with h | h <;> linarith

/-- Normalized degree zero is totally nonnegative at every order. -/
theorem normalizedTruncatedPowerKernel_zero_totallyNonnegative :
    TotallyNonnegative (normalizedTruncatedPowerKernel 0) := by
  intro n rows columns rows_mono columns_mono
  rw [show minorMatrix (normalizedTruncatedPowerKernel 0) rows columns =
      minorMatrix strictOrderKernel rows columns by
    ext i j
    exact normalizedTruncatedPowerKernel_zero_eq_strictOrderKernel _ _]
  exact strictOrderKernel_totallyNonnegative n rows columns rows_mono columns_mono

end Fermionic
end GECPKernelStructure
