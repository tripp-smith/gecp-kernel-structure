import GECPKernelStructure.Fermionic.CompleteRun
import GECPKernelStructure.GECP.Determinant

namespace GECPKernelStructure

namespace GECP

universe u v

private theorem finset_prod_range_eq_list_prod_range (f : ℕ → ℝ) (n : ℕ) :
    (∏ i ∈ Finset.range n, f i) = ((List.range n).map f).prod := by
  induction n with
  | zero => simp
  | succ n ih => simp [Finset.prod_range_succ, List.range_succ, ih]

namespace Run

/-- A successful run stores exactly the first `n` pivots of a residual sequence. -/
def RealizesPivotPrefix {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (residual : ℕ → Kernel α β ℝ)
    (rows : ℕ → α) (columns : ℕ → β) (n : ℕ) : Prop :=
  run.pivots = (List.range n).map fun i => residual i (rows i) (columns i)

/-- The absolute selected-core determinant is the product of the absolute pivots. -/
theorem abs_finSelectedCore_det_eq_prod_abs_pivots
    {α : Type u} {β : Type v} {K : Kernel α β ℝ} (run : Run K) :
    abs run.finSelectedCore.det = (run.pivots.map abs).prod := by
  rw [gecp_core_det_eq_prod_pivots]
  induction run.pivots with
  | nil => simp
  | cons pivot pivots ih => simp [abs_mul, ih]

end Run

end GECP

namespace Fermionic

open GECP

universe u v

/-- Exact complete-pivot magnitudes are antitone for a realized strictly sign-regular sequence. -/
theorem strictSignRegular_gecp_pivotMagnitude_antitone
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (residual : ℕ → Kernel α β ℝ) (rows : ℕ → α) (columns : ℕ → β)
    (realized : ∀ n, ∃ run : Run K, run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn rowDomain columnDomain (residual n) (rows n) (columns n)) :
    Antitone fun n => abs (residual n (rows n) (columns n)) := by
  apply antitone_nat_of_succ_le
  intro n
  obtain ⟨run, run_residual⟩ := realized n
  have coherent :
      PivotCrossProductSignCoherent (residual n) (rows n) (columns n) := by
    rw [← run_residual]
    exact strictSignRegular_pivotCrossProductSignCoherent
      regular run (rows n) (columns n)
  obtain ⟨row_mem, column_mem, pivot_max⟩ := complete n
  obtain ⟨next_row_mem, next_column_mem, _⟩ := complete (n + 1)
  rw [updates n]
  apply residualUpdate_le_of_signCoherentOn coherent
  · exact pivot_max (rows (n + 1)) (columns (n + 1))
      next_row_mem next_column_mem
  · rfl
  · exact pivot_max (rows (n + 1)) (columns n) next_row_mem column_mem
  · exact pivot_max (rows n) (columns (n + 1)) row_mem next_column_mem

/-- The current complete-pivot magnitude to the `n`th power is below the prior pivot product. -/
theorem strictSignRegular_gecp_pivotMagnitude_pow_le_product
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (residual : ℕ → Kernel α β ℝ) (rows : ℕ → α) (columns : ℕ → β)
    (realized : ∀ n, ∃ run : Run K, run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn rowDomain columnDomain (residual n) (rows n) (columns n))
    (n : ℕ) :
    abs (residual n (rows n) (columns n)) ^ n ≤
      ∏ i ∈ Finset.range n, abs (residual i (rows i) (columns i)) := by
  have antitone := strictSignRegular_gecp_pivotMagnitude_antitone regular
    rowDomain columnDomain residual rows columns realized pivot_ne updates complete
  calc
    abs (residual n (rows n) (columns n)) ^ n =
        ∏ _i ∈ Finset.range n, abs (residual n (rows n) (columns n)) := by simp
    _ ≤ ∏ i ∈ Finset.range n, abs (residual i (rows i) (columns i)) := by
      apply Finset.prod_le_prod
      · intro i hi
        exact abs_nonneg _
      · intro i hi
        exact antitone (Nat.le_of_lt (Finset.mem_range.mp hi))

/--
The `n`th power of every rank-`n` domain residual entry is bounded by the
absolute determinant of any successful run recording the first `n` pivots.
-/
theorem strictSignRegular_gecp_error_pow_le_selectedCore_det
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (residual : ℕ → Kernel α β ℝ) (rows : ℕ → α) (columns : ℕ → β)
    (realized : ∀ n, ∃ run : Run K, run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn rowDomain columnDomain (residual n) (rows n) (columns n))
    (n : ℕ) (run : Run K)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (x : α) (y : β) (hx : rowDomain x) (hy : columnDomain y) :
    abs (residual n x y) ^ n ≤ abs run.finSelectedCore.det := by
  have entry_le_pivot := (complete n).2.2 x y hx hy
  have entry_pow_le_pivot_pow :=
    pow_le_pow_left₀ (abs_nonneg _) entry_le_pivot n
  have pivot_pow_le_product :=
    strictSignRegular_gecp_pivotMagnitude_pow_le_product regular
      rowDomain columnDomain residual rows columns realized pivot_ne updates complete n
  calc
    abs (residual n x y) ^ n ≤
        abs (residual n (rows n) (columns n)) ^ n := entry_pow_le_pivot_pow
    _ ≤ ∏ i ∈ Finset.range n, abs (residual i (rows i) (columns i)) :=
      pivot_pow_le_product
    _ = abs run.finSelectedCore.det := by
      rw [Run.abs_finSelectedCore_det_eq_prod_abs_pivots, run_pivots]
      rw [List.map_map]
      exact finset_prod_range_eq_list_prod_range
        (abs ∘ fun i => residual i (rows i) (columns i)) n

end Fermionic
end GECPKernelStructure
