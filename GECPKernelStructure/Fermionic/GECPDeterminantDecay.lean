import GECPKernelStructure.Fermionic.DeterminantDecay
import GECPKernelStructure.Fermionic.GeometricMean
import Mathlib.Data.Nat.Log

namespace GECPKernelStructure

namespace GECP

universe u v

namespace Run

/-- Every row selected by a domain-complete run belongs to the row domain. -/
theorem selectedRow_mem_of_completeOn {α : Type u} {β : Type v}
    {rowDomain : α → Prop} {columnDomain : β → Prop} {K : Kernel α β ℝ}
    (run : Run K) (complete : run.CompleteOn rowDomain columnDomain) :
    ∀ i, rowDomain (run.selectedRow i) := by
  induction run with
  | nil K =>
      intro i
      exact Fin.elim0 i
  | step row column pivot_ne tail ih =>
      rcases complete with ⟨pivot_complete, tail_complete⟩
      intro i
      cases i with
      | inl i => exact pivot_complete.1
      | inr i => exact ih tail_complete i

/-- Every column selected by a domain-complete run belongs to the column domain. -/
theorem selectedColumn_mem_of_completeOn {α : Type u} {β : Type v}
    {rowDomain : α → Prop} {columnDomain : β → Prop} {K : Kernel α β ℝ}
    (run : Run K) (complete : run.CompleteOn rowDomain columnDomain) :
    ∀ i, columnDomain (run.selectedColumn i) := by
  induction run with
  | nil K =>
      intro i
      exact Fin.elim0 i
  | step row column pivot_ne tail ih =>
      rcases complete with ⟨pivot_complete, tail_complete⟩
      intro i
      cases i with
      | inl i => exact pivot_complete.2.1
      | inr i => exact ih tail_complete i

/-- A run's selected core is the original kernel sampled at its recursive pivot coordinates. -/
theorem selectedCore_eq_sampleMatrix {α : Type u} {β : Type v}
    {K : Kernel α β ℝ} (run : Run K) :
    run.selectedCore =
      Fermionic.sampleMatrix K run.selectedRow run.selectedColumn := rfl

end Run
end GECP

namespace Fermionic

open GECP

private theorem card_selectedIndex_eq_prefixLength
    {K : Kernel ℝ ℝ ℝ} (residual : ℕ → Kernel ℝ ℝ ℝ)
    (rows columns : ℕ → ℝ) (n : ℕ) (run : Run K)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n) :
    Fintype.card run.SelectedIndex = n := by
  have hpivot_length : run.pivots.length = n := by
    rw [run_pivots]
    simp
  calc
    Fintype.card run.SelectedIndex = Fintype.card (Fin run.pivots.length) :=
      Fintype.card_congr run.selectedIndexEquivFin
    _ = run.pivots.length := Fintype.card_fin _
    _ = n := hpivot_length

/--
Exact composition of the complete-pivot residual-power inequality with the
surviving-subset determinant estimate for the fermionic separated approximation.
-/
theorem fermionicKernel_gecp_error_pow_le_separatedApprox
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (p s : ℕ) (hp : 0 < p)
    (realized : ∀ m, ∃ run : Run fermionicKernel, run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual m) (rows m) (columns m))
    (n : ℕ) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ)
    (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual n x y| ^ n ≤
      ∑ S : Finset run.SelectedIndex,
        if n ≤ 2 * ((s + 1) * (8 * p)) + S.card then
          n.factorial * (1 / 2 : ℝ) ^ (p * S.card) * 2 ^ (n - S.card)
        else 0 := by
  have hpower := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
    residual rows columns realized pivot_ne updates complete n run run_pivots x y hx hy
  rw [Run.finSelectedCore, Matrix.det_reindex_self] at hpower
  have hrows := run.selectedRow_mem_of_completeOn run_complete
  have hcolumns := run.selectedColumn_mem_of_completeOn run_complete
  have hdet := fermionicKernel_sample_det_le_separatedApprox p s hp
    run.selectedRow run.selectedColumn
    (fun i => (hrows i).1) (fun i => (hrows i).2)
    (fun j => (hcolumns j).1) (fun j => (hcolumns j).2)
  rw [← Run.selectedCore_eq_sampleMatrix] at hdet
  have hcard := card_selectedIndex_eq_prefixLength residual rows columns n run run_pivots
  exact hpower.trans (by simpa [hcard] using hdet)

/-- Closed power-form GECP residual bound obtained by composing Phases T and U. -/
theorem fermionicKernel_gecp_error_pow_le_two_pow
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (p s : ℕ) (hp : 0 < p)
    (realized : ∀ m, ∃ run : Run fermionicKernel, run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual m) (rows m) (columns m))
    (n : ℕ) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (hcard : 2 * ((s + 1) * (8 * p)) ≤ n)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual n x y| ^ n ≤
      (2 : ℝ) ^ n * n.factorial *
        (1 / 2 : ℝ) ^ (p * (n - 2 * ((s + 1) * (8 * p)))) *
          2 ^ (2 * ((s + 1) * (8 * p))) := by
  have hpower := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
    residual rows columns realized pivot_ne updates complete n run run_pivots x y hx hy
  rw [Run.finSelectedCore, Matrix.det_reindex_self] at hpower
  have hrows := run.selectedRow_mem_of_completeOn run_complete
  have hcolumns := run.selectedColumn_mem_of_completeOn run_complete
  have hselectedCard := card_selectedIndex_eq_prefixLength residual rows columns n run run_pivots
  have hdet := fermionicKernel_sample_det_le_two_pow p s hp
    run.selectedRow run.selectedColumn
    (fun i => (hrows i).1) (fun i => (hrows i).2)
    (fun j => (hcolumns j).1) (fun j => (hcolumns j).2)
    (by simpa [hselectedCard] using hcard)
  rw [← Run.selectedCore_eq_sampleMatrix] at hdet
  exact hpower.trans (by simpa [hselectedCard] using hdet)

/--
Selected-core residual-power decay using a smaller approximation scale outside
an explicit finite set of exact selected columns.
-/
theorem fermionicKernel_gecp_error_pow_le_separatedApprox_except
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (p approximationScale domainScale : ℕ) (hp : 0 < p)
    (realized : ∀ m, ∃ run : Run fermionicKernel, run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ domainScale : ℕ) : ℝ) ≤ ω ∧
        ω ≤ (2 ^ domainScale : ℕ))
      (residual m) (rows m) (columns m))
    (n : ℕ) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (exceptional : Finset run.SelectedIndex)
    (selectedRow_nonneg : ∀ i, 0 ≤ run.selectedRow i)
    (selectedRow_le_one : ∀ i, run.selectedRow i ≤ 1)
    (selectedColumn_lower : ∀ j, j ∉ exceptional →
      -((2 ^ approximationScale : ℕ) : ℝ) ≤ run.selectedColumn j)
    (selectedColumn_upper : ∀ j, j ∉ exceptional →
      run.selectedColumn j ≤ (2 ^ approximationScale : ℕ))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ domainScale : ℕ) : ℝ) ≤ y ∧
      y ≤ (2 ^ domainScale : ℕ)) :
    |residual n x y| ^ n ≤
      ∑ S : Finset run.SelectedIndex,
        if n ≤ 2 * ((approximationScale + 1) * (8 * p)) +
            exceptional.card + S.card then
          n.factorial * (1 / 2 : ℝ) ^ (p * S.card) * 2 ^ (n - S.card)
        else 0 := by
  have hpower := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -((2 ^ domainScale : ℕ) : ℝ) ≤ ω ∧
      ω ≤ (2 ^ domainScale : ℕ))
    residual rows columns realized pivot_ne updates complete n run run_pivots x y hx hy
  rw [Run.finSelectedCore, Matrix.det_reindex_self] at hpower
  have hdet := fermionicKernel_sample_det_le_separatedApprox_except
    p approximationScale hp run.selectedRow run.selectedColumn exceptional
    selectedRow_nonneg selectedRow_le_one selectedColumn_lower selectedColumn_upper
  rw [← Run.selectedCore_eq_sampleMatrix] at hdet
  have hcard := card_selectedIndex_eq_prefixLength residual rows columns n run run_pivots
  exact hpower.trans (by simpa [hcard] using hdet)

/--
Closed GECP residual-power decay at a smaller approximation scale with an
exact exceptional-column rank shift.
-/
theorem fermionicKernel_gecp_error_pow_le_two_pow_except
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (p approximationScale domainScale : ℕ) (hp : 0 < p)
    (realized : ∀ m, ∃ run : Run fermionicKernel, run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ domainScale : ℕ) : ℝ) ≤ ω ∧
        ω ≤ (2 ^ domainScale : ℕ))
      (residual m) (rows m) (columns m))
    (n : ℕ) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (exceptional : Finset run.SelectedIndex)
    (selectedRow_nonneg : ∀ i, 0 ≤ run.selectedRow i)
    (selectedRow_le_one : ∀ i, run.selectedRow i ≤ 1)
    (selectedColumn_lower : ∀ j, j ∉ exceptional →
      -((2 ^ approximationScale : ℕ) : ℝ) ≤ run.selectedColumn j)
    (selectedColumn_upper : ∀ j, j ∉ exceptional →
      run.selectedColumn j ≤ (2 ^ approximationScale : ℕ))
    (hcard : 2 * ((approximationScale + 1) * (8 * p)) +
      exceptional.card ≤ n)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ domainScale : ℕ) : ℝ) ≤ y ∧
      y ≤ (2 ^ domainScale : ℕ)) :
    |residual n x y| ^ n ≤
      (2 : ℝ) ^ n * n.factorial *
        (1 / 2 : ℝ) ^
          (p * (n - (2 * ((approximationScale + 1) * (8 * p)) +
            exceptional.card))) *
          2 ^ (2 * ((approximationScale + 1) * (8 * p)) +
            exceptional.card) := by
  have hpower := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -((2 ^ domainScale : ℕ) : ℝ) ≤ ω ∧
      ω ≤ (2 ^ domainScale : ℕ))
    residual rows columns realized pivot_ne updates complete n run run_pivots x y hx hy
  rw [Run.finSelectedCore, Matrix.det_reindex_self] at hpower
  have hselectedCard := card_selectedIndex_eq_prefixLength
    residual rows columns n run run_pivots
  have hdet := fermionicKernel_sample_det_le_two_pow_except
    p approximationScale hp run.selectedRow run.selectedColumn exceptional
    selectedRow_nonneg selectedRow_le_one selectedColumn_lower selectedColumn_upper
    (by simpa [hselectedCard] using hcard)
  rw [← Run.selectedCore_eq_sampleMatrix] at hdet
  exact hpower.trans (by simpa [hselectedCard] using hdet)

/--
Choosing odd approximation order `2m+1` and a block of length
`32(2m+1)(s+1)` makes the dyadic exponent divisible by the block length,
so the residual-power estimate has an explicit root.
-/
theorem fermionicKernel_gecp_error_le_oddBlock
    (m s : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ j, ∃ run : Run fermionicKernel, run.finalResidual = residual j)
    (pivot_ne : ∀ j, residual j (rows j) (columns j) ≠ 0)
    (updates : ∀ j, residual (j + 1) =
      residualUpdate (residual j) (rows j) (columns j) (pivot_ne j))
    (complete : ∀ j, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual j) (rows j) (columns j))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns
      (32 * (2 * m + 1) * (s + 1)))
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (2 * m + 1) * (s + 1)) x y| ≤
      2 * (32 * (2 * m + 1) * (s + 1) : ℕ) * (1 / 2 : ℝ) ^ m := by
  let p := 2 * m + 1
  let n := 32 * p * (s + 1)
  let r := 2 * ((s + 1) * (8 * p))
  change |residual n x y| ≤ 2 * n * (1 / 2 : ℝ) ^ m
  have hp : 0 < p := by simp [p]
  have hn : n ≠ 0 := by simp [n, p]
  have hr : r = 16 * p * (s + 1) := by
    simp [r]
    ring
  have hr_le : r ≤ n := by
    rw [hr]
    simp [n]
    omega
  have hpower := fermionicKernel_gecp_error_pow_le_two_pow
    residual rows columns p s hp realized pivot_ne updates complete n run
    (by simpa [n, p] using run_pivots) run_complete (by simpa [r] using hr_le)
    x y hx hy
  have hfac : (n.factorial : ℝ) ≤ (n : ℝ) ^ n := by
    exact_mod_cast n.factorial_le_pow
  have hnsub : n - r = r := by
    rw [hr]
    simp [n]
    ring_nf
    omega
  have hexp : p * (n - r) = n * m + r := by
    rw [hnsub, hr]
    simp [n, p]
    ring
  have hbound :
      (2 : ℝ) ^ n * n.factorial * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r ≤
        (2 * n * (1 / 2 : ℝ) ^ m) ^ n := by
    calc
      (2 : ℝ) ^ n * n.factorial * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r ≤
          (2 : ℝ) ^ n * (n : ℝ) ^ n *
            (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r := by gcongr
      _ = (2 * n * (1 / 2 : ℝ) ^ m) ^ n := by
        rw [hexp, pow_add]
        have hcancel : (1 / 2 : ℝ) ^ r * 2 ^ r = 1 := by
          rw [← mul_pow]
          norm_num
        calc
          (2 : ℝ) ^ n * (n : ℝ) ^ n *
                ((1 / 2 : ℝ) ^ (n * m) * (1 / 2 : ℝ) ^ r) * 2 ^ r =
              ((2 : ℝ) ^ n * (n : ℝ) ^ n * (1 / 2 : ℝ) ^ (n * m)) *
                ((1 / 2 : ℝ) ^ r * 2 ^ r) := by ring
          _ = (2 : ℝ) ^ n * (n : ℝ) ^ n * (1 / 2 : ℝ) ^ (n * m) := by
            rw [hcancel, mul_one]
          _ = (2 * n * (1 / 2 : ℝ) ^ m) ^ n := by
            rw [← mul_pow]
            have hhalf : (1 / 2 : ℝ) ^ (n * m) = ((1 / 2 : ℝ) ^ m) ^ n := by
              rw [Nat.mul_comm]
              exact pow_mul _ _ _
            rw [hhalf]
            rw [← mul_pow]
  have hpow : |residual n x y| ^ n ≤
      (2 * n * (1 / 2 : ℝ) ^ m) ^ n := by
    change |residual n x y| ^ n ≤
      (2 : ℝ) ^ n * n.factorial * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r at hpower
    apply hpower.trans
    exact hbound
  have hroot := (pow_le_pow_iff_left₀ (abs_nonneg _)
    (by positivity : 0 ≤ (2 : ℝ) * n * (1 / 2 : ℝ) ^ m) hn).mp hpow
  exact hroot

/-- An explicit arithmetic condition making the odd-block bound a half contraction. -/
theorem oddBlock_factor_le_half (m s : ℕ)
    (hscale : 128 * (2 * m + 1) * (s + 1) ≤ 2 ^ m) :
    2 * (32 * (2 * m + 1) * (s + 1) : ℕ) * (1 / 2 : ℝ) ^ m ≤ 1 / 2 := by
  have hscale_real :
      (128 * (2 * m + 1) * (s + 1) : ℝ) ≤ (2 : ℝ) ^ m := by
    exact_mod_cast hscale
  rw [div_pow]
  norm_num only [one_pow]
  rw [show
    2 * (32 * (2 * m + 1) * (s + 1) : ℕ) * (1 / (2 : ℝ) ^ m) =
      (2 * (32 * (2 * m + 1) * (s + 1) : ℕ)) / (2 : ℝ) ^ m by ring]
  rw [div_le_iff₀ (by positivity : 0 < (2 : ℝ) ^ m)]
  push_cast
  linarith

/-- An exact arithmetic condition making the odd-block factor at most `2⁻q`. -/
theorem oddBlock_factor_le_accuracy (m s q : ℕ)
    (hscale : 64 * (2 * m + 1) * (s + 1) * 2 ^ q ≤ 2 ^ m) :
    2 * (32 * (2 * m + 1) * (s + 1) : ℕ) * (1 / 2 : ℝ) ^ m ≤
      (1 / 2 : ℝ) ^ q := by
  have hscale_real :
      (64 * (2 * m + 1) * (s + 1) * 2 ^ q : ℝ) ≤ (2 : ℝ) ^ m := by
    exact_mod_cast hscale
  rw [div_pow, div_pow]
  norm_num only [one_pow]
  rw [show
    2 * (32 * (2 * m + 1) * (s + 1) : ℕ) * (1 / (2 : ℝ) ^ m) =
      (2 * (32 * (2 * m + 1) * (s + 1) : ℕ)) / (2 : ℝ) ^ m by ring]
  rw [div_le_div_iff₀ (by positivity : 0 < (2 : ℝ) ^ m)
    (by positivity : 0 < (2 : ℝ) ^ q)]
  push_cast
  norm_num at hscale_real ⊢
  nlinarith [hscale_real]

/-- The simple choice `m = 16(s+1)` satisfies the odd-block half-contraction condition. -/
theorem oddBlock_scale_quadratic (s : ℕ) :
    128 * (2 * (16 * (s + 1)) + 1) * (s + 1) ≤ 2 ^ (16 * (s + 1)) := by
  induction s with
  | zero => norm_num
  | succ s ih =>
      calc
        128 * (2 * (16 * (s + 1 + 1)) + 1) * (s + 1 + 1) ≤
            65536 * (128 * (2 * (16 * (s + 1)) + 1) * (s + 1)) := by
          nlinarith
        _ ≤ 65536 * 2 ^ (16 * (s + 1)) := Nat.mul_le_mul_left 65536 ih
        _ = 2 ^ (16 * (s + 1 + 1)) := by
          rw [show 16 * (s + 1 + 1) = 16 * (s + 1) + 16 by omega, pow_add]
          norm_num
          rw [Nat.mul_comm]

private theorem linear_le_512_mul_pow_two (ell : ℕ) :
    33 + 4 * ell ≤ 512 * 2 ^ ell := by
  induction ell with
  | zero => norm_num
  | succ ell ih =>
      calc
        33 + 4 * (ell + 1) ≤ 2 * (33 + 4 * ell) := by omega
        _ ≤ 2 * (512 * 2 ^ ell) := Nat.mul_le_mul_left 2 ih
        _ = 512 * 2 ^ (ell + 1) := by rw [pow_succ]; ring

/-- A ceil-logarithmic choice of `m` satisfies the odd-block half-contraction condition. -/
theorem oddBlock_scale_logarithmic (s : ℕ) :
    let m := 16 + 2 * Nat.clog 2 (s + 1)
    128 * (2 * m + 1) * (s + 1) ≤ 2 ^ m := by
  let ell := Nat.clog 2 (s + 1)
  have hq : s + 1 ≤ 2 ^ ell := by
    exact Nat.le_pow_clog (by norm_num) (s + 1)
  have hlin := linear_le_512_mul_pow_two ell
  dsimp only
  change 128 * (2 * (16 + 2 * ell) + 1) * (s + 1) ≤ 2 ^ (16 + 2 * ell)
  rw [show 2 * (16 + 2 * ell) + 1 = 33 + 4 * ell by omega]
  calc
    128 * (33 + 4 * ell) * (s + 1) ≤
        128 * (33 + 4 * ell) * 2 ^ ell := by gcongr
    _ ≤ 128 * (512 * 2 ^ ell) * 2 ^ ell := by gcongr
    _ = 2 ^ (16 + 2 * ell) := by
      rw [pow_add, show 2 * ell = ell + ell by omega, pow_add]
      norm_num
      ring

/-- The affine-in-accuracy choice of `m` satisfies the exact dyadic target condition. -/
theorem oddBlock_scale_accuracy (s q : ℕ) :
    let m := 16 + 2 * Nat.clog 2 (s + 1) + 2 * q
    64 * (2 * m + 1) * (s + 1) * 2 ^ q ≤ 2 ^ m := by
  let ell := Nat.clog 2 (s + 1)
  have hq : s + 1 ≤ 2 ^ ell := by
    exact Nat.le_pow_clog (by norm_num) (s + 1)
  have hlin : 33 + 4 * (ell + q) ≤ 1024 * 2 ^ (ell + q) := by
    exact (linear_le_512_mul_pow_two (ell + q)).trans (by gcongr; omega)
  dsimp only
  change 64 * (2 * (16 + 2 * ell + 2 * q) + 1) * (s + 1) * 2 ^ q ≤
    2 ^ (16 + 2 * ell + 2 * q)
  rw [show 2 * (16 + 2 * ell + 2 * q) + 1 = 33 + 4 * (ell + q) by omega]
  calc
    64 * (33 + 4 * (ell + q)) * (s + 1) * 2 ^ q ≤
        64 * (33 + 4 * (ell + q)) * 2 ^ ell * 2 ^ q := by gcongr
    _ ≤ 64 * (1024 * 2 ^ (ell + q)) * 2 ^ ell * 2 ^ q := by gcongr
    _ = 2 ^ (16 + 2 * ell + 2 * q) := by
      rw [pow_add]
      rw [show 16 + 2 * ell + 2 * q = 16 + ell + q + ell + q by omega]
      simp only [pow_add]
      norm_num
      ring

/-- The odd block contracts the exact complete-pivot residual by one half under `hscale`. -/
theorem fermionicKernel_gecp_error_le_half_of_oddBlock
    (m s : ℕ) (hscale : 128 * (2 * m + 1) * (s + 1) ≤ 2 ^ m)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ j, ∃ run : Run fermionicKernel, run.finalResidual = residual j)
    (pivot_ne : ∀ j, residual j (rows j) (columns j) ≠ 0)
    (updates : ∀ j, residual (j + 1) =
      residualUpdate (residual j) (rows j) (columns j) (pivot_ne j))
    (complete : ∀ j, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual j) (rows j) (columns j))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns
      (32 * (2 * m + 1) * (s + 1)))
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (2 * m + 1) * (s + 1)) x y| ≤ 1 / 2 := by
  exact (fermionicKernel_gecp_error_le_oddBlock m s residual rows columns
    realized pivot_ne updates complete run run_pivots run_complete x y hx hy).trans
      (oddBlock_factor_le_half m s hscale)

/-- Every realized physical-domain exact run reaches error at most one half by a quadratic block. -/
theorem fermionicKernel_gecp_error_le_half_quadraticBlock
    (s : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ j, ∃ run : Run fermionicKernel, run.finalResidual = residual j)
    (pivot_ne : ∀ j, residual j (rows j) (columns j) ≠ 0)
    (updates : ∀ j, residual (j + 1) =
      residualUpdate (residual j) (rows j) (columns j) (pivot_ne j))
    (complete : ∀ j, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual j) (rows j) (columns j))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns
      (32 * (2 * (16 * (s + 1)) + 1) * (s + 1)))
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (2 * (16 * (s + 1)) + 1) * (s + 1)) x y| ≤ 1 / 2 := by
  exact fermionicKernel_gecp_error_le_half_of_oddBlock
    (16 * (s + 1)) s (oddBlock_scale_quadratic s) residual rows columns
    realized pivot_ne updates complete run run_pivots run_complete x y hx hy

/-- Every realized physical-domain exact run reaches error at most one half by a
ceil-logarithmic-width odd block. -/
theorem fermionicKernel_gecp_error_le_half_logarithmicBlock
    (s : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ j, ∃ run : Run fermionicKernel, run.finalResidual = residual j)
    (pivot_ne : ∀ j, residual j (rows j) (columns j) ≠ 0)
    (updates : ∀ j, residual (j + 1) =
      residualUpdate (residual j) (rows j) (columns j) (pivot_ne j))
    (complete : ∀ j, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual j) (rows j) (columns j))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns
      (32 * (2 * (16 + 2 * Nat.clog 2 (s + 1)) + 1) * (s + 1)))
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (2 * (16 + 2 * Nat.clog 2 (s + 1)) + 1) * (s + 1)) x y| ≤
      1 / 2 := by
  exact fermionicKernel_gecp_error_le_half_of_oddBlock
    (16 + 2 * Nat.clog 2 (s + 1)) s (oddBlock_scale_logarithmic s)
    residual rows columns realized pivot_ne updates complete run run_pivots run_complete
    x y hx hy

/-- Every realized physical-domain exact run reaches dyadic accuracy `2⁻q` by
an explicit block linear in `q` after a ceil-logarithmic scale startup. -/
theorem fermionicKernel_gecp_error_le_dyadicAccuracyBlock
    (s q : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ j, ∃ run : Run fermionicKernel, run.finalResidual = residual j)
    (pivot_ne : ∀ j, residual j (rows j) (columns j) ≠ 0)
    (updates : ∀ j, residual (j + 1) =
      residualUpdate (residual j) (rows j) (columns j) (pivot_ne j))
    (complete : ∀ j, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ))
      (residual j) (rows j) (columns j))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns
      (32 * (2 * (16 + 2 * Nat.clog 2 (s + 1) + 2 * q) + 1) * (s + 1)))
    (run_complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual
      (32 * (2 * (16 + 2 * Nat.clog 2 (s + 1) + 2 * q) + 1) * (s + 1)) x y| ≤
      (1 / 2 : ℝ) ^ q := by
  exact (fermionicKernel_gecp_error_le_oddBlock
    (16 + 2 * Nat.clog 2 (s + 1) + 2 * q) s residual rows columns realized pivot_ne
    updates complete run run_pivots run_complete x y hx hy).trans
      (oddBlock_factor_le_accuracy
        (16 + 2 * Nat.clog 2 (s + 1) + 2 * q) s q (oddBlock_scale_accuracy s q))

end Fermionic
end GECPKernelStructure
