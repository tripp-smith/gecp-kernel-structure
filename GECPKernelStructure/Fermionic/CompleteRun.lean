import GECPKernelStructure.Fermionic.ExponentialTotalPositivity

namespace GECPKernelStructure
universe u v

namespace GECP

/-- A pivot belongs to the chosen domain and maximizes the residual there. -/
def CompletePivotOn {α : Type u} {β : Type v} (rowDomain : α → Prop)
    (columnDomain : β → Prop) (R : Kernel α β ℝ) (row : α) (column : β) : Prop :=
  rowDomain row ∧ columnDomain column ∧
    ∀ x y, rowDomain x → columnDomain y → abs (R x y) ≤ abs (R row column)

namespace Run

/-- Every pivot in a successful finite run is complete on the chosen domain. -/
def CompleteOn {α : Type u} {β : Type v} (rowDomain : α → Prop)
    (columnDomain : β → Prop) : {K : Kernel α β ℝ} → Run K → Prop
  | _, .nil _ => True
  | K, .step row column _ tail =>
      CompletePivotOn rowDomain columnDomain K row column ∧
        CompleteOn rowDomain columnDomain tail

end Run

end GECP

namespace Fermionic

open GECP

private theorem abs_sub_le_of_sameSign_of_four_bounds {a b B : ℝ}
    (same_sign : 0 ≤ a * b) (ha : abs a ≤ B * B) (hb : abs b ≤ B * B) :
    abs (a - b) ≤ B * B := by
  rw [abs_le]
  rcases (mul_nonneg_iff.mp same_sign) with ⟨ha0, hb0⟩ | ⟨ha0, hb0⟩
  · rw [abs_of_nonneg ha0] at ha
    rw [abs_of_nonneg hb0] at hb
    constructor <;> linarith
  · rw [abs_of_nonpos ha0] at ha
    rw [abs_of_nonpos hb0] at hb
    constructor <;> linarith

/-- Sign coherence gives a nonexpansive update using only bounds on the four involved points. -/
theorem residualUpdate_le_of_signCoherentOn
    {α : Type u} {β : Type v} {R : Kernel α β ℝ}
    {row x : α} {column y : β}
    (coherent : PivotCrossProductSignCoherent R row column) {B : ℝ}
    (bound_xy : abs (R x y) ≤ B) (bound_pivot : abs (R row column) = B)
    (bound_x_column : abs (R x column) ≤ B)
    (bound_row_y : abs (R row y) ≤ B)
    (pivot_ne : R row column ≠ 0) :
    abs (residualUpdate R row column pivot_ne x y) ≤ B := by
  have B_nonneg : 0 ≤ B := by rw [← bound_pivot]; exact abs_nonneg _
  have first_product : abs (R x y * R row column) ≤ B * B := by
    rw [abs_mul, bound_pivot]
    exact mul_le_mul_of_nonneg_right bound_xy B_nonneg
  have second_product : abs (R x column * R row y) ≤ B * B := by
    rw [abs_mul]
    exact mul_le_mul bound_x_column bound_row_y (abs_nonneg _) B_nonneg
  have difference := abs_sub_le_of_sameSign_of_four_bounds
    (coherent x y) first_product second_product
  unfold residualUpdate
  rw [show R x y - R x column * R row y / R row column =
      (R x y * R row column - R x column * R row y) / R row column by
    field_simp]
  rw [abs_div, bound_pivot]
  have B_pos : 0 < B := by
    rw [← bound_pivot]
    exact abs_pos.mpr pivot_ne
  exact (div_le_iff₀ B_pos).2 difference

/--
For a strictly sign-regular original kernel, an exact complete-pivot residual
sequence cannot exceed its initial domain bound.  `realized` is the essential
link to the original kernel: it supplies the finite prefix whose bordered
minors control the sign of each current residual.
-/
theorem strictSignRegular_gecp_error_nonincreasing
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (residual : ℕ → Kernel α β ℝ) (rows : ℕ → α) (columns : ℕ → β)
    (initialBound : ℝ)
    (initial_bound : ∀ x y, rowDomain x → columnDomain y →
      abs (residual 0 x y) ≤ initialBound)
    (realized : ∀ n, ∃ run : Run K, run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn rowDomain columnDomain (residual n) (rows n) (columns n)) :
    ∀ n x y, rowDomain x → columnDomain y →
      abs (residual n x y) ≤ initialBound := by
  intro n
  induction n with
  | zero => exact initial_bound
  | succ n ih =>
      intro x y hx hy
      obtain ⟨run, run_residual⟩ := realized n
      have coherent :
          PivotCrossProductSignCoherent (residual n) (rows n) (columns n) := by
        rw [← run_residual]
        exact strictSignRegular_pivotCrossProductSignCoherent
          regular run (rows n) (columns n)
      obtain ⟨row_mem, column_mem, pivot_max⟩ := complete n
      have step_bound :
          abs (residualUpdate (residual n) (rows n) (columns n) (pivot_ne n) x y) ≤
            abs (residual n (rows n) (columns n)) := by
        apply residualUpdate_le_of_signCoherentOn coherent
        · exact pivot_max x y hx hy
        · rfl
        · exact pivot_max x (columns n) hx column_mem
        · exact pivot_max (rows n) y row_mem hy
      rw [updates n]
      exact step_bound.trans (ih (rows n) (columns n) row_mem column_mem)

/--
Every exact fermionic complete-pivot residual sequence on the cutoff rectangle
is bounded by the initial cutoff-corner pivot, at every finite rank.
-/
theorem fermionicKernel_gecp_error_le_cutoffCorner
    {Λ : ℝ}
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (realized : ∀ n, ∃ run : Run fermionicKernel, run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n)) :
    ∀ n t ω, (0 ≤ t ∧ t ≤ 1) → (-Λ ≤ ω ∧ ω ≤ Λ) →
      abs (residual n t ω) ≤ fermionicKernel 0 Λ := by
  have initial_bound : ∀ t ω, (0 ≤ t ∧ t ≤ 1) → (-Λ ≤ ω ∧ ω ≤ Λ) →
      abs (residual 0 t ω) ≤ fermionicKernel 0 Λ := by
    intro t ω ht hω
    rw [initial, abs_of_pos (fermionicKernel_pos t ω)]
    exact fermionicKernel_le_cutoffCorner ht.1 ht.2 hω.1 hω.2
  exact strictSignRegular_gecp_error_nonincreasing
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
    residual rows columns (fermionicKernel 0 Λ)
    initial_bound realized pivot_ne updates complete

end Fermionic
end GECPKernelStructure
