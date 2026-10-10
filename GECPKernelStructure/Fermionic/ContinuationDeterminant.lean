import GECPKernelStructure.Fermionic.Condensation
import GECPKernelStructure.Fermionic.CentralPrefixRun
import GECPKernelStructure.Fermionic.GeometricMean

namespace GECPKernelStructure

universe u v

namespace GECP

open Fermionic

/-- Every recursive pivot selected by a run has a sign-coherent residual cross. -/
def Run.SignCoherent {α : Type u} {β : Type v} :
    {K : Kernel α β ℝ} → Run K → Prop
  | _, .nil _ => True
  | K, .step row column _ tail =>
      PivotCrossProductSignCoherent K row column ∧ tail.SignCoherent

/-- Every recursive pivot magnitude of a run is bounded by `B`. -/
def Run.PivotsBounded {α : Type u} {β : Type v} (B : ℝ) :
    {K : Kernel α β ℝ} → Run K → Prop
  | _, .nil _ => True
  | K, .step row column _ tail => abs (K row column) ≤ B ∧ tail.PivotsBounded B

namespace Run

/-- Pointwise pivot bounds give the matching factorial-free selected-core bound. -/
theorem abs_selectedCore_det_le_pow_of_pivotsBounded
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) {B : ℝ} (bounded : run.PivotsBounded B) :
    abs run.selectedCore.det ≤ B ^ Fintype.card run.SelectedIndex := by
  induction run with
  | nil K => simp [selectedCore_nil, SelectedIndex]
  | @step K row column pivot_ne tail ih =>
      rcases bounded with ⟨pivot_bound, tail_bound⟩
      rw [selectedCore_step_det, abs_mul]
      have hmul := mul_le_mul pivot_bound (ih tail_bound) (abs_nonneg _) (by
        exact le_trans (abs_nonneg _) pivot_bound)
      rw [show Fintype.card (Run.step row column pivot_ne tail).SelectedIndex =
        Fintype.card tail.SelectedIndex + 1 by
          simp [SelectedIndex, Nat.add_comm]]
      rw [pow_succ]
      simpa [mul_comm] using hmul

/--
Sign-coherent complete pivoting preserves a uniform envelope and bounds every
pivot in a finite run by that same envelope.
-/
theorem pivotsBounded_of_signCoherent_completeOn
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (run : Run K) {B : ℝ}
    (initial_bound : ∀ x y, rowDomain x → columnDomain y → abs (K x y) ≤ B)
    (coherent : run.SignCoherent)
    (complete : run.CompleteOn rowDomain columnDomain) :
    run.PivotsBounded B := by
  induction run with
  | nil K => trivial
  | @step K row column pivot_ne tail ih =>
      rcases coherent with ⟨pivot_coherent, tail_coherent⟩
      rcases complete with ⟨pivot_complete, tail_complete⟩
      obtain ⟨row_mem, column_mem, pivot_max⟩ := pivot_complete
      constructor
      · exact initial_bound row column row_mem column_mem
      · apply ih
        · intro x y hx hy
          have step_bound :
              abs (residualUpdate K row column pivot_ne x y) ≤
                abs (K row column) := by
            apply residualUpdate_le_of_signCoherentOn pivot_coherent
            · exact pivot_max x y hx hy
            · rfl
            · exact pivot_max x column hx column_mem
            · exact pivot_max row y row_mem hy
          exact step_bound.trans (initial_bound row column row_mem column_mem)
        · exact tail_coherent
        · exact tail_complete

/-- A sign-coherent complete run has no factorial in its determinant envelope. -/
theorem abs_selectedCore_det_le_pow_of_signCoherent_completeOn
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (run : Run K) {B : ℝ}
    (initial_bound : ∀ x y, rowDomain x → columnDomain y → abs (K x y) ≤ B)
    (coherent : run.SignCoherent)
    (complete : run.CompleteOn rowDomain columnDomain) :
    abs run.selectedCore.det ≤ B ^ Fintype.card run.SelectedIndex :=
  run.abs_selectedCore_det_le_pow_of_pivotsBounded
    (run.pivotsBounded_of_signCoherent_completeOn rowDomain columnDomain
      initial_bound coherent complete)

@[simp]
theorem signCoherent_transport
    {α : Type u} {β : Type v} {K L : Kernel α β ℝ}
    (kernel_eq : K = L) (run : Run K) :
    (run.transport kernel_eq).SignCoherent ↔ run.SignCoherent := by
  cases kernel_eq
  rfl

end Run
end GECP

namespace Fermionic

open GECP

/-- Every pivot in a continuation of a strictly sign-regular run has a coherent cross. -/
theorem strictSignRegular_continuation_signCoherent
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (initialRun : Run K)
    (continuation : Run initialRun.finalResidual) :
    continuation.SignCoherent := by
  cases continuation with
  | nil => trivial
  | @step R row column pivot_ne tail =>
      constructor
      · exact strictSignRegular_pivotCrossProductSignCoherent
          regular initialRun row column
      · let oneStep : Run initialRun.finalResidual :=
          Run.step row column pivot_ne (Run.nil _)
        let nextRun : Run K := initialRun.append oneStep
        have next_final : nextRun.finalResidual =
            residualUpdate initialRun.finalResidual row column pivot_ne := by
          change (initialRun.append oneStep).finalResidual = _
          rw [Run.finalResidual_append]
          rfl
        let transportedTail : Run nextRun.finalResidual :=
          tail.transport next_final.symm
        have htransported := strictSignRegular_continuation_signCoherent
          regular nextRun transportedTail
        exact (Run.signCoherent_transport next_final.symm tail).mp htransported
termination_by continuation.pivots.length
decreasing_by
  simp_all [Run.pivots]
  exact le_of_eq (congrArg List.length
    (Run.pivots_transport next_final.symm tail))

/-- Every pivot in a complete strictly sign-regular continuation stays inside its initial envelope. -/
theorem strictSignRegular_completeContinuation_pivotsBounded
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (initialRun : Run K) (continuation : Run initialRun.finalResidual)
    {B : ℝ}
    (initial_bound : ∀ x y, rowDomain x → columnDomain y →
      abs (initialRun.finalResidual x y) ≤ B)
    (complete : continuation.CompleteOn rowDomain columnDomain) :
    continuation.PivotsBounded B :=
  continuation.pivotsBounded_of_signCoherent_completeOn
    rowDomain columnDomain initial_bound
    (strictSignRegular_continuation_signCoherent regular initialRun continuation)
    complete

/--
A complete continuation of a strictly sign-regular run has a factorial-free
selected-core determinant envelope.
-/
theorem strictSignRegular_completeContinuation_selectedCore_det_le_pow
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature)
    (rowDomain : α → Prop) (columnDomain : β → Prop)
    (initialRun : Run K) (continuation : Run initialRun.finalResidual)
    {B : ℝ}
    (initial_bound : ∀ x y, rowDomain x → columnDomain y →
      abs (initialRun.finalResidual x y) ≤ B)
    (complete : continuation.CompleteOn rowDomain columnDomain) :
    abs continuation.selectedCore.det ≤
      B ^ Fintype.card continuation.SelectedIndex :=
  continuation.abs_selectedCore_det_le_pow_of_pivotsBounded
    (strictSignRegular_completeContinuation_pivotsBounded regular
      rowDomain columnDomain initialRun continuation initial_bound complete)

/-- Fermionic complete continuations satisfy the same factorial-free envelope. -/
theorem fermionicKernel_completeContinuation_selectedCore_det_le_pow
    (rowDomain columnDomain : ℝ → Prop)
    (initialRun : Run fermionicKernel)
    (continuation : Run initialRun.finalResidual) {B : ℝ}
    (initial_bound : ∀ t ω, rowDomain t → columnDomain ω →
      abs (initialRun.finalResidual t ω) ≤ B)
    (complete : continuation.CompleteOn rowDomain columnDomain) :
    abs continuation.selectedCore.det ≤
      B ^ Fintype.card continuation.SelectedIndex :=
  strictSignRegular_completeContinuation_selectedCore_det_le_pow
    fermionicKernel_strictSignRegular rowDomain columnDomain initialRun continuation
    initial_bound complete

end Fermionic
end GECPKernelStructure
