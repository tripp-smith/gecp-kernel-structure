import GECPKernelStructure.Fermionic.OuterHalfInvariant
import GECPKernelStructure.Fermionic.GECPDeterminantDecay

namespace GECPKernelStructure
namespace Fermionic

open GECP

/-- Prepend the two prescribed symmetric cutoff-corner pivots to a continuation. -/
noncomputable def symmetricTwoCornerRun {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos)) :
    Run fermionicKernel :=
  .step 0 Λ (fermionicKernel_pos 0 Λ).ne'
    (.step 1 (-Λ) (fermionicKernel_firstPivot_reflected_ne cutoff_pos) continuation)

/-- The two recursive selected indices occupied by the prescribed cutoff corners. -/
noncomputable def symmetricTwoCornerExceptional {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos)) :
    Finset (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex := by
  classical
  exact {Sum.inl 0, Sum.inr (Sum.inl 0)}

/-- The prescribed corner exception set contains exactly two selected indices. -/
theorem symmetricTwoCornerExceptional_card {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos)) :
    (symmetricTwoCornerExceptional cutoff_pos continuation).card = 2 := by
  classical
  unfold symmetricTwoCornerExceptional
  apply Finset.card_pair
  intro h
  cases h

/-- Every row of the full two-corner run lies in the physical time interval. -/
theorem symmetricTwoCornerRun_selectedRow_mem {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos))
    {columnDomain : ℝ → Prop}
    (complete : continuation.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      columnDomain) :
    ∀ i, 0 ≤ (symmetricTwoCornerRun cutoff_pos continuation).selectedRow i ∧
      (symmetricTwoCornerRun cutoff_pos continuation).selectedRow i ≤ 1 := by
  intro i
  change Fin 1 ⊕ (Fin 1 ⊕ continuation.SelectedIndex) at i
  rcases i with first | rest
  · change 0 ≤ (0 : ℝ) ∧ (0 : ℝ) ≤ 1
    norm_num
  · rcases rest with second | tail
    · change 0 ≤ (1 : ℝ) ∧ (1 : ℝ) ≤ 1
      norm_num
    · change 0 ≤ continuation.selectedRow tail ∧ continuation.selectedRow tail ≤ 1
      exact continuation.selectedRow_mem_of_completeOn complete tail

/--
Every selected column outside the two prescribed corner indices belongs to
the continuation's column domain.
-/
theorem symmetricTwoCornerRun_selectedColumn_mem_of_not_exceptional
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos))
    {rowDomain : ℝ → Prop} {columnDomain : ℝ → Prop}
    (complete : continuation.CompleteOn rowDomain columnDomain) :
    ∀ i, i ∉ symmetricTwoCornerExceptional cutoff_pos continuation →
      columnDomain
        ((symmetricTwoCornerRun cutoff_pos continuation).selectedColumn i) := by
  intro i hi
  change Fin 1 ⊕ (Fin 1 ⊕ continuation.SelectedIndex) at i
  rcases i with first | rest
  · exfalso
    apply hi
    have hfirst : first = 0 := Subsingleton.elim _ _
    rw [hfirst]
    unfold symmetricTwoCornerExceptional
    exact Finset.mem_insert_self _ _
  · rcases rest with second | tail
    · exfalso
      apply hi
      have hsecond : second = 0 := Subsingleton.elim _ _
      rw [hsecond]
      unfold symmetricTwoCornerExceptional
      exact Finset.mem_insert_of_mem (Finset.mem_singleton_self _)
    · change columnDomain (continuation.selectedColumn tail)
      exact continuation.selectedColumn_mem_of_completeOn complete tail

/--
The selected-core determinant of a two-corner run with a smaller-band complete
continuation pays exactly two exceptional coordinates.
-/
theorem fermionicKernel_symmetricTwoCornerRun_sample_det_le_two_pow
    {Λ : ℝ} (cutoff_pos : 0 < Λ) (p s : ℕ) (hp : 0 < p)
    (continuation : Run (symmetricTwoCornerResidual Λ cutoff_pos))
    (complete : continuation.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧ ω ≤ (2 ^ s : ℕ)))
    (hcard : 2 * ((s + 1) * (8 * p)) + 2 ≤
      Fintype.card (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex) :
    |(symmetricTwoCornerRun cutoff_pos continuation).selectedCore.det| ≤
      (2 : ℝ) ^ Fintype.card
          (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex *
        (Fintype.card
          (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex).factorial *
        (1 / 2 : ℝ) ^
          (p * (Fintype.card
            (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex -
              (2 * ((s + 1) * (8 * p)) + 2))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by
  let run := symmetricTwoCornerRun cutoff_pos continuation
  let exceptional := symmetricTwoCornerExceptional cutoff_pos continuation
  have hrows := symmetricTwoCornerRun_selectedRow_mem
    cutoff_pos continuation complete
  have hcolumns : ∀ j, j ∉ exceptional →
      -((2 ^ s : ℕ) : ℝ) ≤ run.selectedColumn j ∧
        run.selectedColumn j ≤ (2 ^ s : ℕ) := by
    exact symmetricTwoCornerRun_selectedColumn_mem_of_not_exceptional
      cutoff_pos continuation complete
  have hdet := fermionicKernel_sample_det_le_two_pow_except
    p s hp run.selectedRow run.selectedColumn exceptional
    (fun i => (hrows i).1) (fun i => (hrows i).2)
    (fun j hj => (hcolumns j hj).1) (fun j hj => (hcolumns j hj).2)
    (by simpa [run, exceptional, symmetricTwoCornerExceptional_card] using hcard)
  rw [← Run.selectedCore_eq_sampleMatrix] at hdet
  simpa [run, exceptional, symmetricTwoCornerExceptional_card] using hdet

end Fermionic
end GECPKernelStructure
