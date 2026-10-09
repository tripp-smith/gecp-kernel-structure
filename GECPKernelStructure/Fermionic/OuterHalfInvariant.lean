import GECPKernelStructure.Fermionic.ThirdPivotLocalization

namespace GECPKernelStructure
namespace Fermionic

open GECP

/--
After the two symmetric cutoff corners, every later exact complete-pivot
residual remains bounded by one quarter on both outer frequency half-bands.
-/
theorem fermionicKernel_outerHalfBound_preserved
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (realized : ∀ n, ∃ run : Run fermionicKernel,
      run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n)) :
    ∀ n t ω, (0 ≤ t ∧ t ≤ 1) → (-Λ ≤ ω ∧ ω ≤ Λ) →
      Λ / 2 ≤ abs ω → abs (residual n t ω) ≤ 1 / 4 := by
  have initial_bound : ∀ t ω, (0 ≤ t ∧ t ≤ 1) →
      (-Λ ≤ ω ∧ ω ≤ Λ) → Λ / 2 ≤ abs ω →
      abs (residual 0 t ω) ≤ 1 / 4 := by
    intro t ω ht hω frequency_outer
    rw [initial]
    have frequency_order : -Λ < Λ := by linarith
    have residual_nonneg :
        0 ≤ symmetricTwoCornerResidual Λ cutoff_pos t ω := by
      unfold symmetricTwoCornerResidual
      exact fermionicKernel_twoCornerResidual_nonneg frequency_order
        ht.1 ht.2 hω.1 hω.2
    rw [abs_of_nonneg residual_nonneg]
    exact
      fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_absFrequency
        cutoff_pos ht.1 ht.2 hω.1 hω.2 frequency_outer
  have coherent : ∀ n,
      PivotCrossProductSignCoherent (residual n) (rows n) (columns n) := by
    intro n
    obtain ⟨run, run_residual⟩ := realized n
    rw [← run_residual]
    exact fermionicKernel_pivotCrossProductSignCoherent run (rows n) (columns n)
  exact stripBound_preserved_of_signCoherentCompletePivot
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
    (fun ω : ℝ => Λ / 2 ≤ abs ω)
    residual rows columns (1 / 4)
    initial_bound pivot_ne updates coherent complete

/--
Every later complete pivot whose magnitude is above one quarter must lie in
the central frequency half-band.
-/
theorem fermionicKernel_laterCompletePivot_frequency_lt_halfCutoff_of_quarter_lt
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (realized : ∀ n, ∃ run : Run fermionicKernel,
      run.finalResidual = residual n)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (n : ℕ) (quarter_lt :
      1 / 4 < abs (residual n (rows n) (columns n))) :
    abs (columns n) < Λ / 2 := by
  rcases complete n with ⟨row_mem, column_mem, _⟩
  by_contra frequency_not_inner
  have frequency_outer : Λ / 2 ≤ abs (columns n) :=
    le_of_not_gt frequency_not_inner
  have outer_bound := fermionicKernel_outerHalfBound_preserved
    cutoff_pos residual rows columns initial realized pivot_ne updates complete
    n (rows n) (columns n) row_mem column_mem frequency_outer
  linarith

end Fermionic
end GECPKernelStructure
