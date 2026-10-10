import GECPKernelStructure.Fermionic.TwoCornerRun

namespace GECPKernelStructure

namespace GECP
universe u v w

namespace Run

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

/-- Transport a successful dependent run across equality of its initial kernel. -/
noncomputable def transport {K L : Kernel α β 𝕜} (kernel_eq : K = L)
    (run : Run K) : Run L := by
  rw [← kernel_eq]
  exact run

@[simp]
theorem finalResidual_transport {K L : Kernel α β 𝕜} (kernel_eq : K = L)
    (run : Run K) : (transport kernel_eq run).finalResidual = run.finalResidual := by
  cases kernel_eq
  rfl

@[simp]
theorem pivots_transport {K L : Kernel α β 𝕜} (kernel_eq : K = L)
    (run : Run K) : (transport kernel_eq run).pivots = run.pivots := by
  cases kernel_eq
  rfl

@[simp]
theorem card_selectedIndex_transport {K L : Kernel α β 𝕜}
    (kernel_eq : K = L) (run : Run K) :
    Fintype.card (transport kernel_eq run).SelectedIndex =
      Fintype.card run.SelectedIndex := by
  cases kernel_eq
  rfl

/-- Assemble a finite successful run from an indexed exact residual recurrence. -/
noncomputable def ofResidualSequenceFrom
    (residual : ℕ → Kernel α β 𝕜) (rows : ℕ → α) (columns : ℕ → β)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n)) :
    (start length : ℕ) → Run (residual start)
  | start, 0 => .nil (residual start)
  | start, length + 1 =>
      .step (rows start) (columns start) (pivot_ne start)
        (transport (updates start)
          (ofResidualSequenceFrom residual rows columns pivot_ne updates
            (start + 1) length))

/-- The finite run assembled from a recurrence ends at the indexed residual. -/
@[simp]
theorem finalResidual_ofResidualSequenceFrom
    (residual : ℕ → Kernel α β 𝕜) (rows : ℕ → α) (columns : ℕ → β)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (start length : ℕ) :
    (ofResidualSequenceFrom residual rows columns pivot_ne updates start length).finalResidual =
      residual (start + length) := by
  induction length generalizing start with
  | zero => rfl
  | succ length ih =>
      change
        (transport (updates start)
          (ofResidualSequenceFrom residual rows columns pivot_ne updates
            (start + 1) length)).finalResidual =
          residual (start + (length + 1))
      rw [finalResidual_transport, ih]
      congr 1
      omega

/-- The assembled run stores exactly the requested consecutive pivots. -/
@[simp]
theorem pivots_ofResidualSequenceFrom
    (residual : ℕ → Kernel α β 𝕜) (rows : ℕ → α) (columns : ℕ → β)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (start length : ℕ) :
    (ofResidualSequenceFrom residual rows columns pivot_ne updates start length).pivots =
      (List.range' start length).map fun n =>
        residual n (rows n) (columns n) := by
  induction length generalizing start with
  | zero => rfl
  | succ length ih =>
      change
        residual start (rows start) (columns start) ::
            (transport (updates start)
              (ofResidualSequenceFrom residual rows columns pivot_ne updates
                (start + 1) length)).pivots =
          (List.range' start (length + 1)).map fun n =>
            residual n (rows n) (columns n)
      rw [pivots_transport, ih, List.range'_succ]
      rfl

/-- The assembled run has exactly the requested number of selected indices. -/
@[simp]
theorem card_selectedIndex_ofResidualSequenceFrom
    (residual : ℕ → Kernel α β 𝕜) (rows : ℕ → α) (columns : ℕ → β)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (start length : ℕ) :
    Fintype.card
        (ofResidualSequenceFrom residual rows columns pivot_ne updates
          start length).SelectedIndex = length := by
  let run := ofResidualSequenceFrom residual rows columns pivot_ne updates
    start length
  calc
    Fintype.card run.SelectedIndex = Fintype.card (Fin run.pivots.length) :=
      Fintype.card_congr run.selectedIndexEquivFin
    _ = run.pivots.length := Fintype.card_fin _
    _ = length := by
      rw [pivots_ofResidualSequenceFrom]
      simp

section CompleteOn

variable {rowDomain : α → Prop} {columnDomain : β → Prop}

@[simp]
theorem completeOn_transport {K L : Kernel α β ℝ} (kernel_eq : K = L)
    (run : Run K) :
    CompleteOn rowDomain columnDomain (transport kernel_eq run) ↔
      CompleteOn rowDomain columnDomain run := by
  cases kernel_eq
  rfl

/-- Pointwise complete pivots assemble into completeness of the finite run. -/
theorem completeOn_ofResidualSequenceFrom
    (residual : ℕ → Kernel α β ℝ) (rows : ℕ → α) (columns : ℕ → β)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (start length : ℕ)
    (complete : ∀ j, j < length →
      CompletePivotOn rowDomain columnDomain
        (residual (start + j)) (rows (start + j)) (columns (start + j))) :
    CompleteOn rowDomain columnDomain
      (ofResidualSequenceFrom residual rows columns pivot_ne updates start length) := by
  induction length generalizing start with
  | zero => trivial
  | succ length ih =>
      change
        CompletePivotOn rowDomain columnDomain
            (residual start) (rows start) (columns start) ∧
          CompleteOn rowDomain columnDomain
            (transport (updates start)
              (ofResidualSequenceFrom residual rows columns pivot_ne updates
                (start + 1) length))
      constructor
      · simpa using complete 0 (by omega)
      · rw [completeOn_transport]
        apply ih
        intro j hj
        have pointwise := complete (j + 1) (by omega)
        simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using pointwise

end CompleteOn

end Run
end GECP

namespace Fermionic

open GECP

/--
A later full-band complete pivot above one quarter is also complete on the
central frequency half-band.
-/
theorem fermionicKernel_laterCompletePivot_completeOn_centralHalf_of_quarter_lt
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
    CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -(Λ / 2) ≤ ω ∧ ω ≤ Λ / 2)
      (residual n) (rows n) (columns n) := by
  obtain ⟨row_mem, column_mem, pivot_max⟩ := complete n
  have central :=
    fermionicKernel_laterCompletePivot_frequency_lt_halfCutoff_of_quarter_lt
      cutoff_pos residual rows columns initial realized pivot_ne updates complete
      n quarter_lt
  rw [abs_lt] at central
  refine ⟨row_mem, ⟨central.1.le, central.2.le⟩, ?_⟩
  intro t ω time_mem frequency_mem
  apply pivot_max t ω time_mem
  constructor <;> linarith

/-- A finite continuation prefix, transported to the exact two-corner residual. -/
noncomputable def symmetricTwoCornerContinuationPrefix
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (length : ℕ) : Run (symmetricTwoCornerResidual Λ cutoff_pos) :=
  Run.transport initial
    (Run.ofResidualSequenceFrom residual rows columns pivot_ne updates 0 length)

/-- The transported two-corner continuation prefix ends at `residual length`. -/
@[simp]
theorem symmetricTwoCornerContinuationPrefix_finalResidual
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (length : ℕ) :
    (symmetricTwoCornerContinuationPrefix cutoff_pos residual rows columns
      initial pivot_ne updates length).finalResidual = residual length := by
  unfold symmetricTwoCornerContinuationPrefix
  simp

/-- The transported continuation has exactly its requested prefix length. -/
@[simp]
theorem symmetricTwoCornerContinuationPrefix_card
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (length : ℕ) :
    Fintype.card
        (symmetricTwoCornerContinuationPrefix cutoff_pos residual rows columns
          initial pivot_ne updates length).SelectedIndex = length := by
  unfold symmetricTwoCornerContinuationPrefix
  simp

/-- The full prescribed-corner run has the continuation length plus two steps. -/
@[simp]
theorem symmetricTwoCornerPrefixRun_card
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (length : ℕ) :
    Fintype.card
        (symmetricTwoCornerRun cutoff_pos
          (symmetricTwoCornerContinuationPrefix cutoff_pos residual rows columns
            initial pivot_ne updates length)).SelectedIndex = length + 2 := by
  change Fintype.card
      (Fin 1 ⊕ Fin 1 ⊕
        (symmetricTwoCornerContinuationPrefix cutoff_pos residual rows columns
          initial pivot_ne updates length).SelectedIndex) = length + 2
  rw [Fintype.card_sum, Fintype.card_sum]
  simp
  omega

/-- Every indexed continuation residual is realized by an actual fermionic run. -/
theorem symmetricTwoCornerResidualSequence_realized
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n)) :
    ∀ n, ∃ run : Run fermionicKernel, run.finalResidual = residual n := by
  intro n
  let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos residual
    rows columns initial pivot_ne updates n
  refine ⟨symmetricTwoCornerRun cutoff_pos continuation, ?_⟩
  change continuation.finalResidual = residual n
  exact symmetricTwoCornerContinuationPrefix_finalResidual cutoff_pos residual
    rows columns initial pivot_ne updates n

/-- Every finite above-quarter continuation prefix is complete centrally. -/
theorem fermionicKernel_aboveQuarterPrefix_completeOn_centralHalf
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (length : ℕ)
    (above_quarter : ∀ j, j < length →
      1 / 4 < abs (residual j (rows j) (columns j))) :
    Run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun ω : ℝ => -(Λ / 2) ≤ ω ∧ ω ≤ Λ / 2)
      (symmetricTwoCornerContinuationPrefix cutoff_pos residual rows columns
        initial pivot_ne updates length) := by
  unfold symmetricTwoCornerContinuationPrefix
  rw [Run.completeOn_transport]
  apply Run.completeOn_ofResidualSequenceFrom
  intro j hj
  have realized := symmetricTwoCornerResidualSequence_realized cutoff_pos
    residual rows columns initial pivot_ne updates
  simpa using
    fermionicKernel_laterCompletePivot_completeOn_centralHalf_of_quarter_lt
      cutoff_pos residual rows columns initial realized pivot_ne updates complete
      j (above_quarter j hj)

/--
The actual two-corner run followed by an above-quarter trajectory prefix obeys
the smaller-band selected-core determinant bound, with no continuation-domain
or realization hypothesis left to supply.
-/
theorem fermionicKernel_aboveQuarterPrefix_sample_det_le_two_pow
    {Λ : ℝ} (cutoff_pos : 0 < Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = symmetricTwoCornerResidual Λ cutoff_pos)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (length p s : ℕ) (hp : 0 < p)
    (half_cutoff : Λ / 2 = ((2 ^ s : ℕ) : ℝ))
    (above_quarter : ∀ j, j < length →
      1 / 4 < abs (residual j (rows j) (columns j)))
    (long_enough : 2 * ((s + 1) * (8 * p)) ≤ length) :
    let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos
      residual rows columns initial pivot_ne updates length
    let run := symmetricTwoCornerRun cutoff_pos continuation
    |run.selectedCore.det| ≤
      (2 : ℝ) ^ (length + 2) * (length + 2).factorial *
        (1 / 2 : ℝ) ^
          (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by
  dsimp only
  let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos
    residual rows columns initial pivot_ne updates length
  have continuation_complete :
      continuation.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -((2 ^ s : ℕ) : ℝ) ≤ ω ∧
          ω ≤ (2 ^ s : ℕ)) := by
    have central := fermionicKernel_aboveQuarterPrefix_completeOn_centralHalf
      cutoff_pos residual rows columns initial pivot_ne updates complete length
      above_quarter
    simpa [continuation, half_cutoff] using central
  have card_bound : 2 * ((s + 1) * (8 * p)) + 2 ≤
      Fintype.card (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex := by
    rw [show Fintype.card
        (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex = length + 2 by
      simp [continuation]]
    omega
  have determinant_bound :=
    fermionicKernel_symmetricTwoCornerRun_sample_det_le_two_pow cutoff_pos p s hp
      continuation continuation_complete card_bound
  rw [show Fintype.card
      (symmetricTwoCornerRun cutoff_pos continuation).SelectedIndex = length + 2 by
    simp [continuation]] at determinant_bound
  simpa [Nat.add_sub_add_right] using determinant_bound

end Fermionic
end GECPKernelStructure
