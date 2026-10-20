import GECPKernelStructure.PositiveDefinite.StrictKernel
import GECPKernelStructure.PositiveDefinite.KernelFillDistance
import GECPKernelStructure.Fermionic.CompleteRun

namespace GECPKernelStructure.PositiveDefinite

open Matrix GECP

universe u
variable {X : Type u} {K : Kernel X X ℝ}

/-- A nonzero complete maximum cannot be off diagonal in the strict fresh domain. -/
theorem strict_completePivot_eq {D : X → Prop} (hPSD : KernelPosSemidef K)
    (hstrict : StrictKernelOn K D) (hzero : ∀ x, ¬D x → ∀ y, K x y = 0)
    {p q : X} (hp : K p q ≠ 0)
    (hmax : CompletePivotOn (fun _ => True) (fun _ => True) K p q) : p = q := by
  have hpD : D p := by
    by_contra hn
    exact hp (hzero p hn q)
  have hqD : D q := by
    by_contra hn
    apply hp
    rw [← hPSD.symmetric p q, hzero q hn p]
  by_contra hne
  have hs := hstrict.abs_sq_lt hpD hqD hne
  have hpp := (le_abs_self (K p p)).trans (hmax.2.2 p p trivial trivial)
  have hqq := (le_abs_self (K q q)).trans (hmax.2.2 q q trivial trivial)
  have hprod := mul_le_mul hpp hqq (hPSD.diag_nonneg q) (abs_nonneg (K p q))
  nlinarith

theorem zeroOutside_residualUpdate {D : X → Prop}
    (hzero : ∀ x, ¬D x → ∀ y, K x y = 0) {p : X} (hp : K p p ≠ 0) :
    ∀ x, ¬(D x ∧ x ≠ p) → ∀ y, GECP.residualUpdate K p p hp x y = 0 := by
  intro x hx y
  by_cases hxD : D x
  · have heq : x = p := by
      by_contra hne
      exact hx ⟨hxD, hne⟩
    subst x
    exact GECP.residualUpdate_selected_row K p p hp y
  · exact GECP.residualUpdate_preserves_zero_row K p p hp (hzero x hxD) y

/-- Fresh strictness and diagonal necessity are propagated through the actual run. -/
theorem completeRun_diagonal_and_strict (run : Run K) {D : X → Prop}
    (hPSD : KernelPosSemidef K) (hstrict : StrictKernelOn K D)
    (hzero : ∀ x, ¬D x → ∀ y, K x y = 0)
    (hcomplete : run.CompleteOn (fun _ => True) (fun _ => True)) :
    DiagonalRun run ∧ StrictKernelOn run.finalResidual (fun x => D x ∧ x ∉ run.rows) := by
  induction run generalizing D with
  | nil K =>
    exact ⟨trivial, hstrict.mono (fun _ hx => hx.1)⟩
  | @step K p q hp tail ih =>
    obtain ⟨hmax, htail⟩ := hcomplete
    obtain rfl := strict_completePivot_eq hPSD hstrict hzero hp hmax
    have hpD : D p := by
      by_contra hn
      exact hp (hzero p hn p)
    obtain ⟨hdiag, hfresh⟩ := ih (hPSD.residualUpdate p hp)
      (hstrict.residualUpdate hpD hp) (zeroOutside_residualUpdate hzero hp) htail
    refine ⟨⟨rfl, hdiag⟩, hfresh.mono ?_⟩
    intro x hx
    change D x ∧ x ∉ p :: tail.rows at hx
    have hn : x ≠ p ∧ x ∉ tail.rows := by simpa only [List.mem_cons, not_or] using hx.2
    exact ⟨⟨hx.1, hn.1⟩, hn.2⟩

/-- Every successful complete pivot is diagonal for an SPD original kernel.
No diagonal-pivot preference or residual strictness is assumed. -/
theorem completeRun_diagonal_of_posDef (run : Run K) (hK : Matrix.PosDef K)
    (hcomplete : run.CompleteOn (fun _ => True) (fun _ => True)) : DiagonalRun run :=
  (completeRun_diagonal_and_strict run (KernelPosSemidef.of_matrix hK.posSemidef)
    (strictKernelOn_of_matrix hK (fun _ => True)) (by simp) hcomplete).1

theorem posDef_completeRun_strictOn_fresh (run : Run K) (hK : Matrix.PosDef K)
    (hcomplete : run.CompleteOn (fun _ => True) (fun _ => True)) :
    StrictKernelOn run.finalResidual (fun x => x ∉ run.rows) := by
  have h := (completeRun_diagonal_and_strict run (KernelPosSemidef.of_matrix hK.posSemidef)
    (strictKernelOn_of_matrix hK (fun _ => True)) (by simp) hcomplete).2
  exact h.mono (fun _ hx => ⟨trivial, hx⟩)

/-- The independent diagonal-maximization rule on actual successive residuals. -/
def DiagonalMaxRun : {K : Kernel X X ℝ} → Run K → Prop
  | _, .nil _ => True
  | K, .step p q _ tail => p = q ∧ (∀ x, K x x ≤ K p p) ∧ DiagonalMaxRun tail

theorem diagonalCompleteRun_to_pivotedCholesky (run : Run K)
    (hPSD : KernelPosSemidef K) (hdiag : DiagonalRun run)
    (hcomplete : run.CompleteOn (fun _ => True) (fun _ => True)) : DiagonalMaxRun run := by
  induction run with
  | nil => trivial
  | @step K p q hp tail ih =>
    obtain ⟨rfl, hdiagtail⟩ := hdiag
    obtain ⟨hmax, htail⟩ := hcomplete
    refine ⟨rfl, ?_, ih (hPSD.residualUpdate p hp) hdiagtail htail⟩
    intro x
    simpa only [abs_of_nonneg (hPSD.diag_nonneg x), abs_of_nonneg (hPSD.diag_nonneg p)]
      using hmax.2.2 x x trivial trivial

theorem pivotedCholesky_to_completeRun (run : Run K)
    (hPSD : KernelPosSemidef K) (htrace : DiagonalMaxRun run) :
    run.CompleteOn (fun _ => True) (fun _ => True) := by
  induction run with
  | nil => trivial
  | @step K p q hp tail ih =>
    obtain ⟨rfl, hmax, htail⟩ := htrace
    refine ⟨⟨trivial, trivial, ?_⟩, ih (hPSD.residualUpdate p hp) htail⟩
    intro x y _ _
    rw [abs_of_nonneg (hPSD.diag_nonneg p)]
    exact (hPSD.diagonal_bound_iff _).mp hmax x y

/-- Actual continuous-domain GECP and pivoted Cholesky have the same finite runs. -/
theorem posDef_completeRun_iff_pivotedCholesky (run : Run K) (hK : Matrix.PosDef K) :
    run.CompleteOn (fun _ => True) (fun _ => True) ↔ DiagonalMaxRun run := by
  constructor
  · intro hc
    exact diagonalCompleteRun_to_pivotedCholesky run (KernelPosSemidef.of_matrix hK.posSemidef)
      (completeRun_diagonal_of_posDef run hK hc) hc
  · exact pivotedCholesky_to_completeRun run (KernelPosSemidef.of_matrix hK.posSemidef)

theorem posDef_completeRun_residual_posSemidef (run : Run K) (hK : Matrix.PosDef K)
    (hc : run.CompleteOn (fun _ => True) (fun _ => True)) :
    KernelPosSemidef run.finalResidual :=
  (diagonalRun_residual_and_removed_posSemidef run (completeRun_diagonal_of_posDef run hK hc)
    (KernelPosSemidef.of_matrix hK.posSemidef)).1

theorem posDef_completeRun_residual_abs_le_diag (run : Run K) (hK : Matrix.PosDef K)
    (hc : run.CompleteOn (fun _ => True) (fun _ => True)) (x y : X) :
    |run.finalResidual x y| ^ 2 ≤ run.finalResidual x x * run.finalResidual y y :=
  (posDef_completeRun_residual_posSemidef run hK hc).abs_sq_le x y

theorem posDef_completeRun_exists_diagonal_completePivot [TopologicalSpace X]
    [CompactSpace X] [Nonempty X] (run : Run K) (hK : Matrix.PosDef K)
    (hc : run.CompleteOn (fun _ => True) (fun _ => True))
    (hcont : Continuous (fun z : X × X => K z.1 z.2)) :
    ∃ p, CompletePivotOn (fun _ => True) (fun _ => True) run.finalResidual p p := by
  obtain ⟨p, hp⟩ := diagonalRun_exists_diagonal_completePivot run
    (completeRun_diagonal_of_posDef run hK hc) (KernelPosSemidef.of_matrix hK.posSemidef) hcont
  refine ⟨p, trivial, trivial, ?_⟩
  intro x y _ _
  exact (hp x y).trans (le_abs_self _)

theorem posDef_completeRun_error_le_original_fillDistance [PseudoMetricSpace X]
    (run : Run K) (hK : Matrix.PosDef K)
    (hc : run.CompleteOn (fun _ => True) (fun _ => True)) {L : NNReal} {h : ℝ}
    (hL : ∀ y, LipschitzWith L (fun x => K x y))
    (hfill : ∀ x, ∃ s ∈ run.rows, dist x s ≤ h) :
    ∀ x y, |run.finalResidual x y| ≤ 2 * (L : ℝ) * h :=
  diagonalRun_error_le_original_fillDistance run (completeRun_diagonal_of_posDef run hK hc)
    (KernelPosSemidef.of_matrix hK.posSemidef) hL hfill

end GECPKernelStructure.PositiveDefinite
