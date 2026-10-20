import GECPKernelStructure.PositiveDefinite.KernelGeometry
import Mathlib.Topology.MetricSpace.Lipschitz

namespace GECPKernelStructure.PositiveDefinite

open GECP

universe u

/-- The original kernel's coordinate regularity bounds its squared feature distance. -/
theorem original_energy_le_lipschitz {X : Type u} [PseudoMetricSpace X]
    {K : Kernel X X ℝ} (hK : KernelPosSemidef K) {L : NNReal}
    (hL : ∀ y, LipschitzWith L (fun x => K x y)) (x s : X) :
    K x x - 2 * K x s + K s s ≤ 2 * (L : ℝ) * dist x s := by
  have hx := (hL x).dist_le_mul x s
  have hs := (hL s).dist_le_mul s x
  rw [Real.dist_eq, hK.symmetric x s] at hx
  rw [Real.dist_eq, dist_comm s x] at hs
  have hx' := (le_abs_self (K x x - K x s)).trans hx
  have hs' := (le_abs_self (K s s - K x s)).trans hs
  linarith

/-- Rank-independent control from Lipschitz regularity of the original kernel.

No Lipschitz condition on the residual is assumed, and pivots need not be maximal.
The selected-node cover is the geometric input; the PSD and order invariants are
derived from the actual finite diagonal run.
-/
theorem diagonalRun_error_le_original_fillDistance {X : Type u} [PseudoMetricSpace X]
    {K : Kernel X X ℝ} (run : Run K) (hdiag : DiagonalRun run)
    (hK : KernelPosSemidef K) {L : NNReal} {h : ℝ}
    (hL : ∀ y, LipschitzWith L (fun x => K x y))
    (hfill : ∀ x, ∃ s ∈ run.rows, dist x s ≤ h) :
    ∀ x y, |run.finalResidual x y| ≤ 2 * (L : ℝ) * h := by
  have hR := (diagonalRun_residual_and_removed_posSemidef run hdiag hK).1
  apply (hR.diagonal_bound_iff _).mp
  intro x
  obtain ⟨s, hs, hdist⟩ := hfill x
  exact (diagonalRun_power_le_original_energy run hdiag hK hs x).trans
    ((original_energy_le_lipschitz hK hL x s).trans
      (mul_le_mul_of_nonneg_left hdist (by positivity)))

end GECPKernelStructure.PositiveDefinite
