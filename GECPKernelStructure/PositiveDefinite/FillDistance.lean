import Mathlib.Topology.MetricSpace.Lipschitz

namespace GECPKernelStructure
namespace PositiveDefinite

universe u

/-- Conditional fill-distance control for a Lipschitz residual diagonal. -/
theorem posDef_gecp_error_le_fillDistance {X : Type u} [PseudoMetricSpace X]
    {diagonal : X → ℝ} {selected : Set X} {L : NNReal} {h : ℝ}
    (lipschitz : LipschitzWith L diagonal)
    (interpolates : ∀ point ∈ selected, diagonal point = 0)
    (fillDistance : ∀ x, ∃ point ∈ selected, dist x point ≤ h) :
    ∀ x, |diagonal x| ≤ L * h := by
  intro x
  obtain ⟨point, point_mem, close⟩ := fillDistance x
  have hdist := lipschitz.dist_le_mul x point
  rw [Real.dist_eq, interpolates point point_mem, sub_zero] at hdist
  exact hdist.trans (mul_le_mul_of_nonneg_left close L.coe_nonneg)

end PositiveDefinite
end GECPKernelStructure
