import Mathlib.Topology.MetricSpace.Lipschitz

namespace GECPKernelStructure
namespace Discretization

universe u

/-- A Lipschitz function's continuous supremum is controlled by a grid maximum
plus the grid covering radius. -/
theorem grid_sup_le_max_add_lipschitz {X : Type u} [PseudoMetricSpace X]
    {objective : X → ℝ} {grid : Set X} {L : NNReal} {radius gridMax : ℝ}
    (lipschitz : LipschitzWith L objective)
    (gridBound : ∀ point ∈ grid, objective point ≤ gridMax)
    (cover : ∀ x, ∃ point ∈ grid, dist x point ≤ radius) :
    ∀ x, objective x ≤ gridMax + L * radius := by
  intro x
  obtain ⟨point, point_mem, close⟩ := cover x
  exact (lipschitz.le_add_mul x point).trans <|
    add_le_add (gridBound point point_mem)
      (mul_le_mul_of_nonneg_left close L.coe_nonneg)

end Discretization
end GECPKernelStructure
