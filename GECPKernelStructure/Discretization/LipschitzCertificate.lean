import GECPKernelStructure.Discretization.Grid

namespace GECPKernelStructure
namespace Discretization

/-- A certified upper bound converts a sampled candidate into an
`η`-approximate complete pivot. -/
theorem approxPivot_of_grid_certificate {X : Type*} {objective : X → ℝ}
    {candidate upper η : ℝ} (eta_nonneg : 0 ≤ η)
    (global_upper : ∀ x, objective x ≤ upper)
    (scaled_upper_le_candidate : η * upper ≤ candidate) :
    ∀ x, η * objective x ≤ candidate := by
  intro x
  exact (mul_le_mul_of_nonneg_left (global_upper x) eta_nonneg).trans
    scaled_upper_le_candidate

end Discretization
end GECPKernelStructure
