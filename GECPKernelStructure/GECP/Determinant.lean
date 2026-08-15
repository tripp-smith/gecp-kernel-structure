import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u

variable {𝕜 : Type u} [Field 𝕜]

/-- LDU data produced by exact elimination of a selected finite core. -/
structure CoreDecomposition (n : ℕ) where
  core : Matrix (Fin n) (Fin n) 𝕜
  lower : Matrix (Fin n) (Fin n) 𝕜
  upper : Matrix (Fin n) (Fin n) 𝕜
  pivots : Fin n → 𝕜
  factorization : core = lower * diagonal pivots * upper
  det_lower : lower.det = 1
  det_upper : upper.det = 1

/-- A successful core decomposition has no zero pivot. -/
structure SuccessfulCoreDecomposition (n : ℕ) extends CoreDecomposition (𝕜 := 𝕜) n where
  pivot_ne : ∀ i, pivots i ≠ 0

/-- The ordered selected-core determinant is the product of elimination pivots. -/
theorem gecp_core_det_eq_prod_pivots {n : ℕ} (decomposition : CoreDecomposition (𝕜 := 𝕜) n) :
    decomposition.core.det = ∏ i, decomposition.pivots i := by
  rw [decomposition.factorization, det_mul, det_mul, det_diagonal,
    decomposition.det_lower, decomposition.det_upper]
  simp

/-- A selected core with a successful LDU decomposition is nonsingular. -/
theorem gecp_core_nonsingular {n : ℕ}
    (decomposition : SuccessfulCoreDecomposition (𝕜 := 𝕜) n) :
    decomposition.core.det ≠ 0 := by
  rw [gecp_core_det_eq_prod_pivots decomposition.toCoreDecomposition]
  exact Finset.prod_ne_zero_iff.mpr fun i _ ↦ decomposition.pivot_ne i

end GECP
end GECPKernelStructure
