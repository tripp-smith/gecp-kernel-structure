import GECPKernelStructure.PositiveDefinite.ResidualPSD

namespace GECPKernelStructure
namespace PositiveDefinite

open Matrix

universe u

/-- A diagonal maximizer of a finite real PSD residual is a complete-pivot maximizer.

Together with a diagonal-preferring deterministic tie rule, this is the exact
one-step reduction used by canonical pivoted Cholesky.
-/
theorem gecp_eq_pivotedCholesky_of_posDef {n : Type u} [Fintype n]
    (M : Matrix n n ℝ) (hM : M.PosSemidef) (pivot : n) :
    (∀ i, M i i ≤ M pivot pivot) ↔ ∀ i j, |M i j| ≤ M pivot pivot :=
  posDef_gecp_residual_sup_eq_diag_sup M hM (M pivot pivot)

end PositiveDefinite
end GECPKernelStructure
