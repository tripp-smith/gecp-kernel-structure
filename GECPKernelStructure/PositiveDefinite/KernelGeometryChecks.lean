import GECPKernelStructure.PositiveDefinite.KernelFillDistance

namespace GECPKernelStructure.PositiveDefinite

open GECP

/-- A nontrivial continuous feature-Gram kernel on an infinite domain. -/
def linearFeatureKernel (x y : ℝ) : ℝ := 1 + x * y

theorem linearFeatureKernel_posSemidef : KernelPosSemidef linearFeatureKernel := by
  intro n nodes
  have h1 := Matrix.posSemidef_vecMulVec_self_star (fun _ : Fin n => (1 : ℝ))
  have h2 := Matrix.posSemidef_vecMulVec_self_star nodes
  have heq : (fun i j => linearFeatureKernel (nodes i) (nodes j)) =
      (Matrix.vecMulVec (fun _ : Fin n => (1 : ℝ)) (star (fun _ => (1 : ℝ))) +
        Matrix.vecMulVec nodes (star nodes) : Matrix (Fin n) (Fin n) ℝ) := by
    ext i j
    change 1 + nodes i * nodes j = 1 * star (1 : ℝ) + nodes i * star (nodes j)
    simp
  rw [heq]
  exact h1.add h2

noncomputable def linearFeatureRun : Run linearFeatureKernel :=
  .step 0 0 (by norm_num [linearFeatureKernel]) (.nil _)

theorem linearFeatureRun_diagonal : DiagonalRun linearFeatureRun := ⟨rfl, trivial⟩

/-- The abstract residual theorem computes the actual feature remainder, not a supplied residual. -/
theorem linearFeatureRun_residual (x y : ℝ) :
    linearFeatureRun.finalResidual x y = x * y := by
  simp [linearFeatureRun, Run.finalResidual, GECP.residualUpdate, linearFeatureKernel]

theorem linearFeatureRun_power (x : ℝ) :
    linearFeatureRun.finalResidual x x ≤
      linearFeatureKernel x x - 2 * linearFeatureKernel x 0 + linearFeatureKernel 0 0 := by
  apply diagonalRun_power_le_original_energy linearFeatureRun linearFeatureRun_diagonal
    linearFeatureKernel_posSemidef
  simp [linearFeatureRun, Run.rows]

end GECPKernelStructure.PositiveDefinite
