import GECPKernelStructure.Fermionic.ExponentialTotalPositivity
import GECPKernelStructure.Matrix.AndreiefSymmetrization

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/--
On an ordered intermediate sample, the Andréief-symmetrized exponential /
truncated-power integrand has the exponential kernel's prescribed weak sign.
-/
theorem eighthOrderComposition_symmetrized_nonneg
    {n : ℕ} (rows columns sample : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (sample_mono : StrictMono sample) :
    0 ≤ expKernelSignature n *
      ∑ σ : Equiv.Perm (Fin n), Matrix.compositionSampleIntegrand
        (fun i v => expKernel (rows i) v)
        (fun v j => eighthOrderTruncatedPowerKernel v (columns j))
        (sample ∘ σ) := by
  rw [Matrix.sum_compositionSampleIntegrand_comp_perm]
  have exponential_sign := expKernel_strictSignRegular n rows sample
    rows_mono sample_mono
  have truncated_nonneg := eighthOrderTruncatedPowerKernel_totallyNonnegative n
    sample columns sample_mono columns_mono
  change 0 ≤ expKernelSignature n *
    ((GECP.minorMatrix expKernel rows sample).det *
      (GECP.minorMatrix eighthOrderTruncatedPowerKernel sample columns).det)
  rw [← mul_assoc]
  exact mul_nonneg exponential_sign.le truncated_nonneg

end Fermionic
end GECPKernelStructure
