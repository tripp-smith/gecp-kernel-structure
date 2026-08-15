import GECPKernelStructure.Fermionic.Centering

namespace GECPKernelStructure
namespace Fermionic

/-- The fermionic kernel is jointly continuous on all real coordinates. -/
theorem fermionicKernel_continuous :
    Continuous (fun point : ℝ × ℝ ↦ fermionicKernel point.1 point.2) := by
  unfold fermionicKernel
  apply Continuous.div
  · fun_prop
  · fun_prop
  · intro point
    exact fermionicKernel_denominator_ne point.2

end Fermionic
end GECPKernelStructure
