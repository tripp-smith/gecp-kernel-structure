import GECPKernelStructure.Fermionic.Kernel

namespace GECPKernelStructure
namespace Fermionic

/-- Reflection across `(t, ω) ↦ (1-t, -ω)`. -/
theorem fermionicKernel_reflection (t ω : ℝ) :
    fermionicKernel (1 - t) (-ω) = fermionicKernel t ω := by
  unfold fermionicKernel
  have hexp : Real.exp ((1 - t) * ω) = Real.exp (-t * ω) * Real.exp ω := by
    rw [show (1 - t) * ω = -t * ω + ω by ring, Real.exp_add]
  rw [show -(1 - t) * -ω = (1 - t) * ω by ring,
    show - -ω = ω by ring, hexp, Real.exp_neg]
  field_simp [ne_of_gt (Real.exp_pos ω)]
  ring

end Fermionic
end GECPKernelStructure
