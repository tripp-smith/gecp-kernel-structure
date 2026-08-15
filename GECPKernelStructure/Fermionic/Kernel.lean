import Mathlib.Analysis.SpecialFunctions.Exp

namespace GECPKernelStructure
namespace Fermionic

/-- The dimensionless fermionic Lehmann kernel. -/
noncomputable def fermionicKernel (t ω : ℝ) : ℝ :=
  Real.exp (-t * ω) / (1 + Real.exp (-ω))

theorem fermionicKernel_denominator_pos (ω : ℝ) :
    0 < 1 + Real.exp (-ω) := by positivity

theorem fermionicKernel_denominator_ne (ω : ℝ) :
    1 + Real.exp (-ω) ≠ 0 :=
  ne_of_gt (fermionicKernel_denominator_pos ω)

theorem fermionicKernel_pos (t ω : ℝ) : 0 < fermionicKernel t ω := by
  unfold fermionicKernel
  positivity

end Fermionic
end GECPKernelStructure
