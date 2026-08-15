import GECPKernelStructure.Fermionic.Symmetry
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

namespace GECPKernelStructure
namespace Fermionic

/-- The centered hyperbolic-cosine form of the fermionic kernel. -/
noncomputable def centeredKernel (t ω : ℝ) : ℝ :=
  Real.exp (-(t - (1 / 2 : ℝ)) * ω) / (2 * Real.cosh (ω / 2))

/-- The direct and centered closed forms agree exactly. -/
theorem fermionicKernel_centered (t ω : ℝ) :
    fermionicKernel t ω = centeredKernel t ω := by
  have hnum : Real.exp (-(t - (1 / 2 : ℝ)) * ω) =
      Real.exp (-t * ω) * Real.exp (ω / 2) := by
    rw [show -(t - (1 / 2 : ℝ)) * ω = -t * ω + ω / 2 by ring, Real.exp_add]
  have hnegative : Real.exp (-(ω / 2)) = Real.exp (ω / 2) * Real.exp (-ω) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have hden : 2 * Real.cosh (ω / 2) =
      Real.exp (ω / 2) * (1 + Real.exp (-ω)) := by
    rw [Real.cosh_eq, hnegative]
    ring
  unfold fermionicKernel centeredKernel
  rw [hnum, hden]
  field_simp [ne_of_gt (Real.exp_pos (ω / 2)), fermionicKernel_denominator_ne]

end Fermionic
end GECPKernelStructure
