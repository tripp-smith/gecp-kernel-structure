import GECPKernelStructure.Fermionic.FiniteRunRate
import Mathlib.Topology.Order.Compact

namespace GECPKernelStructure.Fermionic

open GECP

/-- Each finite continuous residual has a complete pivot on the physical
compact rectangle. A zero maximum certifies termination, not another update. -/
theorem fermionicKernel_finiteRun_exists_completePivot (Λ : ℝ) (hΛ : 0 ≤ Λ)
    (run : Run fermionicKernel) :
    ∃ p q : ℝ, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) run.finalResidual p q := by
  have hcompact : IsCompact (Set.Icc (0 : ℝ) 1 ×ˢ Set.Icc (-Λ) Λ) :=
    isCompact_Icc.prod isCompact_Icc
  have hnonempty : (Set.Icc (0 : ℝ) 1 ×ˢ Set.Icc (-Λ) Λ).Nonempty :=
    ⟨(0, 0), ⟨by constructor <;> norm_num, by constructor <;> linarith⟩⟩
  obtain ⟨z, hz, hmax⟩ := hcompact.exists_isMaxOn hnonempty
    (run.continuous_finalResidual fermionicKernel_continuous).abs.continuousOn
  refine ⟨z.1, z.2, hz.1, hz.2, ?_⟩
  intro t omega ht homega
  exact hmax (show (t, omega) ∈ Set.Icc (0 : ℝ) 1 ×ˢ Set.Icc (-Λ) Λ from ⟨ht, homega⟩)

end GECPKernelStructure.Fermionic
