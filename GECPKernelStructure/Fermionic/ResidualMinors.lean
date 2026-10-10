import GECPKernelStructure.GECP.AugmentedDeterminant
import GECPKernelStructure.Fermionic.ExponentialTotalPositivity

namespace GECPKernelStructure
namespace Fermionic

open GECP

universe z

/-- Fermionic residual minors retain the exact orientation-aware shifted sign. -/
theorem fermionicKernel_finalResidual_minor_oriented_pos
    (run : Run fermionicKernel) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows columns : γ → ℝ)
    (rows_injective : Function.Injective (run.augmentedRow γ rows))
    (columns_injective : Function.Injective (run.augmentedColumn γ columns)) :
    0 < expKernelSignature (Fintype.card (run.AugmentedIndex γ)) *
        run.augmentedRowOrientation γ rows *
        run.augmentedColumnOrientation γ columns *
        run.selectedCore.det *
        Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) :=
  strictSignRegular_finalResidual_minor_oriented_pos
    fermionicKernel_strictSignRegular run γ rows columns
    rows_injective columns_injective

/-- Every injectively augmented fermionic run has nonsingular final-residual minors. -/
theorem fermionicKernel_finalResidual_minor_ne_zero
    (run : Run fermionicKernel) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows columns : γ → ℝ)
    (rows_injective : Function.Injective (run.augmentedRow γ rows))
    (columns_injective : Function.Injective (run.augmentedColumn γ columns)) :
    Matrix.det (fun i j => run.finalResidual (rows i) (columns j)) ≠ 0 :=
  strictSignRegular_finalResidual_minor_ne_zero
    fermionicKernel_strictSignRegular run γ rows columns
    rows_injective columns_injective

end Fermionic
end GECPKernelStructure
