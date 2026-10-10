import GECPKernelStructure.Fermionic.CompoundResidual

namespace GECPKernelStructure

open Matrix

universe u v

namespace GECP

/-- Two values represented as a `Fin 2` tuple in the displayed order. -/
def orderedPair {α : Type u} (first second : α) : Fin 2 → α :=
  ![first, second]

/-- Same-sign subtraction costs at most the larger absolute value, not their sum. -/
theorem abs_sub_le_max_of_mul_nonneg {a b : ℝ} (same_sign : 0 ≤ a * b) :
    abs (a - b) ≤ max (abs a) (abs b) := by
  rw [abs_le]
  rcases mul_nonneg_iff.mp same_sign with ⟨a_nonneg, b_nonneg⟩ | ⟨a_nonpos, b_nonpos⟩
  · rw [abs_of_nonneg a_nonneg, abs_of_nonneg b_nonneg]
    constructor
    · have a_le := le_max_left a b
      have b_le := le_max_right a b
      linarith
    · have a_le := le_max_left a b
      have b_le := le_max_right a b
      linarith
  · rw [abs_of_nonpos a_nonpos, abs_of_nonpos b_nonpos]
    constructor
    · have neg_a_le := le_max_left (-a) (-b)
      have neg_b_le := le_max_right (-a) (-b)
      linarith
    · have neg_a_le := le_max_left (-a) (-b)
      have neg_b_le := le_max_right (-a) (-b)
      linarith

namespace Run

/--
Selected-core Desnanot--Jacobi identity for two appended rows and columns,
proved through exact Schur factorization rather than an inverse.
-/
theorem augmentedCore_fin_two_mul_selectedCore_det
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (row₀ row₁ : α) (column₀ column₁ : β) :
    (run.augmentedCore (Fin 2) (orderedPair row₀ row₁)
        (orderedPair column₀ column₁)).det * run.selectedCore.det =
      (run.borderedCore row₀ column₀).det *
          (run.borderedCore row₁ column₁).det -
        (run.borderedCore row₀ column₁).det *
          (run.borderedCore row₁ column₀).det := by
  rw [run.augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor
    (Fin 2) (orderedPair row₀ row₁) (orderedPair column₀ column₁)]
  rw [run.borderedCore_det_eq_selectedCore_det_mul_finalResidual row₀ column₀,
    run.borderedCore_det_eq_selectedCore_det_mul_finalResidual row₁ column₁,
    run.borderedCore_det_eq_selectedCore_det_mul_finalResidual row₀ column₁,
    run.borderedCore_det_eq_selectedCore_det_mul_finalResidual row₁ column₀]
  let residualMinor : Matrix (Fin 2) (Fin 2) ℝ :=
    fun i j => run.finalResidual (orderedPair row₀ row₁ i)
      (orderedPair column₀ column₁ j)
  change (run.selectedCore.det * residualMinor.det) * run.selectedCore.det = _
  rw [Matrix.det_fin_two residualMinor]
  simp only [residualMinor, orderedPair, Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

end Run
end GECP

namespace Fermionic

open GECP

/--
Strict sign regularity turns the exact two-border identity into a factor-one
absolute condensation bound.
-/
theorem strictSignRegular_augmentedCore_fin_two_condensation
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (row₀ row₁ : α) (column₀ column₁ : β) :
    abs ((run.augmentedCore (Fin 2) (orderedPair row₀ row₁)
        (orderedPair column₀ column₁)).det * run.selectedCore.det) ≤
      max
        (abs ((run.borderedCore row₀ column₀).det *
          (run.borderedCore row₁ column₁).det))
        (abs ((run.borderedCore row₀ column₁).det *
          (run.borderedCore row₁ column₀).det)) := by
  rw [run.augmentedCore_fin_two_mul_selectedCore_det]
  apply abs_sub_le_max_of_mul_nonneg
  exact strictSignRegular_borderedMinorSignCoherent regular run
    row₁ column₁ row₀ column₀

/-- Every fermionic run satisfies the factor-one two-border condensation bound. -/
theorem fermionicKernel_augmentedCore_fin_two_condensation
    (run : Run fermionicKernel) (row₀ row₁ column₀ column₁ : ℝ) :
    abs ((run.augmentedCore (Fin 2) (orderedPair row₀ row₁)
        (orderedPair column₀ column₁)).det * run.selectedCore.det) ≤
      max
        (abs ((run.borderedCore row₀ column₀).det *
          (run.borderedCore row₁ column₁).det))
        (abs ((run.borderedCore row₀ column₁).det *
          (run.borderedCore row₁ column₀).det)) :=
  strictSignRegular_augmentedCore_fin_two_condensation
    fermionicKernel_strictSignRegular run row₀ row₁ column₀ column₁

end Fermionic
end GECPKernelStructure
