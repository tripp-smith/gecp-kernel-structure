import GECPKernelStructure.GECP.SignRegularity
import Mathlib.Tactic

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v

variable {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]

/-- A kernel has noncrossing support when two oppositely ordered row-to-column
edges cannot both carry nonzero weight. -/
def HasNoncrossingSupport (K : Kernel α β ℝ) : Prop :=
  ∀ ⦃x₁ x₂ : α⦄ ⦃y₁ y₂ : β⦄, x₁ < x₂ → y₁ < y₂ →
    K x₂ y₁ = 0 ∨ K x₁ y₂ = 0

/-- A nonnegative kernel with noncrossing support is totally nonnegative.
Every nonidentity term in an ordered minor's Leibniz expansion contains a
crossing pair and therefore vanishes. -/
theorem totallyNonnegative_of_nonnegative_of_noncrossing
    {K : Kernel α β ℝ}
    (entry_nonnegative : ∀ x y, 0 ≤ K x y)
    (noncrossing : HasNoncrossingSupport K) :
    TotallyNonnegative K := by
  intro n rows columns rows_mono columns_mono
  rw [Matrix.det_apply']
  classical
  rw [Finset.sum_eq_single (1 : Equiv.Perm (Fin n))]
  · simp only [Equiv.Perm.coe_one, id_eq, Equiv.Perm.sign_one, Units.val_one,
      Int.cast_one, one_mul]
    exact Finset.prod_nonneg fun i _ => entry_nonnegative (rows i) (columns i)
  · intro σ _ σ_ne_one
    have not_strict_mono : ¬StrictMono σ := by
      intro strict_mono
      have : σ = 1 := by
        apply Equiv.ext
        intro i
        simpa using StrictMono.apply_eq strict_mono (x := i)
      exact σ_ne_one this
    rw [StrictMono] at not_strict_mono
    push Not at not_strict_mono
    obtain ⟨i, j, i_lt_j, not_lt⟩ := not_strict_mono
    have sigma_ne : σ i ≠ σ j := σ.injective.ne i_lt_j.ne
    have sigma_reverse : σ j < σ i :=
      lt_of_le_of_ne not_lt sigma_ne.symm
    obtain crossed_zero | crossed_zero :=
      noncrossing (rows_mono sigma_reverse) (columns_mono i_lt_j)
    · rw [Finset.prod_eq_zero (Finset.mem_univ i) crossed_zero, mul_zero]
    · rw [Finset.prod_eq_zero (Finset.mem_univ j) crossed_zero, mul_zero]
  · simp

end GECP
end GECPKernelStructure
