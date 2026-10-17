import GECPKernelStructure.Matrix.WeightedCenteredGroups
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open Finset Matrix

/-- One pivot's positive Taylor feature weight at degree `k`. -/
noncomputable def onePivotTaylorWeight {n : ℕ}
    (x0 y0 : ℝ) (k : Fin n) : ℝ :=
  (x0 * y0) ^ k.1 / (k.1.factorial : ℝ)

/-- The finite positive Taylor kernel through the first `n` features. -/
noncomputable def finitePositiveTaylorKernel (n : ℕ) (x y : ℝ) : ℝ :=
  ∑ k : Fin n, (x * y) ^ k.1 / (k.1.factorial : ℝ)

/-- The left one-pivot bordered coordinate at one Taylor degree. -/
noncomputable def onePivotTaylorLeftCoordinate {n : ℕ}
    (x0 x1 y0 : ℝ) (k : Fin n) : ℝ :=
  finitePositiveTaylorKernel n x0 y0 *
      (x1 ^ k.1 / (k.1.factorial : ℝ)) -
    (x0 ^ k.1 / (k.1.factorial : ℝ)) *
      finitePositiveTaylorKernel n x1 y0

/-- The right one-pivot bordered coordinate at the same Taylor degree. -/
noncomputable def onePivotTaylorRightCoordinate {n : ℕ}
    (x0 y0 y1 : ℝ) (k : Fin n) : ℝ :=
  finitePositiveTaylorKernel n x0 y0 * y1 ^ k.1 -
    finitePositiveTaylorKernel n x0 y1 * y0 ^ k.1

/-- The matched one-pivot bordered product at one Taylor degree. -/
noncomputable def onePivotTaylorProduct {n : ℕ}
    (x0 x1 y0 y1 : ℝ) (k : Fin n) : ℝ :=
  onePivotTaylorLeftCoordinate x0 x1 y0 k *
    onePivotTaylorRightCoordinate x0 y0 y1 k

@[simp]
theorem weightedMass_onePivotTaylorWeight {n : ℕ} (x0 y0 : ℝ) :
    weightedMass Finset.univ (onePivotTaylorWeight (n := n) x0 y0) =
      finitePositiveTaylorKernel n x0 y0 := by
  rfl

theorem weightedMoment_onePivotTaylorWeight_rowRatio
    {n : ℕ} {x0 : ℝ} (x0_ne : x0 ≠ 0) (x1 y0 : ℝ) :
    weightedMoment Finset.univ (onePivotTaylorWeight (n := n) x0 y0)
        (fun k ↦ (x1 / x0) ^ k.1) =
      finitePositiveTaylorKernel n x1 y0 := by
  unfold weightedMoment onePivotTaylorWeight finitePositiveTaylorKernel
  apply sum_congr rfl
  intro k hk
  change
    (x0 * y0) ^ k.1 / (k.1.factorial : ℝ) * (x1 / x0) ^ k.1 =
      (x1 * y0) ^ k.1 / (k.1.factorial : ℝ)
  rw [div_pow]
  field_simp [x0_ne]
  ring

theorem weightedMoment_onePivotTaylorWeight_columnRatio
    {n : ℕ} {y0 : ℝ} (y0_ne : y0 ≠ 0) (x0 y1 : ℝ) :
    weightedMoment Finset.univ (onePivotTaylorWeight (n := n) x0 y0)
        (fun k ↦ (y1 / y0) ^ k.1) =
      finitePositiveTaylorKernel n x0 y1 := by
  unfold weightedMoment onePivotTaylorWeight finitePositiveTaylorKernel
  apply sum_congr rfl
  intro k hk
  change
    (x0 * y0) ^ k.1 / (k.1.factorial : ℝ) * (y1 / y0) ^ k.1 =
      (x0 * y1) ^ k.1 / (k.1.factorial : ℝ)
  rw [div_pow]
  field_simp [y0_ne]
  ring

theorem onePivotTaylorLeftCoordinate_eq_centered
    {n : ℕ} {x0 : ℝ} (x0_ne : x0 ≠ 0)
    (x1 y0 : ℝ) (k : Fin n) :
    onePivotTaylorLeftCoordinate x0 x1 y0 k =
      (x0 ^ k.1 / (k.1.factorial : ℝ)) *
      (weightedMass Finset.univ (onePivotTaylorWeight (n := n) x0 y0) *
            (x1 / x0) ^ k.1 -
          weightedMoment Finset.univ (onePivotTaylorWeight (n := n) x0 y0)
            (fun j : Fin n ↦ (x1 / x0) ^ j.1)) := by
  rw [weightedMass_onePivotTaylorWeight,
    weightedMoment_onePivotTaylorWeight_rowRatio x0_ne]
  unfold onePivotTaylorLeftCoordinate
  rw [div_pow]
  field_simp

theorem onePivotTaylorRightCoordinate_eq_centered
    {n : ℕ} {y0 : ℝ} (y0_ne : y0 ≠ 0)
    (x0 y1 : ℝ) (k : Fin n) :
    onePivotTaylorRightCoordinate x0 y0 y1 k =
      y0 ^ k.1 *
        (weightedMass Finset.univ (onePivotTaylorWeight (n := n) x0 y0) *
            (y1 / y0) ^ k.1 -
          weightedMoment Finset.univ (onePivotTaylorWeight (n := n) x0 y0)
            (fun j : Fin n ↦ (y1 / y0) ^ j.1)) := by
  rw [weightedMass_onePivotTaylorWeight,
    weightedMoment_onePivotTaylorWeight_columnRatio y0_ne]
  unfold onePivotTaylorRightCoordinate
  rw [div_pow]
  field_simp

/-- Every matched Taylor bordered product is exactly the corresponding
weighted centered product of the row and column geometric ratios. -/
theorem onePivotTaylorProduct_eq_geometricCenteredProduct
    {n : ℕ} {x0 y0 : ℝ} (x0_ne : x0 ≠ 0) (y0_ne : y0 ≠ 0)
    (x1 y1 : ℝ) (k : Fin n) :
    onePivotTaylorProduct x0 x1 y0 y1 k =
      geometricCenteredProduct (onePivotTaylorWeight x0 y0)
        (x1 / x0) (y1 / y0) k := by
  rw [onePivotTaylorProduct, onePivotTaylorLeftCoordinate_eq_centered x0_ne,
    onePivotTaylorRightCoordinate_eq_centered y0_ne]
  unfold geometricCenteredProduct onePivotTaylorWeight
  rw [mul_pow]
  ring

theorem onePivotTaylorWeight_pos
    {n : ℕ} {x0 y0 : ℝ} (x0_pos : 0 < x0) (y0_pos : 0 < y0)
    (k : Fin n) :
    0 < onePivotTaylorWeight x0 y0 k := by
  unfold onePivotTaylorWeight
  positivity

/-- Every nonempty degree prefix of an arbitrary finite one-pivot positive
Taylor truncation has a nonnegative grouped bordered pairing. -/
theorem onePivotTaylorProduct_prefix_nonneg
    {n cutoff : ℕ} (n_pos : 0 < n) (cutoff_pos : 0 < cutoff)
    {x0 x1 y0 y1 : ℝ}
    (x0_pos : 0 < x0) (x_order : x0 ≤ x1)
    (y0_pos : 0 < y0) (y_order : y0 ≤ y1) :
    0 ≤ ∑ k ∈ degreePrefix n cutoff,
      onePivotTaylorProduct x0 x1 y0 y1 k := by
  have x0_ne : x0 ≠ 0 := ne_of_gt x0_pos
  have y0_ne : y0 ≠ 0 := ne_of_gt y0_pos
  have rowRatio_ge_one : 1 ≤ x1 / x0 := (one_le_div x0_pos).mpr x_order
  have columnRatio_ge_one : 1 ≤ y1 / y0 := (one_le_div y0_pos).mpr y_order
  rw [sum_congr rfl fun k hk ↦
    onePivotTaylorProduct_eq_geometricCenteredProduct x0_ne y0_ne x1 y1 k]
  exact geometricCenteredProduct_prefix_nonneg n_pos cutoff_pos
    (onePivotTaylorWeight x0 y0) (onePivotTaylorWeight_pos x0_pos y0_pos)
    rowRatio_ge_one columnRatio_ge_one

/-- Every nonempty degree suffix satisfies the symmetric grouped sign. -/
theorem onePivotTaylorProduct_suffix_nonneg
    {n cutoff : ℕ} (cutoff_lt : cutoff < n)
    {x0 x1 y0 y1 : ℝ}
    (x0_pos : 0 < x0) (x_order : x0 ≤ x1)
    (y0_pos : 0 < y0) (y_order : y0 ≤ y1) :
    0 ≤ ∑ k ∈ degreeSuffix n cutoff,
      onePivotTaylorProduct x0 x1 y0 y1 k := by
  have x0_ne : x0 ≠ 0 := ne_of_gt x0_pos
  have y0_ne : y0 ≠ 0 := ne_of_gt y0_pos
  have rowRatio_ge_one : 1 ≤ x1 / x0 := (one_le_div x0_pos).mpr x_order
  have columnRatio_ge_one : 1 ≤ y1 / y0 := (one_le_div y0_pos).mpr y_order
  rw [sum_congr rfl fun k hk ↦
    onePivotTaylorProduct_eq_geometricCenteredProduct x0_ne y0_ne x1 y1 k]
  exact geometricCenteredProduct_suffix_nonneg cutoff_lt
    (onePivotTaylorWeight x0 y0) (onePivotTaylorWeight_pos x0_pos y0_pos)
    rowRatio_ge_one columnRatio_ge_one

/-- Every nontrivial two-level contiguous feature step has a nonnegative
one-pivot bordered pairing for arbitrary finite Taylor order. -/
theorem onePivotTaylorProduct_step_nonneg
    {n cutoff : ℕ} (cutoff_pos : 0 < cutoff) (cutoff_lt : cutoff < n)
    {x0 x1 y0 y1 weight : ℝ}
    (x0_pos : 0 < x0) (x_order : x0 ≤ x1)
    (y0_pos : 0 < y0) (y_order : y0 ≤ y1)
    (weight_nonneg : 0 ≤ weight) :
    0 ≤ (∑ k ∈ degreePrefix n cutoff,
        onePivotTaylorProduct x0 x1 y0 y1 k) +
      weight * (∑ k ∈ degreeSuffix n cutoff,
        onePivotTaylorProduct x0 x1 y0 y1 k) := by
  have prefix_nonneg := onePivotTaylorProduct_prefix_nonneg
    (Nat.zero_lt_of_lt cutoff_lt) cutoff_pos x0_pos x_order y0_pos y_order
  have suffix_nonneg := onePivotTaylorProduct_suffix_nonneg
    cutoff_lt x0_pos x_order y0_pos y_order
  exact add_nonneg prefix_nonneg (mul_nonneg weight_nonneg suffix_nonneg)

end Fermionic
end GECPKernelStructure
