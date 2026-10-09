import GECPKernelStructure.Fermionic.Derivatives
import GECPKernelStructure.GECP.Definitions
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

/--
The fermionic kernel has an exact multiplicative-Monge cross ratio. Its
frequency-only normalization cancels from the four-point ratio.
-/
theorem fermionicKernel_cross_product (t t₀ ω ω₀ : ℝ) :
    fermionicKernel t ω₀ * fermionicKernel t₀ ω =
      Real.exp ((t - t₀) * (ω - ω₀)) *
        (fermionicKernel t₀ ω₀ * fermionicKernel t ω) := by
  unfold fermionicKernel
  field_simp [fermionicKernel_denominator_ne]
  rw [← Real.exp_add, ← Real.exp_add, ← Real.exp_add]
  congr 1
  ring

/--
One exact elimination step factors into the original positive kernel and a
scalar depending only on the pivot displacement product.
-/
theorem fermionicKernel_residualUpdate_factor (t t₀ ω ω₀ : ℝ) :
    residualUpdate fermionicKernel t₀ ω₀ (fermionicKernel_pos t₀ ω₀).ne' t ω =
      fermionicKernel t ω *
        (1 - Real.exp ((t - t₀) * (ω - ω₀))) := by
  have cross := fermionicKernel_cross_product t t₀ ω ω₀
  have ratio :
      fermionicKernel t ω₀ * fermionicKernel t₀ ω / fermionicKernel t₀ ω₀ =
        Real.exp ((t - t₀) * (ω - ω₀)) * fermionicKernel t ω := by
    apply (div_eq_iff (fermionicKernel_pos t₀ ω₀).ne').2
    calc
      fermionicKernel t ω₀ * fermionicKernel t₀ ω =
          Real.exp ((t - t₀) * (ω - ω₀)) *
            (fermionicKernel t₀ ω₀ * fermionicKernel t ω) := cross
      _ = (Real.exp ((t - t₀) * (ω - ω₀)) * fermionicKernel t ω) *
          fermionicKernel t₀ ω₀ := by ring
  unfold residualUpdate
  rw [ratio]
  ring

/-- A fermionic update is nonnegative when its two pivot displacements have opposite signs. -/
theorem fermionicKernel_residualUpdate_nonneg_of_cross_nonpos {t t₀ ω ω₀ : ℝ}
    (cross_nonpos : (t - t₀) * (ω - ω₀) ≤ 0) :
    0 ≤ residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω := by
  rw [fermionicKernel_residualUpdate_factor]
  apply mul_nonneg (fermionicKernel_pos t ω).le
  rw [sub_nonneg, Real.exp_le_one_iff]
  exact cross_nonpos

/--
On the physical time interval, an oppositely ordered update is bounded by the
negated displacement product, its dimensionless time-frequency cross area.
-/
theorem fermionicKernel_residualUpdate_le_crossArea {t t₀ ω ω₀ : ℝ}
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (cross_nonpos : (t - t₀) * (ω - ω₀) ≤ 0) :
    abs (residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω) ≤
        -((t - t₀) * (ω - ω₀)) := by
  let area := -((t - t₀) * (ω - ω₀))
  have area_nonneg : 0 ≤ area := neg_nonneg.mpr cross_nonpos
  have exponent_eq : (t - t₀) * (ω - ω₀) = -area := by
    dsimp [area]
    ring
  have factor_nonneg : 0 ≤ 1 - Real.exp (-area) := by
    rw [sub_nonneg, Real.exp_le_one_iff]
    exact neg_nonpos.mpr area_nonneg
  have factor_le_area : 1 - Real.exp (-area) ≤ area := by
    linarith [Real.one_sub_le_exp_neg area]
  change abs (residualUpdate fermionicKernel t₀ ω₀
    (fermionicKernel_pos t₀ ω₀).ne' t ω) ≤ area
  rw [fermionicKernel_residualUpdate_factor, exponent_eq]
  rw [abs_of_nonneg (mul_nonneg (fermionicKernel_pos t ω).le factor_nonneg)]
  calc
    fermionicKernel t ω * (1 - Real.exp (-area)) ≤
        1 * (1 - Real.exp (-area)) :=
      mul_le_mul_of_nonneg_right
        (fermionicKernel_le_one time_nonneg time_le_one) factor_nonneg
    _ ≤ 1 * area := mul_le_mul_of_nonneg_left factor_le_area zero_le_one
    _ = area := one_mul _

/--
An oppositely ordered update contracts locally by one half whenever its
dimensionless cross area is at most `log 2`.
-/
theorem fermionicKernel_residualUpdate_le_half_of_crossArea {t t₀ ω ω₀ : ℝ}
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (cross_nonpos : (t - t₀) * (ω - ω₀) ≤ 0)
    (area_le : -((t - t₀) * (ω - ω₀)) ≤ Real.log 2) :
    abs (residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω) ≤ 1 / 2 := by
  let area := -((t - t₀) * (ω - ω₀))
  have area_nonneg : 0 ≤ area := neg_nonneg.mpr cross_nonpos
  have exponent_eq : (t - t₀) * (ω - ω₀) = -area := by
    dsimp [area]
    ring
  have exp_neg_log_two : Real.exp (-Real.log 2) = (1 / 2 : ℝ) := by
    rw [Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 2)]
    norm_num
  have half_le_exp : (1 / 2 : ℝ) ≤ Real.exp (-area) := by
    rw [← exp_neg_log_two]
    exact Real.exp_le_exp.mpr (neg_le_neg (show area ≤ Real.log 2 from area_le))
  have factor_nonneg : 0 ≤ 1 - Real.exp (-area) := by
    rw [sub_nonneg, Real.exp_le_one_iff]
    exact neg_nonpos.mpr area_nonneg
  have factor_le_half : 1 - Real.exp (-area) ≤ (1 / 2 : ℝ) := by
    linarith
  rw [fermionicKernel_residualUpdate_factor, exponent_eq]
  rw [abs_of_nonneg (mul_nonneg (fermionicKernel_pos t ω).le factor_nonneg)]
  calc
    fermionicKernel t ω * (1 - Real.exp (-area)) ≤
        1 * (1 - Real.exp (-area)) :=
      mul_le_mul_of_nonneg_right
        (fermionicKernel_le_one time_nonneg time_le_one) factor_nonneg
    _ ≤ 1 * (1 / 2 : ℝ) := mul_le_mul_of_nonneg_left factor_le_half zero_le_one
    _ = 1 / 2 := one_mul _

/-- An update from the northwest corner of an oppositely ordered rectangle is nonnegative. -/
theorem fermionicKernel_cornerResidual_nonneg {t t₀ ω ω₀ : ℝ}
    (time_order : t₀ ≤ t) (frequency_order : ω ≤ ω₀) :
    0 ≤ residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω := by
  apply fermionicKernel_residualUpdate_nonneg_of_cross_nonpos
  exact mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr time_order)
    (sub_nonpos.mpr frequency_order)

/--
On the physical time interval, a northwest-corner update is bounded by the
dimensionless time-frequency area between the pivot and evaluation point.
-/
theorem fermionicKernel_cornerResidual_le_area {t t₀ ω ω₀ : ℝ}
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (time_order : t₀ ≤ t) (frequency_order : ω ≤ ω₀) :
    abs (residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω) ≤
        (t - t₀) * (ω₀ - ω) := by
  have cross_nonpos : (t - t₀) * (ω - ω₀) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr time_order)
      (sub_nonpos.mpr frequency_order)
  have bound := fermionicKernel_residualUpdate_le_crossArea
    time_nonneg time_le_one cross_nonpos
  convert bound using 1
  ring

/--
A northwest-corner update contracts locally by one half whenever its
dimensionless time-frequency area is at most `log 2`.
-/
theorem fermionicKernel_cornerResidual_le_half {t t₀ ω ω₀ : ℝ}
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (time_order : t₀ ≤ t) (frequency_order : ω ≤ ω₀)
    (area_le : (t - t₀) * (ω₀ - ω) ≤ Real.log 2) :
    abs (residualUpdate fermionicKernel t₀ ω₀
      (fermionicKernel_pos t₀ ω₀).ne' t ω) ≤ 1 / 2 := by
  have cross_nonpos : (t - t₀) * (ω - ω₀) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr time_order)
      (sub_nonpos.mpr frequency_order)
  apply fermionicKernel_residualUpdate_le_half_of_crossArea
    time_nonneg time_le_one cross_nonpos
  convert area_le using 1
  ring

end Fermionic
end GECPKernelStructure
