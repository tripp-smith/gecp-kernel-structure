import GECPKernelStructure.Fermionic.CrossRatioLocalization
import Mathlib.Analysis.Convex.SpecificFunctions.Pow
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

/-- The fermionic kernel in the positive coordinate `z = exp (-omega)`. -/
theorem fermionicKernel_eq_rpow (t ω : ℝ) :
    fermionicKernel t ω =
      Real.exp (-ω) ^ t / (1 + Real.exp (-ω)) := by
  unfold fermionicKernel
  congr 1
  rw [Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp]
  congr 1
  ring

/-- Linear interpolation of the concave power `z^t` between positive endpoints. -/
noncomputable def powerSecant (a b t z : ℝ) : ℝ :=
  ((b - z) / (b - a)) * a ^ t + ((z - a) / (b - a)) * b ^ t

/-- For an exponent in `[0,1]`, the power secant lies below the power function. -/
theorem powerSecant_le_rpow {a b t z : ℝ}
    (a_nonneg : 0 ≤ a) (endpoints : a < b)
    (z_lower : a ≤ z) (z_upper : z ≤ b)
    (exponent_nonneg : 0 ≤ t) (exponent_le_one : t ≤ 1) :
    powerSecant a b t z ≤ z ^ t := by
  let leftWeight := (b - z) / (b - a)
  let rightWeight := (z - a) / (b - a)
  have denominator_pos : 0 < b - a := sub_pos.mpr endpoints
  have left_nonneg : 0 ≤ leftWeight := by
    exact div_nonneg (sub_nonneg.mpr z_upper) denominator_pos.le
  have right_nonneg : 0 ≤ rightWeight := by
    exact div_nonneg (sub_nonneg.mpr z_lower) denominator_pos.le
  have weights_sum : leftWeight + rightWeight = 1 := by
    dsimp [leftWeight, rightWeight]
    field_simp
    ring
  have weighted_point : leftWeight * a + rightWeight * b = z := by
    dsimp [leftWeight, rightWeight]
    field_simp
    ring
  have b_nonneg : 0 ≤ b := a_nonneg.trans endpoints.le
  have concavity := (Real.concaveOn_rpow exponent_nonneg exponent_le_one).2
    a_nonneg b_nonneg left_nonneg right_nonneg weights_sum
  simp only [smul_eq_mul] at concavity
  rw [weighted_point] at concavity
  simpa [powerSecant, leftWeight, rightWeight] using concavity

/-- The same power secant dominates the line through the right endpoint and the origin. -/
theorem mul_rpow_sub_one_le_powerSecant {a b t z : ℝ}
    (a_pos : 0 < a) (endpoints : a < b)
    (z_lower : a ≤ z) (z_upper : z ≤ b) (exponent_le_one : t ≤ 1) :
    z * b ^ (t - 1) ≤ powerSecant a b t z := by
  let leftWeight := (b - z) / (b - a)
  let rightWeight := (z - a) / (b - a)
  have denominator_pos : 0 < b - a := sub_pos.mpr endpoints
  have left_nonneg : 0 ≤ leftWeight := by
    exact div_nonneg (sub_nonneg.mpr z_upper) denominator_pos.le
  have right_nonneg : 0 ≤ rightWeight := by
    exact div_nonneg (sub_nonneg.mpr z_lower) denominator_pos.le
  have weighted_point : leftWeight * a + rightWeight * b = z := by
    dsimp [leftWeight, rightWeight]
    field_simp
    ring
  have powers_order : b ^ (t - 1) ≤ a ^ (t - 1) :=
    Real.rpow_le_rpow_of_nonpos a_pos endpoints.le (sub_nonpos.mpr exponent_le_one)
  have a_factor : a * a ^ (t - 1) = a ^ t := by
    calc
      a * a ^ (t - 1) = a ^ 1 * a ^ (t - 1) := by rw [Real.rpow_one]
      _ = a ^ (1 + (t - 1)) := (Real.rpow_add a_pos 1 (t - 1)).symm
      _ = a ^ t := by ring_nf
  have b_pos : 0 < b := a_pos.trans endpoints
  have b_factor : b * b ^ (t - 1) = b ^ t := by
    calc
      b * b ^ (t - 1) = b ^ 1 * b ^ (t - 1) := by rw [Real.rpow_one]
      _ = b ^ (1 + (t - 1)) := (Real.rpow_add b_pos 1 (t - 1)).symm
      _ = b ^ t := by ring_nf
  calc
    z * b ^ (t - 1) =
        (leftWeight * a + rightWeight * b) * b ^ (t - 1) := by rw [weighted_point]
    _ = leftWeight * (a * b ^ (t - 1)) +
        rightWeight * (b * b ^ (t - 1)) := by ring
    _ ≤ leftWeight * (a * a ^ (t - 1)) +
        rightWeight * (b * b ^ (t - 1)) := by
          exact add_le_add
            (mul_le_mul_of_nonneg_left
              (mul_le_mul_of_nonneg_left powers_order a_pos.le) left_nonneg)
            le_rfl
    _ = powerSecant a b t z := by
      rw [a_factor, b_factor]
      rfl

/--
After the upper-frequency corner pivot at time zero, the residual is a
difference of positive-coordinate powers with a common logistic denominator.
-/
theorem fermionicKernel_firstCornerResidual_eq_rpow_sub (t ω ωHigh : ℝ) :
    residualUpdate fermionicKernel 0 ωHigh
      (fermionicKernel_pos 0 ωHigh).ne' t ω =
        (Real.exp (-ω) ^ t - Real.exp (-ωHigh) ^ t) /
          (1 + Real.exp (-ω)) := by
  have product_identity :
      Real.exp (-ω) ^ t * Real.exp ((t - 0) * (ω - ωHigh)) =
        Real.exp (-ωHigh) ^ t := by
    rw [Real.rpow_def_of_pos (Real.exp_pos _),
      Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp, Real.log_exp,
      ← Real.exp_add]
    congr 1
    ring
  rw [fermionicKernel_residualUpdate_factor, fermionicKernel_eq_rpow]
  field_simp
  rw [mul_sub, mul_one, product_identity]

/--
The reflected lower-frequency/time-one corner update has the corresponding
right-endpoint power-line form.
-/
theorem fermionicKernel_secondCornerResidual_eq_rpow_sub (t ω ωLow : ℝ) :
    residualUpdate fermionicKernel 1 ωLow
      (fermionicKernel_pos 1 ωLow).ne' t ω =
        (Real.exp (-ω) ^ t -
          Real.exp (-ω) * Real.exp (-ωLow) ^ (t - 1)) /
          (1 + Real.exp (-ω)) := by
  have product_identity :
      Real.exp (-ω) ^ t * Real.exp ((t - 1) * (ω - ωLow)) =
        Real.exp (-ω) * Real.exp (-ωLow) ^ (t - 1) := by
    rw [Real.rpow_def_of_pos (Real.exp_pos _),
      Real.rpow_def_of_pos (Real.exp_pos _), Real.log_exp, Real.log_exp,
      ← Real.exp_add, ← Real.exp_add]
    congr 1
    ring
  rw [fermionicKernel_residualUpdate_factor, fermionicKernel_eq_rpow]
  field_simp
  rw [mul_sub, mul_one, product_identity]

/-- The lower-frequency endpoint is a positive second pivot after the first corner update. -/
theorem fermionicKernel_firstCorner_secondPivot_pos {ωLow ωHigh : ℝ}
    (frequency_order : ωLow < ωHigh) :
    0 < residualUpdate fermionicKernel 0 ωHigh
      (fermionicKernel_pos 0 ωHigh).ne' 1 ωLow := by
  rw [fermionicKernel_firstCornerResidual_eq_rpow_sub]
  simp only [Real.rpow_one]
  apply div_pos
  · exact sub_pos.mpr (Real.exp_lt_exp.mpr (neg_lt_neg frequency_order))
  · positivity

/--
Two endpoint-time/corner-frequency pivots turn the fermionic residual into the
error of linear interpolation of `z^t`, divided by the logistic denominator.
-/
theorem fermionicKernel_twoCornerResidual_eq_powerSecantError
    {ωLow ωHigh : ℝ} (frequency_order : ωLow < ωHigh) (t ω : ℝ) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω =
        (Real.exp (-ω) ^ t -
          powerSecant (Real.exp (-ωHigh)) (Real.exp (-ωLow)) t
            (Real.exp (-ω))) /
          (1 + Real.exp (-ω)) := by
  let firstResidual := residualUpdate fermionicKernel 0 ωHigh
    (fermionicKernel_pos 0 ωHigh).ne'
  have first_formula (s η : ℝ) :
      firstResidual s η =
        (Real.exp (-η) ^ s - Real.exp (-ωHigh) ^ s) /
          (1 + Real.exp (-η)) :=
    fermionicKernel_firstCornerResidual_eq_rpow_sub s η ωHigh
  have transformed_order : Real.exp (-ωHigh) < Real.exp (-ωLow) :=
    Real.exp_lt_exp.mpr (neg_lt_neg frequency_order)
  have endpoint_difference_ne :
      Real.exp (-ωLow) - Real.exp (-ωHigh) ≠ 0 :=
    (sub_pos.mpr transformed_order).ne'
  change firstResidual t ω -
      firstResidual t ωLow * firstResidual 1 ω / firstResidual 1 ωLow = _
  rw [first_formula, first_formula, first_formula, first_formula]
  simp only [Real.rpow_one]
  unfold powerSecant
  field_simp [endpoint_difference_ne]
  ring

/-- The exact two-corner residual is nonnegative throughout its frequency band. -/
theorem fermionicKernel_twoCornerResidual_nonneg
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : ωLow ≤ ω) (frequency_upper : ω ≤ ωHigh) :
    0 ≤ residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω := by
  rw [fermionicKernel_twoCornerResidual_eq_powerSecantError frequency_order]
  apply div_nonneg
  · apply sub_nonneg.mpr
    apply powerSecant_le_rpow (Real.exp_pos _).le
      (Real.exp_lt_exp.mpr (neg_lt_neg frequency_order))
    · exact Real.exp_le_exp.mpr (neg_le_neg frequency_upper)
    · exact Real.exp_le_exp.mpr (neg_le_neg frequency_lower)
    · exact time_nonneg
    · exact time_le_one
  · positivity

/-- The two-corner residual lies below the upper-corner time-frequency area. -/
theorem fermionicKernel_twoCornerResidual_le_upperCrossArea
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_upper : ω ≤ ωHigh) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ≤
        t * (ωHigh - ω) := by
  let firstResidual := residualUpdate fermionicKernel 0 ωHigh
    (fermionicKernel_pos 0 ωHigh).ne'
  have first_at_lower_nonneg : 0 ≤ firstResidual t ωLow :=
    fermionicKernel_cornerResidual_nonneg time_nonneg frequency_order.le
  have first_at_row_nonneg : 0 ≤ firstResidual 1 ω :=
    fermionicKernel_cornerResidual_nonneg (by norm_num) frequency_upper
  have pivot_pos : 0 < firstResidual 1 ωLow :=
    fermionicKernel_firstCorner_secondPivot_pos frequency_order
  have update_le_first :
      residualUpdate firstResidual 1 ωLow pivot_pos.ne' t ω ≤ firstResidual t ω := by
    unfold residualUpdate
    exact sub_le_self _ (div_nonneg
      (mul_nonneg first_at_lower_nonneg first_at_row_nonneg) pivot_pos.le)
  have first_nonneg : 0 ≤ firstResidual t ω :=
    fermionicKernel_cornerResidual_nonneg time_nonneg frequency_upper
  have first_area := fermionicKernel_cornerResidual_le_area
    time_nonneg time_le_one (show (0 : ℝ) ≤ t from time_nonneg) frequency_upper
  calc
    residualUpdate firstResidual 1 ωLow pivot_pos.ne' t ω ≤ firstResidual t ω := update_le_first
    _ = abs (firstResidual t ω) := (abs_of_nonneg first_nonneg).symm
    _ ≤ (t - 0) * (ωHigh - ω) := first_area
    _ = t * (ωHigh - ω) := by ring

/-- The two-corner residual also lies below the reflected lower-corner cross area. -/
theorem fermionicKernel_twoCornerResidual_le_lowerCrossArea
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : ωLow ≤ ω) (frequency_upper : ω ≤ ωHigh) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ≤
        (1 - t) * (ω - ωLow) := by
  have transformed_order : Real.exp (-ωHigh) < Real.exp (-ωLow) :=
    Real.exp_lt_exp.mpr (neg_lt_neg frequency_order)
  have transformed_lower : Real.exp (-ωHigh) ≤ Real.exp (-ω) :=
    Real.exp_le_exp.mpr (neg_le_neg frequency_upper)
  have transformed_upper : Real.exp (-ω) ≤ Real.exp (-ωLow) :=
    Real.exp_le_exp.mpr (neg_le_neg frequency_lower)
  have secant_lower := mul_rpow_sub_one_le_powerSecant
    (Real.exp_pos (-ωHigh)) transformed_order transformed_lower transformed_upper time_le_one
  have denominator_nonneg : 0 ≤ 1 + Real.exp (-ω) := by positivity
  have residual_le_second :
      residualUpdate
        (residualUpdate fermionicKernel 0 ωHigh
          (fermionicKernel_pos 0 ωHigh).ne')
        1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ≤
          residualUpdate fermionicKernel 1 ωLow
            (fermionicKernel_pos 1 ωLow).ne' t ω := by
    rw [fermionicKernel_twoCornerResidual_eq_powerSecantError frequency_order,
      fermionicKernel_secondCornerResidual_eq_rpow_sub]
    exact div_le_div_of_nonneg_right (sub_le_sub_left secant_lower _) denominator_nonneg
  have second_nonneg := fermionicKernel_residualUpdate_nonneg_of_cross_nonpos
    (mul_nonpos_of_nonpos_of_nonneg (sub_nonpos.mpr time_le_one)
      (sub_nonneg.mpr frequency_lower))
  have second_area := fermionicKernel_residualUpdate_le_crossArea
    time_nonneg time_le_one
      (mul_nonpos_of_nonpos_of_nonneg (sub_nonpos.mpr time_le_one)
        (sub_nonneg.mpr frequency_lower))
  calc
    _ ≤ residualUpdate fermionicKernel 1 ωLow
        (fermionicKernel_pos 1 ωLow).ne' t ω := residual_le_second
    _ = abs (residualUpdate fermionicKernel 1 ωLow
        (fermionicKernel_pos 1 ωLow).ne' t ω) := (abs_of_nonneg second_nonneg).symm
    _ ≤ -((t - 1) * (ω - ωLow)) := second_area
    _ = (1 - t) * (ω - ωLow) := by ring

/--
The exact two-corner residual is a nonnegative tent controlled by its distances
from both transformed corners.
-/
theorem fermionicKernel_twoCornerResidual_mem_Icc_minCrossArea
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : ωLow ≤ ω) (frequency_upper : ω ≤ ωHigh) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ∈
        Set.Icc 0 (min (t * (ωHigh - ω)) ((1 - t) * (ω - ωLow))) := by
  refine ⟨fermionicKernel_twoCornerResidual_nonneg frequency_order time_nonneg time_le_one
    frequency_lower frequency_upper, ?_⟩
  rw [le_min_iff]
  exact ⟨fermionicKernel_twoCornerResidual_le_upperCrossArea frequency_order
      time_nonneg time_le_one frequency_upper,
    fermionicKernel_twoCornerResidual_le_lowerCrossArea frequency_order
      time_nonneg time_le_one frequency_lower frequency_upper⟩

end Fermionic
end GECPKernelStructure
