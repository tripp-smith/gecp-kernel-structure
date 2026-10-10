import GECPKernelStructure.Fermionic.TaylorTailBeta
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

namespace GECPKernelStructure
namespace Fermionic

open Set

/-- The first tilted moment of the beta weight from the eighth-order tail. -/
noncomputable def eighthOrderBetaMomentOne (x : ℝ) : ℝ :=
  ∫ u in (0 : ℝ)..1, u * (1 - u) ^ 7 * Real.exp (-u * x)

/-- The second tilted moment of the beta weight from the eighth-order tail. -/
noncomputable def eighthOrderBetaMomentTwo (x : ℝ) : ℝ :=
  ∫ u in (0 : ℝ)..1, u ^ 2 * (1 - u) ^ 7 * Real.exp (-u * x)

/-- An explicit primitive for the unweighted first-minus-twice-second moment. -/
private noncomputable def eighthOrderMomentGapPrimitive (u : ℝ) : ℝ :=
  (1 / 2 : ℝ) * u ^ 2 - 3 * u ^ 3 + (35 / 4 : ℝ) * u ^ 4 -
    (77 / 5 : ℝ) * u ^ 5 + (35 / 2 : ℝ) * u ^ 6 - 13 * u ^ 7 +
    (49 / 8 : ℝ) * u ^ 8 - (5 / 3 : ℝ) * u ^ 9 + (1 / 5 : ℝ) * u ^ 10

private theorem eighthOrderMomentGapPrimitive_hasDerivAt (u : ℝ) :
    HasDerivAt eighthOrderMomentGapPrimitive
      (u * (1 - 2 * u) * (1 - u) ^ 7) u := by
  unfold eighthOrderMomentGapPrimitive
  have h₂ := (hasDerivAt_pow 2 u).const_mul (1 / 2 : ℝ)
  have h₃ := (hasDerivAt_pow 3 u).const_mul (3 : ℝ)
  have h₄ := (hasDerivAt_pow 4 u).const_mul (35 / 4 : ℝ)
  have h₅ := (hasDerivAt_pow 5 u).const_mul (77 / 5 : ℝ)
  have h₆ := (hasDerivAt_pow 6 u).const_mul (35 / 2 : ℝ)
  have h₇ := (hasDerivAt_pow 7 u).const_mul (13 : ℝ)
  have h₈ := (hasDerivAt_pow 8 u).const_mul (49 / 8 : ℝ)
  have h₉ := (hasDerivAt_pow 9 u).const_mul (5 / 3 : ℝ)
  have h₁₀ := (hasDerivAt_pow 10 u).const_mul (1 / 5 : ℝ)
  have derivative := (((((((h₂.sub h₃).add h₄).sub h₅).add h₆).sub h₇).add h₈).sub h₉).add h₁₀
  have derivative' := derivative.congr_deriv
    (g' := u * (1 - 2 * u) * (1 - u) ^ 7) (by
    norm_num
    ring)
  apply derivative'.congr_of_eventuallyEq
  filter_upwards with z
  simp only [Pi.sub_apply, Pi.add_apply]

/-- The unweighted moment gap is exactly `1/120`. -/
theorem eighthOrder_unweighted_moment_gap :
    (∫ u in (0 : ℝ)..1, u * (1 - 2 * u) * (1 - u) ^ 7) = 1 / 120 := by
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt
    (fun u _ => eighthOrderMomentGapPrimitive_hasDerivAt u)
    ((by fun_prop : Continuous
      (fun u : ℝ => u * (1 - 2 * u) * (1 - u) ^ 7)).intervalIntegrable 0 1)]
  norm_num [eighthOrderMomentGapPrimitive]

private theorem eighthOrder_moment_gap_weight_pointwise
    {x u : ℝ} (hx : 0 ≤ x) (hu : u ∈ Icc (0 : ℝ) 1) :
    u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-x / 2) ≤
      u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-u * x) := by
  by_cases hlo : u ≤ 1 / 2
  · have hgap : 0 ≤ u * (1 - 2 * u) * (1 - u) ^ 7 := by
      exact mul_nonneg (mul_nonneg hu.1 (by linarith)) (pow_nonneg (by linarith [hu.2]) 7)
    apply mul_le_mul_of_nonneg_left _ hgap
    rw [Real.exp_le_exp]
    nlinarith [mul_nonneg hx (sub_nonneg.mpr hlo)]
  · have hhalf : 1 / 2 ≤ u := le_of_not_ge hlo
    have hgap : u * (1 - 2 * u) * (1 - u) ^ 7 ≤ 0 := by
      have : 1 - 2 * u ≤ 0 := by linarith
      exact mul_nonpos_of_nonpos_of_nonneg
        (mul_nonpos_of_nonneg_of_nonpos hu.1 this)
        (pow_nonneg (by linarith [hu.2]) 7)
    apply mul_le_mul_of_nonpos_left _ hgap
    rw [Real.exp_le_exp]
    nlinarith [mul_nonneg hx (sub_nonneg.mpr hhalf)]

private theorem eighthOrder_moment_gap_integral_identity (x : ℝ) :
    eighthOrderBetaMomentOne x - 2 * eighthOrderBetaMomentTwo x =
      ∫ u in (0 : ℝ)..1,
        u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-u * x) := by
  have first_integrable : IntervalIntegrable
      (fun u : ℝ => u * (1 - u) ^ 7 * Real.exp (-u * x))
      MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun u : ℝ => u * (1 - u) ^ 7 * Real.exp (-u * x))).intervalIntegrable 0 1
  have second_integrable : IntervalIntegrable
      (fun u : ℝ => u ^ 2 * (1 - u) ^ 7 * Real.exp (-u * x))
      MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous
      (fun u : ℝ => u ^ 2 * (1 - u) ^ 7 * Real.exp (-u * x))).intervalIntegrable 0 1
  rw [eighthOrderBetaMomentOne, eighthOrderBetaMomentTwo,
    ← intervalIntegral.integral_const_mul,
    ← intervalIntegral.integral_sub first_integrable (second_integrable.const_mul 2)]
  apply intervalIntegral.integral_congr
  intro u hu
  ring

/-- Exponential tilting preserves a strict first-versus-second moment gap. -/
theorem eighthOrderBetaMoment_two_mul_two_lt_one {x : ℝ} (hx : 0 ≤ x) :
    2 * eighthOrderBetaMomentTwo x < eighthOrderBetaMomentOne x := by
  have lower_bound :
      Real.exp (-x / 2) * (1 / 120 : ℝ) ≤
        ∫ u in (0 : ℝ)..1,
          u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-u * x) := by
    rw [← eighthOrder_unweighted_moment_gap,
      ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_mono_on (by norm_num)
      ((by fun_prop : Continuous (fun u : ℝ =>
        Real.exp (-x / 2) * (u * (1 - 2 * u) * (1 - u) ^ 7))).intervalIntegrable 0 1)
      ((by fun_prop : Continuous (fun u : ℝ =>
        u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-u * x))).intervalIntegrable 0 1)
    intro u hu
    simpa [mul_assoc, mul_left_comm, mul_comm] using
      eighthOrder_moment_gap_weight_pointwise hx hu
  have gap_pos :
      0 < ∫ u in (0 : ℝ)..1,
        u * (1 - 2 * u) * (1 - u) ^ 7 * Real.exp (-u * x) :=
    lt_of_lt_of_le (mul_pos (Real.exp_pos _) (by norm_num)) lower_bound
  rw [← eighthOrder_moment_gap_integral_identity] at gap_pos
  linarith

private theorem eighthOrderBetaMomentOne_pos (x : ℝ) :
    0 < eighthOrderBetaMomentOne x := by
  rw [eighthOrderBetaMomentOne]
  apply intervalIntegral.intervalIntegral_pos_of_pos_on
    ((by fun_prop : Continuous (fun u : ℝ =>
      u * (1 - u) ^ 7 * Real.exp (-u * x))).intervalIntegrable 0 1) _ (by norm_num)
  intro u hu
  exact mul_pos (mul_pos hu.1 (pow_pos (by linarith [hu.2]) 7)) (Real.exp_pos _)

private theorem eighthOrderBetaMomentTwo_nonneg (x : ℝ) :
    0 ≤ eighthOrderBetaMomentTwo x := by
  rw [eighthOrderBetaMomentTwo]
  apply intervalIntegral.integral_nonneg (by norm_num)
  intro u hu
  exact mul_nonneg (mul_nonneg (sq_nonneg u) (pow_nonneg (by linarith [hu.2]) 7))
    (Real.exp_pos _).le

/-- The tilted-beta variance numerator has the strict local bound needed by logarithmic elasticity. -/
theorem eighthOrderBetaMoment_variance_numerator_lt
    {x : ℝ} (hx0 : 0 ≤ x) (hx2 : x ≤ 2) :
    x * (eighthOrderBetaMomentTwo x * eighthOrderBetaMoment x -
        eighthOrderBetaMomentOne x ^ 2) <
      eighthOrderBetaMomentOne x * eighthOrderBetaMoment x := by
  have moment_gap := eighthOrderBetaMoment_two_mul_two_lt_one hx0
  have moment0_pos := eighthOrderBetaMoment_pos x
  have moment1_pos := eighthOrderBetaMomentOne_pos x
  have moment2_nonneg := eighthOrderBetaMomentTwo_nonneg x
  have x_moment2_lt : x * eighthOrderBetaMomentTwo x < eighthOrderBetaMomentOne x := by
    calc
      x * eighthOrderBetaMomentTwo x ≤ 2 * eighthOrderBetaMomentTwo x := by
        exact mul_le_mul_of_nonneg_right hx2 moment2_nonneg
      _ < eighthOrderBetaMomentOne x := moment_gap
  have scaled_lt :
      x * eighthOrderBetaMomentTwo x * eighthOrderBetaMoment x <
        eighthOrderBetaMomentOne x * eighthOrderBetaMoment x :=
    mul_lt_mul_of_pos_right x_moment2_lt moment0_pos
  have square_nonneg : 0 ≤ x * eighthOrderBetaMomentOne x ^ 2 := by positivity
  nlinarith

end Fermionic
end GECPKernelStructure
