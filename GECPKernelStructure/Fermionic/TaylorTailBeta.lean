import GECPKernelStructure.Fermionic.ExponentialFamily
import Mathlib.Analysis.Calculus.TaylorIntegral

namespace GECPKernelStructure
namespace Fermionic

open scoped BigOperators
open Finset Set

/-- The unmasked error after the first eight terms of `exp (-x)`. -/
noncomputable def eighthOrderTail (x : ℝ) : ℝ :=
  Real.exp (-x) - expNegTaylor 8 x

/-- The beta-weighted Laplace moment that remains after removing `x^8 / 7!`. -/
noncomputable def eighthOrderBetaMoment (x : ℝ) : ℝ :=
  ∫ u in (0 : ℝ)..1, (1 - u) ^ 7 * Real.exp (-u * x)

private theorem eighthOrder_iteratedFDeriv_apply (x u : ℝ) :
    iteratedFDeriv ℝ 8 (fun y : ℝ => Real.exp (-y)) (u * x) (fun _ => x) =
      x ^ 8 * Real.exp (-(u * x)) := by
  rw [iteratedFDeriv_apply_eq_iteratedDeriv_mul_prod]
  rw [show (fun y : ℝ => Real.exp (-y)) =
    (fun y : ℝ => Real.exp ((-1) * y)) by
      funext y
      congr 1
      ring]
  rw [congrFun (iteratedDeriv_exp_const_mul 8 (-1)) (u * x)]
  simp only [Finset.prod_const, Finset.card_fin, smul_eq_mul]
  norm_num

private theorem eighthOrder_taylor_term (x : ℝ) (k : ℕ) :
    ((k.factorial : ℝ)⁻¹ •
      iteratedFDeriv ℝ k (fun y : ℝ => Real.exp (-y)) 0 (fun _ => x)) =
        (-x) ^ k / (k.factorial : ℝ) := by
  rw [iteratedFDeriv_apply_eq_iteratedDeriv_mul_prod]
  rw [show (fun y : ℝ => Real.exp (-y)) =
    (fun y : ℝ => Real.exp ((-1) * y)) by
      funext y
      congr 1
      ring]
  rw [congrFun (iteratedDeriv_exp_const_mul k (-1)) 0]
  simp only [Finset.prod_const, Finset.card_fin, smul_eq_mul, mul_zero]
  rw [Real.exp_zero, mul_one]
  rw [show (-x) ^ k = (-1 : ℝ) ^ k * x ^ k by ring]
  field_simp

/-- Taylor's integral remainder exposes the eighth-order tail as a beta moment. -/
theorem eighthOrderTail_eq_betaMoment (x : ℝ) :
    eighthOrderTail x = x ^ 8 / (Nat.factorial 7 : ℝ) * eighthOrderBetaMoment x := by
  have taylor := map_add_eq_sum_add_integral_iteratedFDeriv
    (f := fun y : ℝ => Real.exp (-y)) (x := 0) (y := x) (n := 7)
    (by intro u hu; fun_prop)
  simp_rw [eighthOrder_taylor_term] at taylor
  rw [show (∑ k ∈ Finset.range 8, (-x) ^ k / (k.factorial : ℝ)) =
    expNegTaylor 8 x by rfl] at taylor
  simp only [zero_add, smul_eq_mul] at taylor
  norm_num at taylor
  have derivative_rewrite (u : ℝ) :
      iteratedFDeriv ℝ (7 + 1) (fun y : ℝ => Real.exp (-y)) (u * x) (fun _ => x) =
        x ^ 8 * Real.exp (-(u * x)) := by
    simpa only using eighthOrder_iteratedFDeriv_apply x u
  simp_rw [derivative_rewrite] at taylor
  have tail_integral :
      eighthOrderTail x =
        (1 / 5040 : ℝ) * ∫ u in (0 : ℝ)..1,
          (1 - u) ^ 7 * (x ^ 8 * Real.exp (-(u * x))) := by
    rw [eighthOrderTail, taylor]
    ring
  rw [tail_integral, eighthOrderBetaMoment]
  norm_num
  conv_lhs => rw [← intervalIntegral.integral_const_mul]
  conv_rhs => rw [← intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_congr
  intro u hu
  ring

/-- The normalized beta moment is strictly positive for every real argument. -/
theorem eighthOrderBetaMoment_pos (x : ℝ) :
    0 < eighthOrderBetaMoment x := by
  rw [eighthOrderBetaMoment]
  have integrable : IntervalIntegrable
      (fun u : ℝ => (1 - u) ^ 7 * Real.exp (-u * x)) MeasureTheory.volume 0 1 := by
    exact (by fun_prop : Continuous
      (fun u : ℝ => (1 - u) ^ 7 * Real.exp (-u * x))).intervalIntegrable 0 1
  apply intervalIntegral.intervalIntegral_pos_of_pos_on
    integrable _ (by norm_num)
  intro u hu
  have hu1 : 0 < 1 - u := by linarith [hu.2]
  exact mul_pos (pow_pos hu1 7) (Real.exp_pos _)

/-- The unmasked eighth-order tail is strictly positive on positive arguments. -/
theorem eighthOrderTail_pos {x : ℝ} (hx : 0 < x) :
    0 < eighthOrderTail x := by
  rw [eighthOrderTail_eq_betaMoment]
  exact mul_pos (div_pos (pow_pos hx 8) (by positivity))
    (eighthOrderBetaMoment_pos x)

/--
Every order-two tail minor is a positive monomial scaling of the corresponding
beta-moment minor.
-/
theorem eighthOrderTail_fin_two_det_factor (t₀ t₁ ω₀ ω₁ : ℝ) :
    eighthOrderTail (t₀ * ω₀) * eighthOrderTail (t₁ * ω₁) -
        eighthOrderTail (t₀ * ω₁) * eighthOrderTail (t₁ * ω₀) =
      (t₀ * t₁ * ω₀ * ω₁) ^ 8 / (Nat.factorial 7 : ℝ) ^ 2 *
        (eighthOrderBetaMoment (t₀ * ω₀) * eighthOrderBetaMoment (t₁ * ω₁) -
          eighthOrderBetaMoment (t₀ * ω₁) * eighthOrderBetaMoment (t₁ * ω₀)) := by
  rw [eighthOrderTail_eq_betaMoment, eighthOrderTail_eq_betaMoment,
    eighthOrderTail_eq_betaMoment, eighthOrderTail_eq_betaMoment]
  ring

/-- Positive row and column arguments leave the order-two sign entirely in the beta moment. -/
theorem eighthOrderTail_fin_two_det_neg_iff_betaMoment
    {t₀ t₁ ω₀ ω₁ : ℝ} (ht₀ : 0 < t₀) (ht₁ : 0 < t₁)
    (hω₀ : 0 < ω₀) (hω₁ : 0 < ω₁) :
    eighthOrderTail (t₀ * ω₀) * eighthOrderTail (t₁ * ω₁) -
          eighthOrderTail (t₀ * ω₁) * eighthOrderTail (t₁ * ω₀) < 0 ↔
      eighthOrderBetaMoment (t₀ * ω₀) * eighthOrderBetaMoment (t₁ * ω₁) -
          eighthOrderBetaMoment (t₀ * ω₁) * eighthOrderBetaMoment (t₁ * ω₀) < 0 := by
  rw [eighthOrderTail_fin_two_det_factor]
  have coefficient_pos :
      0 < (t₀ * t₁ * ω₀ * ω₁) ^ 8 / (Nat.factorial 7 : ℝ) ^ 2 := by
    positivity
  constructor
  · intro product_neg
    rcases mul_neg_iff.mp product_neg with ⟨_, hmoment⟩ | ⟨hcoefficient, _⟩
    · exact hmoment
    · linarith
  · exact mul_neg_of_pos_of_neg coefficient_pos

end Fermionic
end GECPKernelStructure
