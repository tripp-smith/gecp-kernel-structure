import GECPKernelStructure.Fermionic.BetaMomentElasticity

namespace GECPKernelStructure
namespace Fermionic

open Set

/-- Ratio of the normalized beta moment at two time scales. -/
noncomputable def eighthOrderBetaMomentScaleRatio (t₀ t₁ ω : ℝ) : ℝ :=
  eighthOrderBetaMoment (t₀ * ω) / eighthOrderBetaMoment (t₁ * ω)

/-- Exact derivative of the beta-moment scale ratio. -/
theorem eighthOrderBetaMomentScaleRatio_hasDerivAt (t₀ t₁ ω : ℝ) :
    HasDerivAt (eighthOrderBetaMomentScaleRatio t₀ t₁)
      (((-t₀ * eighthOrderBetaMomentOne (t₀ * ω)) *
            eighthOrderBetaMoment (t₁ * ω) -
          eighthOrderBetaMoment (t₀ * ω) *
            (-t₁ * eighthOrderBetaMomentOne (t₁ * ω))) /
        eighthOrderBetaMoment (t₁ * ω) ^ 2) ω := by
  have hlinear₀ : HasDerivAt (fun x : ℝ => t₀ * x) t₀ ω := by
    simpa only [id_eq, mul_one] using (hasDerivAt_id ω).const_mul t₀
  have hlinear₁ : HasDerivAt (fun x : ℝ => t₁ * x) t₁ ω := by
    simpa only [id_eq, mul_one] using (hasDerivAt_id ω).const_mul t₁
  have hmoment₀raw := (eighthOrderBetaMoment_hasDerivAt (t₀ * ω)).comp ω hlinear₀
  have hmoment₁raw := (eighthOrderBetaMoment_hasDerivAt (t₁ * ω)).comp ω hlinear₁
  have hmoment₀ : HasDerivAt (fun x => eighthOrderBetaMoment (t₀ * x))
      (-t₀ * eighthOrderBetaMomentOne (t₀ * ω)) ω := by
    have h := hmoment₀raw.congr_deriv
      (g' := -t₀ * eighthOrderBetaMomentOne (t₀ * ω)) (by ring)
    exact h.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun x => by
      simp only [Function.comp_apply]
  have hmoment₁ : HasDerivAt (fun x => eighthOrderBetaMoment (t₁ * x))
      (-t₁ * eighthOrderBetaMomentOne (t₁ * ω)) ω := by
    have h := hmoment₁raw.congr_deriv
      (g' := -t₁ * eighthOrderBetaMomentOne (t₁ * ω)) (by ring)
    exact h.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun x => by
      simp only [Function.comp_apply]
  have denominator_ne : eighthOrderBetaMoment (t₁ * ω) ≠ 0 :=
    ne_of_gt (eighthOrderBetaMoment_pos _)
  exact hmoment₀.div hmoment₁ denominator_ne

private theorem eighthOrderBetaMomentScaleRatio_deriv_pos
    {t₀ t₁ ω : ℝ} (ht₀ : 0 < t₀) (ht : t₀ < t₁)
    (hω : 0 < ω) (hupper : t₁ * ω ≤ 2) :
    0 < ((-t₀ * eighthOrderBetaMomentOne (t₀ * ω)) *
            eighthOrderBetaMoment (t₁ * ω) -
          eighthOrderBetaMoment (t₀ * ω) *
            (-t₁ * eighthOrderBetaMomentOne (t₁ * ω))) /
        eighthOrderBetaMoment (t₁ * ω) ^ 2 := by
  have ht₁ : 0 < t₁ := lt_trans ht₀ ht
  have hproduct : t₀ * ω < t₁ * ω := mul_lt_mul_of_pos_right ht hω
  have hsmall_mem : t₀ * ω ∈ Icc (0 : ℝ) 2 :=
    ⟨(mul_pos ht₀ hω).le, le_trans hproduct.le hupper⟩
  have hlarge_mem : t₁ * ω ∈ Icc (0 : ℝ) 2 :=
    ⟨(mul_pos ht₁ hω).le, hupper⟩
  have helasticity := eighthOrderBetaElasticity_strictAntiOn
    hsmall_mem hlarge_mem hproduct
  have hmoment₀ := eighthOrderBetaMoment_pos (t₀ * ω)
  have hmoment₁ := eighthOrderBetaMoment_pos (t₁ * ω)
  unfold eighthOrderBetaElasticity at helasticity
  rw [div_lt_div_iff₀ hmoment₁ hmoment₀] at helasticity
  have numerator_pos :
      0 < (-t₀ * eighthOrderBetaMomentOne (t₀ * ω)) *
            eighthOrderBetaMoment (t₁ * ω) -
          eighthOrderBetaMoment (t₀ * ω) *
            (-t₁ * eighthOrderBetaMomentOne (t₁ * ω)) := by
    ring_nf at helasticity ⊢
    nlinarith
  exact div_pos numerator_pos (sq_pos_of_pos hmoment₁)

/-- For ordered positive time scales, the beta-moment scale ratio strictly increases locally. -/
theorem eighthOrderBetaMomentScaleRatio_strictMonoOn
    {t₀ t₁ ω₀ ω₁ : ℝ} (ht₀ : 0 < t₀) (ht : t₀ < t₁)
    (hω₀ : 0 < ω₀) (hupper : t₁ * ω₁ ≤ 2) :
    StrictMonoOn (eighthOrderBetaMomentScaleRatio t₀ t₁) (Icc ω₀ ω₁) := by
  apply strictMonoOn_of_deriv_pos (convex_Icc ω₀ ω₁)
  · exact fun x _ =>
      (eighthOrderBetaMomentScaleRatio_hasDerivAt t₀ t₁ x).continuousAt.continuousWithinAt
  · intro x hx
    rw [interior_Icc] at hx
    rw [(eighthOrderBetaMomentScaleRatio_hasDerivAt t₀ t₁ x).deriv]
    apply eighthOrderBetaMomentScaleRatio_deriv_pos ht₀ ht
    · exact lt_trans hω₀ hx.1
    · have ht₁ : 0 < t₁ := lt_trans ht₀ ht
      exact le_trans (mul_le_mul_of_nonneg_left hx.2.le ht₁.le) hupper

/-- Every positive product-bounded ordered beta-moment minor has the expected negative sign. -/
theorem eighthOrderBetaMoment_fin_two_det_neg
    {t₀ t₁ ω₀ ω₁ : ℝ} (ht₀ : 0 < t₀) (ht : t₀ < t₁)
    (hω₀ : 0 < ω₀) (hω : ω₀ < ω₁) (hupper : t₁ * ω₁ ≤ 2) :
    eighthOrderBetaMoment (t₀ * ω₀) * eighthOrderBetaMoment (t₁ * ω₁) -
        eighthOrderBetaMoment (t₀ * ω₁) * eighthOrderBetaMoment (t₁ * ω₀) < 0 := by
  have hratio := eighthOrderBetaMomentScaleRatio_strictMonoOn
    ht₀ ht hω₀ hupper (left_mem_Icc.mpr hω.le) (right_mem_Icc.mpr hω.le) hω
  unfold eighthOrderBetaMomentScaleRatio at hratio
  rw [div_lt_div_iff₀ (eighthOrderBetaMoment_pos (t₁ * ω₀))
    (eighthOrderBetaMoment_pos (t₁ * ω₁))] at hratio
  linarith

/-- The unmasked eighth-order Taylor tail has the local order-two negative sign. -/
theorem eighthOrderTail_fin_two_det_neg
    {t₀ t₁ ω₀ ω₁ : ℝ} (ht₀ : 0 < t₀) (ht : t₀ < t₁)
    (hω₀ : 0 < ω₀) (hω : ω₀ < ω₁) (hupper : t₁ * ω₁ ≤ 2) :
    eighthOrderTail (t₀ * ω₀) * eighthOrderTail (t₁ * ω₁) -
        eighthOrderTail (t₀ * ω₁) * eighthOrderTail (t₁ * ω₀) < 0 := by
  have ht₁ : 0 < t₁ := lt_trans ht₀ ht
  have hω₁ : 0 < ω₁ := lt_trans hω₀ hω
  rw [eighthOrderTail_fin_two_det_neg_iff_betaMoment ht₀ ht₁ hω₀ hω₁]
  exact eighthOrderBetaMoment_fin_two_det_neg ht₀ ht hω₀ hω hupper

end Fermionic
end GECPKernelStructure
