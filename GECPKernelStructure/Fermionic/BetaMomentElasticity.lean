import GECPKernelStructure.Fermionic.BetaMomentGap
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.Calculus.ParametricIntervalIntegral

namespace GECPKernelStructure
namespace Fermionic

open Set
open Filter
open scoped Interval Topology

/-- The beta-weighted power moments used to differentiate the normalized tail. -/
private noncomputable def eighthOrderBetaPowerMoment (k : ℕ) (x : ℝ) : ℝ :=
  ∫ u in (0 : ℝ)..1, u ^ k * (1 - u) ^ 7 * Real.exp (-u * x)

private theorem eighthOrderBetaPowerMoment_hasDerivAt (k : ℕ) (x₀ : ℝ) :
    HasDerivAt (eighthOrderBetaPowerMoment k)
      (-eighthOrderBetaPowerMoment (k + 1) x₀) x₀ := by
  let F : ℝ → ℝ → ℝ := fun x u =>
    u ^ k * (1 - u) ^ 7 * Real.exp (-u * x)
  let F' : ℝ → ℝ → ℝ := fun x u =>
    -u ^ (k + 1) * (1 - u) ^ 7 * Real.exp (-u * x)
  let bound : ℝ → ℝ := fun _ => Real.exp (|x₀| + 1)
  have hs : Metric.ball x₀ 1 ∈ 𝓝 x₀ := Metric.ball_mem_nhds x₀ zero_lt_one
  have hF_meas : ∀ᶠ x in 𝓝 x₀,
      MeasureTheory.AEStronglyMeasurable (F x)
        (MeasureTheory.volume.restrict (Ι (0 : ℝ) 1)) := by
    filter_upwards with x
    exact (by fun_prop : Continuous (F x)).aestronglyMeasurable
  have hF_int : IntervalIntegrable (F x₀) MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous (F x₀)).intervalIntegrable 0 1
  have hF'_meas : MeasureTheory.AEStronglyMeasurable (F' x₀)
      (MeasureTheory.volume.restrict (Ι (0 : ℝ) 1)) :=
    (by fun_prop : Continuous (F' x₀)).aestronglyMeasurable
  have h_bound : ∀ᵐ u ∂MeasureTheory.volume, u ∈ Ι (0 : ℝ) 1 →
      ∀ x ∈ Metric.ball x₀ 1, ‖F' x u‖ ≤ bound u := by
    filter_upwards with u
    intro hu x hx
    rw [uIoc_of_le (by norm_num)] at hu
    have hu0 : 0 ≤ u := le_of_lt hu.1
    have hu1 : u ≤ 1 := hu.2
    have hux : -u * x ≤ |x₀| + 1 := by
      have hdist : |x - x₀| < 1 := by simpa [Real.dist_eq] using hx
      have habsx : |x| < |x₀| + 1 := by
        calc
          |x| = |x - x₀ + x₀| := by congr 1; ring
          _ ≤ |x - x₀| + |x₀| := abs_add_le _ _
          _ < |x₀| + 1 := by linarith
      calc
        -u * x ≤ |-u * x| := by
          exact le_abs_self _
        _ = u * |x| := by
          rw [abs_mul, abs_neg, abs_of_nonneg hu0]
        _ ≤ |x| := mul_le_of_le_one_left (abs_nonneg x) hu1
        _ ≤ |x₀| + 1 := habsx.le
    have hexp : Real.exp (-u * x) ≤ Real.exp (|x₀| + 1) :=
      Real.exp_le_exp.mpr hux
    have hu_pow : u ^ (k + 1) ≤ 1 := pow_le_one₀ hu0 hu1
    have hone_pow : (1 - u) ^ 7 ≤ 1 := by
      exact pow_le_one₀ (by linarith) (by linarith)
    simp only [F', bound, Real.norm_eq_abs, abs_mul, abs_neg, Real.abs_exp,
      abs_pow, abs_of_nonneg hu0, abs_of_nonneg (sub_nonneg.mpr hu1)]
    calc
      u ^ (k + 1) * (1 - u) ^ 7 * Real.exp (-u * x)
          ≤ 1 * 1 * Real.exp (|x₀| + 1) := by gcongr
      _ = Real.exp (|x₀| + 1) := by ring
  have bound_integrable : IntervalIntegrable bound MeasureTheory.volume 0 1 :=
    (by fun_prop : Continuous bound).intervalIntegrable 0 1
  have h_diff : ∀ᵐ u ∂MeasureTheory.volume, u ∈ Ι (0 : ℝ) 1 →
      ∀ x ∈ Metric.ball x₀ 1, HasDerivAt (fun x => F x u) (F' x u) x := by
    filter_upwards with u
    intro _ x _
    have hlin : HasDerivAt (fun x : ℝ => -u * x) (-u) x := by
      simpa only [id_eq, mul_one] using (hasDerivAt_id x).const_mul (-u)
    have hexp : HasDerivAt (fun x : ℝ => Real.exp (-u * x))
        (-u * Real.exp (-u * x)) x := by
      have hraw := (Real.hasDerivAt_exp (-u * x)).comp x hlin
      have hraw' := hraw.congr_deriv
        (g' := -u * Real.exp (-u * x)) (by ring)
      exact hraw'.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
        simp only [Function.comp_apply]
    have h := hexp.const_mul (u ^ k * (1 - u) ^ 7)
    have h' : HasDerivAt (fun x => F x u) (F' x u) x := by
      apply h.congr_deriv
      simp only [F', pow_succ]
      ring
    exact h'.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
      simp only [F]
  have differentiated :=
    intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
      hs hF_meas hF_int hF'_meas h_bound bound_integrable h_diff
  have hderiv : (∫ u in (0 : ℝ)..1, F' x₀ u) =
      -eighthOrderBetaPowerMoment (k + 1) x₀ := by
    rw [eighthOrderBetaPowerMoment, ← intervalIntegral.integral_neg]
    apply intervalIntegral.integral_congr
    intro u _
    simp only [F']
    ring
  have h := differentiated.2.congr_deriv hderiv
  exact h.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
    rfl

/-- The derivative of the normalized beta moment is minus its first tilted moment. -/
theorem eighthOrderBetaMoment_hasDerivAt (x : ℝ) :
    HasDerivAt eighthOrderBetaMoment (-eighthOrderBetaMomentOne x) x := by
  have h := eighthOrderBetaPowerMoment_hasDerivAt 0 x
  have hderiv : -eighthOrderBetaPowerMoment 1 x = -eighthOrderBetaMomentOne x := by
    simp [eighthOrderBetaPowerMoment, eighthOrderBetaMomentOne]
  have h' := h.congr_deriv hderiv
  exact h'.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
    simp [eighthOrderBetaPowerMoment, eighthOrderBetaMoment]

/-- The derivative of the first tilted moment is minus the second tilted moment. -/
theorem eighthOrderBetaMomentOne_hasDerivAt (x : ℝ) :
    HasDerivAt eighthOrderBetaMomentOne (-eighthOrderBetaMomentTwo x) x := by
  have h := eighthOrderBetaPowerMoment_hasDerivAt 1 x
  have hderiv : -eighthOrderBetaPowerMoment 2 x = -eighthOrderBetaMomentTwo x := by
    simp [eighthOrderBetaPowerMoment, eighthOrderBetaMomentTwo]
  have h' := h.congr_deriv hderiv
  exact h'.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
    simp [eighthOrderBetaPowerMoment, eighthOrderBetaMomentOne]

/-- Logarithmic elasticity of the normalized beta moment. -/
noncomputable def eighthOrderBetaElasticity (x : ℝ) : ℝ :=
  -x * eighthOrderBetaMomentOne x / eighthOrderBetaMoment x

/-- Exact derivative of the normalized beta moment's logarithmic elasticity. -/
theorem eighthOrderBetaElasticity_hasDerivAt (x : ℝ) :
    HasDerivAt eighthOrderBetaElasticity
      ((x * (eighthOrderBetaMomentTwo x * eighthOrderBetaMoment x -
          eighthOrderBetaMomentOne x ^ 2) -
        eighthOrderBetaMomentOne x * eighthOrderBetaMoment x) /
        eighthOrderBetaMoment x ^ 2) x := by
  have h0 := eighthOrderBetaMoment_hasDerivAt x
  have h1 := eighthOrderBetaMomentOne_hasDerivAt x
  have h0_ne : eighthOrderBetaMoment x ≠ 0 := ne_of_gt (eighthOrderBetaMoment_pos x)
  have h := (((hasDerivAt_id x).neg.mul h1).div h0 h0_ne)
  have hderiv :
      ((-1 * eighthOrderBetaMomentOne x + (-x) * (-eighthOrderBetaMomentTwo x)) *
          eighthOrderBetaMoment x -
        ((-x) * eighthOrderBetaMomentOne x) * (-eighthOrderBetaMomentOne x)) /
          eighthOrderBetaMoment x ^ 2 =
        (x * (eighthOrderBetaMomentTwo x * eighthOrderBetaMoment x -
            eighthOrderBetaMomentOne x ^ 2) -
          eighthOrderBetaMomentOne x * eighthOrderBetaMoment x) /
          eighthOrderBetaMoment x ^ 2 := by
    ring
  have h' := h.congr_deriv hderiv
  exact h'.congr_of_eventuallyEq <| Filter.Eventually.of_forall fun y => by
    rfl

/-- The beta-moment logarithmic elasticity strictly decreases on the local Taylor regime. -/
theorem eighthOrderBetaElasticity_strictAntiOn :
    StrictAntiOn eighthOrderBetaElasticity (Icc (0 : ℝ) 2) := by
  apply strictAntiOn_of_deriv_neg (convex_Icc 0 2)
  · exact fun x _ => (eighthOrderBetaElasticity_hasDerivAt x).continuousAt.continuousWithinAt
  · intro x hx
    rw [interior_Icc] at hx
    rw [(eighthOrderBetaElasticity_hasDerivAt x).deriv]
    apply div_neg_of_neg_of_pos
    · exact sub_neg.mpr
        (eighthOrderBetaMoment_variance_numerator_lt hx.1.le hx.2.le)
    · exact sq_pos_of_pos (eighthOrderBetaMoment_pos x)

end Fermionic
end GECPKernelStructure
