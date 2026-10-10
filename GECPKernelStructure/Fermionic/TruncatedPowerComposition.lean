import GECPKernelStructure.Fermionic.BetaMomentMinor

namespace GECPKernelStructure
namespace Fermionic

open Set

/-- Seventh-degree truncated-power kernel underlying the eighth-order Taylor tail. -/
noncomputable def eighthOrderTruncatedPowerKernel (v ω : ℝ) : ℝ :=
  max (ω - v) 0 ^ 7

/-- Scaling the beta moment exposes a polynomial-weighted exponential integral. -/
theorem eighthOrderBetaMoment_scaled_eq_truncatedPower (t ω : ℝ) :
    ω ^ 8 * eighthOrderBetaMoment (t * ω) =
      ∫ v in (0 : ℝ)..ω, (ω - v) ^ 7 * Real.exp (-t * v) := by
  let g : ℝ → ℝ := fun v => (ω - v) ^ 7 * Real.exp (-t * v)
  have hscale := intervalIntegral.smul_integral_comp_mul_left
    (a := (0 : ℝ)) (b := 1) g ω
  simp only [smul_eq_mul, mul_one] at hscale
  have hintegrand :
      (∫ u in (0 : ℝ)..1, g (ω * u)) =
        ω ^ 7 * eighthOrderBetaMoment (t * ω) := by
    rw [eighthOrderBetaMoment, ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro u _
    simp only [g]
    have hexp : Real.exp (-t * (ω * u)) = Real.exp (-u * (t * ω)) := by
      congr 1
      ring
    rw [hexp]
    ring
  rw [hintegrand] at hscale
  calc
    ω ^ 8 * eighthOrderBetaMoment (t * ω) =
        ω * (ω ^ 7 * eighthOrderBetaMoment (t * ω)) := by ring
    _ = ∫ v in (0 : ℝ)..ω, (ω - v) ^ 7 * Real.exp (-t * v) := by
      simpa only [mul_zero, g] using hscale

private theorem eighthOrderTruncatedPowerKernel_eq_on_Icc
    {v ω : ℝ} (hv : v ∈ Icc (0 : ℝ) ω) :
    eighthOrderTruncatedPowerKernel v ω = (ω - v) ^ 7 := by
  rw [eighthOrderTruncatedPowerKernel, max_eq_left]
  linarith [hv.2]

private theorem eighthOrderTruncatedPowerKernel_eq_zero_of_le
    {v ω : ℝ} (hv : ω ≤ v) :
    eighthOrderTruncatedPowerKernel v ω = 0 := by
  rw [eighthOrderTruncatedPowerKernel, max_eq_right]
  · norm_num
  · linarith

/-- Above its knot, the truncated-power composition integrand vanishes exactly. -/
theorem eighthOrderTruncatedPower_integral_eq_zero
    {t ω Ω : ℝ} (hωΩ : ω ≤ Ω) :
    (∫ v in ω..Ω,
      eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v)) = 0 := by
  calc
    (∫ v in ω..Ω,
        eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v)) =
        ∫ _v in ω..Ω, (0 : ℝ) := by
      apply intervalIntegral.integral_congr
      intro v hv
      rw [uIcc_of_le hωΩ] at hv
      change eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v) = 0
      rw [eighthOrderTruncatedPowerKernel_eq_zero_of_le hv.1, zero_mul]
    _ = 0 := by simp

/-- Any common cutoff above the knot gives the same truncated-power composition integral. -/
theorem eighthOrderTruncatedPower_integral_commonCutoff
    {t ω Ω : ℝ} (hω : 0 ≤ ω) (hωΩ : ω ≤ Ω) :
    (∫ v in (0 : ℝ)..Ω,
        eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v)) =
      ∫ v in (0 : ℝ)..ω, (ω - v) ^ 7 * Real.exp (-t * v) := by
  let f : ℝ → ℝ := fun v =>
    eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v)
  have hf : Continuous f := by
    unfold f eighthOrderTruncatedPowerKernel
    fun_prop
  have hfleft : IntervalIntegrable f MeasureTheory.volume 0 ω :=
    hf.intervalIntegrable 0 ω
  have hfright : IntervalIntegrable f MeasureTheory.volume ω Ω :=
    hf.intervalIntegrable ω Ω
  have hleft : (∫ v in (0 : ℝ)..ω, f v) =
      ∫ v in (0 : ℝ)..ω, (ω - v) ^ 7 * Real.exp (-t * v) := by
    apply intervalIntegral.integral_congr
    intro v hv
    rw [uIcc_of_le hω] at hv
    simp only [f, eighthOrderTruncatedPowerKernel_eq_on_Icc hv]
  have hadd := intervalIntegral.integral_add_adjacent_intervals
    hfleft hfright
  have hright : (∫ v in ω..Ω, f v) = 0 := by
    simpa only [f] using eighthOrderTruncatedPower_integral_eq_zero (t := t) hωΩ
  rw [hright, add_zero] at hadd
  rw [← hadd, hleft]

/-- The unmasked tail is a positive row scaling of an exponential/truncated-power composition. -/
theorem eighthOrderTail_eq_truncatedPowerComposition (t ω : ℝ) :
    eighthOrderTail (t * ω) =
      t ^ 8 / (Nat.factorial 7 : ℝ) *
        ∫ v in (0 : ℝ)..ω, (ω - v) ^ 7 * Real.exp (-t * v) := by
  rw [eighthOrderTail_eq_betaMoment,
    ← eighthOrderBetaMoment_scaled_eq_truncatedPower t ω]
  ring

/-- A shared cutoff makes the tail representation ready for kernel composition. -/
theorem eighthOrderTail_eq_truncatedPowerComposition_commonCutoff
    {t ω Ω : ℝ} (hω : 0 ≤ ω) (hωΩ : ω ≤ Ω) :
    eighthOrderTail (t * ω) =
      t ^ 8 / (Nat.factorial 7 : ℝ) *
        ∫ v in (0 : ℝ)..Ω,
          eighthOrderTruncatedPowerKernel v ω * Real.exp (-t * v) := by
  rw [eighthOrderTail_eq_truncatedPowerComposition,
    eighthOrderTruncatedPower_integral_commonCutoff hω hωΩ]

end Fermionic
end GECPKernelStructure
