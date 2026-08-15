import GECPKernelStructure.Fermionic.Centering
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

namespace GECPKernelStructure
namespace Fermionic

/-- The fermionic kernel is jointly continuous on all real coordinates. -/
theorem fermionicKernel_continuous :
    Continuous (fun point : ℝ × ℝ ↦ fermionicKernel point.1 point.2) := by
  unfold fermionicKernel
  apply Continuous.div
  · fun_prop
  · fun_prop
  · intro point
    exact fermionicKernel_denominator_ne point.2

/-- The logistic factor occurring in the frequency derivative. -/
noncomputable def fermionicLogistic (ω : ℝ) : ℝ :=
  Real.exp (-ω) / (1 + Real.exp (-ω))

theorem fermionicKernel_le_one_of_nonneg {t ω : ℝ} (ht : 0 ≤ t) (hω : 0 ≤ ω) :
    fermionicKernel t ω ≤ 1 := by
  unfold fermionicKernel
  rw [div_le_one (fermionicKernel_denominator_pos ω)]
  simpa only [neg_mul] using
    (Real.exp_le_one_iff.mpr (neg_nonpos.mpr (mul_nonneg ht hω))).trans
      (le_add_of_nonneg_right (le_of_lt (Real.exp_pos (-ω))))

/-- The fermionic kernel is bounded by one on the physical time interval. -/
theorem fermionicKernel_le_one {t ω : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    fermionicKernel t ω ≤ 1 := by
  by_cases hω : 0 ≤ ω
  · exact fermionicKernel_le_one_of_nonneg ht0 hω
  · rw [← fermionicKernel_reflection t ω]
    exact fermionicKernel_le_one_of_nonneg (sub_nonneg.mpr ht1)
      (neg_nonneg.mpr (le_of_not_ge hω))

/-- Exact derivative with respect to imaginary time. -/
theorem fermionicKernel_hasDerivAt_time (t ω : ℝ) :
    HasDerivAt (fun s : ℝ => fermionicKernel s ω)
      (fermionicKernel t ω * (-ω)) t := by
  have h := ((hasDerivAt_id t).mul_const (-ω)).exp.div_const
    (1 + Real.exp (-ω))
  have hf : (fun x : ℝ => Real.exp (id x * (-ω)) / (1 + Real.exp (-ω))) =
      (fun s : ℝ => fermionicKernel s ω) := by
    funext s
    unfold fermionicKernel
    congr 2
    simp only [id_eq]
    ring
  have hd : Real.exp (id t * (-ω)) * (1 * (-ω)) / (1 + Real.exp (-ω)) =
      fermionicKernel t ω * (-ω) := by
    unfold fermionicKernel
    have hexp : Real.exp (id t * (-ω)) = Real.exp (-t * ω) := by
      congr 1
      simp only [id_eq]
      ring
    rw [hexp]
    ring
  rw [hf, hd] at h
  exact h

/-- Exact derivative with respect to frequency. -/
theorem fermionicKernel_hasDerivAt_frequency (t ω : ℝ) :
    HasDerivAt (fun r : ℝ => fermionicKernel t r)
      (fermionicKernel t ω * (fermionicLogistic ω - t)) ω := by
  have hnum := ((hasDerivAt_id ω).mul_const (-t)).exp
  have hden := ((hasDerivAt_id ω).mul_const (-1 : ℝ)).exp.const_add 1
  have hden_ne : 1 + Real.exp (id ω * (-1 : ℝ)) ≠ 0 := by
    simpa using fermionicKernel_denominator_ne ω
  have h := hnum.fun_div hden hden_ne
  have hf : (fun x : ℝ => Real.exp (id x * (-t)) /
      (1 + Real.exp (id x * (-1)))) =
      (fun r : ℝ => fermionicKernel t r) := by
    funext r
    unfold fermionicKernel
    have hnumarg : id r * (-t) = -t * r := by simp [mul_comm]
    have hdenarg : id r * (-1 : ℝ) = -r := by simp
    rw [hnumarg, hdenarg]
  have hd :
      (Real.exp (id ω * (-t)) * (1 * (-t)) *
          (1 + Real.exp (id ω * (-1))) -
        Real.exp (id ω * (-t)) * (Real.exp (id ω * (-1)) * (1 * (-1)))) /
          (1 + Real.exp (id ω * (-1))) ^ 2 =
        fermionicKernel t ω * (fermionicLogistic ω - t) := by
    have hnumarg : id ω * (-t) = -t * ω := by simp [mul_comm]
    have hdenarg : id ω * (-1 : ℝ) = -ω := by simp
    rw [hnumarg, hdenarg]
    unfold fermionicKernel fermionicLogistic
    field_simp [fermionicKernel_denominator_ne ω]
    ring
  rw [hf, hd] at h
  exact h

theorem fermionicLogistic_pos (ω : ℝ) : 0 < fermionicLogistic ω := by
  unfold fermionicLogistic
  positivity

theorem fermionicLogistic_lt_one (ω : ℝ) : fermionicLogistic ω < 1 := by
  unfold fermionicLogistic
  rw [div_lt_one (fermionicKernel_denominator_pos ω)]
  linarith [Real.exp_pos (-ω)]

/-- The frequency logarithmic slope has magnitude at most one. -/
theorem fermionicKernel_frequencySlope_bound {t ω : ℝ}
    (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |fermionicLogistic ω - t| ≤ 1 := by
  rw [abs_le]
  constructor <;> linarith [fermionicLogistic_pos ω, fermionicLogistic_lt_one ω]

/-- Uniform frequency-derivative bound on the physical time interval. -/
theorem fermionicKernel_frequencyDerivative_bound {t ω : ℝ}
    (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |fermionicKernel t ω * (fermionicLogistic ω - t)| ≤ 1 := by
  rw [abs_mul, abs_of_pos (fermionicKernel_pos t ω)]
  calc
    fermionicKernel t ω * |fermionicLogistic ω - t| ≤ 1 * 1 :=
      mul_le_mul (fermionicKernel_le_one ht0 ht1)
        (fermionicKernel_frequencySlope_bound ht0 ht1)
        (abs_nonneg _) zero_le_one
    _ = 1 := one_mul 1

/-- Uniform time-derivative bound on a frequency interval of radius `cutoff`. -/
theorem fermionicKernel_timeDerivative_bound {t ω cutoff : ℝ}
    (ht0 : 0 ≤ t) (ht1 : t ≤ 1) (hω : |ω| ≤ cutoff) :
    |fermionicKernel t ω * (-ω)| ≤ cutoff := by
  rw [abs_mul, abs_neg, abs_of_pos (fermionicKernel_pos t ω)]
  calc
    fermionicKernel t ω * |ω| ≤ 1 * cutoff :=
      mul_le_mul (fermionicKernel_le_one ht0 ht1) hω (abs_nonneg _) zero_le_one
    _ = cutoff := one_mul cutoff

end Fermionic
end GECPKernelStructure
