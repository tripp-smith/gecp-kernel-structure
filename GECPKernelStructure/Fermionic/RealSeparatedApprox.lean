import GECPKernelStructure.Fermionic.RealCutoffRate
import GECPKernelStructure.Fermionic.SeparatedApprox

namespace GECPKernelStructure.Fermionic

/-- D1 with its stated real cutoff/tolerance quantifiers and one explicit
universal constant. The witness is the existing dyadic Taylor construction. -/
theorem fermionicKernel_separatedApprox_real (Λ epsilon : ℝ)
    (hΛ : 1 ≤ Λ) (hepsilon : 0 < epsilon) (hepsilon_one : epsilon < 1) :
    ∃ terms : List SeparatedTerm,
      (terms.length : ℝ) ≤ 128 * Real.log (1 + Λ) * Real.log (128 / epsilon) ∧
      ∀ t omega : ℝ, 0 ≤ t → t ≤ 1 → -Λ ≤ omega → omega ≤ Λ →
        |fermionicKernel t omega - evalSeparated terms t omega| ≤ epsilon := by
  obtain ⟨s, hs, hS⟩ := exists_logarithmic_cutoff_cover hΛ
  have ha : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hratio : 0 < Real.log (1 / epsilon) / Real.log 2 := by
    apply div_pos (Real.log_pos ?_) ha
    exact (lt_div_iff₀ hepsilon).mpr (by simpa using hepsilon_one)
  let p := ⌈Real.log (1 / epsilon) / Real.log 2⌉₊
  have hceil : Real.log (1 / epsilon) / Real.log 2 ≤ (p : ℝ) := Nat.le_ceil _
  have hp : 0 < p := by
    by_contra hnot
    have hz : p = 0 := by omega
    rw [hz, Nat.cast_zero] at hceil
    exact (not_le_of_gt hratio) hceil
  have haccuracy : (1 / 2 : ℝ) ^ p ≤ epsilon := by
    have hh := (div_le_iff₀ ha).mp hceil
    rw [Real.log_div (by norm_num : (1 : ℝ) ≠ 0) hepsilon.ne', Real.log_one, zero_sub] at hh
    calc
      (1 / 2 : ℝ) ^ p = Real.exp (-(Real.log 2) * (p : ℝ)) := by
        rw [mul_comm, Real.exp_nat_mul, Real.exp_neg,
          Real.exp_log (by norm_num : (0 : ℝ) < 2)]
        norm_num
      _ ≤ Real.exp (Real.log epsilon) := Real.exp_le_exp.mpr (by nlinarith)
      _ = epsilon := Real.exp_log hepsilon
  have hP : (p : ℝ) * Real.log 2 ≤ Real.log (2 / epsilon) := by
    have hh : (p : ℝ) < Real.log (1 / epsilon) / Real.log 2 + 1 :=
      Nat.ceil_lt_add_one hratio.le
    have hm := (lt_div_iff₀ ha).mp (by linarith :
      (p : ℝ) - 1 < Real.log (1 / epsilon) / Real.log 2)
    rw [Real.log_div (by norm_num : (1 : ℝ) ≠ 0) hepsilon.ne', Real.log_one,
      zero_sub] at hm
    rw [Real.log_div (by norm_num : (2 : ℝ) ≠ 0) hepsilon.ne']
    nlinarith
  have hL : 0 ≤ Real.log (1 + Λ) := (Real.log_pos (by linarith)).le
  have hbase : (1 / 2 : ℝ) ≤ Real.log 2 := by
    have hh := Real.le_log_one_add_of_nonneg (show (0 : ℝ) ≤ 1 by norm_num)
    norm_num at hh
    linarith
  have hprod : (p : ℝ) * ((s : ℝ) + 1) * (Real.log 2) ^ 2 ≤
      2 * Real.log (1 + Λ) * Real.log (2 / epsilon) := by
    have h1 := mul_le_mul_of_nonneg_left hS (mul_nonneg (Nat.cast_nonneg p) ha.le)
    have h2 := mul_le_mul_of_nonneg_left hP (mul_nonneg (by norm_num : (0 : ℝ) ≤ 2) hL)
    nlinarith
  have hquarter : (1 / 4 : ℝ) ≤ (Real.log 2) ^ 2 := by nlinarith
  have hscaled := mul_le_mul_of_nonneg_left hquarter
    (show 0 ≤ (p : ℝ) * ((s : ℝ) + 1) by positivity)
  have hrank : 16 * (p : ℝ) * ((s : ℝ) + 1) ≤
      128 * Real.log (1 + Λ) * Real.log (2 / epsilon) := by nlinarith
  have hlog : Real.log (2 / epsilon) ≤ Real.log (128 / epsilon) :=
    Real.log_le_log (by positivity) (div_le_div_of_nonneg_right (by norm_num) hepsilon.le)
  refine ⟨fermionicSeparatedTerms p s, ?_, ?_⟩
  · rw [fermionicSeparatedTerms_length]
    push_cast
    have hh := hrank.trans (mul_le_mul_of_nonneg_left hlog
      (show 0 ≤ 128 * Real.log (1 + Λ) by positivity))
    nlinarith
  · intro t omega ht0 ht1 hw0 hw1
    exact (fermionicKernel_separatedApprox_error p s hp ht0 ht1
      (by linarith) (by linarith)).trans haccuracy

end GECPKernelStructure.Fermionic
