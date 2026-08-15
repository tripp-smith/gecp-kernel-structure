import GECPKernelStructure.Fermionic.DyadicPartition
import Mathlib.Analysis.Calculus.Taylor
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Stirling

namespace GECPKernelStructure
namespace Fermionic

open scoped BigOperators
open Finset Set

/-- Natural-number node count for a scale-by-order exponential construction. -/
def expFamilyNodeCount (scales order : ℕ) : ℕ := (scales + 1) * (order + 1)

theorem expFamilyNodeCount_pos (scales order : ℕ) :
    0 < expFamilyNodeCount scales order := by
  simp [expFamilyNodeCount]

theorem expFamilyNodeCount_mono {s₁ s₂ n₁ n₂ : ℕ}
    (hs : s₁ ≤ s₂) (hn : n₁ ≤ n₂) :
    expFamilyNodeCount s₁ n₁ ≤ expFamilyNodeCount s₂ n₂ := by
  exact Nat.mul_le_mul (Nat.add_le_add_right hs 1) (Nat.add_le_add_right hn 1)

/-- The first `terms` terms of the Taylor series for `exp (-x)`. -/
noncomputable def expNegTaylor (terms : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range terms, (-x) ^ k / (k.factorial : ℝ)

lemma expNegTaylor_eq_taylorWithinEval {terms : ℕ} (hterms : 0 < terms)
    {x : ℝ} (hx : 0 < x) :
    expNegTaylor terms x =
      taylorWithinEval (fun y : ℝ => Real.exp (-y)) (terms - 1) (Icc 0 x) 0 x := by
  rw [taylor_within_apply]
  unfold expNegTaylor
  rw [Nat.sub_add_cancel hterms]
  apply Finset.sum_congr rfl
  intro k hk
  rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hx)
    (by fun_prop) (left_mem_Icc.mpr hx.le)]
  rw [show (fun y : ℝ => Real.exp (-y)) = (fun y : ℝ => Real.exp ((-1) * y)) by
    funext y; congr 1; ring]
  rw [congrFun (iteratedDeriv_exp_const_mul k (-1)) 0]
  simp only [smul_eq_mul, sub_zero]
  rw [show (-x) ^ k = (-1 : ℝ) ^ k * x ^ k by
    rw [show -x = (-1 : ℝ) * x by ring, mul_pow]]
  simp [div_eq_inv_mul, Real.exp_zero]
  ring

/-- Lagrange-form remainder for the truncated negative exponential. -/
theorem expNegTaylor_error_le_next {terms : ℕ} (hterms : 0 < terms)
    {x : ℝ} (hx0 : 0 ≤ x) :
    |Real.exp (-x) - expNegTaylor terms x| ≤
      x ^ terms / (terms.factorial : ℝ) := by
  rcases hx0.eq_or_lt with rfl | hx
  · obtain ⟨n, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hterms)
    have hz : expNegTaylor (n + 1) 0 = 1 := by
      unfold expNegTaylor
      rw [Finset.sum_eq_single 0]
      · norm_num
      · intro b hb hne
        simp [hne]
      · simp
    rw [hz]
    simp
  · rw [expNegTaylor_eq_taylorWithinEval hterms hx]
    obtain ⟨point, point_mem, remainder⟩ :=
      taylor_mean_remainder_lagrange_iteratedDeriv
        (f := fun y : ℝ => Real.exp (-y)) (x := x) (x₀ := 0)
        (n := terms - 1) hx.ne (by fun_prop : ContDiffOn ℝ (terms - 1 + 1)
          (fun y : ℝ => Real.exp (-y)) (uIcc 0 x))
    have hsub : terms - 1 + 1 = terms := by omega
    simp only [hsub] at remainder
    simp [Set.uIcc, hx.le] at remainder point_mem
    rw [remainder, abs_div, abs_mul, abs_pow]
    rw [show (fun y : ℝ => Real.exp (-y)) = (fun y : ℝ => Real.exp ((-1) * y)) by
      funext y; congr 1; ring]
    rw [congrFun (iteratedDeriv_exp_const_mul terms (-1)) point]
    simp only [abs_of_nonneg hx0,
      abs_of_nonneg (show (0 : ℝ) ≤ (terms.factorial : ℝ) by positivity)]
    rw [show (-1 : ℝ) * point = -point by ring]
    have hpoint : 0 < point := point_mem.1
    have hexp : Real.exp (-point) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    rw [abs_mul, abs_pow]
    norm_num only [abs_neg, abs_one, one_pow, one_mul]
    rw [abs_of_pos (Real.exp_pos (-point))]
    exact div_le_div_of_nonneg_right (mul_le_of_le_one_left (by positivity) hexp)
      (by positivity)

/-- Eight Taylor terms per accuracy bit suffice uniformly on `0 ≤ x ≤ 2p`. -/
theorem pow_div_factorial_eight_mul_le (p : ℕ) (hp : 0 < p)
    {x : ℝ} (hx0 : 0 ≤ x) (hx : x ≤ 2 * p) :
    x ^ (8 * p) / ((8 * p).factorial : ℝ) ≤ (1 / 2 : ℝ) ^ p := by
  let n := 8 * p
  have hn : 0 < n := by simp [n, hp]
  have hsqrt : 1 ≤ Real.sqrt (2 * Real.pi * n) := by
    have hpi : 3 < Real.pi := Real.pi_gt_three
    have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
    have : (1 : ℝ) ≤ 2 * Real.pi * n := by nlinarith
    exact (Real.le_sqrt (by norm_num : (0 : ℝ) ≤ 1) (by positivity)).mpr
      (by simpa using this)
  have hstirling := Stirling.le_factorial_stirling n
  have hbase : (n / Real.exp 1) ^ n ≤ (n.factorial : ℝ) := by
    calc
      (n / Real.exp 1) ^ n ≤ Real.sqrt (2 * Real.pi * n) * (n / Real.exp 1) ^ n := by
        exact le_mul_of_one_le_left (by positivity) hsqrt
      _ ≤ (n.factorial : ℝ) := hstirling
  have hexp : Real.exp 1 < 3 := Real.exp_one_lt_three
  have hnp : (0 : ℝ) ≤ n := by positivity
  have hthree : ((n : ℝ) / 3) ≤ n / Real.exp 1 := by
    gcongr
  have hthird : ((n : ℝ) / 3) ^ n ≤ (n.factorial : ℝ) :=
    (pow_le_pow_left₀ (by positivity) hthree n).trans hbase
  have hxthird : x ≤ (3 / 4 : ℝ) * ((n : ℝ) / 3) := by
    simp [n] at hx ⊢
    linarith
  have hxpow : x ^ n ≤ (3 / 4 : ℝ) ^ n * ((n : ℝ) / 3) ^ n := by
    rw [← mul_pow]
    exact pow_le_pow_left₀ hx0 hxthird n
  have hratio : x ^ n / (n.factorial : ℝ) ≤ (3 / 4 : ℝ) ^ n := by
    apply (div_le_iff₀ (by positivity : (0 : ℝ) < (n.factorial : ℝ))).mpr
    calc
      x ^ n ≤ (3 / 4 : ℝ) ^ n * ((n : ℝ) / 3) ^ n := hxpow
      _ ≤ (3 / 4 : ℝ) ^ n * (n.factorial : ℝ) := by gcongr
  have hnum : (3 / 4 : ℝ) ^ 8 ≤ (1 / 2 : ℝ) := by norm_num
  have hpows : (3 / 4 : ℝ) ^ n ≤ (1 / 2 : ℝ) ^ p := by
    rw [show n = 8 * p by rfl, pow_mul]
    exact pow_le_pow_left₀ (by positivity) hnum p
  simpa [n] using hratio.trans hpows

/-- One separated term, represented by its time and frequency factors. -/
abbrev SeparatedTerm := (ℝ → ℝ) × (ℝ → ℝ)

/-- Evaluate a finite list of separated terms. -/
noncomputable def evalSeparated (terms : List SeparatedTerm) (t ω : ℝ) : ℝ :=
  (terms.map fun term => term.1 t * term.2 ω).sum

noncomputable def taylorSeparatedTerms (q : ℕ) : List SeparatedTerm :=
  (List.range q).map fun k =>
    (fun t : ℝ => (-t) ^ k / (k.factorial : ℝ), fun ω : ℝ => ω ^ k)

theorem taylorSeparatedTerms_length (q : ℕ) :
    (taylorSeparatedTerms q).length = q := by simp [taylorSeparatedTerms]

theorem eval_taylorSeparatedTerms (q : ℕ) (t ω : ℝ) :
    evalSeparated (taylorSeparatedTerms q) t ω = expNegTaylor q (t * ω) := by
  induction q with
  | zero => simp [evalSeparated, taylorSeparatedTerms, expNegTaylor]
  | succ q ih =>
      calc
        evalSeparated (taylorSeparatedTerms (q + 1)) t ω =
            evalSeparated (taylorSeparatedTerms q) t ω +
              (-t) ^ q / (q.factorial : ℝ) * ω ^ q := by
                simp [evalSeparated, taylorSeparatedTerms, List.range_succ]
        _ = expNegTaylor q (t * ω) + (-t) ^ q / (q.factorial : ℝ) * ω ^ q := by
          rw [ih]
        _ = expNegTaylor (q + 1) (t * ω) := by
          rw [show expNegTaylor (q + 1) (t * ω) =
              expNegTaylor q (t * ω) + (-(t * ω)) ^ q / (q.factorial : ℝ) by
            simp [expNegTaylor, Finset.sum_range_succ]]
          rw [show (-(t * ω)) ^ q = (-t) ^ q * ω ^ q by
            rw [show -(t * ω) = (-t) * ω by ring, mul_pow]]
          ring

noncomputable def maskSeparatedTerms (terms : List SeparatedTerm)
    (leftMask rightMask : ℝ → Prop) [DecidablePred leftMask] [DecidablePred rightMask] :
    List SeparatedTerm :=
  terms.map fun term =>
    (fun t => if leftMask t then term.1 t else 0,
      fun ω => if rightMask ω then term.2 ω else 0)

theorem maskSeparatedTerms_length (terms : List SeparatedTerm)
    (leftMask rightMask : ℝ → Prop) [DecidablePred leftMask] [DecidablePred rightMask] :
    (maskSeparatedTerms terms leftMask rightMask).length = terms.length := by
  simp [maskSeparatedTerms]

theorem eval_maskSeparatedTerms (terms : List SeparatedTerm)
    (leftMask rightMask : ℝ → Prop) [DecidablePred leftMask] [DecidablePred rightMask]
    (t ω : ℝ) :
    evalSeparated (maskSeparatedTerms terms leftMask rightMask) t ω =
      if leftMask t ∧ rightMask ω then evalSeparated terms t ω else 0 := by
  by_cases hleft : leftMask t <;> by_cases hright : rightMask ω <;>
    simp [evalSeparated, maskSeparatedTerms, hleft, hright, List.map_map,
      Function.comp_def]

theorem evalSeparated_append (first second : List SeparatedTerm) (t ω : ℝ) :
    evalSeparated (first ++ second) t ω =
      evalSeparated first t ω + evalSeparated second t ω := by
  simp [evalSeparated]

/-- Dyadic time-cutoff construction. Scale `s` covers frequencies through `2^s`. -/
noncomputable def dyadicSeparatedTerms (q p : ℕ) : ℕ → List SeparatedTerm
  | 0 => taylorSeparatedTerms q
  | s + 1 =>
      maskSeparatedTerms (dyadicSeparatedTerms q p s) (fun _ => True)
          (fun ω => ω ≤ (2 ^ s : ℕ)) ++
        maskSeparatedTerms (taylorSeparatedTerms q)
          (fun t => t * (2 ^ s : ℕ) ≤ p) (fun ω => (2 ^ s : ℕ) < ω)

theorem dyadicSeparatedTerms_length (q p s : ℕ) :
    (dyadicSeparatedTerms q p s).length = (s + 1) * q := by
  induction s with
  | zero => simp [dyadicSeparatedTerms, taylorSeparatedTerms_length]
  | succ s ih =>
      simp [dyadicSeparatedTerms, maskSeparatedTerms_length,
        taylorSeparatedTerms_length, ih]
      ring

theorem eval_dyadicSeparatedTerms_zero (q p : ℕ) (t ω : ℝ) :
    evalSeparated (dyadicSeparatedTerms q p 0) t ω = expNegTaylor q (t * ω) := by
  simp [dyadicSeparatedTerms, eval_taylorSeparatedTerms]

theorem eval_dyadicSeparatedTerms_succ (q p s : ℕ) (t ω : ℝ) :
    evalSeparated (dyadicSeparatedTerms q p (s + 1)) t ω =
      if ω ≤ (2 ^ s : ℕ) then evalSeparated (dyadicSeparatedTerms q p s) t ω
      else if t * (2 ^ s : ℕ) ≤ p then expNegTaylor q (t * ω) else 0 := by
  rw [dyadicSeparatedTerms, evalSeparated_append, eval_maskSeparatedTerms,
    eval_maskSeparatedTerms, eval_taylorSeparatedTerms]
  split_ifs <;> simp_all <;> linarith

theorem exp_neg_nat_le_two_pow_neg (p : ℕ) :
    Real.exp (-(p : ℝ)) ≤ (1 / 2 : ℝ) ^ p := by
  rw [show -(p : ℝ) = (p : ℝ) * (-1) by ring, Real.exp_nat_mul]
  exact pow_le_pow_left₀ (Real.exp_nonneg _) Real.exp_neg_one_lt_half.le p

/-- Pointwise error certificate for the explicit dyadic exponential construction. -/
theorem expFamily_separatedApprox_error (p s : ℕ) (hp : 0 < p)
    {t ω : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1)
    (hω0 : 0 ≤ ω) (hω1 : ω ≤ (2 ^ s : ℕ)) :
    |Real.exp (-t * ω) -
        evalSeparated (dyadicSeparatedTerms (8 * p) p s) t ω| ≤
      (1 / 2 : ℝ) ^ p := by
  have hrewrite : -t * ω = -(t * ω) := by ring
  rw [hrewrite]
  induction s with
  | zero =>
      rw [eval_dyadicSeparatedTerms_zero]
      apply (expNegTaylor_error_le_next (by simp [hp]) (mul_nonneg ht0 hω0)).trans
      apply pow_div_factorial_eight_mul_le p hp (mul_nonneg ht0 hω0)
      norm_num at hω1
      have htw : t * ω ≤ 1 * 1 := mul_le_mul ht1 hω1 hω0 (by norm_num)
      have hp1 : (1 : ℝ) ≤ p := by exact_mod_cast hp
      nlinarith
  | succ s ih =>
      rw [eval_dyadicSeparatedTerms_succ]
      by_cases hlow : ω ≤ (2 ^ s : ℕ)
      · rw [if_pos hlow]
        exact ih hlow
      · rw [if_neg hlow]
        have hhigh : (2 ^ s : ℕ) < ω := lt_of_not_ge hlow
        by_cases hcut : t * (2 ^ s : ℕ) ≤ p
        · rw [if_pos hcut]
          apply (expNegTaylor_error_le_next (by simp [hp]) (mul_nonneg ht0 hω0)).trans
          apply pow_div_factorial_eight_mul_le p hp (mul_nonneg ht0 hω0)
          have hpow : ((2 ^ (s + 1) : ℕ) : ℝ) = 2 * ((2 ^ s : ℕ) : ℝ) := by
            norm_num [pow_succ]
            ring
          rw [hpow] at hω1
          nlinarith
        · rw [if_neg hcut, sub_zero, abs_of_pos (Real.exp_pos (-(t * ω)))]
          have hprod : (p : ℝ) < t * ω := by
            have hbpos : (0 : ℝ) < (2 ^ s : ℕ) := by positivity
            have htpos : 0 < t := by
              by_contra h
              have : t = 0 := le_antisymm (le_of_not_gt h) ht0
              subst t
              simp at hcut
            have hmul : t * (2 ^ s : ℕ) < t * ω :=
              mul_lt_mul_of_pos_left hhigh htpos
            exact (lt_of_not_ge hcut).trans hmul
          exact (Real.exp_le_exp.mpr (by linarith)).trans (exp_neg_nat_le_two_pow_neg p)

/--
An explicit separated approximation of `exp (-tω)` on
`[0,1] × [0,2^s]` with `8p(s+1)` terms and uniform error `2⁻ᵖ`.
-/
theorem expFamily_separatedApprox (p s : ℕ) (hp : 0 < p) :
    ∃ terms : List SeparatedTerm,
      terms.length = (s + 1) * (8 * p) ∧
      ∀ {t ω : ℝ}, 0 ≤ t → t ≤ 1 → 0 ≤ ω → ω ≤ (2 ^ s : ℕ) →
        |Real.exp (-t * ω) - evalSeparated terms t ω| ≤ (1 / 2 : ℝ) ^ p := by
  refine ⟨dyadicSeparatedTerms (8 * p) p s, dyadicSeparatedTerms_length _ _ _, ?_⟩
  intro t ω ht0 ht1 hω0 hω1
  exact expFamily_separatedApprox_error p s hp ht0 ht1 hω0 hω1

end Fermionic
end GECPKernelStructure
