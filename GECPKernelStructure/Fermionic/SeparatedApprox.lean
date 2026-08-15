import GECPKernelStructure.Fermionic.Derivatives
import GECPKernelStructure.Fermionic.ExponentialFamily
import GECPKernelStructure.Fermionic.Symmetry

namespace GECPKernelStructure
namespace Fermionic

/-- Absorb the fermionic denominator into each frequency factor. -/
noncomputable def divideSeparatedRight (terms : List SeparatedTerm) :
    List SeparatedTerm :=
  terms.map fun term =>
    (term.1, fun ω => term.2 ω / (1 + Real.exp (-ω)))

theorem divideSeparatedRight_length (terms : List SeparatedTerm) :
    (divideSeparatedRight terms).length = terms.length := by
  simp [divideSeparatedRight]

theorem eval_divideSeparatedRight (terms : List SeparatedTerm) (t ω : ℝ) :
    evalSeparated (divideSeparatedRight terms) t ω =
      evalSeparated terms t ω / (1 + Real.exp (-ω)) := by
  induction terms with
  | nil => simp [divideSeparatedRight, evalSeparated]
  | cons term terms ih =>
      simp only [divideSeparatedRight, List.map_cons, evalSeparated, List.map_cons,
        List.sum_cons]
      change term.1 t * (term.2 ω / (1 + Real.exp (-ω))) +
          evalSeparated (divideSeparatedRight terms) t ω =
        (term.1 t * term.2 ω + evalSeparated terms t ω) / (1 + Real.exp (-ω))
      rw [ih]
      ring

/-- Reflect separated terms through `(t,ω) ↦ (1-t,-ω)`. -/
noncomputable def reflectSeparatedTerms (terms : List SeparatedTerm) :
    List SeparatedTerm :=
  terms.map fun term => (fun t => term.1 (1 - t), fun ω => term.2 (-ω))

theorem reflectSeparatedTerms_length (terms : List SeparatedTerm) :
    (reflectSeparatedTerms terms).length = terms.length := by
  simp [reflectSeparatedTerms]

theorem eval_reflectSeparatedTerms (terms : List SeparatedTerm) (t ω : ℝ) :
    evalSeparated (reflectSeparatedTerms terms) t ω = evalSeparated terms (1 - t) (-ω) := by
  induction terms with
  | nil => simp [reflectSeparatedTerms, evalSeparated]
  | cons term terms ih =>
      simp only [reflectSeparatedTerms, List.map_cons, evalSeparated, List.map_cons,
        List.sum_cons]
      change term.1 (1 - t) * term.2 (-ω) +
          evalSeparated (reflectSeparatedTerms terms) t ω =
        term.1 (1 - t) * term.2 (-ω) + evalSeparated terms (1 - t) (-ω)
      rw [ih]

/-- Explicit positive/negative-frequency fermionic separated terms. -/
noncomputable def fermionicSeparatedTerms (p s : ℕ) : List SeparatedTerm :=
  let positive := divideSeparatedRight (dyadicSeparatedTerms (8 * p) p s)
  maskSeparatedTerms positive (fun _ => True) (fun ω => 0 ≤ ω) ++
    maskSeparatedTerms (reflectSeparatedTerms positive) (fun _ => True) (fun ω => ω < 0)

theorem fermionicSeparatedTerms_length (p s : ℕ) :
    (fermionicSeparatedTerms p s).length = 2 * ((s + 1) * (8 * p)) := by
  simp [fermionicSeparatedTerms, maskSeparatedTerms_length,
    reflectSeparatedTerms_length, divideSeparatedRight_length,
    dyadicSeparatedTerms_length]
  ring

theorem eval_fermionicSeparatedTerms_nonneg (p s : ℕ) (t ω : ℝ) (hω : 0 ≤ ω) :
    evalSeparated (fermionicSeparatedTerms p s) t ω =
      evalSeparated (dyadicSeparatedTerms (8 * p) p s) t ω /
        (1 + Real.exp (-ω)) := by
  simp [fermionicSeparatedTerms, evalSeparated_append, eval_maskSeparatedTerms,
    eval_divideSeparatedRight, hω]

theorem eval_fermionicSeparatedTerms_neg (p s : ℕ) (t ω : ℝ) (hω : ω < 0) :
    evalSeparated (fermionicSeparatedTerms p s) t ω =
      evalSeparated (dyadicSeparatedTerms (8 * p) p s) (1 - t) (-ω) /
        (1 + Real.exp (-(-ω))) := by
  simp [fermionicSeparatedTerms, evalSeparated_append, eval_maskSeparatedTerms,
    eval_reflectSeparatedTerms, eval_divideSeparatedRight, hω]

theorem fermionicKernel_positive_separatedApprox (p s : ℕ) (hp : 0 < p)
    {t ω : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1)
    (hω0 : 0 ≤ ω) (hω1 : ω ≤ (2 ^ s : ℕ)) :
    |fermionicKernel t ω -
        evalSeparated (dyadicSeparatedTerms (8 * p) p s) t ω /
          (1 + Real.exp (-ω))| ≤ (1 / 2 : ℝ) ^ p := by
  have hnum := expFamily_separatedApprox_error p s hp ht0 ht1 hω0 hω1
  have hden_pos : 0 < 1 + Real.exp (-ω) := by positivity
  have hden_one : 1 ≤ 1 + Real.exp (-ω) := by
    linarith [Real.exp_pos (-ω)]
  have heps : 0 ≤ (1 / 2 : ℝ) ^ p := by positivity
  rw [fermionicKernel]
  rw [← sub_div]
  rw [abs_div, abs_of_pos hden_pos]
  apply (div_le_div_of_nonneg_right hnum hden_pos.le).trans
  apply (div_le_iff₀ hden_pos).mpr
  nlinarith [mul_le_mul_of_nonneg_left hden_one heps]

/-- Pointwise error certificate for the explicit fermionic separated terms. -/
theorem fermionicKernel_separatedApprox_error (p s : ℕ) (hp : 0 < p)
    {t ω : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1)
    (hω_lower : -((2 ^ s : ℕ) : ℝ) ≤ ω)
    (hω_upper : ω ≤ (2 ^ s : ℕ)) :
    |fermionicKernel t ω - evalSeparated (fermionicSeparatedTerms p s) t ω| ≤
      (1 / 2 : ℝ) ^ p := by
  by_cases hω : 0 ≤ ω
  · rw [eval_fermionicSeparatedTerms_nonneg p s t ω hω]
    exact fermionicKernel_positive_separatedApprox p s hp ht0 ht1 hω hω_upper
  · have hωneg : ω < 0 := lt_of_not_ge hω
    rw [eval_fermionicSeparatedTerms_neg p s t ω hωneg]
    have ht0' : 0 ≤ 1 - t := by linarith
    have ht1' : 1 - t ≤ 1 := by linarith
    have hneg0 : 0 ≤ -ω := by linarith
    have hneg1 : -ω ≤ ((2 ^ s : ℕ) : ℝ) := by linarith
    have result := fermionicKernel_positive_separatedApprox p s hp ht0' ht1' hneg0 hneg1
    rw [fermionicKernel_reflection] at result
    exact result

/--
An explicit separated approximation of the fermionic kernel on
`[0,1] × [-2^s,2^s]` with `16p(s+1)` terms and uniform error `2⁻ᵖ`.

This is a low-rank existence theorem. It makes no claim about GECP selecting
these terms or inheriting the same rate.
-/
theorem fermionicKernel_separatedApprox (p s : ℕ) (hp : 0 < p) :
    ∃ terms : List SeparatedTerm,
      terms.length = 2 * ((s + 1) * (8 * p)) ∧
      ∀ {t ω : ℝ}, 0 ≤ t → t ≤ 1 →
        -((2 ^ s : ℕ) : ℝ) ≤ ω → ω ≤ (2 ^ s : ℕ) →
        |fermionicKernel t ω - evalSeparated terms t ω| ≤ (1 / 2 : ℝ) ^ p := by
  refine ⟨fermionicSeparatedTerms p s, fermionicSeparatedTerms_length p s, ?_⟩
  intro t ω ht0 ht1 hω_lower hω_upper
  exact fermionicKernel_separatedApprox_error p s hp ht0 ht1 hω_lower hω_upper

end Fermionic
end GECPKernelStructure
