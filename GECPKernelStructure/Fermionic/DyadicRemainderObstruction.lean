import GECPKernelStructure.Fermionic.ContinuationDeterminant

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- Error of the existing positive-frequency masked dyadic Taylor approximation. -/
noncomputable def dyadicTaylorError (p s : ℕ) (t ω : ℝ) : ℝ :=
  Real.exp (-t * ω) -
    evalSeparated (dyadicSeparatedTerms (8 * p) p s) t ω

/-- The first dyadic mask transition reverses the expected order-two minor sign. -/
theorem eighthOrder_maskTransition_det_pos :
    0 <
      (Real.exp (-1) - expNegTaylor 8 1) * Real.exp (-(3 / 2 : ℝ)) -
        (Real.exp (-(5 / 4 : ℝ)) - expNegTaylor 8 (5 / 4 : ℝ)) *
          (Real.exp (-(6 / 5 : ℝ)) - expNegTaylor 8 (6 / 5 : ℝ)) := by
  let a := Real.exp (-1) - expNegTaylor 8 1
  let b := Real.exp (-(5 / 4 : ℝ)) - expNegTaylor 8 (5 / 4 : ℝ)
  let c := Real.exp (-(6 / 5 : ℝ)) - expNegTaylor 8 (6 / 5 : ℝ)
  let d := Real.exp (-(3 / 2 : ℝ))
  have exp_at_one := Real.exp_bound (x := -(1 : ℝ)) (n := 9) (by norm_num) (by norm_num)
  have sum_nine :
      (∑ m ∈ Finset.range 9, (-(1 : ℝ)) ^ m / (m.factorial : ℝ)) =
        expNegTaylor 8 1 + 1 / (Nat.factorial 8 : ℝ) := by
    rw [Finset.sum_range_succ]
    norm_num [expNegTaylor]
  have a_lower : (1 / 50000 : ℝ) ≤ a := by
    rw [sum_nine] at exp_at_one
    rw [abs_le] at exp_at_one
    dsimp only [a]
    norm_num at exp_at_one ⊢
    linarith
  have b_abs := expNegTaylor_error_le_next (terms := 8) (x := (5 / 4 : ℝ))
    (by norm_num) (by norm_num)
  have b_upper : abs b ≤ (1 / 6700 : ℝ) := by
    dsimp only [b]
    norm_num at b_abs ⊢
    exact b_abs.trans (by norm_num)
  have c_abs := expNegTaylor_error_le_next (terms := 8) (x := (6 / 5 : ℝ))
    (by norm_num) (by norm_num)
  have c_upper : abs c ≤ (1 / 9000 : ℝ) := by
    dsimp only [c]
    norm_num at c_abs ⊢
    exact c_abs.trans (by norm_num)
  have exp_three_quarters := Real.add_one_le_exp (-(3 / 4 : ℝ))
  have d_lower : (1 / 16 : ℝ) ≤ d := by
    have d_factor : d = Real.exp (-(3 / 4 : ℝ)) ^ 2 := by
      dsimp only [d]
      rw [pow_two, ← Real.exp_add]
      congr 1
      ring
    rw [d_factor]
    norm_num at exp_three_quarters ⊢
    nlinarith
  have ad_lower : (1 / 800000 : ℝ) ≤ a * d := by
    calc
      (1 / 800000 : ℝ) = (1 / 50000 : ℝ) * (1 / 16 : ℝ) := by norm_num
      _ ≤ a * (1 / 16 : ℝ) := by gcongr
      _ ≤ a * d := by
        have ha : 0 ≤ a := by linarith
        exact mul_le_mul_of_nonneg_left d_lower ha
  have bc_upper : b * c ≤ (1 / 60300000 : ℝ) := by
    calc
      b * c ≤ abs (b * c) := le_abs_self _
      _ = abs b * abs c := abs_mul b c
      _ ≤ (1 / 6700 : ℝ) * (1 / 9000 : ℝ) :=
        mul_le_mul b_upper c_upper (abs_nonneg _) (by positivity)
      _ = (1 / 60300000 : ℝ) := by norm_num
  change 0 < a * d - b * c
  linarith

/-- The two time samples straddle the first dyadic time-mask transition. -/
noncomputable def maskTransitionRows : Fin 2 → ℝ :=
  orderedPair (1 / 2 : ℝ) (3 / 5 : ℝ)

/-- The two frequency samples straddle the first dyadic frequency threshold. -/
noncomputable def maskTransitionColumns : Fin 2 → ℝ :=
  orderedPair (2 : ℝ) (5 / 2 : ℝ)

theorem maskTransitionRows_strictMono : StrictMono maskTransitionRows := by
  intro i j hij
  fin_cases i <;> fin_cases j <;>
    norm_num [maskTransitionRows, orderedPair, div_eq_mul_inv] at *

theorem maskTransitionColumns_strictMono : StrictMono maskTransitionColumns := by
  intro i j hij
  fin_cases i <;> fin_cases j <;>
    norm_num [maskTransitionColumns, orderedPair, div_eq_mul_inv] at *

@[simp]
theorem dyadicTaylorError_maskTransition_zero_zero :
    dyadicTaylorError 1 2 (maskTransitionRows 0) (maskTransitionColumns 0) =
      Real.exp (-1) - expNegTaylor 8 1 := by
  norm_num [dyadicTaylorError, maskTransitionRows, maskTransitionColumns, orderedPair,
    eval_dyadicSeparatedTerms_succ, eval_dyadicSeparatedTerms_zero]

@[simp]
theorem dyadicTaylorError_maskTransition_zero_one :
    dyadicTaylorError 1 2 (maskTransitionRows 0) (maskTransitionColumns 1) =
      Real.exp (-(5 / 4 : ℝ)) - expNegTaylor 8 (5 / 4 : ℝ) := by
  norm_num [dyadicTaylorError, maskTransitionRows, maskTransitionColumns, orderedPair,
    eval_dyadicSeparatedTerms_succ, eval_dyadicSeparatedTerms_zero]

@[simp]
theorem dyadicTaylorError_maskTransition_one_zero :
    dyadicTaylorError 1 2 (maskTransitionRows 1) (maskTransitionColumns 0) =
      Real.exp (-(6 / 5 : ℝ)) - expNegTaylor 8 (6 / 5 : ℝ) := by
  norm_num [dyadicTaylorError, maskTransitionRows, maskTransitionColumns, orderedPair,
    eval_dyadicSeparatedTerms_succ, eval_dyadicSeparatedTerms_zero]

@[simp]
theorem dyadicTaylorError_maskTransition_one_one :
    dyadicTaylorError 1 2 (maskTransitionRows 1) (maskTransitionColumns 1) =
      Real.exp (-(3 / 2 : ℝ)) := by
  norm_num [dyadicTaylorError, maskTransitionRows, maskTransitionColumns, orderedPair,
    eval_dyadicSeparatedTerms_succ, eval_dyadicSeparatedTerms_zero]

/-- The concrete masked dyadic remainder has a positive ordered `2 × 2` minor. -/
theorem dyadicTaylorError_maskTransition_det_pos :
    0 < (minorMatrix (dyadicTaylorError 1 2)
      maskTransitionRows maskTransitionColumns).det := by
  rw [Matrix.det_fin_two]
  simp only [minorMatrix, dyadicTaylorError_maskTransition_zero_zero,
    dyadicTaylorError_maskTransition_zero_one,
    dyadicTaylorError_maskTransition_one_zero,
    dyadicTaylorError_maskTransition_one_one]
  simpa using eighthOrder_maskTransition_det_pos

/--
The current masked dyadic approximation cannot inherit the exponential
kernel's strict sign-regularity route: its first mask transition reverses the
required sign already at order two.
-/
theorem dyadicTaylorError_not_strictSignRegularAtOrder_two :
    ¬StrictSignRegularAtOrder (dyadicTaylorError 1 2) 2 (expKernelSignature 2) := by
  intro regular
  have signed := regular maskTransitionRows maskTransitionColumns
    maskTransitionRows_strictMono maskTransitionColumns_strictMono
  have positive := dyadicTaylorError_maskTransition_det_pos
  rw [show expKernelSignature 2 = -1 by norm_num [expKernelSignature, Nat.choose]] at signed
  nlinarith

end Fermionic
end GECPKernelStructure
