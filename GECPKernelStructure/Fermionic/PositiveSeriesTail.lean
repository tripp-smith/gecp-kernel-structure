import GECPKernelStructure.Fermionic.TailGECPEnvelope

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- Separated positive monomials for the Taylor polynomial of `exp (u * omega)`. -/
noncomputable def positiveTaylorSeparatedTerms (terms : ℕ) : List SeparatedTerm :=
  (List.range terms).map fun k =>
    (fun u : ℝ => u ^ k / (k.factorial : ℝ), fun ω : ℝ => ω ^ k)

theorem positiveTaylorSeparatedTerms_length (terms : ℕ) :
    (positiveTaylorSeparatedTerms terms).length = terms := by
  simp [positiveTaylorSeparatedTerms]

theorem eval_positiveTaylorSeparatedTerms (terms : ℕ) (u ω : ℝ) :
    evalSeparated (positiveTaylorSeparatedTerms terms) u ω =
      expNegTaylor terms (-(u * ω)) := by
  induction terms with
  | zero => simp [positiveTaylorSeparatedTerms, evalSeparated, expNegTaylor]
  | succ terms ih =>
      calc
        evalSeparated (positiveTaylorSeparatedTerms (terms + 1)) u ω =
            evalSeparated (positiveTaylorSeparatedTerms terms) u ω +
              u ^ terms / (terms.factorial : ℝ) * ω ^ terms := by
                simp [positiveTaylorSeparatedTerms, evalSeparated, List.range_succ]
        _ = expNegTaylor terms (-(u * ω)) +
              u ^ terms / (terms.factorial : ℝ) * ω ^ terms := by rw [ih]
        _ = expNegTaylor (terms + 1) (-(u * ω)) := by
          rw [show expNegTaylor (terms + 1) (-(u * ω)) =
              expNegTaylor terms (-(u * ω)) +
                (-(-(u * ω))) ^ terms / (terms.factorial : ℝ) by
            simp [expNegTaylor, Finset.sum_range_succ]]
          congr 1
          ring

/-- Absorb the positive-coordinate fermionic denominator into each right feature. -/
noncomputable def dividePositiveSeparatedRight (terms : List SeparatedTerm) :
    List SeparatedTerm :=
  terms.map fun term =>
    (term.1, fun ω => term.2 ω / (1 + Real.exp ω))

theorem dividePositiveSeparatedRight_length (terms : List SeparatedTerm) :
    (dividePositiveSeparatedRight terms).length = terms.length := by
  simp [dividePositiveSeparatedRight]

theorem eval_dividePositiveSeparatedRight
    (terms : List SeparatedTerm) (u ω : ℝ) :
    evalSeparated (dividePositiveSeparatedRight terms) u ω =
      evalSeparated terms u ω / (1 + Real.exp ω) := by
  induction terms with
  | nil => simp [dividePositiveSeparatedRight, evalSeparated]
  | cons term terms ih =>
      simp only [dividePositiveSeparatedRight, List.map_cons, evalSeparated,
        List.sum_cons]
      change term.1 u * (term.2 ω / (1 + Real.exp ω)) +
          evalSeparated (dividePositiveSeparatedRight terms) u ω =
        (term.1 u * term.2 ω + evalSeparated terms u ω) /
          (1 + Real.exp ω)
      rw [ih]
      ring

private lemma positiveExpTaylor_eq_taylorWithinEval
    {terms : ℕ} (hterms : 0 < terms) {x : ℝ} (hx : 0 < x) :
    expNegTaylor terms (-x) =
      taylorWithinEval Real.exp (terms - 1) (Set.Icc 0 x) 0 x := by
  rw [taylor_within_apply]
  unfold expNegTaylor
  rw [Nat.sub_add_cancel hterms]
  apply Finset.sum_congr rfl
  intro k _hk
  rw [iteratedDerivWithin_eq_iteratedDeriv (uniqueDiffOn_Icc hx)
    (by fun_prop) (Set.left_mem_Icc.mpr hx.le)]
  rw [show Real.exp = (fun y : ℝ => Real.exp ((1 : ℝ) * y)) by
    funext y
    congr 1
    ring]
  rw [congrFun (iteratedDeriv_exp_const_mul k 1) 0]
  simp only [smul_eq_mul, mul_zero, Real.exp_zero, mul_one, one_pow]
  rw [show (- -x) ^ k = x ^ k by ring]
  ring

/-- Lagrange envelope for the positive exponential Taylor remainder. -/
theorem positiveExpTaylor_error_le_next_mul_exp
    {terms : ℕ} (hterms : 0 < terms) {x : ℝ} (hx0 : 0 ≤ x) :
    abs (Real.exp x - expNegTaylor terms (-x)) ≤
      Real.exp x * x ^ terms / (terms.factorial : ℝ) := by
  rcases hx0.eq_or_lt with rfl | hx
  · obtain ⟨n, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hterms)
    have hz : expNegTaylor (n + 1) (-0) = 1 := by
      unfold expNegTaylor
      rw [Finset.sum_eq_single 0]
      · norm_num
      · intro b _hb hne
        simp [hne]
      · simp
    rw [hz]
    simp
  · rw [positiveExpTaylor_eq_taylorWithinEval hterms hx]
    obtain ⟨point, point_mem, remainder⟩ :=
      taylor_mean_remainder_lagrange_iteratedDeriv
        (f := Real.exp) (x := x) (x₀ := 0)
        (n := terms - 1) hx.ne (by fun_prop : ContDiffOn ℝ (terms - 1 + 1)
          Real.exp (Set.uIcc 0 x))
    have hsub : terms - 1 + 1 = terms := by omega
    simp only [hsub] at remainder
    simp [Set.uIcc, hx.le] at remainder point_mem
    rw [remainder, abs_div, abs_mul]
    rw [show Real.exp = (fun y : ℝ => Real.exp ((1 : ℝ) * y)) by
      funext y
      congr 1
      ring]
    rw [congrFun (iteratedDeriv_exp_const_mul terms 1) point]
    simp only [one_mul, one_pow, abs_pow, abs_of_nonneg hx0,
      abs_of_nonneg (show (0 : ℝ) ≤ (terms.factorial : ℝ) by positivity)]
    rw [abs_of_pos (Real.exp_pos point)]
    have point_le_x : point ≤ x := point_mem.2.le
    have exp_le : Real.exp point ≤ Real.exp x := Real.exp_le_exp.mpr point_le_x
    apply div_le_div_of_nonneg_right _ (by positivity)
    exact mul_le_mul_of_nonneg_right exp_le (pow_nonneg hx0 _)

/-- The eighth-order tail theorem only needs nonzero rows because its row power is even. -/
theorem eighthOrderTail_minor_signed_pos_of_rows_ne_zero
    {n : ℕ} (rows columns : Fin (n + 1) → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_ne : ∀ i, rows i ≠ 0) (columns_pos : ∀ j, 0 < columns j) :
    0 < expKernelSignature (n + 1) *
      Matrix.det (fun i j => eighthOrderTail (rows i * columns j)) := by
  let upper := columns (Fin.last n)
  let composition : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    fun i j => ∫ v in (0 : ℝ)..upper,
      expKernel (rows i) v * eighthOrderTruncatedPowerKernel v (columns j)
  let weights : Fin (n + 1) → ℝ :=
    fun i => rows i ^ 8 / (Nat.factorial 7 : ℝ)
  let scaled : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    fun i j => weights i * composition i j
  let tailMatrix : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    fun i j => eighthOrderTail (rows i * columns j)
  have composition_pos :
      0 < expKernelSignature (n + 1) * composition.det := by
    exact eighthOrderComposition_interval_det_pos
      rows columns rows_mono columns_mono columns_pos
  have tail_matrix_eq : tailMatrix = scaled := by
    ext i j
    dsimp only [tailMatrix, scaled, weights, composition]
    rw [eighthOrderTail_eq_truncatedPowerComposition_commonCutoff
      (columns_pos j).le
      (columns_mono.monotone (Fin.le_last j))]
    congr 1
    apply intervalIntegral.integral_congr
    intro v _
    simp only [expKernel]
    ring
  have weights_product_pos : 0 < ∏ i, weights i := by
    apply Finset.prod_pos
    intro i _
    have row_pow_pos : 0 < rows i ^ 8 :=
      (by norm_num : Even 8).pow_pos (rows_ne i)
    exact div_pos row_pow_pos (by positivity)
  have scaled_det : scaled.det = (∏ i, weights i) * composition.det := by
    exact Matrix.det_mul_column weights composition
  change 0 < expKernelSignature (n + 1) * tailMatrix.det
  rw [tail_matrix_eq, scaled_det]
  calc
    0 < (∏ i, weights i) *
        (expKernelSignature (n + 1) * composition.det) :=
      mul_pos weights_product_pos composition_pos
    _ = expKernelSignature (n + 1) *
        ((∏ i, weights i) * composition.det) := by ring

/-- The positive exponential tail, written through the already proved negative tail. -/
noncomputable def positiveEighthOrderTail (x : ℝ) : ℝ :=
  eighthOrderTail (-x)

/-- The positive-series eighth-order tail kernel on the positive quadrant. -/
noncomputable def positiveSeriesEighthOrderTailKernel :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun u ω => positiveEighthOrderTail ((u : ℝ) * (ω : ℝ))

/-- Positive exponential Taylor tails are strictly totally positive. -/
theorem positiveSeriesEighthOrderTailKernel_strictSignRegular :
    StrictSignRegular positiveSeriesEighthOrderTailKernel (fun _ => 1) := by
  intro order
  cases order with
  | zero =>
      intro rows columns _rows_mono _columns_mono
      simp
  | succ n =>
      intro rows columns rows_mono columns_mono
      let reversedRows : Fin (n + 1) → ℝ :=
        fun i => -(rows i.rev : ℝ)
      have reversedRows_mono : StrictMono reversedRows := by
        intro i j hij
        exact neg_lt_neg (rows_mono (by simpa [reversedRows] using hij))
      have reversedRows_ne : ∀ i, reversedRows i ≠ 0 := by
        intro i
        dsimp only [reversedRows]
        exact neg_ne_zero.mpr (ne_of_gt (rows i.rev).property)
      have signed_reversed := eighthOrderTail_minor_signed_pos_of_rows_ne_zero
        reversedRows (fun j => (columns j : ℝ))
        reversedRows_mono (fun _ _ hij => columns_mono hij)
        reversedRows_ne (fun j => (columns j).property)
      let original : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
        fun i j => positiveSeriesEighthOrderTailKernel (rows i) (columns j)
      have reversed_matrix :
          (fun i j => eighthOrderTail (reversedRows i * (columns j : ℝ))) =
            original.submatrix Fin.revPerm id := by
        ext i j
        change eighthOrderTail (reversedRows i * (columns j : ℝ)) =
          original (Fin.revPerm i) (id j)
        dsimp only [reversedRows, original, positiveSeriesEighthOrderTailKernel,
          positiveEighthOrderTail, id_eq]
        rw [Fin.revPerm_apply]
        congr 1
        ring
      have sign_cast :
          (((Equiv.Perm.sign (@Fin.revPerm (n + 1)) : ℤ)) : ℝ) =
            expKernelSignature (n + 1) := by
        rw [expKernelSignature]
        norm_cast
        exact congrArg Units.val (Fin.sign_revPerm (n + 1))
      rw [reversed_matrix, Matrix.det_permute, sign_cast] at signed_reversed
      have signature_sq :
          expKernelSignature (n + 1) * expKernelSignature (n + 1) = 1 := by
        simp only [expKernelSignature]
        rw [← pow_add, ← two_mul, pow_mul]
        norm_num
      change 0 < (1 : ℝ) * original.det
      nlinarith

/-- The eight positive Taylor terms in the `u * omega` coordinate. -/
noncomputable def positiveSeriesEighthOrderTaylorApprox :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun u ω => expNegTaylor 8 (-((u : ℝ) * (ω : ℝ))) *
    (1 + Real.exp (ω : ℝ))⁻¹

/-- The denominator-scaled positive-series tail. -/
noncomputable def positiveSeriesEighthOrderFermionicTailKernel :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun u ω => positiveSeriesEighthOrderTailKernel u ω *
    (1 + Real.exp (ω : ℝ))⁻¹

/-- The positive-coordinate Taylor component uses exactly eight separated terms. -/
theorem positiveSeriesEighthOrderTaylorApprox_eq_evalSeparated
    (u ω : PositiveArgument) :
    positiveSeriesEighthOrderTaylorApprox u ω =
      evalSeparated
        (dividePositiveSeparatedRight (positiveTaylorSeparatedTerms 8))
        (u : ℝ) (ω : ℝ) := by
  rw [eval_dividePositiveSeparatedRight, eval_positiveTaylorSeparatedTerms]
  simp only [positiveSeriesEighthOrderTaylorApprox, div_eq_mul_inv]

/-- The positive coordinate is exactly the original time reflection. -/
theorem fermionicKernel_one_sub_eq_positive_exp
    (u ω : ℝ) :
    fermionicKernel (1 - u) ω =
      Real.exp (u * ω) / (1 + Real.exp ω) := by
  have reflected := fermionicKernel_reflection u (-ω)
  rw [neg_neg] at reflected
  rw [reflected]
  simp only [fermionicKernel]
  rw [show -u * -ω = u * ω by ring, show - -ω = ω by ring]

/-- Exact rank-eight-plus-tail decomposition in the positive-series coordinate. -/
theorem fermionicKernel_one_sub_eq_positiveSeriesTaylorApprox_add_tail
    (u ω : PositiveArgument) :
    fermionicKernel (1 - (u : ℝ)) (ω : ℝ) =
      positiveSeriesEighthOrderTaylorApprox u ω +
        positiveSeriesEighthOrderFermionicTailKernel u ω := by
  rw [fermionicKernel_one_sub_eq_positive_exp]
  simp only [positiveSeriesEighthOrderTaylorApprox,
    positiveSeriesEighthOrderFermionicTailKernel,
    positiveSeriesEighthOrderTailKernel, positiveEighthOrderTail,
    eighthOrderTail, div_eq_mul_inv]
  ring_nf

/-- Positive denominator scaling preserves strict total positivity of the tail. -/
theorem positiveSeriesEighthOrderFermionicTailKernel_strictSignRegular :
    StrictSignRegular positiveSeriesEighthOrderFermionicTailKernel (fun _ => 1) := by
  intro order
  have scaled := strictSignRegularAtOrder_columnScale
    (positiveSeriesEighthOrderTailKernel_strictSignRegular order)
    (fun ω : PositiveArgument => (1 + Real.exp (ω : ℝ))⁻¹)
    (fun _ => inv_pos.mpr (by positivity))
  change StrictSignRegularAtOrder
    (fun u ω : PositiveArgument => positiveSeriesEighthOrderTailKernel u ω *
      (1 + Real.exp (ω : ℝ))⁻¹) order 1
  exact scaled

/-- The positive-series fermionic tail is at most `1 / 8!` on the positive unit square. -/
theorem positiveSeriesEighthOrderFermionicTailKernel_abs_le_inv_factorial
    {u ω : PositiveArgument} (hu : (u : ℝ) ≤ 1) (hω : (ω : ℝ) ≤ 1) :
    abs (positiveSeriesEighthOrderFermionicTailKernel u ω) ≤
      1 / (Nat.factorial 8 : ℝ) := by
  let x := (u : ℝ) * (ω : ℝ)
  let denominator := 1 + Real.exp (ω : ℝ)
  have x_nonneg : 0 ≤ x := mul_nonneg u.property.le ω.property.le
  have x_le_omega : x ≤ (ω : ℝ) := by
    dsimp only [x]
    nlinarith [mul_le_mul_of_nonneg_right hu ω.property.le]
  have x_le_one : x ≤ 1 := x_le_omega.trans hω
  have denominator_pos : 0 < denominator := by
    dsimp only [denominator]
    positivity
  have exp_le_denominator : Real.exp x ≤ denominator := by
    dsimp only [denominator]
    exact (Real.exp_le_exp.mpr x_le_omega).trans (le_add_of_nonneg_left (by norm_num))
  have remainder := positiveExpTaylor_error_le_next_mul_exp
    (terms := 8) (by norm_num) x_nonneg
  have quotient_bound :
      abs (Real.exp x - expNegTaylor 8 (-x)) / denominator ≤
        x ^ 8 / (Nat.factorial 8 : ℝ) := by
    apply (div_le_iff₀ denominator_pos).2
    calc
      abs (Real.exp x - expNegTaylor 8 (-x)) ≤
          Real.exp x * x ^ 8 / (Nat.factorial 8 : ℝ) := remainder
      _ ≤ denominator * x ^ 8 / (Nat.factorial 8 : ℝ) := by
        gcongr
      _ = x ^ 8 / (Nat.factorial 8 : ℝ) * denominator := by ring
  calc
    abs (positiveSeriesEighthOrderFermionicTailKernel u ω) =
        abs (Real.exp x - expNegTaylor 8 (-x)) / denominator := by
      simp only [positiveSeriesEighthOrderFermionicTailKernel,
        positiveSeriesEighthOrderTailKernel, positiveEighthOrderTail,
        eighthOrderTail, x, denominator, abs_mul, abs_inv,
        abs_of_pos denominator_pos, div_eq_mul_inv, neg_neg]
    _ ≤ x ^ 8 / (Nat.factorial 8 : ℝ) := quotient_bound
    _ ≤ 1 / (Nat.factorial 8 : ℝ) := by
      gcongr
      exact pow_le_one₀ x_nonneg x_le_one

/-- The positive-series eight-term fermionic approximation inherits the tail bound. -/
theorem positiveSeriesEighthOrderTaylorApprox_error_le_inv_factorial
    {u ω : PositiveArgument} (hu : (u : ℝ) ≤ 1) (hω : (ω : ℝ) ≤ 1) :
    abs (fermionicKernel (1 - (u : ℝ)) (ω : ℝ) -
      positiveSeriesEighthOrderTaylorApprox u ω) ≤
        1 / (Nat.factorial 8 : ℝ) := by
  rw [fermionicKernel_one_sub_eq_positiveSeriesTaylorApprox_add_tail]
  ring_nf
  exact positiveSeriesEighthOrderFermionicTailKernel_abs_le_inv_factorial hu hω

/-- Complete pivoting on the positive-series tail has the same factorial-free envelope. -/
theorem positiveSeriesEighthOrderFermionicTailKernel_completeRun_selectedCore_det_le
    (run : Run positiveSeriesEighthOrderFermionicTailKernel)
    (complete : run.CompleteOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)) :
    abs run.selectedCore.det ≤
      (1 / (Nat.factorial 8 : ℝ)) ^ Fintype.card run.SelectedIndex := by
  exact strictSignRegular_completeContinuation_selectedCore_det_le_pow
    positiveSeriesEighthOrderFermionicTailKernel_strictSignRegular
    (fun u : PositiveArgument => (u : ℝ) ≤ 1)
    (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)
    (Run.nil positiveSeriesEighthOrderFermionicTailKernel) run
    (fun _u _ω hu hω =>
      positiveSeriesEighthOrderFermionicTailKernel_abs_le_inv_factorial hu hω)
    complete

end Fermionic
end GECPKernelStructure
