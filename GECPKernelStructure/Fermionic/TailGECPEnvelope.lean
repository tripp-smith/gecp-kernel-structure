import GECPKernelStructure.Fermionic.ContinuationDeterminant
import GECPKernelStructure.Fermionic.SeparatedApprox
import GECPKernelStructure.Fermionic.TruncatedPowerStrictSupport

namespace GECPKernelStructure
namespace Fermionic

open GECP

/-- A positive real argument, used to state the tail theorem as global strict sign regularity. -/
abbrev PositiveArgument := {x : ℝ // 0 < x}

/-- The eighth-order exponential Taylor tail on the positive quadrant. -/
noncomputable def positiveEighthOrderTailKernel :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun t ω => eighthOrderTail ((t : ℝ) * (ω : ℝ))

/--
The positive-quadrant eighth-order tail is a genuinely strictly sign-regular
kernel, rather than merely a family of separately stated ordered-minor facts.
-/
theorem positiveEighthOrderTailKernel_strictSignRegular :
    StrictSignRegular positiveEighthOrderTailKernel expKernelSignature := by
  intro order
  cases order with
  | zero =>
      intro rows columns _rows_mono _columns_mono
      simp [expKernelSignature]
  | succ n =>
      intro rows columns rows_mono columns_mono
      have rows_mono_real : StrictMono (fun i => (rows i : ℝ)) := by
        intro i j hij
        exact rows_mono hij
      have columns_mono_real : StrictMono (fun j => (columns j : ℝ)) := by
        intro i j hij
        exact columns_mono hij
      change 0 < expKernelSignature (n + 1) *
        Matrix.det (fun i j => eighthOrderTail
          ((rows i : ℝ) * (columns j : ℝ)))
      exact eighthOrderTail_minor_signed_pos
        (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ))
        rows_mono_real columns_mono_real
        (fun i => (rows i).property) (fun j => (columns j).property)

/-- The positive-frequency fermionic remainder obtained by denominator scaling. -/
noncomputable def positiveEighthOrderFermionicTailKernel :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun t ω => positiveEighthOrderTailKernel t ω *
    (1 + Real.exp (-(ω : ℝ)))⁻¹

/-- The rank-eight local Taylor component of the positive-frequency fermionic kernel. -/
noncomputable def positiveEighthOrderFermionicTaylorApprox :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun t ω => expNegTaylor 8 ((t : ℝ) * (ω : ℝ)) *
    (1 + Real.exp (-(ω : ℝ)))⁻¹

/-- The local Taylor component is represented by exactly eight separated terms. -/
theorem positiveEighthOrderFermionicTaylorApprox_eq_evalSeparated
    (t ω : PositiveArgument) :
    positiveEighthOrderFermionicTaylorApprox t ω =
      evalSeparated (divideSeparatedRight (taylorSeparatedTerms 8))
        (t : ℝ) (ω : ℝ) := by
  rw [eval_divideSeparatedRight, eval_taylorSeparatedTerms]
  simp only [positiveEighthOrderFermionicTaylorApprox, div_eq_mul_inv]

/-- Exact sign-compatible rank-eight-plus-tail decomposition on the positive quadrant. -/
theorem fermionicKernel_eq_positiveEighthOrderTaylorApprox_add_tail
    (t ω : PositiveArgument) :
    fermionicKernel (t : ℝ) (ω : ℝ) =
      positiveEighthOrderFermionicTaylorApprox t ω +
        positiveEighthOrderFermionicTailKernel t ω := by
  simp only [fermionicKernel, positiveEighthOrderFermionicTaylorApprox,
    positiveEighthOrderFermionicTailKernel,
    positiveEighthOrderTailKernel, eighthOrderTail, div_eq_mul_inv]
  ring_nf

/-- Positive fermionic denominator scaling preserves every strict tail-minor sign. -/
theorem positiveEighthOrderFermionicTailKernel_strictSignRegular :
    StrictSignRegular positiveEighthOrderFermionicTailKernel expKernelSignature := by
  intro order
  have scaled := strictSignRegularAtOrder_columnScale
    (positiveEighthOrderTailKernel_strictSignRegular order)
    (fun ω : PositiveArgument => (1 + Real.exp (-(ω : ℝ)))⁻¹)
    (fun ω => inv_pos.mpr (by positivity))
  change StrictSignRegularAtOrder
    (fun t ω : PositiveArgument => positiveEighthOrderTailKernel t ω *
      (1 + Real.exp (-(ω : ℝ)))⁻¹)
    order (expKernelSignature order)
  exact scaled

/-- The eighth-order Taylor tail is at most `1 / 8!` on the positive unit square. -/
theorem positiveEighthOrderTailKernel_abs_le_inv_factorial
    {t ω : PositiveArgument} (ht : (t : ℝ) ≤ 1) (hω : (ω : ℝ) ≤ 1) :
    abs (positiveEighthOrderTailKernel t ω) ≤
      1 / (Nat.factorial 8 : ℝ) := by
  have product_nonneg : 0 ≤ (t : ℝ) * (ω : ℝ) :=
    mul_nonneg t.property.le ω.property.le
  have product_le_one : (t : ℝ) * (ω : ℝ) ≤ 1 := by
    calc
      (t : ℝ) * (ω : ℝ) ≤ 1 * 1 :=
        mul_le_mul ht hω ω.property.le (by norm_num)
      _ = 1 := by norm_num
  have remainder := expNegTaylor_error_le_next
    (terms := 8) (by norm_num) product_nonneg
  change abs (Real.exp (-((t : ℝ) * (ω : ℝ))) -
    expNegTaylor 8 ((t : ℝ) * (ω : ℝ))) ≤ _
  calc
    abs (Real.exp (-((t : ℝ) * (ω : ℝ))) -
        expNegTaylor 8 ((t : ℝ) * (ω : ℝ))) ≤
        ((t : ℝ) * (ω : ℝ)) ^ 8 / (Nat.factorial 8 : ℝ) := remainder
    _ ≤ 1 / (Nat.factorial 8 : ℝ) := by
      gcongr
      exact pow_le_one₀ product_nonneg product_le_one

/-- The denominator-scaled fermionic tail obeys the same unit-square envelope. -/
theorem positiveEighthOrderFermionicTailKernel_abs_le_inv_factorial
    {t ω : PositiveArgument} (ht : (t : ℝ) ≤ 1) (hω : (ω : ℝ) ≤ 1) :
    abs (positiveEighthOrderFermionicTailKernel t ω) ≤
      1 / (Nat.factorial 8 : ℝ) := by
  have tail_bound := positiveEighthOrderTailKernel_abs_le_inv_factorial ht hω
  have denominator_one : 1 ≤ 1 + Real.exp (-(ω : ℝ)) := by
    linarith [Real.exp_pos (-(ω : ℝ))]
  have inverse_le_one : abs ((1 + Real.exp (-(ω : ℝ)))⁻¹) ≤ 1 := by
    rw [abs_of_pos (inv_pos.mpr (by positivity))]
    exact (inv_le_one₀ (by positivity : 0 < 1 + Real.exp (-(ω : ℝ)))).2
      denominator_one
  rw [positiveEighthOrderFermionicTailKernel, abs_mul]
  calc
    abs (positiveEighthOrderTailKernel t ω) *
        abs ((1 + Real.exp (-(ω : ℝ)))⁻¹) ≤
        (1 / (Nat.factorial 8 : ℝ)) * 1 := by
      gcongr
    _ = 1 / (Nat.factorial 8 : ℝ) := by ring

/-- The explicit eight-term fermionic approximation has the same unit-square error bound. -/
theorem positiveEighthOrderFermionicTaylorApprox_error_le_inv_factorial
    {t ω : PositiveArgument} (ht : (t : ℝ) ≤ 1) (hω : (ω : ℝ) ≤ 1) :
    abs (fermionicKernel (t : ℝ) (ω : ℝ) -
      positiveEighthOrderFermionicTaylorApprox t ω) ≤
        1 / (Nat.factorial 8 : ℝ) := by
  rw [fermionicKernel_eq_positiveEighthOrderTaylorApprox_add_tail]
  ring_nf
  exact positiveEighthOrderFermionicTailKernel_abs_le_inv_factorial ht hω

/--
Every complete GECP run on the positive unit-band Taylor tail has a
factorial-free selected-core determinant envelope `(1 / 8!)^rank`.
-/
theorem positiveEighthOrderTailKernel_completeRun_selectedCore_det_le
    (run : Run positiveEighthOrderTailKernel)
    (complete : run.CompleteOn
      (fun t : PositiveArgument => (t : ℝ) ≤ 1)
      (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)) :
    abs run.selectedCore.det ≤
      (1 / (Nat.factorial 8 : ℝ)) ^ Fintype.card run.SelectedIndex := by
  exact strictSignRegular_completeContinuation_selectedCore_det_le_pow
    positiveEighthOrderTailKernel_strictSignRegular
    (fun t : PositiveArgument => (t : ℝ) ≤ 1)
    (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)
    (Run.nil positiveEighthOrderTailKernel) run
    (fun _t _ω ht hω =>
      positiveEighthOrderTailKernel_abs_le_inv_factorial ht hω)
    complete

/-- The denominator-scaled positive-frequency fermionic tail has the same GECP envelope. -/
theorem positiveEighthOrderFermionicTailKernel_completeRun_selectedCore_det_le
    (run : Run positiveEighthOrderFermionicTailKernel)
    (complete : run.CompleteOn
      (fun t : PositiveArgument => (t : ℝ) ≤ 1)
      (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)) :
    abs run.selectedCore.det ≤
      (1 / (Nat.factorial 8 : ℝ)) ^ Fintype.card run.SelectedIndex := by
  exact strictSignRegular_completeContinuation_selectedCore_det_le_pow
    positiveEighthOrderFermionicTailKernel_strictSignRegular
    (fun t : PositiveArgument => (t : ℝ) ≤ 1)
    (fun ω : PositiveArgument => (ω : ℝ) ≤ 1)
    (Run.nil positiveEighthOrderFermionicTailKernel) run
    (fun _t _ω ht hω =>
      positiveEighthOrderFermionicTailKernel_abs_le_inv_factorial ht hω)
    complete

end Fermionic
end GECPKernelStructure
