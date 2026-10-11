import GECPKernelStructure.Fermionic.TailAmplification
import GECPKernelStructure.Fermionic.ContinuationDeterminant
import GECPKernelStructure.Matrix.TotalPositiveHadamard

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The positive-series fermionic kernel with its eighth-order tail multiplied
by `lambda`. -/
noncomputable def positiveSeriesAmplifiedFermionicKernel (lambda : ℝ) :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun u omega =>
    (lambda * positiveEighthOrderTail ((u : ℝ) * (omega : ℝ)) +
      expNegTaylor 8 (-((u : ℝ) * (omega : ℝ)))) *
        (1 + Real.exp (omega : ℝ))⁻¹

/-- Sampling the amplified kernel gives the matrix pencil used in the rank
amplification theorem. -/
theorem positiveSeriesAmplifiedFermionicKernel_minorMatrix
    {order : ℕ} (lambda : ℝ)
    (rows columns : Fin order → PositiveArgument) :
    minorMatrix (positiveSeriesAmplifiedFermionicKernel lambda) rows columns =
      positiveSeriesAmplifiedFermionicMatrix lambda
        (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ)) := by
  ext i j
  simp only [minorMatrix, positiveSeriesAmplifiedFermionicKernel,
    positiveSeriesAmplifiedFermionicMatrix, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul, positiveEighthOrderTailMatrix_apply,
    positiveEighthOrderTaylorMatrix, expTaylorMatrix_apply]
  rw [expNegTaylor, ← Fin.sum_univ_eq_sum_range]
  congr 2
  apply Finset.sum_congr rfl
  intro k _
  field_simp
  ring

/-- Every nonnegative tail amplification above one preserves strict total
positivity of the positive-series fermionic kernel. -/
theorem positiveSeriesAmplifiedFermionicKernel_strictSignRegular
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    StrictSignRegular (positiveSeriesAmplifiedFermionicKernel lambda)
      (fun _ => 1) := by
  intro order rows columns rows_mono columns_mono
  let realRows : Fin order → ℝ := fun i => (rows i : ℝ)
  let realColumns : Fin order → ℝ := fun j => (columns j : ℝ)
  have realRows_nonneg : ∀ i, 0 ≤ realRows i :=
    fun i => (rows i).property.le
  have realColumns_nonneg : ∀ j, 0 ≤ realColumns j :=
    fun j => (columns j).property.le
  have realRows_mono : StrictMono realRows := by
    intro i j hij
    exact rows_mono hij
  have realColumns_mono : StrictMono realColumns := by
    intro i j hij
    exact columns_mono hij
  let tail := positiveEighthOrderTailMatrix realRows realColumns
  let low := positiveEighthOrderTaylorMatrix realRows realColumns
  have coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree tail low := by
    intro degree
    exact positiveEighthOrderTail_mixedColumnCoefficient_nonneg
      realRows realColumns realRows_nonneg realColumns_nonneg
      realRows_mono realColumns_mono degree
  have raw_mono :
      ((1 : ℝ) • tail + low).det ≤ (lambda • tail + low).det :=
    mixedColumnPolynomial_det_mono tail low coefficient_nonneg
      (by norm_num) lambda_one_le
  have split : tail + low = expMatrix realRows realColumns := by
    simp [tail, low, positiveEighthOrderTailMatrix]
  have raw_pos : 0 < (lambda • tail + low).det := by
    have exp_pos := expMatrix_det_pos realRows realColumns
      realRows_nonneg realColumns_nonneg realRows_mono realColumns_mono
    rw [one_smul, split] at raw_mono
    exact exp_pos.trans_le raw_mono
  let weights : Fin order → ℝ :=
    fun j => (1 + Real.exp (realColumns j))⁻¹
  have weights_product_pos : 0 < ∏ j, weights j := by
    apply Finset.prod_pos
    intro j _
    exact inv_pos.mpr (by positivity)
  have sampled_det :
      (minorMatrix (positiveSeriesAmplifiedFermionicKernel lambda)
          rows columns).det =
        (∏ j, weights j) * (lambda • tail + low).det := by
    rw [positiveSeriesAmplifiedFermionicKernel_minorMatrix]
    rw [show positiveSeriesAmplifiedFermionicMatrix lambda realRows realColumns =
        fun i j => weights j * (lambda • tail + low) i j by
      ext i j
      simp only [positiveSeriesAmplifiedFermionicMatrix, weights, tail, low,
        realRows, realColumns]
      ring]
    exact Matrix.det_mul_row weights (lambda • tail + low)
  rw [sampled_det]
  simpa only [one_mul] using mul_pos weights_product_pos raw_pos

/-- At the balancing weight `8!`, the amplified positive-series kernel stays
inside the fixed envelope `2` on the positive unit band. -/
theorem positiveSeriesAmplifiedFermionicKernel_factorial_abs_le_two
    {u omega : PositiveArgument} (hu : (u : ℝ) ≤ 1)
    (homega : (omega : ℝ) ≤ 1) :
    abs (positiveSeriesAmplifiedFermionicKernel (Nat.factorial 8 : ℝ)
      u omega) ≤ 2 := by
  let x := (u : ℝ) * (omega : ℝ)
  let denominator := 1 + Real.exp (omega : ℝ)
  have x_le_omega : x ≤ (omega : ℝ) := by
    dsimp only [x]
    nlinarith [mul_le_mul_of_nonneg_right hu (omega.property.le)]
  have denominator_pos : 0 < denominator := by
    dsimp only [denominator]
    positivity
  have base_abs_le_one :
      abs (Real.exp x * denominator⁻¹) ≤ 1 := by
    rw [abs_mul, abs_of_pos (Real.exp_pos x),
      abs_of_pos (inv_pos.mpr denominator_pos)]
    rw [← div_eq_mul_inv]
    apply (div_le_one denominator_pos).2
    dsimp only [denominator]
    exact (Real.exp_le_exp.mpr x_le_omega).trans
      (le_add_of_nonneg_left (by norm_num))
  have tail_bound :=
    positiveSeriesEighthOrderFermionicTailKernel_abs_le_inv_factorial
      hu homega
  have decomposition :
      positiveSeriesAmplifiedFermionicKernel (Nat.factorial 8 : ℝ) u omega =
        Real.exp x * denominator⁻¹ +
          ((Nat.factorial 8 : ℝ) - 1) *
            positiveSeriesEighthOrderFermionicTailKernel u omega := by
    simp only [positiveSeriesAmplifiedFermionicKernel,
      positiveSeriesEighthOrderFermionicTailKernel,
      positiveSeriesEighthOrderTailKernel, positiveEighthOrderTail,
      eighthOrderTail, x, denominator]
    ring_nf
  rw [decomposition]
  calc
    abs (Real.exp x * denominator⁻¹ +
        ((Nat.factorial 8 : ℝ) - 1) *
          positiveSeriesEighthOrderFermionicTailKernel u omega) ≤
        abs (Real.exp x * denominator⁻¹) +
          abs (((Nat.factorial 8 : ℝ) - 1) *
            positiveSeriesEighthOrderFermionicTailKernel u omega) :=
      abs_add_le _ _
    _ ≤ 1 + 1 := by
      gcongr
      rw [abs_mul, abs_of_nonneg (by positivity :
        0 ≤ (Nat.factorial 8 : ℝ) - 1)]
      calc
        ((Nat.factorial 8 : ℝ) - 1) *
            abs (positiveSeriesEighthOrderFermionicTailKernel u omega) ≤
            (Nat.factorial 8 : ℝ) *
              (1 / (Nat.factorial 8 : ℝ)) := by
          norm_num at tail_bound ⊢
          nlinarith [abs_nonneg
            (positiveSeriesEighthOrderFermionicTailKernel u omega)]
        _ = 1 := by norm_num
    _ = 2 := by norm_num

/-- Complete pivoting on the `8!`-amplified kernel inherits a factorial-free
`2^rank` selected-core determinant envelope. -/
theorem positiveSeriesAmplifiedFermionicKernel_factorial_completeRun_det_le
    (run : Run (positiveSeriesAmplifiedFermionicKernel
      (Nat.factorial 8 : ℝ)))
    (complete : run.CompleteOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun omega : PositiveArgument => (omega : ℝ) ≤ 1)) :
    abs run.selectedCore.det ≤
      (2 : ℝ) ^ Fintype.card run.SelectedIndex := by
  exact strictSignRegular_completeContinuation_selectedCore_det_le_pow
    (positiveSeriesAmplifiedFermionicKernel_strictSignRegular (by norm_num))
    (fun u : PositiveArgument => (u : ℝ) ≤ 1)
    (fun omega : PositiveArgument => (omega : ℝ) ≤ 1)
    (Run.nil (positiveSeriesAmplifiedFermionicKernel
      (Nat.factorial 8 : ℝ))) run
    (fun _u _omega hu homega =>
      positiveSeriesAmplifiedFermionicKernel_factorial_abs_le_two hu homega)
    complete

/-- Every ordered unit-band sample of the `8!`-amplified kernel has
determinant at most `2^order`, independently of any pivot path. -/
theorem positiveSeriesAmplifiedFermionicKernel_factorial_minor_det_le_two_pow
    {order : ℕ} (rows columns : Fin order → PositiveArgument)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_le_one : ∀ i, (rows i : ℝ) ≤ 1)
    (columns_le_one : ∀ j, (columns j : ℝ) ≤ 1) :
    (minorMatrix
      (positiveSeriesAmplifiedFermionicKernel (Nat.factorial 8 : ℝ))
      rows columns).det ≤ (2 : ℝ) ^ order := by
  apply minorMatrix_det_le_pow_of_strictSignRegular_one
    (positiveSeriesAmplifiedFermionicKernel_strictSignRegular (by norm_num))
    rows columns rows_mono columns_mono
  intro i
  exact (le_abs_self _).trans
    (positiveSeriesAmplifiedFermionicKernel_factorial_abs_le_two
      (rows_le_one i) (columns_le_one i))

/-- Rank-eight amplification plus the total-positive Hadamard envelope gives
an explicit fixed-base determinant-decay comparison for every ordered
positive-series fermionic sample. -/
theorem positiveSeriesFermionic_factorial_amplified_det_le_two_pow
    {order : ℕ} (rows columns : Fin order → PositiveArgument)
    (eight_le : 8 ≤ order)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_le_one : ∀ i, (rows i : ℝ) ≤ 1)
    (columns_le_one : ∀ j, (columns j : ℝ) ≤ 1) :
    (Nat.factorial 8 : ℝ) ^ (order - 8) *
        (positiveSeriesFermionicMatrix
          (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ))).det ≤
      (2 : ℝ) ^ order := by
  have realRows_mono : StrictMono (fun i => (rows i : ℝ)) := by
    intro i j hij
    exact rows_mono hij
  have realColumns_mono : StrictMono (fun j => (columns j : ℝ)) := by
    intro i j hij
    exact columns_mono hij
  have amplification := positiveSeriesFermionic_rank_amplification
    (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ)) eight_le
    (fun i => (rows i).property.le) (fun j => (columns j).property.le)
    realRows_mono realColumns_mono (lambda := (Nat.factorial 8 : ℝ))
    (by norm_num)
  have amplified_bound :=
    positiveSeriesAmplifiedFermionicKernel_factorial_minor_det_le_two_pow
      rows columns rows_mono columns_mono rows_le_one columns_le_one
  rw [positiveSeriesAmplifiedFermionicKernel_minorMatrix] at amplified_bound
  exact amplification.trans amplified_bound

/-- Division by the positive amplification factor exposes the determinant's
fixed-base exponential decay. -/
theorem positiveSeriesFermionic_det_le_two_pow_div_factorial_pow
    {order : ℕ} (rows columns : Fin order → PositiveArgument)
    (eight_le : 8 ≤ order)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_le_one : ∀ i, (rows i : ℝ) ≤ 1)
    (columns_le_one : ∀ j, (columns j : ℝ) ≤ 1) :
    (positiveSeriesFermionicMatrix
        (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ))).det ≤
      (2 : ℝ) ^ order /
        (Nat.factorial 8 : ℝ) ^ (order - 8) := by
  apply (le_div_iff₀ (pow_pos (by norm_num) (order - 8))).2
  simpa only [mul_comm] using
    (positiveSeriesFermionic_factorial_amplified_det_le_two_pow
      rows columns eight_le rows_mono columns_mono rows_le_one columns_le_one)

end Fermionic
end GECPKernelStructure
