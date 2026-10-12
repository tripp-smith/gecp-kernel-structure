import GECPKernelStructure.Fermionic.PositiveSeriesGECP

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The affine positive-coordinate form of the fermionic kernel, with
`z = 1 - omega`. -/
noncomputable def signedUnitShiftedKernel (t z : ℝ) : ℝ :=
  Real.exp (-t) * Real.exp (t * z) * (1 + Real.exp (z - 1))⁻¹

/-- The affine frequency shift identifies the full signed unit band with a
positive exponential product band of width two. -/
theorem fermionicKernel_eq_signedUnitShiftedKernel (t omega : ℝ) :
    fermionicKernel t omega = signedUnitShiftedKernel t (1 - omega) := by
  simp only [fermionicKernel, signedUnitShiftedKernel, div_eq_mul_inv]
  rw [← Real.exp_add]
  congr 1
  · ring_nf
  · ring_nf

/-- A finite sample of the shifted fermionic kernel. -/
noncomputable def signedUnitShiftedMatrix {order : ℕ}
    (rows columns : Fin order → ℝ) : Matrix (Fin order) (Fin order) ℝ :=
  fun i j => signedUnitShiftedKernel (rows i) (columns j)

/-- The shifted sample with its eighth-order positive exponential tail
multiplied by `lambda`. -/
noncomputable def signedUnitShiftedAmplifiedMatrix {order : ℕ}
    (lambda : ℝ) (rows columns : Fin order → ℝ) :
    Matrix (Fin order) (Fin order) ℝ :=
  fun i j =>
    Real.exp (-rows i) *
      (lambda * positiveEighthOrderTail (rows i * columns j) +
        expNegTaylor 8 (-(rows i * columns j))) *
      (1 + Real.exp (columns j - 1))⁻¹

private theorem signedUnitShiftedMatrix_eq_scaled_expMatrix
    {order : ℕ} (rows columns : Fin order → ℝ) :
    signedUnitShiftedMatrix rows columns =
      fun i j => Real.exp (-rows i) *
        ((1 + Real.exp (columns j - 1))⁻¹ * expMatrix rows columns i j) := by
  ext i j
  simp only [signedUnitShiftedMatrix, signedUnitShiftedKernel, expMatrix]
  ring_nf

private theorem signedUnitShiftedAmplifiedMatrix_eq_scaled_raw
    {order : ℕ} (lambda : ℝ) (rows columns : Fin order → ℝ) :
    signedUnitShiftedAmplifiedMatrix lambda rows columns =
      fun i j => Real.exp (-rows i) *
        ((1 + Real.exp (columns j - 1))⁻¹ *
          (lambda • positiveEighthOrderTailMatrix rows columns +
            positiveEighthOrderTaylorMatrix rows columns) i j) := by
  have taylor_apply (i j : Fin order) :
      positiveEighthOrderTaylorMatrix rows columns i j =
        expNegTaylor 8 (-(rows i * columns j)) := by
    simp only [positiveEighthOrderTaylorMatrix, expTaylorMatrix_apply]
    rw [expNegTaylor, ← Fin.sum_univ_eq_sum_range]
    apply Finset.sum_congr rfl
    intro k _
    field_simp
    ring_nf
  ext i j
  simp only [signedUnitShiftedAmplifiedMatrix, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul, positiveEighthOrderTailMatrix_apply,
    taylor_apply]
  ring

private theorem signedUnitShiftedMatrix_det
    {order : ℕ} (rows columns : Fin order → ℝ) :
    (signedUnitShiftedMatrix rows columns).det =
      (∏ i, Real.exp (-rows i)) *
        ((∏ j, (1 + Real.exp (columns j - 1))⁻¹) *
          (expMatrix rows columns).det) := by
  rw [signedUnitShiftedMatrix_eq_scaled_expMatrix]
  let columnScaled : Matrix (Fin order) (Fin order) ℝ :=
    fun i j => (1 + Real.exp (columns j - 1))⁻¹ * expMatrix rows columns i j
  have column_det : columnScaled.det =
      (∏ j, (1 + Real.exp (columns j - 1))⁻¹) *
        (expMatrix rows columns).det := by
    exact Matrix.det_mul_row
      (fun j => (1 + Real.exp (columns j - 1))⁻¹) (expMatrix rows columns)
  change (Matrix.det fun i j => Real.exp (-rows i) * columnScaled i j) = _
  calc
    (Matrix.det fun i j => Real.exp (-rows i) * columnScaled i j) =
        (∏ i, Real.exp (-rows i)) * columnScaled.det :=
      Matrix.det_mul_column (fun i => Real.exp (-rows i)) columnScaled
    _ = _ := by rw [column_det]

private theorem signedUnitShiftedAmplifiedMatrix_det
    {order : ℕ} (lambda : ℝ) (rows columns : Fin order → ℝ) :
    (signedUnitShiftedAmplifiedMatrix lambda rows columns).det =
      (∏ i, Real.exp (-rows i)) *
        ((∏ j, (1 + Real.exp (columns j - 1))⁻¹) *
          (lambda • positiveEighthOrderTailMatrix rows columns +
            positiveEighthOrderTaylorMatrix rows columns).det) := by
  rw [signedUnitShiftedAmplifiedMatrix_eq_scaled_raw]
  let raw := lambda • positiveEighthOrderTailMatrix rows columns +
    positiveEighthOrderTaylorMatrix rows columns
  let columnScaled : Matrix (Fin order) (Fin order) ℝ :=
    fun i j => (1 + Real.exp (columns j - 1))⁻¹ * raw i j
  have column_det : columnScaled.det =
      (∏ j, (1 + Real.exp (columns j - 1))⁻¹) * raw.det := by
    exact Matrix.det_mul_row
      (fun j => (1 + Real.exp (columns j - 1))⁻¹) raw
  change (Matrix.det fun i j => Real.exp (-rows i) * columnScaled i j) = _
  calc
    (Matrix.det fun i j => Real.exp (-rows i) * columnScaled i j) =
        (∏ i, Real.exp (-rows i)) * columnScaled.det :=
      Matrix.det_mul_column (fun i => Real.exp (-rows i)) columnScaled
    _ = _ := by rw [column_det]

/-- Tail amplification on the shifted signed band preserves strict positivity
of every ordered finite sample, including samples at the closed endpoints. -/
theorem signedUnitShiftedAmplifiedMatrix_det_pos
    {order : ℕ} {lambda : ℝ} (lambda_one_le : 1 ≤ lambda)
    (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ j, 0 ≤ columns j)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 < (signedUnitShiftedAmplifiedMatrix lambda rows columns).det := by
  let tail := positiveEighthOrderTailMatrix rows columns
  let low := positiveEighthOrderTaylorMatrix rows columns
  have coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree tail low := by
    intro degree
    exact positiveEighthOrderTail_mixedColumnCoefficient_nonneg
      rows columns rows_nonneg columns_nonneg rows_mono columns_mono degree
  have raw_mono :
      ((1 : ℝ) • tail + low).det ≤ (lambda • tail + low).det :=
    mixedColumnPolynomial_det_mono tail low coefficient_nonneg
      (by norm_num) lambda_one_le
  have split : tail + low = expMatrix rows columns := by
    simp [tail, low, positiveEighthOrderTailMatrix]
  have raw_pos : 0 < (lambda • tail + low).det := by
    have exp_pos := expMatrix_det_pos rows columns rows_nonneg columns_nonneg
      rows_mono columns_mono
    rw [one_smul, split] at raw_mono
    exact exp_pos.trans_le raw_mono
  rw [signedUnitShiftedAmplifiedMatrix_det]
  exact mul_pos
    (Finset.prod_pos fun i _ => Real.exp_pos _)
    (mul_pos
      (Finset.prod_pos fun j _ => inv_pos.mpr (by positivity)) raw_pos)

@[simp]
theorem signedUnitShiftedAmplifiedMatrix_one
    {order : ℕ} (rows columns : Fin order → ℝ) :
    signedUnitShiftedAmplifiedMatrix 1 rows columns =
      signedUnitShiftedMatrix rows columns := by
  ext i j
  simp only [signedUnitShiftedAmplifiedMatrix, signedUnitShiftedMatrix,
    signedUnitShiftedKernel, one_mul, positiveEighthOrderTail,
    eighthOrderTail]
  ring_nf

/-- The unamplified shifted matrix has positive determinant on every strictly
ordered nonnegative sample, including the closed endpoints. -/
theorem signedUnitShiftedMatrix_det_pos
    {order : ℕ} (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ j, 0 ≤ columns j)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 < (signedUnitShiftedMatrix rows columns).det := by
  rw [← signedUnitShiftedAmplifiedMatrix_one]
  exact signedUnitShiftedAmplifiedMatrix_det_pos (by norm_num)
    rows columns rows_nonneg columns_nonneg rows_mono columns_mono

private theorem signedUnitShiftedKernel_eq_fermionicKernel (t z : ℝ) :
    signedUnitShiftedKernel t z = fermionicKernel t (1 - z) := by
  simpa using (fermionicKernel_eq_signedUnitShiftedKernel t (1 - z)).symm

private theorem signedUnitShiftedTail_abs_le
    {t z : ℝ} (t_nonneg : 0 ≤ t) (t_le_one : t ≤ 1)
    (z_nonneg : 0 ≤ z) (z_le_two : z ≤ 2) :
    abs (Real.exp (-t) * positiveEighthOrderTail (t * z) *
      (1 + Real.exp (z - 1))⁻¹) ≤
        (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
  let x := t * z
  let denominator := 1 + Real.exp (z - 1)
  have x_nonneg : 0 ≤ x := mul_nonneg t_nonneg z_nonneg
  have x_le_two : x ≤ 2 := by
    dsimp only [x]
    nlinarith [mul_le_mul_of_nonneg_right t_le_one z_nonneg]
  have denominator_pos : 0 < denominator := by
    dsimp only [denominator]
    positivity
  have weighted_exp_le_denominator :
      Real.exp (-t) * Real.exp x ≤ denominator := by
    rw [← Real.exp_add]
    have exponent_eq : -t + x = t * (z - 1) := by
      dsimp only [x]
      ring
    rw [exponent_eq]
    by_cases hz : z ≤ 1
    · have exponent_nonpos : t * (z - 1) ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos t_nonneg (by linarith)
      exact (Real.exp_le_one_iff.mpr exponent_nonpos).trans
        (le_add_of_nonneg_right (Real.exp_pos (z - 1)).le)
    · have one_le_z : 1 ≤ z := le_of_not_ge hz
      have exponent_le : t * (z - 1) ≤ z - 1 := by
        nlinarith [mul_le_mul_of_nonneg_right t_le_one (sub_nonneg.mpr one_le_z)]
      exact (Real.exp_le_exp.mpr exponent_le).trans
        (le_add_of_nonneg_left (by norm_num))
  have remainder := positiveExpTaylor_error_le_next_mul_exp
    (terms := 8) (by norm_num) x_nonneg
  have tail_remainder :
      abs (positiveEighthOrderTail x) ≤
        Real.exp x * x ^ 8 / (Nat.factorial 8 : ℝ) := by
    simpa [positiveEighthOrderTail, eighthOrderTail] using remainder
  calc
    abs (Real.exp (-t) * positiveEighthOrderTail (t * z) *
        (1 + Real.exp (z - 1))⁻¹) =
        Real.exp (-t) * abs (positiveEighthOrderTail x) / denominator := by
      simp only [x, denominator, abs_mul, abs_of_pos (Real.exp_pos _),
        abs_inv, abs_of_pos denominator_pos, div_eq_mul_inv]
    _ ≤ Real.exp (-t) *
        (Real.exp x * x ^ 8 / (Nat.factorial 8 : ℝ)) / denominator := by
      gcongr
    _ = (Real.exp (-t) * Real.exp x / denominator) *
        (x ^ 8 / (Nat.factorial 8 : ℝ)) := by ring
    _ ≤ 1 * (x ^ 8 / (Nat.factorial 8 : ℝ)) := by
      gcongr
      exact (div_le_one denominator_pos).2 weighted_exp_le_denominator
    _ ≤ (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
      rw [one_mul]
      apply div_le_div_of_nonneg_right
      exact pow_le_pow_left₀ x_nonneg x_le_two 8
      positivity

/-- The weight-`128` amplified shifted kernel stays inside the fixed envelope
`2` on `0 <= t <= 1`, `0 <= z <= 2`. -/
theorem signedUnitShiftedAmplifiedMatrix_one_two_abs_le_two
    {order : ℕ} (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2) (i j : Fin order) :
    abs (signedUnitShiftedAmplifiedMatrix 128 rows columns i j) ≤ 2 := by
  let t := rows i
  let z := columns j
  let tail := Real.exp (-t) * positiveEighthOrderTail (t * z) *
    (1 + Real.exp (z - 1))⁻¹
  have decomposition :
      signedUnitShiftedAmplifiedMatrix 128 rows columns i j =
        signedUnitShiftedKernel t z + 127 * tail := by
    simp only [signedUnitShiftedAmplifiedMatrix, signedUnitShiftedKernel,
      tail, t, z, positiveEighthOrderTail, eighthOrderTail]
    ring_nf
  have base_abs_le_one : abs (signedUnitShiftedKernel t z) ≤ 1 := by
    rw [signedUnitShiftedKernel_eq_fermionicKernel,
      abs_of_pos (fermionicKernel_pos t (1 - z))]
    exact fermionicKernel_le_one (rows_nonneg i) (rows_le_one i)
  have tail_abs_le : abs tail ≤
      (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
    exact signedUnitShiftedTail_abs_le (rows_nonneg i) (rows_le_one i)
      (columns_nonneg j) (columns_le_two j)
  rw [decomposition]
  calc
    abs (signedUnitShiftedKernel t z + 127 * tail) ≤
        abs (signedUnitShiftedKernel t z) + abs (127 * tail) := abs_add_le _ _
    _ ≤ 1 + 127 * ((2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ)) := by
      rw [abs_mul, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 127)]
      gcongr
    _ ≤ 2 := by norm_num [Nat.factorial]

/-- Every ordered shifted-band amplified sample is strictly totally positive. -/
theorem signedUnitShiftedAmplifiedMatrix_strictlyTotallyPositive
    {order : ℕ} {lambda : ℝ} (lambda_one_le : 1 ≤ lambda)
    (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ j, 0 ≤ columns j)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    Matrix.StrictlyTotallyPositive
      (signedUnitShiftedAmplifiedMatrix lambda rows columns) := by
  intro m selectedRows selectedColumns selectedRows_mono selectedColumns_mono
  change 0 < (signedUnitShiftedAmplifiedMatrix lambda
    (rows ∘ selectedRows) (columns ∘ selectedColumns)).det
  exact signedUnitShiftedAmplifiedMatrix_det_pos lambda_one_le
    (rows ∘ selectedRows) (columns ∘ selectedColumns)
    (fun i => rows_nonneg (selectedRows i))
    (fun j => columns_nonneg (selectedColumns j))
    (rows_mono.comp selectedRows_mono) (columns_mono.comp selectedColumns_mono)

/-- Ordered shifted-band determinants acquire one factor `128⁻¹` per sample
beyond the first eight, up to the fixed amplified envelope. -/
theorem signedUnitShiftedMatrix_det_le_two_pow_div_128_pow
    {order : ℕ} (rows columns : Fin order → ℝ) (eight_le : 8 ≤ order)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    (signedUnitShiftedMatrix rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let weight : ℝ :=
    (∏ i, Real.exp (-rows i)) *
      ∏ j, (1 + Real.exp (columns j - 1))⁻¹
  have weight_pos : 0 < weight := by
    dsimp only [weight]
    exact mul_pos (Finset.prod_pos fun i _ => Real.exp_pos _)
      (Finset.prod_pos fun j _ => inv_pos.mpr (by positivity))
  have rank_amplification := positiveEighthOrderTail_rank_amplification
    rows columns eight_le rows_nonneg columns_nonneg rows_mono columns_mono
    (lambda := (128 : ℝ)) (by norm_num)
  have shifted_rank_amplification :
      (128 : ℝ) ^ (order - 8) *
          (signedUnitShiftedMatrix rows columns).det ≤
        (signedUnitShiftedAmplifiedMatrix 128 rows columns).det := by
    rw [signedUnitShiftedMatrix_det,
      signedUnitShiftedAmplifiedMatrix_det]
    calc
      (128 : ℝ) ^ (order - 8) *
          ((∏ i, Real.exp (-rows i)) *
            ((∏ j, (1 + Real.exp (columns j - 1))⁻¹) *
              (expMatrix rows columns).det)) =
        weight * ((128 : ℝ) ^ (order - 8) *
          (expMatrix rows columns).det) := by
            dsimp only [weight]
            ring
      _ ≤ weight *
          ((128 : ℝ) • positiveEighthOrderTailMatrix rows columns +
            positiveEighthOrderTaylorMatrix rows columns).det :=
        mul_le_mul_of_nonneg_left rank_amplification weight_pos.le
      _ = (∏ i, Real.exp (-rows i)) *
          ((∏ j, (1 + Real.exp (columns j - 1))⁻¹) *
            ((128 : ℝ) • positiveEighthOrderTailMatrix rows columns +
              positiveEighthOrderTaylorMatrix rows columns).det) := by
        dsimp only [weight]
        ring
  have positive := signedUnitShiftedAmplifiedMatrix_strictlyTotallyPositive
    (lambda := (128 : ℝ)) (by norm_num) rows columns rows_nonneg
    columns_nonneg rows_mono columns_mono
  have amplified_bound :
      (signedUnitShiftedAmplifiedMatrix 128 rows columns).det ≤
        (2 : ℝ) ^ order := by
    have hadamard :=
      Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive
        (signedUnitShiftedAmplifiedMatrix 128 rows columns) positive
    have diagonal_bound :
        (∏ i, signedUnitShiftedAmplifiedMatrix 128 rows columns i i) ≤
          ∏ _i : Fin order, (2 : ℝ) := by
      apply Finset.prod_le_prod
      · intro i _
        exact (positive.entry_pos i i).le
      · intro i _
        exact (le_abs_self _).trans
          (signedUnitShiftedAmplifiedMatrix_one_two_abs_le_two
            rows columns rows_nonneg rows_le_one columns_nonneg columns_le_two i i)
    exact hadamard.trans (by simpa using diagonal_bound)
  apply (le_div_iff₀ (pow_pos (by norm_num) (order - 8))).2
  simpa only [mul_comm] using shifted_rank_amplification.trans amplified_bound

/-- The signed-band determinant envelope is invariant under independent
permutations of injective row and shifted-frequency samples. -/
theorem signedUnitShiftedMatrix_abs_det_le_two_pow_div_128_pow
    {order : ℕ} (rows columns : Fin order → ℝ) (eight_le : 8 ≤ order)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2) :
    abs (signedUnitShiftedMatrix rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let rowPerm := Tuple.sort rows
  let columnPerm := Tuple.sort columns
  have rows_sorted : StrictMono (rows ∘ rowPerm) :=
    (Tuple.monotone_sort rows).strictMono_of_injective
      (rows_injective.comp rowPerm.injective)
  have columns_sorted : StrictMono (columns ∘ columnPerm) :=
    (Tuple.monotone_sort columns).strictMono_of_injective
      (columns_injective.comp columnPerm.injective)
  let original := signedUnitShiftedMatrix rows columns
  have sorted_matrix :
      signedUnitShiftedMatrix (rows ∘ rowPerm) (columns ∘ columnPerm) =
        original.submatrix rowPerm columnPerm := by
    rfl
  have sorted_pos := signedUnitShiftedMatrix_det_pos
    (rows ∘ rowPerm) (columns ∘ columnPerm)
    (fun i => rows_nonneg (rowPerm i))
    (fun j => columns_nonneg (columnPerm j)) rows_sorted columns_sorted
  have sorted_bound := signedUnitShiftedMatrix_det_le_two_pow_div_128_pow
    (rows ∘ rowPerm) (columns ∘ columnPerm) eight_le
    (fun i => rows_nonneg (rowPerm i)) (fun i => rows_le_one (rowPerm i))
    (fun j => columns_nonneg (columnPerm j))
    (fun j => columns_le_two (columnPerm j)) rows_sorted columns_sorted
  calc
    abs (signedUnitShiftedMatrix rows columns).det = abs original.det := rfl
    _ = abs (original.submatrix rowPerm columnPerm).det := by
      symm
      exact Matrix.abs_det_submatrix_equiv_equiv rowPerm columnPerm original
    _ = abs (signedUnitShiftedMatrix
        (rows ∘ rowPerm) (columns ∘ columnPerm)).det := by rw [sorted_matrix]
    _ = (signedUnitShiftedMatrix
        (rows ∘ rowPerm) (columns ∘ columnPerm)).det := abs_of_pos sorted_pos
    _ ≤ (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := sorted_bound

/-- Every injective fermionic sample on the full signed unit band has the
fixed-base determinant decay exposed by the affine frequency shift. -/
theorem fermionicKernel_signedUnit_abs_sample_det_le
    {order : ℕ} (rows columns : Fin order → ℝ) (eight_le : 8 ≤ order)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_neg_one_le : ∀ j, -1 ≤ columns j)
    (columns_le_one : ∀ j, columns j ≤ 1) :
    abs (sampleMatrix fermionicKernel rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let shiftedColumns : Fin order → ℝ := fun j => 1 - columns j
  have shifted_injective : Function.Injective shiftedColumns := by
    intro i j hij
    apply columns_injective
    dsimp only [shiftedColumns] at hij
    linarith
  have shifted_nonneg : ∀ j, 0 ≤ shiftedColumns j := by
    intro j
    dsimp only [shiftedColumns]
    linarith [columns_le_one j]
  have shifted_le_two : ∀ j, shiftedColumns j ≤ 2 := by
    intro j
    dsimp only [shiftedColumns]
    linarith [columns_neg_one_le j]
  have matrix_eq :
      sampleMatrix fermionicKernel rows columns =
        signedUnitShiftedMatrix rows shiftedColumns := by
    ext i j
    exact fermionicKernel_eq_signedUnitShiftedKernel (rows i) (columns j)
  rw [matrix_eq]
  exact signedUnitShiftedMatrix_abs_det_le_two_pow_div_128_pow
    rows shiftedColumns eight_le rows_injective shifted_injective
    rows_nonneg rows_le_one shifted_nonneg shifted_le_two

/-- Complete-pivot residual powers on the full signed unit band inherit the
fixed-base determinant decay from the affine shifted amplification. -/
theorem fermionicKernel_signedUnit_gecp_error_pow_le
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ m, ∃ run : Run fermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -1 ≤ omega ∧ omega ≤ 1)
      (residual m) (rows m) (columns m))
    (n : ℕ) (eight_le : 8 ≤ n) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (run_complete : run.CompleteOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -1 ≤ omega ∧ omega ≤ 1))
    (t omega : ℝ) (t_mem : 0 ≤ t ∧ t ≤ 1)
    (omega_mem : -1 ≤ omega ∧ omega ≤ 1) :
    abs (residual n t omega) ^ n ≤
      (2 : ℝ) ^ n / (128 : ℝ) ^ (n - 8) := by
  have residual_power := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun omega : ℝ => -1 ≤ omega ∧ omega ≤ 1)
    residual rows columns realized pivot_ne updates complete n run run_pivots
    t omega t_mem omega_mem
  have pivot_length : run.pivots.length = n := by
    rw [run_pivots]
    simp
  let finiteRows : Fin run.pivots.length → ℝ :=
    run.selectedRow ∘ run.selectedIndexEquivFin.symm
  let finiteColumns : Fin run.pivots.length → ℝ :=
    run.selectedColumn ∘ run.selectedIndexEquivFin.symm
  have finiteRows_injective : Function.Injective finiteRows :=
    run.selectedRow_injective.comp run.selectedIndexEquivFin.symm.injective
  have finiteColumns_injective : Function.Injective finiteColumns :=
    run.selectedColumn_injective.comp run.selectedIndexEquivFin.symm.injective
  have finiteRows_nonneg : ∀ i, 0 ≤ finiteRows i := by
    intro i
    exact (run.selectedRow_mem_of_completeOn run_complete _).1
  have finiteRows_le_one : ∀ i, finiteRows i ≤ 1 := by
    intro i
    exact (run.selectedRow_mem_of_completeOn run_complete _).2
  have finiteColumns_neg_one_le : ∀ j, -1 ≤ finiteColumns j := by
    intro j
    exact (run.selectedColumn_mem_of_completeOn run_complete _).1
  have finiteColumns_le_one : ∀ j, finiteColumns j ≤ 1 := by
    intro j
    exact (run.selectedColumn_mem_of_completeOn run_complete _).2
  have eight_le_length : 8 ≤ run.pivots.length := by
    simpa [pivot_length] using eight_le
  have determinant_bound := fermionicKernel_signedUnit_abs_sample_det_le
    finiteRows finiteColumns eight_le_length finiteRows_injective
    finiteColumns_injective finiteRows_nonneg finiteRows_le_one
    finiteColumns_neg_one_le finiteColumns_le_one
  have finite_sample :
      sampleMatrix fermionicKernel finiteRows finiteColumns =
        run.finSelectedCore := by
    rfl
  rw [finite_sample] at determinant_bound
  exact residual_power.trans (by simpa [pivot_length] using determinant_bound)

/-- Twelve exact complete pivots contract every full signed-unit-band residual
entry by at least one half. -/
theorem fermionicKernel_signedUnit_gecp_error_le_half_twelve
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (realized : ∀ m, ∃ run : Run fermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -1 ≤ omega ∧ omega ≤ 1)
      (residual m) (rows m) (columns m))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns 12)
    (run_complete : run.CompleteOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -1 ≤ omega ∧ omega ≤ 1))
    (t omega : ℝ) (t_mem : 0 ≤ t ∧ t ≤ 1)
    (omega_mem : -1 ≤ omega ∧ omega ≤ 1) :
    abs (residual 12 t omega) ≤ (1 / 2 : ℝ) := by
  have power_bound := fermionicKernel_signedUnit_gecp_error_pow_le
    residual rows columns realized pivot_ne updates complete 12 (by norm_num)
    run run_pivots run_complete t omega t_mem omega_mem
  apply le_of_pow_le_pow_left₀ (by norm_num : 12 ≠ 0)
    (by norm_num : (0 : ℝ) ≤ 1 / 2)
  exact power_bound.trans (by norm_num)

end Fermionic
end GECPKernelStructure
