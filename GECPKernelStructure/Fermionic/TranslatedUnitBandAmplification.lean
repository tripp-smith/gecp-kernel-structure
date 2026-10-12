import GECPKernelStructure.Fermionic.SignedUnitBandAmplification

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The positive-coordinate fermionic kernel on the width-two band centered at
`center`, with `z = center + 1 - omega`. -/
noncomputable def translatedUnitShiftedKernel
    (center t z : ℝ) : ℝ :=
  Real.exp (-t * (center + 1)) * Real.exp (t * z) *
    (1 + Real.exp (z - (center + 1)))⁻¹

/-- Every translated width-two frequency band has the same positive
exponential coordinate interval `0 <= z <= 2`. -/
theorem fermionicKernel_eq_translatedUnitShiftedKernel
    (center t omega : ℝ) :
    fermionicKernel t omega =
      translatedUnitShiftedKernel center t (center + 1 - omega) := by
  simp only [fermionicKernel, translatedUnitShiftedKernel, div_eq_mul_inv]
  rw [← Real.exp_add]
  congr 1 <;> ring_nf

/-- A finite sample of the center-parametric shifted kernel. -/
noncomputable def translatedUnitShiftedMatrix {order : ℕ}
    (center : ℝ) (rows columns : Fin order → ℝ) :
    Matrix (Fin order) (Fin order) ℝ :=
  fun i j => translatedUnitShiftedKernel center (rows i) (columns j)

/-- The center-parametric shifted sample with its positive exponential tail
multiplied by `lambda`. -/
noncomputable def translatedUnitShiftedAmplifiedMatrix {order : ℕ}
    (center lambda : ℝ) (rows columns : Fin order → ℝ) :
    Matrix (Fin order) (Fin order) ℝ :=
  fun i j =>
    Real.exp (-rows i * (center + 1)) *
      (lambda * positiveEighthOrderTail (rows i * columns j) +
        expNegTaylor 8 (-(rows i * columns j))) *
      (1 + Real.exp (columns j - (center + 1)))⁻¹

private theorem translatedUnitShiftedMatrix_det
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ) :
    (translatedUnitShiftedMatrix center rows columns).det =
      (∏ i, Real.exp (-rows i * (center + 1))) *
        ((∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) *
          (expMatrix rows columns).det) := by
  let raw := expMatrix rows columns
  let columnScaled : Matrix (Fin order) (Fin order) ℝ :=
    fun i j =>
      (1 + Real.exp (columns j - (center + 1)))⁻¹ * raw i j
  have matrix_eq : translatedUnitShiftedMatrix center rows columns =
      fun i j => Real.exp (-rows i * (center + 1)) * columnScaled i j := by
    ext i j
    simp only [translatedUnitShiftedMatrix, translatedUnitShiftedKernel,
      columnScaled, raw, expMatrix]
    ring
  have column_det : columnScaled.det =
      (∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) * raw.det := by
    exact Matrix.det_mul_row
      (fun j => (1 + Real.exp (columns j - (center + 1)))⁻¹) raw
  rw [matrix_eq]
  calc
    (Matrix.det fun i j =>
        Real.exp (-rows i * (center + 1)) * columnScaled i j) =
        (∏ i, Real.exp (-rows i * (center + 1))) * columnScaled.det :=
      Matrix.det_mul_column
        (fun i => Real.exp (-rows i * (center + 1))) columnScaled
    _ = _ := by rw [column_det]

private theorem translatedUnitShiftedAmplifiedMatrix_det
    {order : ℕ} (center lambda : ℝ)
    (rows columns : Fin order → ℝ) :
    (translatedUnitShiftedAmplifiedMatrix center lambda rows columns).det =
      (∏ i, Real.exp (-rows i * (center + 1))) *
        ((∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) *
          (lambda • positiveEighthOrderTailMatrix rows columns +
            positiveEighthOrderTaylorMatrix rows columns).det) := by
  let raw := lambda • positiveEighthOrderTailMatrix rows columns +
    positiveEighthOrderTaylorMatrix rows columns
  let columnScaled : Matrix (Fin order) (Fin order) ℝ :=
    fun i j =>
      (1 + Real.exp (columns j - (center + 1)))⁻¹ * raw i j
  have taylor_apply (i j : Fin order) :
      positiveEighthOrderTaylorMatrix rows columns i j =
        expNegTaylor 8 (-(rows i * columns j)) := by
    simp only [positiveEighthOrderTaylorMatrix, expTaylorMatrix_apply]
    rw [expNegTaylor, ← Fin.sum_univ_eq_sum_range]
    apply Finset.sum_congr rfl
    intro k _
    field_simp
    ring_nf
  have matrix_eq :
      translatedUnitShiftedAmplifiedMatrix center lambda rows columns =
        fun i j => Real.exp (-rows i * (center + 1)) * columnScaled i j := by
    ext i j
    simp only [translatedUnitShiftedAmplifiedMatrix, columnScaled, raw,
      Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
      positiveEighthOrderTailMatrix_apply, taylor_apply]
    ring
  have column_det : columnScaled.det =
      (∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) * raw.det := by
    exact Matrix.det_mul_row
      (fun j => (1 + Real.exp (columns j - (center + 1)))⁻¹) raw
  rw [matrix_eq]
  calc
    (Matrix.det fun i j =>
        Real.exp (-rows i * (center + 1)) * columnScaled i j) =
        (∏ i, Real.exp (-rows i * (center + 1))) * columnScaled.det :=
      Matrix.det_mul_column
        (fun i => Real.exp (-rows i * (center + 1))) columnScaled
    _ = _ := by rw [column_det]

/-- Tail amplification preserves strict positivity on every translated
width-two band, uniformly in its center. -/
theorem translatedUnitShiftedAmplifiedMatrix_det_pos
    {order : ℕ} (center : ℝ) {lambda : ℝ} (lambda_one_le : 1 ≤ lambda)
    (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ j, 0 ≤ columns j)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 < (translatedUnitShiftedAmplifiedMatrix
      center lambda rows columns).det := by
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
  rw [translatedUnitShiftedAmplifiedMatrix_det]
  exact mul_pos
    (Finset.prod_pos fun i _ => Real.exp_pos _)
    (mul_pos
      (Finset.prod_pos fun j _ => inv_pos.mpr (by positivity)) raw_pos)

@[simp]
theorem translatedUnitShiftedAmplifiedMatrix_one
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ) :
    translatedUnitShiftedAmplifiedMatrix center 1 rows columns =
      translatedUnitShiftedMatrix center rows columns := by
  ext i j
  simp only [translatedUnitShiftedAmplifiedMatrix,
    translatedUnitShiftedMatrix, translatedUnitShiftedKernel, one_mul,
    positiveEighthOrderTail, eighthOrderTail]
  ring_nf

private theorem translatedUnitShiftedKernel_eq_fermionicKernel
    (center t z : ℝ) :
    translatedUnitShiftedKernel center t z =
      fermionicKernel t (center + 1 - z) := by
  simpa using
    (fermionicKernel_eq_translatedUnitShiftedKernel
      center t (center + 1 - z)).symm

private theorem translatedUnitShiftedTail_abs_le
    {center t z : ℝ} (t_nonneg : 0 ≤ t) (t_le_one : t ≤ 1)
    (z_nonneg : 0 ≤ z) (z_le_two : z ≤ 2) :
    abs (Real.exp (-t * (center + 1)) *
      positiveEighthOrderTail (t * z) *
      (1 + Real.exp (z - (center + 1)))⁻¹) ≤
        (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
  let x := t * z
  let denominator := 1 + Real.exp (z - (center + 1))
  have x_nonneg : 0 ≤ x := mul_nonneg t_nonneg z_nonneg
  have x_le_two : x ≤ 2 := by
    dsimp only [x]
    nlinarith [mul_le_mul_of_nonneg_right t_le_one z_nonneg]
  have denominator_pos : 0 < denominator := by
    dsimp only [denominator]
    positivity
  have weighted_exp_le_one :
      Real.exp (-t * (center + 1)) * Real.exp x / denominator ≤ 1 := by
    have kernel_le : fermionicKernel t (center + 1 - z) ≤ 1 :=
      fermionicKernel_le_one t_nonneg t_le_one
    rw [← translatedUnitShiftedKernel_eq_fermionicKernel] at kernel_le
    simpa [translatedUnitShiftedKernel, x, denominator, div_eq_mul_inv]
      using kernel_le
  have remainder := positiveExpTaylor_error_le_next_mul_exp
    (terms := 8) (by norm_num) x_nonneg
  have tail_remainder :
      abs (positiveEighthOrderTail x) ≤
        Real.exp x * x ^ 8 / (Nat.factorial 8 : ℝ) := by
    simpa [positiveEighthOrderTail, eighthOrderTail] using remainder
  calc
    abs (Real.exp (-t * (center + 1)) *
        positiveEighthOrderTail (t * z) *
        (1 + Real.exp (z - (center + 1)))⁻¹) =
        Real.exp (-t * (center + 1)) *
          abs (positiveEighthOrderTail x) / denominator := by
      simp only [x, denominator, abs_mul, abs_of_pos (Real.exp_pos _),
        abs_inv, abs_of_pos denominator_pos, div_eq_mul_inv]
    _ ≤ Real.exp (-t * (center + 1)) *
        (Real.exp x * x ^ 8 / (Nat.factorial 8 : ℝ)) / denominator := by
      gcongr
    _ = (Real.exp (-t * (center + 1)) * Real.exp x / denominator) *
        (x ^ 8 / (Nat.factorial 8 : ℝ)) := by ring
    _ ≤ 1 * (x ^ 8 / (Nat.factorial 8 : ℝ)) := by
      gcongr
    _ ≤ (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
      rw [one_mul]
      apply div_le_div_of_nonneg_right
      exact pow_le_pow_left₀ x_nonneg x_le_two 8
      positivity

/-- The weight-`128` amplified kernel stays inside the fixed envelope `2` on
every translated width-two band, independently of its center. -/
theorem translatedUnitShiftedAmplifiedMatrix_abs_le_two
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2) (i j : Fin order) :
    abs (translatedUnitShiftedAmplifiedMatrix
      center 128 rows columns i j) ≤ 2 := by
  let t := rows i
  let z := columns j
  let tail := Real.exp (-t * (center + 1)) *
    positiveEighthOrderTail (t * z) *
    (1 + Real.exp (z - (center + 1)))⁻¹
  have decomposition :
      translatedUnitShiftedAmplifiedMatrix center 128 rows columns i j =
        translatedUnitShiftedKernel center t z + 127 * tail := by
    simp only [translatedUnitShiftedAmplifiedMatrix,
      translatedUnitShiftedKernel, tail, t, z, positiveEighthOrderTail,
      eighthOrderTail]
    ring_nf
  have base_abs_le_one :
      abs (translatedUnitShiftedKernel center t z) ≤ 1 := by
    rw [translatedUnitShiftedKernel_eq_fermionicKernel,
      abs_of_pos (fermionicKernel_pos t (center + 1 - z))]
    exact fermionicKernel_le_one (rows_nonneg i) (rows_le_one i)
  have tail_abs_le : abs tail ≤
      (2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ) := by
    exact translatedUnitShiftedTail_abs_le
      (rows_nonneg i) (rows_le_one i) (columns_nonneg j) (columns_le_two j)
  rw [decomposition]
  calc
    abs (translatedUnitShiftedKernel center t z + 127 * tail) ≤
        abs (translatedUnitShiftedKernel center t z) + abs (127 * tail) :=
      abs_add_le _ _
    _ ≤ 1 + 127 * ((2 : ℝ) ^ 8 / (Nat.factorial 8 : ℝ)) := by
      rw [abs_mul, abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 127)]
      gcongr
    _ ≤ 2 := by norm_num [Nat.factorial]

/-- Every ordered amplified sample from a translated width-two band is
strictly totally positive. -/
theorem translatedUnitShiftedAmplifiedMatrix_strictlyTotallyPositive
    {order : ℕ} (center : ℝ) {lambda : ℝ} (lambda_one_le : 1 ≤ lambda)
    (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ j, 0 ≤ columns j)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    Matrix.StrictlyTotallyPositive
      (translatedUnitShiftedAmplifiedMatrix center lambda rows columns) := by
  intro m selectedRows selectedColumns selectedRows_mono selectedColumns_mono
  change 0 < (translatedUnitShiftedAmplifiedMatrix center lambda
    (rows ∘ selectedRows) (columns ∘ selectedColumns)).det
  exact translatedUnitShiftedAmplifiedMatrix_det_pos center lambda_one_le
    (rows ∘ selectedRows) (columns ∘ selectedColumns)
    (fun i => rows_nonneg (selectedRows i))
    (fun j => columns_nonneg (selectedColumns j))
    (rows_mono.comp selectedRows_mono) (columns_mono.comp selectedColumns_mono)

/-- Ordered translated-band determinants acquire one factor `128⁻¹` per
sample beyond the first eight, with constants independent of the center. -/
theorem translatedUnitShiftedMatrix_det_le_two_pow_div_128_pow
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ)
    (eight_le : 8 ≤ order)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    (translatedUnitShiftedMatrix center rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let weight : ℝ :=
    (∏ i, Real.exp (-rows i * (center + 1))) *
      ∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹
  have weight_pos : 0 < weight := by
    dsimp only [weight]
    exact mul_pos (Finset.prod_pos fun i _ => Real.exp_pos _)
      (Finset.prod_pos fun j _ => inv_pos.mpr (by positivity))
  have rank_amplification := positiveEighthOrderTail_rank_amplification
    rows columns eight_le rows_nonneg columns_nonneg rows_mono columns_mono
    (lambda := (128 : ℝ)) (by norm_num)
  have shifted_rank_amplification :
      (128 : ℝ) ^ (order - 8) *
          (translatedUnitShiftedMatrix center rows columns).det ≤
        (translatedUnitShiftedAmplifiedMatrix
          center 128 rows columns).det := by
    rw [translatedUnitShiftedMatrix_det,
      translatedUnitShiftedAmplifiedMatrix_det]
    calc
      (128 : ℝ) ^ (order - 8) *
          ((∏ i, Real.exp (-rows i * (center + 1))) *
            ((∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) *
              (expMatrix rows columns).det)) =
        weight * ((128 : ℝ) ^ (order - 8) *
          (expMatrix rows columns).det) := by
            dsimp only [weight]
            ring
      _ ≤ weight *
          ((128 : ℝ) • positiveEighthOrderTailMatrix rows columns +
            positiveEighthOrderTaylorMatrix rows columns).det :=
        mul_le_mul_of_nonneg_left rank_amplification weight_pos.le
      _ = (∏ i, Real.exp (-rows i * (center + 1))) *
          ((∏ j, (1 + Real.exp (columns j - (center + 1)))⁻¹) *
            ((128 : ℝ) • positiveEighthOrderTailMatrix rows columns +
              positiveEighthOrderTaylorMatrix rows columns).det) := by
        dsimp only [weight]
        ring
  have positive :=
    translatedUnitShiftedAmplifiedMatrix_strictlyTotallyPositive
      center (lambda := (128 : ℝ)) (by norm_num) rows columns rows_nonneg
      columns_nonneg rows_mono columns_mono
  have amplified_bound :
      (translatedUnitShiftedAmplifiedMatrix
        center 128 rows columns).det ≤ (2 : ℝ) ^ order := by
    have hadamard := Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive
      (translatedUnitShiftedAmplifiedMatrix center 128 rows columns) positive
    have diagonal_bound :
        (∏ i, translatedUnitShiftedAmplifiedMatrix
          center 128 rows columns i i) ≤ ∏ _i : Fin order, (2 : ℝ) := by
      apply Finset.prod_le_prod
      · intro i _
        exact (positive.entry_pos i i).le
      · intro i _
        exact (le_abs_self _).trans
          (translatedUnitShiftedAmplifiedMatrix_abs_le_two center rows columns
            rows_nonneg rows_le_one columns_nonneg columns_le_two i i)
    exact hadamard.trans (by simpa using diagonal_bound)
  apply (le_div_iff₀ (pow_pos (by norm_num) (order - 8))).2
  simpa only [mul_comm] using shifted_rank_amplification.trans amplified_bound

/-- The translated-band determinant envelope is invariant under independent
permutations of injective row and shifted-frequency samples. -/
theorem translatedUnitShiftedMatrix_abs_det_le_two_pow_div_128_pow
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ)
    (eight_le : 8 ≤ order)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_nonneg : ∀ j, 0 ≤ columns j)
    (columns_le_two : ∀ j, columns j ≤ 2) :
    abs (translatedUnitShiftedMatrix center rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let rowPerm := Tuple.sort rows
  let columnPerm := Tuple.sort columns
  have rows_sorted : StrictMono (rows ∘ rowPerm) :=
    (Tuple.monotone_sort rows).strictMono_of_injective
      (rows_injective.comp rowPerm.injective)
  have columns_sorted : StrictMono (columns ∘ columnPerm) :=
    (Tuple.monotone_sort columns).strictMono_of_injective
      (columns_injective.comp columnPerm.injective)
  let original := translatedUnitShiftedMatrix center rows columns
  have sorted_matrix :
      translatedUnitShiftedMatrix center
        (rows ∘ rowPerm) (columns ∘ columnPerm) =
          original.submatrix rowPerm columnPerm := by
    rfl
  have sorted_pos := translatedUnitShiftedAmplifiedMatrix_det_pos
    center (lambda := (1 : ℝ)) (by norm_num)
    (rows ∘ rowPerm) (columns ∘ columnPerm)
    (fun i => rows_nonneg (rowPerm i))
    (fun j => columns_nonneg (columnPerm j)) rows_sorted columns_sorted
  rw [translatedUnitShiftedAmplifiedMatrix_one] at sorted_pos
  have sorted_bound :=
    translatedUnitShiftedMatrix_det_le_two_pow_div_128_pow center
      (rows ∘ rowPerm) (columns ∘ columnPerm) eight_le
      (fun i => rows_nonneg (rowPerm i)) (fun i => rows_le_one (rowPerm i))
      (fun j => columns_nonneg (columnPerm j))
      (fun j => columns_le_two (columnPerm j)) rows_sorted columns_sorted
  calc
    abs (translatedUnitShiftedMatrix center rows columns).det =
        abs original.det := rfl
    _ = abs (original.submatrix rowPerm columnPerm).det := by
      symm
      exact Matrix.abs_det_submatrix_equiv_equiv rowPerm columnPerm original
    _ = abs (translatedUnitShiftedMatrix center
        (rows ∘ rowPerm) (columns ∘ columnPerm)).det := by rw [sorted_matrix]
    _ = (translatedUnitShiftedMatrix center
        (rows ∘ rowPerm) (columns ∘ columnPerm)).det := abs_of_pos sorted_pos
    _ ≤ (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := sorted_bound

/-- Every injective fermionic sample on a translated width-two frequency band
has the same fixed-base determinant decay as the centered signed unit band. -/
theorem fermionicKernel_translatedUnit_abs_sample_det_le
    {order : ℕ} (center : ℝ) (rows columns : Fin order → ℝ)
    (eight_le : 8 ≤ order)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_le_one : ∀ i, rows i ≤ 1)
    (columns_lower : ∀ j, center - 1 ≤ columns j)
    (columns_upper : ∀ j, columns j ≤ center + 1) :
    abs (sampleMatrix fermionicKernel rows columns).det ≤
      (2 : ℝ) ^ order / (128 : ℝ) ^ (order - 8) := by
  let shiftedColumns : Fin order → ℝ :=
    fun j => center + 1 - columns j
  have shifted_injective : Function.Injective shiftedColumns := by
    intro i j hij
    apply columns_injective
    dsimp only [shiftedColumns] at hij
    linarith
  have shifted_nonneg : ∀ j, 0 ≤ shiftedColumns j := by
    intro j
    dsimp only [shiftedColumns]
    linarith [columns_upper j]
  have shifted_le_two : ∀ j, shiftedColumns j ≤ 2 := by
    intro j
    dsimp only [shiftedColumns]
    linarith [columns_lower j]
  have matrix_eq :
      sampleMatrix fermionicKernel rows columns =
        translatedUnitShiftedMatrix center rows shiftedColumns := by
    ext i j
    exact fermionicKernel_eq_translatedUnitShiftedKernel
      center (rows i) (columns j)
  rw [matrix_eq]
  exact translatedUnitShiftedMatrix_abs_det_le_two_pow_div_128_pow
    center rows shiftedColumns eight_le rows_injective shifted_injective
    rows_nonneg rows_le_one shifted_nonneg shifted_le_two

/-- Complete-pivot residual powers on any translated width-two band inherit
the same center-independent fixed-base determinant decay. -/
theorem fermionicKernel_translatedUnit_gecp_error_pow_le
    (center : ℝ) (residual : ℕ → Kernel ℝ ℝ ℝ)
    (rows columns : ℕ → ℝ)
    (realized : ∀ m, ∃ run : Run fermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => center - 1 ≤ omega ∧ omega ≤ center + 1)
      (residual m) (rows m) (columns m))
    (n : ℕ) (eight_le : 8 ≤ n) (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (run_complete : run.CompleteOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => center - 1 ≤ omega ∧ omega ≤ center + 1))
    (t omega : ℝ) (t_mem : 0 ≤ t ∧ t ≤ 1)
    (omega_mem : center - 1 ≤ omega ∧ omega ≤ center + 1) :
    abs (residual n t omega) ^ n ≤
      (2 : ℝ) ^ n / (128 : ℝ) ^ (n - 8) := by
  have residual_power := strictSignRegular_gecp_error_pow_le_selectedCore_det
    fermionicKernel_strictSignRegular
    (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun omega : ℝ => center - 1 ≤ omega ∧ omega ≤ center + 1)
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
  have finiteColumns_lower : ∀ j, center - 1 ≤ finiteColumns j := by
    intro j
    exact (run.selectedColumn_mem_of_completeOn run_complete _).1
  have finiteColumns_upper : ∀ j, finiteColumns j ≤ center + 1 := by
    intro j
    exact (run.selectedColumn_mem_of_completeOn run_complete _).2
  have eight_le_length : 8 ≤ run.pivots.length := by
    simpa [pivot_length] using eight_le
  have determinant_bound := fermionicKernel_translatedUnit_abs_sample_det_le
    center finiteRows finiteColumns eight_le_length finiteRows_injective
    finiteColumns_injective finiteRows_nonneg finiteRows_le_one
    finiteColumns_lower finiteColumns_upper
  have finite_sample :
      sampleMatrix fermionicKernel finiteRows finiteColumns =
        run.finSelectedCore := by
    rfl
  rw [finite_sample] at determinant_bound
  exact residual_power.trans (by simpa [pivot_length] using determinant_bound)

/-- Twelve exact complete pivots contract every residual entry by one half on
every translated width-two frequency band, with no dependence on its center. -/
theorem fermionicKernel_translatedUnit_gecp_error_le_half_twelve
    (center : ℝ) (residual : ℕ → Kernel ℝ ℝ ℝ)
    (rows columns : ℕ → ℝ)
    (realized : ∀ m, ∃ run : Run fermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => center - 1 ≤ omega ∧ omega ≤ center + 1)
      (residual m) (rows m) (columns m))
    (run : Run fermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns 12)
    (run_complete : run.CompleteOn
      (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => center - 1 ≤ omega ∧ omega ≤ center + 1))
    (t omega : ℝ) (t_mem : 0 ≤ t ∧ t ≤ 1)
    (omega_mem : center - 1 ≤ omega ∧ omega ≤ center + 1) :
    abs (residual 12 t omega) ≤ (1 / 2 : ℝ) := by
  have power_bound := fermionicKernel_translatedUnit_gecp_error_pow_le
    center residual rows columns realized pivot_ne updates complete 12
    (by norm_num) run run_pivots run_complete t omega t_mem omega_mem
  apply le_of_pow_le_pow_left₀ (by norm_num : 12 ≠ 0)
    (by norm_num : (0 : ℝ) ≤ 1 / 2)
  exact power_bound.trans (by norm_num)

end Fermionic
end GECPKernelStructure
