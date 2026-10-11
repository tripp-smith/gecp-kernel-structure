import GECPKernelStructure.Fermionic.PositiveSeriesMixedCoefficients
import GECPKernelStructure.Matrix.LowRankDeterminant

namespace Matrix

open Finset Polynomial

variable {n k : Type*} [Fintype n] [DecidableEq n]
  [Fintype k] [DecidableEq k]

theorem mixedColumnCoefficient_eq_zero_of_rank
    (degree : ℕ) (E : Matrix n n ℝ) (L : Matrix n k ℝ)
    (R : Matrix k n ℝ)
    (hcard : Fintype.card k + degree < Fintype.card n) :
    mixedColumnCoefficient degree E (L * R) = 0 := by
  unfold mixedColumnCoefficient
  apply Finset.sum_eq_zero
  intro s _
  by_cases hs : s.card = degree
  · rw [if_pos hs]
    apply det_columnChoice_eq_zero_of_card_add_lt
    simpa [hs] using hcard
  · rw [if_neg hs]

/-- A determinant pencil with nonnegative mixed-column coefficients is monotone
on nonnegative scalar weights. -/
theorem mixedColumnPolynomial_det_mono
    (E B : Matrix n n ℝ)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E B)
    {mu lambda : ℝ} (mu_nonneg : 0 ≤ mu) (mu_le_lambda : mu ≤ lambda) :
    (mu • E + B).det ≤ (lambda • E + B).det := by
  let p := mixedColumnPolynomial E B
  have p_coefficient_nonneg (degree : ℕ) : 0 ≤ p.coeff degree := by
    rw [mixedColumnPolynomial_coeff]
    exact coefficient_nonneg degree
  rw [← mixedColumnPolynomial_eval mu E B,
    ← mixedColumnPolynomial_eval lambda E B]
  rw [Polynomial.eval_eq_sum, Polynomial.eval_eq_sum]
  apply Finset.sum_le_sum
  intro degree _
  exact mul_le_mul_of_nonneg_left
    (pow_le_pow_left₀ mu_nonneg mu_le_lambda degree)
    (p_coefficient_nonneg degree)

/-- Low-rank vanishing and nonnegative coefficients amplify every tail weight above one. -/
theorem mixedColumnPolynomial_rank_amplification
    (E : Matrix n n ℝ) (L : Matrix n k ℝ) (R : Matrix k n ℝ)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_le : Fintype.card k ≤ Fintype.card n)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card n - Fintype.card k) * (E + L * R).det ≤
      (lambda • E + L * R).det := by
  let p := mixedColumnPolynomial E (L * R)
  have p_coefficient_nonneg (degree : ℕ) : 0 ≤ p.coeff degree := by
    rw [mixedColumnPolynomial_coeff]
    exact coefficient_nonneg degree
  have p_coefficient_zero (degree : ℕ)
      (hdegree : Fintype.card k + degree < Fintype.card n) :
      p.coeff degree = 0 := by
    rw [mixedColumnPolynomial_coeff]
    exact mixedColumnCoefficient_eq_zero_of_rank degree E L R hdegree
  have eval_one := mixedColumnPolynomial_eval 1 E (L * R)
  simp only [one_smul] at eval_one
  rw [← eval_one, ← mixedColumnPolynomial_eval lambda E (L * R)]
  rw [Polynomial.eval_eq_sum, Polynomial.eval_eq_sum]
  simp only [one_pow, mul_one]
  change lambda ^ (Fintype.card n - Fintype.card k) *
      (∑ degree ∈ p.support, p.coeff degree) ≤
    ∑ degree ∈ p.support, p.coeff degree * lambda ^ degree
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro degree _
  by_cases hdegree : Fintype.card k + degree < Fintype.card n
  · rw [p_coefficient_zero degree hdegree]
    simp
  · have exponent_le : Fintype.card n - Fintype.card k ≤ degree := by
      omega
    calc
      lambda ^ (Fintype.card n - Fintype.card k) * p.coeff degree =
          p.coeff degree * lambda ^ (Fintype.card n - Fintype.card k) := by ring
      _ ≤ p.coeff degree * lambda ^ degree :=
        mul_le_mul_of_nonneg_left
          (pow_le_pow_right₀ lambda_one_le exponent_le)
          (p_coefficient_nonneg degree)

end Matrix

namespace GECPKernelStructure
namespace Fermionic

open Matrix

theorem positiveEighthOrderTail_rank_amplification
    {order : ℕ} (rows columns : Fin order → ℝ)
    (eight_le : 8 ≤ order)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (order - 8) * (expMatrix rows columns).det ≤
      (lambda • positiveEighthOrderTailMatrix rows columns +
        positiveEighthOrderTaylorMatrix rows columns).det := by
  have amplification := mixedColumnPolynomial_rank_amplification
    (n := Fin order) (k := Fin 8)
    (positiveEighthOrderTailMatrix rows columns)
    (expTaylorLeft (terms := 8) rows)
    (expTaylorRight (terms := 8) columns)
    (positiveEighthOrderTail_mixedColumnCoefficient_nonneg
      rows columns rows_nonneg columns_nonneg rows_mono columns_mono)
    (by simpa using eight_le)
    lambda_one_le
  rw [show expTaylorLeft (terms := 8) rows *
      expTaylorRight (terms := 8) columns =
        positiveEighthOrderTaylorMatrix rows columns by rfl] at amplification
  have split :
      positiveEighthOrderTailMatrix rows columns +
          positiveEighthOrderTaylorMatrix rows columns =
        expMatrix rows columns := by
    simp [positiveEighthOrderTailMatrix]
  simpa [split] using amplification

/-- Reflected positive-series fermionic sample matrix. -/
noncomputable def positiveSeriesFermionicMatrix {order : ℕ}
    (rows columns : Fin order → ℝ) : Matrix (Fin order) (Fin order) ℝ :=
  fun i j => expMatrix rows columns i j * (1 + Real.exp (columns j))⁻¹

/-- Tail-amplified reflected fermionic sample matrix. -/
noncomputable def positiveSeriesAmplifiedFermionicMatrix {order : ℕ}
    (lambda : ℝ) (rows columns : Fin order → ℝ) :
    Matrix (Fin order) (Fin order) ℝ :=
  fun i j =>
    (lambda • positiveEighthOrderTailMatrix rows columns +
      positiveEighthOrderTaylorMatrix rows columns) i j *
        (1 + Real.exp (columns j))⁻¹

theorem positiveSeriesFermionic_rank_amplification
    {order : ℕ} (rows columns : Fin order → ℝ)
    (eight_le : 8 ≤ order)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (order - 8) * (positiveSeriesFermionicMatrix rows columns).det ≤
      (positiveSeriesAmplifiedFermionicMatrix lambda rows columns).det := by
  let weights : Fin order → ℝ :=
    fun j => (1 + Real.exp (columns j))⁻¹
  have weights_pos : 0 < ∏ j, weights j := by
    apply Finset.prod_pos
    intro j _
    exact inv_pos.mpr (by positivity)
  have original_det :
      (positiveSeriesFermionicMatrix rows columns).det =
        (∏ j, weights j) * (expMatrix rows columns).det := by
    rw [show positiveSeriesFermionicMatrix rows columns =
        fun i j => weights j * expMatrix rows columns i j by
          ext i j
          simp only [positiveSeriesFermionicMatrix, weights]
          ring]
    exact Matrix.det_mul_row weights (expMatrix rows columns)
  let amplified :=
    lambda • positiveEighthOrderTailMatrix rows columns +
      positiveEighthOrderTaylorMatrix rows columns
  have amplified_det :
      (positiveSeriesAmplifiedFermionicMatrix lambda rows columns).det =
        (∏ j, weights j) * amplified.det := by
    rw [show positiveSeriesAmplifiedFermionicMatrix lambda rows columns =
        fun i j => weights j * amplified i j by
          ext i j
          simp only [positiveSeriesAmplifiedFermionicMatrix, weights, amplified]
          ring]
    exact Matrix.det_mul_row weights amplified
  rw [original_det, amplified_det]
  have base := positiveEighthOrderTail_rank_amplification rows columns
    eight_le rows_nonneg columns_nonneg rows_mono columns_mono lambda_one_le
  nlinarith [mul_le_mul_of_nonneg_left base weights_pos.le]

end Fermionic
end GECPKernelStructure
