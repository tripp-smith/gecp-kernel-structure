import GECPKernelStructure.Fermionic.ResidualMinors
import GECPKernelStructure.Fermionic.TailAmplification
import Mathlib.Tactic

namespace Matrix

open Finset

variable {n k : Type*} [Fintype n] [DecidableEq n]
  [Fintype k] [DecidableEq k]

/-- Nonnegative mixed coefficients and a strictly rank-deficient low block
make both the original and amplified determinants nonnegative, so the generic
rank-amplification theorem is available in absolute value. -/
theorem mixedColumnPolynomial_rank_amplification_abs
    (E : Matrix n n ℝ) (L : Matrix n k ℝ) (R : Matrix k n ℝ)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_lt : Fintype.card k < Fintype.card n)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card n - Fintype.card k) * abs (E + L * R).det ≤
      abs (lambda • E + L * R).det := by
  have low_det_zero : (L * R).det = 0 :=
    det_mul_rect_eq_zero_of_card_lt L R rank_lt
  have original_nonneg : 0 ≤ (E + L * R).det := by
    have mono := mixedColumnPolynomial_det_mono E (L * R)
      coefficient_nonneg (mu := (0 : ℝ)) (lambda := (1 : ℝ))
      (by norm_num) (by norm_num)
    simpa [low_det_zero] using mono
  have amplified_nonneg : 0 ≤ (lambda • E + L * R).det := by
    have mono := mixedColumnPolynomial_det_mono E (L * R)
      coefficient_nonneg (mu := (0 : ℝ)) (lambda := lambda)
      (by norm_num) (by linarith)
    simpa [low_det_zero] using mono
  rw [abs_of_nonneg original_nonneg, abs_of_nonneg amplified_nonneg]
  exact mixedColumnPolynomial_rank_amplification E L R coefficient_nonneg
    rank_lt.le lambda_one_le

end Matrix

namespace GECPKernelStructure

namespace GECP.Run

open Matrix

universe u v z w

/-- A compatible low-rank-plus-tail split made after exact Schur conditioning
has an amplification exponent depending only on the fresh sample order and
the residual low rank, not on the length of the selected prefix. -/
theorem finalResidual_rank_amplification
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows : γ → α) (columns : γ → β)
    (E : Matrix γ γ ℝ) (L : Matrix γ κ ℝ) (R : Matrix κ γ ℝ)
    (split : (fun i j => run.finalResidual (rows i) (columns j)) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      abs (lambda • E + L * R).det := by
  rw [split]
  exact mixedColumnPolynomial_rank_amplification_abs E L R
    coefficient_nonneg rank_lt lambda_one_le

/-- Phase AJ's arbitrary-border factorization composes with residual tail
amplification without charging the selected prefix against the amplification
exponent.  The prefix appears only through its exact selected-core factor. -/
theorem augmentedCore_rank_amplification
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows : γ → α) (columns : γ → β)
    (E : Matrix γ γ ℝ) (L : Matrix γ κ ℝ) (R : Matrix κ γ ℝ)
    (split : (fun i j => run.finalResidual (rows i) (columns j)) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (run.augmentedCore γ rows columns).det ≤
      abs run.selectedCore.det * abs (lambda • E + L * R).det := by
  have residual_amplification := run.finalResidual_rank_amplification γ κ
    rows columns E L R split coefficient_nonneg rank_lt lambda_one_le
  rw [augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor
      (𝕜 := ℝ) γ run rows columns, abs_mul]
  calc
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        (abs run.selectedCore.det *
          abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j))) =
      abs run.selectedCore.det *
        (lambda ^ (Fintype.card γ - Fintype.card κ) *
          abs (Matrix.det fun i j =>
            run.finalResidual (rows i) (columns j))) := by ring
    _ ≤ abs run.selectedCore.det * abs (lambda • E + L * R).det :=
      mul_le_mul_of_nonneg_left residual_amplification (abs_nonneg _)

end GECP.Run

namespace Fermionic

open GECP Matrix

/-- Exact rational row nodes for the frozen-column amplification obstruction. -/
def frozenColumnRows : Fin 4 → ℚ := ![0, 1 / 4, 1 / 2, 1]

/-- Exact rational column nodes, with the last column kept unmodified. -/
def frozenColumnColumns : Fin 4 → ℚ := ![1 / 4, 1 / 2, 1, 2]

/-- A finite positive exponential Taylor sum over the rationals. -/
def finitePositiveTaylorRat (terms : ℕ) (x : ℚ) : ℚ :=
  ∑ degree ∈ Finset.range terms,
    x ^ degree / (degree.factorial : ℚ)

/-- The high Taylor block on the three fresh columns, zero on the old column. -/
def frozenColumnTail : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  if j = 3 then 0
  else
    finitePositiveTaylorRat 6 (frozenColumnRows i * frozenColumnColumns j) - 1

/-- The rank-two background: three constant columns and one exact old column. -/
def frozenColumnLow : Matrix (Fin 4) (Fin 4) ℚ := fun i j =>
  finitePositiveTaylorRat 6 (frozenColumnRows i * frozenColumnColumns j) -
    frozenColumnTail i j

/-- Explicit two-feature factor for the frozen-column background. -/
def frozenColumnLowLeft : Matrix (Fin 4) (Fin 2) ℚ := fun i k =>
  if k = 0 then 1
  else finitePositiveTaylorRat 6 (frozenColumnRows i * 2)

/-- The matching feature-to-column factor for the rank-two background. -/
def frozenColumnLowRight : Matrix (Fin 2) (Fin 4) ℚ := fun k j =>
  if k = 0 then if j = 3 then 0 else 1
  else if j = 3 then 1 else 0

theorem frozenColumnLow_eq_mul :
    frozenColumnLow = frozenColumnLowLeft * frozenColumnLowRight := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [frozenColumnLow, frozenColumnTail, frozenColumnLowLeft,
      frozenColumnLowRight, finitePositiveTaylorRat, frozenColumnRows,
      frozenColumnColumns, Matrix.mul_apply, Fin.sum_univ_succ]

theorem frozenColumnTaylor_original_det :
    (frozenColumnTail + frozenColumnLow).det =
      (91063777 / 85899345920 : ℚ) := by
  have matrix_eq : frozenColumnTail + frozenColumnLow =
      !![(1 : ℚ), 1, 1, 1;
        133944401 / 125829120, 4455721 / 3932160, 157781 / 122880,
          6331 / 3840;
        4455721 / 3932160, 157781 / 122880, 6331 / 3840, 163 / 60;
        157781 / 122880, 6331 / 3840, 163 / 60, 109 / 15] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [frozenColumnLow, frozenColumnTail, finitePositiveTaylorRat,
        frozenColumnRows, frozenColumnColumns] <;> norm_num
  let lower : Matrix (Fin 4) (Fin 4) ℚ :=
    !![1, 0, 0, 0;
      133944401 / 125829120, 1, 0, 0;
      4455721 / 3932160, 6328224 / 2879557, 1, 0;
      157781 / 122880, 15295488 / 2879557, 1571932976 / 195205593, 1]
  let upper : Matrix (Fin 4) (Fin 4) ℚ :=
    !![1, 1, 1, 1;
      0, 2879557 / 41943040, 9207781 / 41943040, 24503269 / 41943040;
      0, 0, 195205593 / 5897332736, 1767138569 / 5897332736;
      0, 0, 0, 91063777 / 195205593]
  have factorization : lower * upper = frozenColumnTail + frozenColumnLow := by
    rw [matrix_eq]
    ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [lower, upper, Matrix.mul_apply, Fin.sum_univ_succ]
  have lower_triangular : lower.IsLowerTriangular := by
    intro i j hij
    change i < j at hij
    fin_cases i <;> fin_cases j
    all_goals simp at hij
    all_goals norm_num [lower]
  have upper_triangular : upper.IsUpperTriangular := by
    intro i j hij
    change j < i at hij
    fin_cases i <;> fin_cases j
    all_goals simp at hij
    all_goals norm_num [upper]
  rw [← factorization, Matrix.det_mul,
    Matrix.det_of_isLowerTriangular lower lower_triangular,
    Matrix.det_of_isUpperTriangular upper_triangular]
  norm_num [Fin.prod_univ_succ, lower, upper]

theorem frozenColumnTaylor_amplified_det :
    ((128 : ℚ) • frozenColumnTail + frozenColumnLow).det =
      (-32176499 / 983040 : ℚ) := by
  have matrix_eq : (128 : ℚ) • frozenColumnTail + frozenColumnLow =
      !![(1 : ℚ), 1, 1, 1;
        9098321 / 983040, 554281 / 30720, 35861 / 960, 6331 / 3840;
        554281 / 30720, 35861 / 960, 2521 / 30, 163 / 60;
        35861 / 960, 2521 / 30, 3311 / 15, 109 / 15] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [frozenColumnLow, frozenColumnTail, finitePositiveTaylorRat,
        frozenColumnRows, frozenColumnColumns] <;> norm_num
  let lower : Matrix (Fin 4) (Fin 4) ℚ :=
    !![1, 0, 0, 0;
      9098321 / 983040, 1, 0, 0;
      554281 / 30720, 6328224 / 2879557, 1, 0;
      35861 / 960, 15295488 / 2879557, 1571932976 / 195205593, 1]
  let upper : Matrix (Fin 4) (Fin 4) ℚ :=
    !![1, 1, 1, 1;
      0, 2879557 / 327680, 9207781 / 327680, -1495517 / 196608;
      0, 0, 195205593 / 46072912, 192152113 / 138218736;
      0, 0, 0, -514823984 / 585616779]
  have factorization : lower * upper =
      (128 : ℚ) • frozenColumnTail + frozenColumnLow := by
    rw [matrix_eq]
    ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [lower, upper, Matrix.mul_apply, Fin.sum_univ_succ]
  have lower_triangular : lower.IsLowerTriangular := by
    intro i j hij
    change i < j at hij
    fin_cases i <;> fin_cases j
    all_goals simp at hij
    all_goals norm_num [lower]
  have upper_triangular : upper.IsUpperTriangular := by
    intro i j hij
    change j < i at hij
    fin_cases i <;> fin_cases j
    all_goals simp at hij
    all_goals norm_num [upper]
  rw [← factorization, Matrix.det_mul,
    Matrix.det_of_isLowerTriangular lower lower_triangular,
    Matrix.det_of_isUpperTriangular upper_triangular]
  norm_num [Fin.prod_univ_succ, lower, upper]

/-- Freezing only an old column destroys the positive-tail amplification
mechanism: despite a rank-two low block and a positive original determinant,
weight `128` makes the determinant negative. -/
theorem frozenColumnTaylor_rankAmplification_fails :
    ((128 : ℚ) • frozenColumnTail + frozenColumnLow).det <
      128 ^ (4 - 2) * (frozenColumnTail + frozenColumnLow).det := by
  rw [frozenColumnTaylor_original_det, frozenColumnTaylor_amplified_det]
  norm_num

/-- Fermionic final residuals use the generic prefix-independent bridge once
the residual sample itself has a compatible low-rank-plus-tail split. -/
theorem fermionicKernel_finalResidual_rank_amplification
    (run : Run fermionicKernel)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows columns : γ → ℝ)
    (E : Matrix γ γ ℝ) (L : Matrix γ κ ℝ) (R : Matrix κ γ ℝ)
    (split : (fun i j => run.finalResidual (rows i) (columns j)) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      abs (lambda • E + L * R).det :=
  run.finalResidual_rank_amplification γ κ rows columns E L R split
    coefficient_nonneg rank_lt lambda_one_le

/-- The corresponding fermionic augmented determinant pays for the existing
prefix only through its selected-core determinant, not in the tail exponent. -/
theorem fermionicKernel_augmentedCore_rank_amplification
    (run : Run fermionicKernel)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows columns : γ → ℝ)
    (E : Matrix γ γ ℝ) (L : Matrix γ κ ℝ) (R : Matrix κ γ ℝ)
    (split : (fun i j => run.finalResidual (rows i) (columns j)) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree E (L * R))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (run.augmentedCore γ rows columns).det ≤
      abs run.selectedCore.det * abs (lambda • E + L * R).det :=
  run.augmentedCore_rank_amplification γ κ rows columns E L R split
    coefficient_nonneg rank_lt lambda_one_le

end Fermionic
end GECPKernelStructure
