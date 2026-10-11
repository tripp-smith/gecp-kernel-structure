import GECPKernelStructure.Matrix.BlockFeatureCauchyBinet
import GECPKernelStructure.Matrix.LowRankDeterminant
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.LinearAlgebra.Matrix.Polynomial

namespace Matrix

open Finset Polynomial

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Sum of mixed-column determinants using exactly `degree` columns from `E`. -/
def mixedColumnCoefficient (degree : ℕ) (E B : Matrix n n ℝ) : ℝ :=
  ∑ s : Finset n, if s.card = degree then (columnChoice s E B).det else 0

/-- Determinant multilinearity packaged as a polynomial by error-column count. -/
noncomputable def mixedColumnPolynomial (E B : Matrix n n ℝ) : ℝ[X] :=
  ∑ s : Finset n, C (columnChoice s E B).det * X ^ s.card

theorem mixedColumnPolynomial_coeff (degree : ℕ) (E B : Matrix n n ℝ) :
    (mixedColumnPolynomial E B).coeff degree =
      mixedColumnCoefficient degree E B := by
  simp [mixedColumnPolynomial, mixedColumnCoefficient,
    Polynomial.coeff_C_mul, Polynomial.coeff_X_pow, eq_comm]

theorem det_columnChoice_smul (s : Finset n) (lambda : ℝ)
    (E B : Matrix n n ℝ) :
    (columnChoice s (lambda • E) B).det =
      lambda ^ s.card * (columnChoice s E B).det := by
  have matrix_eq :
      columnChoice s (lambda • E) B =
        Matrix.of fun i j => (if j ∈ s then lambda else 1) *
          columnChoice s E B i j := by
    ext i j
    by_cases hj : j ∈ s <;> simp [columnChoice, hj]
  rw [matrix_eq, Matrix.det_mul_row]
  congr 1
  simp

theorem mixedColumnPolynomial_eval (lambda : ℝ) (E B : Matrix n n ℝ) :
    (mixedColumnPolynomial E B).eval lambda = (lambda • E + B).det := by
  rw [mixedColumnPolynomial, Polynomial.eval_finsetSum]
  simp_rw [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow,
    Polynomial.eval_X]
  rw [Matrix.det_add_eq_sum_columnChoices]
  apply Finset.sum_congr rfl
  intro s _
  rw [det_columnChoice_smul]
  ring

variable {a b : ℕ}

/-- The unmarked block of a rectangular feature matrix. -/
noncomputable def lowFeatureLeft (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) : Matrix (Fin a) (Fin b) ℝ :=
  fun i k => if high k then 0 else A i k

/-- The marked block of a rectangular feature matrix. -/
noncomputable def highFeatureLeft (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) : Matrix (Fin a) (Fin b) ℝ :=
  fun i k => if high k then A i k else 0

theorem weightedFeatureLeft_eq_low_add_smul_high
    (high : Fin b → Prop) [DecidablePred high] (lambda : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ) :
    weightedFeatureLeft high lambda A =
      lowFeatureLeft high A + lambda • highFeatureLeft high A := by
  ext i k
  by_cases hk : high k <;>
    simp [weightedFeatureLeft, featureBlockWeight, lowFeatureLeft,
      highFeatureLeft, hk]

theorem weightedFeatureProduct_eq_low_add_smul_high
    (high : Fin b → Prop) [DecidablePred high] (lambda : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ) :
    weightedFeatureLeft high lambda A * B =
      lowFeatureLeft high A * B +
        lambda • (highFeatureLeft high A * B) := by
  rw [weightedFeatureLeft_eq_low_add_smul_high, Matrix.add_mul,
    Matrix.smul_mul]

/-- The feature-count determinant expansion as a polynomial. -/
noncomputable def blockFeaturePolynomial
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ) : ℝ[X] :=
  ∑ q : {q : Fin a → Fin b // StrictMono q},
    C (signedOrderedCauchyBinetTerm signature A B q) *
      X ^ highFeatureCount high q.1

theorem blockFeaturePolynomial_coeff
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ)
    (degree : ℕ) :
    (blockFeaturePolynomial signature high A B).coeff degree =
      blockFeatureCoefficient signature high A B degree := by
  simp [blockFeaturePolynomial, blockFeatureCoefficient,
    Polynomial.coeff_C_mul, Polynomial.coeff_X_pow, eq_comm]

theorem blockFeaturePolynomial_eval
    (high : Fin b → Prop) [DecidablePred high] (lambda : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ) :
    (blockFeaturePolynomial 1 high A B).eval lambda =
      (weightedFeatureLeft high lambda A * B).det := by
  rw [blockFeaturePolynomial, Polynomial.eval_finsetSum]
  simp_rw [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow,
    Polynomial.eval_X]
  simpa [mul_comm] using
    (signature_mul_det_weightedFeatureProduct_eq_sum
      1 high lambda A B).symm

theorem blockFeaturePolynomial_eq_mixedColumnPolynomial
    (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ) :
    blockFeaturePolynomial 1 high A B =
      mixedColumnPolynomial (highFeatureLeft high A * B)
        (lowFeatureLeft high A * B) := by
  apply Polynomial.funext
  intro lambda
  rw [blockFeaturePolynomial_eval, mixedColumnPolynomial_eval,
    weightedFeatureProduct_eq_low_add_smul_high]
  simp [add_comm]

theorem blockFeatureCoefficient_eq_mixedColumnCoefficient
    (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ)
    (degree : ℕ) :
    blockFeatureCoefficient 1 high A B degree =
      mixedColumnCoefficient degree (highFeatureLeft high A * B)
        (lowFeatureLeft high A * B) := by
  rw [← blockFeaturePolynomial_coeff, ← mixedColumnPolynomial_coeff,
    blockFeaturePolynomial_eq_mixedColumnPolynomial]

theorem mixedColumnCoefficient_nonneg_of_orderedTerms
    (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ)
    (base_nonneg : ∀ q : {q : Fin a → Fin b // StrictMono q},
      0 ≤ signedOrderedCauchyBinetTerm 1 A B q)
    (degree : ℕ) :
    0 ≤ mixedColumnCoefficient degree (highFeatureLeft high A * B)
      (lowFeatureLeft high A * B) := by
  rw [← blockFeatureCoefficient_eq_mixedColumnCoefficient]
  exact blockFeatureCoefficient_nonneg 1 high A B base_nonneg degree

end Matrix

namespace GECPKernelStructure
namespace Fermionic

open Matrix

theorem expTaylor_mixedColumnCoefficient_nonneg
    {order terms : ℕ} (lowTerms : ℕ) (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (degree : ℕ) :
    0 ≤ mixedColumnCoefficient degree
      (highFeatureLeft (fun k : Fin terms => lowTerms ≤ k.1)
        (expTaylorLeft (terms := terms) rows) *
          expTaylorRight (terms := terms) columns)
      (lowFeatureLeft (fun k : Fin terms => lowTerms ≤ k.1)
        (expTaylorLeft (terms := terms) rows) *
          expTaylorRight (terms := terms) columns) := by
  apply mixedColumnCoefficient_nonneg_of_orderedTerms
  intro q
  exact expTaylor_signedOrderedCauchyBinetTerm_nonneg rows columns
    rows_nonneg columns_nonneg rows_mono columns_mono q

end Fermionic
end GECPKernelStructure
