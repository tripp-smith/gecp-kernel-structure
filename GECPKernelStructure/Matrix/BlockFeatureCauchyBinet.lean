import GECPKernelStructure.Fermionic.ExponentialTotalPositivity

namespace Matrix

open Finset

/-- Scalar weight that marks one finite feature block by `lambda`. -/
def featureBlockWeight {b : ℕ} (high : Fin b → Prop) [DecidablePred high]
    (lambda : ℝ) (k : Fin b) : ℝ :=
  if high k then lambda else 1

/-- Scale exactly the marked columns of a rectangular left feature matrix. -/
def weightedFeatureLeft {a b : ℕ} (high : Fin b → Prop) [DecidablePred high]
    (lambda : ℝ) (A : Matrix (Fin a) (Fin b) ℝ) :
    Matrix (Fin a) (Fin b) ℝ :=
  fun i k => featureBlockWeight high lambda k * A i k

/-- Number of marked features used by one ordered Cauchy--Binet selection. -/
def highFeatureCount {a b : ℕ} (high : Fin b → Prop) [DecidablePred high]
    (q : Fin a → Fin b) : ℕ :=
  (Finset.univ.filter fun i => high (q i)).card

theorem prod_featureBlockWeight_comp {a b : ℕ}
    (high : Fin b → Prop) [DecidablePred high] (lambda : ℝ)
    (q : Fin a → Fin b) :
    (∏ i, featureBlockWeight high lambda (q i)) =
      lambda ^ highFeatureCount high q := by
  simp [featureBlockWeight, highFeatureCount, Finset.prod_ite]

/-- Scaling a feature block scales each selected left minor by its marked count. -/
theorem det_weightedFeatureLeft_submatrix {a b : ℕ}
    (high : Fin b → Prop) [DecidablePred high] (lambda : ℝ)
    (A : Matrix (Fin a) (Fin b) ℝ) (q : Fin a → Fin b) :
    ((weightedFeatureLeft high lambda A).submatrix id q).det =
      lambda ^ highFeatureCount high q * (A.submatrix id q).det := by
  have matrix_eq :
      (weightedFeatureLeft high lambda A).submatrix id q =
        Matrix.of fun i j => featureBlockWeight high lambda (q j) *
          (A.submatrix id q) i j := by
    ext i j
    rfl
  rw [matrix_eq, Matrix.det_mul_row, prod_featureBlockWeight_comp]

/--
Exact block-feature Cauchy--Binet expansion. Each ordered selection receives
one power of `lambda` for every marked feature it uses.
-/
theorem signature_mul_det_weightedFeatureProduct_eq_sum {a b : ℕ}
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (lambda : ℝ) (A : Matrix (Fin a) (Fin b) ℝ)
    (B : Matrix (Fin b) (Fin a) ℝ) :
    signature * (weightedFeatureLeft high lambda A * B).det =
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        lambda ^ highFeatureCount high q.1 *
          signedOrderedCauchyBinetTerm signature A B q := by
  rw [signature_mul_det_mul_rect_eq_sum]
  apply Finset.sum_congr rfl
  intro q _
  rw [signedOrderedCauchyBinetTerm, det_weightedFeatureLeft_submatrix]
  rw [signedOrderedCauchyBinetTerm]
  ring

/-- Fixed-degree coefficient in the block-feature determinant expansion. -/
def blockFeatureCoefficient {a b : ℕ}
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ)
    (degree : ℕ) : ℝ :=
  ∑ q : {q : Fin a → Fin b // StrictMono q},
    if highFeatureCount high q.1 = degree then
      signedOrderedCauchyBinetTerm signature A B q else 0

/-- Every fixed block-degree coefficient is nonnegative when every base term is. -/
theorem blockFeatureCoefficient_nonneg {a b : ℕ}
    (signature : ℝ) (high : Fin b → Prop) [DecidablePred high]
    (A : Matrix (Fin a) (Fin b) ℝ) (B : Matrix (Fin b) (Fin a) ℝ)
    (base_nonneg : ∀ q : {q : Fin a → Fin b // StrictMono q},
      0 ≤ signedOrderedCauchyBinetTerm signature A B q)
    (degree : ℕ) :
    0 ≤ blockFeatureCoefficient signature high A B degree := by
  apply Finset.sum_nonneg
  intro q _
  by_cases hdegree : highFeatureCount high q.1 = degree
  · rw [if_pos hdegree]
    exact base_nonneg q
  · rw [if_neg hdegree]

end Matrix

namespace GECPKernelStructure
namespace Fermionic

open Matrix

/-- Positive Taylor feature selections have nonnegative Cauchy--Binet terms. -/
theorem expTaylor_signedOrderedCauchyBinetTerm_nonneg
    {order terms : ℕ} (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (q : {q : Fin order → Fin terms // StrictMono q}) :
    0 ≤ signedOrderedCauchyBinetTerm 1
      (expTaylorLeft (terms := terms) rows)
      (expTaylorRight (terms := terms) columns) q := by
  let exponents : Fin order → ℕ := fun i => (q.1 i).1
  have exponents_mono : StrictMono exponents :=
    Fin.val_strictMono.comp q.2
  have left_factor :
      (expTaylorLeft (terms := terms) rows).submatrix id q.1 =
        generalizedVandermonde rows exponents *
          Matrix.diagonal (fun j : Fin order =>
            ((exponents j).factorial : ℝ)⁻¹) := by
    ext i j
    simp [expTaylorLeft, generalizedVandermonde, exponents]
  have left_nonneg :
      0 ≤ ((expTaylorLeft (terms := terms) rows).submatrix id q.1).det := by
    rw [left_factor, Matrix.det_mul, Matrix.det_diagonal]
    exact mul_nonneg
      (generalizedVandermonde_det_nonneg rows exponents rows_nonneg
        rows_mono exponents_mono)
      (Finset.prod_nonneg fun j _ => inv_nonneg.mpr (Nat.cast_nonneg _))
  have right_transpose :
      ((expTaylorRight (terms := terms) columns).submatrix q.1 id)ᵀ =
        generalizedVandermonde columns exponents := by
    ext i j
    rfl
  have right_nonneg :
      0 ≤ ((expTaylorRight (terms := terms) columns).submatrix q.1 id).det := by
    rw [← Matrix.det_transpose, right_transpose]
    exact generalizedVandermonde_det_nonneg columns exponents columns_nonneg
      columns_mono exponents_mono
  simpa [signedOrderedCauchyBinetTerm] using mul_nonneg left_nonneg right_nonneg

/-- Positive finite Taylor block coefficients are nonnegative at every degree. -/
theorem expTaylor_blockFeatureCoefficient_nonneg
    {order terms : ℕ} (lowTerms : ℕ) (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (degree : ℕ) :
    0 ≤ blockFeatureCoefficient 1 (fun k : Fin terms => lowTerms ≤ k.1)
      (expTaylorLeft (terms := terms) rows)
      (expTaylorRight (terms := terms) columns) degree := by
  apply blockFeatureCoefficient_nonneg
  intro q
  exact expTaylor_signedOrderedCauchyBinetTerm_nonneg rows columns
    rows_nonneg columns_nonneg rows_mono columns_mono q

/-- Exact finite positive-Taylor determinant expansion by high-feature count. -/
theorem expTaylor_blockWeighted_det_eq_sum
    {order terms : ℕ} (lowTerms : ℕ) (lambda : ℝ)
    (rows columns : Fin order → ℝ) :
    (weightedFeatureLeft (fun k : Fin terms => lowTerms ≤ k.1) lambda
        (expTaylorLeft (terms := terms) rows) *
      expTaylorRight (terms := terms) columns).det =
      ∑ q : {q : Fin order → Fin terms // StrictMono q},
        lambda ^ highFeatureCount (fun k : Fin terms => lowTerms ≤ k.1) q.1 *
          signedOrderedCauchyBinetTerm 1
            (expTaylorLeft (terms := terms) rows)
            (expTaylorRight (terms := terms) columns) q := by
  simpa using signature_mul_det_weightedFeatureProduct_eq_sum
    1 (fun k : Fin terms => lowTerms ≤ k.1) lambda
    (expTaylorLeft (terms := terms) rows)
    (expTaylorRight (terms := terms) columns)

end Fermionic
end GECPKernelStructure
