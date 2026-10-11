import GECPKernelStructure.Matrix.AndreiefIntegral
import Mathlib.Data.Fin.Tuple.Sort

namespace Matrix

open Function
open MeasureTheory

variable {n : ℕ} {X : Type*}

/--
Applying one permutation to the sampled columns of the left matrix and the
sampled rows of the right matrix leaves the product of their determinants
unchanged: the two permutation signs cancel.
-/
theorem det_mul_det_comp_perm
    (left : Fin n → X → ℝ) (right : X → Fin n → ℝ)
    (sample : Fin n → X) (σ : Equiv.Perm (Fin n)) :
    Matrix.det (fun i j => left i ((sample ∘ σ) j)) *
        Matrix.det (fun i j => right ((sample ∘ σ) i) j) =
      Matrix.det (fun i j => left i (sample j)) *
        Matrix.det (fun i j => right (sample i) j) := by
  let leftMatrix : Matrix (Fin n) (Fin n) ℝ := fun i j => left i (sample j)
  let rightMatrix : Matrix (Fin n) (Fin n) ℝ := fun i j => right (sample i) j
  change (leftMatrix.submatrix id σ).det *
      (rightMatrix.submatrix σ id).det = leftMatrix.det * rightMatrix.det
  rw [Matrix.det_permute', Matrix.det_permute]
  have sign_square :
      (((σ.sign : ℤ) : ℝ) * ((σ.sign : ℤ) : ℝ)) = 1 := by
    norm_num [← Int.cast_mul, ← Units.val_mul, Int.units_mul_self]
  calc
    ((σ.sign : ℤ) : ℝ) * leftMatrix.det *
          (((σ.sign : ℤ) : ℝ) * rightMatrix.det) =
        (((σ.sign : ℤ) : ℝ) * ((σ.sign : ℤ) : ℝ)) *
          (leftMatrix.det * rightMatrix.det) := by ring
    _ = leftMatrix.det * rightMatrix.det := by rw [sign_square, one_mul]

/--
If a product of sampled determinants has a prescribed weak sign on strictly
ordered samples, then it has that sign on every finite sample.  Injective
samples are sorted by one common permutation; repeated samples make the right
determinant vanish.
-/
theorem signed_det_mul_det_nonneg_of_strictMono
    [LinearOrder X] (signature : ℝ)
    (left : Fin n → X → ℝ) (right : X → Fin n → ℝ)
    (ordered_nonneg : ∀ ordered : Fin n → X, StrictMono ordered →
      0 ≤ signature *
        (Matrix.det (fun i j => left i (ordered j)) *
          Matrix.det (fun i j => right (ordered i) j)))
    (sample : Fin n → X) :
    0 ≤ signature *
      (Matrix.det (fun i j => left i (sample j)) *
        Matrix.det (fun i j => right (sample i) j)) := by
  by_cases sample_injective : Injective sample
  · let σ := Tuple.sort sample
    have sorted_strict : StrictMono (sample ∘ σ) :=
      (Tuple.monotone_sort sample).strictMono_of_injective
        (sample_injective.comp σ.injective)
    rw [← det_mul_det_comp_perm left right sample σ]
    exact ordered_nonneg (sample ∘ σ) sorted_strict
  · rw [not_injective_iff] at sample_injective
    obtain ⟨i, j, sample_eq, hij⟩ := sample_injective
    have right_det_zero :
        Matrix.det (fun row column => right (sample row) column) = 0 := by
      apply Matrix.det_zero_of_row_eq hij
      funext column
      rw [sample_eq]
    rw [right_det_zero, mul_zero, mul_zero]

/-- Pointwise ordered-sample sign control transfers to a product-measure integral. -/
theorem signature_mul_integral_det_mul_det_nonneg
    [LinearOrder X] [MeasurableSpace X] (μ : Measure X)
    (signature : ℝ) (left : Fin n → X → ℝ) (right : X → Fin n → ℝ)
    (ordered_nonneg : ∀ ordered : Fin n → X, StrictMono ordered →
      0 ≤ signature *
        (Matrix.det (fun i j => left i (ordered j)) *
          Matrix.det (fun i j => right (ordered i) j))) :
    0 ≤ signature *
      ∫ sample,
        Matrix.det (fun i j => left i (sample j)) *
          Matrix.det (fun i j => right (sample i) j)
        ∂Measure.pi (fun _ : Fin n => μ) := by
  rw [← MeasureTheory.integral_const_mul]
  exact MeasureTheory.integral_nonneg_of_ae
    (Filter.Eventually.of_forall fun sample =>
      signed_det_mul_det_nonneg_of_strictMono signature left right
        ordered_nonneg sample)

/--
Weak-sign closure for a finite integral-kernel composition.  This is the
nonnegative consequence of Andréief after the common sample is sorted.
-/
theorem signature_mul_det_integral_nonneg
    [LinearOrder X] [MeasurableSpace X] (μ : Measure X) [SigmaFinite μ]
    (signature : ℝ) (left : Fin n → X → ℝ) (right : X → Fin n → ℝ)
    (pair_integrable : ∀ i j, Integrable (fun x => left i x * right x j) μ)
    (ordered_nonneg : ∀ ordered : Fin n → X, StrictMono ordered →
      0 ≤ signature *
        (Matrix.det (fun i j => left i (ordered j)) *
          Matrix.det (fun i j => right (ordered i) j))) :
    0 ≤ signature *
      Matrix.det (fun i j => ∫ x, left i x * right x j ∂μ) := by
  have signed_integral_nonneg :=
    signature_mul_integral_det_mul_det_nonneg μ signature left right ordered_nonneg
  rw [← factorial_mul_det_integral_eq_integral_det_mul_det μ left right
    pair_integrable] at signed_integral_nonneg
  have factorial_pos : 0 < ((Nat.factorial n : ℕ) : ℝ) := by positivity
  apply nonneg_of_mul_nonneg_left _ factorial_pos
  simpa [Fintype.card_fin, mul_assoc, mul_left_comm, mul_comm] using
    signed_integral_nonneg

end Matrix
