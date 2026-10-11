import GECPKernelStructure.Fermionic.TruncatedPowerTotalNonnegative

namespace Matrix

open scoped BigOperators
open Finset

universe u v

variable {ι : Type u} [Fintype ι] [DecidableEq ι]
variable {X : Type v} {R : Type*} [CommRing R]

/-- The row-wise product integrand occurring before Andréief symmetrization. -/
def compositionSampleIntegrand
    (left : ι → X → R) (right : X → ι → R) (sample : ι → X) : R :=
  (∏ i, left i (sample i)) * Matrix.det (fun i j => right (sample i) j)

/--
Summing the row-wise composition integrand over all coordinate permutations
equals the product of the two sampled determinants.
-/
theorem sum_compositionSampleIntegrand_comp_perm
    (left : ι → X → R) (right : X → ι → R) (sample : ι → X) :
    ∑ σ : Equiv.Perm ι, compositionSampleIntegrand left right (sample ∘ σ) =
      Matrix.det (fun i j => left i (sample j)) *
        Matrix.det (fun i j => right (sample i) j) := by
  let leftMatrix : Matrix ι ι R := fun i j => left i (sample j)
  let rightMatrix : Matrix ι ι R := fun i j => right (sample i) j
  have right_permute (σ : Equiv.Perm ι) :
      Matrix.det (fun i j => right ((sample ∘ σ) i) j) =
        (((σ.sign : ℤ) : R) * rightMatrix.det) := by
    change (rightMatrix.submatrix σ id).det = _
    exact Matrix.det_permute σ rightMatrix
  change (∑ σ : Equiv.Perm ι,
      compositionSampleIntegrand left right (sample ∘ σ)) =
    leftMatrix.det * rightMatrix.det
  rw [← Matrix.det_transpose leftMatrix, Matrix.det_apply']
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro σ _
  unfold compositionSampleIntegrand
  rw [show Matrix.det (fun i j => right ((sample ∘ σ) i) j) =
    (((σ.sign : ℤ) : R) * rightMatrix.det) from right_permute σ]
  simp only [Function.comp_apply, transpose_apply]
  simp [leftMatrix, mul_comm, mul_left_comm]

end Matrix
