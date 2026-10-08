import GECPKernelStructure.Fermionic.SignRegularity
import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- An increasing real tuple has a strictly positive Vandermonde determinant. -/
theorem vandermonde_det_pos_of_strictMono {n : ℕ} (values : Fin n → ℝ)
    (values_mono : StrictMono values) :
    0 < (Matrix.vandermonde values).det := by
  rw [Matrix.det_vandermonde]
  refine Finset.prod_pos fun i _ => ?_
  refine Finset.prod_pos fun j hj => ?_
  exact sub_pos.mpr (values_mono (Finset.mem_Ioi.mp hj))

/--
The square principal block of the exponential series, using powers
`0, ..., n - 1` and the positive weights `1 / k!`.
-/
noncomputable def expTaylorPrincipalMatrix {n : ℕ}
    (rows columns : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.vandermonde rows *
    (Matrix.diagonal (fun k : Fin n => ((k.1.factorial : ℝ)⁻¹)) *
      (Matrix.vandermonde columns)ᵀ)

/-- The principal block evaluates the first `n` terms of `exp(x*y)`. -/
theorem expTaylorPrincipalMatrix_apply {n : ℕ}
    (rows columns : Fin n → ℝ) (i j : Fin n) :
    expTaylorPrincipalMatrix rows columns i j =
      ∑ k : Fin n, rows i ^ k.1 * (k.1.factorial : ℝ)⁻¹ * columns j ^ k.1 := by
  rw [expTaylorPrincipalMatrix, Matrix.mul_apply]
  simp only [Matrix.diagonal_mul, Matrix.transpose_apply, Matrix.vandermonde_apply]
  congr 1
  funext k
  ring

/--
The principal exponential-series block has positive determinant on increasing
row and column tuples. This certifies one strictly positive Cauchy--Binet term;
it does not yet control the signs of the remaining generalized terms.
-/
theorem expTaylorPrincipalMatrix_det_pos {n : ℕ}
    (rows columns : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 < (expTaylorPrincipalMatrix rows columns).det := by
  have rows_pos := vandermonde_det_pos_of_strictMono rows rows_mono
  have columns_pos := vandermonde_det_pos_of_strictMono columns columns_mono
  have weights_pos :
      0 < (Matrix.diagonal
        (fun k : Fin n => ((k.1.factorial : ℝ)⁻¹))).det := by
    rw [Matrix.det_diagonal]
    exact Finset.prod_pos fun k _ =>
      inv_pos.mpr (Nat.cast_pos.mpr k.1.factorial_pos)
  simp only [expTaylorPrincipalMatrix, Matrix.det_mul, Matrix.det_transpose]
  positivity

end Fermionic
end GECPKernelStructure
