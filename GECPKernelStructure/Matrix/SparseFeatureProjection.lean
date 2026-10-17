import GECPKernelStructure.Matrix.CauchyBinet

namespace Matrix

open Finset

/-- The one-hot inclusion of an ordered `a`-tuple into a `b`-dimensional
feature space. -/
def orderedFeatureInclusion {a b : ℕ} (q : Fin a → Fin b) :
    Matrix (Fin b) (Fin a) ℝ :=
  fun k i => if k = q i then 1 else 0

@[simp]
theorem mul_orderedFeatureInclusion {a b : ℕ}
    (L : Matrix (Fin a) (Fin b) ℝ) (q : Fin a → Fin b) :
    L * orderedFeatureInclusion q = L.submatrix id q := by
  ext i j
  simp [Matrix.mul_apply, orderedFeatureInclusion]

@[simp]
theorem orderedFeatureInclusion_transpose_mul {a b : ℕ}
    (q : Fin a → Fin b) (R : Matrix (Fin b) (Fin a) ℝ) :
    (orderedFeatureInclusion q)ᵀ * R = R.submatrix q id := by
  ext i j
  simp [Matrix.mul_apply, orderedFeatureInclusion]

/-- The feature-space projection supported on one ordered coordinate tuple. -/
def sparseFeatureProjection {a b : ℕ} (q : Fin a → Fin b) :
    Matrix (Fin b) (Fin b) ℝ :=
  orderedFeatureInclusion q * (orderedFeatureInclusion q)ᵀ

/-- Every entry of a sparse feature projection is nonnegative. -/
theorem sparseFeatureProjection_nonneg {a b : ℕ} (q : Fin a → Fin b)
    (i j : Fin b) :
    0 ≤ sparseFeatureProjection q i j := by
  classical
  rw [sparseFeatureProjection, Matrix.mul_apply]
  apply Finset.sum_nonneg
  intro k _
  simp only [Matrix.transpose_apply, orderedFeatureInclusion]
  split_ifs <;> norm_num

/-- For an injective feature tuple, the sparse feature projection has zero
off-diagonal entries; together with `sparseFeatureProjection_nonneg`, this
makes it a nonnegative diagonal action. -/
theorem sparseFeatureProjection_offDiagonal_eq_zero {a b : ℕ}
    (q : Fin a → Fin b)
    {i j : Fin b} (hij : i ≠ j) :
    sparseFeatureProjection q i j = 0 := by
  classical
  rw [sparseFeatureProjection, Matrix.mul_apply]
  apply Finset.sum_eq_zero
  intro k _
  by_cases hik : i = q k
  · have hjk : j ≠ q k := by
      intro hjk
      exact hij (hik.trans hjk.symm)
    simp [orderedFeatureInclusion, hik, hjk]
  · simp [orderedFeatureInclusion, hik]

/-- Inserting a sparse feature projection isolates exactly one matched pair
of rectangular minors. -/
theorem det_mul_sparseFeatureProjection_mul {a b : ℕ}
    (L : Matrix (Fin a) (Fin b) ℝ) (R : Matrix (Fin b) (Fin a) ℝ)
    (q : Fin a → Fin b) :
    (L * sparseFeatureProjection q * R).det =
      (L.submatrix id q).det * (R.submatrix q id).det := by
  have factorization :
      L * sparseFeatureProjection q * R =
        (L * orderedFeatureInclusion q) *
          ((orderedFeatureInclusion q)ᵀ * R) := by
    simp only [sparseFeatureProjection, Matrix.mul_assoc]
  rw [factorization, mul_orderedFeatureInclusion,
    orderedFeatureInclusion_transpose_mul, Matrix.det_mul]

/-- A negative matched Plucker coordinate is therefore a counterexample to
sign preservation by arbitrary nonnegative diagonal feature actions. -/
theorem mul_sparseFeatureProjection_mul_det_neg {a b : ℕ}
    (signature : ℝ)
    (L : Matrix (Fin a) (Fin b) ℝ) (R : Matrix (Fin b) (Fin a) ℝ)
    (q : Fin a → Fin b)
    (negative : signature *
      ((L.submatrix id q).det * (R.submatrix q id).det) < 0) :
    signature * (L * sparseFeatureProjection q * R).det < 0 := by
  rw [det_mul_sparseFeatureProjection_mul]
  exact negative

end Matrix
