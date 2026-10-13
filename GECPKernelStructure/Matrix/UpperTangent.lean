import GECPKernelStructure.Matrix.DoubleBorderedCompanion

namespace Matrix

universe u v w z

variable {row : Type u} {column : Type v} {core : Type w} {fresh : Type z}

/-- The upper-triangular two-layer, or tangent, lift of a base matrix and a
companion matrix. -/
def upperTangentMatrix {𝕜 : Type*} [Zero 𝕜]
    (K H : Matrix row column 𝕜) :
    Matrix (row ⊕ row) (column ⊕ column) 𝕜 :=
  Matrix.fromBlocks K H 0 K

/-- Rows selecting two copies of a core tuple and a fresh tuple in the upper
layer. -/
def upperTangentDoubleBorderRow
    (selectedRows : core → row) (freshRows : fresh → row) :
    (core ⊕ core) ⊕ fresh → row ⊕ row
  | Sum.inl (Sum.inl i) => Sum.inl (selectedRows i)
  | Sum.inl (Sum.inr i) => Sum.inr (selectedRows i)
  | Sum.inr i => Sum.inl (freshRows i)

/-- Columns selecting two copies of a core tuple and a fresh tuple in the
lower layer. -/
def upperTangentDoubleBorderColumn
    (selectedColumns : core → column) (freshColumns : fresh → column) :
    (core ⊕ core) ⊕ fresh → column ⊕ column
  | Sum.inl (Sum.inl j) => Sum.inl (selectedColumns j)
  | Sum.inl (Sum.inr j) => Sum.inr (selectedColumns j)
  | Sum.inr j => Sum.inr (freshColumns j)

/-- The double-bordered companion matrix is one distinguished square minor of
the upper-tangent lift. -/
theorem submatrix_upperTangentMatrix_eq_doubleBorderedCompanionMatrix
    {𝕜 : Type*} [Zero 𝕜]
    (K H : Matrix row column 𝕜)
    (selectedRows : core → row) (selectedColumns : core → column)
    (freshRows : fresh → row) (freshColumns : fresh → column) :
    Matrix.submatrix (upperTangentMatrix K H)
        (upperTangentDoubleBorderRow selectedRows freshRows)
        (upperTangentDoubleBorderColumn selectedColumns freshColumns) =
      doubleBorderedCompanionMatrix
        (Matrix.submatrix K selectedRows selectedColumns)
        (Matrix.submatrix K freshRows selectedColumns)
        (Matrix.submatrix K selectedRows freshColumns)
        (Matrix.submatrix H selectedRows selectedColumns)
        (Matrix.submatrix H freshRows selectedColumns)
        (Matrix.submatrix H selectedRows freshColumns)
        (Matrix.submatrix H freshRows freshColumns) := by
  ext (i | i) (j | j)
  · rcases i with i | i <;> rcases j with j | j <;> rfl
  · rcases i with i | i <;> rfl
  · rcases j with j | j <;> rfl
  · rfl

section Factorization

variable {feature : Type*} [Fintype feature] [DecidableEq feature]
  {𝕜 : Type*} [CommRing 𝕜]

/-- Tangent lifting respects a shared left/right feature factorization. The
companion action is isolated in one upper-unitriangular middle block. -/
theorem upperTangentMatrix_mul_factorization
    (A : Matrix row feature 𝕜) (W : Matrix feature feature 𝕜)
    (B : Matrix feature column 𝕜) :
    upperTangentMatrix (A * B) (A * W * B) =
      Matrix.fromBlocks A 0 0 A *
        upperTangentMatrix (1 : Matrix feature feature 𝕜) W *
          Matrix.fromBlocks B 0 0 B := by
  simp only [upperTangentMatrix, Matrix.fromBlocks_multiply,
    Matrix.mul_one, Matrix.mul_zero, Matrix.zero_mul,
    add_zero, zero_add, Matrix.mul_assoc]

end Factorization

end Matrix
