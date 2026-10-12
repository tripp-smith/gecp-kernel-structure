import GECPKernelStructure.GECP.AugmentedDeterminant
import Mathlib.Order.Fin.InsertNth

namespace Matrix

open Finset
open GECPKernelStructure.GECP

private def finOneSumEquiv (m : ℕ) : Fin 1 ⊕ Fin m ≃ Fin (m + 1) where
  toFun := Sum.elim (fun _ => 0) Fin.succ
  invFun := Fin.cases (Sum.inl 0) Sum.inr
  left_inv x := by
    rcases x with i | i
    · exact Fin.cases rfl (fun j => Fin.elim0 j) i
    · rfl
  right_inv i := by
    exact Fin.cases rfl (fun _ => rfl) i

/-- Every ordered square subminor of a finite real matrix is strictly positive. -/
def StrictlyTotallyPositive {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  ∀ (m : ℕ) (rows columns : Fin m → Fin n),
    StrictMono rows → StrictMono columns →
      0 < (A.submatrix rows columns).det

theorem StrictlyTotallyPositive.entry_pos {n : ℕ}
    {A : Matrix (Fin n) (Fin n) ℝ} (positive : StrictlyTotallyPositive A)
    (i j : Fin n) : 0 < A i j := by
  have minor := positive 1 (fun _ => i) (fun _ => j)
    (by simp [Fin.strictMono_iff_lt_succ])
    (by simp [Fin.strictMono_iff_lt_succ])
  simpa using minor

/-- The scalar northwest Schur complement after eliminating entry `(0,0)`. -/
noncomputable def northwestSchur {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (pivot_ne : A 0 0 ≠ 0) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => residualUpdate A 0 0 pivot_ne i.succ j.succ

@[simp]
theorem northwestSchur_apply {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (pivot_ne : A 0 0 ≠ 0) (i j : Fin n) :
    northwestSchur A pivot_ne i j =
      A i.succ j.succ - A i.succ 0 * A 0 j.succ / A 0 0 := by
  rfl

set_option maxHeartbeats 800000 in
/-- One scalar Schur step factors every minor formed from the northwest pivot
and an arbitrary subminor of the remaining rows and columns. -/
theorem det_cons_succ_eq_pivot_mul_northwestSchur_det
    {n m : ℕ} (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (pivot_ne : A 0 0 ≠ 0) (rows columns : Fin m → Fin n) :
    (A.submatrix (Fin.cons 0 fun i => (rows i).succ)
        (Fin.cons 0 fun j => (columns j).succ)).det =
      A 0 0 * ((northwestSchur A pivot_ne).submatrix rows columns).det := by
  classical
  let residual : GECPKernelStructure.Kernel (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    residualUpdate A 0 0 pivot_ne
  let tail : Run residual := Run.nil residual
  have factor := Run.augmentedCore_step_det (γ := Fin m)
    0 0 pivot_ne tail (fun i => (rows i).succ) (fun j => (columns j).succ)
  let e : Fin 1 ⊕ Fin m ≃ Fin (m + 1) := finOneSumEquiv m
  have reindexed :
      Matrix.reindex e e
          ((Run.step 0 0 pivot_ne tail).augmentedCore
            (Fin m) (fun i => (rows i).succ)
            (fun j => (columns j).succ)) =
        A.submatrix (Fin.cons 0 fun i => (rows i).succ)
          (Fin.cons 0 fun j => (columns j).succ) := by
    ext i j
    refine Fin.cases ?_ (fun i => ?_) i <;>
      refine Fin.cases ?_ (fun j => ?_) j <;>
      rfl
  calc
    (A.submatrix (Fin.cons 0 fun i => (rows i).succ)
        (Fin.cons 0 fun j => (columns j).succ)).det =
        (Matrix.reindex e e
          ((Run.step 0 0 pivot_ne tail).augmentedCore
            (Fin m) (fun i => (rows i).succ)
              (fun j => (columns j).succ))).det :=
      congrArg Matrix.det reindexed.symm
    _ = @Matrix.det
          (Run.AugmentedIndex (Fin m) (Run.step 0 0 pivot_ne tail))
          (Run.augmentedIndexDecidableEq (Fin m) (Run.step 0 0 pivot_ne tail))
          (Run.augmentedIndexFintype (Fin m) (Run.step 0 0 pivot_ne tail))
          ℝ Real.commRing
          ((Run.step 0 0 pivot_ne tail).augmentedCore
            (Fin m) (fun i => (rows i).succ)
              (fun j => (columns j).succ)) :=
      @Matrix.det_reindex_self
        (Run.AugmentedIndex (Fin m) (Run.step 0 0 pivot_ne tail))
        (Fin (m + 1)) inferInstance inferInstance
        (Run.augmentedIndexDecidableEq (Fin m) (Run.step 0 0 pivot_ne tail))
        (Run.augmentedIndexFintype (Fin m) (Run.step 0 0 pivot_ne tail))
        ℝ Real.commRing e _
    _ = A 0 0 * (tail.augmentedCore
          (Fin m) (fun i => (rows i).succ)
            (fun j => (columns j).succ)).det := factor
    _ = A 0 0 * ((northwestSchur A pivot_ne).submatrix rows columns).det := rfl

/-- Strict total positivity is preserved by the scalar northwest Schur step. -/
theorem StrictlyTotallyPositive.northwestSchur {n : ℕ}
    {A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ}
    (positive : StrictlyTotallyPositive A) :
    StrictlyTotallyPositive
      (northwestSchur A (ne_of_gt (positive.entry_pos 0 0))) := by
  let pivot_ne : A 0 0 ≠ 0 := ne_of_gt (positive.entry_pos 0 0)
  intro m rows columns rows_mono columns_mono
  have augmented_rows_mono :
      StrictMono (Fin.cons 0 fun i => (rows i).succ) := by
    rw [Fin.strictMono_cons]
    exact ⟨fun i => Fin.succ_pos (rows i), Fin.strictMono_succ.comp rows_mono⟩
  have augmented_columns_mono :
      StrictMono (Fin.cons 0 fun j => (columns j).succ) := by
    rw [Fin.strictMono_cons]
    exact ⟨fun j => Fin.succ_pos (columns j),
      Fin.strictMono_succ.comp columns_mono⟩
  have augmented_pos := positive (m + 1)
    (Fin.cons 0 fun i => (rows i).succ)
    (Fin.cons 0 fun j => (columns j).succ)
    augmented_rows_mono augmented_columns_mono
  have factor := det_cons_succ_eq_pivot_mul_northwestSchur_det
    A pivot_ne rows columns
  rw [factor] at augmented_pos
  exact (mul_pos_iff_of_pos_left (positive.entry_pos 0 0)).mp augmented_pos

/-- Every Schur-complement diagonal entry is at most the corresponding
original diagonal entry. -/
theorem northwestSchur_diagonal_le {n : ℕ}
    {A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ}
    (positive : StrictlyTotallyPositive A) (i : Fin n) :
    northwestSchur A (ne_of_gt (positive.entry_pos 0 0)) i i ≤
      A i.succ i.succ := by
  rw [northwestSchur_apply]
  exact sub_le_self _ (div_nonneg
    (mul_nonneg (positive.entry_pos i.succ 0).le
      (positive.entry_pos 0 i.succ).le)
    (positive.entry_pos 0 0).le)

/-- Hadamard's diagonal-product inequality for strictly totally positive
finite matrices. -/
theorem det_le_prod_diagonal_of_strictlyTotallyPositive :
    ∀ {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ),
      StrictlyTotallyPositive A → A.det ≤ ∏ i, A i i
  | 0, A, _ => by simp
  | n + 1, A, positive => by
      let pivot_ne : A 0 0 ≠ 0 := ne_of_gt (positive.entry_pos 0 0)
      let S := northwestSchur A pivot_ne
      have schur_positive : StrictlyTotallyPositive S := by
        exact positive.northwestSchur
      have schur_bound : S.det ≤ ∏ i, S i i :=
        det_le_prod_diagonal_of_strictlyTotallyPositive S schur_positive
      have diagonal_product_le : (∏ i, S i i) ≤
          ∏ i : Fin n, A i.succ i.succ := by
        apply Finset.prod_le_prod
        · intro i _
          exact (schur_positive.entry_pos i i).le
        · intro i _
          exact northwestSchur_diagonal_le positive i
      have determinant_factor : A.det = A 0 0 * S.det := by
        dsimp only [S]
        have factor := det_cons_succ_eq_pivot_mul_northwestSchur_det
          A pivot_ne id id
        simpa only [id_eq, Fin.cons_zero_succ, Matrix.submatrix_id_id] using factor
      calc
        A.det = A 0 0 * S.det := determinant_factor
        _ ≤ A 0 0 * (∏ i, S i i) :=
          mul_le_mul_of_nonneg_left schur_bound (positive.entry_pos 0 0).le
        _ ≤ A 0 0 * (∏ i : Fin n, A i.succ i.succ) :=
          mul_le_mul_of_nonneg_left diagonal_product_le
            (positive.entry_pos 0 0).le
        _ = ∏ i, A i i := by rw [Fin.prod_univ_succ]

end Matrix

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v

/-- An ordered sample of a kernel with positive minors is a strictly totally
positive finite matrix. -/
theorem strictlyTotallyPositive_minorMatrix_of_strictSignRegular_one
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ}
    (regular : StrictSignRegular K (fun _ => 1))
    {n : ℕ} (rows : Fin n → α) (columns : Fin n → β)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    Matrix.StrictlyTotallyPositive (minorMatrix K rows columns) := by
  intro m selectedRows selectedColumns selectedRows_mono selectedColumns_mono
  have composedRows_mono : StrictMono (rows ∘ selectedRows) :=
    rows_mono.comp selectedRows_mono
  have composedColumns_mono : StrictMono (columns ∘ selectedColumns) :=
    columns_mono.comp selectedColumns_mono
  have minor_pos := regular m (rows ∘ selectedRows) (columns ∘ selectedColumns)
    composedRows_mono composedColumns_mono
  simp only [one_mul] at minor_pos
  change 0 < Matrix.det (fun i j =>
    K (rows (selectedRows i)) (columns (selectedColumns j))) at minor_pos ⊢
  exact minor_pos

/-- Every ordered sample of a positive-minor kernel is bounded by the product
of any common upper bounds on its paired diagonal entries. -/
theorem minorMatrix_det_le_pow_of_strictSignRegular_one
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ}
    (regular : StrictSignRegular K (fun _ => 1))
    {n : ℕ} (rows : Fin n → α) (columns : Fin n → β)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    {B : ℝ} (diagonal_bound : ∀ i, K (rows i) (columns i) ≤ B) :
    (minorMatrix K rows columns).det ≤ B ^ n := by
  let A := minorMatrix K rows columns
  have positive : Matrix.StrictlyTotallyPositive A :=
    strictlyTotallyPositive_minorMatrix_of_strictSignRegular_one
      regular rows columns rows_mono columns_mono
  have determinant_bound :=
    Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive A positive
  have diagonal_product_bound : (∏ i, A i i) ≤ ∏ _i : Fin n, B := by
    apply Finset.prod_le_prod
    · intro i _
      exact (positive.entry_pos i i).le
    · intro i _
      exact diagonal_bound i
  exact determinant_bound.trans (by simpa using diagonal_product_bound)

/-- An ordered determinant envelope for a positive-minor kernel controls the
absolute determinant of every injective sample, independently of the tuple
order. -/
theorem abs_minorMatrix_det_le_pow_of_strictSignRegular_one
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ}
    (regular : StrictSignRegular K (fun _ => 1))
    {n : ℕ} (rows : Fin n → α) (columns : Fin n → β)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    {B : ℝ} (entry_bound : ∀ i j, K (rows i) (columns j) ≤ B) :
    abs (minorMatrix K rows columns).det ≤ B ^ n := by
  let rowPerm := Tuple.sort rows
  let columnPerm := Tuple.sort columns
  have rows_sorted : StrictMono (rows ∘ rowPerm) :=
    (Tuple.monotone_sort rows).strictMono_of_injective
      (rows_injective.comp rowPerm.injective)
  have columns_sorted : StrictMono (columns ∘ columnPerm) :=
    (Tuple.monotone_sort columns).strictMono_of_injective
      (columns_injective.comp columnPerm.injective)
  have sorted_pos := regular n (rows ∘ rowPerm) (columns ∘ columnPerm)
    rows_sorted columns_sorted
  simp only [one_mul] at sorted_pos
  have sorted_bound := minorMatrix_det_le_pow_of_strictSignRegular_one
    regular (rows ∘ rowPerm) (columns ∘ columnPerm)
    rows_sorted columns_sorted
    (fun i => entry_bound (rowPerm i) (columnPerm i))
  let original := minorMatrix K rows columns
  have sorted_matrix :
      minorMatrix K (rows ∘ rowPerm) (columns ∘ columnPerm) =
        original.submatrix rowPerm columnPerm := by
    rfl
  have det_sorted :
      (original.submatrix rowPerm columnPerm).det =
        tupleOrientation rows * tupleOrientation columns * original.det := by
    calc
      (original.submatrix rowPerm columnPerm).det =
          ((original.submatrix rowPerm id).submatrix id columnPerm).det := by
            rfl
      _ = tupleOrientation columns * (original.submatrix rowPerm id).det := by
        rw [Matrix.det_permute']
        rfl
      _ = tupleOrientation columns * (tupleOrientation rows * original.det) := by
        rw [Matrix.det_permute]
        rfl
      _ = tupleOrientation rows * tupleOrientation columns * original.det := by
        ring
  have row_orientation_abs : abs (tupleOrientation rows) = 1 := by
    have square := tupleOrientation_sq rows
    rw [← sq_abs] at square
    nlinarith [abs_nonneg (tupleOrientation rows)]
  have column_orientation_abs : abs (tupleOrientation columns) = 1 := by
    have square := tupleOrientation_sq columns
    rw [← sq_abs] at square
    nlinarith [abs_nonneg (tupleOrientation columns)]
  have abs_eq : abs original.det =
      (minorMatrix K (rows ∘ rowPerm) (columns ∘ columnPerm)).det := by
    calc
      abs original.det =
          abs (tupleOrientation rows * tupleOrientation columns * original.det) := by
            rw [abs_mul, abs_mul, row_orientation_abs, column_orientation_abs]
            ring
      _ = abs (original.submatrix rowPerm columnPerm).det :=
        congrArg abs det_sorted.symm
      _ = abs (minorMatrix K (rows ∘ rowPerm) (columns ∘ columnPerm)).det := by
        rw [sorted_matrix]
      _ = (minorMatrix K (rows ∘ rowPerm) (columns ∘ columnPerm)).det :=
        abs_of_pos sorted_pos
  rw [abs_eq]
  exact sorted_bound

end GECP
end GECPKernelStructure
