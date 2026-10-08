import GECPKernelStructure.GECP.BorderedDeterminant
import Mathlib.Data.Fin.Tuple.Sort
import Mathlib.Tactic

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v

variable {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]

/-- A square kernel minor evaluated on finite row and column tuples. -/
def minorMatrix (K : Kernel α β ℝ) {n : ℕ}
    (rows : Fin n → α) (columns : Fin n → β) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => K (rows i) (columns j)

/--
All `n x n` minors on strictly increasing row and column tuples have the
prescribed strict sign. A signature of `1` means positive ordered minors and a
signature of `-1` means negative ordered minors.
-/
def StrictSignRegularAtOrder (K : Kernel α β ℝ) (n : ℕ) (signature : ℝ) : Prop :=
  ∀ (rows : Fin n → α) (columns : Fin n → β),
    StrictMono rows → StrictMono columns →
      0 < signature * (minorMatrix K rows columns).det

/-- A kernel is strictly sign regular when every order has a prescribed sign. -/
def StrictSignRegular (K : Kernel α β ℝ) (signature : ℕ → ℝ) : Prop :=
  ∀ n, StrictSignRegularAtOrder K n (signature n)

/-- The orientation sign of a finite tuple, computed by its sorting permutation. -/
noncomputable def tupleOrientation {γ : Type*} [LinearOrder γ] {n : ℕ}
    (values : Fin n → γ) : ℝ :=
  (((Tuple.sort values).sign : ℤ) : ℝ)

@[simp]
theorem tupleOrientation_sq {γ : Type*} [LinearOrder γ] {n : ℕ}
    (values : Fin n → γ) : tupleOrientation values ^ 2 = 1 := by
  unfold tupleOrientation
  rw [pow_two]
  norm_cast
  exact congrArg Units.val (Int.units_mul_self (Tuple.sort values).sign)

omit [LinearOrder α] [LinearOrder β] in
theorem det_kernel_eq_zero_of_not_injective_left
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (K : Kernel α β ℝ) (rows : ι → α) (columns : ι → β)
    (not_injective : ¬Function.Injective rows) :
    Matrix.det (fun i j => K (rows i) (columns j)) = 0 := by
  rw [Function.Injective] at not_injective
  push Not at not_injective
  obtain ⟨i, j, hij, hne⟩ := not_injective
  apply Matrix.det_zero_of_row_eq hne
  funext k
  exact congrArg (fun value => K value (columns k)) hij

omit [LinearOrder α] [LinearOrder β] in
theorem det_kernel_eq_zero_of_not_injective_right
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (K : Kernel α β ℝ) (rows : ι → α) (columns : ι → β)
    (not_injective : ¬Function.Injective columns) :
    Matrix.det (fun i j => K (rows i) (columns j)) = 0 := by
  rw [Function.Injective] at not_injective
  push Not at not_injective
  obtain ⟨i, j, hij, hne⟩ := not_injective
  apply Matrix.det_zero_of_column_eq hne
  intro k
  exact congrArg (fun value => K (rows k) value) hij

/--
Strict sign regularity on sorted tuples determines the oriented sign of every
minor whose row and column tuples are injective.
-/
theorem oriented_minor_pos_of_strictSignRegularAtOrder
    {K : Kernel α β ℝ} {n : ℕ} {signature : ℝ}
    (regular : StrictSignRegularAtOrder K n signature)
    (rows : Fin n → α) (columns : Fin n → β)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns) :
    0 < signature * tupleOrientation rows * tupleOrientation columns *
      (minorMatrix K rows columns).det := by
  let rowPerm := Tuple.sort rows
  let columnPerm := Tuple.sort columns
  have rows_sorted : StrictMono (rows ∘ rowPerm) :=
    (Tuple.monotone_sort rows).strictMono_of_injective
      (rows_injective.comp rowPerm.injective)
  have columns_sorted : StrictMono (columns ∘ columnPerm) :=
    (Tuple.monotone_sort columns).strictMono_of_injective
      (columns_injective.comp columnPerm.injective)
  have sorted_pos := regular (rows ∘ rowPerm) (columns ∘ columnPerm)
    rows_sorted columns_sorted
  let original : Matrix (Fin n) (Fin n) ℝ :=
    minorMatrix K rows columns
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
      _ = tupleOrientation rows * tupleOrientation columns * original.det := by ring
  rw [sorted_matrix, det_sorted] at sorted_pos
  simpa [original, minorMatrix, mul_assoc] using sorted_pos

/-- Positive column scaling preserves every prescribed ordered-minor sign. -/
theorem strictSignRegularAtOrder_columnScale
    {K : Kernel α β ℝ} {n : ℕ} {signature : ℝ}
    (regular : StrictSignRegularAtOrder K n signature)
    (weight : β → ℝ) (weight_pos : ∀ y, 0 < weight y) :
    StrictSignRegularAtOrder (fun x y => K x y * weight y) n signature := by
  intro rows columns rows_mono columns_mono
  have base := regular rows columns rows_mono columns_mono
  let weights : Fin n → ℝ := fun j => weight (columns j)
  let core : Matrix (Fin n) (Fin n) ℝ := fun i j => K (rows i) (columns j)
  have weights_product_pos : 0 < ∏ j, weights j :=
    Finset.prod_pos fun j _ => weight_pos (columns j)
  have scaled_det :
      (minorMatrix (fun x y => K x y * weight y) rows columns).det =
        (∏ j, weights j) * core.det := by
    rw [show minorMatrix (fun x y => K x y * weight y) rows columns =
      fun i j => weights j * core i j by
        ext i j
        change K (rows i) (columns j) * weight (columns j) =
          weight (columns j) * K (rows i) (columns j)
        ring]
    exact Matrix.det_mul_row weights core
  rw [scaled_det]
  calc
    0 < (∏ j, weights j) * (signature * core.det) :=
      mul_pos weights_product_pos base
    _ = signature * ((∏ j, weights j) * core.det) := by ring

/-- Positive row scaling preserves every prescribed ordered-minor sign. -/
theorem strictSignRegularAtOrder_rowScale
    {K : Kernel α β ℝ} {n : ℕ} {signature : ℝ}
    (regular : StrictSignRegularAtOrder K n signature)
    (weight : α → ℝ) (weight_pos : ∀ x, 0 < weight x) :
    StrictSignRegularAtOrder (fun x y => weight x * K x y) n signature := by
  intro rows columns rows_mono columns_mono
  have base := regular rows columns rows_mono columns_mono
  let weights : Fin n → ℝ := fun i => weight (rows i)
  let core : Matrix (Fin n) (Fin n) ℝ := fun i j => K (rows i) (columns j)
  have weights_product_pos : 0 < ∏ i, weights i :=
    Finset.prod_pos fun i _ => weight_pos (rows i)
  have scaled_det :
      (minorMatrix (fun x y => weight x * K x y) rows columns).det =
        (∏ i, weights i) * core.det := by
    rw [show minorMatrix (fun x y => weight x * K x y) rows columns =
      fun i j => weights i * core i j by rfl]
    exact Matrix.det_mul_column weights core
  rw [scaled_det]
  calc
    0 < (∏ i, weights i) * (signature * core.det) :=
      mul_pos weights_product_pos base
    _ = signature * ((∏ i, weights i) * core.det) := by ring

end GECP
end GECPKernelStructure
