import GECPKernelStructure.GECP.SignRegularity
import GECPKernelStructure.Matrix.DeterminantIntegral

namespace GECPKernelStructure
namespace GECP

open Matrix MeasureTheory Set

universe u v

variable {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]

/-- Total nonnegativity also controls weakly ordered tuples: repetitions make the minor zero. -/
theorem TotallyNonnegativeAtOrder.det_minorMatrix_nonneg_of_monotone
    {K : Kernel α β ℝ} {n : ℕ} (regular : TotallyNonnegativeAtOrder K n)
    (rows : Fin n → α) (columns : Fin n → β)
    (rows_mono : Monotone rows) (columns_mono : Monotone columns) :
    0 ≤ (minorMatrix K rows columns).det := by
  by_cases rows_injective : Function.Injective rows
  · by_cases columns_injective : Function.Injective columns
    · exact regular rows columns
        (rows_mono.strictMono_of_injective rows_injective)
        (columns_mono.strictMono_of_injective columns_injective)
    · have hzero := det_kernel_eq_zero_of_not_injective_right
        K rows columns columns_injective
      change 0 ≤ Matrix.det (fun i j => K (rows i) (columns j))
      rw [hzero]
  · have hzero := det_kernel_eq_zero_of_not_injective_left
      K rows columns rows_injective
    change 0 ≤ Matrix.det (fun i j => K (rows i) (columns j))
    rw [hzero]

section IteratedRows

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Integrate exactly the rows in `integrated`; sample every other row at `z`. -/
noncomputable def partialRowIntervalIntegral
    (integrated : Finset ι) (f : ι → ℝ → ι → ℝ)
    (lower upper z : ι → ℝ) : Matrix ι ι ℝ :=
  fun i j =>
    if i ∈ integrated then (∫ x in lower i..upper i, f i x) j
    else f i (z i) j

/--
Finite iteration of one-row determinant integration. If every simultaneous
choice of one row from each interval has nonnegative determinant, then the
matrix of row integrals has nonnegative determinant.
-/
theorem det_rowIntervalIntegral_nonneg
    (f : ι → ℝ → ι → ℝ) (lower upper : ι → ℝ)
    (bounds : ∀ i, lower i ≤ upper i)
    (integrable : ∀ i, IntervalIntegrable (f i) volume (lower i) (upper i))
    (pointwise : ∀ z : ι → ℝ,
      (∀ i, z i ∈ Icc (lower i) (upper i)) →
        0 ≤ Matrix.det (fun i j => f i (z i) j)) :
    0 ≤ Matrix.det (fun i j => (∫ x in lower i..upper i, f i x) j) := by
  have partial_nonneg : ∀ integrated : Finset ι, ∀ z : ι → ℝ,
      (∀ i, i ∉ integrated → z i ∈ Icc (lower i) (upper i)) →
        0 ≤ (partialRowIntervalIntegral integrated f lower upper z).det := by
    intro integrated
    induction integrated using Finset.induction_on with
    | empty =>
        intro z hz
        have empty_matrix :
            partialRowIntervalIntegral ∅ f lower upper z =
              fun i j => f i (z i) j := by
          ext i j
          simp [partialRowIntervalIntegral]
        rw [empty_matrix]
        exact pointwise z (fun i => hz i (by simp))
    | @insert i integrated hi ih =>
        intro z hz
        let M := partialRowIntervalIntegral integrated f lower upper z
        have target_matrix :
            partialRowIntervalIntegral (insert i integrated) f lower upper z =
              M.updateRow i (∫ x in lower i..upper i, f i x) := by
          ext r c
          by_cases hri : r = i
          · subst r
            simp [partialRowIntervalIntegral, M, hi]
          · simp [partialRowIntervalIntegral, M, hri]
        rw [target_matrix]
        apply Matrix.det_updateRow_intervalIntegral_nonneg M i (f i)
          (bounds i) (integrable i)
        intro x hx
        have sampled_nonneg := ih (Function.update z i x) (fun r hr => by
          by_cases hri : r = i
          · subst r
            simpa using hx
          · rw [Function.update_of_ne hri]
            exact hz r (by simp [hr, hri]))
        have sampled_matrix :
            M.updateRow i (f i x) =
              partialRowIntervalIntegral integrated f lower upper
                (Function.update z i x) := by
          ext r c
          by_cases hri : r = i
          · subst r
            simp [partialRowIntervalIntegral, M, hi]
          · simp [partialRowIntervalIntegral, M, hri]
        rwa [sampled_matrix]
  have all_nonneg := partial_nonneg Finset.univ (fun _ => 0) (by simp)
  have all_matrix :
      partialRowIntervalIntegral Finset.univ f lower upper (fun _ => 0) =
        fun i j => (∫ x in lower i..upper i, f i x) j := by
    ext i j
    simp [partialRowIntervalIntegral]
  rwa [all_matrix] at all_nonneg

end IteratedRows

/--
Ordered disjoint row intervals preserve a totally nonnegative minor after
integrating the kernel over every row interval.
-/
theorem det_rowIntervalIntegral_kernel_nonneg
    {K : Kernel ℝ β ℝ} {n : ℕ} (regular : TotallyNonnegativeAtOrder K n)
    (lower upper : Fin n → ℝ) (columns : Fin n → β)
    (bounds : ∀ i, lower i ≤ upper i)
    (intervals_ordered : ∀ ⦃i j : Fin n⦄, i < j → upper i ≤ lower j)
    (columns_mono : StrictMono columns)
    (integrable : ∀ i, IntervalIntegrable
      (fun x => fun j => K x (columns j)) volume (lower i) (upper i)) :
    0 ≤ Matrix.det (fun i j => ∫ x in lower i..upper i, K x (columns j)) := by
  have vector_integral_nonneg := det_rowIntervalIntegral_nonneg
    (fun _i x j => K x (columns j)) lower upper bounds integrable
    (fun z hz => by
      apply regular.det_minorMatrix_nonneg_of_monotone z columns
      · intro i j hij
        rcases eq_or_lt_of_le hij with rfl | hlt
        · exact le_rfl
        · exact (hz i).2.trans ((intervals_ordered hlt).trans (hz j).1)
      · exact columns_mono.monotone)
  have integral_matrix :
      (fun i j => ∫ x in lower i..upper i, K x (columns j)) =
        (fun i j => (∫ x in lower i..upper i,
          fun j => K x (columns j)) j) := by
    funext i j
    simpa using
      ((ContinuousLinearMap.proj j :
        (Fin n → ℝ) →L[ℝ] ℝ).intervalIntegral_comp_comm (integrable i))
  rw [integral_matrix]
  exact vector_integral_nonneg

end GECP
end GECPKernelStructure
