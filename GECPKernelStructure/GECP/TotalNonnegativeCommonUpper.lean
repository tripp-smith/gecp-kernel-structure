import GECPKernelStructure.Matrix.AdjacentRowDifference

namespace GECPKernelStructure
namespace GECP

open Matrix MeasureTheory

universe v

variable {β : Type v} [LinearOrder β]

/--
Common-upper row integrals of a totally nonnegative kernel have nonnegative
ordered minors. Successive row differences reduce them to ordered disjoint
intervals.
-/
theorem det_commonUpperRowIntegral_kernel_nonneg
    {K : Kernel ℝ β ℝ} {n : ℕ}
    (regular : TotallyNonnegativeAtOrder K (n + 1))
    (knots : Fin (n + 1) → ℝ) (columns : Fin (n + 1) → β) (upper : ℝ)
    (knots_mono : StrictMono knots) (columns_mono : StrictMono columns)
    (last_le_upper : knots (Fin.last n) ≤ upper)
    (integrable : ∀ i, IntervalIntegrable
      (fun x => fun j => K x (columns j)) volume (knots i) upper) :
    0 ≤ Matrix.det (fun i j => ∫ x in knots i..upper, K x (columns j)) := by
  let intervalUpper : Fin (n + 1) → ℝ :=
    Fin.lastCases upper (fun i => knots i.succ)
  have bounds : ∀ i, knots i ≤ intervalUpper i := by
    intro i
    refine Fin.lastCases (motive := fun i => knots i ≤ intervalUpper i) ?_
      (fun k => ?_) i
    · simpa [intervalUpper] using last_le_upper
    · simpa [intervalUpper] using (knots_mono Fin.castSucc_lt_succ).le
  have intervals_ordered : ∀ ⦃i j : Fin (n + 1)⦄,
      i < j → intervalUpper i ≤ knots j := by
    intro i j
    refine Fin.lastCases
      (motive := fun i => i < j → intervalUpper i ≤ knots j) ?_
      (fun k => ?_) i <;> intro hij
    · exact False.elim ((not_lt_of_ge (Fin.le_last j)) hij)
    · simpa [intervalUpper] using
        knots_mono.monotone (Fin.castSucc_lt_iff_succ_le.mp hij)
  have adjacent_integrable : ∀ i, IntervalIntegrable
      (fun x => fun j => K x (columns j)) volume (knots i) (intervalUpper i) := by
    intro i
    refine Fin.lastCases (motive := fun i => IntervalIntegrable
      (fun x => fun j => K x (columns j)) volume (knots i) (intervalUpper i))
      ?_ (fun k => ?_) i
    · simpa [intervalUpper] using integrable (Fin.last n)
    · simpa [intervalUpper] using
        (integrable k.castSucc).trans (integrable k.succ).symm
  have adjacent_nonneg := det_rowIntervalIntegral_kernel_nonneg regular
    knots intervalUpper columns bounds intervals_ordered columns_mono adjacent_integrable
  have scalar_integrable : ∀ i j, IntervalIntegrable
      (fun x => K x (columns j)) volume (knots i) upper := by
    intro i j
    constructor
    · simpa [Function.comp_def] using
        (ContinuousLinearMap.proj j :
          (Fin (n + 1) → ℝ) →L[ℝ] ℝ).integrableOn_comp (integrable i).1
    · simpa [Function.comp_def] using
        (ContinuousLinearMap.proj j :
          (Fin (n + 1) → ℝ) →L[ℝ] ℝ).integrableOn_comp (integrable i).2
  rw [Matrix.det_commonUpperRowIntegral_eq_adjacent
    (fun x j => K x (columns j)) knots upper scalar_integrable]
  have adjacent_matrix :
      (fun i j => Fin.lastCases
        (∫ x in knots (Fin.last n)..upper, K x (columns j))
        (fun k => ∫ x in knots k.castSucc..knots k.succ,
          K x (columns j)) i) =
      (fun i j => ∫ x in knots i..intervalUpper i, K x (columns j)) := by
    funext i j
    refine Fin.lastCases ?_ (fun k => ?_) i <;> simp [intervalUpper]
  rw [adjacent_matrix]
  exact adjacent_nonneg

end GECP
end GECPKernelStructure
