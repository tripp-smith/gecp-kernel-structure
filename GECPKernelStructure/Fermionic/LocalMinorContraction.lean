import GECPKernelStructure.Fermionic.PositiveSeriesGECP
import GECPKernelStructure.Matrix.TotalPositiveFischer

namespace GECPKernelStructure
namespace GECP

open Matrix

/-- One column-maximal Schur step preserves a local minor budget at every
order. The old row tuple deletes the pivot column's sorted rank; it need not
delete the pivot row. The columns are exactly the fresh local columns. -/
theorem localMinor_schur_le_of_strictTotalPositive
    {α β : Type*} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} (regular : StrictSignRegular K (fun _ ↦ 1))
    {n : ℕ} (rows : Fin (n + 1) → α) (columns : Fin (n + 1) → β)
    (hr : StrictMono rows) (hc : StrictMono columns) (p q : Fin (n + 1))
    (columnMax : ∀ i, K (rows i) (columns q) ≤ K (rows p) (columns q))
    (hne : K (rows p) (columns q) ≠ 0) :
    |(minorMatrix (residualUpdate K (rows p) (columns q) hne)
        (rows ∘ p.succAbove) (columns ∘ q.succAbove)).det| ≤
      (minorMatrix K (rows ∘ q.succAbove) (columns ∘ q.succAbove)).det := by
  have positive := strictlyTotallyPositive_minorMatrix_of_strictSignRegular_one
    regular rows columns hr hc
  exact Matrix.abs_det_pivotSchur_le_matchedComplement_of_columnMax
    positive p q columnMax

/-- An initial envelope for the same local columns bounds the new minor;
no pivot count or factorial loss is introduced by this one step. -/
theorem localMinor_schur_le_budget_of_strictTotalPositive
    {α β : Type*} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} (regular : StrictSignRegular K (fun _ ↦ 1))
    {n : ℕ} (rows : Fin (n + 1) → α) (columns : Fin (n + 1) → β)
    (hr : StrictMono rows) (hc : StrictMono columns) (p q : Fin (n + 1))
    (columnMax : ∀ i, K (rows i) (columns q) ≤ K (rows p) (columns q))
    (hne : K (rows p) (columns q) ≠ 0) {B : ℝ}
    (budget : ∀ r : Fin n → α,
      |(minorMatrix K r (columns ∘ q.succAbove)).det| ≤ B) :
    |(minorMatrix (residualUpdate K (rows p) (columns q) hne)
        (rows ∘ p.succAbove) (columns ∘ q.succAbove)).det| ≤ B := by
  exact (localMinor_schur_le_of_strictTotalPositive
    regular rows columns hr hc p q columnMax hne).trans
      ((le_abs_self _).trans (budget _))

end GECP

namespace Fermionic

open GECP

/-- The local-minor contraction bridge for the actual positive-series
fermionic kernel, not a finite Taylor surrogate. No unit-band restriction is
needed for this algebraic one-pivot theorem. -/
theorem positiveSeriesFermionicKernel_localMinor_schur_le
    {n : ℕ} (rows columns : Fin (n + 1) → PositiveArgument)
    (hr : StrictMono rows) (hc : StrictMono columns) (p q : Fin (n + 1))
    (columnMax : ∀ i, positiveSeriesFermionicKernel (rows i) (columns q) ≤
      positiveSeriesFermionicKernel (rows p) (columns q))
    (hne : positiveSeriesFermionicKernel (rows p) (columns q) ≠ 0) :
    |(minorMatrix
        (residualUpdate positiveSeriesFermionicKernel (rows p) (columns q) hne)
        (rows ∘ p.succAbove) (columns ∘ q.succAbove)).det| ≤
      (minorMatrix positiveSeriesFermionicKernel
        (rows ∘ q.succAbove) (columns ∘ q.succAbove)).det :=
  localMinor_schur_le_of_strictTotalPositive
    positiveSeriesFermionicKernel_strictSignRegular rows columns hr hc p q columnMax hne

end Fermionic
end GECPKernelStructure
