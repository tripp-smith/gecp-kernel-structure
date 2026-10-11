import GECPKernelStructure.Fermionic.OrderKernelTotalNonnegative
import GECPKernelStructure.Fermionic.TruncatedPowerConvolution
import GECPKernelStructure.GECP.TotalNonnegativeCommonUpper
import Mathlib.MeasureTheory.SpecificCodomains.Pi

namespace GECPKernelStructure

namespace GECP

open Matrix

universe u v

variable {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]

/-- Multiplication by a nonnegative constant preserves every ordered minor sign. -/
theorem TotallyNonnegativeAtOrder.const_mul
    {K : Kernel α β ℝ} {n : ℕ} (regular : TotallyNonnegativeAtOrder K n)
    (constant : ℝ) (constant_nonneg : 0 ≤ constant) :
    TotallyNonnegativeAtOrder (fun x y => constant * K x y) n := by
  intro rows columns rows_mono columns_mono
  have base_nonneg := regular rows columns rows_mono columns_mono
  have matrix_eq :
      minorMatrix (fun x y => constant * K x y) rows columns =
        constant • minorMatrix K rows columns := by
    ext i j
    rfl
  rw [matrix_eq, Matrix.det_smul]
  exact mul_nonneg (pow_nonneg constant_nonneg _) base_nonneg

end GECP

namespace Fermionic

open GECP Matrix MeasureTheory

/-- Total nonnegativity propagates from one normalized truncated-power degree to the next. -/
theorem normalizedTruncatedPowerKernel_succ_totallyNonnegativeAtOrder
    (degree : ℕ) {order : ℕ}
    (regular : TotallyNonnegativeAtOrder
      (normalizedTruncatedPowerKernel degree) order) :
    TotallyNonnegativeAtOrder
      (normalizedTruncatedPowerKernel (degree + 1)) order := by
  cases order with
  | zero =>
      intro rows columns _rows_mono _columns_mono
      simp
  | succ n =>
      intro rows columns rows_mono columns_mono
      let upper := max (rows (Fin.last n)) (columns (Fin.last n))
      have row_le_upper (i : Fin (n + 1)) : rows i ≤ upper :=
        (rows_mono.monotone (Fin.le_last i)).trans (le_max_left _ _)
      have column_le_upper (j : Fin (n + 1)) : columns j ≤ upper :=
        (columns_mono.monotone (Fin.le_last j)).trans (le_max_right _ _)
      have scalar_integrable (i j : Fin (n + 1)) : IntervalIntegrable
          (fun z => normalizedTruncatedPowerKernel degree z (columns j))
          volume (rows i) upper :=
        normalizedTruncatedPowerKernel_intervalIntegrable_to_commonUpper degree
          (row_le_upper i) (column_le_upper j)
      have vector_integrable (i : Fin (n + 1)) : IntervalIntegrable
          (fun z => fun j => normalizedTruncatedPowerKernel degree z (columns j))
          volume (rows i) upper := by
        rw [intervalIntegrable_iff]
        apply Integrable.of_eval
        intro j
        exact intervalIntegrable_iff.mp (scalar_integrable i j)
      have tails_nonneg := det_commonUpperRowIntegral_kernel_nonneg regular
        rows columns upper rows_mono columns_mono
        (le_max_left _ _) vector_integrable
      have matrix_eq :
          minorMatrix (normalizedTruncatedPowerKernel (degree + 1)) rows columns =
            fun i j => ∫ z in rows i..upper,
              normalizedTruncatedPowerKernel degree z (columns j) := by
        ext i j
        exact normalizedTruncatedPowerKernel_succ_integral_to_commonUpper degree
          (row_le_upper i) (column_le_upper j)
      rw [matrix_eq]
      exact tails_nonneg

/-- Every factorially normalized truncated-power kernel is totally nonnegative. -/
theorem normalizedTruncatedPowerKernel_totallyNonnegative (degree : ℕ) :
    TotallyNonnegative (normalizedTruncatedPowerKernel degree) := by
  induction degree with
  | zero => exact normalizedTruncatedPowerKernel_zero_totallyNonnegative
  | succ degree ih =>
      intro order
      exact normalizedTruncatedPowerKernel_succ_totallyNonnegativeAtOrder degree (ih order)

/-- The seventh-degree truncated-power factor has nonnegative ordered minors of every order. -/
theorem eighthOrderTruncatedPowerKernel_totallyNonnegative :
    TotallyNonnegative eighthOrderTruncatedPowerKernel := by
  intro order
  have normalized := normalizedTruncatedPowerKernel_totallyNonnegative 7 order
  have scaled := normalized.const_mul (Nat.factorial 7 : ℝ) (by positivity)
  have kernel_eq : eighthOrderTruncatedPowerKernel =
      fun v ω => (Nat.factorial 7 : ℝ) * normalizedTruncatedPowerKernel 7 v ω := by
    funext v ω
    exact eighthOrderTruncatedPowerKernel_eq_factorial_mul_normalized v ω
  rw [kernel_eq]
  exact scaled

end Fermionic
end GECPKernelStructure
