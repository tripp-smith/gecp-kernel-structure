import GECPKernelStructure.Matrix.CauchyBinet
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Topology.Algebra.Module.FiniteDimension

namespace Matrix

open MeasureTheory Set

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- With every other row fixed, the determinant is a continuous linear functional of one row. -/
noncomputable def detUpdateRowCLM (M : Matrix n n ℝ) (i : n) :
    (n → ℝ) →L[ℝ] ℝ :=
  LinearMap.toContinuousLinearMap
    { toFun := fun row => (M.updateRow i row).det
      map_add' := fun left right => det_updateRow_add M i left right
      map_smul' := fun scalar row => by
        simpa only [RingHom.id_apply, smul_eq_mul] using
          det_updateRow_smul M i scalar row }

@[simp]
theorem detUpdateRowCLM_apply (M : Matrix n n ℝ) (i : n) (row : n → ℝ) :
    detUpdateRowCLM M i row = (M.updateRow i row).det :=
  rfl

/-- Updating one row by an interval integral commutes exactly with the determinant. -/
theorem det_updateRow_intervalIntegral (M : Matrix n n ℝ) (i : n)
    (f : ℝ → n → ℝ) {a b : ℝ}
    (hf : IntervalIntegrable f volume a b) :
    (M.updateRow i (∫ x in a..b, f x)).det =
      ∫ x in a..b, (M.updateRow i (f x)).det := by
  symm
  simpa only [detUpdateRowCLM_apply] using
    (detUpdateRowCLM M i).intervalIntegral_comp_comm hf

/-- Pointwise determinant nonnegativity transfers through one row integral. -/
theorem det_updateRow_intervalIntegral_nonneg (M : Matrix n n ℝ) (i : n)
    (f : ℝ → n → ℝ) {a b : ℝ} (hab : a ≤ b)
    (hf : IntervalIntegrable f volume a b)
    (hdet : ∀ x ∈ Icc a b, 0 ≤ (M.updateRow i (f x)).det) :
    0 ≤ (M.updateRow i (∫ x in a..b, f x)).det := by
  rw [det_updateRow_intervalIntegral M i f hf]
  exact intervalIntegral.integral_nonneg hab hdet

/-- A continuous nonnegative determinant path that is positive somewhere transfers strictness. -/
theorem det_updateRow_intervalIntegral_pos (M : Matrix n n ℝ) (i : n)
    (f : ℝ → n → ℝ) {a b : ℝ} (hab : a < b) (hf : Continuous f)
    (hdet : ∀ x ∈ Ioc a b, 0 ≤ (M.updateRow i (f x)).det)
    (hpos : ∃ x ∈ Icc a b, 0 < (M.updateRow i (f x)).det) :
    0 < (M.updateRow i (∫ x in a..b, f x)).det := by
  have hfint : IntervalIntegrable f volume a b := hf.intervalIntegrable a b
  rw [det_updateRow_intervalIntegral M i f hfint]
  apply intervalIntegral.integral_pos hab
  · exact (detUpdateRowCLM M i).continuous.comp hf |>.continuousOn
  · exact hdet
  · exact hpos

end Matrix
