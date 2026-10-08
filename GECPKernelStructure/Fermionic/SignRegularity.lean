import GECPKernelStructure.Fermionic.GECPBounds

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The exponential interaction underlying the positively column-scaled fermionic kernel. -/
noncomputable def expKernel (t ω : ℝ) : ℝ :=
  Real.exp (-t * ω)

/-- The predicted sign of an ordered exponential minor of order `n`. -/
def expKernelSignature (n : ℕ) : ℝ :=
  (-1 : ℝ) ^ n.choose 2

/-- The exponential interaction has positive ordered minors at order one. -/
theorem expKernel_strictSignRegularAtOrder_one :
    StrictSignRegularAtOrder expKernel 1 (expKernelSignature 1) := by
  change StrictSignRegularAtOrder expKernel 1 1
  intro rows columns _rows_mono _columns_mono
  rw [Matrix.det_fin_one]
  change 0 < 1 * Real.exp (-rows 0 * columns 0)
  positivity

/-- Every ordered `2 x 2` minor of `exp(-t*ω)` is strictly negative. -/
theorem expKernel_strictSignRegularAtOrder_two :
    StrictSignRegularAtOrder expKernel 2 (expKernelSignature 2) := by
  rw [show expKernelSignature 2 = -1 by norm_num [expKernelSignature, Nat.choose]]
  intro rows columns rows_mono columns_mono
  have rows_lt : rows 0 < rows 1 := rows_mono (by decide)
  have columns_lt : columns 0 < columns 1 := columns_mono (by decide)
  rw [Matrix.det_fin_two]
  change 0 < -1 *
    (Real.exp (-rows 0 * columns 0) * Real.exp (-rows 1 * columns 1) -
      Real.exp (-rows 0 * columns 1) * Real.exp (-rows 1 * columns 0))
  rw [← Real.exp_add, ← Real.exp_add]
  have exponent_lt :
      -rows 0 * columns 0 + -rows 1 * columns 1 <
        -rows 0 * columns 1 + -rows 1 * columns 0 := by
    nlinarith [mul_pos (sub_pos.mpr rows_lt) (sub_pos.mpr columns_lt)]
  have exponential_lt := Real.exp_lt_exp.mpr exponent_lt
  nlinarith

/-- The fermionic kernel is a positive frequency-only scaling of `expKernel`. -/
theorem fermionicKernel_eq_expKernel_mul_weight (t ω : ℝ) :
    fermionicKernel t ω = expKernel t ω * (1 + Real.exp (-ω))⁻¹ := by
  simp [fermionicKernel, expKernel, div_eq_mul_inv]

/-- Positive denominator scaling transfers every exponential ordered-minor sign. -/
theorem fermionicKernel_strictSignRegularAtOrder_of_expKernel
    {n : ℕ} {signature : ℝ}
    (regular : StrictSignRegularAtOrder expKernel n signature) :
    StrictSignRegularAtOrder fermionicKernel n signature := by
  have scaled := strictSignRegularAtOrder_columnScale regular
    (fun ω => (1 + Real.exp (-ω))⁻¹)
    (fun ω => inv_pos.mpr (fermionicKernel_denominator_pos ω))
  intro rows columns rows_mono columns_mono
  have result := scaled rows columns rows_mono columns_mono
  change 0 < signature * Matrix.det
    (fun i j => fermionicKernel (rows i) (columns j))
  change 0 < signature * Matrix.det
    (fun i j => expKernel (rows i) (columns j) *
      (1 + Real.exp (-columns j))⁻¹) at result
  simpa only [fermionicKernel, expKernel, div_eq_mul_inv] using result

/-- The fermionic kernel has positive ordered minors at order one. -/
theorem fermionicKernel_strictSignRegularAtOrder_one :
    StrictSignRegularAtOrder fermionicKernel 1 (expKernelSignature 1) :=
  fermionicKernel_strictSignRegularAtOrder_of_expKernel
    expKernel_strictSignRegularAtOrder_one

/-- Every ordered `2 x 2` fermionic minor is strictly negative. -/
theorem fermionicKernel_strictSignRegularAtOrder_two :
    StrictSignRegularAtOrder fermionicKernel 2 (expKernelSignature 2) :=
  fermionicKernel_strictSignRegularAtOrder_of_expKernel
    expKernel_strictSignRegularAtOrder_two

/-- All-orders exponential strict sign regularity transfers to the fermionic kernel. -/
theorem fermionicKernel_strictSignRegular_of_expKernel
    {signature : ℕ → ℝ} (regular : StrictSignRegular expKernel signature) :
    StrictSignRegular fermionicKernel signature := by
  intro n
  exact fermionicKernel_strictSignRegularAtOrder_of_expKernel (regular n)

/--
The remaining analytic exponential sign theorem is sufficient for exact
selected-pivot sign coherence of every fermionic GECP run.
-/
theorem fermionicKernel_pivotCrossProductSignCoherent_of_expKernel_strictSignRegular
    (regular : StrictSignRegular expKernel expKernelSignature)
    (run : Run fermionicKernel) (row column : ℝ) :
    PivotCrossProductSignCoherent run.finalResidual row column :=
  strictSignRegular_pivotCrossProductSignCoherent
    (fermionicKernel_strictSignRegular_of_expKernel regular) run row column

end Fermionic
end GECPKernelStructure
