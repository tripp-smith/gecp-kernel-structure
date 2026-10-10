import GECPKernelStructure.Fermionic.ResidualMinors
import Mathlib.Data.Real.Sign

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

universe u v z

/-- The anchor-column entry used as a row scaling weight. -/
def anchoredRowWeight {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (anchorColumn : β) (row : α) : ℝ :=
  R row anchorColumn

/-- The anchor-row entry, corrected by the anchor sign, used as a column weight. -/
def anchoredColumnWeight {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (anchorRow : α) (anchorColumn : β) (column : β) : ℝ :=
  R anchorRow column * R anchorRow anchorColumn

/-- Raw diagonal scaling determined by one anchor row and column. -/
def anchoredGaugeKernel {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (anchorRow : α) (anchorColumn : β) : Kernel α β ℝ :=
  fun row column =>
    anchoredRowWeight R anchorColumn row * R row column *
      anchoredColumnWeight R anchorRow anchorColumn column

/-- Magnitude-one version of the anchor scaling. -/
noncomputable def anchoredSignGaugeKernel {α : Type u} {β : Type v}
    (R : Kernel α β ℝ) (anchorRow : α) (anchorColumn : β) : Kernel α β ℝ :=
  fun row column =>
    Real.sign (anchoredRowWeight R anchorColumn row) * R row column *
      Real.sign (anchoredColumnWeight R anchorRow anchorColumn column)

/-- Cross-product sign coherence makes the raw anchored scaling positive. -/
theorem anchoredGaugeKernel_pos_of_crossProductSignCoherent
    {α : Type u} {β : Type v} {R : Kernel α β ℝ}
    (coherent : CrossProductSignCoherent R)
    {anchorRow row : α} {anchorColumn column : β}
    (entry_ne : R row column ≠ 0)
    (row_anchor_ne : R row anchorColumn ≠ 0)
    (anchor_column_ne : R anchorRow column ≠ 0)
    (anchor_ne : R anchorRow anchorColumn ≠ 0) :
    0 < anchoredGaugeKernel R anchorRow anchorColumn row column := by
  have same_sign := coherent anchorRow anchorColumn row column
  have product_ne :
      (R row column * R anchorRow anchorColumn) *
          (R row anchorColumn * R anchorRow column) ≠ 0 :=
    mul_ne_zero (mul_ne_zero entry_ne anchor_ne)
      (mul_ne_zero row_anchor_ne anchor_column_ne)
  rw [show anchoredGaugeKernel R anchorRow anchorColumn row column =
      (R row column * R anchorRow anchorColumn) *
        (R row anchorColumn * R anchorRow column) by
      simp only [anchoredGaugeKernel, anchoredRowWeight, anchoredColumnWeight]
      ring]
  exact lt_of_le_of_ne same_sign product_ne.symm

private theorem real_sign_sq_of_ne_zero (x : ℝ) (x_ne : x ≠ 0) :
    Real.sign x ^ 2 = 1 := by
  rcases Real.sign_apply_eq_of_ne_zero x x_ne with sign_neg | sign_pos
  · rw [sign_neg]
    norm_num
  · rw [sign_pos]
    norm_num

/-- The magnitude-one anchor gauge is positive wherever all four entries are nonzero. -/
theorem anchoredSignGaugeKernel_pos_of_crossProductSignCoherent
    {α : Type u} {β : Type v} {R : Kernel α β ℝ}
    (coherent : CrossProductSignCoherent R)
    {anchorRow row : α} {anchorColumn column : β}
    (entry_ne : R row column ≠ 0)
    (row_anchor_ne : R row anchorColumn ≠ 0)
    (anchor_column_ne : R anchorRow column ≠ 0)
    (anchor_ne : R anchorRow anchorColumn ≠ 0) :
    0 < anchoredSignGaugeKernel R anchorRow anchorColumn row column := by
  let rowWeight := anchoredRowWeight R anchorColumn row
  let columnWeight := anchoredColumnWeight R anchorRow anchorColumn column
  have rowWeight_ne : rowWeight ≠ 0 := row_anchor_ne
  have columnWeight_ne : columnWeight ≠ 0 :=
    mul_ne_zero anchor_column_ne anchor_ne
  have raw_pos := anchoredGaugeKernel_pos_of_crossProductSignCoherent coherent
    entry_ne row_anchor_ne anchor_column_ne anchor_ne
  have row_factor_pos : 0 < Real.sign rowWeight * rowWeight :=
    Real.sign_mul_pos_of_ne_zero rowWeight rowWeight_ne
  have column_factor_pos : 0 < Real.sign columnWeight * columnWeight :=
    Real.sign_mul_pos_of_ne_zero columnWeight columnWeight_ne
  have factorization :
      (Real.sign rowWeight * rowWeight) *
          (Real.sign columnWeight * columnWeight) *
            anchoredSignGaugeKernel R anchorRow anchorColumn row column =
        anchoredGaugeKernel R anchorRow anchorColumn row column := by
    calc
      _ = Real.sign rowWeight ^ 2 * Real.sign columnWeight ^ 2 *
          (rowWeight * R row column * columnWeight) := by
        simp only [anchoredSignGaugeKernel]
        ring
      _ = rowWeight * R row column * columnWeight := by
        rw [real_sign_sq_of_ne_zero rowWeight rowWeight_ne,
          real_sign_sq_of_ne_zero columnWeight columnWeight_ne]
        ring
      _ = anchoredGaugeKernel R anchorRow anchorColumn row column := by
        rfl
  rw [← factorization] at raw_pos
  exact (mul_pos_iff_of_pos_left (mul_pos row_factor_pos column_factor_pos)).mp raw_pos

/-- Sign gauging changes no entry magnitude. -/
theorem abs_anchoredSignGaugeKernel
    {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (anchorRow row : α) (anchorColumn column : β)
    (row_anchor_ne : R row anchorColumn ≠ 0)
    (anchor_column_ne : R anchorRow column ≠ 0)
    (anchor_ne : R anchorRow anchorColumn ≠ 0) :
    abs (anchoredSignGaugeKernel R anchorRow anchorColumn row column) =
      abs (R row column) := by
  have rowWeight_ne : anchoredRowWeight R anchorColumn row ≠ 0 := row_anchor_ne
  have columnWeight_ne :
      anchoredColumnWeight R anchorRow anchorColumn column ≠ 0 :=
    mul_ne_zero anchor_column_ne anchor_ne
  have sign_abs (x : ℝ) (x_ne : x ≠ 0) : abs (Real.sign x) = 1 := by
    rcases Real.sign_apply_eq_of_ne_zero x x_ne with sign_neg | sign_pos
    · rw [sign_neg]
      norm_num
    · rw [sign_pos]
      norm_num
  simp only [anchoredSignGaugeKernel, abs_mul,
    sign_abs _ rowWeight_ne, sign_abs _ columnWeight_ne, one_mul, mul_one]

/-- Every finite minor receives exactly the products of its row and column gauges. -/
theorem anchoredSignGaugeKernel_minor_det
    {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (anchorRow : α) (anchorColumn : β)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows : γ → α) (columns : γ → β) :
    Matrix.det
        (fun i j => anchoredSignGaugeKernel R anchorRow anchorColumn
          (rows i) (columns j)) =
      (∏ i, Real.sign (anchoredRowWeight R anchorColumn (rows i))) *
        (∏ j, Real.sign
          (anchoredColumnWeight R anchorRow anchorColumn (columns j))) *
        Matrix.det (fun i j => R (rows i) (columns j)) := by
  let rowWeights : γ → ℝ :=
    fun i => Real.sign (anchoredRowWeight R anchorColumn (rows i))
  let columnWeights : γ → ℝ :=
    fun j => Real.sign
      (anchoredColumnWeight R anchorRow anchorColumn (columns j))
  let core : Matrix γ γ ℝ := fun i j => R (rows i) (columns j)
  let columnScaled : Matrix γ γ ℝ := fun i j => columnWeights j * core i j
  have matrix_eq :
      (fun i j => anchoredSignGaugeKernel R anchorRow anchorColumn
        (rows i) (columns j)) =
      (fun i j => rowWeights i * columnScaled i j) := by
    ext i j
    simp only [anchoredSignGaugeKernel, rowWeights, columnWeights, columnScaled, core]
    ring
  have columnScaled_det :
      columnScaled.det = (∏ j, columnWeights j) * core.det := by
    change Matrix.det (fun i j : γ => columnWeights j * core i j) =
      (∏ j, columnWeights j) * core.det
    exact Matrix.det_mul_row columnWeights core
  rw [matrix_eq]
  calc
    Matrix.det (fun i j : γ => rowWeights i * columnScaled i j) =
        (∏ i, rowWeights i) * columnScaled.det :=
      Matrix.det_mul_column rowWeights columnScaled
    _ = (∏ i, rowWeights i) * ((∏ j, columnWeights j) * core.det) := by
      rw [columnScaled_det]
    _ = _ := by
      simp only [rowWeights, columnWeights, core]
      ring

end Fermionic
namespace GECP.Run

/-- A proposed row is distinct from every row already selected by the run. -/
def FreshRow {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (row : α) : Prop :=
  Function.Injective (run.augmentedRow (Fin 1) (fun _ => row))

/-- A proposed column is distinct from every column already selected by the run. -/
def FreshColumn {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (column : β) : Prop :=
  Function.Injective (run.augmentedColumn (Fin 1) (fun _ => column))

end GECP.Run

namespace Fermionic

open GECP Matrix

/-- A fresh row and column give a nonzero final-residual entry. -/
theorem strictSignRegular_finalResidual_ne_zero_of_fresh
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    {row : α} {column : β}
    (row_fresh : run.FreshRow row) (column_fresh : run.FreshColumn column) :
    run.finalResidual row column ≠ 0 := by
  let minor : Matrix (Fin 1) (Fin 1) ℝ :=
    fun _ _ => run.finalResidual row column
  have minor_ne : minor.det ≠ 0 := by
    simpa only [minor] using
      (strictSignRegular_finalResidual_minor_ne_zero regular run
        (Fin 1) (fun _ => row) (fun _ => column) row_fresh column_fresh)
  rw [Matrix.det_fin_one minor] at minor_ne
  exact minor_ne

/--
Every strictly sign-regular final residual has a positive, magnitude-preserving
anchor gauge on fresh coordinates.
-/
theorem strictSignRegular_finalResidual_anchoredSignGauge_pos
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    {anchorRow row : α} {anchorColumn column : β}
    (anchorRow_fresh : run.FreshRow anchorRow)
    (row_fresh : run.FreshRow row)
    (anchorColumn_fresh : run.FreshColumn anchorColumn)
    (column_fresh : run.FreshColumn column) :
    0 < anchoredSignGaugeKernel run.finalResidual anchorRow anchorColumn row column := by
  apply anchoredSignGaugeKernel_pos_of_crossProductSignCoherent
    (fun pivotRow pivotColumn =>
      strictSignRegular_pivotCrossProductSignCoherent
        regular run pivotRow pivotColumn)
  · exact strictSignRegular_finalResidual_ne_zero_of_fresh
      regular run row_fresh column_fresh
  · exact strictSignRegular_finalResidual_ne_zero_of_fresh
      regular run row_fresh anchorColumn_fresh
  · exact strictSignRegular_finalResidual_ne_zero_of_fresh
      regular run anchorRow_fresh column_fresh
  · exact strictSignRegular_finalResidual_ne_zero_of_fresh
      regular run anchorRow_fresh anchorColumn_fresh

/-- Fermionic final residuals have the positive anchor gauge on all fresh entries. -/
theorem fermionicKernel_finalResidual_anchoredSignGauge_pos
    (run : Run fermionicKernel)
    {anchorRow row anchorColumn column : ℝ}
    (anchorRow_fresh : run.FreshRow anchorRow)
    (row_fresh : run.FreshRow row)
    (anchorColumn_fresh : run.FreshColumn anchorColumn)
    (column_fresh : run.FreshColumn column) :
    0 < anchoredSignGaugeKernel run.finalResidual anchorRow anchorColumn row column :=
  strictSignRegular_finalResidual_anchoredSignGauge_pos
    fermionicKernel_strictSignRegular run
    anchorRow_fresh row_fresh anchorColumn_fresh column_fresh

end Fermionic
end GECPKernelStructure
