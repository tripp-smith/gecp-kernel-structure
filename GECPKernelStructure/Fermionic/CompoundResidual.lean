import GECPKernelStructure.Fermionic.ResidualGauge

namespace GECPKernelStructure

open Matrix

universe u v z

namespace GECP.Run

/-- A finite row tuple that remains injective after all selected rows are prepended. -/
structure FreshRows {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) where
  values : γ → α
  injective : Function.Injective (run.augmentedRow γ values)

/-- A finite column tuple that remains injective after all selected columns are prepended. -/
structure FreshColumns {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) where
  values : γ → β
  injective : Function.Injective (run.augmentedColumn γ values)

/-- The fixed-order compound kernel of final-residual minors on fresh tuples. -/
noncomputable def finalResidualMinorKernel {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) [Fintype γ] [DecidableEq γ] :
    Kernel (run.FreshRows γ) (run.FreshColumns γ) ℝ :=
  fun rows columns =>
    Matrix.det (fun i j => run.finalResidual (rows.values i) (columns.values j))

end GECP.Run

namespace Fermionic

open GECP

/-- Every fresh compound entry of a strictly sign-regular final residual is nonzero. -/
theorem strictSignRegular_finalResidualMinorKernel_ne_zero
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows : run.FreshRows γ) (columns : run.FreshColumns γ) :
    run.finalResidualMinorKernel γ rows columns ≠ 0 :=
  strictSignRegular_finalResidual_minor_ne_zero regular run γ
    rows.values columns.values rows.injective columns.injective

/--
At every fixed order, the compound kernel of fresh residual minors has a
balanced bipartite signing.
-/
theorem strictSignRegular_finalResidualMinorKernel_crossProductSignCoherent
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ] :
    CrossProductSignCoherent (run.finalResidualMinorKernel γ) := by
  intro anchorRows anchorColumns rows columns
  let order := Fintype.card (run.AugmentedIndex γ)
  let sign := signature order
  let core := run.selectedCore.det
  let rowFactor : run.FreshRows γ → ℝ :=
    fun tuple => run.augmentedRowOrientation γ tuple.values
  let columnFactor : run.FreshColumns γ → ℝ :=
    fun tuple => run.augmentedColumnOrientation γ tuple.values
  let minor : run.FreshRows γ → run.FreshColumns γ → ℝ :=
    run.finalResidualMinorKernel γ
  have hxy := strictSignRegular_finalResidual_minor_oriented_pos regular run γ
    rows.values columns.values rows.injective columns.injective
  have haa := strictSignRegular_finalResidual_minor_oriented_pos regular run γ
    anchorRows.values anchorColumns.values anchorRows.injective anchorColumns.injective
  have hxa := strictSignRegular_finalResidual_minor_oriented_pos regular run γ
    rows.values anchorColumns.values rows.injective anchorColumns.injective
  have hay := strictSignRegular_finalResidual_minor_oriented_pos regular run γ
    anchorRows.values columns.values anchorRows.injective columns.injective
  change 0 < sign * rowFactor rows * columnFactor columns * core *
    minor rows columns at hxy
  change 0 < sign * rowFactor anchorRows * columnFactor anchorColumns * core *
    minor anchorRows anchorColumns at haa
  change 0 < sign * rowFactor rows * columnFactor anchorColumns * core *
    minor rows anchorColumns at hxa
  change 0 < sign * rowFactor anchorRows * columnFactor columns * core *
    minor anchorRows columns at hay
  have rowFactor_sq (tuple : run.FreshRows γ) : rowFactor tuple ^ 2 = 1 := by
    exact tupleOrientation_sq _
  have columnFactor_sq (tuple : run.FreshColumns γ) : columnFactor tuple ^ 2 = 1 := by
    exact tupleOrientation_sq _
  have sign_ne : sign ≠ 0 := by
    intro sign_zero
    simp [sign_zero] at hxy
  have core_ne : core ≠ 0 := run.selectedCore_det_ne
  have combined := mul_pos (mul_pos hxy haa) (mul_pos hxa hay)
  have factorization :
      (sign * rowFactor rows * columnFactor columns * core * minor rows columns) *
          (sign * rowFactor anchorRows * columnFactor anchorColumns * core *
            minor anchorRows anchorColumns) *
        ((sign * rowFactor rows * columnFactor anchorColumns * core *
            minor rows anchorColumns) *
          (sign * rowFactor anchorRows * columnFactor columns * core *
            minor anchorRows columns)) =
        sign ^ 4 * core ^ 4 *
          ((minor rows columns * minor anchorRows anchorColumns) *
            (minor rows anchorColumns * minor anchorRows columns)) := by
    calc
      _ = sign ^ 4 * core ^ 4 *
          rowFactor rows ^ 2 * rowFactor anchorRows ^ 2 *
          columnFactor columns ^ 2 * columnFactor anchorColumns ^ 2 *
          ((minor rows columns * minor anchorRows anchorColumns) *
            (minor rows anchorColumns * minor anchorRows columns)) := by
        ring
      _ = sign ^ 4 * core ^ 4 *
          ((minor rows columns * minor anchorRows anchorColumns) *
            (minor rows anchorColumns * minor anchorRows columns)) := by
        rw [rowFactor_sq rows, rowFactor_sq anchorRows,
          columnFactor_sq columns, columnFactor_sq anchorColumns]
        ring
  rw [factorization] at combined
  have coefficient_pos : 0 < sign ^ 4 * core ^ 4 :=
    mul_pos ((show Even 4 by decide).pow_pos sign_ne)
      ((show Even 4 by decide).pow_pos core_ne)
  exact ((mul_pos_iff_of_pos_left coefficient_pos).mp combined).le

/-- Every compound kernel has a strictly positive magnitude-one anchor gauge. -/
theorem strictSignRegular_finalResidualMinorKernel_anchoredSignGauge_pos
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (anchorRows rows : run.FreshRows γ)
    (anchorColumns columns : run.FreshColumns γ) :
    0 < anchoredSignGaugeKernel (run.finalResidualMinorKernel γ)
      anchorRows anchorColumns rows columns := by
  apply anchoredSignGaugeKernel_pos_of_crossProductSignCoherent
    (strictSignRegular_finalResidualMinorKernel_crossProductSignCoherent
      regular run γ)
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ rows columns
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ rows anchorColumns
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ anchorRows columns
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ anchorRows anchorColumns

/-- The compound anchor gauge preserves the absolute value of every residual minor. -/
theorem abs_strictSignRegular_finalResidualMinorKernel_anchoredSignGauge
    {α : Type u} {β : Type v} [LinearOrder α] [LinearOrder β]
    {K : Kernel α β ℝ} {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) (run : Run K)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (anchorRows rows : run.FreshRows γ)
    (anchorColumns columns : run.FreshColumns γ) :
    abs (anchoredSignGaugeKernel (run.finalResidualMinorKernel γ)
      anchorRows anchorColumns rows columns) =
      abs (run.finalResidualMinorKernel γ rows columns) := by
  apply abs_anchoredSignGaugeKernel
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ rows anchorColumns
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ anchorRows columns
  · exact strictSignRegular_finalResidualMinorKernel_ne_zero
      regular run γ anchorRows anchorColumns

/-- Every fermionic residual compound kernel has balanced cross-product signs. -/
theorem fermionicKernel_finalResidualMinorKernel_crossProductSignCoherent
    (run : Run fermionicKernel) (γ : Type z) [Fintype γ] [DecidableEq γ] :
    CrossProductSignCoherent (run.finalResidualMinorKernel γ) :=
  strictSignRegular_finalResidualMinorKernel_crossProductSignCoherent
    fermionicKernel_strictSignRegular run γ

/-- Every fermionic residual compound kernel has a positive anchor gauge. -/
theorem fermionicKernel_finalResidualMinorKernel_anchoredSignGauge_pos
    (run : Run fermionicKernel) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (anchorRows rows : run.FreshRows γ)
    (anchorColumns columns : run.FreshColumns γ) :
    0 < anchoredSignGaugeKernel (run.finalResidualMinorKernel γ)
      anchorRows anchorColumns rows columns :=
  strictSignRegular_finalResidualMinorKernel_anchoredSignGauge_pos
    fermionicKernel_strictSignRegular run γ
    anchorRows rows anchorColumns columns

end Fermionic
end GECPKernelStructure
