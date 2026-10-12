import GECPKernelStructure.GECP.CompanionConditioning
import GECPKernelStructure.Fermionic.PositiveSeriesGECP
import GECPKernelStructure.Fermionic.ResidualTailAmplification

namespace GECPKernelStructure
namespace GECP.Run

open Matrix

universe u v z w

/-- Phase BQ applies canonically after transporting an initial feature split
through the fixed row operations of the selected prefix.  The residual low
block still factors through the original feature type `κ`. -/
theorem conditionRows_rank_amplification
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows : γ → α) (columns : γ → β)
    (E : Matrix α β ℝ) (L : Matrix α κ ℝ) (R : Matrix κ β ℝ)
    (split : (K : Matrix α β ℝ) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree
        (Matrix.submatrix (run.conditionRows E) rows columns)
        (Matrix.submatrix (run.conditionRows L) rows id *
          Matrix.submatrix R id columns))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      abs (lambda • Matrix.submatrix (run.conditionRows E) rows columns +
        Matrix.submatrix (run.conditionRows L) rows id *
          Matrix.submatrix R id columns).det := by
  let conditionedTail : Matrix γ γ ℝ :=
    Matrix.submatrix (run.conditionRows E) rows columns
  let conditionedLeft : Matrix γ κ ℝ :=
    Matrix.submatrix (run.conditionRows L) rows id
  let sampledRight : Matrix κ γ ℝ := Matrix.submatrix R id columns
  have residual_split :
      (fun i j => run.finalResidual (rows i) (columns j)) =
        conditionedTail + conditionedLeft * sampledRight := by
    have sampled := run.finalResidual_sample_eq_conditionRows_add_mul
      rows columns E L R split
    ext i j
    have entry := congrFun (congrFun sampled i) j
    change run.finalResidual (rows i) (columns j) = _ at entry
    change run.finalResidual (rows i) (columns j) =
      (conditionedTail + conditionedLeft * sampledRight) i j
    exact entry
  exact run.finalResidual_rank_amplification γ κ rows columns
    conditionedTail conditionedLeft sampledRight residual_split
    coefficient_nonneg rank_lt lambda_one_le

/-- The selected prefix enters the corresponding augmented determinant only
through its exact core determinant; its length does not consume feature rank
or reduce the fresh-sample amplification exponent. -/
theorem conditionRows_augmentedCore_rank_amplification
    {α : Type u} {β : Type v} {K : Kernel α β ℝ}
    (run : Run K) (γ : Type z) [Fintype γ] [DecidableEq γ]
    (κ : Type w) [Fintype κ] [DecidableEq κ]
    (rows : γ → α) (columns : γ → β)
    (E : Matrix α β ℝ) (L : Matrix α κ ℝ) (R : Matrix κ β ℝ)
    (split : (K : Matrix α β ℝ) = E + L * R)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree
        (Matrix.submatrix (run.conditionRows E) rows columns)
        (Matrix.submatrix (run.conditionRows L) rows id *
          Matrix.submatrix R id columns))
    (rank_lt : Fintype.card κ < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - Fintype.card κ) *
        abs (run.augmentedCore γ rows columns).det ≤
      abs run.selectedCore.det *
        abs (lambda • Matrix.submatrix (run.conditionRows E) rows columns +
          Matrix.submatrix (run.conditionRows L) rows id *
            Matrix.submatrix R id columns).det := by
  let conditionedTail : Matrix γ γ ℝ :=
    Matrix.submatrix (run.conditionRows E) rows columns
  let conditionedLeft : Matrix γ κ ℝ :=
    Matrix.submatrix (run.conditionRows L) rows id
  let sampledRight : Matrix κ γ ℝ := Matrix.submatrix R id columns
  have residual_split :
      (fun i j => run.finalResidual (rows i) (columns j)) =
        conditionedTail + conditionedLeft * sampledRight := by
    have sampled := run.finalResidual_sample_eq_conditionRows_add_mul
      rows columns E L R split
    ext i j
    have entry := congrFun (congrFun sampled i) j
    change run.finalResidual (rows i) (columns j) = _ at entry
    change run.finalResidual (rows i) (columns j) =
      (conditionedTail + conditionedLeft * sampledRight) i j
    exact entry
  exact run.augmentedCore_rank_amplification γ κ rows columns
    conditionedTail conditionedLeft sampledRight residual_split
    coefficient_nonneg rank_lt lambda_one_le

end GECP.Run

namespace Fermionic

open GECP Matrix

/-- Matrix-typed view of the positive-series fermionic Taylor tail. -/
noncomputable def positiveSeriesTailMatrix :
    Matrix PositiveArgument PositiveArgument ℝ :=
  positiveSeriesEighthOrderFermionicTailKernel

/-- The eight left monomial features in the positive-series fermionic split. -/
noncomputable def positiveSeriesLowLeft :
    Matrix PositiveArgument (Fin 8) ℝ :=
  fun u degree =>
    (u : ℝ) ^ (degree : ℕ) / ((degree : ℕ).factorial : ℝ)

/-- The matching denominator-scaled right monomial features. -/
noncomputable def positiveSeriesLowRight :
    Matrix (Fin 8) PositiveArgument ℝ :=
  fun degree omega =>
    (omega : ℝ) ^ (degree : ℕ) / (1 + Real.exp (omega : ℝ))

theorem positiveSeriesLow_mul :
    positiveSeriesLowLeft * positiveSeriesLowRight =
      positiveSeriesEighthOrderTaylorApprox := by
  ext u omega
  rw [Matrix.mul_apply]
  simp only [positiveSeriesLowLeft, positiveSeriesLowRight,
    positiveSeriesEighthOrderTaylorApprox, div_eq_mul_inv, expNegTaylor,
    ← Fin.sum_univ_eq_sum_range]
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro degree _
  ring

/-- The reflected fermionic kernel has a global rank-eight-plus-tail feature
split before conditioning. -/
theorem positiveSeriesFermionicKernel_eq_tail_add_low :
    (positiveSeriesFermionicKernel :
      Matrix PositiveArgument PositiveArgument ℝ) =
      positiveSeriesTailMatrix +
        positiveSeriesLowLeft * positiveSeriesLowRight := by
  ext u omega
  have split :=
    fermionicKernel_one_sub_eq_positiveSeriesTaylorApprox_add_tail u omega
  rw [fermionicKernel_one_sub_eq_positive_exp] at split
  change positiveSeriesFermionicKernel u omega =
    positiveSeriesTailMatrix u omega +
      (positiveSeriesLowLeft * positiveSeriesLowRight) u omega
  simp only [positiveSeriesTailMatrix]
  rw [positiveSeriesLow_mul]
  change Real.exp ((u : ℝ) * (omega : ℝ)) /
      (1 + Real.exp (omega : ℝ)) = _
  linarith

/-- Every exact positive-series fermionic prefix canonically retains the
original eight-dimensional right feature space. -/
theorem positiveSeriesFermionicKernel_finalResidual_featureSplit
    (run : Run positiveSeriesFermionicKernel) :
    (run.finalResidual : Matrix PositiveArgument PositiveArgument ℝ) =
      run.conditionRows positiveSeriesTailMatrix +
        run.conditionRows positiveSeriesLowLeft * positiveSeriesLowRight :=
  run.finalResidual_eq_conditionRows_add_mul
    positiveSeriesTailMatrix
    positiveSeriesLowLeft positiveSeriesLowRight
    positiveSeriesFermionicKernel_eq_tail_add_low

/-- With the remaining conditioned mixed-coefficient signs, the canonical
fermionic feature split supplies Phase BQ directly. -/
theorem positiveSeriesFermionicKernel_conditioned_rank_amplification
    (run : Run positiveSeriesFermionicKernel)
    (γ : Type z) [Fintype γ] [DecidableEq γ]
    (rows columns : γ → PositiveArgument)
    (coefficient_nonneg : ∀ degree,
      0 ≤ mixedColumnCoefficient degree
        (Matrix.submatrix
          (run.conditionRows positiveSeriesTailMatrix)
          rows columns)
        (Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix positiveSeriesLowRight id columns))
    (eight_lt : 8 < Fintype.card γ)
    {lambda : ℝ} (lambda_one_le : 1 ≤ lambda) :
    lambda ^ (Fintype.card γ - 8) *
        abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      abs (lambda • Matrix.submatrix
          (run.conditionRows positiveSeriesTailMatrix)
          rows columns +
        Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix positiveSeriesLowRight id columns).det := by
  simpa using run.conditionRows_rank_amplification γ (Fin 8) rows columns
    positiveSeriesTailMatrix
    positiveSeriesLowLeft positiveSeriesLowRight
    positiveSeriesFermionicKernel_eq_tail_add_low coefficient_nonneg
    (by simpa using eight_lt) lambda_one_le

end Fermionic
end GECPKernelStructure
