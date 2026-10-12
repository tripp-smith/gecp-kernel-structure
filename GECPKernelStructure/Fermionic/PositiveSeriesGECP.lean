import GECPKernelStructure.Fermionic.AmplifiedKernelEnvelope
import GECPKernelStructure.Fermionic.GeometricMean

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The reflected positive-series presentation of the fermionic kernel. -/
noncomputable def positiveSeriesFermionicKernel :
    Kernel PositiveArgument PositiveArgument ℝ :=
  fun u omega =>
    Real.exp ((u : ℝ) * (omega : ℝ)) / (1 + Real.exp (omega : ℝ))

/-- Unit tail weight recovers the unamplified positive-series kernel. -/
theorem positiveSeriesFermionicKernel_eq_amplified_one :
    positiveSeriesFermionicKernel =
      positiveSeriesAmplifiedFermionicKernel 1 := by
  funext u omega
  simp only [positiveSeriesFermionicKernel,
    positiveSeriesAmplifiedFermionicKernel, positiveEighthOrderTail,
    eighthOrderTail, neg_neg, one_mul, div_eq_mul_inv]
  ring

/-- The reflected positive-series kernel is strictly totally positive. -/
theorem positiveSeriesFermionicKernel_strictSignRegular :
    StrictSignRegular positiveSeriesFermionicKernel (fun _ => 1) := by
  rw [positiveSeriesFermionicKernel_eq_amplified_one]
  exact positiveSeriesAmplifiedFermionicKernel_strictSignRegular (by norm_num)

/-- Sampling the positive-series kernel gives its canonical finite matrix. -/
theorem positiveSeriesFermionicKernel_minorMatrix
    {order : ℕ} (rows columns : Fin order → PositiveArgument) :
    minorMatrix positiveSeriesFermionicKernel rows columns =
      positiveSeriesFermionicMatrix
        (fun i => (rows i : ℝ)) (fun j => (columns j : ℝ)) := by
  ext i j
  simp only [minorMatrix, positiveSeriesFermionicKernel,
    positiveSeriesFermionicMatrix, expMatrix, div_eq_mul_inv]

/-- Phase BM's ordered determinant decay is invariant under independent
permutations of injective row and column samples. -/
theorem positiveSeriesFermionicKernel_abs_minor_det_le
    {order : ℕ} (rows columns : Fin order → PositiveArgument)
    (eight_le : 8 ≤ order)
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns)
    (rows_le_one : ∀ i, (rows i : ℝ) ≤ 1)
    (columns_le_one : ∀ j, (columns j : ℝ) ≤ 1) :
    abs (minorMatrix positiveSeriesFermionicKernel rows columns).det ≤
      (2 : ℝ) ^ order /
        (Nat.factorial 8 : ℝ) ^ (order - 8) := by
  let rowPerm := Tuple.sort rows
  let columnPerm := Tuple.sort columns
  have rows_sorted : StrictMono (rows ∘ rowPerm) :=
    (Tuple.monotone_sort rows).strictMono_of_injective
      (rows_injective.comp rowPerm.injective)
  have columns_sorted : StrictMono (columns ∘ columnPerm) :=
    (Tuple.monotone_sort columns).strictMono_of_injective
      (columns_injective.comp columnPerm.injective)
  have ordered_bound := positiveSeriesFermionic_det_le_two_pow_div_factorial_pow
    (rows ∘ rowPerm) (columns ∘ columnPerm) eight_le
    rows_sorted columns_sorted
    (fun i => rows_le_one (rowPerm i))
    (fun j => columns_le_one (columnPerm j))
  have ordered_pos := positiveSeriesFermionicKernel_strictSignRegular order
    (rows ∘ rowPerm) (columns ∘ columnPerm) rows_sorted columns_sorted
  simp only [one_mul] at ordered_pos
  let original := minorMatrix positiveSeriesFermionicKernel rows columns
  have sorted_matrix :
      minorMatrix positiveSeriesFermionicKernel
          (rows ∘ rowPerm) (columns ∘ columnPerm) =
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
      _ = tupleOrientation rows * tupleOrientation columns * original.det := by
        ring
  have row_orientation_abs : abs (tupleOrientation rows) = 1 := by
    have square := tupleOrientation_sq rows
    rw [← sq_abs] at square
    nlinarith [abs_nonneg (tupleOrientation rows)]
  have column_orientation_abs : abs (tupleOrientation columns) = 1 := by
    have square := tupleOrientation_sq columns
    rw [← sq_abs] at square
    nlinarith [abs_nonneg (tupleOrientation columns)]
  calc
    abs (minorMatrix positiveSeriesFermionicKernel rows columns).det =
        abs original.det := rfl
    _ = abs (tupleOrientation rows * tupleOrientation columns * original.det) := by
      rw [abs_mul, abs_mul, row_orientation_abs, column_orientation_abs]
      ring
    _ = abs (original.submatrix rowPerm columnPerm).det :=
      congrArg abs det_sorted.symm
    _ = abs (minorMatrix positiveSeriesFermionicKernel
        (rows ∘ rowPerm) (columns ∘ columnPerm)).det := by
      rw [sorted_matrix]
    _ = (positiveSeriesFermionicMatrix
        (fun i => ((rows ∘ rowPerm) i : ℝ))
        (fun j => ((columns ∘ columnPerm) j : ℝ))).det := by
      rw [positiveSeriesFermionicKernel_minorMatrix, abs_of_pos]
      simpa only [positiveSeriesFermionicKernel_minorMatrix] using ordered_pos
    _ ≤ (2 : ℝ) ^ order /
        (Nat.factorial 8 : ℝ) ^ (order - 8) := ordered_bound

private theorem card_selectedIndex_eq_prefixLength
    {K : Kernel PositiveArgument PositiveArgument ℝ}
    (residual : ℕ → Kernel PositiveArgument PositiveArgument ℝ)
    (rows columns : ℕ → PositiveArgument) (n : ℕ) (run : Run K)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n) :
    Fintype.card run.SelectedIndex = n := by
  have hpivot_length : run.pivots.length = n := by
    rw [run_pivots]
    simp
  calc
    Fintype.card run.SelectedIndex = Fintype.card (Fin run.pivots.length) :=
      Fintype.card_congr run.selectedIndexEquivFin
    _ = run.pivots.length := Fintype.card_fin _
    _ = n := hpivot_length

/-- Complete-pivot residual powers on the positive unit band inherit the
fixed-base determinant decay from Phase BM. -/
theorem positiveSeriesFermionicKernel_gecp_error_pow_le
    (residual : ℕ → Kernel PositiveArgument PositiveArgument ℝ)
    (rows columns : ℕ → PositiveArgument)
    (realized : ∀ m, ∃ run : Run positiveSeriesFermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun omega : PositiveArgument => (omega : ℝ) ≤ 1)
      (residual m) (rows m) (columns m))
    (n : ℕ) (eight_le : 8 ≤ n) (run : Run positiveSeriesFermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns n)
    (run_complete : run.CompleteOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun omega : PositiveArgument => (omega : ℝ) ≤ 1))
    (u omega : PositiveArgument) (hu : (u : ℝ) ≤ 1)
    (homega : (omega : ℝ) ≤ 1) :
    abs (residual n u omega) ^ n ≤
      (2 : ℝ) ^ n / (Nat.factorial 8 : ℝ) ^ (n - 8) := by
  have residual_power := strictSignRegular_gecp_error_pow_le_selectedCore_det
    positiveSeriesFermionicKernel_strictSignRegular
    (fun z : PositiveArgument => (z : ℝ) ≤ 1)
    (fun z : PositiveArgument => (z : ℝ) ≤ 1)
    residual rows columns realized pivot_ne updates complete n run run_pivots
    u omega hu homega
  have pivot_length : run.pivots.length = n := by
    rw [run_pivots]
    simp
  let finiteRows : Fin run.pivots.length → PositiveArgument :=
    run.selectedRow ∘ run.selectedIndexEquivFin.symm
  let finiteColumns : Fin run.pivots.length → PositiveArgument :=
    run.selectedColumn ∘ run.selectedIndexEquivFin.symm
  have finiteRows_injective : Function.Injective finiteRows :=
    run.selectedRow_injective.comp run.selectedIndexEquivFin.symm.injective
  have finiteColumns_injective : Function.Injective finiteColumns :=
    run.selectedColumn_injective.comp run.selectedIndexEquivFin.symm.injective
  have finiteRows_bound : ∀ i, (finiteRows i : ℝ) ≤ 1 := by
    intro i
    exact run.selectedRow_mem_of_completeOn run_complete _
  have finiteColumns_bound : ∀ j, (finiteColumns j : ℝ) ≤ 1 := by
    intro j
    exact run.selectedColumn_mem_of_completeOn run_complete _
  have eight_le_length : 8 ≤ run.pivots.length := by
    simpa [pivot_length] using eight_le
  have determinant_bound := positiveSeriesFermionicKernel_abs_minor_det_le
    finiteRows finiteColumns eight_le_length
    finiteRows_injective finiteColumns_injective
    finiteRows_bound finiteColumns_bound
  have finite_sample :
      minorMatrix positiveSeriesFermionicKernel finiteRows finiteColumns =
        run.finSelectedCore := by
    rfl
  rw [finite_sample] at determinant_bound
  exact residual_power.trans (by simpa [pivot_length] using determinant_bound)

/-- Ten exact complete pivots contract every positive-unit-band residual entry
by at least one half. -/
theorem positiveSeriesFermionicKernel_gecp_error_le_half_ten
    (residual : ℕ → Kernel PositiveArgument PositiveArgument ℝ)
    (rows columns : ℕ → PositiveArgument)
    (realized : ∀ m, ∃ run : Run positiveSeriesFermionicKernel,
      run.finalResidual = residual m)
    (pivot_ne : ∀ m, residual m (rows m) (columns m) ≠ 0)
    (updates : ∀ m, residual (m + 1) =
      residualUpdate (residual m) (rows m) (columns m) (pivot_ne m))
    (complete : ∀ m, CompletePivotOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun omega : PositiveArgument => (omega : ℝ) ≤ 1)
      (residual m) (rows m) (columns m))
    (run : Run positiveSeriesFermionicKernel)
    (run_pivots : run.RealizesPivotPrefix residual rows columns 10)
    (run_complete : run.CompleteOn
      (fun u : PositiveArgument => (u : ℝ) ≤ 1)
      (fun omega : PositiveArgument => (omega : ℝ) ≤ 1))
    (u omega : PositiveArgument) (hu : (u : ℝ) ≤ 1)
    (homega : (omega : ℝ) ≤ 1) :
    abs (residual 10 u omega) ≤ (1 / 2 : ℝ) := by
  have power_bound := positiveSeriesFermionicKernel_gecp_error_pow_le
    residual rows columns realized pivot_ne updates complete 10 (by norm_num)
    run run_pivots run_complete u omega hu homega
  apply le_of_pow_le_pow_left₀ (by norm_num : 10 ≠ 0) (by norm_num : (0 : ℝ) ≤ 1 / 2)
  exact power_bound.trans (by norm_num [Nat.factorial])

end Fermionic
end GECPKernelStructure
