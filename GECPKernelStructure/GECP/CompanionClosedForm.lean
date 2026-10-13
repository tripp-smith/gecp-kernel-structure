import GECPKernelStructure.GECP.CompanionPlucker
import GECPKernelStructure.GECP.Interpolation
import GECPKernelStructure.Matrix.DoubleBorderedCompanion

namespace GECPKernelStructure
namespace GECP

open Matrix

universe u v w

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

/-- A fixed left companion update commutes with an arbitrary fixed right
companion update. -/
theorem companionRowUpdate_companionColumnUpdate
    (K G : Kernel α β 𝕜) (H : Matrix α β 𝕜)
    (leftRow : α) (leftColumn : β)
    (left_ne : K leftRow leftColumn ≠ 0)
    (rightRow : α) (rightColumn : β)
    (right_ne : G rightRow rightColumn ≠ 0) :
    companionRowUpdate K
        (companionColumnUpdate G H rightRow rightColumn right_ne)
        leftRow leftColumn left_ne =
      companionColumnUpdate G
        (companionRowUpdate K H leftRow leftColumn left_ne)
        rightRow rightColumn right_ne := by
  ext x y
  simp only [companionRowUpdate, companionColumnUpdate]
  field_simp
  ring

namespace Run

variable {K : Kernel α β 𝕜}

/-- Original-kernel columns at the coordinates selected by a run. -/
def selectedKernelColumns (run : Run K) : Matrix α run.SelectedIndex 𝕜 :=
  fun x i => K x (run.selectedColumn i)

/-- Original-kernel rows at the coordinates selected by a run. -/
def selectedKernelRows (run : Run K) : Matrix run.SelectedIndex β 𝕜 :=
  fun i y => K (run.selectedRow i) y

/-- Companion rows at the row coordinates selected by a run. -/
def selectedCompanionRows {δ : Type*} (run : Run K)
    (H : Matrix α δ 𝕜) : Matrix run.SelectedIndex δ 𝕜 :=
  fun i y => H (run.selectedRow i) y

/-- Companion columns at the column coordinates selected by a run. -/
def selectedCompanionColumns {δ : Type*} (run : Run K)
    (H : Matrix δ β 𝕜) : Matrix δ run.SelectedIndex 𝕜 :=
  fun x i => H x (run.selectedColumn i)

/-- A companion sampled on both selected coordinate tuples. -/
def selectedCompanionCore (run : Run K) (H : Matrix α β 𝕜) :
    Matrix run.SelectedIndex run.SelectedIndex 𝕜 :=
  fun i j => H (run.selectedRow i) (run.selectedColumn j)

theorem selectedRow_mem_rows (run : Run K) (i : run.SelectedIndex) :
    run.selectedRow i ∈ run.rows := by
  induction run with
  | nil initial => exact Fin.elim0 i
  | step row column pivot_ne tail ih =>
      cases i with
      | inl i => simp [selectedRow, rows]
      | inr i =>
          change tail.selectedRow i ∈ row :: tail.rows
          exact List.mem_cons_of_mem row (ih i)

theorem selectedColumn_mem_columns (run : Run K) (i : run.SelectedIndex) :
    run.selectedColumn i ∈ run.columns := by
  induction run with
  | nil initial => exact Fin.elim0 i
  | step row column pivot_ne tail ih =>
      cases i with
      | inl i => simp [selectedColumn, columns]
      | inr i =>
          change tail.selectedColumn i ∈ column :: tail.columns
          exact List.mem_cons_of_mem column (ih i)

theorem companionRowUpdate_submatrix_right {δ ε : Type*}
    (K : Kernel α β 𝕜) (H : Matrix α δ 𝕜) (f : ε → δ)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionRowUpdate K (Matrix.submatrix H id f) row column pivot_ne =
      Matrix.submatrix (companionRowUpdate K H row column pivot_ne) id f := by
  rfl

theorem conditionRows_submatrix_right {δ ε : Type*}
    (run : Run K) (H : Matrix α δ 𝕜) (f : ε → δ) :
    run.conditionRows (Matrix.submatrix H id f) =
      Matrix.submatrix (run.conditionRows H) id f := by
  induction run generalizing H with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionRows_step, companionRowUpdate_submatrix_right, ih]

theorem companionColumnUpdate_submatrix_left {δ ε : Type*}
    (K : Kernel α β 𝕜) (H : Matrix δ β 𝕜) (f : ε → δ)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionColumnUpdate K (Matrix.submatrix H f id) row column pivot_ne =
      Matrix.submatrix (companionColumnUpdate K H row column pivot_ne) f id := by
  rfl

theorem conditionColumns_submatrix_left {δ ε : Type*}
    (run : Run K) (H : Matrix δ β 𝕜) (f : ε → δ) :
    run.conditionColumns (Matrix.submatrix H f id) =
      Matrix.submatrix (run.conditionColumns H) f id := by
  induction run generalizing H with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionColumns_step, companionColumnUpdate_submatrix_left, ih]

@[simp]
theorem companionColumnUpdate_self (K : Kernel α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionColumnUpdate K K row column pivot_ne =
      residualUpdate K row column pivot_ne :=
  rfl

/-- Conditioning the original kernel from the right gives the same final
residual as ordinary row elimination. -/
theorem conditionColumns_self (run : Run K) :
    run.conditionColumns K = run.finalResidual := by
  induction run with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      change tail.conditionColumns
          (residualUpdate _ row column pivot_ne) = tail.finalResidual
      exact ih

theorem conditionColumns_add {δ : Type*} (run : Run K)
    (H G : Matrix δ β 𝕜) :
    run.conditionColumns (H + G) =
      run.conditionColumns H + run.conditionColumns G := by
  induction run generalizing H G with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionColumns_step, companionColumnUpdate_add, ih]

/-- Fixed column elimination commutes with a fixed left feature factor. -/
theorem conditionColumns_mul {δ : Type*} {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (run : Run K) (L : Matrix δ κ 𝕜) (R : Matrix κ β 𝕜) :
    run.conditionColumns (L * R) = L * run.conditionColumns R := by
  induction run generalizing R with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionColumns_step, companionColumnUpdate_mul, ih]

theorem conditionRows_selectedKernelColumns (run : Run K) :
    run.conditionRows run.selectedKernelColumns = 0 := by
  calc
    run.conditionRows run.selectedKernelColumns =
        Matrix.submatrix
          (run.conditionRows (fun x y => K x y : Matrix α β 𝕜))
          id run.selectedColumn :=
      conditionRows_submatrix_right run
        (fun x y => K x y : Matrix α β 𝕜) run.selectedColumn
    _ = Matrix.submatrix run.finalResidual id run.selectedColumn := by
      rw [conditionRows_self]
    _ = 0 := by
      ext x i
      exact run.gecp_interpolates_selected_cols
        (run.selectedColumn_mem_columns i) x

theorem conditionColumns_selectedKernelRows (run : Run K) :
    run.conditionColumns run.selectedKernelRows = 0 := by
  calc
    run.conditionColumns run.selectedKernelRows =
        Matrix.submatrix
          (run.conditionColumns (fun x y => K x y : Matrix α β 𝕜))
          run.selectedRow id :=
      conditionColumns_submatrix_left run
        (fun x y => K x y : Matrix α β 𝕜) run.selectedRow
    _ = Matrix.submatrix run.finalResidual run.selectedRow id := by
      rw [conditionColumns_self]
    _ = 0 := by
      ext i y
      exact run.gecp_interpolates_selected_rows
        (run.selectedRow_mem_rows i) y

theorem conditionRows_eq_self_of_selectedRows_zero {δ : Type*}
    (run : Run K) (H : Matrix α δ 𝕜)
    (zero_on_selected : ∀ i y, H (run.selectedRow i) y = 0) :
    run.conditionRows H = H := by
  induction run generalizing H with
  | nil => rfl
  | @step K' row column pivot_ne tail ih =>
      have head_zero (y : δ) : H row y = 0 :=
        zero_on_selected (Sum.inl 0) y
      have update_eq : companionRowUpdate K' H row column pivot_ne = H := by
        ext x y
        simp [companionRowUpdate, head_zero]
      rw [conditionRows_step, update_eq]
      apply ih
      intro i y
      exact zero_on_selected (Sum.inr i) y

theorem conditionColumns_eq_self_of_selectedColumns_zero {δ : Type*}
    (run : Run K) (H : Matrix δ β 𝕜)
    (zero_on_selected : ∀ x i, H x (run.selectedColumn i) = 0) :
    run.conditionColumns H = H := by
  induction run generalizing H with
  | nil => rfl
  | @step K' row column pivot_ne tail ih =>
      have head_zero (x : δ) : H x column = 0 :=
        zero_on_selected x (Sum.inl 0)
      have update_eq : companionColumnUpdate K' H row column pivot_ne = H := by
        ext x y
        simp [companionColumnUpdate, head_zero]
      rw [conditionColumns_step, update_eq]
      apply ih
      intro x i
      exact zero_on_selected x (Sum.inr i)

set_option linter.style.haveILetI false in
/-- Closed selected-core solve for every fixed-run row conditioning. -/
theorem conditionRows_eq_selectedCoreSchur {δ : Type*}
    (run : Run K) (H : Matrix α δ 𝕜) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    run.conditionRows H =
      H - run.selectedKernelColumns * ⅟run.selectedCore *
        run.selectedCompanionRows H := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  let correction := run.selectedKernelColumns * ⅟run.selectedCore *
    run.selectedCompanionRows H
  let remainder := H - correction
  have remainder_zero (i : run.SelectedIndex) (y : δ) :
      remainder (run.selectedRow i) y = 0 := by
    have correction_entry :
        correction (run.selectedRow i) y = H (run.selectedRow i) y := by
      change ((run.selectedCore * ⅟run.selectedCore) *
        run.selectedCompanionRows H) i y = _
      rw [mul_invOf_self, Matrix.one_mul]
      rfl
    simp [remainder, correction_entry]
  have conditioned_remainder : run.conditionRows remainder = remainder :=
    run.conditionRows_eq_self_of_selectedRows_zero remainder remainder_zero
  have decomposition : H = remainder + correction := by
    simp [remainder]
  rw [decomposition, conditionRows_add, conditioned_remainder]
  rw [show correction = run.selectedKernelColumns *
      (⅟run.selectedCore * run.selectedCompanionRows H) by
        simp [correction, Matrix.mul_assoc]]
  rw [conditionRows_mul, conditionRows_selectedKernelColumns]
  simp [remainder, correction, Matrix.mul_assoc]

set_option linter.style.haveILetI false in
/-- Closed selected-core solve for every fixed-run column conditioning. -/
theorem conditionColumns_eq_selectedCoreSchur {δ : Type*}
    (run : Run K) (H : Matrix δ β 𝕜) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    run.conditionColumns H =
      H - run.selectedCompanionColumns H * ⅟run.selectedCore *
        run.selectedKernelRows := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  let correction := run.selectedCompanionColumns H * ⅟run.selectedCore *
    run.selectedKernelRows
  let remainder := H - correction
  have remainder_zero (x : δ) (i : run.SelectedIndex) :
      remainder x (run.selectedColumn i) = 0 := by
    have correction_entry :
        correction x (run.selectedColumn i) = H x (run.selectedColumn i) := by
      change ((run.selectedCompanionColumns H * ⅟run.selectedCore) *
        run.selectedCore) x i = _
      rw [Matrix.mul_assoc, invOf_mul_self, Matrix.mul_one]
      rfl
    simp [remainder, correction_entry]
  have conditioned_remainder : run.conditionColumns remainder = remainder :=
    run.conditionColumns_eq_self_of_selectedColumns_zero remainder remainder_zero
  have decomposition : H = remainder + correction := by
    simp [remainder]
  rw [decomposition, conditionColumns_add, conditioned_remainder]
  rw [show correction =
      (run.selectedCompanionColumns H * ⅟run.selectedCore) *
        run.selectedKernelRows by simp [correction]]
  rw [conditionColumns_mul, conditionColumns_selectedKernelRows]
  simp [remainder, correction, Matrix.mul_assoc]

/-- All row projectors of a run commute past one fixed right projector. -/
theorem conditionRows_companionColumnUpdate
    (run : Run K) (G : Kernel α β 𝕜) (H : Matrix α β 𝕜)
    (row : α) (column : β) (pivot_ne : G row column ≠ 0) :
    run.conditionRows
        (companionColumnUpdate G H row column pivot_ne) =
      companionColumnUpdate G (run.conditionRows H) row column pivot_ne := by
  induction run generalizing H with
  | nil => rfl
  | step selectedRow selectedColumn selected_ne tail ih =>
      simp only [conditionRows_step]
      rw [companionRowUpdate_companionColumnUpdate, ih]

/-- Iterated two-sided conditioning is the composition of all left projectors
followed by all right projectors. -/
theorem conditionBoth_eq_conditionColumns_conditionRows
    (run : Run K) (H : Matrix α β 𝕜) :
    run.conditionBoth H = run.conditionColumns (run.conditionRows H) := by
  induction run generalizing H with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionBoth_step, conditionRows_step, conditionColumns_step]
      rw [ih]
      simp only [companionBiUpdate]
      rw [conditionRows_companionColumnUpdate]

set_option linter.style.haveILetI false in
/-- Sampling the selected columns after row conditioning gives the matching
selected-core correction. -/
theorem selectedCompanionColumns_conditionRows (run : Run K)
    (H : Matrix α β 𝕜) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    run.selectedCompanionColumns (run.conditionRows H) =
      run.selectedCompanionColumns H -
        run.selectedKernelColumns * ⅟run.selectedCore *
          run.selectedCompanionCore H := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  rw [run.conditionRows_eq_selectedCoreSchur H]
  rfl

set_option linter.style.haveILetI false in
/-- Every recursively conditioned two-sided companion is exactly the global
four-term Schur projection through the original selected core. -/
theorem conditionBoth_eq_twoSidedCompanionSchur (run : Run K)
    (H : Matrix α β 𝕜) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    run.conditionBoth H =
      Matrix.twoSidedCompanionSchur run.selectedCore
        run.selectedKernelColumns run.selectedKernelRows
        (run.selectedCompanionCore H)
        (run.selectedCompanionColumns H)
        (run.selectedCompanionRows H) H := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  rw [conditionBoth_eq_conditionColumns_conditionRows,
    conditionColumns_eq_selectedCoreSchur,
    selectedCompanionColumns_conditionRows,
    conditionRows_eq_selectedCoreSchur]
  simp only [Matrix.twoSidedCompanionSchur, Matrix.sub_mul, Matrix.mul_assoc]
  noncomm_ring

set_option linter.style.haveILetI false in
/-- A finite sample of recursive two-sided conditioning is the corresponding
sampled selected-core Schur expression. -/
theorem conditionBoth_sample_eq_twoSidedCompanionSchur
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    (run : Run K) (H : Matrix α β 𝕜)
    (rows : γ → α) (columns : γ → β) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    Matrix.submatrix (run.conditionBoth H) rows columns =
      Matrix.twoSidedCompanionSchur run.selectedCore
        (Matrix.submatrix run.selectedKernelColumns rows id)
        (Matrix.submatrix run.selectedKernelRows id columns)
        (run.selectedCompanionCore H)
        (Matrix.submatrix (run.selectedCompanionColumns H) rows id)
        (Matrix.submatrix (run.selectedCompanionRows H) id columns)
        (Matrix.submatrix H rows columns) := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  rw [run.conditionBoth_eq_twoSidedCompanionSchur H]
  ext i j
  simp only [Matrix.twoSidedCompanionSchur, Matrix.submatrix_apply,
    Matrix.sub_apply, Matrix.add_apply, Matrix.mul_apply, id_eq]

set_option linter.style.haveILetI false in
/-- The Phase BX double border specialized to any finite sample of a
successful real GECP run. Its determinant is the squared selected-core
determinant times the sampled recursively conditioned determinant. -/
theorem doubleBorderedCompanionMatrix_det_eq_conditionBoth_sample
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    {Kℝ : Kernel α β ℝ} (run : Run Kℝ) (H : Matrix α β ℝ)
    (rows : γ → α) (columns : γ → β) :
    letI : Invertible run.selectedCore :=
      Matrix.invertibleOfIsUnitDet
        run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
    (Matrix.doubleBorderedCompanionMatrix run.selectedCore
        (Matrix.submatrix run.selectedKernelColumns rows id)
        (Matrix.submatrix run.selectedKernelRows id columns)
        (run.selectedCompanionCore H)
        (Matrix.submatrix (run.selectedCompanionColumns H) rows id)
        (Matrix.submatrix (run.selectedCompanionRows H) id columns)
        (Matrix.submatrix H rows columns)).det =
      run.selectedCore.det ^ 2 *
        (Matrix.submatrix (run.conditionBoth H) rows columns).det := by
  letI : Invertible run.selectedCore :=
    Matrix.invertibleOfIsUnitDet
      run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  rw [Matrix.doubleBorderedCompanionMatrix_det]
  rw [← run.conditionBoth_sample_eq_twoSidedCompanionSchur H rows columns]

end Run
end GECP
end GECPKernelStructure
