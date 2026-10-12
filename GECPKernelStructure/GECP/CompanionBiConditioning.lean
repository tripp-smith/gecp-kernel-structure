import GECPKernelStructure.GECP.CompanionConditioning
import Mathlib.Tactic

namespace GECPKernelStructure
namespace GECP

universe u v w z

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

/-- Apply the right projector of one fixed GECP pivot to a companion matrix. -/
def companionColumnUpdate {δ : Type z} (K : Kernel α β 𝕜)
    (H : Matrix δ β 𝕜) (row : α) (column : β)
    (_pivot_ne : K row column ≠ 0) : Matrix δ β 𝕜 :=
  fun x y => H x y - H x column * K row y / K row column

/-- Apply both fixed oblique pivot projectors to a companion matrix.  Unlike
ordinary Schur elimination, the pivot data always come from `K`; the map is
therefore linear in `H`. -/
def companionBiUpdate (K : Kernel α β 𝕜) (H : Matrix α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    Matrix α β 𝕜 :=
  companionColumnUpdate K
    (companionRowUpdate K H row column pivot_ne)
    row column pivot_ne

theorem companionBiUpdate_apply (K : Kernel α β 𝕜) (H : Matrix α β 𝕜)
    (row x : α) (column y : β) (pivot_ne : K row column ≠ 0) :
    companionBiUpdate K H row column pivot_ne x y =
      H x y - K x column * H row y / K row column -
        H x column * K row y / K row column +
          K x column * H row column * K row y /
            (K row column) ^ 2 := by
  simp only [companionBiUpdate, companionColumnUpdate, companionRowUpdate]
  field_simp
  ring

@[simp]
theorem companionColumnUpdate_selectedColumn {δ : Type z}
    (K : Kernel α β 𝕜) (H : Matrix δ β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) (x : δ) :
    companionColumnUpdate K H row column pivot_ne x column = 0 := by
  simp only [companionColumnUpdate]
  field_simp
  ring

@[simp]
theorem companionBiUpdate_selectedRow (K : Kernel α β 𝕜)
    (H : Matrix α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) (y : β) :
    companionBiUpdate K H row column pivot_ne row y = 0 := by
  rw [companionBiUpdate_apply]
  field_simp
  ring

@[simp]
theorem companionBiUpdate_selectedColumn (K : Kernel α β 𝕜)
    (H : Matrix α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) (x : α) :
    companionBiUpdate K H row column pivot_ne x column = 0 :=
  companionColumnUpdate_selectedColumn K
    (companionRowUpdate K H row column pivot_ne)
    row column pivot_ne x

@[simp]
theorem companionBiUpdate_self (K : Kernel α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
  companionBiUpdate K K row column pivot_ne =
      residualUpdate K row column pivot_ne := by
  ext x y
  simp only [companionBiUpdate, companionColumnUpdate,
    companionRowUpdate, residualUpdate]
  field_simp
  ring

theorem companionColumnUpdate_add {δ : Type z} (K : Kernel α β 𝕜)
    (H G : Matrix δ β 𝕜) (row : α) (column : β)
    (pivot_ne : K row column ≠ 0) :
    companionColumnUpdate K (H + G) row column pivot_ne =
      companionColumnUpdate K H row column pivot_ne +
        companionColumnUpdate K G row column pivot_ne := by
  ext x y
  simp only [companionColumnUpdate, Matrix.add_apply]
  ring

theorem companionColumnUpdate_smul {δ : Type z} (K : Kernel α β 𝕜)
    (c : 𝕜) (H : Matrix δ β 𝕜) (row : α) (column : β)
    (pivot_ne : K row column ≠ 0) :
    companionColumnUpdate K (c • H) row column pivot_ne =
      c • companionColumnUpdate K H row column pivot_ne := by
  ext x y
  simp only [companionColumnUpdate, Matrix.smul_apply, smul_eq_mul]
  ring

theorem companionColumnUpdate_mul {δ : Type z} {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (K : Kernel α β 𝕜) (L : Matrix δ κ 𝕜) (R : Matrix κ β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionColumnUpdate K (L * R) row column pivot_ne =
      L * companionColumnUpdate K R row column pivot_ne := by
  ext x y
  simp only [companionColumnUpdate, Matrix.mul_apply, div_eq_mul_inv]
  rw [mul_assoc, Finset.sum_mul]
  simp_rw [mul_sub]
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro k _
  ring

theorem companionBiUpdate_add (K : Kernel α β 𝕜)
    (H G : Matrix α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionBiUpdate K (H + G) row column pivot_ne =
      companionBiUpdate K H row column pivot_ne +
        companionBiUpdate K G row column pivot_ne := by
  rw [companionBiUpdate, companionRowUpdate_add,
    companionColumnUpdate_add]
  rfl

theorem companionBiUpdate_smul (K : Kernel α β 𝕜) (c : 𝕜)
    (H : Matrix α β 𝕜) (row : α) (column : β)
    (pivot_ne : K row column ≠ 0) :
    companionBiUpdate K (c • H) row column pivot_ne =
      c • companionBiUpdate K H row column pivot_ne := by
  rw [companionBiUpdate, companionRowUpdate_smul,
    companionColumnUpdate_smul]
  rfl

/-- Two-sided fixed-pivot projection transports a low-rank product by
conditioning its left and right feature factors separately. -/
theorem companionBiUpdate_mul {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (K : Kernel α β 𝕜) (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionBiUpdate K (L * R) row column pivot_ne =
      companionRowUpdate K L row column pivot_ne *
        companionColumnUpdate K R row column pivot_ne := by
  change companionColumnUpdate K
      (companionRowUpdate K (L * R) row column pivot_ne)
      row column pivot_ne = _
  rw [companionRowUpdate_mul, companionColumnUpdate_mul]

namespace Run

variable {K : Kernel α β 𝕜}

/-- Transport a right feature factor through all right projectors of a fixed
GECP run. -/
noncomputable def conditionColumns {δ : Type z} (run : Run K)
    (R : Matrix δ β 𝕜) : Matrix δ β 𝕜 :=
  Run.rec (motive := fun _ _ => Matrix δ β 𝕜 → Matrix δ β 𝕜)
    (fun _ R => R)
    (fun {K} row column pivot_ne _tail recurse R =>
      recurse (companionColumnUpdate K R row column pivot_ne))
    run R

/-- Transport a companion kernel through both projectors of every fixed pivot
in a GECP run. -/
noncomputable def conditionBoth (run : Run K)
    (H : Matrix α β 𝕜) : Matrix α β 𝕜 :=
  Run.rec (motive := fun _ _ => Matrix α β 𝕜 → Matrix α β 𝕜)
    (fun _ H => H)
    (fun {K} row column pivot_ne _tail recurse H =>
      recurse (companionBiUpdate K H row column pivot_ne))
    run H

@[simp]
theorem conditionColumns_nil {δ : Type z} (K : Kernel α β 𝕜)
    (R : Matrix δ β 𝕜) :
    (Run.nil K).conditionColumns R = R :=
  rfl

@[simp]
theorem conditionBoth_nil (K : Kernel α β 𝕜) (H : Matrix α β 𝕜) :
    (Run.nil K).conditionBoth H = H :=
  rfl

@[simp]
theorem conditionColumns_step {δ : Type z} {K : Kernel α β 𝕜}
    (row : α) (column : β) (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (R : Matrix δ β 𝕜) :
    (Run.step row column pivot_ne tail).conditionColumns R =
      tail.conditionColumns
        (companionColumnUpdate K R row column pivot_ne) :=
  rfl

@[simp]
theorem conditionBoth_step {K : Kernel α β 𝕜}
    (row : α) (column : β) (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (H : Matrix α β 𝕜) :
    (Run.step row column pivot_ne tail).conditionBoth H =
      tail.conditionBoth (companionBiUpdate K H row column pivot_ne) :=
  rfl

theorem conditionBoth_self (run : Run K) :
    run.conditionBoth K = run.finalResidual := by
  induction run with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      change tail.conditionBoth
        (companionBiUpdate _ _ row column pivot_ne) = tail.finalResidual
      rw [companionBiUpdate_self]
      exact ih

theorem conditionBoth_add (run : Run K) (H G : Matrix α β 𝕜) :
    run.conditionBoth (H + G) =
      run.conditionBoth H + run.conditionBoth G := by
  induction run generalizing H G with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionBoth_step, companionBiUpdate_add, ih]

theorem conditionBoth_smul (run : Run K) (c : 𝕜)
    (H : Matrix α β 𝕜) :
    run.conditionBoth (c • H) = c • run.conditionBoth H := by
  induction run generalizing H with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionBoth_step, companionBiUpdate_smul, ih]

/-- Full two-sided conditioning preserves a fixed feature dimension while
transporting both feature factors. -/
theorem conditionBoth_mul {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (run : Run K) (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜) :
    run.conditionBoth (L * R) =
      run.conditionRows L * run.conditionColumns R := by
  induction run generalizing L R with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      change tail.conditionBoth
        (companionBiUpdate _ (L * R) row column pivot_ne) = _
      rw [companionBiUpdate_mul, ih]
      rfl

theorem finalResidual_eq_conditionBoth_add_mul {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (run : Run K) (E : Matrix α β 𝕜)
    (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜)
    (split : (K : Matrix α β 𝕜) = E + L * R) :
    (run.finalResidual : Matrix α β 𝕜) =
      run.conditionBoth E +
      run.conditionRows L * run.conditionColumns R := by
  subst K
  rw [← run.conditionBoth_self, run.conditionBoth_add,
    run.conditionBoth_mul]

theorem finalResidual_sample_eq_conditionBoth_add_mul
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    {κ : Type*} [Fintype κ] [DecidableEq κ]
    (run : Run K) (rows : γ → α) (columns : γ → β)
    (E : Matrix α β 𝕜) (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜)
    (split : (K : Matrix α β 𝕜) = E + L * R) :
    Matrix.submatrix run.finalResidual rows columns =
      Matrix.submatrix (run.conditionBoth E) rows columns +
        Matrix.submatrix (run.conditionRows L) rows id *
          Matrix.submatrix (run.conditionColumns R) id columns := by
  have conditioned := run.finalResidual_eq_conditionBoth_add_mul E L R split
  ext i j
  have entry := congrFun (congrFun conditioned (rows i)) (columns j)
  change run.finalResidual (rows i) (columns j) =
    run.conditionBoth E (rows i) (columns j) +
      ∑ k, run.conditionRows L (rows i) k *
        run.conditionColumns R k (columns j)
  simpa only [Matrix.submatrix_apply, id_eq, Matrix.add_apply,
    Matrix.mul_apply] using entry

end Run
end GECP
end GECPKernelStructure
