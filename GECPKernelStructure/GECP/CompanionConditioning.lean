import GECPKernelStructure.GECP.Residual
import Mathlib.Data.Matrix.Mul

namespace GECPKernelStructure
namespace GECP

universe u v w z

variable {α : Type u} {β : Type v} {𝕜 : Type w} [Field 𝕜]

/-- Apply one fixed GECP row-elimination operation to a companion matrix.  The
elimination coefficient comes from `K`, while the selected row comes from
`H`; thus the operation is linear in `H`. -/
def companionRowUpdate {δ : Type z} (K : Kernel α β 𝕜)
    (H : Matrix α δ 𝕜) (row : α) (column : β)
    (_pivot_ne : K row column ≠ 0) : Matrix α δ 𝕜 :=
  fun x y => H x y - K x column * H row y / K row column

@[simp]
theorem companionRowUpdate_self (K : Kernel α β 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionRowUpdate K K row column pivot_ne =
      residualUpdate K row column pivot_ne :=
  rfl

theorem companionRowUpdate_add {δ : Type z} (K : Kernel α β 𝕜)
    (H G : Matrix α δ 𝕜) (row : α) (column : β)
    (pivot_ne : K row column ≠ 0) :
    companionRowUpdate K (H + G) row column pivot_ne =
      companionRowUpdate K H row column pivot_ne +
        companionRowUpdate K G row column pivot_ne := by
  ext x y
  simp only [companionRowUpdate, Matrix.add_apply]
  ring

theorem companionRowUpdate_smul {δ : Type z} (K : Kernel α β 𝕜)
    (c : 𝕜) (H : Matrix α δ 𝕜) (row : α) (column : β)
    (pivot_ne : K row column ≠ 0) :
    companionRowUpdate K (c • H) row column pivot_ne =
      c • companionRowUpdate K H row column pivot_ne := by
  ext x y
  simp only [companionRowUpdate, Matrix.smul_apply, smul_eq_mul]
  ring

theorem companionRowUpdate_mul {δ : Type z} {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (K : Kernel α β 𝕜) (L : Matrix α κ 𝕜) (R : Matrix κ δ 𝕜)
    (row : α) (column : β) (pivot_ne : K row column ≠ 0) :
    companionRowUpdate K (L * R) row column pivot_ne =
      companionRowUpdate K L row column pivot_ne * R := by
  ext x y
  simp only [companionRowUpdate, Matrix.mul_apply, div_eq_mul_inv]
  rw [Finset.mul_sum, Finset.sum_mul]
  simp_rw [sub_mul]
  rw [Finset.sum_sub_distrib]
  congr 1
  apply Finset.sum_congr rfl
  intro k _
  ring

namespace Run

variable {K : Kernel α β 𝕜}

/-- Transport any companion matrix through the fixed row operations of a GECP
run.  No pivots are recomputed from the companion. -/
noncomputable def conditionRows {δ : Type z} (run : Run K)
    (H : Matrix α δ 𝕜) : Matrix α δ 𝕜 :=
  Run.rec (motive := fun _ _ => Matrix α δ 𝕜 → Matrix α δ 𝕜)
    (fun _ H => H)
    (fun {K} row column pivot_ne _tail recurse H =>
      recurse (companionRowUpdate K H row column pivot_ne))
    run H

@[simp]
theorem conditionRows_nil {δ : Type z} (K : Kernel α β 𝕜)
    (H : Matrix α δ 𝕜) :
    (Run.nil K).conditionRows H = H :=
  rfl

@[simp]
theorem conditionRows_step {δ : Type z} {K : Kernel α β 𝕜}
    (row : α) (column : β) (pivot_ne : K row column ≠ 0)
    (tail : Run (residualUpdate K row column pivot_ne))
    (H : Matrix α δ 𝕜) :
    (Run.step row column pivot_ne tail).conditionRows H =
      tail.conditionRows (companionRowUpdate K H row column pivot_ne) :=
  rfl

/-- Conditioning the original kernel itself recovers the run's final
residual. -/
theorem conditionRows_self (run : Run K) :
    run.conditionRows K = run.finalResidual := by
  induction run with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      change tail.conditionRows (residualUpdate _ row column pivot_ne) =
        tail.finalResidual
      exact ih

theorem conditionRows_add {δ : Type z} (run : Run K)
    (H G : Matrix α δ 𝕜) :
    run.conditionRows (H + G) =
      run.conditionRows H + run.conditionRows G := by
  induction run generalizing H G with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionRows_step, companionRowUpdate_add, ih]

theorem conditionRows_smul {δ : Type z} (run : Run K)
    (c : 𝕜) (H : Matrix α δ 𝕜) :
    run.conditionRows (c • H) = c • run.conditionRows H := by
  induction run generalizing H with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionRows_step, companionRowUpdate_smul, ih]

/-- Fixed row elimination commutes with a fixed right feature factor. -/
theorem conditionRows_mul {δ : Type z} {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (run : Run K) (L : Matrix α κ 𝕜) (R : Matrix κ δ 𝕜) :
    run.conditionRows (L * R) = run.conditionRows L * R := by
  induction run generalizing L with
  | nil => rfl
  | step row column pivot_ne tail ih =>
      simp only [conditionRows_step, companionRowUpdate_mul, ih]

/-- An initial low-rank-plus-tail split survives any exact GECP prefix with
the same right feature dimension. -/
theorem finalResidual_eq_conditionRows_add_mul {κ : Type*}
    [Fintype κ] [DecidableEq κ]
    (run : Run K) (E : Matrix α β 𝕜)
    (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜)
    (split : (K : Matrix α β 𝕜) = E + L * R) :
    (run.finalResidual : Matrix α β 𝕜) =
      run.conditionRows E + run.conditionRows L * R := by
  subst K
  rw [← run.conditionRows_self, run.conditionRows_add,
    run.conditionRows_mul]

/-- Sampling the canonical conditioned split yields exactly the finite matrix
required by residual determinant amplification. -/
theorem finalResidual_sample_eq_conditionRows_add_mul
    {γ : Type*} [Fintype γ] [DecidableEq γ]
    {κ : Type*} [Fintype κ] [DecidableEq κ]
    (run : Run K) (rows : γ → α) (columns : γ → β)
    (E : Matrix α β 𝕜) (L : Matrix α κ 𝕜) (R : Matrix κ β 𝕜)
    (split : (K : Matrix α β 𝕜) = E + L * R) :
    Matrix.submatrix run.finalResidual rows columns =
      Matrix.submatrix (run.conditionRows E) rows columns +
        Matrix.submatrix (run.conditionRows L) rows id *
          Matrix.submatrix R id columns := by
  have conditioned := run.finalResidual_eq_conditionRows_add_mul E L R split
  ext i j
  have entry := congrFun (congrFun conditioned (rows i)) (columns j)
  change run.finalResidual (rows i) (columns j) =
    run.conditionRows E (rows i) (columns j) +
      ∑ k, run.conditionRows L (rows i) k * R k (columns j)
  simpa only [Matrix.submatrix_apply, id_eq, Matrix.add_apply,
    Matrix.mul_apply] using entry

end Run
end GECP
end GECPKernelStructure
