import GECPKernelStructure.GECP.Definitions
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

universe u v

/-- The only classes of rigorous outcome allowed to close structural Milestone E. -/
inductive GECPResearchOutcome where
  | targetRate
  | weakerImprovedRate
  | containingStructuralClass
  | certifiedObstruction
  deriving DecidableEq, Repr

/--
A selected cross controls the next residual relative to a uniform bound `B`.
This is the quantitative refinement that survives the pivot-order obstruction.
-/
def CrossRatioControl {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (row : α) (column : β) (θ B : ℝ) : Prop :=
  ∀ x y, abs (R x y * R row column - R x column * R row y) ≤
    θ * abs (R row column) * B

/-- Cross-ratio control gives the corresponding one-step residual contraction. -/
theorem residualUpdate_le_of_crossRatioControl {α : Type u} {β : Type v}
    {R : Kernel α β ℝ} {row : α} {column : β} (pivot_ne : R row column ≠ 0)
    {θ B : ℝ} (control : CrossRatioControl R row column θ B) (x : α) (y : β) :
    abs (residualUpdate R row column pivot_ne x y) ≤ θ * B := by
  unfold residualUpdate
  rw [show R x y - R x column * R row y / R row column =
      (R x y * R row column - R x column * R row y) / R row column by
    field_simp]
  rw [abs_div]
  apply (div_le_div_of_nonneg_right (control x y) (abs_nonneg _)).trans_eq
  field_simp [abs_ne_zero.mpr pivot_ne]

/--
A quantitative cross-ratio hypothesis at every selected pivot is sufficient
for geometric GECP residual decay. The theorem deliberately does not assert
that the fermionic kernel satisfies the hypothesis.
-/
theorem gecp_error_le_geometric_of_crossRatioControl
    {α : Type u} {β : Type v} (residual : ℕ → Kernel α β ℝ)
    (rows : ℕ → α) (columns : ℕ → β) (θ initialBound : ℝ)
    (theta_nonneg : 0 ≤ θ) (initial_nonneg : 0 ≤ initialBound)
    (initial_bound : ∀ x y, abs (residual 0 x y) ≤ initialBound)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (cross_ratio : ∀ n B, 0 ≤ B →
      (∀ x y, abs (residual n x y) ≤ B) →
      CrossRatioControl (residual n) (rows n) (columns n) θ B) :
    ∀ n x y, abs (residual n x y) ≤ θ ^ n * initialBound := by
  intro n
  induction n with
  | zero => simpa using initial_bound
  | succ n ih =>
      intro x y
      rw [show n + 1 = Nat.succ n by omega, updates n]
      have bound_nonneg : 0 ≤ θ ^ n * initialBound :=
        mul_nonneg (pow_nonneg theta_nonneg n) initial_nonneg
      have step := residualUpdate_le_of_crossRatioControl (pivot_ne n)
        (cross_ratio n (θ ^ n * initialBound) bound_nonneg ih) x y
      calc
        abs (residualUpdate (residual n) (rows n) (columns n) (pivot_ne n) x y) ≤
            θ * (θ ^ n * initialBound) := step
        _ = θ ^ (n + 1) * initialBound := by rw [pow_succ]; ring

end Fermionic
end GECPKernelStructure
