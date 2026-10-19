import GECPKernelStructure.GECP.Definitions
import Mathlib.Data.Real.Basic

namespace GECPKernelStructure.GECP.Run

private def castKernel {K L : Kernel ℝ ℝ ℝ} (h : K = L) (run : Run K) : Run L :=
  h ▸ run

private theorem finalResidual_castKernel {K L : Kernel ℝ ℝ ℝ} (h : K = L) (run : Run K) :
    (castKernel h run).finalResidual = run.finalResidual := by
  subst h
  rfl

private theorem finalResidual_step_eq {K : Kernel ℝ ℝ ℝ}
    (p q : ℝ) (hne : K p q ≠ 0) (tail : Run (residualUpdate K p q hne)) :
    (Run.step p q hne tail).finalResidual = tail.finalResidual := rfl

/-- Transport an exact run under the frequency-order reversing bijection.
The same time nodes and negated frequency nodes produce the reflected run. -/
noncomputable def negateColumns {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    Run (fun x y ↦ K x (-y)) := by
  induction run with
  | nil initial => exact Run.nil _
  | @step K p q hne tail ih =>
    have hp : (fun x y ↦ K x (-y)) p (-q) ≠ 0 := by simpa only [neg_neg] using hne
    have hh : (fun x y ↦ residualUpdate K p q hne x (-y)) =
        residualUpdate (fun x y ↦ K x (-y)) p (-q) hp := by
      funext x y
      simp only [residualUpdate, neg_neg]
    exact Run.step p (-q) hp (castKernel hh ih)

theorem finalResidual_negateColumns {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    run.negateColumns.finalResidual = fun x y ↦ run.finalResidual x (-y) := by
  induction run with
  | nil initial => rfl
  | step p q hne tail ih =>
    simpa only [negateColumns, finalResidual_step_eq, finalResidual_castKernel] using ih

end GECPKernelStructure.GECP.Run
