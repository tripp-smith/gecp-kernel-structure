import GECPKernelStructure.GECP.ResidualTotalPositivity
import GECPKernelStructure.GECP.FrequencyReflection

namespace GECPKernelStructure.Fermionic

open GECP Matrix

/-- Reverse frequency order, not time: this presentation is defined on all
real coordinates and includes both signed frequency halves and time endpoints. -/
noncomputable def frequencyReflectedFermionicKernel : Kernel ℝ ℝ ℝ :=
  fun t y ↦ fermionicKernel t (-y)

theorem frequencyReflectedFermionicKernel_eq (t y : ℝ) :
    frequencyReflectedFermionicKernel t y = Real.exp (t * y) / (1 + Real.exp y) := by
  simp [frequencyReflectedFermionicKernel, fermionicKernel]

/-- The frequency-reversed fermionic kernel has positive ordered minors on
the full real coordinate domains, with no strict interior or cutoff premise. -/
theorem frequencyReflectedFermionicKernel_strictSignRegular :
    StrictSignRegular frequencyReflectedFermionicKernel (fun _ ↦ 1) := by
  have base : StrictSignRegular (fun t y : ℝ ↦ Real.exp (t * y)) (fun _ ↦ 1) := by
    intro n rows columns hr hc
    simp only [one_mul]
    exact expMatrix_det_pos_of_strictMono rows columns hr hc
  intro n
  have scaled := strictSignRegularAtOrder_columnScale (base n)
    (fun y : ℝ ↦ (1 + Real.exp y)⁻¹) (fun _ ↦ inv_pos.mpr (by positivity))
  have heq : frequencyReflectedFermionicKernel =
      fun t y : ℝ ↦ Real.exp (t * y) * (1 + Real.exp y)⁻¹ := by
    funext t y
    rw [frequencyReflectedFermionicKernel_eq, div_eq_mul_inv]
  rw [heq]
  exact scaled

/-- The actual fermionic family inherits the common all-orders gauge through
arbitrary successful prefixes, including boundary times and signed frequencies. -/
theorem frequencyReflectedFermionicKernel_gaugedFinalResidual_strictSignRegular
    (run : Run frequencyReflectedFermionicKernel) :
    StrictSignRegular run.gaugedFinalResidual (fun _ ↦ 1) :=
  Run.gaugedFinalResidual_strictSignRegular
    frequencyReflectedFermionicKernel_strictSignRegular run

/-- Any run of the original kernel transports to the positive-minor
presentation without reselecting pivots. -/
theorem fermionicKernel_reflectedRun_gaugedFinalResidual_strictSignRegular
    (run : Run fermionicKernel) :
    StrictSignRegular run.negateColumns.gaugedFinalResidual (fun _ ↦ 1) :=
  Run.gaugedFinalResidual_strictSignRegular
    frequencyReflectedFermionicKernel_strictSignRegular run.negateColumns

/-- The transported final residual is exactly the original final residual,
with the inverse frequency reflection. -/
theorem fermionicKernel_reflectedRun_finalResidual
    (run : Run fermionicKernel) (t omega : ℝ) :
    run.negateColumns.finalResidual t (-omega) = run.finalResidual t omega := by
  rw [Run.finalResidual_negateColumns]
  simp only [neg_neg]

end GECPKernelStructure.Fermionic
