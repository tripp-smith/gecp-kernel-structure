import GECPKernelStructure.GECP.FiniteLocalMinorWitness
import GECPKernelStructure.GreenFunction

namespace GECPKernelStructure.GECP

/-- A zero-valued complete maximum certifies exact domain termination. -/
theorem residual_eq_zero_of_complete_zero {R : Kernel ℝ ℝ ℝ} {D C : ℝ → Prop}
    {p q : ℝ} (complete : CompletePivotOn D C R p q) (hzero : R p q = 0)
    (x y : ℝ) (hx : D x) (hy : C y) : R x y = 0 := by
  have hh := complete.2.2 x y hx hy
  rw [hzero, abs_zero] at hh
  exact abs_eq_zero.mp (le_antisymm hh (abs_nonneg _))

namespace Run

variable {K : Kernel ℝ ℝ ℝ}

/-- The algebraic selected-core cross. Numerical evaluation uses solves. -/
noncomputable def crossApproximation (run : Run K) : Kernel ℝ ℝ ℝ :=
  letI : Invertible run.selectedCore := Matrix.invertibleOfIsUnitDet
    run.selectedCore (isUnit_iff_ne_zero.mpr run.selectedCore_det_ne)
  run.selectedKernelColumns * ⅟run.selectedCore * run.selectedKernelRows

/-- The iterative finite residual is exactly the selected-core cross error. -/
theorem finalResidual_eq_sub_crossApproximation (run : Run K) :
    run.finalResidual = fun x y ↦ K x y - run.crossApproximation x y := by
  have hh := run.conditionColumns_eq_selectedCoreSchur K
  rw [run.conditionColumns_self] at hh
  exact hh

theorem kernel_sub_crossApproximation_eq_finalResidual (run : Run K) (x y : ℝ) :
    K x y - run.crossApproximation x y = run.finalResidual x y := by
  exact (congrFun (congrFun run.finalResidual_eq_sub_crossApproximation x) y).symm

theorem crossApproximation_eq_sub_finalResidual (run : Run K) :
    run.crossApproximation = fun x y ↦ K x y - run.finalResidual x y := by
  funext x y
  have hh := run.kernel_sub_crossApproximation_eq_finalResidual x y
  linarith

theorem crossApproximation_eq_kernel_of_complete_zero (run : Run K) {D C : ℝ → Prop}
    {p q : ℝ} (complete : CompletePivotOn D C run.finalResidual p q)
    (hzero : run.finalResidual p q = 0) (x y : ℝ) (hx : D x) (hy : C y) :
    run.crossApproximation x y = K x y := by
  have hz := residual_eq_zero_of_complete_zero complete hzero x y hx hy
  have hh := run.kernel_sub_crossApproximation_eq_finalResidual x y
  linarith

/-- Every successful finite scalar Schur run preserves joint continuity. -/
theorem continuous_finalResidual (run : Run K)
    (hK : Continuous (fun z : ℝ × ℝ ↦ K z.1 z.2)) :
    Continuous (fun z : ℝ × ℝ ↦ run.finalResidual z.1 z.2) := by
  induction run with
  | nil => exact hK
  | @step K p q hp tail ih =>
    apply ih
    exact hK.sub (((hK.comp (continuous_fst.prodMk continuous_const)).mul
      (hK.comp (continuous_const.prodMk continuous_snd))).div_const (K p q))

theorem continuous_crossApproximation (run : Run K)
    (hK : Continuous (fun z : ℝ × ℝ ↦ K z.1 z.2)) :
    Continuous (fun z : ℝ × ℝ ↦ run.crossApproximation z.1 z.2) := by
  rw [run.crossApproximation_eq_sub_finalResidual]
  exact hK.sub (run.continuous_finalResidual hK)

/-- The actual difference of Green integrals is bounded by a finite-run
uniform cross error. Integrability of both terms is derived, not presumed. -/
theorem greenError_le_uniformCrossError (run : Run K)
    (hK : Continuous (fun z : ℝ × ℝ ↦ K z.1 z.2))
    (μ : MeasureTheory.Measure ℝ) (rho : ℝ → ℝ) (hrho : MeasureTheory.Integrable rho μ)
    (x B : ℝ) (hKbound : ∀ᵐ y ∂μ, |K x y| ≤ 1)
    (herror : ∀ᵐ y ∂μ, |run.finalResidual x y| ≤ B) :
    |(∫ y, rho y * K x y ∂μ) - (∫ y, rho y * run.crossApproximation x y ∂μ)| ≤
      B * ∫ y, |rho y| ∂μ := by
  have hKcont : Continuous (fun y : ℝ ↦ K x y) :=
    hK.comp (continuous_const.prodMk continuous_id)
  have hAcont : Continuous (fun y : ℝ ↦ run.crossApproximation x y) :=
    (run.continuous_crossApproximation hK).comp (continuous_const.prodMk continuous_id)
  have hiK : MeasureTheory.Integrable (fun y ↦ rho y * K x y) μ := by
    apply hrho.norm.mono' (hrho.aestronglyMeasurable.mul hKcont.aestronglyMeasurable)
    filter_upwards [hKbound] with y hy
    simp only [Pi.mul_apply, Real.norm_eq_abs, abs_mul]
    nlinarith [abs_nonneg (rho y)]
  have hiA : MeasureTheory.Integrable (fun y ↦ rho y * run.crossApproximation x y) μ := by
    apply (hrho.norm.const_mul (1 + B)).mono'
      (hrho.aestronglyMeasurable.mul hAcont.aestronglyMeasurable)
    filter_upwards [hKbound, herror] with y hy he
    have ha : |run.crossApproximation x y| ≤ 1 + B := by
      rw [run.crossApproximation_eq_sub_finalResidual]
      exact (abs_sub _ _).trans (add_le_add hy he)
    simp only [Pi.mul_apply, Real.norm_eq_abs, abs_mul]
    nlinarith [abs_nonneg (rho y)]
  rw [← MeasureTheory.integral_sub hiK hiA]
  have heq : (fun y ↦ rho y * K x y - rho y * run.crossApproximation x y) =
      fun y ↦ rho y * run.finalResidual x y := by
    funext y
    rw [← mul_sub, run.kernel_sub_crossApproximation_eq_finalResidual]
  rw [heq]
  exact greenError_le_kernelError_mul_l1 μ rho (fun y ↦ run.finalResidual x y) B hrho herror

end Run
end GECPKernelStructure.GECP
