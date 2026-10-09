import GECPKernelStructure.Fermionic.PowerSecantResidual
import GECPKernelStructure.Fermionic.TwoCorner
import GECPKernelStructure.Fermionic.CompleteRun
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

/-- The actual residual after the two symmetric cutoff-corner pivots. -/
noncomputable def symmetricTwoCornerResidual (Λ : ℝ) (cutoff_pos : 0 < Λ) :
    Kernel ℝ ℝ ℝ :=
  residualUpdate
    (residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne')
    1 (-Λ) (fermionicKernel_firstPivot_reflected_ne cutoff_pos)

/-- A two-corner residual is bounded above by its first upper-corner residual. -/
theorem fermionicKernel_twoCornerResidual_le_firstCornerResidual
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_nonneg : 0 ≤ t) (frequency_upper : ω ≤ ωHigh) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ≤
        residualUpdate fermionicKernel 0 ωHigh
          (fermionicKernel_pos 0 ωHigh).ne' t ω := by
  let firstResidual := residualUpdate fermionicKernel 0 ωHigh
    (fermionicKernel_pos 0 ωHigh).ne'
  have first_at_lower_nonneg : 0 ≤ firstResidual t ωLow :=
    fermionicKernel_cornerResidual_nonneg time_nonneg frequency_order.le
  have first_at_row_nonneg : 0 ≤ firstResidual 1 ω :=
    fermionicKernel_cornerResidual_nonneg (by norm_num) frequency_upper
  have pivot_pos : 0 < firstResidual 1 ωLow :=
    fermionicKernel_firstCorner_secondPivot_pos frequency_order
  unfold residualUpdate
  exact sub_le_self _ (div_nonneg
    (mul_nonneg first_at_lower_nonneg first_at_row_nonneg) pivot_pos.le)

/-- A two-corner residual is also bounded by the reflected lower-corner residual. -/
theorem fermionicKernel_twoCornerResidual_le_secondCornerResidual
    {ωLow ωHigh t ω : ℝ} (frequency_order : ωLow < ωHigh)
    (time_le_one : t ≤ 1) (frequency_lower : ωLow ≤ ω)
    (frequency_upper : ω ≤ ωHigh) :
    residualUpdate
      (residualUpdate fermionicKernel 0 ωHigh
        (fermionicKernel_pos 0 ωHigh).ne')
      1 ωLow (fermionicKernel_firstCorner_secondPivot_pos frequency_order).ne' t ω ≤
        residualUpdate fermionicKernel 1 ωLow
          (fermionicKernel_pos 1 ωLow).ne' t ω := by
  have transformed_order : Real.exp (-ωHigh) < Real.exp (-ωLow) :=
    Real.exp_lt_exp.mpr (neg_lt_neg frequency_order)
  have transformed_lower : Real.exp (-ωHigh) ≤ Real.exp (-ω) :=
    Real.exp_le_exp.mpr (neg_le_neg frequency_upper)
  have transformed_upper : Real.exp (-ω) ≤ Real.exp (-ωLow) :=
    Real.exp_le_exp.mpr (neg_le_neg frequency_lower)
  have secant_lower := mul_rpow_sub_one_le_powerSecant
    (Real.exp_pos (-ωHigh)) transformed_order transformed_lower transformed_upper
      time_le_one
  rw [fermionicKernel_twoCornerResidual_eq_powerSecantError frequency_order,
    fermionicKernel_secondCornerResidual_eq_rpow_sub]
  exact div_le_div_of_nonneg_right (sub_le_sub_left secant_lower _) (by positivity)

/-- The lower-corner residual is the reflection of the upper-corner residual. -/
theorem fermionicKernel_secondCornerResidual_reflection (Λ t ω : ℝ) :
    residualUpdate fermionicKernel 1 (-Λ)
      (fermionicKernel_pos 1 (-Λ)).ne' t ω =
        residualUpdate fermionicKernel 0 Λ
          (fermionicKernel_pos 0 Λ).ne' (1 - t) (-ω) := by
  have kernel_reflection : fermionicKernel (1 - t) (-ω) = fermionicKernel t ω :=
    fermionicKernel_reflection t ω
  have upper_column_reflection :
      fermionicKernel (1 - t) Λ = fermionicKernel t (-Λ) := by
    simpa using fermionicKernel_reflection t (-Λ)
  have row_reflection : fermionicKernel 0 (-ω) = fermionicKernel 1 ω := by
    simpa using fermionicKernel_reflection 1 ω
  have pivot_reflection : fermionicKernel 0 Λ = fermionicKernel 1 (-Λ) := by
    simpa using fermionicKernel_reflection 1 (-Λ)
  unfold residualUpdate
  rw [kernel_reflection, upper_column_reflection, row_reflection, pivot_reflection]

/-- The upper-corner residual is at most one quarter on the positive outer half-band. -/
theorem fermionicKernel_firstCornerResidual_le_quarter_of_halfCutoff_le_frequency
    {Λ t ω : ℝ} (cutoff_pos : 0 < Λ) (time_nonneg : 0 ≤ t)
    (frequency_half : Λ / 2 ≤ ω) (frequency_upper : ω ≤ Λ) :
    residualUpdate fermionicKernel 0 Λ
      (fermionicKernel_pos 0 Λ).ne' t ω ≤ 1 / 4 := by
  let x := Real.exp (-t * (Λ / 2))
  have x_pos : 0 < x := Real.exp_pos _
  have x_le_one : x ≤ 1 := by
    dsimp [x]
    rw [Real.exp_le_one_iff]
    nlinarith [mul_nonneg time_nonneg cutoff_pos.le]
  have exp_frequency_le : Real.exp (-t * ω) ≤ x := by
    dsimp [x]
    exact Real.exp_le_exp.mpr (by nlinarith)
  have exp_cutoff_eq : Real.exp (-t * Λ) = x * x := by
    dsimp [x]
    rw [← Real.exp_add]
    congr 1
    ring
  rw [fermionicKernel_firstPivot_residual]
  have numerator_nonneg :
      0 ≤ Real.exp (-t * ω) - Real.exp (-t * Λ) := by
    exact sub_nonneg.mpr (Real.exp_le_exp.mpr (by nlinarith))
  calc
    (Real.exp (-t * ω) - Real.exp (-t * Λ)) / (1 + Real.exp (-ω)) ≤
        Real.exp (-t * ω) - Real.exp (-t * Λ) := by
      apply (div_le_iff₀ (fermionicKernel_denominator_pos ω)).2
      nlinarith [Real.exp_pos (-ω)]
    _ ≤ x - x * x := by rw [exp_cutoff_eq]; linarith
    _ ≤ 1 / 4 := by nlinarith [sq_nonneg (x - 1 / 2)]

/-- The symmetric two-corner residual is at most one quarter on the positive outer half. -/
theorem fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_frequency
    {Λ t ω : ℝ} (cutoff_pos : 0 < Λ) (time_nonneg : 0 ≤ t)
    (frequency_half : Λ / 2 ≤ ω) (frequency_upper : ω ≤ Λ) :
    symmetricTwoCornerResidual Λ cutoff_pos t ω ≤ 1 / 4 := by
  have frequency_order : -Λ < Λ := by linarith
  calc
    symmetricTwoCornerResidual Λ cutoff_pos t ω ≤
        residualUpdate fermionicKernel 0 Λ
          (fermionicKernel_pos 0 Λ).ne' t ω := by
      unfold symmetricTwoCornerResidual
      exact fermionicKernel_twoCornerResidual_le_firstCornerResidual
        frequency_order time_nonneg frequency_upper
    _ ≤ 1 / 4 :=
      fermionicKernel_firstCornerResidual_le_quarter_of_halfCutoff_le_frequency
        cutoff_pos time_nonneg frequency_half frequency_upper

/-- The symmetric two-corner residual is at most one quarter on the negative outer half. -/
theorem fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_frequency_le_negHalfCutoff
    {Λ t ω : ℝ} (cutoff_pos : 0 < Λ) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_half : ω ≤ -(Λ / 2)) :
    symmetricTwoCornerResidual Λ cutoff_pos t ω ≤ 1 / 4 := by
  have frequency_order : -Λ < Λ := by linarith
  calc
    symmetricTwoCornerResidual Λ cutoff_pos t ω ≤
        residualUpdate fermionicKernel 1 (-Λ)
          (fermionicKernel_pos 1 (-Λ)).ne' t ω := by
      unfold symmetricTwoCornerResidual
      exact fermionicKernel_twoCornerResidual_le_secondCornerResidual
        frequency_order time_le_one frequency_lower (by linarith)
    _ = residualUpdate fermionicKernel 0 Λ
        (fermionicKernel_pos 0 Λ).ne' (1 - t) (-ω) :=
      fermionicKernel_secondCornerResidual_reflection Λ t ω
    _ ≤ 1 / 4 :=
      fermionicKernel_firstCornerResidual_le_quarter_of_halfCutoff_le_frequency
        cutoff_pos (by linarith) (by linarith) (by linarith)

/-- Both outer frequency half-bands obey the uniform quarter bound. -/
theorem fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_absFrequency
    {Λ t ω : ℝ} (cutoff_pos : 0 < Λ)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ)
    (frequency_outer : Λ / 2 ≤ abs ω) :
    symmetricTwoCornerResidual Λ cutoff_pos t ω ≤ 1 / 4 := by
  by_cases frequency_nonneg : 0 ≤ ω
  · rw [abs_of_nonneg frequency_nonneg] at frequency_outer
    exact fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_frequency
      cutoff_pos time_nonneg frequency_outer frequency_upper
  · have frequency_nonpos : ω ≤ 0 := le_of_not_ge frequency_nonneg
    rw [abs_of_nonpos frequency_nonpos] at frequency_outer
    apply fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_frequency_le_negHalfCutoff
      cutoff_pos time_le_one frequency_lower
    linarith

/-- Exact value of the symmetric two-corner residual at the fixed center. -/
theorem fermionicKernel_symmetricTwoCornerResidual_center {Λ : ℝ}
    (cutoff_pos : 0 < Λ) :
    symmetricTwoCornerResidual Λ cutoff_pos (1 / 2) 0 =
      1 / 2 - 1 / (2 * Real.cosh (Λ / 2)) := by
  unfold symmetricTwoCornerResidual
  rw [fermionicKernel_twoCornerResidual_eq_sub_interpolate cutoff_pos]
  unfold twoCornerInterpolate cornerLeftWeight cornerRightWeight fermionicKernel
  norm_num
  have half_pos : 0 < Λ / 2 := by linarith
  have sinh_half_ne : Real.sinh (Λ / 2) ≠ 0 :=
    (Real.sinh_pos_iff.mpr half_pos).ne'
  have sinh_cutoff_ne : Real.sinh Λ ≠ 0 :=
    (Real.sinh_pos_iff.mpr cutoff_pos).ne'
  have double_angle : Real.sinh Λ =
      2 * Real.sinh (Λ / 2) * Real.cosh (Λ / 2) := by
    simpa only [show 2 * (Λ / 2) = Λ by ring] using
      Real.sinh_two_mul (Λ / 2)
  rw [double_angle]
  field_simp [sinh_half_ne, sinh_cutoff_ne]
  norm_num

/-- Above the explicit threshold, the center residual is strictly larger than one quarter. -/
theorem fermionicKernel_symmetricTwoCornerResidual_center_gt_quarter {Λ : ℝ}
    (cutoff_pos : 0 < Λ) (large_cutoff : 2 * Real.log 4 ≤ Λ) :
    1 / 4 < symmetricTwoCornerResidual Λ cutoff_pos (1 / 2) 0 := by
  rw [fermionicKernel_symmetricTwoCornerResidual_center cutoff_pos]
  have log_four_le_half : Real.log 4 ≤ Λ / 2 := by linarith
  have four_le_exp : (4 : ℝ) ≤ Real.exp (Λ / 2) := by
    rw [← Real.exp_log (by norm_num : (0 : ℝ) < 4)]
    exact Real.exp_le_exp.mpr log_four_le_half
  have cosh_gt_two : 2 < Real.cosh (Λ / 2) := by
    rw [Real.cosh_eq]
    calc
      (2 : ℝ) = 4 / 2 := by norm_num
      _ ≤ Real.exp (Λ / 2) / 2 := by gcongr
      _ < (Real.exp (Λ / 2) + Real.exp (-(Λ / 2))) / 2 := by
        apply div_lt_div_of_pos_right _ (by norm_num)
        exact lt_add_of_pos_right _ (Real.exp_pos _)
  have denominator_pos : 0 < 2 * Real.cosh (Λ / 2) := by
    nlinarith
  have inverse_lt_quarter : 1 / (2 * Real.cosh (Λ / 2)) < (1 / 4 : ℝ) := by
    apply (div_lt_iff₀ denominator_pos).2
    nlinarith
  linarith

/--
Every third complete pivot after the two symmetric cutoff corners has frequency
strictly inside the central half-band once `2 * log 4 ≤ Λ`.
-/
theorem fermionicKernel_thirdCompletePivot_frequency_lt_halfCutoff
    {Λ t ω : ℝ} (cutoff_pos : 0 < Λ) (large_cutoff : 2 * Real.log 4 ≤ Λ)
    (complete_pivot : CompletePivotOn (fun s : ℝ => 0 ≤ s ∧ s ≤ 1)
      (fun η : ℝ => -Λ ≤ η ∧ η ≤ Λ)
      (symmetricTwoCornerResidual Λ cutoff_pos) t ω) :
    abs ω < Λ / 2 := by
  rcases complete_pivot with
    ⟨⟨time_nonneg, time_le_one⟩, ⟨frequency_lower, frequency_upper⟩, maximal⟩
  by_contra frequency_not_inner
  have frequency_outer : Λ / 2 ≤ abs ω := le_of_not_gt frequency_not_inner
  have selected_upper :=
    fermionicKernel_symmetricTwoCornerResidual_le_quarter_of_halfCutoff_le_absFrequency
      cutoff_pos time_nonneg time_le_one frequency_lower frequency_upper frequency_outer
  have frequency_order : -Λ < Λ := by linarith
  have selected_nonneg : 0 ≤ symmetricTwoCornerResidual Λ cutoff_pos t ω := by
    unfold symmetricTwoCornerResidual
    exact fermionicKernel_twoCornerResidual_nonneg frequency_order time_nonneg time_le_one
      frequency_lower frequency_upper
  have center_dominates := maximal (1 / 2) 0
    ⟨by norm_num, by norm_num⟩ ⟨by linarith, by linarith⟩
  have center_pos : 0 < symmetricTwoCornerResidual Λ cutoff_pos (1 / 2) 0 :=
    (fermionicKernel_symmetricTwoCornerResidual_center_gt_quarter
      cutoff_pos large_cutoff).trans' (by norm_num)
  rw [abs_of_pos center_pos, abs_of_nonneg selected_nonneg] at center_dominates
  linarith [fermionicKernel_symmetricTwoCornerResidual_center_gt_quarter
    cutoff_pos large_cutoff]

end Fermionic
end GECPKernelStructure
