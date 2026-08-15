import Mathlib.MeasureTheory.Integral.Bochner.Basic

namespace GECPKernelStructure

open MeasureTheory

universe u

/-- Uniform kernel error transfers to an `L¹`-weighted Green-function error. -/
theorem greenError_le_kernelError_mul_l1 {Ω : Type u} [MeasurableSpace Ω]
    (μ : Measure Ω) (ρ kernelError : Ω → ℝ) (ε : ℝ)
    (rho_integrable : Integrable ρ μ)
    (kernel_error_le : ∀ᵐ ω ∂μ, |kernelError ω| ≤ ε) :
    |∫ ω, ρ ω * kernelError ω ∂μ| ≤ ε * ∫ ω, |ρ ω| ∂μ := by
  have dominating_integrable : Integrable (fun ω ↦ ε * |ρ ω|) μ := by
    simpa [Real.norm_eq_abs] using rho_integrable.norm.const_mul ε
  have pointwise : ∀ᵐ ω ∂μ, ‖ρ ω * kernelError ω‖ ≤ ε * |ρ ω| :=
    kernel_error_le.mono fun ω error_le ↦ by
      rw [Real.norm_eq_abs, abs_mul]
      exact mul_le_mul_of_nonneg_left error_le (abs_nonneg (ρ ω)) |>.trans_eq
        (mul_comm |ρ ω| ε)
  calc
    |∫ ω, ρ ω * kernelError ω ∂μ| = ‖∫ ω, ρ ω * kernelError ω ∂μ‖ :=
      (Real.norm_eq_abs _).symm
    _ ≤ ∫ ω, ε * |ρ ω| ∂μ :=
      norm_integral_le_of_norm_le dominating_integrable pointwise
    _ = ε * ∫ ω, |ρ ω| ∂μ := integral_const_mul ε _
