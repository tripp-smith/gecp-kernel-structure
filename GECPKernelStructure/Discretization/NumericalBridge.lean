import GECPKernelStructure.Discretization.LipschitzCertificate
import Mathlib.Analysis.Complex.Norm

namespace GECPKernelStructure.Discretization

/-- A uniform discrepancy transfers the quality of a rounded candidate to the
true residual on exactly the same domain. In a finite-grid application the
domain remains that grid; this does not supply a continuous-domain bound. -/
theorem approxPivot_of_uniform_discrepancy {X : Type*} (R Rhat : X → ℝ)
    (p : X) {delta L U : ℝ} (delta_nonneg : 0 ≤ delta)
    (delta_lt : delta < L) (discrepancy : ∀ x, |R x - Rhat x| ≤ delta)
    (candidate_lower : L ≤ |Rhat p|) (global_upper : ∀ x, |Rhat x| ≤ U) :
    0 < (L - delta) / (U + delta) ∧
      (L - delta) / (U + delta) ≤ 1 ∧
      ∀ x, (L - delta) / (U + delta) * |R x| ≤ |R p| := by
  have L_le_U : L ≤ U := candidate_lower.trans (global_upper p)
  have den_pos : 0 < U + delta := by linarith
  have num_pos : 0 < L - delta := sub_pos.mpr delta_lt
  have eta_pos : 0 < (L - delta) / (U + delta) := div_pos num_pos den_pos
  refine ⟨eta_pos, (div_le_one den_pos).mpr (by linarith), ?_⟩
  have candidate_true : L - delta ≤ |R p| := by
    have h := abs_sub_abs_le_abs_sub (Rhat p) (R p)
    rw [abs_sub_comm] at h
    linarith [discrepancy p]
  intro x
  have true_upper : |R x| ≤ U + delta := by
    have h := abs_sub_abs_le_abs_sub (R x) (Rhat x)
    linarith [discrepancy x, global_upper x]
  calc
    (L - delta) / (U + delta) * |R x|
      ≤ (L - delta) / (U + delta) * (U + delta) :=
        mul_le_mul_of_nonneg_left true_upper eta_pos.le
    _ = L - delta := div_mul_cancel₀ _ den_pos.ne'
    _ ≤ |R p| := candidate_true

/-- A computed complex dot product inherits an error budget from its final
arithmetic, its rounded design entries and its computed coefficients. The ideal
coefficients can be exact selected-core or least-squares solutions. No stability
property of the algorithm producing the rounded values is assumed. -/
theorem dot_error_le_of_component_bounds {ι : Type*} (s : Finset ι)
    (a ahat c chat : ι → ℂ) (computed : ℂ)
    (design_error coefficient_error : ι → ℝ) (dot_error : ℝ)
    (design_bound : ∀ i ∈ s, ‖ahat i - a i‖ ≤ design_error i)
    (coefficient_bound : ∀ i ∈ s, ‖chat i - c i‖ ≤ coefficient_error i)
    (dot_bound : ‖computed - ∑ i ∈ s, ahat i * chat i‖ ≤ dot_error) :
    ‖computed - ∑ i ∈ s, a i * c i‖ ≤ dot_error +
      ∑ i ∈ s, (design_error i * ‖chat i‖ + ‖a i‖ * coefficient_error i) := by
  have entry_bound (i : ι) (hi : i ∈ s) :
      ‖ahat i * chat i - a i * c i‖ ≤
        design_error i * ‖chat i‖ + ‖a i‖ * coefficient_error i := by
    have entry_eq : ahat i * chat i - a i * c i =
        (ahat i - a i) * chat i + a i * (chat i - c i) := by ring
    rw [entry_eq]
    calc
      ‖(ahat i - a i) * chat i + a i * (chat i - c i)‖
        ≤ ‖(ahat i - a i) * chat i‖ + ‖a i * (chat i - c i)‖ := norm_add_le _ _
      _ = ‖ahat i - a i‖ * ‖chat i‖ + ‖a i‖ * ‖chat i - c i‖ := by
        rw [Complex.norm_mul, Complex.norm_mul]
      _ ≤ design_error i * ‖chat i‖ + ‖a i‖ * coefficient_error i :=
        add_le_add (mul_le_mul_of_nonneg_right (design_bound i hi) (norm_nonneg _))
          (mul_le_mul_of_nonneg_left (coefficient_bound i hi) (norm_nonneg _))
  calc
    ‖computed - ∑ i ∈ s, a i * c i‖
      = ‖(computed - ∑ i ∈ s, ahat i * chat i) +
          ((∑ i ∈ s, ahat i * chat i) - ∑ i ∈ s, a i * c i)‖ := by
        congr 1
        abel
    _ ≤ ‖computed - ∑ i ∈ s, ahat i * chat i‖ +
        ‖(∑ i ∈ s, ahat i * chat i) - ∑ i ∈ s, a i * c i‖ := norm_add_le _ _
    _ ≤ dot_error + ∑ i ∈ s,
        (design_error i * ‖chat i‖ + ‖a i‖ * coefficient_error i) := by
      apply add_le_add dot_bound
      rw [← Finset.sum_sub_distrib]
      exact (norm_sum_le _ _).trans (Finset.sum_le_sum entry_bound)

end GECPKernelStructure.Discretization
