import GECPKernelStructure.Fermionic.TailCompositionSymmetrization
import GECPKernelStructure.Matrix.OrderedSampleSign

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix MeasureTheory Set

/--
The exponential/truncated-power sampled determinant product has the prescribed
weak sign for every sample, not merely for samples already in strict order.
-/
theorem eighthOrderComposition_sampledDeterminants_nonneg
    {n : ℕ} (rows columns sample : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 ≤ expKernelSignature n *
      (Matrix.det (fun i j => expKernel (rows i) (sample j)) *
        Matrix.det (fun i j =>
          eighthOrderTruncatedPowerKernel (sample i) (columns j))) := by
  apply Matrix.signed_det_mul_det_nonneg_of_strictMono
    (expKernelSignature n)
    (fun i v => expKernel (rows i) v)
    (fun v j => eighthOrderTruncatedPowerKernel v (columns j))
  intro ordered ordered_mono
  have symmetrized := eighthOrderComposition_symmetrized_nonneg
    rows columns ordered rows_mono columns_mono ordered_mono
  rw [Matrix.sum_compositionSampleIntegrand_comp_perm] at symmetrized
  exact symmetrized

/--
Weak-sign closure for the exponential/truncated-power integral composition
over any common sigma-finite measure on the intermediate variable.
-/
theorem eighthOrderComposition_det_integral_nonneg
    {n : ℕ} (μ : Measure ℝ) [SigmaFinite μ]
    (rows columns : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (pair_integrable : ∀ i j, Integrable
      (fun v => expKernel (rows i) v *
        eighthOrderTruncatedPowerKernel v (columns j)) μ) :
    0 ≤ expKernelSignature n *
      Matrix.det (fun i j => ∫ v,
        expKernel (rows i) v *
          eighthOrderTruncatedPowerKernel v (columns j) ∂μ) := by
  apply Matrix.signature_mul_det_integral_nonneg μ
    (expKernelSignature n)
    (fun i v => expKernel (rows i) v)
    (fun v j => eighthOrderTruncatedPowerKernel v (columns j))
    pair_integrable
  intro ordered ordered_mono
  exact eighthOrderComposition_sampledDeterminants_nonneg
    rows columns ordered rows_mono columns_mono

/--
The common-cutoff exponential/truncated-power composition has the prescribed
weak sign at every order.
-/
theorem eighthOrderComposition_interval_det_nonneg
    {n : ℕ} {upper : ℝ} (upper_nonneg : 0 ≤ upper)
    (rows columns : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 ≤ expKernelSignature n *
      Matrix.det (fun i j => ∫ v in (0 : ℝ)..upper,
        expKernel (rows i) v *
          eighthOrderTruncatedPowerKernel v (columns j)) := by
  have pair_integrable (i j : Fin n) : Integrable
      (fun v => expKernel (rows i) v *
        eighthOrderTruncatedPowerKernel v (columns j))
      (volume.restrict (Icc (0 : ℝ) upper)) := by
    apply Continuous.integrableOn_Icc
    unfold expKernel eighthOrderTruncatedPowerKernel
    fun_prop
  have result := eighthOrderComposition_det_integral_nonneg
    (volume.restrict (Icc (0 : ℝ) upper)) rows columns
    rows_mono columns_mono pair_integrable
  simpa only [intervalIntegral.integral_of_le upper_nonneg,
    ← MeasureTheory.integral_Icc_eq_integral_Ioc] using result

/--
Every ordered minor of the eighth-order Taylor-tail kernel on positive times
and nonnegative frequencies has the exponential kernel's prescribed weak sign.
-/
theorem eighthOrderTail_minor_signed_nonneg
    {n : ℕ} (rows columns : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_pos : ∀ i, 0 < rows i) (columns_nonneg : ∀ j, 0 ≤ columns j) :
    0 ≤ expKernelSignature n *
      Matrix.det (fun i j => eighthOrderTail (rows i * columns j)) := by
  cases n with
  | zero =>
      let tailMatrix : Matrix (Fin 0) (Fin 0) ℝ :=
        fun i j => eighthOrderTail (rows i * columns j)
      change 0 ≤ expKernelSignature 0 * tailMatrix.det
      rw [Matrix.det_fin_zero]
      norm_num [expKernelSignature]
  | succ n =>
      let upper := columns (Fin.last n)
      let composition : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
        fun i j => ∫ v in (0 : ℝ)..upper,
          expKernel (rows i) v *
            eighthOrderTruncatedPowerKernel v (columns j)
      let weights : Fin (n + 1) → ℝ :=
        fun i => rows i ^ 8 / (Nat.factorial 7 : ℝ)
      let scaled : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
        fun i j => weights i * composition i j
      let tailMatrix : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
        fun i j => eighthOrderTail (rows i * columns j)
      have upper_nonneg : 0 ≤ upper := columns_nonneg (Fin.last n)
      have composition_nonneg :
          0 ≤ expKernelSignature (n + 1) * composition.det := by
        exact eighthOrderComposition_interval_det_nonneg upper_nonneg
          rows columns rows_mono columns_mono
      have tail_matrix_eq : tailMatrix = scaled := by
        ext i j
        dsimp only [tailMatrix, scaled, weights, composition]
        rw [eighthOrderTail_eq_truncatedPowerComposition_commonCutoff
          (columns_nonneg j)
          (columns_mono.monotone (Fin.le_last j))]
        congr 1
        apply intervalIntegral.integral_congr
        intro v _
        simp only [expKernel]
        ring
      have weights_product_pos : 0 < ∏ i, weights i := by
        apply Finset.prod_pos
        intro i _
        exact div_pos (pow_pos (rows_pos i) _) (by positivity)
      have scaled_det : scaled.det = (∏ i, weights i) * composition.det := by
        exact Matrix.det_mul_column weights composition
      change 0 ≤ expKernelSignature (n + 1) * tailMatrix.det
      rw [tail_matrix_eq, scaled_det]
      calc
        0 ≤ (∏ i, weights i) *
            (expKernelSignature (n + 1) * composition.det) :=
          mul_nonneg weights_product_pos.le composition_nonneg
        _ = expKernelSignature (n + 1) *
            ((∏ i, weights i) * composition.det) := by ring

end Fermionic
end GECPKernelStructure
