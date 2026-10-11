import GECPKernelStructure.Matrix.AndreiefSymmetrization
import Mathlib.MeasureTheory.Integral.Pi

namespace Matrix

open scoped BigOperators
open Finset MeasureTheory

universe u v

variable {ι : Type u} [Fintype ι]
variable {X : Type v} [MeasurableSpace X]

/-- A finite product of one common measure is invariant under coordinate permutations. -/
theorem integral_comp_perm_pi
    (μ : Measure X) [SigmaFinite μ] (σ : Equiv.Perm ι)
    (f : (ι → X) → ℝ) :
    (∫ sample, f (sample ∘ σ) ∂Measure.pi (fun _ : ι => μ)) =
      ∫ sample, f sample ∂Measure.pi (fun _ : ι => μ) := by
  let e : (ι → X) ≃ᵐ (ι → X) :=
    MeasurableEquiv.piCongrLeft (fun _ : ι => X) σ.symm
  have preserving : MeasurePreserving e
      (Measure.pi (fun _ : ι => μ)) (Measure.pi (fun _ : ι => μ)) := by
    simpa only using
      (measurePreserving_piCongrLeft (fun _ : ι => μ) σ.symm)
  have integral_eq := preserving.integral_comp' f
  have e_apply (sample : ι → X) : e sample = sample ∘ σ := by
    funext i
    simpa [e, Function.comp_def] using
      (MeasurableEquiv.piCongrLeft_apply_apply
        (β := fun _ : ι => X) σ.symm sample (σ i))
  simpa only [e_apply] using integral_eq

/-- Integrability on a finite power is invariant under coordinate permutations. -/
theorem integrable_comp_perm_pi
    (μ : Measure X) [SigmaFinite μ] (σ : Equiv.Perm ι)
    (f : (ι → X) → ℝ) :
    Integrable (fun sample => f (sample ∘ σ)) (Measure.pi (fun _ : ι => μ)) ↔
      Integrable f (Measure.pi (fun _ : ι => μ)) := by
  let e : (ι → X) ≃ᵐ (ι → X) :=
    MeasurableEquiv.piCongrLeft (fun _ : ι => X) σ.symm
  have preserving : MeasurePreserving e
      (Measure.pi (fun _ : ι => μ)) (Measure.pi (fun _ : ι => μ)) := by
    simpa only using
      (measurePreserving_piCongrLeft (fun _ : ι => μ) σ.symm)
  have integrable_iff := preserving.integrable_comp_emb e.measurableEmbedding (g := f)
  have e_apply (sample : ι → X) : e sample = sample ∘ σ := by
    funext i
    simpa [e, Function.comp_def] using
      (MeasurableEquiv.piCongrLeft_apply_apply
        (β := fun _ : ι => X) σ.symm sample (σ i))
  have comp_e : f ∘ e = fun sample => f (sample ∘ σ) := by
    funext sample
    rw [Function.comp_apply, e_apply]
  rw [← comp_e]
  exact integrable_iff

/--
The integrated form of finite permutation symmetrization.  The factorial is
the exact normalization relating one row-wise composition integrand to the
product of the two sampled determinants.
-/
theorem factorial_mul_integral_compositionSampleIntegrand
    [DecidableEq ι] (μ : Measure X) [SigmaFinite μ]
    (left : ι → X → ℝ) (right : X → ι → ℝ)
    (hintegrable : Integrable (compositionSampleIntegrand left right)
      (Measure.pi (fun _ : ι => μ))) :
    ((Fintype.card ι).factorial : ℝ) *
        ∫ sample, compositionSampleIntegrand left right sample
          ∂Measure.pi (fun _ : ι => μ) =
      ∫ sample,
        Matrix.det (fun i j => left i (sample j)) *
          Matrix.det (fun i j => right (sample i) j)
        ∂Measure.pi (fun _ : ι => μ) := by
  let integrand := compositionSampleIntegrand left right
  have hintegrable_perm (σ : Equiv.Perm ι) :
      Integrable (fun sample => integrand (sample ∘ σ))
        (Measure.pi (fun _ : ι => μ)) :=
    (integrable_comp_perm_pi μ σ integrand).2 hintegrable
  calc
    ((Fintype.card ι).factorial : ℝ) *
          ∫ sample, integrand sample ∂Measure.pi (fun _ : ι => μ) =
        ∑ σ : Equiv.Perm ι,
          ∫ sample, integrand sample ∂Measure.pi (fun _ : ι => μ) := by
            simp [Fintype.card_perm, nsmul_eq_mul]
    _ = ∑ σ : Equiv.Perm ι,
          ∫ sample, integrand (sample ∘ σ) ∂Measure.pi (fun _ : ι => μ) := by
            apply Finset.sum_congr rfl
            intro σ _
            exact (integral_comp_perm_pi μ σ integrand).symm
    _ = ∫ sample, ∑ σ : Equiv.Perm ι, integrand (sample ∘ σ)
          ∂Measure.pi (fun _ : ι => μ) := by
            symm
            exact MeasureTheory.integral_finsetSum Finset.univ
              (fun σ _ => hintegrable_perm σ)
    _ = ∫ sample,
          Matrix.det (fun i j => left i (sample j)) *
            Matrix.det (fun i j => right (sample i) j)
          ∂Measure.pi (fun _ : ι => μ) := by
            apply MeasureTheory.integral_congr_ae
            filter_upwards with sample
            exact sum_compositionSampleIntegrand_comp_perm left right sample

omit [MeasurableSpace X] in
/-- The row-wise composition integrand expanded into coordinatewise products. -/
theorem compositionSampleIntegrand_eq_sum_prod
    [DecidableEq ι] (left : ι → X → ℝ) (right : X → ι → ℝ)
    (sample : ι → X) :
    compositionSampleIntegrand left right sample =
      ∑ σ : Equiv.Perm ι, (((σ.sign : ℤ) : ℝ) *
        ∏ i, left i (sample i) * right (sample i) (σ i)) := by
  let rightMatrix : Matrix ι ι ℝ := fun i j => right (sample i) j
  unfold compositionSampleIntegrand
  change (∏ i, left i (sample i)) * rightMatrix.det = _
  rw [← Matrix.det_transpose rightMatrix, Matrix.det_apply', Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro σ _
  change (∏ i, left i (sample i)) *
      (((σ.sign : ℤ) : ℝ) * ∏ i, right (sample i) (σ i)) = _
  rw [Finset.prod_mul_distrib]
  ring

/--
Integrating the row-wise composition integrand over a finite product measure
equals the determinant of the matrix of one-dimensional pair integrals.
-/
theorem integral_compositionSampleIntegrand_eq_det_integral
    [DecidableEq ι] (μ : Measure X) [SigmaFinite μ]
    (left : ι → X → ℝ) (right : X → ι → ℝ)
    (pair_integrable : ∀ i j, Integrable (fun x => left i x * right x j) μ) :
    (∫ sample, compositionSampleIntegrand left right sample
        ∂Measure.pi (fun _ : ι => μ)) =
      Matrix.det (fun i j => ∫ x, left i x * right x j ∂μ) := by
  let term : Equiv.Perm ι → (ι → X) → ℝ := fun σ sample =>
    (((σ.sign : ℤ) : ℝ) *
      ∏ i, left i (sample i) * right (sample i) (σ i))
  have product_integrable (σ : Equiv.Perm ι) :
      Integrable
        (fun sample : ι → X =>
          ∏ i, left i (sample i) * right (sample i) (σ i))
        (Measure.pi (fun _ : ι => μ)) :=
    Integrable.fintype_prod (fun i => pair_integrable i (σ i))
  have term_integrable (σ : Equiv.Perm ι) :
      Integrable (term σ) (Measure.pi (fun _ : ι => μ)) :=
    (product_integrable σ).const_mul _
  calc
    (∫ sample, compositionSampleIntegrand left right sample
        ∂Measure.pi (fun _ : ι => μ)) =
        ∫ sample, ∑ σ : Equiv.Perm ι, term σ sample
          ∂Measure.pi (fun _ : ι => μ) := by
            apply MeasureTheory.integral_congr_ae
            filter_upwards with sample
            exact compositionSampleIntegrand_eq_sum_prod left right sample
    _ = ∑ σ : Equiv.Perm ι,
          ∫ sample, term σ sample ∂Measure.pi (fun _ : ι => μ) := by
            exact MeasureTheory.integral_finsetSum Finset.univ
              (fun σ _ => term_integrable σ)
    _ = ∑ σ : Equiv.Perm ι, (((σ.sign : ℤ) : ℝ) *
          ∏ i, ∫ x, left i x * right x (σ i) ∂μ) := by
            apply Finset.sum_congr rfl
            intro σ _
            rw [MeasureTheory.integral_const_mul]
            congr 1
            exact MeasureTheory.integral_fintype_prod_eq_prod
              (fun i x => left i x * right x (σ i))
    _ = Matrix.det (fun i j => ∫ x, left i x * right x j ∂μ) := by
            let pairMatrix : Matrix ι ι ℝ :=
              fun i j => ∫ x, left i x * right x j ∂μ
            change _ = pairMatrix.det
            rw [← Matrix.det_transpose pairMatrix, Matrix.det_apply']
            rfl

/--
Finite-order Andréief identity over one common measure, stated without division:
the determinant of pair integrals times `n!` equals the product-measure integral
of the two sampled determinants.
-/
theorem factorial_mul_det_integral_eq_integral_det_mul_det
    [DecidableEq ι] (μ : Measure X) [SigmaFinite μ]
    (left : ι → X → ℝ) (right : X → ι → ℝ)
    (pair_integrable : ∀ i j, Integrable (fun x => left i x * right x j) μ) :
    ((Fintype.card ι).factorial : ℝ) *
        Matrix.det (fun i j => ∫ x, left i x * right x j ∂μ) =
      ∫ sample,
        Matrix.det (fun i j => left i (sample j)) *
          Matrix.det (fun i j => right (sample i) j)
        ∂Measure.pi (fun _ : ι => μ) := by
  have composition_integrable :
      Integrable (compositionSampleIntegrand left right)
        (Measure.pi (fun _ : ι => μ)) := by
    rw [show compositionSampleIntegrand left right =
      fun sample => ∑ σ : Equiv.Perm ι,
        (((σ.sign : ℤ) : ℝ) *
          ∏ i, left i (sample i) * right (sample i) (σ i)) by
      funext sample
      exact compositionSampleIntegrand_eq_sum_prod left right sample]
    apply integrable_finsetSum
    intro σ _
    exact (Integrable.fintype_prod
      (fun i => pair_integrable i (σ i))).const_mul _
  rw [← integral_compositionSampleIntegrand_eq_det_integral μ left right
    pair_integrable]
  exact factorial_mul_integral_compositionSampleIntegrand μ left right
    composition_integrable

end Matrix
