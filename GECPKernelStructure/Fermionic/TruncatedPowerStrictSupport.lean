import GECPKernelStructure.Fermionic.TailCompositionSign
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix MeasureTheory Set

/--
An explicit Schoenberg--Whitney sample: the first point is halfway from zero
to the first knot, and each later point is halfway between adjacent knots.
-/
noncomputable def interlacingSample {n : ℕ}
    (columns : Fin (n + 1) → ℝ) : Fin (n + 1) → ℝ :=
  Fin.cases (columns 0 / 2)
    (fun i => (columns i.castSucc + columns i.succ) / 2)

@[simp]
theorem interlacingSample_zero {n : ℕ} (columns : Fin (n + 1) → ℝ) :
    interlacingSample columns 0 = columns 0 / 2 := by
  rfl

@[simp]
theorem interlacingSample_succ {n : ℕ} (columns : Fin (n + 1) → ℝ)
    (i : Fin n) :
    interlacingSample columns i.succ =
      (columns i.castSucc + columns i.succ) / 2 := by
  rfl

/-- Every explicit interlacing sample lies strictly above zero. -/
theorem interlacingSample_pos {n : ℕ} {columns : Fin (n + 1) → ℝ}
    (columns_pos : ∀ j, 0 < columns j) :
    ∀ i, 0 < interlacingSample columns i := by
  intro i
  refine Fin.cases ?_ (fun k => ?_) i
  · simp only [interlacingSample_zero]
    linarith [columns_pos 0]
  · simp only [interlacingSample_succ]
    linarith [columns_pos k.castSucc, columns_pos k.succ]

/-- Each explicit interlacing sample lies below its matching knot. -/
theorem interlacingSample_lt_column {n : ℕ} {columns : Fin (n + 1) → ℝ}
    (columns_mono : StrictMono columns) (columns_pos : 0 < columns 0) :
    ∀ i, interlacingSample columns i < columns i := by
  intro i
  refine Fin.cases ?_ (fun k => ?_) i
  · simp only [interlacingSample_zero]
    linarith
  · simp only [interlacingSample_succ]
    have adjacent := columns_mono k.castSucc_lt_succ
    linarith

/-- A predecessor knot lies below the next explicit interlacing sample. -/
theorem column_lt_interlacingSample_succ {n : ℕ}
    {columns : Fin (n + 1) → ℝ} (columns_mono : StrictMono columns)
    (i : Fin n) :
    columns i.castSucc < interlacingSample columns i.succ := by
  simp only [interlacingSample_succ]
  have adjacent := columns_mono i.castSucc_lt_succ
  linarith

/-- The explicit interlacing sample is strictly increasing. -/
theorem interlacingSample_strictMono {n : ℕ}
    {columns : Fin (n + 1) → ℝ} (columns_mono : StrictMono columns)
    (columns_pos : 0 < columns 0) :
    StrictMono (interlacingSample columns) := by
  rw [Fin.strictMono_iff_lt_succ]
  intro i
  exact (interlacingSample_lt_column columns_mono columns_pos i.castSucc).trans
    (column_lt_interlacingSample_succ columns_mono i)

/-- Every earlier knot lies below a later explicit interlacing sample. -/
theorem column_lt_interlacingSample_of_lt {n : ℕ}
    {columns : Fin (n + 1) → ℝ} (columns_mono : StrictMono columns)
    {i j : Fin (n + 1)} (hji : j < i) :
    columns j < interlacingSample columns i := by
  have i_ne_zero : i ≠ 0 := by
    intro hi
    subst i
    exact (Fin.not_lt_zero j) hji
  let k : Fin n := i.pred i_ne_zero
  have k_succ : k.succ = i := Fin.succ_pred i i_ne_zero
  have hji' : j < k.succ := by simpa only [k_succ] using hji
  have j_le : j ≤ k.castSucc := by
    apply Fin.le_iff_val_le_val.mpr
    apply Nat.le_of_lt_succ
    exact hji'
  rw [← k_succ]
  exact (columns_mono.monotone j_le).trans_lt
    (column_lt_interlacingSample_succ columns_mono k)

/--
At the explicit interlacing sample the truncated-power matrix is upper
triangular with a positive diagonal, hence has positive determinant at every
order.
-/
theorem eighthOrderTruncatedPowerKernel_interlacing_det_pos
    {n : ℕ} (columns : Fin (n + 1) → ℝ)
    (columns_mono : StrictMono columns) (columns_pos : ∀ j, 0 < columns j) :
    0 < Matrix.det (fun i j =>
      eighthOrderTruncatedPowerKernel (interlacingSample columns i) (columns j)) := by
  let sampled : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    fun i j => eighthOrderTruncatedPowerKernel
      (interlacingSample columns i) (columns j)
  have upper_triangular : sampled.IsUpperTriangular := by
    intro i j hji
    have knot_lt := column_lt_interlacingSample_of_lt columns_mono
      (show j < i by simpa only [id_eq] using hji)
    change eighthOrderTruncatedPowerKernel
      (interlacingSample columns i) (columns j) = 0
    rw [eighthOrderTruncatedPowerKernel,
      max_eq_right (sub_nonpos.mpr knot_lt.le)]
    norm_num
  have diagonal_pos (i : Fin (n + 1)) : 0 < sampled i i := by
    have sample_lt := interlacingSample_lt_column columns_mono (columns_pos 0) i
    simp only [sampled, eighthOrderTruncatedPowerKernel]
    rw [max_eq_left (sub_nonneg.mpr sample_lt.le)]
    positivity
  change 0 < sampled.det
  rw [Matrix.det_of_isUpperTriangular upper_triangular]
  exact Finset.prod_pos fun i _ => diagonal_pos i

/--
At the explicit interlacing sample, the signed exponential/truncated-power
determinant product is strictly positive.
-/
theorem eighthOrderComposition_interlacing_pos
    {n : ℕ} (rows columns : Fin (n + 1) → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (columns_pos : ∀ j, 0 < columns j) :
    0 < expKernelSignature (n + 1) *
      (Matrix.det (fun i j =>
          expKernel (rows i) (interlacingSample columns j)) *
        Matrix.det (fun i j => eighthOrderTruncatedPowerKernel
          (interlacingSample columns i) (columns j))) := by
  have sample_mono := interlacingSample_strictMono
    columns_mono (columns_pos 0)
  have exponential_pos := expKernel_strictSignRegular (n + 1)
    rows (interlacingSample columns) rows_mono sample_mono
  have truncated_pos := eighthOrderTruncatedPowerKernel_interlacing_det_pos
    columns columns_mono columns_pos
  calc
    0 < (expKernelSignature (n + 1) *
          Matrix.det (fun i j => expKernel (rows i) (interlacingSample columns j))) *
        Matrix.det (fun i j => eighthOrderTruncatedPowerKernel
          (interlacingSample columns i) (columns j)) :=
      mul_pos exponential_pos truncated_pos
    _ = expKernelSignature (n + 1) *
        (Matrix.det (fun i j => expKernel (rows i) (interlacingSample columns j)) *
          Matrix.det (fun i j => eighthOrderTruncatedPowerKernel
            (interlacingSample columns i) (columns j))) := by ring

/--
The signed sampled determinant product has a strictly positive integral over
the common cutoff cube.  A compactly supported bump turns the explicit
interlacing witness into a positive-measure support certificate.
-/
theorem eighthOrderComposition_signed_cubeIntegral_pos
    {n : ℕ} (rows columns : Fin (n + 1) → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (columns_pos : ∀ j, 0 < columns j) :
    let upper := columns (Fin.last n)
    0 < ∫ sample in Icc (fun _ : Fin (n + 1) => (0 : ℝ)) (fun _ => upper),
      expKernelSignature (n + 1) *
        (Matrix.det (fun i j => expKernel (rows i) (sample j)) *
          Matrix.det (fun i j =>
            eighthOrderTruncatedPowerKernel (sample i) (columns j))) := by
  let upper := columns (Fin.last n)
  let cube : Set (Fin (n + 1) → ℝ) :=
    Icc (fun _ => (0 : ℝ)) (fun _ => upper)
  let signedIntegrand : (Fin (n + 1) → ℝ) → ℝ := fun sample =>
    expKernelSignature (n + 1) *
      (Matrix.det (fun i j => expKernel (rows i) (sample j)) *
        Matrix.det (fun i j =>
          eighthOrderTruncatedPowerKernel (sample i) (columns j)))
  let bump : (Fin (n + 1) → ℝ) → ℝ := fun sample =>
    ∏ i, max (sample i) 0 * max (upper - sample i) 0
  let localized : (Fin (n + 1) → ℝ) → ℝ := fun sample =>
    bump sample * signedIntegrand sample
  have signed_continuous : Continuous signedIntegrand := by
    dsimp only [signedIntegrand]
    apply Continuous.mul continuous_const
    apply Continuous.mul
    · apply Continuous.matrix_det
      unfold expKernel
      fun_prop
    · apply Continuous.matrix_det
      unfold eighthOrderTruncatedPowerKernel
      fun_prop
  have bump_continuous : Continuous bump := by
    dsimp only [bump]
    fun_prop
  have localized_continuous : Continuous localized := by
    exact bump_continuous.mul signed_continuous
  have signed_nonneg (sample : Fin (n + 1) → ℝ) :
      0 ≤ signedIntegrand sample := by
    exact eighthOrderComposition_sampledDeterminants_nonneg
      rows columns sample rows_mono columns_mono
  have bump_nonneg (sample : Fin (n + 1) → ℝ) : 0 ≤ bump sample := by
    dsimp only [bump]
    exact Finset.prod_nonneg fun _ _ =>
      mul_nonneg (le_max_right _ _) (le_max_right _ _)
  have localized_nonneg : 0 ≤ localized := by
    intro sample
    exact mul_nonneg (bump_nonneg sample) (signed_nonneg sample)
  have sample_pos : ∀ i, 0 < interlacingSample columns i :=
    interlacingSample_pos columns_pos
  have sample_lt_upper : ∀ i, interlacingSample columns i < upper := by
    intro i
    exact (interlacingSample_lt_column columns_mono (columns_pos 0) i).trans_le
      (columns_mono.monotone (Fin.le_last i))
  have bump_witness_pos : 0 < bump (interlacingSample columns) := by
    dsimp only [bump]
    apply Finset.prod_pos
    intro i _
    have lower := sample_pos i
    have upper_bound := sample_lt_upper i
    rw [max_eq_left lower.le, max_eq_left (sub_nonneg.mpr upper_bound.le)]
    exact mul_pos lower (sub_pos.mpr upper_bound)
  have signed_witness_pos : 0 < signedIntegrand (interlacingSample columns) :=
    eighthOrderComposition_interlacing_pos rows columns
      rows_mono columns_mono columns_pos
  have localized_witness_pos : 0 < localized (interlacingSample columns) :=
    mul_pos bump_witness_pos signed_witness_pos
  have localized_zero_outside (sample : Fin (n + 1) → ℝ)
      (sample_outside : sample ∉ cube) : localized sample = 0 := by
    have bump_zero : bump sample = 0 := by
      by_contra bump_ne
      have factor_ne (i : Fin (n + 1)) :
          max (sample i) 0 * max (upper - sample i) 0 ≠ 0 :=
        Finset.prod_ne_zero_iff.mp bump_ne i (Finset.mem_univ i)
      have lower_bound (i : Fin (n + 1)) : 0 ≤ sample i := by
        by_contra h
        apply factor_ne i
        rw [max_eq_right (le_of_not_ge h), zero_mul]
      have upper_bound (i : Fin (n + 1)) : sample i ≤ upper := by
        by_contra h
        apply factor_ne i
        rw [max_eq_right (sub_nonpos.mpr (le_of_not_ge h)), mul_zero]
      apply sample_outside
      exact ⟨lower_bound, upper_bound⟩
    change bump sample * signedIntegrand sample = 0
    rw [bump_zero, zero_mul]
  have localized_compact : HasCompactSupport localized :=
    HasCompactSupport.intro isCompact_Icc localized_zero_outside
  have localized_integrable : Integrable localized volume :=
    localized_continuous.integrable_of_hasCompactSupport localized_compact
  have localized_integral_pos : 0 < ∫ sample, localized sample :=
    MeasureTheory.integral_pos_of_integrable_nonneg_nonzero
      localized_continuous localized_integrable localized_nonneg
      localized_witness_pos.ne'
  have localized_support_pos :
      0 < (volume : Measure (Fin (n + 1) → ℝ))
        (Function.support localized) :=
    (MeasureTheory.integral_pos_iff_support_of_nonneg
      localized_nonneg localized_integrable).mp localized_integral_pos
  have support_subset : Function.support localized ⊆
      Function.support signedIntegrand ∩ cube := by
    intro sample sample_mem
    constructor
    · intro signed_zero
      apply sample_mem
      simp [localized, signed_zero]
    · by_contra sample_outside
      exact sample_mem (localized_zero_outside sample sample_outside)
  have signed_support_pos :
      0 < (volume : Measure (Fin (n + 1) → ℝ))
        (Function.support signedIntegrand ∩ cube) :=
    localized_support_pos.trans_le (measure_mono support_subset)
  have signed_integrable_on : IntegrableOn signedIntegrand cube volume :=
    signed_continuous.continuousOn.integrableOn_compact isCompact_Icc
  change 0 < ∫ sample in cube, signedIntegrand sample
  exact (MeasureTheory.setIntegral_pos_iff_support_of_nonneg_ae
    (Filter.Eventually.of_forall fun sample => signed_nonneg sample)
    signed_integrable_on).mpr signed_support_pos

/--
Strict-sign closure for the common-cutoff exponential/truncated-power
composition at every positive order.
-/
theorem eighthOrderComposition_interval_det_pos
    {n : ℕ} (rows columns : Fin (n + 1) → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (columns_pos : ∀ j, 0 < columns j) :
    let upper := columns (Fin.last n)
    0 < expKernelSignature (n + 1) *
      Matrix.det (fun i j => ∫ v in (0 : ℝ)..upper,
        expKernel (rows i) v *
          eighthOrderTruncatedPowerKernel v (columns j)) := by
  let upper := columns (Fin.last n)
  let interval : Set ℝ := Icc (0 : ℝ) upper
  let cube : Set (Fin (n + 1) → ℝ) :=
    Icc (fun _ => (0 : ℝ)) (fun _ => upper)
  let left : Fin (n + 1) → ℝ → ℝ := fun i v => expKernel (rows i) v
  let right : ℝ → Fin (n + 1) → ℝ := fun v j =>
    eighthOrderTruncatedPowerKernel v (columns j)
  have upper_pos : 0 < upper := columns_pos (Fin.last n)
  have pair_integrable (i j : Fin (n + 1)) : Integrable
      (fun v => left i v * right v j) (volume.restrict interval) := by
    apply Continuous.integrableOn_Icc
    dsimp only [left, right]
    unfold expKernel eighthOrderTruncatedPowerKernel
    fun_prop
  have andreief := Matrix.factorial_mul_det_integral_eq_integral_det_mul_det
    (volume.restrict interval) left right pair_integrable
  have product_measure_eq :
      Measure.pi (fun _ : Fin (n + 1) => volume.restrict interval) =
        volume.restrict cube := by
    calc
      Measure.pi (fun _ : Fin (n + 1) => volume.restrict interval) =
          (Measure.pi (fun _ : Fin (n + 1) => volume)).restrict
            (Set.univ.pi fun _ : Fin (n + 1) => interval) :=
        (Measure.restrict_pi_pi (μ := fun _ : Fin (n + 1) => volume)
          (fun _ : Fin (n + 1) => interval)).symm
      _ = volume.restrict cube := by
        rw [← MeasureTheory.volume_pi, Set.pi_univ_Icc]
  have signed_product_integral_pos :
      0 < expKernelSignature (n + 1) *
        ∫ sample,
          Matrix.det (fun i j => left i (sample j)) *
            Matrix.det (fun i j => right (sample i) j)
          ∂Measure.pi (fun _ : Fin (n + 1) => volume.restrict interval) := by
    rw [← MeasureTheory.integral_const_mul, product_measure_eq]
    exact eighthOrderComposition_signed_cubeIntegral_pos
      rows columns rows_mono columns_mono columns_pos
  rw [← andreief] at signed_product_integral_pos
  have factorial_pos : 0 < ((Nat.factorial (n + 1) : ℕ) : ℝ) := by positivity
  have signed_set_det_pos :
      0 < expKernelSignature (n + 1) *
        Matrix.det (fun i j => ∫ v in interval,
          expKernel (rows i) v *
            eighthOrderTruncatedPowerKernel v (columns j)) := by
    apply pos_of_mul_pos_right _ factorial_pos.le
    simpa only [Fintype.card_fin, left, right, interval,
      mul_assoc, mul_left_comm, mul_comm] using signed_product_integral_pos
  simpa only [intervalIntegral.integral_of_le upper_pos.le,
    ← MeasureTheory.integral_Icc_eq_integral_Ioc, interval, upper] using signed_set_det_pos

/--
Every positive-order Taylor-tail minor on strictly ordered positive times and
strictly ordered positive frequencies has the exponential signature strictly.
-/
theorem eighthOrderTail_minor_signed_pos
    {n : ℕ} (rows columns : Fin (n + 1) → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (rows_pos : ∀ i, 0 < rows i) (columns_pos : ∀ j, 0 < columns j) :
    0 < expKernelSignature (n + 1) *
      Matrix.det (fun i j => eighthOrderTail (rows i * columns j)) := by
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
  have composition_pos :
      0 < expKernelSignature (n + 1) * composition.det := by
    exact eighthOrderComposition_interval_det_pos
      rows columns rows_mono columns_mono columns_pos
  have tail_matrix_eq : tailMatrix = scaled := by
    ext i j
    dsimp only [tailMatrix, scaled, weights, composition]
    rw [eighthOrderTail_eq_truncatedPowerComposition_commonCutoff
      (columns_pos j).le
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
  change 0 < expKernelSignature (n + 1) * tailMatrix.det
  rw [tail_matrix_eq, scaled_det]
  calc
    0 < (∏ i, weights i) *
        (expKernelSignature (n + 1) * composition.det) :=
      mul_pos weights_product_pos composition_pos
    _ = expKernelSignature (n + 1) *
        ((∏ i, weights i) * composition.det) := by ring

end Fermionic
end GECPKernelStructure
