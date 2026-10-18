import Mathlib.Algebra.Order.Chebyshev
import Mathlib.Tactic

namespace Matrix

open Finset

variable {ι : Type*}

/-- Total weight on a finite index set. -/
noncomputable def weightedMass (s : Finset ι) (w : ι → ℝ) : ℝ :=
  ∑ i ∈ s, w i

/-- An unnormalized first moment on a finite index set. -/
noncomputable def weightedMoment (s : Finset ι) (w f : ι → ℝ) : ℝ :=
  ∑ i ∈ s, w i * f i

/-- The numerator of the weighted covariance: mass times the mixed moment
minus the product of the two first moments. -/
noncomputable def weightedCovarianceNumerator
    (s : Finset ι) (w f g : ι → ℝ) : ℝ :=
  weightedMass s w * (∑ i ∈ s, w i * f i * g i) -
    weightedMoment s w f * weightedMoment s w g

/-- Weighted covariance is exactly half a pairwise sum of same-direction
increments. -/
theorem two_mul_weightedCovarianceNumerator
    (s : Finset ι) (w f g : ι → ℝ) :
    2 * weightedCovarianceNumerator s w f g =
      ∑ i ∈ s, ∑ j ∈ s,
        w i * w j * (f i - f j) * (g i - g j) := by
  have diagonal_swap :
      (∑ i ∈ s, ∑ j ∈ s, w i * w j * f j * g j) =
        ∑ i ∈ s, ∑ j ∈ s, w i * w j * f i * g i := by
    rw [sum_comm]
    apply sum_congr rfl
    intro i hi
    apply sum_congr rfl
    intro j hj
    ring
  have cross_swap :
      (∑ i ∈ s, ∑ j ∈ s, w i * w j * f i * g j) =
        ∑ i ∈ s, ∑ j ∈ s, w i * w j * f j * g i := by
    rw [sum_comm]
    apply sum_congr rfl
    intro i hi
    apply sum_congr rfl
    intro j hj
    ring
  have left_diagonal_comm :
      (∑ i ∈ s, ∑ j ∈ s, w j * w i * f i * g i) =
        ∑ i ∈ s, ∑ j ∈ s, w i * w j * f i * g i := by
    apply sum_congr rfl
    intro i hi
    apply sum_congr rfl
    intro j hj
    ring
  have left_cross_comm :
      (∑ i ∈ s, ∑ j ∈ s, w j * f j * w i * g i) =
        ∑ i ∈ s, ∑ j ∈ s, w i * w j * f j * g i := by
    apply sum_congr rfl
    intro i hi
    apply sum_congr rfl
    intro j hj
    ring
  unfold weightedCovarianceNumerator weightedMass weightedMoment
  simp only [sum_mul, mul_sum]
  ring_nf
  simp only [sum_sub_distrib, sum_add_distrib]
  rw [diagonal_swap, cross_swap, left_diagonal_comm, left_cross_comm]
  ring

/-- Weighted Chebyshev inequality in covariance-numerator form. -/
theorem weightedCovarianceNumerator_nonneg
    (s : Finset ι) (w f g : ι → ℝ)
    (weight_nonneg : ∀ i ∈ s, 0 ≤ w i)
    (monovary : MonovaryOn f g s) :
    0 ≤ weightedCovarianceNumerator s w f g := by
  have pairwise_nonneg : ∀ i ∈ s, ∀ j ∈ s,
      0 ≤ w i * w j * (f i - f j) * (g i - g j) := by
    intro i hi j hj
    have weights_nonneg : 0 ≤ w i * w j :=
      mul_nonneg (weight_nonneg i hi) (weight_nonneg j hj)
    have increments_nonneg : 0 ≤ (f i - f j) * (g i - g j) := by
      calc
        0 ≤ (f j - f i) * (g j - g i) :=
          monovary.sub_mul_sub_nonneg hi hj
        _ = (f i - f j) * (g i - g j) := by ring
    simpa only [mul_assoc] using mul_nonneg weights_nonneg increments_nonneg
  have double_sum_nonneg :
      0 ≤ ∑ i ∈ s, ∑ j ∈ s,
        w i * w j * (f i - f j) * (g i - g j) := by
    exact sum_nonneg fun i hi ↦ sum_nonneg fun j hj ↦ pairwise_nonneg i hi j hj
  rw [← two_mul_weightedCovarianceNumerator] at double_sum_nonneg
  linarith

/-- The sum over `group` of feature products centered and scaled by moments
over the ambient `universe`.  The scaling avoids all division by the ambient
mass. -/
noncomputable def weightedCenteredGroup
    (ambient group : Finset ι) (w f g : ι → ℝ) : ℝ :=
  ∑ i ∈ group, w i *
    (weightedMass ambient w * f i - weightedMoment ambient w f) *
    (weightedMass ambient w * g i - weightedMoment ambient w g)

theorem weightedCenteredGroup_eq
    (ambient group : Finset ι) (w f g : ι → ℝ) :
    weightedCenteredGroup ambient group w f g =
      weightedMass ambient w ^ 2 * (∑ i ∈ group, w i * f i * g i) -
        weightedMass ambient w * weightedMoment ambient w f *
          weightedMoment group w g -
        weightedMass ambient w * weightedMoment ambient w g *
          weightedMoment group w f +
        weightedMass group w * weightedMoment ambient w f *
          weightedMoment ambient w g := by
  unfold weightedCenteredGroup
  calc
    _ = ∑ i ∈ group, (
        weightedMass ambient w ^ 2 * (w i * f i * g i) -
          weightedMass ambient w * weightedMoment ambient w f * (w i * g i) -
          weightedMass ambient w * weightedMoment ambient w g * (w i * f i) +
          w i * weightedMoment ambient w f * weightedMoment ambient w g) := by
      apply sum_congr rfl
      intro i hi
      ring
    _ = _ := by
      simp only [sum_sub_distrib, sum_add_distrib]
      simp_rw [← mul_sum, ← sum_mul]
      change
        weightedMass ambient w ^ 2 * (∑ i ∈ group, w i * f i * g i) -
            weightedMass ambient w * weightedMoment ambient w f *
              weightedMoment group w g -
            weightedMass ambient w * weightedMoment ambient w g *
              weightedMoment group w f +
            weightedMass group w * weightedMoment ambient w f *
              weightedMoment ambient w g = _
      ring

/-- Conditional-covariance decomposition for a group centered at ambient
weighted means.  The second term is the product of the two group-mean shifts. -/
theorem weightedMass_mul_weightedCenteredGroup
    (ambient group : Finset ι) (w f g : ι → ℝ) :
    weightedMass group w * weightedCenteredGroup ambient group w f g =
      weightedMass ambient w ^ 2 *
          weightedCovarianceNumerator group w f g +
        (weightedMoment group w f * weightedMass ambient w -
            weightedMass group w * weightedMoment ambient w f) *
          (weightedMoment group w g * weightedMass ambient w -
            weightedMass group w * weightedMoment ambient w g) := by
  rw [weightedCenteredGroup_eq]
  unfold weightedCovarianceNumerator
  ring

/-- A positive-mass group has a nonnegative globally centered product sum once
its local covariance and its two mean shifts have the same sign. -/
theorem weightedCenteredGroup_nonneg
    (ambient group : Finset ι) (w f g : ι → ℝ)
    (group_mass_pos : 0 < weightedMass group w)
    (weight_nonneg : ∀ i ∈ group, 0 ≤ w i)
    (monovary : MonovaryOn f g group)
    (mean_shifts_nonneg : 0 ≤
      (weightedMoment group w f * weightedMass ambient w -
          weightedMass group w * weightedMoment ambient w f) *
        (weightedMoment group w g * weightedMass ambient w -
          weightedMass group w * weightedMoment ambient w g)) :
    0 ≤ weightedCenteredGroup ambient group w f g := by
  have covariance_nonneg :=
    weightedCovarianceNumerator_nonneg group w f g weight_nonneg monovary
  have covariance_term_nonneg :
      0 ≤ weightedMass ambient w ^ 2 *
        weightedCovarianceNumerator group w f g :=
    mul_nonneg (sq_nonneg _) covariance_nonneg
  have scaled_nonneg :
      0 ≤ weightedMass group w *
        weightedCenteredGroup ambient group w f g := by
    rw [weightedMass_mul_weightedCenteredGroup]
    exact add_nonneg covariance_term_nonneg mean_shifts_nonneg
  exact nonneg_of_mul_nonneg_left (by simpa [mul_comm] using scaled_nonneg)
    group_mass_pos

/-- The difference between two weighted block means, before normalization, is
the pairwise sum of their cross-block increments. -/
theorem weightedMoment_mul_mass_sub_mass_mul_moment_eq_sum
    (first second : Finset ι) (w f : ι → ℝ) :
    weightedMoment first w f * weightedMass second w -
        weightedMass first w * weightedMoment second w f =
      ∑ i ∈ first, ∑ j ∈ second, w i * w j * (f i - f j) := by
  unfold weightedMoment weightedMass
  calc
    _ = ∑ i ∈ first,
        ((∑ j ∈ second, w i * w j * f i) -
          ∑ j ∈ second, w i * w j * f j) := by
      rw [sum_mul, sum_mul]
      simp_rw [mul_sum]
      rw [← sum_sub_distrib]
      apply sum_congr rfl
      intro i hi
      congr 1
      · apply sum_congr rfl
        intro j hj
        ring
      · apply sum_congr rfl
        intro j hj
        ring
    _ = _ := by
      apply sum_congr rfl
      intro i hi
      rw [← sum_sub_distrib]
      apply sum_congr rfl
      intro j hj
      ring

variable [DecidableEq ι]

theorem weightedMoment_mul_mass_sub_mass_mul_moment_nonpos_of_cross_le
    (first second : Finset ι) (w f : ι → ℝ)
    (weight_nonneg : ∀ i ∈ first ∪ second, 0 ≤ w i)
    (cross_le : ∀ i ∈ first, ∀ j ∈ second, f i ≤ f j) :
    weightedMoment first w f * weightedMass second w -
      weightedMass first w * weightedMoment second w f ≤ 0 := by
  rw [weightedMoment_mul_mass_sub_mass_mul_moment_eq_sum]
  apply sum_nonpos
  intro i hi
  apply sum_nonpos
  intro j hj
  have weights_nonneg : 0 ≤ w i * w j :=
    mul_nonneg
      (weight_nonneg i (mem_union_left second hi))
      (weight_nonneg j (mem_union_right first hj))
  have increment_nonpos : f i - f j ≤ 0 := sub_nonpos.mpr (cross_le i hi j hj)
  exact mul_nonpos_of_nonneg_of_nonpos weights_nonneg increment_nonpos

theorem weightedMoment_mul_mass_sub_mass_mul_moment_nonneg_of_cross_ge
    (first second : Finset ι) (w f : ι → ℝ)
    (weight_nonneg : ∀ i ∈ first ∪ second, 0 ≤ w i)
    (cross_ge : ∀ i ∈ first, ∀ j ∈ second, f j ≤ f i) :
    0 ≤ weightedMoment first w f * weightedMass second w -
      weightedMass first w * weightedMoment second w f := by
  rw [weightedMoment_mul_mass_sub_mass_mul_moment_eq_sum]
  apply sum_nonneg
  intro i hi
  apply sum_nonneg
  intro j hj
  have weights_nonneg : 0 ≤ w i * w j :=
    mul_nonneg
      (weight_nonneg i (mem_union_left second hi))
      (weight_nonneg j (mem_union_right first hj))
  have increment_nonneg : 0 ≤ f i - f j := sub_nonneg.mpr (cross_ge i hi j hj)
  exact mul_nonneg weights_nonneg increment_nonneg

/-- The weighted sum of products of increments from a lower block to an upper
block.  This is the mixed coefficient of the rank-two determinant model. -/
noncomputable def weightedCrossIncrement
    (lower upper : Finset ι) (w f g : ι → ℝ) : ℝ :=
  ∑ i ∈ lower, ∑ j ∈ upper,
    w i * w j * (f j - f i) * (g j - g i)

omit [DecidableEq ι] in
/-- Expansion of the cross-block increment sum into masses and moments. -/
theorem weightedCrossIncrement_eq
    (lower upper : Finset ι) (w f g : ι → ℝ) :
    weightedCrossIncrement lower upper w f g =
      weightedMass lower w * (∑ j ∈ upper, w j * f j * g j) +
        weightedMass upper w * (∑ i ∈ lower, w i * f i * g i) -
        weightedMoment lower w f * weightedMoment upper w g -
        weightedMoment upper w f * weightedMoment lower w g := by
  have upper_diagonal_swap :
      (∑ i ∈ lower, ∑ j ∈ upper, w i * (w j * f j * g j)) =
        ∑ j ∈ upper, ∑ i ∈ lower, w i * (w j * f j * g j) := by
    rw [sum_comm]
  have lower_upper_cross_swap :
      (∑ i ∈ lower, ∑ j ∈ upper, w i * f i * (w j * g j)) =
        ∑ j ∈ upper, ∑ i ∈ lower, w i * f i * (w j * g j) := by
    rw [sum_comm]
  unfold weightedCrossIncrement weightedMass weightedMoment
  calc
    _ = ∑ i ∈ lower, ∑ j ∈ upper,
        (w i * (w j * f j * g j) + w j * (w i * f i * g i) -
          (w i * f i) * (w j * g j) - (w j * f j) * (w i * g i)) := by
      apply sum_congr rfl
      intro i hi
      apply sum_congr rfl
      intro j hj
      ring
    _ = _ := by
      simp only [sum_sub_distrib, sum_add_distrib, sum_mul, mul_sum]
      rw [upper_diagonal_swap, lower_upper_cross_swap]

omit [DecidableEq ι] in
/-- Division-free two-block law of total covariance.  The cross increment is
the sum of the two rescaled within-block covariances and the product of the
two block-mean shifts. -/
theorem weightedMass_mul_weightedMass_mul_crossIncrement
    (lower upper : Finset ι) (w f g : ι → ℝ) :
    weightedMass lower w * weightedMass upper w *
        weightedCrossIncrement lower upper w f g =
      weightedMass upper w ^ 2 *
          weightedCovarianceNumerator lower w f g +
        weightedMass lower w ^ 2 *
          weightedCovarianceNumerator upper w f g +
        (weightedMass lower w * weightedMoment upper w f -
            weightedMass upper w * weightedMoment lower w f) *
          (weightedMass lower w * weightedMoment upper w g -
            weightedMass upper w * weightedMoment lower w g) := by
  rw [weightedCrossIncrement_eq]
  unfold weightedCovarianceNumerator
  ring

/-- Ordered blocks have a nonnegative cross-increment sum. -/
theorem weightedCrossIncrement_nonneg_of_cross_le
    (lower upper : Finset ι) (w f g : ι → ℝ)
    (weight_nonneg : ∀ i ∈ lower ∪ upper, 0 ≤ w i)
    (f_cross : ∀ i ∈ lower, ∀ j ∈ upper, f i ≤ f j)
    (g_cross : ∀ i ∈ lower, ∀ j ∈ upper, g i ≤ g j) :
    0 ≤ weightedCrossIncrement lower upper w f g := by
  unfold weightedCrossIncrement
  apply sum_nonneg
  intro i hi
  apply sum_nonneg
  intro j hj
  have wi_nonneg := weight_nonneg i (mem_union_left upper hi)
  have wj_nonneg := weight_nonneg j (mem_union_right lower hj)
  exact mul_nonneg
    (mul_nonneg (mul_nonneg wi_nonneg wj_nonneg)
      (sub_nonneg.mpr (f_cross i hi j hj)))
    (sub_nonneg.mpr (g_cross i hi j hj))

/-- The cross-block mean-shift numerator is nonnegative for an ordered lower
and upper block. -/
theorem weightedMass_mul_moment_sub_mass_mul_moment_nonneg_of_cross_le
    (lower upper : Finset ι) (w f : ι → ℝ)
    (weight_nonneg : ∀ i ∈ lower ∪ upper, 0 ≤ w i)
    (cross_le : ∀ i ∈ lower, ∀ j ∈ upper, f i ≤ f j) :
    0 ≤ weightedMass lower w * weightedMoment upper w f -
      weightedMass upper w * weightedMoment lower w f := by
  have nonpos :=
    weightedMoment_mul_mass_sub_mass_mul_moment_nonpos_of_cross_le
      lower upper w f weight_nonneg cross_le
  linarith

omit [DecidableEq ι] in
/-- The sharp two-block reverse covariance inequality needs only the aggregate
product of the two block-mean shifts.  In particular, it does not require
every cross-block pair of coordinates to have the same orientation. -/
theorem four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq_of_meanShiftProduct_nonneg
    (lower upper : Finset ι) (w f g : ι → ℝ)
    (lower_mass_pos : 0 < weightedMass lower w)
    (upper_mass_pos : 0 < weightedMass upper w)
    (lower_cov_nonneg :
      0 ≤ weightedCovarianceNumerator lower w f g)
    (upper_cov_nonneg :
      0 ≤ weightedCovarianceNumerator upper w f g)
    (mean_shift_product_nonneg : 0 ≤
      (weightedMass lower w * weightedMoment upper w f -
          weightedMass upper w * weightedMoment lower w f) *
        (weightedMass lower w * weightedMoment upper w g -
          weightedMass upper w * weightedMoment lower w g)) :
    4 * weightedCovarianceNumerator lower w f g *
        weightedCovarianceNumerator upper w f g ≤
      weightedCrossIncrement lower upper w f g ^ 2 := by
  let lowerCov := weightedCovarianceNumerator lower w f g
  let upperCov := weightedCovarianceNumerator upper w f g
  let cross := weightedCrossIncrement lower upper w f g
  let lowerMass := weightedMass lower w
  let upperMass := weightedMass upper w
  have decomposition :=
    weightedMass_mul_weightedMass_mul_crossIncrement lower upper w f g
  change lowerMass * upperMass * cross =
      upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov + _ at decomposition
  have covariance_terms_nonneg :
      0 ≤ upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov :=
    add_nonneg (mul_nonneg (sq_nonneg _) lower_cov_nonneg)
      (mul_nonneg (sq_nonneg _) upper_cov_nonneg)
  have scaled_lower_bound :
      upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov ≤
        lowerMass * upperMass * cross := by
    rw [decomposition]
    exact le_add_of_nonneg_right mean_shift_product_nonneg
  have scaled_cross_nonneg : 0 ≤ lowerMass * upperMass * cross := by
    rw [decomposition]
    exact add_nonneg covariance_terms_nonneg mean_shift_product_nonneg
  have squared_lower_bound :
      (upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov) ^ 2 ≤
        (lowerMass * upperMass * cross) ^ 2 :=
    (sq_le_sq₀ covariance_terms_nonneg scaled_cross_nonneg).2 scaled_lower_bound
  have am_gm :
      4 * (lowerMass * upperMass) ^ 2 * lowerCov * upperCov ≤
        (upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov) ^ 2 := by
    nlinarith [sq_nonneg (upperMass ^ 2 * lowerCov -
      lowerMass ^ 2 * upperCov)]
  have mass_product_sq_pos : 0 < (lowerMass * upperMass) ^ 2 :=
    sq_pos_of_pos (mul_pos lower_mass_pos upper_mass_pos)
  apply le_of_mul_le_mul_left _ mass_product_sq_pos
  calc
    (lowerMass * upperMass) ^ 2 *
        (4 * weightedCovarianceNumerator lower w f g *
          weightedCovarianceNumerator upper w f g) =
      4 * (lowerMass * upperMass) ^ 2 * lowerCov * upperCov := by
        simp only [lowerCov, upperCov]
        ring
    _ ≤ (upperMass ^ 2 * lowerCov + lowerMass ^ 2 * upperCov) ^ 2 := am_gm
    _ ≤ (lowerMass * upperMass * cross) ^ 2 := squared_lower_bound
    _ = (lowerMass * upperMass) ^ 2 *
        weightedCrossIncrement lower upper w f g ^ 2 := by
      simp only [cross]
      ring

/-- Sharp two-block reverse covariance inequality.  It is the
fresh-order-two Lorentzian discriminant inequality for pair weights generated
by two comonotone ratio sequences. -/
theorem four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq
    (lower upper : Finset ι) (w f g : ι → ℝ)
    (lower_mass_pos : 0 < weightedMass lower w)
    (upper_mass_pos : 0 < weightedMass upper w)
    (weight_nonneg : ∀ i ∈ lower ∪ upper, 0 ≤ w i)
    (lower_monovary : MonovaryOn f g lower)
    (upper_monovary : MonovaryOn f g upper)
    (f_cross : ∀ i ∈ lower, ∀ j ∈ upper, f i ≤ f j)
    (g_cross : ∀ i ∈ lower, ∀ j ∈ upper, g i ≤ g j) :
    4 * weightedCovarianceNumerator lower w f g *
        weightedCovarianceNumerator upper w f g ≤
      weightedCrossIncrement lower upper w f g ^ 2 := by
  apply
    four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq_of_meanShiftProduct_nonneg
      lower upper w f g lower_mass_pos upper_mass_pos
  · apply weightedCovarianceNumerator_nonneg lower w f g
    · intro i hi
      exact weight_nonneg i (mem_union_left upper hi)
    · exact lower_monovary
  · apply weightedCovarianceNumerator_nonneg upper w f g
    · intro i hi
      exact weight_nonneg i (mem_union_right lower hi)
    · exact upper_monovary
  · apply mul_nonneg
    · exact weightedMass_mul_moment_sub_mass_mul_moment_nonneg_of_cross_le
        lower upper w f weight_nonneg f_cross
    · exact weightedMass_mul_moment_sub_mass_mul_moment_nonneg_of_cross_le
        lower upper w g weight_nonneg g_cross

theorem weightedBlockMeanShift_eq
    (first second : Finset ι) (w f : ι → ℝ)
    (disjoint : Disjoint first second) :
    weightedMoment first w f * weightedMass (first ∪ second) w -
        weightedMass first w * weightedMoment (first ∪ second) w f =
      weightedMoment first w f * weightedMass second w -
        weightedMass first w * weightedMoment second w f := by
  simp only [weightedMass, weightedMoment, sum_union disjoint]
  ring

/-- A lower ordered block has a nonnegative centered product sum relative to
the union of the lower and upper blocks. -/
theorem weightedCenteredLowerGroup_nonneg
    (lower upper : Finset ι) (w f g : ι → ℝ)
    (disjoint : Disjoint lower upper)
    (lower_mass_pos : 0 < weightedMass lower w)
    (weight_nonneg : ∀ i ∈ lower ∪ upper, 0 ≤ w i)
    (monovary : MonovaryOn f g lower)
    (f_cross : ∀ i ∈ lower, ∀ j ∈ upper, f i ≤ f j)
    (g_cross : ∀ i ∈ lower, ∀ j ∈ upper, g i ≤ g j) :
    0 ≤ weightedCenteredGroup (lower ∪ upper) lower w f g := by
  have f_shift_nonpos :
      weightedMoment lower w f * weightedMass (lower ∪ upper) w -
          weightedMass lower w * weightedMoment (lower ∪ upper) w f ≤ 0 := by
    rw [weightedBlockMeanShift_eq lower upper w f disjoint]
    exact weightedMoment_mul_mass_sub_mass_mul_moment_nonpos_of_cross_le
      lower upper w f weight_nonneg f_cross
  have g_shift_nonpos :
      weightedMoment lower w g * weightedMass (lower ∪ upper) w -
          weightedMass lower w * weightedMoment (lower ∪ upper) w g ≤ 0 := by
    rw [weightedBlockMeanShift_eq lower upper w g disjoint]
    exact weightedMoment_mul_mass_sub_mass_mul_moment_nonpos_of_cross_le
      lower upper w g weight_nonneg g_cross
  apply weightedCenteredGroup_nonneg (lower ∪ upper) lower w f g
    lower_mass_pos
  · intro i hi
    exact weight_nonneg i (mem_union_left upper hi)
  · exact monovary
  · exact mul_nonneg_of_nonpos_of_nonpos f_shift_nonpos g_shift_nonpos

/-- An upper ordered block satisfies the symmetric centered-product sign. -/
theorem weightedCenteredUpperGroup_nonneg
    (lower upper : Finset ι) (w f g : ι → ℝ)
    (disjoint : Disjoint lower upper)
    (upper_mass_pos : 0 < weightedMass upper w)
    (weight_nonneg : ∀ i ∈ lower ∪ upper, 0 ≤ w i)
    (monovary : MonovaryOn f g upper)
    (f_cross : ∀ i ∈ lower, ∀ j ∈ upper, f i ≤ f j)
    (g_cross : ∀ i ∈ lower, ∀ j ∈ upper, g i ≤ g j) :
    0 ≤ weightedCenteredGroup (lower ∪ upper) upper w f g := by
  have f_shift_nonneg : 0 ≤
      weightedMoment upper w f * weightedMass (lower ∪ upper) w -
        weightedMass upper w * weightedMoment (lower ∪ upper) w f := by
    rw [union_comm, weightedBlockMeanShift_eq upper lower w f disjoint.symm]
    exact weightedMoment_mul_mass_sub_mass_mul_moment_nonneg_of_cross_ge
      upper lower w f (by simpa [union_comm] using weight_nonneg)
        (fun i hi j hj ↦ f_cross j hj i hi)
  have g_shift_nonneg : 0 ≤
      weightedMoment upper w g * weightedMass (lower ∪ upper) w -
        weightedMass upper w * weightedMoment (lower ∪ upper) w g := by
    rw [union_comm, weightedBlockMeanShift_eq upper lower w g disjoint.symm]
    exact weightedMoment_mul_mass_sub_mass_mul_moment_nonneg_of_cross_ge
      upper lower w g (by simpa [union_comm] using weight_nonneg)
        (fun i hi j hj ↦ g_cross j hj i hi)
  apply weightedCenteredGroup_nonneg (lower ∪ upper) upper w f g
    upper_mass_pos
  · intro i hi
    exact weight_nonneg i (mem_union_right lower hi)
  · exact monovary
  · exact mul_nonneg f_shift_nonneg g_shift_nonneg

/-- Degrees below a cutoff. -/
def degreePrefix (n cutoff : ℕ) : Finset (Fin n) :=
  Finset.univ.filter fun k ↦ k.1 < cutoff

/-- Degrees at or above a cutoff. -/
def degreeSuffix (n cutoff : ℕ) : Finset (Fin n) :=
  Finset.univ.filter fun k ↦ cutoff ≤ k.1

theorem degreePrefix_disjoint_degreeSuffix (n cutoff : ℕ) :
    Disjoint (degreePrefix n cutoff) (degreeSuffix n cutoff) := by
  rw [Finset.disjoint_left]
  intro k hk_low hk_high
  simp only [degreePrefix, mem_filter, mem_univ, true_and] at hk_low
  simp only [degreeSuffix, mem_filter, mem_univ, true_and] at hk_high
  omega

theorem degreePrefix_union_degreeSuffix (n cutoff : ℕ) :
    degreePrefix n cutoff ∪ degreeSuffix n cutoff = Finset.univ := by
  ext k
  by_cases hk : k.1 < cutoff
  · simp [degreePrefix, degreeSuffix, hk]
  · have hk' : cutoff ≤ k.1 := Nat.le_of_not_gt hk
    simp [degreePrefix, degreeSuffix, hk, hk']

omit [DecidableEq ι] in theorem weightedMass_pos_of_pos
    (s : Finset ι) (w : ι → ℝ) (nonempty : s.Nonempty)
    (weight_pos : ∀ i ∈ s, 0 < w i) :
    0 < weightedMass s w := by
  unfold weightedMass
  exact sum_pos weight_pos nonempty

/-- A weighted product centered at the global geometric moments. -/
noncomputable def geometricCenteredProduct {n : ℕ}
    (w : Fin n → ℝ) (rowRatio columnRatio : ℝ) (degree : Fin n) : ℝ :=
  w degree *
    (weightedMass Finset.univ w * rowRatio ^ degree.1 -
      weightedMoment Finset.univ w fun k ↦ rowRatio ^ k.1) *
    (weightedMass Finset.univ w * columnRatio ^ degree.1 -
      weightedMoment Finset.univ w fun k ↦ columnRatio ^ k.1)

theorem sum_geometricCenteredProduct_eq_weightedCenteredGroup
    {n : ℕ} (group : Finset (Fin n))
    (w : Fin n → ℝ) (rowRatio columnRatio : ℝ) :
    (∑ k ∈ group, geometricCenteredProduct w rowRatio columnRatio k) =
      weightedCenteredGroup Finset.univ group w
        (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1) := by
  rfl

theorem geometricCenteredProduct_prefix_nonneg
    {n cutoff : ℕ} (n_pos : 0 < n) (cutoff_pos : 0 < cutoff)
    (w : Fin n → ℝ) (weight_pos : ∀ k, 0 < w k)
    {rowRatio columnRatio : ℝ}
    (rowRatio_ge_one : 1 ≤ rowRatio)
    (columnRatio_ge_one : 1 ≤ columnRatio) :
    0 ≤ ∑ k ∈ degreePrefix n cutoff,
      geometricCenteredProduct w rowRatio columnRatio k := by
  let lower := degreePrefix n cutoff
  let upper := degreeSuffix n cutoff
  have lower_nonempty : lower.Nonempty := by
    refine ⟨⟨0, n_pos⟩, ?_⟩
    simp [lower, degreePrefix, cutoff_pos]
  have lower_mass_pos : 0 < weightedMass lower w :=
    weightedMass_pos_of_pos lower w lower_nonempty fun i hi ↦ weight_pos i
  have row_monotone : Monotone (fun k : Fin n ↦ rowRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ rowRatio_ge_one (Fin.mk_le_mk.mp hij)
  have column_monotone : Monotone (fun k : Fin n ↦ columnRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ columnRatio_ge_one (Fin.mk_le_mk.mp hij)
  have cross_index : ∀ i ∈ lower, ∀ j ∈ upper, i ≤ j := by
    intro i hi j hj
    simp only [lower, degreePrefix, mem_filter, mem_univ, true_and] at hi
    simp only [upper, degreeSuffix, mem_filter, mem_univ, true_and] at hj
    exact Fin.mk_le_mk.mpr (Nat.le_trans (Nat.le_of_lt hi) hj)
  rw [sum_geometricCenteredProduct_eq_weightedCenteredGroup]
  rw [← degreePrefix_union_degreeSuffix n cutoff]
  apply weightedCenteredLowerGroup_nonneg lower upper w
    (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1)
    (degreePrefix_disjoint_degreeSuffix n cutoff) lower_mass_pos
  · intro k hk
    exact (weight_pos k).le
  · exact row_monotone.monovary column_monotone |>.monovaryOn lower
  · intro i hi j hj
    exact row_monotone (cross_index i hi j hj)
  · intro i hi j hj
    exact column_monotone (cross_index i hi j hj)

theorem geometricCenteredProduct_suffix_nonneg
    {n cutoff : ℕ} (cutoff_lt : cutoff < n)
    (w : Fin n → ℝ) (weight_pos : ∀ k, 0 < w k)
    {rowRatio columnRatio : ℝ}
    (rowRatio_ge_one : 1 ≤ rowRatio)
    (columnRatio_ge_one : 1 ≤ columnRatio) :
    0 ≤ ∑ k ∈ degreeSuffix n cutoff,
      geometricCenteredProduct w rowRatio columnRatio k := by
  let lower := degreePrefix n cutoff
  let upper := degreeSuffix n cutoff
  have upper_nonempty : upper.Nonempty := by
    refine ⟨⟨cutoff, cutoff_lt⟩, ?_⟩
    simp [upper, degreeSuffix]
  have upper_mass_pos : 0 < weightedMass upper w :=
    weightedMass_pos_of_pos upper w upper_nonempty fun i hi ↦ weight_pos i
  have row_monotone : Monotone (fun k : Fin n ↦ rowRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ rowRatio_ge_one (Fin.mk_le_mk.mp hij)
  have column_monotone : Monotone (fun k : Fin n ↦ columnRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ columnRatio_ge_one (Fin.mk_le_mk.mp hij)
  have cross_index : ∀ i ∈ lower, ∀ j ∈ upper, i ≤ j := by
    intro i hi j hj
    simp only [lower, degreePrefix, mem_filter, mem_univ, true_and] at hi
    simp only [upper, degreeSuffix, mem_filter, mem_univ, true_and] at hj
    exact Fin.mk_le_mk.mpr (Nat.le_trans (Nat.le_of_lt hi) hj)
  rw [sum_geometricCenteredProduct_eq_weightedCenteredGroup]
  rw [← degreePrefix_union_degreeSuffix n cutoff]
  apply weightedCenteredUpperGroup_nonneg lower upper w
    (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1)
    (degreePrefix_disjoint_degreeSuffix n cutoff) upper_mass_pos
  · intro k hk
    exact (weight_pos k).le
  · exact row_monotone.monovary column_monotone |>.monovaryOn upper
  · intro i hi j hj
    exact row_monotone (cross_index i hi j hj)
  · intro i hi j hj
    exact column_monotone (cross_index i hi j hj)

/-- The degree-prefix/suffix split of two geometric ratio coordinates obeys
the sharp Lorentzian discriminant inequality.  This is the ratio-coordinate
form consumed by conditioned rank-two feature factorizations. -/
theorem geometricCrossIncrement_lorentzian
    {n cutoff : ℕ} (cutoff_pos : 0 < cutoff) (cutoff_lt : cutoff < n)
    (w : Fin n → ℝ) (weight_pos : ∀ k, 0 < w k)
    {rowRatio columnRatio : ℝ}
    (rowRatio_ge_one : 1 ≤ rowRatio)
    (columnRatio_ge_one : 1 ≤ columnRatio) :
    4 * weightedCovarianceNumerator (degreePrefix n cutoff) w
          (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1) *
        weightedCovarianceNumerator (degreeSuffix n cutoff) w
          (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1) ≤
      weightedCrossIncrement (degreePrefix n cutoff) (degreeSuffix n cutoff) w
          (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1) ^ 2 := by
  let lower := degreePrefix n cutoff
  let upper := degreeSuffix n cutoff
  have n_pos : 0 < n := lt_trans cutoff_pos cutoff_lt
  have lower_nonempty : lower.Nonempty := by
    refine ⟨⟨0, n_pos⟩, ?_⟩
    simp [lower, degreePrefix, cutoff_pos]
  have upper_nonempty : upper.Nonempty := by
    refine ⟨⟨cutoff, cutoff_lt⟩, ?_⟩
    simp [upper, degreeSuffix]
  have lower_mass_pos : 0 < weightedMass lower w :=
    weightedMass_pos_of_pos lower w lower_nonempty fun i hi ↦ weight_pos i
  have upper_mass_pos : 0 < weightedMass upper w :=
    weightedMass_pos_of_pos upper w upper_nonempty fun i hi ↦ weight_pos i
  have row_monotone : Monotone (fun k : Fin n ↦ rowRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ rowRatio_ge_one (Fin.mk_le_mk.mp hij)
  have column_monotone : Monotone (fun k : Fin n ↦ columnRatio ^ k.1) := by
    intro i j hij
    exact pow_le_pow_right₀ columnRatio_ge_one (Fin.mk_le_mk.mp hij)
  have cross_index : ∀ i ∈ lower, ∀ j ∈ upper, i ≤ j := by
    intro i hi j hj
    simp only [lower, degreePrefix, mem_filter, mem_univ, true_and] at hi
    simp only [upper, degreeSuffix, mem_filter, mem_univ, true_and] at hj
    exact Fin.mk_le_mk.mpr (Nat.le_trans (Nat.le_of_lt hi) hj)
  apply four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq
    lower upper w (fun k ↦ rowRatio ^ k.1) (fun k ↦ columnRatio ^ k.1)
    lower_mass_pos upper_mass_pos
  · intro k hk
    exact (weight_pos k).le
  · exact row_monotone.monovary column_monotone |>.monovaryOn lower
  · exact row_monotone.monovary column_monotone |>.monovaryOn upper
  · intro i hi j hj
    exact row_monotone (cross_index i hi j hj)
  · intro i hi j hj
    exact column_monotone (cross_index i hi j hj)

end Matrix
