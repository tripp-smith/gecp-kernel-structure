import GECPKernelStructure.Fermionic.Centering
import GECPKernelStructure.Fermionic.GECPBounds
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

/-- The second-order remainder of the real exponential at zero. -/
noncomputable def expRemainder (z : ℝ) : ℝ :=
  Real.exp z - 1 - z

/-- On the unit interval, the exponential remainder lies between zero and `z²`. -/
theorem expRemainder_mem_Icc {z : ℝ} (hz : abs z ≤ 1) :
    expRemainder z ∈ Set.Icc 0 (z ^ 2) := by
  constructor
  · unfold expRemainder
    linarith [Real.add_one_le_exp z]
  · exact (le_abs_self (expRemainder z)).trans (by
      simpa only [expRemainder] using Real.abs_exp_sub_one_sub_id_le hz)

/-- Linear interpolation of `exp (-x ω)` between the frequencies `-Λ` and `Λ`. -/
noncomputable def centeredExpSecant (Λ x ω : ℝ) : ℝ :=
  ((Λ - ω) / (2 * Λ)) * Real.exp (x * Λ) +
    ((Λ + ω) / (2 * Λ)) * Real.exp (-x * Λ)

/-- The rank-two centered secant approximation to the fermionic kernel. -/
noncomputable def twoCornerApproximation (Λ t ω : ℝ) : ℝ :=
  centeredExpSecant Λ (t - (1 / 2 : ℝ)) ω /
    (2 * Real.cosh (ω / 2))

private theorem abs_sub_le_of_Icc {a b M : ℝ}
    (ha : a ∈ Set.Icc 0 M) (hb : b ∈ Set.Icc 0 M) :
    abs (a - b) ≤ M := by
  rw [abs_le]
  constructor <;> linarith [ha.1, ha.2, hb.1, hb.2]

/-- The centered exponential secant has a uniform quadratic error on the base scale. -/
theorem centeredExpSecant_error_le_quarter {Λ x ω : ℝ}
    (cutoff_pos : 0 < Λ) (cutoff_le_one : Λ ≤ 1)
    (centered_time : abs x ≤ (1 / 2 : ℝ))
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (Real.exp (-x * ω) - centeredExpSecant Λ x ω) ≤ 1 / 4 := by
  let a : ℝ := (Λ - ω) / (2 * Λ)
  let b : ℝ := (Λ + ω) / (2 * Λ)
  let M : ℝ := (x * Λ) ^ 2
  have cutoff_nonneg : 0 ≤ Λ := cutoff_pos.le
  have denominator_pos : 0 < 2 * Λ := by positivity
  have a_nonneg : 0 ≤ a := by
    dsimp [a]
    exact div_nonneg (sub_nonneg.mpr frequency_upper) denominator_pos.le
  have b_nonneg : 0 ≤ b := by
    dsimp [b]
    exact div_nonneg (by linarith) denominator_pos.le
  have weights_sum : a + b = 1 := by
    dsimp [a, b]
    field_simp
    ring
  have frequency_abs : abs ω ≤ Λ := by
    rw [abs_le]
    exact ⟨frequency_lower, frequency_upper⟩
  have centered_cutoff_abs : abs (x * Λ) ≤ 1 := by
    rw [abs_mul, abs_of_nonneg cutoff_nonneg]
    calc
      abs x * Λ ≤ (1 / 2 : ℝ) * 1 :=
        mul_le_mul centered_time cutoff_le_one cutoff_nonneg (by norm_num)
      _ ≤ 1 := by norm_num
  have centered_frequency_abs : abs (-x * ω) ≤ 1 := by
    rw [abs_mul, abs_neg]
    calc
      abs x * abs ω ≤ (1 / 2 : ℝ) * Λ :=
        mul_le_mul centered_time frequency_abs (abs_nonneg _) (by norm_num)
      _ ≤ (1 / 2 : ℝ) * 1 := mul_le_mul_of_nonneg_left cutoff_le_one (by norm_num)
      _ ≤ 1 := by norm_num
  have remainder_pos := expRemainder_mem_Icc centered_cutoff_abs
  have negative_centered_cutoff_abs : abs (-x * Λ) ≤ 1 := by
    simpa only [abs_mul, abs_neg] using centered_cutoff_abs
  have remainder_neg : expRemainder (-x * Λ) ∈ Set.Icc 0 M := by
    have raw := expRemainder_mem_Icc negative_centered_cutoff_abs
    refine ⟨raw.1, raw.2.trans_eq ?_⟩
    dsimp [M]
    ring
  have remainder_frequency := expRemainder_mem_Icc centered_frequency_abs
  have frequency_square_le : (-x * ω) ^ 2 ≤ M := by
    have product_abs_le : abs (-x * ω) ≤ abs (x * Λ) := by
      rw [abs_mul, abs_neg, abs_mul, abs_of_nonneg cutoff_nonneg]
      exact mul_le_mul_of_nonneg_left frequency_abs (abs_nonneg x)
    have h := mul_self_le_mul_self (abs_nonneg _) product_abs_le
    calc
      (-x * ω) ^ 2 = abs (-x * ω) ^ 2 := (sq_abs (-x * ω)).symm
      _ ≤ abs (x * Λ) ^ 2 := by simpa only [pow_two] using h
      _ = (x * Λ) ^ 2 := sq_abs (x * Λ)
      _ = M := by rfl
  have frequency_remainder_mem : expRemainder (-x * ω) ∈ Set.Icc 0 M :=
    ⟨remainder_frequency.1, remainder_frequency.2.trans frequency_square_le⟩
  have weighted_remainder_mem :
      a * expRemainder (x * Λ) + b * expRemainder (-x * Λ) ∈ Set.Icc 0 M := by
    have remainder_pos_upper : expRemainder (x * Λ) ≤ M := by
      simpa only [M] using remainder_pos.2
    have remainder_neg_upper : expRemainder (-x * Λ) ≤ M := remainder_neg.2
    constructor
    · exact add_nonneg (mul_nonneg a_nonneg remainder_pos.1)
        (mul_nonneg b_nonneg remainder_neg.1)
    · calc
        a * expRemainder (x * Λ) + b * expRemainder (-x * Λ) ≤
            a * M + b * M := by
              exact add_le_add (mul_le_mul_of_nonneg_left remainder_pos_upper a_nonneg)
                (mul_le_mul_of_nonneg_left remainder_neg_upper b_nonneg)
        _ = M := by rw [← add_mul, weights_sum, one_mul]
  have remainder_difference :
      Real.exp (-x * ω) - centeredExpSecant Λ x ω =
        expRemainder (-x * ω) -
          (a * expRemainder (x * Λ) + b * expRemainder (-x * Λ)) := by
    unfold centeredExpSecant expRemainder
    dsimp only [a, b]
    field_simp [cutoff_pos.ne']
    ring
  rw [remainder_difference]
  exact (abs_sub_le_of_Icc frequency_remainder_mem weighted_remainder_mem).trans (by
    dsimp [M]
    have product_abs : abs (x * Λ) ≤ (1 / 2 : ℝ) := by
      rw [abs_mul, abs_of_nonneg cutoff_nonneg]
      exact (mul_le_mul centered_time cutoff_le_one cutoff_nonneg (by norm_num)).trans_eq
        (by norm_num)
    have h := mul_self_le_mul_self (abs_nonneg (x * Λ)) product_abs
    nlinarith [h, sq_abs (x * Λ)])

/-- The explicit rank-two approximation has fermionic-kernel error at most `1/8`. -/
theorem twoCornerApproximation_error_le_eighth {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ) (cutoff_le_one : Λ ≤ 1)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (fermionicKernel t ω - twoCornerApproximation Λ t ω) ≤ 1 / 8 := by
  have centered_time : abs (t - (1 / 2 : ℝ)) ≤ 1 / 2 := by
    rw [abs_le]
    constructor <;> linarith
  rw [fermionicKernel_centered]
  unfold centeredKernel twoCornerApproximation
  rw [← sub_div]
  rw [abs_div]
  have numerator_bound := centeredExpSecant_error_le_quarter cutoff_pos cutoff_le_one
    centered_time frequency_lower frequency_upper
  have denominator_lower : (2 : ℝ) ≤ abs (2 * Real.cosh (ω / 2)) := by
    rw [abs_of_pos (by positivity)]
    nlinarith [Real.one_le_cosh (ω / 2)]
  calc
    abs (Real.exp (-(t - 1 / 2) * ω) -
        centeredExpSecant Λ (t - 1 / 2) ω) /
        abs (2 * Real.cosh (ω / 2)) ≤ (1 / 4 : ℝ) / 2 := by
          exact div_le_div₀ (by norm_num) numerator_bound (by norm_num) denominator_lower
    _ = 1 / 8 := by norm_num

/-- Left endpoint cardinal function in the two cutoff-column span. -/
noncomputable def cornerLeftWeight (Λ t : ℝ) : ℝ :=
  Real.sinh (Λ * (1 - t)) / Real.sinh Λ

/-- Right endpoint cardinal function in the two cutoff-column span. -/
noncomputable def cornerRightWeight (Λ t : ℝ) : ℝ :=
  Real.sinh (Λ * t) / Real.sinh Λ

/-- Interpolation at the time endpoints using the cutoff-column cardinal functions. -/
noncomputable def twoCornerInterpolate (Λ : ℝ) (g : ℝ → ℝ → ℝ) (t ω : ℝ) : ℝ :=
  cornerLeftWeight Λ t * g 0 ω + cornerRightWeight Λ t * g 1 ω

/-- The two endpoint cardinal weights are nonnegative and have sum at most one. -/
theorem cornerWeights_nonnegative_sum_le_one {Λ t : ℝ}
    (cutoff_pos : 0 < Λ) (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1) :
    0 ≤ cornerLeftWeight Λ t ∧ 0 ≤ cornerRightWeight Λ t ∧
      cornerLeftWeight Λ t + cornerRightWeight Λ t ≤ 1 := by
  have sinh_cutoff_pos : 0 < Real.sinh Λ := Real.sinh_pos_iff.mpr cutoff_pos
  have left_argument_nonneg : 0 ≤ Λ * (1 - t) :=
    mul_nonneg cutoff_pos.le (sub_nonneg.mpr time_le_one)
  have right_argument_nonneg : 0 ≤ Λ * t := mul_nonneg cutoff_pos.le time_nonneg
  have left_sinh_nonneg : 0 ≤ Real.sinh (Λ * (1 - t)) :=
    Real.sinh_nonneg_iff.mpr left_argument_nonneg
  have right_sinh_nonneg : 0 ≤ Real.sinh (Λ * t) :=
    Real.sinh_nonneg_iff.mpr right_argument_nonneg
  have left_nonneg : 0 ≤ cornerLeftWeight Λ t := by
    exact div_nonneg left_sinh_nonneg sinh_cutoff_pos.le
  have right_nonneg : 0 ≤ cornerRightWeight Λ t := by
    exact div_nonneg right_sinh_nonneg sinh_cutoff_pos.le
  refine ⟨left_nonneg, right_nonneg, ?_⟩
  rw [cornerLeftWeight, cornerRightWeight, ← add_div]
  apply (div_le_one sinh_cutoff_pos).mpr
  have left_cosh_bound :
      Real.sinh (Λ * (1 - t)) ≤
        Real.sinh (Λ * (1 - t)) * Real.cosh (Λ * t) := by
    nlinarith [mul_le_mul_of_nonneg_left (Real.one_le_cosh (Λ * t)) left_sinh_nonneg]
  have right_cosh_bound :
      Real.sinh (Λ * t) ≤
        Real.cosh (Λ * (1 - t)) * Real.sinh (Λ * t) := by
    nlinarith [mul_le_mul_of_nonneg_right
      (Real.one_le_cosh (Λ * (1 - t))) right_sinh_nonneg]
  calc
    Real.sinh (Λ * (1 - t)) + Real.sinh (Λ * t) ≤
        Real.sinh (Λ * (1 - t)) * Real.cosh (Λ * t) +
          Real.cosh (Λ * (1 - t)) * Real.sinh (Λ * t) :=
            add_le_add left_cosh_bound right_cosh_bound
    _ = Real.sinh Λ := by
      rw [← Real.sinh_add]
      congr 1
      ring

private theorem cornerWeights_reconstruct_exp {Λ t : ℝ} (cutoff_pos : 0 < Λ) :
    cornerLeftWeight Λ t * Real.exp (-Λ / 2) +
        cornerRightWeight Λ t * Real.exp (Λ / 2) =
      Real.exp ((t - 1 / 2) * Λ) := by
  have sinh_ne : Real.sinh Λ ≠ 0 := (Real.sinh_pos_iff.mpr cutoff_pos).ne'
  have h₁ : Real.exp (Λ * (1 - t)) * Real.exp (-Λ / 2) =
      Real.exp (Λ / 2 - Λ * t) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₂ : Real.exp (-(Λ * (1 - t))) * Real.exp (-Λ / 2) =
      Real.exp (Λ * t - 3 * Λ / 2) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₃ : Real.exp (Λ * t) * Real.exp (Λ / 2) =
      Real.exp (Λ * t + Λ / 2) := by rw [← Real.exp_add]
  have h₄ : Real.exp (-(Λ * t)) * Real.exp (Λ / 2) =
      Real.exp (Λ / 2 - Λ * t) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₅ : Real.exp Λ * Real.exp ((t - 1 / 2) * Λ) =
      Real.exp (Λ * t + Λ / 2) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₆ : Real.exp (-Λ) * Real.exp ((t - 1 / 2) * Λ) =
      Real.exp (Λ * t - 3 * Λ / 2) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have h₅' : Real.exp ((t - 1 / 2) * Λ) * Real.exp Λ =
      Real.exp (Λ * t + Λ / 2) := by simpa [mul_comm] using h₅
  have h₆' : Real.exp ((t - 1 / 2) * Λ) * Real.exp (-Λ) =
      Real.exp (Λ * t - 3 * Λ / 2) := by simpa [mul_comm] using h₆
  unfold cornerLeftWeight cornerRightWeight
  rw [div_mul_eq_mul_div, div_mul_eq_mul_div, ← add_div]
  apply (div_eq_iff sinh_ne).mpr
  rw [Real.sinh_eq, Real.sinh_eq, Real.sinh_eq]
  field_simp
  have rhs_exp : Real.exp (Λ * (t * 2 - 1) / 2) =
      Real.exp ((t - 1 / 2) * Λ) := by
    congr 1
    ring
  rw [rhs_exp]
  rw [sub_mul, sub_mul]
  rw [show -(Λ / 2) = -Λ / 2 by ring]
  rw [h₁, h₂, h₃, h₄]
  rw [mul_sub, h₅', h₆']
  ring

private theorem cornerWeights_reconstruct_exp_neg {Λ t : ℝ} (cutoff_pos : 0 < Λ) :
    cornerLeftWeight Λ t * Real.exp (Λ / 2) +
        cornerRightWeight Λ t * Real.exp (-Λ / 2) =
      Real.exp (-(t - 1 / 2) * Λ) := by
  have reflected := cornerWeights_reconstruct_exp (Λ := Λ) (t := 1 - t) cutoff_pos
  have left_reflection : cornerLeftWeight Λ (1 - t) = cornerRightWeight Λ t := by
    unfold cornerLeftWeight cornerRightWeight
    congr 2
    ring
  have right_reflection : cornerRightWeight Λ (1 - t) = cornerLeftWeight Λ t := by
    unfold cornerLeftWeight cornerRightWeight
    congr 2
  rw [left_reflection, right_reflection] at reflected
  rw [add_comm] at reflected
  calc
    cornerLeftWeight Λ t * Real.exp (Λ / 2) +
        cornerRightWeight Λ t * Real.exp (-Λ / 2) =
      Real.exp ((1 - t - 1 / 2) * Λ) := reflected
    _ = Real.exp (-(t - 1 / 2) * Λ) := by
      congr 1
      ring

/-- The endpoint operator fixes the explicit rank-two secant approximation. -/
theorem twoCornerInterpolate_approximation {Λ t ω : ℝ} (cutoff_pos : 0 < Λ) :
    twoCornerInterpolate Λ (twoCornerApproximation Λ) t ω =
      twoCornerApproximation Λ t ω := by
  have positive_reconstruction := cornerWeights_reconstruct_exp (t := t) cutoff_pos
  have negative_reconstruction := cornerWeights_reconstruct_exp_neg (t := t) cutoff_pos
  unfold twoCornerInterpolate twoCornerApproximation centeredExpSecant
  rw [show (0 - (1 / 2 : ℝ)) * Λ = -Λ / 2 by ring,
    show -(0 - (1 / 2 : ℝ)) * Λ = Λ / 2 by ring,
    show (1 - (1 / 2 : ℝ)) * Λ = Λ / 2 by ring,
    show -(1 - (1 / 2 : ℝ)) * Λ = -Λ / 2 by ring,
    ← positive_reconstruction, ← negative_reconstruction]
  ring

/-- Endpoint interpolation amplifies a uniform approximation error by at most one. -/
theorem twoCornerInterpolate_error_le {Λ t ω δ : ℝ}
    (cutoff_pos : 0 < Λ) (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    {g a : ℝ → ℝ → ℝ}
    (error : ∀ s, s = 0 ∨ s = 1 → abs (g s ω - a s ω) ≤ δ) :
    abs (twoCornerInterpolate Λ g t ω - twoCornerInterpolate Λ a t ω) ≤ δ := by
  rcases cornerWeights_nonnegative_sum_le_one cutoff_pos time_nonneg time_le_one with
    ⟨left_nonneg, right_nonneg, weight_sum⟩
  have δ_nonneg : 0 ≤ δ := (abs_nonneg _).trans (error 0 (Or.inl rfl))
  unfold twoCornerInterpolate
  rw [show cornerLeftWeight Λ t * g 0 ω + cornerRightWeight Λ t * g 1 ω -
      (cornerLeftWeight Λ t * a 0 ω + cornerRightWeight Λ t * a 1 ω) =
      cornerLeftWeight Λ t * (g 0 ω - a 0 ω) +
        cornerRightWeight Λ t * (g 1 ω - a 1 ω) by ring]
  calc
    abs (cornerLeftWeight Λ t * (g 0 ω - a 0 ω) +
        cornerRightWeight Λ t * (g 1 ω - a 1 ω)) ≤
        cornerLeftWeight Λ t * abs (g 0 ω - a 0 ω) +
          cornerRightWeight Λ t * abs (g 1 ω - a 1 ω) := by
            calc
              _ ≤ abs (cornerLeftWeight Λ t * (g 0 ω - a 0 ω)) +
                  abs (cornerRightWeight Λ t * (g 1 ω - a 1 ω)) := abs_add_le _ _
              _ = _ := by
                simp only [abs_mul, abs_of_nonneg left_nonneg,
                  abs_of_nonneg right_nonneg]
    _ ≤ cornerLeftWeight Λ t * δ + cornerRightWeight Λ t * δ := by
      exact add_le_add (mul_le_mul_of_nonneg_left (error 0 (Or.inl rfl)) left_nonneg)
        (mul_le_mul_of_nonneg_left (error 1 (Or.inr rfl)) right_nonneg)
    _ = (cornerLeftWeight Λ t + cornerRightWeight Λ t) * δ := by ring
    _ ≤ δ := by nlinarith

/-- The endpoint cross approximation is uniformly within `1/4` of the kernel. -/
theorem fermionicKernel_twoCornerInterpolate_error_le_quarter {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ) (cutoff_le_one : Λ ≤ 1)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (fermionicKernel t ω -
      twoCornerInterpolate Λ fermionicKernel t ω) ≤ 1 / 4 := by
  have point_error := twoCornerApproximation_error_le_eighth cutoff_pos cutoff_le_one
    time_nonneg time_le_one frequency_lower frequency_upper
  have endpoint_error : ∀ s : ℝ, s = 0 ∨ s = 1 →
      abs (fermionicKernel s ω - twoCornerApproximation Λ s ω) ≤ 1 / 8 := by
    intro s hs
    rcases hs with rfl | rfl
    · exact twoCornerApproximation_error_le_eighth cutoff_pos cutoff_le_one
        (by norm_num) (by norm_num) frequency_lower frequency_upper
    · exact twoCornerApproximation_error_le_eighth cutoff_pos cutoff_le_one
        (by norm_num) (by norm_num) frequency_lower frequency_upper
  have interpolated_error := twoCornerInterpolate_error_le cutoff_pos time_nonneg
    time_le_one endpoint_error
  have approximation_fixed := twoCornerInterpolate_approximation
    (Λ := Λ) (t := t) (ω := ω) cutoff_pos
  calc
    abs (fermionicKernel t ω - twoCornerInterpolate Λ fermionicKernel t ω) =
        abs ((fermionicKernel t ω - twoCornerApproximation Λ t ω) +
          (twoCornerInterpolate Λ (twoCornerApproximation Λ) t ω -
            twoCornerInterpolate Λ fermionicKernel t ω)) := by
              rw [approximation_fixed]
              congr 1
              ring
    _ ≤ abs (fermionicKernel t ω - twoCornerApproximation Λ t ω) +
        abs (twoCornerInterpolate Λ (twoCornerApproximation Λ) t ω -
          twoCornerInterpolate Λ fermionicKernel t ω) := abs_add_le _ _
    _ ≤ (1 / 8 : ℝ) + 1 / 8 := by
      gcongr
      simpa only [abs_sub_comm] using interpolated_error
    _ = 1 / 4 := by norm_num

/-- The reflected corner is a nonzero second pivot whenever the cutoff is positive. -/
theorem fermionicKernel_firstPivot_reflected_ne {Λ : ℝ} (cutoff_pos : 0 < Λ) :
    residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
      1 (-Λ) ≠ 0 := by
  rw [fermionicKernel_firstPivot_reflected_residual]
  exact sub_ne_zero.mpr (ne_of_gt (Real.exp_lt_one_iff.mpr (by linarith)))

private theorem cornerRightWeight_eq_firstPivot_ratio {Λ t : ℝ}
    (cutoff_pos : 0 < Λ) :
    cornerRightWeight Λ t =
      residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' t (-Λ) /
        residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) := by
  rw [fermionicKernel_firstPivot_residual,
    fermionicKernel_firstPivot_reflected_residual]
  unfold cornerRightWeight
  rw [Real.sinh_eq, Real.sinh_eq]
  have pivot_pos : 0 < 1 - Real.exp (-Λ) :=
    sub_pos.mpr (Real.exp_lt_one_iff.mpr (by linarith))
  have exp_diff_pos : 0 < Real.exp Λ - Real.exp (-Λ) :=
    sub_pos.mpr (Real.exp_lt_exp.mpr (by linarith))
  have factor_identity :
      (1 + Real.exp Λ) * (1 - Real.exp (-Λ)) =
        Real.exp Λ - Real.exp (-Λ) := by
    have exp_mul : Real.exp Λ * Real.exp (-Λ) = 1 := by
      rw [← Real.exp_add]
      simp
    nlinarith
  field_simp [pivot_pos.ne', exp_diff_pos.ne', Real.exp_ne_zero]
  calc
    (Real.exp (Λ * t) - Real.exp (-(Λ * t))) * (1 + Real.exp Λ) *
        (1 - Real.exp (-Λ)) =
      (Real.exp (Λ * t) - Real.exp (-(Λ * t))) *
        ((1 + Real.exp Λ) * (1 - Real.exp (-Λ))) := by ring
    _ = (Real.exp (Λ * t) - Real.exp (-(Λ * t))) *
        (Real.exp Λ - Real.exp (-Λ)) := by rw [factor_identity]

private theorem cornerLeftWeight_eq_sequentialCoefficient {Λ t : ℝ}
    (cutoff_pos : 0 < Λ) :
    cornerLeftWeight Λ t =
      fermionicKernel t Λ / fermionicKernel 0 Λ -
        cornerRightWeight Λ t * (fermionicKernel 1 Λ / fermionicKernel 0 Λ) := by
  have sinh_ne : Real.sinh Λ ≠ 0 := (Real.sinh_pos_iff.mpr cutoff_pos).ne'
  have exp_diff_pos : 0 < Real.exp Λ - Real.exp (-Λ) :=
    sub_pos.mpr (Real.exp_lt_exp.mpr (by linarith))
  unfold cornerLeftWeight cornerRightWeight fermionicKernel
  simp only [zero_mul, neg_zero, Real.exp_zero, Real.sinh_eq]
  field_simp [sinh_ne, exp_diff_pos.ne', fermionicKernel_denominator_ne,
    Real.exp_ne_zero]
  have first_product : Real.exp (Λ - Λ * t) =
      Real.exp Λ * Real.exp (-(Λ * t)) := by
    rw [← Real.exp_add]
    congr 1
  have second_product : Real.exp (-Λ + Λ * t) =
      Real.exp (-Λ) * Real.exp (Λ * t) := by
    rw [← Real.exp_add]
  rw [show Λ * (1 - t) = Λ - Λ * t by ring]
  rw [show -(Λ - Λ * t) = -Λ + Λ * t by ring]
  rw [first_product, second_product]
  ring

/-- The recursive two-corner GECP residual is exactly the endpoint interpolation error. -/
theorem fermionicKernel_twoCornerResidual_eq_sub_interpolate {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ) :
    residualUpdate
        (residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne')
        1 (-Λ) (fermionicKernel_firstPivot_reflected_ne cutoff_pos) t ω =
      fermionicKernel t ω - twoCornerInterpolate Λ fermionicKernel t ω := by
  let R : Kernel ℝ ℝ ℝ :=
    residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
  have pivot_ne : R 1 (-Λ) ≠ 0 := by
    simpa only [R] using fermionicKernel_firstPivot_reflected_ne cutoff_pos
  have residual_formula (s ν : ℝ) :
      R s ν = fermionicKernel s ν -
        fermionicKernel s Λ * fermionicKernel 0 ν / fermionicKernel 0 Λ := by
    rfl
  have right_ratio : cornerRightWeight Λ t = R t (-Λ) / R 1 (-Λ) := by
    simpa only [R] using cornerRightWeight_eq_firstPivot_ratio (t := t) cutoff_pos
  have left_coefficient := cornerLeftWeight_eq_sequentialCoefficient (t := t) cutoff_pos
  change residualUpdate R 1 (-Λ) pivot_ne t ω = _
  unfold residualUpdate
  have factor_reassociation :
      R t (-Λ) * R 1 ω / R 1 (-Λ) =
        (R t (-Λ) / R 1 (-Λ)) * R 1 ω := by
    field_simp [pivot_ne]
  rw [factor_reassociation, ← right_ratio, residual_formula t ω,
    residual_formula 1 ω]
  unfold twoCornerInterpolate
  rw [left_coefficient]
  ring

/--
For `0 < Λ ≤ 1`, the first two complete fermionic GECP pivots contract the
continuous residual to at most one quarter in absolute norm.
-/
theorem fermionicKernel_twoCornerResidual_abs_le_quarter {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ) (cutoff_le_one : Λ ≤ 1)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (residualUpdate
      (residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne')
      1 (-Λ) (fermionicKernel_firstPivot_reflected_ne cutoff_pos) t ω) ≤ 1 / 4 := by
  rw [fermionicKernel_twoCornerResidual_eq_sub_interpolate cutoff_pos]
  exact fermionicKernel_twoCornerInterpolate_error_le_quarter cutoff_pos cutoff_le_one
    time_nonneg time_le_one frequency_lower frequency_upper

/--
On the base cutoff scale, the actual two-step GECP residual is at most half of
the initial complete-pivot magnitude. This is the first proved block-contraction
case of the cutoff-dependent conjecture.
-/
theorem fermionicKernel_twoCornerResidual_le_half_initial {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ) (cutoff_le_one : Λ ≤ 1)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (residualUpdate
      (residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne')
      1 (-Λ) (fermionicKernel_firstPivot_reflected_ne cutoff_pos) t ω) ≤
        (1 / 2 : ℝ) * fermionicKernel 0 Λ := by
  have residual_bound := fermionicKernel_twoCornerResidual_abs_le_quarter cutoff_pos
    cutoff_le_one time_nonneg time_le_one frequency_lower frequency_upper
  have initial_lower : (1 / 2 : ℝ) ≤ fermionicKernel 0 Λ := by
    unfold fermionicKernel
    simp only [zero_mul, neg_zero, Real.exp_zero]
    have exp_le_one : Real.exp (-Λ) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith)
    apply (le_div_iff₀ (fermionicKernel_denominator_pos Λ)).mpr
    nlinarith
  exact residual_bound.trans (by nlinarith)

end Fermionic
end GECPKernelStructure
