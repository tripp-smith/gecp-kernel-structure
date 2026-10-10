import GECPKernelStructure.Fermionic.TruncatedPowerComposition
import GECPKernelStructure.GECP.SignRegularity

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The seventh-degree truncated-power cross product is nonnegative on ordered pairs. -/
theorem eighthOrderTruncatedPowerKernel_cross_nonneg
    {v₀ v₁ ω₀ ω₁ : ℝ} (hv : v₀ ≤ v₁) (hω : ω₀ ≤ ω₁) :
    0 ≤
      eighthOrderTruncatedPowerKernel v₀ ω₀ *
          eighthOrderTruncatedPowerKernel v₁ ω₁ -
        eighthOrderTruncatedPowerKernel v₀ ω₁ *
          eighthOrderTruncatedPowerKernel v₁ ω₀ := by
  unfold eighthOrderTruncatedPowerKernel
  by_cases hfull : v₁ < ω₀
  · have h₀₀ : 0 ≤ ω₀ - v₀ := by linarith
    have h₁₁ : 0 ≤ ω₁ - v₁ := by linarith
    have h₀₁ : 0 ≤ ω₁ - v₀ := by linarith
    have h₁₀ : 0 ≤ ω₀ - v₁ := by linarith
    rw [max_eq_left h₀₀, max_eq_left h₁₁,
      max_eq_left h₀₁, max_eq_left h₁₀]
    rw [← mul_pow, ← mul_pow]
    apply sub_nonneg.mpr
    apply pow_le_pow_left₀ (mul_nonneg h₀₁ h₁₀)
    nlinarith
  · have h₁₀ : ω₀ - v₁ ≤ 0 := by linarith
    rw [max_eq_right h₁₀]
    norm_num
    positivity

/--
The order-two minor is strictly positive under the paired support conditions;
this is the first Schoenberg--Whitney region for the truncated-power kernel.
-/
theorem eighthOrderTruncatedPowerKernel_cross_pos
    {v₀ v₁ ω₀ ω₁ : ℝ} (hv : v₀ < v₁) (hω : ω₀ < ω₁)
    (hdiag₀ : v₀ < ω₀) (hdiag₁ : v₁ < ω₁) :
    0 <
      eighthOrderTruncatedPowerKernel v₀ ω₀ *
          eighthOrderTruncatedPowerKernel v₁ ω₁ -
        eighthOrderTruncatedPowerKernel v₀ ω₁ *
          eighthOrderTruncatedPowerKernel v₁ ω₀ := by
  unfold eighthOrderTruncatedPowerKernel
  by_cases hfull : v₁ < ω₀
  · have h₀₀ : 0 ≤ ω₀ - v₀ := by linarith
    have h₁₁ : 0 ≤ ω₁ - v₁ := by linarith
    have h₀₁ : 0 ≤ ω₁ - v₀ := by linarith
    have h₁₀ : 0 ≤ ω₀ - v₁ := by linarith
    rw [max_eq_left h₀₀, max_eq_left h₁₁,
      max_eq_left h₀₁, max_eq_left h₁₀]
    rw [← mul_pow, ← mul_pow]
    apply sub_pos.mpr
    apply pow_lt_pow_left₀
    · nlinarith
    · exact mul_nonneg h₀₁ h₁₀
    · norm_num
  · have h₁₀ : ω₀ - v₁ ≤ 0 := by linarith
    have h₀₀ : 0 < ω₀ - v₀ := by linarith
    have h₁₁ : 0 < ω₁ - v₁ := by linarith
    rw [max_eq_right h₁₀, max_eq_left h₀₀.le, max_eq_left h₁₁.le]
    norm_num
    positivity

/-- On strict row and column pairs, the paired support conditions exactly characterize strictness. -/
theorem eighthOrderTruncatedPowerKernel_cross_pos_iff
    {v₀ v₁ ω₀ ω₁ : ℝ} (hv : v₀ < v₁) (hω : ω₀ < ω₁) :
    0 <
        eighthOrderTruncatedPowerKernel v₀ ω₀ *
            eighthOrderTruncatedPowerKernel v₁ ω₁ -
          eighthOrderTruncatedPowerKernel v₀ ω₁ *
            eighthOrderTruncatedPowerKernel v₁ ω₀ ↔
      v₀ < ω₀ ∧ v₁ < ω₁ := by
  constructor
  · intro hpos
    constructor
    · by_contra hnot
      have h₀₀ : ω₀ - v₀ ≤ 0 := by linarith
      have h₁₀ : ω₀ - v₁ ≤ 0 := by linarith
      unfold eighthOrderTruncatedPowerKernel at hpos
      rw [max_eq_right h₀₀, max_eq_right h₁₀] at hpos
      norm_num at hpos
    · by_contra hnot
      have h₁₁ : ω₁ - v₁ ≤ 0 := by linarith
      have h₁₀ : ω₀ - v₁ ≤ 0 := by linarith
      unfold eighthOrderTruncatedPowerKernel at hpos
      rw [max_eq_right h₁₁, max_eq_right h₁₀] at hpos
      norm_num at hpos
  · rintro ⟨hdiag₀, hdiag₁⟩
    exact eighthOrderTruncatedPowerKernel_cross_pos hv hω hdiag₀ hdiag₁

/-- The truncated-power kernel is nonnegative at order one. -/
theorem eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_one :
    TotallyNonnegativeAtOrder eighthOrderTruncatedPowerKernel 1 := by
  intro rows columns _ _
  rw [Matrix.det_fin_one]
  exact pow_nonneg (le_max_right _ _) _

/-- The truncated-power kernel is totally nonnegative at the first nontrivial order. -/
theorem eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_two :
    TotallyNonnegativeAtOrder eighthOrderTruncatedPowerKernel 2 := by
  intro rows columns rows_mono columns_mono
  rw [Matrix.det_fin_two]
  exact eighthOrderTruncatedPowerKernel_cross_nonneg
    (rows_mono (show (0 : Fin 2) < 1 by decide)).le
    (columns_mono (show (0 : Fin 2) < 1 by decide)).le

/-- The paired support conditions make every ordered `2 x 2` minor strictly positive. -/
theorem eighthOrderTruncatedPowerKernel_fin_two_det_pos
    (rows columns : Fin 2 → ℝ) (rows_mono : StrictMono rows)
    (columns_mono : StrictMono columns)
    (support₀ : rows 0 < columns 0) (support₁ : rows 1 < columns 1) :
    0 < (minorMatrix eighthOrderTruncatedPowerKernel rows columns).det := by
  rw [Matrix.det_fin_two]
  exact eighthOrderTruncatedPowerKernel_cross_pos
    (rows_mono (show (0 : Fin 2) < 1 by decide))
    (columns_mono (show (0 : Fin 2) < 1 by decide)) support₀ support₁

/-- The order-two Schoenberg--Whitney support conditions are also necessary. -/
theorem eighthOrderTruncatedPowerKernel_fin_two_det_pos_iff
    (rows columns : Fin 2 → ℝ) (rows_mono : StrictMono rows)
    (columns_mono : StrictMono columns) :
    0 < (minorMatrix eighthOrderTruncatedPowerKernel rows columns).det ↔
      rows 0 < columns 0 ∧ rows 1 < columns 1 := by
  rw [Matrix.det_fin_two]
  exact eighthOrderTruncatedPowerKernel_cross_pos_iff
    (rows_mono (show (0 : Fin 2) < 1 by decide))
    (columns_mono (show (0 : Fin 2) < 1 by decide))

end Fermionic
end GECPKernelStructure
