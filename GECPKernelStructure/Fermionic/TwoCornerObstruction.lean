import GECPKernelStructure.Fermionic.PowerSecantResidual
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

/--
At the rational-coordinate witness `a = 24^-3`, `b = 24^3`,
`z = (5/4)^3`, and `t = 2/3`, the normalized power-secant error is an
explicit rational number.
-/
theorem powerSecant_twoThirds_counterexample_value :
    (((125 / 64 : ℝ) ^ (2 / 3 : ℝ) -
      powerSecant (1 / 13824) 13824 (2 / 3) (125 / 64)) /
      (1 + 125 / 64)) = 4495348 / 8973531 := by
  norm_num [powerSecant]

/-- The exact excess over half the initial pivot at the rational witness. -/
theorem powerSecant_twoThirds_counterexample_margin :
    (((125 / 64 : ℝ) ^ (2 / 3 : ℝ) -
      powerSecant (1 / 13824) 13824 (2 / 3) (125 / 64)) /
      (1 + 125 / 64)) - (1 / 2) * (13824 / 13825) =
        222676 / 224338275 := by
  rw [powerSecant_twoThirds_counterexample_value]
  norm_num

/-- The rational-coordinate witness lies strictly above half the initial pivot. -/
theorem powerSecant_twoThirds_counterexample_gt_half :
    (1 / 2 : ℝ) * (13824 / 13825) <
      (((125 / 64 : ℝ) ^ (2 / 3 : ℝ) -
        powerSecant (1 / 13824) 13824 (2 / 3) (125 / 64)) /
        (1 + 125 / 64)) := by
  have margin_pos : (0 : ℝ) < 222676 / 224338275 := by norm_num
  linarith [powerSecant_twoThirds_counterexample_margin]

private theorem counterexample_frequency_order :
    -Real.log (13824 : ℝ) < Real.log 13824 := by
  have log_pos : 0 < Real.log (13824 : ℝ) := Real.log_pos (by norm_num)
  linarith

/--
The rational secant witness is realized by the actual two-corner fermionic
Schur residual on the symmetric band with cutoff `log 13824`.
-/
theorem fermionicKernel_twoCornerResidual_counterexample_value :
    residualUpdate
      (residualUpdate fermionicKernel 0 (Real.log 13824)
        (fermionicKernel_pos 0 (Real.log 13824)).ne')
      1 (-Real.log 13824)
        (fermionicKernel_firstCorner_secondPivot_pos
          counterexample_frequency_order).ne'
      (2 / 3) (-Real.log (125 / 64)) = 4495348 / 8973531 := by
  rw [fermionicKernel_twoCornerResidual_eq_powerSecantError
    counterexample_frequency_order]
  rw [show Real.exp (-(-Real.log (125 / 64))) = (125 / 64 : ℝ) by
        rw [neg_neg, Real.exp_log (by norm_num : (0 : ℝ) < 125 / 64)],
      show Real.exp (-Real.log 13824) = (1 / 13824 : ℝ) by
        rw [Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 13824)]
        norm_num,
      show Real.exp (-(-Real.log 13824)) = (13824 : ℝ) by
        rw [neg_neg, Real.exp_log (by norm_num : (0 : ℝ) < 13824)]]
  exact powerSecant_twoThirds_counterexample_value

/--
Two symmetric corner pivots do not contract the fermionic residual by one half
at the explicit interior witness.
-/
theorem fermionicKernel_twoCornerResidual_counterexample_gt_half_initial :
    (1 / 2 : ℝ) * fermionicKernel 0 (Real.log 13824) <
      residualUpdate
        (residualUpdate fermionicKernel 0 (Real.log 13824)
          (fermionicKernel_pos 0 (Real.log 13824)).ne')
        1 (-Real.log 13824)
          (fermionicKernel_firstCorner_secondPivot_pos
            counterexample_frequency_order).ne'
        (2 / 3) (-Real.log (125 / 64)) := by
  rw [fermionicKernel_twoCornerResidual_counterexample_value,
    fermionicKernel_eq_rpow]
  rw [show Real.exp (-Real.log 13824) = (1 / 13824 : ℝ) by
        rw [Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 13824)]
        norm_num]
  norm_num

/--
At cutoff `log 13824`, the two symmetric corner pivots fail the universal
half-contraction property on the physical time-frequency rectangle.
-/
theorem fermionicKernel_twoCornerResidual_not_half_contraction :
    ¬ ∀ t ω : ℝ,
      0 ≤ t → t ≤ 1 → -Real.log 13824 ≤ ω → ω ≤ Real.log 13824 →
        abs (residualUpdate
          (residualUpdate fermionicKernel 0 (Real.log 13824)
            (fermionicKernel_pos 0 (Real.log 13824)).ne')
          1 (-Real.log 13824)
            (fermionicKernel_firstCorner_secondPivot_pos
              counterexample_frequency_order).ne' t ω) ≤
          (1 / 2 : ℝ) * fermionicKernel 0 (Real.log 13824) := by
  intro half_contraction
  have frequency_lower :
      -Real.log (13824 : ℝ) ≤ -Real.log (125 / 64) := by
    apply neg_le_neg
    exact Real.log_le_log (by norm_num) (by norm_num)
  have frequency_upper :
      -Real.log (125 / 64 : ℝ) ≤ Real.log 13824 := by
    have interior_log_pos : 0 < Real.log (125 / 64 : ℝ) :=
      Real.log_pos (by norm_num)
    have cutoff_log_pos : 0 < Real.log (13824 : ℝ) :=
      Real.log_pos (by norm_num)
    linarith
  have bound := half_contraction (2 / 3) (-Real.log (125 / 64))
    (by norm_num) (by norm_num) frequency_lower frequency_upper
  have residual_pos :
      0 < residualUpdate
        (residualUpdate fermionicKernel 0 (Real.log 13824)
          (fermionicKernel_pos 0 (Real.log 13824)).ne')
        1 (-Real.log 13824)
          (fermionicKernel_firstCorner_secondPivot_pos
            counterexample_frequency_order).ne'
        (2 / 3) (-Real.log (125 / 64)) := by
    rw [fermionicKernel_twoCornerResidual_counterexample_value]
    norm_num
  rw [abs_of_pos residual_pos] at bound
  exact (not_le_of_gt
    fermionicKernel_twoCornerResidual_counterexample_gt_half_initial) bound

end Fermionic
end GECPKernelStructure
