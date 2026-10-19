import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Tactic

/-! A three-node sign criterion obtained from Cauchy's mean value theorem.
The derivative denominator need not be positive: two negative increments
also have the positive product required by the criterion. -/

namespace Matrix

open Set

def threePointCross (f g : ℝ → ℝ) (x y z : ℝ) : ℝ :=
  (f y - f x) * (g z - g y) - (g y - g x) * (f z - f y)

/-- Cauchy's mean value theorem converts an increasing derivative ratio into
the sign of a three-node generalized Vandermonde determinant. -/
theorem threePointCross_pos_of_derivativeRatio
    (f g f' g' : ℝ → ℝ) {x y z : ℝ}
    (hxy : x < y) (hyz : y < z)
    (fcont : ContinuousOn f (Icc x z)) (gcont : ContinuousOn g (Icc x z))
    (fderiv : ∀ t ∈ Ioo x z, HasDerivAt f (f' t) t)
    (gderiv : ∀ t ∈ Ioo x z, HasDerivAt g (g' t) t)
    (fderiv_ne : ∀ t ∈ Ioo x z, f' t ≠ 0)
    (ratio_mono : StrictMonoOn (fun t ↦ g' t / f' t) (Ioo x z))
    (increments_pos : 0 < (f y - f x) * (f z - f y)) :
    0 < threePointCross f g x y z := by
  obtain ⟨u, hu, eu⟩ := exists_ratio_hasDerivAt_eq_ratio_slope
    f f' hxy
    (fcont.mono (Icc_subset_Icc le_rfl hyz.le))
    (fun t ht ↦ fderiv t ⟨ht.1, ht.2.trans hyz⟩)
    g g'
    (gcont.mono (Icc_subset_Icc le_rfl hyz.le))
    (fun t ht ↦ gderiv t ⟨ht.1, ht.2.trans hyz⟩)
  obtain ⟨v, hv, ev⟩ := exists_ratio_hasDerivAt_eq_ratio_slope
    f f' hyz
    (fcont.mono (Icc_subset_Icc hxy.le le_rfl))
    (fun t ht ↦ fderiv t ⟨hxy.trans ht.1, ht.2⟩)
    g g'
    (gcont.mono (Icc_subset_Icc hxy.le le_rfl))
    (fun t ht ↦ gderiv t ⟨hxy.trans ht.1, ht.2⟩)
  have hu' : u ∈ Ioo x z := ⟨hu.1, hu.2.trans hyz⟩
  have hv' : v ∈ Ioo x z := ⟨hxy.trans hv.1, hv.2⟩
  have gu : g y - g x = (f y - f x) * (g' u / f' u) := by
    calc
      _ = ((f y - f x) * g' u) / f' u :=
        (eq_div_iff (fderiv_ne u hu')).mpr eu
      _ = _ := by ring
  have gv : g z - g y = (f z - f y) * (g' v / f' v) := by
    calc
      _ = ((f z - f y) * g' v) / f' v :=
        (eq_div_iff (fderiv_ne v hv')).mpr ev
      _ = _ := by ring
  have ratio_lt := ratio_mono hu' hv' (hu.2.trans hv.1)
  calc
    0 < ((f y - f x) * (f z - f y)) * (g' v / f' v - g' u / f' u) :=
      mul_pos increments_pos (sub_pos.mpr ratio_lt)
    _ = threePointCross f g x y z := by
      unfold threePointCross
      rw [gu, gv]
      ring

theorem threePointCross_neg_right (f g : ℝ → ℝ) (x y z : ℝ) :
    threePointCross f (fun t ↦ -g t) x y z = -threePointCross f g x y z := by
  unfold threePointCross
  ring

end Matrix
