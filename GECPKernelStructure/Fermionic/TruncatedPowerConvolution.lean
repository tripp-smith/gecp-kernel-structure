import GECPKernelStructure.Fermionic.TruncatedPowerComposition
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

namespace GECPKernelStructure
namespace Fermionic

open Set

/-- Factorially normalized truncated powers, with the support boundary excluded. -/
noncomputable def normalizedTruncatedPowerKernel (degree : ℕ) (v ω : ℝ) : ℝ :=
  if v < ω then (ω - v) ^ degree / (Nat.factorial degree : ℝ) else 0

/-- Degree zero is the elementary strict order kernel. -/
theorem normalizedTruncatedPowerKernel_zero (v ω : ℝ) :
    normalizedTruncatedPowerKernel 0 v ω = if v < ω then 1 else 0 := by
  simp [normalizedTruncatedPowerKernel]

/-- The existing seventh-degree kernel is the unnormalized member of the family. -/
theorem eighthOrderTruncatedPowerKernel_eq_factorial_mul_normalized
    (v ω : ℝ) :
    eighthOrderTruncatedPowerKernel v ω =
      (Nat.factorial 7 : ℝ) * normalizedTruncatedPowerKernel 7 v ω := by
  by_cases h : v < ω
  · rw [eighthOrderTruncatedPowerKernel, max_eq_left (by linarith)]
    simp only [normalizedTruncatedPowerKernel, if_pos h]
    field_simp
  · have hnonpos : ω - v ≤ 0 := by linarith
    rw [eighthOrderTruncatedPowerKernel, max_eq_right hnonpos]
    simp [normalizedTruncatedPowerKernel, h]

/-- Integrating one normalized truncated power raises its degree by one. -/
theorem normalizedTruncatedPowerKernel_succ_integral
    (degree : ℕ) {v ω : ℝ} (hvo : v ≤ ω) :
    normalizedTruncatedPowerKernel (degree + 1) v ω =
      ∫ z in v..ω, normalizedTruncatedPowerKernel degree z ω := by
  by_cases hlt : v < ω
  · rw [intervalIntegral.integral_congr_Ioo_of_le hvo]
    · simp only [normalizedTruncatedPowerKernel, if_pos hlt]
      rw [intervalIntegral.integral_div (Nat.factorial degree : ℝ)
        (fun z : ℝ => (ω - z) ^ degree)]
      rw [intervalIntegral.integral_comp_sub_left (fun x : ℝ => x ^ degree) ω]
      simp only [sub_self, integral_pow, zero_pow (Nat.succ_ne_zero degree), sub_zero]
      rw [Nat.factorial_succ]
      push_cast
      field_simp
    · intro z hz
      simp [normalizedTruncatedPowerKernel, hz.2]
  · have heq : v = ω := le_antisymm hvo (le_of_not_gt hlt)
    subst ω
    simp [normalizedTruncatedPowerKernel]

/--
The degree-zero order kernel makes the successor recurrence valid on any
shared interval containing the two endpoints.
-/
theorem normalizedTruncatedPowerKernel_succ_integral_commonDomain
    (degree : ℕ) {lower v ω upper : ℝ} (hlv : lower ≤ v)
    (hvo : v ≤ ω) (hωu : ω ≤ upper) :
    normalizedTruncatedPowerKernel (degree + 1) v ω =
      ∫ z in lower..upper,
        normalizedTruncatedPowerKernel 0 v z *
          normalizedTruncatedPowerKernel degree z ω := by
  let f : ℝ → ℝ := fun z =>
    normalizedTruncatedPowerKernel 0 v z *
      normalizedTruncatedPowerKernel degree z ω
  let g : ℝ → ℝ := fun z =>
    (ω - z) ^ degree / (Nat.factorial degree : ℝ)
  have hfleft : IntervalIntegrable f MeasureTheory.volume lower v := by
    apply IntervalIntegrable.zero.congr_uIoo
    intro z hz
    rw [uIoo_of_le hlv] at hz
    simp [f, normalizedTruncatedPowerKernel, not_lt_of_ge hz.2.le]
  have hfmiddle : IntervalIntegrable f MeasureTheory.volume v ω := by
    apply ((show Continuous g by unfold g; fun_prop).intervalIntegrable v ω).congr_uIoo
    intro z hz
    rw [uIoo_of_le hvo] at hz
    simp [f, g, normalizedTruncatedPowerKernel, hz.1, hz.2]
  have hfright : IntervalIntegrable f MeasureTheory.volume ω upper := by
    apply IntervalIntegrable.zero.congr_uIoo
    intro z hz
    rw [uIoo_of_le hωu] at hz
    simp [f, normalizedTruncatedPowerKernel, not_lt_of_ge hz.1.le]
  have hleft : (∫ z in lower..v, f z) = 0 := by
    calc
      (∫ z in lower..v, f z) = ∫ _z in lower..v, (0 : ℝ) := by
        apply intervalIntegral.integral_congr_Ioo_of_le hlv
        intro z hz
        simp [f, normalizedTruncatedPowerKernel, not_lt_of_ge hz.2.le]
      _ = 0 := by simp
  have hmiddle :
      (∫ z in v..ω, f z) =
        ∫ z in v..ω, normalizedTruncatedPowerKernel degree z ω := by
    rw [intervalIntegral.integral_congr_Ioo_of_le hvo]
    intro z hz
    simp [f, normalizedTruncatedPowerKernel, hz.1]
  have hright : (∫ z in ω..upper, f z) = 0 := by
    calc
      (∫ z in ω..upper, f z) = ∫ _z in ω..upper, (0 : ℝ) := by
        apply intervalIntegral.integral_congr_Ioo_of_le hωu
        intro z hz
        simp [f, normalizedTruncatedPowerKernel, not_lt_of_ge hz.1.le]
      _ = 0 := by simp
  have haddleft := intervalIntegral.integral_add_adjacent_intervals hfleft hfmiddle
  have haddright := intervalIntegral.integral_add_adjacent_intervals
    (hfleft.trans hfmiddle) hfright
  calc
    normalizedTruncatedPowerKernel (degree + 1) v ω =
        ∫ z in v..ω, normalizedTruncatedPowerKernel degree z ω :=
      normalizedTruncatedPowerKernel_succ_integral degree hvo
    _ = ∫ z in v..ω, f z := hmiddle.symm
    _ = ∫ z in lower..ω, f z := by rw [← haddleft, hleft, zero_add]
    _ = ∫ z in lower..upper, f z := by rw [← haddright, hright, add_zero]

/-- A normalized truncated power is interval-integrable up to any common upper bound. -/
theorem normalizedTruncatedPowerKernel_intervalIntegrable_to_commonUpper
    (degree : ℕ) {lower ω upper : ℝ} (hlower : lower ≤ upper)
    (hω : ω ≤ upper) :
    IntervalIntegrable
      (fun z => normalizedTruncatedPowerKernel degree z ω)
      MeasureTheory.volume lower upper := by
  let polynomial : ℝ → ℝ := fun z =>
    (ω - z) ^ degree / (Nat.factorial degree : ℝ)
  by_cases hlω : lower ≤ ω
  · have left : IntervalIntegrable
        (fun z => normalizedTruncatedPowerKernel degree z ω)
        MeasureTheory.volume lower ω := by
      apply ((show Continuous polynomial by unfold polynomial; fun_prop).intervalIntegrable
        lower ω).congr_uIoo
      intro z hz
      rw [uIoo_of_le hlω] at hz
      simp [normalizedTruncatedPowerKernel, polynomial, hz.2]
    have right : IntervalIntegrable
        (fun z => normalizedTruncatedPowerKernel degree z ω)
        MeasureTheory.volume ω upper := by
      apply IntervalIntegrable.zero.congr_uIoo
      intro z hz
      rw [uIoo_of_le hω] at hz
      simp [normalizedTruncatedPowerKernel, not_lt_of_ge hz.1.le]
    exact left.trans right
  · have hωlower : ω ≤ lower := le_of_not_ge hlω
    apply IntervalIntegrable.zero.congr_uIoo
    intro z hz
    rw [uIoo_of_le hlower] at hz
    simp [normalizedTruncatedPowerKernel,
      not_lt_of_ge (hωlower.trans hz.1.le)]

/-- The successor recurrence may end at any common upper bound beyond the column knot. -/
theorem normalizedTruncatedPowerKernel_succ_integral_to_commonUpper
    (degree : ℕ) {v ω upper : ℝ} (hv : v ≤ upper) (hω : ω ≤ upper) :
    normalizedTruncatedPowerKernel (degree + 1) v ω =
      ∫ z in v..upper, normalizedTruncatedPowerKernel degree z ω := by
  by_cases hvo : v ≤ ω
  · rw [normalizedTruncatedPowerKernel_succ_integral_commonDomain
      degree le_rfl hvo hω]
    apply intervalIntegral.integral_congr_Ioo_of_le hv
    intro z hz
    simp [normalizedTruncatedPowerKernel, hz.1]
  · have hωv : ω ≤ v := le_of_not_ge hvo
    have left_zero : normalizedTruncatedPowerKernel (degree + 1) v ω = 0 := by
      simp [normalizedTruncatedPowerKernel, not_lt_of_ge hωv]
    rw [left_zero]
    symm
    calc
      (∫ z in v..upper, normalizedTruncatedPowerKernel degree z ω) =
          ∫ _z in v..upper, (0 : ℝ) := by
        apply intervalIntegral.integral_congr_Ioo_of_le hv
        intro z hz
        simp [normalizedTruncatedPowerKernel,
          not_lt_of_ge (hωv.trans hz.1.le)]
      _ = 0 := by simp

end Fermionic
end GECPKernelStructure
