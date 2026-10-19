import GECPKernelStructure.Matrix.ThreePointChebyshev
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

/-! Signs of actual-exponential/monomial augmented minors on the unit band.

After dividing rows by their positive exponential factors, explicit ratios
of derivatives give the three-node sign. The degree-zero branch changes
sign and must be handled separately. These results concern one exponential
column and two monomial columns, not arbitrary conditioned pivot prefixes. -/

namespace GECPKernelStructure.Fermionic

open Set Matrix

noncomputable def expNormalizedPower (rate : ℝ) (degree : ℕ) (t : ℝ) : ℝ :=
  t ^ degree * Real.exp (-rate * t)

noncomputable def expNormalizedPowerSuccDeriv (rate : ℝ) (m : ℕ) (t : ℝ) : ℝ :=
  t ^ m * Real.exp (-rate * t) * ((m + 1 : ℝ) - rate * t)

theorem hasDerivAt_expNormalizedPower_succ (rate : ℝ) (m : ℕ) (t : ℝ) :
    HasDerivAt (expNormalizedPower rate (m + 1))
      (expNormalizedPowerSuccDeriv rate m t) t := by
  have hp : HasDerivAt (fun x : ℝ ↦ x ^ (m + 1))
      ((m + 1 : ℝ) * t ^ m) t := by
    simpa using hasDerivAt_pow (m + 1) t
  have he : HasDerivAt (fun x : ℝ ↦ Real.exp (-rate * x))
      (Real.exp (-rate * t) * (-rate)) t := by
    simpa using ((hasDerivAt_id t).const_mul (-rate)).exp
  apply (hp.mul he).congr_deriv
  unfold expNormalizedPowerSuccDeriv
  rw [pow_succ]
  ring

noncomputable def expPositiveDegreeRatio (rate : ℝ) (m n : ℕ) (t : ℝ) : ℝ :=
  t ^ (n + 1) * ((m + n + 2 : ℝ) - rate * t) /
    ((m + 1 : ℝ) - rate * t)

noncomputable def expPositiveDegreeRatioDeriv (rate : ℝ) (m n : ℕ) (t : ℝ) : ℝ :=
  (n + 1 : ℝ) * t ^ n *
      (((m + 1 : ℝ) * (m + n + 2 : ℝ)) -
        (2 * m + n + 2 : ℝ) * (rate * t) + (rate * t) ^ 2) /
    ((m + 1 : ℝ) - rate * t) ^ 2

theorem hasDerivAt_expPositiveDegreeRatio
    (rate : ℝ) (m n : ℕ) (t : ℝ)
    (denom_ne : (m + 1 : ℝ) - rate * t ≠ 0) :
    HasDerivAt (expPositiveDegreeRatio rate m n)
      (expPositiveDegreeRatioDeriv rate m n t) t := by
  have hnum := ((hasDerivAt_id t).pow (n + 1)).mul
    ((hasDerivAt_const t (m + n + 2 : ℝ)).sub
      ((hasDerivAt_id t).const_mul rate))
  have hden := (hasDerivAt_const t (m + 1 : ℝ)).sub
    ((hasDerivAt_id t).const_mul rate)
  apply (hnum.div hden denom_ne).congr_deriv
  simp only [expPositiveDegreeRatioDeriv, Nat.add_sub_cancel,
    Nat.cast_add, Nat.cast_one, id_eq, mul_one, zero_sub,
    Pi.sub_apply, Pi.mul_apply, Pi.pow_apply]
  rw [pow_succ]
  field_simp
  ring

theorem expPositiveDegreeRatioDeriv_pos
    {rate t : ℝ} (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (m n : ℕ) (ht : t ∈ Ioo (0 : ℝ) 1) :
    0 < expPositiveDegreeRatioDeriv rate m n t := by
  have zpos : 0 < rate * t := mul_pos rate_pos ht.1
  have zlt : rate * t < 1 := by
    calc
      rate * t ≤ 1 * t := mul_le_mul_of_nonneg_right rate_le ht.1.le
      _ < 1 := by simpa using ht.2
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg _
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg _
  have one_minus_z : 0 ≤ 1 - rate * t := by linarith
  have extra : 0 ≤ (2 * m + n : ℝ) * (1 - rate * t) :=
    mul_nonneg (by positivity) one_minus_z
  have bracket_pos : 0 < ((m + 1 : ℝ) * (m + n + 2 : ℝ)) -
      (2 * m + n + 2 : ℝ) * (rate * t) + (rate * t) ^ 2 := by
    have hmn : 0 ≤ (m : ℝ) * (m + n + 1 : ℝ) := by positivity
    nlinarith [sq_nonneg (1 - rate * t)]
  unfold expPositiveDegreeRatioDeriv
  apply div_pos
  · exact mul_pos (mul_pos (by positivity) (pow_pos ht.1 n)) bracket_pos
  · exact sq_pos_of_pos (by linarith : 0 < (m + 1 : ℝ) - rate * t)

theorem expNormalizedPowerSuccDeriv_pos
    {rate t : ℝ} (rate_le : rate ≤ 1) (m : ℕ)
    (ht : t ∈ Ioo (0 : ℝ) 1) :
    0 < expNormalizedPowerSuccDeriv rate m t := by
  have zlt : rate * t < 1 := by
    calc
      rate * t ≤ 1 * t := mul_le_mul_of_nonneg_right rate_le ht.1.le
      _ < 1 := by simpa using ht.2
  unfold expNormalizedPowerSuccDeriv
  apply mul_pos (mul_pos (pow_pos ht.1 _) (Real.exp_pos _))
  have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg _
  linarith

theorem expNormalizedPowerSuccDeriv_ratio (rate : ℝ) (m n : ℕ) (t : ℝ)
    (hd : (m + 1 : ℝ) - rate * t ≠ 0) :
    expNormalizedPowerSuccDeriv rate (m + n + 1) t =
      expNormalizedPowerSuccDeriv rate m t * expPositiveDegreeRatio rate m n t := by
  unfold expNormalizedPowerSuccDeriv expPositiveDegreeRatio
  simp only [Nat.cast_add, Nat.cast_one]
  field_simp
  rw [show m + n + 1 = m + (n + 1) by omega, pow_add]
  ring

theorem expNormalizedPower_succ_strictMonoOn {rate : ℝ} (rate_le : rate ≤ 1) (m : ℕ) :
    StrictMonoOn (expNormalizedPower rate (m + 1)) (Icc 0 1) := by
  apply strictMonoOn_of_deriv_pos (convex_Icc 0 1)
  · unfold expNormalizedPower
    fun_prop
  · intro t ht
    rw [(hasDerivAt_expNormalizedPower_succ rate m t).deriv]
    exact expNormalizedPowerSuccDeriv_pos rate_le m (by simpa using ht)

theorem expPositiveDegreeRatio_strictMonoOn {rate : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1) (m n : ℕ) :
    StrictMonoOn (expPositiveDegreeRatio rate m n) (Ioo 0 1) := by
  have hd : ∀ t ∈ Ioo (0 : ℝ) 1, (m + 1 : ℝ) - rate * t ≠ 0 := by
    intro t ht
    have hz : rate * t < 1 := (mul_le_mul_of_nonneg_right rate_le ht.1.le).trans_lt
      (by simpa using ht.2)
    have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg _
    linarith
  apply strictMonoOn_of_deriv_pos (convex_Ioo 0 1)
  · exact fun t ht ↦ (hasDerivAt_expPositiveDegreeRatio rate m n t (hd t ht)).continuousAt.continuousWithinAt
  · intro t ht
    have ht' : t ∈ Ioo (0 : ℝ) 1 := by simpa using ht
    rw [(hasDerivAt_expPositiveDegreeRatio rate m n t (hd t ht')).deriv]
    exact expPositiveDegreeRatioDeriv_pos rate_pos rate_le m n ht'

theorem expNormalizedPower_succ_cross_pos {rate x y z : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1) (m n : ℕ) :
    0 < threePointCross (expNormalizedPower rate (m + 1))
      (expNormalizedPower rate (m + n + 2)) x y z := by
  have inside : ∀ t ∈ Ioo x z, t ∈ Ioo (0 : ℝ) 1 :=
    fun t ht ↦ ⟨hx.trans_lt ht.1, ht.2.trans_le hz⟩
  have deriv_pos := fun t ht ↦ expNormalizedPowerSuccDeriv_pos rate_le m (inside t ht)
  have ratio_eq : ∀ t ∈ Ioo x z,
      expNormalizedPowerSuccDeriv rate (m + n + 1) t /
        expNormalizedPowerSuccDeriv rate m t = expPositiveDegreeRatio rate m n t := by
    intro t ht
    have hm : (0 : ℝ) ≤ m := Nat.cast_nonneg _
    have ht' := inside t ht
    have hrt : rate * t < 1 :=
      (mul_le_mul_of_nonneg_right rate_le ht'.1.le).trans_lt (by simpa using ht'.2)
    rw [expNormalizedPowerSuccDeriv_ratio rate m n t (by linarith)]
    exact mul_div_cancel_left₀ _ (ne_of_gt (deriv_pos t ht))
  have mono := expNormalizedPower_succ_strictMonoOn rate_le m
  have mx : x ∈ Icc (0 : ℝ) 1 := ⟨hx, (hxy.trans hyz).le.trans hz⟩
  have my : y ∈ Icc (0 : ℝ) 1 := ⟨hx.trans hxy.le, hyz.le.trans hz⟩
  have mz : z ∈ Icc (0 : ℝ) 1 := ⟨hx.trans (hxy.trans hyz).le, hz⟩
  apply threePointCross_pos_of_derivativeRatio _ _
    (expNormalizedPowerSuccDeriv rate m)
    (expNormalizedPowerSuccDeriv rate (m + n + 1)) hxy hyz
  · unfold expNormalizedPower; fun_prop
  · unfold expNormalizedPower; fun_prop
  · exact fun t _ ↦ hasDerivAt_expNormalizedPower_succ rate m t
  · intro t _
    simpa [show m + n + 1 + 1 = m + n + 2 by omega] using
      hasDerivAt_expNormalizedPower_succ rate (m + n + 1) t
  · exact fun t ht ↦ ne_of_gt (deriv_pos t ht)
  · intro u hu v hv huv
    dsimp only
    rw [ratio_eq u hu, ratio_eq v hv]
    exact expPositiveDegreeRatio_strictMonoOn rate_pos rate_le m n
      (inside u hu) (inside v hv) huv
  · exact mul_pos (sub_pos.mpr (mono mx my hxy)) (sub_pos.mpr (mono my mz hyz))

theorem hasDerivAt_expNormalizedPower_zero (rate t : ℝ) :
    HasDerivAt (expNormalizedPower rate 0) (-rate * Real.exp (-rate * t)) t := by
  unfold expNormalizedPower
  simp only [pow_zero, one_mul]
  simpa only [mul_comm, id_eq, mul_one] using
    ((hasDerivAt_id t).const_mul (-rate)).exp

theorem expNormalizedPower_zero_strictAnti {rate : ℝ} (rate_pos : 0 < rate) :
    StrictAnti (expNormalizedPower rate 0) := by
  apply strictAnti_of_deriv_neg
  intro t
  rw [(hasDerivAt_expNormalizedPower_zero rate t).deriv]
  exact mul_neg_of_neg_of_pos (neg_neg_of_pos rate_pos) (Real.exp_pos _)

theorem expNormalizedPower_zero_cross_pos_of_ratio {rate x y z : ℝ}
    (rate_pos : 0 < rate) (hxy : x < y) (hyz : y < z)
    (g g' ratio : ℝ → ℝ) (gcont : ContinuousOn g (Icc x z))
    (gderiv : ∀ t ∈ Ioo x z, HasDerivAt g (g' t) t)
    (factor : ∀ t ∈ Ioo x z, g' t = (-rate * Real.exp (-rate * t)) * ratio t)
    (ratio_mono : StrictMonoOn ratio (Ioo x z)) :
    0 < threePointCross (expNormalizedPower rate 0) g x y z := by
  apply threePointCross_pos_of_derivativeRatio _ _
    (fun t ↦ -rate * Real.exp (-rate * t)) g' hxy hyz
  · unfold expNormalizedPower; fun_prop
  · exact gcont
  · exact fun t _ ↦ hasDerivAt_expNormalizedPower_zero rate t
  · exact gderiv
  · intro t _
    exact mul_ne_zero (neg_ne_zero.mpr (ne_of_gt rate_pos)) (Real.exp_ne_zero _)
  · intro u hu v hv huv
    dsimp only
    rw [factor u hu, factor v hv, mul_div_cancel_left₀, mul_div_cancel_left₀]
    · exact ratio_mono hu hv huv
    all_goals exact mul_ne_zero (neg_ne_zero.mpr (ne_of_gt rate_pos)) (Real.exp_ne_zero _)
  · have anti := expNormalizedPower_zero_strictAnti rate_pos
    exact mul_pos_of_neg_of_neg (sub_neg.mpr (anti hxy)) (sub_neg.mpr (anti hyz))

theorem expNormalizedPower_zero_one_cross_pos {rate x y z : ℝ}
    (rate_pos : 0 < rate) (hxy : x < y) (hyz : y < z) :
    0 < threePointCross (expNormalizedPower rate 0) (expNormalizedPower rate 1) x y z := by
  apply expNormalizedPower_zero_cross_pos_of_ratio rate_pos hxy hyz
    (expNormalizedPower rate 1) (expNormalizedPowerSuccDeriv rate 0)
    (fun t ↦ t - 1 / rate)
  · unfold expNormalizedPower; fun_prop
  · exact fun t _ ↦ hasDerivAt_expNormalizedPower_succ rate 0 t
  · intro t _
    unfold expNormalizedPowerSuccDeriv
    field_simp
    ring
  · intro u _ v _ huv
    dsimp only
    linarith

noncomputable def expZeroDegreeRatio (rate : ℝ) (n : ℕ) (t : ℝ) : ℝ :=
  t ^ (n + 1) * ((n + 2 : ℝ) - rate * t) / rate

theorem hasDerivAt_expZeroDegreeRatio (rate t : ℝ) (n : ℕ) :
    HasDerivAt (expZeroDegreeRatio rate n)
      (((n + 2 : ℝ) * t ^ n * ((n + 1 : ℝ) - rate * t)) / rate) t := by
  apply ((((hasDerivAt_id t).pow (n + 1)).mul
    ((hasDerivAt_const t (n + 2 : ℝ)).sub
      ((hasDerivAt_id t).const_mul rate))).div_const rate).congr_deriv
  simp only [Nat.add_sub_cancel, Nat.cast_add, Nat.cast_one, id_eq, mul_one, zero_sub,
    Pi.sub_apply, Pi.pow_apply]
  rw [pow_succ]
  ring

theorem expZeroDegreeRatio_strictMonoOn {rate : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1) (n : ℕ) :
    StrictMonoOn (expZeroDegreeRatio rate n) (Ioo 0 1) := by
  apply strictMonoOn_of_deriv_pos (convex_Ioo 0 1)
  · exact fun t _ ↦ (hasDerivAt_expZeroDegreeRatio rate t n).continuousAt.continuousWithinAt
  · intro t ht
    have ht' : t ∈ Ioo (0 : ℝ) 1 := by simpa using ht
    rw [(hasDerivAt_expZeroDegreeRatio rate t n).deriv]
    have hz : rate * t < 1 := (mul_le_mul_of_nonneg_right rate_le ht'.1.le).trans_lt
      (by simpa using ht'.2)
    have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg _
    exact div_pos (mul_pos (mul_pos (by positivity) (pow_pos ht'.1 n))
      (by linarith)) rate_pos

theorem expNormalizedPower_zero_high_cross_neg {rate x y z : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1) (n : ℕ) :
    threePointCross (expNormalizedPower rate 0)
      (expNormalizedPower rate (n + 2)) x y z < 0 := by
  have h := expNormalizedPower_zero_cross_pos_of_ratio rate_pos hxy hyz
    (fun t ↦ -expNormalizedPower rate (n + 2) t)
    (fun t ↦ -expNormalizedPowerSuccDeriv rate (n + 1) t)
    (expZeroDegreeRatio rate n)
    (by unfold expNormalizedPower; fun_prop)
    (fun t _ ↦ (hasDerivAt_expNormalizedPower_succ rate (n + 1) t).neg)
    (by
      intro t _
      unfold expNormalizedPowerSuccDeriv expZeroDegreeRatio
      simp only [Nat.cast_add, Nat.cast_one]
      field_simp
      ring)
    (fun u hu v hv huv ↦ expZeroDegreeRatio_strictMonoOn rate_pos rate_le n
      ⟨hx.trans_lt hu.1, hu.2.trans_le hz⟩
      ⟨hx.trans_lt hv.1, hv.2.trans_le hz⟩ huv)
  rw [threePointCross_neg_right] at h
  linarith

noncomputable def expMonomialMinor (rate x y z : ℝ) (i j : ℕ) : ℝ :=
  Matrix.det (Matrix.of ![![Real.exp (rate * x), x ^ i, x ^ j],
    ![Real.exp (rate * y), y ^ i, y ^ j],
    ![Real.exp (rate * z), z ^ i, z ^ j]])

theorem expMonomialMinor_eq_normalized (rate x y z : ℝ) (i j : ℕ) :
    expMonomialMinor rate x y z i j =
      (Real.exp (rate * x) * Real.exp (rate * y) * Real.exp (rate * z)) *
        threePointCross (expNormalizedPower rate i) (expNormalizedPower rate j) x y z := by
  unfold expMonomialMinor
  rw [Matrix.det_fin_three]
  simp [threePointCross, expNormalizedPower, Real.exp_neg]
  field_simp
  ring

theorem expMonomialMinor_pos_of_pos_degrees {rate x y z : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1)
    {i j : ℕ} (hi : 0 < i) (hij : i < j) :
    0 < expMonomialMinor rate x y z i j := by
  obtain ⟨m, rfl⟩ : ∃ m, i = m + 1 := ⟨i - 1, by omega⟩
  obtain ⟨n, rfl⟩ : ∃ n, j = m + n + 2 := ⟨j - m - 2, by omega⟩
  rw [expMonomialMinor_eq_normalized]
  exact mul_pos (by positivity)
    (expNormalizedPower_succ_cross_pos rate_pos rate_le hx hxy hyz hz m n)

theorem expMonomialMinor_zero_one_pos {rate x y z : ℝ}
    (rate_pos : 0 < rate) (hxy : x < y) (hyz : y < z) :
    0 < expMonomialMinor rate x y z 0 1 := by
  rw [expMonomialMinor_eq_normalized]
  exact mul_pos (by positivity) (expNormalizedPower_zero_one_cross_pos rate_pos hxy hyz)

theorem expMonomialMinor_neg_of_zero_left {rate x y z : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1)
    {j : ℕ} (hj : 2 ≤ j) :
    expMonomialMinor rate x y z 0 j < 0 := by
  obtain ⟨n, rfl⟩ : ∃ n, j = n + 2 := ⟨j - 2, by omega⟩
  rw [expMonomialMinor_eq_normalized]
  exact mul_neg_of_pos_of_neg (by positivity)
    (expNormalizedPower_zero_high_cross_neg rate_pos rate_le hx hxy hyz hz n)

/-- The common sign depends on degrees, not on the ordered node triple. -/
noncomputable def expMonomialSign (i j : ℕ) : ℝ :=
  if i = 0 ∧ 2 ≤ j then -1 else 1

theorem expMonomialMinor_signed_pos {rate x y z : ℝ}
    (rate_pos : 0 < rate) (rate_le : rate ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1)
    {i j : ℕ} (hij : i < j) :
    0 < expMonomialSign i j * expMonomialMinor rate x y z i j := by
  by_cases hi : i = 0
  · subst i
    by_cases hj : j = 1
    · subst j
      simpa [expMonomialSign] using expMonomialMinor_zero_one_pos rate_pos hxy hyz
    · have hj' : 2 ≤ j := by omega
      have h := expMonomialMinor_neg_of_zero_left rate_pos rate_le hx hxy hyz hz hj'
      simpa [expMonomialSign, hj'] using neg_pos.mpr h
  · have h := expMonomialMinor_pos_of_pos_degrees rate_pos rate_le hx hxy hyz hz
      (Nat.pos_of_ne_zero hi) hij
    simpa [expMonomialSign, hi] using h

/-- Independently sampled one-pivot augmented minors have matching signs.
The rates may differ, but both must belong to the positive unit band. -/
theorem expMonomialMinor_matched_mul_pos {a b x y z u v w : ℝ}
    (a_pos : 0 < a) (a_le : a ≤ 1) (b_pos : 0 < b) (b_le : b ≤ 1)
    (hx : 0 ≤ x) (hxy : x < y) (hyz : y < z) (hz : z ≤ 1)
    (hu : 0 ≤ u) (huv : u < v) (hvw : v < w) (hw : w ≤ 1)
    {i j : ℕ} (hij : i < j) :
    0 < expMonomialMinor a x y z i j * expMonomialMinor b u v w i j := by
  have h := mul_pos
    (expMonomialMinor_signed_pos a_pos a_le hx hxy hyz hz hij)
    (expMonomialMinor_signed_pos b_pos b_le hu huv hvw hw hij)
  unfold expMonomialSign at h
  split_ifs at h <;> nlinarith

end GECPKernelStructure.Fermionic
