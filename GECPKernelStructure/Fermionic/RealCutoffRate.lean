import GECPKernelStructure.Fermionic.LocalPanelRate
import GECPKernelStructure.Fermionic.CentralPrefixRun
import Mathlib.Analysis.SpecialFunctions.Log.Base

namespace GECPKernelStructure.GECP

/-- Initialization and successful exact updates already supply every finite run. -/
theorem realized_of_initial_updates {K : Kernel ℝ ℝ ℝ}
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = K)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i)) :
    ∀ i, ∃ run : Run K, run.finalResidual = residual i := by
  intro i
  refine ⟨Run.transport initial (Run.ofResidualSequenceFrom residual rows columns pivot_ne updates 0 i), ?_⟩
  simp only [Run.finalResidual_transport, Run.finalResidual_ofResidualSequenceFrom, zero_add]

end GECPKernelStructure.GECP

namespace GECPKernelStructure.Fermionic

open GECP Matrix

/-- The panel count uses a dyadic *cover*, while maxima stay on the actual domain. -/
theorem frequencyReflectedFermionicKernel_gecp_error_le_realCutoffBlock
    (Λ : ℝ) (s q : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ))
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = frequencyReflectedFermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual (32 * (s + 1) * (33 + 4 * q)) x y| ≤ (1 / 2 : ℝ) ^ q := by
  classical
  have realized := realized_of_initial_updates residual rows columns initial pivot_ne updates
  let n := 16 * (33 + 4 * q)
  let N := 32 * (s + 1) * (33 + 4 * q)
  let epsilon : ℝ := (1 / 2 : ℝ) ^ q
  let D : ℝ → Prop := fun t ↦ 0 ≤ t ∧ t ≤ 1
  let C : ℝ → Prop := fun omega ↦ -Λ ≤ omega ∧ omega ≤ Λ
  change |residual N x y| ≤ epsilon
  by_contra hnot
  have hlarge : epsilon < |residual N x y| := lt_of_not_ge hnot
  have antitone := strictSignRegular_gecp_pivotMagnitude_antitone
    frequencyReflectedFermionicKernel_strictSignRegular D C residual rows columns
    realized pivot_ne updates complete
  have hfinal : epsilon < |residual N (rows N) (columns N)| :=
    hlarge.trans_le ((complete N).2.2 x y hx hy)
  have hall (i : ℕ) (hi : i < N) : epsilon < |residual i (rows i) (columns i)| :=
    hfinal.trans_le (antitone hi.le)
  have hneg (i : ℕ) : -((2 ^ s : ℕ) : ℝ) ≤ -columns i ∧ -columns i ≤ (2 ^ s : ℕ) := by
    obtain ⟨hl, hu⟩ := (complete i).2.1
    constructor <;> linarith
  let label : Fin N → Bool × Fin (s + 1) :=
    fun i ↦ Classical.choose (fermionicPanelDomain_cover s (hneg i.val))
  have label_valid (i : Fin N) :
      fermionicPanelDomain (label i).1 (label i).2.val (-columns i.val) :=
    Classical.choose_spec (fermionicPanelDomain_cover s (hneg i.val))
  have hsize : Fintype.card (Bool × Fin (s + 1)) * n = Fintype.card (Fin N) := by
    rw [fermionicPanel_count, Fintype.card_fin]
    dsimp [n, N]
    ring
  obtain ⟨panel, hfiber⟩ := Fintype.exists_le_card_fiber_of_mul_le_card label hsize.le
  let fiber : Finset (Fin N) := Finset.univ.filter (fun i ↦ label i = panel)
  let S : Finset ℕ := fiber.map Fin.valEmbedding
  have hSN : S ⊆ Finset.range N := by
    intro i hi
    obtain ⟨j, _, rfl⟩ := Finset.mem_map.mp hi
    exact Finset.mem_range.mpr j.isLt
  have hScard : n ≤ S.card := by simpa [S, fiber] using hfiber
  let B : ℝ → Prop := fun omega ↦ fermionicPanelDomain panel.1 panel.2.val (-omega)
  have hSB : ∀ i ∈ S, B (columns i) := by
    intro i hi
    obtain ⟨j, hj, rfl⟩ := Finset.mem_map.mp hi
    have heq : label j = panel := (Finset.mem_filter.mp hj).2
    have hh := label_valid j
    change fermionicPanelDomain panel.1 panel.2.val (-columns j.val)
    simpa only [heq] using hh
  have budget : ∀ rs cs : Fin n → ℝ, (∀ i, D (rs i)) → (∀ j, B (cs j)) →
      |(minorMatrix frequencyReflectedFermionicKernel rs cs).det| ≤ epsilon ^ n := by
    intro rs cs hD hB
    change |(sampleMatrix fermionicKernel rs (fun j ↦ -cs j)).det| ≤
      ((1 / 2 : ℝ) ^ q) ^ (16 * (33 + 4 * q))
    exact fermionicKernel_panel_minor_det_le_dyadicBudget q panel.2.val panel.1 rs
      (fun j ↦ -cs j) (fun i ↦ (hD i).1) (fun i ↦ (hD i).2) hB
  have hsmall := panel_pivot_count_lt_of_initialBudget
    frequencyReflectedFermionicKernel_strictSignRegular residual rows columns initial
    realized pivot_ne updates D B (fun i ↦ (complete i).1)
    (fun i z hz ↦ (complete i).2.2 z (columns i) hz (complete i).2.1)
    n (by dsimp [n]; positivity) epsilon (by dsimp [epsilon]; positivity)
    budget N S hSN hSB (fun i hi ↦ hall i (Finset.mem_range.mp (hSN hi)))
  exact (Nat.not_lt_of_ge hScard) hsmall

/-- Actual-cutoff dyadic accuracy at every rank beyond the explicit block.
There is no independent realized-run or residual-budget hypothesis. -/
theorem fermionicKernel_gecp_error_le_realCutoffBlock
    (Λ : ℝ) (s q k : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ))
    (hcount : 32 * (s + 1) * (33 + 4 * q) ≤ k)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual k x y| ≤ (1 / 2 : ℝ) ^ q := by
  let R : ℕ → Kernel ℝ ℝ ℝ := fun i t omega ↦ residual i t (-omega)
  let cs : ℕ → ℝ := fun i ↦ -columns i
  have hinit : R 0 = frequencyReflectedFermionicKernel := by
    funext t omega
    change residual 0 t (-omega) = fermionicKernel t (-omega)
    rw [initial]
  have hp : ∀ i, R i (rows i) (cs i) ≠ 0 := by
    intro i
    simpa [R, cs] using pivot_ne i
  have hu : ∀ i, R (i + 1) = residualUpdate (R i) (rows i) (cs i) (hp i) := by
    intro i
    funext t omega
    have hh := congrFun (congrFun (updates i) t) (-omega)
    simpa only [R, cs, residualUpdate, neg_neg] using hh
  have hc : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (R i) (rows i) (cs i) := by
    intro i
    obtain ⟨hrow, hcol, hmax⟩ := complete i
    refine ⟨hrow, ?_, ?_⟩
    · dsimp [cs]
      constructor <;> linarith
    · intro t omega ht homega
      have hneg : -Λ ≤ -omega ∧ -omega ≤ Λ := by constructor <;> linarith
      simpa only [R, cs, neg_neg] using hmax t (-omega) ht hneg
  let N := 32 * (s + 1) * (33 + 4 * q)
  have hb := frequencyReflectedFermionicKernel_gecp_error_le_realCutoffBlock Λ s q hΛ
    R rows cs hinit hp hu hc (rows N) (cs N) (hc N).1 (hc N).2.1
  have realized := realized_of_initial_updates R rows cs hinit hp hu
  have ha := strictSignRegular_gecp_pivotMagnitude_antitone
    frequencyReflectedFermionicKernel_strictSignRegular
    (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1) (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ)
    R rows cs realized hp hu hc
  have hy' : -Λ ≤ -y ∧ -y ≤ Λ := by constructor <;> linarith
  have hh := ((hc k).2.2 x (-y) hx hy').trans ((ha hcount).trans hb)
  simpa only [R, neg_neg] using hh

/-- Convert the explicit natural block bounds to an every-rank exponential
bound. The constants cover the initial ranks as well as all block endpoints. -/
theorem localPanelBlocks_le_exponential (s k : ℕ) (E : ℝ) (hone : E ≤ 1)
    (hblock : ∀ q : ℕ, 32 * (s + 1) * (33 + 4 * q) ≤ k → E ≤ (1 / 2 : ℝ) ^ q) :
    E ≤ 1024 * Real.exp (-(Real.log 2) * k / (128 * (s + 1))) := by
  let j := k / (128 * (s + 1))
  let q := j - 9
  have hqj : j ≤ q + 9 := by dsimp [q]; omega
  have hE : E ≤ (1 / 2 : ℝ) ^ q := by
    by_cases hj : 9 ≤ j
    · apply hblock
      have heq : q + 9 = j := Nat.sub_add_cancel hj
      have hfirst : 32 * (s + 1) * (33 + 4 * q) ≤ 128 * (s + 1) * j := by
        rw [← heq]
        nlinarith
      exact hfirst.trans (Nat.mul_div_le k (128 * (s + 1)))
    · have hq : q = 0 := by dsimp [q]; omega
      simpa only [hq, pow_zero] using hone
  have hk : (k : ℝ) < 128 * ((s : ℝ) + 1) * ((q : ℝ) + 10) := by
    have hd := Nat.lt_mul_div_succ k (by positivity : 0 < 128 * (s + 1))
    have hn : k < 128 * (s + 1) * (q + 10) :=
      hd.trans_le (Nat.mul_le_mul_left _ (by dsimp [j] at hqj ⊢; omega))
    exact_mod_cast hn
  have hlog : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hden : (0 : ℝ) < 128 * ((s : ℝ) + 1) := by positivity
  have hexp : -(Real.log 2) * (q : ℝ) ≤
      10 * Real.log 2 - Real.log 2 * k / (128 * (s + 1)) := by
    have hh : (k : ℝ) / (128 * ((s : ℝ) + 1)) < (q : ℝ) + 10 :=
      (div_lt_iff₀ hden).mpr (by nlinarith [hk])
    have hm := mul_le_mul_of_nonneg_left hh.le hlog.le
    rw [← mul_div_assoc] at hm
    nlinarith
  have hhalf : (1 / 2 : ℝ) ^ q = Real.exp (-(Real.log 2) * (q : ℝ)) := by
    rw [mul_comm, Real.exp_nat_mul, Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 2)]
    norm_num
  have hconstant : Real.exp (10 * Real.log 2) = (1024 : ℝ) := by
    rw [show (10 : ℝ) = ((10 : ℕ) : ℝ) by norm_num, Real.exp_nat_mul,
      Real.exp_log (by norm_num : (0 : ℝ) < 2)]
    norm_num
  calc
    E ≤ Real.exp (-(Real.log 2) * (q : ℝ)) := by rwa [← hhalf]
    _ ≤ Real.exp (10 * Real.log 2 - Real.log 2 * k / (128 * (s + 1))) :=
      Real.exp_le_exp.mpr hexp
    _ = _ := by rw [sub_eq_add_neg, Real.exp_add, hconstant]; congr 2; ring

/-- Every-rank exponential decay, with the dyadic cover only in the rate. -/
theorem fermionicKernel_gecp_error_le_coverExponential
    (Λ : ℝ) (s k : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ))
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual k x y| ≤ 1024 * Real.exp (-(Real.log 2) * k / (128 * (s + 1))) := by
  have realized := realized_of_initial_updates residual rows columns initial pivot_ne updates
  have hone := (fermionicKernel_gecp_error_le_cutoffCorner residual rows columns initial
    realized pivot_ne updates complete k x y hx hy).trans
      (fermionicKernel_le_one (t := 0) (ω := Λ) (by norm_num) (by norm_num))
  exact localPanelBlocks_le_exponential s k _ hone (fun q hq ↦
    fermionicKernel_gecp_error_le_realCutoffBlock Λ s q k hΛ hq residual rows columns
      initial pivot_ne updates complete x y hx hy)

/-- A ceiling-logarithmic cover, with a bound suitable for the large-cutoff
`log(1+Lambda)` rate. No numerical rounding determines the cover. -/
theorem exists_logarithmic_cutoff_cover {Λ : ℝ} (hΛ : 1 ≤ Λ) :
    ∃ s : ℕ, Λ ≤ (2 ^ s : ℕ) ∧
      ((s : ℝ) + 1) * Real.log 2 ≤ 2 * Real.log (1 + Λ) := by
  have hl : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hpos : 0 < Λ := by linarith
  have hnonneg : 0 ≤ Real.log Λ / Real.log 2 :=
    div_nonneg (Real.log_nonneg hΛ) hl.le
  let s := ⌈Real.log Λ / Real.log 2⌉₊
  refine ⟨s, ?_, ?_⟩
  · have hh : Real.log Λ ≤ (s : ℝ) * Real.log 2 :=
      (div_le_iff₀ hl).mp (Nat.le_ceil _)
    have he := Real.exp_le_exp.mpr hh
    simpa only [Real.exp_log hpos, Real.exp_nat_mul,
      Real.exp_log (by norm_num : (0 : ℝ) < 2), Nat.cast_pow, Nat.cast_ofNat] using he
  · have hs : (s : ℝ) < Real.log Λ / Real.log 2 + 1 := Nat.ceil_lt_add_one hnonneg
    have hs' : ((s : ℝ) + 1) * Real.log 2 ≤ Real.log Λ + 2 * Real.log 2 := by
      have hh := (lt_div_iff₀ hl).mp (by linarith : (s : ℝ) - 1 < Real.log Λ / Real.log 2)
      nlinarith
    have hpoly : 4 * Λ ≤ (1 + Λ) ^ 2 := by nlinarith [sq_nonneg (Λ - 1)]
    have hlogs := Real.log_le_log (by positivity : (0 : ℝ) < 4 * Λ) hpoly
    rw [Real.log_mul (by norm_num : (4 : ℝ) ≠ 0) hpos.ne', Real.log_pow] at hlogs
    have hfour : Real.log 4 = 2 * Real.log 2 := by
      convert Real.log_pow (2 : ℝ) 2 using 1 <;> norm_num
    rw [hfour] at hlogs
    norm_num only [Nat.cast_ofNat] at hlogs
    exact hs'.trans (by linarith)

/-- A logarithmic cover converts the explicit rate without changing its pivots. -/
theorem localPanel_exponent_le_logarithmic (s k : ℕ) (L : ℝ) (hL : 0 < L)
    (hcover : ((s : ℝ) + 1) * Real.log 2 ≤ 2 * L) :
    -(Real.log 2) * k / (128 * (s + 1)) ≤
      -((Real.log 2) ^ 2 / 256) * k / L := by
  have hl : 0 ≤ Real.log 2 := (Real.log_pos (by norm_num : (1 : ℝ) < 2)).le
  have hd : (0 : ℝ) < 128 * ((s : ℝ) + 1) := by positivity
  apply (div_le_div_iff₀ hd hL).mpr
  have hh := mul_le_mul_of_nonneg_right hcover (mul_nonneg hl (Nat.cast_nonneg k : (0 : ℝ) ≤ k))
  nlinarith

/-- Explicit G1-shaped rate on its large-cutoff range, with universal
`C=1024` and `c=(log 2)^2/256`. -/
theorem fermionicKernel_gecp_error_le_logarithmic
    (Λ : ℝ) (hΛ : 1 ≤ Λ) (k : ℕ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual k x y| ≤ 1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * k / Real.log (1 + Λ)) := by
  obtain ⟨s, hs, hlog⟩ := exists_logarithmic_cutoff_cover hΛ
  have hh := fermionicKernel_gecp_error_le_coverExponential Λ s k hs residual rows columns
    initial pivot_ne updates complete x y hx hy
  refine hh.trans (mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by norm_num))
  exact localPanel_exponent_le_logarithmic s k _ (Real.log_pos (by linarith)) hlog

/-- Every positive cutoff has a logarithmic cover with `log(2+Lambda)`. -/
theorem exists_allCutoff_logarithmic_cover {Λ : ℝ} (hΛ : 0 < Λ) :
    ∃ s : ℕ, Λ ≤ (2 ^ s : ℕ) ∧
      ((s : ℝ) + 1) * Real.log 2 ≤ 2 * Real.log (2 + Λ) := by
  obtain ⟨s, hs, hlog⟩ := exists_logarithmic_cutoff_cover (le_max_left 1 Λ)
  refine ⟨s, (le_max_right 1 Λ).trans hs, hlog.trans ?_⟩
  have hm : 1 + max 1 Λ ≤ 2 + Λ := by
    rw [max_def]
    split_ifs <;> linarith
  exact mul_le_mul_of_nonneg_left (Real.log_le_log (by positivity) hm) (by norm_num)

/-- A cutoff-uniform version valid down to arbitrarily small positive cutoffs.
The denominator is deliberately `log(2+Lambda)`, not a rewritten SPEC claim. -/
theorem fermionicKernel_gecp_error_le_allCutoff_logarithmic
    (Λ : ℝ) (hΛ : 0 < Λ) (k : ℕ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual k x y| ≤ 1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * k / Real.log (2 + Λ)) := by
  obtain ⟨s, hs, hlog⟩ := exists_allCutoff_logarithmic_cover hΛ
  have hh := fermionicKernel_gecp_error_le_coverExponential Λ s k hs residual rows columns
    initial pivot_ne updates complete x y hx hy
  refine hh.trans (mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by norm_num))
  exact localPanel_exponent_le_logarithmic s k _ (Real.log_pos (by linarith)) hlog

/-- An explicit real rank threshold for any positive requested tolerance. -/
theorem logarithmic_envelope_le_tolerance (L epsilon : ℝ) (hL : 0 < L)
    (hepsilon : 0 < epsilon) (k : ℕ)
    (hcount : 256 / (Real.log 2) ^ 2 * L * Real.log (1024 / epsilon) ≤ k) :
    1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * k / L) ≤ epsilon := by
  have hl : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have ha : (0 : ℝ) < (Real.log 2) ^ 2 := sq_pos_of_pos hl
  have hm := mul_le_mul_of_nonneg_left hcount
    (show 0 ≤ (Real.log 2) ^ 2 / (256 * L) by positivity)
  have hid : (Real.log 2) ^ 2 / (256 * L) *
      (256 / (Real.log 2) ^ 2 * L * Real.log (1024 / epsilon)) =
      Real.log (1024 / epsilon) := by field_simp
  rw [hid] at hm
  have hexp : -((Real.log 2) ^ 2 / 256) * (k : ℝ) / L ≤
      -Real.log (1024 / epsilon) := by
    have heq : (Real.log 2) ^ 2 / (256 * L) * (k : ℝ) =
        ((Real.log 2) ^ 2 / 256) * k / L := by ring
    rw [heq] at hm
    simpa only [neg_div, neg_mul] using neg_le_neg hm
  have he := Real.exp_le_exp.mpr hexp
  have heval : 1024 * Real.exp (-Real.log (1024 / epsilon)) = epsilon := by
    rw [Real.exp_neg, Real.exp_log (by positivity : (0 : ℝ) < 1024 / epsilon)]
    field_simp
  exact (mul_le_mul_of_nonneg_left he (by norm_num : (0 : ℝ) ≤ 1024)).trans_eq heval

/-- Real accuracy with an explicit ceiling rank, on all positive cutoffs.
Using `log(2+Lambda)` states the small-cutoff range without ambiguity. -/
theorem fermionicKernel_gecp_error_le_realTolerance
    (Λ epsilon : ℝ) (hΛ : 0 < Λ) (hepsilon : 0 < epsilon) (k : ℕ)
    (hcount : ⌈256 / (Real.log 2) ^ 2 * Real.log (2 + Λ) * Real.log (1024 / epsilon)⌉₊ ≤ k)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |residual k x y| ≤ epsilon := by
  have hh := fermionicKernel_gecp_error_le_allCutoff_logarithmic Λ hΛ k residual rows columns
    initial pivot_ne updates complete x y hx hy
  exact hh.trans (logarithmic_envelope_le_tolerance _ epsilon (Real.log_pos (by linarith))
    hepsilon k ((Nat.le_ceil _).trans (Nat.cast_le.mpr hcount)))

/-- The literal G1 denominator cannot hold uniformly as `Lambda` tends to
zero: an actual first complete pivot already contradicts any proposed constants. -/
theorem fermionicKernel_smallCutoff_logarithmic_obstruction
    (C c : ℝ) (hC : 0 < C) (hc : 0 < c) :
    ∃ Λ : ℝ, 0 < Λ ∧ Λ ≤ 1 ∧
      CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
        (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) fermionicKernel 0 Λ ∧
      C * Real.exp (-c / Real.log (1 + Λ)) <
        |residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' 1 (-Λ)| := by
  let Λ := min 1 (c ^ 2 / (16 * C))
  have hp : 0 < Λ := lt_min (by norm_num) (by positivity)
  have hsmall : Λ ≤ 1 := min_le_left _ _
  have hchoice : 16 * C * Λ ≤ c ^ 2 := by
    have hh := (le_div_iff₀ (by positivity : (0 : ℝ) < 16 * C)).mp (min_le_right 1 (c ^ 2 / (16 * C)))
    nlinarith
  let L := Real.log (1 + Λ)
  have hL : 0 < L := Real.log_pos (by linarith)
  have hLupper : L ≤ Λ := by
    have hh := Real.log_le_sub_one_of_pos (by linarith : 0 < 1 + Λ)
    dsimp [L]
    linarith
  have hquadratic : (c / L) ^ 2 / 2 ≤ Real.exp (c / L) := by
    have hh := Real.quadratic_le_exp_of_nonneg (by positivity : 0 ≤ c / L)
    have hnonneg : 0 ≤ c / L := by positivity
    linarith
  have hupper : C * Real.exp (-c / L) ≤ 2 * C * L ^ 2 / c ^ 2 := by
    rw [neg_div, Real.exp_neg, ← div_eq_mul_inv]
    apply (div_le_iff₀ (Real.exp_pos _)).mpr
    calc
      C = (2 * C * L ^ 2 / c ^ 2) * ((c / L) ^ 2 / 2) := by field_simp
      _ ≤ _ := mul_le_mul_of_nonneg_left hquadratic (by positivity)
  have hupper' : C * Real.exp (-c / L) ≤ Λ / 4 := by
    refine hupper.trans ?_
    apply (div_le_iff₀ (sq_pos_of_pos hc)).mpr
    have hsq : L ^ 2 ≤ Λ ^ 2 := pow_le_pow_left₀ hL.le hLupper 2
    have hm := mul_le_mul_of_nonneg_right hchoice hp.le
    nlinarith [mul_le_mul_of_nonneg_left hsq hC.le]
  have hinv : Real.exp (-Λ) ≤ 1 / (1 + Λ) := by
    rw [Real.exp_neg, ← one_div]
    apply one_div_le_one_div_of_le (by linarith)
    linarith [Real.add_one_le_exp Λ]
  have hlower : Λ / 2 ≤ 1 - Real.exp (-Λ) := by
    have hf : Λ / 2 ≤ Λ / (1 + Λ) :=
      (div_le_div_iff₀ (by norm_num : (0 : ℝ) < 2) (by linarith)).mpr (by nlinarith)
    have heq : Λ / (1 + Λ) = 1 - 1 / (1 + Λ) := by field_simp; ring
    rw [heq] at hf
    linarith
  refine ⟨Λ, hp, hsmall, ?_, ?_⟩
  · refine ⟨⟨le_rfl, zero_le_one⟩, ⟨by linarith, le_rfl⟩, ?_⟩
    intro t omega ht homega
    rw [abs_of_pos (fermionicKernel_pos t omega), abs_of_pos (fermionicKernel_pos 0 Λ)]
    exact fermionicKernel_le_cutoffCorner ht.1 ht.2 homega.1 homega.2
  · rw [fermionicKernel_firstPivot_reflected_residual, abs_of_nonneg (by linarith)]
    change C * Real.exp (-c / L) < 1 - Real.exp (-Λ)
    linarith

end GECPKernelStructure.Fermionic
