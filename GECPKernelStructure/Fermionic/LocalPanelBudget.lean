import GECPKernelStructure.Fermionic.DeterminantDecay
import GECPKernelStructure.Fermionic.GECPDeterminantDecay

namespace GECPKernelStructure.Fermionic

open Matrix

/-- Positive central panel followed by disjoint right-closed dyadic panels. -/
def positivePanelDomain : ℕ → ℝ → Prop
  | 0, omega => 0 ≤ omega ∧ omega ≤ 1
  | j + 1, omega => ((2 ^ j : ℕ) : ℝ) < omega ∧ omega ≤ (2 ^ (j + 1) : ℕ)

/-- A signed panel; zero belongs only to the positive central panel. -/
def fermionicPanelDomain : Bool → ℕ → ℝ → Prop
  | false, j, omega => positivePanelDomain j omega
  | true, j, omega => omega < 0 ∧ positivePanelDomain j (-omega)

/-- A single positive panel needs only `8p` features, regardless of its scale. -/
noncomputable def positivePanelTerms (p : ℕ) : ℕ → List SeparatedTerm
  | 0 => divideSeparatedRight (taylorSeparatedTerms (8 * p))
  | j + 1 => divideSeparatedRight
    (maskSeparatedTerms (taylorSeparatedTerms (8 * p))
      (fun t ↦ t * (2 ^ j : ℕ) ≤ p) (fun _ ↦ True))

theorem positivePanelTerms_length (p j : ℕ) : (positivePanelTerms p j).length = 8 * p := by
  cases j <;> simp [positivePanelTerms, divideSeparatedRight_length,
    maskSeparatedTerms_length, taylorSeparatedTerms_length]

theorem positivePanelDomain_bounds {j : ℕ} {omega : ℝ}
    (h : positivePanelDomain j omega) : 0 ≤ omega ∧ omega ≤ (2 ^ j : ℕ) := by
  cases j with
  | zero => simpa [positivePanelDomain] using h
  | succ j =>
    exact ⟨(show (0 : ℝ) ≤ (2 ^ j : ℕ) by positivity).trans h.1.le, h.2⟩

theorem eval_positivePanelTerms_eq_dyadic (p j : ℕ) (t omega : ℝ)
    (h : positivePanelDomain j omega) :
    evalSeparated (positivePanelTerms p j) t omega =
      evalSeparated (dyadicSeparatedTerms (8 * p) p j) t omega /
        (1 + Real.exp (-omega)) := by
  cases j with
  | zero => simp [positivePanelTerms, eval_divideSeparatedRight, dyadicSeparatedTerms]
  | succ j =>
    rw [positivePanelTerms, eval_divideSeparatedRight, eval_maskSeparatedTerms,
      eval_taylorSeparatedTerms, eval_dyadicSeparatedTerms_succ, if_neg h.1.not_ge]
    simp only [and_true]

theorem positivePanelTerms_error (p j : ℕ) (hp : 0 < p)
    {t omega : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1) (h : positivePanelDomain j omega) :
    |fermionicKernel t omega - evalSeparated (positivePanelTerms p j) t omega| ≤
      (1 / 2 : ℝ) ^ p := by
  rw [eval_positivePanelTerms_eq_dyadic p j t omega h]
  exact fermionicKernel_positive_separatedApprox p j hp ht0 ht1
    (positivePanelDomain_bounds h).1 (positivePanelDomain_bounds h).2

noncomputable def fermionicPanelTerms (p : ℕ) : Bool → ℕ → List SeparatedTerm
  | false, j => positivePanelTerms p j
  | true, j => reflectSeparatedTerms (positivePanelTerms p j)

theorem fermionicPanelTerms_length (p j : ℕ) (negative : Bool) :
    (fermionicPanelTerms p negative j).length = 8 * p := by
  cases negative <;> simp [fermionicPanelTerms, reflectSeparatedTerms_length,
    positivePanelTerms_length]

/-- The actual signed fermionic kernel has the same rank/error budget on
every panel, including the central and boundary cases. -/
theorem fermionicPanelTerms_error (p j : ℕ) (negative : Bool) (hp : 0 < p)
    {t omega : ℝ} (ht0 : 0 ≤ t) (ht1 : t ≤ 1)
    (h : fermionicPanelDomain negative j omega) :
    |fermionicKernel t omega - evalSeparated (fermionicPanelTerms p negative j) t omega| ≤
      (1 / 2 : ℝ) ^ p := by
  cases negative with
  | false => exact positivePanelTerms_error p j hp ht0 ht1 h
  | true =>
    have hh := positivePanelTerms_error p j hp (by linarith : 0 ≤ 1 - t)
      (by linarith : 1 - t ≤ 1) h.2
    rw [fermionicKernel_reflection] at hh
    simpa [fermionicPanelTerms, eval_reflectSeparatedTerms] using hh

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Initial determinant decay on one signed panel. The effective rank `8p`
and all constants are independent of the panel index and any global cutoff. -/
theorem fermionicKernel_panel_sample_det_le_two_pow (p j : ℕ) (negative : Bool)
    (hp : 0 < p) (rows columns : ι → ℝ)
    (ht0 : ∀ i, 0 ≤ rows i) (ht1 : ∀ i, rows i ≤ 1)
    (hpanel : ∀ k, fermionicPanelDomain negative j (columns k))
    (hcard : 8 * p ≤ Fintype.card ι) :
    |(sampleMatrix fermionicKernel rows columns).det| ≤
      (2 : ℝ) ^ Fintype.card ι * (Fintype.card ι).factorial *
        (1 / 2 : ℝ) ^ (p * (Fintype.card ι - 8 * p)) * 2 ^ (8 * p) := by
  let terms := fermionicPanelTerms p negative j
  let L := separatedLeft terms rows
  let R := separatedRight terms columns
  have hLR : L * R = sampleMatrix (evalSeparated terms) rows columns :=
    separatedLeft_mul_separatedRight terms rows columns
  have happrox : ∀ i k,
      |sampleMatrix fermionicKernel rows columns i k - (L * R) i k| ≤
        (1 / 2 : ℝ) ^ p := by
    intro i k
    rw [hLR]
    exact fermionicPanelTerms_error p j negative hp (ht0 i) (ht1 i) (hpanel k)
  have heps : (1 / 2 : ℝ) ^ p ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
  have hfactor : ∀ i k, |(L * R) i k| ≤ (2 : ℝ) := by
    intro i k
    rw [hLR]
    change |evalSeparated terms (rows i) (columns k)| ≤ 2
    have hpos := fermionicKernel_pos (rows i) (columns k)
    have hle := fermionicKernel_le_one (ht0 i) (ht1 i) (ω := columns k)
    have herr := fermionicPanelTerms_error p j negative hp (ht0 i) (ht1 i) (hpanel k)
    change |fermionicKernel (rows i) (columns k) - evalSeparated terms (rows i) (columns k)| ≤
      (1 / 2 : ℝ) ^ p at herr
    rcases abs_le.mp herr with ⟨hl, hu⟩
    apply abs_le.mpr
    constructor <;> linarith
  have h := Matrix.abs_det_le_two_pow_mul_of_factors_approx
    (sampleMatrix fermionicKernel rows columns) L R ((1 / 2 : ℝ) ^ p) 2
    happrox hfactor (by simpa [terms, fermionicPanelTerms_length] using hcard)
    (by positivity) heps (by norm_num)
  simpa [terms, fermionicPanelTerms_length, pow_mul] using h

/-- Odd approximation order makes the determinant exponent divisible by the
local sample order. There is no global-scale factor in this order. -/
theorem fermionicKernel_panel_minor_det_le_oddBudget (m j : ℕ) (negative : Bool)
    (rows columns : Fin (16 * (2 * m + 1)) → ℝ)
    (ht0 : ∀ i, 0 ≤ rows i) (ht1 : ∀ i, rows i ≤ 1)
    (hpanel : ∀ k, fermionicPanelDomain negative j (columns k)) :
    |(sampleMatrix fermionicKernel rows columns).det| ≤
      ((2 : ℝ) * (16 * (2 * m + 1) : ℕ) * (1 / 2 : ℝ) ^ m) ^
        (16 * (2 * m + 1)) := by
  let p := 2 * m + 1
  let n := 16 * p
  let r := 8 * p
  have hp : 0 < p := by omega
  have hcard : 8 * p ≤ Fintype.card (Fin n) := by simp only [Fintype.card_fin]; dsimp [n]; omega
  have hdet := fermionicKernel_panel_sample_det_le_two_pow p j negative hp rows columns
    ht0 ht1 hpanel hcard
  simp only [Fintype.card_fin] at hdet
  change |(sampleMatrix fermionicKernel rows columns).det| ≤
    (2 : ℝ) ^ n * n.factorial * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r at hdet
  change |(sampleMatrix fermionicKernel rows columns).det| ≤
    ((2 : ℝ) * n * (1 / 2 : ℝ) ^ m) ^ n
  have hfac : (n.factorial : ℝ) ≤ (n : ℝ) ^ n := by exact_mod_cast n.factorial_le_pow
  have hnsub : n - r = r := by dsimp [n, r]; omega
  have hexp : p * (n - r) = n * m + r := by rw [hnsub]; dsimp [n, r, p]; ring
  refine hdet.trans ?_
  calc
    (2 : ℝ) ^ n * n.factorial * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r ≤
        (2 : ℝ) ^ n * (n : ℝ) ^ n * (1 / 2 : ℝ) ^ (p * (n - r)) * 2 ^ r := by gcongr
    _ = ((2 : ℝ) * n * (1 / 2 : ℝ) ^ m) ^ n := by
      rw [hexp, pow_add]
      have hcancel : (1 / 2 : ℝ) ^ r * 2 ^ r = 1 := by rw [← mul_pow]; norm_num
      have hhalf : (1 / 2 : ℝ) ^ (n * m) = ((1 / 2 : ℝ) ^ m) ^ n := by
        rw [Nat.mul_comm n m, pow_mul]
      calc
        _ = ((2 : ℝ) ^ n * (n : ℝ) ^ n * (1 / 2 : ℝ) ^ (n * m)) *
            ((1 / 2 : ℝ) ^ r * 2 ^ r) := by ring
        _ = (2 : ℝ) ^ n * (n : ℝ) ^ n * (1 / 2 : ℝ) ^ (n * m) := by rw [hcancel, mul_one]
        _ = _ := by rw [hhalf, ← mul_pow, ← mul_pow]

/-- The local odd-block coefficient reaches dyadic accuracy with a parameter
choice independent of panel scale. -/
theorem panel_oddBudget_factor_le_accuracy (q : ℕ) :
    (2 : ℝ) * (16 * (33 + 4 * q) : ℕ) * (1 / 2 : ℝ) ^ (16 + 2 * q) ≤
      (1 / 2 : ℝ) ^ q := by
  have hscale : 64 * (2 * (16 + 2 * q) + 1) * (0 + 1) * 2 ^ q ≤
      2 ^ (16 + 2 * q) := by simpa using oddBlock_scale_accuracy 0 q
  have hh := oddBlock_factor_le_accuracy (16 + 2 * q) 0 q hscale
  have hm : 2 * (16 + 2 * q) + 1 = 33 + 4 * q := by omega
  rw [hm] at hh
  have hn : (16 * (33 + 4 * q) : ℕ) ≤ 32 * (33 + 4 * q) * (0 + 1) := by omega
  have hnr : ((16 * (33 + 4 * q) : ℕ) : ℝ) ≤
      ((32 * (33 + 4 * q) * (0 + 1) : ℕ) : ℝ) := by exact_mod_cast hn
  calc
    _ ≤ (2 : ℝ) * (32 * (33 + 4 * q) * (0 + 1) : ℕ) *
        (1 / 2 : ℝ) ^ (16 + 2 * q) :=
      mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hnr (by norm_num)) (by positivity)
    _ ≤ _ := hh

/-- A scale-free initial local determinant budget at every dyadic accuracy.
All rows/columns may be unordered or repeated; only panel membership and
the physical time interval are required. -/
theorem fermionicKernel_panel_minor_det_le_dyadicBudget (q j : ℕ) (negative : Bool)
    (rows columns : Fin (16 * (33 + 4 * q)) → ℝ)
    (ht0 : ∀ i, 0 ≤ rows i) (ht1 : ∀ i, rows i ≤ 1)
    (hpanel : ∀ k, fermionicPanelDomain negative j (columns k)) :
    |(sampleMatrix fermionicKernel rows columns).det| ≤
      ((1 / 2 : ℝ) ^ q) ^ (16 * (33 + 4 * q)) := by
  have hm : 2 * (16 + 2 * q) + 1 = 33 + 4 * q := by omega
  have hbudget := fermionicKernel_panel_minor_det_le_oddBudget (16 + 2 * q) j negative
  rw [hm] at hbudget
  exact (hbudget rows columns ht0 ht1 hpanel).trans
    (pow_le_pow_left₀ (by positivity) (panel_oddBudget_factor_le_accuracy q) (16 * (33 + 4 * q)))

end GECPKernelStructure.Fermionic
