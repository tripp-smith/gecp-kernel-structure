import GECPKernelStructure.GECP.LocalMinorWitness
import GECPKernelStructure.Fermionic.LocalPanelBudget
import GECPKernelStructure.Fermionic.CommonResidualGauge
import GECPKernelStructure.Fermionic.GeometricMean
import Mathlib.Combinatorics.Pigeonhole

namespace GECPKernelStructure.Fermionic

open GECP Matrix

/-- The central/dyadic positive panels cover the entire finite cutoff interval. -/
theorem positivePanelDomain_cover (s : ℕ) {omega : ℝ}
    (h0 : 0 ≤ omega) (hupper : omega ≤ (2 ^ s : ℕ)) :
    ∃ j : Fin (s + 1), positivePanelDomain j.val omega := by
  induction s with
  | zero =>
    refine ⟨0, ?_⟩
    exact ⟨h0, by simpa using hupper⟩
  | succ s ih =>
    by_cases hlow : omega ≤ (2 ^ s : ℕ)
    · obtain ⟨j, hj⟩ := ih hlow
      exact ⟨j.castSucc, hj⟩
    · exact ⟨Fin.last (s + 1), lt_of_not_ge hlow, hupper⟩

theorem positivePanelDomain_unique {j k : ℕ} {omega : ℝ}
    (hj : positivePanelDomain j omega) (hk : positivePanelDomain k omega) : j = k := by
  have not_lt {a b : ℕ} (ha : positivePanelDomain a omega)
      (hb : positivePanelDomain b omega) : ¬a < b := by
    intro hab
    cases b with
    | zero => omega
    | succ b =>
      have hpow : ((2 ^ a : ℕ) : ℝ) ≤ ((2 ^ b : ℕ) : ℝ) := by
        exact_mod_cast (pow_le_pow_right₀ (by norm_num : 1 ≤ (2 : ℕ)) (by omega : a ≤ b))
      have hu := (positivePanelDomain_bounds ha).2
      exact (not_lt_of_ge (hu.trans hpow)) hb.1
  exact le_antisymm (Nat.le_of_not_gt (not_lt hk hj)) (Nat.le_of_not_gt (not_lt hj hk))

/-- Every signed-cutoff frequency belongs to one of `2(s+1)` signed panels. -/
theorem fermionicPanelDomain_cover (s : ℕ) {omega : ℝ}
    (h : -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ)) :
    ∃ panel : Bool × Fin (s + 1), fermionicPanelDomain panel.1 panel.2.val omega := by
  by_cases hn : omega < 0
  · obtain ⟨j, hj⟩ := positivePanelDomain_cover s (by linarith : 0 ≤ -omega)
      (by linarith : -omega ≤ (2 ^ s : ℕ))
    exact ⟨(true, j), hn, hj⟩
  · obtain ⟨j, hj⟩ := positivePanelDomain_cover s (le_of_not_gt hn) h.2
    exact ⟨(false, j), hj⟩

theorem fermionicPanelDomain_unique {s : ℕ} {omega : ℝ}
    {a b : Bool × Fin (s + 1)} (ha : fermionicPanelDomain a.1 a.2.val omega)
    (hb : fermionicPanelDomain b.1 b.2.val omega) : a = b := by
  rcases a with ⟨a, j⟩
  rcases b with ⟨b, k⟩
  cases a <;> cases b
  · exact Prod.ext rfl (Fin.ext (positivePanelDomain_unique ha hb))
  · exact ((positivePanelDomain_bounds ha).1.not_gt hb.1).elim
  · exact ((positivePanelDomain_bounds hb).1.not_gt ha.1).elim
  · exact Prod.ext rfl (Fin.ext (positivePanelDomain_unique ha.2 hb.2))

theorem fermionicPanel_count (s : ℕ) : Fintype.card (Bool × Fin (s + 1)) = 2 * (s + 1) := by
  simp

/-- The positive-minor presentation reaches dyadic accuracy in a number of
actual complete pivots linear in both the dyadic cutoff and accuracy orders. -/
theorem frequencyReflectedFermionicKernel_gecp_error_le_localPanelBlock
    (s q : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = frequencyReflectedFermionicKernel)
    (realized : ∀ i, ∃ run : Run frequencyReflectedFermionicKernel, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ))
      (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (s + 1) * (33 + 4 * q)) x y| ≤ (1 / 2 : ℝ) ^ q := by
  classical
  let n := 16 * (33 + 4 * q)
  let N := 32 * (s + 1) * (33 + 4 * q)
  let epsilon : ℝ := (1 / 2 : ℝ) ^ q
  let D : ℝ → Prop := fun t ↦ 0 ≤ t ∧ t ≤ 1
  let C : ℝ → Prop := fun omega ↦
    -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ)
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

/-- Explicit global dyadic-accuracy rate for the original fermionic kernel.
The algorithm hypotheses are only initialization, exact successful updates,
realized finite prefixes, and complete maxima on the physical rectangle. -/
theorem fermionicKernel_gecp_error_le_localPanelBlock
    (s q : ℕ) (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (realized : ∀ i, ∃ run : Run fermionicKernel, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ))
      (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual (32 * (s + 1) * (33 + 4 * q)) x y| ≤ (1 / 2 : ℝ) ^ q := by
  let R : ℕ → Kernel ℝ ℝ ℝ := fun i t omega ↦ residual i t (-omega)
  let cs : ℕ → ℝ := fun i ↦ -columns i
  have hinit : R 0 = frequencyReflectedFermionicKernel := by
    funext t omega
    change residual 0 t (-omega) = fermionicKernel t (-omega)
    rw [initial]
  have hrealized : ∀ i, ∃ run : Run frequencyReflectedFermionicKernel, run.finalResidual = R i := by
    intro i
    obtain ⟨run, hr⟩ := realized i
    refine ⟨run.negateColumns, ?_⟩
    refine (run.finalResidual_negateColumns).trans ?_
    funext t omega
    change run.finalResidual t (-omega) = residual i t (-omega)
    rw [hr]
  have hp : ∀ i, R i (rows i) (cs i) ≠ 0 := by
    intro i
    simpa [R, cs] using pivot_ne i
  have hu : ∀ i, R (i + 1) = residualUpdate (R i) (rows i) (cs i) (hp i) := by
    intro i
    funext t omega
    have hh := congrFun (congrFun (updates i) t) (-omega)
    simpa only [R, cs, residualUpdate, neg_neg] using hh
  have hc : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ))
      (R i) (rows i) (cs i) := by
    intro i
    obtain ⟨hrow, hcol, hmax⟩ := complete i
    refine ⟨hrow, ?_, ?_⟩
    · dsimp [cs]
      constructor <;> linarith
    · intro t omega ht homega
      have hneg : -((2 ^ s : ℕ) : ℝ) ≤ -omega ∧ -omega ≤ (2 ^ s : ℕ) := by
        constructor <;> linarith
      simpa only [R, cs, neg_neg] using hmax t (-omega) ht hneg
  have hy' : -((2 ^ s : ℕ) : ℝ) ≤ -y ∧ -y ≤ (2 ^ s : ℕ) := by
    constructor <;> linarith
  have hh := frequencyReflectedFermionicKernel_gecp_error_le_localPanelBlock s q R rows cs
    hinit hrealized hp hu hc x (-y) hx hy'
  simpa only [R, neg_neg] using hh

/-- The dyadic-accuracy guarantee holds at every later rank, not just at the
named block endpoint. -/
theorem fermionicKernel_gecp_error_le_of_localPanelBlock
    (s q k : ℕ) (hcount : 32 * (s + 1) * (33 + 4 * q) ≤ k)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 = fermionicKernel)
    (realized : ∀ i, ∃ run : Run fermionicKernel, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (rows i) (columns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) = residualUpdate (residual i) (rows i) (columns i) (pivot_ne i))
    (complete : ∀ i, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ))
      (residual i) (rows i) (columns i))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1)
    (hy : -((2 ^ s : ℕ) : ℝ) ≤ y ∧ y ≤ (2 ^ s : ℕ)) :
    |residual k x y| ≤ (1 / 2 : ℝ) ^ q := by
  let N := 32 * (s + 1) * (33 + 4 * q)
  have ha := strictSignRegular_gecp_pivotMagnitude_antitone fermionicKernel_strictSignRegular
    (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
    (fun omega : ℝ ↦ -((2 ^ s : ℕ) : ℝ) ≤ omega ∧ omega ≤ (2 ^ s : ℕ))
    residual rows columns realized pivot_ne updates complete
  have hb := fermionicKernel_gecp_error_le_localPanelBlock s q residual rows columns
    initial realized pivot_ne updates complete (rows N) (columns N)
    (complete N).1 (complete N).2.1
  exact ((complete k).2.2 x y hx hy).trans ((ha hcount).trans hb)

end GECPKernelStructure.Fermionic
