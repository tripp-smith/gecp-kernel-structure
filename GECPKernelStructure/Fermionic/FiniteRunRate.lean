import GECPKernelStructure.GECP.RunCrossApproximation
import GECPKernelStructure.Fermionic.RealCutoffRate

namespace GECPKernelStructure.GECP.Run

private theorem length_rows_cast {K L : Kernel ℝ ℝ ℝ} (h : K = L) (run : Run K) :
    (h ▸ run).rows.length = run.rows.length := by
  cases h
  rfl

private theorem completeOn_cast {K L : Kernel ℝ ℝ ℝ} (h : K = L) (run : Run K)
    (D C : ℝ → Prop) : (h ▸ run).CompleteOn D C ↔ run.CompleteOn D C := by
  cases h
  rfl

/-- Frequency reflection preserves the actual finite successful count. -/
theorem length_rows_negateColumns {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    run.negateColumns.rows.length = run.rows.length := by
  induction run with
  | nil => rfl
  | @step K p q hp tail ih =>
    have hp' : (fun x y ↦ K x (-y)) p (-q) ≠ 0 := by simpa only [neg_neg] using hp
    have hh : (fun x y ↦ residualUpdate K p q hp x (-y)) =
        residualUpdate (fun x y ↦ K x (-y)) p (-q) hp' := by
      funext x y
      simp only [residualUpdate, neg_neg]
    change (hh ▸ tail.negateColumns).rows.length + 1 = tail.rows.length + 1
    rw [length_rows_cast, ih]

/-- The same stored pivots remain complete on the reflected physical rectangle. -/
theorem completeOn_negateColumns {K : Kernel ℝ ℝ ℝ} (run : Run K) (Λ : ℝ)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ)) :
    run.negateColumns.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) := by
  induction run with
  | nil => trivial
  | @step K p q hp tail ih =>
    have hp' : (fun x y ↦ K x (-y)) p (-q) ≠ 0 := by simpa only [neg_neg] using hp
    have hh : (fun x y ↦ residualUpdate K p q hp x (-y)) =
        residualUpdate (fun x y ↦ K x (-y)) p (-q) hp' := by
      funext x y
      simp only [residualUpdate, neg_neg]
    change CompletePivotOn _ _ (fun x y ↦ K x (-y)) p (-q) ∧ _
    refine ⟨?_, ?_⟩
    · obtain ⟨hrow, hcol, hmax⟩ := complete.1
      refine ⟨hrow, ⟨by linarith, by linarith⟩, ?_⟩
      intro t omega ht hw
      have hn : -Λ ≤ -omega ∧ -omega ≤ Λ := by constructor <;> linarith
      simpa only [neg_neg] using hmax t (-omega) ht hn
    · change (hh ▸ tail.negateColumns).CompleteOn _ _
      exact (completeOn_cast hh _ _ _).2 (ih complete.2)

end GECPKernelStructure.GECP.Run

namespace GECPKernelStructure.Fermionic

open GECP Matrix

/-- Finite local-panel count on the positive-minor presentation. No pivot
after the stored successful run is required. -/
theorem frequencyReflectedFermionicKernel_finiteRun_error_le_localPanelBlock
    (Λ : ℝ) (s q : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ))
    (run : Run frequencyReflectedFermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (hcount : 32 * (s + 1) * (33 + 4 * q) ≤ run.rows.length)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤ (1 / 2 : ℝ) ^ q := by
  classical
  let n := 16 * (33 + 4 * q)
  let N := 32 * (s + 1) * (33 + 4 * q)
  let epsilon : ℝ := (1 / 2 : ℝ) ^ q
  let D : ℝ → Prop := fun t ↦ 0 ≤ t ∧ t ≤ 1
  let C : ℝ → Prop := fun omega ↦ -Λ ≤ omega ∧ omega ≤ Λ
  by_contra hnot
  have hlarge : epsilon < |run.finalResidual x y| := lt_of_not_ge hnot
  have hall (i : ℕ) (hi : i < run.rows.length) :
      epsilon < |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)| :=
    hlarge.trans_le (run.finalResidual_abs_le_pivot
      frequencyReflectedFermionicKernel_strictSignRegular complete i hi x y hx hy)
  have hneg (i : Fin N) :
      -((2 ^ s : ℕ) : ℝ) ≤ -run.pivotColumnAt i.val ∧
        -run.pivotColumnAt i.val ≤ (2 ^ s : ℕ) := by
    have hc := (run.completePivotAt complete i.val (i.isLt.trans_le hcount)).2.1
    obtain ⟨hl, hu⟩ := hc
    constructor <;> linarith
  let label : Fin N → Bool × Fin (s + 1) :=
    fun i ↦ Classical.choose (fermionicPanelDomain_cover s (hneg i))
  have label_valid (i : Fin N) :
      fermionicPanelDomain (label i).1 (label i).2.val (-run.pivotColumnAt i.val) :=
    Classical.choose_spec (fermionicPanelDomain_cover s (hneg i))
  have hsize : Fintype.card (Bool × Fin (s + 1)) * n = Fintype.card (Fin N) := by
    rw [fermionicPanel_count, Fintype.card_fin]
    dsimp [n, N]
    ring
  obtain ⟨panel, hfiber⟩ := Fintype.exists_le_card_fiber_of_mul_le_card label hsize.le
  let fiber : Finset (Fin N) := Finset.univ.filter (fun i ↦ label i = panel)
  let S : Finset ℕ := fiber.map Fin.valEmbedding
  have hS : S ⊆ Finset.range run.rows.length := by
    intro i hi
    obtain ⟨j, _, rfl⟩ := Finset.mem_map.mp hi
    exact Finset.mem_range.mpr (j.isLt.trans_le hcount)
  have hScard : n ≤ S.card := by simpa [S, fiber] using hfiber
  let B : ℝ → Prop := fun omega ↦ fermionicPanelDomain panel.1 panel.2.val (-omega)
  have hSB : ∀ i ∈ S, B (run.pivotColumnAt i) := by
    intro i hi
    obtain ⟨j, hj, rfl⟩ := Finset.mem_map.mp hi
    have heq : label j = panel := (Finset.mem_filter.mp hj).2
    have hh := label_valid j
    change fermionicPanelDomain panel.1 panel.2.val (-run.pivotColumnAt j.val)
    simpa only [heq] using hh
  have budget : ∀ rs cs : Fin n → ℝ, (∀ i, D (rs i)) → (∀ j, B (cs j)) →
      |(minorMatrix frequencyReflectedFermionicKernel rs cs).det| ≤ epsilon ^ n := by
    intro rs cs hD hB
    change |(sampleMatrix fermionicKernel rs (fun j ↦ -cs j)).det| ≤
      ((1 / 2 : ℝ) ^ q) ^ (16 * (33 + 4 * q))
    exact fermionicKernel_panel_minor_det_le_dyadicBudget q panel.2.val panel.1 rs
      (fun j ↦ -cs j) (fun i ↦ (hD i).1) (fun i ↦ (hD i).2) hB
  have hsmall := run.panel_pivot_count_lt_of_initialBudget
    frequencyReflectedFermionicKernel_strictSignRegular D B
    (fun i hi ↦ (run.completePivotAt complete i hi).1)
    (fun i hi z hz ↦ (run.completePivotAt complete i hi).2.2 z _ hz
      (run.completePivotAt complete i hi).2.1)
    n (by dsimp [n]; positivity) epsilon (by dsimp [epsilon]; positivity)
    budget S hS hSB (fun i hi ↦ hall i (Finset.mem_range.mp (hS hi)))
  exact (Nat.not_lt_of_ge hScard) hsmall

/-- The actual fermionic finite run satisfies the same explicit panel bound. -/
theorem fermionicKernel_finiteRun_error_le_localPanelBlock
    (Λ : ℝ) (s q : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ)) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (hcount : 32 * (s + 1) * (33 + 4 * q) ≤ run.rows.length)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤ (1 / 2 : ℝ) ^ q := by
  have hc := run.completeOn_negateColumns Λ complete
  have hn : 32 * (s + 1) * (33 + 4 * q) ≤ run.negateColumns.rows.length := by
    rwa [run.length_rows_negateColumns]
  have hy' : -Λ ≤ -y ∧ -y ≤ Λ := by constructor <;> linarith
  have hh := frequencyReflectedFermionicKernel_finiteRun_error_le_localPanelBlock Λ s q hΛ
    run.negateColumns hc hn x (-y) hx hy'
  have heq := congrFun (congrFun run.finalResidual_negateColumns x) (-y)
  have heq' : run.negateColumns.finalResidual x (-y) = run.finalResidual x y := by
    simpa only [neg_neg] using heq
  calc
    |run.finalResidual x y| = |run.negateColumns.finalResidual x (-y)| := congrArg abs heq'.symm
    _ ≤ _ := hh

theorem fermionicKernel_finiteRun_error_le_coverExponential
    (Λ : ℝ) (s : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ)) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤
      1024 * Real.exp (-(Real.log 2) * run.rows.length / (128 * (s + 1))) := by
  have hone := run.finalResidual_abs_le_initialBound fermionicKernel_strictSignRegular
    complete 1 (fun t omega ht _ ↦ by
      rw [abs_of_pos (fermionicKernel_pos t omega)]
      exact fermionicKernel_le_one ht.1 ht.2) x y hx hy
  exact localPanelBlocks_le_exponential s run.rows.length _ hone (fun q hq ↦
    fermionicKernel_finiteRun_error_le_localPanelBlock Λ s q hΛ run complete hq x y hx hy)

/-- Explicit G1-shaped bound for any finite exact complete run. There is no
assumption about a successful or complete pivot after its final residual. -/
theorem fermionicKernel_finiteRun_error_le_exp
    (Λ : ℝ) (hΛ : 1 ≤ Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ)) := by
  obtain ⟨s, hs, hlog⟩ := exists_logarithmic_cutoff_cover hΛ
  have hh := fermionicKernel_finiteRun_error_le_coverExponential Λ s hs run complete x y hx hy
  refine hh.trans (mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by norm_num))
  exact localPanel_exponent_le_logarithmic s run.rows.length _ (Real.log_pos (by linarith)) hlog

theorem fermionicKernel_finiteRun_error_le_allCutoff_exp
    (Λ : ℝ) (hΛ : 0 < Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (2 + Λ)) := by
  obtain ⟨s, hs, hlog⟩ := exists_allCutoff_logarithmic_cover hΛ
  have hh := fermionicKernel_finiteRun_error_le_coverExponential Λ s hs run complete x y hx hy
  refine hh.trans (mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) (by norm_num))
  exact localPanel_exponent_le_logarithmic s run.rows.length _ (Real.log_pos (by linarith)) hlog

/-- Exact termination is a zero-valued complete maximum, not a successful
zero-pivot update. It permits stopping before the sufficient rank budget. -/
theorem fermionicKernel_finiteRun_error_le_of_stoppedOrBudget
    (Λ : ℝ) (s q : ℕ) (hΛ : Λ ≤ (2 ^ s : ℕ)) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (hstop : 32 * (s + 1) * (33 + 4 * q) ≤ run.rows.length ∨
      ∃ p q : ℝ, CompletePivotOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
        (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ) run.finalResidual p q ∧
        run.finalResidual p q = 0)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤ (1 / 2 : ℝ) ^ q := by
  rcases hstop with hcount | ⟨p, q, hc, hz⟩
  · exact fermionicKernel_finiteRun_error_le_localPanelBlock Λ s q hΛ run complete hcount x y hx hy
  · rw [residual_eq_zero_of_complete_zero hc hz x y hx hy, abs_zero]
    positivity

theorem fermionicKernel_finiteRun_error_le_realTolerance
    (Λ epsilon : ℝ) (hΛ : 0 < Λ) (hepsilon : 0 < epsilon) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (hcount : ⌈256 / (Real.log 2) ^ 2 * Real.log (2 + Λ) * Real.log (1024 / epsilon)⌉₊ ≤
      run.rows.length)
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |run.finalResidual x y| ≤ epsilon := by
  exact (fermionicKernel_finiteRun_error_le_allCutoff_exp Λ hΛ run complete x y hx hy).trans
    (logarithmic_envelope_le_tolerance _ epsilon (Real.log_pos (by linarith)) hepsilon
      run.rows.length ((Nat.le_ceil _).trans (Nat.cast_le.mpr hcount)))

/-- The rate is a bound on the original selected-core cross approximation. -/
theorem fermionicKernel_finiteRun_cross_error_le_exp
    (Λ : ℝ) (hΛ : 1 ≤ Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (x y : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) (hy : -Λ ≤ y ∧ y ≤ Λ) :
    |fermionicKernel x y - run.crossApproximation x y| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ)) := by
  rw [run.kernel_sub_crossApproximation_eq_finalResidual]
  exact fermionicKernel_finiteRun_error_le_exp Λ hΛ run complete x y hx hy

/-- The actual difference of Green integrals inherits the finite-run rate. -/
theorem fermionicKernel_finiteRun_greenError_le_exp
    (Λ : ℝ) (hΛ : 1 ≤ Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (μ : MeasureTheory.Measure ℝ) (rho : ℝ → ℝ) (hrho : MeasureTheory.Integrable rho μ)
    (hsupport : ∀ᵐ y ∂μ, -Λ ≤ y ∧ y ≤ Λ) (x : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) :
    |(∫ y, rho y * fermionicKernel x y ∂μ) -
        (∫ y, rho y * run.crossApproximation x y ∂μ)| ≤
      (1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ))) *
        ∫ y, |rho y| ∂μ := by
  apply run.greenError_le_uniformCrossError fermionicKernel_continuous μ rho hrho
  · exact Filter.Eventually.of_forall (fun y ↦ by
      rw [abs_of_pos (fermionicKernel_pos x y)]
      exact fermionicKernel_le_one hx.1 hx.2)
  · exact hsupport.mono (fun y hy ↦ fermionicKernel_finiteRun_error_le_exp
      Λ hΛ run complete x y hx hy)

/-- The physical interval-integral form needs no separate support premise. -/
theorem fermionicKernel_finiteRun_greenError_le_exp_interval
    (Λ : ℝ) (hΛ : 1 ≤ Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ ↦ 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ ↦ -Λ ≤ omega ∧ omega ≤ Λ))
    (μ : MeasureTheory.Measure ℝ) (rho : ℝ → ℝ)
    (hrho : MeasureTheory.Integrable rho (μ.restrict (Set.Icc (-Λ) Λ)))
    (x : ℝ) (hx : 0 ≤ x ∧ x ≤ 1) :
    |(∫ y, rho y * fermionicKernel x y ∂μ.restrict (Set.Icc (-Λ) Λ)) -
        (∫ y, rho y * run.crossApproximation x y ∂μ.restrict (Set.Icc (-Λ) Λ))| ≤
      (1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ))) *
        ∫ y, |rho y| ∂μ.restrict (Set.Icc (-Λ) Λ) := by
  exact fermionicKernel_finiteRun_greenError_le_exp Λ hΛ run complete
    (μ.restrict (Set.Icc (-Λ) Λ)) rho hrho (MeasureTheory.ae_restrict_mem measurableSet_Icc) x hx

end GECPKernelStructure.Fermionic
