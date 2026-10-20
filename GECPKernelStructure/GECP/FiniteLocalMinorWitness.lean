import GECPKernelStructure.GECP.LocalMinorWitness

namespace GECPKernelStructure.GECP.Run

variable {K : Kernel ℝ ℝ ℝ}

/-- Extract at most the first `n` successful steps. No continuation is added. -/
noncomputable def take : {K : Kernel ℝ ℝ ℝ} → Run K → ℕ → Run K
  | K, _, 0 => .nil K
  | _, .nil K, _ + 1 => .nil K
  | _, .step p q hp tail, n + 1 => .step p q hp (take tail n)

@[simp] theorem take_zero (run : Run K) : run.take 0 = .nil K := by
  cases run <;> rfl

/-- The residual at a finite prefix; after the last step it stays unchanged. -/
noncomputable def residualAt (run : Run K) (n : ℕ) : Kernel ℝ ℝ ℝ :=
  (run.take n).finalResidual

@[simp] theorem residualAt_zero (run : Run K) : run.residualAt 0 = K := by
  simp only [residualAt, take_zero, finalResidual]

/-- Pivot coordinates are queried only below the stored successful length. -/
noncomputable def pivotRowAt : {K : Kernel ℝ ℝ ℝ} → Run K → ℕ → ℝ
  | _, .nil _, _ => 0
  | _, .step p _ _ _, 0 => p
  | _, .step _ _ _ tail, n + 1 => pivotRowAt tail n

noncomputable def pivotColumnAt : {K : Kernel ℝ ℝ ℝ} → Run K → ℕ → ℝ
  | _, .nil _, _ => 0
  | _, .step _ q _ _, 0 => q
  | _, .step _ _ _ tail, n + 1 => pivotColumnAt tail n

theorem take_eq_self (run : Run K) (n : ℕ) (h : run.rows.length ≤ n) :
    run.take n = run := by
  induction run generalizing n with
  | nil => cases n <;> rfl
  | step p q hp tail ih =>
    cases n with
    | zero => simp [rows] at h
    | succ n =>
      change Run.step p q hp (tail.take n) = Run.step p q hp tail
      congr 1
      exact ih n (by change tail.rows.length + 1 ≤ n + 1 at h; omega)

@[simp] theorem residualAt_length (run : Run K) :
    run.residualAt run.rows.length = run.finalResidual := by
  unfold residualAt
  rw [take_eq_self run _ le_rfl]

theorem length_rows_take (run : Run K) (n : ℕ) :
    (run.take n).rows.length = min n run.rows.length := by
  induction run generalizing n with
  | nil => cases n <;> simp [take, rows]
  | step p q hp tail ih =>
    cases n with
    | zero => simp [take, rows]
    | succ n =>
      change (tail.take n).rows.length + 1 = min (n + 1) (tail.rows.length + 1)
      have hh := ih n
      omega

theorem completeOn_take (run : Run K) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (n : ℕ) : (run.take n).CompleteOn D C := by
  induction run generalizing n with
  | nil => cases n <;> trivial
  | step p q hp tail ih =>
    cases n with
    | zero => trivial
    | succ n => exact ⟨complete.1, ih complete.2 n⟩

theorem pivotAt_ne (run : Run K) (i : ℕ) (hi : i < run.rows.length) :
    run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i) ≠ 0 := by
  induction run generalizing i with
  | nil => simp [rows] at hi
  | step p q hp tail ih =>
    cases i with
    | zero => exact hp
    | succ i =>
      exact ih i (by change i + 1 < tail.rows.length + 1 at hi; omega)

theorem residualAt_succ (run : Run K) (i : ℕ) (hi : i < run.rows.length) :
    run.residualAt (i + 1) = residualUpdate (run.residualAt i)
      (run.pivotRowAt i) (run.pivotColumnAt i) (run.pivotAt_ne i hi) := by
  induction run generalizing i with
  | nil => simp [rows] at hi
  | step p q hp tail ih =>
    cases i with
    | zero => simp only [residualAt, take, take_zero, finalResidual, pivotRowAt, pivotColumnAt]
    | succ i =>
      exact ih i (by change i + 1 < tail.rows.length + 1 at hi; omega)

theorem completePivotAt (run : Run K) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (i : ℕ) (hi : i < run.rows.length) :
    CompletePivotOn D C (run.residualAt i) (run.pivotRowAt i) (run.pivotColumnAt i) := by
  induction run generalizing i with
  | nil => simp [rows] at hi
  | step p q hp tail ih =>
    cases i with
    | zero => exact complete.1
    | succ i =>
      exact ih complete.2 i (by change i + 1 < tail.rows.length + 1 at hi; omega)

/-- The next residual is controlled by this stored pivot, without a next pivot. -/
theorem abs_residualAt_succ_le_pivot (run : Run K) {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (i : ℕ) (hi : i < run.rows.length)
    (x y : ℝ) (hx : D x) (hy : C y) :
    |run.residualAt (i + 1) x y| ≤
      |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)| := by
  have hc := run.completePivotAt complete i hi
  have coherent := Fermionic.strictSignRegular_pivotCrossProductSignCoherent regular
    (run.take i) (run.pivotRowAt i) (run.pivotColumnAt i)
  rw [run.residualAt_succ i hi]
  exact Fermionic.residualUpdate_le_of_signCoherentOn coherent
    (hc.2.2 x y hx hy) rfl (hc.2.2 x _ hx hc.2.1)
    (hc.2.2 _ y hc.1 hy) (run.pivotAt_ne i hi)

theorem abs_residualAt_add_succ_le_pivot (run : Run K) {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (i d : ℕ) (hi : i + d < run.rows.length)
    (x y : ℝ) (hx : D x) (hy : C y) :
    |run.residualAt (i + d + 1) x y| ≤
      |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)| := by
  induction d generalizing i with
  | zero => exact run.abs_residualAt_succ_le_pivot regular complete i hi x y hx hy
  | succ d ih =>
    have hn : i + 1 < run.rows.length := by omega
    have hc := run.completePivotAt complete (i + 1) hn
    have ht := ih (i + 1) (by omega)
    have heq : i + 1 + d + 1 = i + (d + 1) + 1 := by omega
    rw [heq] at ht
    exact ht.trans (run.abs_residualAt_succ_le_pivot regular complete i (by omega)
      _ _ hc.1 hc.2.1)

/-- Every stored pivot dominates the final residual on the physical domain. -/
theorem finalResidual_abs_le_pivot (run : Run K) {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (i : ℕ) (hi : i < run.rows.length)
    (x y : ℝ) (hx : D x) (hy : C y) :
    |run.finalResidual x y| ≤
      |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)| := by
  have hh := run.abs_residualAt_add_succ_le_pivot regular complete i
    (run.rows.length - i - 1) (by omega) x y hx hy
  have heq : i + (run.rows.length - i - 1) + 1 = run.rows.length := by omega
  simpa only [heq, residualAt_length] using hh

theorem finalResidual_abs_le_initialBound (run : Run K) {signature : ℕ → ℝ}
    (regular : StrictSignRegular K signature) {D C : ℝ → Prop}
    (complete : run.CompleteOn D C) (B : ℝ)
    (hB : ∀ x y, D x → C y → |K x y| ≤ B)
    (x y : ℝ) (hx : D x) (hy : C y) : |run.finalResidual x y| ≤ B := by
  by_cases hn : run.rows.length = 0
  · rw [← run.residualAt_length, hn]
    simpa only [residualAt_zero] using hB x y hx hy
  · have hi : 0 < run.rows.length := Nat.pos_of_ne_zero hn
    have hc := run.completePivotAt complete 0 hi
    exact (run.finalResidual_abs_le_pivot regular complete 0 hi x y hx hy).trans
      (by simpa only [residualAt_zero] using hB _ _ hc.1 hc.2.1)

theorem markedPivotProduct_pos (run : Run K) (mark : ℕ → Bool) (start length : ℕ)
    (hbound : start + length ≤ run.rows.length) :
    0 < GECP.markedPivotProduct
      (fun i ↦ run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)) mark start length := by
  apply Finset.prod_pos
  intro i hi
  have hii : start + i < run.rows.length := by
    have hh := Finset.mem_range.mp hi
    omega
  split_ifs
  · exact abs_pos.mpr (run.pivotAt_ne _ hii)
  · norm_num

/-- Backwards witnesses require hypotheses only inside this finite stored run. -/
theorem exists_marked_minor_witness (run : Run K)
    (regular : StrictSignRegular K (fun _ ↦ 1))
    (D B : ℝ → Prop)
    (pivot_mem : ∀ i, i < run.rows.length → D (run.pivotRowAt i))
    (columnMax : ∀ i, i < run.rows.length → ∀ x, D x →
      |run.residualAt i x (run.pivotColumnAt i)| ≤
        |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)|)
    (mark : ℕ → Bool)
    (marked_mem : ∀ i, i < run.rows.length → mark i = true → B (run.pivotColumnAt i))
    (start length : ℕ) (hbound : start + length ≤ run.rows.length) :
    ∃ rows columns : Fin (markedCount mark start length) → ℝ,
      (∀ i, D (rows i)) ∧ (∀ j, B (columns j)) ∧
      GECP.markedPivotProduct
        (fun i ↦ run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)) mark start length ≤
        |(minorMatrix (run.residualAt start) rows columns).det| := by
  induction length generalizing start with
  | zero =>
    refine ⟨Fin.elim0, Fin.elim0, ?_, ?_, ?_⟩
    · intro i; exact Fin.elim0 i
    · intro i; exact Fin.elim0 i
    · change (1 : ℝ) ≤ |(minorMatrix (run.residualAt start)
        (Fin.elim0 : Fin 0 → ℝ) Fin.elim0).det|
      rw [Matrix.det_fin_zero]
      norm_num
  | succ length ih =>
    have hi : start < run.rows.length := by omega
    obtain ⟨rows, columns, hD, hB, hbound'⟩ := ih (start + 1) (by omega)
    have hpos := run.markedPivotProduct_pos mark (start + 1) length (by omega)
    have hminor : (minorMatrix (run.residualAt (start + 1)) rows columns).det ≠ 0 := by
      intro hz
      rw [hz, abs_zero] at hbound'
      linarith
    rw [markedCount_succ, GECP.markedPivotProduct_succ]
    cases hm : mark start with
    | false =>
      simp only [Bool.false_eq_true, if_false, mul_one]
      have hminor' := hminor
      rw [run.residualAt_succ start hi] at hminor'
      obtain ⟨oldRows, oldColumns, hoD, hoB, hb⟩ := exists_minor_witness_before_step
        regular (run.take start) D B (run.pivotRowAt start) (run.pivotColumnAt start)
          (run.pivotAt_ne start hi) (pivot_mem start hi) (columnMax start hi)
          rows columns hD hB hminor'
      refine ⟨oldRows, oldColumns, hoD, hoB, hbound'.trans ?_⟩
      change |(minorMatrix (residualUpdate (run.residualAt start) _ _ _) rows columns).det| ≤
        |(minorMatrix (run.residualAt start) oldRows oldColumns).det| at hb
      rwa [← run.residualAt_succ start hi] at hb
    | true =>
      refine ⟨Fin.cons (run.pivotRowAt start) rows, Fin.cons (run.pivotColumnAt start) columns, ?_, ?_, ?_⟩
      · exact Fin.cases (pivot_mem start hi) hD
      · exact Fin.cases (marked_mem start hi hm) hB
      · rw [minor_cons_det_eq_pivot_mul_update _ _ (run.pivotAt_ne start hi), abs_mul,
          ← run.residualAt_succ start hi]
        exact (mul_le_mul_of_nonneg_right hbound' (abs_nonneg _)).trans_eq (mul_comm _ _)

/-- A fixed-order initial local budget bounds marked pivots of a finite run. -/
theorem panel_pivot_count_lt_of_initialBudget (run : Run K)
    (regular : StrictSignRegular K (fun _ ↦ 1))
    (D B : ℝ → Prop)
    (pivot_mem : ∀ i, i < run.rows.length → D (run.pivotRowAt i))
    (columnMax : ∀ i, i < run.rows.length → ∀ x, D x →
      |run.residualAt i x (run.pivotColumnAt i)| ≤
        |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)|)
    (n : ℕ) (hn : 0 < n) (epsilon : ℝ) (hepsilon : 0 < epsilon)
    (budget : ∀ rows columns : Fin n → ℝ, (∀ i, D (rows i)) →
      (∀ j, B (columns j)) → |(minorMatrix K rows columns).det| ≤ epsilon ^ n)
    (S : Finset ℕ) (hS : S ⊆ Finset.range run.rows.length)
    (hSB : ∀ i ∈ S, B (run.pivotColumnAt i))
    (hSlarge : ∀ i ∈ S, epsilon <
      |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)|) : S.card < n := by
  classical
  by_contra hnot
  obtain ⟨T, hTS, hTcard⟩ := Finset.exists_subset_card_eq (Nat.le_of_not_gt hnot)
  have hT : T ⊆ Finset.range run.rows.length := hTS.trans hS
  have hm : ∀ i, i < run.rows.length → decide (i ∈ T) = true → B (run.pivotColumnAt i) := by
    intro i _ hi
    exact hSB i (hTS (of_decide_eq_true hi))
  have hh := run.exists_marked_minor_witness regular D B pivot_mem columnMax
    (fun i ↦ decide (i ∈ T)) hm 0 run.rows.length (by omega)
  rw [markedCount_of_finset T _ hT, markedPivotProduct_of_finset _ T _ hT, hTcard] at hh
  obtain ⟨rows, columns, hD, hB, hprod⟩ := hh
  rw [residualAt_zero] at hprod
  have hneT : T.Nonempty := Finset.card_pos.mp (by omega)
  have hstrict : epsilon ^ n < ∏ i ∈ T,
      |run.residualAt i (run.pivotRowAt i) (run.pivotColumnAt i)| := by
    rw [← hTcard, ← Finset.prod_const]
    exact Finset.prod_lt_prod_of_nonempty (fun _ _ ↦ hepsilon)
      (fun i hi ↦ hSlarge i (hTS hi)) hneT
  exact (not_lt_of_ge (hprod.trans (budget rows columns hD hB))) hstrict

end GECPKernelStructure.GECP.Run
