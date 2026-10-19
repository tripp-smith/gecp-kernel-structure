import GECPKernelStructure.GECP.ResidualTotalPositivity
import GECPKernelStructure.GECP.CompanionClosedForm

namespace GECPKernelStructure.GECP

open Matrix

namespace Run

/-- Any nonzero residual entry certifies freshness of its row. -/
theorem freshRow_of_finalResidual_ne_zero {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {x y : ℝ} (hne : run.finalResidual x y ≠ 0) : run.FreshRow x := by
  have avoid (i : run.SelectedIndex) : run.selectedRow i ≠ x := by
    intro heq
    apply hne
    rw [← heq]
    exact run.gecp_interpolates_selected_rows (run.selectedRow_mem_rows i) y
  have hh : Function.Injective (run.augmentedRow (Fin 1) (fun _ ↦ x) ∘
      run.selectedBorderEquiv (Fin 1)) := by
    rw [run.augmentedRow_comp_selectedBorderEquiv]
    exact run.selectedRow_injective.sumElim
      (fun _ _ _ ↦ Subsingleton.elim _ _) (fun i _ ↦ avoid i)
  intro i j heq
  apply (run.selectedBorderEquiv (Fin 1)).symm.injective
  apply hh
  simpa only [Function.comp_apply, Equiv.apply_symm_apply] using heq

/-- Any nonzero residual entry certifies freshness of its column. -/
theorem freshColumn_of_finalResidual_ne_zero {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {x y : ℝ} (hne : run.finalResidual x y ≠ 0) : run.FreshColumn y := by
  have avoid (i : run.SelectedIndex) : run.selectedColumn i ≠ y := by
    intro heq
    apply hne
    rw [← heq]
    exact run.gecp_interpolates_selected_cols (run.selectedColumn_mem_columns i) x
  have hh : Function.Injective (run.augmentedColumn (Fin 1) (fun _ ↦ y) ∘
      run.selectedBorderEquiv (Fin 1)) := by
    rw [run.augmentedColumn_comp_selectedBorderEquiv]
    exact run.selectedColumn_injective.sumElim
      (fun _ _ _ ↦ Subsingleton.elim _ _) (fun i _ ↦ avoid i)
  intro i j heq
  apply (run.selectedBorderEquiv (Fin 1)).symm.injective
  apply hh
  simpa only [Function.comp_apply, Equiv.apply_symm_apply] using heq

theorem freshRows_of_minor_ne_zero {K : Kernel ℝ ℝ ℝ} (run : Run K) {n : ℕ}
    (rows columns : Fin n → ℝ) (hne : (minorMatrix run.finalResidual rows columns).det ≠ 0) :
    ∀ i, run.FreshRow (rows i) := by
  intro i
  have hh : ∃ j, run.finalResidual (rows i) (columns j) ≠ 0 := by
    by_contra h
    push Not at h
    exact hne (Matrix.det_eq_zero_of_row_eq_zero i h)
  obtain ⟨j, hj⟩ := hh
  exact run.freshRow_of_finalResidual_ne_zero hj

theorem freshColumns_of_minor_ne_zero {K : Kernel ℝ ℝ ℝ} (run : Run K) {n : ℕ}
    (rows columns : Fin n → ℝ) (hne : (minorMatrix run.finalResidual rows columns).det ≠ 0) :
    ∀ j, run.FreshColumn (columns j) := by
  intro j
  have hh : ∃ i, run.finalResidual (rows i) (columns j) ≠ 0 := by
    by_contra h
    push Not at h
    exact hne (Matrix.det_eq_zero_of_column_eq_zero j h)
  obtain ⟨i, hi⟩ := hh
  exact run.freshColumn_of_finalResidual_ne_zero hi

end Run

/-- One pivot factors the minor on its row/column plus any finite border. -/
theorem minor_cons_det_eq_pivot_mul_update {R : Kernel ℝ ℝ ℝ} (p q : ℝ)
    (hp : R p q ≠ 0) {n : ℕ} (rows columns : Fin n → ℝ) :
    (minorMatrix R (Fin.cons p rows) (Fin.cons q columns)).det =
      R p q * (minorMatrix (residualUpdate R p q hp) rows columns).det := by
  let A := minorMatrix R (Fin.cons p rows) (Fin.cons q columns)
  have hA : A 0 0 ≠ 0 := hp
  have hs : Matrix.northwestSchur A hA =
      minorMatrix (residualUpdate R p q hp) rows columns := by
    ext i j
    simp [Matrix.northwestSchur, residualUpdate, A, minorMatrix]
  have hh := Matrix.det_eq_northwest_pivot_mul_schur A hA
  rw [hs] at hh
  exact hh

/-- An unmarked column-maximal step can be crossed backwards without losing
a nonzero local minor. The replacement row tuple may change, while every
replacement column still belongs to the original local column predicate. -/
theorem exists_minor_witness_before_step {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1)) (run : Run K)
    (D B : ℝ → Prop) (p q : ℝ) (hp : run.finalResidual p q ≠ 0)
    (hDp : D p) (columnMax : ∀ x, D x → |run.finalResidual x q| ≤ |run.finalResidual p q|)
    {n : ℕ} (rows columns : Fin n → ℝ)
    (hD : ∀ i, D (rows i)) (hB : ∀ j, B (columns j))
    (hne : (minorMatrix (residualUpdate run.finalResidual p q hp) rows columns).det ≠ 0) :
    ∃ oldRows oldColumns : Fin n → ℝ,
      (∀ i, D (oldRows i)) ∧ (∀ j, B (oldColumns j)) ∧
      |(minorMatrix (residualUpdate run.finalResidual p q hp) rows columns).det| ≤
        |(minorMatrix run.finalResidual oldRows oldColumns).det| := by
  classical
  let ar : Fin (n + 1) → ℝ := Fin.cons p rows
  let ac : Fin (n + 1) → ℝ := Fin.cons q columns
  have haug : (minorMatrix run.finalResidual ar ac).det ≠ 0 := by
    rw [minor_cons_det_eq_pivot_mul_update p q hp]
    exact mul_ne_zero hp hne
  have hrinj : Function.Injective ar := by
    by_contra h
    exact haug (det_kernel_eq_zero_of_not_injective_left run.finalResidual ar ac h)
  have hcinj : Function.Injective ac := by
    by_contra h
    exact haug (det_kernel_eq_zero_of_not_injective_right run.finalResidual ar ac h)
  let rf : Fin (n + 1) → {x // run.FreshRow x} :=
    fun i ↦ ⟨ar i, run.freshRows_of_minor_ne_zero ar ac haug i⟩
  let cf : Fin (n + 1) → {y // run.FreshColumn y} :=
    fun j ↦ ⟨ac j, run.freshColumns_of_minor_ne_zero ar ac haug j⟩
  have hrf : Function.Injective rf := fun i j h ↦ hrinj (congrArg Subtype.val h)
  have hcf : Function.Injective cf := fun i j h ↦ hcinj (congrArg Subtype.val h)
  let rp := Tuple.sort rf
  let cp := Tuple.sort cf
  let r := rf ∘ rp
  let c := cf ∘ cp
  have hr : StrictMono r :=
    (Tuple.monotone_sort rf).strictMono_of_injective (hrf.comp rp.injective)
  have hc : StrictMono c :=
    (Tuple.monotone_sort cf).strictMono_of_injective (hcf.comp cp.injective)
  let A := minorMatrix run.gaugedFinalResidual r c
  have hA : Matrix.StrictlyTotallyPositive A :=
    strictlyTotallyPositive_minorMatrix_of_strictSignRegular_one
      (run.gaugedFinalResidual_strictSignRegular regular) r c hr hc
  let k := cp.symm 0
  have hck : (c k : ℝ) = q := by simp [c, cf, k, ac]
  have hDom (i : Fin (n + 1)) : D (rf i : ℝ) := by
    exact Fin.cases hDp (fun j ↦ hD j) i
  have hEntry (i j : Fin (n + 1)) : A i j = |run.finalResidual (r i) (c j)| := by
    calc
      A i j = |A i j| := (abs_of_pos (hA.entry_pos i j)).symm
      _ = _ := run.abs_gaugedFinalResidual (r i) (c j)
  have hdiag : A k k ≤ |run.finalResidual p q| := by
    rw [hEntry, hck]
    exact columnMax _ (hDom (rp k))
  let oldRows : Fin n → ℝ := fun i ↦ (r (k.succAbove i) : ℝ)
  let oldColumns : Fin n → ℝ := fun j ↦ (c (k.succAbove j) : ℝ)
  have oldD : ∀ i, D (oldRows i) := fun i ↦ hDom (rp (k.succAbove i))
  have oldB : ∀ j, B (oldColumns j) := by
    intro j
    have hne0 : cp (k.succAbove j) ≠ 0 := by
      intro hz
      have heq : k.succAbove j = k := cp.injective (by simpa [k] using hz)
      exact Fin.succAbove_ne k j heq
    have htail : ∀ l : Fin (n + 1), l ≠ 0 → B (cf l : ℝ) := by
      refine Fin.cases (fun h ↦ (h rfl).elim) (fun l _ ↦ ?_)
      exact hB l
    exact htail _ hne0
  have hminor : 0 < (A.submatrix k.succAbove k.succAbove).det :=
    hA n _ _ (Fin.strictMono_succAbove _) (Fin.strictMono_succAbove _)
  have holdeq : (A.submatrix k.succAbove k.succAbove).det =
      |(minorMatrix run.finalResidual oldRows oldColumns).det| := by
    calc
      _ = |(A.submatrix k.succAbove k.succAbove).det| := (abs_of_pos hminor).symm
      _ = _ := run.abs_gaugedFinalResidual_minor_det (r ∘ k.succAbove) (c ∘ k.succAbove)
  have hAdetpos : 0 < A.det := by simpa using hA (n + 1) id id strictMono_id strictMono_id
  have hAdet : A.det = |run.finalResidual p q| *
      |(minorMatrix (residualUpdate run.finalResidual p q hp) rows columns).det| := by
    calc
      A.det = |A.det| := (abs_of_pos hAdetpos).symm
      _ = |(minorMatrix run.gaugedFinalResidual rf cf).det| :=
        Matrix.abs_det_submatrix_equiv_equiv rp cp (minorMatrix run.gaugedFinalResidual rf cf)
      _ = |(minorMatrix run.finalResidual ar ac).det| :=
        run.abs_gaugedFinalResidual_minor_det rf cf
      _ = _ := by rw [minor_cons_det_eq_pivot_mul_update p q hp, abs_mul]
  have hf := Matrix.det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive A hA k
  have hprod : |run.finalResidual p q| *
      |(minorMatrix (residualUpdate run.finalResidual p q hp) rows columns).det| ≤
      |run.finalResidual p q| * |(minorMatrix run.finalResidual oldRows oldColumns).det| := by
    calc
      _ = A.det := hAdet.symm
      _ ≤ A k k * (A.submatrix k.succAbove k.succAbove).det := hf
      _ ≤ |run.finalResidual p q| * (A.submatrix k.succAbove k.succAbove).det :=
        mul_le_mul_of_nonneg_right hdiag hminor.le
      _ = _ := by rw [holdeq]
  exact ⟨oldRows, oldColumns, oldD, oldB,
    (mul_le_mul_iff_right₀ (abs_pos.mpr hp)).mp hprod⟩

/-- Count marked steps in an interval of a pivot sequence. -/
def markedCount (mark : ℕ → Bool) (start length : ℕ) : ℕ :=
  ∑ i ∈ Finset.range length, if mark (start + i) then 1 else 0

/-- Product of the selected absolute pivots, with unmarked steps contributing one. -/
noncomputable def markedPivotProduct (pivot : ℕ → ℝ) (mark : ℕ → Bool)
    (start length : ℕ) : ℝ :=
  ∏ i ∈ Finset.range length, if mark (start + i) then |pivot (start + i)| else 1

theorem markedCount_succ (mark : ℕ → Bool) (start length : ℕ) :
    markedCount mark start (length + 1) =
      markedCount mark (start + 1) length + if mark start then 1 else 0 := by
  unfold markedCount
  rw [Finset.sum_range_succ']
  simp only [Nat.add_zero]
  congr 1
  apply Finset.sum_congr rfl
  intro i _
  rw [show start + (i + 1) = start + 1 + i by omega]

theorem markedPivotProduct_succ (pivot : ℕ → ℝ) (mark : ℕ → Bool) (start length : ℕ) :
    markedPivotProduct pivot mark start (length + 1) =
      markedPivotProduct pivot mark (start + 1) length *
        if mark start then |pivot start| else 1 := by
  unfold markedPivotProduct
  rw [Finset.prod_range_succ']
  simp only [Nat.add_zero]
  congr 1
  apply Finset.prod_congr rfl
  intro i _
  rw [show start + (i + 1) = start + 1 + i by omega]

theorem markedPivotProduct_pos (pivot : ℕ → ℝ) (mark : ℕ → Bool)
    (hne : ∀ i, pivot i ≠ 0) (start length : ℕ) :
    0 < markedPivotProduct pivot mark start length := by
  apply Finset.prod_pos
  intro i _
  split_ifs
  · exact abs_pos.mpr (hne _)
  · norm_num

/-- Backwards induction constructs an actual minor dominating the marked
pivot product. No residual-envelope premise is assumed. -/
theorem exists_marked_minor_witness {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1))
    (residual : ℕ → Kernel ℝ ℝ ℝ) (pivotRows pivotColumns : ℕ → ℝ)
    (realized : ∀ i, ∃ run : Run K, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (pivotRows i) (pivotColumns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) =
      residualUpdate (residual i) (pivotRows i) (pivotColumns i) (pivot_ne i))
    (D B : ℝ → Prop) (pivot_mem : ∀ i, D (pivotRows i))
    (columnMax : ∀ i x, D x → |residual i x (pivotColumns i)| ≤
      |residual i (pivotRows i) (pivotColumns i)|)
    (mark : ℕ → Bool) (marked_mem : ∀ i, mark i = true → B (pivotColumns i))
    (start length : ℕ) :
    ∃ rows columns : Fin (markedCount mark start length) → ℝ,
      (∀ i, D (rows i)) ∧ (∀ j, B (columns j)) ∧
      markedPivotProduct (fun i ↦ residual i (pivotRows i) (pivotColumns i)) mark start length ≤
        |(minorMatrix (residual start) rows columns).det| := by
  induction length generalizing start with
  | zero =>
    refine ⟨Fin.elim0, Fin.elim0, ?_, ?_, ?_⟩
    · intro i; exact Fin.elim0 i
    · intro i; exact Fin.elim0 i
    · change (1 : ℝ) ≤ |(minorMatrix (residual start)
        (Fin.elim0 : Fin 0 → ℝ) Fin.elim0).det|
      rw [Matrix.det_fin_zero]
      norm_num
  | succ length ih =>
    obtain ⟨rows, columns, hD, hB, hbound⟩ := ih (start + 1)
    have hpos := markedPivotProduct_pos
      (fun i ↦ residual i (pivotRows i) (pivotColumns i)) mark pivot_ne (start + 1) length
    have hminor : (minorMatrix (residual (start + 1)) rows columns).det ≠ 0 := by
      intro hz
      rw [hz, abs_zero] at hbound
      linarith
    rw [markedCount_succ, markedPivotProduct_succ]
    cases hm : mark start with
    | false =>
      simp only [Bool.false_eq_true, if_false, mul_one]
      obtain ⟨run, hr⟩ := realized start
      have hp : run.finalResidual (pivotRows start) (pivotColumns start) ≠ 0 := by
        simpa only [hr] using pivot_ne start
      have hnext : residual (start + 1) =
          residualUpdate run.finalResidual (pivotRows start) (pivotColumns start) hp := by
        rw [updates start]
        funext x y
        simp only [residualUpdate, hr]
      have hminor' : (minorMatrix (residualUpdate run.finalResidual
          (pivotRows start) (pivotColumns start) hp) rows columns).det ≠ 0 := by
        simpa only [← hnext] using hminor
      have hmax : ∀ x, D x → |run.finalResidual x (pivotColumns start)| ≤
          |run.finalResidual (pivotRows start) (pivotColumns start)| := by
        simpa only [hr] using columnMax start
      obtain ⟨oldRows, oldColumns, hoD, hoB, hb⟩ := exists_minor_witness_before_step
        regular run D B (pivotRows start) (pivotColumns start) hp (pivot_mem start)
          hmax rows columns hD hB hminor'
      refine ⟨oldRows, oldColumns, hoD, hoB, hbound.trans ?_⟩
      simpa only [hr, ← hnext] using hb
    | true =>
      refine ⟨Fin.cons (pivotRows start) rows, Fin.cons (pivotColumns start) columns, ?_, ?_, ?_⟩
      · exact Fin.cases (pivot_mem start) hD
      · exact Fin.cases (marked_mem start hm) hB
      · rw [minor_cons_det_eq_pivot_mul_update (pivotRows start) (pivotColumns start)
          (pivot_ne start), abs_mul, ← updates start]
        calc
          _ ≤ |(minorMatrix (residual (start + 1)) rows columns).det| *
              |residual start (pivotRows start) (pivotColumns start)| :=
            mul_le_mul_of_nonneg_right hbound (abs_nonneg _)
          _ = _ := mul_comm _ _

/-- A marked pivot product is dominated by an original-kernel minor on the
same local column predicate and the physical row domain. -/
theorem markedPivotProduct_le_initialMinor {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1))
    (residual : ℕ → Kernel ℝ ℝ ℝ) (pivotRows pivotColumns : ℕ → ℝ)
    (initial : residual 0 = K)
    (realized : ∀ i, ∃ run : Run K, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (pivotRows i) (pivotColumns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) =
      residualUpdate (residual i) (pivotRows i) (pivotColumns i) (pivot_ne i))
    (D B : ℝ → Prop) (pivot_mem : ∀ i, D (pivotRows i))
    (columnMax : ∀ i x, D x → |residual i x (pivotColumns i)| ≤
      |residual i (pivotRows i) (pivotColumns i)|)
    (mark : ℕ → Bool) (marked_mem : ∀ i, mark i = true → B (pivotColumns i)) (length : ℕ) :
    ∃ rows columns : Fin (markedCount mark 0 length) → ℝ,
      (∀ i, D (rows i)) ∧ (∀ j, B (columns j)) ∧
      markedPivotProduct (fun i ↦ residual i (pivotRows i) (pivotColumns i)) mark 0 length ≤
        |(minorMatrix K rows columns).det| := by
  have hh := exists_marked_minor_witness regular residual pivotRows pivotColumns realized
    pivot_ne updates D B pivot_mem columnMax mark marked_mem 0 length
  simpa only [initial] using hh

theorem markedCount_of_finset (S : Finset ℕ) (length : ℕ)
    (hS : S ⊆ Finset.range length) :
    markedCount (fun i ↦ decide (i ∈ S)) 0 length = S.card := by
  unfold markedCount
  simp only [Nat.zero_add, decide_eq_true_eq]
  exact (Finset.card_eq_sum_ite hS).symm

theorem markedPivotProduct_of_finset (pivot : ℕ → ℝ) (S : Finset ℕ) (length : ℕ)
    (hS : S ⊆ Finset.range length) :
    markedPivotProduct pivot (fun i ↦ decide (i ∈ S)) 0 length = ∏ i ∈ S, |pivot i| := by
  unfold markedPivotProduct
  simp only [Nat.zero_add, decide_eq_true_eq]
  rw [Finset.prod_ite]
  simp only [Finset.filter_mem_eq_inter, Finset.prod_const_one, mul_one,
    Finset.inter_eq_right.mpr hS]

/-- A fixed-order initial local budget bounds the number of above-threshold
pivots in that panel. Exactly the budget order is marked in the contradiction. -/
theorem panel_pivot_count_lt_of_initialBudget {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1))
    (residual : ℕ → Kernel ℝ ℝ ℝ) (pivotRows pivotColumns : ℕ → ℝ)
    (initial : residual 0 = K)
    (realized : ∀ i, ∃ run : Run K, run.finalResidual = residual i)
    (pivot_ne : ∀ i, residual i (pivotRows i) (pivotColumns i) ≠ 0)
    (updates : ∀ i, residual (i + 1) =
      residualUpdate (residual i) (pivotRows i) (pivotColumns i) (pivot_ne i))
    (D B : ℝ → Prop) (pivot_mem : ∀ i, D (pivotRows i))
    (columnMax : ∀ i x, D x → |residual i x (pivotColumns i)| ≤
      |residual i (pivotRows i) (pivotColumns i)|)
    (n : ℕ) (hn : 0 < n) (epsilon : ℝ) (hepsilon : 0 < epsilon)
    (budget : ∀ rows columns : Fin n → ℝ, (∀ i, D (rows i)) →
      (∀ j, B (columns j)) → |(minorMatrix K rows columns).det| ≤ epsilon ^ n)
    (length : ℕ) (S : Finset ℕ) (hS : S ⊆ Finset.range length)
    (hSB : ∀ i ∈ S, B (pivotColumns i))
    (hSlarge : ∀ i ∈ S, epsilon < |residual i (pivotRows i) (pivotColumns i)|) :
    S.card < n := by
  classical
  by_contra hnot
  obtain ⟨T, hTS, hTcard⟩ := Finset.exists_subset_card_eq (Nat.le_of_not_gt hnot)
  have hT : T ⊆ Finset.range length := hTS.trans hS
  have hm : ∀ i, decide (i ∈ T) = true → B (pivotColumns i) := by
    intro i hi
    exact hSB i (hTS (of_decide_eq_true hi))
  have hh := markedPivotProduct_le_initialMinor regular residual pivotRows pivotColumns
    initial realized pivot_ne updates D B pivot_mem columnMax (fun i ↦ decide (i ∈ T)) hm length
  rw [markedCount_of_finset T length hT, markedPivotProduct_of_finset _ T length hT,
    hTcard] at hh
  obtain ⟨rows, columns, hD, hB, hprod⟩ := hh
  have hneT : T.Nonempty := Finset.card_pos.mp (by omega)
  have hstrict : epsilon ^ n < ∏ i ∈ T, |residual i (pivotRows i) (pivotColumns i)| := by
    rw [← hTcard, ← Finset.prod_const]
    exact Finset.prod_lt_prod_of_nonempty (fun _ _ ↦ hepsilon)
      (fun i hi ↦ hSlarge i (hTS hi)) hneT
  exact (not_lt_of_ge (hprod.trans (budget rows columns hD hB))) hstrict

end GECPKernelStructure.GECP
