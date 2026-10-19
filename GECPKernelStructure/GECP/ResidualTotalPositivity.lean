import GECPKernelStructure.Matrix.VandermondeBorder
import GECPKernelStructure.Fermionic.ResidualGauge
import GECPKernelStructure.Fermionic.LocalMinorContraction

namespace GECPKernelStructure.GECP

open Matrix

namespace Run

/-- Separate the selected positions from an arbitrary finite border without
changing their individual order. -/
noncomputable def selectedBorderEquiv {K : Kernel ℝ ℝ ℝ} :
    (run : Run K) → (γ : Type) → run.SelectedIndex ⊕ γ ≃ run.AugmentedIndex γ
  | .nil _, γ => Equiv.emptySum (Fin 0) γ
  | .step _ _ _ tail, γ =>
    (Equiv.sumAssoc _ _ _).trans
      (Equiv.sumCongr (Equiv.refl _) (tail.selectedBorderEquiv γ))

theorem augmentedRow_comp_selectedBorderEquiv {K : Kernel ℝ ℝ ℝ}
    (run : Run K) (γ : Type) (rows : γ → ℝ) :
    run.augmentedRow γ rows ∘ run.selectedBorderEquiv γ =
      Sum.elim run.selectedRow rows := by
  induction run with
  | nil =>
    funext i
    rcases i with i | i
    · exact Fin.elim0 i
    · rfl
  | step row column hne tail ih =>
    funext i
    rcases i with (i | i) | i
    · rfl
    · exact congrFun ih (Sum.inl i)
    · exact congrFun ih (Sum.inr i)

theorem augmentedColumn_comp_selectedBorderEquiv {K : Kernel ℝ ℝ ℝ}
    (run : Run K) (γ : Type) (columns : γ → ℝ) :
    run.augmentedColumn γ columns ∘ run.selectedBorderEquiv γ =
      Sum.elim run.selectedColumn columns := by
  induction run with
  | nil =>
    funext i
    rcases i with i | i
    · exact Fin.elim0 i
    · rfl
  | step row column hne tail ih =>
    funext i
    rcases i with (i | i) | i
    · rfl
    · exact congrFun ih (Sum.inl i)
    · exact congrFun ih (Sum.inr i)

/-- A canonical Fin enumeration of the selected rows. -/
noncomputable def selectedRowsFin {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    Fin run.pivots.length → ℝ := run.selectedRow ∘ run.selectedIndexEquivFin.symm

/-- A canonical Fin enumeration of the selected columns. -/
noncomputable def selectedColumnsFin {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    Fin run.pivots.length → ℝ := run.selectedColumn ∘ run.selectedIndexEquivFin.symm

noncomputable def augmentedFinEquiv {K : Kernel ℝ ℝ ℝ} (run : Run K) (n : ℕ) :
    Fin (run.pivots.length + n) ≃ run.AugmentedIndex (Fin n) :=
  finSumFinEquiv.symm.trans
    ((Equiv.sumCongr run.selectedIndexEquivFin.symm (Equiv.refl _)).trans
      (run.selectedBorderEquiv (Fin n)))

theorem augmentedRow_comp_augmentedFinEquiv {K : Kernel ℝ ℝ ℝ}
    (run : Run K) {n : ℕ} (rows : Fin n → ℝ) :
    run.augmentedRow (Fin n) rows ∘ run.augmentedFinEquiv n =
      Fin.append run.selectedRowsFin rows := by
  funext i
  refine Fin.addCases (fun i ↦ ?_) (fun i ↦ ?_) i
  · simpa [augmentedFinEquiv, Function.comp_def, Fin.append_left, selectedRowsFin]
      using congrFun (run.augmentedRow_comp_selectedBorderEquiv (Fin n) rows)
        (Sum.inl (run.selectedIndexEquivFin.symm i))
  · simpa [augmentedFinEquiv, Function.comp_def, Fin.append_right]
      using congrFun (run.augmentedRow_comp_selectedBorderEquiv (Fin n) rows) (Sum.inr i)

theorem augmentedColumn_comp_augmentedFinEquiv {K : Kernel ℝ ℝ ℝ}
    (run : Run K) {n : ℕ} (columns : Fin n → ℝ) :
    run.augmentedColumn (Fin n) columns ∘ run.augmentedFinEquiv n =
      Fin.append run.selectedColumnsFin columns := by
  funext i
  refine Fin.addCases (fun i ↦ ?_) (fun i ↦ ?_) i
  · simpa [augmentedFinEquiv, Function.comp_def, Fin.append_left, selectedColumnsFin]
      using congrFun (run.augmentedColumn_comp_selectedBorderEquiv (Fin n) columns)
        (Sum.inl (run.selectedIndexEquivFin.symm i))
  · simpa [augmentedFinEquiv, Function.comp_def, Fin.append_right]
      using congrFun (run.augmentedColumn_comp_selectedBorderEquiv (Fin n) columns) (Sum.inr i)

/-- The multi-bordered identity in explicit selected-then-fresh Fin order. -/
theorem appended_minor_det {K : Kernel ℝ ℝ ℝ} (run : Run K) {n : ℕ}
    (rows columns : Fin n → ℝ) :
    (minorMatrix K (Fin.append run.selectedRowsFin rows)
      (Fin.append run.selectedColumnsFin columns)).det =
      run.selectedCore.det * (minorMatrix run.finalResidual rows columns).det := by
  have heq : minorMatrix K (Fin.append run.selectedRowsFin rows)
      (Fin.append run.selectedColumnsFin columns) =
        (run.augmentedCore (Fin n) rows columns).submatrix
          (run.augmentedFinEquiv n) (run.augmentedFinEquiv n) := by
    ext i j
    rw [← run.augmentedRow_comp_augmentedFinEquiv rows,
      ← run.augmentedColumn_comp_augmentedFinEquiv columns]
    rfl
  rw [heq, Matrix.det_submatrix_equiv_self]
  exact augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor (Fin n) run rows columns

/-- The polynomial whose roots are the selected row coordinates. -/
noncomputable def rowPolynomialWeight {K : Kernel ℝ ℝ ℝ} (run : Run K) (x : ℝ) : ℝ :=
  ∏ i, (x - run.selectedRowsFin i)

/-- The polynomial whose roots are the selected column coordinates. -/
noncomputable def columnPolynomialWeight {K : Kernel ℝ ℝ ℝ} (run : Run K) (y : ℝ) : ℝ :=
  ∏ i, (y - run.selectedColumnsFin i)

theorem selectedRow_ne_of_fresh {K : Kernel ℝ ℝ ℝ} (run : Run K) {x : ℝ}
    (hx : run.FreshRow x) (i : run.SelectedIndex) : run.selectedRow i ≠ x := by
  have hh : Function.Injective
      (Sum.elim run.selectedRow (fun _ : Fin 1 ↦ x)) := by
    rw [← run.augmentedRow_comp_selectedBorderEquiv (Fin 1) (fun _ ↦ x)]
    exact hx.comp (run.selectedBorderEquiv (Fin 1)).injective
  intro heq
  have bad := hh (a₁ := Sum.inl i) (a₂ := Sum.inr 0) heq
  cases bad

theorem selectedColumn_ne_of_fresh {K : Kernel ℝ ℝ ℝ} (run : Run K) {y : ℝ}
    (hy : run.FreshColumn y) (i : run.SelectedIndex) : run.selectedColumn i ≠ y := by
  have hh : Function.Injective
      (Sum.elim run.selectedColumn (fun _ : Fin 1 ↦ y)) := by
    rw [← run.augmentedColumn_comp_selectedBorderEquiv (Fin 1) (fun _ ↦ y)]
    exact hy.comp (run.selectedBorderEquiv (Fin 1)).injective
  intro heq
  have bad := hh (a₁ := Sum.inl i) (a₂ := Sum.inr 0) heq
  cases bad

theorem rowPolynomialWeight_ne_zero_of_fresh {K : Kernel ℝ ℝ ℝ}
    (run : Run K) {x : ℝ} (hx : run.FreshRow x) : run.rowPolynomialWeight x ≠ 0 := by
  apply Finset.prod_ne_zero_iff.mpr
  intro i _
  exact sub_ne_zero.mpr (run.selectedRow_ne_of_fresh hx (run.selectedIndexEquivFin.symm i)).symm

theorem columnPolynomialWeight_ne_zero_of_fresh {K : Kernel ℝ ℝ ℝ}
    (run : Run K) {y : ℝ} (hy : run.FreshColumn y) : run.columnPolynomialWeight y ≠ 0 := by
  apply Finset.prod_ne_zero_iff.mpr
  intro i _
  exact sub_ne_zero.mpr (run.selectedColumn_ne_of_fresh hy (run.selectedIndexEquivFin.symm i)).symm

theorem augmentedRow_injective_of_fresh {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {n : ℕ} (rows : Fin n → ℝ) (hr : Function.Injective rows)
    (hf : ∀ i, run.FreshRow (rows i)) :
    Function.Injective (run.augmentedRow (Fin n) rows) := by
  have hh : Function.Injective (run.augmentedRow (Fin n) rows ∘
      run.selectedBorderEquiv (Fin n)) := by
    rw [run.augmentedRow_comp_selectedBorderEquiv]
    exact run.selectedRow_injective.sumElim hr
      (fun i j ↦ run.selectedRow_ne_of_fresh (hf j) i)
  intro i j heq
  apply (run.selectedBorderEquiv (Fin n)).symm.injective
  apply hh
  simpa only [Function.comp_apply, Equiv.apply_symm_apply] using heq

theorem augmentedColumn_injective_of_fresh {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {n : ℕ} (columns : Fin n → ℝ) (hc : Function.Injective columns)
    (hf : ∀ i, run.FreshColumn (columns i)) :
    Function.Injective (run.augmentedColumn (Fin n) columns) := by
  have hh : Function.Injective (run.augmentedColumn (Fin n) columns ∘
      run.selectedBorderEquiv (Fin n)) := by
    rw [run.augmentedColumn_comp_selectedBorderEquiv]
    exact run.selectedColumn_injective.sumElim hc
      (fun i j ↦ run.selectedColumn_ne_of_fresh (hf j) i)
  intro i j heq
  apply (run.selectedBorderEquiv (Fin n)).symm.injective
  apply hh
  simpa only [Function.comp_apply, Equiv.apply_symm_apply] using heq

/-- Raw polynomial gauging makes every ordered fresh minor positive, using
one pair of row/column weights independent of the minor order. -/
theorem polynomialGaugedFinalResidual_minor_pos {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1)) (run : Run K) {n : ℕ}
    (rows columns : Fin n → ℝ) (hr : StrictMono rows) (hc : StrictMono columns)
    (hfr : ∀ i, run.FreshRow (rows i)) (hfc : ∀ i, run.FreshColumn (columns i)) :
    0 < (∏ i, run.rowPolynomialWeight (rows i)) *
      (∏ j, run.columnPolynomialWeight (columns j)) *
        (minorMatrix run.finalResidual rows columns).det := by
  have hir := run.augmentedRow_injective_of_fresh rows hr.injective hfr
  have hic := run.augmentedColumn_injective_of_fresh columns hc.injective hfc
  have har : Function.Injective (Fin.append run.selectedRowsFin rows) := by
    rw [← run.augmentedRow_comp_augmentedFinEquiv rows]
    exact hir.comp (run.augmentedFinEquiv n).injective
  have hac : Function.Injective (Fin.append run.selectedColumnsFin columns) := by
    rw [← run.augmentedColumn_comp_augmentedFinEquiv columns]
    exact hic.comp (run.augmentedFinEquiv n).injective
  have ha := minor_mul_vandermonde_pos_of_strictTotalPositive regular
    (Fin.append run.selectedRowsFin rows) (Fin.append run.selectedColumnsFin columns) har hac
  rw [run.appended_minor_det, Matrix.det_vandermonde_append,
    Matrix.det_vandermonde_append] at ha
  have hs := minor_mul_vandermonde_pos_of_strictTotalPositive regular
    run.selectedRowsFin run.selectedColumnsFin
    (run.selectedRow_injective.comp run.selectedIndexEquivFin.symm.injective)
    (run.selectedColumn_injective.comp run.selectedIndexEquivFin.symm.injective)
  change 0 < (run.selectedCore.submatrix run.selectedIndexEquivFin.symm
    run.selectedIndexEquivFin.symm).det * _ * _ at hs
  rw [Matrix.det_submatrix_equiv_self] at hs
  have hf := mul_pos (Fermionic.vandermonde_det_pos_of_strictMono rows hr)
    (Fermionic.vandermonde_det_pos_of_strictMono columns hc)
  have heq :
      run.selectedCore.det * (minorMatrix run.finalResidual rows columns).det *
        ((vandermonde run.selectedRowsFin).det * (∏ j, ∏ i, (rows j - run.selectedRowsFin i)) *
          (vandermonde rows).det) *
        ((vandermonde run.selectedColumnsFin).det * (∏ j, ∏ i, (columns j - run.selectedColumnsFin i)) *
          (vandermonde columns).det) =
      (run.selectedCore.det * (vandermonde run.selectedRowsFin).det *
          (vandermonde run.selectedColumnsFin).det) *
        ((∏ i, run.rowPolynomialWeight (rows i)) *
          (∏ j, run.columnPolynomialWeight (columns j)) *
            (minorMatrix run.finalResidual rows columns).det) *
        ((vandermonde rows).det * (vandermonde columns).det) := by
    unfold rowPolynomialWeight columnPolynomialWeight
    ring
  rw [heq] at ha
  exact (mul_pos_iff_of_pos_left hs).mp ((mul_pos_iff_of_pos_right hf).mp ha)

/-- The same magnitude-one row/column signing is used at every minor order. -/
noncomputable def gaugedFinalResidual {K : Kernel ℝ ℝ ℝ} (run : Run K) :
    Kernel {x // run.FreshRow x} {y // run.FreshColumn y} ℝ :=
  fun x y ↦ Real.sign (run.rowPolynomialWeight x) * run.finalResidual x y *
    Real.sign (run.columnPolynomialWeight y)

private theorem det_diagonalScale {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ)
    (a b : Fin n → ℝ) :
    (Matrix.of fun i j ↦ a i * A i j * b j).det =
      (∏ i, a i) * (∏ j, b j) * A.det := by
  have heq : (Matrix.of fun i j ↦ a i * A i j * b j) =
      (Matrix.of fun i j ↦ a i * (b j * A i j)) := by
    ext i j
    simp only [Matrix.of_apply]
    ring
  calc
    _ = (∏ i, a i) * (Matrix.of fun i j ↦ b j * A i j).det := by
      rw [heq]
      exact Matrix.det_mul_column a (Matrix.of fun i j ↦ b j * A i j)
    _ = _ := by rw [Matrix.det_mul_row b A]; ring

private theorem abs_mul_sign (x : ℝ) : |x| * Real.sign x = x := by
  rcases lt_trichotomy x 0 with hn | rfl | hp
  · rw [abs_of_neg hn, Real.sign_of_neg hn]; ring
  · simp
  · rw [abs_of_pos hp, Real.sign_of_pos hp]; ring

private theorem abs_sign_of_ne (x : ℝ) (hx : x ≠ 0) : |Real.sign x| = 1 := by
  rcases Real.sign_apply_eq_of_ne_zero x hx with hn | hp
  · rw [hn]; norm_num
  · rw [hp]; norm_num

/-- One fixed magnitude-preserving signing makes the residual strictly totally
positive at all orders on fresh coordinates, after any successful prefix. -/
theorem gaugedFinalResidual_strictSignRegular {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1)) (run : Run K) :
    StrictSignRegular run.gaugedFinalResidual (fun _ ↦ 1) := by
  intro n rows columns hr hc
  simp only [one_mul]
  let r : Fin n → ℝ := fun i ↦ (rows i : ℝ)
  let c : Fin n → ℝ := fun j ↦ (columns j : ℝ)
  have hr' : StrictMono r := hr
  have hc' : StrictMono c := hc
  have hraw := run.polynomialGaugedFinalResidual_minor_pos regular r c hr' hc'
    (fun i ↦ (rows i).property) (fun j ↦ (columns j).property)
  let a : Fin n → ℝ := fun i ↦ |run.rowPolynomialWeight (r i)|
  let b : Fin n → ℝ := fun j ↦ |run.columnPolynomialWeight (c j)|
  have ha : 0 < ∏ i, a i := Finset.prod_pos fun i _ ↦
    abs_pos.mpr (run.rowPolynomialWeight_ne_zero_of_fresh (rows i).property)
  have hb : 0 < ∏ j, b j := Finset.prod_pos fun j _ ↦
    abs_pos.mpr (run.columnPolynomialWeight_ne_zero_of_fresh (columns j).property)
  have heq : (Matrix.of fun i j ↦ a i *
      (minorMatrix run.gaugedFinalResidual rows columns) i j * b j) =
      (Matrix.of fun i j ↦ run.rowPolynomialWeight (r i) *
        run.finalResidual (r i) (c j) * run.columnPolynomialWeight (c j)) := by
    ext i j
    dsimp [a, b, minorMatrix, gaugedFinalResidual, r, c]
    calc
      _ = (|run.rowPolynomialWeight (rows i)| *
          Real.sign (run.rowPolynomialWeight (rows i))) *
        run.finalResidual (rows i) (columns j) *
        (|run.columnPolynomialWeight (columns j)| *
          Real.sign (run.columnPolynomialWeight (columns j))) := by ring
      _ = _ := by rw [abs_mul_sign, abs_mul_sign]
  have hd := congrArg Matrix.det heq
  have hrawdet :
      (Matrix.of fun i j ↦ run.rowPolynomialWeight (r i) *
        run.finalResidual (r i) (c j) * run.columnPolynomialWeight (c j)).det =
      (∏ i, run.rowPolynomialWeight (r i)) *
        (∏ j, run.columnPolynomialWeight (c j)) *
          (minorMatrix run.finalResidual r c).det :=
    det_diagonalScale (minorMatrix run.finalResidual r c) _ _
  rw [det_diagonalScale, hrawdet] at hd
  rw [← hd] at hraw
  exact (mul_pos_iff_of_pos_left (mul_pos ha hb)).mp hraw

theorem gaugedFinalResidual_minor_det {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {n : ℕ} (rows : Fin n → {x // run.FreshRow x})
    (columns : Fin n → {y // run.FreshColumn y}) :
    (minorMatrix run.gaugedFinalResidual rows columns).det =
      (∏ i, Real.sign (run.rowPolynomialWeight (rows i))) *
        (∏ j, Real.sign (run.columnPolynomialWeight (columns j))) *
          (minorMatrix run.finalResidual (fun i ↦ (rows i : ℝ))
            (fun j ↦ (columns j : ℝ))).det :=
  det_diagonalScale _ _ _

/-- Every entry magnitude is unchanged; no complete-pivot assumption is needed. -/
theorem abs_gaugedFinalResidual {K : Kernel ℝ ℝ ℝ} (run : Run K)
    (x : {x // run.FreshRow x}) (y : {y // run.FreshColumn y}) :
    |run.gaugedFinalResidual x y| = |run.finalResidual x y| := by
  simp only [gaugedFinalResidual, abs_mul,
    abs_sign_of_ne _ (run.rowPolynomialWeight_ne_zero_of_fresh x.property),
    abs_sign_of_ne _ (run.columnPolynomialWeight_ne_zero_of_fresh y.property),
    one_mul, mul_one]

/-- The common signing preserves absolute minors at every order, including
unordered tuples, repeated coordinates, and the empty minor. -/
theorem abs_gaugedFinalResidual_minor_det {K : Kernel ℝ ℝ ℝ} (run : Run K)
    {n : ℕ} (rows : Fin n → {x // run.FreshRow x})
    (columns : Fin n → {y // run.FreshColumn y}) :
    |(minorMatrix run.gaugedFinalResidual rows columns).det| =
      |(minorMatrix run.finalResidual (fun i ↦ (rows i : ℝ))
        (fun j ↦ (columns j : ℝ))).det| := by
  rw [gaugedFinalResidual_minor_det, abs_mul, abs_mul]
  have ha : |∏ i, Real.sign (run.rowPolynomialWeight (rows i))| = 1 := by
    rw [Finset.abs_prod]
    exact Finset.prod_eq_one fun i _ ↦
      abs_sign_of_ne _ (run.rowPolynomialWeight_ne_zero_of_fresh (rows i).property)
  have hb : |∏ j, Real.sign (run.columnPolynomialWeight (columns j))| = 1 := by
    rw [Finset.abs_prod]
    exact Finset.prod_eq_one fun j _ ↦
      abs_sign_of_ne _ (run.columnPolynomialWeight_ne_zero_of_fresh (columns j).property)
  rw [ha, hb, one_mul, one_mul]

private theorem abs_det_diagonalScale {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ)
    (a b : Fin n → ℝ) (ha : ∀ i, |a i| = 1) (hb : ∀ j, |b j| = 1) :
    |(Matrix.of fun i j ↦ a i * A i j * b j).det| = |A.det| := by
  rw [det_diagonalScale, abs_mul, abs_mul, Finset.abs_prod, Finset.abs_prod]
  simp only [ha, hb, Finset.prod_const_one, one_mul]

private theorem residualUpdate_diagonalScale {α β : Type*} (R : Kernel α β ℝ)
    (a : α → ℝ) (b : β → ℝ) (p : α) (q : β) (hne : R p q ≠ 0)
    (hscaled : a p * R p q * b q ≠ 0) (ha : a p ≠ 0) (hb : b q ≠ 0)
    (x : α) (y : β) :
    residualUpdate (fun u v ↦ a u * R u v * b v) p q hscaled x y =
      a x * residualUpdate R p q hne x y * b y := by
  dsimp [residualUpdate]
  field_simp

/-- Phase DC's same-column minor bound applies to the actual next residual
after any successful prefix. Only the next pivot must be column-maximal;
none of the preceding pivots is assumed complete or confined to a band. -/
theorem finalResidual_localMinor_schur_le {K : Kernel ℝ ℝ ℝ}
    (regular : StrictSignRegular K (fun _ ↦ 1)) (run : Run K) {n : ℕ}
    (rows : Fin (n + 1) → {x // run.FreshRow x})
    (columns : Fin (n + 1) → {y // run.FreshColumn y})
    (hr : StrictMono rows) (hc : StrictMono columns) (p q : Fin (n + 1))
    (columnMax : ∀ i, |run.finalResidual (rows i) (columns q)| ≤
      |run.finalResidual (rows p) (columns q)|)
    (hne : run.finalResidual (rows p) (columns q) ≠ 0) :
    |(minorMatrix (residualUpdate run.finalResidual (rows p) (columns q) hne)
        (fun i ↦ (rows (p.succAbove i) : ℝ))
        (fun j ↦ (columns (q.succAbove j) : ℝ))).det| ≤
      |(minorMatrix run.finalResidual (fun i ↦ (rows (q.succAbove i) : ℝ))
        (fun j ↦ (columns (q.succAbove j) : ℝ))).det| := by
  have hg := run.gaugedFinalResidual_strictSignRegular regular
  have hpos (i j : Fin (n + 1)) : 0 < run.gaugedFinalResidual (rows i) (columns j) := by
    have hh := hg 1 (fun _ ↦ rows i) (fun _ ↦ columns j)
      (by simp [Fin.strictMono_iff_lt_succ]) (by simp [Fin.strictMono_iff_lt_succ])
    simpa [minorMatrix] using hh
  have hvalue (i j : Fin (n + 1)) :
      run.gaugedFinalResidual (rows i) (columns j) =
        |run.finalResidual (rows i) (columns j)| := by
    rw [← run.abs_gaugedFinalResidual (rows i) (columns j), abs_of_pos (hpos i j)]
  have hmax : ∀ i, run.gaugedFinalResidual (rows i) (columns q) ≤
      run.gaugedFinalResidual (rows p) (columns q) := by simpa only [hvalue] using columnMax
  have htne := (hpos p q).ne'
  have bound := localMinor_schur_le_of_strictTotalPositive hg rows columns hr hc p q hmax htne
  have hqp := Fin.strictMono_succAbove q
  have hold := hg n (rows ∘ q.succAbove) (columns ∘ q.succAbove) (hr.comp hqp) (hc.comp hqp)
  simp only [one_mul] at hold
  have holdeq : (minorMatrix run.gaugedFinalResidual
      (rows ∘ q.succAbove) (columns ∘ q.succAbove)).det =
        |(minorMatrix run.finalResidual (fun i ↦ (rows (q.succAbove i) : ℝ))
          (fun j ↦ (columns (q.succAbove j) : ℝ))).det| := by
    exact (abs_of_pos hold).symm.trans
      (run.abs_gaugedFinalResidual_minor_det (rows ∘ q.succAbove) (columns ∘ q.succAbove))
  have hp : Real.sign (run.rowPolynomialWeight (rows p)) ≠ 0 :=
    Real.sign_eq_zero_iff.not.mpr (run.rowPolynomialWeight_ne_zero_of_fresh (rows p).property)
  have hq : Real.sign (run.columnPolynomialWeight (columns q)) ≠ 0 :=
    Real.sign_eq_zero_iff.not.mpr (run.columnPolynomialWeight_ne_zero_of_fresh (columns q).property)
  have heq : minorMatrix (residualUpdate run.gaugedFinalResidual (rows p) (columns q) htne)
      (rows ∘ p.succAbove) (columns ∘ q.succAbove) =
      (Matrix.of fun i j ↦ Real.sign (run.rowPolynomialWeight (rows (p.succAbove i))) *
        (minorMatrix (residualUpdate run.finalResidual (rows p) (columns q) hne)
          (fun i ↦ (rows (p.succAbove i) : ℝ))
          (fun j ↦ (columns (q.succAbove j) : ℝ))) i j *
        Real.sign (run.columnPolynomialWeight (columns (q.succAbove j)))) := by
    ext i j
    exact residualUpdate_diagonalScale
      (fun x y ↦ run.finalResidual (x : ℝ) (y : ℝ))
      (fun x ↦ Real.sign (run.rowPolynomialWeight x))
      (fun y ↦ Real.sign (run.columnPolynomialWeight y))
      (rows p) (columns q) hne htne hp hq (rows (p.succAbove i)) (columns (q.succAbove j))
  rw [heq, abs_det_diagonalScale _ _ _
    (fun i ↦ abs_sign_of_ne _ (run.rowPolynomialWeight_ne_zero_of_fresh (rows (p.succAbove i)).property))
    (fun j ↦ abs_sign_of_ne _ (run.columnPolynomialWeight_ne_zero_of_fresh (columns (q.succAbove j)).property)),
    holdeq] at bound
  exact bound

end Run
end GECPKernelStructure.GECP
