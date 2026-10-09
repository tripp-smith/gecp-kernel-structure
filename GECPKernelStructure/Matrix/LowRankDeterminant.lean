import GECPKernelStructure.Matrix.CauchyBinet
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.AbsoluteValue

namespace Matrix

open Equiv Finset Function

variable {n k : Type*} [Fintype n] [DecidableEq n]
  [Fintype k] [DecidableEq k]

/-- Choose the columns of `E` in `s` and the columns of `B` outside `s`. -/
def columnChoice (s : Finset n) (E B : Matrix n n ℝ) : Matrix n n ℝ :=
  fun i j => if j ∈ s then E i j else B i j

/-- A square matrix factoring through a strictly smaller finite type has zero determinant. -/
theorem det_mul_rect_eq_zero_of_card_lt (L : Matrix n k ℝ) (R : Matrix k n ℝ)
    (hcard : Fintype.card k < Fintype.card n) : (L * R).det = 0 := by
  rw [det_mul_rect_eq_sum_injective]
  apply Finset.sum_eq_zero
  intro p hp
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hp
  exact (Fintype.not_injective_of_card_lt p hcard hp).elim

/-- Multilinearity expands `B + E` over all choices of error columns. -/
theorem det_add_eq_sum_columnChoices (B E : Matrix n n ℝ) :
    (E + B).det = ∑ s : Finset n, (columnChoice s E B).det := by
  rw [Matrix.det_apply']
  simp only [Matrix.add_apply, Fintype.prod_add]
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro s _
  rw [Matrix.det_apply']
  apply Finset.sum_congr rfl
  intro σ _
  simp only [columnChoice]
  rw [Finset.prod_ite]
  rw [Finset.filter_mem_eq_inter, Finset.filter_notMem_eq_sdiff]
  rw [Finset.univ_inter, ← Finset.compl_eq_univ_sdiff]

/-- A determinant bound with an individual absolute bound for each column. -/
theorem abs_det_le_factorial_mul_prod_columnBounds (A : Matrix n n ℝ)
    (bound : n → ℝ) (hbound : ∀ i j, |A i j| ≤ bound j) :
    |A.det| ≤ (Fintype.card n).factorial * ∏ j, bound j := by
  rw [Matrix.det_apply']
  calc
    |∑ σ : Perm n, Equiv.Perm.sign σ * ∏ i, A (σ i) i| ≤
        ∑ σ : Perm n, |Equiv.Perm.sign σ * ∏ i, A (σ i) i| :=
      abs_sum_le_sum_abs _ _
    _ = ∑ _σ : Perm n, ∏ i, |A (_σ i) i| := by
      apply Finset.sum_congr rfl
      intro σ _
      rw [abs_mul, Finset.abs_prod, ← Int.cast_abs, Equiv.Perm.sign_abs,
        Int.cast_one, one_mul]
    _ ≤ ∑ _σ : Perm n, ∏ i, bound i := by
      apply Finset.sum_le_sum
      intro σ _
      exact Finset.prod_le_prod (fun _ _ => abs_nonneg _) fun i _ => hbound (σ i) i
    _ = (Fintype.card n).factorial * ∏ j, bound j := by
      simp [Fintype.card_perm, nsmul_eq_mul]

/-- The left factor that combines the low-rank coordinates with the chosen error columns. -/
def mixedLeft (s : Finset n) (L : Matrix n k ℝ) (E : Matrix n n ℝ) :
    Matrix n (k ⊕ s) ℝ
  | i, Sum.inl a => L i a
  | i, Sum.inr j => E i j.1

/-- The right factor that selects low-rank columns outside `s` and error columns in `s`. -/
def mixedRight (s : Finset n) (R : Matrix k n ℝ) : Matrix (k ⊕ s) n ℝ
  | Sum.inl a, j => if j ∈ s then 0 else R a j
  | Sum.inr a, j => if a.1 = j then 1 else 0

omit [Fintype n] in
private theorem sum_subtype_ite_eq (s : Finset n) (f : n → ℝ) (j : n) :
    (∑ a : s, if a.1 = j then f a.1 else 0) = if j ∈ s then f j else 0 := by
  by_cases hj : j ∈ s
  · rw [Fintype.sum_eq_single (⟨j, hj⟩ : s)]
    · simp [hj]
    · intro b hb
      rw [if_neg]
      intro hbj
      exact hb (Subtype.ext hbj)
  · rw [if_neg hj]
    apply Finset.sum_eq_zero
    intro a _
    rw [if_neg]
    intro haj
    apply hj
    simpa [haj] using a.property

omit [Fintype n] [DecidableEq k] in
/-- Every mixed column choice factors through the low-rank coordinates plus its error columns. -/
theorem mixedLeft_mul_mixedRight (s : Finset n) (L : Matrix n k ℝ)
    (R : Matrix k n ℝ) (E : Matrix n n ℝ) :
    mixedLeft s L E * mixedRight s R = columnChoice s E (L * R) := by
  ext i j
  rw [Matrix.mul_apply, Fintype.sum_sum_type]
  by_cases hj : j ∈ s
  · simp [mixedLeft, mixedRight, columnChoice, hj]
    simpa [hj] using sum_subtype_ite_eq s (E i) j
  · simp [mixedLeft, mixedRight, columnChoice, Matrix.mul_apply, hj]
    simpa [hj] using sum_subtype_ite_eq s (E i) j

/-- A mixed choice with more ambient columns than available coordinates is singular. -/
theorem det_columnChoice_eq_zero_of_card_add_lt (s : Finset n)
    (L : Matrix n k ℝ) (R : Matrix k n ℝ) (E : Matrix n n ℝ)
    (hcard : Fintype.card k + s.card < Fintype.card n) :
    (columnChoice s E (L * R)).det = 0 := by
  rw [← mixedLeft_mul_mixedRight s L R E]
  apply det_mul_rect_eq_zero_of_card_lt
  simpa using hcard

/-- Columnwise error and background bounds give a cardinality-sensitive determinant bound. -/
theorem abs_det_columnChoice_le (s : Finset n) (E B : Matrix n n ℝ)
    (ε C : ℝ) (hE : ∀ i j, |E i j| ≤ ε) (hB : ∀ i j, |B i j| ≤ C) :
    |(columnChoice s E B).det| ≤
      (Fintype.card n).factorial * ε ^ s.card * C ^ (Fintype.card n - s.card) := by
  have h := abs_det_le_factorial_mul_prod_columnBounds (columnChoice s E B)
    (fun j => if j ∈ s then ε else C) (by
      intro i j
      by_cases hj : j ∈ s
      · simpa [columnChoice, hj] using hE i j
      · simpa [columnChoice, hj] using hB i j)
  have hcard : (Finset.univ \ s).card = Fintype.card n - s.card := by
    rw [Finset.card_sdiff]
    simp
  simpa [Finset.prod_ite, Finset.filter_mem_eq_inter,
    Finset.filter_notMem_eq_sdiff, hcard, mul_assoc] using h

/--
If `A` is entrywise close to a matrix factoring through `k`, its determinant is
bounded by the mixed terms containing at least `card n - card k` error columns.
This is the finite-dimensional exterior-power obstruction behind the estimate.
-/
theorem abs_det_le_of_factors_approx (A : Matrix n n ℝ) (L : Matrix n k ℝ)
    (R : Matrix k n ℝ) (ε C : ℝ)
    (happrox : ∀ i j, |A i j - (L * R) i j| ≤ ε)
    (hfactor : ∀ i j, |(L * R) i j| ≤ C) :
    |A.det| ≤ ∑ s : Finset n,
      if Fintype.card n ≤ Fintype.card k + s.card then
        (Fintype.card n).factorial * ε ^ s.card * C ^ (Fintype.card n - s.card)
      else 0 := by
  let E : Matrix n n ℝ := A - L * R
  have hA : A = E + L * R := by
    ext i j
    simp [E]
  have hE : ∀ i j, |E i j| ≤ ε := by
    intro i j
    simpa [E] using happrox i j
  rw [hA, det_add_eq_sum_columnChoices]
  calc
    |∑ s : Finset n, (columnChoice s E (L * R)).det| ≤
        ∑ s : Finset n, |(columnChoice s E (L * R)).det| :=
      abs_sum_le_sum_abs _ _
    _ ≤ ∑ s : Finset n,
        if Fintype.card n ≤ Fintype.card k + s.card then
          (Fintype.card n).factorial * ε ^ s.card * C ^ (Fintype.card n - s.card)
        else 0 := by
      apply Finset.sum_le_sum
      intro s _
      by_cases hs : Fintype.card n ≤ Fintype.card k + s.card
      · rw [if_pos hs]
        exact abs_det_columnChoice_le s E (L * R) ε C hE hfactor
      · rw [if_neg hs, det_columnChoice_eq_zero_of_card_add_lt s L R E (by omega), abs_zero]

/--
A finite set of exact columns augments the effective factorization rank by its
cardinality while every other column retains the approximation error bound.
-/
theorem abs_det_le_of_factors_approx_except (A : Matrix n n ℝ)
    (L : Matrix n k ℝ) (R : Matrix k n ℝ) (exceptional : Finset n)
    (ε C : ℝ)
    (happrox : ∀ i j, j ∉ exceptional → |A i j - (L * R) i j| ≤ ε)
    (hA : ∀ i j, j ∈ exceptional → |A i j| ≤ C)
    (hfactor : ∀ i j, j ∉ exceptional → |(L * R) i j| ≤ C)
    (hε0 : 0 ≤ ε) :
    |A.det| ≤ ∑ s : Finset n,
      if Fintype.card n ≤ Fintype.card k + exceptional.card + s.card then
        (Fintype.card n).factorial * ε ^ s.card * C ^ (Fintype.card n - s.card)
      else 0 := by
  let L' := mixedLeft exceptional L A
  let R' := mixedRight exceptional R
  have hproduct : L' * R' = columnChoice exceptional A (L * R) := by
    exact mixedLeft_mul_mixedRight exceptional L R A
  have happrox' : ∀ i j, |A i j - (L' * R') i j| ≤ ε := by
    intro i j
    rw [hproduct]
    by_cases hj : j ∈ exceptional
    · simpa [columnChoice, hj] using hε0
    · simpa [columnChoice, hj] using happrox i j hj
  have hfactor' : ∀ i j, |(L' * R') i j| ≤ C := by
    intro i j
    rw [hproduct]
    by_cases hj : j ∈ exceptional
    · simpa [columnChoice, hj] using hA i j hj
    · simpa [columnChoice, hj] using hfactor i j hj
  have h := abs_det_le_of_factors_approx A L' R' ε C happrox' hfactor'
  simpa [L', R', Fintype.card_sum, Fintype.card_coe] using h

/-- A coarse closed form of `abs_det_le_of_factors_approx`. -/
theorem abs_det_le_two_pow_mul_of_factors_approx (A : Matrix n n ℝ)
    (L : Matrix n k ℝ) (R : Matrix k n ℝ) (ε C : ℝ)
    (happrox : ∀ i j, |A i j - (L * R) i j| ≤ ε)
    (hfactor : ∀ i j, |(L * R) i j| ≤ C)
    (hcard : Fintype.card k ≤ Fintype.card n)
    (hε0 : 0 ≤ ε) (hε1 : ε ≤ 1) (hC1 : 1 ≤ C) :
    |A.det| ≤ (2 : ℝ) ^ Fintype.card n * (Fintype.card n).factorial *
      ε ^ (Fintype.card n - Fintype.card k) * C ^ Fintype.card k := by
  refine (abs_det_le_of_factors_approx A L R ε C happrox hfactor).trans ?_
  calc
    (∑ s : Finset n,
        if Fintype.card n ≤ Fintype.card k + s.card then
          (Fintype.card n).factorial * ε ^ s.card * C ^ (Fintype.card n - s.card)
        else 0) ≤
        ∑ _s : Finset n, (Fintype.card n).factorial *
          ε ^ (Fintype.card n - Fintype.card k) * C ^ Fintype.card k := by
      apply Finset.sum_le_sum
      intro s _
      by_cases hs : Fintype.card n ≤ Fintype.card k + s.card
      · rw [if_pos hs]
        have hεpow : ε ^ s.card ≤ ε ^ (Fintype.card n - Fintype.card k) :=
          pow_le_pow_of_le_one hε0 hε1 (by omega)
        have hCpow : C ^ (Fintype.card n - s.card) ≤ C ^ Fintype.card k :=
          pow_le_pow_right₀ hC1 (by omega)
        gcongr
      · rw [if_neg hs]
        positivity
    _ = (2 : ℝ) ^ Fintype.card n * (Fintype.card n).factorial *
        ε ^ (Fintype.card n - Fintype.card k) * C ^ Fintype.card k := by
      simp [Fintype.card_finset, nsmul_eq_mul, mul_assoc]

/-- Closed determinant bound with a finite set of exactly represented columns. -/
theorem abs_det_le_two_pow_mul_of_factors_approx_except (A : Matrix n n ℝ)
    (L : Matrix n k ℝ) (R : Matrix k n ℝ) (exceptional : Finset n)
    (ε C : ℝ)
    (happrox : ∀ i j, j ∉ exceptional → |A i j - (L * R) i j| ≤ ε)
    (hA : ∀ i j, j ∈ exceptional → |A i j| ≤ C)
    (hfactor : ∀ i j, j ∉ exceptional → |(L * R) i j| ≤ C)
    (hcard : Fintype.card k + exceptional.card ≤ Fintype.card n)
    (hε0 : 0 ≤ ε) (hε1 : ε ≤ 1) (hC1 : 1 ≤ C) :
    |A.det| ≤ (2 : ℝ) ^ Fintype.card n * (Fintype.card n).factorial *
      ε ^ (Fintype.card n - (Fintype.card k + exceptional.card)) *
        C ^ (Fintype.card k + exceptional.card) := by
  let L' := mixedLeft exceptional L A
  let R' := mixedRight exceptional R
  have hproduct : L' * R' = columnChoice exceptional A (L * R) := by
    exact mixedLeft_mul_mixedRight exceptional L R A
  have happrox' : ∀ i j, |A i j - (L' * R') i j| ≤ ε := by
    intro i j
    rw [hproduct]
    by_cases hj : j ∈ exceptional
    · simpa [columnChoice, hj] using hε0
    · simpa [columnChoice, hj] using happrox i j hj
  have hfactor' : ∀ i j, |(L' * R') i j| ≤ C := by
    intro i j
    rw [hproduct]
    by_cases hj : j ∈ exceptional
    · simpa [columnChoice, hj] using hA i j hj
    · simpa [columnChoice, hj] using hfactor i j hj
  have h := abs_det_le_two_pow_mul_of_factors_approx
    A L' R' ε C happrox' hfactor'
    (by simpa [L', Fintype.card_sum, Fintype.card_coe] using hcard)
    hε0 hε1 hC1
  simpa [L', R', Fintype.card_sum, Fintype.card_coe] using h

end Matrix
