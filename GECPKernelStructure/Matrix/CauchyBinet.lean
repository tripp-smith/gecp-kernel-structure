import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Data.Fin.Tuple.Sort

namespace Matrix

open Equiv Equiv.Perm Finset Function

variable {R n m : Type*} [CommRing R] [Fintype n] [DecidableEq n]
  [Fintype m] [DecidableEq m]

local notation "ε " σ:arg => ((Equiv.Perm.sign σ : ℤ) : R)

omit [Fintype m] [DecidableEq m] in
private theorem det_mul_rect_noninjective_aux (A : Matrix n m R) (B : Matrix m n R)
    {p : n → m} (hp : ¬Injective p) :
    (∑ σ : Perm n, ε σ * ∏ i, A (σ i) (p i) * B (p i) i) = 0 := by
  obtain ⟨i, j, hpij, hij⟩ : ∃ i j, p i = p j ∧ i ≠ j := by
    exact Function.not_injective_iff.mp hp
  exact
    sum_involution (fun σ _ => σ * Equiv.swap i j)
      (fun σ _ => by
        have hProduct : (∏ x, A (σ x) (p x)) =
            ∏ x, A ((σ * Equiv.swap i j) x) (p x) :=
          Fintype.prod_equiv (Equiv.swap i j) _ _ (by simp [Equiv.apply_swap_eq_self hpij])
        simp [hProduct, sign_swap hij, -sign_swap', Finset.prod_mul_distrib])
      (fun σ _ _ => (not_congr mul_swap_eq_iff).mpr hij)
      (fun _ _ => Finset.mem_univ _) fun σ _ => Equiv.mul_swap_involutive i j σ

/--
The determinant of a rectangular matrix product, expanded over injective
choices of intermediate indices. Grouping injections with the same range
gives the usual subset-indexed Cauchy--Binet formula.
-/
theorem det_mul_rect_eq_sum_injective (A : Matrix n m R) (B : Matrix m n R) :
    (A * B).det =
      ∑ p : n → m with Injective p,
        (A.submatrix id p).det * ∏ i, B (p i) i := by
  calc
    (A * B).det =
        ∑ p : n → m, ∑ σ : Perm n,
          ε σ * ∏ i, A (σ i) (p i) * B (p i) i := by
      simp only [Matrix.det_apply', Matrix.mul_apply, Finset.prod_univ_sum,
        Finset.mul_sum, Fintype.piFinset_univ]
      rw [Finset.sum_comm]
    _ = ∑ p : n → m with Injective p, ∑ σ : Perm n,
          ε σ * ∏ i, A (σ i) (p i) * B (p i) i := by
      refine (Finset.sum_subset (Finset.filter_subset _ _) fun p _ hp => ?_).symm
      apply det_mul_rect_noninjective_aux A B
      simpa only [Finset.mem_filter, Finset.mem_univ, true_and] using hp
    _ = ∑ p : n → m with Injective p,
          (A.submatrix id p).det * ∏ i, B (p i) i := by
      refine Finset.sum_congr rfl fun p hp => ?_
      rw [Matrix.det_apply']
      simp only [Matrix.submatrix_apply, id_eq, Finset.prod_mul_distrib]
      rw [Finset.sum_mul]
      refine Finset.sum_congr rfl fun σ _ => ?_
      ring

/--
An injective finite tuple is uniquely its increasing rearrangement followed by
a permutation. This is the combinatorial grouping used by Cauchy--Binet.
-/
noncomputable def injectiveFunctionEquivStrictMonoPerm {a b : ℕ} :
    {p : Fin a → Fin b // Injective p} ≃
      {q : Fin a → Fin b // StrictMono q} × Equiv.Perm (Fin a) where
  toFun p :=
    let σ := Tuple.sort p.1
    (⟨p.1 ∘ σ, (Tuple.monotone_sort p.1).strictMono_of_injective
      (p.2.comp σ.injective)⟩, σ⁻¹)
  invFun qσ := ⟨qσ.1.1 ∘ qσ.2, qσ.1.2.injective.comp qσ.2.injective⟩
  left_inv p := by
    apply Subtype.ext
    funext i
    simp [Function.comp_apply]
  right_inv qσ := by
    rcases qσ with ⟨q, τ⟩
    apply Prod.ext
    · apply Subtype.ext
      exact Tuple.comp_perm_comp_sort_eq_comp_sort.trans <| by
        rw [Tuple.sort_eq_refl_iff_monotone.mpr q.2.monotone]
        rfl
    · have hSort : Tuple.sort (q.1 ∘ τ) = τ⁻¹ := by
        apply (Tuple.eq_sort_iff.mpr ⟨?_, ?_⟩).symm
        · simpa [Function.comp_assoc] using q.2.monotone
        · intro i j hij hEq
          exact (hij.ne (q.2.injective (by
            simpa [Function.comp_apply] using hEq))).elim
      simp [hSort]

/-- Cauchy--Binet indexed by increasing intermediate-index tuples. -/
theorem det_mul_rect {a b : ℕ} (A : Matrix (Fin a) (Fin b) R)
    (B : Matrix (Fin b) (Fin a) R) :
    (A * B).det =
      ∑ q : {q : Fin a → Fin b // StrictMono q},
        (A.submatrix id q.1).det * (B.submatrix q.1 id).det := by
  rw [det_mul_rect_eq_sum_injective]
  rw [← Finset.sum_subtype_eq_sum_filter]
  simp only [Finset.subtype_univ]
  let e := injectiveFunctionEquivStrictMonoPerm (a := a) (b := b)
  calc
    ∑ p : {p : Fin a → Fin b // Injective p},
        (A.submatrix id p.1).det * ∏ i, B (p.1 i) i =
        ∑ qτ : {q : Fin a → Fin b // StrictMono q} × Equiv.Perm (Fin a),
          (A.submatrix id ((e.symm qτ).1)).det *
            ∏ i, B ((e.symm qτ).1 i) i := by
      exact Fintype.sum_equiv e _ _ fun p => by simp
    _ = ∑ q : {q : Fin a → Fin b // StrictMono q}, ∑ τ : Equiv.Perm (Fin a),
          (A.submatrix id (q.1 ∘ τ)).det * ∏ i, B (q.1 (τ i)) i := by
      rw [Fintype.sum_prod_type]
      rfl
    _ = ∑ q : {q : Fin a → Fin b // StrictMono q},
          (A.submatrix id q.1).det * (B.submatrix q.1 id).det := by
      refine Finset.sum_congr rfl fun q _ => ?_
      rw [Matrix.det_apply' (B.submatrix q.1 id), Finset.mul_sum]
      refine Finset.sum_congr rfl fun τ _ => ?_
      have hPermute : A.submatrix id (q.1 ∘ τ) =
          (A.submatrix id q.1).submatrix id τ := by
        ext i j
        rfl
      rw [hPermute, Matrix.det_permute']
      simp only [Matrix.submatrix_apply, id_eq]
      ring

end Matrix
