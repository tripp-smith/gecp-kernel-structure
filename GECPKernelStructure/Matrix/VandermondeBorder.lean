import GECPKernelStructure.Fermionic.ExponentialTotalPositivity

namespace Matrix

variable {R : Type*} [CommRing R]

theorem det_vandermonde_cast {n m : ℕ} (h : n = m) (v : Fin m → R) :
    (vandermonde (v ∘ Fin.cast h)).det = (vandermonde v).det := by
  subst h
  rfl

/-- Prepending one node factors off all of its cross differences. -/
theorem det_vandermonde_cons {n : ℕ} (a : R) (v : Fin n → R) :
    (vandermonde (Fin.cons a v)).det =
      (∏ i, (v i - a)) * (vandermonde v).det := by
  rw [det_vandermonde, Fin.prod_univ_succ]
  simp only [Fin.Ioi_zero_eq_map, Finset.prod_map, Fin.coe_succEmb,
    Fin.cons_zero, Fin.cons_succ, Fin.prod_Ioi_succ]
  rw [det_vandermonde]

/-- The exact selected/fresh block factorization. Neither tuple must be
sorted, and no nonzero assumption is needed for this polynomial identity. -/
theorem det_vandermonde_append {n m : ℕ} (u : Fin n → R) (v : Fin m → R) :
    (vandermonde (Fin.append u v)).det =
      (vandermonde u).det * (∏ j, ∏ i, (v j - u i)) * (vandermonde v).det := by
  induction n with
  | zero =>
    have hu : u = Fin.elim0 := funext (fun i ↦ Fin.elim0 i)
    subst u
    rw [Fin.elim0_append, det_vandermonde_cast]
    simp
  | succ n ih =>
    rw [← Fin.cons_self_tail u, Fin.append_cons, det_vandermonde_cast,
      det_vandermonde_cons, ih, det_vandermonde_cons]
    rw [Fin.prod_univ_add]
    simp only [Fin.append_left, Fin.append_right]
    simp_rw [Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, Finset.prod_mul_distrib]
    ring

end Matrix

namespace GECPKernelStructure.GECP

/-- A Vandermonde determinant records exactly the orientation of an injective
real tuple. This does not assume the tuple is already sorted. -/
theorem tupleOrientation_mul_vandermonde_pos {n : ℕ} (v : Fin n → ℝ)
    (hinj : Function.Injective v) :
    0 < tupleOrientation v * (Matrix.vandermonde v).det := by
  let p := Tuple.sort v
  have hm : StrictMono (v ∘ p) :=
    (Tuple.monotone_sort v).strictMono_of_injective (hinj.comp p.injective)
  have hp := Fermionic.vandermonde_det_pos_of_strictMono (v ∘ p) hm
  have heq : Matrix.vandermonde (v ∘ p) =
      (Matrix.vandermonde v).submatrix p id := rfl
  rw [heq, Matrix.det_permute] at hp
  exact hp

/-- Multiplying a positive-minor kernel determinant by its two Vandermonde
determinants removes all tuple-order signs. -/
theorem minor_mul_vandermonde_pos_of_strictTotalPositive
    {K : Kernel ℝ ℝ ℝ} (regular : StrictSignRegular K (fun _ ↦ 1))
    {n : ℕ} (rows columns : Fin n → ℝ)
    (hr : Function.Injective rows) (hc : Function.Injective columns) :
    0 < (minorMatrix K rows columns).det *
      (Matrix.vandermonde rows).det * (Matrix.vandermonde columns).det := by
  have hk := oriented_minor_pos_of_strictSignRegularAtOrder
    (regular n) rows columns hr hc
  simp only [one_mul] at hk
  have hh := mul_pos (mul_pos hk (tupleOrientation_mul_vandermonde_pos rows hr))
    (tupleOrientation_mul_vandermonde_pos columns hc)
  have heq :
      (tupleOrientation rows * tupleOrientation columns * (minorMatrix K rows columns).det) *
          (tupleOrientation rows * (Matrix.vandermonde rows).det) *
        (tupleOrientation columns * (Matrix.vandermonde columns).det) =
      tupleOrientation rows ^ 2 * tupleOrientation columns ^ 2 *
        ((minorMatrix K rows columns).det *
          (Matrix.vandermonde rows).det * (Matrix.vandermonde columns).det) := by ring
  rw [heq, tupleOrientation_sq, tupleOrientation_sq, one_mul, one_mul] at hh
  exact hh

end GECPKernelStructure.GECP
