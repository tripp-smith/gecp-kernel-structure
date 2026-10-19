import GECPKernelStructure.Matrix.TotalPositiveHadamard
import Mathlib.Data.Fin.Rev

namespace Matrix

theorem StrictlyTotallyPositive.submatrix {n m : ℕ}
    {A : Matrix (Fin n) (Fin n) ℝ} (positive : StrictlyTotallyPositive A)
    (rows columns : Fin m → Fin n) (hr : StrictMono rows) (hc : StrictMono columns) :
    StrictlyTotallyPositive (A.submatrix rows columns) := by
  intro k r c hrr hcc
  exact positive k (rows ∘ r) (columns ∘ c) (hr.comp hrr) (hc.comp hcc)

theorem StrictlyTotallyPositive.reverse {n : ℕ}
    {A : Matrix (Fin n) (Fin n) ℝ} (positive : StrictlyTotallyPositive A) :
    StrictlyTotallyPositive (A.submatrix Fin.rev Fin.rev) := by
  intro m rows columns hr hc
  let r : Fin m → Fin n := fun i ↦ Fin.rev (rows (Fin.rev i))
  let c : Fin m → Fin n := fun i ↦ Fin.rev (columns (Fin.rev i))
  have hrr : StrictMono r := by
    intro i j hij
    exact Fin.rev_lt_rev.mpr (hr (Fin.rev_lt_rev.mpr hij))
  have hcc : StrictMono c := by
    intro i j hij
    exact Fin.rev_lt_rev.mpr (hc (Fin.rev_lt_rev.mpr hij))
  have hpos := positive m r c hrr hcc
  have heq : A.submatrix r c =
      ((A.submatrix Fin.rev Fin.rev).submatrix rows columns).submatrix
        Fin.revPerm Fin.revPerm := by rfl
  rw [heq, Matrix.det_submatrix_equiv_self] at hpos
  exact hpos

theorem det_eq_northwest_pivot_mul_schur {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ) (hne : A 0 0 ≠ 0) :
    A.det = A 0 0 * (northwestSchur A hne).det := by
  simpa only [id_eq, Fin.cons_zero_succ, Matrix.submatrix_id_id] using
    det_cons_succ_eq_pivot_mul_northwestSchur_det A hne id id

private theorem fischer_positive_index_step {n : ℕ}
    (smaller : ∀ B : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ,
      StrictlyTotallyPositive B → ∀ i : Fin (n + 1),
      B.det ≤ B i i * (B.submatrix i.succAbove i.succAbove).det)
    (A : Matrix (Fin (n + 2)) (Fin (n + 2)) ℝ)
    (positive : StrictlyTotallyPositive A) (k : Fin (n + 1)) :
    A.det ≤ A k.succ k.succ *
      (A.submatrix k.succ.succAbove k.succ.succAbove).det := by
  let hp : A 0 0 ≠ 0 := (positive.entry_pos 0 0).ne'
  let S := northwestSchur A hp
  let C := A.submatrix k.succ.succAbove k.succ.succAbove
  have hS : StrictlyTotallyPositive S := positive.northwestSchur
  have hC : StrictlyTotallyPositive C :=
    positive.submatrix _ _ (Fin.strictMono_succAbove _) (Fin.strictMono_succAbove _)
  have c00 : C 0 0 = A 0 0 := by simp [C]
  have hcne : C 0 0 ≠ 0 := (hC.entry_pos 0 0).ne'
  have schur_eq : northwestSchur C hcne = S.submatrix k.succAbove k.succAbove := by
    ext i j
    simp [northwestSchur, GECPKernelStructure.GECP.residualUpdate, C, S,
      ← Fin.succ_succAbove_succ]
  have factorC : C.det = A 0 0 * (S.submatrix k.succAbove k.succAbove).det := by
    rw [det_eq_northwest_pivot_mul_schur C hcne, c00, schur_eq]
  have hminor : 0 ≤ (S.submatrix k.succAbove k.succAbove).det :=
    (hS n _ _ (Fin.strictMono_succAbove _) (Fin.strictMono_succAbove _)).le
  calc
    A.det = A 0 0 * S.det := det_eq_northwest_pivot_mul_schur A hp
    _ ≤ A 0 0 * (S k k * (S.submatrix k.succAbove k.succAbove).det) :=
      mul_le_mul_of_nonneg_left (smaller S hS k) (positive.entry_pos 0 0).le
    _ ≤ A 0 0 * (A k.succ k.succ * (S.submatrix k.succAbove k.succAbove).det) := by
      apply mul_le_mul_of_nonneg_left _ (positive.entry_pos 0 0).le
      exact mul_le_mul_of_nonneg_right (northwestSchur_diagonal_le positive k) hminor
    _ = A k.succ k.succ * C.det := by rw [factorC]; ring

/-- Singleton Hadamard--Fischer for strictly totally positive matrices.
No symmetry or positive-definite hypothesis is required. -/
theorem det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive :
    ∀ {n : ℕ} (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ),
      StrictlyTotallyPositive A → ∀ i : Fin (n + 1),
      A.det ≤ A i i * (A.submatrix i.succAbove i.succAbove).det
  | 0, A, _, i => by fin_cases i; simp
  | n + 1, A, positive, i => by
    have step := fischer_positive_index_step (n := n)
      (fun B hB j ↦ det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive B hB j)
    refine Fin.cases ?_ (fun k ↦ step A positive k) i
    let B := A.submatrix Fin.rev Fin.rev
    have hB : StrictlyTotallyPositive B := positive.reverse
    have hb := step B hB (Fin.last n)
    have hdet : B.det = A.det := Matrix.det_submatrix_equiv_self Fin.revPerm A
    have hdiag : B (Fin.last n).succ (Fin.last n).succ = A 0 0 := by
      simp [B]
    have hcomp :
        (B.submatrix (Fin.last n).succ.succAbove (Fin.last n).succ.succAbove).det =
          (A.submatrix (0 : Fin (n + 2)).succAbove (0 : Fin (n + 2)).succAbove).det := by
      have hm : B.submatrix (Fin.last n).succ.succAbove (Fin.last n).succ.succAbove =
          (A.submatrix Fin.succ Fin.succ).submatrix Fin.revPerm Fin.revPerm := by
        ext j k
        simp [B, Fin.rev_castSucc]
      rw [hm, Matrix.det_submatrix_equiv_self]
      rfl
    rwa [hdet, hdiag, hcomp] at hb

/-- Arbitrary-pivot scalar Schur complement, retaining the ordered complementary
rows and columns. -/
noncomputable def pivotSchur {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (p q : Fin (n + 1)) (hne : A p q ≠ 0) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j ↦ GECPKernelStructure.GECP.residualUpdate A p q hne
    (p.succAbove i) (q.succAbove j)

/-- The absolute determinant identity does not require matching row/column
pivot ranks. -/
theorem abs_det_eq_abs_pivot_mul_pivotSchur {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (p q : Fin (n + 1)) (hne : A p q ≠ 0) :
    |A.det| = |A p q| * |(pivotSchur A p q hne).det| := by
  let B := A.submatrix p.cycleRange.symm q.cycleRange.symm
  have hb : B 0 0 ≠ 0 := by simpa [B] using hne
  have hs : northwestSchur B hb = pivotSchur A p q hne := by
    ext i j
    simp [northwestSchur, pivotSchur, GECPKernelStructure.GECP.residualUpdate, B]
  calc
    |A.det| = |B.det| := (abs_det_submatrix_equiv_equiv _ _ A).symm
    _ = |B 0 0| * |(northwestSchur B hb).det| := by
      rw [det_eq_northwest_pivot_mul_schur B hb, abs_mul]
    _ = |A p q| * |(pivotSchur A p q hne).det| := by simp [B, hs]

/-- A column-maximal pivot cannot enlarge the largest local minor of any
order: delete the row with the pivot column's sorted rank, not necessarily
the pivot row. This distinction is essential for off-diagonal pivots. -/
theorem abs_det_pivotSchur_le_matchedComplement_of_columnMax {n : ℕ}
    {A : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ}
    (positive : StrictlyTotallyPositive A) (p q : Fin (n + 1))
    (columnMax : ∀ i, A i q ≤ A p q) :
    |(pivotSchur A p q (positive.entry_pos p q).ne').det| ≤
      (A.submatrix q.succAbove q.succAbove).det := by
  have hdet : 0 < A.det := by
    simpa using positive (n + 1) id id strictMono_id strictMono_id
  have hminor : 0 ≤ (A.submatrix q.succAbove q.succAbove).det :=
    (positive n _ _ (Fin.strictMono_succAbove _) (Fin.strictMono_succAbove _)).le
  have hf := det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive
    A positive q
  have hfactor := abs_det_eq_abs_pivot_mul_pivotSchur A p q
    (positive.entry_pos p q).ne'
  rw [abs_of_pos hdet, abs_of_pos (positive.entry_pos p q)] at hfactor
  apply (mul_le_mul_iff_right₀ (positive.entry_pos p q)).mp
  calc
    A p q * |(pivotSchur A p q (positive.entry_pos p q).ne').det| = A.det :=
      hfactor.symm
    _ ≤ A q q * (A.submatrix q.succAbove q.succAbove).det := hf
    _ ≤ A p q * (A.submatrix q.succAbove q.succAbove).det :=
      mul_le_mul_of_nonneg_right (columnMax q) hminor

end Matrix
