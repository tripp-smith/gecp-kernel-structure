import GECPKernelStructure.PositiveDefinite.KernelGeometry

namespace GECPKernelStructure.PositiveDefinite

open Matrix GECP

universe u
variable {X : Type u} {K : Kernel X X ℝ}

/-- Strict positivity is required only for injective samples in the fresh domain. -/
def StrictKernelOn (K : Kernel X X ℝ) (D : X → Prop) : Prop :=
  ∀ (n : ℕ) (nodes : Fin n → X), Function.Injective nodes → (∀ i, D (nodes i)) →
    Matrix.PosDef (fun i j => K (nodes i) (nodes j))

theorem strictKernelOn_of_matrix (hK : Matrix.PosDef K) (D : X → Prop) :
    StrictKernelOn K D := fun _ _ hinj _ => hK.submatrix hinj

theorem StrictKernelOn.mono {D E : X → Prop} (hK : StrictKernelOn K D)
    (hED : ∀ x, E x → D x) : StrictKernelOn K E :=
  fun n nodes hinj hnodes => hK n nodes hinj (fun i => hED _ (hnodes i))

/-- Strict scalar Schur preservation on fresh coordinates uses an injective rectangular map. -/
theorem posDef_diagonalSchur_fresh {n : ℕ} (M : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ)
    (hM : M.PosDef) : Matrix.PosDef (fun i j : Fin n =>
      M i.succ j.succ - M i.succ 0 * M 0 j.succ / M 0 0) := by
  classical
  let E : Matrix (Fin (n + 1)) (Fin n) ℝ :=
    Fin.cons (fun j => -M 0 j.succ / M 0 0) (fun i j => if i = j then 1 else 0)
  have hE : Function.Injective E.mulVec := by
    intro v w hvw
    funext i
    have h := congrFun hvw i.succ
    simpa [Matrix.mulVec, dotProduct, E] using h
  let R : Matrix (Fin (n + 1)) (Fin n) ℝ :=
    fun i j => M i j.succ - M i 0 * M 0 j.succ / M 0 0
  have hME : M * E = R := by
    ext i j
    change (∑ k : Fin (n + 1), M i k * E k j) = R i j
    simp [Fin.sum_univ_succ, E, R, mul_ite]
    ring
  have hzero : ∀ j, R 0 j = 0 := by
    intro j
    simp [R, ne_of_gt (hM.diag_pos (i := 0))]
  have hfactor : Eᴴ * M * E = (fun i j : Fin n =>
      M i.succ j.succ - M i.succ 0 * M 0 j.succ / M 0 0) := by
    rw [Matrix.mul_assoc, hME]
    ext i j
    change (∑ k : Fin (n + 1), star (E k i) * R k j) = _
    simp only [Fin.sum_univ_succ, hzero, mul_zero, zero_add]
    simp [E, ite_mul, R]
  rw [← hfactor]
  exact hM.conjTranspose_mul_mul_same hE

/-- The fresh-domain strict invariant is derived from the current original samples. -/
theorem StrictKernelOn.residualUpdate {D : X → Prop} (hK : StrictKernelOn K D)
    {p : X} (hpD : D p) (hp : K p p ≠ 0) :
    StrictKernelOn (GECP.residualUpdate K p p hp) (fun x => D x ∧ x ≠ p) := by
  intro n nodes hinj hnodes
  have hfresh : p ∉ Set.range nodes := by
    rintro ⟨i, hi⟩
    exact (hnodes i).2 hi
  have hsample := hK (n + 1) (Fin.cons p nodes)
    (Fin.cons_injective_of_injective hfresh hinj)
    (by intro i; refine Fin.cases hpD (fun j => (hnodes j).1) i)
  exact posDef_diagonalSchur_fresh _ hsample

/-- On a strict two-point sample, the off-diagonal square is strictly smaller. -/
theorem StrictKernelOn.abs_sq_lt {D : X → Prop} (hK : StrictKernelOn K D)
    {p q : X} (hp : D p) (hq : D q) (hne : p ≠ q) :
    |K p q| ^ 2 < K p p * K q q := by
  have hinj : Function.Injective (![p, q] : Fin 2 → X) := by
    intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all
  have hM := hK 2 ![p, q] hinj (by intro i; fin_cases i <;> assumption)
  have hdet := hM.det_pos
  have hsym : K q p = K p q := by simpa using hM.isHermitian.apply 0 1
  rw [Matrix.det_fin_two (fun i j : Fin 2 => K (![p, q] i) (![p, q] j))] at hdet
  simpa [hsym, sq_abs, pow_two] using sub_pos.mp hdet

end GECPKernelStructure.PositiveDefinite
