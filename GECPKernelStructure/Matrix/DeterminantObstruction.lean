import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.HadamardMatrix
import Mathlib.Tactic

namespace GECPKernelStructure

open scoped Matrix Kronecker

/-- Recursive index type for the order-`2^(k+1)` Sylvester--Hadamard matrix. -/
def SylvesterIndex : ℕ → Type
  | 0 => Fin 2
  | k + 1 => SylvesterIndex k × Fin 2

@[instance_reducible] private def sylvesterIndexFintype : (k : ℕ) → Fintype (SylvesterIndex k)
  | 0 => inferInstanceAs (Fintype (Fin 2))
  | k + 1 => @instFintypeProd (SylvesterIndex k) (Fin 2) (sylvesterIndexFintype k) inferInstance

@[instance_reducible] private def sylvesterIndexDecidableEq : (k : ℕ) → DecidableEq (SylvesterIndex k)
  | 0 => inferInstanceAs (DecidableEq (Fin 2))
  | k + 1 => @instDecidableEqProd (SylvesterIndex k) (Fin 2)
      (sylvesterIndexDecidableEq k) inferInstance

attribute [local instance] sylvesterIndexFintype sylvesterIndexDecidableEq

/-- The real Hadamard matrix of order two. -/
def hadamardTwo : Matrix (Fin 2) (Fin 2) ℝ := !![1, 1; 1, -1]

theorem hadamardTwo_isHadamard : hadamardTwo.IsHadamard := by
  constructor
  · intro i j
    fin_cases i <;> fin_cases j <;> norm_num [hadamardTwo, Unitary.mem_iff]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [hadamardTwo, Matrix.mul_apply, Matrix.conjTranspose_apply,
        Fin.sum_univ_two, Matrix.one_apply]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [hadamardTwo, Matrix.mul_apply, Matrix.conjTranspose_apply,
        Fin.sum_univ_two, Matrix.one_apply]

/-- The Sylvester--Hadamard matrix obtained by repeated Kronecker products. -/
noncomputable def sylvesterHadamard :
    (k : ℕ) → Matrix (SylvesterIndex k) (SylvesterIndex k) ℝ
  | 0 => hadamardTwo
  | k + 1 => sylvesterHadamard k ⊗ₖ hadamardTwo

theorem sylvesterHadamard_isHadamard (k : ℕ) :
    (sylvesterHadamard k).IsHadamard := by
  induction k with
  | zero => exact hadamardTwo_isHadamard
  | succ k ih =>
      rw [sylvesterHadamard]
      exact ih.kronecker hadamardTwo_isHadamard

@[simp] theorem sylvesterIndex_card (k : ℕ) :
    Fintype.card (SylvesterIndex k) = 2 ^ (k + 1) := by
  induction k with
  | zero => exact Fintype.card_fin 2
  | succ k ih =>
      change Fintype.card (SylvesterIndex k × Fin 2) = 2 ^ (k + 1 + 1)
      rw [Fintype.card_prod, ih, Fintype.card_fin]
      simp only [pow_succ]

/-- Exact squared determinant of the order-`2^(k+1)` Sylvester--Hadamard matrix. -/
theorem sylvesterHadamard_det_square (k : ℕ) :
    (sylvesterHadamard k).det ^ 2 =
      ((2 ^ (k + 1) : ℕ) : ℝ) ^ (2 ^ (k + 1)) := by
  have hdet := (sylvesterHadamard_isHadamard k).det_mul_star_det
  simpa [pow_two] using hdet

theorem sylvesterHadamard_abs_apply (k : ℕ) (i j : SylvesterIndex k) :
    |sylvesterHadamard k i j| = 1 := by
  have hunit := (sylvesterHadamard_isHadamard k).apply_mem i j
  rw [Unitary.mem_iff] at hunit
  have hsquare : (sylvesterHadamard k i j) ^ 2 = 1 := by
    simpa [pow_two] using hunit.1
  nlinarith [sq_abs (sylvesterHadamard k i j), abs_nonneg (sylvesterHadamard k i j)]

/-- Bounded entries alone cannot support a determinant bound with a universal constant base. -/
theorem entrywise_determinant_constant_base_obstruction (C : ℕ) :
    ∃ k : ℕ, ∃ A : Matrix (SylvesterIndex k) (SylvesterIndex k) ℝ,
      (∀ i j, |A i j| = 1) ∧
      (C : ℝ) ^ Fintype.card (SylvesterIndex k) < |A.det| := by
  let k := 2 * C + 1
  let N := 2 ^ (k + 1)
  have hC : C < 2 ^ C := C.lt_two_pow_self
  have hbase : C ^ 2 < N := by
    dsimp [N, k]
    rw [show 2 * C + 1 + 1 = C + C + 2 by omega, pow_add, pow_add]
    norm_num
    nlinarith
  have hNpos : 0 < N := by positivity
  have hpow_nat : (C ^ 2) ^ N < N ^ N :=
    Nat.pow_lt_pow_left hbase hNpos.ne'
  have hpow_real : ((C : ℝ) ^ 2) ^ N < (N : ℝ) ^ N := by
    exact_mod_cast hpow_nat
  have hdet := sylvesterHadamard_det_square k
  have hsq : ((C : ℝ) ^ N) ^ 2 < |(sylvesterHadamard k).det| ^ 2 := by
    calc
      ((C : ℝ) ^ N) ^ 2 = ((C : ℝ) ^ 2) ^ N := by
        rw [← pow_mul, ← pow_mul]
        congr 1
        omega
      _ < (N : ℝ) ^ N := hpow_real
      _ = (sylvesterHadamard k).det ^ 2 := by simpa [N] using hdet.symm
      _ = |(sylvesterHadamard k).det| ^ 2 := by rw [sq_abs]
  refine ⟨k, sylvesterHadamard k, sylvesterHadamard_abs_apply k, ?_⟩
  rw [sylvesterIndex_card]
  change (C : ℝ) ^ N < |(sylvesterHadamard k).det|
  exact (sq_lt_sq₀ (by positivity) (abs_nonneg _)).mp hsq

end GECPKernelStructure
