import GECPKernelStructure.Fermionic.SchurFeatureElimination
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open Finset Polynomial Matrix GECP

/-- Nine positive tail weights forming the Phase BS degree-eight norming set. -/
def rankEightAmplificationNode : Fin 9 → ℝ :=
  ![56, 65, 74, 83, 92, 101, 110, 119, 128]

/-- Exact Lagrange coefficients which recover evaluation at one from the
Phase BS norming set. -/
noncomputable def rankEightLagrangeWeight : Fin 9 → ℝ :=
  ![1271136797320 / 387420489,
    -8739065481575 / 387420489,
    26815762573600 / 387420489,
    -47745138240800 / 387420489,
    53778864502000 / 387420489,
    -39151013357456 / 387420489,
    17959180439200 / 387420489,
    -4739832125600 / 387420489,
    550492313800 / 387420489]

theorem rankEightLagrange_moment (degree : Fin 9) :
    ∑ j, rankEightLagrangeWeight j *
      rankEightAmplificationNode j ^ (degree : ℕ) = 1 := by
  fin_cases degree <;>
    norm_num [rankEightLagrangeWeight, rankEightAmplificationNode,
      Fin.sum_univ_succ]

theorem rankEightLagrange_abs_sum_le :
    ∑ j, abs (rankEightLagrangeWeight j) ≤ (518173 : ℝ) := by
  norm_num [rankEightLagrangeWeight, Fin.sum_univ_succ, abs_of_nonneg,
    abs_of_nonpos]

/-- Evaluation at one is reconstructed exactly from the nine norming values
for every polynomial of degree at most eight. -/
theorem rankEightLagrange_eval_one (q : ℝ[X]) (degree_le : q.natDegree ≤ 8) :
    q.eval 1 = ∑ j, rankEightLagrangeWeight j *
      q.eval (rankEightAmplificationNode j) := by
  have degree_lt : q.natDegree < 9 := by omega
  have eval_node (j : Fin 9) :
      q.eval (rankEightAmplificationNode j) =
        ∑ degree : Fin 9, q.coeff degree *
          rankEightAmplificationNode j ^ (degree : ℕ) := by
    rw [eval_eq_sum_range' degree_lt]
    exact (Fin.sum_univ_eq_sum_range (fun degree : ℕ =>
      q.coeff degree *
        rankEightAmplificationNode j ^ degree) 9).symm
  rw [eval_eq_sum_range' degree_lt, Finset.sum_range]
  simp only [one_pow, mul_one]
  simp_rw [eval_node]
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro degree _
  calc
    q.coeff degree = q.coeff degree * 1 := by ring
    _ = q.coeff degree *
        (∑ j, rankEightLagrangeWeight j *
          rankEightAmplificationNode j ^ (degree : ℕ)) := by
      rw [rankEightLagrange_moment]
    _ = ∑ j, rankEightLagrangeWeight j *
        (q.coeff degree *
          rankEightAmplificationNode j ^ (degree : ℕ)) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro j _
      ring

theorem rankEightAmplificationNode_ge (j : Fin 9) :
    (56 : ℝ) ≤ rankEightAmplificationNode j := by
  fin_cases j <;> norm_num [rankEightAmplificationNode]

/-- A uniform bound at the nine norming nodes controls evaluation at one
without any coefficient-sign or root-location assumption. -/
theorem rankEightLagrange_norming_bound (q : ℝ[X])
    (degree_le : q.natDegree ≤ 8) (m : ℕ) {D : ℝ} (D_nonneg : 0 ≤ D)
    (node_bound : ∀ j, 56 ^ m *
      abs (q.eval (rankEightAmplificationNode j)) ≤ D) :
    56 ^ m * abs (q.eval 1) ≤ 518173 * D := by
  rw [rankEightLagrange_eval_one q degree_le]
  calc
    56 ^ m * abs (∑ j, rankEightLagrangeWeight j *
        q.eval (rankEightAmplificationNode j)) ≤
        56 ^ m * ∑ j, abs (rankEightLagrangeWeight j *
          q.eval (rankEightAmplificationNode j)) := by
      gcongr
      exact Finset.abs_sum_le_sum_abs _ _
    _ = ∑ j, abs (rankEightLagrangeWeight j) *
        (56 ^ m * abs (q.eval (rankEightAmplificationNode j))) := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro j _
      rw [abs_mul]
      ring
    _ ≤ ∑ j, abs (rankEightLagrangeWeight j) * D := by
      apply Finset.sum_le_sum
      intro j _
      exact mul_le_mul_of_nonneg_left (node_bound j) (abs_nonneg _)
    _ = (∑ j, abs (rankEightLagrangeWeight j)) * D := by
      rw [Finset.sum_mul]
    _ ≤ 518173 * D :=
      mul_le_mul_of_nonneg_right rankEightLagrange_abs_sum_le D_nonneg

end Fermionic
end GECPKernelStructure

namespace Matrix

open Finset Polynomial

variable {n : Type*} [Fintype n] [DecidableEq n]

theorem mixedColumnPolynomial_natDegree_le_card (E B : Matrix n n ℝ) :
    (mixedColumnPolynomial E B).natDegree ≤ Fintype.card n := by
  unfold mixedColumnPolynomial
  apply Polynomial.natDegree_sum_le_of_forall_le
  intro s _
  exact (Polynomial.natDegree_C_mul_X_pow_le _ _).trans
    (Finset.card_le_univ s)

/-- Rank-eight vanishing plus nine bounded positive evaluations gives a
sign-free determinant amplification estimate. -/
theorem mixedColumnPolynomial_rankEight_finiteWeight_bound
    (E : Matrix n n ℝ) (L : Matrix n (Fin 8) ℝ)
    (R : Matrix (Fin 8) n ℝ) (eight_le : 8 ≤ Fintype.card n)
    {D : ℝ} (D_nonneg : 0 ≤ D)
    (node_det_bound : ∀ j,
      abs (GECPKernelStructure.Fermionic.rankEightAmplificationNode j • E +
        L * R).det ≤ D) :
    56 ^ (Fintype.card n - 8) * abs (E + L * R).det ≤ 518173 * D := by
  let p := mixedColumnPolynomial E (L * R)
  let m := Fintype.card n - 8
  have p_low_zero (degree : ℕ) (degree_lt : degree < m) :
      p.coeff degree = 0 := by
    rw [mixedColumnPolynomial_coeff]
    apply mixedColumnCoefficient_eq_zero_of_rank degree E L R
    simp only [Fintype.card_fin]
    omega
  have X_pow_dvd : X ^ m ∣ p :=
    Polynomial.X_pow_dvd_iff.mpr p_low_zero
  obtain ⟨q, p_eq⟩ := X_pow_dvd
  have q_degree_le : q.natDegree ≤ 8 := by
    by_cases q_zero : q = 0
    · simp [q_zero]
    · have p_degree_eq : p.natDegree = q.natDegree + m := by
        rw [p_eq, Polynomial.natDegree_X_pow_mul m q_zero]
      have p_degree_le : p.natDegree ≤ Fintype.card n :=
        mixedColumnPolynomial_natDegree_le_card E (L * R)
      omega
  have p_eval (x : ℝ) :
      p.eval x = x ^ m * q.eval x := by
    rw [p_eq, Polynomial.eval_mul, Polynomial.eval_pow,
      Polynomial.eval_X]
  have q_node_bound (j : Fin 9) :
      56 ^ m * abs (q.eval
        (GECPKernelStructure.Fermionic.rankEightAmplificationNode j)) ≤ D := by
    have node_nonneg :
        0 ≤ GECPKernelStructure.Fermionic.rankEightAmplificationNode j :=
      (by linarith
        [GECPKernelStructure.Fermionic.rankEightAmplificationNode_ge j])
    have pow_le :
        (56 : ℝ) ^ m ≤
          GECPKernelStructure.Fermionic.rankEightAmplificationNode j ^ m :=
      pow_le_pow_left₀ (by norm_num)
        (GECPKernelStructure.Fermionic.rankEightAmplificationNode_ge j) m
    calc
      56 ^ m * abs (q.eval
          (GECPKernelStructure.Fermionic.rankEightAmplificationNode j)) ≤
          GECPKernelStructure.Fermionic.rankEightAmplificationNode j ^ m *
            abs (q.eval
              (GECPKernelStructure.Fermionic.rankEightAmplificationNode j)) := by
        gcongr
      _ = abs (p.eval
          (GECPKernelStructure.Fermionic.rankEightAmplificationNode j)) := by
        rw [p_eval, abs_mul, abs_of_nonneg (pow_nonneg node_nonneg m)]
      _ = abs
          (GECPKernelStructure.Fermionic.rankEightAmplificationNode j • E +
            L * R).det := by
        rw [mixedColumnPolynomial_eval]
      _ ≤ D := node_det_bound j
  have norming :=
    GECPKernelStructure.Fermionic.rankEightLagrange_norming_bound
      q q_degree_le m D_nonneg q_node_bound
  have eval_one : q.eval 1 = (E + L * R).det := by
    have at_one := p_eval 1
    rw [mixedColumnPolynomial_eval] at at_one
    simpa using at_one.symm
  simpa [m, eval_one] using norming

/-- At fresh order eighteen, the Phase BS norming constants turn nine
`2^18` amplified determinant envelopes into the half-power bound needed by a
complete-pivot contraction block. -/
theorem rankEight_finiteWeight_det_le_halfPow
    (E : Matrix (Fin 18) (Fin 18) ℝ)
    (L : Matrix (Fin 18) (Fin 8) ℝ)
    (R : Matrix (Fin 8) (Fin 18) ℝ)
    (node_det_bound : ∀ j,
      abs (GECPKernelStructure.Fermionic.rankEightAmplificationNode j • E +
        L * R).det ≤ (2 : ℝ) ^ 18) :
    abs (E + L * R).det ≤ ((1 / 2 : ℝ) ^ 18) := by
  have amplified := mixedColumnPolynomial_rankEight_finiteWeight_bound
    E L R (by norm_num) (D := (2 : ℝ) ^ 18) (by positivity) node_det_bound
  norm_num at amplified ⊢
  linarith [abs_nonneg (E + L * R).det]

end Matrix

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- The canonical Phase BR residual split needs only nine fixed companion
determinant envelopes—not coefficientwise positivity—to yield an order-18
half-power determinant bound. -/
theorem positiveSeriesFermionicKernel_conditioned_det_le_halfPow
    (run : Run positiveSeriesFermionicKernel)
    (rows columns : Fin 18 → PositiveArgument)
    (node_det_bound : ∀ j,
      abs (rankEightAmplificationNode j •
          Matrix.submatrix (run.conditionRows positiveSeriesTailMatrix)
            rows columns +
        Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id *
          Matrix.submatrix positiveSeriesLowRight id columns).det ≤
        (2 : ℝ) ^ 18) :
    abs (Matrix.det fun i j => run.finalResidual (rows i) (columns j)) ≤
      ((1 / 2 : ℝ) ^ 18) := by
  let conditionedTail : Matrix (Fin 18) (Fin 18) ℝ :=
    Matrix.submatrix (run.conditionRows positiveSeriesTailMatrix) rows columns
  let conditionedLeft : Matrix (Fin 18) (Fin 8) ℝ :=
    Matrix.submatrix (run.conditionRows positiveSeriesLowLeft) rows id
  let sampledRight : Matrix (Fin 8) (Fin 18) ℝ :=
    Matrix.submatrix positiveSeriesLowRight id columns
  have determinant_bound := Matrix.rankEight_finiteWeight_det_le_halfPow
    conditionedTail conditionedLeft sampledRight node_det_bound
  have residual_split := run.finalResidual_sample_eq_conditionRows_add_mul
    rows columns positiveSeriesTailMatrix positiveSeriesLowLeft
    positiveSeriesLowRight positiveSeriesFermionicKernel_eq_tail_add_low
  change abs (Matrix.submatrix run.finalResidual rows columns).det ≤ _
  rw [residual_split]
  exact determinant_bound

end Fermionic
end GECPKernelStructure
