import GECPKernelStructure.Fermionic.SignRegularity
import GECPKernelStructure.Matrix.CauchyBinet
import Mathlib.Algebra.Polynomial.RuleOfSigns
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.Topology.Instances.Matrix
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- A nonzero real polynomial has fewer sign variations than nonzero coefficients. -/
theorem signVariations_lt_card_support (P : Polynomial ℝ) (hP : P ≠ 0) :
    P.signVariations < P.support.card := by
  generalize hc : P.support.card = c
  induction c using Nat.strong_induction_on generalizing P with
  | h c ih =>
      by_cases hErase : P.eraseLead = 0
      · have hCard : P.support.card = 1 :=
          P.card_support_eq_one_of_eraseLead_eq_zero hP hErase
        have hMonomial : Polynomial.monomial P.natDegree P.leadingCoeff = P := by
          simpa only [hErase, zero_add] using
            P.eraseLead_add_monomial_natDegree_leadingCoeff
        rw [← hMonomial, Polynomial.signVariations_monomial]
        omega
      · have hCardLt : P.eraseLead.support.card < c := by
          rw [← hc]
          exact P.eraseLead_support_card_lt hP
        have hInduction := ih _ hCardLt P.eraseLead hErase rfl
        have hVariations := P.signVariations_le_eraseLead_succ
        have hCardStep := P.card_support_eraseLead_add_one hP
        omega

/-- Descartes' rule bounds positive roots by one less than the support size. -/
theorem positiveRoots_lt_card_support (P : Polynomial ℝ) (hP : P ≠ 0) :
    P.roots.countP (0 < ·) < P.support.card :=
  (P.roots_countP_pos_le_signVariations).trans_lt
    (signVariations_lt_card_support P hP)

/-- The generalized Vandermonde matrix with prescribed natural exponents. -/
def generalizedVandermonde {n : ℕ} (nodes : Fin n → ℝ)
    (exponents : Fin n → ℕ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => nodes i ^ exponents j

/-- Distinct positive nodes and distinct exponents give a nonsingular power matrix. -/
theorem generalizedVandermonde_det_ne_zero {n : ℕ} (nodes : Fin n → ℝ)
    (exponents : Fin n → ℕ) (nodes_pos : ∀ i, 0 < nodes i)
    (nodes_inj : Function.Injective nodes)
    (exponents_inj : Function.Injective exponents) :
    (generalizedVandermonde nodes exponents).det ≠ 0 := by
  intro hDet
  obtain ⟨coefficients, coefficients_ne, hKernel⟩ :=
    Matrix.exists_mulVec_eq_zero_iff.mpr hDet
  let P : Polynomial ℝ :=
    ∑ j : Fin n, Polynomial.monomial (exponents j) (coefficients j)
  have hCoeff (j : Fin n) : P.coeff (exponents j) = coefficients j := by
    simp [P, Polynomial.coeff_monomial, exponents_inj.eq_iff]
  have hP : P ≠ 0 := by
    intro hPZero
    apply coefficients_ne
    funext j
    simpa [hPZero] using (hCoeff j).symm
  have hSupport : P.support ⊆ Finset.univ.image exponents := by
    intro degree hDegree
    by_contra hDegreeImage
    have hAll : ∀ j : Fin n, degree ≠ exponents j := by
      intro j hEq
      apply hDegreeImage
      exact Finset.mem_image.mpr ⟨j, Finset.mem_univ j, hEq.symm⟩
    have hCoeffZero : P.coeff degree = 0 := by
      have hAll' : ∀ j : Fin n, exponents j ≠ degree := fun j => (hAll j).symm
      simp [P, Polynomial.coeff_monomial, hAll']
    exact (Polynomial.mem_support_iff.mp hDegree) hCoeffZero
  have hSupportCard : P.support.card ≤ n := by
    calc
      P.support.card ≤ (Finset.univ.image exponents).card :=
        Finset.card_le_card hSupport
      _ ≤ Finset.univ.card := Finset.card_image_le
      _ = n := Fintype.card_fin n
  have hEval (i : Fin n) : P.eval (nodes i) = 0 := by
    have hRow := congrFun hKernel i
    simpa only [P, Polynomial.eval_finsetSum, Polynomial.eval_monomial,
      generalizedVandermonde, Matrix.mulVec, dotProduct, Pi.zero_apply, mul_comm] using hRow
  let positiveRoots := (P.roots.filter (0 < ·)).toFinset
  have hNodesSubset : Finset.univ.image nodes ⊆ positiveRoots := by
    intro x hx
    obtain ⟨i, _, rfl⟩ := Finset.mem_image.mp hx
    change nodes i ∈ (P.roots.filter (0 < ·)).toFinset
    rw [Multiset.mem_toFinset, Multiset.mem_filter]
    exact ⟨(P.mem_roots hP).mpr (hEval i), nodes_pos i⟩
  have hNodeCard : (Finset.univ.image nodes).card = n := by
    rw [Finset.card_image_of_injective _ nodes_inj, Finset.card_univ, Fintype.card_fin]
  have hRootLower : n ≤ P.roots.countP (0 < ·) := by
    rw [Multiset.countP_eq_card_filter]
    calc
      n = (Finset.univ.image nodes).card := hNodeCard.symm
      _ ≤ positiveRoots.card := Finset.card_le_card hNodesSubset
      _ ≤ (P.roots.filter (0 < ·)).card := by
        simpa only [positiveRoots] using
          (Multiset.toFinset_card_le (P.roots.filter (0 < ·)))
  have hRootUpper := positiveRoots_lt_card_support P hP
  omega

/-- An increasing real tuple has a strictly positive Vandermonde determinant. -/
theorem vandermonde_det_pos_of_strictMono {n : ℕ} (values : Fin n → ℝ)
    (values_mono : StrictMono values) :
    0 < (Matrix.vandermonde values).det := by
  rw [Matrix.det_vandermonde]
  refine Finset.prod_pos fun i _ => ?_
  refine Finset.prod_pos fun j hj => ?_
  exact sub_pos.mpr (values_mono (Finset.mem_Ioi.mp hj))

/-- The generalized Vandermonde determinant is positive at geometric nodes. -/
theorem geometric_generalizedVandermonde_det_pos {n : ℕ}
    (exponents : Fin n → ℕ) (exponents_mono : StrictMono exponents) :
    0 < (generalizedVandermonde (fun i : Fin n => (2 : ℝ) ^ i.1) exponents).det := by
  let values : Fin n → ℝ := fun j => (2 : ℝ) ^ exponents j
  have values_mono : StrictMono values :=
    (pow_right_strictMono₀ (by norm_num : (1 : ℝ) < 2)).comp exponents_mono
  have hTranspose :
      (generalizedVandermonde (fun i : Fin n => (2 : ℝ) ^ i.1) exponents)ᵀ =
        Matrix.vandermonde values := by
    ext i j
    simp only [generalizedVandermonde, Matrix.transpose_apply, Matrix.vandermonde_apply,
      values]
    rw [← pow_mul, ← pow_mul, Nat.mul_comm]
  have hPositive := vandermonde_det_pos_of_strictMono values values_mono
  rw [← hTranspose, Matrix.det_transpose] at hPositive
  exact hPositive

/--
The generalized Vandermonde determinant is positive on positive, strictly
ordered nodes with strictly ordered natural exponents.
-/
theorem generalizedVandermonde_det_pos {n : ℕ} (nodes : Fin n → ℝ)
    (exponents : Fin n → ℕ) (nodes_pos : ∀ i, 0 < nodes i)
    (nodes_mono : StrictMono nodes) (exponents_mono : StrictMono exponents) :
    0 < (generalizedVandermonde nodes exponents).det := by
  let geometricNodes : Fin n → ℝ := fun i => (2 : ℝ) ^ i.1
  let nodePath : ℝ → Fin n → ℝ := fun t i =>
    (1 - t) * nodes i + t * geometricNodes i
  let determinantPath : ℝ → ℝ := fun t =>
    (generalizedVandermonde (nodePath t) exponents).det
  have geometricNodes_pos : ∀ i, 0 < geometricNodes i := by
    intro i
    simp only [geometricNodes]
    positivity
  have geometricNodes_mono : StrictMono geometricNodes := by
    exact (pow_right_strictMono₀ (by norm_num : (1 : ℝ) < 2)).comp Fin.val_strictMono
  have nodePath_pos (t : ℝ) (ht : t ∈ Set.Icc (0 : ℝ) 1) :
      ∀ i, 0 < nodePath t i := by
    intro i
    dsimp only [nodePath]
    by_cases htZero : t = 0
    · simp [htZero, nodes_pos i]
    · have htPos : 0 < t := lt_of_le_of_ne ht.1 (Ne.symm htZero)
      exact add_pos_of_nonneg_of_pos
        (mul_nonneg (sub_nonneg.mpr ht.2) (nodes_pos i).le)
        (mul_pos htPos (geometricNodes_pos i))
  have nodePath_mono (t : ℝ) (ht : t ∈ Set.Icc (0 : ℝ) 1) :
      StrictMono (nodePath t) := by
    intro i j hij
    dsimp only [nodePath]
    have hNodes := nodes_mono hij
    have hGeometric := geometricNodes_mono hij
    by_cases htZero : t = 0
    · simpa [htZero] using hNodes
    · have htPos : 0 < t := lt_of_le_of_ne ht.1 (Ne.symm htZero)
      have hNodeTerm : 0 ≤ (1 - t) * (nodes j - nodes i) :=
        mul_nonneg (sub_nonneg.mpr ht.2) (sub_nonneg.mpr hNodes.le)
      have hGeometricTerm : 0 < t * (geometricNodes j - geometricNodes i) :=
        mul_pos htPos (sub_pos.mpr hGeometric)
      nlinarith
  have determinantPath_ne_zero (t : ℝ) (ht : t ∈ Set.Icc (0 : ℝ) 1) :
      determinantPath t ≠ 0 := by
    exact generalizedVandermonde_det_ne_zero (nodePath t) exponents
      (nodePath_pos t ht) (nodePath_mono t ht).injective exponents_mono.injective
  have determinantPath_continuous : Continuous determinantPath := by
    apply Continuous.matrix_det
    apply continuous_pi
    intro i
    apply continuous_pi
    intro j
    dsimp only [determinantPath, generalizedVandermonde, nodePath, geometricNodes]
    fun_prop
  have determinantPath_zero :
      determinantPath 0 = (generalizedVandermonde nodes exponents).det := by
    simp [determinantPath, nodePath]
  have determinantPath_one :
      determinantPath 1 = (generalizedVandermonde geometricNodes exponents).det := by
    simp [determinantPath, nodePath]
  have determinantPath_one_pos : 0 < determinantPath 1 := by
    rw [determinantPath_one]
    exact geometric_generalizedVandermonde_det_pos exponents exponents_mono
  rw [← determinantPath_zero]
  by_contra hNotPositive
  have determinantPath_zero_ne : determinantPath 0 ≠ 0 :=
    determinantPath_ne_zero 0 (by constructor <;> norm_num)
  have determinantPath_zero_neg : determinantPath 0 < 0 := by
    exact lt_of_le_of_ne (le_of_not_gt hNotPositive) determinantPath_zero_ne
  have hZeroBetween : 0 ∈ Set.Icc (determinantPath 0) (determinantPath 1) :=
    ⟨determinantPath_zero_neg.le, determinantPath_one_pos.le⟩
  obtain ⟨t, ht, hPathZero⟩ :=
    (intermediate_value_Icc (by norm_num : (0 : ℝ) ≤ 1)
      determinantPath_continuous.continuousOn) hZeroBetween
  exact determinantPath_ne_zero t ht hPathZero

/--
The generalized Vandermonde determinant remains nonnegative when the least
strictly ordered node is allowed to be zero.
-/
theorem generalizedVandermonde_det_nonneg {n : ℕ} (nodes : Fin n → ℝ)
    (exponents : Fin n → ℕ) (nodes_nonneg : ∀ i, 0 ≤ nodes i)
    (nodes_mono : StrictMono nodes) (exponents_mono : StrictMono exponents) :
    0 ≤ (generalizedVandermonde nodes exponents).det := by
  let shiftedNodes : ℕ → Fin n → ℝ := fun k i =>
    nodes i + ((k + 1 : ℕ) : ℝ)⁻¹
  have shiftedNodes_pos (k : ℕ) : ∀ i, 0 < shiftedNodes k i := by
    intro i
    dsimp only [shiftedNodes]
    exact add_pos_of_nonneg_of_pos (nodes_nonneg i) (by positivity)
  have shiftedNodes_mono (k : ℕ) : StrictMono (shiftedNodes k) := by
    intro i j hij
    dsimp only [shiftedNodes]
    simpa only [add_comm] using
      (add_lt_add_right (nodes_mono hij) (((k + 1 : ℕ) : ℝ)⁻¹))
  have shiftedDet_nonneg (k : ℕ) :
      0 ≤ (generalizedVandermonde (shiftedNodes k) exponents).det :=
    (generalizedVandermonde_det_pos (shiftedNodes k) exponents
      (shiftedNodes_pos k) (shiftedNodes_mono k) exponents_mono).le
  have shiftedNodes_tendsto : Filter.Tendsto shiftedNodes Filter.atTop (nhds nodes) := by
    rw [tendsto_pi_nhds]
    intro i
    simpa only [shiftedNodes, Nat.cast_add, Nat.cast_one, one_div, add_zero] using
      (tendsto_const_nhds.add tendsto_one_div_add_atTop_nhds_zero_nat :
        Filter.Tendsto (fun k : ℕ => nodes i + 1 / ((k : ℝ) + 1))
          Filter.atTop (nhds (nodes i + 0)))
  have determinant_continuous : Continuous fun values : Fin n → ℝ =>
      (generalizedVandermonde values exponents).det := by
    apply Continuous.matrix_det
    apply continuous_pi
    intro i
    apply continuous_pi
    intro j
    exact (continuous_apply i).pow (exponents j)
  have shiftedDet_tendsto :
      Filter.Tendsto
        (fun k => (generalizedVandermonde (shiftedNodes k) exponents).det)
        Filter.atTop (nhds (generalizedVandermonde nodes exponents).det) :=
    determinant_continuous.continuousAt.tendsto.comp shiftedNodes_tendsto
  exact ge_of_tendsto shiftedDet_tendsto (Filter.Eventually.of_forall shiftedDet_nonneg)

/-- The row-side weighted monomial features in a finite exponential series. -/
noncomputable def expTaylorLeft {n terms : ℕ} (rows : Fin n → ℝ) :
    Matrix (Fin n) (Fin terms) ℝ :=
  fun i k => rows i ^ k.1 * (k.1.factorial : ℝ)⁻¹

/-- The column-side monomial features in a finite exponential series. -/
noncomputable def expTaylorRight {n terms : ℕ} (columns : Fin n → ℝ) :
    Matrix (Fin terms) (Fin n) ℝ :=
  fun k j => columns j ^ k.1

/-- A finite Taylor truncation of the kernel `exp(row * column)`. -/
noncomputable def expTaylorMatrix {n terms : ℕ} (rows columns : Fin n → ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  expTaylorLeft (terms := terms) rows * expTaylorRight (terms := terms) columns

/-- The finite feature product evaluates the truncated exponential series. -/
theorem expTaylorMatrix_apply {n terms : ℕ} (rows columns : Fin n → ℝ)
    (i j : Fin n) :
    expTaylorMatrix (terms := terms) rows columns i j =
      ∑ k : Fin terms,
        rows i ^ k.1 * (k.1.factorial : ℝ)⁻¹ * columns j ^ k.1 := by
  simp only [expTaylorMatrix, expTaylorLeft, expTaylorRight, Matrix.mul_apply]

private theorem expTaylorLeftMinor_det_nonneg {n terms : ℕ}
    (rows : Fin n → ℝ) (q : Fin n → Fin terms)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (rows_mono : StrictMono rows)
    (q_mono : StrictMono q) :
    0 ≤ ((expTaylorLeft (terms := terms) rows).submatrix id q).det := by
  let exponents : Fin n → ℕ := fun j => (q j).1
  have exponents_mono : StrictMono exponents := Fin.val_strictMono.comp q_mono
  have hFactor :
      (expTaylorLeft (terms := terms) rows).submatrix id q =
        generalizedVandermonde rows exponents *
          Matrix.diagonal (fun j : Fin n => ((exponents j).factorial : ℝ)⁻¹) := by
    ext i j
    simp [expTaylorLeft, generalizedVandermonde, exponents]
  rw [hFactor, Matrix.det_mul, Matrix.det_diagonal]
  exact mul_nonneg
    (generalizedVandermonde_det_nonneg rows exponents rows_nonneg rows_mono
      exponents_mono)
    (Finset.prod_nonneg fun j _ => inv_nonneg.mpr (Nat.cast_nonneg _))

private theorem expTaylorRightMinor_det_nonneg {n terms : ℕ}
    (columns : Fin n → ℝ) (q : Fin n → Fin terms)
    (columns_nonneg : ∀ i, 0 ≤ columns i) (columns_mono : StrictMono columns)
    (q_mono : StrictMono q) :
    0 ≤ ((expTaylorRight (terms := terms) columns).submatrix q id).det := by
  let exponents : Fin n → ℕ := fun i => (q i).1
  have exponents_mono : StrictMono exponents := Fin.val_strictMono.comp q_mono
  have hTranspose :
      ((expTaylorRight (terms := terms) columns).submatrix q id)ᵀ =
        generalizedVandermonde columns exponents := by
    ext i j
    rfl
  calc
    0 ≤ (((expTaylorRight (terms := terms) columns).submatrix q id)ᵀ).det := by
      rw [hTranspose]
      exact generalizedVandermonde_det_nonneg columns exponents columns_nonneg
        columns_mono exponents_mono
    _ = ((expTaylorRight (terms := terms) columns).submatrix q id).det :=
      Matrix.det_transpose _

/--
A finite exponential Taylor feature matrix has positive determinant once it
contains at least the first `n` powers. The proof is an ordered
Cauchy--Binet expansion: every generalized Vandermonde summand is nonnegative,
and the summand indexed by powers `0, ..., n - 1` is strictly positive.
-/
theorem expTaylorMatrix_det_pos {n terms : ℕ} (rows columns : Fin n → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (hterms : n ≤ terms) :
    0 < (expTaylorMatrix (terms := terms) rows columns).det := by
  let principal : Fin n → Fin terms := fun i =>
    ⟨i.1, lt_of_lt_of_le i.2 hterms⟩
  have principal_mono : StrictMono principal := by
    intro i j hij
    exact hij
  have leftFactor :
      (expTaylorLeft (terms := terms) rows).submatrix id principal =
        Matrix.vandermonde rows *
          Matrix.diagonal (fun j : Fin n => ((j.1.factorial : ℝ)⁻¹)) := by
    ext i j
    simp [expTaylorLeft, principal, Matrix.vandermonde_apply]
  have leftPrincipal_pos :
      0 < ((expTaylorLeft (terms := terms) rows).submatrix id principal).det := by
    rw [leftFactor, Matrix.det_mul, Matrix.det_diagonal]
    exact mul_pos (vandermonde_det_pos_of_strictMono rows rows_mono)
      (Finset.prod_pos fun j _ => inv_pos.mpr (Nat.cast_pos.mpr j.1.factorial_pos))
  have rightTranspose :
      ((expTaylorRight (terms := terms) columns).submatrix principal id)ᵀ =
        Matrix.vandermonde columns := by
    ext i j
    rfl
  have rightPrincipal_pos :
      0 < ((expTaylorRight (terms := terms) columns).submatrix principal id).det := by
    calc
      0 < (((expTaylorRight (terms := terms) columns).submatrix principal id)ᵀ).det := by
        rw [rightTranspose]
        exact vandermonde_det_pos_of_strictMono columns columns_mono
      _ = ((expTaylorRight (terms := terms) columns).submatrix principal id).det :=
        Matrix.det_transpose _
  rw [expTaylorMatrix, Matrix.det_mul_rect]
  apply Finset.sum_pos'
  · intro q _
    exact mul_nonneg
      (expTaylorLeftMinor_det_nonneg rows q.1 rows_nonneg rows_mono q.2)
      (expTaylorRightMinor_det_nonneg columns q.1 columns_nonneg columns_mono q.2)
  · exact ⟨⟨principal, principal_mono⟩, Finset.mem_univ _,
      mul_pos leftPrincipal_pos rightPrincipal_pos⟩

/--
The square principal block of the exponential series, using powers
`0, ..., n - 1` and the positive weights `1 / k!`.
-/
noncomputable def expTaylorPrincipalMatrix {n : ℕ}
    (rows columns : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  Matrix.vandermonde rows *
    (Matrix.diagonal (fun k : Fin n => ((k.1.factorial : ℝ)⁻¹)) *
      (Matrix.vandermonde columns)ᵀ)

/-- The principal block evaluates the first `n` terms of `exp(x*y)`. -/
theorem expTaylorPrincipalMatrix_apply {n : ℕ}
    (rows columns : Fin n → ℝ) (i j : Fin n) :
    expTaylorPrincipalMatrix rows columns i j =
      ∑ k : Fin n, rows i ^ k.1 * (k.1.factorial : ℝ)⁻¹ * columns j ^ k.1 := by
  rw [expTaylorPrincipalMatrix, Matrix.mul_apply]
  simp only [Matrix.diagonal_mul, Matrix.transpose_apply, Matrix.vandermonde_apply]
  congr 1
  funext k
  ring

/--
The principal exponential-series block has positive determinant on increasing
row and column tuples. This certifies one strictly positive Cauchy--Binet term;
it does not yet control the signs of the remaining generalized terms.
-/
theorem expTaylorPrincipalMatrix_det_pos {n : ℕ}
    (rows columns : Fin n → ℝ)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns) :
    0 < (expTaylorPrincipalMatrix rows columns).det := by
  have rows_pos := vandermonde_det_pos_of_strictMono rows rows_mono
  have columns_pos := vandermonde_det_pos_of_strictMono columns columns_mono
  have weights_pos :
      0 < (Matrix.diagonal
        (fun k : Fin n => ((k.1.factorial : ℝ)⁻¹))).det := by
    rw [Matrix.det_diagonal]
    exact Finset.prod_pos fun k _ =>
      inv_pos.mpr (Nat.cast_pos.mpr k.1.factorial_pos)
  simp only [expTaylorPrincipalMatrix, Matrix.det_mul, Matrix.det_transpose]
  positivity

end Fermionic
end GECPKernelStructure
