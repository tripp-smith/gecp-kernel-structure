import GECPKernelStructure.Fermionic.SeparatedApprox
import GECPKernelStructure.Matrix.LowRankDeterminant

namespace GECPKernelStructure
namespace Fermionic

open Matrix

variable {n : Type*}

/-- Sample a bivariate function on prescribed row and column coordinates. -/
noncomputable def sampleMatrix (f : ℝ → ℝ → ℝ) (t ω : n → ℝ) : Matrix n n ℝ :=
  fun i j => f (t i) (ω j)

/-- Row factor of a finite separated representation. -/
noncomputable def separatedLeft (terms : List SeparatedTerm) (t : n → ℝ) :
    Matrix n (Fin terms.length) ℝ :=
  fun i a => (terms.get a).1 (t i)

/-- Column factor of a finite separated representation. -/
noncomputable def separatedRight (terms : List SeparatedTerm) (ω : n → ℝ) :
    Matrix (Fin terms.length) n ℝ :=
  fun a j => (terms.get a).2 (ω j)

/-- Matrix multiplication realizes evaluation of a separated representation. -/
theorem separatedLeft_mul_separatedRight (terms : List SeparatedTerm) (t ω : n → ℝ) :
    separatedLeft terms t * separatedRight terms ω =
      sampleMatrix (evalSeparated terms) t ω := by
  ext i j
  rw [Matrix.mul_apply]
  simp only [separatedLeft, separatedRight, sampleMatrix]
  rw [evalSeparated, ← Fin.sum_univ_fun_getElem]
  simp

variable [Fintype n] [DecidableEq n]

/--
Every fermionic sample determinant is controlled by the explicit separated
approximation terms that contain enough error columns to exceed its rank.
-/
theorem fermionicKernel_sample_det_le_separatedApprox (p s : ℕ) (hp : 0 < p)
    (t ω : n → ℝ) (ht0 : ∀ i, 0 ≤ t i) (ht1 : ∀ i, t i ≤ 1)
    (hω_lower : ∀ j, -((2 ^ s : ℕ) : ℝ) ≤ ω j)
    (hω_upper : ∀ j, ω j ≤ (2 ^ s : ℕ)) :
    |(sampleMatrix fermionicKernel t ω).det| ≤
      ∑ S : Finset n,
        if Fintype.card n ≤ 2 * ((s + 1) * (8 * p)) + S.card then
          (Fintype.card n).factorial * (1 / 2 : ℝ) ^ (p * S.card) *
            2 ^ (Fintype.card n - S.card)
        else 0 := by
  let terms := fermionicSeparatedTerms p s
  let L := separatedLeft terms t
  let R := separatedRight terms ω
  have hLR : L * R = sampleMatrix (evalSeparated terms) t ω :=
    separatedLeft_mul_separatedRight terms t ω
  have happrox : ∀ i j,
      |sampleMatrix fermionicKernel t ω i j - (L * R) i j| ≤ (1 / 2 : ℝ) ^ p := by
    intro i j
    rw [hLR]
    exact fermionicKernel_separatedApprox_error p s hp (ht0 i) (ht1 i)
      (hω_lower j) (hω_upper j)
  have heps : (1 / 2 : ℝ) ^ p ≤ 1 := by
    exact pow_le_one₀ (by norm_num) (by norm_num)
  have hfactor : ∀ i j, |(L * R) i j| ≤ (2 : ℝ) := by
    intro i j
    rw [hLR]
    change |evalSeparated terms (t i) (ω j)| ≤ (2 : ℝ)
    have hkpos := fermionicKernel_pos (t i) (ω j)
    have hkle := fermionicKernel_le_one (t := t i) (ω := ω j) (ht0 i) (ht1 i)
    have herr := fermionicKernel_separatedApprox_error p s hp (ht0 i) (ht1 i)
      (hω_lower j) (hω_upper j)
    have htri : |evalSeparated terms (t i) (ω j)| ≤
        |evalSeparated terms (t i) (ω j) - fermionicKernel (t i) (ω j)| +
          |fermionicKernel (t i) (ω j)| := by
      calc
        |evalSeparated terms (t i) (ω j)| =
            |(evalSeparated terms (t i) (ω j) - fermionicKernel (t i) (ω j)) +
              fermionicKernel (t i) (ω j)| := by ring_nf
        _ ≤ _ := abs_add_le _ _
    rw [abs_sub_comm] at htri
    rw [abs_of_pos hkpos] at htri
    linarith
  have h := Matrix.abs_det_le_of_factors_approx
    (sampleMatrix fermionicKernel t ω) L R ((1 / 2 : ℝ) ^ p) 2 happrox hfactor
  simpa [terms, L, R, fermionicSeparatedTerms_length, pow_mul, mul_assoc] using h

/--
A closed determinant-decay bound. Once the sample size exceeds the separated
rank, every additional column contributes another factor `2⁻ᵖ`.
-/
theorem fermionicKernel_sample_det_le_two_pow (p s : ℕ) (hp : 0 < p)
    (t ω : n → ℝ) (ht0 : ∀ i, 0 ≤ t i) (ht1 : ∀ i, t i ≤ 1)
    (hω_lower : ∀ j, -((2 ^ s : ℕ) : ℝ) ≤ ω j)
    (hω_upper : ∀ j, ω j ≤ (2 ^ s : ℕ))
    (hcard : 2 * ((s + 1) * (8 * p)) ≤ Fintype.card n) :
    |(sampleMatrix fermionicKernel t ω).det| ≤
      (2 : ℝ) ^ Fintype.card n * (Fintype.card n).factorial *
        (1 / 2 : ℝ) ^ (p * (Fintype.card n - 2 * ((s + 1) * (8 * p)))) *
          2 ^ (2 * ((s + 1) * (8 * p))) := by
  refine (fermionicKernel_sample_det_le_separatedApprox p s hp t ω ht0 ht1
    hω_lower hω_upper).trans ?_
  calc
    (∑ S : Finset n,
        if Fintype.card n ≤ 2 * ((s + 1) * (8 * p)) + S.card then
          (Fintype.card n).factorial * (1 / 2 : ℝ) ^ (p * S.card) *
            2 ^ (Fintype.card n - S.card)
        else 0) ≤
      ∑ _S : Finset n, (Fintype.card n).factorial *
        (1 / 2 : ℝ) ^ (p * (Fintype.card n - 2 * ((s + 1) * (8 * p)))) *
          2 ^ (2 * ((s + 1) * (8 * p))) := by
      apply Finset.sum_le_sum
      intro S _
      by_cases hS : Fintype.card n ≤ 2 * ((s + 1) * (8 * p)) + S.card
      · rw [if_pos hS]
        have heps : (1 / 2 : ℝ) ^ (p * S.card) ≤
            (1 / 2 : ℝ) ^ (p * (Fintype.card n - 2 * ((s + 1) * (8 * p)))) :=
          pow_le_pow_of_le_one (by norm_num) (by norm_num)
            (Nat.mul_le_mul_left p (by omega))
        have htwo : (2 : ℝ) ^ (Fintype.card n - S.card) ≤
            2 ^ (2 * ((s + 1) * (8 * p))) :=
          pow_le_pow_right₀ (by norm_num) (by omega)
        gcongr
      · rw [if_neg hS]
        positivity
    _ = (2 : ℝ) ^ Fintype.card n * (Fintype.card n).factorial *
        (1 / 2 : ℝ) ^ (p * (Fintype.card n - 2 * ((s + 1) * (8 * p)))) *
          2 ^ (2 * ((s + 1) * (8 * p))) := by
      simp [Fintype.card_finset, nsmul_eq_mul, mul_assoc]

end Fermionic
end GECPKernelStructure
