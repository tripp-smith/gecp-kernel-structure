import GECPKernelStructure.Fermionic.Symmetry
import GECPKernelStructure.GECP.Definitions
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

open GECP

universe u v

/-- The only classes of rigorous outcome allowed to close structural Milestone E. -/
inductive GECPResearchOutcome where
  | targetRate
  | weakerImprovedRate
  | containingStructuralClass
  | certifiedObstruction
  deriving DecidableEq, Repr

/--
A selected cross controls the next residual relative to a uniform bound `B`.
This is the quantitative refinement that survives the pivot-order obstruction.
-/
def CrossRatioControl {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (row : α) (column : β) (θ B : ℝ) : Prop :=
  ∀ x y, abs (R x y * R row column - R x column * R row y) ≤
    θ * abs (R row column) * B

/-- The two products through a proposed pivot have compatible signs. -/
def PivotCrossProductSignCoherent {α : Type u} {β : Type v} (R : Kernel α β ℝ)
    (row : α) (column : β) : Prop :=
  ∀ x y,
    0 ≤ (R x y * R row column) * (R x column * R row y)

/-- Every cross of a residual has compatible product signs. -/
def CrossProductSignCoherent {α : Type u} {β : Type v} (R : Kernel α β ℝ) : Prop :=
  ∀ row column, PivotCrossProductSignCoherent R row column

private theorem abs_sub_le_of_sameSign_of_abs_le {a b B : ℝ}
    (same_sign : 0 ≤ a * b) (ha : abs a ≤ B) (hb : abs b ≤ B) :
    abs (a - b) ≤ B := by
  rw [abs_le]
  rcases (mul_nonneg_iff.mp same_sign) with ⟨ha0, hb0⟩ | ⟨ha0, hb0⟩
  · rw [abs_of_nonneg ha0] at ha
    rw [abs_of_nonneg hb0] at hb
    constructor <;> linarith
  · rw [abs_of_nonpos ha0] at ha
    rw [abs_of_nonpos hb0] at hb
    constructor <;> linarith

/--
Sign-coherent crosses and an exact complete pivot imply factor-one
cross-ratio control. Thus this structure rules out residual element growth,
although it does not by itself give strict contraction.
-/
theorem crossRatioControl_one_of_signCoherent
    {α : Type u} {β : Type v} {R : Kernel α β ℝ}
    {row : α} {column : β} (coherent : PivotCrossProductSignCoherent R row column)
    {B : ℝ}
    (bound : ∀ x y, abs (R x y) ≤ B) (pivot_max : abs (R row column) = B) :
    CrossRatioControl R row column 1 B := by
  intro x y
  have first_product : abs (R x y * R row column) ≤ B * B := by
    rw [abs_mul]
    exact mul_le_mul (bound x y) (bound row column) (abs_nonneg _) (by
      exact (abs_nonneg _).trans (bound x y))
  have second_product : abs (R x column * R row y) ≤ B * B := by
    rw [abs_mul]
    exact mul_le_mul (bound x column) (bound row y) (abs_nonneg _) (by
      exact (abs_nonneg _).trans (bound x column))
  have difference := abs_sub_le_of_sameSign_of_abs_le
    (coherent x y) first_product second_product
  rw [pivot_max, one_mul]
  exact difference

/-- Cross-ratio control gives the corresponding one-step residual contraction. -/
theorem residualUpdate_le_of_crossRatioControl {α : Type u} {β : Type v}
    {R : Kernel α β ℝ} {row : α} {column : β} (pivot_ne : R row column ≠ 0)
    {θ B : ℝ} (control : CrossRatioControl R row column θ B) (x : α) (y : β) :
    abs (residualUpdate R row column pivot_ne x y) ≤ θ * B := by
  unfold residualUpdate
  rw [show R x y - R x column * R row y / R row column =
      (R x y * R row column - R x column * R row y) / R row column by
    field_simp]
  rw [abs_div]
  apply (div_le_div_of_nonneg_right (control x y) (abs_nonneg _)).trans_eq
  field_simp [abs_ne_zero.mpr pivot_ne]

/-- Exact complete pivoting is nonexpansive on a sign-coherent residual. -/
theorem residualUpdate_le_of_signCoherent {α : Type u} {β : Type v}
    {R : Kernel α β ℝ} {row : α} {column : β}
    (coherent : PivotCrossProductSignCoherent R row column) {B : ℝ}
    (bound : ∀ x y, abs (R x y) ≤ B) (pivot_max : abs (R row column) = B)
    (pivot_ne : R row column ≠ 0) (x : α) (y : β) :
    abs (residualUpdate R row column pivot_ne x y) ≤ B := by
  simpa using residualUpdate_le_of_crossRatioControl pivot_ne
    (crossRatioControl_one_of_signCoherent coherent bound pivot_max) x y

/--
A quantitative cross-ratio hypothesis at every selected pivot is sufficient
for geometric GECP residual decay. The theorem deliberately does not assert
that the fermionic kernel satisfies the hypothesis.
-/
theorem gecp_error_le_geometric_of_crossRatioControl
    {α : Type u} {β : Type v} (residual : ℕ → Kernel α β ℝ)
    (rows : ℕ → α) (columns : ℕ → β) (θ initialBound : ℝ)
    (theta_nonneg : 0 ≤ θ) (initial_nonneg : 0 ≤ initialBound)
    (initial_bound : ∀ x y, abs (residual 0 x y) ≤ initialBound)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (cross_ratio : ∀ n B, 0 ≤ B →
      (∀ x y, abs (residual n x y) ≤ B) →
      CrossRatioControl (residual n) (rows n) (columns n) θ B) :
    ∀ n x y, abs (residual n x y) ≤ θ ^ n * initialBound := by
  intro n
  induction n with
  | zero => simpa using initial_bound
  | succ n ih =>
      intro x y
      rw [show n + 1 = Nat.succ n by omega, updates n]
      have bound_nonneg : 0 ≤ θ ^ n * initialBound :=
        mul_nonneg (pow_nonneg theta_nonneg n) initial_nonneg
      have step := residualUpdate_le_of_crossRatioControl (pivot_ne n)
        (cross_ratio n (θ ^ n * initialBound) bound_nonneg ih) x y
      calc
        abs (residualUpdate (residual n) (rows n) (columns n) (pivot_ne n) x y) ≤
            θ * (θ ^ n * initialBound) := step
        _ = θ ^ (n + 1) * initialBound := by rw [pow_succ]; ring

/--
Cross-ratio factors may vary from step to step. This is the form needed for
cutoff-dependent or blockwise contraction: individual factors need not be
bounded by one, provided their products over the relevant blocks contract.
-/
theorem gecp_error_le_product_of_crossRatioControl
    {α : Type u} {β : Type v} (residual : ℕ → Kernel α β ℝ)
    (rows : ℕ → α) (columns : ℕ → β) (θ : ℕ → ℝ) (initialBound : ℝ)
    (theta_nonneg : ∀ n, 0 ≤ θ n) (initial_nonneg : 0 ≤ initialBound)
    (initial_bound : ∀ x y, abs (residual 0 x y) ≤ initialBound)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (cross_ratio : ∀ n B, 0 ≤ B →
      (∀ x y, abs (residual n x y) ≤ B) →
      CrossRatioControl (residual n) (rows n) (columns n) (θ n) B) :
    ∀ n x y, abs (residual n x y) ≤
      (∏ i ∈ Finset.range n, θ i) * initialBound := by
  intro n
  induction n with
  | zero => simpa using initial_bound
  | succ n ih =>
      intro x y
      rw [updates n]
      have product_nonneg : 0 ≤ ∏ i ∈ Finset.range n, θ i :=
        Finset.prod_nonneg fun i hi ↦ theta_nonneg i
      have bound_nonneg : 0 ≤ (∏ i ∈ Finset.range n, θ i) * initialBound :=
        mul_nonneg product_nonneg initial_nonneg
      have step := residualUpdate_le_of_crossRatioControl (pivot_ne n)
        (cross_ratio n ((∏ i ∈ Finset.range n, θ i) * initialBound)
          bound_nonneg ih) x y
      calc
        abs (residualUpdate (residual n) (rows n) (columns n) (pivot_ne n) x y) ≤
            θ n * ((∏ i ∈ Finset.range n, θ i) * initialBound) := step
        _ = (∏ i ∈ Finset.range (Nat.succ n), θ i) * initialBound := by
          rw [Finset.prod_range_succ]
          ring

/--
A dyadically contracting product over cutoff-sized blocks gives the natural
`block * accuracyBits` rank scale. This isolates the exact quantitative lemma
still required for the fermionic kernel.
-/
theorem gecp_error_le_dyadic_of_crossRatioProduct
    {α : Type u} {β : Type v} (residual : ℕ → Kernel α β ℝ)
    (rows : ℕ → α) (columns : ℕ → β) (θ : ℕ → ℝ)
    (initialBound : ℝ) (block : ℕ)
    (theta_nonneg : ∀ n, 0 ≤ θ n) (initial_nonneg : 0 ≤ initialBound)
    (initial_bound : ∀ x y, abs (residual 0 x y) ≤ initialBound)
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (cross_ratio : ∀ n B, 0 ≤ B →
      (∀ x y, abs (residual n x y) ≤ B) →
      CrossRatioControl (residual n) (rows n) (columns n) (θ n) B)
    (block_product : ∀ p,
      (∏ i ∈ Finset.range (block * p), θ i) ≤ (1 / 2 : ℝ) ^ p) :
    ∀ p x y, abs (residual (block * p) x y) ≤
      (1 / 2 : ℝ) ^ p * initialBound := by
  intro p x y
  exact (gecp_error_le_product_of_crossRatioControl residual rows columns θ
    initialBound theta_nonneg initial_nonneg initial_bound pivot_ne updates
    cross_ratio (block * p) x y).trans
      (mul_le_mul_of_nonneg_right (block_product p) initial_nonneg)

/-- On the nonnegative-frequency half-domain, the positive cutoff corner is maximal. -/
theorem fermionicKernel_le_positiveCutoffCorner {Λ t ω : ℝ}
    (ht : 0 ≤ t) (hω : 0 ≤ ω) (hωΛ : ω ≤ Λ) :
    fermionicKernel t ω ≤ fermionicKernel 0 Λ := by
  have numerator_le : Real.exp (-t * ω) ≤ 1 :=
    Real.exp_le_one_iff.mpr (by nlinarith)
  have denominator_mono : 1 + Real.exp (-Λ) ≤ 1 + Real.exp (-ω) := by
    gcongr
  unfold fermionicKernel
  calc
    Real.exp (-t * ω) / (1 + Real.exp (-ω)) ≤
        1 / (1 + Real.exp (-ω)) := by gcongr
    _ ≤ 1 / (1 + Real.exp (-Λ)) := by gcongr
    _ = Real.exp (-0 * Λ) / (1 + Real.exp (-Λ)) := by norm_num

/-- Both reflected cutoff corners are complete pivots on the cutoff rectangle. -/
theorem fermionicKernel_le_cutoffCorner {Λ t ω : ℝ}
    (ht0 : 0 ≤ t) (ht1 : t ≤ 1) (hω0 : -Λ ≤ ω) (hω1 : ω ≤ Λ) :
    fermionicKernel t ω ≤ fermionicKernel 0 Λ := by
  by_cases hω : 0 ≤ ω
  · exact fermionicKernel_le_positiveCutoffCorner ht0 hω hω1
  · rw [← fermionicKernel_reflection]
    have reflected_time : 0 ≤ 1 - t := by linarith
    have reflected_frequency : 0 ≤ -ω := by linarith
    have reflected_cutoff : -ω ≤ Λ := by linarith
    exact fermionicKernel_le_positiveCutoffCorner reflected_time
      reflected_frequency reflected_cutoff

/-- The first complete pivot `(0, Λ)` leaves this exact reflected-corner residual. -/
theorem fermionicKernel_firstPivot_reflected_residual (Λ : ℝ) :
    residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
      1 (-Λ) = 1 - Real.exp (-Λ) := by
  unfold residualUpdate fermionicKernel
  simp only [zero_mul, neg_zero, Real.exp_zero]
  rw [Real.exp_neg]
  field_simp [Real.exp_ne_zero]
  have exp_mul : Real.exp (-Λ) * Real.exp Λ = 1 := by
    rw [← Real.exp_add]
    simp
  nlinarith [exp_mul]

/-- Closed form of the residual after selecting the positive cutoff corner. -/
theorem fermionicKernel_firstPivot_residual (Λ t ω : ℝ) :
    residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
      t ω = (Real.exp (-t * ω) - Real.exp (-t * Λ)) /
        (1 + Real.exp (-ω)) := by
  unfold residualUpdate
  have cross_term :
      fermionicKernel t Λ * fermionicKernel 0 ω / fermionicKernel 0 Λ =
        Real.exp (-t * Λ) / (1 + Real.exp (-ω)) := by
    unfold fermionicKernel
    simp only [zero_mul, neg_zero, Real.exp_zero]
    field_simp [fermionicKernel_denominator_ne]
  rw [cross_term]
  unfold fermionicKernel
  ring

/--
For every positive cutoff, the first residual is nonnegative on the cutoff
rectangle and its complete maximum is the reflected corner `(1, -Λ)`. Thus the
first two canonical fermionic GECP pivots are rigorously identified.
-/
theorem fermionicKernel_firstPivot_abs_le_reflectedCorner {Λ t ω : ℝ}
    (cutoff_pos : 0 < Λ)
    (time_nonneg : 0 ≤ t) (time_le_one : t ≤ 1)
    (frequency_lower : -Λ ≤ ω) (frequency_upper : ω ≤ Λ) :
    abs (residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' t ω) ≤
      residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) := by
  rw [fermionicKernel_firstPivot_residual,
    fermionicKernel_firstPivot_reflected_residual]
  have numerator_nonneg :
      0 ≤ Real.exp (-t * ω) - Real.exp (-t * Λ) := by
    exact sub_nonneg.mpr (Real.exp_le_exp.mpr (by nlinarith))
  rw [abs_of_nonneg (div_nonneg numerator_nonneg
    (fermionicKernel_denominator_pos ω).le)]
  have pivot_nonneg : 0 ≤ 1 - Real.exp (-Λ) := by
    exact sub_nonneg.mpr (Real.exp_le_one_iff.mpr (by linarith))
  by_cases frequency_nonneg : 0 ≤ ω
  · have exp_time_frequency_le_one : Real.exp (-t * ω) ≤ 1 :=
      Real.exp_le_one_iff.mpr (by nlinarith)
    have exp_cutoff_le_time : Real.exp (-Λ) ≤ Real.exp (-t * Λ) :=
      Real.exp_le_exp.mpr (by nlinarith)
    have numerator_le_pivot :
        Real.exp (-t * ω) - Real.exp (-t * Λ) ≤ 1 - Real.exp (-Λ) := by
      linarith
    apply (div_le_iff₀ (fermionicKernel_denominator_pos ω)).mpr
    nlinarith [Real.exp_pos (-ω)]
  · have frequency_nonpos : ω ≤ 0 := le_of_not_ge frequency_nonneg
    have first_exp_le : Real.exp (-t * ω) ≤ Real.exp (-ω) :=
      Real.exp_le_exp.mpr (by nlinarith)
    have second_exp_le : Real.exp (-Λ) ≤ Real.exp (-t * Λ) :=
      Real.exp_le_exp.mpr (by nlinarith)
    have numerator_endpoint :
        Real.exp (-t * ω) - Real.exp (-t * Λ) ≤
          Real.exp (-ω) - Real.exp (-Λ) := by linarith
    have first_bound :
        (Real.exp (-t * ω) - Real.exp (-t * Λ)) / (1 + Real.exp (-ω)) ≤
          (Real.exp (-ω) - Real.exp (-Λ)) / (1 + Real.exp (-ω)) := by
      gcongr
    have frequency_exp_le : Real.exp (-ω) ≤ Real.exp Λ :=
      Real.exp_le_exp.mpr (by linarith)
    have second_bound :
        (Real.exp (-ω) - Real.exp (-Λ)) / (1 + Real.exp (-ω)) ≤
          (Real.exp Λ - Real.exp (-Λ)) / (1 + Real.exp Λ) := by
      have denominator_pos := fermionicKernel_denominator_pos ω
      have cutoff_denominator_pos : 0 < 1 + Real.exp Λ := by positivity
      apply (div_le_div_iff₀ denominator_pos cutoff_denominator_pos).mpr
      nlinarith [Real.exp_pos (-Λ), Real.exp_pos Λ]
    calc
      (Real.exp (-t * ω) - Real.exp (-t * Λ)) / (1 + Real.exp (-ω)) ≤
          (Real.exp (-ω) - Real.exp (-Λ)) / (1 + Real.exp (-ω)) := first_bound
      _ ≤ (Real.exp Λ - Real.exp (-Λ)) / (1 + Real.exp Λ) := second_bound
      _ = 1 - Real.exp (-Λ) := by
        rw [Real.exp_neg]
        field_simp [Real.exp_ne_zero]
        ring

/-- The reflected-corner residual divided by the initial complete-pivot value. -/
theorem fermionicKernel_firstPivot_reflected_ratio (Λ : ℝ) :
    residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
      1 (-Λ) / fermionicKernel 0 Λ = 1 - Real.exp (-(2 * Λ)) := by
  rw [fermionicKernel_firstPivot_reflected_residual]
  unfold fermionicKernel
  rw [show -(2 * Λ) = -Λ + -Λ by ring, Real.exp_add]
  field_simp [fermionicKernel_denominator_ne]
  norm_num
  ring

/--
No contraction factor below one controls even the first reflected-corner
residual uniformly in the cutoff. The witness is explicit, not asymptotic.
-/
theorem fermionicKernel_no_uniform_firstStep_contraction
    {θ : ℝ} (theta_nonneg : 0 ≤ θ) (theta_lt_one : θ < 1) :
    ∃ Λ : ℝ, 0 < Λ ∧
      θ < residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne'
        1 (-Λ) / fermionicKernel 0 Λ := by
  let q : ℝ := (1 - θ) / 2
  have q_pos : 0 < q := by dsimp [q]; linarith
  have q_lt_one : q < 1 := by dsimp [q]; linarith
  refine ⟨-(Real.log q) / 2, ?_, ?_⟩
  · have log_q_neg : Real.log q < 0 := Real.log_neg q_pos q_lt_one
    linarith
  · rw [fermionicKernel_firstPivot_reflected_ratio]
    have exponent_identity : -(2 * (-(Real.log q) / 2)) = Real.log q := by ring
    rw [exponent_identity, Real.exp_log q_pos]
    dsimp [q]
    linarith

end Fermionic
end GECPKernelStructure
