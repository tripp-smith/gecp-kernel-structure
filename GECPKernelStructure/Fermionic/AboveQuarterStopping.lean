import GECPKernelStructure.Fermionic.CentralPrefixRun
import Mathlib.Analysis.Complex.ExponentialBounds

namespace GECPKernelStructure

namespace GECP
universe u v

namespace Run

/-- Reindexing the selected core does not change its absolute pivot product. -/
theorem abs_selectedCore_det_eq_prod_abs_pivots
    {α : Type u} {β : Type v} {K : Kernel α β ℝ} (run : Run K) :
    |run.selectedCore.det| = (run.pivots.map abs).prod := by
  simpa [finSelectedCore, Matrix.det_reindex_self] using
    run.abs_finSelectedCore_det_eq_prod_abs_pivots

end Run
end GECP

namespace Fermionic

open GECP

/-- The positive cutoff-corner pivot is strictly larger than one half. -/
theorem fermionicKernel_cutoffCorner_gt_half {Λ : ℝ} (cutoff_pos : 0 < Λ) :
    1 / 2 < fermionicKernel 0 Λ := by
  have denominator_pos : 0 < 1 + Real.exp (-Λ) := by positivity
  have exponential_lt_one : Real.exp (-Λ) < 1 :=
    Real.exp_lt_one_iff.mpr (by linarith)
  unfold fermionicKernel
  simp only [zero_mul, neg_zero, Real.exp_zero]
  rw [lt_div_iff₀ denominator_pos]
  nlinarith

/-- The reflected second corner pivot is larger than one half for cutoff at least one. -/
theorem fermionicKernel_firstPivot_reflected_gt_half {Λ : ℝ}
    (cutoff_ge_one : 1 ≤ Λ) :
    1 / 2 < residualUpdate fermionicKernel 0 Λ
      (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) := by
  rw [fermionicKernel_firstPivot_reflected_residual]
  have exponential_le : Real.exp (-Λ) ≤ Real.exp (-1) := by
    exact Real.exp_le_exp.mpr (by linarith)
  linarith [Real.exp_neg_one_lt_half]

private theorem quarter_pow_le_prod_abs_range (f : ℕ → ℝ) (length : ℕ)
    (above_quarter : ∀ j, j < length → 1 / 4 < abs (f j)) :
    (1 / 4 : ℝ) ^ length ≤
      ((List.range length).map fun j => abs (f j)).prod := by
  induction length with
  | zero => simp
  | succ length ih =>
      rw [List.range_succ, List.map_append, List.prod_append]
      simp only [List.map_singleton, List.prod_singleton, pow_succ]
      have prefix_bound := ih fun j hj => above_quarter j (by omega)
      apply mul_le_mul prefix_bound
      · exact (above_quarter length (by omega)).le
      · positivity
      · exact (by positivity : 0 ≤ (1 / 4 : ℝ) ^ length).trans prefix_bound

/--
An all-above-quarter continuation forces a strict quarter-power lower bound on
the selected-core determinant of the full prescribed-corner run.
-/
theorem fermionicKernel_aboveQuarterPrefix_sample_det_gt_quarter_pow
    {Λ : ℝ} (cutoff_ge_one : 1 ≤ Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 =
      symmetricTwoCornerResidual Λ (lt_of_lt_of_le zero_lt_one cutoff_ge_one))
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (length : ℕ)
    (above_quarter : ∀ j, j < length →
      1 / 4 < abs (residual j (rows j) (columns j))) :
    let cutoff_pos := lt_of_lt_of_le zero_lt_one cutoff_ge_one
    let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos
      residual rows columns initial pivot_ne updates length
    let run := symmetricTwoCornerRun cutoff_pos continuation
    (1 / 4 : ℝ) ^ (length + 1) < |run.selectedCore.det| := by
  dsimp only
  let cutoff_pos : 0 < Λ := lt_of_lt_of_le zero_lt_one cutoff_ge_one
  let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos
    residual rows columns initial pivot_ne updates length
  let run := symmetricTwoCornerRun cutoff_pos continuation
  have first_gt : 1 / 2 < abs (fermionicKernel 0 Λ) := by
    rw [abs_of_pos (fermionicKernel_pos 0 Λ)]
    exact fermionicKernel_cutoffCorner_gt_half cutoff_pos
  have second_pos : 0 < residualUpdate fermionicKernel 0 Λ
      (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) :=
    (fermionicKernel_firstPivot_reflected_gt_half cutoff_ge_one).trans'
      (by norm_num)
  have second_gt : 1 / 2 < abs (residualUpdate fermionicKernel 0 Λ
      (fermionicKernel_pos 0 Λ).ne' 1 (-Λ)) := by
    rw [abs_of_pos second_pos]
    exact fermionicKernel_firstPivot_reflected_gt_half cutoff_ge_one
  have corners_gt : (1 / 4 : ℝ) <
      abs (fermionicKernel 0 Λ) *
        abs (residualUpdate fermionicKernel 0 Λ
          (fermionicKernel_pos 0 Λ).ne' 1 (-Λ)) := by
    nlinarith [mul_pos (sub_pos.mpr first_gt) (sub_pos.mpr second_gt)]
  have continuation_bound : (1 / 4 : ℝ) ^ length ≤
      ((List.range length).map fun j =>
        abs (residual j (rows j) (columns j))).prod :=
    quarter_pow_le_prod_abs_range
      (fun j => residual j (rows j) (columns j)) length above_quarter
  have pivot_list : run.pivots =
      fermionicKernel 0 Λ ::
        residualUpdate fermionicKernel 0 Λ
          (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) ::
        (List.range length).map (fun j =>
          residual j (rows j) (columns j)) := by
    change
      fermionicKernel 0 Λ ::
          residualUpdate fermionicKernel 0 Λ
            (fermionicKernel_pos 0 Λ).ne' 1 (-Λ) ::
          (Run.transport initial
            (Run.ofResidualSequenceFrom residual rows columns pivot_ne updates
              0 length)).pivots = _
    rw [Run.pivots_transport, Run.pivots_ofResidualSequenceFrom]
    simp [← List.range_eq_range']
  rw [run.abs_selectedCore_det_eq_prod_abs_pivots, pivot_list]
  simp only [List.map_cons, List.prod_cons]
  rw [List.map_map]
  calc
    (1 / 4 : ℝ) ^ (length + 1) =
        (1 / 4) * (1 / 4 : ℝ) ^ length := by rw [pow_succ]; ring
    _ < (abs (fermionicKernel 0 Λ) *
          abs (residualUpdate fermionicKernel 0 Λ
            (fermionicKernel_pos 0 Λ).ne' 1 (-Λ))) *
        (1 / 4 : ℝ) ^ length :=
      mul_lt_mul_of_pos_right corners_gt (by positivity)
    _ ≤ (abs (fermionicKernel 0 Λ) *
          abs (residualUpdate fermionicKernel 0 Λ
            (fermionicKernel_pos 0 Λ).ne' 1 (-Λ))) *
        ((List.range length).map fun j =>
          abs (residual j (rows j) (columns j))).prod := by
      gcongr
    _ = abs (fermionicKernel 0 Λ) *
          (abs (residualUpdate fermionicKernel 0 Λ
            (fermionicKernel_pos 0 Λ).ne' 1 (-Λ)) *
          ((List.range length).map fun j =>
            abs (residual j (rows j) (columns j))).prod) := by
      ring

/--
A binary size bound and an exponent budget absorb the determinant prefactor
into the quarter-power threshold used by the stopping contradiction.
-/
theorem aboveQuarter_determinant_upper_le_quarter_pow
    (length p s q : ℕ)
    (size_le : length + 2 ≤ 2 ^ q)
    (exponent_budget :
      (q + 3) * (length + 2) + 2 * ((s + 1) * (8 * p)) ≤
        p * (length - 2 * ((s + 1) * (8 * p)))) :
    (2 : ℝ) ^ (length + 2) * (length + 2).factorial *
          (1 / 2 : ℝ) ^
            (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) ≤
      (1 / 4 : ℝ) ^ (length + 1) := by
  let n := length + 2
  let r := 2 * ((s + 1) * (8 * p))
  have factorial_bound : ((length + 2).factorial : ℝ) ≤
      (length + 2 : ℝ) ^ (length + 2) := by
    exact_mod_cast (length + 2).factorial_le_pow
  have size_real : (length + 2 : ℝ) ≤ (2 : ℝ) ^ q := by
    exact_mod_cast size_le
  have size_power : (length + 2 : ℝ) ^ (length + 2) ≤
      ((2 : ℝ) ^ q) ^ (length + 2) :=
    pow_le_pow_left₀ (by positivity) size_real _
  have half_power :
      (1 / 2 : ℝ) ^
          (p * (length - 2 * ((s + 1) * (8 * p)))) ≤
        (1 / 2 : ℝ) ^
          ((q + 3) * (length + 2) + 2 * ((s + 1) * (8 * p))) :=
    pow_le_pow_of_le_one (by norm_num) (by norm_num) exponent_budget
  calc
    (2 : ℝ) ^ (length + 2) * (length + 2).factorial *
          (1 / 2 : ℝ) ^
            (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) ≤
      (2 : ℝ) ^ (length + 2) * (length + 2 : ℝ) ^ (length + 2) *
          (1 / 2 : ℝ) ^
            (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by gcongr
    _ ≤ (2 : ℝ) ^ (length + 2) * ((2 : ℝ) ^ q) ^ (length + 2) *
          (1 / 2 : ℝ) ^
            (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by gcongr
    _ ≤ (2 : ℝ) ^ (length + 2) * ((2 : ℝ) ^ q) ^ (length + 2) *
          (1 / 2 : ℝ) ^
            ((q + 3) * (length + 2) + 2 * ((s + 1) * (8 * p))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by gcongr
    _ = (1 / 4 : ℝ) ^ (length + 1) := by
      change (2 : ℝ) ^ n * ((2 : ℝ) ^ q) ^ n *
          (1 / 2 : ℝ) ^ ((q + 3) * n + r) * 2 ^ (r + 2) =
        (1 / 4 : ℝ) ^ (length + 1)
      have exponent_eq : (q + 3) * n + r =
          (n + q * n + r + 2) + 2 * (length + 1) := by
        dsimp [n]
        ring_nf
      rw [exponent_eq, pow_add (1 / 2 : ℝ)]
      have powers_of_two :
          (2 : ℝ) ^ n * ((2 : ℝ) ^ q) ^ n * 2 ^ (r + 2) =
            2 ^ (n + q * n + r + 2) := by
        have nested : ((2 : ℝ) ^ q) ^ n = 2 ^ (q * n) :=
          (pow_mul (2 : ℝ) q n).symm
        have combine_first : (2 : ℝ) ^ n * 2 ^ (q * n) =
            2 ^ (n + q * n) := (pow_add (2 : ℝ) n (q * n)).symm
        have combine_second : (2 : ℝ) ^ (n + q * n) * 2 ^ (r + 2) =
            2 ^ ((n + q * n) + (r + 2)) :=
          (pow_add (2 : ℝ) (n + q * n) (r + 2)).symm
        calc
          (2 : ℝ) ^ n * ((2 : ℝ) ^ q) ^ n * 2 ^ (r + 2) =
              2 ^ n * 2 ^ (q * n) * 2 ^ (r + 2) := by rw [nested]
          _ = 2 ^ (n + q * n) * 2 ^ (r + 2) := by rw [combine_first]
          _ = 2 ^ ((n + q * n) + (r + 2)) := combine_second
          _ = 2 ^ (n + q * n + r + 2) := by congr 1
      calc
        (2 : ℝ) ^ n * ((2 : ℝ) ^ q) ^ n *
              ((1 / 2 : ℝ) ^ (n + q * n + r + 2) *
                (1 / 2 : ℝ) ^ (2 * (length + 1))) *
            2 ^ (r + 2) =
          ((2 : ℝ) ^ n * ((2 : ℝ) ^ q) ^ n * 2 ^ (r + 2)) *
            (1 / 2 : ℝ) ^ (n + q * n + r + 2) *
            (1 / 2 : ℝ) ^ (2 * (length + 1)) := by ring
        _ = (2 : ℝ) ^ (n + q * n + r + 2) *
            (1 / 2 : ℝ) ^ (n + q * n + r + 2) *
            (1 / 2 : ℝ) ^ (2 * (length + 1)) := by rw [powers_of_two]
        _ = (1 / 2 : ℝ) ^ (2 * (length + 1)) := by
          rw [← mul_pow]
          norm_num
        _ = (1 / 4 : ℝ) ^ (length + 1) := by
          rw [pow_mul]
          norm_num

/--
Any parameter choice satisfying the determinant exponent budget forces a
continuation pivot at or below one quarter.
-/
theorem fermionicKernel_exists_pivot_le_quarter_of_exponent_budget
    {Λ : ℝ} (cutoff_ge_one : 1 ≤ Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 =
      symmetricTwoCornerResidual Λ (lt_of_lt_of_le zero_lt_one cutoff_ge_one))
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (length p s q : ℕ) (hp : 0 < p)
    (half_cutoff : Λ / 2 = ((2 ^ s : ℕ) : ℝ))
    (rank_le : 2 * ((s + 1) * (8 * p)) ≤ length)
    (size_le : length + 2 ≤ 2 ^ q)
    (exponent_budget :
      (q + 3) * (length + 2) + 2 * ((s + 1) * (8 * p)) ≤
        p * (length - 2 * ((s + 1) * (8 * p)))) :
    ∃ j, j < length ∧
      abs (residual j (rows j) (columns j)) ≤ 1 / 4 := by
  let cutoff_pos : 0 < Λ := lt_of_lt_of_le zero_lt_one cutoff_ge_one
  by_contra no_small_pivot
  push Not at no_small_pivot
  have above_quarter : ∀ j, j < length →
      1 / 4 < abs (residual j (rows j) (columns j)) := by
    intro j hj
    exact no_small_pivot j hj
  let continuation := symmetricTwoCornerContinuationPrefix cutoff_pos
    residual rows columns initial pivot_ne updates length
  let run := symmetricTwoCornerRun cutoff_pos continuation
  have determinant_lower : (1 / 4 : ℝ) ^ (length + 1) <
      |run.selectedCore.det| := by
    simpa [cutoff_pos, continuation, run] using
      fermionicKernel_aboveQuarterPrefix_sample_det_gt_quarter_pow
        cutoff_ge_one residual rows columns initial pivot_ne updates length
        above_quarter
  have determinant_upper : |run.selectedCore.det| ≤
      (2 : ℝ) ^ (length + 2) * (length + 2).factorial *
          (1 / 2 : ℝ) ^
            (p * (length - 2 * ((s + 1) * (8 * p)))) *
        2 ^ (2 * ((s + 1) * (8 * p)) + 2) := by
    simpa [cutoff_pos, continuation, run] using
      fermionicKernel_aboveQuarterPrefix_sample_det_le_two_pow cutoff_pos
        residual rows columns initial pivot_ne updates complete length p s hp
        half_cutoff above_quarter rank_le
  have upper_le_quarter := aboveQuarter_determinant_upper_le_quarter_pow
    length p s q size_le exponent_budget
  exact (not_le_of_gt determinant_lower)
    (determinant_upper.trans upper_le_quarter)

/-- Conservative separated-approximation order for the explicit stopping block. -/
def aboveQuarterStoppingOrder (s : ℕ) : ℕ := 24 * (s + 1)

/-- Conservative continuation length for the explicit stopping block. -/
def aboveQuarterStoppingLength (s : ℕ) : ℕ := 2048 * (s + 1) ^ 2

/-- Binary exponent controlling the factorial through `n! <= n^n`. -/
def aboveQuarterStoppingBits (s : ℕ) : ℕ := 2 * (s + 1) + 10

theorem aboveQuarterStoppingOrder_pos (s : ℕ) :
    0 < aboveQuarterStoppingOrder s := by
  simp [aboveQuarterStoppingOrder]

theorem aboveQuarterStopping_rank_le (s : ℕ) :
    2 * ((s + 1) * (8 * aboveQuarterStoppingOrder s)) ≤
      aboveQuarterStoppingLength s := by
  unfold aboveQuarterStoppingOrder aboveQuarterStoppingLength
  nlinarith

theorem aboveQuarterStopping_size_le (s : ℕ) :
    aboveQuarterStoppingLength s + 2 ≤ 2 ^ aboveQuarterStoppingBits s := by
  let scale := s + 1
  have base := Nat.two_mul_sq_add_one_le_two_pow_two_mul scale
  calc
    aboveQuarterStoppingLength s + 2 = 2048 * scale ^ 2 + 2 := by
      simp [aboveQuarterStoppingLength, scale]
    _ ≤ 1024 * (2 * scale ^ 2 + 1) := by omega
    _ ≤ 1024 * 2 ^ (2 * scale) := Nat.mul_le_mul_left 1024 base
    _ = 2 ^ (2 * scale + 10) := by rw [pow_add]; norm_num; ring
    _ = 2 ^ aboveQuarterStoppingBits s := by
      simp [aboveQuarterStoppingBits, scale]

theorem aboveQuarterStopping_exponent_budget (s : ℕ) :
    (aboveQuarterStoppingBits s + 3) *
          (aboveQuarterStoppingLength s + 2) +
        2 * ((s + 1) * (8 * aboveQuarterStoppingOrder s)) ≤
      aboveQuarterStoppingOrder s *
        (aboveQuarterStoppingLength s -
          2 * ((s + 1) * (8 * aboveQuarterStoppingOrder s))) := by
  let scale := s + 1
  have scale_pos : 0 < scale := by omega
  have rank_eq :
      2 * (scale * (8 * (24 * scale))) = 384 * scale ^ 2 := by ring
  have diff_eq :
      2048 * scale ^ 2 - 384 * scale ^ 2 = 1664 * scale ^ 2 := by omega
  change ((2 * scale + 10) + 3) * (2048 * scale ^ 2 + 2) +
      2 * (scale * (8 * (24 * scale))) ≤
    (24 * scale) *
      (2048 * scale ^ 2 - 2 * (scale * (8 * (24 * scale))))
  rw [rank_eq, diff_eq]
  have square_le_cube : scale ^ 2 ≤ scale ^ 3 := by nlinarith
  have linear_le_cube : scale ≤ scale ^ 3 := by nlinarith
  have one_le_cube : 1 ≤ scale ^ 3 := by nlinarith
  ring_nf
  omega

/--
At a dyadic half-cutoff, some exact complete continuation pivot is at most one
quarter within the explicit quadratic stopping block.
-/
theorem fermionicKernel_exists_pivot_le_quarter_explicit
    {Λ : ℝ} (cutoff_ge_one : 1 ≤ Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 =
      symmetricTwoCornerResidual Λ (lt_of_lt_of_le zero_lt_one cutoff_ge_one))
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (s : ℕ) (half_cutoff : Λ / 2 = ((2 ^ s : ℕ) : ℝ)) :
    ∃ j, j < aboveQuarterStoppingLength s ∧
      abs (residual j (rows j) (columns j)) ≤ 1 / 4 := by
  exact fermionicKernel_exists_pivot_le_quarter_of_exponent_budget
    cutoff_ge_one residual rows columns initial pivot_ne updates complete
    (aboveQuarterStoppingLength s) (aboveQuarterStoppingOrder s)
    s (aboveQuarterStoppingBits s) (aboveQuarterStoppingOrder_pos s)
    half_cutoff (aboveQuarterStopping_rank_le s)
    (aboveQuarterStopping_size_le s) (aboveQuarterStopping_exponent_budget s)

/-- The stopping pivot certifies a global quarter bound for its residual. -/
theorem fermionicKernel_exists_residual_le_quarter_explicit
    {Λ : ℝ} (cutoff_ge_one : 1 ≤ Λ)
    (residual : ℕ → Kernel ℝ ℝ ℝ) (rows columns : ℕ → ℝ)
    (initial : residual 0 =
      symmetricTwoCornerResidual Λ (lt_of_lt_of_le zero_lt_one cutoff_ge_one))
    (pivot_ne : ∀ n, residual n (rows n) (columns n) ≠ 0)
    (updates : ∀ n, residual (n + 1) =
      residualUpdate (residual n) (rows n) (columns n) (pivot_ne n))
    (complete : ∀ n,
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun ω : ℝ => -Λ ≤ ω ∧ ω ≤ Λ)
        (residual n) (rows n) (columns n))
    (s : ℕ) (half_cutoff : Λ / 2 = ((2 ^ s : ℕ) : ℝ)) :
    ∃ j, j < aboveQuarterStoppingLength s ∧
      ∀ t ω, (0 ≤ t ∧ t ≤ 1) → (-Λ ≤ ω ∧ ω ≤ Λ) →
        abs (residual j t ω) ≤ 1 / 4 := by
  obtain ⟨j, index_lt, pivot_le⟩ :=
    fermionicKernel_exists_pivot_le_quarter_explicit cutoff_ge_one residual
      rows columns initial pivot_ne updates complete s half_cutoff
  refine ⟨j, index_lt, ?_⟩
  intro t ω time_mem frequency_mem
  exact ((complete j).2.2 t ω time_mem frequency_mem).trans pivot_le

end Fermionic
end GECPKernelStructure
