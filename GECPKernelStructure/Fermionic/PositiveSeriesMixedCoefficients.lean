import GECPKernelStructure.Matrix.MixedColumnCoefficients
import GECPKernelStructure.Fermionic.PositiveSeriesTail

namespace Matrix

open Filter Finset Topology

variable {n : Type*} [Fintype n] [DecidableEq n]

set_option maxHeartbeats 800000 in
theorem mixedColumnCoefficient_tendsto {ι : Type*} {l : Filter ι}
    (degree : ℕ) (E B : ι → Matrix n n ℝ) (E₀ B₀ : Matrix n n ℝ)
    (hE : Tendsto E l (𝓝 E₀)) (hB : Tendsto B l (𝓝 B₀)) :
    Tendsto (fun k => mixedColumnCoefficient degree (E k) (B k)) l
      (𝓝 (mixedColumnCoefficient degree E₀ B₀)) := by
  have continuous_coefficient : Continuous
      (fun pair : Matrix n n ℝ × Matrix n n ℝ =>
        mixedColumnCoefficient degree pair.1 pair.2) := by
    unfold mixedColumnCoefficient
    apply continuous_finsetSum
    intro s _
    by_cases hcard : s.card = degree
    · simp only [if_pos hcard]
      apply Continuous.matrix_det
      apply continuous_pi
      intro i
      apply continuous_pi
      intro j
      by_cases hj : j ∈ s
      · simp [columnChoice, hj]
        fun_prop
      · simp [columnChoice, hj]
        fun_prop
    · simp [hcard]
      fun_prop
  have hp : Tendsto (fun k => (E k, B k)) l (𝓝 (E₀, B₀)) := hE.prodMk_nhds hB
  exact continuous_coefficient.continuousAt.tendsto.comp hp

theorem mixedColumnCoefficient_nonneg_of_tendsto {ι : Type*} {l : Filter ι}
    [NeBot l] (degree : ℕ) (E B : ι → Matrix n n ℝ)
    (E₀ B₀ : Matrix n n ℝ)
    (hE : Tendsto E l (𝓝 E₀)) (hB : Tendsto B l (𝓝 B₀))
    (hnonneg : ∀ᶠ k in l, 0 ≤ mixedColumnCoefficient degree (E k) (B k)) :
    0 ≤ mixedColumnCoefficient degree E₀ B₀ := by
  exact ge_of_tendsto
    (mixedColumnCoefficient_tendsto degree E B E₀ B₀ hE hB) hnonneg

end Matrix

namespace GECPKernelStructure
namespace Fermionic

open Filter Matrix Topology

/-- The first eight positive exponential Taylor features as a square sample matrix. -/
noncomputable def positiveEighthOrderTaylorMatrix {order : ℕ}
    (rows columns : Fin order → ℝ) : Matrix (Fin order) (Fin order) ℝ :=
  expTaylorMatrix (terms := 8) rows columns

/-- A finite truncation of the positive exponential tail beyond degree seven. -/
noncomputable def positiveEighthOrderFiniteTailMatrix {order : ℕ}
    (terms : ℕ) (rows columns : Fin order → ℝ) :
    Matrix (Fin order) (Fin order) ℝ :=
  expTaylorMatrix (terms := terms + 8) rows columns -
    positiveEighthOrderTaylorMatrix rows columns

/-- The exact positive exponential eighth-order tail sample matrix. -/
noncomputable def positiveEighthOrderTailMatrix {order : ℕ}
    (rows columns : Fin order → ℝ) : Matrix (Fin order) (Fin order) ℝ :=
  expMatrix rows columns - positiveEighthOrderTaylorMatrix rows columns

theorem positiveEighthOrderTailMatrix_apply {order : ℕ}
    (rows columns : Fin order → ℝ) (i j : Fin order) :
    positiveEighthOrderTailMatrix rows columns i j =
      positiveEighthOrderTail (rows i * columns j) := by
  simp only [positiveEighthOrderTailMatrix, Matrix.sub_apply, expMatrix,
    positiveEighthOrderTaylorMatrix, expTaylorMatrix_apply,
    positiveEighthOrderTail, eighthOrderTail]
  rw [show -(-(rows i * columns j)) = rows i * columns j by ring]
  congr 1
  rw [expNegTaylor, ← Fin.sum_univ_eq_sum_range]
  apply Finset.sum_congr rfl
  intro k hk
  field_simp
  ring

theorem positiveEighthOrderFiniteTailMatrix_tendsto {order : ℕ}
    (rows columns : Fin order → ℝ) :
    Tendsto (fun terms => positiveEighthOrderFiniteTailMatrix terms rows columns)
      atTop (𝓝 (positiveEighthOrderTailMatrix rows columns)) := by
  exact ((expTaylorMatrix_tendsto_expMatrix rows columns).comp
    (tendsto_add_atTop_nat 8)).sub_const
      (positiveEighthOrderTaylorMatrix rows columns)

theorem lowFeatureProduct_eight_apply {order terms : ℕ}
    (rows columns : Fin order → ℝ) (i j : Fin order) :
    (lowFeatureLeft (fun k : Fin (terms + 8) => 8 ≤ k.1)
      (expTaylorLeft (terms := terms + 8) rows) *
        expTaylorRight (terms := terms + 8) columns) i j =
      positiveEighthOrderTaylorMatrix rows columns i j := by
  rw [Matrix.mul_apply, positiveEighthOrderTaylorMatrix,
    expTaylorMatrix_apply]
  simp only [lowFeatureLeft, expTaylorLeft, expTaylorRight, ite_mul, zero_mul]
  let f : ℕ → ℝ := fun k =>
    rows i ^ k * (k.factorial : ℝ)⁻¹ * columns j ^ k
  change (∑ k : Fin (terms + 8), if 8 ≤ k.1 then 0 else f k.1) =
    ∑ k : Fin 8, f k.1
  rw [Fin.sum_univ_eq_sum_range (fun k => if 8 ≤ k then 0 else f k),
    Fin.sum_univ_eq_sum_range f]
  change (∑ k ∈ Finset.range (terms + 8), if 8 ≤ k then 0 else f k) =
    ∑ k ∈ Finset.range 8, f k
  calc
    (∑ k ∈ Finset.range (terms + 8), if 8 ≤ k then 0 else f k) =
        ∑ k ∈ (Finset.range (terms + 8)).filter (fun k => k < 8), f k := by
      rw [Finset.sum_filter]
      apply Finset.sum_congr rfl
      intro k hk
      by_cases h : k < 8
      · simp [h, Nat.not_le.mpr h]
      · simp [h, Nat.le_of_not_gt h]
    _ = ∑ k ∈ Finset.range 8, f k := by
      congr 1
      ext k
      simp
      omega

theorem lowFeatureProduct_eight {order terms : ℕ}
    (rows columns : Fin order → ℝ) :
    lowFeatureLeft (fun k : Fin (terms + 8) => 8 ≤ k.1)
      (expTaylorLeft (terms := terms + 8) rows) *
        expTaylorRight (terms := terms + 8) columns =
      positiveEighthOrderTaylorMatrix rows columns := by
  ext i j
  exact lowFeatureProduct_eight_apply rows columns i j

theorem highFeatureProduct_eight {order terms : ℕ}
    (rows columns : Fin order → ℝ) :
    highFeatureLeft (fun k : Fin (terms + 8) => 8 ≤ k.1)
      (expTaylorLeft (terms := terms + 8) rows) *
        expTaylorRight (terms := terms + 8) columns =
      positiveEighthOrderFiniteTailMatrix terms rows columns := by
  have split := weightedFeatureProduct_eq_low_add_smul_high
    (fun k : Fin (terms + 8) => 8 ≤ k.1) 1
    (expTaylorLeft (terms := terms + 8) rows)
    (expTaylorRight (terms := terms + 8) columns)
  have weighted_eq :
      weightedFeatureLeft (fun k : Fin (terms + 8) => 8 ≤ k.1) 1
          (expTaylorLeft (terms := terms + 8) rows) *
        expTaylorRight (terms := terms + 8) columns =
          expTaylorMatrix (terms := terms + 8) rows columns := by
    have left_eq :
        weightedFeatureLeft (fun k : Fin (terms + 8) => 8 ≤ k.1) 1
            (expTaylorLeft (terms := terms + 8) rows) =
          expTaylorLeft (terms := terms + 8) rows := by
      ext i k
      simp [weightedFeatureLeft, featureBlockWeight]
    rw [left_eq]
    rfl
  rw [weighted_eq, one_smul, lowFeatureProduct_eight] at split
  rw [positiveEighthOrderFiniteTailMatrix]
  exact eq_sub_of_add_eq (by simpa [add_comm] using split.symm)

theorem positiveEighthOrderFiniteTail_mixedColumnCoefficient_nonneg
    {order : ℕ} (terms : ℕ) (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (degree : ℕ) :
    0 ≤ mixedColumnCoefficient degree
      (positiveEighthOrderFiniteTailMatrix terms rows columns)
      (positiveEighthOrderTaylorMatrix rows columns) := by
  simpa [highFeatureProduct_eight, lowFeatureProduct_eight] using
    (expTaylor_mixedColumnCoefficient_nonneg (terms := terms + 8) 8
      rows columns rows_nonneg columns_nonneg rows_mono columns_mono degree)

theorem positiveEighthOrderTail_mixedColumnCoefficient_nonneg
    {order : ℕ} (rows columns : Fin order → ℝ)
    (rows_nonneg : ∀ i, 0 ≤ rows i) (columns_nonneg : ∀ i, 0 ≤ columns i)
    (rows_mono : StrictMono rows) (columns_mono : StrictMono columns)
    (degree : ℕ) :
    0 ≤ mixedColumnCoefficient degree
      (positiveEighthOrderTailMatrix rows columns)
      (positiveEighthOrderTaylorMatrix rows columns) := by
  apply mixedColumnCoefficient_nonneg_of_tendsto degree
    (fun terms => positiveEighthOrderFiniteTailMatrix terms rows columns)
    (fun _ => positiveEighthOrderTaylorMatrix rows columns)
    (positiveEighthOrderTailMatrix rows columns)
    (positiveEighthOrderTaylorMatrix rows columns)
    (positiveEighthOrderFiniteTailMatrix_tendsto rows columns)
    tendsto_const_nhds
  exact Filter.Eventually.of_forall fun terms =>
    positiveEighthOrderFiniteTail_mixedColumnCoefficient_nonneg terms
      rows columns rows_nonneg columns_nonneg rows_mono columns_mono degree

end Fermionic
end GECPKernelStructure
