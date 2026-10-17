import GECPKernelStructure.GECP.SignRegularity
import GECPKernelStructure.Matrix.NoncrossingTotalNonnegative
import GECPKernelStructure.Matrix.UpperTangentCauchyBinet
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Logic.Equiv.Prod
import Mathlib.Order.Lex
import Mathlib.Tactic

namespace Matrix

open GECPKernelStructure GECP

/-- Identify the two tangent layers with a layer bit paired with the feature
index. -/
def layerFeatureEquiv (b : ℕ) :
    (Fin b ⊕ Fin b) ≃ (Fin 2 × Fin b) :=
  ((Equiv.prodCongr finTwoEquiv (Equiv.refl (Fin b))).trans
    (Equiv.boolProdEquivSum (Fin b))).symm

@[simp]
theorem layerFeatureEquiv_inl {b : ℕ} (k : Fin b) :
    layerFeatureEquiv b (Sum.inl k) = (0, k) := by
  rfl

@[simp]
theorem layerFeatureEquiv_inr {b : ℕ} (k : Fin b) :
    layerFeatureEquiv b (Sum.inr k) = (1, k) := by
  rfl

/-- Change from the layer-major finite ordering used by `upperTangentMatrix`
to feature-major ordering, with the two copies of each feature adjacent. -/
def featureInterleaveEquiv (b : ℕ) :
    Fin (b + b) ≃ (Fin b ×ₗ Fin 2) :=
  finSumFinEquiv.symm |>.trans (layerFeatureEquiv b) |>.trans
    (Equiv.prodComm (Fin 2) (Fin b)) |>.trans toLex

/-- One weighted two-layer tangent cell. -/
def tangentCell {R : Type*} [Zero R] [One R] (weight : R) :
    Matrix (Fin 2) (Fin 2) R :=
  fun i j =>
    if i = 0 then
      if j = 0 then 1 else weight
    else if j = 0 then 0 else 1

@[simp]
theorem tangentCell_zero_zero {R : Type*} [Zero R] [One R] (weight : R) :
    tangentCell weight 0 0 = 1 := by
  simp [tangentCell]

@[simp]
theorem tangentCell_zero_one {R : Type*} [Zero R] [One R] (weight : R) :
    tangentCell weight 0 1 = weight := by
  simp [tangentCell]

@[simp]
theorem tangentCell_one_zero {R : Type*} [Zero R] [One R] (weight : R) :
    tangentCell weight 1 0 = 0 := by
  simp [tangentCell]

@[simp]
theorem tangentCell_one_one {R : Type*} [Zero R] [One R] (weight : R) :
    tangentCell weight 1 1 = 1 := by
  simp [tangentCell]

theorem tangentCell_entry_nonneg (weight : ℝ) (weight_nonneg : 0 ≤ weight)
    (i j : Fin 2) : 0 ≤ tangentCell weight i j := by
  fin_cases i <;> fin_cases j <;> simp [weight_nonneg]

@[simp]
theorem tangentCell_det {R : Type*} [CommRing R] (weight : R) :
    (tangentCell weight).det = 1 := by
  rw [Matrix.det_fin_two]
  simp

/-- Every ordered square minor of one weighted tangent cell is nonnegative. -/
theorem tangentCell_totallyNonnegative (weight : ℝ)
    (weight_nonneg : 0 ≤ weight) :
    TotallyNonnegative (fun i j => tangentCell weight i j) := by
  intro n rows columns rows_mono columns_mono
  have order_le : n ≤ 2 := by
    simpa using Fintype.card_le_of_injective rows rows_mono.injective
  interval_cases n
  · simp
  · rw [Matrix.det_fin_one]
    exact tangentCell_entry_nonneg weight weight_nonneg (rows 0) (columns 0)
  · rw [rows_mono.eq_id, columns_mono.eq_id]
    change 0 ≤ (tangentCell weight).det
    simp

/-- In layer/feature product coordinates, a diagonal tangent middle matrix is
exactly a block-diagonal collection of independent weighted cells. -/
theorem reindex_upperTangent_diagonal_eq_blockDiagonal {b : ℕ}
    (weight : Fin b → ℝ) :
    Matrix.reindex (layerFeatureEquiv b) (layerFeatureEquiv b)
        (upperTangentMatrix (1 : Matrix (Fin b) (Fin b) ℝ)
          (Matrix.diagonal weight)) =
      Matrix.blockDiagonal (fun k => tangentCell (weight k)) := by
  ext ⟨layer, k⟩ ⟨layer', k'⟩
  fin_cases layer <;> fin_cases layer' <;>
    simp [Matrix.reindex_apply, layerFeatureEquiv, upperTangentMatrix,
      tangentCell, Matrix.blockDiagonal_apply, finTwoEquiv,
      Equiv.boolProdEquivSum, Matrix.one_apply, Matrix.diagonal_apply]

/-- The feature-major version of the diagonal tangent middle matrix. -/
def interleavedTangentMiddle {b : ℕ} (weight : Fin b → ℝ) :
    Matrix (Fin b ×ₗ Fin 2) (Fin b ×ₗ Fin 2) ℝ :=
  Matrix.reindex ((Equiv.prodComm (Fin 2) (Fin b)).trans toLex)
    ((Equiv.prodComm (Fin 2) (Fin b)).trans toLex)
    (Matrix.blockDiagonal (fun k => tangentCell (weight k)))

@[simp]
theorem interleavedTangentMiddle_apply {b : ℕ} (weight : Fin b → ℝ)
    (i j : Fin b ×ₗ Fin 2) :
    interleavedTangentMiddle weight i j =
      if (ofLex i).1 = (ofLex j).1 then
        tangentCell (weight (ofLex i).1) (ofLex i).2 (ofLex j).2
      else 0 := by
  rfl

/-- The feature-major middle matrix is the canonical reindexing of the
layer-major finite tangent factor used by Phase CA. -/
theorem interleavedTangentMiddle_eq_reindex_upperTangentFeatureMiddle
    {b : ℕ} (weight : Fin b → ℝ) :
    interleavedTangentMiddle weight =
      Matrix.reindex (featureInterleaveEquiv b) (featureInterleaveEquiv b)
        (upperTangentFeatureMiddle (Matrix.diagonal weight)) := by
  rw [interleavedTangentMiddle, upperTangentFeatureMiddle]
  rw [← reindex_upperTangent_diagonal_eq_blockDiagonal]
  ext i j
  simp [Matrix.reindex_apply, featureInterleaveEquiv, layerFeatureEquiv,
    finTwoEquiv, Equiv.boolProdEquivSum]

@[simp]
theorem interleavedTangentMiddle_det {b : ℕ} (weight : Fin b → ℝ) :
    (interleavedTangentMiddle weight).det = 1 := by
  rw [interleavedTangentMiddle, Matrix.det_reindex_self,
    Matrix.det_blockDiagonal]
  simp

theorem interleavedTangentMiddle_entry_nonneg {b : ℕ}
    (weight : Fin b → ℝ) (weight_nonneg : ∀ k, 0 ≤ weight k)
    (i j : Fin b ×ₗ Fin 2) :
    0 ≤ interleavedTangentMiddle weight i j := by
  rw [interleavedTangentMiddle_apply]
  split_ifs
  · exact tangentCell_entry_nonneg _ (weight_nonneg _) _ _
  · exact le_rfl

theorem interleavedTangentMiddle_hasNoncrossingSupport {b : ℕ}
    (weight : Fin b → ℝ) :
    HasNoncrossingSupport
      (fun i j => interleavedTangentMiddle weight i j) := by
  intro x₁ x₂ y₁ y₂ x_lt y_lt
  by_cases crossed_left : (ofLex x₂).1 = (ofLex y₁).1
  · by_cases crossed_right : (ofLex x₁).1 = (ofLex y₂).1
    · have x_feature_le : (ofLex x₁).1 ≤ (ofLex x₂).1 :=
        Prod.Lex.monotone_fst_ofLex x_lt.le
      have y_feature_le : (ofLex y₁).1 ≤ (ofLex y₂).1 :=
        Prod.Lex.monotone_fst_ofLex y_lt.le
      have x_feature_eq : (ofLex x₁).1 = (ofLex x₂).1 := by
        apply le_antisymm x_feature_le
        calc
          (ofLex x₂).1 = (ofLex y₁).1 := crossed_left
          _ ≤ (ofLex y₂).1 := y_feature_le
          _ = (ofLex x₁).1 := crossed_right.symm
      have y_feature_eq : (ofLex y₁).1 = (ofLex y₂).1 := by
        calc
          (ofLex y₁).1 = (ofLex x₂).1 := crossed_left.symm
          _ = (ofLex x₁).1 := x_feature_eq.symm
          _ = (ofLex y₂).1 := crossed_right
      have x_layer_lt : (ofLex x₁).2 < (ofLex x₂).2 :=
        (Prod.Lex.lt_iff'.mp x_lt).2 x_feature_eq
      have y_layer_lt : (ofLex y₁).2 < (ofLex y₂).2 :=
        (Prod.Lex.lt_iff'.mp y_lt).2 y_feature_eq
      left
      change interleavedTangentMiddle weight x₂ y₁ = 0
      rw [interleavedTangentMiddle_apply, if_pos crossed_left]
      have x₂_layer_one : (ofLex x₂).2 = 1 := by omega
      have y₁_layer_zero : (ofLex y₁).2 = 0 := by omega
      simp [x₂_layer_one, y₁_layer_zero]
    · right
      change interleavedTangentMiddle weight x₁ y₂ = 0
      rw [interleavedTangentMiddle_apply, if_neg crossed_right]
  · left
    change interleavedTangentMiddle weight x₂ y₁ = 0
    rw [interleavedTangentMiddle_apply, if_neg crossed_left]

/-- The feature-major diagonal tangent middle factor is totally nonnegative
at every order for arbitrary nonnegative feature weights. -/
theorem interleavedTangentMiddle_totallyNonnegative {b : ℕ}
    (weight : Fin b → ℝ) (weight_nonneg : ∀ k, 0 ≤ weight k) :
    TotallyNonnegative (fun i j => interleavedTangentMiddle weight i j) :=
  totallyNonnegative_of_nonnegative_of_noncrossing
    (interleavedTangentMiddle_entry_nonneg weight weight_nonneg)
    (interleavedTangentMiddle_hasNoncrossingSupport weight)

/-- Total nonnegativity in feature-major order supplies exactly the shuffle
orientation required for every layer-major middle flag minor. -/
theorem upperTangentFeatureMiddle_oriented_minor_nonneg
    {a b : ℕ} (weight : Fin b → ℝ)
    (interleaved_nonnegative : TotallyNonnegativeAtOrder
      (fun i j => interleavedTangentMiddle weight i j) a)
    (rows columns : Fin a → Fin (b + b))
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns) :
    0 ≤ tupleOrientation (featureInterleaveEquiv b ∘ rows) *
        tupleOrientation (featureInterleaveEquiv b ∘ columns) *
        ((upperTangentFeatureMiddle (Matrix.diagonal weight)).submatrix
          rows columns).det := by
  have oriented := oriented_minor_nonneg_of_totallyNonnegativeAtOrder
    interleaved_nonnegative (featureInterleaveEquiv b ∘ rows)
      (featureInterleaveEquiv b ∘ columns)
      ((featureInterleaveEquiv b).injective.comp rows_injective)
      ((featureInterleaveEquiv b).injective.comp columns_injective)
  have sampled :
      minorMatrix (fun i j => interleavedTangentMiddle weight i j)
          (featureInterleaveEquiv b ∘ rows)
          (featureInterleaveEquiv b ∘ columns) =
        (upperTangentFeatureMiddle (Matrix.diagonal weight)).submatrix
          rows columns := by
    rw [interleavedTangentMiddle_eq_reindex_upperTangentFeatureMiddle]
    ext i j
    simp [minorMatrix, Matrix.reindex_apply, Function.comp_apply]
  rw [sampled] at oriented
  exact oriented

/-- Nonnegative diagonal weights unconditionally supply the shuffle-oriented
middle-minor sign required by the upper-tangent Cauchy--Binet factorization. -/
theorem upperTangentFeatureMiddle_diagonal_oriented_minor_nonneg
    {a b : ℕ} (weight : Fin b → ℝ)
    (weight_nonneg : ∀ k, 0 ≤ weight k)
    (rows columns : Fin a → Fin (b + b))
    (rows_injective : Function.Injective rows)
    (columns_injective : Function.Injective columns) :
    0 ≤ tupleOrientation (featureInterleaveEquiv b ∘ rows) *
        tupleOrientation (featureInterleaveEquiv b ∘ columns) *
        ((upperTangentFeatureMiddle (Matrix.diagonal weight)).submatrix
          rows columns).det :=
  upperTangentFeatureMiddle_oriented_minor_nonneg weight
    (interleavedTangentMiddle_totallyNonnegative weight weight_nonneg a)
    rows columns rows_injective columns_injective

end Matrix
