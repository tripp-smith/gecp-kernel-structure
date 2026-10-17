import GECPKernelStructure.GECP.SignRegularity
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

end Matrix
