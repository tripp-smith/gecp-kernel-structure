import GECPKernelStructure.Fermionic.AdaptiveCoverageObstructionCore

namespace GECPKernelStructure.Fermionic

open GECP Matrix

set_option maxRecDepth 16384

def phaseDBStep (n : ℕ) : Fin 7 := ⟨n % 7, Nat.mod_lt _ (by decide)⟩
def phaseDBFreshRows : Fin 2 → Fin 9 := ![6, 4]
def phaseDBFreshColumns : Fin 2 → Fin 9 := ![1, 2]
def phaseDBFeatureScale : ℤ := phaseDBDivisors 7

def phaseDBRawLeftQ : Matrix (Fin 9) (Fin 13) ℚ := fun i k ↦ phaseDBRawLeft i k
def phaseDBRawRightQ : Matrix (Fin 13) (Fin 9) ℚ := fun k j ↦ phaseDBRawRight k j

def phaseDBLeftProjectorQ : ℕ → Matrix (Fin 9) (Fin 9) ℚ
  | 0 => 1
  | n + 1 => companionRowUpdate (phaseDBResidualQ (phaseDBStep n).castSucc)
      (phaseDBLeftProjectorQ n) (phaseDBPivotRows (phaseDBStep n))
      (phaseDBPivotColumns (phaseDBStep n)) (phaseDB_pivot_ne (phaseDBStep n))

def phaseDBRightProjectorQ : ℕ → Matrix (Fin 9) (Fin 9) ℚ
  | 0 => 1
  | n + 1 => companionColumnUpdate (phaseDBResidualQ (phaseDBStep n).castSucc)
      (phaseDBRightProjectorQ n) (phaseDBPivotRows (phaseDBStep n))
      (phaseDBPivotColumns (phaseDBStep n)) (phaseDB_pivot_ne (phaseDBStep n))

def phaseDBLeftFeaturesQ : ℕ → Matrix (Fin 9) (Fin 13) ℚ
  | 0 => phaseDBRawLeftQ
  | n + 1 => companionRowUpdate (phaseDBResidualQ (phaseDBStep n).castSucc)
      (phaseDBLeftFeaturesQ n) (phaseDBPivotRows (phaseDBStep n))
      (phaseDBPivotColumns (phaseDBStep n)) (phaseDB_pivot_ne (phaseDBStep n))

def phaseDBRightFeaturesQ : ℕ → Matrix (Fin 13) (Fin 9) ℚ
  | 0 => phaseDBRawRightQ
  | n + 1 => companionColumnUpdate (phaseDBResidualQ (phaseDBStep n).castSucc)
      (phaseDBRightFeaturesQ n) (phaseDBPivotRows (phaseDBStep n))
      (phaseDBPivotColumns (phaseDBStep n)) (phaseDB_pivot_ne (phaseDBStep n))

theorem phaseDB_leftFeatures_projector (n : ℕ) :
    phaseDBLeftFeaturesQ n = phaseDBLeftProjectorQ n * phaseDBRawLeftQ := by
  induction n with
  | zero => simp [phaseDBLeftFeaturesQ, phaseDBLeftProjectorQ]
  | succ n ih =>
    simp only [phaseDBLeftFeaturesQ, phaseDBLeftProjectorQ, ih,
      companionRowUpdate_mul]

theorem phaseDB_rightFeatures_projector (n : ℕ) :
    phaseDBRightFeaturesQ n = phaseDBRawRightQ * phaseDBRightProjectorQ n := by
  induction n with
  | zero => simp [phaseDBRightFeaturesQ, phaseDBRightProjectorQ]
  | succ n ih =>
    simp only [phaseDBRightFeaturesQ, phaseDBRightProjectorQ, ih,
      companionColumnUpdate_mul]

set_option maxHeartbeats 20000000 in
theorem phaseDB_left_projector_certificate (i : Fin 2) (j : Fin 9) :
    (phaseDBFeatureScale : ℚ) * phaseDBLeftProjectorQ 7 (phaseDBFreshRows i) j =
      (phaseDBLeftProjectorNumerators i j : ℚ) := by
  revert i j
  decide +kernel

set_option maxHeartbeats 20000000 in
theorem phaseDB_right_projector_certificate (i : Fin 9) (j : Fin 2) :
    (phaseDBFeatureScale : ℚ) * phaseDBRightProjectorQ 7 i (phaseDBFreshColumns j) =
      (phaseDBRightProjectorNumerators i j : ℚ) := by
  revert i j
  decide +kernel

theorem phaseDB_left_numerator_product (i : Fin 2) (k : Fin 13) :
    phaseDBLeftNumerators i k = ∑ j : Fin 9,
      phaseDBLeftProjectorNumerators i j * phaseDBRawLeft j k := by
  fin_cases i <;> fin_cases k <;> decide +kernel

theorem phaseDB_right_numerator_product (k : Fin 13) (j : Fin 2) :
    phaseDBRightNumerators k j = ∑ i : Fin 9,
      phaseDBRawRight k i * phaseDBRightProjectorNumerators i j := by
  fin_cases k <;> fin_cases j <;> decide +kernel

theorem phaseDB_left_feature_certificate (i : Fin 2) (k : Fin 13) :
    (phaseDBFeatureScale : ℚ) * phaseDBLeftFeaturesQ 7 (phaseDBFreshRows i) k =
      (phaseDBLeftNumerators i k : ℚ) := by
  rw [phaseDB_leftFeatures_projector]
  simp only [Matrix.mul_apply, Finset.mul_sum, ← mul_assoc,
    phaseDB_left_projector_certificate, phaseDBRawLeftQ]
  exact_mod_cast (phaseDB_left_numerator_product i k).symm

theorem phaseDB_right_feature_certificate (k : Fin 13) (j : Fin 2) :
    (phaseDBFeatureScale : ℚ) * phaseDBRightFeaturesQ 7 k (phaseDBFreshColumns j) =
      (phaseDBRightNumerators k j : ℚ) := by
  rw [phaseDB_rightFeatures_projector]
  simp only [Matrix.mul_apply, Finset.mul_sum, phaseDBRawRightQ]
  calc
    _ = ∑ i : Fin 9, (phaseDBRawRight k i : ℚ) *
        ((phaseDBFeatureScale : ℚ) * phaseDBRightProjectorQ 7 i (phaseDBFreshColumns j)) := by
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ = _ := by
      simp only [phaseDB_right_projector_certificate]
      exact_mod_cast (phaseDB_right_numerator_product k j).symm

def phaseDBBlockQ (degrees : Finset (Fin 13)) : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j ↦ ∑ k ∈ degrees,
    (phaseDBLeftNumerators i k : ℚ) * (phaseDBRightNumerators k j : ℚ)

def phaseDBActualBlockQ (degrees : Finset (Fin 13)) : Matrix (Fin 2) (Fin 2) ℚ :=
  fun i j ↦ ∑ k ∈ degrees,
    phaseDBLeftFeaturesQ 7 (phaseDBFreshRows i) k *
      phaseDBRightFeaturesQ 7 k (phaseDBFreshColumns j)

theorem phaseDB_block_certificate (degrees : Finset (Fin 13)) :
    phaseDBBlockQ degrees = (phaseDBFeatureScale : ℚ) ^ 2 • phaseDBActualBlockQ degrees := by
  ext i j
  simp only [phaseDBBlockQ, phaseDBActualBlockQ, Matrix.smul_apply, smul_eq_mul,
    Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k hk
  rw [← phaseDB_left_feature_certificate, ← phaseDB_right_feature_certificate]
  ring

def phaseDBMixedQ (low high : Matrix (Fin 2) (Fin 2) ℚ) : ℚ :=
  low 0 0 * high 1 1 + high 0 0 * low 1 1 -
    low 0 1 * high 1 0 - high 0 1 * low 1 0

def phaseDBTransferQ (low moved upper : Matrix (Fin 2) (Fin 2) ℚ) : ℚ :=
  low.det * phaseDBMixedQ moved upper ^ 2 -
    phaseDBMixedQ low (moved + upper) * phaseDBMixedQ low moved *
      phaseDBMixedQ moved upper +
    (moved + upper).det * phaseDBMixedQ low moved ^ 2

theorem phaseDBTransferQ_cast (low moved upper : Matrix (Fin 2) (Fin 2) ℚ) :
    (phaseDBTransferQ low moved upper : ℝ) =
      quadraticTransferNodeNumerator
        (low.map (fun x ↦ (x : ℝ))).det
        (finTwoPencilMixed (low.map (fun x ↦ (x : ℝ)))
          ((moved + upper).map (fun x ↦ (x : ℝ))))
        ((moved + upper).map (fun x ↦ (x : ℝ))).det
        (finTwoPencilMixed (low.map (fun x ↦ (x : ℝ)))
          (moved.map (fun x ↦ (x : ℝ))))
        (finTwoPencilMixed (moved.map (fun x ↦ (x : ℝ)))
          (upper.map (fun x ↦ (x : ℝ)))) := by
  simp only [phaseDBTransferQ, phaseDBMixedQ, quadraticTransferNodeNumerator,
    finTwoPencilMixed, Matrix.det_fin_two, Matrix.map_apply]
  push_cast
  rfl

def phaseDBLowerSeven : Finset (Fin 13) := {0, 1, 2, 3, 4, 5, 6}
def phaseDBUpperSeven : Finset (Fin 13) := {8, 9, 10, 11, 12}
def phaseDBLowerEight : Finset (Fin 13) := {0, 1, 2, 3, 4, 5, 6, 7}
def phaseDBUpperEight : Finset (Fin 13) := {9, 10, 11, 12}

set_option maxHeartbeats 20000000 in
theorem phaseDB_adjacent_transferNodes_pos :
    0 < phaseDBTransferQ (phaseDBBlockQ phaseDBLowerSeven) (phaseDBBlockQ {7})
      (phaseDBBlockQ phaseDBUpperSeven) ∧
    0 < phaseDBTransferQ (phaseDBBlockQ phaseDBLowerEight) (phaseDBBlockQ {8})
      (phaseDBBlockQ phaseDBUpperEight) := by
  decide +kernel

set_option maxHeartbeats 20000000 in
theorem phaseDB_cutoff_discriminant_pos :
    0 < phaseDBMixedQ (phaseDBBlockQ phaseDBLowerEight)
          (phaseDBBlockQ ({8} ∪ phaseDBUpperEight)) ^ 2 -
      4 * (phaseDBBlockQ phaseDBLowerEight).det *
        (phaseDBBlockQ ({8} ∪ phaseDBUpperEight)).det := by
  decide +kernel

theorem phaseDB_numerators_reconstruct (i j : Fin 2) :
    (∑ k : Fin 13, phaseDBLeftNumerators i k * phaseDBRightNumerators k j) =
      phaseDBFeatureScale * phaseDBState 7 (phaseDBFreshRows i) (phaseDBFreshColumns j) := by
  fin_cases i <;> fin_cases j <;> decide +kernel

theorem phaseDB_conditioned_features_reconstruct :
    phaseDBActualBlockQ Finset.univ =
      Matrix.of (fun i j ↦ phaseDBResidualQ 7 (phaseDBFreshRows i) (phaseDBFreshColumns j)) := by
  ext i j
  have hs : (phaseDBFeatureScale : ℚ) ≠ 0 := by
    exact_mod_cast (phaseDB_divisors_pos 7).ne'
  have hc : (phaseDBFeatureScale : ℚ) ^ 2 * phaseDBActualBlockQ Finset.univ i j =
      (phaseDBFeatureScale : ℚ) *
        (phaseDBState 7 (phaseDBFreshRows i) (phaseDBFreshColumns j) : ℚ) := by
    calc
      _ = phaseDBBlockQ Finset.univ i j := by rw [phaseDB_block_certificate]; rfl
      _ = _ := by
        simp only [phaseDBBlockQ]
        exact_mod_cast phaseDB_numerators_reconstruct i j
  apply mul_left_cancel₀ (pow_ne_zero 2 hs)
  rw [hc]
  change (phaseDBFeatureScale : ℚ) *
      (phaseDBState 7 (phaseDBFreshRows i) (phaseDBFreshColumns j) : ℚ) =
    (phaseDBFeatureScale : ℚ) ^ 2 *
      ((phaseDBState 7 (phaseDBFreshRows i) (phaseDBFreshColumns j) : ℚ) /
        (phaseDBFeatureScale : ℚ))
  field_simp

theorem phaseDBTransferQ_smul (s : ℚ) (low moved upper : Matrix (Fin 2) (Fin 2) ℚ) :
    phaseDBTransferQ (s • low) (s • moved) (s • upper) =
      s ^ 6 * phaseDBTransferQ low moved upper := by
  simp only [phaseDBTransferQ, phaseDBMixedQ, Matrix.det_fin_two,
    Matrix.smul_apply, Matrix.add_apply, smul_eq_mul]
  ring

theorem phaseDB_actual_adjacent_transferNodes_pos :
    0 < phaseDBTransferQ (phaseDBActualBlockQ phaseDBLowerSeven)
      (phaseDBActualBlockQ {7}) (phaseDBActualBlockQ phaseDBUpperSeven) ∧
    0 < phaseDBTransferQ (phaseDBActualBlockQ phaseDBLowerEight)
      (phaseDBActualBlockQ {8}) (phaseDBActualBlockQ phaseDBUpperEight) := by
  have h := phaseDB_adjacent_transferNodes_pos
  simp only [phaseDB_block_certificate, phaseDBTransferQ_smul] at h
  have hp : (0 : ℚ) < ((phaseDBFeatureScale : ℚ) ^ 2) ^ 6 := by
    have hs : (0 : ℚ) < (phaseDBFeatureScale : ℚ) := by
      exact_mod_cast phaseDB_divisors_pos 7
    positivity
  exact ⟨(mul_pos_iff_of_pos_left hp).mp h.1, (mul_pos_iff_of_pos_left hp).mp h.2⟩

theorem phaseDB_not_adaptive_coverage :
    ¬ (phaseDBTransferQ (phaseDBActualBlockQ phaseDBLowerSeven)
        (phaseDBActualBlockQ {7}) (phaseDBActualBlockQ phaseDBUpperSeven) ≤ 0 ∨
      phaseDBTransferQ (phaseDBActualBlockQ phaseDBLowerEight)
        (phaseDBActualBlockQ {8}) (phaseDBActualBlockQ phaseDBUpperEight) ≤ 0) := by
  obtain ⟨h7, h8⟩ := phaseDB_actual_adjacent_transferNodes_pos
  exact fun h ↦ h.elim (not_le_of_gt h7) (not_le_of_gt h8)

set_option maxHeartbeats 20000000 in
theorem phaseDB_endpoints_and_transfers_pos :
    0 < (phaseDBBlockQ phaseDBLowerSeven).det ∧
    0 < (phaseDBBlockQ ({7} ∪ phaseDBUpperSeven)).det ∧
    0 < (phaseDBBlockQ phaseDBLowerEight).det ∧
    0 < (phaseDBBlockQ ({8} ∪ phaseDBUpperEight)).det ∧
    0 < (phaseDBBlockQ (phaseDBLowerEight ∪ {8})).det ∧
    0 < (phaseDBBlockQ phaseDBUpperEight).det ∧
    0 < phaseDBMixedQ (phaseDBBlockQ phaseDBLowerSeven) (phaseDBBlockQ {7}) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ {7}) (phaseDBBlockQ phaseDBUpperSeven) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ phaseDBLowerEight) (phaseDBBlockQ {8}) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ {8}) (phaseDBBlockQ phaseDBUpperEight) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ phaseDBLowerSeven)
      (phaseDBBlockQ ({7} ∪ phaseDBUpperSeven)) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ phaseDBLowerEight)
      (phaseDBBlockQ ({8} ∪ phaseDBUpperEight)) ∧
    0 < phaseDBMixedQ (phaseDBBlockQ (phaseDBLowerEight ∪ {8}))
      (phaseDBBlockQ phaseDBUpperEight) := by
  decide +kernel

theorem phaseDBMixedQ_smul (s : ℚ) (low high : Matrix (Fin 2) (Fin 2) ℚ) :
    phaseDBMixedQ (s • low) (s • high) = s ^ 2 * phaseDBMixedQ low high := by
  simp only [phaseDBMixedQ, Matrix.smul_apply, smul_eq_mul]
  ring

theorem phaseDB_actual_cutoff_discriminant_pos :
    0 < phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerEight)
          (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight)) ^ 2 -
      4 * (phaseDBActualBlockQ phaseDBLowerEight).det *
        (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight)).det := by
  have h := phaseDB_cutoff_discriminant_pos
  simp only [phaseDB_block_certificate, phaseDBMixedQ_smul, Matrix.det_smul,
    Fintype.card_fin] at h
  have hp : (0 : ℚ) < (phaseDBFeatureScale : ℚ) ^ 8 := by
    have hs : (0 : ℚ) < (phaseDBFeatureScale : ℚ) := by
      exact_mod_cast phaseDB_divisors_pos 7
    positivity
  apply (mul_pos_iff_of_pos_left hp).mp
  calc
    0 < ((((phaseDBFeatureScale : ℚ) ^ 2) ^ 2) *
        phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerEight)
          (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight))) ^ 2 -
        4 * ((((phaseDBFeatureScale : ℚ) ^ 2) ^ 2) *
          (phaseDBActualBlockQ phaseDBLowerEight).det) *
          ((((phaseDBFeatureScale : ℚ) ^ 2) ^ 2) *
            (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight)).det) := h
    _ = _ := by ring

theorem phaseDB_actual_common_cutoff :
    phaseDBActualBlockQ phaseDBLowerSeven + phaseDBActualBlockQ {7} =
        phaseDBActualBlockQ phaseDBLowerEight ∧
      phaseDBActualBlockQ phaseDBUpperSeven =
        phaseDBActualBlockQ {8} + phaseDBActualBlockQ phaseDBUpperEight := by
  have add_blocks (s t : Finset (Fin 13)) (hd : Disjoint s t) :
      phaseDBActualBlockQ (s ∪ t) = phaseDBActualBlockQ s + phaseDBActualBlockQ t := by
    ext i j
    exact Finset.sum_union hd
  constructor
  · rw [show phaseDBLowerEight = phaseDBLowerSeven ∪ {7} by decide +kernel]
    exact (add_blocks _ _ (by decide +kernel)).symm
  · rw [show phaseDBUpperSeven = {8} ∪ phaseDBUpperEight by decide +kernel]
    exact add_blocks _ _ (by decide +kernel)

theorem phaseDB_block_det_pos_iff (s : Finset (Fin 13)) :
    0 < (phaseDBBlockQ s).det ↔ 0 < (phaseDBActualBlockQ s).det := by
  rw [phaseDB_block_certificate, Matrix.det_smul]
  have hp : (0 : ℚ) < ((phaseDBFeatureScale : ℚ) ^ 2) ^ Fintype.card (Fin 2) := by
    have hs : (0 : ℚ) < (phaseDBFeatureScale : ℚ) := by
      exact_mod_cast phaseDB_divisors_pos 7
    positivity
  exact mul_pos_iff_of_pos_left hp

theorem phaseDB_block_mixed_pos_iff (s t : Finset (Fin 13)) :
    0 < phaseDBMixedQ (phaseDBBlockQ s) (phaseDBBlockQ t) ↔
      0 < phaseDBMixedQ (phaseDBActualBlockQ s) (phaseDBActualBlockQ t) := by
  rw [phaseDB_block_certificate, phaseDB_block_certificate, phaseDBMixedQ_smul]
  have hp : (0 : ℚ) < ((phaseDBFeatureScale : ℚ) ^ 2) ^ 2 := by
    have hs : (0 : ℚ) < (phaseDBFeatureScale : ℚ) := by
      exact_mod_cast phaseDB_divisors_pos 7
    positivity
  exact mul_pos_iff_of_pos_left hp

theorem phaseDB_actual_endpoints_and_transfers_pos :
    0 < (phaseDBActualBlockQ phaseDBLowerSeven).det ∧
    0 < (phaseDBActualBlockQ ({7} ∪ phaseDBUpperSeven)).det ∧
    0 < (phaseDBActualBlockQ phaseDBLowerEight).det ∧
    0 < (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight)).det ∧
    0 < (phaseDBActualBlockQ (phaseDBLowerEight ∪ {8})).det ∧
    0 < (phaseDBActualBlockQ phaseDBUpperEight).det ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerSeven) (phaseDBActualBlockQ {7}) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ {7}) (phaseDBActualBlockQ phaseDBUpperSeven) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerEight) (phaseDBActualBlockQ {8}) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ {8}) (phaseDBActualBlockQ phaseDBUpperEight) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerSeven)
      (phaseDBActualBlockQ ({7} ∪ phaseDBUpperSeven)) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ phaseDBLowerEight)
      (phaseDBActualBlockQ ({8} ∪ phaseDBUpperEight)) ∧
    0 < phaseDBMixedQ (phaseDBActualBlockQ (phaseDBLowerEight ∪ {8}))
      (phaseDBActualBlockQ phaseDBUpperEight) := by
  simpa only [phaseDB_block_det_pos_iff, phaseDB_block_mixed_pos_iff] using
    phaseDB_endpoints_and_transfers_pos

end GECPKernelStructure.Fermionic
