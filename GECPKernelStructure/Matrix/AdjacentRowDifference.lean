import GECPKernelStructure.GECP.TotalNonnegativeIntegral

namespace Matrix

open MeasureTheory

variable {R : Type*} [CommRing R]

/-- Replace each nonfinal row by its difference with the next row. -/
def successiveRowDifference {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) R) :
    Matrix (Fin (n + 1)) (Fin (n + 1)) R :=
  fun i => Fin.lastCases (A (Fin.last n))
    (fun k => A k.castSucc - A k.succ) i

@[simp]
theorem successiveRowDifference_last {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) R) :
    successiveRowDifference A (Fin.last n) = A (Fin.last n) := by
  simp [successiveRowDifference]

@[simp]
theorem successiveRowDifference_castSucc {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) R) (i : Fin n) :
    successiveRowDifference A i.castSucc = A i.castSucc - A i.succ := by
  simp [successiveRowDifference]

section

variable [IsDomain R] [CharZero R]

/-- Successive adjacent row differences preserve the determinant. -/
theorem det_successiveRowDifference {n : ℕ}
    (A : Matrix (Fin (n + 1)) (Fin (n + 1)) R) :
    (successiveRowDifference A).det = A.det := by
  let reverseA := A.submatrix Fin.rev id
  let reverseDifference := (successiveRowDifference A).submatrix Fin.rev id
  have reversed_det : reverseA.det = reverseDifference.det := by
    apply det_eq_of_forall_row_eq_smul_add_pred (fun _ => 1)
    · intro j
      change A (Fin.rev 0) j = successiveRowDifference A (Fin.rev 0) j
      rw [Fin.rev_zero, successiveRowDifference_last]
    · intro i j
      change A (Fin.rev i.succ) j =
        successiveRowDifference A (Fin.rev i.succ) j +
          1 * A (Fin.rev i.castSucc) j
      rw [Fin.rev_succ, Fin.rev_castSucc, successiveRowDifference_castSucc]
      simp
  let rev : Equiv.Perm (Fin (n + 1)) := Fin.revPerm
  have permutedA := det_permute rev A
  have permutedDifference := det_permute rev (successiveRowDifference A)
  have signed_eq :
      (((rev.sign : ℤ) : R) * A.det) =
        ((rev.sign : ℤ) : R) * (successiveRowDifference A).det := by
    rw [← permutedA, ← permutedDifference]
    exact reversed_det
  have sign_ne : (((rev.sign : ℤ) : R)) ≠ 0 :=
    Int.cast_ne_zero.mpr (Units.ne_zero rev.sign)
  have original_eq_difference : A.det = (successiveRowDifference A).det :=
    mul_left_cancel₀ sign_ne signed_eq
  exact original_eq_difference.symm

end

/--
Common-upper row integrals become adjacent-interval integrals after successive
row differences, so their determinants agree.
-/
theorem det_commonUpperRowIntegral_eq_adjacent {n : ℕ}
    (f : ℝ → Fin (n + 1) → ℝ) (knots : Fin (n + 1) → ℝ) (upper : ℝ)
    (integrable : ∀ i j,
      IntervalIntegrable (fun x => f x j) volume (knots i) upper) :
    Matrix.det (fun i j => ∫ x in knots i..upper, f x j) =
      Matrix.det (fun i j => Fin.lastCases
        (∫ x in knots (Fin.last n)..upper, f x j)
        (fun k => ∫ x in knots k.castSucc..knots k.succ, f x j) i) := by
  let tails : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ :=
    fun i j => ∫ x in knots i..upper, f x j
  have difference_eq : successiveRowDifference tails =
      fun i j => Fin.lastCases
        (∫ x in knots (Fin.last n)..upper, f x j)
        (fun k => ∫ x in knots k.castSucc..knots k.succ, f x j) i := by
    ext i j
    refine Fin.lastCases ?_ (fun k => ?_) i
    · rw [successiveRowDifference_last]
      simp [tails]
    · rw [successiveRowDifference_castSucc]
      change tails k.castSucc j - tails k.succ j = _
      rw [Fin.lastCases_castSucc]
      have adjacent_integrable :
          IntervalIntegrable (fun x => f x j) volume
            (knots k.castSucc) (knots k.succ) :=
        (integrable k.castSucc j).trans (integrable k.succ j).symm
      change (∫ x in knots k.castSucc..upper, f x j) -
          (∫ x in knots k.succ..upper, f x j) =
        ∫ x in knots k.castSucc..knots k.succ, f x j
      have add_tails := intervalIntegral.integral_add_adjacent_intervals
        adjacent_integrable (integrable k.succ j)
      linarith
  change tails.det = _
  rw [← difference_eq]
  exact (det_successiveRowDifference tails).symm

end Matrix
