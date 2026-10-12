import GECPKernelStructure.Matrix.MixedColumnCoefficients

namespace Matrix

open Finset

/-- The first `p` coordinates of a combined fixed-plus-feature index. -/
def isFixedCoordinate {p k : ℕ} (index : Fin (p + k)) : Prop :=
  index.1 < p

instance {p k : ℕ} : DecidablePred (isFixedCoordinate (p := p) (k := k)) :=
  fun index => inferInstanceAs (Decidable (index.1 < p))

@[simp]
theorem isFixedCoordinate_castAdd {p k : ℕ} (index : Fin p) :
    isFixedCoordinate (p := p) (k := k) (Fin.castAdd k index) := by
  exact index.2

@[simp]
theorem not_isFixedCoordinate_natAdd {p k : ℕ} (index : Fin k) :
    ¬isFixedCoordinate (p := p) (k := k) (Fin.natAdd p index) := by
  simp [isFixedCoordinate]

/-- Adjoin a fixed coordinate block `U` before an ordinary left feature
matrix `A`. -/
def combinedFeatureLeft {p n k : ℕ}
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ) :
    Matrix (Fin (p + n)) (Fin (p + k)) ℝ :=
  fun i => Fin.addCases (fun j => U i j) (fun j => A i j)

/-- Adjoin a fixed coordinate block `V` above an ordinary right feature
matrix `B`. -/
def combinedFeatureRight {p n k : ℕ}
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    Matrix (Fin (p + k)) (Fin (p + n)) ℝ :=
  Fin.addCases (fun i j => V i j) (fun i j => B i j)

/-- The marked part of the combined factorization is exactly the fixed
rank-`p` product. -/
theorem highFeatureLeft_combined_mul {p n k : ℕ}
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    highFeatureLeft (isFixedCoordinate (p := p) (k := k))
        (combinedFeatureLeft U A) * combinedFeatureRight V B = U * V := by
  ext i j
  rw [Matrix.mul_apply, Fin.sum_univ_add, Matrix.mul_apply]
  simp [highFeatureLeft, combinedFeatureLeft, combinedFeatureRight]

/-- The unmarked part of the combined factorization is exactly the ordinary
feature product. -/
theorem lowFeatureLeft_combined_mul {p n k : ℕ}
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    lowFeatureLeft (isFixedCoordinate (p := p) (k := k))
        (combinedFeatureLeft U A) * combinedFeatureRight V B = A * B := by
  ext i j
  rw [Matrix.mul_apply, Fin.sum_univ_add, Matrix.mul_apply]
  simp [lowFeatureLeft, combinedFeatureLeft, combinedFeatureRight]

/-- Weighting the fixed coordinate block produces the physical determinant
pencil `mu • (U*V) + A*B`. -/
theorem weightedCombinedFeatureProduct {p n k : ℕ} (mu : ℝ)
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    weightedFeatureLeft (isFixedCoordinate (p := p) (k := k)) mu
        (combinedFeatureLeft U A) * combinedFeatureRight V B =
      mu • (U * V) + A * B := by
  rw [weightedFeatureProduct_eq_low_add_smul_high,
    highFeatureLeft_combined_mul, lowFeatureLeft_combined_mul]
  exact add_comm _ _

/-- The graded Plücker pairing at the full fixed-coordinate degree is the
physical mixed-column coefficient of the fixed product against the feature
product. -/
theorem fixedBlockPluckerPairing_eq_mixedColumnCoefficient
    {p n k : ℕ}
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    blockFeatureCoefficient 1 (isFixedCoordinate (p := p) (k := k))
        (combinedFeatureLeft U A) (combinedFeatureRight V B) p =
      mixedColumnCoefficient p (U * V) (A * B) := by
  rw [blockFeatureCoefficient_eq_mixedColumnCoefficient,
    highFeatureLeft_combined_mul, lowFeatureLeft_combined_mul]

/-- Expanded form of the fixed-degree identity.  Individual summands may
have either sign; the theorem deliberately preserves their signed sum. -/
theorem fixedBlockPluckerPairing_eq_sum {p n k : ℕ}
    (U : Matrix (Fin (p + n)) (Fin p) ℝ)
    (V : Matrix (Fin p) (Fin (p + n)) ℝ)
    (A : Matrix (Fin (p + n)) (Fin k) ℝ)
    (B : Matrix (Fin k) (Fin (p + n)) ℝ) :
    mixedColumnCoefficient p (U * V) (A * B) =
      ∑ q : {q : Fin (p + n) → Fin (p + k) // StrictMono q},
        if highFeatureCount (isFixedCoordinate (p := p) (k := k)) q.1 = p then
          signedOrderedCauchyBinetTerm 1
            (combinedFeatureLeft U A) (combinedFeatureRight V B) q
        else 0 := by
  rw [← fixedBlockPluckerPairing_eq_mixedColumnCoefficient]
  rfl

end Matrix
