import GECPKernelStructure.Fermionic.DyadicPartition

namespace GECPKernelStructure
namespace Fermionic

/-- Natural-number node count for a scale-by-order exponential construction. -/
def expFamilyNodeCount (scales order : ℕ) : ℕ := (scales + 1) * (order + 1)

theorem expFamilyNodeCount_pos (scales order : ℕ) :
    0 < expFamilyNodeCount scales order := by
  simp [expFamilyNodeCount]

theorem expFamilyNodeCount_mono {s₁ s₂ n₁ n₂ : ℕ}
    (hs : s₁ ≤ s₂) (hn : n₁ ≤ n₂) :
    expFamilyNodeCount s₁ n₁ ≤ expFamilyNodeCount s₂ n₂ := by
  exact Nat.mul_le_mul (Nat.add_le_add_right hs 1) (Nat.add_le_add_right hn 1)

end Fermionic
end GECPKernelStructure
