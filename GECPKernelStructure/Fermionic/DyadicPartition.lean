import Mathlib.Data.Nat.Log

namespace GECPKernelStructure
namespace Fermionic

/-- Explicit dyadic scale count used by the constructive research contract. -/
def dyadicScaleCount (cutoff : ℕ) : ℕ := Nat.log 2 cutoff + 1

theorem dyadicScaleCount_pos (cutoff : ℕ) : 0 < dyadicScaleCount cutoff := by
  simp [dyadicScaleCount]

end Fermionic
end GECPKernelStructure
