import Mathlib.Data.Finset.Max
import Mathlib.Data.Real.Basic

namespace GECPKernelStructure

/-- A finite-grid complete pivot maximizes absolute value over an explicit set. -/
def IsFiniteCompletePivot {ι : Type*} (grid : Finset ι) (value : ι → ℝ) (pivot : ι) : Prop :=
  pivot ∈ grid ∧ ∀ point ∈ grid, abs (value point) ≤ abs (value pivot)

end GECPKernelStructure
