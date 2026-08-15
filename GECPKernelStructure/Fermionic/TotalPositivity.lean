import GECPKernelStructure.SmallInstanceChecks
import Mathlib.Tactic.FinCases

namespace GECPKernelStructure
namespace Fermionic

/-- Exact two-point sign check used to seed, but not establish, the total-positivity hypothesis. -/
theorem geometricSurrogate_twoPoint_sign {q : ℚ} (q_pos : 0 < q) (q_lt_one : q < 1) :
    0 < 1 - q ^ 2 :=
  geometricTwoByTwo_det_pos q_pos q_lt_one

/-- Exact complete-pivot predicate for the minimal rational obstruction witnesses. -/
def IsCompletePivotRat (R : Kernel (Fin 2) (Fin 2) ℚ) (row column : Fin 2) : Prop :=
  ∀ i j, |R i j| ≤ |R row column|

/-- Every entry is positive and the ordered determinant has negative sign. -/
def StrictSignRegularTwo (R : Kernel (Fin 2) (Fin 2) ℚ) : Prop :=
  (∀ i j, 0 < R i j) ∧ R 0 0 * R 1 1 - R 0 1 * R 1 0 < 0

def signRegularLeft : Kernel (Fin 2) (Fin 2) ℚ :=
  fun i j => !![(4 : ℚ), 3; 2, 1] i j

def signRegularRight : Kernel (Fin 2) (Fin 2) ℚ :=
  fun i j => !![(1 : ℚ), 2; 3, 4] i j

/--
Two exact `2 × 2` kernels have the same strict minor-sign pattern but no common
complete pivot. Thus minor signs alone do not determine GECP pivot locations.
This is the minimized obstruction replacing the failed stronger conjecture.
-/
theorem signRegular_not_sufficient_for_parameterIndependent_pivots :
    StrictSignRegularTwo signRegularLeft ∧
    StrictSignRegularTwo signRegularRight ∧
    IsCompletePivotRat signRegularLeft 0 0 ∧
    IsCompletePivotRat signRegularRight 1 1 ∧
    ¬ ∃ row column,
      IsCompletePivotRat signRegularLeft row column ∧
      IsCompletePivotRat signRegularRight row column := by
  constructor
  · constructor
    · intro i j
      fin_cases i <;> fin_cases j <;> norm_num [signRegularLeft]
    · norm_num [signRegularLeft]
  constructor
  · constructor
    · intro i j
      fin_cases i <;> fin_cases j <;> norm_num [signRegularRight]
    · norm_num [signRegularRight]
  constructor
  · intro i j
    fin_cases i <;> fin_cases j <;> norm_num [signRegularLeft]
  constructor
  · intro i j
    fin_cases i <;> fin_cases j <;> norm_num [signRegularRight]
  · rintro ⟨row, column, leftMax, rightMax⟩
    fin_cases row <;> fin_cases column
    · have impossible := rightMax 1 1
      norm_num [signRegularRight] at impossible
    · have impossible := leftMax 0 0
      norm_num [signRegularLeft] at impossible
    · have impossible := leftMax 0 0
      norm_num [signRegularLeft] at impossible
    · have impossible := leftMax 0 0
      norm_num [signRegularLeft] at impossible

end Fermionic
end GECPKernelStructure
