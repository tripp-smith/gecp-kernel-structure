import GECPKernelStructure.GECP.FactorizedTangentCollapse
import GECPKernelStructure.Fermionic.ExponentialTotalPositivity

namespace GECPKernelStructure
namespace Fermionic

open GECP Matrix

/-- Finite positive-Taylor tangent minors collapse to one determinant of the
conditioned monomial flags and the diagonal degree weights.  This is the
finite-dimensional cancellation-preserving interface for the nonnegative-
weight problem left by Phase CE. -/
theorem expTaylor_upperTangent_diagonal_distinguishedMinor_det
    {n terms : ℕ} {fresh : Type*}
    [Fintype fresh] [DecidableEq fresh]
    (rowNodes columnNodes : Fin n → ℝ)
    (weight : Fin terms → ℝ)
    (run : GECP.Run (fun i j =>
      (expTaylorLeft (terms := terms) rowNodes *
        expTaylorRight (terms := terms) columnNodes) i j))
    (rows columns : fresh → Fin n) :
    (Matrix.submatrix
        (Matrix.upperTangentMatrix
          (expTaylorLeft (terms := terms) rowNodes *
            expTaylorRight (terms := terms) columnNodes)
          (expTaylorLeft (terms := terms) rowNodes *
            Matrix.diagonal weight *
              expTaylorRight (terms := terms) columnNodes))
        (Matrix.upperTangentDoubleBorderRow run.selectedRow rows)
        (Matrix.upperTangentDoubleBorderColumn run.selectedColumn columns)).det =
      run.selectedCore.det ^ 2 *
        (Matrix.submatrix
          (run.conditionRows (expTaylorLeft (terms := terms) rowNodes) *
            Matrix.diagonal weight *
              run.conditionColumns
                (expTaylorRight (terms := terms) columnNodes))
          rows columns).det := by
  exact run.upperTangentMatrix_factorized_distinguishedMinor_det
    (expTaylorLeft (terms := terms) rowNodes)
    (Matrix.diagonal weight)
    (expTaylorRight (terms := terms) columnNodes) rows columns

end Fermionic
end GECPKernelStructure
