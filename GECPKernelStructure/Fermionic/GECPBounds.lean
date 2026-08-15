namespace GECPKernelStructure
namespace Fermionic

/-- The only classes of rigorous outcome allowed to close structural Milestone E. -/
inductive GECPResearchOutcome where
  | targetRate
  | weakerImprovedRate
  | containingStructuralClass
  | certifiedObstruction
  deriving DecidableEq, Repr

end Fermionic
end GECPKernelStructure
