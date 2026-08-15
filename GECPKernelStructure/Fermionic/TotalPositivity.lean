import GECPKernelStructure.SmallInstanceChecks

namespace GECPKernelStructure
namespace Fermionic

/-- Exact two-point sign check used to seed, but not establish, the total-positivity hypothesis. -/
theorem geometricSurrogate_twoPoint_sign {q : ℚ} (q_pos : 0 < q) (q_lt_one : q < 1) :
    0 < 1 - q ^ 2 :=
  geometricTwoByTwo_det_pos q_pos q_lt_one

end Fermionic
end GECPKernelStructure
