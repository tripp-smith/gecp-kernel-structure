namespace GECPKernelStructure

universe u v w

/-- A two-variable kernel with row domain `α`, column domain `β`, and values in `𝕜`. -/
abbrev Kernel (α : Type u) (β : Type v) (𝕜 : Type w) := α → β → 𝕜

end GECPKernelStructure
