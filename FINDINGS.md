# Findings and claim ledger

## Proved

- A successful dependent GECP run leaves every selected row and column zero in
  its final residual.
- The finite core obtained by evaluating the original kernel at the ordered
  rows and columns of any successful dependent `GECP.Run` has determinant
  equal to its actual residual-pivot product. The stored nonzero witnesses
  therefore prove that core is nonsingular; no external LDU assumption is
  required.
- A positive diagonal GECP/Cholesky update on a finite real PSD matrix remains
  PSD. The proof constructs an exact elimination map `P` and verifies the
  residual factorization `R = Pᴴ M P`.
- Every entry of a finite real PSD residual obeys
  `|Rᵢⱼ|² ≤ Rᵢᵢ Rⱼⱼ`; bounding all diagonal entries is equivalent to bounding
  all absolute entries. Hence a diagonal maximizer is a complete-pivot
  maximizer for the canonical diagonal-preferring rule.
- Recursive traces using the least diagonal maximizer are equivalent in both
  directions for complete-pivot GECP and diagonal-pivoted Cholesky.
- A Lipschitz residual diagonal has the stated conditional fill-distance bound.
- The fermionic kernel has the exact reflection and centered forms and is
  jointly continuous.
- A Lipschitz grid cover yields an explicit continuous supremum bound, and a
  certified upper bound implies an approximate-pivot inequality.
- Uniform kernel error transfers to the Green-function error with the
  `L¹` norm of the spectral density as weight.

## Observed

- All minors enumerated in the 21 exact matrices `q^(ij)` for sizes 2–8 and
  `q ∈ {1/2, 2/3, 3/4}` are nonzero with the predicted sign.
- Exact complete-pivot paths vary with `q` by size eight. This rules out a
  single `q`-independent pivot order for this surrogate family.
- The 21-case float64 finite-grid fermionic census converged on its sampled
  grids with ranks 6–35. Ranks plateau for the three largest cutoffs, which is
  potentially a float64/grid-resolution effect and is not promoted to a
  structural conclusion.
- The fixed two-delta fixture is recovered below `1e-8` using two atoms.

## Conjectured

The target continuous fermionic GECP rate remains a research objective as
specified in `SPEC.md`.

The strongest current structural hypothesis is sign regularity combined with
a quantitative cross-ratio or near-max-volume estimate. Exact minor signs
alone do not control GECP residual decay.

## Not claimed

- A solution of Simons Problem 4.2.
- A GECP rate inferred from an empirical decay curve.
- Diagonal pivot selection for every possible PSD tie-breaking rule.
- The planned explicit dyadic/Chebyshev uniform approximation theorem. Only its
  natural-number counting scaffold and an independently tested numerical
  piecewise interpolant are present.
- A rigorous 128-bit continuous/adaptive pivot census. The committed census is
  finite-grid float64 data.
