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
- For every positive natural `p` and natural scale `s`, an explicit dyadic
  construction approximates `exp(-tω)` on `[0,1] × [0,2ˢ]` with
  `8p(s+1)` separated terms and uniform error at most `2⁻ᵖ`.
- Reflection and the positive fermionic denominator give an explicit
  `16p(s+1)`-term approximation of the fermionic kernel on
  `[0,1] × [-2ˢ,2ˢ]` with the same error. This is the machine-checked
  `O(log Λ log(1/ε))` low-rank theorem, not a GECP theorem.
- Strict minor-sign data alone do not determine a complete-pivot location:
  `signRegularLeft` and `signRegularRight` are exact `2 × 2` rational kernels
  with the same strict sign-regular pattern and no common complete pivot.
- If every selected residual cross satisfies `CrossRatioControl` with factor
  `θ`, the exact GECP update sequence has residual bound `θⁿ` times its initial
  uniform bound. For `0 ≤ θ < 1` this is geometric decay.
- Variable cross-ratio factors give the exact cumulative-product residual
  bound. A product bounded by `2⁻ᵖ` at rank `block * p` gives the corresponding
  dyadic residual rate.
- The positive cutoff corner is a complete pivot on the fermionic cutoff
  rectangle. Its reflected-corner first-update ratio is exactly
  `1 - exp(-2Λ)`, so no factor below one contracts every first step uniformly
  in the cutoff.
- For every positive cutoff, the residual after the first corner pivot is
  nonnegative and is maximized in absolute value at the reflected corner
  `(1,-Λ)`. Thus the first two canonical complete pivots are proved exactly,
  not inferred from the census.
- For `0 < Λ ≤ 1`, the rank-two secant approximation has kernel error at most
  `1/8`; its endpoint cardinal weights are nonnegative and sum to at most one.
  The resulting recursive two-corner GECP residual is at most `1/4`, hence at
  most one half of the initial complete-pivot magnitude. This is the first
  formally proved strict block-contraction case for the fermionic trajectory.
- A residual that is `PivotCrossProductSignCoherent` at its selected pivot cannot grow under an exact complete
  pivot: compatible cross-product signs sharpen the generic factor-two update
  estimate to factor one.
- Its time and frequency derivatives are proved exactly. On
  `[0,1] × [-Λ,Λ]`, the kernel is at most one and the coordinate derivative
  magnitudes are bounded by `Λ` and one, respectively.
- A Lipschitz grid cover yields an explicit continuous supremum bound, and a
  certified upper bound implies an approximate-pivot inequality.
- Uniform kernel error transfers to the Green-function error with the
  `L¹` norm of the spectral density as weight.

## Observed

- All minors enumerated in the 21 exact matrices `q^(ij)` for sizes 2–8 and
  `q ∈ {1/2, 2/3, 3/4}` are nonzero with the predicted sign.
- Exact complete-pivot paths vary with `q` by size eight. This rules out a
  single `q`-independent pivot order for this surrogate family.
- The endpoint-resolved 128-bit finite-grid fermionic census converged in all
  21 cases. At tolerance `1e-10`, sampled ranks grow from 8 at cutoff 1 to 124
  at cutoff `1e6`. Two complete executions produced byte-identical JSONL.
- The fixed two-delta fixture is recovered below `1e-8` using two atoms.
- One 128-bit, 12-pivot fermionic GECP basis at cutoff 8 compresses the
  2,001-node validation frequency grid by a factor of 166.75. For the
  normalized Hubbard-like three-peak density, the held-out kernel error is
  `5.404231354739705e-9` and the Green-function error is
  `4.162832301091157e-10`, or 7.70% of the discrete `L¹` transfer bound.
- The same universal basis applied to a normalized gapped two-band density
  gives Green-function error `2.909588125987739e-10`, or 5.38% of the same
  transfer bound. This demonstrates reuse across qualitatively different
  spectra rather than density-specific basis fitting.
- A four-entry known transition library recovers a synthetic
  quasiparticle/satellite spectrum with four atoms and
  `3.3306690738754696e-16` held-out error. A separate dense blind scan of 1,604
  candidates with deterministic `1e-8`-scale input noise uses eight effective
  atoms and reaches `8.559113029438237e-9` error against the noiseless held-out
  Green function.
- On a 72-site clustered Gaussian covariance, complete-pivot GECP and
  diagonal-pivoted Cholesky select the same 31 landmarks. The selected cross
  has maximum error `7.782740553130552e-7` and relative Frobenius error
  `1.1432201168664416e-7`.
- The canonical synthetic-application JSON is byte-identical across repeated
  executions and has SHA-256
  `f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2`.
- Using blocks of `2(s+1)` pivots for the least `s` with `Λ ≤ 2ˢ`, every
  complete block in the strictest-tolerance 128-bit census reduces the sampled
  residual by less than one half. The largest observed complete-block ratio is
  about `0.07741` at cutoff one; the largest among cutoffs `10` through `10⁶`
  is about `0.002475`.
- Outward-rounded interval residual bounds certify a continuous two-step ratio
  below `0.078` at cutoff one along the returned certified approximate-pivot
  trajectory.

## Conjectured

The target continuous fermionic GECP rate remains a research objective as
specified in `SPEC.md`.

The cutoff-one base case of the `2(s+1)` block hypothesis is now proved for
the actual continuous GECP trajectory. The conjectured extension is block
contraction after `2(s+1)` pivots for every `Λ ≤ 2ˢ`. The formal
`PivotCrossProductSignCoherent` condition
would supply nonexpansiveness between strict contractions and is expected to
follow from a stronger residual sign-regularity theorem. Exact minor signs
alone do not control GECP pivot locations or residual decay.

## Not claimed

- A solution of Simons Problem 4.2.
- A GECP rate inferred from an empirical decay curve.
- Diagonal pivot selection for every possible PSD tie-breaking rule.
- Formalization of the Gimbutas–Marshall–Rokhlin selected-exponential
  Chebyshev construction. The delivered formal theorem instead uses an
  explicit dyadic truncated-Taylor construction with the same asymptotic rank
  scale.
- A continuous-domain GECP convergence conclusion from the 128-bit census. The
  committed dataset is endpoint-resolved finite-grid evidence; adaptive
  interval pivot certificates are validated separately on bounded cases.
- Proof that the fermionic residual sequence satisfies `CrossRatioControl`
  with a cutoff-uniform contraction factor below one.
- Proof that fermionic residuals are `PivotCrossProductSignCoherent` at every
  selected pivot, or that the cutoff-scaled half-contraction extends from the
  proved base regime `0 < Λ ≤ 1` to all dyadic cutoff scales.
- Material-specific validation of the Hubbard-like or gapped fixtures. They
  are stylized synthetic densities, not outputs fitted to experiment, a named
  compound, DMFT, or quantum Monte Carlo.
- Unique recovery of the generating frequencies in the noisy blind scan. Its
  eight atoms are an accurate effective representation of the held-out Green
  function, not an identifiable reconstruction of the four true atoms.
- An optimal sensor-placement theorem for the PSD covariance landmarks.

## Run provenance

Agent model/mode metadata, elapsed time, aggregate token usage, and the
non-billing cost estimate are recorded separately in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They are
operational provenance and are not mathematical or numerical evidence for any
claim in this ledger.
