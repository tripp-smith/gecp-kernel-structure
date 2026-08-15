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

## Conjectured

The target continuous fermionic GECP rate remains a research objective as
specified in `SPEC.md`.

The strongest current structural hypothesis is the formal
`CrossRatioControl`, potentially derived from sign regularity plus a
near-max-volume estimate. Exact minor signs alone do not control GECP pivot
locations or residual decay.

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
