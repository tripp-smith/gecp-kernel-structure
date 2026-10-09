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
- For every successful dependent GECP run, appending an arbitrary row and
  column to its ordered selected core gives determinant
  `det(core) * finalResidual(row,column)`. Because the selected core is
  nonsingular, `PivotCrossProductSignCoherent` is equivalent to a
  `BorderedMinorSignCoherent` four-determinant product condition on the
  original kernel.
- `StrictSignRegularAtOrder` and `StrictSignRegular` encode prescribed signs
  of ordered kernel minors. All-orders strict sign regularity implies the
  four-bordered-minor condition, including proposed points that duplicate a
  selected row or column, and therefore implies exact selected-residual
  cross-product sign coherence for every successful run.
- Positive column scaling preserves every prescribed ordered-minor sign. The
  fermionic kernel is such a scaling of `exp(-t*omega)`, so strict sign
  regularity transfers order by order and at all orders. Unconditionally, the
  expected signs are proved for both kernels at orders one and two.
- Positive row scaling likewise preserves every ordered-minor signature. For
  every pair of increasing real `n`-tuples, the ordinary Vandermonde
  determinants are positive, and the square exponential Taylor feature block
  with exponents `0, ..., n-1` and weights `1/k!` has strictly positive
  determinant. This proves one positive term of the planned Cauchy--Binet
  expansion, not positivity of the full exponential determinant.
- For every positive strictly increasing real node tuple and strictly
  increasing natural exponent tuple, the generalized Vandermonde determinant
  `det(x_i ^ m_j)` is strictly positive. The proof combines a sparse-polynomial
  positive-root bound from Descartes' rule, nonsingularity, positivity at
  geometric nodes, and a zero-avoiding homotopy. This strict result does not
  include a zero node.
- For nonnegative strictly increasing nodes, the same generalized
  Vandermonde determinant is nonnegative. Positive uniform shifts reduce the
  result to the interior theorem, and determinant continuity closes the
  zero-node boundary.
- Rectangular matrix products over finite ordered index types satisfy the
  ordered Cauchy--Binet expansion `Matrix.det_mul_rect`. Applying it to finite
  exponential Taylor features shows that every generalized-Vandermonde
  summand is nonnegative and the principal summand is positive. Consequently,
  `expTaylorMatrix_det_pos` proves strict positivity at every order whenever
  the term count is at least the matrix order and both node tuples are
  nonnegative and strictly increasing.
- The finite determinants retain the fixed determinant of the principal
  `0, ..., n-1` feature block as a lower bound. Their entries converge to
  `exp(row_i * column_j)`, determinant continuity preserves that lower bound,
  and `expMatrix_det_pos` therefore proves strict positivity of the full
  exponential determinant for nonnegative strictly increasing row and column
  tuples.
- Positive exponential row and column scalings remove arbitrary tuple offsets,
  so `exp(x*y)` has positive determinants for every pair of strictly increasing
  real tuples. Reversing the negated frequency tuple contributes exactly
  `(-1)^(n.choose 2)`, proving `expKernel_strictSignRegular` and, by positive
  fermionic denominator scaling, `fermionicKernel_strictSignRegular` at every
  order. Consequently every selected fermionic residual cross is
  `PivotCrossProductSignCoherent` and exact complete-pivot updates are
  nonexpansive.
- Successful dependent runs compose with `Run.append`. On restricted row and
  column domains, `GECP.CompletePivotOn` records membership and exact residual
  maximization, while `residualUpdate_le_of_signCoherentOn` needs bounds only
  at the four points in the selected cross. Therefore
  `strictSignRegular_gecp_error_nonincreasing` propagates an initial domain
  bound through every realized complete-pivot residual sequence. In
  particular, every exact fermionic sequence on
  `[0,1] × [-Λ,Λ]` remains bounded by the initial cutoff-corner pivot.
- Along every realized strictly sign-regular exact complete-pivot sequence,
  the selected pivot magnitudes are antitone. If a successful length-`n` run
  records the first `n` sequence pivots, then for every domain point
  `(x,y)`, the rank-`n` residual satisfies
  `|R_n(x,y)|^n ≤ |det(K[X_n,Y_n])|`. The proof identifies the absolute
  selected-core determinant with the product of the absolute pivots. This is
  a geometric-mean reduction of residual decay to determinant decay, not a
  determinant decay theorem.
- If an `n x n` matrix is entrywise within `epsilon` of a matrix factoring
  through `r` coordinates, its determinant expansion has no term with fewer
  than `n-r` error columns. Consequently, for `r ≤ n`, `0 ≤ epsilon ≤ 1`,
  and approximant entries bounded by `C ≥ 1`,
  `|det A| ≤ 2^n n! epsilon^(n-r) C^r`. The proof is an exact column-choice
  expansion and finite-dimensional factorization, not a singular-value or
  floating-point argument.
- For arbitrary fermionic sample points in
  `[0,1] x [-2^s,2^s]`, the explicit rank `r = 16p(s+1)` separated
  approximation gives
  `|det K_n| ≤ 2^n n! 2^(-p(n-r)) 2^r` whenever `r ≤ n`. The theorem
  converts low-rank existence into sampled determinant decay; the next two
  results record its GECP composition and explicit parameter choice.
- The sampled determinant estimate now composes with every realized
  physical-domain exact complete-pivot run. At
  `n = 32(2m+1)(s+1)`, the rank-`n` residual obeys
  `|R_n(x,y)| <= 2n 2^-m`. The proof chooses an odd approximation order so the
  determinant exponent is exactly divisible by `n`, then uses
  `n! <= n^n` and extracts the positive `n`th root.
- The explicit choice `m = 16(s+1)` satisfies the required arithmetic
  inequality and proves `|R_n(x,y)| <= 1/2` by
  `n = 32(32(s+1)+1)(s+1)`. This is a rigorous quadratic-in-log-cutoff strict
  contraction bound, not the conjectured linear `2(s+1)` block.
- The sharper choice `m = 16 + 2 ceil(log_2(s+1))` also satisfies the exact
  arithmetic condition. It proves `|R_n(x,y)| <= 1/2` by
  `n = 32(33 + 4 ceil(log_2(s+1)))(s+1)`, improving the certified block to
  `O((s+1) log(s+1))` without changing the complete-pivot hypotheses.
- More generally, for every dyadic accuracy order `q`, the choice
  `m = 16 + 2 ceil(log_2(s+1)) + 2q` proves
  `|R_n(x,y)| <= 2^-q` by
  `n = 32(33 + 4 ceil(log_2(s+1)) + 4q)(s+1)`. Thus the determinant route
  yields rank `O((s+1)(log(s+1)+q))`, or
  `O(log(1+Λ)(log log(1+Λ)+log(1/epsilon)))` after choosing a dyadic target.
- Bounded entries alone cannot improve the determinant prefactor to `C^n`
  for any universal constant base `C`. For every natural `C`, an exact
  Sylvester--Hadamard matrix with unit-modulus entries satisfies
  `C^n < |det A|`; equivalently, this family has determinant root `sqrt(n)`.
  Thus the remaining GECP startup loss cannot be removed by a generic
  entrywise determinant theorem and needs fermionic or residual-selected
  structure.
- The fermionic kernel has the exact multiplicative cross ratio that Phase Y's
  arbitrary matrices lack. An update at `(t0, omega0)` factors as
  `K(t,omega) * (1-exp((t-t0)*(omega-omega0)))`. On the oppositely ordered
  corner region `t0 <= t` and `omega <= omega0`, this residual is nonnegative,
  is at most the cross area `(t-t0)*(omega0-omega)`, and is at most one half
  whenever that area is at most `log 2`. The cross-area formulation covers
  the reflected southeast orientation as well.
- After pivots at `(0,omegaHigh)` and `(1,omegaLow)` for any
  `omegaLow < omegaHigh`, the exact nested fermionic residual is
  `(z^t-S_t(z))/(1+z)`, where `z=exp(-omega)` and `S_t` is the linear secant
  of the concave power function between the transformed band endpoints. The
  residual is nonnegative on the whole band and obeys the two-sided tent bound
  `min(t*(omegaHigh-omega), (1-t)*(omega-omegaLow))`. In particular, it
  vanishes on both time edges and both pivot-frequency edges without a
  symmetric-band or small-cutoff assumption.
- The two symmetric cutoff-corner pivots do not yield a half contraction
  uniformly in cutoff. At `Lambda = log 13824`, `t = 2/3`, and
  `omega = -log(125/64)`, the exact two-corner residual is
  `4495348/8973531`, while half the initial pivot is
  `(1/2)*(13824/13825)`. The residual exceeds that target by the exact positive
  margin `222676/224338275`. The witness lies inside the physical cutoff
  rectangle, so the universal two-pivot half-contraction statement is formally
  false even though the first two symmetric corners are complete pivots.
- The obstruction still forces the next complete pivot inward. After the two
  symmetric cutoff-corner pivots, the residual is at most `1/4` whenever
  `|omega| >= Lambda/2`: on the positive side it is dominated by
  `x*(1-x)` for `x=exp(-t*Lambda/2)`, and the negative side follows by exact
  reflection. At `(t,omega)=(1/2,0)` the residual is exactly
  `1/2 - 1/(2*cosh(Lambda/2))`, which is strictly greater than `1/4` when
  `2*log 4 <= Lambda`. Consequently every actual third complete pivot has
  `|omega| < Lambda/2` in that regime. This is the first proved dyadic
  frequency-localization step for the canonical fermionic GECP trajectory.
- This outer-half exclusion persists through every later exact complete-pivot
  update. Sign coherence makes the two Schur cross products have a common
  sign; the factors evaluated in the outer strip remain bounded by `1/4`,
  while complete-pivot maximality controls the other factor by the pivot.
  Consequently the updated outer-strip residual is still at most `1/4` in
  absolute value. Iterating this mixed local/global estimate proves that every
  later selected pivot with magnitude greater than `1/4` has
  `|omega| < Lambda/2`.
- Low-rank determinant decay survives finitely many arbitrary columns outside
  the approximation band. Representing an exceptional set exactly augments
  the factorization rank by precisely its cardinality; every remaining column
  beyond that augmented rank still contributes the original approximation-
  error power. For fermionic sample matrices, cutoff hypotheses are therefore
  required only on nonexceptional frequencies. In particular, the two cutoff
  corners cost two factor coordinates while every Phase AD-localized
  continuation column can use the smaller central-band approximation.
- The exceptional-column estimate now composes with the actual selected core
  of a realized fermionic GECP prefix. Complete pivoting and residual
  evaluation may use a larger dyadic rectangle, while the determinant bound
  uses an independent smaller approximation scale on every nonexceptional
  selected column. The resulting residual-power exponent charges exactly
  `16*p*(s+1) + exceptional.card` effective coordinates.
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
- Every selected cross in those 21 exact geometric-surrogate runs is sign
  coherent. The diagnostic returns the first negative cross-product witness;
  on the generic matrix `[[-4,-1],[-3,4]]` it exactly reports step zero,
  pivot `(0,0)`, point `(1,1)`, and product `-48`, so the criterion is not
  automatic for arbitrary matrices.
- For the `q = 2/3` geometric surrogate at sizes two through four, exhaustive
  exact enumeration of every row and column permutation confirms that each
  determinant sign is the ordered-minor sign times the two permutation signs.
- Exact rational Taylor-feature matrices at sizes one through six have
  determinant equal to the product of the row Vandermonde, reciprocal
  factorial weights, and column Vandermonde, and that product is positive.
- Exact rational generalized Vandermonde matrices at sizes one through five,
  over every increasing exponent tuple selected from `0, ..., n+2`, have
  positive determinant.
- With the least node set to zero, the same exact family is nonnegative and is
  strictly positive exactly when its least exponent is zero.
- Exact rational finite Taylor matrices at sizes one through five and term
  counts from the matrix order through three additional terms have strictly
  positive determinants, including node tuples that start at zero.
- At 100 decimal digits, full exponential matrices at sizes one through five
  have positive determinants and retain the exact principal Vandermonde-
  factorial lower bound.
- At 100 decimal digits, ordered `exp(-t*omega)` determinants through size six
  are nonzero with sign `(-1)^(n.choose 2)` for tuples spanning negative and
  positive values.
- A 70-decimal finite-grid scan through 24 fermionic pivots at cutoffs
  `1,2,10,100` found no substantive sign-coherence failure. Tiny normalized
  negatives between about `1e-89` and `1e-131` occurred only on rows or
  columns that are algebraically zero after selection and are treated as
  cancellation evidence, not a proof of nonnegativity.
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

The cutoff-one base case of the `2(s+1)` block hypothesis, all-orders sign
regularity, complete-run nonexpansiveness, the determinant geometric-mean
reduction, sampled-core determinant decay, and their exact composition are now
proved. The composition now yields arbitrary dyadic accuracy with rank linear
in the requested accuracy order after an additive logarithmic scale startup.
It remains open whether that startup factor and the constants can be removed,
or the argument localized enough to prove the conjectured linear block and the
cutoff-uniform form of Conjecture G1.
The generic bounded-entry route is now formally ruled out; total positivity,
divided differences, or GECP-selected residual structure remain viable because
the obstruction family does not satisfy the fermionic hypotheses.
The exact corner-area law and its two-corner secant tent now supply a local,
scale-invariant mechanism, but the exact Phase AB obstruction rules out using
only the two outer corners at every scale. Phase AC forces the third complete
pivot into the central half-band, and Phase AD proves that the original
outer-half quarter exclusion survives every later update. It remains
conjectural whether the Phase AE exceptional-column determinant bound can be
instantiated with exactly the two recursive cutoff-corner indices and Phase
AD's localized continuation to restart the argument on the central
restriction and iterate shrinking bands using only linearly many pivots.

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
- Proof that the cutoff-scaled half-contraction extends from the proved base
  regime `0 < Λ ≤ 1` to all dyadic cutoff scales.
- A cutoff-uniform strict contraction factor or proof of the dyadic block
  localization required by Conjecture G1. All-orders sign regularity gives
  nonexpansiveness, not strict decay.
- A proved cutoff-uniform GECP rate obtained by combining the new sampled-core
  determinant bound with the formal geometric-mean inequality at the
  conjectured linear-in-`s` rank scale. The delivered composition instead has
  an explicit arbitrary-accuracy rate with an additive
  `(s+1) log(s+1)` startup cost.
- A counterexample to a fermionic-specific or totally-positive determinant
  estimate. The Sylvester--Hadamard obstruction concerns arbitrary
  bounded-entry matrices and is deliberately outside those structural
  classes.
- A global GECP contraction obtained by applying the corner-area theorem to
  later residuals. The proved exponential cross ratio belongs to the original
  fermionic kernel and is not asserted to be preserved verbatim by Schur
  updates.
- Identification of the two prescribed corners with the first two complete
  pivots on every asymmetric subband, or control of the actual third complete
  pivot by the proved tent. The Phase AA theorem concerns the exact residual
  after those prescribed nonzero pivots.
- A counterexample to Conjecture G1 or to cutoff-dependent multi-pivot block
  contraction. Phase AB rules out only the stronger shortcut asserting that
  the first two symmetric corner pivots halve the residual for every cutoff.
- A shrinking-band iteration of the central-half localization theorem. Phase
  AD preserves the fixed outer-half quarter bound through all later Schur
  updates, while Phases AE–AF permit exact exceptional columns in the actual
  GECP selected-core estimate. The first two recursive indices have not yet
  been packaged as the exceptional set, and no new relative central-half
  bound follows yet.
- Material-specific validation of the Hubbard-like or gapped fixtures. They
  are stylized synthetic densities, not outputs fitted to experiment, a named
  compound, DMFT, or quantum Monte Carlo.
- Unique recovery of the generating frequencies in the noisy blind scan. Its
  eight atoms are an accurate effective representation of the held-out Green
  function, not an identifiable reconstruction of the four true atoms.
- An optimal sensor-placement theorem for the PSD covariance landmarks.

## Verification evidence

- On 2026-10-07 the Phase J root command passed: Lean build, public axiom
  audit, placeholder rejection, Ruff, formatting, mypy, 48 Python tests, and
  dependency audit with no known vulnerabilities; the same change passed CI
  and merged as PR #14. Phase H separately verified the root entry point in
  maintained and fresh frozen environments.
- On 2026-10-07 the Phase K root command passed the same gate with 49 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #15.
- On 2026-10-08 the Phase L root command passed the same gate with 50 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #16.
- On 2026-10-08 the Phase M root command passed the same gate with 51 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #17.
- On 2026-10-08 the Phase N root command passed the same gate with 52 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #18.
- On 2026-10-08 the Phase O root command passed the same gate with 53 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #19.
- On 2026-10-08 the Phase P root command passed the same gate with 54 Python
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #20.
- On 2026-10-08 the Phase Q root command passed the same gate with 54 Python
  tests and no known dependency vulnerabilities; its four new structural
  results use only the permitted Lean axioms. The change passed CI and merged
  as PR #21.
- On 2026-10-08 the Phase T root command passed the same gate with 54 Python
  tests and no known dependency vulnerabilities. Its four new structural
  results use only the permitted Lean axioms. The change passed CI and merged
  as PR #22.
- On 2026-10-08 the Phase U root command passed the same gate with 55 Python
  tests and no known dependency vulnerabilities. Its eleven audited public
  results use only the permitted Lean axioms; the exact rational low-rank
  perturbation census and Lean rank-one `3 x 3` check also pass. The change
  passed CI and merged as PR #23.
- On 2026-10-08 the Phase V root command passed the same gate with 56 Python
  tests and no known dependency vulnerabilities. Its ten new public results
  use only the permitted Lean axioms; the independent exact integer scale
  regression also passes. The change passed CI and merged as PR #24.
- On 2026-10-09 the Phase W root command passed the same gate with 57 Python
  tests and no known dependency vulnerabilities. Its two new public results
  use only the permitted Lean axioms; the exact integer regression covers
  ordinary scales and dyadic transition boundaries through `2^64`. The change
  passed CI and merged as PR #25.
- On 2026-10-09 the Phase X root command passed the same gate with 58 Python
  tests and no known dependency vulnerabilities. Its three new public results
  use only the permitted Lean axioms; the exact integer regression covers 65
  accuracy orders and dyadic transition boundaries through `2^64`. The change
  passed CI and merged as PR #26.
- On 2026-10-09 the Phase Y root command passed the same gate with 59 Python
  tests and no known dependency vulnerabilities. Its four new public results
  use only the permitted Lean axioms; the exact Sylvester regression covers
  orders two through 64 and obstruction witnesses for bases zero through 256.
  The change passed CI and merged as PR #27.
- On 2026-10-09 the Phase Z root command passed the same gate with 60 Python
  tests and no known dependency vulnerabilities. Its eight new public results
  use only the permitted Lean axioms; the independent 100-decimal regression
  checks four arbitrary cross-ratio factorizations and 4,590 admissible
  oriented corner tiles. The change passed CI and merged as PR #28.
- On 2026-10-09 the Phase AA root command passed the same gate with 61 Python
  tests and no known dependency vulnerabilities. Its eleven new public
  results use only the permitted Lean axioms; an independent 100-decimal
  regression checks 1,156 points across asymmetric and reflected bands for
  the exact nested-update identity, sign, endpoint zeros, and two-sided tent
  bound. The change passed CI and merged as PR #29.
- On 2026-10-09 the Phase AB root command passed the same gate with 62 Python
  tests and no known dependency vulnerabilities. Its six new public results
  use only the permitted Lean axioms; exact `Fraction` arithmetic verifies the
  witness and sweeps endpoint bases 24 through 64, while a direct 100-decimal
  nested update reproduces the formal residual and margin. The change passed
  CI and merged as PR #30.
- On 2026-10-09 the Phase AC root command passed the same gate with 63 Python
  tests and no known dependency vulnerabilities. Its ten new public results
  use only the permitted Lean axioms; a 100-decimal regression checks the
  exact center formula, both outer-half quarter bounds, and sampled-maximizer
  localization across five cutoffs from the threshold through 16. The change
  passed CI and merged as PR #31.
- On 2026-10-09 the Phase AD root command passed the same gate with 64 Python
  tests and no known dependency vulnerabilities. Its four new public results
  use only the permitted Lean axioms; exact rational sign-orientation cases
  check the mixed-bound cancellation, and an 80-decimal regression follows
  five further sampled complete-pivot updates at cutoffs 4, 8, and 16. The
  change passed CI and merged as PR #32.
- On 2026-10-09 the Phase AE root command passed the same gate with 65 Python
  tests and no known dependency vulnerabilities. Its four new public results
  use only the permitted Lean axioms; an exact `Fraction` regression checks
  determinant expansion, the effective-rank vanishing threshold, and the
  coarse error power with zero and two exceptional columns. The change passed
  CI and merged as PR #33.
- On 2026-10-09 the Phase AF root command passed the same gate with 66 Python
  tests and no known dependency vulnerabilities. Its two new public results
  use only the permitted Lean axioms; an exact `Fraction` regression verifies
  a full-rank GECP selected core, its pivot product, and the two-exception
  effective-rank determinant bound. The change passed CI and merged as PR #34.
- The verification refresh changes no proved, observed, conjectured, or
  not-claimed mathematical statement above.

## Run provenance

Agent model/mode metadata, elapsed time, aggregate token usage, and the
non-billing cost estimate are recorded separately in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They are
operational provenance and are not mathematical or numerical evidence for any
claim in this ledger.
