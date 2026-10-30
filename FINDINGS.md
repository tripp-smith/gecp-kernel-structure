# Findings and claim ledger

## Current original-SPEC audit (Phase DQ verified; merged PR #123)

The source-level [final audit](FINAL_AUDIT.md) supports the exact original-kernel
finite-run theorem: for every real `Λ≥1`, every actual finite complete GECP run
on `[0,1]×[-Λ,Λ]`, and every physical point,
`|R_k| ≤ 1024 exp(-((log 2)^2/256) k/log(1+Λ))`.
The selected-core cross and actual Green-integral error obey the corresponding
bound. Real-parameter separated approximation, arbitrary-domain strict-SPD
diagonal necessity/Cholesky equivalence, and original-kernel `2*L*h`
fill-distance control are also proved. No convergence conclusion is supplied
as a premise. Literal all-positive-cutoff `log(1+Λ)` wording is refuted; the
all-positive `log(2+Λ)` variant is separate.

**The required original-SPEC scope is delivered and verified at its stated boundaries.**
The explicit-item evidence in [FINAL_AUDIT.md](FINAL_AUDIT.md) inspects every
original section, theorem hypotheses, diagnostics, tests, six comparison
methods and proposed imports. The required scope is accounted for without
promoting optional proof routes or approximate-pivot extensions into proved
results. DQ also reproduces and repairs genuine conditional pivot-certificate
rounding and dyadic-parameter coverage defects. Its new theorem consumer uses
only the public import and checks 22 non-tautological existing interfaces.
No new formal result is asserted. Exact DQ code `f2ffd5e` passes the 757-test
root, 3,634 Lean jobs, 876 allowed axiom reports, 22 public-import examples,
style/type/security gates and push/PR CI `40050118292`/`40050129457`.
Clean rebuilt archives and both outside-checkout installed checks pass.
Final branch `0a6d363` repeats the root (757 tests in 149.97 seconds), clean
rebuilt archives and both installed environments; CI `40050605960`/`40050610230`
pass all four checks. PR #123 merges at 2026-10-24 01:41:41 UTC as
`d9978607b88a45d98d917f44318bcaf99b4c3239`. Exact main `e8c35ed` passed the root
(757 tests in 150.31 seconds), all proof/type/style/security gates,
[CI 40051138962](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40051138962),
archive inspection and both installed environments. The original goal was
marked complete after these gates. The completion estimate is 100% of required scope,
not a claim that every optional research question has been resolved positively.

DO reopened R25
because the fixed order-24 census and three order-32 reruns did not supply the
required increasingly dense residual-adaptive experiment. DP now records all
seven cutoffs and 28 tolerance prefixes passing on genuinely refined nested
grids, with two fresh clean-source runs producing identical canonical bytes.
These are observed finite-sample results, not a new theorem. The normalized
plot has been rendered, visually inspected and repeated byte-for-byte.
Independent post-run MP192 checking passes all 21 selected-core checks and
28 full final crosses. Final branch head
`3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes the exact local root,
source/wheel and both outside-checkout install checks, plus green push/PR CI.
[PR #122](https://github.com/tripp-smith/gecp-kernel-structure/pull/122) merged
on 2026-10-23 at 20:20:48 UTC as
`3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`. Exact pushed-main closure
metadata/gates subsequently pass at main `5a5fef0`, including CI `40018905634`.
DQ supplies the subsequent original-scope audit and final repairs; its
completion decision follows original requirements and actual gates, not
the historical forty-group count or earlier 99% effort estimate.

DO's numerical repairs and historical verification remain valid: `f0d082d`
passes the 569-test root, 3,634 Lean jobs, all 876 axiom reports, exact-head
push/PR CI and isolated wheel checks with SciPy 1.18 and minimum 1.15.0.
PR #121 merged on 2026-10-23. Historical artifacts retain their original source
hashes; DP adds evidence without overwriting those fixed-grid datasets.
Certificates concern the stored rounded residual model, not the ideal exact
run or a separately reconstructed cross. No arbitrary-η convergence theorem,
continuous numerical census certificate, optimal sparse support, or literature
priority is claimed.

The phase entries below are a chronological development record. Their
“remaining”/“open” language describes the named phase's delivery, not current
status; later proved results and this audit supersede earlier route limits.

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
- The canonical two-corner prefix is now represented as an actual dependent
  elimination run. Its recursive selected-index type has a proved two-element
  exceptional set consisting exactly of the `(0,Lambda)` and `(1,-Lambda)`
  steps. Every other selected column belongs to the continuation, so a
  continuation complete on the smaller band gives the closed determinant
  exponent `16*p*(s+1)+2` without positional list assumptions.
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
- The historical endpoint-resolved 128-bit finite-grid fermionic census
  converged in all 21 cases. At tolerance `1e-10`, sampled ranks grow from 8
  at cutoff 1 to 124 at cutoff `1e6`. Two complete executions produced
  byte-identical JSONL.
- The DP residual-adaptive MP128 census records seven successful nested-grid
  trajectories and all 28 requested tolerance prefixes passing separate
  training/probe/surplus and final MP192 sampled checks. The `1e-12` chosen
  ranks are `9,17,36,64,97,124,151`; two clean-source runs produce identical
  JSON bytes. Independent MP192 checking passes 21 selected-core checks and
  all 28 full final crosses within their requested tolerances. DQ's source
  closure audit is verified; no continuous supremum bound is inferred.
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
  `3.3306690738754696e-16` held-out error. A separate dense scan of 1,604
  candidates (including explicitly added generating frequencies, so not
  support-blind) with deterministic `1e-8`-scale input noise uses eight effective
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
conjectural whether Phase AD's pointwise above-quarter localization can be
lifted to `Run.CompleteOn` for a finite continuation and then iterated to
restart the argument on the central restriction using only linearly many
pivots.

## Not claimed

- A solution of Simons Problem 4.2.
- A GECP rate inferred from an empirical decay curve.
- Diagonal pivot selection for every possible PSD tie-breaking rule.
- Formalization of the Gimbutas–Marshall–Rokhlin selected-exponential
  Chebyshev construction. The delivered formal theorem instead uses an
  explicit dyadic truncated-Taylor construction with the same asymptotic rank
  scale.
- A continuous-domain GECP convergence conclusion from either 128-bit census.
  The historical fixed-grid and DP residual-adaptive datasets are finite-grid
  evidence, including separately labeled finite validation samples. Adaptive
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
  updates, while Phases AE–AG permit the exact two corner columns in the
  smaller-band selected-core estimate. A finite continuation has not yet been
  proved `CompleteOn` that smaller band from the above-quarter condition, and
  no new relative central-half bound follows yet.
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
- On 2026-10-09 the Phase AG root command passed the same gate with 67 Python
  tests and no known dependency vulnerabilities. Its five new public results
  use only the permitted Lean axioms; an exact rational prescribed-corner run
  verifies the two exceptional indices, central continuation columns,
  selected-core determinant, and complete pivot product. The change passed CI
  and merged as PR #35.
- On 2026-10-09 the Phase AH root command passed the same gate with 68 Python
  tests and no known dependency vulnerabilities. Its trajectory constructor
  and central-prefix results use only the permitted Lean axioms; an exact
  rational recurrence independently reconstructs three above-quarter complete
  pivots, their central columns, and the final residual. The change passed CI
  and merged as PR #36.
- On 2026-10-09 the Phase AI root command passed the same gate with 69 Python
  tests and no known dependency vulnerabilities. Its determinant-sandwich and
  explicit stopping results use only the permitted Lean axioms; exact integer
  regression verifies the effective-rank, binary-size, and exponent budgets
  across ordinary and power-of-two-adjacent scales. The change passed CI and
  merged as PR #37.
- On 2026-10-10 the Phase AJ root command passed the same gate with 70 Python
  tests and no known dependency vulnerabilities. Its arbitrary-border
  determinant factorization, orientation-aware strict-sign transfer, and
  fermionic residual-minor nonsingularity results use only the permitted Lean
  axioms; exact `Fraction` arithmetic independently checks border sizes zero
  through three. The change passed CI and merged as PR #38.
- On 2026-10-10 the Phase AK root command passed the same gate with 71 Python
  tests and no known dependency vulnerabilities. Its balanced-sign theorem
  turns residual cross-product coherence into an explicit magnitude-one
  diagonal gauge, proves positivity on every fresh fermionic entry, preserves
  all entry magnitudes, and records exact finite-minor scaling. An exact
  `Fraction` regression checks every unselected anchor and entry after three
  geometric-surrogate pivots. The change passed CI and merged as PR #39.
- On 2026-10-10 the Phase AL root command passed the same gate with 72 Python
  tests and no known dependency vulnerabilities. Fresh residual minors at
  every finite order now form a cross-product-sign-coherent compound kernel,
  and an anchor-based magnitude-one gauge makes every compound entry positive.
  Exact `Fraction` regressions exhaust every anchor and entry at compound
  orders one through three after three geometric-surrogate pivots. The change
  passed CI and merged as PR #40.
- On 2026-10-10 the Phase AM root command passed the same gate with 73 Python
  tests and no known dependency vulnerabilities. The exact selected-core
  Desnanot--Jacobi identity is now formalized for two arbitrary border rows and
  columns, and strict sign regularity turns it into a factor-one maximum bound
  rather than the generic two-term triangle bound. Exact `Fraction` arithmetic
  exhausts all 256 ordered fresh two-border choices after a three-pivot
  geometric-surrogate run. This is a verified subtraction-sensitive gain, but
  it propagates one-border magnitude estimates and does not by itself generate
  determinant decay, prove Conjecture G1, or solve Problem 4.2. The change
  passed CI and merged as PR #41.
- On 2026-10-10 the Phase AN root command passed the same gate with 75 Python
  tests and no known dependency vulnerabilities. Every finite complete
  continuation of a strictly sign-regular run now inherits recursive selected-
  cross coherence, keeps all pivots inside its initial uniform envelope `B`,
  and has selected-core determinant at most `B^m` with no factorial. Exact
  `Fraction` arithmetic checks every suffix of an eight-point geometric run.
  A separate exact coefficient `6243374306 / 4938550965 > 1` rules out
  automatic factor-one core dominance. Applying the factorial-free mechanism
  to the dyadic approximation still requires sign-regularity and mixed-column
  compatibility for its remainder; Conjecture G1 and Problem 4.2 remain open.
  The change passed CI and merged as PR #42.
- Phase AO proves that the current hard-masked dyadic Taylor remainder is not
  strictly sign regular even at order two. At `p = 1`, `s = 2`, rows
  `(1/2, 3/5)`, and columns `(2, 5/2)`, its ordered determinant is strictly
  positive, whereas `expKernelSignature 2 = -1` requires it to be negative.
  Lean derives the four entries from the actual dyadic evaluator and proves
  the sign with rational exponential bounds. A 100-decimal independent check
  places the determinant in `(4.96e-6, 4.97e-6)` and a 36-point neighborhood
  sweep shows the same reversal. Thus the Phase AN factorial-free mechanism
  cannot use this masked remainder directly; smooth overlaps, unmasked local
  tails, or a different compatible decomposition remain viable. Conjecture G1
  and Problem 4.2 remain open. The root command passes with 76 tests and no
  known dependency vulnerabilities. The change passed CI and merged as PR #43.
- Phase AP identifies the exact structure that survives Phase AO's mask
  obstruction. The unmasked eighth-order tail is
  `x^8 / 7!` times the strictly positive beta-weighted Laplace moment
  `integral_0^1 (1-u)^7 exp(-ux) du`. Consequently, every order-two tail
  minor is a positive monomial row/column scaling of the corresponding beta-
  moment minor: no sign is lost in normalization. At 100 decimal digits, the
  integral, convergent series tail, direct subtraction, and repository
  implementation agree throughout `0.001 <= x <= 2`; all 315 ordered minors
  on the independent positive grid have the expected negative sign, with
  normalized margin above `4.3e-6`. This is evidence and an exact reduction,
  not yet a proof of the beta-moment sign inequality or all-orders strict sign
  regularity. The root command passes with 77 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #44.
- Phase AQ proves the first analytic inequality needed for the beta-moment
  order-two sign. The exact unweighted gap is `1/120`, and monotonicity of the
  exponential tilt around its single sign change at `u = 1/2` gives
  `2 J_2(x) < J_1(x)` for every `x >= 0`. On the local approximation regime
  `0 <= x <= 2`, Lean then proves
  `x (J_2(x) M(x) - J_1(x)^2) < J_1(x) M(x)`, the strict variance-numerator
  inequality governing logarithmic elasticity. Exact `Fraction` arithmetic
  independently recomputes `1/120`; a 100-decimal endpoint/logarithmic/dense
  sweep keeps the moment gap above `0.0062` and the variance margin above
  `0.0008`. Differentiation and the final ordered-minor theorem remain to be
  formalized. The root command passes with 78 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #45.
- Phase AR supplies the analytic hinge left open by Phase AQ. A single
  dominated-parametric-integral lemma proves
  `M'(x) = -J_1(x)` and `J_1'(x) = -J_2(x)` for every real `x`. The exact
  quotient rule then identifies the derivative of
  `E(x) = -x J_1(x) / M(x)` with the Phase AQ variance numerator divided by
  `M(x)^2`, so `eighthOrderBetaElasticity_strictAntiOn` proves that `E` is
  strictly decreasing on `0 <= x <= 2`. At 100 decimal digits, independently
  differentiated integrals agree at endpoint and logarithmic test points; a
  dense local sweep keeps `-E'(x)` above `0.078`. The two-variable ordered-
  product transfer and final beta-moment minor theorem remain separate tasks.
  The root command passes with 79 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #46.
- Phase AS completes the local order-two sign chain for the unmasked Taylor
  tail. For `0 < t₀ < t₁`, the exact derivative of
  `M(t₀ ω) / M(t₁ ω)` is positive on positive intervals with
  `t₁ ω <= 2`, because Phase AR's decreasing elasticity cross-multiplies to
  precisely the quotient-rule numerator. Therefore every ordered normalized
  beta-moment minor is strictly negative, and Phase AP's positive monomial
  factor gives the same sign for the eighth-order tail. A 100-decimal
  independent check validates the ratio derivative on asymmetric scales and
  every minor of a product-bounded grid. This is local order-two sign
  regularity, not yet an all-orders or global smooth-remainder theorem. The
  root command passes with 80 tests and no known dependency vulnerabilities.
  The change passed CI and merged as PR #47.
- Phase AT uncovers a different route to all orders. The exact substitution
  `v = ωu` rewrites the normalized beta moment as
  `ω^8 M(tω) = integral_0^ω (ω-v)^7 exp(-tv) dv`, and the tail becomes the
  positive row scaling `t^8 / 7!` of this integral. Lean also extends the
  integral to any shared cutoff using the nonnegative truncated-power kernel
  `max (ω-v) 0 ^ 7`, putting every column on the same integration domain.
  This removes Phase AP's product-argument obstacle to a basic composition
  formula. Numerically, the expected signature persists through orders two to
  five on a product-bounded grid. Total positivity of the truncated-power
  factor and continuous Cauchy--Binet strictness remain to be formalized. The
  root command passes with 81 tests and no known dependency vulnerabilities.
  The change passed CI and merged as PR #48.
- Phase AU isolates the finite strictness mechanism required by that route.
  The signed ordered Cauchy--Binet determinant is nonnegative when every
  factor-minor product has the common orientation, every selected summand is
  an exact lower bound for the signed product determinant, and one positive
  summand makes the determinant strictly signed. An exact `2 x 3` by `3 x 2`
  regression independently obtains summands `1, 4, 1` and determinant `6`.
  This closes finite algebra only: truncated-power total nonnegativity and the
  passage to integrals remain open. The root command passes with 82 tests and
  no known dependency vulnerabilities. The change passed CI and merged as PR
  #49.
- Phase AV begins the missing factor theorem without erasing its support
  geometry. New total-nonnegativity predicates separate weak minor signs from
  strict sign regularity. The seventh-degree truncated-power kernel is proved
  totally nonnegative at orders one and two, and its ordered `2 x 2` minor is
  positive exactly when both paired row knots lie below their corresponding
  column knots. Exact arithmetic checks all 784 ordered pair configurations
  on `[-3,4]`, including boundary zeros and the full-support determinant
  `14197`. The all-orders theorem and continuous composition remain open. The
  root command passes with 83 tests and no known dependency vulnerabilities.
- Phase AW turns the classical truncated-power clue into an exact internal
  construction. `normalizedTruncatedPowerKernel_succ_integral` raises the
  degree by integration on the natural interval, while
  `normalizedTruncatedPowerKernel_succ_integral_commonDomain` inserts the
  degree-zero strict order kernel and proves the same identity on any shared
  containing interval. The existing seventh-degree kernel is exactly `7!`
  times normalized degree seven. Exact rational polynomial integration checks
  degrees zero through seven, including coincident endpoints. This establishes
  the convolution ladder, not total nonnegativity or continuous
  Cauchy--Binet. The root command passes with 84 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #51.
- Phase AX proves the all-orders base of that ladder. Every ordered square
  minor of `strictOrderKernel` is exactly zero or one: the first column is
  either zero, creates duplicate all-one rows, or has a single leading one and
  reduces with positive sign to the tail minor. Consequently both the named
  kernel and `normalizedTruncatedPowerKernel 0` are totally nonnegative at
  every order. Exact arithmetic checks all 3,381 ordered row/column tuple
  pairs on `[-3,3]` through order five. Higher-degree total nonnegativity and
  continuous Cauchy--Binet remain open. The root command passes with 85 tests
  and no known dependency vulnerabilities. The change passed CI and merged as
  PR #52.
- Phase AY formalizes the analytic atom for iterating that result.
  `Matrix.det_updateRow_intervalIntegral` proves that integrating one row
  commutes exactly with the determinant because the fixed-row determinant is
  a continuous linear functional. Its nonnegative and strict corollaries
  transfer pointwise signs through the integral. An exact affine-row example
  independently gives determinant integral `6` on both sides. Multi-row
  iteration and higher-degree total nonnegativity remain open. The root
  command passes with 86 tests and no known dependency vulnerabilities. The
  change passed CI and merged as PR #53.
- Phase AZ iterates that atom over every row. A finite-set induction freezes
  unintegrated rows, integrates one selected row, and repeats; no product
  measure is needed. Total nonnegativity also controls weakly ordered samples,
  because any repeated row or column makes the determinant zero. Consequently,
  integrating a totally nonnegative kernel over ordered disjoint row intervals
  preserves the minor's nonnegative sign. The adjacent-row reduction needed to
  turn common-upper truncated-power integrals into these disjoint intervals,
  and hence higher-degree total nonnegativity, remain open.
  Exact arithmetic checks the positive staircase and a collapsed-interval zero
  boundary. The root command passes with 87 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #54.
- Phase BA closes the all-orders truncated-power factor theorem. Reversing row
  order turns subtraction of each next row into mathlib's existing
  predecessor-row determinant invariance; interval additivity then converts
  common-upper tails into ordered adjacent intervals. This proves that total
  nonnegativity propagates from normalized degree `d` to degree `d+1`.
  Starting from the strict order kernel, every normalized degree is therefore
  totally nonnegative, and positive `7!` scaling gives the repository's
  seventh-degree factor at every order. Exact arithmetic checks 7,088 minors
  across degrees zero through seven and orders one through four. Exact strict
  support and all-orders tail-composition signs remain open. The root command
  passes with 88 tests and no known dependency vulnerabilities. The change
  passed CI and merged as PR #55.
- Phase BB isolates the exact algebra behind continuous Cauchy--Binet. Summing
  the row-wise composition integrand over all permutations of an intermediate
  sample equals the product of the sampled exponential and truncated-power
  determinants. On strictly ordered samples this product has the exponential
  signature's weak sign at every order, by exponential strict sign regularity
  and the new all-orders truncated-power theorem. An exact rational `3 x 3`
  regression obtains `-2205` both from six permutation terms and from the
  determinant product. Product-measure permutation invariance, integration,
  and a strict positive-measure witness remain open. The root command passes
  with 89 tests and no known dependency vulnerabilities. The change passed CI
  and merged as PR #56.
- Phase BC completes the finite-product analytic identity behind Andreief.
  Coordinate permutations preserve both integrability and integration for a
  product of one common sigma-finite measure. Finite-sum integration turns the
  Phase BB symmetrization into the exact factorial normalization, while a
  separate determinant expansion and finite-product Fubini theorem identify
  the unsymmetrized integral with the determinant of pair integrals. Together
  these give a division-free finite-order Andreief theorem. An independent
  exact three-monomial unit-cube expansion gives `1/360` on both sides.
  Extending the ordered exponential/truncated-power sign to arbitrary samples
  by a common sorting permutation, and then proving strict support, remain the
  next analytic steps. The root command passes with 90 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #57.
- Phase BD removes the ordered-intermediate-sample restriction. Every
  injective finite real sample has a monotone sorting permutation, and the
  permutation signs of the sampled exponential columns and truncated-power
  rows cancel in their determinant product. If the sample is not injective,
  the repeated row makes the truncated-power determinant zero. This pointwise
  result integrates through Andreief to give weak-sign closure for finite
  integral-kernel composition. The common-cutoff representation and positive
  `t^8/7!` row scaling therefore prove the prescribed weak sign for every
  ordered eighth-order Taylor-tail minor on positive times and nonnegative
  frequencies. An exact rational regression checks all six permutations and
  the repeated-sample boundary. Strict support and hence strict all-orders
  tail sign regularity remain open. The root command passes with 91 tests and
  no known dependency vulnerabilities. The change passed CI and merged as
  PR #58.
- Phase BE closes the all-orders strict-support gap for the eighth-order
  Taylor tail. The explicit intermediate tuple starts halfway between zero
  and the first positive frequency knot, then uses adjacent-knot midpoints.
  Its truncated-power matrix is upper triangular with positive diagonal, so
  its determinant is positive at every order; the exponential determinant has
  the matching strict signed orientation. A continuous bump supported inside
  the common-cutoff cube turns this point witness into positive product
  measure, and the exact Andreief identity makes the composed determinant
  strict. Positive row scaling transfers strictness to every ordered tail
  minor on positive times and frequencies. Exact arithmetic checks the
  triangular witness through order 12, including orders beyond the polynomial
  degree. The remaining problem is quantitative: exploit this compatible tail
  sign in a cutoff-uniform GECP contraction. The root command passes with 92
  tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #59.
- The verification refresh changes no proved, observed, conjectured, or
  not-claimed mathematical statement above.

- Phase BF packages the positive-quadrant eighth-order Taylor tail as a global
  strictly sign-regular kernel and preserves that structure after the
  positive fermionic denominator scaling. The positive-frequency fermionic
  kernel is exactly an eight-separated-term Taylor component plus this tail.
  On the positive unit band the tail has uniform magnitude at most `1 / 8!`,
  so every exact complete-pivot run on the tail has selected-core determinant
  at most `(1 / 8!)^rank`, with no generic determinant factorial. A 120-digit
  eight-step regression independently checks the decomposition, pivots, and
  cumulative determinant envelope. This does not yet control determinants of
  the low-rank-plus-tail sum or extend the result to the `8p` dyadic orders.
  The root command passes with 93 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #60.

- Phase BG identifies and repairs a mixed-sign obstruction in the local
  decomposition. For the direct negative-exponential Taylor split, the
  one-tail-column coefficient of a `2 x 2` determinant can have the opposite
  sign from the full exponential determinant. Reflecting time instead gives
  `K(1-u,omega)=exp(u*omega)/(1+exp(omega))`: all eight low features are
  positive monomials, and the high tail is proved strictly totally positive.
  The denominator-scaled tail retains the unit-band `1 / 8!` error and
  `(1 / 8!)^rank` complete-GECP determinant bounds. A 140-digit deterministic
  sweep finds all mixed coefficients positive through order six; that
  coefficientwise statement remains experimental pending a finite
  Cauchy--Binet theorem. The root command passes with 94 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #61.

- Phase BH proves the finite coefficientwise statement selected by Phase BG.
  Scaling an arbitrary marked block of rectangular features by `lambda`
  weights each ordered Cauchy--Binet term by exactly the number of marked
  features it selects. Grouping terms by that count produces nonnegative
  coefficients whenever the base ordered terms are nonnegative. Generalized
  Vandermonde positivity discharges this hypothesis for finite positive
  exponential Taylor features at every order and cutoff. An independent exact
  rational `3 x 3` calculation agrees between determinant multilinearity and
  Cauchy--Binet grouping and finds every coefficient strictly positive. The
  infinite-tail limit and a quantitative mixed-determinant bound remain open.
  The root command passes with 95 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #62.

- Phase BI identifies each grouped feature-count coefficient with the sum of
  physical mixed-column determinants using exactly that many tail columns.
  Both expansions are packaged as determinant polynomials in a scalar tail
  weight; equality at every real weight and polynomial uniqueness give the
  coefficient identity. Consequently every fixed-cardinality mixed-column sum
  for a finite positive exponential Taylor split is nonnegative. The exact
  rational regression independently compares both coefficient constructions.
  Passing these signs to the exact infinite tail is now a finite determinant-
  continuity problem. The root command passes with 95 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #63.

- Phase BJ passes the finite positive-Taylor mixed-column signs to the exact
  infinite eighth-order tail. A generic continuity theorem handles each
  fixed-cardinality sum of column-choice determinants. Exact low/high feature
  identities connect the finite truncations to Phase BI, and entrywise Taylor
  convergence identifies the limit with `positiveEighthOrderTail`. Therefore
  every exact mixed-column coefficient sum is nonnegative at every order and
  tail-column count on nonnegative strictly ordered nodes. A 120-digit `3 x 3`
  regression sees monotone convergence through cutoff 64. A quantitative bound
  on the sum is still needed for GECP contraction. The root command passes with
  96 tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #64.

- Phase BK turns the exact coefficient signs into a quantitative scaling law.
  A low block factoring through `k` coordinates forces every mixed coefficient
  below degree `n-k` to vanish. Hence, for every `lambda >= 1`, nonnegative
  coefficients give the division-free comparison
  `lambda^(n-k) det(E+B) <= det(lambda E+B)`. The exact positive exponential
  split has `k=8`, and positive fermionic denominator scaling preserves the
  inequality. A 350-digit sweep verifies strict slack at orders 9--16 and
  weights 2, 8, and `8!`. The remaining quantitative target is a fixed-base
  upper envelope for the amplified determinant. The root command passes with
  97 tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #65.

- Phase BL proves that the amplified target itself has the needed qualitative
  structure. Nonnegative mixed-column coefficients make its determinant
  monotone in the tail weight, so every weight at least one preserves strict
  total positivity. At the balancing weight `8!`, the existing tail remainder
  estimate and the bounded base kernel give a uniform entry bound of two; a
  complete run of this amplified kernel consequently has selected-core
  determinant at most `2^n`. A 250-digit regression checks the entry envelope
  and positive ordered determinants through order eight. The outstanding gap
  is not positivity or entry control, but transfer of the `2^n` envelope to
  arbitrary ordered sample minors selected by the original GECP run. The root
  command passes with 98 tests and no known dependency vulnerabilities. The
  change passed CI and merged as PR #66.

- Phase BM removes the amplified-kernel pivot-path mismatch. A generic
  northwest-Schur induction proves that every strictly totally positive finite
  matrix has determinant at most the product of its diagonal entries. The
  kernel transfer bounds every ordered amplified positive-series sample by
  `2^n`, including samples selected by the original kernel, and Phase BK then
  yields `(8!)^(n-8) det K_n <= 2^n`, equivalently
  `det K_n <= 2^n / (8!)^(n-8)`, on the positive unit band for `n >= 8`.
  Exact arithmetic checks all 12,869 nonempty ordered subminors at order eight
  and every Schur stage. GECP composition, band/reflection assembly,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 99
  tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #67.

- Phase BN composes the new arbitrary-sample determinant decay with exact
  complete pivoting. Independent sorting of the selected row and column tuples
  preserves determinant magnitude, so the algorithmic pivot order needs no
  monotonicity assumption. On the reflected positive-series kernel over the
  positive unit band, every rank-`n >= 8` residual satisfies the Phase BM
  determinant-power envelope, and at `n = 10` its magnitude is uniformly at
  most `1/2`. A 180-digit three-grid regression independently checks the
  selected-core determinant and residual bounds. Endpoint closure, mixed-sign
  frequency assembly, dyadic rescaling, Conjecture G1, and Problem 4.2 remain
  open. The root command passes with 100 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #68.

- Phase BO eliminates the mixed-sign Laplace-expansion obligation. The affine
  change `z = 1 - omega` rewrites the original fermionic kernel on
  `[0,1] x [-1,1]` as positive row and column weights times `exp(t*z)`, with
  `z` in the closed interval `[0,2]`. Weighting the eighth-order positive tail
  by `128` preserves strict total positivity and stays within the pointwise
  envelope `2`, yielding the arbitrary-sample bound
  `|det K_n| <= 2^n / 128^(n-8)`. Exact complete-pivot composition gives a
  uniform one-half residual bound after twelve pivots on the full closed
  signed unit band. A 180-digit three-grid regression checks the determinant
  and residual bounds. Dyadic rescaling/continuation, Conjecture G1, and
  Problem 4.2 remain open. The root command passes with 101 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #69.

- Phase BP removes absolute frequency location from the local contraction
  theorem. For every real `center`, the affine coordinate
  `z = center + 1 - omega` maps `[center-1,center+1]` to `[0,2]`. The weighted
  positive-exponential term is exactly the original fermionic kernel and is
  therefore at most one on the full time strip, uniformly in `center`.
  Weight-`128` tail amplification consequently gives the same
  `|det K_n| <= 2^n / 128^(n-8)` envelope and the same twelve-pivot one-half
  residual bound on every translated width-two band. A 180-digit regression
  checks centers `-8`, `-2`, `3`, and `9`. This is an original-kernel local
  block, not yet a restart theorem for a residual conditioned by earlier
  pivots. The root command passes with 102 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #70.

- Phase BQ proves the exact algebraic restart bridge and rules out the
  misleading one-sided substitute. A rational finite positive-Taylor matrix
  with a rank-two frozen background has positive original determinant
  `91063777 / 85899345920`, while weight-`128` amplification of only the three
  fresh-column tails makes the determinant `-32176499 / 983040`. In contrast,
  if the actual final Schur residual sample admits `E + L*R` with nonnegative
  mixed coefficients and residual low rank `k`, then its determinant magnitude
  amplifies by `lambda^(n-k)`, independently of the earlier prefix length.
  Phase AJ transports this to the augmented original-kernel determinant with
  the selected core only as a multiplicative factor. The remaining analytic
  target is now exact and unavoidable: construct such a compatible split for
  arbitrary fermionic final residuals. The root command passes with 103 tests
  and no known dependency vulnerabilities. The change passed CI and merged as
  PR #71.

- Phase BR constructs the compatible split left open by Phase BQ. A fixed
  GECP run acts linearly on companion matrices, commutes with multiplication
  by a fixed right feature factor, and therefore transports
  `K = E + L*R` exactly to
  `finalResidual = conditionRows(E) + conditionRows(L)*R`. Applied to the
  positive-series fermionic decomposition, this preserves the rank-eight
  feature block through an arbitrarily long selected prefix and feeds
  directly into the Phase BQ amplification theorem. Exact three-pivot
  arithmetic independently verifies the transported split and rank. It also
  reveals the next obstruction: conditioned-tail mixed coefficients need not
  remain nonnegative—the first example has signs `0, 0, +, -, +`. The next
  useful invariant should therefore exploit alternating signs, real-rootedness,
  or a direct two-weight determinant comparison rather than coefficientwise
  positivity. The root command passes with 104 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #72.

- Phase BS removes the false analytic premise exposed by Phase BR. After the
  rank-eight forced zero is factored out, every conditioned determinant pencil
  has degree at most eight, regardless of its coefficient signs or roots.
  Exact Lagrange interpolation at the nine weights `56, 65, ..., 128` gives a
  norming constant below `518173`. Lean combines this with rank vanishing to
  prove a sign-free determinant amplification theorem and the sharp usable
  corollary: at fresh order eighteen, a `2^18` determinant envelope at each
  norming weight forces the physical determinant below `2^-18`. Exact
  arithmetic simultaneously exhibits failures of coefficient positivity,
  real-rootedness, and single-weight `+128` domination, confirming that the
  finite norming set is not cosmetic. The remaining target is the nine
  constant-base companion envelopes. The root command passes with 105 tests
  and no known dependency vulnerabilities. The change passed CI and merged as
  PR #73.

- Phase BT replaces asymmetric companion row conditioning by the natural
  two-sided fixed-pivot oblique projection. The map is linear, annihilates the
  selected row and column, recovers the ordinary Schur update on the pivoting
  matrix, and transports `L*R` as separately conditioned left and right
  factors through an arbitrary prefix. The positive-series fermionic residual
  consequently retains an exact two-sided rank-eight split whose low-block
  determinant is an ordered Cauchy--Binet sum of matched conditioned feature
  minors. In exact arithmetic, the earlier one-sided sign witness
  `0, 0, +, -, +` becomes `0, 0, +, +, +`; all 2,760 exhaustively enumerated
  samples from deterministic randomized prefixes have a common nonzero
  coefficient orientation. This evidence identifies matched-minor sign
  coherence as the next theorem but does not prove it. The root command passes
  a 3,568-job Lean build with 106 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #74.

- Phase BU identifies the exact Plücker coordinates behind Phase BT's sign
  evidence. A left conditioned feature minor equals an original selected-core
  determinant bordered by left feature columns, divided algebraically by the
  selected core; the right feature minor has the transposed identity.
  Multiplying them makes the selected-core determinant a square, so every
  conditioned matched-minor sign is equivalent to the sign of a product of
  two original bordered determinants. The fermionic specialization exposes
  this equivalence term by term in the ordered Cauchy--Binet expansion. Exact
  arithmetic verifies 2,382 left, right, and product identities. The remaining
  factors are confluent exponential determinants: `t^d/d!` and the
  denominator-scaled `omega^d` features arise from derivatives at exponential
  parameter zero. The root command passes a 3,570-job Lean build with 106
  tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #75.

- Phase BV shows that the termwise route from Phase BU is false but that its
  invariant sum survives. In the exact positive-unit finite-Taylor witness,
  three complete pivots select rows `[6,0,3]` and columns `[6,0,2]`; among the
  56 order-three feature subsets, exactly two products of left and right
  bordered determinants have the opposite sign. Adjoining the fixed pivot
  blocks `U,V` to arbitrary feature factors `A,B` instead packages all of
  those products into one graded Cauchy--Binet coefficient. Lean proves that
  weighting the fixed coordinates gives `mu • (U*V) + A*B` and that the
  full-fixed-degree Plücker pairing is exactly the corresponding mixed-column
  coefficient. The witness verifies the signed sum identity and a positive
  total without discarding the two negative coordinates. The root command
  passes a 3,571-job Lean build with 107 tests and no known dependency
  vulnerabilities. Positivity of the full coefficient for arbitrary
  fermionic prefixes remains the analytic target. The change passed CI and
  merged as PR #76.

- Phase BW keeps the Phase BV cancellations inside one matrix. The bordered
  pencil with blocks `1`, `V`, `-mu • U`, and `B` has determinant exactly
  `det(mu • (U*V) + B)` for arbitrary finite real matrices, with no
  invertibility premise. It therefore linearizes both the mixed-column
  polynomial and the combined fixed-plus-feature block polynomial at every
  weight. Exact arithmetic verifies the identity on the Phase BV witness at
  weight one and all nine Phase BS weights. The same witness exposes a
  negative signed four-cycle crossing the four blocks. Since that cycle sign
  is invariant under independent row and column sign changes, no sign gauge
  can make the bordered matrix entrywise nonnegative, ruling out a direct
  totally-nonnegative-matrix proof. The root command passes a 3,572-job Lean
  build with 107 tests and no known dependency vulnerabilities. The change
  passed CI and merged as PR #77.

- Phase BX aligns the bordered approach with the actual two-sided companion
  envelope. For an invertible selected core `C`, the closed companion block is
  `H - K_L C⁻¹ H_R - H_L C⁻¹ K_R + K_L C⁻¹ H_C C⁻¹ K_R`.
  A single double-bordered matrix has upper-left block
  `[[C,H_C],[0,C]]` and this expression as its Schur complement. Lean proves
  its determinant and absolute determinant are `det(C)^2` and
  `|det(C)|^2`, respectively, times the companion determinant. Exact
  solve-based rational checks show that this closed expression agrees
  entrywise with three iterative two-sided companion pivots at weight one and
  all nine Phase BS weights. Negative order-two minors persist in the large
  matrix, so the factorization exposes rather than removes the cancellation.
  The root command passes a 3,573-job Lean build with 107 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #78.

- Phase BY removes the last algebraic ambiguity in that representation. For
  every successful run, Lean proves that fixed row and column projectors
  commute, that each iterated one-sided operator is its selected-core
  correction, and that `conditionBoth` is exactly the global four-term Schur
  expression. Arbitrary finite samples therefore satisfy Phase BX's
  double-border identity with determinant
  `det(selectedCore)^2 * det(sample(conditionBoth H))`. The inverse instance
  comes only from the already proved nonzero selected-core determinant. An
  independent exact regression checks 144 deterministic randomized cases
  across prefix lengths one through three, four companion weights, and
  varying fresh blocks using rational solves rather than explicit inverses.
  This identifies the nine Phase BS companion premises exactly with
  normalized double-border bounds, but does not prove those bounds or the
  resulting half contraction. The root command passes a 3,574-job Lean build,
  the public axiom audit, and 107 tests with no known dependency
  vulnerabilities. The change passed CI and merged as PR #79.

- Phase BZ recognizes the double border as a tangent/dual-number lift. The
  matrix `[[K,H],[0,K]]` contains every selected-core double border as one
  distinguished minor, and when `K=A*B`, `H=A*W*B` it factors through
  `diag(A,A)`, the single upper-unitriangular block `[[1,W],[0,1]]`, and
  `diag(B,B)`. Lean specializes this to the positive-series amplified kernel
  and proves that a normalized `2^18` envelope for the nine distinguished
  minors implies the Phase BS order-eighteen half-power bound. Exact rational
  tests verify 24 arbitrary feature factorizations and 360 complete-prefix
  determinant polynomials. Although the full lift still has negative
  order-two minors, every tested polynomial coefficient has the explicit
  cross-inversion orientation predicted by sorting the selected and fresh
  row/column tuples. This is strong evidence for a flag-minor or two-layer
  planar-network theorem, not a proof of the missing envelope. The root
  command passes a 3,577-job Lean build, the public axiom audit, and 108 tests
  with no known dependency vulnerabilities. The change passed CI and merged
  as PR #80.

- Phase CA makes the tangent sign problem factorwise. An ordered
  three-factor Cauchy--Binet theorem expands each sampled tangent determinant
  over a pair of feature flags, and a modular Lean theorem composes separate
  left, middle, and right orientations. Twenty-four exact shuffled
  rectangular checks reproduce the expansion. The middle matrix is not
  totally nonnegative in layer-major order, but interleaving each feature's
  two layers turns it into contiguous `2 x 2` unitriangular cells; the
  induced shuffle parity is therefore the most promising next exact lemma.
  The root command passes a 3,579-job Lean build, the public axiom audit, and
  108 tests. The change passed CI and merged as PR #81.

- Phase CB makes the middle-factor sign concrete. Feature interleaving turns
  diagonal `[[1,W],[0,1]]` into independent `[[1,w_k],[0,1]]` cells. Lean
  proves every cell totally nonnegative, the full determinant equal to one,
  and the exact tuple-orientation transfer from interleaved TN back to the
  layer-major minors used in Phase CA. Exhaustive rational arithmetic checks
  198,626 feature-major minors and 198,626 oriented layer-major minors through
  five cells; negative unoriented layer-major minors still occur. The single
  remaining middle lemma is all-orders TN of the sparse upper-bidiagonal
  direct sum. The root command passes a 3,580-job Lean build, the public axiom
  audit, and 109 tests. The change passed CI and merged as PR #82.

- Phase CC closes the middle-factor theorem at every order. The reusable
  `totallyNonnegative_of_nonnegative_of_noncrossing` criterion shows that an
  entrywise nonnegative matrix with planar, noncrossing support is totally
  nonnegative because every nonidentity determinant matching vanishes.
  Applied to the interleaved tangent cells, it removes Phase CA's conditional
  middle-factor premise for arbitrary nonnegative diagonal weights. The root
  command passes a 3,581-job Lean build and 110 tests; the change passed CI
  and merged as PR #83.

- Phase CD identifies a decisive limitation of the factorwise route. In the
  first positive-monomial witness with mixed signs—two physical nodes and
  three features—the ordered triple Cauchy--Binet expansion has 12 nonzero
  terms, including the formally certified term `-1/192`, even though the
  complete determinant is the positive value `24929/27648`. Lean therefore
  proves that no separate left/middle/right orientation functions can satisfy
  Phase CA's termwise interface on this witness. The two-feature case has only
  the positive terms `1/2` and `1/36`, so three features are minimal in this
  family. The surviving route must preserve cancellation globally, for
  example through a confluent planar network, compound/exterior algebra, or a
  sign-reversing involution; proving the final envelope, Conjecture G1, and
  Problem 4.2 remains open. The root command passes a 3,582-job Lean build and
  111 tests; the change passed CI and merged as PR #84.

- Phase CE performs the cancellation before estimating signs. For every
  shared factorization `K=A*B`, `H=A*W*B`, the complete distinguished tangent
  determinant is `det(selectedCore)^2` times the determinant of
  `conditionRows(A) * W * conditionColumns(B)`, and Lean identifies the full
  ordered double-flag sum with exactly this quantity. The Phase CD negative
  flags therefore live inside a globally positive conditioned determinant.
  In 5,000 deterministic exact generalized-monomial cases, arbitrary
  nonnegative diagonal weights retain the selected/fresh cross-inversion
  orientation (3,964 strict and 1,036 zero), while the signed weight pair
  `(-2,0)` already gives `-4/169`. The dense census initially pointed toward a
  nonnegative diagonal-action theorem, but Phase CF's sparse adversarial test
  shows that inference is false. The root command passes a 3,584-job Lean
  build and 112 tests.
  The change passed CI and merged as PR #85.

- Phase CF isolates individual matched coordinates with a one-hot inclusion
  matrix. The induced projection is entrywise nonnegative and diagonal, yet
  its inserted determinant is exactly one left/right minor product. Thus the
  two negative Phase BV coordinates, supported on `(0,1,5)` and `(0,2,5)`,
  formally obstruct positivity on the whole nonnegative diagonal cone. The
  earlier dense random search missed this sparse boundary. The physical
  family survives the correction: all 90 combinations of contiguous cutoffs
  `0,...,8` and weights `1,56,65,...,128` are strictly positive in the same
  exact witness. The next target is therefore positivity of grouped tail-count
  coefficients for contiguous step weights. The root command passes a
  3,585-job Lean build and 112 tests.
  The change passed CI and merged as PR #86.

- Phase CG identifies the cancellation unit that survives Phase CF. After
  multiplying by the selected-core square, each conditioned coefficient of a
  contiguous step-weight determinant is exactly the sum of original bordered
  Plucker products having one fixed number of tail indices. The full step
  determinant weights these groups by the corresponding powers of `lambda`.
  This statement keeps the two sparse negative coordinates but pairs them
  inside their tail-count groups. In the exact Phase BV witness every group is
  nonnegative at every cutoff, and all 90 cutoff/node determinants equal the
  grouped polynomial evaluation and are positive. The remaining theorem is
  now the analytic sign of these grouped bordered pairings for arbitrary
  positive-Taylor prefixes. The root command passes a 3,586-job Lean build and
  112 tests.
  The change passed CI and merged as PR #87.

- Phase CH isolates the smallest observed cancellation inside one grouped
  coefficient. At one pivot with three consecutive Taylor features and exact
  nodes `x=(1/8,1/4)`, `y=(6,22)`, the bordered products are `561/64`,
  `-7/256`, and `801/32`: coordinatewise positivity fails, but every proper
  nonempty prefix and suffix group is positive. Lean proves both nontrivial
  contiguous-step determinants strictly positive for every nonnegative
  weight. Exact arithmetic exhausts 6,468 ordered grid cases without any
  two-feature cancellation and finds 800 three-feature cancellations, with
  this witness first lexicographically. This points to an adjacent-degree
  variation-diminishing inequality as the next analytic target; it does not
  establish the arbitrary-prefix grouped theorem. The root command passes a
  3,587-job Lean build and 113 tests.
  The change passed CI and merged as PR #88.

- Phase CI turns the Phase CH witness into a full-domain theorem. Each of the
  three generic one-pivot Taylor coordinates has an exact closed factorization,
  and both adjacent two-coordinate sums factor as the product of the ordered
  row gap, ordered column gap, and a polynomial with positive coefficients.
  Hence both groups are strictly positive for arbitrary ordered nonnegative
  nodes—even beyond the physical row interval—and both nontrivial step
  determinants remain positive for every nonnegative scalar weight under the
  sharp boundary hypotheses. Exact arithmetic verifies 2,198 boundary,
  interior, super-physical, and randomized rational cases. The structure now
  suggests a general conditional-covariance proof: normalize the pivot
  feature weights, center two increasing exponential sequences, and decompose
  each prefix or suffix into conditional covariance plus same-sign mean
  shifts. The root command passes a 3,588-job Lean build and 114 tests.
  The change passed CI and merged as PR #89.

- Phase CJ removes the three-feature restriction at one pivot. A weighted
  covariance numerator is half a sum of pairwise same-direction increments.
  A two-block identity then decomposes every globally centered prefix or
  suffix into a local covariance term and the product of two same-direction
  conditional-mean shifts. Lean proves both terms nonnegative for arbitrary
  positive weights and increasing geometric sequences. The finite-Taylor
  bordered products match this model exactly, so every nonempty prefix,
  suffix, and nonnegative two-level step is nonnegative for arbitrary feature
  count and ordered positive nodes. Exact checks cover 4,644 weighted groups
  and 6,782 Taylor groups through 12 features, while retaining 218 negative
  individual coordinates. The next obstacle is genuinely compound: multiple
  fresh rows/columns require association on a feature-subset lattice, likely
  via an `MTP2`/FKG or biorthogonal-ensemble argument. The root command passes
  a 3,592-job Lean build and 115 tests.
  The change passed CI and merged as PR #90.

- Phase CK identifies the surviving higher-order invariant more sharply.
  Absolute-value `MTP2` fails on the Phase BV ordered-subset weights in 372 of
  3,136 comparisons, so a direct FKG proof is unavailable. Nevertheless all
  nine exact threshold polynomials are coefficientwise nonnegative,
  log-concave, and real-rooted by exact quadratic/cubic discriminants. The
  same signatures hold for 1,481 oriented threshold polynomials from 200
  deterministic exact compound Taylor cases through fresh order three, even
  though 147 cases violate `MTP2`. Lean proves a division-free common-factor
  determinant-pencil formula and nonnegative/strict-positive consequences
  from nonnegative diagonal generalized-spectrum coordinates. This is a
  reduction, not the missing compound theorem: simultaneous diagonalization
  of arbitrary conditioned fermionic low/high pairs, the infinite-tail
  passage, normalized tangent envelope, order-eighteen contraction, dyadic
  recurrence, Conjecture G1, and Problem 4.2 remain open. The root command
  passes a 3,593-job Lean build and 116 tests. The change passed CI and merged
  as PR #91.

- Phase CL makes the fresh-order-two case exact. Lean expands every `2 x 2`
  determinant pencil into its endpoint determinants and mixed coefficient,
  identifies the discriminant inequality `4*c0*c2 ≤ c1^2`, and constructively
  factors every nonnegative Lorentzian quadratic into two affine factors with
  nonnegative offsets. This yields an all-nonnegative-weight determinant
  certificate without an inverse. In 500 deterministic exact finite-Taylor
  runs, all 3,868 oriented threshold pencils have nonnegative coefficients and
  satisfy the Lorentzian inequality. Only 1,030 satisfy either tempting
  single row- or column-cross-minor sign; 2,838 failures show that the next
  proof must retain the full mixed compound sum, as in an
  Alexandrov--Fenchel or Hodge--Riemann argument. The fermionic Lorentzian
  inequality itself, higher fresh order, the infinite-tail passage,
  normalized tangent envelope, order-eighteen contraction, dyadic recurrence,
  Conjecture G1, and Problem 4.2 remain open. The root command passes a
  3,594-job Lean build and 117 tests. The change passed CI and merged as PR
  #92.

- Phase CM supplies the missing abstract Lorentzian inequality behind the
  Phase CL audit. For two ordered blocks, the cross-increment coefficient
  obeys a division-free law of total covariance: after multiplying by the
  two block masses it is the sum of two rescaled within-block covariance
  numerators and a product of same-direction mean shifts. Every term is
  nonnegative, and the first two terms dominate the geometric mean by one
  explicit square, proving `4*c0*c2 ≤ c1^2`. Lean proves the generic finite
  weighted result and a geometric prefix/suffix specialization. Exact
  arithmetic checks 1,764 weighted splits and all 11,024 oriented
  two-feature products from the Phase CL conditioned-Taylor runs, with no
  failure. The remaining bridge is representational: arbitrary conditioned
  Taylor rank-two factors must be rewritten as positive pair weights times
  two comonotone ratio differences. Higher fresh order and the infinite-tail
  passage remain open. The root command passes a 3,594-job Lean build and 117
  tests. The change passed CI and merged as PR #93.

- Phase CN turns the Phase CM inequality into a direct rank-two matrix-pencil
  interface. Affine chart factorizations rewrite each matched Plucker product
  as a positive pair weight times two ratio increments; the within-block and
  cross-block sums are exactly the weighted covariance and cross increment,
  while scalar-row Cauchy--Binet identifies them with the endpoint
  determinants and mixed coefficient of the original `2 x 2` pencil. Lean
  therefore proves the Lorentzian discriminant from a division-free chart
  certificate. The obvious global chart is genuinely insufficient: among
  500 exact conditioned-Taylor cases, only 42 are degree-monotone and 276 of
  2,868 nontrivial splits have opposite aggregate mean shifts. Projectively,
  however, every left/right pair has the same cyclic order and every degree
  prefix is a circular interval in all 2,868 splits. The next target is now
  exact: prove this pivot-depth cyclic order, then choose a cutoff-specific
  point at infinity. Higher fresh order and the infinite-tail passage remain
  open. The root command passes a 3,595-job Lean build and 117 tests.
  The change passed CI and merged as PR #94.

- Phase CO makes the Phase CN certificate genuinely projective. The mixed
  coefficient under `M ↦ P*M*Q` scales by `det(P)det(Q)`, exactly like both
  endpoint determinants, so the discriminant scales by
  `(det(P)det(Q))^2`. Lean proves Lorentzianity is preserved even through a
  singular degeneration and is reflected by every invertible basis change.
  The projective-affine composition theorem now allows each degree cutoff to
  choose its own point at infinity and transfer the resulting ordered-chart
  proof back to the original pencil. Exact rational checks cover 512
  invertible changes, 271 with negative determinant product, and a singular
  boundary case. The only remaining finite fresh-order-two analytic input is
  the conditioned Taylor cyclic-prefix theorem and its separating-chart
  construction. The root command passes a 3,595-job Lean build and 118 tests.
  The change passed CI and merged as PR #95.

- Phase CP constructs the cutoff-specific chart explicitly. The pole basis
  `[[-c,1],[1,0]]` always has determinant `-1` and sends each homogeneous
  feature `a_k(1,f_k)` to coordinates
  `(a_k(f_k-c),a_k)`. Away from a collision, the new affine ratio is exactly
  `1/(f_k-c)`. Lean proves the two-sided block transformation and composes it
  with the Phase CO invariance theorem, reducing every original rank-two
  Lorentzian pencil to positive transformed weights and reciprocal-ratio
  block order. Exact rational checks exercise three exterior/interior pole
  pairs, all five nontrivial cutoffs of a six-feature example, and one rejected
  collision. The remaining finite fresh-order-two input is now only the
  conditioned-Taylor cyclic-prefix/separating-pole theorem; CP does not prove
  that analytic invariant. The root command passes a 3,595-job Lean build and
  119 tests.
  The change passed CI and merged as PR #96.

- Phase CQ separates the genuinely projective data from chart bookkeeping.
  An oriented pole basis uses the ratio `σ/(f-c)`; an independent common scale
  flips all transformed singleton weights without changing that ratio. Lean
  proves the exact determinant and block actions, shows that the matched
  coordinate scales by `scale²*σ`, and derives reciprocal `MonovaryOn`
  automatically from positive transformed weights and oriented nonnegative
  matched minors. The final theorem leaves only cross-ordering separators and
  projective minor signs to the conditioned-Taylor argument. Exact rational
  gap search constructs the full certificate on the first candidate pair for
  all 2,868 nontrivial Phase CN cutoffs; both negative orientations and
  negative raw chart scales occur, so neither parameter is redundant. This
  does not prove the arbitrary-pivot Taylor invariant. The root command passes
  a 3,595-job Lean build and 119 tests.
  The change passed CI and merged as PR #97.

- Phase CR clears the conditioning denominators from the remaining rank-two
  sign problem. For an arbitrary successful run, the selected-core determinant
  squared times a conditioned matched coordinate is exactly the product of
  the corresponding original left and right feature-augmented core
  determinants. Since the core square is strictly positive, multiplying by
  any orientation scalar preserves a two-way nonnegativity equivalence. Exact
  rational reconstruction verifies the identity for all 11,024 feature pairs
  in the 500 arbitrary-pivot conditioned-Taylor cases. The next analytic target
  is therefore a pure sign/cyclic-interval theorem for augmented generalized-
  Vandermonde minors, with no Schur quotients left in the statement. The root
  command passes a 3,596-job Lean build and 119 tests.
  The change passed CI and merged as PR #98.

- Phase CS identifies what does and does not survive before the augmented-
  determinant product is formed. A generic Lean lemma proves that two scalar
  Plucker coordinates with independent reorientations and a shared nonzero
  sign schedule yield the desired oriented product sign; composing it with
  Phase CR gives the conditioned matched-coordinate theorem. The exact audit
  confirms shared nonzero schedules for all 11,024 augmented pairs and sees
  both schedule signs. It also refutes the tempting stronger rule that the raw
  schedule is determined only by pivot depth and the degree pair. Thus the next
  theorem must establish left/right oriented-matroid synchronization, not a
  false universal confluent sign formula. The root command passes a 3,596-job
  Lean build and 119 tests.
  The change passed CI and merged as PR #99.

- Phase CT converts the remaining left/right sign comparison into a single
  topological obstruction. A continuous scalar path that never vanishes has
  same-sign endpoints, and the determinant specialization derives continuity
  from entrywise continuity. After positive denominator and factorial
  normalization, the two conditioned Taylor determinants become endpoints of
  one canonical confluent exponential determinant path. Exact arithmetic
  verifies both endpoint identities, midpoint nonvanishing, and common sign
  for all 11,024 augmented pairs. This finite audit is evidence for, not a
  proof of, nonvanishing at every path parameter. The root command passes a
  3,597-job Lean build and 119 tests.
  The change passed CI and merged as PR #100.

- Phase CU rules out the overbroad ordered-chamber completion of the Phase CT
  homotopy. With row tuples `(1,2,3,4)/9` and `(5,6,7,8)/9`, logarithmic
  column nodes `9 log(12/11) < 9 log(11/10)` in `(0,1)`, and appended columns
  `1,x^3`, positive exponential scaling reduces the endpoint determinants to
  the exact rationals `83/194059800` and `-1/194059800`. Lean proves the
  straight interpolation remains strictly ordered in the open unit chamber,
  so continuity forces a singular member. This does not refute the physical
  paired-subtuple path: actual endpoints additionally satisfy
  `X_sel ⊂ X_aug` and `Y_sel ⊂ Y_aug`. The root command passes a 3,598-job
  Lean build and 120 tests.
  The change passed CI and merged as PR #101.

- Phase CV shows that complete-pivot coupling does not repair the chirotope
  claim. The exact kernel `K(i,j)=(101/100)^((i+1)(j+21))` is a genuine
  exponential kernel on increasing nodes in `(0,1)`. Lean certifies `(3,3)`
  as its first complete pivot and `(0,0)` as the second complete pivot of the
  Schur residual. Those selections are the two endpoints on both axes, yet
  the crossed degree-`(0,3)` augmented determinants have opposite signs.
  Thus neither ordered-chamber homotopy, paired-subtuple containment, endpoint
  pattern, nor complete-pivot coupling yields universal common chirotopes.
  Aggregate cancellation or an absolute determinant bound remains possible.
  The root command passes a 3,600-job Lean build and 121 tests.
  The change passed CI and merged as PR #102.

- Phase CW removes coordinatewise chirotope signs from the
  fresh-order-two critical path. The exact law of total covariance shows that
  positive block masses, nonnegative within-block covariance numerators, and
  one nonnegative product of aggregate block-mean shifts imply the sharp
  Lorentzian discriminant inequality. The affine rank-two and matrix-pencil
  lifts require neither division nor an inverse. The Phase CV witness retains
  its negative degree-`(0,3)` coordinate, yet all 42 of its tested aggregate
  pencils remain Lorentzian; 40 use the new certificate and two have a
  rank-one endpoint. A broader deterministic exact corpus verifies 16,424
  conditioned pencils with minimum normalized discriminant margin above
  `0.37`. Proving the aggregate mean-shift sign for arbitrary infinite-tail
  fermionic prefixes remains open.
  The change passed CI and merged as PR #103.

- Phase CX rules out a tempting global Hodge--Riemann completion while
  preserving exactly the cone used by degree cutoffs. The rational kernel
  `K(i,j)=(21/20)^((i+1)(j+1))` is an exponential kernel on nodes strictly
  inside the physical open unit square, and `(3,3)` is its first complete
  pivot. After conditioning the first four normalized Taylor features, the
  matched-coordinate Hessian has an explicit two-dimensional positive-
  definite restriction: the vectors `(10,177,686,706)` and
  `(999,-38,-6,0)` violate reverse Cauchy--Schwarz. Nevertheless, each of the
  three nontrivial prefix/suffix indicator planes satisfies the sharp reverse
  Cauchy--Schwarz inequality. Thus global stability or a one-positive-
  eigenvalue theorem cannot be the missing invariant; the next proof must use
  the one-sign-change/prefix cone, plausibly through variation diminishing or
  interlacing of cumulative mixed determinants.
  The change passed CI and merged as PR #104.

- Phase CY proves a local mechanism for neighboring degree cutoffs. Moving a
  rank-one feature changes `A+B*t+C*t^2` by `(1-t)*(l+r*t)`, so the two
  pencils share the value at `t=-l/r`. Clearing the denominator gives the
  invariant `T=A*r^2-B*l*r+C*l^2`. With positive transfer sums and nonnegative
  endpoint determinants, `T≤0` implies both adjacent discriminants are
  nonnegative. More specifically, Lean proves
  `T=-projected_column_cross*projected_row_cross` for an outer-product moved
  feature, and derives the two matrix inequalities from that aggregate sign.
  The exact corpus checks the identity, including the determinant orientation,
  on 13,924 transfers. All 9,170 strict interior transitions have `T<0`;
  positive numerators occur next to rank-one endpoints. The Phase CV witness
  has one such exception among 35 transitions and passes all 21 strict
  interior transitions. The next missing theorem is sign agreement of the
  cumulative projected crosses for arbitrary interior conditioned Taylor
  blocks. Adjacent-root ordering, the infinite-tail sign, and the final
  determinant envelope remain open.
  Local root verification passes a 3,602-job Lean build and 123 tests.
  The change passed CI and merged as PR #105.

- Phase CZ refutes the universal interior sign proposed after Phase CY. The
  exact nine-term Taylor kernel on rows `(17,18,19,20)/21` and columns
  `(15,16,17,19)/20` has complete pivots `(3,3)` followed by `(0,0)`.
  All nine features reconstruct the pivot kernel. After those two pivots,
  moving degree three gives a strictly positive canonical-node numerator even
  though all four neighboring endpoint determinants and both transfer sums
  are positive. Both cutoff discriminants remain positive. Moving degrees
  two and four instead supplies negative-node certificates for the affected
  cutoffs. This is a strict interior failure with full-rank blocks, rather
  than the boundary exception seen in Phase CY.
  The bounded exact audit checks every remaining fresh row/column pair in 100
  deterministic finite-Taylor cases: 54,845 transfers, four failed fixed
  nodes, and coverage of all 47,496 interior cutoffs by one of their two
  neighboring nodes. A 500-case run observes 339,950 transfers, six failed
  nodes, and coverage of all 292,076 cutoffs. These finite observations leave
  adaptive coverage open for arbitrary runs and the infinite tail.
  An independent actual-exponential regression at 80 and 160 decimal digits
  agrees on the positive normalized degree-three numerator
  `0.00010354916980416281785667342918876537`, with negative neighboring nodes.
  This is numerical evidence only; all complete-pivot claims in this phase
  concern the specified finite node sets.
  Root verification passes a 3,603-job Lean build and 126 tests.
  The change passed CI and merged as PR #106.

- Phase DA proves an analytic one-pivot augmented-minor sign theorem for the
  actual exponential. For `0<a≤1`, `0≤x<y<z≤1`, and natural `i<j`,
  `det[exp(a*t),t^i,t^j]` is negative precisely when `i=0,j≥2` and positive
  otherwise. `expMonomialMinor_signed_pos` makes that schedule explicit;
  `expMonomialMinor_matched_mul_pos` proves that independent ordered triples,
  even at independent rates in `(0,1]`, have a strictly positive product for
  the same degree pair. The proof normalizes rows and uses a monotone ratio
  of derivatives through Cauchy's mean value theorem. Boundary nodes zero
  and one are included. This is a proved infinite-family base case, not a
  claim of sign synchronization for arbitrary prefixes; Phase CV still
  rules out that latter claim.
  Independent 80/160-decimal regressions check 324 determinants; a separate
  100-decimal check covers 66 paired products. An exact depth-ten
  coherent conditioned-Taylor fixture has a bad
  degree-ten node, refuting the proposed `depth+1` localization shortcut in
  that setting. The deeper 400-case exploratory audit covers 5,964 cutoffs
  among 6,364 transfers, with sixteen positive nodes and ten outside the
  proposed exceptional degree. Adaptive coverage remains a conjecture.
  Local root verification passes a 3,605-job Lean build and 129 tests; the
  new public results use only the accepted Lean foundations.
  The change passed push and PR CI and merged as PR #107.

Phase DB refutes the neighboring canonical-node coverage conjecture suggested
after Phase CZ. A coherent 13-term Taylor kernel on nine open-unit rows and
nine open-unit columns has seven certified complete pivots, after which
degrees seven and eight both have positive canonical numerators at the same
interior cutoff. The left and right transfer sums are all positive, the three
neighboring pencils have positive coefficients, and the middle discriminant
is positive. `phaseDB_not_adaptive_coverage` therefore rules out that proof
mechanism without refuting the surviving Lorentzian inequality.
The generated integer certificate checks every pivot update and maximum,
Taylor reconstruction, both feature projectors, and transfer of the signs
back from normalized integers to the actual conditioned factors. Independent
160/320-decimal actual-exponential checks agree on both failed signs to 45
significant digits; they are numerical evidence, not a formal infinite-tail
obstruction. A separate exact weighted outer-prefix regression also loses
positive transfer sums while preserving its local positive discriminant.
All complete-pivot claims here concern the specified finite node sets, not
maximization over the full continuous domain. At DB delivery the proposed
next routes were a direct cumulative-pencil bound or a sign-free companion
envelope; DD later establishes a different local-minor route.
Local root verification passes a 3,608-job Lean build and 133 tests, using
only the accepted public theorem axioms.
The change passed push and PR CI and merged as PR #108.

Phase DC proves a new all-orders one-step bridge, rather than DC's original
arbitrary-prefix companion envelope. For a strictly totally positive ordered
`(n+1)`-sample, a column-maximal pivot at ranks `(p,q)` satisfies
`abs(det(Schur[p,q])) ≤ det(A[delete q,delete q])`.
The old row set changes but the local columns do not. In particular the
replacement old row set may contain the pivot row; deleting `(p,q)` on both
sides would be a different claim. The proof establishes singleton
Hadamard--Fischer without symmetry, factors the arbitrary-pivot determinant
using independent row/column permutations, and applies column maximality.
`positiveSeriesFermionicKernel_localMinor_schur_le` applies this result to the
actual exponential/fermionic kernel on positive arguments with no unit-band
restriction. Every minor order, including zero, is covered.

Independent exact rational checks cover 2,520 column-maximal sample steps
through fresh order four. A separate 1,376-check finite complete-prefix audit
finds contraction of the same-column maximum over old row tuples and positive
minors after a single cumulative row/column sign gauge. The latter is evidence
for the next theorem, not a Lean proof of arbitrary-prefix closure. Phase AL's
independent compound-order gauges are not silently promoted to that property.

The direct companion audit uses the actual exponential, both conditioned
feature factors, original-kernel pivots, eighteen fresh local columns, and all
nine weights `56,65,...,128`. Nine configurations through depth 64 and cutoff
256 agree at 320/640 decimal digits on pivots and 45-significant-digit
determinant ratios. No sampled envelope fails. At depth 64/cutoff 256, the
largest normalized determinant is approximately `1.0326550577e-1234`; this
quantity is stored as a decimal string, not an underflowed machine zero.
These are finite-grid, non-interval-certified numerical observations. At DC delivery the
continuous/arbitrary-prefix companion envelope, residual gauge closure,
initial dyadic local budget, G1, and Problem 4.2 remain open.
The root gate passes a 3,610-job Lean build and 140 tests with only the accepted
public axioms and no known dependency vulnerabilities. Push and PR CI pass;
the change merged as PR #109 on 2026-10-19.

Phase DD proves the common gauge that Phase DC was missing. For any real
kernel with positive ordered minors and any successful exact prefix, let
`P(x)=prod_i(x-selectedRow_i)` and `Q(y)=prod_i(y-selectedColumn_i)`. Then
`sign(P(x))*R(x,y)*sign(Q(y))` has strictly positive ordered minors at every
order on the fresh coordinate sets. The same two weights are used at every
order. `Run.gaugedFinalResidual_strictSignRegular` proves this directly from
the multi-bordered determinant identity and the selected/fresh factorization
of Vandermonde determinants; it does not infer it from Phase AL's separate
compound gauges. The entry and absolute-minor preservation theorems include
unordered/repeated samples and empty minors where applicable.

`Run.finalResidual_localMinor_schur_le` now applies the Fischer bridge to the
**actual next residual** after an arbitrary successful prefix. Only the next
pivot must be absolute-column-maximal on the augmented sample. Previous
pivots need not be complete or local, and the matched old row can differ from
the newly pivot-deleted row. This is not yet the iterated initial-budget or
band-pivot counting theorem.

The specialization covers `exp(t*y)/(1+exp(y))` on all real coordinates.
`Run.negateColumns` transports any run of the original fermionic kernel
without reselecting pivots, and its final residual is exactly the original
residual at the inverse frequency reflection. Thus neither signed frequencies
nor time endpoints are excluded. No canonical tie-policy equivalence is
claimed: the chosen pivot sequence is transported.

Independent regressions check 25 Vandermonde block identities, including
unsorted/repeated nodes; 3,825 exact fresh minors along complete, smallest-
pivot, and interior-pivot prefixes; and 350 exact next-residual matched-rank
bounds. The matrix entries are exact values of the actual fermionic family
at `t=i/5`, `y=5*j*log(2)`, including zero/one times and both frequency signs.
Nineteen residual minors agree with these exact values at independent
160/320-digit precision runs, to 45 significant digits. These decimal checks
are evidence only; the all-real structural claims are proved in Lean.
At DD delivery the initial variable-degree local determinant budget, its iteration/pivot
counting, a cutoff-uniform G1 bound, and Problem 4.2 remain open.
Local root verification passes a 3,614-job Lean build and 146 tests, with
only the accepted public theorem axioms and no known dependency vulnerabilities.
Push and PR CI pass; the change merged as PR #110 on 2026-10-19.

Phase DE proves the scale-free **initial** local determinant budget. Each
positive central/dyadic panel has an approximation using `8p` separated
terms and uniform error `2^(-p)` on physical times. Negative panels follow by
the existing fermionic reflection; zero is assigned to the positive central
panel. The identity with one active branch of `dyadicSeparatedTerms` reuses
the actual-exponential error theorem rather than introducing a new remainder
estimate. The number of features does not grow with the panel index.

For any order `n>=8p`, `fermionicKernel_panel_sample_det_le_two_pow` proves
`abs(det K_n) <= 2^n * n! * 2^(-p*(n-8p)) * 2^(8p)` for arbitrary sample rows
in `[0,1]` and columns in that one signed panel. The generic determinant
bound retains its factorial; it now depends only on local minor order, not
the global cutoff scale. Odd-order exponent cancellation and the explicit
choice `m=16+2q` prove
`fermionicKernel_panel_minor_det_le_dyadicBudget`: at `n=16*(33+4q)`, every
such minor is bounded by `(2^(-q))^n`, uniformly over all panel indices and
both signs. Tuples need not be sorted or injective. This is a theorem for
the actual kernel, not for its Taylor approximant.

Independent checks cover 345 exact masked-feature identities with rational
column weights, 435 actual-exponential approximation samples at 160/320
decimal digits, and twelve small-order actual determinant samples at 320/640
digits. Canonical quantities in the precision comparisons are decimal strings
with 45 significant digits. Exact integer tests check 65 exponent identities,
1,001 scalar accuracy conditions, and seventeen full factorial-budget
inequalities. A finite boundary-membership regression checks both frequency
signs and zero without claiming the universal partition theorem. Large
target-order floating-point determinants are not used to certify the theorem.
Backward transport of nonzero minors, the panel-local pivot count, the global
partition/count assembly, G1, and Problem 4.2 remain open.
Local root verification passes a 3,615-job Lean build and 151 tests, using
only the accepted public theorem axioms and finding no known dependency
vulnerabilities.
Push and PR CI pass; the change merged as PR #111 on 2026-10-19.

## Phase DF — backward witnesses and the global finite rate

Lean proves `fermionicKernel_gecp_error_le_localPanelBlock`: on
`[0,1] x [-2^s,2^s]`, every exact successful complete-pivot sequence starting
from the actual kernel has error at most `2^(-q)` after
`32*(s+1)*(33+4*q)` pivots. The same guarantee holds at every later rank by
`fermionicKernel_gecp_error_le_of_localPanelBlock`. The extra logarithmic
startup overhead of Phase X is absent. Boundary times, signed frequencies,
and arbitrary choices among complete ties are included.

The mechanism is a backwards minor witness, not positivity of separate
conditioned features. `GECP.markedPivotProduct_le_initialMinor` constructs an
actual initial minor dominating any marked pivot product. Unmarked steps
use a matched-rank row deletion that keeps the local columns and physical
row domain; marked steps gain the actual pivot factor. Selecting exactly
`n=16*(33+4*q)` marks allows DE's fixed-order budget to rule out `n` large
pivots in a panel. A proved signed partition has `2*(s+1)` panels, completing
the count. The generic pivot-count theorem assumes an initial budget, but
the actual fermionic rate discharges it with DE; no residual envelope or
convergence conclusion is assumed.

Independent tests check 64 exact mark sets, 135 independently enumerated
initial budgets, signed boundaries, 780 count identities, and three actual
exponential determinants at independent 160/320-digit precision. This is
finite evidence, not a runtime continuous-domain certificate. Root verification
passes 3,619 Lean jobs, accepted public axioms, 155 tests, and no known
dependency vulnerabilities. Push and PR CI pass; the change merged as PR #112
on 2026-10-19. The real-cutoff/accuracy corollaries and full G1 /
Problem 4.2 completion audit remain open. In particular, dyadic-domain
completeness cannot simply be assumed for a trace on a smaller rectangle.

## Phase DG — actual real cutoffs and the G1 parameter boundary

`fermionicKernel_gecp_error_le_logarithmic` proves the actual original-kernel
bound, for every `Lambda>=1` and natural rank `k`,

`abs(R_k(t,omega)) <= 1024 * exp(-((log 2)^2/256)*k/log(1+Lambda))`.

The complete maxima are on the actual cutoff rectangle, not a larger dyadic
domain. Initialization and exact successful updates construct the finite runs;
there is no independent realized-prefix, determinant-envelope, or convergence
assumption. The current interface assumes an exact successful complete-pivot
sequence for all natural indices. Both signs, physical endpoints, and choices
among complete ties are included. Constants are explicit, not sharp.

`fermionicKernel_gecp_error_le_allCutoff_logarithmic` replaces the denominator
by `log(2+Lambda)` for every positive cutoff. The real-tolerance corollary gives
the explicit ceiling rank
`ceil(256/(log 2)^2 * log(2+Lambda) * log(1024/epsilon))`, with natural ceiling
truncated at zero, for any positive tolerance.

The literal all-positive-cutoff interpretation of G1 is **refuted**, not
silently repaired. For every `C,c>0`, the formal obstruction chooses
`Lambda=min(1,c^2/(16*C))` and proves that the actual first complete pivot
`(0,Lambda)` leaves residual at `(1,-Lambda)` strictly larger than
`C*exp(-c/log(1+Lambda))`. The first-pivot residual is at least `Lambda/2`,
whereas that proposed bound is at most `Lambda/4`. This does not contradict
the proved `Lambda>=1` theorem or the intended large-cutoff scaling.

Independent evidence includes 6,400 exact count cases, 53 exact cutoff-cover
choices, and 160/320-digit comparisons of 212 envelope values, 28 tolerance
ceilings, and nine actual first-pivot witnesses. Zero/nonfinite witness values
are rejected; numerical evidence is not a continuous-domain certificate.
Root verification passes 3,622 Lean jobs, the permitted public axioms, all
159 tests, and no known dependency vulnerabilities. Push and PR CI pass;
PR #113 merged on 2026-10-20 after retrying one toolchain-download TLS failure.
The separate full-SPEC completion audit remains pending.

## Phase DH — finite GECP runs and physical cross error

The rate no longer assumes an infinite successful sequence.
`fermionicKernel_finiteRun_error_le_exp` takes only a finite actual
`Run fermionicKernel` complete on its physical rectangle, `Lambda>=1`, and
physical query coordinates. Its successful count is `run.rows.length`; the
bound is `1024*exp(-((log 2)^2/256)*count/log(1+Lambda))`. No residual
sequence, independently realized prefixes, later nonzero pivot, or future
completeness is assumed. The dyadic block, all-positive-cutoff `log(2+Lambda)`,
and real-tolerance corollaries retain DG's explicit counts/constants.

Finite prefix extraction and bounded-index backwards witnesses provide the
proof. The last successful update controls the final residual, even when no
further successful step is available. A zero-valued actual complete maximum
on the final residual certifies exact domain termination. This is not a
floating-point-zero or certification-budget heuristic, and it never divides
by a zero pivot.

`Run.finalResidual_eq_sub_crossApproximation` identifies iterative error with
the original selected-core cross. The new Green theorem bounds the actual
difference of the original and approximated integrals. Integrability of both
terms is derived from continuity, `rho` integrability, and the bounds; it is
not an additional application premise. Signed densities and any measure
supported almost everywhere on the cutoff interval are allowed.
The interval-integral corollary discharges the support premise with the
restricted physical measure.

Independent exact checks cover all 127 prefix/mark witnesses (642 backwards
steps), 504 selected-cross entries, 84 signed Green-sum identities, and two
rank-deficient exact stops. Another 258 actual-exponential determinant/cross
values agree at 160/320 digits as 45-significant-digit strings. All numerical
Green sums are also recomputed directly at both precisions (42 per precision),
against iterative residual sums and exact references. All numerical
crosses use solves; nonfinite/zero core determinants are rejected. Finite-grid
evidence is not a continuous pivot certificate. Root verification passes 3,625
Lean jobs, the permitted public axioms, 163 tests, and no known dependency
vulnerabilities. Push and PR CI pass; PR #114 merged on 2026-10-20. The
separate full-SPEC audit remains pending. Approximate-pivot robustness and
sharper constants are not claimed.

## Phase DI — original-scope audit and bounded closures

The [completion audit](COMPLETION_AUDIT.md) records 33 supported and seven
incomplete requirement groups. The main finite exact-GECP rate remains proved;
the full original SPEC is **not complete**. Finite PSD canonical traces do not
by themselves prove the continuous SPD baseline, and a Lipschitz residual-
diagonal premise does not establish regularity from the original kernel.
The required census depth/geometry/normalized curves, result metadata outside
GECP, and listed sparse-method comparisons also remain incomplete.

DI closes D1's real quantifiers with an actual separated witness of length
`<=128*log(1+Lambda)*log(128/epsilon)` for `Lambda>=1`, `0<epsilon<1`.
It also derives compact complete-maximum existence for each finite residual,
and certifies a concrete exact two-pivot grid and its minor/pivot signs.
These are proved corollaries, not status-only aliases.

The axiom gate now fails closed on disallowed/missing/malformed reports and
preserves a failed Lean process exit. Clean packaging excludes the local
duplicate note files; a rebuilt wheel passes isolated import/GECP smoke.
No artifact was published and original note files were not changed.

Known continuous pivoted-Cholesky/P-greedy convergence is prior work, including
[Jeong and Townsend's historical v1 Lipschitz-kernel result](https://arxiv.org/abs/2509.13582v1).
The repository's conditional residual-Lipschitz lemma is not claimed to
formalize that full published result. Root verification passes 3,627 Lean
jobs, 838 enforced axiom reports, 177 tests and no known dependency
vulnerabilities. Push and PR CI pass; PR #115 merged on 2026-10-20. The seven
gaps remain required work.

## Phase DJ: original-kernel geometry

**Proved:** `KernelPosSemidef` tests mathlib PSD on every finite sample,
including repeats, with no finite ambient-domain assumption. Every actual
successful diagonal run preserves PSD of both `R` and `K-R` through
`diagonalRun_residual_and_removed_posSemidef`. The latter invariant is not
assumed: each update removes a positive rank-one kernel. It yields
`diagonalRun_power_le_original_energy` at every selected node `s`:
`R(x,x) <= K(x,x)-2*K(x,s)+K(s,s)`.

`diagonalRun_error_le_original_fillDistance` then proves
`abs R(x,y) <= 2*L*h` from coordinatewise Lipschitz regularity of the original
kernel and a selected-node cover of radius `h`. There is no residual-Lipschitz
premise, no maximality assumption, and no rank-dependent constant.
`diagonalRun_exists_diagonal_completePivot` proves availability of a diagonal
complete maximum on every nonempty compact continuous kernel domain, including
the zero-maximum stopping case. `linearFeatureRun_residual` computes `x*y`
after the actual pivot at zero of `K(x,y)=1+x*y` on the infinite real domain.

**Exact independent evidence:** two arbitrary diagonal orders, six prefixes
each, 510 selected-node energy comparisons, all principal minors of repeated
five-node samples of residual/removed kernels, interpolation, exact rank-five
termination, and original-kernel Lipschitz/fill inequalities. A constant
rank-one PSD kernel permits nonzero off-diagonal complete ties. A two-node
triangular-kernel example gives `R(1/4,1/4)=7/16>1/4`, so the factor two in
the generic `2*L*h` theorem cannot simply be dropped to one.

**Not claimed:** every weak-PSD maximizer is diagonal (the tie example refutes
that); strict-SPD fresh-residual necessity (R08 still incomplete); published
asymptotic P-greedy rates or literature priority. DJ closes R09, so the
40-group audit has 34 supported and six incomplete groups: 85% checklist
coverage, not a weighted effort estimate. Full verification/delivery evidence
is maintained in README and RESEARCH: 3,630 Lean jobs, 857 enforced axiom
reports, 180 tests and no known dependency vulnerabilities; green push and
PR CI, merged PR #116 on 2026-10-20.

## Phase DK: strict SPD through actual complete runs

**Proved:** `posDef_diagonalSchur_fresh` expresses the fresh Schur matrix as
`Eᴴ*M*E` with a rectangular elimination map whose lower block is the identity.
Its `mulVec` is injective, so mathlib transports strict positive definiteness.
`StrictKernelOn.residualUpdate` lifts this to every injective fresh sample.
Strict two-point determinant positivity rules out an off-diagonal nonzero
complete maximum. Previously selected rows are zero, so successful later
pivots must be fresh; neither fact is assumed as a residual invariant.

`completeRun_diagonal_of_posDef` therefore proves diagonal necessity for every
actual finite exact complete run of an SPD kernel on any ambient type, with
no finite-type or diagonal tie-preference assumption. Mathlib's
`Matrix.PosDef K` is the original-kernel strict finite-support quadratic-form
condition. `posDef_completeRun_iff_pivotedCholesky` equates those runs with an
independently defined recursive diagonal-maximization rule. Residual PSD,
strictness on fresh samples, squared diagonal domination, compact continuous
complete-maximum existence and original-kernel `2*L*h` fill control follow.
The fill corollary no longer takes diagonality as a supplied premise.

**Independent exact evidence:** all complete-pivot branches and prefixes of a
five-node rational feature-Gram sample and a three-node identity sample;
every fresh principal minor, strict two-point inequalities, rectangular
congruence, interpolation and exact stopping. All six identity-kernel tie
orders are retained. Lean additionally certifies the existing two-pivot real
kernel as an injective Gram factor and applies the general equivalence to it.

**Not claimed:** weak PSD forces diagonal necessity (a signed rank-one
off-diagonal tie refutes that), diagonal maximizers are unique, approximate
floating-point pivots satisfy exact maximality, or the full published
asymptotic P-greedy rates have been formalized. R08 is now supported at its
original staged scope. The audit has 35 supported and five incomplete groups,
87.5% checklist coverage rather than weighted remaining effort.
Full verification passes 3,633 Lean jobs, 876 enforced axiom reports, 183 tests
and no known dependency vulnerabilities. Push and PR CI pass; merged as
PR #117 on 2026-10-20. The five computational gaps remain required work.

## Phase DL: complete result metadata and numerical failure honesty

**Implemented, not a new Lean theorem:** typed numerical metadata accompanies
Cholesky/sparse results and schema-2 synthetic summaries; GECP exposes a view
of recorded evidence and approximation constructions report non-iterative
scopes. Actual precision, requested precision, conditioning stage/status,
ties, certification and error scopes/settings are explicit. No estimate is
silently replaced by a made-up finite value or high-precision claim.

Cholesky uses the actual returned reconstruction to decide convergence, with
an explicit roundoff/symmetry-limited stop. Stable symmetrization avoids
overflow for finite large SPD input. Sparse supports no repeated exhausted
atoms, returns zero atoms for an already-small signal, labels training error,
and separates optimizer termination from returned-fit convergence. Worsening,
nonfinite or failed refinement retains the prior fit. Empty/non-applicable,
rank-deficient, failed-SVD and budget paths have independent tests.

Independent review found and corrected a false `estimated` status for empty
conditioning, an incorrect float64 GECP default tie threshold in the metadata
view, and overflow in post-update Cholesky symmetrization. Tests independently
recompute returned errors and singular values, exercise optimizer failure and
rejection, and enforce deterministic finite JSON. Historical application
artifacts and all Lean theorems are unchanged. R24 is supported; the 40-group
audit has 36 supported and four incomplete groups (90% checklist coverage,
not weighted effort). Census depth/diagnostics/plots and sparse comparisons
remain required work; no numerical evidence is relabeled a formal proof.
The full gate passes 3,633 Lean jobs, 876 enforced axiom reports, 232 tests
and the clean dependency audit after targeted certifi/pillow lock updates.
PNG and isolated wheel smoke checks pass; push and PR CI pass; PR #118 merged
on 2026-10-22. No artifact was published.

## Phase DM: full census evidence (verified; merged PR #119)

**Implemented numerical interface, not a new theorem:** schema-2 records add
complete sampled residual/core histories and physical boundary, signed dyadic
panel and nearest-reflection diagnostics. Prefix metadata uses the requested
tolerance while retaining the source trajectory's strict tolerance. Failed
adaptive certificates cannot imply convergence. Normalized plotting rejects
explicit adaptive/non-sampled or unsuccessful runs instead of mislabeling them.

Independent exact Fraction geometry checks and manually specified prefix,
zero-rank, budget, malformed-history and failure cases pass. Review corrected
floating-point frequency endpoints that could overshoot non-integer physical
cutoffs; byte-level regressions preserve all seven canonical integer grids.
Schema-2 block analysis uses actual sampled maxima, not near-maximal pivots,
and rejects negative residuals rather than silently taking their magnitude.

**Primary observations:** all 28 schema-2 cases converge at the seven cutoffs
`1,10,100,1000,10000,100000,1000000` and four tolerances down to `1e-12`.
The strictest-tolerance ranks are respectively `9,17,37,63,91,119,146`.
Two full runs from clean tracked generator revision
`c0e8c5093134f369b93781c413750b7f1de46036`, configuration hash
`c0ac1bbaa7ddb39a2bb4a1c306985386f89f4103ff75cb4e21f99257e1e9f31b`,
are byte-identical. The JSONL SHA-256 is
`c6854330fcb74007f8f9814fba48c314d34fe26c1cef5461d908ea5182b66bda`.
All 21 overlapping historical prefixes preserve pivots and ranks exactly;
legacy files are unchanged. The new rank and normalized-residual plots have
been rendered and inspected. Normalized decay curves have related scales,
not an exact common exponential law; this dataset has no zero residuals.

**Independent numerical checks:** 192-bit selected-core recomputation at
rank one, a middle prefix and the terminal prefix checks 21 cores. Relative
differences in determinants, signed pivot products and smallest singular
values are below `9e-27`. Representative stored crosses at cutoffs
`1,1000,1000000` have held-out maximum errors approximately
`6.31e-14,6.33e-13,8.61e-13` on `36,3960,15120` point pairs. Elimination uses
128 bits, but grid coordinates originate in float64; both scopes are recorded.
These checks use point arithmetic and factored solves, not interval bounds.

**Representative resolution observations:** all three independently rerun
order-32 trajectories converge at `1e-12`. At cutoffs `1,1000,1000000`, their
ranks are `9,63,146`, equal to the order-24 ranks; their sampled terminal
residuals are approximately `1.78383e-13,9.59729e-13,9.71664e-13`. Grid shapes
are `63×63,621×683,1279×1303`. These records have clean generation revision
`ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61`, configuration hash
`bf49e21ce208a1336bfc83d3e2ec955adf627e64e8056c8be9986d855ee2d3be`,
and artifact SHA-256
`7251cbb528179cc5cbfa6eb5635b8fa280e087d1e04db1348b56a9065026d588`.
The paths differ: exact decimal selected-node common prefixes have lengths
`3,0,0`, and all three order-32 terminal sampled residuals are larger. Both
generation revisions/source hashes remain in the combined report. Equal ranks
do not establish grid independence; only three cutoffs and two grid orders
were compared.

The combined report passes all 30 selected-core checks at 192 bits, including
nine resolution cores. Relative determinant, signed-pivot-product and
singular-value differences remain below `9e-27`. The full branch root gate
passes 3,633 Lean jobs, 876 enforced axiom reports and 430 tests (139.64
seconds), with clean formatting/mypy/dependency audit and green push/PR CI.
PR #119 merged on 2026-10-22; DM supports R25–R27. Exact merged/pushed
closure-head gates passed at `8fb6b03`. At DM closure the audit recorded
39 supported and one incomplete group, 97.5% checklist coverage rather than
weighted effort. DN subsequently verifies R31's comparisons; the final
forty-group audit remains required. No
continuous-error certificate, exact-reflection theorem, grid independence,
optimal tie behavior or new approximate-pivot rate is claimed.

## DN: six-method spectral comparisons (verified, merged PR #120)

The R31 implementation adds the Matsubara rational bridge, SciPy AAA-derived
physical atoms, genuine iterative vector-fitting pole relocation, and
data-ranked strided-Hankel Prony/ESPRIT. The shared-observation comparison
retains every tested budget, including complex/out-of-band pole rejection and
failed fits; no generating order or exact rational truth enters an adapter.
Primary dictionaries are fixed independently of off-grid truth. The
known-library control is separately prior-informed. The frozen 192-fit JSON
and reviewed PNG repeat bytes from clean `072ce10`; every measured error is
independently reconstructed. At `1e-8`, the clean four-atom case needs four
tested atoms for AAA/VF/Prony/NLS versus universal rank 12. VF/Prony/NLS also
pass with four on the noisy case and six on each fixed-quadrature continuous
density reference. AAA fails the other three physical fits; those failures
and nonphysical poles remain explicit. No method is supplied the true order.
These are smallest tested supports, not optimality/uniqueness or continuous
error claims. R31 is supported by the actual comparisons, independent error
reconstruction, 513-test full root gate, green exact-head push/PR CI and
installed-wheel checks (including minimum SciPy 1.15). PR #120 is merged;
the final independent original-SPEC audit remains required.

## Attribution update (DO)

The [workshop problem list, v3](https://arxiv.org/html/2602.05394v3) records an
August 2026 claim that Venkata Siddharth Pendyala's *Local-to-Global Convergence
of Greedy Cross Approximation for Totally Positive Kernels*
([DOI](https://doi.org/10.5281/zenodo.21863274)) solves Problem 4.2. This audit
read the workshop's update, not the inaccessible underlying paper, and does
not independently certify that claim. No priority for the repository's rate
is asserted. The [current v3 Cholesky paper](https://arxiv.org/abs/2509.13582v3)
is by Sungwoo Jeong, Alex Townsend and Christopher Wang; older phase records
refer specifically to the two-author v1. The repository formalizes the staged
fill-distance baseline, not all of that paper's stronger asymptotic results.

## Phase DP: genuine residual-adaptive census (verified; merged PR #122)

**Observed, not proved:** the new explicit-axis
[`adaptive_census.json`](experiments/data/adaptive_census.json) contains all
seven cutoffs and four tolerances `1e-6,1e-8,1e-10,1e-12`. Every case and all
28 prefixes report convergence. Refinement responds to odd-eighth residual
and interpolation-surplus probes, selectively bisects marked panels together
with reflection mates, and retains every previous binary64 coordinate.
Two successive strictest-tolerance stage passes are required. This is actual
adaptive-grid evidence, not a relabeling of the historical fixed-grid census
or the separate interval-pivot API.

| Cutoff | Actual stages | Initial → final training grid | Chosen ranks at the four tolerances | Final sampled error at `1e-12` |
|---:|---:|---:|---|---:|
| 1 | 3 | 9×9 → 25×25 | 6, 7, 8, 9 | `3.488e-13` |
| 10 | 2 | 33×41 → 41×57 | 12, 13, 15, 17 | `4.528e-13` |
| 100 | 3 | 57×65 → 109×89 | 24, 28, 33, 36 | `6.418e-13` |
| 1,000 | 4 | 81×89 → 122×121 | 38, 48, 55, 64 | `3.659e-13` |
| 10,000 | 4 | 113×121 → 170×169 | 53, 66, 80, 97 | `9.792e-14` |
| 100,000 | 5 | 137×145 → 238×217 | 68, 87, 103, 124 | `2.932e-13` |
| 1,000,000 | 6 | 161×169 → 258×257 | 83, 103, 126, 151 | `2.175e-13` |

The table's errors are rounded summaries of the controller's MP192 final
selected-cross checks on six seventh-fraction points per final leaf panel.
Those coordinates are untouched until every refinement and prefix choice is
complete. They do not tune a grid or rank. The first training-passing rank
can fail the probe/surplus test: at cutoff `1e6`, the strict prefix advances
from 147 to 151. The chosen ranks are first jointly passing tested prefixes,
not globally minimum ranks or a claim of grid independence. Full canonical
decimal strings, stage axes, panels, scores, decisions, budgets, pivot/core
histories and physical geometry are retained in the JSON.

Both fresh CLI generations use clean revision
`3cee4186fddf990abdd33dbe24c46ff908b1923f`, configuration SHA-256
`ac47ec5146f5dc50916513917518d2bbfcb31da622b1f659225568cc48d641a1`,
and produce identical bytes with data SHA-256
`730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c`.
The frozen study uses nine stages maximum and generation margin `0.01`
(`1e-14` trajectory tolerance), following explicitly labeled failed pilots;
library defaults remain five stages and margin `0.25`. Coordinates have
53-bit repr-decimal semantics, elimination uses 128 bits, and sampled cross
checks use 192-bit factorizations/solves, never an explicit inverse.
Source/lock hashes and runtime versions are recorded. A fresh run on a later
commit records that later revision and is not promised to reproduce the
historical JSON bytes, even if the numerical paths agree.

The normalized plot has been rendered and visually inspected at 2240×1440,
with seven curves and 28 final-validation markers. One sampled zero is omitted
from logarithmic coordinates without substituting a positive floor. Its
byte-identical repeat has SHA-256
`dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6`.
The independent MP192
[`adaptive_census_checks.json`](experiments/data/adaptive_census_checks.json)
passes full seven-cutoff/28-prefix coverage, all revision-blob provenance
checks, 21 selected-core determinant/pivot-product/singular-value checks and
all 28 complete final sampled crosses. Every final cross is within its
requested tolerance. The report has SHA-256
`4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03`.
These are independent higher-precision point/factored-solve checks, not
interval enclosures.

At final branch head `3c38183e5c2a7d2859867bb7acff1353ee98fb76`, the
exact full root gate passes **713 tests in 150.98 seconds**, 3,634 Lean jobs,
all 876 axiom reports, 95 formatted files, mypy on 21 sources and the clean
dependency audit. Source/wheel archives include all three DP artifacts and
tests while excluding protected note copies, caches and scratch; all 21 wheel
Python files hash-match source and the packaged README is checked. Actual
outside-checkout installed-wheel checks pass with pinned Python 3.12.12,
NumPy 2.5.2, SciPy 1.18.0 and mpmath 1.4.1, and with minimum SciPy 1.15.0 /
NumPy 2.4.6. Each environment exercises all 45 exports, adaptive success and
budget failure, and 48 six-method fits; final archives and both environments
are rebuilt and rechecked at that branch head. Exact push CI `40017937659`
and PR CI `40017945103` pass all four checks. Earlier development and
code/artifact-head roots passed the same 713 tests in 151.22 and 150.40
seconds respectively; those results are not relabeled as the final branch
run. PR #122 merged at 2026-10-23 20:20:48 UTC as
`3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`, and local main is fast-forwarded
to origin. Canonical generation revision `3cee4186` is confirmed as an
ancestor, not replaced by the merge revision. Exact pushed-main closure
metadata/gates repeat before handoff. The next fresh full original-SPEC audit
remains required; whole-goal completion is not claimed. No sampled maximum
or surplus is a continuous supremum bound, and no approximate-pivot rate
theorem is added.

## DQ completed exact-main closure and DR documentation boundary

The DQ closure revision `e8c35eda7662566639c5eac6748d3032df87075f` passed
757 tests in 150.31 seconds, the 3,634-job Lean build, 876 permitted axiom
reports, 22 public-import examples, style/type/security and CI `40051138962`.
Clean source/wheel inspection and both outside-checkout installed environments
passed; local and origin main were equal with no open PRs. The original goal
was marked complete after these checks. This supersedes earlier prospective
closure language in the dated DP/DQ history, not any theorem's assumptions.
DR synchronizes documentation/provenance and adds documentation regressions;
it does not add a mathematical result, change an algorithm, regenerate
historical data or reopen the completed research scope.
Exact DR source `51fa6b3` passed 759 tests in 144.90 seconds, all established
proof/type/style/security gates, CI `40054146369`/`40054160333`, clean archive
inspection and both installed environments. PR #124 merged on 2026-10-23
at 22:30:40 EDT as `a421f7c8ed8b1b8754f38116939e26181636fe82`.

## DT practical-application baseline

The matched [DLR/IR pilot](experiments/data/dlr_benchmark_pilot.json) is a new
application experiment, not a reopening of the completed original theorem scope.
GECP, pydlr and sparse-ir receive the same physical forward observations and
common cached SVD fitting. Eight of nine method/case records meet the sampled
`1e-8` target; float64 GECP at cutoff 1,000 retains its explicit arithmetic-range
failure. Independent 100-digit physical pole sums verify every successful
workload/error flag. [APPLICATION.md](APPLICATION.md#recorded-dt-pilot) records
the frozen generation, ranks, limitations and next decision.

The observation supports investigating cached-basis reuse and numerical-range
handling before a Rust rewrite. It does not establish a fastest algorithm,
minimal rank, stable spectral recovery, continuous certificate or production
materials-solver result. Timings and process RSS are single-machine pilot
measurements, not deterministic numerical artifacts. Final verification and
delivery are recorded in [PR #126](https://github.com/tripp-smith/gecp-kernel-structure/pull/126).

## DU reusable physical basis and arithmetic policy

`FermionicBasis`, `FermionicFit` and `LeastSquaresPlan` promote reusable
solve-based operations into the installed package without changing the existing
research API defaults. Their arrays own immutable data; bounded per-instance
caches retain a core factor, evaluation designs and SVD fit plans. Complex
matrix batches preserve the caller's observation axis. The physical Fourier
transform uses `exp(+i*wn*tau)`, `wn=(2*n+1)*pi/beta`, and ordinary integer `n`.
Independent analytic pole sums test both beta scaling and the Green-function
minus sign, not just consistency with another package helper.

The valid-input cutoff-1,000 float64 failure is not a language-speed problem.
The new explicit auto policy records the arithmetic failure and retries at
reference precision; the original `gecp` behavior is unchanged. A beta-20,
wmax-50, order-six, `tol=1e-8`, maximum-rank-96 check selects rank 45 with
sampled residual `8.295898986e-9`. Reference selection, kernel evaluation and
cached core solves use 128 bits; physical coordinates, design exports, SVD and
final products remain binary64. Nodal cardinal error is `1.68e-31` before the
final coefficient products. Independent 160-bit asymmetric three-pole truth
checks the exported complex matrix fit at 117 time points and eight positive/
negative Matsubara indices: maximum errors `3.82e-9` and `1.23e-9`. This also
passes under NumPy's strict traps without changing the caller's settings.
These are sampled fixture observations, not a certificate for arbitrary inputs.

The smaller runnable tutorial gives time/frequency errors approximately
`2.54e-12`/`2.85e-13`, but its double core has estimated condition `2.31e10`
and nodal cardinal error about `3.04e-7`. Thus a benign smooth-function fit
does not prove arbitrary coefficient stability. Neither a finite-grid stopping
criterion nor high-precision construction supplies an end-to-end accuracy
guarantee. The finite-precision/inexact-pivot bridge remains an open product
requirement, including coefficient, transform and solve error propagation.

Versioned NPZ/JSON persistence validates conventions, numerical shapes, byte
budgets and reconstruction work without pickle or pivot reselection. Loaded
construction metadata is unverified history and cannot grant certification;
training diagnostics are not restored as fresh evidence. Distribution gates
inspect the actual wheel/sdist and exercise installed APIs and typing outside
the checkout. Cross-platform CI results, publication, production comparative
value, a solver loop and the measured Rust decision remain separate gates in
`PRODUCTIZATION.md`. No new mathematical theorem or speedup is claimed here.

## DV pilot: sampling observability is not solver convergence

The new physical comparison remains in progress. Its initial solver pilot
withheld several low Matsubara rows and reached a small observed fixed-point
defect while giving unseen time/frequency errors around `1e-5`. All three
providers exhibited the issue. More iterations, real-coefficient fitting,
particle-hole symmetrization, truncation changes and endpoint constraints did
not repair the missing information. The low-frequency-hole witness is retained
as an executable regression, not relabeled as successful validation.

Keeping the informative low modes makes the GECP projection design well
conditioned on a scout (estimated condition roughly `13` rather than `1e10`)
and gives independent physical errors around `5e-10`. This illustrates the
connection to observability/experimental design: a validation split can remove
essential information from an inverse projection. It does not establish a new
GECP theorem or justify accepting a small sampled defect as global accuracy.
The final protocol retains low-frequency observations, separately checks their
physical truth after the loop, and uses disjoint time/middle-high-frequency
holdouts. Forward time fitting also has fresh unseen low-frequency checks.
The old pilot points are now calibration evidence; final time fractions and
frequency lists are frozen anew before recorded generation. Detailed failed
hypotheses, conventions, acceptance criteria and the numerical-bridge planning
checkpoint are in `RESEARCH.md`. Actual results and delivery gates will follow.

Frozen cohort A now provides actual evidence:11/12 records pass the shared
final1e-8 criterion, and all three physical solvers pass. The retained sparse-ir
cutoff 10 failure has frequency error1.208646949e-8 despite passing calibration;
an independent 95-digit pole calculation confirms it. Native truncation eps and
finite calibration samples do not establish an external absolute-error bound.
Both runs agree on numerical fields, including the failure. Protocol B strengthens
calibration and reserves new validation points; no target is loosened and no
claim of fully matched success is made from A alone.

Cohort B now supplies the matched success: all 12 records,36 complex matrix
forward workloads and three physical solver loops pass the unchanged1e-8
target on fresh validation. Independent 95-digit checking confirms every forward
error/flag, all solver frequencies and selected time integrals. Both runs match
numerical fields, not machine-dependent timings. Maximum forward time/frequency
errors are 5.45e-10/7.97e-10; solver time error is at most 5.74e-10. Details,
source/measurement scopes and rank/native-eps differences are in
`APPLICATION.md` and `RESEARCH.md`.

This establishes bounded practical comparative evidence, not a performance
advantage: the cutoff 1,000 GECP reference constructor takes about19.7 seconds
in the recorded run, versus0.0194/1.337 seconds for pydlr/sparse-ir. The methods
meet one external ceiling but do not have identical errors or native arithmetic.
The next optimization must be measured, and any Rust prototype must preserve
numerical behavior. A standalone native-Sparse-IR recipe and concrete upstream
proposal were present at DV's freeze; external adoption and publication had
not yet occurred. The later DY record below supplies actual publication and
outreach evidence.

## Run provenance

Agent model/mode metadata, elapsed time, aggregate token usage, and the
non-billing cost estimate are recorded separately in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They are
operational provenance and are not mathematical or numerical evidence for any
claim in this ledger.
## DW — numerical discrepancy and observation design

Two conditional Lean results now connect explicit numerical budgets to useful
claims: `Discretization.approxPivot_of_uniform_discrepancy` transfers pivot
quality on the same domain; `Discretization.dot_error_le_of_component_bounds`
combines complex design, coefficient and final-dot errors. A sharp exact
two-point witness attains quality `765/1027`. These theorems do not verify an
interval arithmetic engine or extend a finite grid to a continuous domain.

The implementation independently encloses exact physical selected-core and
full-rank time least-squares problems, then compares them to actual public
fits and outputs. This supplies actual coefficient/core/solve budgets, rather
than assuming the perturbation hypothesis of a scalar lemma. Independent
rational, QR and bordered-determinant checks pass in development. Canonical
frozen generation and final delivery remain the phase's acceptance gates.

A large-cutoff development probe exposed a useful distinction. At beta 20,
wmax 50, tolerance `1e-9`, the rank-49 basis fitted to 61 uniform time points
retains only 39 singular directions in the actual SVD. At 256 interval bits,
the Gram pivot enclosure is inconclusive. At 512 bits it encloses the full-rank
reference, but that reference has enormous coefficient/frequency discrepancy
from the deliberately truncated fit. More arithmetic precision cannot turn
this ill-conditioned full-rank problem into a sensible fitting target.

Resolving endpoint information changes the problem conditioning: the shared
endpoint mesh retains all 49 directions with an estimated design condition
about 277, versus about `1.69e16` on the uniform mesh. The selected-node default
has condition about one. Both alternatives yield useful small numerical budgets
at 512 interval bits. This connects the bridge to experimental design and
regularization: a small physical error, a small equation defect and fidelity to
an exact unregularized least-squares reference answer different questions. The
frozen six-case protocol retains the inconclusive and large-discrepancy cases
alongside the useful alternatives; it does not hide them as failed tests.
The frozen repeat confirms these distinctions. Both complete JSON artifacts
match exactly. All five completed enclosures pass the independent coefficient
and prediction audit; one solve remains inconclusive. Large-cutoff endpoint
sampling gives numerical budgets below `2.50e-14` while physical fixture error
still reaches `3.38e-8`, above a `1e-8` ceiling. Selected-node fitting gives
numerical budgets below `9.65e-16` and physical errors below `4.48e-10`.
The uniform-time full-rank reference remains a valid large-discrepancy case.
No finite-query arithmetic budget is promoted into a physical or continuous
error theorem. Source/hash and final delivery evidence are in `RESEARCH.md`
and [PR #129](https://github.com/tripp-smith/gecp-kernel-structure/pull/129).

Cross-platform CI caught an audit portability error: re-executing BLAS cannot
authenticate another machine's rounded result byte for byte. Cohort C instead
stores those actual results and audits their mathematical discrepancy directly.
Its complete repeat agrees, and every prior scientific/numerical result is
unchanged from B. Historical B artifacts and failed audit/CI evidence remain
retained. This is a data/provenance correction, not a relaxed accuracy check.
## DX — measured optional native evaluation

The application backend now isolates high-precision time-kernel evaluation and
cached triangular substitution behind an explicit optional Rust interface.
Selection, frequency evaluation, normalization/export policy and existing
NumPy fitting/products keep their prior implementations. Python remains the
default, and archive history cannot activate native code automatically.

Development tests cover MPFR state restoration, concurrency, resource/token
limits and independent numerical parity. Review found two concrete private-API
boundary defects before source freeze: reentrant numeric conversion could
mutate input lists and create an inconsistent handle; oversized Unicode tokens
could allocate a UTF-8 copy before the byte limit. Exact builtin input checks
and pre-copy character limits fix both, with retained regressions. These are
engineering corrections, not new mathematical results.

The production cohort frozen at `e5f6e796893d31837545e6c7eac0918b641e3fda`
retains the backend under its predeclared 20% application/process-wall gate.
Three alternating fresh trials per backend give median application times
42.7959 s Python and 21.3320 s Rust (50.15% reduction), and process-wall times
44.7031 s and 23.3508 s (47.76% reduction). Rank is 55; all six physical
stream errors are below `2.57e-11`, against an unchanged `1e-8` target.
Independent 512-bit physical-coordinate QR checks pass; maximum paired
Python/Rust discrepancy is `4.26e-15`, below `1e-12`, but outputs are not
bitwise equal. Each backend's complete numerical arrays and QR values repeat
exactly across its three trials.

The eight-query prefix improves only 8.92%; warm cached calls differ by 0.32%
in this cohort and make zero new design calls. Construction remains about
19 s. These are one-machine observations for changing time meshes, not a
universal speedup, a physical continuous certificate or a new theorem.
Complete arrays, timings, scopes, hashes and failures are retained in
`experiments/data/rust_reference_benchmark.json.gz`; RESEARCH owns provenance.

## DY — published practical release

DX's exact pushed-main root/source/pure/native installed gates pass at
`67e4b14`; [PR #130](https://github.com/tripp-smith/gecp-kernel-structure/pull/130#issuecomment-6400942719)
owns that closure. DY delivers
[kernelgecp 1.1.0 on PyPI](https://pypi.org/project/kernelgecp/1.1.0/) and the
[GitHub release](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0)
from immutable tag/source `9a679930db326bda2456bfa91b3256f37359e4ce`.
The saved Trusted Publisher uploaded the same tested artifact bytes on
2026-10-25 at 23:23:48–23:23:53 UTC. METHODS and CITATION
keep exact-run theorems, independent numerical observations and finite-query
arithmetic enclosures distinct. This is release engineering, not a new theorem.

The [actual upstream offer](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
contains a normalized two-pole complex matrix recipe with disjoint time checks,
both frequency signs and native Sparse-IR fitting. Its sampled time/frequency
errors are `2.76169e-11` and `4.94085e-11` against `1e-8`; canonical spectral
normalization and independent Fourier quadrature tests pass. The earlier
single-atom recipe remains historical development evidence. This records the
initial outreach; EA below records the invited notebook submission and links
the live upstream review/merge state.

[Publication workflow 40244470229](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40244470229)
attempt 2 passes all six Linux/macOS/Windows × Python 3.12/3.13 clean public
wheel installs and the separate public sdist installation. Actual native
execution is checked for Linux x86_64 and macOS arm64; Windows uses the
universal Python fallback. Public JSON and simple-index file hashes match all
four original tested distributions. Complete native source, notices and the
locked dependency inventory were publicly available on GitHub before upload.
The [release receipt](RELEASE_RECEIPT.md) links exact-main verification and
records the source and file hashes.

The initial publication attempt is retained: upload, the public sdist and one
public wheel lane passed, but five postchecks exhausted the short metadata
visibility window with HTTP 404. Only failed verification jobs were rerun;
upload and source publication were not repeated. This distinguishes package
availability from JSON metadata readiness. DZ records that outcome and a
bounded readiness correction; its final repository acceptance is in
[PR #132](https://github.com/tripp-smith/gecp-kernel-structure/pull/132).
The controlled application and publication requirements are delivered at
their recorded workload, arithmetic and platform boundaries. Subsequent
upstream review/merge decisions, broader platform support and a continuous
inexact-trajectory theorem have separate evidence requirements.

## EA — invited upstream notebook submission

The [initial physical-accuracy offer](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
received a [maintainer invitation](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6405100346)
to add the example after the conditioning plot. The resulting
[Sparse-IR documentation PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10)
submits commit `bc4659e14c1a5f319fd6dcc7fa42796c7b29e2bd` against
`c94bf6ae3d3953ca9e2081a8be9998db348cdd24`. One notebook changes by 87 inserted
lines and three new cells; its original eight cells, metadata and cleared
source-output policy are preserved. The notebook SHA-256 is
`c93d98bcea0d0e1ef0b695931a78eda1001f6284d9986db20e1b2ecd94f0e16f`.

The upstream locked environment uses Python 3.13.2, NumPy 2.3.5, sparse-ir 2.1.4
and pylibsparseir 0.9.0. Isolated and full-notebook executions agree exactly:
rank 14, time error `2.7616964271004463e-11` and frequency error
`4.94080361244419e-11`, below the fixed `1e-8` sampled ceiling. Independent
Fourier quadrature, spectral normalization, matrix/index conventions and
read-only content review pass. The full book build exits zero with one
unchanged README-toctree warning. An initial local harness error came from
unnormalized list-valued notebook sources; `nbformat.reads` repaired the
validator without changing the submitted notebook.

This delivers a tested upstream submission, not upstream acceptance or a new
IR method. At submission, upstream workflow 40248110146 was `action_required`
with no jobs; [the live PR](https://github.com/SpM-lab/sparse-ir-doc/pull/10)
records later CI, review and merge decisions. The detailed contribution scope
is in [UPSTREAM_PROPOSAL.md](UPSTREAM_PROPOSAL.md), and local final delivery is
tracked in [PR #133](https://github.com/tripp-smith/gecp-kernel-structure/pull/133).
The published kernelgecp 1.1.0 tag/assets and frozen research evidence are
unchanged.

## EB production-runtime baseline

The optional provider baseline has actual Linux/macOS TRIQS, CT-HYB, MaxEnt
and CPU GP calculations, plus native Linux interacting H2 GF2 and coarse MgB2
QE/EPW calculations. All six execution lanes pass at `d5f58a1`; final-source
and merged-main acceptance is recorded in PR #134. Executed environments have
source/artifact identity records and actual-platform conda locks, with a fresh
Mac replay and independent native resolvent check. The observed-notice catalog
retains missing-license-file limitations rather than inventing classifications.

These establish runtime availability, not interacting DMFT, molecular
compression, EPW restart, covariance-aware continuation, large-GP accuracy,
whole-solver acceleration or publication of 1.2.0. Those are separate EC–EJ
gates. Setup/build/native startup failures and receipt/exporter omissions are
retained in the environment guides and phase receipt; no failed attempt was
counted as a successful provider calculation.

## EC development: observation domain and interoperability boundaries

The new source API fits supplied ordinary-index Matsubara observations with
the same cached complex SVD machinery as time observations, retaining explicit
observation-domain metadata. Coefficients remain cardinal-time coefficients.
Format-v2 archives retain that domain and the reader preserves v1 predictions;
an actual public 1.1.0 wheel generated the frozen compatibility fixture.
Source implementation and domain/provider tests are distinct from pending
phase acceptance and publication of 1.2.0.

Positive-only unconstrained complex least squares does not imply a reliable
time-domain continuation or negative-frequency prediction. A pilot design
with condition number around 1e12 fits positive observations well while its
time prediction misses by around 1e-3. The fitter therefore preserves explicit
observations and exposes rank/conditioning rather than silently enforcing a
symmetry. Real-provider exports from a validated time fit remain a separate
operation. Neither a training residual nor the existing time-fit precision
certificate becomes a frequency-fit physical certificate.

## ED development: interacting transforms and validation isolation

The deterministic scalar Bethe IPT implementation uses the reflected product
`U**2 * G0(tau)**2 * G0(beta - tau)`. An asymmetric finite-pole fixture differs
from the naive cube by more than one, while direct thermal/Fourier checks
agree near binary64 precision. Exact interacting finite-bath data differ from
IPT, as expected for an approximate impurity model; that disagreement is not
a compression error. Static Hartree cancellation and moment assumptions are
explicit in the native reference.

Controller tests separate a small mixed step, an observed unmixed defect,
and full-grid accuracy. In particular, a subset can converge while its
unobserved full-grid defect remains large. This is retained as a regression,
and final predictions are therefore evaluated from fitted coefficients on
untouched queries. Two consecutive observed defects below tolerance remain
necessary. Cached plan identities and numerical results have owning immutable
state, including failure and terminal-bath semantics.

Unfrozen pilots found passing calibration choices for GECP, pydlr and sparse-ir;
they do not establish final acceptance or a speed advantage. Independent review
identified asymmetric import timing, replay-inclusive RSS, empty-replay false
passes and pilot holdout reuse. The campaign distinguishes those scopes and
reserves new outside-primary frequencies and off-grid final times. Actual
native tests and exact finite-model unit checks pass; the canonical frozen
campaign, root command and exact-source delivery remain required.

The frozen `a3a64ba` campaign now passes all 27 compressed and nine native
closed loops across three fresh repetitions, with byte-identical numerical
arrays per method/U. Independent artifact auditing recomputes actual equations
and final errors. Maximum independent G/Sigma error is 6.87e-8 and common-bath
replay error is 8.15e-8; native refinement changes by at most 7.15e-11.
Actual finite fixtures also pass. [The report](INTERACTING_DMFT.md) retains
unfavorable timing: native Fourier is faster than all three compressed methods
on this small beta=10 workload, and Python reference GECP setup dominates.
No universal or interacting-workflow speed advantage follows. Full root,
exact-revision CI, merge and main verification remain separate delivery gates.

## EG calibration: frequency-native real-parity cross

Direct time-cardinal frequency fits of authentic coarse MgB2 output missed the
Z/energy-shift budgets and remain failed evidence. An explicit real Nambu
parity reduction gives three Green channels in squared Matsubara frequency.
For this Cauchy kernel, exact Schur-product identities permit log-magnitude
pivot selection without subtracting nearly equal residuals. High-precision
cardinal construction followed by guarded binary64 evaluation avoids the
ill-conditioned time-cardinal continuation observed in the failed pilots.

The rank-60 candidate meets the training/calibration margins on actual EPW
data. Independent direct-product arithmetic and an analytic finite-pole BCS
eigensystem check the representation. Archive reload, complex symmetry,
ownership, corruption and unsafe-conditioning tests pass. These checks do not
close untouched final-row or native-restart acceptance. The method explicitly
uses real Eliashberg parity, preserves native solver fields, and adds no generic
complex-linearity, continuous certificate or universal speedup claim.

The subsequent actual EPW run reads both native exports successfully but its
default-order real-axis Padé outputs differ far beyond the fixed budgets.
Small imaginary-axis approximation errors therefore do not establish stable
downstream continuation. A parser defect had hidden this numerical failure;
correcting header whitespace preserves the rejection. The failed inputs,
logs and complete comparison are retained. A bounded native conditioning
study must pass independent reproduction and omitted-frequency checks before
any fresh validation. Fitted-restart acceptance remains open.

## EI development: column observations and independent covariance checks

Column-oracle Cholesky separates scalable factor state from a dense covariance
matrix. Consistency checks on selected columns cannot establish positivity or
symmetry of unobserved entries, so those remain explicit callback assumptions.
An immutable factor's convergence uses its actual returned diagonal
reconstruction error, including small negative values before clipping.

Tiny diagonal noise can destroy a subtraction-based Woodbury solution even
for a well-behaved shifted system. A temporary implicit Householder/SVD path
avoids that subtraction without retaining a second full factor or building a
complete orthogonal matrix. The actual approximate-system residual, including
final RHS rescaling, must pass; this does not establish the original-kernel
CG residual. Independent dense prediction and true-kernel CG gates remain
separate from factor stopping.

Pre-validation UCI pilots show why a representation budget belongs in model
calibration: a likelihood-favored narrow kernel can miss dense predictions
at the allowed rank, while a slightly smoother candidate can meet the same
fixed error budget. Rejected candidates are retained. The accepted runtime,
all six real-data sizes, three fresh repeats, provider parity and delivery
gates were pending at that pilot; default-environment results close no EI gate.

The subsequent locked-runtime generation `203e58a` now passes all 18 real-data
trials, actual Mac GPyTorch parity and independent source/blob/array auditing.
At n16384 the method retains 69.37 MiB and uses at most 321.75 MiB fresh RSS;
the original-kernel CG residual is independently recomputed, separate from
the approximate solve. Maximum small-case normalized mean RMS error is
1.25e-4 and observation-standard-deviation RMS error 1.46e-5. [The report](COVARIANCE.md)
retains slower 512/1024 timing cases and large physical prediction error;
representation fidelity does not establish a competitive forecaster. Full
native/root/CI and delivery remain separate open obligations.

The EF molecular comparison at frozen source `56c3211` reports passing
same-solution sampled compression gates on the two H2 geometries and both
native grids: 36 fresh trials have maximum component error 1.726e-8. Its
selection requires both primary and doubled-cutoff controls; the failed
`6050ddd` run and its consumed validation points remain recorded in
[RESEARCH.md](RESEARCH.md). Independent QR/Fourier artifact auditing passes all eight declared gates;
final verification and delivery remain open.
These sampled results neither certify the continuous domain nor establish GF2
accuracy against the exact molecular problem.


### EG native continuation checkpoint

The 17 K pure fit passes held-out imaginary-axis checks, but its native Padé
restart fails. The fixed seven-order study is inconclusive with all candidates
failing, despite passing baseline and compiled-helper rechecks. A lossless
archive and truth-assisted sensitivity diagnostic are development work, not a
proof of fitted restart interoperability; no native sensitivity result or
20/30 K fitted-field score exists. See the [evidence](EPW_INTEROPERABILITY.md#development-checkpoint-native-continuation-and-archival).

The EG 20 K archive checkpoint exposes an arithmetic API boundary: equivalent
real Matsubara formulas need not give the same binary64 coordinates. Exact-node
replay is safe here while a three-ULP neighbor fails the unchanged conditioning
guard. The verifier repair preserves this rejection and does not establish
continuous-frequency stability or native solver acceptance. See
[EPW_INTEROPERABILITY.md](EPW_INTEROPERABILITY.md).
