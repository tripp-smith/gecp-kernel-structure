# Final Handoff Report

## Project summary

The current Lean library proves an exact-GECP rate for the original fermionic
kernel, not just low-rank existence. For every real `Λ≥1`, every actual finite
complete run on `[0,1]×[-Λ,Λ]`, and every physical point, with `k` successful
pivots,

`|R_k| ≤ 1024 exp(-((log 2)^2/256) k/log(1+Λ))`.

`fermionicKernel_finiteRun_error_le_exp` assumes neither a convergence envelope
nor future successful pivots. The original selected-core cross has the same
error, and the actual Green-integral difference is bounded by this error times
the density's L1 norm. Exact-zero stopping, real-parameter logarithmic separated
rank, arbitrary-domain strict-SPD GECP/Cholesky equivalence and original-kernel
`2*L*h` fill-distance control are also proved. An all-positive-cutoff
`log(2+Λ)` rate is separate; the literal all-positive `log(1+Λ)` version is
refuted. No arbitrary-η numerical-pivot rate or literature priority is claimed.

**The required original-SPEC scope is complete at its stated boundaries.**
Exact pushed-main verification passed at `e8c35ed`: 757 tests in 150.31 seconds,
3,634 Lean jobs, 876 allowed axiom reports, 22 public-import examples and
[CI 40051138962](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40051138962),
plus archive inspection and both installed-wheel environments. This is measured
closure evidence, not a status-label inference. The audit rereads all original sections and explicit lists,
checks actual proof hypotheses through a 22-example public consumer, and
recomputes numerical artifact checks instead of relying on inventory labels.
It repairs real conditional Lipschitz-certificate rounding and dyadic-parameter
coverage errors, including one-ULP/subnormal cases and standalone precision
metadata. A live dependency audit also requires a development-only filelock
update; numerical dependencies and historical data manifests are preserved.
Exact repaired-source and final documentation-head root, package/install and
push/PR CI pass; PR #123 and its exact-main closure are delivered. No new
independent multiagent review of the mathematics is
claimed: three attempted reviewers hit an account limit before contributing;
the primary completes this audit locally. The updated estimate is 100% of
required scope, not a theorem probability, a count-based verdict or resolution
of every optional research extension.

DP supplies the genuine
seven-cutoff residual-adaptive experiment identified as missing by the
[DO audit](FINAL_AUDIT.md). Both clean-source executions converge in all seven
cases and 28 tolerance prefixes through `1e-12` and repeat bytes; the normalized
PNG is inspected and byte-repeated. Independent MP192 checking passes all 21
selected-core and 28 final-cross checks. DP's 713-test root and clean package/
installed-wheel checks and final branch CI pass; PR #122 is merged and exact
pushed-main closure subsequently passes at `5a5fef0`, including CI `40018905634`.
DO's certificate/numerical repairs are merged in PR #121. The earlier 40/40
label count did not establish full completion; the DQ original-scope source
audit and final delivery evidence, not that count, govern the current decision.

## Practical application follow-on

The later practical implementation goal is tracked separately in
[PRODUCTIZATION.md](PRODUCTIZATION.md). DU–DX deliver the reusable physical
basis, complex batches/archives, matched DLR/IR comparison and small physical
solver, finite-query precision bridge, and a measured optional Rust backend.
DX is merged at `67e4b14`; its complete exact-main root/platform/install closure
is in [PR #130](https://github.com/tripp-smith/gecp-kernel-structure/pull/130#issuecomment-6400942719).
DY published [kernelgecp 1.1.0 on PyPI](https://pypi.org/project/kernelgecp/1.1.0/)
and [GitHub release v1.1.0](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0)
from immutable source/tag `9a679930db326bda2456bfa91b3256f37359e4ce` on
2026-10-25. The universal Python wheel, Linux x86_64/macOS arm64 native wheels
and source distribution are the unchanged exact-main tested artifacts.
[Publication workflow 40244470229, attempt 2](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40244470229/attempts/2)
is successful: all six Linux/macOS/Windows × Python 3.12/3.13 clean public-index
wheel installs, native execution where built, the Windows Python fallback and
a separate public-sdist install pass. Complete corresponding sources, notices,
dependency inventory and the manifest were public on GitHub before the OIDC
upload. The [publication receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/131#issuecomment-6404917491)
records the exact hashes and checks.

The upload occurred once, at 23:23:48–23:23:53 UTC. Five initial postchecks
exhausted their short readiness window while PyPI JSON still returned HTTP 404;
the upload, public-sdist install and Windows 3.13 verification had succeeded.
Once the public JSON and simple-index hashes matched the tested files, only
those five failed verification jobs were rerun. The first attempt remains in
Actions history; attempt 2 does not imply an initially all-green run. Final DZ
documentation and verification leave the released tag and artifact bytes
unchanged. The functional practical release/publication objective is complete.
[Actual upstream outreach](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
led to a [maintainer invitation](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6405100346)
and [Sparse-IR documentation PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10),
submitted as `bc4659e14c1a5f319fd6dcc7fa42796c7b29e2bd`. It adds three
executable cells after the conditioning example while preserving the original
eight cells and metadata. With the upstream lock, isolated and full-notebook
execution agree at rank 14: sampled time/frequency errors are
`2.7616964271004463e-11` / `4.94080361244419e-11` against `1e-8`. The full book
build and independent quadrature/content review pass. The
[contribution record](UPSTREAM_PROPOSAL.md) gives exact notebook/dependency
identities and preserves the corrected disposable-harness failure.
Submission-time upstream CI required permission and executed no jobs; the live
PR records subsequent review/check/merge state. Local checks and submission
do not establish maintainer acceptance. Local final delivery is recorded in
[PR #133](https://github.com/tripp-smith/gecp-kernel-structure/pull/133), separately
from the immutable 1.1.0 publication.
The original research's 100% scope estimate and this later product's measured
release gates retain their separate scopes.

## Historical implementation and research chronology

The following entries preserve what was available at each named phase. Their
“remaining”, “open” and obstruction-only statements are historical, not a
replacement for the current summary and final audit. v1.0.0 originally
delivered a `sorry`-free library and the `kernelgecp` research package, including
the dyadic `16p(s+1)` separated-rank construction. Its original E contract
closed with a sign-pattern obstruction and a conditional `CrossRatioControl`
rate; DG/DH later supplied the unconditional exact-GECP theorem above by the
local-minor/panel route, without proving that abandoned cross-ratio hypothesis.

The post-v1 Python application phase demonstrates how those components work
together on realistic synthetic data: one universal fermionic basis compresses
two qualitatively different continuous spectra, sparse routines represent
clean and noisy Green functions under explicit dictionary assumptions, and
the PSD baseline selects spatial covariance landmarks. These are application
observations with held-out errors, not new formal convergence or inverse-
problem claims.

Post-v1 Phase G adds a rigorous rate obstruction and the first strict
fermionic block theorem. Lean proves that the first two complete pivots are
the two reflected cutoff corners for every positive cutoff, that no universal
factor below one can contract the first step, and that for `0<Λ≤1` the actual
two-step residual is at most one half of the initial complete-pivot magnitude.
This verifies the base block but does not supply the dyadic scaling lemma
needed for Conjecture G1.

Post-v1 Phase BL closes the qualitative and pointwise parts of the amplified-
determinant strategy. Lean proves that nonnegative mixed coefficients make the
tail-weight pencil determinant-monotone, that every amplification weight at
least one preserves strict total positivity, and that the `8!`-amplified
positive-series kernel is bounded by two on the unit band. Complete runs of
that amplified kernel therefore have a `2^n` selected-core determinant bound.
The remaining rate gap is the arbitrary ordered-sample determinant envelope:
the original GECP nodes are not assumed to be complete pivots of the amplified
kernel.

Post-v1 Phase BM removes that pivot-path restriction. A generic northwest-
Schur induction proves that every strictly totally positive finite matrix has
determinant at most its diagonal product, so an entry bound applies to every
ordered sample, not only a complete run of the bounded kernel. Combining the
arbitrary-sample `2^n` amplified envelope with the Phase BK rank-eight
amplification gives
`det K_n <= 2^n / (8!)^(n-8)` for positive unit-band samples and `n >= 8`.
The remaining rate gap is now composition with the original GECP selected-
core/geometric-mean theorem and assembly across reflections and dyadic bands.

Post-v1 Phase BN completes the first exact GECP composition of that decay.
Independent sorting removes the mismatch between monotone determinant
theorems and algorithmic pivot order, after which the selected-core
geometric-mean inequality gives a residual-power bound. Exact arithmetic at
rank ten proves a uniform one-half contraction for complete pivoting on the
reflected positive-series kernel over the positive unit band. The remaining
global step is to assemble positive and negative frequency columns, close the
zero endpoints, and transport the unit-band block across dyadic scales.

Post-v1 Phase BO closes the full signed-unit-band and endpoint obligations in
one step. The affine coordinate `z = 1 - omega` rewrites the original kernel
as positive diagonal weights times `exp(t*z)` on `[0,1] x [0,2]`, avoiding a
mixed-column Laplace expansion. A weight-`128` eighth-order tail amplification,
strict total positivity, and the diagonal-product envelope prove
`|det K_n| <= 2^n / 128^(n-8)` for arbitrary injective samples. Composing with
the exact selected-core residual-power theorem proves a uniform half
contraction after twelve complete pivots on `[0,1] x [-1,1]`. The remaining
global step is dyadic rescaling/continuation across larger cutoffs.

Post-v1 Phase BP removes absolute band location from that local theorem. For
every real center, the coordinate `z = center + 1 - omega` maps the closed
width-two band to `[0,2]`. The weight-`128` amplified kernel remains strictly
totally positive and bounded by two uniformly in the center, so arbitrary
injective samples satisfy the same `2^n / 128^(n-8)` determinant envelope and
twelve exact complete pivots give the same one-half residual contraction. The
remaining global step is now specifically a conditioned-residual restart
across the selected outer-scale prefix, not an affine recentering problem.

Post-v1 Phase BQ identifies the exact restart interface. It disproves
one-sided frozen-column amplification by a kernel-checked rational
positive-Taylor witness, then proves that an exact low-rank-plus-tail split of
the actual final Schur residual retains the full fresh-order amplification
exponent. Phase AJ's bordered factorization composes this with the selected
core without charging prefix length against that exponent. What remains is no
longer an algebraic composition ambiguity: it is the analytic construction of
a compatible split, with nonnegative mixed coefficients, for arbitrary
fermionic final residuals.

Post-v1 Phase BR constructs that compatible split exactly. The pivot
coefficients of a fixed run define a linear companion row-conditioning
operator which commutes with a fixed right feature factor. Consequently the
positive-series fermionic decomposition remains a conditioned tail plus a
rank-eight block after every prefix, with no prefix-length penalty in the
Phase BQ exponent. Exact arithmetic also shows that conditioning does not
preserve coefficientwise nonnegativity: a mixed coefficient becomes
negative. The remaining analytic target is therefore a signed or real-rooted
conditioned determinant-pencil comparison, not the previously anticipated
all-nonnegative coefficient theorem.

Post-v1 Phase BS removes that failed sign theorem from the critical path.
Rank-eight vanishing leaves a polynomial of degree at most eight, so exact
Lagrange interpolation at nine weights between `56` and `128` controls its
value at one without assumptions on coefficients or roots. The formal
order-eighteen corollary converts nine `2^18` companion determinant envelopes
into the `2^-18` physical determinant bound required for a half-contraction
block. The remaining restart problem is thus a finite constant-base envelope
for those companion residuals, not a coefficientwise positivity theorem.

Post-v1 Phase BT gives that finite-envelope problem its natural two-sided
coordinates. Each fixed pivot applies a left and right oblique projector to a
companion, annihilating both selected directions while preserving linearity.
The projector transports `L*R` by conditioning the two factors separately, so
every fermionic prefix retains an exact two-sided rank-eight split. The
conditioned low determinant is now an ordered Cauchy--Binet sum of matched
left/right feature minors. Exact finite searches consistently recover a
common coefficient orientation, but the corresponding all-orders
matched-minor sign theorem remains open.

Post-v1 Phase BU removes the conditioning from that sign theorem. Exact left
and right feature-augmented determinant identities express each conditioned
feature minor as a Schur quotient by the selected core. In a matched product
the core is squared, so pivot order and core sign cannot affect the result.
For the fermionic rank-eight block, every conditioned Cauchy--Binet term now
has the sign of two original bordered feature determinants. Those borders are
confluent exponential systems: the polynomial features are derivatives of
the exponential family at parameter zero. Proving their oriented signs is the
next analytic step.

Post-v1 Phase BV corrects that proposed next step. An exact three-pivot
finite-Taylor witness has two negative products among its 56 order-three
bordered-feature coordinates, so coordinatewise confluent positivity is
false. The complete signed sum nevertheless has a canonical invariant form:
adjoin the selected kernel blocks `U,V` to the two feature factors `A,B`, mark
the fixed coordinates, and take their full-degree block coefficient. Lean
proves this is exactly the mixed-column coefficient of
`det(mu • (U*V) + A*B)`. In the witness the coefficient equals the signed sum
and remains positive. The open analytic problem is positivity of this
aggregate graded pairing for arbitrary fermionic prefixes, not positivity of
each Plücker coordinate. The root command passes a 3,571-job Lean build with
107 tests and no known dependency vulnerabilities; the change passed CI and
merged as PR #76.

Post-v1 Phase BW linearizes the whole Phase BV determinant pencil before any
coefficient expansion. The block matrix with corners `1`, `V`, `-mu • U`,
and `B` has determinant `det(mu • (U*V) + B)` without invertibility
assumptions, so every Phase BS finite weight is now one bordered determinant.
Exact rational checks cover weight one and all nine norming weights. A
negative four-cycle across the four blocks is unchanged by row/column sign
gauges, so this bordered matrix cannot be made entrywise nonnegative and a
naive total-positivity proof is impossible. The remaining route must control
the global determinant through cancellation-preserving structure. The root
command passes a 3,572-job Lean build with 107 tests and no known dependency
vulnerabilities; the change passed CI and merged as PR #77.

Post-v1 Phase BX redirects the bordered method to the actual Phase BS
companion premise. A fixed two-sided projection has the closed four-term Schur
formula determined by the selected core, the original-kernel cross blocks,
and the companion cross blocks. Lean packages it as one double-bordered
matrix and proves its determinant is the squared selected-core determinant
times the projected companion determinant. Exact arithmetic confirms that
three iterative pivots agree entrywise with this formula at all ten relevant
weights. The large matrix still has negative order-two minors, so the next
step is a normalized double-border bound, not total nonnegativity. The root
command passes a 3,573-job Lean build with 107 tests and no known dependency
vulnerabilities; the change passed CI and merged as PR #78.

Post-v1 Phase BY proves that this double border is not merely a finite witness
construction. For every successful run, the iterated row and column
projectors commute and equal their original selected-core corrections, so
`Run.conditionBoth` is exactly the four-term companion Schur expression.
Every finite sample consequently satisfies the Phase BX determinant identity
with the actual recursive companion determinant. Exact rational solves check
144 deterministic randomized combinations of prefix length, companion weight,
and fresh block. The remaining task is now the normalized double-border
envelope itself, with no recursive-versus-closed-form gap. The root command
passes a 3,574-job Lean build with 107 tests and no known dependency
vulnerabilities; the change passed CI and merged as PR #79.

Post-v1 Phase BZ places the double border inside the upper-tangent kernel
`[[K,H],[0,K]]` and factors the latter through the common feature matrices and
the middle block `[[1,W],[0,1]]`. The positive-series specialization proves
that normalized bounds for its nine distinguished minors imply the Phase BS
half-power determinant estimate. Exact rational checks find negative ordinary
order-two minors but a consistent cross-inversion orientation for all 360
tested complete-prefix coefficient polynomials. This isolates a narrower
flag-minor or two-layer planar-network theorem as the next target. The root
command passes a 3,577-job Lean build with 108 tests and no known dependency
vulnerabilities; the change passed CI and merged as PR #80.

Post-v1 Phase CA expands every finite sampled upper-tangent factorization over
two ordered feature flags. Each summand is now the product of an explicit
left block-diagonal feature minor, a middle unitriangular tangent minor, and a
right block-diagonal feature minor. Lean proves that independent orientation
certificates for those three factors compose to the required determinant
orientation, without assuming the conclusion. Exact rational arithmetic
checks 24 shuffled rectangular triple expansions and retains the negative
ordinary-minor guard. This localizes the remaining theorem to a middle
shuffle classification and two confluent Vandermonde flag-minor signs. The
root command passes a 3,579-job Lean build and 108 tests with no known
dependency vulnerabilities; the change passed CI and merged as PR #81.

Post-v1 Phase CB resolves the source of the middle-factor shuffle sign. After
interleaving the two copies of each feature, a diagonal tangent middle matrix
is exactly a block diagonal family of `[[1,w_k],[0,1]]` cells. Lean proves
every local cell totally nonnegative, computes the complete determinant as
one, and proves that total nonnegativity of the interleaved direct sum yields
the exact tuple-orientation certificate needed by Phase CA. An exhaustive
exact sweep checks 198,626 feature-major minors and their 198,626 oriented
layer-major counterparts through five cells, while retaining a negative
unoriented witness. The remaining middle obligation is the generic
all-orders TN closure for this sparse upper-bidiagonal direct sum. The root
command passes a 3,580-job Lean build and 109 tests with no known dependency
vulnerabilities; the change passed CI and merged as PR #82.

Post-v1 Phase CC closes that middle-factor obligation. A reusable noncrossing
support criterion proves total nonnegativity by eliminating every
nonidentity term in the Leibniz expansion. The feature-major tangent middle
matrix satisfies the criterion for arbitrary nonnegative weights, including
zero weights, and therefore supplies the exact layer-major shuffle-oriented
minor sign required by Phase CA without a remaining hypothesis. Exact
arithmetic checks 225 crossing pairs and 12,403 nonidentity matching terms in
addition to the Phase CB minor census. The two outer flag-minor families are
now the only missing sign inputs to the three-factor expansion. The root
command passes a 3,581-job Lean build and 110 tests with no known dependency
vulnerabilities; the change passed CI and merged as PR #83.

Post-v1 Phase CD shows that the remaining obstacle cannot be split into
independent outer-factor signs. In an exact two-node, three-feature positive
monomial tangent factorization, one ordered double-flag term is `-1/192` and
Lean proves that no separate left/middle/right orientation functions can make
all terms nonnegative with external signature one. The full distinguished
determinant is nevertheless `24929/27648 > 0`. Exact enumeration confirms 10
positive and 2 negative nonzero terms, while the two-feature predecessor has
only two positive terms. The proof route must therefore preserve cancellation
across flags, most plausibly through a global confluent/jet determinant,
planar-network, or exterior-algebra argument.

The root command passes a 3,582-job Lean build and 111 tests with no known
dependency vulnerabilities; the change passed CI and merged as PR #84.

Post-v1 Phase H restores verification integrity after the local editable
pytest entry point became stale and the development lock acquired dependency
audit findings. The root and CI workflows now launch pytest through the
synchronized interpreter, the audited transitive lock resolves `urllib3`
2.8.0, and the complete root command passes in maintained and fresh frozen
environments. This is maintenance evidence and changes no mathematical claim.

Post-v1 Phase I isolates the next exact analytic obligation. Lean proves a
bordered selected-core determinant identity and uses core nonsingularity to
show that selected-residual cross-product sign coherence is equivalent to a
four-bordered-minor sign condition on the original kernel. Exact surrogate
checks retain first-failure witnesses. This is an algebraic reduction, not a
proof that the fermionic kernel satisfies the determinant-sign condition.

Post-v1 Phase J closes the algebraic part of that obligation. Lean formalizes
ordered strict sign regularity, transports it through arbitrary tuple
orientations and positive column scaling, and proves that an all-orders
exponential sign theorem would force selected-residual sign coherence for
every fermionic GECP run. Orders one and two are proved unconditionally. The
all-orders exponential theorem itself remains the precise analytic gap.

Post-v1 Phase K begins the constructive all-orders proof. Lean proves positive
row scaling, strict positivity of ordered ordinary Vandermonde determinants,
and strict positivity of the principal square exponential Taylor feature
block at every order. This certifies the principal Cauchy--Binet term while
leaving generalized Vandermonde positivity and the infinite-series limit open.

Post-v1 Phase L proves the positive-node generalized Vandermonde theorem at
every finite order. Lean derives a sparse-polynomial positive-root bound from
Descartes' rule, uses it to prove power-matrix nonsingularity, anchors the sign
at geometric nodes, and transports positivity through a zero-avoiding
homotopy. The zero-node boundary case and exponential-series determinant limit
remain open.

Post-v1 Phase M closes the zero-node boundary: nonnegative strictly increasing
nodes and strictly increasing natural exponents give a nonnegative generalized
Vandermonde determinant. This supplies the missing termwise sign input for a
finite exponential Cauchy--Binet expansion.

Post-v1 Phase N assembles that finite argument. A reusable rectangular
Cauchy--Binet identity groups injective intermediate tuples into increasing
tuples and permutations. Applied to finite exponential Taylor features, every
minor product is nonnegative and the principal product is strictly positive,
so every truncation with at least as many terms as the matrix order has
positive determinant. Passing to the infinite-series determinant while
retaining strictness remains open.

Post-v1 Phase O passes to that limit without losing strictness. Each finite
Taylor determinant dominates the same positive principal determinant; the
feature matrices converge entrywise to `exp(row_i * column_j)`, and
determinant continuity preserves the lower bound. Thus the full exponential
matrix has positive determinant on nonnegative strictly increasing node
tuples. Arbitrary-node scaling and the frequency-reversal orientation remain
before the all-orders sign-regularity theorem.

Post-v1 Phase P completes the all-orders sign theorem. Positive exponential
row and column scalings remove arbitrary node offsets, finite reversal has
sign `(-1)^(n.choose 2)`, and negated reversed frequencies convert positive
`exp(x*y)` determinants into the prescribed `exp(-t*omega)` signature. The
fermionic kernel inherits strict sign regularity at every order, so every
selected residual cross is sign coherent and exact complete-pivot updates are
nonexpansive. This still does not supply strict dyadic block contraction.

Post-v1 Phase CX identifies the remaining fresh-order-two geometry more
sharply. An exact complete-pivot exponential witness has two independent
positive directions in its normalized matched-coordinate Hessian, so a global
Lorentzian or Hodge--Riemann signature theorem is false. The three contiguous
degree-prefix planes in the same witness still satisfy the sharp reverse
Cauchy--Schwarz inequality. The surviving target is therefore a prefix-cone
variation-diminishing or cumulative-interlacing theorem, not global stability.

Post-v1 Phase CY supplies the adjacent-cutoff algebra. A shared value at the
canonical negative node factors into a product of two projected cross
determinants. Aggregate sign agreement, positive transfer sums, and
nonnegative endpoints then imply both neighboring matrix discriminants. The
finite-Taylor corpus verifies 13,924 identities and 9,170 strict interior
signs. The analytic sign for arbitrary conditioned blocks remains open.

Post-v1 Phase CZ refutes that universal interior fixed-node sign for finite
Taylor kernels. A four-node, nine-term kernel has two certified complete
pivots, a positive degree-three canonical numerator, and positive adjacent
discriminants. The neighboring degree-two and degree-four nodes provide
alternate certificates. The exact 100-case audit retains four failed fixed
nodes while covering all 47,496 interior cutoffs. A proof of adaptive coverage
or a different spectral representation is now the next target.
The sampled actual-exponential calculation agrees at 80 and 160 decimal
digits, but is numerical evidence rather than an infinite-tail proof. The
formal pivot claims are on the four-node sample, not the whole continuous
domain.

Post-v1 Phase DA proves an analytic base case. For `0<a≤1` and
`0≤x<y<z≤1`, the actual-exponential determinant
`det[exp(a*t),t^i,t^j]` has a node-independent sign for every natural `i<j`:
negative for `(0,j)` with `j≥2`, positive otherwise. Independent ordered
triples, even at different positive unit-band rates, therefore have a positive
matched product. The proof uses row normalization, explicit derivative
ratios, and Cauchy's mean value theorem; zero and one boundary nodes are
included. This does not extend the sign rule past one pivot. A separate exact
depth-ten conditioned-feature regression rejects the proposed `depth+1`
localization of failed transfer nodes. Arbitrary-prefix cumulative coverage,
the higher-fresh-order envelope, and the final dyadic rate remain open.

Post-v1 Phase DB rules out neighboring canonical-node coverage as a general
finite-prefix invariant. A coherent thirteen-term Taylor kernel on nine
open-unit rows and columns has seven certified complete pivots, then two
consecutive strictly positive transfer-node numerators at cutoff eight.
All transfer sums and neighboring endpoint/mixed coefficients are positive;
the cutoff discriminant is positive. The integer certificate verifies the
actual residual recurrence, both feature projectors, Taylor reconstruction,
and sign transport back to the actual conditioned blocks. The failure does
not refute the Lorentzian inequality or a continuous-domain complete-pivot
theorem. Actual-exponential 160/320-decimal calculations agree on the two
failed signs to 45 significant digits, but are numerical evidence only.

Post-v1 Phase DC proves an explicitly weaker but all-orders local-minor
contraction bridge. A column-maximal scalar Schur step in a strictly totally
positive ordered sample is bounded by an old minor on the same fresh columns,
with a different row tuple chosen by the pivot column's sorted rank. The
nonsymmetric singleton Fischer inequality, arbitrary-pivot absolute determinant
factorization, generic kernel bridge, and actual positive-series fermionic
specialization are kernel-checked. The original arbitrary-prefix companion
target is still open. Exact prefix/gauge tests and 320/640-digit depth-64
companion stress evidence are not promoted to formal closure. No new public
Python API or release is introduced.

Post-v1 Phase DD proves a common magnitude-one row/column gauge for every
minor order of an actual residual after an arbitrary successful positive-minor
kernel prefix. It preserves absolute minors, and a Schur-scaling argument
applies Phase DC's matched-rank bound to the actual next residual. A verified
frequency-reflection transport covers original fermionic runs, signed
frequencies, and boundary times without reselecting pivots. The initial local
quantitative budget and its iteration/pivot count remain open; this does not
solve G1 or Problem 4.2.

Post-v1 Phase DE proves the initial scale-free local panel budget. A signed
central/dyadic panel needs `8p` separated terms for error `2^(-p)`, independent
of its scale. Generic rank-sensitive determinant decay then gives
`abs(det K_n) <= (2^(-q))^n` at `n=16*(33+4q)` on every panel. The theorem
concerns the actual kernel on closed physical time, with unordered/repeated
samples permitted. Backwards minor-witness iteration and global panel-pivot
counting remain open; this is not yet G1 or Problem 4.2.

Post-v1 Phase DF proves backwards marked-pivot witnesses and the global finite
rate. On `[0,1] x [-2^s,2^s]`, the actual fermionic residual is bounded by
`2^(-q)` after `32*(s+1)*(33+4*q)` exact complete pivots and at every later
rank. All structural and local-budget premises are discharged by DD/DE;
remaining assumptions describe initialization and a successful exact realized
complete-pivot sequence. Independent finite exact/high-precision checks
supplement, but do not prove, this universal result. The real-parameter
corollaries and full completion audit remain open.

Post-v1 Phase DG proves the actual real-cutoff logarithmic rate. For
`Lambda>=1`, the bound is
`abs(R_k) <= 1024*exp(-((log 2)^2/256)*k/log(1+Lambda))`.
For all positive cutoffs the same constants work with `log(2+Lambda)`;
the real-tolerance theorem gives an explicit ceiling rank. Initialization
and exact updates supply realized prefixes without a duplicate premise.
The current public interface describes an exact successful complete-pivot
sequence on all natural indices; it is not an approximate-pivot certificate.

DG also proves a parameter obstruction: literal G1 with `log(1+Lambda)`
cannot hold uniformly down to zero. For any `C,c>0`, the actual first complete
pivot at `Lambda=min(1,c^2/(16*C))` violates that bound. The SPEC conjecture is
preserved, with this boundary recorded as an annotation. The large-cutoff
theorem remains valid. Full completion is still subject to the separate
finite-prefix/termination, cross-approximation, and specification audit.

Post-v1 Phase DH closes the finite-run and physical-error interfaces. The
same explicit logarithmic rate now takes an actual finite complete
`Run fermionicKernel`, with no residual-sequence, future-success, next-pivot,
or infinite-continuation premise. Finite backwards witnesses and the last
successful update supply the bound. Exact stopping is a zero-valued complete
maximum, not a zero-pivot update or a floating-point heuristic.

The selected-core `Run.crossApproximation` is proved to have error exactly
equal to the iterative final residual. The Green consequence bounds the actual
difference of weighted kernel integrals and derives both integrability facts.
Its interval form uses the measure restricted to `[-Lambda,Lambda]`, so no
separate support premise is needed. Numerical checks use solves and compare
exact finite algebra with independent actual-exponential precision. Approximate
pivot robustness, runtime continuous-domain certification and full-SPEC
completion are not claimed by this phase.

Post-v1 Phase DI audits the original specification rather than relying on
historical bounded-phase labels. [COMPLETION_AUDIT.md](COMPLETION_AUDIT.md)
records 33 supported and seven incomplete requirement groups (about 80%
checklist coverage, not weighted effort). The full goal is not complete.
The gaps are continuous SPD/original-kernel Lipschitz baseline, result metadata
outside GECP, census tolerance/diagnostic/normalized-curve scope, and listed
Phase-II sparse comparisons.

DI supplies the real-parameter D1 corollary, compact complete-maximum existence,
exact small-grid certificates, fail-closed axiom-report enforcement and clean
packaging checks. A temporary local source archive containing copies of two
untracked notes was removed; rebuilt archives exclude them and their originals
are unchanged. An isolated wheel smoke passes; nothing was published.

Post-v1 Phase DJ closes R09 with original-kernel, not residual, regularity:
`diagonalRun_error_le_original_fillDistance` proves `abs R(x,y) <= 2*L*h`
for actual finite diagonal runs on arbitrary metric domains. Residual and
removed-kernel PSD, original-energy comparison, arbitrary-domain continuity,
and compact diagonal complete-maximizer availability are derived. An actual
infinite-domain linear-feature run and exact rational repeated-sample/prefix
checks exercise these results. Strict-SPD necessity for every successful
complete pivot remains open in R08. The updated audit has 34/40 supported
groups (85% checklist coverage, not weighted effort) and six incomplete groups.
Verification passes 3,630 Lean jobs, 857 enforced axiom reports, 180 tests and
no known dependency vulnerabilities; push and PR CI pass, and PR #116 is
merged. No new numerical behavior or release is claimed.

Post-v1 Phase DK closes R08's staged continuous baseline. From the original
kernel's strict finite-support `Matrix.PosDef` condition on any ambient type,
every actual finite exact complete run is derived to be diagonal. Strict
positivity on fresh samples is preserved by an injective rectangular Schur
congruence. `posDef_completeRun_iff_pivotedCholesky` equates actual complete
runs with independently defined recursive diagonal maximization, without a
diagonal-preference premise. Residual PSD, diagonal domination, compact
continuous maximum existence and original-kernel fill control follow; the
latter no longer requires supplied diagonality. Exact branch/prefix and Lean
Gram-factor checks accompany the theorem. Weak-PSD off-diagonal ties remain
excluded from necessity; no unique-tie, approximate-pivot or full asymptotic
P-greedy theorem is claimed. The audit is now 35/40 supported groups (87.5%
checklist coverage, not weighted effort) with five numerical/comparison gaps.
Full verification passes 3,633 Lean jobs, 876 enforced axiom reports,
183 tests and no known dependency vulnerabilities; push and PR CI pass;
merged as PR #117.

Post-v1 Phase DL completes R24 with typed numerical evidence metadata.
Cholesky/sparse results expose actual precision, conditioning histories,
tie/settings and certification/error scopes; GECP exposes a view of existing
evidence. Fixed constructions explicitly report non-iterative scopes.
Schema-2 synthetic summaries preserve metadata, while committed historical
schema-1 JSON/PNG artifacts and hashes remain unchanged. Returned-error
convergence, candidate exhaustion, rank-zero/non-applicable states, optimizer
budget/failure/rejection and diagnostic failures are independently tested.
No new Lean theorem, high-precision Cholesky/sparse execution or continuous
floating-point certificate is claimed. The audit has 36 supported and four
incomplete groups (90% checklist coverage, not weighted effort). Full gates
pass 3,633 Lean jobs, 876 enforced axiom reports, 232 tests and a clean
dependency audit after targeted certifi/pillow lock updates. PNG and isolated
wheel smoke checks and push/PR CI pass; merged as PR #118.

Post-v1 Phase DM is verified and merged in PR #119. The full schema-2 census
contains 28 converged sampled-grid cases through `1e-12`, with complete
residual/core and geometry histories. Two full runs from clean generator
`c0e8c5093134f369b93781c413750b7f1de46036` match bytes, SHA-256
`c6854330fcb74007f8f9814fba48c314d34fe26c1cef5461d908ea5182b66bda`;
observed concurrent wall times are 1496.40 and 1508.07 seconds. Strictest ranks
are `9,17,37,63,91,119,146` across the seven cutoffs. Both new plots have been
rendered, inspected and byte-repeated. The legacy artifacts and all 21 shared
pivot/rank traces are preserved. Preliminary 192-bit primary-core checks and
panel-spanning held-out samples pass their stated numerical checks. Actual
order-32 GECP at cutoffs `1,1000,1000000` also converges with ranks `9,63,146`,
matching the corresponding order-24 ranks. The denser dataset is generated
from clean `ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61`, SHA-256
`7251cbb528179cc5cbfa6eb5635b8fa280e087d1e04db1348b56a9065026d588`,
in 2004.01 wall seconds. This is a representative two-resolution study, not
an all-cutoff resolution sweep or grid-independent-pivot claim. The combined
report passes all 30 selected-core comparisons at 192 bits, SHA-256
`b7192e8bc57a72bf64802dcb9334d4c1a4761be5a2a2876368fd72aeda40247d`.
The full local root gate passes 3,633 Lean jobs, 876 axiom reports, 430 tests
in 139.64 seconds at the final branch head, plus the clean dependency audit.
Source/wheel inspection and an isolated installed-wheel smoke pass. Exact-head
push/PR CI pass at `4f99db923e89773e1d248eecde872f2350e0d853`; PR #119 merged
at 16:44:20 UTC on 2026-10-22 as
`b068e1070745d1308ee23101096e1197c345f057`. The pushed main closure head
repeats root/CI gates before final handoff. The audit now supports 39/40 groups
(97.5% checklist coverage, not weighted effort). None of these observations
is a continuous certificate or new Lean rate theorem; R31 comparisons and the
final full-SPEC audit remain required.

## DN six-method comparison (verified, merged PR #120)

DN implements all six SPEC §20 methods on identical time observations,
with original-kernel least-squares weights and separately scoped native,
training and held-out convergence. Genuine VF eigenvalue pole relocation,
SciPy AAA and data-ranked strided-Hankel Prony are not renamed nonlinear LS.
All 192 fits are retained, including nonphysical pole rejections and failures.
The known-library control is separate from the off-grid primary cohort.

The canonical JSON comes from clean
`072ce1004de11422f05ef1fbe3e738093d6f2f13`, configuration SHA-256
`665e158dd17b5722ecd1e34b2d549cbe0dcd413613dd7f5d9c33a3af54cdeeb7`;
two complete runs match bytes, JSON SHA-256
`9f85509c929bc0511b0f34f25af5188420902a165ddd3697631607d16fd5ec31`.
The reviewed PNG repeats bytes, SHA-256
`0368a73f27aa4b15ed9c313cb73a7326157d706b6aea22ca5c798cd75fe834d2`.
The final schema preserves the universal generator's MP decimal histories as
well as labeled float summaries. Independent tests reconstruct every reported
time error and smallest-tested support.

At cutoff eight and `1e-8`, universal kernel rank is 12. AAA/VF/Prony/NLS use
four tested atoms for the clean off-grid fixture. VF/Prony/NLS also pass with
four on the noisy fixture and six on both fixed-quadrature density references.
AAA fails the other three physical fits. No support optimum, uniqueness,
continuous density-integral error or universal method ranking is claimed.
The full root at `256779028bd5a07736730bf3ff70587f80b91c75` passes 513 tests
in 141.91 seconds, 3,633 Lean jobs, all 876 axiom reports and the clean
format/type/dependency gates. Exact-head push CI `39831341046` and PR CI
`39831351890` pass. Source/wheel inspections and installed-wheel smokes outside
the checkout pass with both SciPy 1.18 and the declared minimum 1.15. Nothing
is published. Final branch `eec016d0b4040941e9e41f825bc0a3ac5f62f8bf` repeats
the full root gate (513 tests in 140.52 seconds); push CI `39832373652` and PR
CI `39832383406` pass. PR #120 merged on 2026-10-22 at 18:37:17 UTC as
`a3186a89b27d92677ef0e3338b90f9048a28fdb5`. Exact pushed-main closure gates
repeat before handoff; the final independent whole-SPEC audit remains next.

## DO original-SPEC audit and numerical correctness repairs

The [independent final audit](FINAL_AUDIT.md) maps all thirty original sections
to source/artifact evidence rather than assuming the prior forty labels imply
completion. Its formal assumption audit supports the finite exact rate,
strict-SPD baseline, real separated witness and physical cross/Green claims.
The required all-cutoff residual-adaptive experiment is still absent; fixed-grid
diagnostics and normalized curves remain useful evidence at their actual scope.

DO corrects a counterexample where rounded point cancellation produced a
claimed lower bound above the interval upper bound. Candidate lower magnitudes
now come from singleton interval enclosures; the interval callbacks must bound
the same target function. Adaptive GECP certificates concern the residual model
defined by stored rounded coefficients, not an ideal exact run or a new
float64 core solve. Observed candidate ties do not count all global maximizers.
A representable large-finite rank-one matrix receives a stable update/failure
path, and the centered formula has an independent regression. Historical DN
source hashes remain tied to clean generation revision `072ce10`, not relabeled
after source repairs. Generic completeness definitions move below the
fermionic layer, preserving all theorem interfaces and axiom obligations.

The final implementation `f0d082dfbbc91e4213dea0d2cad06a36c4a78826` passes
569 tests in 142.82 seconds, 3,634 Lean jobs, all 876 axiom reports, 83-file
formatting, mypy on 19 sources and the clean dependency audit. Push CI
`39961732232` and PR CI `39961740682` pass all checks. Source/wheel archives
are clean; all nineteen wheel Python files hash-match that implementation.
Installed-wheel checks outside the checkout pass both pinned SciPy 1.18/
NumPy 2.5.2 and minimum SciPy 1.15.0/NumPy 2.4.6, including all 45 exports,
the repaired numerical paths and 48 six-method fits each. The companion block
diagnostic preserves binary coordinates and the same frozen-model scope;
float64 and MP ties exclude zero/old coordinates. These remain numerical
models, not exact-run certificates. PR #121 merged at 2026-10-23 13:34:53 UTC as
`b041c723add1b9a7a0395cfffbdf4b793c0cfed9`. The final pushed-main metadata
head repeats root/CI and archive checks before handoff. R25 remains incomplete:
39/40 supported groups is 97.5% checklist coverage, with a rough 97% overall
effort estimate, not full completion or a theorem-success probability.

The workshop's [v3 problem-list update](https://arxiv.org/html/2602.05394v3)
attributes an August 2026 claimed solution of Problem 4.2 to Venkata Siddharth
Pendyala's *Local-to-Global Convergence of Greedy Cross Approximation for Totally
Positive Kernels* ([DOI](https://doi.org/10.5281/zenodo.21863274)). This audit
could inspect that update but not the underlying paper. It records the prior
claim without independently verifying it or asserting priority here. The
[Cholesky v1](https://arxiv.org/abs/2509.13582v1) citations are historical;
[v3](https://arxiv.org/abs/2509.13582v3) adds Christopher Wang to Jeong and
Townsend. The staged fill theorem does not claim all published asymptotic rates.

## DP residual-adaptive census

DP delivers `kernelgecp/adaptive_grids.py`, `adaptive_grid.py`, the schema-3
explicit-axis runner, independent checker and stage-aware normalized plot.
Residual/surplus scores drive selective mirrored panel refinement; every old
coordinate is retained. Two consecutive jointly passing strict stages are
required. Requested prefixes advance beyond the first training-passing rank
when probes demand it. Six seventh-fraction samples per final leaf are used
only after all refinement and rank choices; all 28 final checks pass.

The two fresh clean-source executions at generation revision
`3cee4186fddf990abdd33dbe24c46ff908b1923f` yield identical
`experiments/data/adaptive_census.json`, SHA-256
`730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c`.
At cutoffs `1,10,100,1000,10000,100000,1000000`, actual stage counts are
`3,2,3,4,4,5,6`, final axes are respectively `25×25,41×57,109×89,122×121,
170×169,238×217,258×257`, and `1e-12` jointly chosen ranks are
`9,17,36,64,97,124,151`. Strict untouched sampled maxima range from
`9.792e-14` to `6.418e-13`. Complete per-tolerance histories, geometry, initial
and final axes and decisions remain machine-readable. These ranks are not
global minima, and nonmonotone stage behavior is preserved.

The inspected/repeated `adaptive_census.png` has SHA-256
`dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6`.
It plots all seven normalized trajectories, untouched validation markers,
training/joint ranks and stage growth. One zero is explicitly omitted from log
axes without a positive display floor. The independent MP192 post-run report
passes all 21 selected-core and 28 final-cross checks and binds the input and
historical source blobs. Every recomputed final sampled maximum meets its
tolerance. Report SHA-256 is
`4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03`.
The full root passes 713 tests in 151.22 seconds, 3,634 Lean jobs, 876 allowed
axiom reports, 95 formatted files, mypy on 21 sources and the clean dependency
audit. Source/wheel archives are inspected; all 21 wheel Python files and
package README match. Source includes the canonical artifacts/tests and excludes
owner note copies, caches and scratch. Outside-checkout installed-wheel smokes
pass with Python 3.12.12/mpmath 1.4.1 on pinned NumPy 2.5.2/SciPy 1.18.0 and
minimum NumPy 2.4.6/SciPy 1.15.0. Both load all 45 exports, exercise adaptive
success and budget failure, and complete 48 six-method fits. The exact
artifact head `a6488cc0c91400bab3f327d93db153c390ef20f4` repeats the root with
713 tests in 150.40 seconds and passes push CI `40017052450` and PR CI
`40017056920`, all four checks. Both rebuilt installed environments pass again.
The final documentation head `3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes
the root again (713 tests, 150.98 seconds), push CI `40017937659` and PR CI
`40017945103`; rebuilt archives and both installed smokes pass. PR #122 merges
at 2026-10-23 20:20:48 UTC as `3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`.
The generation commit remains an ancestor of main. DP is verified and merged;
exact pushed-main closure repeats the gates before handoff. Archives are rebuilt after final
documentation synchronization rather than assuming an earlier package matches.

Generation uses frozen nine-stage / `0.01` margin controls, not the public
five-stage / `0.25` defaults. Actual environment is Python 3.12.12, NumPy 2.5.2,
SciPy 1.18.0, mpmath 1.4.1 (Python backend), with recorded config/lock/source
hashes. After both immutable runs, the release runner's post-work untracked
experiment-code check is repaired and dirty noncanonical provenance propagates
to every record. Actual generation had no such source change. Preserve its
real `3cee4186` manifest rather than relabeling it with this later hygiene fix.
No numerical core or frozen controls change. All errors are finite sampled
evidence, not continuous, interval, exact-run or separately certified-cross
bounds. Full original-SPEC closure follows this phase.

## DQ original-scope audit and boundary repairs

The explicit-item tables in `FINAL_AUDIT.md` cover all thirty original
sections, eleven §14 diagnostics, eleven §25 numerical test categories,
six §20 methods, seven proposed imports, six milestones and six minimum
outcomes. Actual original-kernel finite-run proofs, arbitrary-domain strict-SPD
types and the literal small-cutoff obstruction are checked through source
dependencies and the new public-import consumer. Optional research routes
and approximate-pivot extensions remain explicitly unclaimed; illustrative
API syntax is distinguished from delivered residual/matrix semantics.

The Lipschitz branch-and-bound repair uses exact binary-rational radii from
the actual sample, outward upper conversion, exact gap checks and downward
eta. Exact fixed-function callback values and valid whole-cell analytic
constants remain required. Unresolved coordinate resolution is uncertified;
out-of-range upper bounds fail explicitly. Standalone certificates expose
point/enclosure/export precision and tie/stopping/certification metadata.
The caller's exact interval bit precision is restored on success and failure.
The Taylor factory/rank helpers now select dyadic parameters by exact integer
ratios, not rounded logarithms, covering adjacent powers and subnormal
tolerances; the DLR count conservatively bounds its central term. None of
these numerical fixes is a new exact-run theorem or rounding certificate for
the Taylor evaluator. Regressions reproduce the old defects before repair.

The full independent adaptive MP192 report and both adaptive/spectral PNGs
are freshly regenerated and compare byte-for-byte. The retained second DP
canonical generation also still equals the committed JSON. The initial
753-test DQ root correctly fails newly reported filelock advisories; updating
only the development cache-lock dependency/constraint to 4.2.1 restores the
audit while preserving NumPy 2.5.2, SciPy 1.18.0 and mpmath 1.4.1. Historical
generation-lock hashes remain checked against their actual generation Git
blobs, not relabeled as the current lock. The repaired preliminary root then
passes 753 tests in 151.32 seconds, 3,634 Lean jobs, all 876 allowed axiom
reports, the 22-example consumer and style/type/security gates. Later dyadic
boundary changes require their own exact-head verification; those preliminary
timings are not final-head evidence.

Exact DQ code head `f2ffd5e344818e2c06ff5152a7f182ef7a31e3a8` subsequently
passes the full root: **757 tests in 150.01 seconds, 3,634 Lean jobs, 876
permitted axiom reports, 22 public-import examples, 97 formatted files, mypy
on 21 sources and clean dependency audit**. The preceding development-tree
root took 148.83 seconds. Push CI `40050118292` and PR CI `40050129457` pass
all four checks at that exact code commit. Clean rebuilt archives hash-match
21 wheel Python sources, all 153 owned Lean sources, claim docs and lock;
they include new consumers/tests and exclude local notes/caches/scratch.
Reinstalled wheels pass outside the checkout in pinned Python 3.12.12/NumPy
2.5.2/SciPy 1.18.0/mpmath 1.4.1 and NumPy 2.4.6/minimum SciPy 1.15.0. Each
checks all 45 exports, adaptive success/budget failure, 48 six-method fits,
the new boundary/metadata paths and exact interval-context restoration.
Draft [PR #123](https://github.com/tripp-smith/gecp-kernel-structure/pull/123)
is attached to the task. Code-head success is not relabeled as a later run.

Final documentation branch `0a6d3635529279beac9fb198c2d9de147c5ee2e5` repeats
the full root with 757 tests in 149.97 seconds and all established Lean,
consumer, axiom, style, type and security gates. Push CI `40050605960` and PR
CI `40050610230` pass all four checks at that exact head. Clean archives are
rebuilt and both reinstalled outside-checkout environments pass again.
PR #123 merges at 2026-10-24 01:41:41 UTC as
`d9978607b88a45d98d917f44318bcaf99b4c3239`; local main fast-forwards and clean
generation `3cee4186` remains an ancestor. Closure main
`e8c35eda7662566639c5eac6748d3032df87075f` passed the exact root with 757 tests
in 150.31 seconds, all 3,634 Lean jobs, 876 permitted axiom reports, 22 consumer
examples and style/type/security gates. Exact-head CI `40051138962` completed
successfully on 2026-10-23 at 21:52:41 EDT (2026-10-24 01:52:41 UTC), with
both Lean and Python successful. Rebuilt source/wheel archives hash-matched
21 Python sources, 153 Lean sources, claim docs and lock; both reinstalled
outside-checkout environments passed. Final Git checks confirmed local/main
and origin/main equal, no tracked changes and no open PRs. The original goal
was marked complete after those actual gates, not from the inventory count.

## DR verified documentation finalization

The user-requested final documentation/hygiene phase updates all nine project
claim documents, records the completed DQ main result and terminal goal
accounting, labels historical snapshots, aligns closed-phase states and makes
exploratory command examples write ignored scratch outputs. Two documentation
regressions check the terminal record and local Markdown link targets. A
read-only reviewer found no factual/scope/provenance defect; the requested
phase-table alignment is included. AGENTS and the preserved autonomous
workflow reference were inspected and left unchanged.

Exact source `51fa6b3bdc28b02b7d2a45cbabaab71fce56e9bd` passed the full root:
759 tests in 144.90 seconds, 3,634 Lean jobs, 876 permitted axiom reports,
22 consumer examples, 97 formatted files, mypy on 21 sources and clean
dependency audit. Push/PR CI `40054146369`/`40054160333` passed all four checks.
Clean source/wheel inspection matched 21 Python/153 Lean sources, claim docs,
lock, documentation tests, AGENTS and CI. Both outside-checkout installed
environments passed all 45 exports, adaptive success/budget failure, 48 six-
method fits and boundary/metadata checks. Original SPEC prefix bytes and all
historical experiment data/manifests are unchanged.

[PR #124](https://github.com/tripp-smith/gecp-kernel-structure/pull/124) merged
on 2026-10-23 at 22:30:40 EDT as
`a421f7c8ed8b1b8754f38116939e26181636fe82`. Closing metadata is delivered on
main. Local version-1.0.0 archives are build snapshots, not a new publication;
old versioned build artifacts/caches remain ignored and protected user notes
remain untouched and untracked. No new goal or phase-only accounting is
invented; the completed original research scope is not reopened.

## Implementation run metadata

The following table is the historical v1 implementation record, captured at
handoff on 2026-08-15 at 17:39 EDT (UTC-04:00).

| Field | Recorded value |
| --- | --- |
| Codex surface | Codex desktop app |
| Declared model identity | OpenAI Codex based on the GPT-5 family |
| Exact serving model/version | Not exposed to the task runtime; no model slug or snapshot can be verified |
| Collaboration mode | Default mode; primary/root agent |
| Reasoning effort | Not exposed to the task runtime |
| Workspace execution mode | Unrestricted local filesystem/network access; approval policy `never` |
| Goal state | Complete |
| Goal-accounting usage at completion | 1,801,132 aggregate tokens and 10,747 seconds (2:59:07) |
| Usage scope | The implementation goal through v1 completion; this post-completion metadata update is not included because its goal counter is no longer available |

The aggregate token counter does not split input, cached input, reasoning, or
output tokens and is not an OpenAI invoice. The runtime also does not expose
the exact serving model. For an order-of-magnitude comparison only, the
historical [OpenAI API pricing](https://openai.com/api/pricing/) comparison used standard
GPT-5.6 Sol text rates of $5 per million input tokens and $30 per million
output tokens. Applying those rates to the 1.801132-million-token aggregate
gives:

- an illustrative 80% input / 20% output estimate of **$18.01**;
- an uncached standard-rate sensitivity envelope of **$9.01–$54.03** if the
  whole aggregate were priced at the input or output rate, respectively.

These figures are API-equivalent estimates, not the actual Codex subscription,
credit, or internal compute cost. They exclude cache discounts, long-context
premiums, service-tier adjustments, tool-call charges, and any usage that was
not represented by the completed goal counter.

At the historical DL snapshot, the post-v1 research goal was incomplete. A runtime snapshot on
2026-10-22 at 09:38 UTC returned 17,290,841 aggregate goal tokens and 90,256
goal-accounting seconds with tracker status `usageLimited`. These cumulative
counters are not Phase DL-only usage, have not advanced between the observed
DL queries, are not a complete measure of subsequent DL work, and do not provide input/cache/reasoning/
output categories or account billing. The current runtime instruction
identifies Codex based on GPT-5; an earlier research snapshot identified the
GPT-6 family. Neither exposes a verifiable exact serving snapshot or
reasoning-effort setting. Collaboration remains Default mode with one primary
agent and unrestricted local filesystem/network access. DL used two temporary
subagents for metadata tests and read-only review; exact serving versions are
not exposed. No new pricing estimate is
attached to these research counters.
That historical snapshot was not an `active` or `complete` tracker report.
The DL audit records four incomplete original-SPEC requirement groups.

The DM continuation's runtime snapshot (2026-10-22, 16:44 UTC) instead
reports tracker status `active`, 18,595,692 aggregate goal tokens and 94,848
goal-accounting seconds. These are cumulative, not DM-only counters, and do
not expose token categories or billing. Default mode, unrestricted local
filesystem/network execution and approval policy `never` remain exposed.
The declared family is Codex based on GPT-5; exact model/version and reasoning
effort remain unavailable. DM uses one primary agent and three temporary
subagents for plotting/validation tests, independent review/geometry tests,
and independent numerical checking. No serving-version attribution or new
API-equivalent price calculation is inferred. The repeated primary numerical
artifact, denser-grid evidence and combined checks are verified; branch root
and CI gates pass and PR #119 is merged. Exact pushed-main closure gates are
repeated before handoff. The original-SPEC goal is not complete: R31 and the
final whole-SPEC audit remain required.

DN's runtime snapshot on 2026-10-22 at 18:37 UTC reports tracker status `active`,
19,536,447 aggregate goal tokens and 97,600 goal-accounting seconds, cumulative
rather than DN-only. Exact serving model/version, reasoning effort, token
categories and billing remain unavailable. Default mode, declared Codex/GPT-5
family and unrestricted execution are unchanged. Three temporary adapter
workers contributed partial implementation before reporting an account usage
limit; the primary agent continues integration and verification locally. No
new cost calculation or serving-version attribution is inferred. DN's full
local and push/PR CI gates pass; R31 is supported and PR #120 is merged. Exact
pushed-main closure gates repeat before handoff. DN's historical working ledger
reported forty of forty supported; DO's independent audit subsequently reopens
the required residual-adaptive census, so that historical count is not current.

DO's runtime snapshot on 2026-10-22 at 18:59 UTC reports tracker status
`active`, 20,178,930 aggregate goal tokens and 98,878 goal-accounting seconds.
These cumulative counters are not DO-only, do not expose token categories,
and are not billing. The exposed declaration remains Codex based on GPT-5,
Default mode, unrestricted local filesystem/network and approval policy
`never`. Exact serving model/version and reasoning effort remain unavailable.
One primary agent and three temporary agents perform independent baseline,
rate and numerical audits, then bounded repairs/review; no exact serving
version or new price calculation is inferred. The required adaptive experiment
remains open, and the full goal remains active rather than complete.

The resumed DO snapshot on 2026-10-23 at 13:16 UTC reports `active`,
21,175,795 cumulative aggregate goal tokens and 99,657 goal-accounting seconds.
It does not include an input/cache/reasoning/output split, exact serving model,
reasoning effort or billing. Root's exposed family declaration remains GPT-5
based Codex in Default mode; no serving identity for temporary workers is
inferred. Two workers reported account usage limits after contributing audit
and repair evidence; the primary continued local integration. These counters
are not DO-only usage and no new API-equivalent cost is calculated.

The DO merge snapshot on 2026-10-23 at 13:35 UTC reports `active`,
21,505,132 cumulative aggregate goal tokens and 100,758 goal-accounting seconds.
All identity/category/billing limitations above remain: declared root family
Codex/GPT-5, Default mode and unrestricted local execution are exposed; exact
serving version, reasoning effort, token splits and account billing are not.
No new cost estimate is inferred. The verified phase is merged; the full goal
stays active because the adaptive experiment and its closure audit are absent.

DP's runtime snapshot on 2026-10-23 at 20:00:53 UTC reports `active`,
23,678,436 cumulative aggregate goal tokens and 109,728 goal-accounting seconds.
These are not DP-only counters, token-category splits or billing. Exposed
declaration remains GPT-5-based Codex, Default mode, unrestricted local
filesystem/network and approval policy `never`; exact serving version and
reasoning effort remain unavailable. One primary agent and three temporary
workers contribute numerical/controller, geometry/test and independent-checker
work and reviews. No exact serving identities or new API-equivalent price are
inferred. Both full sampled runs pass and repeat bytes; the goal remains active
pending shipping and whole-scope closure; independent checking subsequently passes.

The DP merge snapshot on 2026-10-23 at 20:20:53 UTC reports `active`,
23,852,175 cumulative aggregate goal tokens and 110,928 goal-accounting seconds.
These are not phase-only usage, token-category splits or billing. Declared
root identity remains GPT-5-based Codex, Default mode, unrestricted local
filesystem/network and approval policy `never`; exact serving version and
reasoning effort are not exposed. The three temporary workers contribute
controller/numerical QA, geometry/tests and independent checks/documentation.
No new price or serving-version inference is made. DP is verified and merged;
the original-SPEC goal remains active pending its fresh closure audit.

DQ's source-audit snapshot on 2026-10-24 at 01:25:33 UTC reports `active`,
24,576,774 cumulative aggregate goal tokens and 114,384 goal-accounting seconds.
These are not DQ-only counters and expose neither token categories nor billing.
Root is declared GPT-5-based Codex in Default mode, with unrestricted local
filesystem/network and approval policy `never`; exact serving model/version
and reasoning effort remain unavailable. Three attempted temporary reviewers
fail at the account usage limit without contributing DQ work. The primary
continues the source-level audit/repairs and release gates locally. No worker
identity, new price estimate or external peer review is inferred. The goal
remains active until the actual DQ release and merged-main gates pass.

DQ's code-head shipping snapshot on 2026-10-24 at 01:35:20 UTC reports
`active`, 24,675,942 cumulative aggregate goal tokens and 114,959
goal-accounting seconds. These counters are not phase-only, do not expose
input/cache/reasoning/output categories and are not account billing. Declared
root family, Default mode and unrestricted/never-approval execution remain as
above; exact serving version and reasoning effort are unavailable. No new
serving identity or API-equivalent cost is inferred. Code-head root/CI/archive/
install gates pass, while final documentation-head and merged-main closure
remain. The goal is still active at this point-in-time snapshot.

The DQ merge snapshot on 2026-10-24 at 01:41:43 UTC reports `active`,
24,701,811 cumulative aggregate goal tokens and 115,358 goal-accounting seconds.
It is captured after the verified PR merge and before exact-main closure;
it is not the final completion counter or DQ-only usage. No input/cache/
reasoning/output split, account billing, exact serving version or reasoning
effort is exposed. Declared root identity remains GPT-5-based Codex, Default
mode, unrestricted local filesystem/network and approval policy `never`.
The primary finishes local release checks; the three attempted DQ reviewers
contribute no work after their account-limit failures. No new serving-model
identity or API-equivalent price is inferred. The tracker is marked complete
only after the final pushed-main checks, which are outside this point-in-time
snapshot; these operational counters do not prove a theorem or scope decision.

### Completed original-SPEC research goal

The terminal goal record in the preceding completion response was updated on
2026-10-23 at 21:54 EDT (2026-10-24 01:54 UTC). It exposed the following frozen
accounting record, which supersedes the
earlier active snapshots for final goal state, not their historical usage.

| Field | Completed original-goal record |
| --- | --- |
| Goal state at completion | `complete` |
| Verified closure revision | `e8c35eda7662566639c5eac6748d3032df87075f` |
| Exact-main CI | [40051138962](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40051138962), successful Lean and Python jobs |
| Terminal aggregate usage | 24,738,403 cumulative aggregate goal tokens |
| Terminal goal-accounting time | 116,095 seconds (32:14:55) |
| Declared root model identity | GPT-5-based Codex; exact serving model/version unavailable |
| Execution mode | Default collaboration; unrestricted filesystem/network; approval policy `never` |
| Unavailable fields | Reasoning effort, input/cache/reasoning/output token categories and account billing |
| Scope of counters | Cumulative completed research goal, not DQ-only usage, total calendar duration or this later documentation task |

DR is a post-completion documentation/finalization task. Its live goal query
returns null; the original completion result above is preserved from the
actual prior tool response, not fabricated from a current tracker. No DR-only
usage, new cost estimate or account bill is inferred. The historical v1
API-equivalent calculation above remains explicitly historical and is not
applied to this terminal aggregate. DR uses the primary agent and a read-only
documentation reviewer; that review does not imply new independent theorem
or numerical verification, an exact worker model identity or external peer review.

### Later practical application development: runtime provenance

On 2026-10-25 UTC, the current runtime describes the primary as Codex based on
GPT-6 in Default collaboration mode, with unrestricted filesystem/network and
approval policy `never`. It does not expose a verifiable exact serving model
snapshot, reasoning-effort setting, input/cache/reasoning/output token counts
or account billing. The current goal query returns null, so no new aggregate
goal usage, elapsed-goal total or phase-only usage is available. The frozen
completed original-goal record above is preserved; it excludes this later
application development and is not extrapolated or repriced.

The same provenance limits apply to the completed 1.1.0 publication and DZ
repository handoff. Publication timestamps and individual command/CI durations
are observed operational records, not new aggregate implementation accounting.
Successful OIDC upload, public artifact hashes and clean installations establish
publication; no continuing authenticated PyPI browser session is claimed.

Temporary agents contribute implementation, native-boundary tests, packaging
review, benchmark work and upstream recipe/methods preparation. Earlier
account-limit failures are not counted as reviews; later successful reviews
and tests are identified in the phase PR receipts. Exact worker serving
versions and billing are unavailable. Local/CI durations are recorded only for
the measured commands in those receipts, not treated as full-project calendar
duration. No new API-equivalent price or actual account cost is inferred.

## Phase status

States track closed phase contracts. Delivery counts and then-open research
gaps in older rows are historical snapshots; the current scope decision and
later results above take precedence.

| Phase | State | Delivery |
| --- | --- | --- |
| 0 — bootstrap | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| A — exact GECP core | complete | [PR #4](https://github.com/tripp-smith/gecp-kernel-structure/pull/4) |
| B — PSD baseline | complete | [PR #5](https://github.com/tripp-smith/gecp-kernel-structure/pull/5) |
| C — fermionic structure/census | complete | [PR #6](https://github.com/tripp-smith/gecp-kernel-structure/pull/6) |
| D — separated approximation | complete | [PR #7](https://github.com/tripp-smith/gecp-kernel-structure/pull/7) |
| E — structural GECP | complete | [PR #8](https://github.com/tripp-smith/gecp-kernel-structure/pull/8) |
| F — sparse spectral application | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| R — release readiness | complete | [PR #9](https://github.com/tripp-smith/gecp-kernel-structure/pull/9), [closure PR #10](https://github.com/tripp-smith/gecp-kernel-structure/pull/10) |
| S — realistic synthetic applications | complete | implementation `6404311`; canonical evidence delivered directly to `main` |
| G — fermionic GECP rate | complete; bounded historical contract | merged [PR #11](https://github.com/tripp-smith/gecp-kernel-structure/pull/11) |
| BC — product-measure Andréief identity | complete; green CI | merged [PR #57](https://github.com/tripp-smith/gecp-kernel-structure/pull/57) |
| BD — arbitrary-sample and Taylor-tail weak signs | complete; green CI | merged [PR #58](https://github.com/tripp-smith/gecp-kernel-structure/pull/58) |
| BE — strict all-orders Taylor-tail signs | complete; green CI | merged [PR #59](https://github.com/tripp-smith/gecp-kernel-structure/pull/59) |
| BF — factorial-free Taylor-tail GECP envelope | complete; green CI | merged [PR #60](https://github.com/tripp-smith/gecp-kernel-structure/pull/60) |
| BG — positive-series compatible fermionic split | complete; green CI | merged [PR #61](https://github.com/tripp-smith/gecp-kernel-structure/pull/61) |
| BH — block-feature coefficientwise Cauchy--Binet | complete; green CI | merged [PR #62](https://github.com/tripp-smith/gecp-kernel-structure/pull/62) |
| BI — mixed-column coefficient bridge | complete; green CI | merged [PR #63](https://github.com/tripp-smith/gecp-kernel-structure/pull/63) |
| BJ — infinite-tail mixed coefficients | complete; green CI | merged [PR #64](https://github.com/tripp-smith/gecp-kernel-structure/pull/64) |
| BK — low-rank tail amplification | complete; green CI | merged [PR #65](https://github.com/tripp-smith/gecp-kernel-structure/pull/65) |
| BL — amplified-kernel envelope | complete; green CI | merged [PR #66](https://github.com/tripp-smith/gecp-kernel-structure/pull/66) |
| BM — total-positive Hadamard envelope | complete; green CI | merged [PR #67](https://github.com/tripp-smith/gecp-kernel-structure/pull/67) |
| BN — positive-series GECP composition | complete; green CI | merged [PR #68](https://github.com/tripp-smith/gecp-kernel-structure/pull/68) |
| BO — signed unit-band amplification | complete; green CI | merged [PR #69](https://github.com/tripp-smith/gecp-kernel-structure/pull/69) |
| BP — translated unit-band amplification | complete; green CI | merged [PR #70](https://github.com/tripp-smith/gecp-kernel-structure/pull/70) |
| CX — prefix-cone Hodge obstruction | complete; green CI | merged [PR #104](https://github.com/tripp-smith/gecp-kernel-structure/pull/104) |
| CY — adjacent-cutoff discriminant bridge | complete; green CI | merged [PR #105](https://github.com/tripp-smith/gecp-kernel-structure/pull/105); 3,602 Lean jobs and 123 tests |
| CZ — interior transfer-node obstruction | complete; green CI | merged [PR #106](https://github.com/tripp-smith/gecp-kernel-structure/pull/106); 3,603 Lean jobs and 126 tests |
| DA — one-pivot exponential Chebyshev sign | complete; green CI | merged [PR #107](https://github.com/tripp-smith/gecp-kernel-structure/pull/107); 3,605 Lean jobs and 129 tests |
| DB — adaptive transfer coverage obstruction | complete; green CI | merged [PR #108](https://github.com/tripp-smith/gecp-kernel-structure/pull/108); 3,608 Lean jobs and 133 tests |
| DC — all-orders local-minor contraction bridge | complete; green CI | merged [PR #109](https://github.com/tripp-smith/gecp-kernel-structure/pull/109); 3,610 Lean jobs and 140 tests; original companion target open |
| DD — common all-orders residual sign gauge | complete; green CI | merged [PR #110](https://github.com/tripp-smith/gecp-kernel-structure/pull/110); common gauge, actual Schur bound, original-run reflection; 3,614 Lean jobs and 146 tests |
| DE — cutoff-uniform local panel determinant budget | complete; green CI | merged [PR #111](https://github.com/tripp-smith/gecp-kernel-structure/pull/111); actual-kernel local order `16(33+4q)`; 3,615 Lean jobs and 151 tests |
| DF — backwards witnesses and global finite rate | complete; green CI | merged [PR #112](https://github.com/tripp-smith/gecp-kernel-structure/pull/112); actual-kernel bound `2^(-q)` after `32(s+1)(33+4q)` complete pivots; 3,619 Lean jobs and 155 tests |
| DG — actual real-cutoff rate and parameter boundary | complete; green CI | merged [PR #113](https://github.com/tripp-smith/gecp-kernel-structure/pull/113); explicit logarithmic rate / tolerance and actual first-pivot obstruction; 3,622 Lean jobs and 159 tests |
| DH — finite-run and physical-error closure | complete; green CI | merged [PR #114](https://github.com/tripp-smith/gecp-kernel-structure/pull/114); no future-pivot premise; exact stop, selected cross and actual Green-integral difference; 3,625 Lean jobs and 163 tests |
| DI — original-SPEC audit and bounded closures | complete; green CI | merged [PR #115](https://github.com/tripp-smith/gecp-kernel-structure/pull/115); 40 groups / seven incomplete; real D1, compact maxima, exact certificates; 3,627 Lean jobs, 838 audited declarations, 177 tests and wheel smoke |
| DJ — original-kernel fill-distance control | complete; green CI | merged [PR #116](https://github.com/tripp-smith/gecp-kernel-structure/pull/116); arbitrary-domain PSD/order and compact diagonal maxima; original-kernel `2*L*h` bound; R09 closed; 3,630 Lean jobs, 857 audited declarations and 180 tests |
| DK — strict-SPD complete-pivot reduction | complete; green CI | merged [PR #117](https://github.com/tripp-smith/gecp-kernel-structure/pull/117); fresh strict Schur congruence; every complete pivot diagonal; independently defined actual-run Cholesky equivalence; compact/fill corollaries; 3,633 Lean jobs, 876 audited declarations and 183 tests |
| DL — complete numerical result metadata | complete; green CI | merged [PR #118](https://github.com/tripp-smith/gecp-kernel-structure/pull/118); actual/requested precision, conditioning/ties/certification scopes, returned-error convergence and optimizer/failure diagnostics; schema-2 summaries; 232 tests; clean dependency audit; isolated wheel smoke |
| DM — complete census diagnostics and normalized curves | complete; green CI | merged [PR #119](https://github.com/tripp-smith/gecp-kernel-structure/pull/119); byte-repeated 28-case MP128 census, complete diagnostics, normalized plots, actual representative resolution ladder, 30 MP192 core checks, 430 tests and isolated wheel smoke |
| DN — fair six-method spectral comparisons | complete; green CI | merged [PR #120](https://github.com/tripp-smith/gecp-kernel-structure/pull/120); byte-repeated 192-fit JSON/PNG, independent error reconstruction, 513-test full root gate, installed-wheel smokes including SciPy 1.15 |
| DO — independent final SPEC audit | complete; green CI | merged [PR #121](https://github.com/tripp-smith/gecp-kernel-structure/pull/121); 569-test root, 3,634 Lean jobs, 876 axiom reports, installed-wheel checks including minimum SciPy 1.15.0; R25 still incomplete |
| DP — residual-adaptive seven-cutoff census | complete; green exact-main CI | merged [PR #122](https://github.com/tripp-smith/gecp-kernel-structure/pull/122); repeated seven-case/28-prefix data, 21 independent core/28 cross checks, 713-test root, clean archives/two installed environments; main closure 5a5fef0 and CI 40018905634 passed |
| DQ — original-SPEC closure audit | complete; green exact-main CI | merged [PR #123](https://github.com/tripp-smith/gecp-kernel-structure/pull/123); exact main e8c35ed, 757 tests in 150.31 seconds, successful CI 40051138962, clean archives and two installed environments |
| DR — final documentation and repository handoff | verified; green exact-source CI | merged [PR #124](https://github.com/tripp-smith/gecp-kernel-structure/pull/124); 759-test root, CI 40054146369/40054160333, clean archives/two installed environments; synchronized completed-run evidence and no mathematical/numerical changes |
| H — verification integrity refresh | complete locally; green CI | merged [PR #12](https://github.com/tripp-smith/gecp-kernel-structure/pull/12) |
| I — determinantal sign-coherence bridge | complete; green CI | merged [PR #13](https://github.com/tripp-smith/gecp-kernel-structure/pull/13) |
| J — strict-sign-regularity transfer | complete; green CI | merged [PR #14](https://github.com/tripp-smith/gecp-kernel-structure/pull/14) |
| K — exponential total-positivity scaffolding | complete; green CI | merged [PR #15](https://github.com/tripp-smith/gecp-kernel-structure/pull/15) |

## Completed tasks

- T-001–T-003: Codex-native skills, pinned toolchains, repository skeleton,
  CI, axiom audit, and root verification command.
- T-004–T-005: exact GECP residual/run core plus finite-grid and
  arbitrary-precision numerical GECP.
- T-006: positive-semidefinite residual preservation and canonical recursive
  GECP/pivoted-Cholesky equivalence.
- T-007–T-008: fermionic identities, derivatives, endpoint-resolved grids,
  interval pivot certification, and deterministic census.
- T-009–T-010: explicit exponential-family and fermionic separated
  approximation theorems with an independent high-precision evaluator.
- T-011–T-012: exact structural census, minimized sign-regular obstruction,
  and cross-ratio sufficient condition.
- T-013–T-015: discrete-to-continuous certificate bridge, Green-function
  error transfer, and sparse two-atom application.
- T-016: v1.0.0 metadata, release smoke configuration, artifacts, status
  audit, and final handoff.
- T-017: typed synthetic-application API, Hubbard-like and gapped Green-
  function compression, clean/noisy sparse recovery, PSD covariance landmark
  selection, deterministic JSON, summary plot, and application documentation.
- T-018: variable-factor/block contraction infrastructure, exact one-step
  obstruction, exact first-two-pivot localization for every positive cutoff,
  rank-two secant/interpolation analysis, and the proved two-step half
  contraction on `0<Λ≤1`.
- T-019: interpreter-bound local/CI pytest entry point, refreshed audited
  development lock, and maintained-plus-fresh-environment root verification.
- T-020: bordered selected-core determinant recurrence and residual identity,
  equivalent four-minor sign criterion, exact selected-cross diagnostics, and
  a generic minimized failure witness.
- T-021: ordered-minor strict-sign-regularity definitions, orientation-aware
  transfer to bordered-minor and selected-residual coherence, positive column
  scaling, fermionic reduction, order-one/order-two proofs, and exact
  permutation regressions.
- T-022: positive row scaling, ordered Vandermonde determinant positivity,
  positive determinant of the principal exponential Taylor block, and an exact
  rational factorization regression.
- T-023: sparse-polynomial positive-root bounds, generalized Vandermonde
  nonsingularity, and positive-node determinant positivity.
- T-024: generalized Vandermonde nonnegativity at a zero boundary node by
  uniform positive shifts and determinant continuity.
- T-025: rectangular Cauchy--Binet expansion, finite exponential Taylor
  determinant positivity, and exact rational multi-term regressions.
- T-026: fixed principal determinant lower bound, entrywise exponential-series
  convergence, determinant-limit positivity, and a 100-decimal regression.
- T-027: arbitrary-node exponential determinant positivity, finite-reversal
  sign, all-orders exponential and fermionic strict sign regularity, and
  unconditional selected-residual sign coherence.

## Deferred research

No v1 implementation task is blocked. DN verifies R31's six-method comparison;
DP supplies the seven-cutoff residual-adaptive experiment reopened by DO.
DQ performs and verifies the subsequent original-SPEC closure audit and merges
its real boundary repairs in PR #123. No required original-scope research task
remains at the final verified handoff.
DG/DH have since closed the exact complete-pivot rate by a different route.
The historical v1 follow-ups below are optional alternative proofs or extensions:

- as an alternative proof, transport the proved `Λ≤1` two-corner half
  contraction through dyadic frequency layers, using the now-unconditional selected-pivot
  nonexpansiveness invariant between strict contractions, to obtain the
  `O(log Λ)` block length required by Conjecture G1;
- formalize the Gimbutas–Marshall–Rokhlin selected-exponential Chebyshev
  construction as an alternate low-rank theorem;
- formalize stronger published pivoted-Cholesky Lipschitz rates;
- extend sparse recovery to multidimensional kernels; this remains optional.

## Public Lean results

The principal completed-goal interfaces are exported through
`GECPKernelStructure/Theorems.lean` and consumed by
`scripts/CompletionContract.lean`:

- `fermionicKernel_finiteRun_error_le_exp`, the actual finite exact-GECP rate
  for real `Λ≥1`, without a convergence or future-pivot premise;
- `fermionicKernel_finiteRun_cross_error_le_exp` and
  `fermionicKernel_finiteRun_greenError_le_exp_interval`, the original-core
  cross and physical integral consequences;
- `fermionicKernel_separatedApprox_real`, constructive real-parameter D1;
- `fermionicKernel_finiteRun_error_le_allCutoff_exp` and
  `fermionicKernel_smallCutoff_logarithmic_obstruction`, the separate
  all-positive-cutoff variant and literal-parameter counterexample;
- `posDef_completeRun_iff_pivotedCholesky` and
  `posDef_completeRun_error_le_original_fillDistance`, arbitrary-domain
  strict-SPD complete-run equivalence and original-kernel fill control.

Additional public results from the historical phases include:

- exact GECP: `gecp_interpolates_selected_rows`,
  `gecp_interpolates_selected_cols`, `gecp_core_nonsingular`, and
  `gecp_core_det_eq_prod_pivots`;
- PSD baseline: `posDef_gecp_residual_posSemidefinite`,
  `posDef_gecp_residual_abs_le_diag`,
  `posDef_gecp_residual_sup_eq_diag_sup`,
  `gecp_eq_pivotedCholesky_of_posDef`, and
  `posDef_gecp_error_le_fillDistance`;
- fermionic structure: `fermionicKernel_reflection`,
  `fermionicKernel_centered`, `fermionicKernel_continuous`, and the exact
  derivative/bound theorems;
- low rank: `expFamily_separatedApprox` and
  `fermionicKernel_separatedApprox`;
- structural outcome:
  `signRegular_not_sufficient_for_parameterIndependent_pivots`,
  `residualUpdate_le_of_crossRatioControl`, and
  `gecp_error_le_geometric_of_crossRatioControl`;
- Phase G structure:
  `gecp_error_le_product_of_crossRatioControl`,
  `gecp_error_le_dyadic_of_crossRatioProduct`,
  `fermionicKernel_no_uniform_firstStep_contraction`,
  `fermionicKernel_firstPivot_abs_le_reflectedCorner`, and
  `fermionicKernel_twoCornerResidual_le_half_initial`;
- Phase I structure: `Run.borderedCore`,
  `Run.borderedCore_det_eq_selectedCore_det_mul_finalResidual`,
  `BorderedMinorSignCoherent`, and
  `pivotCrossProductSignCoherent_iff_borderedMinorSignCoherent`;
- Phase J structure: `StrictSignRegularAtOrder`, `StrictSignRegular`,
  `oriented_minor_pos_of_strictSignRegularAtOrder`,
  `strictSignRegular_borderedMinorSignCoherent`,
  `strictSignRegular_pivotCrossProductSignCoherent`,
  `strictSignRegularAtOrder_columnScale`,
  `expKernel_strictSignRegularAtOrder_two`, and
  `fermionicKernel_strictSignRegularAtOrder_two`;
- Phase K structure: `strictSignRegularAtOrder_rowScale`,
  `vandermonde_det_pos_of_strictMono`, `expTaylorPrincipalMatrix_apply`, and
  `expTaylorPrincipalMatrix_det_pos`;
- Phase L structure: `signVariations_lt_card_support`,
  `positiveRoots_lt_card_support`, `generalizedVandermonde_det_ne_zero`,
  `geometric_generalizedVandermonde_det_pos`, and
  `generalizedVandermonde_det_pos`;
- Phase M structure: `generalizedVandermonde_det_nonneg`;
- Phase N structure: `Matrix.det_mul_rect_eq_sum_injective`,
  `Matrix.det_mul_rect`, `expTaylorMatrix_apply`, and
  `expTaylorMatrix_det_pos`;
- Phase O structure: `expTaylorMatrix_det_ge_principal`,
  `expTaylorMatrix_tendsto_expMatrix`,
  `expTaylorMatrix_det_tendsto_expMatrix_det`, and `expMatrix_det_pos`;
- Phase P structure: `Fin.sign_revPerm`,
  `expMatrix_det_pos_of_strictMono`, `expKernel_strictSignRegular`,
  `fermionicKernel_strictSignRegular`, and
  `fermionicKernel_pivotCrossProductSignCoherent`;
- application: `grid_sup_le_max_add_lipschitz`,
  `approxPivot_of_grid_certificate`, and
  `greenError_le_kernelError_mul_l1`.

The axiom audit reports only Lean defaults (`propext`, `Quot.sound`, and where
needed `Classical.choice`).

## Python and API results

The package exports stable `FermionicKernel` evaluation, `gecp`,
`gecp_matrix`, pivoted Cholesky, solve-based cross approximation, interval
pivot certificates, deterministic census tools, sparse spectral recovery,
piecewise Chebyshev interpolation, and the Lean-matched
`fermionic_dyadic_taylor_approximation`.

The research API additionally exports exact-form decimal evaluators
`fermionic_first_residual` and `fermionic_two_corner_residual`, census block
analysis, and focused interval-certified trajectories. The two-corner
evaluator independently checks the Lean base theorem; it is not used as the
source of that theorem.

It also exports `SyntheticApplicationConfig`,
`run_synthetic_application_suite`, typed result records,
`synthetic_spectral_density`, `fermionic_green_from_density`, quadrature
weights, and deterministic sensor-point generation. The application runner
combines existing public algorithms; it does not bypass their convergence,
conditioning, or stop metadata.

The research submodule `kernelgecp.adaptive_grid` additionally provides
`AdaptiveGridConfig` and `run_adaptive_grid` without changing the established
45-name top-level export surface. Its explicit nested axes, stage decisions,
prefix diagnostics and untouched validation are documented in APPLICATION.
The controller is separate from interval-certified adaptive pivoting.

Normal non-convergence is represented by `converged` and `stop_reason`.
Certified-pivot budget exhaustion is never reported as certification, and
high-precision runs retain canonical decimal strings.

## Historical experiment commands

These commands describe earlier recorded runs, not fresh DR artifact
generation. Historical data retain their generation revision and manifests;
use the frozen revision for byte reproduction, or separate scratch output
paths for exploratory runs from a later checkout (see APPLICATION).

```bash
./scripts/verify.sh
uv build
uv run python -m kernelgecp.census \
  --config experiments/release-smoke.toml \
  --output /tmp/gecp_release_smoke.jsonl
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/data/synthetic_applications.json \
  --plot experiments/data/synthetic_applications.png
git status --short
```

## Verification results

The completed original-goal baseline is exact main `e8c35ed`: 757 tests in
150.31 seconds, all 3,634 Lean jobs, 876 permitted axiom reports, 22 consumer
examples, style/type/security gates, CI `40051138962`, clean archives and both
installed environments passed. DR's later documentation checks are recorded
in the phase record and do not relabel that historical result.

### Historical verification snapshot (126 tests)

- Lean build and public axiom audit pass; no `sorry`, `admit`, or custom axiom
  is present.
- Ruff, formatting, strict mypy, and all 126 pytest cases pass.
- `pip-audit` reports no known vulnerabilities; the local editable project is
  correctly skipped because it is not a PyPI dependency.
- The canonical 21-case census is byte-reproducible with SHA-256
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
- The exact-surrogate dataset has SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The 128-bit release smoke run converges at tolerance `1e-6` with rank 6 for
  cutoff 1 and rank 82 for cutoff `10⁶`. Its two-record JSONL has SHA-256
  `d7b5f53d61eb2ea49cd1ebde944b47a06f931b9d82687c8fa4458ac4c09907d5`.
- The synthetic-application JSON is byte-identical across repeated runs and
  has SHA-256
  `f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2`;
  the reviewed plot has SHA-256
  `0529263affd2fe10773f71784f06007650bff833374fbbebbb16f2507b2ab169`.
- At cutoff 8 and tolerance `1e-8`, a 12-pivot universal GECP basis gives
  Green errors `4.16e-10` and `2.91e-10` for the Hubbard-like and gapped
  densities, below the common `5.40e-9` discrete transfer bound.
- The known-transition fixture is recovered with four atoms to `3.33e-16`
  held-out error. The historical noisy support-augmented scan uses eight
  effective atoms and reaches `8.56e-9` held-out function error. Its dictionary
  includes the generating frequencies, so it is not a blind benchmark or
  evidence for atom identification.
- GECP and pivoted Cholesky agree on 31 of 72 covariance landmarks; the
  selected cross has maximum error `7.78e-7` and relative Frobenius error
  `1.14e-7`.

## Known limitations and non-claims

No required original-SPEC deliverable remains at the verified final handoff.
DQ's source-level audit and real boundary repairs are verified and merged;
exact-main `e8c35ed` passed all gates, including CI `40051138962`.
DP's seven-cutoff data pass all
28 sampled tolerance cases; that is not a continuous certificate. Phase-specific proof-route
non-claims below preserve their historical contracts; they are not additional
missing original-SPEC requirements, and later DG/DH rate theorems do not depend
on completing those abandoned routes.

- No unqualified all-positive-cutoff `log(1+Λ)` solution or workshop-priority
  adjudication is claimed. The exact `Λ≥1` theorem and the literal small-cutoff
  counterexample above are the actual mathematical outcomes.
- The low-rank theorem is not a GECP theorem.
- Census curves and pivot certificates are evidence, not a continuous GECP
  convergence proof.
- The proved `0<Λ≤1` two-step half contraction alone is not extrapolated to
  larger cutoffs. DG/DH separately prove the exact complete-pivot logarithmic
  rate for `Λ≥1` and the all-positive-cutoff `log(2+Λ)` variant. Extending that
  rate to arbitrary certified approximate pivots remains a conditional research
  extension, not a delivered practical convergence guarantee.
- Phase BR constructs the exact compatible low-rank-plus-tail split for every
  fermionic prefix and preserves the rank-eight feature factor. It does not
  prove the conditioned-tail amplification inequality: exact arithmetic shows
  that the unconditioned coefficientwise-nonnegative premise can fail after
  conditioning.
- Phase BS proves the sign-free finite-weight amplification reduction, but it
  does not yet prove the nine amplified companion determinant envelopes. Those
  envelopes are the remaining premise for the order-eighteen half contraction.
- Phase BT proves the exact two-sided companion factorization and exposes the
  matched-minor interface, but the observed common orientation of those
  conditioned minor products is still numerical evidence rather than a
  theorem.
- Phase BU proves the exact bordered-determinant sign equivalence, but it does
  not make the resulting products sign coherent; Phase BV gives an exact
  counterexample with two negative products among 56 coordinates.
- Phase BV proves the aggregate graded-pairing identity, but the observed
  positivity of that mixed-column coefficient is not yet a theorem for
  arbitrary fermionic prefixes.
- Phase BW proves a single-determinant linearization at every required weight,
  but its exact signed four-cycle obstruction rules out making the bordered
  matrix totally nonnegative by row and column sign changes.
- Phase BY formally identifies arbitrary recursive `Run.conditionBoth` with
  the Phase BX selected-core Schur expression and specializes the
  double-border factorization to arbitrary finite samples. The normalized
  double-border envelope remains open.
- Phase BZ proves that those double borders are distinguished upper-tangent
  minors and that their normalized bounds feed Phase BS directly. The tested
  cross-inversion coefficient orientation is not yet a theorem, and the
  normalized envelope remains open.
- Phase CA proves the exact double-flag Cauchy--Binet expansion and the
  factorwise orientation-composition interface. It does not prove the three
  exponential flag-minor sign families or the normalized envelope.
- Phase CB proves the tangent middle factor's exact interleaved cell
  decomposition and the conditional shuffle-sign transfer.
- Phase CC proves all-orders TN of that interleaved middle matrix and removes
  the conditional premise. Both outer flag families remain open.
- Phase CD proves that those outer families cannot be oriented independently
  to satisfy the Phase CA termwise criterion: a negative double-flag term
  occurs even though the complete tangent determinant is positive. This is a
  proof-interface obstruction, not a counterexample to the target envelope.
- Phase CX proves that the full conditioned matched-coordinate Hessian need
  not have one positive direction. Any completion must exploit the contiguous
  degree-prefix cone; a global Lorentzian-polynomial argument is unavailable.
- Phase CY reduces the adjacent-cutoff sign to cumulative projected cross
  determinants. Their sign agreement is verified on 9,170 finite interior
  samples. Phase CZ gives an exact strict-interior finite-Taylor failure, so
  the universal fixed-node premise cannot complete this route. Adaptive
  neighboring-node coverage and the infinite-tail statement remain open.
- All-orders exponential and fermionic minor signs, and selected-residual
  nonexpansiveness, do not imply a strict contraction factor below one.
- Full exponential determinant positivity is currently stated for
  nonnegative strictly increasing node tuples. The translation/scaling and
  frequency-reversal orientation needed for arbitrary ordered tuples and
  `expKernelSignature` are not yet assembled.
- Arbitrary PSD tie-breaking is not claimed equivalent to diagonal pivoting;
  the formal theorem uses the documented diagonal-preferring canonical rule.
- The formal separated construction is dyadic truncated Taylor, not the
  selected-exponential Chebyshev variant.
- The post-v1 spectral fixtures are synthetic analogues, not material-
  specific calculations, experimental fits, or uncertainty-quantified
  analytic continuation.
- Accurate held-out Green-function recovery in the noisy dense scan does not
  identify the four generating atoms uniquely; the eight recovered atoms are
  an effective representation.
- Matching covariance pivots do not prove optimal sensor placement.

## Numerical stability and performance

- Kernel evaluation uses sign-dependent overflow-safe formulas.
- GECP supports float64, exact rational surrogates, and arbitrary precision;
  the canonical census uses 128-bit arithmetic.
- Numerical cross approximation uses solves/factorizations rather than an
  explicit inverse.
- Every GECP result records pivots, residual history, conditioning, ties,
  convergence, and stopping metadata.
- The 28-case schema-2 census converges in every tracked case; sampled rank at
  tolerance `1e-12` grows from 9 at cutoff 1 to 146 at cutoff `10⁶`. Grid
  coordinates originate at 53 bits and elimination uses 128 bits. The
  historical 21-case `1e-10` rank range 8–124 remains preserved.

## Security and dependency notes

The Lean revision and Python dependency graph are pinned. Runtime algorithms
need no network access, and the repository stores no credentials, personal
data, or telemetry. Generated environments, caches, and package builds remain
ignored.

## Files and artifacts delivered

- Lean library: `GECPKernelStructure`.
- Python package: [kernelgecp v1.1.0](https://pypi.org/project/kernelgecp/1.1.0/).
- Canonical data and plot under `experiments/data/`.
- DP explicit-axis adaptive JSON, independent checker report and normalized/
  stage PNG under `experiments/data/`; independent report passes.
- Synthetic application configuration, canonical JSON, and reviewed summary
  plot under `experiments/` and `experiments/data/`.
- Public release artifacts: `kernelgecp-1.1.0-py3-none-any.whl`,
  `kernelgecp-1.1.0-cp312-abi3-manylinux_2_28_x86_64.whl`,
  `kernelgecp-1.1.0-cp312-abi3-macosx_11_0_arm64.whl` and
  `kernelgecp-1.1.0.tar.gz`; each matches the
  [release manifest](https://github.com/tripp-smith/gecp-kernel-structure/releases/download/v1.1.0/release-manifest.json).
  Its SHA-256 is `2d2d9b2cb7cb1c35f7762cc079aa3722ad4e335434d57cbe3c98ba4894124420`.
  The same GitHub release includes complete locked native sources, notices,
  dependency inventory, source receipt, methods, citation and release notes.
- Claim and research documentation: `README.md`, `FINDINGS.md`, `RESEARCH.md`,
  `APPLICATION.md`, `MATHLIB.md`, and `SPEC.md`; practical release documentation
  adds `PACKAGE_README.md`, `METHODS.md`, `CITATION.cff` and `RELEASE_NOTES.md`.
- Preserved workflow input: `autonomous-implementation.md`, explicitly labeled
  as a reference rather than current status and excluded from package builds.

## How to run locally

```bash
uv sync --all-extras
lake exe cache get
uv run python -c "import kernelgecp; print(kernelgecp.__version__)"
uv run python experiments/synthetic_applications.py \
  --output experiments/scratch/synthetic_rerun.json \
  --plot experiments/scratch/synthetic_rerun.png
```

See `APPLICATION.md` for GECP, certified-pivot, low-rank, and sparse-recovery
examples.

## How to test

```bash
./scripts/verify.sh
```

Focused commands are documented in `AGENTS.md`; a phase is considered verified
only after the root command passes.

## Final closure and optional next research

DO's independent audit and DQ's primary-agent closure in `FINAL_AUDIT.md`
supersede the earlier forty-supported-label
assumption in `COMPLETION_AUDIT.md`; DP addresses the R25 gap it identified.
DJ/DK prove original-kernel Lipschitz fill control, arbitrary-domain PSD
geometry and strict-SPD complete-pivot diagonal necessity through fresh
residual prefixes. DL completes result metadata; DM completes the census
fixed-grid depth/diagnostics/normalized curves and representative resolution checks,
not the required residual-adaptive grids.
DN verifies R31's six-method comparisons with shared observations and retained
failures. DO verifies correctness repairs and merges PR #121. DP's genuine
all-seven-cutoff residual-driven nested-grid experiment through `1e-12` passes
all sampled cases; data and reviewed plots repeat bytes. DQ performs the
fresh original-SPEC audit and repairs actual certificate/parameter boundary
defects. Final-branch root, package/install and exact-head CI pass; PR #123 is
merged and exact-main closure `e8c35ed` passed all gates, including CI
`40051138962`. No additional
required original-SPEC task remains. Approximate-pivot rate extensions remain
explicitly qualified research, not an inferred consequence of exact-pivot
theorems. The original-scope source/explicit-item audit, not DP's inventory
mapping, is the final scope evidence;
the final 100% estimate concerns that required scope only and is neither a
theorem-success probability nor a substitute for actual release gates.

Earlier companion/conditioned-feature routes remain recorded as alternative
research, including their exact obstructions. They are not additional
prerequisites for this different local-minor proof route and do not negate its
proved dyadic rate. No literature-priority claim is made.
