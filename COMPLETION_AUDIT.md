# Original SPEC completion audit

DO audit update (2026-10-22): the independent [FINAL_AUDIT.md](FINAL_AUDIT.md)
reads all thirty original sections and checks actual sources/artifacts. It
overturns the prior 40/40 conclusion: §14 requires residual-adaptive increasingly
dense grids, not merely the fixed-grid census and three resolution reruns.
R25 is reopened. DO also repairs certificate lower bounds, sampled adaptive
ties, a large-finite float64 failure path, centered-formula coverage and stale
current documentation; fresh verification is recorded at phase closure.
Historical phase counts below describe those phases, not the current decision.

Audit date: 2026-10-20. DI entry state: merged DH, `82dd80b7be7d6f65546112412c0623dd6093694f`.
DI additions were verified by the full root gate (3,627 Lean
jobs, 838 enforced axiom reports, 177 tests, no known dependency vulnerabilities)
and green push/PR CI, and merged as [PR #115](https://github.com/tripp-smith/gecp-kernel-structure/pull/115).
DJ adds the original-kernel fill-distance proof and arbitrary-domain PSD
geometry, verified by the full root gate (3,630 Lean jobs, 857 enforced axiom
reports, 180 tests, no known dependency vulnerabilities), green push/PR CI,
and merged as [PR #116](https://github.com/tripp-smith/gecp-kernel-structure/pull/116).
DK adds the strict-SPD complete-run reduction and independently defined
Cholesky equivalence, verified by the full root gate (3,633 Lean jobs,
876 enforced axiom reports, 183 tests, no known dependency vulnerabilities),
green push/PR CI, and merged as [PR #117](https://github.com/tripp-smith/gecp-kernel-structure/pull/117).
DL adds complete numerical result metadata and related failure-path repairs,
verified by the full root gate (3,633 Lean jobs, 876 enforced axiom reports,
232 tests, clean dependency audit), green push/PR CI, an isolated wheel smoke,
and merged as [PR #118](https://github.com/tripp-smith/gecp-kernel-structure/pull/118).
DM delivers in [PR #119](https://github.com/tripp-smith/gecp-kernel-structure/pull/119): the byte-repeated
28-case schema-2 census, actual representative order-32 runs, normalized
plots, and 30 passing 192-bit core checks. The full 430-test local root gate
passes; exact-head push and PR CI pass, and PR #119 merged on 2026-10-22.
The pushed main closure head repeats root/CI gates before final handoff.

DN verifies all six R31 methods, the byte-repeated 192-fit JSON/PNG and
independent reconstruction of every reported error. Final branch `eec016d`
passes the 513-test root gate and exact-head push/PR CI. Installed-wheel smokes
pass with both pinned SciPy 1.18 and minimum 1.15. PR #120 merged as
`a3186a89b27d92677ef0e3338b90f9048a28fdb5` on 2026-10-22 at 18:37:17 UTC.
Exact pushed-main closure gates repeat before handoff; final audit is next.

## Decision

**The full goal is not complete.** The principal exact fermionic GECP rate is
proved, including its finite-run cross and Green consequences. That does not
prove all of the original specification's other requirements.

The working ledger has 40 concrete requirement groups: **39 supported and 1 incomplete**.
R25 remains missing. DO's final implementation head passes the 569-test root
gate, exact-head push/PR CI and isolated package checks; PR #121 is merged.
Merged/pushed-main closure gates repeat before handoff. The inventory is 97.5%
checklist coverage, not 100%. A rough overall effort
estimate is **97%**, not a probability of theorem truth or an exact workload
fraction. The principal theorem is proved, but the residual-adaptive experiment
and subsequent final closure remain. Green tests and status labels cannot
substitute for those missing outputs.

“Supported” below means the named source was inspected at the relevant scope,
with its tests/artifacts/gates, not merely that a similarly named file exists.
Research-route questions and explicitly proposed interfaces are distinguished
from required mathematical/application outcomes. Original SPEC text is not
rewritten.

## Requirement ledger

| ID | Requirement and SPEC section | Status | Authoritative evidence / limitation |
| --- | --- | --- | --- |
| R01 | Original selected-core cross and iterative residual agree (§1,7,8) | Supported | `GECP/RunCrossApproximation.lean`: `Run.finalResidual_eq_sub_crossApproximation`; selected core invertibility is derived from stored successful pivots. `test_finite_run_closure.py` checks 504 entries with exact solves and actual-exponential precision checks. |
| R02 | Exact finite continuous GECP, complete maxima and exact stop (§7) | Supported | `GECP/Definitions.lean`, `Fermionic/CompleteRun.lean`, `Run.CompleteOn`; `CompletePivotExistence.lean` derives a maximum on the compact physical rectangle from finite-run continuity. `residual_eq_zero_of_complete_zero` handles a zero maximum without an update. Numerical grid/adaptive modes are separate. |
| R03 | Selected-row interpolation (§8,27A,22) | Supported | `GECP/Interpolation.lean`: `gecp_interpolates_selected_rows`; `test_gecp.py`, exact cross regressions. |
| R04 | Selected-column interpolation (§8,27A,22) | Supported | Same module: `gecp_interpolates_selected_cols`; selected column and cross checks. |
| R05 | Actual original-core determinant equals pivot product; nonsingularity (§8,27A,22) | Supported | `GECP/Determinant.lean`: `gecp_core_det_eq_prod_pivots`, `gecp_core_nonsingular`; no supplied LDU/factorization premise; exact two-pivot instance and rational determinant checks. |
| R06 | Finite PSD Schur preservation and diagonal domination (§9) | Supported | `PositiveDefinite/ResidualPSD.lean` constructs `Pᴴ M P`, proves the square bound and equivalence of diagonal/all-entry bounds; no residual PSD premise is substituted for preservation. |
| R07 | Canonical finite diagonal-preferring GECP/Cholesky trace equivalence (§9,27B) | Supported | `PositiveDefinite/PivotedCholesky.lean`: both trace conversions and `gecp_eq_pivotedCholesky_of_posDef`; `test_cholesky.py` covers Gaussian PSD examples and diagonal ties. This does not itself close R08. |
| R08 | Continuous-domain SPD baseline at the original kernel scope (§3A,9,27B) | Supported | DK's `completeRun_diagonal_of_posDef` starts from mathlib strict SPD of the original kernel on an arbitrary ambient type and actual `Run.CompleteOn`, with no diagonal preference or residual strictness premise. Injective rectangular congruence derives strict positivity on fresh samples through updates. `posDef_completeRun_iff_pivotedCholesky` equates actual GECP runs with independently defined recursive diagonal maximization. Residual PSD, squared diagonal domination, fresh strictness and compact continuous diagonal complete-maximizer existence are exported. Exact rational all-branch/prefix checks and a Lean Gram-factor-certified two-pivot instance exercise the theorem. Weak-PSD off-diagonal ties are explicitly excluded from necessity. |
| R09 | Fill-distance error derived from Lipschitz **original kernel**, not assumed residual regularity (§10) | Supported | DJ's `diagonalRun_error_le_original_fillDistance` derives `abs R(x,y) <= 2*L*h` for every actual finite diagonal run from original-kernel coordinatewise Lipschitz regularity, finite-sample PSD, and the selected-node cover. `K-R` PSD and the original-energy bound are proved, not supplied as premises. DK's `posDef_completeRun_error_le_original_fillDistance` derives the needed diagonality from original SPD and actual complete pivoting, removing that supplied premise. Exact rational prefix/repeated-sample checks and Lean concrete instances exercise the interfaces. |
| R10 | Fermionic reflection (§4.2,27C,22) | Supported | `Symmetry.lean`: `fermionicKernel_reflection`; float64/mp tests include signs and boundaries. |
| R11 | Centered representation and positive column weight (§4.1,27C,22) | Supported | `Centering.lean`: `fermionicKernel_centered`; denominator positivity and stable kernel evaluation. |
| R12 | Continuity, derivatives and physical bounds (§4,27C) | Supported | `Derivatives.lean`: joint continuity, exact coordinate derivatives, time/frequency bounds; `test_kernels.py` compares float64 and arbitrary precision, including extreme cutoffs. |
| R13 | All-orders ordered-minor signs and scaling (§4.3,13) | Supported | `ExponentialTotalPositivity.lean`, `SignRegularity.lean`, ordered Vandermonde/Taylor/limit modules; `expKernel_strictSignRegular`, `fermionicKernel_strictSignRegular`. Exact and high-precision minor regressions supplement the proof. |
| R14 | A structural invariant along actual GECP prefixes (§12,13,30) | Supported | DD's single all-orders residual gauge, singleton Fischer bound, DF/DH finite backwards marked-minor witnesses and local pivot count. Separate conditioned-feature sign conjectures are not assumed. |
| R15 | D1: all real `Λ≥1`, `0<ε<1`, explicit logarithmic separated rank (§5,6,27D) | Supported | DI's `fermionicKernel_separatedApprox_real`: actual Taylor-list witness, rank `≤128 log(1+Λ) log(128/ε)`, uniform error `≤ε`. It closes the real quantifier/rank interface absent from the earlier natural dyadic theorem. `test_real_separated_approx.py` independently checks parameter/rank arithmetic. |
| R16 | Near-optimal exact GECP cutoff dependence (§1,11,27E,30) | Supported | `FiniteRunRate.lean`: every finite exact complete run on the actual rectangle has error `≤1024 exp(-((log 2)^2/256) k/log(1+Λ))` for `Λ≥1`. No residual envelope, future pivot, independent realization or convergence premise is assumed. |
| R17 | Literal small-cutoff G1 quantifiers (§11,28) | Supported | `RealCutoffRate.lean` refutes the all-positive-cutoff interpretation with an actual first complete pivot for every proposed `C,c>0`. The `Λ≥1` theorem and all-positive-cutoff `log(2+Λ)` variant are separate statements. This is resolution by proved parameter boundary, not a claim the literal conjecture is true. |
| R18 | Lean exact small-grid sequence, minors and pivot pattern (§15) | Supported | DI's `exactTwoPivotRealRun_complete`, positive entries/full determinant, and exact pivot list `[2,1/2]`; original rational determinant example, all-orders exponential sign theorem, and `test_surrogates.py`. |
| R19 | Lipschitz discrete-to-continuous upper bound (§16) | Supported | `Discretization/Grid.lean`: `grid_sup_le_max_add_lipschitz` derives the cover correction from the stated Lipschitz/cover assumptions. It does not assert a numerically guessed constant is rigorous. |
| R20 | Multiplicative approximate-pivot certificate (§17) | Supported | `approxPivot_of_grid_certificate` derives the pointwise `η` inequality from a certified upper bound; `PivotCertificate.eta` and interval enclosure API expose the gap. See the separate scope note below: this is not an arbitrary-η version of R16. |
| R21 | Adaptive maximizer success and budget-failure paths (§16,24) | Supported | Rounded point values could exceed interval upper bounds under cancellation. DO's tested repair uses singleton interval lower magnitudes, conservative enclosure arithmetic and explicit failure handling, including nonzero endpoint underflow. Certificates require valid callbacks and concern the stored rounded residual model, not an ideal exact GECP run or separate float64 cross. |
| R22 | float64, arbitrary precision and exact surrogates (§15,24) | Supported | 53-bit matrix paths, `mpmath` bit settings and decimal-string `HighPrecisionData`, rational `Fraction` surrogates; sub-float64 pivot regression and independent precision checks. |
| R23 | Cross evaluation by solves, not explicit numerical inversion (agent rules) | Supported | `gecp.py`, `applications.py`, exact/MP cross regressions use linear solves. Source scan finds no numerical inverse calls. Formal inverse notation in selected-core identities is algebraic, not a numerical implementation. |
| R24 | Full result metadata: convergence, stopping, conditioning, precision, ties, certification (agent rules) | Supported | DL's typed precision/conditioning/optimizer/convergence metadata remain. DO replaces invented adaptive tie counts with observed sampled candidate counts and preserves rounded-residual scope through census serialization. Bounded-mantissa/exponent Schur updates avoid representable overflow and mixed-scale underflow, preserve selected interpolation zeros, and fail explicitly outside arithmetic range. No floating-point output is relabeled a proof. |
| R25 | Seven-cutoff residual-adaptive increasingly dense census through at least `1e-12` (§14,26,30) | Incomplete | DM's repeated 28-case MP128 data are fixed `pivot="grid"`, order 24. Order-32 reruns cover only `1,1000,1000000` and are not residual-adaptive. These valid sampled artifacts do not satisfy the adaptive-grid requirement. A separate genuine residual-driven nested-refinement experiment at all seven cutoffs remains required. |
| R26 | All specified census diagnostics (§14) | Supported | Schema-2 `CensusRecord` serializes complete sampled residual and signed original-core determinant histories, singular values, signed pivots/magnitudes, all four physical boundary distances, dyadic panels, reflected coordinates and nearest selected partners/distances/tolerance. MP strings, 53-bit coordinate versus 128-bit elimination precision, per-prefix metadata and clean revision/config/lock/runtime provenance are explicit. Thirty independent 192-bit core checks pass; full phase gates pass and PR #119 is merged. |
| R27 | Normalized `k/log(1+Λ)` residual-curve plot (§14) | Supported | `gecp_census_normalized.png` contains seven deepest compatible trajectories, rank zero included. Both new PNGs are rendered, inspected and byte-repeated; numerical-zero omission is explicit and no zeros occur here. No universal fitted law is inferred. Full phase gates pass; merged PR #119. |
| R28 | Deterministic machine-readable experiment dataset (§26) | Supported | `census.py` records config hash/revision and canonical JSONL; small configurations and both full 28-case canonical runs repeat bytes. DM records dependency/lock/backend and clean-generation provenance; its combined report binds primary/resolution inputs and checker/grid sources. The original 21-case data/hash remains separately labeled historical evidence. |
| R29 | Particular Green functions use fewer spectral atoms with explicit error (§19,27F) | Supported | Committed synthetic suite at cutoff 8: clean four-atom library recovery has held-out error `3.33e-16`, versus the same kernel's universal GECP rank 12 at `1e-8`; noisy case and two-delta tests are separately qualified. No sparse uniqueness or optimality claim. |
| R30 | Sparse exponential/spectral interface, OMP and nonlinear least squares (§19,20) | Supported | `SparseResult` frequencies/weights represent sums of `K(·,ω)`; `sparse.py` implements OMP and bounded nonlinear refinement; synthetic fixtures/tests expose training and held-out errors. |
| R31 | Listed Phase-II method comparisons (§20) | Supported | DN implements shared-observation AAA-derived, genuine vector-fitting and Prony adapters against universal GECP, bounded variable-projection nonlinear least squares and unrefined OMP. All 192 fits and reviewed PNG repeat bytes from clean `072ce10`; independent tests reconstruct every measured error. Full root passes 513 tests/3,633 Lean jobs/876 axiom reports; exact-head push/PR CI and installed-wheel checks, including SciPy 1.15, pass. PR #120 is merged. All methods are compared, including failed physical AAA fits; no requirement that every method win is substituted. A formal AAA implementation remains outside initial scope; no sparse optimality or support-uniqueness claim. |
| R32 | Public Python interface (§23) | Supported | All seven proposed imports exist; all 45 `__all__` exports load after DL's three metadata exports. The result's actual residual interface is `residual_matrix`/histories plus `evaluate_residual`, documented in APPLICATION, rather than the illustrative `result.residual` spelling. |
| R33 | Deterministic ties and visible failures (§24) | Supported | Lexicographic selection and candidate tie counts are qualified as sampled, not global multiplicities. DO tests certificate cancellation, exhausted budgets, large-finite and mixed-scale Schur updates, arithmetic-range failures and exclusion of already annihilated coordinates instead of a raw `IndexError` or duplicate-core false convergence. |
| R34 | Enumerated numerical coverage (§25) | Supported | Existing reflection/interpolation/product/cross/PSD/Cholesky/precision/rank/exact tests remain. DO adds independent mpmath centered-expression comparisons at both signs, zero, boundaries and frequency magnitude 1000; comparing the centered alias with itself was not evidence. |
| R35 | Trusted Lean/root/CI gate (§25, agent rules) | Supported | DI enforces the exact axiom allow-list and expected report coverage, rejects empty/malformed/custom/native-reduction reports, preserves Lean's failing exit status, and rejects source placeholders/custom axioms. `lake build`, Ruff/format/mypy/pytest/pip-audit remain root gates. |
| R36 | Build/install handoff and local-note hygiene | Supported | DI rebuilt wheel/sdist, inspected archive member names, excluded the two untracked duplicate note files and local caches, and ran the wheel in an isolated Python 3.12 environment outside this checkout. The first local archive containing note copies was removed; originals untouched; nothing published. |
| R37 | Pinned setup, license and prior attribution (§2,9,10,28) | Supported | `lean-toolchain`, `lake-manifest.json`, `uv.lock`, Apache-2.0 LICENSE, durable commands. The baseline is attributed below; this repo does not claim to originate pivoted Cholesky/P-greedy convergence. |
| R38 | Required documents and claim separation (§28,29) | Supported | FINAL_AUDIT independently maps all original sections. Current summaries foreground the exact finite rate and reopened adaptive experiment, distinguish numerical residual models from exact crosses, correct stale deferred-method text and mark older proof-route claims historical. Original SPEC mathematics is unchanged. |
| R39 | Honest implementation-run provenance (agent rules) | Supported | Canonical FINAL_HANDOFF records exposed model family/mode, unavailable exact model/effort/categories/billing, cumulative goal counters, and explicitly historical API-equivalent cost assumptions; other documents link to it. |
| R40 | Phase merge/push verification (active objective) | Supported | DO's final `f0d082d` passes the 569-test root and push/PR CI `39961732232`/`39961740682`; PR #121 merged as `b041c723add1b9a7a0395cfffbdf4b793c0cfed9` on 2026-10-23. Merged/pushed-main closure checks repeat before handoff. Earlier phase evidence remains historical; Git cleanliness/CI do not close R25. |

## Original section and outcome coverage

| SPEC sections | Interpretation / ledger mapping |
| --- | --- |
| 1 | R01, R16, R17: main rate and application error, with honest parameter range. |
| 2 | Context and prior-work attribution, R37; not a new proof of every published convergence theorem. |
| 3 | Track A's staged structural/fill-distance baseline is proved in R06–R09. Track B's exact-rate outcome is R16/R17. |
| 4.1–4.3 | R10–R14: identities, weights, exponential signs and actual-prefix structure. |
| 5 | R15: D1's real quantifiers and explicit constant. |
| 6 | Candidate construction, not a mandate to implement both Taylor and Chebyshev. The allowed Taylor route is used. |
| 7 | R02 and numerical modes R20/R21: finite successful updates, complete maximum existence, zero stop. |
| 8 | R01, R03–R05: interpolation, cross, core determinant and nonsingularity. |
| 9 | R06/R07 preserve finite interfaces; DJ/DK prove arbitrary-domain PSD geometry and strict-SPD actual complete-run reduction (R08). |
| 10 | R09: original-kernel Lipschitz fill control proved by DJ. Published full asymptotic rates are not claimed by this staged baseline. |
| 11 | R16/R17: intended large-cutoff rate proved, literal arbitrarily-small-cutoff reading refuted. |
| 12 | A proposed approximation-to-greedy route. The successful local-minor/panel mechanism fulfills the structural objective; a separate general best-rank comparison is not claimed. |
| 13 | Research questions, not six independently promised positive answers. Actual all-orders signs/gauge are proved; pivot-location hypotheses and counterexamples remain explicitly distinguished. |
| 14 | **R25 incomplete:** all-cutoff residual-adaptive increasingly dense experiments are absent. R26/R27 diagnostics/normalized plots are supported for the fixed-grid artifacts only. The bounded interval-adaptive API smoke is not the requested census. |
| 15 | R18/R22: exact small certificates and rational/high-precision evidence. |
| 16 | R19/R21: cover theorem and adaptive certified maximizer with failure metadata. |
| 17 | R20 plus the conditional-extension note below; no arbitrary-η near-optimal claim. |
| 18 | R01 and DH's actual Green-integral difference, including the physical restricted measure. |
| 19 | R29/R30, computational sparse outcome, not a minimum-support theorem. |
| 20 | R30 delivered families; DN verifies R31's six comparisons with repeated artifacts and independent errors; PR #120 is merged. Formal AAA remains excluded as explicitly stated. |
| 21 | Suggested roles are implemented. DO moves generic completeness definitions below the fermionic layer so the SPD baseline no longer imports exponential total positivity for them; identifiers/statements are preserved. |
| 22 | See named-interface mapping below; aliases are not substituted for missing mathematics. |
| 23 | R32 plus the actual documented result spelling. |
| 24 | R22/R24/R33: precision and ties are exposed; DL adds explicit actual/requested precision and full metadata outside GECP. Cholesky/sparse remain float64 only, honestly labeled. |
| 25 | R34/R35: source gates and numerical tests, including tested failure paths. |
| 26 | R26–R28 support deterministic fixed-grid data, repeated canonical runs and bound provenance. R25's required adaptive scope remains open. |
| 27 A–F | A: R03–R05; B: R06–R09; C: R10–R14; D: R15; E: R16/R17; F: R29. |
| 28 | R17/R37/R38: no unqualified “solved Problem 4.2”, no “optimal” theorem naming, no literature-priority claim. |
| 29 | R38 and this explicit gap ledger. |
| 30 | Strong exact-rate layer is now available, but the full original scope is not replaced by the easier minimum-publishable threshold. |

## Named interfaces and artifacts

The requested Lean algebraic names are exported and audited in their actual
namespaces: `GECP.Run.gecp_interpolates_selected_rows`,
`gecp_interpolates_selected_cols`, `gecp_core_det_eq_prod_pivots`;
`PositiveDefinite.gecp_eq_pivotedCholesky_of_posDef`,
`posDef_gecp_residual_sup_eq_diag_sup` (finite canonical scope);
`Fermionic.fermionicKernel_reflection`, `fermionicKernel_centered`,
`fermionicKernel_continuous`, `fermionicKernel_separatedApprox`; and
`greenError_le_kernelError_mul_l1`. The stronger actual physical consequence
is `fermionicKernel_finiteRun_greenError_le_exp_interval`.
The proposed `fermionicKernel_gecp_error_le_exp` role is delivered by the more
explicit `fermionicKernel_finiteRun_error_le_exp`, with its cutoff range visible
in the signature; the literal all-positive-cutoff version is not aliased as true.

The seven proposed Python names are `FermionicKernel`, `gecp`,
`pivoted_cholesky`, `cross_approximation`, `evaluate_residual`,
`certified_pivot`, `dlr_rank_bound`. They are present in `kernelgecp.__all__`.
The example `pivoted_cholesky(kernel, domain, ...)` is realized by its documented
sampled matrix interface, not a second unimplemented callable signature.

Required maintained files exist: AGENTS, README, SPEC, FINDINGS, RESEARCH,
APPLICATION, MATHLIB, FINAL_HANDOFF, the Lean public re-export/audit files,
Python package/tests, pinned toolchain/locks, CI, census TOML/JSONL/PNG, exact
surrogate JSON, and synthetic application TOML/JSON/PNG. The canonical census
includes historical 21-case evidence at revision `efe9b0041cf5c05b7ce665fcebc4b41450fd40ce`,
SHA-256 `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
DM's 28-case schema-2 census is generated from clean
`c0e8c5093134f369b93781c413750b7f1de46036`, SHA-256
`c6854330fcb74007f8f9814fba48c314d34fe26c1cef5461d908ea5182b66bda`;
all 21 shared legacy pivot/rank traces remain identical. Denser data uses clean
`ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61`, SHA-256
`7251cbb528179cc5cbfa6eb5635b8fa280e087d1e04db1348b56a9065026d588`.
The independent report binds both datasets, its checker/grid sources and lock;
all 30 selected-core comparisons pass fixed thresholds. Matching resolution
ranks do not imply matching pivot paths: common selected prefixes are `3,0,0`.
The synthetic suite records four atoms versus universal rank 12 and its
documented historical revision/hash; no artifact is relabeled as a fresh run.

Durable commands remain `uv sync --all-extras`, `lake exe cache get`, and
`./scripts/verify.sh`. Focused checks do not replace the root command.
DI's source/wheel builds are local audit artifacts, not a published release.

## Conditional research extensions and attribution

Section 17 asks for an arbitrary-η structural-rate extension “where feasible”.
The η certificate/definition bridge is implemented; a near-optimal rate for
arbitrary inexact pivots is not. Exact local witness transport uses column
maximality, which an arbitrary η-complete pivot need not satisfy. It would be
incorrect to apply R16 to certified η<1 pivots by renaming them exact.
This is an explicit research extension, not a claimed practical convergence
certificate. If a final claim includes the practical approximate-pivot rate,
this extension becomes a required unresolved theorem rather than a solved item.

The baseline convergence literature is prior work. In particular,
[Sungwoo Jeong and Alex Townsend, historical v1, *Convergence of Pivoted Cholesky Algorithm
for Lipschitz Kernels*](https://arxiv.org/abs/2509.13582v1) proves original-kernel
Lipschitz fill-distance control and corresponding complete-pivot convergence
results. DJ formalizes original-kernel fill-distance control, not that paper's
full asymptotic rates. P-greedy/pivoted-Cholesky theory is
the control baseline, not this repository's claimed invention.

The [current v3](https://arxiv.org/abs/2509.13582v3) adds Christopher Wang.
The [workshop v3 update](https://arxiv.org/html/2602.05394v3) also records an
August 2026 claimed solution of Problem 4.2 by Venkata Siddharth Pendyala,
*Local-to-Global Convergence of Greedy Cross Approximation for Totally Positive
Kernels* ([DOI](https://doi.org/10.5281/zenodo.21863274)). The workshop update
was inspected; the underlying paper was unavailable to this audit, so its
claim is attributed rather than independently verified. No priority for the
repository's theorem is asserted.

## Next exact targets, in priority order

DJ/DK close the staged continuous baseline, not all published asymptotic
P-greedy theory. DN closes the six-method computational comparison. Required
work after DO's repair/verification phase is:

1. **Residual-adaptive seven-cutoff census (R25).** Implement deterministic
   nested grids whose refinement choices respond to measured residuals, rerun
   high-precision GECP at every cutoff through `1e-12`, retain stages, budgets,
   geometry, precision/tie/conditioning metadata and all failures, and repeat
   the canonical output. Independently check cores/held-out errors and render
   adaptive normalized curves. Sampled maxima are not continuous bounds.
2. **Final closure.** Reaudit the actual adaptive outputs and original scope;
   repeat root/package/install/exact-head CI, merge/push main and provenance
   checks before deciding whether the whole goal is achieved. An unsuccessful
   adaptive budget must remain visible and cannot be closed by a status label.

No missing requirement is resolved by a status annotation alone.

## Run provenance

Exposed model/mode, cumulative goal counters and unavailable billing/version
fields are recorded in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They are
operational provenance, not evidence for a theorem or a completed requirement.
