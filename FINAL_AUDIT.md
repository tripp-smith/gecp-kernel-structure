# Independent final audit against the original SPEC

Audit opened: 2026-10-22; DP numerical delivery updated: 2026-10-23.
Phase DO starts from synchronized main
`5355c72f2508eeb604b5515d789f48d70b42278d`. The original text in SPEC sections
1–30, rather than the subsequent phase annotations or the forty labels in
`COMPLETION_AUDIT.md`, defines the required scope.

**Decision: NOT COMPLETE.** The principal exact fermionic GECP rate is
supported, including its finite-run selected-core cross and actual Green
integral consequences. DP now supplies the previously missing seven-cutoff
residual-adaptive census: all seven cases and twenty-eight tolerance records
report protocol convergence, two full generations repeat bytes, and the
adaptive plot is rendered, inspected and byte-repeated. Independent MP192
checking passes all twenty-one representative selected-core checks and
twenty-eight full final-cross checks, with every sampled maximum within its
requested tolerance. The subsequent full original-scope closure audit
follows DP delivery. DP is **verified and merged** in
[PR #122](https://github.com/tripp-smith/gecp-kernel-structure/pull/122) at
2026-10-23 20:20:48 UTC as `3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`.
Final branch `3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes its exact local
root and push/PR CI. Exact pushed-main closure metadata/gates remain required;
the preliminary mapping of forty supported groups is not
whole-goal completion. A separate rough overall effort estimate is **99%**,
not a checklist-derived percentage or a completion decision.
DO's verified repairs and merged PR #121 remain historical evidence, not
verification of the new DP experiment.

This audit inspected original requirements, Lean definitions and proof
dependencies, Python source and tests, experiment configurations, and committed
JSON/JSONL evidence. Independent baseline and numerical reviewers supplied
additional source-level checks and counterexamples. The parent run supplies
the fresh root, package, installation, CI and Git evidence recorded below;
this document does not imply those commands were independently repeated by
the source-audit reviewer. Earlier phase verification remains historical.

## Historical DO verification and delivery

The parent run reports `./scripts/verify.sh` passing at code head
`f0d082dfbbc91e4213dea0d2cad06a36c4a78826`: **569 tests in 142.82 seconds,
3,634 Lean jobs, 876 axiom reports, 83 formatted files, mypy on 19 sources,
and a clean dependency audit**. The generic complete-pivot import repair and
all numerical/documentation repair regressions are included in this root
gate; those code-head repairs are no longer pending implementation or root
verification.

The source archive and wheel were inspected and contain `FINAL_AUDIT.md`,
the new generic Lean module and new tests where appropriate, while excluding
local note copies and caches. All nineteen wheel Python files match the
source hashes. Installed-wheel smokes outside the checkout pass with both:

- pinned Python 3.12.12, NumPy 2.5.2, SciPy 1.18.0 and mpmath 1.4.1;
- the declared-minimum SciPy 1.15.0 environment with NumPy 2.4.6.

Each installed environment exercises all 45 exports, cancellation/underflow,
adaptive tie scope, extreme and mixed-scale arithmetic, duplicate-core
handling, and 48 six-method fits. Independent baseline QA also passes 112
focused tests, fifteen signed interval cases, eight mixed-scale cases and
400 coordinate-reuse checks. These complement the full root gate rather than
replacing it.

Final push CI `39961732232` and PR CI `39961740682` pass all four checks.
The companion interval diagnostic preserves binary-float coordinates as
mpmath values and the same frozen-residual scope; its non-binary cutoff `1.3`
callback-capture regression is included in the root and installed checks.
Both float64 and MP grid paths exclude old/zero coordinates from near-tie
selection and preserve known interpolation zeros. Historical implementation
`e094130` passed 567 tests before those final consistency additions; that
earlier evidence is not relabeled as a later run. The archives and both
installed environments were rebuilt/repeated at `f0d082d`.

PR #121 merged at 2026-10-23 13:34:53 UTC as
`b041c723add1b9a7a0395cfffbdf4b793c0cfed9`. The final pushed-main metadata
head repeats root/CI and archive checks before handoff; exact Git/CI evidence
is checked against that head, not assumed from the branch result. No package
is published. At DO closure the required adaptive experiment remained absent;
the DP delivery below is a later, separate source revision and artifact.

## DP adaptive numerical delivery and independent verification

Direct inspection of `experiments/data/adaptive_census.json` finds seven
actual residual-driven nested-grid trajectories and twenty-eight tolerance
records, all reporting protocol convergence and final-sample passage. The
generation revision is clean `3cee4186fddf990abdd33dbe24c46ff908b1923f`;
configuration SHA-256 is
`ac47ec5146f5dc50916513917518d2bbfcb31da622b1f659225568cc48d641a1`.
The 4,617,990-byte canonical JSON has SHA-256
`730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c`.
The parent run reports both full generations exit successfully and compare
byte-for-byte. These data are additional to the preserved historical fixed-grid
artifacts; their generation provenance is not relabeled as current source.

Elimination uses MP128 arithmetic on explicit 53-bit repr-decimal coordinates.
MP192 factored-cross point checks drive odd-eighth residual/surplus probes;
final validation uses all six seventh-fraction points per final leaf panel
after every refinement and tolerance-prefix choice. Canonical controls are a
nine-stage budget, 512 nodes per axis, rank budget 256, four highest marked
panels before reflection mates, and generation margin `0.01`. Successful
termination requires real refinement and two consecutive passing strictest
stages, not merely a small terminal training residual or passing final sample.

| Cutoff | Actual stages | Final training shape | Chosen ranks at `1e-6,1e-8,1e-10,1e-12` |
| --- | --- | --- | --- |
| `1` | 3 | `25×25` | `6,7,8,9` |
| `10` | 2 | `41×57` | `12,13,15,17` |
| `100` | 3 | `109×89` | `24,28,33,36` |
| `1000` | 4 | `122×121` | `38,48,55,64` |
| `10000` | 4 | `170×169` | `53,66,80,97` |
| `100000` | 5 | `238×217` | `68,87,103,124` |
| `1000000` | 6 | `258×257` | `83,103,126,151` |

Every selected prefix retains the requested pivot/magnitude, sampled residual,
original-core determinant/minimum singular-value and physical-geometry
histories. First training-passing ranks and chosen jointly passing tested
ranks are distinct fields; no global minimum-rank claim follows. Full stage
axes, panels, scores, reflected marks, budgets, precision, conditioning, ties
and source/lock/environment provenance remain in the artifact.

The parent run reports the 2240-by-1440 `adaptive_census.png` rendered,
visually inspected and byte-repeated, SHA-256
`dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6`.
It shows seven normalized `k/log(1+Λ)` training curves, twenty-eight untouched
validation markers, training/chosen rank comparisons, stage errors and node
growth. One zero sample is omitted from logarithmic axes without a positive
floor; no unsuccessful cases are removed or present in this canonical run.
No universal fitted law or continuous supremum bound is inferred.

The independent MP192 checker completes successfully and its report confirms
full seven-case/twenty-eight-record coverage, every case's protocol convergence
and the recorded source revision blobs. All twenty-one representative original
selected-core determinant/pivot-product/SVD checks and twenty-eight full
final-cross recomputations pass; every recomputed final sample maximum meets
its requested tolerance. The 37,811-byte `adaptive_census_checks.json` has
SHA-256 `4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03`.
It binds canonical input SHA-256
`730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c`
and checker source SHA-256
`711ac118e26cdab256ed43a3918f876fc7d00e3b0245885dfe994fb761af1316`.
These are independent finite-sample checks, not continuous-error certificates.
The parent run reports DP artifact/code head
`a6488cc0c91400bab3f327d93db153c390ef20f4` passing the exact local
`./scripts/verify.sh`:
**713 tests in 150.40 seconds, 3,634 Lean jobs, all 876 allowed axiom reports,
95 formatted files, mypy on 21 sources and a clean dependency audit**.
The final source/wheel inspector finds clean archives, all three adaptive
artifacts and new source/tests included, and local notes/caches/scratch
excluded. All 21 wheel Python files and the packaged APPLICATION/README
match source bytes.

Outside-checkout installed-wheel checks pass with pinned Python 3.12.12,
NumPy 2.5.2, SciPy 1.18 and mpmath 1.4.1, and with the declared-minimum
NumPy 2.4.6/SciPy 1.15.0 environment. Each exercises all 45 exports, adaptive
success/budget-failure paths and 48 six-method fits. These are current DP
checks, separate from the preserved nineteen-source DO archive evidence.
Both rebuilt wheels repeat these installed checks at `a6488cc`. The earlier
development-tree root passed 713 tests in 151.22 seconds; it is retained as
historical evidence rather than relabeled as the exact-head run.

Exact-head push CI `40017052450` and PR CI `40017056920` pass all four checks
at `a6488cc`; that artifact/code-head evidence remains historical.

Final branch `3c38183e5c2a7d2859867bb7acff1353ee98fb76` repeats the full local
root with 713 tests in 150.98 seconds. Exact-head push CI `40017937659` and
PR CI `40017945103` pass all four checks, and final source/wheel archives and
both outside-checkout installed-wheel smokes are rebuilt and pass again.
PR #122 merges at 2026-10-23 20:20:48 UTC as
`3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`; local main fast-forwards to
origin and generation revision `3cee4186fddf990abdd33dbe24c46ff908b1923f` is
confirmed an ancestor. DP is verified and merged on this evidence. Exact
pushed-main closure metadata/gates repeat before handoff and remain pending;
no pushed-main CI identifier is asserted yet. The next decision is
the full original-SPEC closure audit after DP, not automatic completion from
this numerical delivery or the inventory count.

## The actual headline mathematical result

`Fermionic/FiniteRunRate.lean:166–176` proves, for every real `Λ ≥ 1`, every
finite actual `GECP.Run fermionicKernel` complete on the physical rectangle,
and every `(t,ω) ∈ [0,1] × [-Λ,Λ]`,

\[
 |R_k(t,\omega)| \le 1024\exp\!\left(
 -\frac{(\log 2)^2}{256}\frac{k}{\log(1+\Lambda)}\right),
 \qquad k=\texttt{run.rows.length}.
\]

Lean paths in this document are relative to `GECPKernelStructure/` unless
prefixed otherwise. `GECP/Definitions.lean:12–21` defines a run by exact
successive Schur updates with nonzero stored pivots. `CompletePivotOn` means
domain membership and an absolute global maximum; recursive `Run.CompleteOn`
requires this at every stored step. DO moves these generic predicates from
`Fermionic/CompleteRun.lean` into `GECP/CompletePivot.lean`, preserving their
interfaces (`CompletePivotOn`:11, `Run.CompleteOn`:19) and removing the SPD
baseline's dependency on the fermionic theorem layer. The unchanged predicate
interfaces and the import repair pass the fresh root build/axiom gate.

The rate assumes no convergence estimate, residual determinant envelope,
independent realization certificate, or successful pivot after the stored
run. Its substantive chain is:

- `Fermionic/ExponentialTotalPositivity.lean:599–634` proves all-orders signs
  for the actual exponential and fermionic kernels.
- `Fermionic/CommonResidualGauge.lean:19–33` derives the positive-minor
  frequency-reflected presentation. `GECP/ResidualTotalPositivity.lean:177–222,
  256–260` derives one magnitude-preserving sign gauge at every minor order
  after arbitrary successful prefixes.
- `GECP/LocalMinorWitness.lean:85–181` crosses an unmarked step backwards
  using column maximality and the singleton Fischer inequality, retaining
  the local column predicate. `GECP/FiniteLocalMinorWitness.lean:174–262`
  constructs finite marked-minor witnesses and the resulting pivot count.
- `Fermionic/LocalPanelBudget.lean:183–193` proves the initial local determinant
  budget for original-kernel samples, both frequency signs, physical time
  endpoints, and unordered or repeated tuples. `FiniteRunRate.lean:71–126`
  combines it with the finite panel cover to give error `≤2^(-q)` after
  `32*(s+1)*(33+4*q)` successful complete pivots whenever `Λ≤2^s`.
- `Fermionic/RealCutoffRate.lean:146–246` turns that explicit natural bound
  into an exponential envelope and a logarithmic real-cutoff cover, without
  enlarging the domain on which pivots must be complete.

`FiniteRunRate.lean:178–188` separately proves the `log(2+Λ)` version for all
positive cutoffs. Its `207–217` gives the explicit sufficient rank

\[
 \left\lceil\frac{256}{(\log2)^2}\log(2+\Lambda)
                   \log(1024/\varepsilon)\right\rceil_+
\]

for every positive real tolerance. Literal all-positive-cutoff G1 with
`log(1+Λ)` is false: `RealCutoffRate.lean:340–392` constructs, for any proposed
positive `C,c`, `Λ=min(1,c²/(16C))` and an actual first complete pivot whose
rank-one residual contradicts that bound. This does not refute the `Λ≥1`
result and does not license rewriting the original conjecture.

`GECP/RunCrossApproximation.lean:19–29` uses the original selected kernel rows,
columns and proved-invertible core, and identifies the iterative residual
with the actual selected-core cross error. `FiniteRunRate.lean:220–261`
transfers the rate to that cross and the actual difference of two weighted
Green integrals for every real-valued integrable density. The restricted
measure form on `[-Λ,Λ]` discharges support. Integrability of both integrands
is derived in `RunCrossApproximation.lean:68–99`. A zero complete maximum
implies zero residual on the domain (`RunCrossApproximation.lean:7–12`), and
the stopped-or-budget result at `FiniteRunRate.lean:192–205` requires no
zero-pivot update. Compact attainment for every finite fermionic residual is
proved in `Fermionic/CompletePivotExistence.lean:10–22`.

For D1, `Fermionic/RealSeparatedApprox.lean:8–13,62–70` constructs an actual
finite separated list for all real `Λ≥1`, `0<ε<1`, with
`r≤128*log(1+Λ)*log(128/ε)` and uniform error `≤ε`. The factors are the
dyadic Taylor/mask/reflection construction in `ExponentialFamily.lean` and
`SeparatedApprox.lean`; an approximation hypothesis is not supplied by the
caller. These headline results are public imports in `Theorems.lean:1–10`
and explicit audit targets in `scripts/AxiomAudit.lean:43–98,105–132`.

## Original thirty-section checklist

“Supported” is limited to the scope in the evidence column. “Pending” means a
current repair or final gate has not been verified. “Incomplete” means a
required artifact or outcome is absent. Context, proposed routes and examples
are not silently promoted into mandatory additional theorems.

| SPEC section | Original requirement and independent finding | Evidence and scope |
| --- | --- | --- |
| 1 | Supported exact mathematical objective; full project incomplete | Finite exact rate, original cross and Green consequence above close the large-cutoff linear-to-logarithmic gap. The small-cutoff boundary is explicit. |
| 2 | Supported baseline role and attribution | SPD/Cholesky and fill-distance results below reproduce the control setting. Updated references distinguish the historical two-author Cholesky v1 from v3 and attribute the workshop's separate solution claim without verifying its inaccessible underlying paper or claiming priority. |
| 3 | Supported two mathematical tracks | Track A: `PositiveDefinite/ContinuousCholesky.lean:68–148`; Track B: actual finite-run fermionic rate, rather than a structural hypothesis merely assumed to contain it. |
| 4 | Supported identities and structural infrastructure | `Fermionic/Symmetry.lean:7`, `Centering.lean:12`, `Derivatives.lean:9,40,65,113,125`, all-orders signs and actual-prefix gauge above. No universal pivot-localization or residual-section-monotonicity claim. |
| 5 | Supported D1 at its original real quantifiers | `RealSeparatedApprox.lean:8–13`: real `Λ≥1`, `0<ε<1`, explicit list/rank/error. |
| 6 | Supported allowed construction route | `ExponentialFamily.lean:131–170,191–260`, `SeparatedApprox.lean:54–133`: Taylor, dyadic masks and reflection. The candidate section does not require both Taylor and Chebyshev. |
| 7 | Supported formal model and root-verified numerical repairs | Actual nonzero recursive updates, complete finite runs, compact maxima and exact stop. Finite-grid and adaptive Python modes are distinct in `kernelgecp/gecp.py`; certificate and finite-arithmetic repairs pass the fresh root gate. |
| 8 | Supported algebraic invariants | `GECP/Interpolation.lean:13,25`, `Determinant.lean:53,134,144`, `RunCrossApproximation.lean:25`; original selected core, no supplied factorization conclusion. |
| 9 | Supported continuous strict-SPD baseline | `ContinuousCholesky.lean:68,111,119,125,130`: actual complete runs must be diagonal, equivalent to independently defined diagonal maximization, retain PSD and diagonal domination. Weak PSD may have off-diagonal ties. |
| 10 | Supported staged original-kernel fill bound | `PositiveDefinite/KernelFillDistance.lean:11,29` and `ContinuousCholesky.lean:141`: derive `abs R(x,y) ≤ 2Lh` from original coordinate Lipschitz regularity and a selected-node cover. Published full `k`-dependent rates are not claimed or required for this staged target. |
| 11 | Supported rate and rigorous parameter obstruction | `FiniteRunRate.lean:166,178`; `RealCutoffRate.lean:340`. Do not call the literal all-positive-cutoff conjecture true. |
| 12 | Proposed route; structural objective supported | The question about a general `C(k)σ_(αk)` comparison is not a separately mandatory theorem. Actual-prefix gauge/local-minor/panel counting proves the desired original-kernel rate without that comparison. |
| 13 | Research questions investigated, not six promised positive answers | All-orders signs, gauge and backward witnesses are proved; numerous exact obstruction modules distinguish false stronger sign/localization hypotheses. No universal nested one-dimensional pivot rule is claimed. |
| 14 | Supported verified/merged DP adaptive delivery and independent checks | `adaptive_census.json` supplies seven actual residual-driven nested-grid trajectories and twenty-eight converged sampled tolerance records from clean `3cee4186fddf990abdd33dbe24c46ff908b1923f`; both full generations repeat bytes. Schema-3 records retain all requested diagnostics; the adaptive normalized/stage plot is rendered, inspected and repeated. MP192 structure/provenance, twenty-one representative original-core and twenty-eight full final-cross checks pass. PR #122 is verified and merged; pushed-main closure and the fresh original-scope audit remain pending. Historical fixed-grid artifacts remain separate. |
| 15 | Supported exact small-instance layer | `SmallInstanceChecks.lean:38–107` certifies a concrete complete two-pivot real/rational grid, minors and `[2,1/2]` pivots. Larger rational exponential cases are exact Python evidence in `tests/test_surrogates.py` and `experiments/data/exact_surrogates.json`, not larger Lean enumerations. |
| 16 | Supported formal cover bound and root-verified numerical repairs | `Discretization/Grid.lean:10` derives `gridMax+L*radius` from explicit cover/Lipschitz assumptions. `kernelgecp/certified.py` has tested adaptive interval and conditional analytic-Lipschitz APIs; cancellation, underflow, rounded-midpoint, invalid-enclosure and budget-failure paths pass the fresh root gate. |
| 17 | Supported η-certificate bridge; structural extension unresolved | `Discretization/LipschitzCertificate.lean:8` transports an upper bound to the η inequality. Original wording says strengthen rates “where feasible”; no arbitrary-η near-optimal rate is proved. Exact rates cannot be applied to certified η<1 pivots. |
| 18 | Supported actual Green consequence | `FiniteRunRate.lean:231–261`, `RunCrossApproximation.lean:68–99`; actual difference of integrals, real `L¹` density, physical support and every physical time point. |
| 19 | Supported computational sparse application | `kernelgecp/sparse.py`, `applications.py`, `experiments/data/synthetic_applications.json`: clean four-atom recovery with held-out error about `3.33e-16` versus universal rank 12 at `1e-8`; dictionary/noise/validation scopes remain explicit. No minimum-support or uniqueness theorem. |
| 20 | Supported measured six-method comparison | `kernelgecp/comparison.py:469–558` passes common observations to GECP, AAA-derived, vector-fitting, Prony, nonlinear least-squares and greedy pursuit methods. `spectral_comparison.json` retains 192 fits and failures; `tests/test_spectral_artifacts.py` independently reconstructs errors. Formalizing AAA was explicitly excluded initially. |
| 21 | Supported architecture roles and verified import cleanup | Exact core, SPD, fermionic, discretization, Green and public audit modules exist. Layout was suggested, not a filename contract. DO extracts unchanged generic complete-pivot definitions from the fermionic layer; fresh Lean/axiom checks pass. |
| 22 | Supported semantic public theorem interfaces | Public algebraic, SPD, fermionic and Green identifiers appear in `Theorems.lean`/`scripts/AxiomAudit.lean`. The conservative `fermionicKernel_finiteRun_error_le_exp` delivers the proposed rate role with its range visible; no “optimal” theorem or false all-cutoff alias. |
| 23 | Supported scoped API adaptations | Seven proposed names are exported in `kernelgecp/__init__.py:53–99`. `GECPConfig` supplies controls; result residual access is `residual_matrix`, histories and `evaluate_residual`; Cholesky consumes sampled matrices. These document the proposed example's roles, not literal compatibility with every example spelling/signature. |
| 24 | Supported precision paths and root-verified tie/arithmetic repairs | Float64, mpmath with canonical decimal histories and `Fraction` surrogates exist. DO records observed candidate ties, scopes interval certificates, and handles representable extreme/mixed-scale updates and unrepresentable failures. MP128 census elimination still uses 53-bit grid coordinates; LAPACK conditioning diagnostics remain estimates. |
| 25 | DP exact-head root, package/install and CI verified; PR merged | Final branch `3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes the full local root with 713 tests in 150.98 seconds; exact-head push CI `40017937659` and PR CI `40017945103` pass all four checks. Clean archives and both rebuilt installed-wheel environments pass again; PR #122 is merged. The established 3,634 Lean jobs, 876 allowed axiom reports, 95-file formatting, mypy on 21 sources and clean dependency gates remain documented above. Pushed-main closure is pending; earlier artifact/code-head and DO checks remain historical evidence. |
| 26 | Supported adaptive and historical datasets | Direct inspection finds seven schema-3 adaptive cases and twenty-eight converged tolerance records through `1e-12`, 128-bit elimination, 53-bit explicit coordinates and MP192 probe/final sample metadata. Revision/config/source/lock/environment provenance and exact full-run byte repeats are recorded above; the passing independent report pins its input and checker source and verifies revision blobs. Historical canonical JSONL retains its original provenance; final DP shipping gates remain separate. |
| 27 | A–E formal outcomes supported; F computational outcome supported | A: interpolation/core product; B: actual strict-SPD equivalence/fill bound; C: identities/sign infrastructure; D: real D1; E: actual exact logarithmic rate; F: sparse application with explicit held-out error. DP supplies the independently checked §14 artifact, but full completion still requires final delivery and the subsequent original-scope audit. |
| 28 | Supported claim boundaries and updated attribution | Current summaries preserve small-cutoff obstruction, exact-vs-inexact and proof-vs-observation distinctions. No unqualified all-parameter G1, practical inexact rate, sparse optimality, or novel-priority claim. |
| 29 | Required documents maintained; pushed-main closure synchronization pending | README, SPEC, FINDINGS, RESEARCH, APPLICATION, MATHLIB and canonical FINAL_HANDOFF preserve the theorem and historical verification. Current audits record the adaptive delivery, passing independent checks and verified PR #122 merge. Original mathematical goals are unchanged; exact pushed-main closure metadata/gates must be synchronized before handoff. |
| 30 | Minimum mathematical outcomes supported; full task incomplete | Reusable core, SPD equivalence, fermionic algebra, constructive low rank and an actual exact rate exceed the minimum formal threshold. DP now adds the independently checked seven-cutoff adaptive census across six orders of magnitude; whole-goal closure still requires remaining delivery gates and a fresh original-scope audit. |

## Independent disposition of the forty inventory groups

These identifiers retain compatibility with the older inventory, but their
statuses follow inspected scope rather than its previous “40/40” count.
`S` means supported at the stated scope; `P` means repair or final verification
pending; `I` means a required deliverable is missing. No percentage of overall
completion is inferred from this table. DP's preliminary mapping is now
**40 supported groups and 0 incomplete groups**, using the stated evidence
scopes and preserving pending final gates. This is an input to the
subsequent original-scope closure audit, not that audit's completion decision.

| ID | Disposition | Evidence / boundary |
| --- | --- | --- |
| R01 | S | Original selected-core cross identity: `GECP/RunCrossApproximation.lean:19–29`; exact and two-precision regressions in `tests/test_finite_run_closure.py`. |
| R02 | S | Exact runs, compact fermionic maxima, finite prefixes and zero termination, as audited above. |
| R03 | S | Selected-row interpolation: `GECP/Interpolation.lean:13–22`. |
| R04 | S | Selected-column interpolation: `GECP/Interpolation.lean:25–34`. |
| R05 | S | Original-core pivot product and nonsingularity: `GECP/Determinant.lean:134–151`. |
| R06 | S | PSD Schur preservation/diagonal domination: `PositiveDefinite/ResidualPSD.lean`; arbitrary-domain `KernelGeometry.lean`. |
| R07 | S | Canonical finite diagonal-preferring trace equivalence: `PositiveDefinite/PivotedCholesky.lean`, `tests/test_cholesky.py`. |
| R08 | S | Original strict-SPD actual complete-run reduction/equivalence: `ContinuousCholesky.lean:68–139`; fresh strictness derived in `StrictKernel.lean`. |
| R09 | S | Original-kernel Lipschitz `2Lh` bound: `KernelFillDistance.lean:11–40`, `ContinuousCholesky.lean:141–148`. |
| R10 | S | Exact reflection: `Fermionic/Symmetry.lean:7–15`. |
| R11 | S | Exact centered form/positive weight: `Centering.lean:12–27`, `Kernel.lean:10–19`; independent high-precision centered-formula regression passes under R34. |
| R12 | S | Actual joint continuity, coordinate derivatives and physical bounds: `Derivatives.lean:9–132`. |
| R13 | S | Actual all-orders exponential/fermionic signs: `ExponentialTotalPositivity.lean:599–634`; signature and positive scaling in `SignRegularity.lean`. |
| R14 | S | Common residual sign gauge and finite backwards marked-minor witnesses; no conditioned-feature sign conjecture assumed. |
| R15 | S | Real-quantifier D1 and explicit list witness: `RealSeparatedApprox.lean:8–70`. |
| R16 | S | Actual finite exact complete-run logarithmic cutoff rate: `FiniteRunRate.lean:166–176`; no future-success/convergence premise. |
| R17 | S | Actual first-pivot counterexample resolves the literal small-cutoff boundary; all-positive `log(2+Λ)` variant remains separate. |
| R18 | S | Lean two-by-two complete-grid/minor/pivot certificate and exact Python exponential-surrogate checks; scopes distinguished. |
| R19 | S | Generic explicit Lipschitz/cover transfer: `Discretization/Grid.lean:10–20`. |
| R20 | S | Formal η arithmetic bridge and root-verified numerical lower/upper consistency repairs; no arbitrary-η structural rate. |
| R21 | S | Adaptive interval and conditional analytic-Lipschitz APIs, repaired singleton lower bounds, underflow/rounded-midpoint/invalid-enclosure guards, and success/budget-failure paths pass the fresh root gate. Certificates remain conditional on valid callbacks and their stated residual model. |
| R22 | S | Float64/mpmath/exact-rational paths and sub-float64 pivot regression; precision/coordinate scopes explicit. |
| R23 | S | Cross evaluation uses solves (`kernelgecp/gecp.py` and `applications.py`); formal inverse notation is algebraic. |
| R24 | S | Metadata, frozen-residual certificate scope, observed candidate ties and float64 arithmetic-range/duplicate-core failure paths pass the fresh root gate and installed-wheel smokes. Extreme-scale determinant/singular-value outputs remain numerical estimates. |
| R25 | S, preliminary DP mapping | Actual clean-`3cee4186` canonical adaptive JSON has all seven increasingly dense nested trajectories and twenty-eight converged tolerance records; both full generations repeat bytes. Stage counts `3,2,3,4,4,5,6` and the final largest `258×257` grid are recorded above. Independent MP192 structure/provenance, twenty-one representative original-core and twenty-eight full final-cross checks pass. DP is verified and merged in PR #122; pushed-main closure and the fresh original-scope audit remain pending. |
| R26 | S | DP schema-3 prefixes retain signed pivots/magnitudes, sampled residual/original-core determinant/minimum singular-value histories and all specified physical geometry. Complete stage axes/panels/marks, training/chosen ranks, precision/ties/conditioning/budget/failure scopes are explicit. All twenty-one independent representative core checks and twenty-eight final-cross recomputations pass; historical fixed-grid diagnostics remain intact. |
| R27 | S | The 2240-by-1440 adaptive normalized/stage PNG is rendered, inspected and byte-repeated, SHA-256 `dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6`. Seven curves and twenty-eight validation markers remain sampled evidence; one zero is omitted without a positive floor, and no universal law or continuous-rate evidence is inferred. Historical normalized plots are preserved. |
| R28 | S | Canonical adaptive JSON SHA-256 `730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c` repeats bytes and records exact clean generation revision/config/source/lock/environment provenance. Passing independent report SHA-256 `4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03` binds that input and checker source and verifies full coverage and revision blobs. Historical artifacts retain their original generation metadata. |
| R29 | S | Particular synthetic Green functions have fewer atoms than the universal basis with measured held-out errors; `synthetic_applications.json`, `tests/test_applications.py`. |
| R30 | S | Explicit spectral-frequency/weight interface, OMP and bounded nonlinear refinement in `kernelgecp/sparse.py`. |
| R31 | S | Shared-observation six-method source and 192 retained artifact fits; independent error reconstruction and historical source-manifest compatibility checks pass the fresh root gate. Each installed environment also passes 48 six-method fits; historical artifact generation is preserved. |
| R32 | S | Seven proposed public names exported; documented config/result/matrix-interface adaptations, not literal example signature compatibility. |
| R33 | S | Deterministic finite-grid ties, observed adaptive sample ties, no selected-coordinate reuse, explicit arithmetic-range/duplicate-core/diagnostic failures and exhausted certificate budgets pass root and installed-wheel checks. |
| R34 | S | Enumerated numerical categories plus independent centered-formula, cancellation, signed underflow, rounded-midpoint, invalid-gradient, mixed-scale, coordinate-reuse, duplicate-core and companion-diagnostic regressions pass the final 569-test root gate. |
| R35 | S | DP final branch `3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes the exact local root: 713 tests in 150.98 seconds, with the established Lean/axiom, formatting, mypy and dependency gates. Exact-head push CI `40017937659` and PR CI `40017945103` pass all four checks. Earlier artifact/code-head root and CI evidence remain recorded above; pushed-main closure remains under R40. |
| R36 | S | DP's final source/wheel archives are inspected clean, include all three adaptive artifacts and new source/tests, and exclude local notes/caches/scratch. All 21 wheel Python sources and packaged APPLICATION/README hash-match source. Outside-checkout installed-wheel smokes pass pinned Python 3.12.12/NumPy 2.5.2/SciPy 1.18/mpmath 1.4.1 and minimum NumPy 2.4.6/SciPy 1.15.0; each checks 45 exports, adaptive success/budget failure and 48 six-method fits. Historical nineteen-source DO evidence is preserved separately. |
| R37 | S | Pinned Lean/Python setup, locks and Apache-2.0 license exist; current references distinguish historical and updated Cholesky authorship, attribute the workshop's separate claim, and assert no literature priority. |
| R38 | S | Current audits record the actual finite theorem, completed adaptive generation/repeat/plot, passing independent checks and numerical scope boundaries; historical phase/route text is distinguished. Final DP documentation/gate synchronization remains under R40. |
| R39 | S | `FINAL_HANDOFF.md#implementation-run-metadata` remains canonical; exposed model/mode/time/usage fields and unavailable exact-version/token-category/billing fields are distinguished. Final closure counters/evidence are operational metadata, not completion proof. |
| R40 | S; DP verified/merged, pushed-main closure pending | DP final branch `3c38183e5c2a7d2859867bb7acff1353ee98fb76` passes its exact root, clean rebuilt archives/installed wheels and push/PR CI `40017937659`/`40017945103`. PR #122 merged at 2026-10-23 20:20:48 UTC as `3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`; local main fast-forwards to origin and generation revision `3cee4186` is confirmed an ancestor. Exact pushed-main closure metadata/gates repeat before handoff and remain pending. Artifact/code-head and DO gates are preserved as historical evidence; verified merge and the forty-group mapping do not establish whole-goal completion. |

## Reproduced defects, certificate scope and remaining work

The numerical reviewer reproduced a certificate reporting
`lower_bound=2.220446049250313e-16` greater than
`upper_bound≈1.300000756e-16` while `certified=True`, from 53-bit cancellation
in `(1+1.3e-16)-1`. An enclosure must use a conservative point lower bound
and maintain `0≤lower≤upper`; ordinary rounded point arithmetic cannot be
treated as an exact lower bound. Adaptive ties were not globally measured yet
were represented by compatibility counts of one. A finite rank-one float64
input of magnitude `1e308` also overflowed in an unscaled outer product.
The centered method delegated to the direct evaluator, and the suite lacked
an independent centered exponential/cosh identity check. These reproduced
defects now have code-head fixes and regressions included in the passing root
gate and installed-wheel checks. The former passing suite alone did not
establish that these cases worked.

The final float64 Schur update uses bounded binary mantissas (`frexp`) and
recombines exponents only at the final product (`ldexp`), avoiding both a
product-first overflow and normalization-first underflow on representable
mixed-scale terms (`kernelgecp/gecp.py`, `gecp_matrix`). Selected-row/column
zeros are explicitly preserved and old coordinates are masked from near-tie
selection in both arithmetic paths. Original selected cores with exact duplicate
rows or columns are rejected: roundoff cannot promote their
spurious fresh pivot into a usable core. These are exact duplicate checks on
the stored input, not a general exact-rank test. Genuinely unrepresentable
updates and failed/nonfinite diagnostics raise explicit breakdown errors.

`tests/test_final_numerical_audit.py` independently exercises
mixed-scale products, selected-coordinate reuse, tiny unrepresentable terms,
duplicate cores in both orientations, diagnostic failure and 35 independent
centered-expression cases. Additional independent QA covers the signed
interval, scale and coordinate-reuse cases listed above. At extreme dynamic
range, a LAPACK determinant or smallest-singular-value estimate can still be
zero for a nonsingular mixed-scale core. Such estimates are not proofs of
rank, conditioning or singularity and are not used as exact rejection
criteria. This limitation remains explicit rather than being hidden by the
arithmetic repairs.

Any interval certificate must identify its object: the interval evaluator
encloses a residual defined using stored, rounded pivot/coefficient constants.
It is not automatically a certificate of the exact Lean run or of a separately
evaluated selected-core solve. Caller-supplied analytic Lipschitz constants
are conditional assumptions, and budget exhaustion remains uncertified.
The exact finite-run theorem therefore does not certify numerical η<1 pivot
trajectories or prove convergence of the practical adaptive DLR construction.
Section 17's “where feasible” extension remains explicit research work.

To complete the original scope:

1. Complete exact pushed-main closure metadata/root/CI/provenance checks at the
   final main head before handoff. DP is verified and merged in PR #122;
   DO's verified
   repairs and prior DN evidence remain historical rather than substitutes
   for these new delivery gates. Preserve sampled scope and all failures; the
   passing independent numerical checks do not establish exact complete-pivot
   behavior or continuous supremum bounds.
2. Perform the subsequent full original-SPEC closure audit against all thirty
   original sections and the completed DP evidence. Reassess this **NOT
   COMPLETE** decision only after that audit and its gates; the preliminary
   forty-supported-group mapping cannot make that decision by itself.

No novel-priority conclusion follows from this repository audit. Mathematical
support is the checked content and scope of the results; literature priority,
full published P-greedy asymptotics, arbitrary-η near-optimal rates and sparse
optimality are not claimed. Runtime counters and any explicitly labeled
API-equivalent pricing estimates are operational provenance, not theorem or
completion evidence.

Canonical exposed model/mode/time/usage and unavailable billing/version fields:
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata).
