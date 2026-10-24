# Application guide

The public package supplies stable fermionic-kernel evaluation, deterministic
finite-grid GECP, solve-based cross approximation, diagonal-pivoted Cholesky,
conditional Lipschitz pivot certificates, numerical separated interpolation,
and sparse spectral recovery. The separate residual-adaptive experiment
controller now has seven-cutoff, 28-prefix sampled evidence and passing
independent MP192 checks. DQ's original-scope source audit and boundary repairs
are verified and merged in PR #123. Exact-main closure `e8c35ed` passed 757
tests in 150.31 seconds, all proof/type/style/security gates,
[CI 40051138962](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40051138962),
clean archive inspection and both outside-checkout installed environments.

```python
import numpy as np

from kernelgecp import FermionicKernel, GECPConfig, gecp, evaluate_residual

kernel = FermionicKernel(cutoff=100)
result = gecp(
    kernel,
    config=GECPConfig(pivot="grid", tol=1e-10, precision_bits=128),
)
residual = evaluate_residual(
    kernel,
    result,
    np.linspace(0, 1, 101),
    np.linspace(-100, 100, 401),
)
print(result.rank, result.converged, result.stop_reason, abs(residual).max())
```

`cross_approximation` solves with the selected core; it never forms an explicit
inverse. Results include residual history, determinants, smallest singular
values, near-tie counts, convergence, stop reason, and actual precision.
Runs above 53 bits retain canonical decimal strings in `result.high_precision`;
the NumPy compatibility views remain float64.
Determinants and smallest singular values are non-certifying numerical
estimates: unscaled LAPACK can report zero for a nonsingular mixed-scale core,
even when its true determinant is representable. A zero estimate alone is not
a singularity proof. Float64 updates use bounded binary mantissas/exponents,
preserve known interpolation zeros and exclude old coordinates from selection;
unrepresentable updates and exact duplicate rows/columns in the selected core
raise `NumericalBreakdownError`. Requests below arithmetic resolution do not
give a rigorous accuracy guarantee or exact numerical-rank determination.
For PSD matrices, `pivoted_cholesky` and `gecp_matrix` use the same
diagonal-preferring lexicographic tie rule; the recursive equivalence of those
pivot policies is separately proved in Lean.

## Certified pivots

`certified_pivot` performs deterministic float64 cell subdivision from a
supplied coordinatewise Lipschitz bound. Its result is conditional on that
analytic bound and the supplied point evaluations; their rounding is not
automatically enclosed. `interval_certified_pivot` instead uses `mpmath`
interval arithmetic, outward-rounded float endpoints, singleton-interval
candidate lower bounds, and an optional interval-gradient mean-value enclosure.
Its rigor is conditional on the supplied interval callbacks enclosing the
same function throughout the domain. A rounded point value alone is not a
certified lower bound.

`gecp(FermionicKernel(...), config=GECPConfig(pivot="adaptive", ...))` applies
the interval method to every continuous residual pivot and retains each
successful `PivotCertificate`. Cell-budget exhaustion returns
`stop_reason="uncertified_pivot"`; it never silently becomes convergence.
These certificates audit pivot quality for the residual model defined by the
stored rounded pivot coefficients. They do not enclose all coefficient error
relative to an ideal exact GECP run, nor the error of the separate float64 core
solve used by `evaluate_residual`. Adaptive candidate-tie counts describe only
the observed sample set, not the multiplicity of global continuous maximizers.
The Lean derivative theorems prove the underlying kernel formulas and domain
bounds, but Python execution is not itself a formal proof of a continuous GECP
rate. In particular, the exact-complete-pivot rate theorem is not automatically
applicable to a numerical or certified inexact run.

DO's certificate, tie, extreme/mixed-scale arithmetic, exact-duplicate-core and
centered-expression regressions pass the 569-test root gate and installed-wheel
checks outside the checkout, with both SciPy 1.18 and the declared minimum
1.15.0. [PR #121](https://github.com/tripp-smith/gecp-kernel-structure/pull/121)
is merged. DO's repairs alone did not supply the all-seven-cutoff adaptive
experiment reopened in [FINAL_AUDIT.md](FINAL_AUDIT.md). DP now supplies the
actual adaptive-grid generations, passing independent post-run and final
branch root/CI/package/install checks described below. PR #122 is merged;
DP main closure passed at `5a5fef0`, CI `40018905634`; DQ's later original-
scope closure passed at `e8c35ed`, CI `40051138962`.

## Separated interpolation

`fermionic_separated_approximation` uses piecewise barycentric interpolation on
dyadic Chebyshev frequency bands. `dlr_rank_bound` and
`exp_family_rank_bound` report explicit arithmetic counts for research
comparison.

`fermionic_dyadic_taylor_approximation` is the independent high-precision
implementation of the Lean-verified construction. For cutoff `Λ`, it chooses
`s = ceil(log₂ Λ)` and `p = ceil(log₂(1/tolerance))`; its public `rank` is
`16p(s+1)` and `error_bound` is `2⁻ᵖ`.

```python
from kernelgecp import fermionic_dyadic_taylor_approximation

approximation = fermionic_dyadic_taylor_approximation(1e6, 1e-8)
values = approximation.evaluate([0.0, 0.5, 1.0], [1e6, 0.0, -1e6])
print(approximation.rank, approximation.error_bound, values)
```

This constructive bound certifies low-rank existence. It does not certify a
GECP convergence rate.

## Sparse spectral recovery

`sparse_representation` performs orthogonal matching pursuit and optional
bounded nonlinear frequency refinement. It returns `converged=False` when the
requested tolerance is not met. The fixed two-delta regression reaches
`1e-8` with two atoms. The separate shared-observation comparison suite also
implements AAA-derived physical poles, iterative vector fitting and
Prony/ESPRIT; see the six-method comparison below. Failed or nonphysical fits
remain visible rather than being replaced by another method.

## Realistic synthetic application suite

`run_synthetic_application_suite` connects the public algorithms in three
workloads: continuous fermionic spectral compression, clean/noisy sparse
spectral recovery, and PSD covariance landmark selection.

```python
from kernelgecp import (
    SyntheticApplicationConfig,
    run_synthetic_application_suite,
)

suite = run_synthetic_application_suite(
    SyntheticApplicationConfig(),
    git_commit="local",
)

for record in suite.green_compression:
    print(
        record.fixture,
        record.gecp_rank,
        record.green_validation_error,
        record.transfer_bound,
    )

for record in suite.sparse_recovery:
    print(
        record.fixture,
        record.candidate_strategy,
        record.recovered_atom_count,
        record.held_out_error,
        record.stop_reason,
    )

print(
    suite.psd_landmarks.selected_rank,
    suite.psd_landmarks.pivot_paths_agree,
    suite.psd_landmarks.cross_max_error,
)
```

The two continuous densities are normalized on `[-8,8]` and share one
12-pivot, 128-bit fermionic GECP basis. `GreenCompressionResult` reports the
held-out kernel and Green errors, the discrete `L¹` transfer bound, selected-
core conditioning, and both the observed and theorem-matched ranks. The
Hubbard-like and gapped fixtures reach Green errors `4.16e-10` and `2.91e-10`,
respectively, below the common `5.40e-9` transfer bound.

The sparse results deliberately name their candidate strategy. The clean
quasiparticle/satellite fixture uses a four-entry known transition library and
recovers it to `3.33e-16` held-out error. The historical noisy fixture scans
1,604 frequencies, including the four generating frequencies explicitly
added to its dense dictionary; it is not support-blind. Its eight effective atoms reproduce the
underlying noiseless Green function to `8.56e-9`. The latter is a compact
function representation, not evidence that the four generating frequencies
were uniquely recovered. New method comparisons must use a fixed dictionary
independent of the generating support, with known-library cases separate.

The PSD fixture generates 72 clustered synthetic sensor locations and a
smooth Gaussian covariance. Complete-pivot GECP and diagonal-pivoted Cholesky
choose the same 31 landmarks. The solve-based selected cross has maximum error
`7.78e-7` and relative Frobenius error `1.14e-7`.

All fixtures are synthetic. They are motivated by DLR Green-function
compression, correlated spectral functions, noisy imaginary-time inversion,
and covariance landmark selection, but they are not material-specific fits,
uncertainty-quantified analytic continuation, or optimal sensor-placement
claims.

## Reproducible experiments

These are exploratory reruns from a clean current checkout. They write ignored
scratch outputs, leaving the committed historical evidence and its generation
manifests intact. Byte-for-byte canonical reproduction additionally requires
the artifact's frozen source revision and recorded environment; a later Git
revision correctly produces different provenance.

```bash
uv run python -m kernelgecp.census \
  --config experiments/census.toml \
  --output experiments/scratch/gecp_census_v2.jsonl
uv run python experiments/plot_census.py \
  experiments/scratch/gecp_census_v2.jsonl \
  experiments/scratch/gecp_census_normalized.png --mode normalized
uv run python experiments/plot_census.py \
  experiments/scratch/gecp_census_v2.jsonl \
  experiments/scratch/gecp_census_v2.png
uv run python -m kernelgecp.census \
  --config experiments/census_resolution.toml \
  --output experiments/scratch/gecp_census_resolution.jsonl
uv run python -m experiments.check_census \
  --input experiments/scratch/gecp_census_v2.jsonl \
  --resolution-input experiments/scratch/gecp_census_resolution.jsonl \
  --output experiments/scratch/gecp_census_checks.json
uv run python experiments/exact_surrogates.py \
  --output experiments/scratch/exact_surrogates.json
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/scratch/synthetic_applications.json \
  --plot experiments/scratch/synthetic_applications.png
```

Runtime algorithms require no network access. The historical fixed-grid census
evaluates endpoint-refined grids at 128-bit precision and records decimal strings. Its
four tolerance records per cutoff are prefixes of one strict trajectory down
to `1e-12`. The canonical config requires a clean tracked revision before
generation; exploratory configurations may opt out and record their dirty
state. The historical schema-1 `gecp_census.jsonl` and PNG are preserved;
new schema-2 output uses separate names. Adaptive interval certification is a
separate, more expensive continuous-domain mode.

Schema 2 includes the actual sampled residual maxima (rank + 1 values), signed
original-core determinants, singular values, signed pivots and magnitudes,
boundary distances, dyadic scales/panels, and reflection diagnostics at every
prefix. Each tolerance record recomputes its nearest selected reflection
partners, including fixed-point self-partners, rather than slicing partner
indices from a later prefix. Distances use
`max(abs(dt), abs(domega)/cutoff)`; a match within the recorded tolerance is
numerical proximity, not an exact symmetry theorem. The central panel has
`abs(omega) <= 1`; scale `s > 0` means
`2^(s-1) < abs(omega) <= 2^s`, with sign recorded separately.

Elimination is 128-bit but the composite Chebyshev grid coordinates originate
in float64 (53 bits); both precisions are exposed. Revision, dirty state,
configuration hash, grid shape, Python/dependency versions, mpmath backend
and lock digest describe provenance. No timestamps or runtime-dependent order
are included in byte-reproducible JSONL. Prefix metadata distinguishes the
requested tolerance from the source trajectory's strict tolerance. An
uncertified adaptive terminal upper bound cannot establish convergence.

Normalized plotting uses the deepest compatible trajectory once per cutoff
and every `k/log(1+cutoff)` residual value, including rank zero. Exact numerical
zeros are omitted from the log display and disclosed, not replaced with a
positive floor. Both plotting modes reject explicitly non-grid algorithms,
non-sampled convergence scopes and unsuccessful records, rather than labeling
an adaptive upper bound or a budget-limited rank as converged grid evidence.
Legacy records with absent optional metadata remain accepted. The independent
check script recomputes selected cores at
higher precision and evaluates fixed selected crosses on denser held-out
panel samples using a factored solve. This is a resolution-sensitivity check,
not a denser GECP trajectory, a continuous supremum certificate, or proof of
grid independence.

The separate resolution configuration does run GECP again at order 32 for
cutoffs `1,1000,1000000`, at the same 128-bit precision and `1e-12` tolerance.
The report compares those actual trajectories with order 24, preserving ranks,
residuals, selected-path/tie changes and both generation revisions. A sampled
rank need not decrease when the grid becomes denser. This representative
two-resolution ladder is not residual-adaptive grid refinement or an
exhaustive all-cutoff resolution study. The adaptive interval GECP API and
its success/budget-failure tests remain separate evidence.

### DM primary numerical results (verified; merged PR #119)

The schema-2 primary artifact
[`gecp_census_v2.jsonl`](experiments/data/gecp_census_v2.jsonl) contains 28
converged records at seven cutoffs and tolerances `1e-6,1e-8,1e-10,1e-12`.
Two complete runs from clean generator revision
`c0e8c5093134f369b93781c413750b7f1de46036`, with configuration SHA-256
`c0ac1bbaa7ddb39a2bb4a1c306985386f89f4103ff75cb4e21f99257e1e9f31b`,
produce identical JSONL bytes with SHA-256
`c6854330fcb74007f8f9814fba48c314d34fe26c1cef5461d908ea5182b66bda`.
The `1e-12` ranks are `9,17,37,63,91,119,146` in increasing cutoff order.
All 21 prefixes shared with the preserved legacy artifact have exactly the
same pivots and ranks. The newly rendered rank and normalized-residual plots
show sampled behavior; the normalized curves do not collapse exactly onto
one exponential curve, and this artifact has no numerical zero residuals.

Independent 192-bit point arithmetic recomputes rank-one, middle and terminal
selected cores for each strictest trajectory (21 cores). Relative differences
in stored determinants, signed pivot products and smallest singular values
are all below `9e-27`. At cutoffs `1,1000,1000000`, the stored selected crosses
are also evaluated on `36,3960,15120` held-out point pairs, with maximum errors
approximately `6.31e-14,6.33e-13,8.61e-13`. These are pointwise numerical
checks, not outward-rounded bounds, all-grid maxima or continuous guarantees.
The representative denser GECP runs and combined higher-precision checks pass.
The full branch root gate passes 3,633 Lean jobs, 876 enforced axiom reports
and 430 tests, with clean dependency audit and green push/PR CI. DM is merged
in PR #119 and supports fixed-grid diagnostics/curves (R26/R27), not the
residual-adaptive R25 requirement reopened by DO. Exact merged/pushed closure
gates are recorded in the phase handoff. Numerical results are not continuous
certificates.

### DM representative resolution results

The separate order-32 artifact has three clean, converged `1e-12` records:

| Cutoff | Grid shape | Order-24 / order-32 rank | Order-32 sampled residual |
| --- | --- | --- | --- |
| 1 | 63 × 63 | 9 / 9 | `1.78383e-13` |
| 1,000 | 621 × 683 | 63 / 63 | `9.59729e-13` |
| 1,000,000 | 1,279 × 1,303 | 146 / 146 | `9.71664e-13` |

[`gecp_census_resolution.jsonl`](experiments/data/gecp_census_resolution.jsonl)
records clean generator revision `ccd5c8c9ad07ee4548d77c5e03d3ccefa50eef61`,
configuration hash
`bf49e21ce208a1336bfc83d3e2ec955adf627e64e8056c8be9986d855ee2d3be`,
and artifact SHA-256
`7251cbb528179cc5cbfa6eb5635b8fa280e087d1e04db1348b56a9065026d588`.
Matching ranks are representative numerical observations, not evidence of
unchanged pivot paths, exhaustive resolution convergence or grid independence.
Actual common selected prefixes have lengths `3,0,0`; all three order-32
sampled terminal residuals are larger than their order-24 counterparts.
The [combined report](experiments/data/gecp_census_checks.json) preserves both
generation revisions/source hashes and passes 30 selected-core checks at 192
bits (21 primary and 9 resolution). All stored-determinant, signed-pivot-product
and singular-value relative differences are below `9e-27`. This comparison
does not confuse held-out evaluation of a fixed cross with rerunning denser
GECP, or equate point arithmetic with interval certification.

### Historical diagnostics and current proof scope

These tools were developed while the exact rate was open. DG/DH now prove
the original-kernel finite exact-complete-run rate; the diagnostics below are
not substitutes for that theorem and do not extend it to inexact pivots.
`analyze_census_file` computes
cutoff-normalized block ratios from the stored decimal strings, while
`certify_fermionic_block_contraction` uses outward-rounded interval residuals
along a returned approximate-pivot trajectory:

```python
from decimal import Decimal
from pathlib import Path

from kernelgecp import (
    analyze_census_file,
    certify_fermionic_block_contraction,
    fermionic_two_corner_residual,
)

blocks = analyze_census_file(Path("experiments/data/gecp_census.jsonl"))
cutoff_one = certify_fermionic_block_contraction(1.0, 2)
assert cutoff_one.half_contraction_certified
base_value = fermionic_two_corner_residual(Decimal("1"), Decimal("0.5"), Decimal("0"))
assert abs(base_value) <= Decimal(1) / 4
```

`BlockContractionEvidence.continuous_certified` is always false because the
canonical census is finite-grid evidence. A certified block result applies to
the frozen rounded residual model of the returned approximate-pivot trajectory; it
does not establish a cutoff-uniform exact-GECP theorem or resolve arbitrary
tie choices.

The two-corner evaluator is different: it independently implements the exact
closed form formalized in Lean. Lean proves that the first two corner pivots
are complete for every positive cutoff and that this residual contracts by
one half for `0<Λ≤1`. Neither the evaluator nor that base theorem supplies a
dyadic scaling argument for larger cutoffs. The later exact finite-run rate
uses the independent local-minor/panel-count proof, not that historical route.

For exact structural experiments, `first_exact_gecp_sign_coherence_failure`
checks every selected residual cross using `Fraction` arithmetic and returns
the earliest step, pivot, coordinate, and negative product. Returning `None`
means only that the supplied finite matrix passed this exact diagnostic. Lean
separately proves that this residual condition is equivalent to a
four-bordered-minor sign condition on the original kernel, proves the required
all-orders fermionic sign theorem, and propagates the resulting
nonexpansiveness through realized exact complete-pivot sequences. The finite
diagnostic remains an independent regression, not the proof.

The exact surrogate regression also enumerates every row and column
permutation at sizes two through four for `q = 2/3`. It checks the determinant
orientation convention used by the Lean strict-sign-regularity transfer. This
is a finite exact cross-check rather than the proof of the now-delivered
all-orders exponential minor theorem.

Lean now also proves a functional geometric-mean reduction for realized exact
complete-pivot sequences: every rank-`n` domain residual magnitude raised to
the `n`th power is bounded by the absolute determinant of the selected
`n × n` core. The exact two-pivot rational example independently exercises
this identity. This adds no Python API and does not turn a finite numerical
determinant history into a certified convergence rate. The required quantitative
determinant budgets and their actual-GECP composition are now proved in Lean;
the numerical determinant history is not that proof.

A second exact regression constructs the first `n` exponential Taylor features
for rational row and column nodes through size six. Its determinant agrees
exactly with the two Vandermonde products times `prod(1/k!)`. This exercises the
Lean principal-term theorem but does not numerically stand in for the
generalized-term or infinite-series arguments, subsequently proved in the
all-orders exponential total-positivity modules.

For repository verification, pytest is launched through the synchronized
project interpreter rather than a standalone editable-install console script:

```bash
uv sync --all-extras --frozen
./scripts/verify.sh
```

The root command also performs formatting, typing, Lean axiom, placeholder,
and dependency-audit checks. This developer workflow does not alter the public
Python API or any numerical certification boundary described above.

Phase BF adds no public Python API. Its 120-digit complete-pivot regression is
an independent check of the formal positive-unit-band tail envelope; it does
not certify mixed low-rank-plus-tail determinants or a global GECP rate.

Phase BG also adds no public Python API. Its 140-digit mixed-coefficient sweep
selects the next formal theorem; the observed common sign is not exposed as a
runtime certificate or claimed as a proved cutoff-uniform rate.

Phase BH also adds no public Python API. It formally proves the corresponding
finite positive-Taylor coefficient signs and retains the exact rational check
as independent evidence. It does not expose infinite-tail coefficients or a
new runtime convergence certificate.

Phase BI also adds no public Python API. Its determinant-polynomial bridge is
formal infrastructure for the forthcoming infinite-tail argument; the exact
rational regression remains evidence rather than a runtime certificate.

Phase BJ also adds no public Python API. The 120-digit truncation experiment
independently checks the proved infinite-tail coefficient limit; it does not
yet provide a quantitative runtime convergence certificate for GECP.

Phase BK also adds no public Python API. Its 350-digit sweep independently
checks the formal determinant-pencil amplification law; it is not a runtime
certificate and does not yet expose a completed GECP rate.

Phase BL also adds no public Python API. Its 250-digit unit-band sweep
independently checks the formal amplified-kernel envelope and ordered-minor
positivity; it is research evidence rather than a runtime certificate, and it
does not identify original and amplified complete-pivot paths.

Phase BM also adds no public Python API. Its exact rational exhaustive-minor
check independently exercises the formal total-positive Hadamard argument.
The resulting fixed-base determinant decay is a Lean theorem, not a new
runtime GECP convergence certificate; GECP composition and global band
assembly remain research obligations.

Phase BN also adds no public Python API. Its 180-digit finite-grid GECP sweep
independently checks the formal positive-unit-band half contraction, but is
not a runtime certificate for arbitrary continuous pivot traces or the full
signed-frequency domain.

Phase CK adds no public Python API. Its exact rational compound-pencil audit
and Lean generalized-spectrum certificate are research infrastructure, not a
runtime real-rootedness or convergence certificate.

Phase CL also adds no public Python API. Its exact fresh-order-two discriminant
audit and formal quadratic factorization narrow a proof obligation; they do
not expose a runtime spectral or convergence certificate.

Phase CM also adds no public Python API. Its exact weighted-split and
conditioned-coordinate audits independently check a formal reverse-covariance
theorem; they do not certify that an arbitrary runtime conditioned factor has
the required positive ratio-coordinate representation.

Phase CN also adds no public Python API. Its rank-two chart certificate is a
formal matrix-pencil interface, while the exact cyclic-order audit selects
the next theorem; neither is exposed as a runtime spectral certificate.

Phase CO also adds no public Python API. Its two-sided basis-change and
discriminant identities justify cutoff-specific formal charts; they do not
construct a runtime chart or certify conditioned numerical pencils.

Phase CP also adds no public Python API. Its explicit pole basis and
reciprocal-ratio certificate are formal proof infrastructure. The exact
rational regression checks the chart identities and collision boundary, but
does not discover or certify a separating pole for arbitrary runtime data.

Phase CQ also adds no public Python API. Its oriented/scaled separator theorem
is formal proof infrastructure. The exact rational search constructs
certificates only for the deterministic research corpus and is not exposed as
a runtime separator or arbitrary-data certification routine.

Phase CR also adds no public Python API. Its selected-core-square identity is
formal proof infrastructure, and the exact augmented-core reconstruction is a
deterministic regression rather than a runtime sign or separator certificate.

Phase CS also adds no public Python API. Its common-schedule theorem is formal
proof infrastructure, while the exact sign audit records both an invariant and
a failed stronger conjecture; it does not certify arbitrary runtime factors.

Phase CT also adds no public Python API. Its endpoint theorem is formal
topological infrastructure, and the exact midpoint sweep is a deterministic
research regression. It does not certify determinant nonvanishing along an
arbitrary runtime path or expose a numerical homotopy checker.

Phase CU also adds no public Python API. Its exact regression independently
reconstructs the two rational endpoint determinants and ordered path samples;
the delivered result is a formal obstruction to an overbroad proof strategy,
not a runtime singularity detector or certificate for user-supplied nodes.

Phase CV also adds no public Python API. Its deterministic exact regression
certifies one research counterexample—kernel entries, Schur residual, two
complete pivots, and crossed determinant signs. It is not exposed as a
runtime pivot or chirotope diagnostic.

Phase CW also adds no public Python API. Its exact regressions test a formal
aggregate determinant-pencil criterion on conditioned research fixtures; they
do not provide a runtime Lorentzian or convergence certificate for arbitrary
user data.

Phase CZ adds a standalone exact research audit, without changing the stable
`kernelgecp` API:

```bash
uv run python experiments/transfer_nodes.py --cases 100
```

The JSON records the revision, seed, exact-rational arithmetic, pivot tie
policy, stopping status, and every failed fixed-node certificate. The default
audit checks all fresh row/column pairs in 100 deterministic Taylor cases;
`--cases 500` reproduces the broader research run. Its neighboring-node
coverage result certifies only the finite samples, not arbitrary GECP runs
or the infinite Taylor tail.

The synthetic-application JSON excludes timestamps, records the implementation
Git commit and configuration hash, and is byte-identical across repeated runs.
The PNG is a human-facing summary and is not used as the source of numerical
truth.

Phase DA changes no public Python API. Its actual-exponential decimal
regressions independently check a formal one-pivot sign theorem; the separate
exact depth-ten fixture rejects a research localization hypothesis. Neither
is a runtime convergence or arbitrary-prefix sign certificate.

Phase DB adds a research certificate generator without changing `kernelgecp`:

```bash
uv run python scripts/generate_coverage_certificate.py --check
```

The command checks byte-for-byte regeneration of the committed Lean integer
data and reports exact arithmetic, pivot/tie policy, conditioning scope,
revision, and stopping metadata. Without `--check` it emits the Lean source
to stdout; it does not write or overwrite a file. Lean independently checks
the arithmetic and its connection to the actual conditioned Taylor factors.
This is a fixed counterexample certificate, not a runtime convergence or
arbitrary-input sign certificate. The actual-exponential 160/320-decimal
regressions remain numerical evidence only.

Phase DC changes no public Python API. Its reproducible research-only stress
audit is:

```bash
uv run python experiments/companion_envelope.py
uv run python experiments/companion_envelope.py --depth 64 --cutoff 256
```

The default run checks nine configurations at 320 and 640 decimal digits.
Pivots come from the original kernel; the two-sided conditioned companion
uses both transported feature factors. JSON records the revision/dirty state,
precision, pivot/tie diagnostics, stop reason, and decimal-string determinant
ratios. Precision disagreement is explicitly uncertified and exits nonzero;
zero pivots or determinants are rejected rather than interpreted as success.
No condition number is estimated; precision agreement is not an interval or
continuous-domain certificate. The exact finite local-minor regression is
separate from this numerical audit and from the new Lean theorem.

Phase DD changes no public Python API or package release. Its independent
regressions can be run with:

```bash
uv run python -m pytest tests/test_common_residual_gauge.py
```

The exact tests include all fresh minor orders, non-complete prefixes,
boundary times, signed frequency labels, repeated samples, and transport of
the original chosen pivots under frequency reflection. Separate 160/320-digit
checks compare actual exponential evaluations to exact rational values.
The new Lean theorem is an all-real structural result, but the package does
not expose a runtime sign-gauge or convergence certificate. The dyadic budget,
its iteration, and the final convergence bound are still research tasks.

Phase DE changes no public Python API or release. Its focused independent
checks are:

```bash
uv run python -m pytest tests/test_local_panel_budget.py
```

They separate exact finite-feature/exponent arithmetic from actual-exponential
precision checks. Zero or non-finite sampled determinants are rejected as
uncertified; agreement at two precisions is not interval certification. The
large target-order local budget is proved in Lean, not inferred from a dense
floating-point determinant. No runtime arbitrary-prefix convergence or panel
certificate is exposed; pivot-witness iteration and counting remain open.

Phase DF adds formal rate theorems but no public Python API or new release.
Its focused checks are:

```bash
uv run python -m pytest tests/test_local_panel_rate.py
```

Exact rational arithmetic independently constructs backwards witnesses through
all subsets of a six-step complete trace, including ignored in-panel pivots,
outside-panel pivots, and matched-rank row changes. Initial local budgets are
enumerated separately, not inferred from residual behavior. Three actual
exponential witness determinants are compared at 160/320 digits; zero or
nonfinite values are rejected, and precision agreement is not an interval
certificate. The all-orders finite dyadic convergence theorem is in Lean.
This does not certify an approximate numerical pivot as a continuous complete
maximum or expose a runtime convergence certificate. General real-cutoff
and accuracy corollaries still need formalization.

Phase DG adds no Python API or package release. Focused checks are:

```bash
uv run python -m pytest tests/test_real_cutoff_rate.py
```

Cutoff covers and block ranks use exact rational/integer arithmetic, including
dyadic boundaries. Envelope, real-tolerance, and actual first-pivot witness
values are compared at 160/320 digits as decimal strings. Witness residuals
must be positive and finite. The finite grid observation that a first corner
is maximal is not a runtime continuous-domain certificate; the corresponding
continuous statement is proved in Lean. The new rate does not certify finite-
precision approximate pivots, expose a runtime rate certificate, or change
stopping/conditioning/tie metadata. Its small-cutoff obstruction explicitly
separates the `Lambda>=1` theorem from literal all-positive-cutoff G1 wording.

Phase DH changes no public Python API. Focused independent checks are:

```bash
uv run python -m pytest tests/test_finite_run_closure.py
```

Each prefix test receives only its stored finite pivots and states; no later
pivot is supplied to justify that prefix's error. Exact selected-core solves
match iterative residuals and signed Green sums. Rank-zero/rank-one examples
stop at an exact complete zero without a zero-pivot update. Actual exponential
crosses and nonzero core determinants are independently compared at 160/320
digits; no explicit numerical inverse is formed. Precision agreement does
not constitute interval or continuous-domain certification. The new Lean
finite-run and exact-stop results do not turn floating-point zeros, approximate
pivots, or certification-budget exhaustion into exact mathematical hypotheses.
Stopping, conditioning, precision, tie, and certification metadata are unchanged.

At Phase DI's [audit](COMPLETION_AUDIT.md) entry, the historical census stopped
at `1e-10`, lacked the full requested geometry/normalized curves, and the
non-GECP result types needed fuller status metadata. DL supplies result
metadata; DM's primary/resolution data and PR #119 support fixed-grid
diagnostics/curves (R26/R27). At DM closure the historical audit recorded
39/40 supported groups; DN then delivered R31. DO's independent audit rejects
the resulting 40/40 assumption and reopens the adaptive-grid R25 scope. At DO
closure its 39/40 count was 97.5% checklist coverage, not weighted effort.
DP now records the actual adaptive generations and passing independent
post-run checks described below. DQ's original-scope audit and boundary repairs
are described below; final branch gates pass and PR #123 is merged.
Those algorithms are not behavior silently supplied by the exact theorem.

The public numerical API is unchanged in DI. Packaging excludes the two local
duplicate working-note files; wheel/sdist were rebuilt and the wheel was
smoke-tested in isolation outside the checkout. No release was published.
The axiom checker now enforces permitted dependencies and declaration coverage;
it also propagates Lean failure. New focused checks are:

```bash
uv run python -m pytest tests/test_axiom_output.py tests/test_real_separated_approx.py tests/test_docs.py
```

## Original-kernel fill-distance theorem (DJ)

For a PSD kernel on any metric domain, the Lean interface
`PositiveDefinite.diagonalRun_error_le_original_fillDistance` accepts an actual
finite successful diagonal `GECP.Run`, original-kernel coordinatewise
`LipschitzWith L`, and a selected-node radius-`h` cover. It proves the absolute
residual bound `2*L*h` everywhere. It does not require a Lipschitz residual,
nor does it require each diagonal pivot to be maximal. Positivity of `R` and
`K-R` is derived from the run. Continuous kernels on nonempty compact domains
have a diagonal complete maximum; a zero maximum is termination.

This is a formal exact-arithmetic interface, not an automatic certificate for
floating-point Cholesky or a derivation of unrecorded numerical Lipschitz
constants/fill distances. It does not apply to the nonsymmetric fermionic
kernel by calling it PSD. Weak PSD allows off-diagonal maximizing ties; strict
SPD diagonal-necessity remains the separate R08 target. Numerical APIs and
metadata are unchanged in DJ; R24 remains incomplete. Focused independent
checks: `uv run python -m pytest tests/test_kernel_geometry.py`.

## Strict-SPD exact complete runs (DK)

`PositiveDefinite.posDef_completeRun_iff_pivotedCholesky` now works with the
actual dependent `GECP.Run` on arbitrary domains, not just finite canonical
trace definitions. From original-kernel `Matrix.PosDef K`, every nonzero exact
complete pivot is derived to be diagonal. The independent `DiagonalMaxRun`
rule then describes exactly the same successful finite runs. Diagonal ties
are allowed; a deterministic implementation can refine them without changing
this set-valued theorem. Weak PSD still permits off-diagonal maximizing ties.

`posDef_completeRun_error_le_original_fillDistance` supplies the original-kernel
`2*L*h` bound without requiring a separate diagonal-run proof. On a nonempty
compact continuous domain, `posDef_completeRun_exists_diagonal_completePivot`
provides an actual next complete maximum, possibly zero (stop). Domain subsets
can be represented by the ambient subtype. No automatic positivity, exact
maximality, Lipschitz constant, fill-distance or rounding certificate is added
to floating-point results. Numerical metadata remains an open R24 requirement;
Python APIs are unchanged. Focused checks:
`uv run python -m pytest tests/test_strict_kernel.py`.

## Numerical metadata (DL)

`GECPResult`, `CholeskyResult`, and `SparseResult` expose `result.metadata` with
typed `NumericalMetadata`, `ConditioningDiagnostic`, and
`RefinementDiagnostic` records. Use `metadata.as_dict()` or
`metadata.canonical_json()` for deterministic JSON; the latter rejects NaN and
infinity. Empty/unavailable estimates use explicit statuses and null values,
not fabricated condition numbers. Legacy manually constructed result records
default to `not_recorded` metadata.
Metadata snapshots include convergence and stop reason too, so application
summaries do not lose an underlying algorithm's distinct stopping condition.

| Path | Actual arithmetic / scope |
| --- | --- |
| GECP grid | float64 or configured mpmath precision; sampled-grid convergence, no continuous certificate |
| GECP interval adaptive | mpmath; interval upper-bound scope; `certified` only on certified-tolerance convergence; near-tie multiplicity not estimated |
| Cholesky | float64, even if the shared config requests 128 bits; actual/requested bits differ explicitly |
| Sparse OMP/refinement | float64; training-sample relative infinity error only, no held-out/continuous/support certificate |

Conditioning histories identify their stage: selected original Cholesky core,
OMP design, candidate refinement Jacobian, or returned sparse design. Each
estimate records numerical rank, smallest singular value and condition number
when applicable. A rank-deficient or failed diagnostic is not a successful
conditioning certificate. GECP preserves its canonical high-precision singular
value strings and explicitly reports that its metadata view estimates only
the smallest singular value, not an uncomputed condition number.

Cholesky recomputes `original - factor @ factor.T`; `residual_matrix` is that
actual error and `converged` requires its infinity norm to meet `tol`.
Roundoff or accepted near-symmetry can instead give
`roundoff_or_symmetry_limit`. Materially invalid PSD residuals raise
`NumericalBreakdownError`; PSD itself is not certified by floating-point code.
`residual_history` retains the working diagonal maxima; metadata's
`returned_error` is independently recomputed from the returned factor.
Zero tie tolerance means exact ties, and selected zero pivots cannot be reused.

Sparse `training_error` aliases the legacy `validation_error` (which has always
been measured on training samples). Candidate exhaustion is an explicit stop.
Nonlinear refinement is accepted only if finite and nonworsening in training
infinity error. Optimizer status/success/evaluation budget are separate from
representation convergence: an exhausted optimizer can return a useful fit,
and a successful optimizer can still miss the requested tolerance. Failed or
worsening refinement retains the OMP fit with an explicit diagnostic.
`iterations` counts OMP selections; optimizer evaluations are recorded separately.
`tie_rtol` and `refinement_max_nfev` are optional sparse controls. Candidate,
time and target fingerprints plus settings support deterministic repetition.

Fixed-order interpolants/Taylor constructions report non-iterative scopes:
interpolation conditioning is not estimated; the Taylor evaluator uses mpmath
internally but returns float64 arrays. Neither claims an outward-rounded
numerical error certificate. Barycentric construction uses a fixed random
state. Synthetic application outputs now use schema version 2 and retain
nested algorithm metadata. Historical schema-1 committed JSON/PNG hashes and
revision remain valid historical evidence, not silently rewritten output.

Focused checks: `uv run python -m pytest tests/test_result_metadata.py`.

## Six-method spectral comparison (DN, verified; merged PR #120)

`ComparisonConfig` and `run_spectral_comparison` live in
`kernelgecp.comparison`. The runner is
`uv run python -m experiments.spectral_comparison`; plotting uses
`uv run python -m experiments.plot_spectral_comparison`. Canonical runs reject
a dirty tracked tree, record revision/configuration/lock/source hashes and
versions, reject nonfinite JSON, and atomically replace output. The SciPy
minimum is now 1.15 because its public `interpolate.AAA` implementation is used.

The primary cohort gives GECP, AAA, genuine pole-relocation vector fitting,
strided-Hankel ESPRIT/Prony, bounded variable-projection nonlinear least
squares, and unrefined greedy pursuit the same observed uniform time samples.
All returned weights are least-squares fitted on those observations. GECP's
universal nodes depend on the kernel alone; OMP's fixed dictionary is
independent and the four-atom primary spectra are off-grid. A separate
known-transition-library control is prior-informed, not a blind comparison.

AAA and vector fitting see a cubic-spline/Gauss-quadrature Matsubara transform
of those samples. Quadrature and source-resolution discrepancies are empirical
diagnostics, not rigorous interpolation/noise bounds. Ground-truth rational
data is used only to assess transform error, never supplied to the fit.
Nearly-real poles may be projected only within a recorded tolerance;
nonphysical poles and budget failures remain visible. Prony permits nodes
larger than one for negative physical frequencies. No stability-half-plane
flip or hidden nonlinear refinement is applied to any other algorithm.

Held-out infinity errors use disjoint noiseless time coordinates. Passes are
absolute `error <= epsilon`; genuine relative errors use the held-out signal
norm (null for zero signal). Training convergence, native rational/Prony
convergence, optimizer termination and held-out success are separate fields.
The smallest **tested** passing support is compared with the sampled universal
kernel rank; it is not the globally minimal support, a support-identification
certificate, or a continuous-error guarantee. Results are 53-bit fits; only
the separately labeled universal node generator uses the configured precision.

The frozen artifact contains 192 fits (four fixtures, eight budgets, six
methods), and the JSON and reviewed plot repeat bytes. Universal kernel ranks
at cutoff eight are 11/12 at `1e-6`/`1e-8`. At `1e-8`, AAA/VF/Prony/NLS reach
four tested atoms on the clean off-grid case; VF/Prony/NLS also reach four on
the noisy case and six on each continuous-density reference. AAA's other
physical fits fail and remain in the plot/data. Density references are fixed
normalized 2001-node quadratures; reported errors are not against the exact
continuous integral. Full root (513 tests), exact-head push/PR CI and
installed-wheel smokes pass, including the minimum supported SciPy 1.15.
Those are DN's historical verification results and PR #120 is merged. DO's
subsequent independent audit, 569-test verification and merged PR #121 are
recorded above. DP's actual adaptive-grid generations are now available as
described below, with independent MP192 checks passing. DQ subsequently
performs the original-scope audit and verifies its repaired-source delivery.
No package is published.

## DP residual-adaptive sampled grids (verified; merged PR #122)

`kernelgecp.adaptive_grid.run_adaptive_grid(cutoff, config=AdaptiveGridConfig(...))`
is a separate experiment controller, not the interval `pivot="adaptive"` mode.
It runs MP GECP on explicit nested quarter-lattices, probes odd-eighth points,
marks influential panels and their reflection mates, and bisects selectively.
Every previous binary64 coordinate is retained exactly. Elimination uses
MP128 by default; grid coordinates remain 53-bit repr-decimal values. MP192
factored selected-cross checks are higher-precision point arithmetic, not
interval enclosures or continuous supremum bounds.

Import `AdaptiveGridConfig` and `run_adaptive_grid` from
`kernelgecp.adaptive_grid`. The public configuration defaults are five maximum
stages and generation margin `0.25`; they are bounded research defaults, not a
promise that every cutoff converges. The frozen seven-cutoff study instead
uses `experiments/adaptive_census.toml`: nine stages, margin `0.01`, MP128
elimination, MP192 checking, maximum rank 256, maximum 512 nodes per axis,
and four highest marked panels per axis before reflection mates. Its four
tolerances extend through `1e-12`, and its source trajectory tolerance is
`1e-14`. These controls were fixed after labeled pilot failures, before
canonical final validation.

The returned `AdaptiveGridCase` preserves every stage, leaf panel, explicit
axis, mark, generation tolerance, budget, pivot/core history and probe score.
Two consecutive strictest-tolerance stage passes are required. Requested
tolerance prefixes advance from the first training-passing rank until a tested
prefix also passes its residual/surplus probes; this is not a global minimum
rank claim. Only after all refinement and prefix choices, six seventh-fraction
points per final panel supply independent validation. They never select a
grid or rank. `fit_passed` and `controller_passed` remain separate: failed
rank/stage/node budgets cannot be reported as convergence even if a local fit
passes. Final validation failures remain visible and are not reused for tuning.

`uv run python -m experiments.adaptive_census` reads the frozen TOML controls and
writes schema-3 explicit-axis JSON. Canonical generation requires clean tracked
source, snapshots source/lock/config hashes and environment versions, and
rejects source/revision/lock/environment changes during a run. Use a fresh CLI
process from the frozen checkout for canonical generation: disk guards cannot
attest stale imported functions in a long-lived API process. Canonical `root`
must match the executing runner, and untracked experiment-package bootstrap
code is included in the clean-source gate. The independent
`experiments.check_adaptive_census` checks geometry, nesting, marks, source
blobs, selected cores and final crosses; `experiments.plot_adaptive_census`
shows normalized sampled curves, independent validation and stage growth.
These tools preserve failures and distinguish zero samples from positive
log-plot values. The historical fixed-grid artifacts are not overwritten.

The canonical [`adaptive_census.json`](experiments/data/adaptive_census.json)
records all seven cases and all 28 tolerance prefixes passing the controller
and its separate final sampled checks. For cutoffs
`1,10,100,1000,10000,100000,1000000`, actual stage counts are
`3,2,3,4,4,5,6`, final training grids are
`25×25,41×57,109×89,122×121,170×169,238×217,258×257`, and chosen `1e-12`
ranks are `9,17,36,64,97,124,151`. Every stage preserves canonical MP pivot,
sampled residual, original selected-core determinant and smallest-singular-
value histories. Final-prefix records also contain boundary distances,
dyadic scales, nearest reflection partners and their stated matching metric.

Two fresh CLI generations from clean revision
`3cee4186fddf990abdd33dbe24c46ff908b1923f` produce identical bytes, SHA-256
`730fc70a9fb3d12889361e6b57d1bc4cce63e07723452bb5b8709f83b5bca93c`.
The JSON records configuration hash
`ac47ec5146f5dc50916513917518d2bbfcb31da622b1f659225568cc48d641a1`,
all source/lock hashes, dirty state and environment. To reproduce the
historical bytes, use that clean frozen revision and its recorded environment
in a fresh CLI process. A run from a later checkout correctly records a
different `git_commit`; byte identity across different revisions is not the
reproducibility claim. Write exploratory reruns to a separate output path,
without replacing historical evidence:

```bash
uv run python -m experiments.adaptive_census \
  --config experiments/adaptive_census.toml \
  --output /tmp/gecp_adaptive_rerun.json
```

The normalized PNG has been rendered, visually inspected at 2240×1440 and
repeated byte-for-byte, SHA-256
`dd1d9f08a00e1ddbd9f4d5ac0224b62241e91d42781857b4537f0c81f352f3d6`.
It shows seven curves and 28 validation markers; one sampled zero is omitted
from the logarithmic view without inventing a floor. The independent MP192
[`adaptive_census_checks.json`](experiments/data/adaptive_census_checks.json)
passes full seven-cutoff/28-prefix coverage and revision-blob checks, all 21
selected-core determinant/pivot-product/singular-value checks, and all 28
complete final sampled crosses within their requested tolerances. Its SHA-256
is `4949243cb779e5fd1132df2768a696465573046f4f0210d872ad2911b4c2be03`.
At final branch head `3c38183e5c2a7d2859867bb7acff1353ee98fb76`, the DP
exact full root gate passes 713 tests in 150.98 seconds, 3,634 Lean jobs,
all 876 axiom reports, formatting of 95 files, mypy on 21 sources and the
clean dependency audit. Source/wheel audits include the three DP artifacts
and tests, exclude protected notes/caches/scratch, and verify all 21 packaged
Python hashes plus README. Installed-wheel checks run outside the checkout
and pass in both pinned Python 3.12.12 / NumPy 2.5.2 / SciPy 1.18.0 /
mpmath 1.4.1 and minimum SciPy 1.15.0 / NumPy 2.4.6 environments. Each checks
all 45 exports, actual adaptive success and budget failure, and 48 six-method
fits; final archives and both environments pass again after rebuilding that
exact branch head. Push CI `40017937659` and PR CI `40017945103` pass all four
checks. Earlier development/code-artifact roots took 151.22/150.40 seconds;
they are not the final branch timing.
[PR #122](https://github.com/tripp-smith/gecp-kernel-structure/pull/122) merged
at 2026-10-23 20:20:48 UTC as
`3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f`, with local main fast-forwarded to
origin. The original clean generation revision `3cee4186` is confirmed as an
ancestor; its provenance remains unchanged. Exact pushed-main closure
metadata/gates subsequently pass at `5a5fef0`; at DP closure the next required
phase was the fresh original-SPEC audit and the effort estimate was about 99%,
not a checklist count or completion verdict. DQ below delivers that audit and
the final repairs. Successful
sampled validation is not a continuous residual-supremum bound, an exact
complete-pivot trajectory, or a new Lean rate theorem. No package publication
is claimed.

## DQ boundary safety and standalone certificate metadata

`certified_pivot` is conditional on **exact callback point values of one fixed
function** and coordinatewise Lipschitz constants valid throughout each cell,
not just at its sample. It does not enclose arbitrary callback evaluation error.
DQ computes radii/upper bounds and gap comparisons using exact binary-rational
arithmetic before outward float64 export, and rounds the lower/upper eta ratio
down. One-ULP cells and subnormal uncertainty no longer disappear through
midpoint, width or product rounding. A non-splittable unresolved cell returns
`coordinate_resolution` with `certified=False` and its upper bound retained;
an unrepresentable upper bound raises `NumericalBreakdownError`.

Standalone `PivotCertificate.metadata` now records arithmetic, actual point,
interval enclosure and exported-bound/coordinate precision, distinct sampled
tie counts, stop reason, objective/certification scope and scalar conditioning
as not applicable. Interval certificates still require their separately
documented valid interval callback. They restore the caller's exact interval
bit precision, including failure paths; high internal precision does not make
53-bit exported coordinates or bounds high precision.

The Taylor factory and `dyadic_taylor_rank_bound` choose the smallest exact
integer parameters covering the **actual binary64 cutoff/tolerance**, including
values immediately adjacent to powers of two and subnormal tolerances.
`exp_family_rank_bound` computes its published integer floor count without
floating logarithms. `dlr_rank_bound` uses a conservative binary-log envelope
for the central count, so it may be larger than the earlier rounded natural-log
helper; neither it nor a sampled GECP rank is the Lean construction's rank.
Exported Taylor error bounds never silently underflow to exact zero. These
fixes do not certify rounded Taylor evaluation or an inexact numerical GECP
trajectory with the ideal exact-run theorem.

DQ also updates the development-only cache-lock dependency to filelock 4.2.1
after the live dependency gate reports advisories. The numerical stack and all
historical experiment manifests remain unchanged. DQ's source-level audit is
recorded in [FINAL_AUDIT.md](FINAL_AUDIT.md). Exact code head `f2ffd5e` passes
757 tests in 150.01 seconds and all Lean/consumer/axiom/style/type/security
gates. Push/PR CI `40050118292`/`40050129457` pass all four checks. Rebuilt
archives match all 21 Python/153 Lean sources and claim docs; installed wheels
pass both pinned and minimum-SciPy environments, 45 exports, adaptive success/
budget failure, 48 six-method fits and the new boundary/metadata checks each.
Final documentation branch `0a6d363` repeats the full root (757 tests in
149.97 seconds), rebuilt archives and both installed environments; push/PR CI
`40050605960`/`40050610230` pass all four checks. PR #123 merges at
2026-10-24 01:41:41 UTC as `d9978607b88a45d98d917f44318bcaf99b4c3239`.
Exact pushed-main `e8c35ed` passed the full root (757 tests in 150.31 seconds),
CI `40051138962`, archive inspection and both installed environments. Required original-
SPEC scope is complete at these explicit boundaries; no package publication,
external peer review or ideal/inexact-run numerical theorem is implied.

## Release provenance

The implementation agent's exposed model family and mode, goal-accounting
time and token snapshot, and the assumptions behind the API-equivalent cost
estimate are documented in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They describe
the historical v1 implementation, completed original-SPEC research goal and
later documentation-only task boundary, not package runtime requirements,
phase-only costs or benchmark performance.
