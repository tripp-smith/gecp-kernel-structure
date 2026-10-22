# Application guide

The public package supplies stable fermionic-kernel evaluation, deterministic
finite-grid GECP, solve-based cross approximation, diagonal-pivoted Cholesky,
conditional Lipschitz pivot certificates, numerical separated interpolation,
and sparse spectral recovery.

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
For PSD matrices, `pivoted_cholesky` and `gecp_matrix` use the same
diagonal-preferring lexicographic tie rule; the recursive equivalence of those
pivot policies is separately proved in Lean.

## Certified pivots

`certified_pivot` performs deterministic float64 cell subdivision from a
supplied coordinatewise Lipschitz bound. Its result is conditional on that
analytic bound. `interval_certified_pivot` instead uses `mpmath` interval
arithmetic, outward-rounded float endpoints, and an optional interval-gradient
mean-value enclosure.

`gecp(FermionicKernel(...), config=GECPConfig(pivot="adaptive", ...))` applies
the interval method to every continuous residual pivot and retains each
successful `PivotCertificate`. Cell-budget exhaustion returns
`stop_reason="uncertified_pivot"`; it never silently becomes convergence.
These certificates audit numerical pivot quality. The Lean derivative theorems
prove the underlying kernel formulas and domain bounds, but Python execution is
not itself a formal proof of a continuous GECP rate.

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
`1e-8` with two atoms. AAA and vector fitting are deferred extension methods,
not placeholder APIs.

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

```bash
uv run python -m kernelgecp.census \
  --config experiments/census.toml \
  --output experiments/data/gecp_census_v2.jsonl
uv run python experiments/plot_census.py \
  experiments/data/gecp_census_v2.jsonl \
  experiments/data/gecp_census_normalized.png --mode normalized
uv run python experiments/plot_census.py \
  experiments/data/gecp_census_v2.jsonl \
  experiments/data/gecp_census_v2.png
uv run python -m kernelgecp.census \
  --config experiments/census_resolution.toml \
  --output experiments/data/gecp_census_resolution.jsonl
uv run python -m experiments.check_census \
  --input experiments/data/gecp_census_v2.jsonl \
  --resolution-input experiments/data/gecp_census_resolution.jsonl \
  --output experiments/data/gecp_census_checks.json
uv run python experiments/exact_surrogates.py
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/data/synthetic_applications.json \
  --plot experiments/data/synthetic_applications.png
```

Runtime algorithms require no network access. The canonical census evaluates
endpoint-refined grids at 128-bit precision and records decimal strings. Its
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
in PR #119 and supports R25–R27; exact merged/pushed closure-head gates will
be repeated before handoff. Numerical results are not continuous certificates.

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

For the open GECP-rate question, `analyze_census_file` computes
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
the explicitly returned interval-certified approximate-pivot trajectory; it
does not establish a cutoff-uniform exact-GECP theorem or resolve arbitrary
tie choices.

The two-corner evaluator is different: it independently implements the exact
closed form formalized in Lean. Lean proves that the first two corner pivots
are complete for every positive cutoff and that this residual contracts by
one half for `0<Λ≤1`. Neither the evaluator nor that base theorem supplies the
still-missing dyadic scaling argument for larger cutoffs.

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
determinant history into a certified convergence rate; a quantitative
determinant-decay theorem is still required.

A second exact regression constructs the first `n` exponential Taylor features
for rational row and column nodes through size six. Its determinant agrees
exactly with the two Vandermonde products times `prod(1/k!)`. This exercises the
new Lean principal-term theorem but does not numerically stand in for the
missing generalized-term or infinite-series arguments.

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
metadata; DM's verified primary and resolution results above support R25–R27
and are merged in PR #119. The audit records 39/40 supported groups (97.5%
checklist coverage, not weighted effort). Listed AAA/vector-fitting/Prony
comparisons remain required in R31, followed by the final forty-group audit;
they are not behavior silently supplied by the new exact theorem.

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

## Six-method spectral comparison (DN, in progress)

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
continuous integral. Full phase gates and merge are pending.

## Release provenance

The implementation agent's exposed model family and mode, goal-accounting
time and token snapshot, and the assumptions behind the API-equivalent cost
estimate are documented in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They describe
the historical v1 implementation and active research runs, not package runtime requirements or benchmark
performance.
