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
recovers it to `3.33e-16` held-out error. The noisy fixture performs a blind
scan over 1,604 frequencies; its eight effective atoms reproduce the
underlying noiseless Green function to `8.56e-9`. The latter is a compact
function representation, not evidence that the four generating frequencies
were uniquely recovered.

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
  --output experiments/data/gecp_census.jsonl
uv run python experiments/exact_surrogates.py
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/data/synthetic_applications.json \
  --plot experiments/data/synthetic_applications.png
```

Runtime algorithms require no network access. The canonical census evaluates
endpoint-refined grids at 128-bit precision and records decimal strings. Its
three tolerance records per cutoff are exact prefixes of one strict trajectory;
the checked dataset is byte-reproducible. Adaptive interval certification is a
separate, more expensive continuous-domain mode.

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

## Release provenance

The implementation agent's exposed model family and mode, goal-accounting
time and token snapshot, and the assumptions behind the API-equivalent cost
estimate are documented in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They describe
the historical v1 implementation and active research runs, not package runtime requirements or benchmark
performance.
