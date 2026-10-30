# Application guide

The [production provider environment guide](environments/README.md) describes
the optional EB–EJ development runtimes. Actual provider availability does not
make the planned domain interfaces part of the public 1.1.0 package below.

[kernelgecp 1.1.0 is available on PyPI](https://pypi.org/project/kernelgecp/1.1.0/):

```bash
python -m pip install kernelgecp==1.1.0
```

Python 3.12/3.13 public installs are verified on Linux, macOS and Windows;
supported Linux x86_64/macOS arm64 wheels include the optional native backend.
Other targets use the universal Python wheel. See the
[package quickstart](PACKAGE_README.md) and [release verification](#versioned-practical-release).

## Reusable physical Green functions

The practical facade uses inverse temperature `beta`, spectral energy cutoff
`wmax`, and physical imaginary time `tau` in the closed interval `[0,beta]`.
The dimensionless research coordinates are `t=tau/beta` and
`Omega=beta*energy`; `beta*wmax` must be at least one. For a positive spectral
measure, the physical fermionic Green function has a minus sign:
`G(tau) = -integral exp(-tau*energy)/(1+exp(-beta*energy)) d rho(energy)`.

```python
import numpy as np
from kernelgecp import FermionicBasis, FermionicFit

beta, energy = 4.0, 0.7
basis = FermionicBasis(beta=beta, wmax=3.0, tol=1e-9)
tau = np.linspace(0, beta, 101)
values = -np.exp(-tau * energy) / (1 + np.exp(-beta * energy))
fit = basis.fit_tau(values, tau)
print(fit.evaluate_tau([0.1, 1.2, 3.9]))
n = np.array([-3, -1, 0, 2])  # ordinary integers, not reduced odd indices
print(fit.evaluate_matsubara(n))
# Convention: integral exp(+i*wn*tau)*G(tau) d tau, wn=(2*n+1)*pi/beta.
# For this single pole the answer is 1/(1j*wn-energy).
fit.save("green.npz")
restored = FermionicFit.load("green.npz")
```

The instance caches one selected-core factorization, bounded evaluation designs
and bounded SVD fit plans. Repeated `fit_tau` calls on identical coordinates and
`rcond` reuse the SVD. `prepare_tau` exposes the reusable `LeastSquaresPlan`;
`design_tau` and `design_matsubara` expose immutable design arrays. The cache
capacity is eight entries per cache by default, not a byte limit; use
`cache_size`, inspect `cache_info`, or release designs/plans with `clear_cache`.
Core factors remain until the basis is released. An instance's mutable caches
are not a concurrent execution API.

Values can be complex and have any number of batch or matrix axes. `axis=0`
is the default observation axis; fitting replaces it with the coefficient
axis, and evaluation replaces that same axis with query coordinates. For
example, values of shape `(ntau, n_orbital, n_orbital)` produce coefficients of
shape `(rank, n_orbital, n_orbital)`. Use `axis=-1` for trailing observations.
An independent `LeastSquaresPlan` also supports complex frequency-domain designs.
SVD truncation uses `s > rcond*smax`, default `eps*max(design.shape)`; rank loss,
untruncated/retained conditioning, arithmetic and training residual scopes are
visible in diagnostics. An arbitrary function can fit its training samples
well and still fail elsewhere: validate independent time and frequency data.

Arithmetic is explicit. `precision="float64"` is a strict fast path;
`precision="reference"` uses `precision_bits=128` by default for selection,
kernel evaluations and cached cardinal solves, then exports doubles. `"auto"`
first tries float64 selection and retries reference only on a numerical
arithmetic breakdown, preserving the failed attempt in metadata. It does not
retry invalid input or a rank-budget failure. Reference precision is supported
from 64 through 4096 bits. Physical node coordinates, exported designs, SVD
fitting and final matrix products are always binary64/complex128. Thus reference
mode is still mixed precision; high-precision construction alone is not an
end-to-end guarantee. `tol` is a finite-grid construction stopping criterion,
not an accuracy contract for arbitrary spectral inputs, unsampled points,
coefficient sensitivity, or Fourier transforms. No facade result is a
continuous or rounding certificate. The pre-existing research API and its
adaptive 128-bit defaults are unchanged.

Basis and fit archives use a versioned, validated NPZ/JSON schema, never pickle.
Loading reconstructs core factors without choosing new pivots and retains
construction diagnostics as unverified historical reports. Training diagnostics
are not restored as fresh evidence. Archives reject incompatible conventions,
objects, malformed numerical payloads and oversized data (64 MiB uncompressed).
Restoration additionally limits rank to 256 and `rank**2 * precision_bits` to
`2**24` before allocating a core, independently of claimed archive settings.
Arrays are immutable owned copies. These are portable numerical artifacts, not
signed scientific certificates; validate accuracy again after transfer.
The runnable [matrix-valued tutorial](examples/reusable_basis.py) checks time
and frequency signs, reuse and persistence on an independent analytic pole.

## Frequency fitting and TRIQS interoperability (unreleased EC source)

These additions target 1.2.0 and currently require this checkout; the immutable
1.1.0 public release does not include them. `prepare_matsubara(indices,
rcond=None)` caches a `LeastSquaresPlan` just as `prepare_tau` does.
`fit_matsubara(values, indices, axis=0, rcond=None)` uses the supplied ordinary
integer indices, preserving order, complex values and every other array axis.

```python
n = np.arange(-24, 24)
observations = 1 / (1j * (2 * n + 1) * np.pi / beta - energy)
frequency_fit = basis.fit_matsubara(observations, n)
assert frequency_fit.observation_domain == "matsubara"
print(frequency_fit.evaluate_tau([0.2, 1.0, 3.8]))
print(frequency_fit.diagnostics.least_squares.retained_condition_number)
```

Coefficients retain the same cardinal-time basis meaning in both fitting
domains; `observation_domain` describes their source observations. Both fit
paths use binary64/complex128 truncated SVD and expose rank/conditioning and
training residual. An underdetermined or truncated fit is a minimum-norm fit,
not an independent accuracy claim. Positive-only frequency observations are
used exactly as supplied: no negative-frequency completion, transpose,
conjugation, moment repair or static-term removal is inferred. Such designs
can be very ill-conditioned and fit their positive observations while giving
poor time/negative-frequency predictions. Validate both domains independently.
The existing finite-query time-fit certificate does not become a frequency-fit
certificate through this API addition.

New source archives use format version 2 and retain the fit observation domain.
Readers accept v1/v2; v1 fits receive the historical default `tau` domain.
Training diagnostics are still not restored as fresh evidence. Loading and
resaving a v1 fit preserves coefficients/predictions and writes v2. Public
1.1.0 readers reject v2, so preserve original files when exchanging with 1.1.0.
The regression fixture was produced by an actual installed public 1.1.0 wheel,
with its wheel/archive hashes recorded in `tests/data`.

The optional `kernelgecp.integrations.triqs` module supplies
`fit_triqs_tau(basis, gf, rcond=None)`, `to_triqs_tau(fit, mesh)` and
`to_triqs_iw(fit, mesh)`. It imports the pinned TRIQS 4.0.2 provider only when
called; ordinary kernelgecp imports need no TRIQS. It reads actual mesh
coordinates, validates fermionic statistics and beta, preserves 0+/beta-
endpoints and native frequency storage order, and returns owning Gf objects.
TRIQS 4's actual modules are `triqs.gfs` and `triqs.mesh`.

```python
from triqs.gfs import BlockGf
from kernelgecp.integrations.triqs import fit_triqs_tau, to_triqs_iw

# source_tau is a BlockGf; iw_mesh is an actual fermionic MeshImFreq.
names = [name for name, _ in source_tau]
fits = [fit_triqs_tau(basis, block) for _, block in source_tau]
converted = BlockGf(
    name_list=names,
    block_list=[to_triqs_iw(fit, iw_mesh) for fit in fits],
    make_copies=True,
)
```

The public unit is a single Gf; this explicit recipe preserves ordered block
names and unequal matrix sizes. The adapter does not infer spectral positivity,
normalization, physical symmetry or high-frequency tails. Consult
[the isolated runtime instructions](environments/TRIQS.md). Actual provider
tests and the installed-wheel smoke gate are separate from optional skips in
the ordinary base environment. No interacting DMFT claim follows from EC.

## DLR/IR interoperability and a physical solver loop

From a checkout with optional providers, run the same complex matrix workload
and small scalar bath self-consistency loop with all three bases:

```bash
uv run --with pydlr==1.0.1 --with sparse-ir==2.1.6 \
  python -m examples.dlr_interoperability
```

The [executable example](examples/dlr_interoperability.py) fits identical time
observations using `LeastSquaresPlan` and provider design matrices, then checks
physical time and Matsubara predictions. It also iterates the half-filled,
noninteracting Bethe equation `G(z)=1/(z-hopping**2*G(z))`, with cached projection
factors, an explicit fixed-point defect and independent iteration budget.
Its reference is the causal semicircle resolvent and separate spectral
quadrature, not another run of the same iteration. This is a small physical
self-consistency workload, not an interacting impurity solver or production
DMFT package. The Bethe bath relation follows the
[TRIQS tutorial](https://triqs.github.io/triqs/latest/userguide/python/dmft_one_page.html);
here the interaction is explicitly zero.

Adapters expose ordinary integer `n` consistently. Sparse-ir's native `uhat`
instead expects the reduced odd index `2*n+1`, as documented in its
[Fourier convention](https://sparse-ir.readthedocs.io/en/latest/basis.html#sparse_ir.FiniteTempBasis.uhat).
Pydlr's native frequency helper has a transpose/conjugation convention for
matrix inputs. The adapter constructs its frequency design from a real identity
and applies coefficients linearly, so general complex batches are not silently
conjugated. These are interoperability adapters, not replacements for the
providers' native optimized fitting APIs.

The separate matched-accuracy application benchmark uses a predeclared native
eps list, shared calibration values in both domains, then a fresh measured
worker for the selected constructor and reuse. Equal eps parameters alone are
not claimed to mean equal accuracy. A clean tracked revision is required:

```bash
uv run --with pydlr==1.0.1 --with sparse-ir==2.1.6 \
  python -m experiments.dlr_application_benchmark \
  --config experiments/dlr_application_benchmark.toml \
  --output experiments/scratch/dlr_application_benchmark.json
```

The report separates construction, design/SVD preparation, cached fitting,
cached/uncached evaluation and solver setup/loop times, plus rank, native-coordinate
conditioning and fresh-worker process peak memory. Memory includes imports,
the harness and native libraries; it is not basis-only allocation. Calibration
has its own worker and cost rather than contaminating the selected worker's
memory measurement. Timing fields are nondeterministic, not a universal speedup
claim. The tutorial's fixed eps is likewise not a calibrated performance result.

The first recorded [cohort A](experiments/data/dlr_applications_cohort_a.json.gz)
retains a real calibration miss: sparse-ir at cutoff 10 passes its calibration
but reaches final frequency error1.21e-8 against target 1e-8. Eleven of twelve
records pass; all three physical solvers pass their separate checks. Two runs
have identical numerical fields and an independent 95-digit physical audit
confirms the failure. Its source/provenance stays frozen, not rewritten to
reflect later improvements. The follow-up protocol uses denser calibration and
a 1/10 calibration guardband while preserving the final1e-8 target and choosing
new untouched final points. A guardband is not a uniform error theorem.

`scripts/audit_dlr_applications.py` is a separate95-digit verifier: it imports
no benchmark truth helper, checks immutable generation-era Git blobs and
container/payload hashes, recomputes all forward physical errors and checks
solver frequencies plus selected direct energy integrals. Its CLI requires the
expected full source SHA and artifact SHA256, and accepts JSON or compressed
JSON. It distinguishes a correctly reported scientific failure from an invalid
artifact. Exact environment/provenance checking requires the recorded optional
versions in `experiments/dlr_benchmark_requirements.txt` and the recorded platform;
it does not claim a new platform reproduces historical performance. Gzip only
reduces storage size and does not alter the underlying generated JSON payload.

### Recorded DV matched cohort B

The [cohort B artifact](experiments/data/dlr_applications_cohort_b.json.gz) is
generated from frozen source `8a4dd6bec0d606764667a1721e06469286e82cef`.
All 12 method/case records pass the shared final1e-8 criterion, with all 36 complex
matrix forward workloads independently checked in time and Matsubara frequency.
Maximum forward errors are 5.45e-10 (time) and 7.97e-10 (frequency). All three
Bethe loops converge and pass: maximum time error5.74e-10, unseen frequency
error4.23e-11, and independently checked observed-grid error4.49e-10. This is
finite sampled evidence for the fixed physical fixtures, not a uniform bound.

| Dimensionless cutoff | GECP native eps / rank | pydlr native eps / rank | sparse-ir native eps / rank |
| --- | --- | --- | --- |
| 10 | 1e-8 / 14 | 1e-8 / 13 | 1e-10 / 14 |
| 100 | 1e-8 / 28 | 1e-10 / 31 | 1e-10 / 30 |
| 1,000 | 1e-10 / 55 | 1e-10 / 51 | 1e-10 / 52 |
| 20 (solver case) | 1e-8 / 17 | 1e-10 / 19 | 1e-10 / 18 |

These are the first passing choices on a coarse predeclared candidate list,
not proved minimal ranks or optimally tuned runtimes. GECP uses128-bit
construction/application at cutoff 1,000 after recorded float64 breakdown;
other GECP cases use53-bit arithmetic. All designs/fit products export doubles.
Provider construction arithmetic is not inferred from its exported dtype.

In the first recorded run, cutoff 1,000 constructor times are 19.67s(GECP),
0.0194s(pydlr), and 1.337s(sparse-ir). Fresh selected-worker high-water RSS is
approximately98.2/91.9/178.0MB, respectively, including harness/import/native
state before output JSON encoding. The complete raw measurements, setup and
reuse timings and scopes remain in the artifact. This is a single-machine
application comparison, not a universal speedup claim; GECP's reference
construction is an obvious performance target rather than an established
advantage. Mature native numerical routines are already present beneath the
Python APIs. Profile the actual workload before deciding what Rust can improve.

Both complete runs have byte-identical numerical fields after removing all
measurements. The [independent B audit](experiments/data/dlr_applications_cohort_b_audit.json)
checks43frozen source blobs,838 array hashes, all forward errors, all solver
frequencies/defects and 18 selected95-digit direct time energy integrals. Maximum
stored forward truth discrepancy is 1.31e-16 and reported-error discrepancy
5.95e-17. [A's audit](experiments/data/dlr_applications_cohort_a_audit.json) remains
alongside its failed record. Reproduction requires the recorded generation
revision/configuration and optional environment; a current-revision rerun
properly produces different provenance. Source, payload and repeat hashes are
centralized in `RESEARCH.md`.

The separately tested [standalone Sparse-IR recipe](examples/sparse_ir_physical_recipe.py)
requires no kernelgecp import and uses native sampling. Following a
[maintainer invitation](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6405100346),
[upstream PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10) adds three
executable cells to its sparse-sampling notebook. With the upstream lock
(Python 3.13.2, NumPy 2.3.5, sparse-ir 2.1.4, pylibsparseir 0.9.0), isolated
and full-notebook execution agree: rank 14, sampled time/frequency errors
`2.76170e-11` / `4.94081e-11` against `1e-8`. The full book builds successfully;
independent Fourier quadrature and content/convention review pass. The
[contribution record](UPSTREAM_PROPOSAL.md) preserves exact source/hash,
validation scope and the corrected local harness failure. Submission-time
upstream CI required permission and ran no jobs; the live PR owns later
checks, review and merge status. These local results do not establish upstream
acceptance or change the separate package-publication evidence.

Sampling matters. An explicitly retained pilot removes important low Matsubara
rows and reaches a tiny sampled defect while predicting unseen values poorly.
Reality/symmetry constraints and extra iterations did not repair that missing
information. The final solver preserves its low-frequency observations and
checks independent truth on the complete solved grid (not heldout), plus disjoint
middle/high-frequency and time holdouts. Forward fitting separately has unseen
low-frequency checks. Convergence, each accuracy check and certification remain
different fields. All accuracy evidence is sampled; `tol`, singular values and
an observed defect do not certify continuous or roundoff error.

## Research interfaces

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

## Column-oracle covariance factors (EI implementation)

`KernelCholeskyFactor.from_columns(diagonal, column, max_rank=256, atol=1e-10,
rtol=0)` approximates a finite real symmetric PSD covariance from its diagonal
and one full length-n callback column per selected pivot. It stores an owning,
immutable float64 root `factor` with shape `(n, rank)` and the original and
residual diagonals; it never constructs a full covariance or residual matrix.
The callback is not retained. The existing dense `pivoted_cholesky` is unchanged.

```python
import numpy as np

from kernelgecp import KernelCholeskyFactor

points = np.linspace(-2, 2, 101)
factor = KernelCholeskyFactor.from_columns(
    np.ones(points.size),
    lambda j: np.exp(-((points - points[j]) ** 2)),
    max_rank=20,
)
system = factor.shifted(0.01)
solution = system.solve_with_diagnostics(np.sin(points))
print(factor.rank, factor.stop_reason, solution.relative_residual)
```

`factor.indices` gives the pivot order. `extend(cross_columns, diagonal=None)`
expects only `K[new, selected pivots]`, with exactly that column order, and
returns immutable `features` plus an optional checked `residual_diagonal`.
It uses a triangular solve with the selected root. An approximate cross
covariance is `extension.features @ factor.factor.T`; no training/test
covariance is built implicitly. Supply a new diagonal to check the extension's
PSD residual constraint. A caller-supplied kernel remains responsible for
consistent physical units, row order and hyperparameters.

`as_linear_operator(shift=0)` implements `L @ L.T + diag(shift)` as a SciPy
operator. Positive scalar/vector noise is prepared once with `shifted(shift)`;
its `solve`, `solve_with_diagnostics` and `as_preconditioner` reuse a small
rank-by-rank Cholesky system. The preconditioner is passed as `M` to SciPy CG
on the original kernel operator. A successful approximate solve does not
establish the original-kernel CG residual. Woodbury refinement is checked
against the approximate system; a temporary implicit Householder QR/SVD
fallback handles subtraction loss for tiny noise without forming a complete
orthogonal matrix. Unresolved or unrepresentable solves raise
`NumericalBreakdownError`, including loss during final RHS rescaling.

`to_gpytorch(shift=0)` copies the identical frozen root into a CPU float64
`linear_operator` object usable by GPyTorch. It imports optional providers
only on that call and owns its tensor storage. Pivot selection is outside
autograd. The supported runtime starts with EB's pinned torch 2.14.1,
GPyTorch 1.15.2 and linear_operator 0.6.1; no GPU or alternate runtime is
claimed by the small actual-provider check.

`factor.diagnostics` and its fresh `metadata` snapshot report convergence,
stop reason, selected-root SVD conditioning, tie counts, roundoff slack and
actual float64 precision. `storage_nbytes` counts retained numerical arrays
and scalar histories, excluding Python overhead; the shifted object's count
includes the shared factor once. Global symmetry/PSD is an explicit callback
assumption: only requested columns are checked. Diagonal convergence and
small synthetic prediction errors are numerical evidence, not a continuous
certificate. The [frozen attributed UCI campaign](COVARIANCE.md) passes its
18 local trials, actual Mac provider and independent artifact audit. Linux
provider CI, final native/root and delivery gates remain pending under
PRODUCTION_INTEGRATIONS.md; numerical acceptance alone does not close EI.

The [real-data campaign guide](experiments/COVARIANCE.md) defines attributed
data hashes, calibration isolation and independent receipts. Its predictive
standard deviation includes the declared observation-noise variance and is
normalized by the training target standard deviation. The same convention is
used for dense and approximate comparisons; the acceptance thresholds are
unchanged. This is an assumed Gaussian observation model, not a calibrated
environmental forecasting or uncertainty-coverage claim.

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

## DT matched application benchmark

The first practical-application extension compares GECP, `pydlr` and `sparse-ir`
bases for the same physical negative fermionic Green functions. Three small
normalized atomic spectra are closed-form forward workloads, not support-recovery
or material-specific results. All bases receive identical endpoint-resolved
training times/observations and fixed disjoint midpoint holdouts. Physical
conventions are `Lambda=beta*wmax`, `t=tau/beta` and `Omega=beta*omega`.

`experiments.dlr_benchmark` uses one cached SVD least-squares protocol for every
basis. It separates provider import, actual adapter construction (including its
package defaults), training-design evaluation, fit-plan factorization, repeated
fitting, uncached evaluation and cached-grid evaluation. It does not compare the
libraries' native fitting algorithms; those are separately smoke-tested. Native
tolerance meanings remain distinct and every method faces the same sampled
absolute-error acceptance test. Raw design conditioning depends on coefficient
coordinates and is not a basis-independent quality ranking.

Each case/method runs in a fresh worker with six disclosed thread settings fixed
to one. Peak process RSS includes runtime, imports and native libraries, not just
algorithm arrays. Timings/RSS are machine-dependent pilot measurements, separated
from numerical records; no byte-repeatable timing, speedup or certification claim
is made. Unavailable methods, numerical failures, timeouts and inaccurate fits
are retained. Workers must match frozen parent source/configuration and input
fingerprints before an output is accepted. Use a fresh CLI from clean tracked
source; imported API functions cannot attest that stale code matches disk.

Install optional packages without changing the core dependency lock:

```bash
uv run --isolated \
  --with-requirements experiments/dlr_benchmark_requirements.txt \
  python -m experiments.dlr_benchmark \
  --config experiments/dlr_benchmark.toml \
  --output experiments/scratch/dlr_benchmark_pilot.json
```

The public package, default precision/pivot behavior and original research
artifacts are unchanged. Rust optimization, a public reusable-basis API,
Matsubara round trips, noisy-data fitting and physical solver integration are
subsequent targets, not features delivered by this benchmark harness.

### Recorded DT pilot

The [pilot artifact](experiments/data/dlr_benchmark_pilot.json) was generated
from clean source `c0dda11ba91ad32f011eb7364a2c9ec2f9622d8e` with the pinned
optional requirements. It contains nine fresh-worker cases: three methods at
cutoffs 10, 100 and 1,000. Eight meet the common `1e-8` heldout target for all
three workloads; float64 GECP at cutoff 1,000 returns an explicit arithmetic-
range failure. Successful ranks are:

| Cutoff | GECP | pydlr | sparse-ir |
| --- | --- | --- | --- |
| 10 | 13 | 13 | 12 |
| 100 | 27 | 26 | 26 |
| 1,000 | numerical-range failure | 43 | 43 |

The largest successful sampled error is `4.3812592443e-9`. Independent literal
100-digit pole sums check all 24 successful workloads/3,498 stored prediction
points: maximum stored-truth discrepancy `1.40e-16`, error-report discrepancy
`5.97e-17`, and all pass flags agree. Two fresh full runs match numerical and
environment/provenance fields exactly after excluding `measurements`; their
timings and memory are not expected to repeat. Artifact SHA-256 is
`70ce453dddba309b488ec72c7740a8cd122ba6a6cc22dbb3ff0ecc152cfaa7ab`.

Constructor measurements include each adapter's native initialization features;
they are single-machine pilot samples, not a statistically established speed
ranking. Provider-import times are incremental after the common harness, not
matched cold-start measurements. The six pinned optional distributions and
their installed transitive environment pass a dependency audit. The immediate
decision is to investigate the float64 range failure and design a reusable
cached basis with an explicit precision policy before committing to Rust.

Current phase/root/CI/delivery evidence is maintained in
[PR #126](https://github.com/tripp-smith/gecp-kernel-structure/pull/126), separate
from the frozen artifact's generation source and original-goal accounting.

## Release provenance

DR's documentation source `51fa6b3` passed 759 tests in 144.90 seconds,
CI `40054146369`/`40054160333`, clean source/wheel inspection and both
outside-checkout installed environments. PR #124 is merged. These are later
documentation checks, not relabeled DQ experiments or new runtime behavior.
Exploratory commands above use ignored scratch outputs; committed artifacts
retain their original generation manifests. Local build snapshots are not a
new package publication.

The implementation agent's exposed model family and mode, goal-accounting
time and token snapshot, and the assumptions behind the API-equivalent cost
estimate are documented in
[FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). They describe
the historical v1 implementation, completed original-SPEC research goal and
later documentation-only task boundary, not package runtime requirements,
phase-only costs or benchmark performance.
## Finite-query numerical certificates

`certify_fit_precision(fit, values, tau, query_tau=..., matsubara_indices=...)`
independently encloses the difference between an actual `FermionicFit` and the
unique full-column-rank least-squares fit on the supplied time samples. Every
stored binary64 input denotes its exact binary rational value; the ideal kernel
uses exact physical `tau/beta`. This encloses normalization, kernel/core solve,
coefficient fit, design export and final output arithmetic effects. It does not
assume that a working precision or a small backward residual guarantees a small
forward error.

```python
from kernelgecp import certify_fit_precision, certify_grid_pivot_precision

report = certify_fit_precision(
    fit, values, tau, query_tau=queries, matsubara_indices=[-1, 0, 4]
)
if report.certified:
    assert report.coefficients is not None
    assert report.time is not None and report.matsubara is not None
    print(report.coefficients.max_error)
    print(report.time.max_error, report.matsubara.max_error)
else:
    print(report.termination_reason)

pivot = certify_grid_pivot_precision(basis, queries, [-2.3, -0.3, 1.7])
print(pivot.certified, pivot.eta, pivot.objective_scope)
```

An enclosed report is not an accuracy-tolerance pass. A truncated SVD fit can
differ greatly from the ideal full-rank fit, especially when observations omit
endpoint information. Compare the reported budget to the application's error
allowance. Physical spectral-model and observation uncertainty are separate;
the exact-data least-squares reference need not be a good physical model.
Supplying `tau=None` uses the selected physical nodes, matching `fit_tau`.
The observation axis is `fit.axis`; complex matrix batches are preserved.

`report.core` distinguishes exact physical versus rounded normalization and a
reconstructed binary64 kernel proxy. The proxy is not the internal reference
factorization; its error already includes normalization and is not an additional
disjoint stage. Training/query design bounds compare the actual facade exports
to independent interval core solves and include true-core backward residuals.
Coefficient enclosures use interval Gram solves solely for certification; the
ordinary public fitter still uses its cached SVD. Query reports contain a direct
prediction enclosure and the separately valid coefficient/design/final-dot
triangle budget. The final BLAS result is compared to an exact rational dot of
the actual exported design and computed coefficients, so no vendor-specific
reduction-order assumption is needed.

Reports and bound arrays are immutable and hash all numerical inputs plus the
actual designs/outputs. Default interval precision is 256 bits, with binary64
outward bound exports; rank/work budgets are explicit. Zero-containing interval
pivots, underdetermined data, export overflow and budget exhaustion produce
uncertified reports with missing bounds, not zero error or a singularity proof.
Higher precision can resolve interval overestimation but cannot repair an
ill-conditioned exact fitting problem. Successful pivot gates establish
regularity, not a rigorous norm condition number.

The grid certificate bounds the discrepancy between the actual rounded selected
cross and the exact physical selected cross on exactly the supplied Cartesian
grid. It transfers the first row-major rounded maximum into a positive true
pivot quality when the discrepancy is small enough; ties are explicit. It
does not cover points between grid samples or identify an inexact trajectory
with an exact complete-pivot run. An exactly zero enclosed grid is labeled
separately and has trivial quality one.

The conditional pivot-transfer and complex finite-dot budget theorems are
kernel-checked Lean results. The interval LU/Gram implementation trusts the
mpmath interval engine and is independently checked with rational/Cramer,
high-precision QR and bordered-determinant formulations; the numerical program
itself has not been verified in Lean. mpmath documents interval containment and
the experimental scope of its engine in its [context documentation](https://mpmath.org/doc/current/contexts.html).

Run `uv run python -m examples.precision_bridge` for a complex matrix example
with an independent QR audit. The frozen experiment uses
`experiments/precision_bridge.toml` and the module
`experiments.precision_bridge`; its source must be committed and tracked-clean
before generation. Final canonical results and phase delivery are recorded in
the current README receipt.

DW's canonical six-case artifact is
`experiments/data/precision_bridge.json.gz`, generated twice from `a4995b4`
with identical complete JSON bytes. Five enclosures complete and one remains
inconclusive; a completed enclosure can have a large valid discrepancy.
At cutoff 1000 the endpoint-resolved numerical budget is below `2.50e-14`,
but its sampled physical error reaches `3.38e-8`. Selected-node fitting gives
a numerical budget below `9.65e-16` and physical error below `4.48e-10`.
The uniform-time truncated fit has a frequency discrepancy around `1.10e15`
from the ill-conditioned exact full-rank target. These outcomes are retained;
they illustrate why an arithmetic certificate and physical accuracy need
separate acceptance criteria. See `RESEARCH.md` for source/container/payload
hashes and [PR #129](https://github.com/tripp-smith/gecp-kernel-structure/pull/129)
for final revision, platform/install and authorized delivery checks.

Recheck the immutable data independently, without importing the implementation:

```bash
uv run python scripts/audit_precision_bridge.py \
  --source experiments/data/precision_bridge.json.gz \
  --output experiments/scratch/precision_bridge_audit.json \
  --generation-commit a4995b471506fe043f3102114d42620c3c645d76 \
  --expected-sha256 b88c47735487433a482620a60577485cb1749e652ccfde201db2145112f9d09a \
  --expected-payload-sha256 b77e38e42aad25082d14fbafb9d767525f61c302942b7ecca607c54df5631782
```

This checks the generation-era Git/config/dependency manifests, complete arrays,
all grid points by 1280-bit bordered determinants, exact selected-node
coefficients, selected-node predictions and all physical pole checks. It
preserves inconclusive/large-error cases and states which other LS checks were
performed by the generation-time QR audit rather than this separate tool.
Schema 2 stores all 90 complete numerical arrays, including the actual rounded
grid residuals and binary64 core proxies. Auditing checks these observed values
directly and does not require a different machine to reproduce BLAS or exp
rounding bit for bit. A separate small guard checks physical fixture identity;
the certificate error bounds retain their strict checks. Historical cohort B
and its generation-platform audit remain in the `precision_bridge_cohort_b*`
files. CI fetches the Git history required by the provenance test.
## Optional Rust reference time backend

The physical facade accepts `reference_backend="rust"` when a compatible native
wheel is installed. The default remains `"python"` and the universal wheel
needs no Rust compiler or native module:

```python
from kernelgecp import FermionicBasis, FermionicFit

basis = FermionicBasis(
    beta=20.0,
    wmax=50.0,
    tol=1e-10,
    precision="reference",
    precision_bits=128,
    reference_backend="rust",
)
print(basis.diagnostics.time_evaluation_backend)
print(basis.diagnostics.native_backend_info)
```

This backend evaluates reference time-kernel rows and solves with the existing
cached LU using 128-bit MPFR arithmetic. Physical normalization stays binary64
in Python. GECP selection and Matsubara evaluation remain mpmath; fits and
final products remain NumPy/SciPy. Diagnostics identify every stage. Float64
construction uses the existing SciPy path, so a reference-backend request does
not make that path native. Other reference precisions and ranks above the
native limit of 128 are not supported by this backend.

A missing/incompatible module or native arithmetic failure raises a clear
error; it does not silently substitute a different implementation or precision.
The checked native handle owns its factors and releases the GIL during work.
The public basis caches themselves retain their existing concurrency contract.
No new rounding or continuous physical guarantee follows from choosing Rust;
the explicit finite-query certificate can audit actual native output too.

Archives remain portable. `FermionicBasis.load(path)` and
`FermionicFit.load(path)` restore Python execution by default, regardless of
saved backend history. Explicit `load(path, reference_backend="rust")` rebuilds
a validated native handle. No native objects or executable state are loaded
from an archive.

The measured opportunity is repeated evaluation on changing time meshes after
construction. DX's frozen three-trial cohort gives median application time
42.80 s Python versus 21.33 s Rust for rank 55, construction/fit/setup and 80
distinct 64-row evaluations with cache capacity eight. Full-process reduction
is 47.76%; the eight-query prefix improves only 8.92%, and warm cache hits show
little difference. All physical stream errors are below `2.57e-11` against
`1e-8`; paired differences are below `4.26e-15`, not bitwise equal. This is
single-machine sampled evidence, not a universal speed or accuracy guarantee.
The complete artifact is `experiments/data/rust_reference_benchmark.json.gz`;
run `python -m experiments.rust_reference_benchmark --help` for the frozen
configuration runner. Existing research GECP defaults are unchanged.

Official native wheels are built from the inspected corresponding-source asset
using `scripts/build_native_wheel.py`. They carry complete dependency notices,
an explicit composite SPDX expression, a sanitized dependency inventory and
the matching source-asset identity. The Python/source distributions retain the
project's Apache license. Native distribution metadata is specific to the code
it contains, consistent with the [PyPA metadata specification](https://packaging.python.org/en/latest/specifications/core-metadata/#license-expression).
The complete locked third-party sources and offline rebuild instructions are
prepared with `scripts/prepare_corresponding_source.py`; see
`native/SOURCE_BUILD.md`. DX's actual native build/install coverage is Linux
x86_64 (glibc 2.28+) and macOS arm64 (11+) on Python 3.12/3.13. Windows, Intel
macOS, Linux ARM and other unbuilt native targets use the universal Python
fallback; no native support is inferred from that fallback.

## Interacting Bethe IPT (unreleased, targeting 1.2.0)

`kernelgecp.dmft` supplies a deterministic scalar half-filled paramagnetic
Bethe IPT loop. It cancels the Hartree term analytically and fits the dynamic
self-energy separately. IPT is an approximate impurity solver; compression
accuracy is checked against the same IPT model rather than exact Hubbard
physics. The canonical supported cases are t=1, beta=10 and U=0/1/2.

```python
from kernelgecp import FermionicBasis
from kernelgecp.dmft import IPTConfig, IPTPlan

cfg = IPTConfig(interaction=2.0)
green = FermionicBasis(10, 4, 1e-8, precision="reference", reference_backend="python")
sigma = FermionicBasis(10, 12, 1e-8, precision="reference", reference_backend="python")
plan = IPTPlan(green, sigma, config=cfg)
result = plan.solve()
assert result.converged, (result.stop_reason, result.failure_detail)
```

Run `uv run python -m examples.interacting_dmft`. `IPTPlan.update(bath_iw)`
exposes a single unmixed update. `solve(initial_bath=None,
record_trajectory=False)` returns immutable arrays, coefficients, projection
and SVD diagnostics, iteration history, and a stop reason. A numerical failure
retains the last valid candidate and failed bath; first-update failure has no
substitute numerical solution. Iteration-budget exhaustion is unconverged.

The reflected diagram is `U**2 * G0(tau)**2 * G0(beta - tau)`. Its cube is
equivalent only under the declared particle-hole symmetry. This is not a
doped IPT solver. Green and self-energy bases must share beta, with the latter
cutoff at least three times the former. That finite-basis requirement is not
a compact-support theorem about the infinite interacting spectrum. Designs
and SVDs are cached; fitted complex arithmetic is binary64 even when GECP
construction/evaluation uses reference precision. Generic bases implementing
`PhysicalBasis` must retain stable parameters throughout a plan's lifetime.

The stop criterion is an unmixed Bethe equation defect below tolerance for
two consecutive updates. `frequency_fit_indices` can select observations;
omitted modes influence neither the SVD nor stopping. Full-grid defects remain
separate diagnostics. `g_iw` is the raw Dyson candidate; `g_tau` and
`g_coefficients` describe the fitted representation. Predict unseen modes from
`green.design_matsubara(indices) @ result.g_coefficients`. The returned bath is
the actual input that produced the final candidate, not an unused mixed step.

`kernelgecp.integrations.triqs_dmft.DenseIPTPlan` is the optional native
TRIQS 4.0.2 Fourier reference. It applies explicit half-filled model moments;
the compressed loop does not repair moments or impose symmetry. Native
interpolation requires independent grid refinement. `kernelgecp.lehmann`
contains bounded exact atom/one-bath Fock references and finite-pole IPT
triples with independent matrix-exponential/Fourier tests.

The optional source campaign is `python -m experiments.interacting_dmft
--output-dir <new-directory>` inside the pinned TRIQS environment with
`environments/locks/dlr-physics-pip.txt` installed. Its TOML freezes all
observations, cutoffs, native tolerances, reference refinement and three fresh
process repetitions. Pilots use `--pilot` and cannot pass canonical acceptance.
Calibration, independent closed-loop validation and common-trajectory replay
are distinct records. Final frequencies lie outside the primary mesh; final
time queries use an off-grid Gauss rule. Interpolated points are evaluation
queries, not additional observations. Run the separate artifact auditor before
interpreting a campaign. Tail-window fits are diagnostics without interacting
remainder certificates. Loop RSS excludes later validation/replay; full-worker
wall/RSS includes them. Coefficient and output-array sizes have distinct scopes
and do not imply a total-memory or universal performance advantage.

## Covariance-aware MaxEnt continuation (1.2 development)

`kernelgecp.integrations.maxent` is an optional, lazy TRIQS MaxEnt 4.0.1
adapter for real scalar or diagonal positive spectral functions. Ordinary
`import kernelgecp` requires no TRIQS installation. Calling the provider fails
clearly when it is unavailable; no sparse-recovery surrogate is substituted.

`SpectralDataset(beta, tau, values, covariance)` owns immutable observations
and their complete covariance. It preserves the distinct zero-plus and
beta-minus endpoints and rejects implicit conversion of complex observations.
Covariance eigenmode selection has declared relative/absolute cutoffs, retained
rank and discarded-mode diagnostics. `diagonal_spectral_datasets` extracts
declared marginals from a time-major joint covariance; it does not construct a
joint posterior or support off-diagonal spectral continuation.

Use `prepare_spectral_projection(basis, source_tau, target_tau, rcond=...)`
and `projection.apply(dataset)` for compression. The result carries the full
linear observation map P and the original source mesh. Continuation sends
P G, P C P-transpose and P K to the actual `TauMaxEnt` provider before its
covariance rotation. New interpolated points therefore remain correlated.
For a same-mesh projection, retained SVD vectors form U-r U-r-transpose
directly, avoiding the cancellation in ill-scaled coefficient maps.

`continue_tau_maxent(dataset, MaxEntConfig(...))` retains each alpha trial,
its convergence/iteration status and chi-square. The selected converged alpha
minimizes the declared discrepancy to retained noise rank. Explicit optional
mass normalization retains both the raw mass and the recomputed normalized
training defect. With normalization enabled, alpha selection compares the
normalized physical forward model against training data, including the complete
projection and retained covariance modes. Provider chi-square values remain
separate; diagnostics also retain each trial mass and selection discrepancy.
Without normalization, selection uses the provider chi-square. Normalization is
explicit postprocessing, not a constraint inside native optimization. `MaxEntSpectrum.evaluate_tau` performs physical forward
evaluation; it is separate from the projected observation model.
`save_spectral_dataset`/`load_spectral_dataset` preserve full covariance,
projection, source mesh and a versioned hash-bound application manifest with
caller-declared units/provider/state ordering/settings.

The [synthetic campaign](SPECTROSCOPY.md) uses the actual provider, independent
held-out forward checks and a separate arithmetic audit. Measured CT-HYB
validation and full phase delivery remain required. Bootstrap variation and
default-model sensitivity do not establish unique peak recovery or complete
uncertainty coverage. Existing sparse-recovery experiments remain separately
labeled comparisons.

Provider failures retain the native exception type, message and completed
trial history. No partial alpha family substitutes for a successful result,
and continuation restores the caller's NumPy error policy even on failure.
The experimental Linux-only fixed-bath acquisition is a separate research
workflow, not a replacement DMFT solver or a published 1.1.0 interface. It
requires a committed source freeze and passed independent calibration before
new measurements. Its 256-chain acceptance and the measured spectroscopy
phase remain open; see the [current evidence](SPECTROSCOPY.md).

## EPW anisotropic frequency interoperability (1.2 development)

`kernelgecp.integrations.epw` reads the tested qe-7.5 / EPW 6.0 anisotropic
imaginary-axis format. An explicit application context supplies temperature,
eV units, ordinal state ordering, reference Fermi energy and chemical-potential
offset. Native printed frequencies and state ordinals are preserved. A changed
native chemical potential requires its measured offset; it is never inferred
from the fitted Green function. Lossless parser/writer replay preserves rounding.

`kernelgecp.integrations.epw_frequency.fit_epw_frequency_cross(data,
real_nambu_parity=True, fit_rows=..., rank=60)` fits the actual channels
`[Re G11, Im G11/omega, Re G12]` using a shared Cauchy cross at `omega**2`.
The caller explicitly requests the real Eliashberg Nambu parity needed to
reconstruct both diagonal and anomalous entries. This real-linear interface
has a distinct representation schema from the generic complex-linear
`FermionicFit`; it cannot represent arbitrary complex Nambu matrices.

The returned immutable `EPWFrequencyCrossFit` evaluates ordinary signed
Matsubara indices, reconstructs native Z/gap/energy-shift fields, and saves or
loads a bounded numeric NPZ plus a hash-bound JSON manifest. The archive retains
selected nodes, poles, binary64 cardinal weights, coefficients, native layout,
channel units and numerical diagnostics. Its shared version-2 transport leaves
the existing version-1 and version-2 `FermionicFit` readers unchanged. Loading
reuses the saved basis rather than selecting new pivots. Evaluation and native
field recovery reject unsafe conditioning instead of substituting another fit.

The 13 eV pole dictionary is a numerical approximation parameter, not a
certified spectral cutoff. No tail correction or self-energy subtraction is
implicit. The historical fit-only workflow is `experiments.epw_frequency` with
`experiments.audit_epw_frequency`; the supported archive workflow is
`experiments.epw_restart_roundtrip`, `experiments.audit_epw_restart_roundtrip`
and `scripts.providers.epw_restart`, using
the isolated locked `environments/epw-science` runtime. Frozen held-out
Green/field accuracy passes. The supported restart path explicitly combines the
model with lossless native XOR corrections and exports the restored fields.
Unmodified EPW passes the 20 K raw/restored comparison; this does not establish
a pure fitted restart, which failed at 17 K under unchanged budgets. Coarse MgB2 cases provide
interoperability evidence, not converged experimental materials predictions or
EPW solver acceleration. EPW's specialized IR-object format is unsupported.

## Published release evidence

Version 1.1.0 is public on [PyPI](https://pypi.org/project/kernelgecp/1.1.0/) and
[GitHub](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0),
with a short [package quickstart](PACKAGE_README.md),
[computational methods](METHODS.md), [citation metadata](CITATION.cff) and
[release notes](RELEASE_NOTES.md). The immutable source/tag is
`9a679930db326bda2456bfa91b3256f37359e4ce`. The release workflow promoted the
unchanged exact-main tested universal wheel, two supported native wheels and
sdist; complete native sources, notices, dependency inventory and manifest were
public on GitHub before the Trusted Publisher's OIDC upload on 2026-10-25 at
23:23:48–23:23:53 UTC.

[Workflow 40244470229, attempt 2](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40244470229/attempts/2)
passes all ten jobs, including clean public-index wheel installations on
Linux/macOS/Windows with Python 3.12/3.13, actual native execution where built,
the Windows universal fallback and a separate public-sdist installation.
Public PyPI JSON and simple-index file hashes match the four tested files in
the [release manifest](https://github.com/tripp-smith/gecp-kernel-structure/releases/download/v1.1.0/release-manifest.json).
The [publication receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/131#issuecomment-6404917491)
preserves the first attempt: five postchecks exhausted a short PyPI-JSON HTTP
404 readiness window after successful upload. Only those five verification
jobs were rerun once the exact public hashes were visible; no upload or source
job was repeated. These release checks leave the finite-fixture, conditional
certificate and exact-theorem scopes above unchanged. Upstream outreach now
includes the maintainer-invited [notebook PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10);
its live page owns subsequent review, check and merge status.

## Optional actual CT-HYB measured-data workflow

The source module `kernelgecp.integrations.cthyb_dmft` requires the isolated
TRIQS 4.0.2/CT-HYB 4.0.1 runtime; importing kernelgecp does not import either
provider. `CTHYBConfig` declares a bounded single-band paramagnetic Bethe
problem. `CTHYBPlan(config, basis, fit_tau_indices=...)` accepts the same physical
time/frequency design protocol as IPTPlan. `sample(bath_iw, cycles=..., seed=...)`
creates and thermalizes a fresh actual native chain and returns immutable
spin-resolved time/frequency data, input bath, native density-matrix Green/Sigma
moments, occupations, double occupancy, sign and solve status.

`represent(measurement)` fits only the declared measured time rows, evaluates
G(iw) and applies the same native tail-window replacement as the dense path.
It reports retained rank, conditioning, training/excluded-row residuals and
native Fourier discrepancy. `bath_update(measurement, representation)` uses
the represented frequency values in the mixed Bethe update; omitting the
representation selects native Fourier values. Spin averaging applies only to
the declared paramagnetic bath. Static mu=U/2 remains explicit in native G0;
no finite-Monte-Carlo density or particle-hole symmetry is silently imposed.

`observables` returns total density, double occupancy, and real/imaginary
spin-mean Green values at ordinary signed integer indices. Pass the optional
representation to inspect the actual compressed output. `estimate_replicas`
uses sample covariance across at least two independent fixed-bath chain means;
the caller must establish independence with noncolliding RNG states and fresh chains. In the pinned native default RNG, even seed `2k` and odd seed `2k+1` map to the same state; distinct seed labels alone are insufficient. Use distinct odd seeds and independently reject repeated raw measurements.
Time-bin scatter and loop-iteration scatter are not standard-error estimates.
The returned covariance is the sample covariance of chain means; covariance
of the final ensemble average is that matrix divided by the replica count.
The independent saved-chain audit stores both explicitly for later EH use.

The source campaign `python -m experiments.cthyb_dmft all --output PATH` uses
the frozen protocol in `experiments/cthyb_dmft.toml`; the fresh scientific
campaign uses `experiments/cthyb_fresh_odd.toml`. Calibration/final seeds
and Green query sets are disjoint. It retains hash-bound numeric archives,
refuses changed source/config/runtime on resume, freezes the required cycle
ladder before final data, and returns nonzero for inconclusive uncertainty,
representation, self-consistency or refinement. The fast actual-provider CI
probe proves wiring, not canonical statistical/physical acceptance. No CT-HYB
speed or exact interacting-model accuracy claim follows from these interfaces.

## Molecular Green GF2 checkpoints

The 1.2.0 candidate source contains an optional molecular reader for
Green/green-mbtools 1.0.0 restricted GF2 output in a declared symmetric Lowdin
basis. Published 1.1.0 files remain unchanged. Install the candidate checkout
with `python -m pip install '.[green]'` and import it explicitly:

```python
from kernelgecp.integrations.green import read_green_tau

data = read_green_tau("input.h5", "sim.h5", "bound-manifest.json")
assert data.g_tau.shape[1:] == data.sigma_tau.shape[1:]
# Axes: tau, spin, molecular k=1, orbital row, orbital column.
print(data.iteration, data.mu, data.energies.total_hartree)
```

The required `kernelgecp.green-molecular-run.v1` JSON manifest binds the actual
HDF5 files, native convergence log, source/build provenance and provider receipt
by SHA256. It declares units, geometry and method settings that Green's output
does not store. For retained native provider outputs, run
`python scripts/providers/bind_green_manifest.py NATIVE_DIRECTORY PROVIDER_RECEIPT`
to enrich their original explicit run settings. It refuses files absent from
the source provider receipt. See the attributed fixture in
`tests/data/green_h2_074` and the [official output contract](https://green-phys.org/docs/components/data-formats/output-files/).

The reader validates versions, dimensions, finite values, overlap, native log
and independently recomputed final iteration-pair deltas. It returns owned,
immutable arrays without collapsing matrix/spin axes, unpickling `Cell`,
conjugating data or repairing symmetry. `/iter` must identify a complete
checkpoint matching the converged manifest; incomplete output never triggers
an implicit older fallback. Explicit historical iterations retain their source
identity and have `converged=False`.

`data.sigma_static`, `data.h_core`, `data.overlap` and `data.mu` remain explicit.
Only decaying `g_tau`/dynamic `sigma_tau` belong in the fermionic basis. The
dynamic time-domain self-energy has Hartree-squared units; energy and static
self-energy use Hartree, imaginary time inverse Hartree, and geometry Angstrom.
The total energy is `Energy_HF + Energy_2b`; `Energy_1b` must not be added again.
The matrix archive recipe is `python -m examples.molecular_green --help`.

EF scientific acceptance remains in progress. Frozen comparison `56c3211`
passes the declared same-solution sampled compression gates across both
geometries and grids, with independent QR/Fourier audit; exact-source CI remains pending.
This is not a continuous-domain or GF2 model-accuracy guarantee. The separate scientific environment
and native commands are documented in `environments/green-science/README.md`.


EE source campaigns use exact loop-state checkpoints. Normal segmentation
commits the next bath, controller state and content-bound head before stopping;
completed loops retain their exact final bath. Legacy folders without that state
and partial/orphan writes fail clearly. New canonical Linux runs use
experiments/cthyb_exact_state.toml; the historical failed Linux artifact remains
separate. This workflow does not relax exact sampling input checks or establish
physical convergence from a checkpoint alone.

### EPW restart development checkpoint

An explicit model-plus-XOR restart archive and 20 K roundtrip driver/auditor
are under development. The archive validates its baseline hash and reconstructs
original field bits without silent fallback or refitting. This lossless path is
separate from the required pure fitted native restart, which has failed. The
[EPW evidence](EPW_INTEROPERABILITY.md#development-checkpoint-native-continuation-and-archival)
records size, test scope and remaining native and cross-platform gates.

The candidate EPW frequency-cross interface uses ordinary integer Matsubara
indices with declared binary64 beta. Its archived coordinates and arithmetic
order are part of the numerical contract: nearby floating frequency values
can fail conditioning even when the exact native-index query passes. The
explicit restart archive retains native correction information separately
from pure-model predictions; actual native readback passes for the declared
20 K workflow. Fit-only restart stability remains unsupported.

### EG supported route and scope

The Cauchy-frequency GECP fit shares the cached `LeastSquaresPlan`, immutable
helpers and safe NPZ primitives. `EPWRestartArchive` explicitly retains the
model and lossless native corrections; it does not silently substitute a fit.
The earlier generic `FermionicBasis`/`FermionicFit` route and a pure-fit native
restart are unsupported historical approaches, not accepted interfaces for
this EPW task. Original accuracy, noise and native comparison budgets remain
unchanged. See [EPW evidence](EPW_INTEROPERABILITY.md#supported-route-acceptance)
for measured storage and provenance; this is finite-case interoperability, not
continuous-frequency stability or a universal performance claim.

For a saved supported bundle, the approximate model and exact native data have
separate public operations:

```python
from kernelgecp.integrations.epw_restart_archive import EPWRestartArchive

archive = EPWRestartArchive.load("bundle.json")
green = archive.fit.evaluate_matsubara([0, 1, 2])
native = archive.reconstruct_native()
manifest = archive.write_native_restart("MgB2.imag_aniso_020.00")
```

The bundle carries its declared temperature, units, state order and chemical
potential offset. The fitted Green values are model predictions; native export
uses the explicitly stored corrections. This example requires the unreleased
1.2 implementation or its source checkout, not the published 1.1.0 package.


The source-only prerelease `LeastSquaresPlan.observation_projection` method
has moved into the optional spectroscopy adapter. Use public
`prepare_spectral_projection` for the complete observation map. Its cached SVD,
rank cutoff, ownership and covariance semantics are unchanged. This development
identifier was never part of the published 1.1.0 API.

### Explicit MaxEnt initialization

The optional `MaxEntConfig.minimizer_mu0` sets the stock Levenberg minimizer's
initial positive finite damping. It is not a lower bound. The default preserves
TRIQS/MaxEnt 4.0.1's 1e-18; scientific candidate campaigns explicitly record
1e-8 for every method. Full production acceptance is pending. Singular or
exhausted native solves still raise an explicit numerical failure with their
cause and completed trial history. The package does not retry another optimizer
or silently discard alpha trials. Archived datasets keep their original
calculation settings; read historical evidence with its recorded source.
