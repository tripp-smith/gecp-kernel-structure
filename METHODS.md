# Computational methods and reproducibility

`kernelgecp` studies Gaussian elimination with complete pivoting (GECP) for the
fermionic kernel and supplies a reusable physical Green-function basis. Formal
exact-run results, sampled application comparisons and finite-query numerical
certificates have different scopes. This document describes the retained
application experiments; [SPEC.md](SPEC.md) and [FINAL_AUDIT.md](FINAL_AUDIT.md)
define the original mathematical scope.

## Physical convention and numerical representation

For inverse temperature `beta`, spectral cutoff `wmax`, physical time
`tau in [0,beta]`, dimensionless frequency `Omega=beta*energy` and
`t=tau/beta`, the positive kernel is

$
K(t,\Omega)=\frac{e^{-t\Omega}}{1+e^{-\Omega}},\qquad
G(\tau)=-\int K(\tau/\beta,\beta\epsilon)\,d\rho(\epsilon).
$

The Fourier transform uses `exp(+i*wn*tau)`, with ordinary integer indices
`n` and `wn=(2*n+1)*pi/beta`. A unit pole at energy `epsilon` therefore has
`G(i*wn)=1/(i*wn-epsilon)`. Sparse-ir's reduced odd index is explicitly
`2*n+1`; the pydlr adapter builds a complex-linear frequency design from a
real identity rather than applying its matrix conjugation/transpose convention
to general complex batches.

A selected core `C[i,j]=K(tau_i/beta,Omega_j)` defines cardinal rows through
`C.T @ D(tau).T = k(tau).T`. Implementations cache LU factors and solve this
system; they do not explicitly invert a core. A cached thin SVD fits supplied
observations and supports arbitrary complex batch/matrix axes. Singular values
are retained when `s > rcond*smax`. `tol` controls finite construction-grid
stopping and is not a physical error guarantee for arbitrary observations.

Physical coordinates and exported arrays are binary64/complex128. Reference
selection and cardinal evaluation use the reported working precision, normally
128 bits. Least squares and final products remain NumPy arithmetic. The optional
Rust backend changes only 128-bit reference time-kernel evaluation and cached
triangular substitution; Matsubara evaluation remains Python/mpmath. Python is
the default backend. Archive loading reconstructs factors and does not turn
saved historical diagnostics into new certification evidence.

For finite physical checks we report
`max(abs(predicted-truth))` over every query and matrix entry, using the complex
modulus entrywise. This is not a matrix spectral or Frobenius norm. Relative
errors, where reported, divide that maximum by `max(abs(truth))`; acceptance in
the experiments below uses the declared absolute target. Conditioning refers
to the reported exported matrix and coordinate system, not an invariant quality
ordering between different basis parameterizations.

## Matched DLR/IR comparison: DV

[The configuration](experiments/dlr_application_benchmark.toml) compares GECP,
pydlr and sparse-ir at cutoffs 10, 100 and 1000, plus a cutoff-20 solver case.
Every provider receives identical complex matrix observations and the same
cached SVD fitter. A predeclared native-eps sequence is selected only on common
time/frequency calibration data. The final target is `1e-8`; cohort B uses a
`0.1` calibration guardband and untouched final points. Equal native eps values
are not assumed to mean equal accuracy. The first passing tested rank is not
claimed minimal, and the common-fitter result is not a contest between providers'
native optimized fitting implementations.

Separate fresh workers perform calibration and selected measurements. Reports
separate construction, design/SVD setup, repeated fitting, cached/uncached
evaluation and solver work. Process high-water RSS includes imports, harness
and native libraries; it is not incremental basis allocation. Fixed thread
environment settings are recorded, not inferred operating-system thread counts.

The solver is the half-filled **noninteracting** Bethe equation
`G(z)=1/(z-hopping**2*G(z))`, solved with cached projection and damping. Stopping
requires the actual unprojected equation defect on the observation grid.
Independent time/frequency truth and observed-grid truth are separate checks;
observed-grid truth is explicitly not held out. This is not an interacting
impurity solver or a production materials calculation.

Cohort B passes all 12 method/case records. Its standalone 95-decimal-digit audit
uses finite-pole formulas, the causal semicircle resolvent and selected direct
spectral integrals. It imports no benchmark truth helper. Cohort A's sparse-ir
cutoff-10 frequency miss (`1.20865e-8`) remains retained. The earlier missing-low-
frequency observation failure likewise motivated a new sampling protocol rather
than being hidden. All results remain finite-fixture, single-machine evidence.

## Finite-query arithmetic bridge: DW

`certify_fit_precision` encloses the exact physical selected core and the unique
full-column-rank least-squares fit to the supplied exact binary64 observations.
Owned interval contexts enclose physical coordinate division, kernels, core
solves and normal equations. The actual fit may use truncation or unrelated
coefficients; its discrepancy from this particular reference can legitimately
be large. The certificate includes coefficient/design errors and exact binary
rational checks of final dot-product rounding. It does not certify unknown
physical truth or all points in a continuous domain. A zero-containing pivot,
insufficient precision or exhausted budget is inconclusive, not a zero error.

`certify_grid_pivot_precision` compares the actual frozen rounded residual to
the exact selected cross on the supplied finite grid. If the discrepancy is at
most `delta`, the rounded candidate is at least `L>delta`, and the rounded global
bound is `U`, its true multiplicative quality is at least
`(L-delta)/(U+delta)`. The Lean theorems
`Discretization.approxPivot_of_uniform_discrepancy` and
`Discretization.dot_error_le_of_component_bounds` prove the conditional transfer
and complex error-budget algebra. They do not formally verify the interval
engine or turn an inexact trajectory into an exact complete-pivot run.

The six-case cohort preserves five completed enclosures and one inconclusive
256-bit uniform-grid fit. A 512-bit uniform-grid reference exposes the mismatch
between a severely ill-conditioned full-rank problem and the deliberately
truncated actual fit. Endpoint resolution and selected-node fitting improve
conditioning. Nevertheless the resolved large-cutoff physical errors reach
`1.29e-8` in time and `3.37e-8` in frequency despite tiny arithmetic discrepancies.
A small arithmetic certificate is therefore not a physical accuracy pass.

The independent auditor checks stored binary values rather than replaying
another architecture's BLAS/libm byte patterns. It uses 1280-bit centered-kernel
bordered determinants, exact selected-node coefficient reproduction and selected
independent predictions. Other LS checks are the retained generation-time QR
checks, not a second independent recomputation by that auditor. Its separate
64-epsilon fixture-identity guard never relaxes certificate error bounds.

## Optional native acceleration: DX

The frozen comparison uses three alternating fresh processes per backend at
`beta=20`, `wmax=50`, rank 55, actual 128-bit application arithmetic and default
cache capacity eight. One selected-node complex matrix fit precedes 80 distinct
64-row time meshes. The application timer includes construction, observations,
fit and native loading/initialization. Full-process wall also includes startup,
provenance checks, independent validation and serialization. Actual underlying
time-design calls are counted; the 100 repeated final-grid calls make none.

The predeclared retention rule requires at least 20% median application and
full-process reduction, paired discrepancies at most `1e-12`, physical errors
at most `1e-8`, and the stated independent QR checks. Median application times
are 42.795850 s Python and 21.332023 s Rust (50.15% reduction); full-process times
are 44.703104 s and 23.350768 s (47.76%). The eight-query prefix improves only
8.92%, and warm cache times differ by 0.32%. The gain concerns repeated changing
meshes, not every workload. Construction remains expensive.

Maximum paired discrepancy is `4.25522e-15`: production outputs are **not
bitwise equal**. Physical stream errors are below `2.57e-11`. Independent
512-bit centered-kernel QR checks use exact physical coordinates and disjoint
queries, applied to the actual frozen coefficients. Within each backend, all
three trials reproduce complete numerical arrays and oracle values. Replaying
saved-record validation is a consistency check, not another independent oracle.
Timing describes the recorded macOS arm64 build and machine; it does not promise
the same gain for another wheel, platform or input.

## Frozen data and reproduction

| Experiment | Canonical data | Generation commit |
| --- | --- | --- |
| DV matched cohort B | [dlr_applications_cohort_b.json.gz](experiments/data/dlr_applications_cohort_b.json.gz) | `8a4dd6bec0d606764667a1721e06469286e82cef` |
| DW portable cohort C | [precision_bridge.json.gz](experiments/data/precision_bridge.json.gz) | `a4995b471506fe043f3102114d42620c3c645d76` |
| DX native cohort | [rust_reference_benchmark.json.gz](experiments/data/rust_reference_benchmark.json.gz) | `e5f6e796893d31837545e6c7eac0918b641e3fda` |

[RESEARCH.md](RESEARCH.md) records the container/payload/source/configuration
hashes, measured environments and preserved failures. The DV and DW independent
receipts are adjacent `*_audit.json` files. Later documentation, wheel-path
repairs or lockfiles do not rewrite these historical generation environments.
Use a separate detached checkout of the appropriate commit and the recorded
Python/dependency/native build identities for historical reruns. Generators
require a committed, tracked-clean source and fresh workers.

```sh
# DV: use its frozen checkout and exact optional-provider requirements.
uv run --with-requirements experiments/dlr_benchmark_requirements.txt \
  python -m experiments.dlr_application_benchmark \
  --config experiments/dlr_application_benchmark.toml \
  --output experiments/scratch/dv-repeat.json

# DW: regenerate at its source revision; retain failed/inconclusive cases.
uv run python -m experiments.precision_bridge \
  --config experiments/precision_bridge.toml \
  --output experiments/scratch/dw-repeat.json

# DX: in its frozen checkout, with Rust 1.92.0 and the documented C tools.
uv sync --python 3.12.12 --all-extras --frozen
export RUSTUP_TOOLCHAIN=1.92.0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RAYON_NUM_THREADS=1
export KERNELGECP_SOURCE_REVISION="$(git rev-parse HEAD)"
export PYO3_PYTHON="$(uv python find)"
export CARGO_TARGET_DIR="$PWD/native/target"
export GMP_MPFR_SYS_CACHE=_
cargo build --manifest-path native/Cargo.toml --release --locked \
  --features extension-module
case "$(uname -s)-$(uname -m)" in
  Linux-x86_64) native_library="$CARGO_TARGET_DIR/release/lib_native.so" ;;
  Darwin-arm64) native_library="$CARGO_TARGET_DIR/release/lib_native.dylib" ;;
  *) echo "This native experiment targets Linux x86_64 or macOS arm64" >&2; exit 1 ;;
esac
cp "$native_library" kernelgecp/_native.abi3.so
uv run python -m experiments.rust_reference_benchmark \
  --config experiments/rust_reference_benchmark.toml \
  --native-binary "$PWD/kernelgecp/_native.abi3.so" \
  --output experiments/scratch/dx-repeat.json
```

The DX commands compile a local experiment module, not a distribution wheel.
Official wheel/source/notices builds use the inspected corresponding-source
procedure in [native/SOURCE_BUILD.md](native/SOURCE_BUILD.md). Choose a new output
name for each generation. A generator's zero exit means evidence was recorded;
inspect its status, certification and retention fields.

For independent rechecking use `scripts/audit_dlr_applications.py` and
`scripts/audit_precision_bridge.py` with the immutable generation SHA and
container/payload hashes shown in [APPLICATION.md](APPLICATION.md). These audit
commands and `./scripts/verify.sh` have different roles: scientific evidence
validation does not replace Lean/Python/archive/install/CI release checks.
Gzip storage leaves the original JSON payload unchanged.

## Citation and related implementations

Cite the software release or commit used via [CITATION.cff](CITATION.cff), and
cite relevant upstream methods/libraries independently. The official
[libdlr citation guidance](https://libdlr.readthedocs.io/en/v1.0.0/background.html#citation-information)
and [sparse-ir citation guidance](https://sparse-ir.readthedocs.io/en/latest/#license-and-citation)
identify their own papers. This repository's experiments and formalization do
not establish literature priority or replace those citations.
