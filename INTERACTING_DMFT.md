# Deterministic interacting DMFT evidence

This is the ED numerical/model validation, targeting kernelgecp 1.2.0. It
does not establish exact Hubbard physics, a continuous error certificate, or
an interacting-workflow speed advantage. The public solver contract and
example are in [APPLICATION.md](APPLICATION.md#interacting-bethe-ipt-unreleased-targeting-120).

## Frozen calculation

Generation source is `a3a64baa4f48a7627a7f7bd9edf2410c327c92f7`, with a clean
tracked tree and the committed [configuration](experiments/interacting_dmft.toml).
The [full report](experiments/data/interacting_dmft_a3a64ba.json.gz) retains
13 calibration trials, nine doubled-cutoff controls, 27 independent compressed
closed loops, nine native closed loops and common-native-trajectory replays.
All final arrays repeat identically across three fresh processes for each
method/U pair. Timings and other nondeterministic diagnostics are separate.
The [independent audit](experiments/data/interacting_dmft_a3a64ba_audit.json)
checks the saved arrays, source/config bindings, physical diagnostics, actual
equations, final residuals, complete replay counts and repetition identity.

The half-filled paramagnetic Bethe IPT model uses beta=10, t=1, U=0/1/2,
mixing 0.5, at most 200 iterations, and two consecutive unmixed defects below
1e-8. Hartree and chemical-potential constants cancel analytically. Green
wmax=4 and dynamic-Sigma wmax=12 pass calibration; controls double both.
GECP uses 128-bit Python reference construction/evaluation with binary64
design export and complex SVD fitting. DLR/IR native truncation tolerances
have their own meanings; they are selected by measured common error budgets.

Actual TRIQS 4.0.2 is the dense Fourier reference, with declared PH moments.
Primary native meshes have 1024 signed ordinary frequencies and 6145 time
rows; refinement uses 2048/12289. Compressed Green fits omit eight reserved
calibration frequencies; all methods use the same compressed observations.
Final frequency queries are outside the entire primary grid, and 64 final
Gauss-Legendre times overlap neither primary nor refined uniform rows.
Raw Dyson values at omitted modes do not become held-out predictions.

## Results and limitations

| Required quantity | Observed maximum |
| --- | ---: |
| Independent closed-loop G/Sigma error, either domain | 6.87e-8 |
| Common-native-trajectory replay error, either domain | 8.15e-8 |
| Native grid-refinement change | 7.15e-11 |
| Per-spin density error from 0.5 | 8.14e-8 |
| Separate asymmetric-pole/atomic representation fixture error | 8.85e-7 |
| Independently recomputed scalar equation/product residual | 8.33e-16 |

Density, sampled causality, PH symmetry, endpoint mass, SVD conditioning,
projection residuals and tail fits are recorded separately. Interacting tail
fits have two-window sensitivity diagnostics without a remainder certificate.
An independent 70-digit finite-pole quadrature verifies signs/normalization;
its six-term Laurent remainder 1.416e-8 is below its analytic finite-measure
bound 2.345e-7. That bound applies only to that finite fixture.

At U=2, the following are medians of three fresh-process repetitions on the
recorded macOS arm64 runtime, using one thread per process. Application time
includes construction, preparation and the loop after provider imports;
loop RSS is the process high-water mark before final truth/replay evaluation.

| Method | Green/Sigma rank | Application seconds | Loop seconds | Loop RSS MiB |
| --- | --- | ---: | ---: | ---: |
| Native TRIQS Fourier | full grids | 0.0218 | 0.0176 | 112.6 |
| GECP Python reference | 18/25 | 20.59 | 0.0473 | 195.8 |
| pydlr | 20/27 | 0.0823 | 0.0544 | 121.0 |
| sparse-ir | 20/27 | 1.149 | 0.0663 | 134.6 |

Native Fourier is faster on this small beta=10 workload. Reference GECP setup
dominates its complete workflow; these results support interoperability and
accuracy, not a performance advantage. The table does not measure a reused
basis campaign, the optional Rust backend or lower-temperature regimes.
Full-worker wall/RSS also includes validation and compressed-only native
replay, so those broader metrics must not be substituted into this comparison.
Coefficient bytes and dense output bytes have different storage scopes.

## Reproduction and failure evidence

Use the actual locked TRIQS environment in [environments/README.md](environments/README.md)
and install the hash-locked `environments/locks/dlr-physics-pip.txt` without
replacing its conda NumPy/SciPy. From a clean committed checkout:

```sh
python -m experiments.interacting_dmft --output-dir artifacts/interacting-dmft
python -m experiments.audit_interacting_dmft artifacts/interacting-dmft --output artifacts/dmft-audit.json
```

Provider CI reruns these commands and retains complete raw arrays, requests,
logs, source hashes, fit archives, environment freeze and audit receipts as
workflow artifacts. The small attributed actual-native fixture under
`tests/data/dmft_native_ipt_v1.*` also runs in ordinary package CI.

Unfrozen pilots remain tuning evidence. They exposed complex-scalar JSON
serialization and review gaps in import timing, replay/RSS scope, empty replay
and validation reuse. A deliberate negative test also caught a CLI that could
exit zero when only the repetition gate failed; it now returns failure.
Actual artifact auditing caught two oracle integration mismatches: dense and
compressed solvers use different stop-reason strings, and NumPy versions round
independently generated Gauss coordinates by a few ulps. The auditor now
recognizes the actual native status and permits only coordinate-rounding
slack; saved prediction/reference coordinate arrays remain exactly equal.
Numerical error budgets were not relaxed. Exact-head/root/main delivery
receipts in the phase PR determine final acceptance beyond this source run.
