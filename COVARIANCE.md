# EI real-data covariance results

The frozen UCI campaign passes all 18 numerical trials, the actual CPU
GPyTorch comparison and independent artifact audit at generation
`203e58ad168053b428c59af4ac5623ddc77c5ce8`. Calibration source is `00f4356`;
every scientific file/lock hash is bound to the actual Git blobs at both
revisions. [Calibration](data/covariance/calibration.json),
[full summary and raw-artifact hashes](data/covariance/summary.json),
[independent audit](data/covariance/audit.json), and the
[reproduction guide](experiments/COVARIANCE.md) retain scope and provenance.
Final native/root, exact-head CI and delivery remain separate gates.

The data are the attributed [UCI Beijing Multi-Site Air Quality dataset](https://doi.org/10.24432/C5RK5G),
Aotizhongxin station, under CC BY 4.0. Five meteorological features predict
PM2.5. Complete cases are chronological; training-prefix normalization,
calibration and test observations remain separate. Calibration selects RBF
lengthscale 3 and observation-noise variance 0.8 in normalized units, retaining
all rejected candidates. Final test data do not change these choices.

| Training size | Selected rank | Mean RMS error | Mean max error | Observation-std RMS error | Retained state | Maximum fresh RSS |
| --- | --- | --- | --- | --- | --- | --- |
| 512 | 128 | 1.25e-4 | 7.30e-4 | 1.45e-5 | 2.54 MiB | 136.59 MiB |
| 1024 | 256 | 1.19e-5 | 7.37e-5 | 3.82e-7 | 4.46 MiB | 151.05 MiB |
| 2048 | 256 | 5.58e-5 | 3.93e-4 | 5.00e-7 | 6.56 MiB | 251.59 MiB |
| 4096 | 512 | not computed | not computed | not computed | 20.24 MiB | 189.03 MiB |
| 8192 | 512 | not computed | not computed | not computed | 36.62 MiB | 261.58 MiB |
| 16384 | 512 | not computed | not computed | not computed | 69.37 MiB | 321.75 MiB |

Errors compare the same kernel/noise model with independent dense SciPy at
n<=2048, in units of the training target standard deviation. Observation
standard deviation includes the declared noise variance. Its maximum error
is 7.37e-5; both RMS and maximum targets (1e-3/1e-2) pass. Dense larger systems
were deliberately not formed. The small chosen rank is the first passing
tested calibration budget, not a proved minimal rank.

Across all sizes, the approximate-system solve residual is at most 1.62e-11
against the 1e-10 gate. Original-kernel preconditioned CG takes two or three
iterations, with maximum independently recomputed residual 1.03e-10 against
the 1e-8 gate. The auditor recomputes the complete original-kernel action
using coordinatewise distance sums, separately from the runner's cdist path.
All numerical arrays agree exactly across three fresh processes per size.
At 16384 points, numerical state is 69.37 MiB against 128 MiB and worst fresh
RSS is 321.75 MiB against 4 GiB. Source parsing, normalized arrays and retained
results are included in the state accounting; RSS covers the whole process.

The actual locked Mac CPU provider uses NumPy 2.2.6, SciPy 1.18.1, torch
2.14.1, GPyTorch 1.15.2 and linear_operator 0.6.1. On the real n512 frozen
factor, maximum mean/std parity error is 3.36e-13, solve difference 4.17e-15
and log-probability difference 7.96e-13. PR CI independently reruns the small
real case and actual Linux provider; final exact-source acceptance is pending.

| Training size | Median application time | Median bounded dense-reference time | Median fresh-process wall time |
| --- | --- | --- | --- |
| 512 | 0.0122 s | 0.00961 s | 1.41 s |
| 1024 | 0.0406 s | 0.0293 s | 1.20 s |
| 2048 | 0.0493 s | 0.181 s | 1.46 s |
| 4096 | 0.269 s | not computed | 1.82 s |
| 8192 | 0.488 s | not computed | 3.23 s |
| 16384 | 1.07 s | not computed | 7.68 s |

These are three fresh single-thread process repetitions on the recorded Mac
runtime. Application time includes normalization, factorization and prediction;
dense time measures its already-normalized covariance/solve/prediction call.
The separate compute/validation timer includes dense/CG checks and NPZ writing.
Parent-measured process wall time additionally covers interpreter startup,
data loading and final receipts. Small 512/1024 cases are slower than the
dense call and remain in the report. No general speed advantage or dense
16384 comparison is claimed.

This is faithful representation of one declared GP model, not a claim of
competitive air-quality forecasting or calibrated uncertainty. Test predictive
RMSE ranges from 58.4 to 77.2 micrograms/m^3 across the cohorts; that physical
prediction error is distinct from the much smaller representation discrepancy.
Kernel positivity remains an assumed column-callback property, not a new
continuous certificate. Large raw factors/logs stay outside Git; the full
scientific workflow uploads its own rerun artifacts. Small numerical summaries
and all local canonical raw-artifact hashes are committed here.
