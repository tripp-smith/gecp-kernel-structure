# Reproducing the EI covariance campaign

This is a finite real-data representation and linear-algebra experiment.
It does not establish a competitive PM2.5 forecaster, calibrated uncertainty,
optimal sensor placement, a continuous PSD certificate or universal speedup.
The full acceptance contract remains in PRODUCTION_INTEGRATIONS.md.

The source is [UCI Beijing Multi-Site Air Quality](https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data),
attributed to Chen, S. (2017), DOI
[10.24432/C5RK5G](https://doi.org/10.24432/C5RK5G), under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The 8,192,212-byte outer archive SHA256 is
`b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4`;
the Aotizhongxin station CSV SHA256 is
`a94fbcfc71708b6ffdc033360163d91efb7f675e496f6b5860929aa96b96351b`.
Both hashes are checked before parsing. Downloaded archives and large frozen
factor artifacts belong outside Git; the committed calibration/summary names
below are populated only after their actual generation succeeds.

The features are TEMP/PRES/DEWP/RAIN/WSPM, in Celsius/hPa/Celsius/mm/m/s,
with PM2.5 in micrograms per cubic meter as target. Keep complete cases only
in those six columns and enforce strictly chronological station rows. The
first 16,384 complete rows form the training pool, the next 512 calibration,
and the next 512 test. Cohorts 512/1024/2048/4096/8192/16384 use training
prefixes and fit all feature/target normalizers on their own prefix only.
The source has 34,119 complete cases and excludes 945 missing rows.

Calibration evaluates lengthscales 0.7/1.5/3 and observation-noise variances
0.05/0.2/0.8, in standardized units. Choose calibration predictive likelihood
subject to rank512 meeting tenfold stricter dense-accuracy margins at n2048.
Then evaluate ranks 32/64/128/256/512 at every cohort. Small cohorts choose the
smallest passing calibration rank; larger cohorts explicitly use 512 because
the independent dense oracle is bounded to n2048. Retain all rejected and
slower cases. No accuracy against a dense n16384 solve is claimed.

The reported standard deviation is the posterior predictive **observation**
standard deviation, including the declared observation-noise variance, divided
by the training target standard deviation. Both methods use this same
convention. Raw physical predictive RMSE is reported separately. The kernel
amplitude is one in normalized units; no hyperparameter or rank is retuned
using test errors. After a test-driven decision, select a new untouched test
block before presenting another final-validation claim.

Use the same locked CPU environment for every scientific command:

```bash
uv sync --project environments/gp --frozen --python 3.12
export PYTHONPATH="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
runtime_python="$PWD/environments/gp/.venv/bin/python"
```

Default development environments can have a different NumPy version and are
rejected for scientific generation. NumPy 2.2.6 and SciPy 1.18.1 must match
calibration, all campaign processes and the actual GPyTorch provider. Record
actual Python/platform versions; selection is deterministic and uses no
stochastic seed. Commit the implementation before calibration and leave
scientific sources/locks frozen during calibration and final generation.

```bash
"$runtime_python" -m experiments.covariance_benchmark calibrate \
  --data /tmp/uci-501.zip --output /tmp/ei-calibration.json
```

An inconclusive calibration exits nonzero and retains rejected candidates.
After a successful calibration, retain it as data/covariance/calibration.json
and commit that metadata without changing any scientific source/lock hash.
Git-blob binding permits this distinct metadata revision while preserving
the original calibration source identity.

```bash
"$runtime_python" -m experiments.covariance_benchmark campaign \
  --data /tmp/uci-501.zip --frozen data/covariance/calibration.json \
  --output /tmp/ei-campaign
"$runtime_python" -m experiments.covariance_benchmark provider \
  --data /tmp/uci-501.zip --frozen data/covariance/calibration.json \
  --output /tmp/ei-provider.json
"$runtime_python" -m scripts.audit_covariance \
  --manifest /tmp/ei-campaign/campaign.json --provider /tmp/ei-provider.json \
  --data /tmp/uci-501.zip --frozen data/covariance/calibration.json \
  --source-revision "$(git rev-parse HEAD)" --output /tmp/ei-audit.json
```

The campaign launches each selected cohort in three fresh processes. Retain
all NPZ factors, JSON receipts and process logs. The auditor verifies hashes,
actual source blobs at each declared revision, all 18 exits/cases, numerical
repeat equality, original dataset/normalization, low-rank system residuals,
independent full true-kernel CG, and independently assembled small dense
predictions. It also requires the actual identical-factor GPyTorch run.
Its coordinate-difference kernel computation differs from the runner's cdist.
Both implementations block large products in two dimensions.

Small-case RMS/max mean and observation-standard-deviation discrepancies must
be <=1e-3/1e-2 in normalized units, approximate solve residual <=1e-10, original
kernel CG residual <=1e-8, retained numerical state <=128 MiB and fresh-process
RSS <=4 GiB. Memory components include the factor/solver, source observations,
normalization split and result arrays; RSS includes the complete process.
Slower and failed cases remain evidence. Timing is reported per construction,
prediction, dense reference, true-kernel CG and complete workflow, rather than
treated as an automatic speedup claim.

The dedicated covariance workflow runs a small attributed case and actual CPU
provider on PRs. Explicit workflow dispatch can rerun the complete campaign
within the 330-minute job budget. A small-case CI pass is not the full domain
gate. Phase closure also requires ./scripts/verify.sh, exact-source CI,
reviewed documentation, authorized merge and confirmed pushed-main state.
