# Synthetic MaxEnt spectroscopy evidence

The frozen synthetic campaign at clean source revision
`dc97680a0fa01789cefd1b699a3fb4c0f00a2699` used actual TRIQS MaxEnt
4.0.1 `TauMaxEnt` on macOS arm64. Calibration has 24 passing cases and the
final manifest has 72 passing fresh-process children. The independent audit
reports 24 independently audited calibration cases and 24 independently
audited final cases among the 72 children; it passes both receipts. A later
merge or documentation update does not make the consumed final data blind
again. Actual CT-HYB spectroscopy, full root verification, Linux CI and
delivery remain open. This is finite synthetic evidence, not proof of unique
peak recovery, complete uncertainty coverage, or a universal method advantage.

The frozen configuration compares raw, KernelGECP, pydlr and sparse-ir on
independent and correlated off-grid atomic spectra and a correlated smooth
spectrum, with flat and broad Gaussian defaults. Calibration seeds are
101/102/103 and final seeds 1001/1002/1003; calibration and final query orders
are 24 and 36 respectively. Controls include epsilon 1e-8, 33 alpha points,
401 frequency points and three process repeats. Eight bootstrap continuations
are performed only on repeat zero; the other two repeats are not equivalent
bootstrap workloads. The provider and numerical dependencies, source hashes,
thread settings and actual runtime paths are recorded in the frozen receipt.

The full observation covariance is retained with each dataset and transformed
as `P covariance Pᵀ`; the forward kernel is transformed as `P K`. The archived
projection, source-time coordinates and dataset permit independent checking.
The former `6fe118e` calibration's independent audit failed on projection
symmetry for pydlr. Its candidate data and failed audit remain in the evidence
archive. The corrected same-mesh projection uses retained SVD factors
`U_r U_rᵀ`; it did not silently relax the tolerance.

Across the final children, positive/normalized spectrum, selected-provider
convergence, held-out whitened residual and coefficient/projection agreement
gates pass. Effective ranks are 65 raw, 21 KernelGECP, 20 pydlr and 20
sparse-ir; the largest recorded projection-SVD condition number is about
1.10e8. CDF integrated absolute error reaches 0.06356 and absolute window
weight error reaches 0.02826 across the finite set. The largest reported
bootstrap window standard deviation is 0.06707; this empirical variation is
not a coverage guarantee. The largest observed difference between default-model
window weights for matching case/method/repeat is 0.000296. Individual CDF,
window, bootstrap, conditioning, per-gate and default-specific values are in
the [72-child inventory](data/spectroscopy/summary.json).

Process wall-time ranges in seconds are raw 4.85–122.60, KernelGECP 3.41–80.64,
pydlr 3.18–73.34 and sparse-ir 3.42–77.90; they include the uneven bootstrap
workload. Peak RSS ranges in bytes are 124,583,936–156,286,976,
123,535,360–163,987,456, 112,918,528–156,499,968 and
129,269,760–163,397,632 respectively. These are process measurements, not
comparative speedup claims. The manifest records each of the three fresh
processes, input, output, log hashes and status.

The [campaign archive](experiments/data/spectroscopy-dc97680.tar.gz) retains
all calibration and final case datasets, projections, spectra, logs, fit
archives, manifests and the original audit receipts, plus the earlier failed
calibration and failed audit. Its SHA-256 is
`fdf3d4508d373abef4e12dce5927a0986e1ba723bd4650abf0d53a62de636dbb`.
The [freeze](data/spectroscopy/freeze.json),
[calibration audit](data/spectroscopy/calibration-audit.json),
[final audit](data/spectroscopy/campaign-audit.json),
[failed audit](data/spectroscopy/failed-calibration-audit.json) and portable
summary are also available separately. Historical runtime paths in receipts
are provenance, not portable paths. The audit and manifest hashes bind the
original files; the archive does not rewrite their content.
