# MaxEnt spectroscopy evidence

Current EH acceptance remains open in [PR #142](https://github.com/tripp-smith/gecp-kernel-structure/pull/142). The measured final failure and normalized-model repair below supersede earlier preparation status; the synthetic results in this section are historical evidence for their frozen selector.

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

The first measured-data calibration at `d507050` is retained separately as
[failed evidence](data/spectroscopy/native-baseline-failure.json). Raw and GECP
pass both defaults; pydlr/flat has held-out RMS 2.058574 and fails the limit 2.
The next common calibration candidate tightens basis tolerance and SVD cutoff;
no noise allowance or mandatory provider is removed. Actual measured-data
acceptance and untouched final scoring remain open.

## Measured final failure and normalization repair

All eight measured calibration cases at `edca95c` and their independent audit
passed. The untouched dense final cohort then failed its first raw/flat case:
held-out whitened RMS 6.495244 exceeds the unchanged limit 2. The selected
spectrum had mass 1.010949. Independent adaptive quadrature gives training RMS
1.06817 before unit-mass rescaling and 5.04594 after it; held-out RMS is 2.80434
before rescaling and 6.49524 after it. Rescaling worsens the result, but does
not alone explain the failed held-out budget. These are post-hoc diagnostics,
not an alpha selection based on held-out data.

The [failure receipt](data/spectroscopy/native-final-failure.json),
[diagnosis](data/spectroscopy/native-final-normalization-diagnosis.json),
[failed final archive](experiments/data/cthyb-spectroscopy-edca95c-failed.tar.gz)
and [original calibration archive](experiments/data/cthyb-spectroscopy-edca95c-calibration.tar.gz)
preserve original output, worker logs, source/runtime identities and hashes.
The previous native-chi-square alpha selector is historical; its passing
calibration does not validate the repaired selector.

The repair keeps actual TauMaxEnt's full trial family and original chi-square
values. When normalization is requested, it selects among converged positive
finite-mass trials using the normalized physical forward model's training
chi-square per retained covariance rank. It records every trial mass and
selection discrepancy, applies the complete projection, and uses no held-out
values for selection. This is explicit postprocessing, not a normalization
constraint inside the provider's optimization. Both scientific auditors
independently integrate the normalized trial spectra before checking selection.

Recalibration retains the original kernelgecp-final observations. The consumed
`9d28053` dense cohort is rejected for future final acceptance; the prospective
final cohort is exclusively the disjoint `7200ec0` Linux dense-final streams,
seed base 30100000, after complete EE scientific acceptance and independent
binding. Synthetic final seeds become 2001/2002/2003 with 48 untouched quadrature
queries, replacing the consumed 1001/1002/1003 and 36-query validation. Numerical,
noise, covariance and provider gates remain unchanged. New calibration, audits,
final runs, full verification, exact-source CI and delivery remain mandatory.

## Normalized-selector evidence and second measured final failure

At `1452a8a`, all 72 fresh Mac synthetic cases and their independent audit
pass. Their [archive](experiments/data/spectroscopy-normalized-1452.tar.gz)
has SHA-256 `5c8b6752f3f6ed57b4c13a6e586ce5941026d30c58c45e771d928b5881105c77`.
The eight repaired-selector measured calibration cases and their independent
audit also pass; [archive](experiments/data/cthyb-spectroscopy-normalized-1452.tar.gz)
SHA-256 `1d9b502818f6818841fabb30780eb956dd4772d9a9e710137287afe3465672e6`.
The [checkpoint](data/spectroscopy/normalized-1452-checkpoint.json) retains
their original identities. Synthetic Linux run 40978733084 failed inside the
native provider before completing all alpha trials. An isolated diagnostic
replay passed, so its cause remains unresolved. The complete failed artifact
is preserved; no partial trial family is accepted and no provider controls
were changed. New failure instrumentation retains the input, exception cause
and completed trial history and restores the caller's NumPy error policy.

The new `7200ec0` dense-final observations have now been scored and are
consumed. Raw/flat and raw/Gaussian held-out RMS are 1.631636 and 1.596542.
GECP/flat fails at 2.023200438825109 against the unchanged limit 2 despite
training RMS 0.999109. Positivity, normalization and selected convergence pass.
The [failure receipt](data/spectroscopy/native-final-7200-failure.json) and
[complete partial archive](experiments/data/cthyb-spectroscopy-1452a8a-7200-failed.tar.gz),
SHA-256 `a742063f8edbf19ef49d1b5e02f4261aece3f338d0c84de9b596b6ebbfbda5f9`,
retain all worker outputs and the terminal parent log. Independent diagnostic
recomputation confirms the inspected arithmetic, but is not a passing final
audit. Finite-sample inverse-covariance variability is a plausible contributor,
not an established cause or uncertainty-coverage guarantee.

The [bounded next contract](RESEARCH.md#eh-bounded-fixed-bath-covariance-estimation-contract)
uses 256 new fixed-bath chains with a predeclared 128/128 split, the accepted
EE bath and unchanged physics and scientific gates. It requires a separately
bound acquisition and independent raw-chain audit, not a relabeled full EE
campaign. No new observations may be generated before source/protocol freeze;
no automatic retries, enlargement or covariance inflation follow a failed
final case. All synthetic, measured, full-root, CI and delivery gates remain.


The revised CI separates calibration from the single canonical acquisition.
Automatic runs cannot regenerate final chains. The explicit manual acquisition
runs 80/80/96 chains with strict checkpoints, then independently audits the raw
cohort and all 24 final workers. Root records its canonical run before data
acquisition; no workflow rerun or replacement dispatch after it begins is
permitted by this protocol. Each hosted job is bounded at 330 minutes, with
240-minute acquisition segments. Calibration-only CI success leaves measured
acceptance and release blocked. The Linux runtime-freeze preflight, complete
new final and every original scientific budget remain mandatory.


The final pre-acquisition resource plan uses actual Linux timing evidence,
not the earlier Mac estimate: approximately 10 hours of serial measurements
plus setup, continuation and audit. Four 64-chain segments replace 80/80/96
for scheduling headroom; the 256 samples, split, seeds and scientific budgets
are identical. Small actual-provider CI remains automatic; complete scientific
calibration and acquisition require explicit dispatch after source verification.
Calibration-only or ordinary CI success does not close the measured gate.
See [resource evidence](data/spectroscopy/fixed-bath-linux-resource-plan.json).
