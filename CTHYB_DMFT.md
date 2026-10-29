# CT-HYB Bethe DMFT evidence

This is a bounded local scientific result, not Phase EE delivery. Actual TRIQS
4.0.2/CT-HYB 4.0.1 on Darwin arm64 produced the finite-grid stochastic campaign
at generation revision `9d28053060329cd20647dc03e038a383c351adf1`. An independent
audit at revision `e83d195cb11ec9e2657c24327e95489edff1c5d6` passed the seven
reported gates. The original EE requirement includes a real Linux execution and
independent audit in CI; that and exact-head delivery remain open. The local
root verification at `e83d195` reported 2,345 passed and 39 skipped. This local
result is not a continuous-domain certificate, exact Hubbard-model result,
performance claim, or evidence of trajectory-error bounds.

## Protocol and independent observations

The paramagnetic single-band Bethe model uses beta 10, U 2, hopping 1 and
mu=U/2. Dense and GECP paths use measured CT-HYB feedback, native moments and
explicit tail handling; the GECP path reconstructs the measured imaginary-time
Green function. The calibration used eight fresh distinct odd-stride streams
per budget. The independent audit found maximum standard errors 0.0008387404,
0.0003650166 and 0.0002455173 at 100,000, 400,000 and 1,600,000 cycles:
the first two fail the unchanged 0.0003 bound, and the third passes. The
frozen choice is 1,600,000 cycles and cutoff 8. Final dense and GECP estimates
use 64 fresh fixed-bath replicas per path, from disjoint seed families; the
audit counts 208 independent fixed-bath chains across calibration, final and
refinement. Unique seed mapping and non-identical fingerprints do not prove
probabilistic independence. The paths are not paired after their baths diverge.

The audit independently recomputed covariance of chain means and pivoted-QR
representation replays. Its maximum final observable standard error is
8.4563e-5 (limit 3e-4). The largest dense/GECP observable difference is
9.5774e-5; all 14 observables meet their individual 1e-3 plus three combined
standard-error limits. These standard errors describe fixed-bath chains, not
uncertainty in the stochastic DMFT trajectory.

Both stochastic controllers stop after 12 iterations; their finite thresholds
are not deterministic convergence certificates. Independently checked maximum
absolute fixed-bath equation defects are 1.8049e-4 for dense native,
2.0141e-4 for GECP native and 2.2221e-4 for GECP represented feedback; each
component passes its own acceptance limit. Native-grid refinement passes with
maximum observable difference 1.2946e-4. Independent low-frequency quadrature
has maximum discrepancy 1.7479e-5. All three frozen-measurement provider replays
and both cutoff checks pass; individual QR discrepancies are in the audit and
summary. These are separate representation, refinement and physical checks,
not replacements for the stochastic controller or uncertainty checks. Native
density-matrix Green and self-energy moments are recorded in the raw archives;
they are not inferred from a noisy tail fit.

## Durable evidence and verification

- [Campaign archive](experiments/data/cthyb_9d28053_campaign.tar.gz) contains
  the complete frozen campaign tree, including raw per-measurement arrays,
  checkpoints, failed calibration budgets, freeze, manifest and result. Its
  SHA-256 is `c50e2db4464001e6bc30a948e14a0da9176d3d005d3ae3c7cb43d13c3267bfd3`.
  The transient empty lock file was excluded.
- [Independent audit](data/cthyb/independent-audit-e83d195.json),
  [binding](data/cthyb/audit-binding-e83d195.json), and
  [portable summary](data/cthyb/summary.json) record hashes, complete gates,
  values and revisions. The binding retains the actual historical `/tmp`
  command and runtime path as provenance; those paths are not portable inputs.
  The summary uses repository-relative artifact paths. Verify binding hashes
  against the extracted campaign tree; no paths or original receipts have been
  rewritten to imply a different execution.
- The original [45d42c7 inconclusive campaign](experiments/data/cthyb_45d42c7_inconclusive.tar.gz)
  remains in the repository. Subsequent independent investigation found native
  adjacent even/odd seeds collide (`seed | 1U`), invalidating nominal
  independence in both 45d42c7 and the consumed 8818c3a revalidation. The
  [failed independence audit](data/cthyb/invalid-independence-audit.json)
  records rejection. The latter raw directory currently exists at
  `/tmp/gecp-ee-revalidation-8818c3a`; its compressed archive was 35 MB, above
  the 30 MB packaging threshold. Its SHA-256 and size are in the summary; a
  durable release asset is still needed for that invalid historical record.
  These failed or invalid cohorts are never pooled into the corrected campaign.

To inspect the raw evidence, extract the campaign archive and compare each
`measurement_archives` SHA-256 in the binding against its relative file. To
rerun the independent audit, use the pinned provider runtime and the audit
revision's `experiments.audit_cthyb_dmft` against the extracted tree; the
binding preserves the original command and source/runtime identity. No new
Monte Carlo execution is needed for such an audit.

## Linux checkpoint failure and new canonical run

Run 40853707061 completed calibration and its first validation segment, then
failed resuming the compressed loop on another runner. The complete
[failed checkpoint](experiments/data/cthyb-linux-e83d195-incomplete.zip) and
[receipt](data/cthyb/linux-checkpoint-failure.json) are retained. Its 164
measurements are historical evidence and are not imported into the next run.
Read-only cross-platform replay found compressed bath differences around
8.8e-16; dense replay was exact. The failed log does not identify its precise
operand, so that replay is evidence of a failure mechanism rather than proof
of the exact failing operation.

New runs persist each exact next bath, controller history and consecutive-pass
counter in a predecessor-bound state chain. A separate atomic head commits the
step before the normal segment-stop notification. Completed loops retain their
exact final averaged bath and report. Resume checks configuration, initial bath,
all measurement hashes/seeds/cycles, every state link and the final mean's
binary64 rounding allowance; it returns the stored bits without replaying prior
representation updates. The sampling bath equality check is unchanged.

An arbitrary interruption between measurement, state and head writes remains a
fail-closed orphan; this is not general crash recovery. Legacy checkpoints
without state are rejected. Hashes provide content integrity, not signatures.
`experiments/cthyb_exact_state.toml` declares a fresh disjoint seed family with
base 30100000 and unchanged physics, cycle ladder, ensemble sizes and numerical
budgets. A new complete Linux campaign and its independent audit remain required.
