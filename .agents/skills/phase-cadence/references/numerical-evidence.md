# Numerical certificates and canonical experiments

Read for changes to numerical certification or canonical experiment generation.
Apply the relevant sections, not every case to every numerical edit.

## Certificate scope and arithmetic attacks

State the exact objective being bounded: an ideal kernel/run, a frozen rounded
model, or a selected-core solve. State callback and regularity assumptions and
which quantities are interval-enclosed. Higher working precision is not itself
a certificate. Report actual component precision, including exported binary64
coordinates/bounds; one headline bit count can hide mixed-precision boundaries.

Before declaring a bound/stopping repair sound, attack the affected boundaries:

- Adjacent representable endpoints and nonrepresentable midpoints; dyadic
  cutoffs/tolerances and their immediate neighbors; zero and subnormal scales.
- Cancellation, overflow, mixed-scale underflow and transposed inputs. Repairing
  one overflow example must not erase another representable contribution.
- Repeated coordinates, ill-conditioned cores, ties/oversized tie tolerances,
  and all relevant budget or coordinate-resolution exhaustion paths.
- Caller precision, rounding and trap settings. Set required contexts explicitly
  and restore their actual prior state, not a rounded proxy for it.

Use exact rational values for representable inputs where possible, plus an
independent formulation/solve for nontrivial checks. Require conservative lower
and upper bounds, nonnegative ordered gaps, outward error accounting and stopping
thresholds that do not round in favor of success. Underflow is not proof of zero
error. If assumptions, budget, conditioning or resolution prevent certification,
return an explicit uncertified/failure status, not a zero gap by default. Ask an
independent reviewer to attack the proposed fix when available; label any missing
review honestly. Retain the smallest failing input as a regression.

## Pilot, freeze, generate, validate

1. Use labeled pilots for tuning; preserve failed/budget-exhausted cases.
2. Freeze the source revision, configuration, seeds, dependencies and precision
   before canonical generation. Launch a fresh process so stale imported code
   cannot disagree with on-disk provenance. Do not alter the source mid-run.
3. Record generation revision, source/config/lock hashes and environment. Detect
   unexpected changes before replacing canonical outputs. A later maintenance
   lockfile does not retroactively become the historical generation environment;
   verify against the recorded Git blobs rather than rewriting old manifests.
4. For a deterministic-output claim, repeat the frozen configuration and compare
   canonical bytes, with timing or other nondeterministic diagnostics separate.
   Retain the run's process/session, log, output location and configuration so a
   long job can be resumed/observed after a tool timeout rather than duplicated.
   File existence alone does not establish a successful process exit.
5. Keep final validation points untouched by fitting, refinement, rank selection
   and stopping. If validation changes a decision, those samples are now tuning
   data: disclose this and use a new untouched final holdout. Independently check
   representative results using a different formulation/factorization, not an
   alias of the same implementation. Distinguish training pass, controller
   convergence, heldout pass and certification in the artifact and report.

Mechanism-specific acceptance is essential: adaptive sampling must actually
refine from observed residual/surplus information and preserve its refinement
history. A fixed dense grid is not evidence of that mechanism. For comparisons,
give methods the same permitted observations and budget; label truth-assisted
controls separately. A sampled/heldout pass is not a continuous-domain certificate,
and the smallest tested passing rank is not a proved minimum.
