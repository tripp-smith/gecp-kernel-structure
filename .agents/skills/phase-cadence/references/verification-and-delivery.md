# Verification gates and revision-specific delivery

Read when changing a gate, closing a release/final handoff, or recording delivery
evidence. User authorization and repository policy still govern external actions.

## Test that gates reject bad inputs

When a critical gate is introduced or changed, exercise intentional rejection
fixtures in an isolated temporary directory or captured-output test harness,
not in normal Lean library modules. Existing gates need a fresh negative pass
at release closure if their rejection behavior lacks current evidence.

- Check process exit status, not only printed warnings. Cover a forbidden axiom,
  missing audited declaration, empty/malformed/truncated output and hidden Lean
  or pipeline failure where relevant. Require expected audit coverage as well
  as an allowed-axiom list; missing output must not silently pass.
- Consumers importing only the public interface should type-check the original
  theorem contracts. Review their assumptions separately for circularity.
- For packaging changes or release closure, inspect wheel/sdist members and test
  the installed public interface outside the checkout. Confirm imported paths
  point to the installed package, not source accidentally found on `sys.path`.
  Inspect for missing sources/docs and included secrets, caches or local notes.

Reuse repository helpers and fixtures. If a temporary verification procedure is
repeated across releases, promote it to a parameterized repository helper rather
than rely on an old `/tmp` path. Test its failure behavior as well as success.
Do not build a new automation framework for a one-off check.

## Keep a canonical verification receipt

Record the following in the PR or a linked canonical evidence record:

- Checked source SHA, source-tree state and scope of the checks. Pre-commit
  verification can be attributed to a commit only after confirming its content
  matches the tested tree. Explain excluded unrelated local files without
  reading, staging or deleting protected notes.
- Root command outcome and logs; exact-SHA CI/check URLs and conclusions;
  relevant independent checks, artifact hashes and installed-package checks.
- Original obligations closed or still open, explicit limitations, and missing
  or failed checks. A worker that never ran is not an independent review.
- PR/merge result and final local/remote SHAs when those actions are authorized.

Freeze content before final verification. A file cannot contain its own final
commit SHA: keep a receipt external to that checked commit, or identify its
checked ancestor and label subsequent changes accurately. Never claim older CI
verified a new head. Any content change requires the applicable gates again;
this convention does not waive the phase's full root command or exact-head CI.

Use README for current status and links; keep detailed evidence in one canonical
place and label historical snapshots as historical. Original-SPEC completion,
phase verification and PR delivery are distinct states. At authorized shipping,
open the draft PR, record its URL, finish the gates on the final branch revision,
then mark ready and merge only under policy. Confirm main's expected contents,
its required CI, successful push and local/`origin/main` agreement before claiming
remote completion. Preserve user changes and never force-push to resolve drift.

At release/final handoff, keep implementation-run provenance canonical in
`FINAL_HANDOFF.md#implementation-run-metadata` and synchronize linked summaries.
Record only exposed model/mode/time/usage fields. Distinguish cumulative goal
accounting from phase-only usage and wall time. Mark unavailable model versions,
token categories and billing as unavailable; any price calculation must be an
explicit API-equivalent estimate with stated rates/assumptions, not an invoice.
Post-completion maintenance must not rewrite frozen original-goal accounting.
