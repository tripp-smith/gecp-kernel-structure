---
name: phase-cadence
description: Execute one bounded GECP Kernel Structure phase from specification through implementation, verification, documentation, and a draft pull request. Use when asked to execute or close a phase, milestone, theorem contract, research obstruction, or release-readiness step in this repository.
---

# GECP phase cadence

Deliver exactly one scoped phase from `SPEC.md`. Do not reopen a completed
phase unless the user explicitly asks. Keep theorem claims, experiments, and
research hypotheses visibly distinct.

This workflow does not grant authority to change files or publish work. For a
read-only assessment, use the relevant checks and report recommendations;
write contracts, push branches, open PRs or merge only within user authorization.

## 1. Specify

Before implementation, write the phase contract in `RESEARCH.md` and set the
README phase to `in progress`.

- Map the SPEC milestone to exact Lean identifiers and Python artifacts.
- Record statements, files, independent checks, and explicit non-claims.
- Identify dependencies on earlier verified results.
- Maintain the affected original-SPEC obligations in `RESEARCH.md`, linking
  existing audit records rather than creating a competing ledger. For each,
  retain the original scope/quantifiers, assumptions, delivered evidence and
  acceptance test. Distinguish mandatory outcomes from proposed proof routes
  and conditional extensions; a weaker result does not close the original item.
- At closure, reread the original SPEC text, not only later status labels.
  Mechanism-specific requirements need evidence of that mechanism: fixed-grid
  runs, for example, do not establish residual-driven adaptivity. Progress
  estimates describe remaining obligations and delivery gates, not counts of
  lemmas, tests or PRs, nor confidence that a conjecture is true.
- Use branch `codex/phase-<id>-<slug>` from synchronized `main`.
- Do not silently broaden a theorem because a planned proof route failed.

When an open proof route remains uncertain, use `$lean-research-loop`, including
its route-level checkpoint before investing in more auxiliary abstractions.

## 2. Implement

- Keep every library module free of `sorry` and undeclared axioms.
- Prefer exact finite algebra before analytic or asymptotic corollaries.
- Re-export public Lean results through `GECPKernelStructure/Theorems.lean`.
- Keep numerical factorization solve-based; do not form explicit inverses.
- Return convergence, stop reason, conditioning, and certification status.
- Store canonical high-precision quantities as decimal strings.
- Preserve deterministic tie-breaking and record near-tie diagnostics.

For numerical certificates or canonical experiments, read
[numerical evidence](references/numerical-evidence.md) before implementation.
It covers arithmetic boundary attacks, validation isolation and frozen provenance;
it is not a checklist for unrelated documentation or Lean-only changes.

If the planned theorem is not proved, deliver only a precisely named weaker
theorem or a certified obstruction allowed by `SPEC.md`. Never create a
tautological certificate wrapper merely to claim a milestone name.

## 3. Verify

Run the phase-focused checks while developing, then the full root command:

```bash
./scripts/verify.sh
```

For each public structural theorem, add it to `scripts/AxiomAudit.lean` and
confirm that only Lean defaults (`propext`, `Classical.choice`, `Quot.sound`)
appear. Independently recompute claimed closed forms, explicit constructions,
and counterexample witnesses. Certification-budget exhaustion is an
uncertified result, not success.

Execute the fail-closed wrapper `scripts/check_axioms.sh`. As original-SPEC public
theorem interfaces become available, add expected-type consumers to
`scripts/CompletionContract.lean`
importing only `GECPKernelStructure.Theorems`. Check the required quantifiers and
absence of extra assumptions, not just that the exported name compiles.

When adding/changing verification gates or closing a release, read
[verification and delivery evidence](references/verification-and-delivery.md).
Require independent computations to differ meaningfully from the implementation;
an alias or the same algorithm at higher precision is a consistency check, not
an independent oracle. Distinguish primary review, independent review and checks
that could not run; failed worker calls are not review evidence.

## 4. Document

Update the same logical change in all affected documents:

- `README.md`: current phase, state, verification, claim level, outputs, PR.
- `FINDINGS.md`: proved results, observations, conjectures, and non-claims.
- `RESEARCH.md`: contract, failed hypotheses, witnesses, next decision.
- `APPLICATION.md`: stable public Python behavior and limitations.
- `MATHLIB.md`: genuinely reusable abstractions only.
- `SPEC.md`: status annotations and delivered names only; never rewrite goals.

Remove stale status text. Set the phase to `verified` only after the final
content passes the full local command and exact-revision CI. Record the checked
revision explicitly; a later metadata commit is not that same checked head.
Mark delivery `complete` only after confirmed merge and remote-state checks;
recording this at the next phase opening still requires that evidence.
Do not reopen completed mathematical scope for workflow-only maintenance.

## 5. Ship

Keep commits logical and intentional. Push the phase branch and open a draft
PR against `main`. Record the PR in README. Merge only when the pre-merge checks
below pass, repository policy allows it and the user authorized merging.
Keep one canonical revision-specific verification receipt in the PR or a linked
record. Summaries reference it instead of copying changing hashes and counts.
Use the delivery reference for exact-head checks and post-merge synchronization;
a receipt never substitutes for a required gate.

## Definition of done

Pre-merge acceptance:

- The contract has named outputs and non-claims.
- Original-SPEC obligations are matched to actual evidence; public expected-type
  consumers cover changed theorem interfaces where applicable.
- Required Lean and Python artifacts are real, tested, and documented.
- `./scripts/verify.sh` passes with no `sorry`.
- Public theorem axioms satisfy the repository policy.
- Exact, high-precision, failure, and certificate paths relevant to the phase
  have evidence.
- README and supporting documents match delivered claims.
- For authorized PR delivery, the branch is pushed and a draft PR exists.

Delivery completion, when merging is authorized, additionally requires confirmed
merge, required main checks and local/`origin/main` agreement. These post-merge
checks are not prerequisites for performing the merge itself.
