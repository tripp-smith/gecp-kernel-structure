---
name: phase-cadence
description: Execute one bounded GECP Kernel Structure phase from specification through implementation, verification, documentation, and a draft pull request. Use when asked to execute or close a phase, milestone, theorem contract, research obstruction, or release-readiness step in this repository.
---

# GECP phase cadence

Deliver exactly one scoped phase from `SPEC.md`. Do not reopen a completed
phase unless the user explicitly asks. Keep theorem claims, experiments, and
research hypotheses visibly distinct.

## 1. Specify

Before implementation, write the phase contract in `RESEARCH.md` and set the
README phase to `in progress`.

- Map the SPEC milestone to exact Lean identifiers and Python artifacts.
- Record statements, files, independent checks, and explicit non-claims.
- Identify dependencies on earlier verified results.
- Use branch `codex/phase-<id>-<slug>` from synchronized `main`.
- Do not silently broaden a theorem because a planned proof route failed.

## 2. Implement

- Keep every library module free of `sorry` and undeclared axioms.
- Prefer exact finite algebra before analytic or asymptotic corollaries.
- Re-export public Lean results through `GECPKernelStructure/Theorems.lean`.
- Keep numerical factorization solve-based; do not form explicit inverses.
- Return convergence, stop reason, conditioning, and certification status.
- Store canonical high-precision quantities as decimal strings.
- Preserve deterministic tie-breaking and record near-tie diagnostics.

If the planned theorem is not proved, deliver only a precisely named weaker
theorem or a certified obstruction allowed by `SPEC.md`. Never create a
tautological certificate wrapper merely to claim a milestone name.

## 3. Verify

Run the phase-focused checks while developing, then the full root command:

```bash
./scripts/verify.sh
```

For each public structural theorem, add it to `scripts/check_axioms.sh` and
confirm that only Lean defaults (`propext`, `Classical.choice`, `Quot.sound`)
appear. Independently recompute claimed closed forms, explicit constructions,
and counterexample witnesses. Certification-budget exhaustion is an
uncertified result, not success.

## 4. Document

Update the same logical change in all affected documents:

- `README.md`: current phase, state, verification, claim level, outputs, PR.
- `FINDINGS.md`: proved results, observations, conjectures, and non-claims.
- `RESEARCH.md`: contract, failed hypotheses, witnesses, next decision.
- `APPLICATION.md`: stable public Python behavior and limitations.
- `MATHLIB.md`: genuinely reusable abstractions only.
- `SPEC.md`: status annotations and delivered names only; never rewrite goals.

Remove stale status text. Set the phase to `verified` only after the final
commit passes the full local command and CI. Mark it `complete` at the start of
the next phase after merge.

## 5. Ship

Keep commits logical and intentional. Push the phase branch and open a draft
PR against `main`. Record the PR in README. Merge only when the definition of
done is satisfied and repository policy allows it.

## Definition of done

- The contract has named outputs and non-claims.
- Required Lean and Python artifacts are real, tested, and documented.
- `./scripts/verify.sh` passes with no `sorry`.
- Public theorem axioms satisfy the repository policy.
- Exact, high-precision, failure, and certificate paths relevant to the phase
  have evidence.
- README and supporting documents match delivered claims.
- The branch is pushed and a draft PR exists.
