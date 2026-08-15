---
name: lean-research-loop
description: Run a disciplined proof-or-obstruction loop for an open Lean theorem in this repository. Use when a phase is blocked on a formal claim, a structural invariant needs exact testing, a proof attempt repeatedly fails, or Milestone D or E needs a rigorous next decision.
---

# Lean research loop

Use this inside a `$phase-cadence` implementation stage when the mathematical
route is uncertain. Work on one exact statement at a time.

## 1. Freeze the claim

- Write the exact proposed theorem, assumptions, namespace, and non-claims in
  `RESEARCH.md`.
- Search mathlib for the needed abstraction before proving a local duplicate.
- Separate algebraic, analytic, and asymptotic obligations.

## 2. Attack the smallest exact case

- Prove the finite algebraic core first.
- Test rational surrogates with exact arithmetic, never float output.
- If the claim fails, minimize the witness and independently recompute it.
- Record failed hypotheses; do not delete them from research history.

## 3. Escalate proof search deliberately

Try, in order:

1. simplification and existing named lemmas;
2. `ring`, `linarith`, `nlinarith`, or `positivity` for closed algebra;
3. a helper lemma exposing the missing invariant;
4. a stronger reusable abstraction already present in mathlib;
5. a revised theorem with explicit additional hypotheses.

After three failures with the same missing fact, stop tactic churn and restate
that missing fact as the next research contract.

## 4. Audit the result

- Build the focused module.
- Add the public result to `scripts/AxiomAudit.lean`.
- Run the exact or high-precision independent check appropriate to the claim.
- Inspect assumptions for circularity or a hidden restatement of the conclusion.

## 5. Decide

Return exactly one of:

- proved theorem with exact public identifier;
- weaker theorem with the lost strength stated explicitly;
- certified counterexample or obstruction;
- `blocked (research)` with the precise missing lemma and evidence gathered.

Experiments can choose the next theorem; they cannot close the theorem.
