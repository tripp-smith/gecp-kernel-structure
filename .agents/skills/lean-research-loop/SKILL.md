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
- Preserve the original SPEC target and its expected public type alongside the
  proposed lemma. Record its downstream consumer: which assumption, missing
  bridge or explicit bound will this result discharge? Added hypotheses or a
  weaker conclusion must leave the original obligation visibly open.

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

### Check the route, not just the tactic

At phase planning, reconsider a route when the same central bridge survives
successive phases or the next abstraction has no concrete target consumer.
Successful auxiliary lemmas alone do not establish progress on the original
goal. Record the bottleneck and compare credible alternatives before extending
the route; this is a review point, not automatic abandonment after a fixed count.

- Search existing repository results as well as mathlib for ingredients before
  building new machinery. Check applicability against the actual hypotheses.
- Consider a change of representation, normalization/sign gauge, local-to-global
  counting or a backward finite witness when it removes the difficult bridge.
  Here, common residual sign gauges, local panel determinant budgets and
  backwards-transported minor witnesses supplied the eventual rate route; these
  are examples of a search strategy, not requirements for every theorem.
- Give an exploratory route an observable success/exit criterion. Prefer a
  discriminating exact case or critical lemma over a long speculative chain.
  Keep reusable lemmas and minimized obstructions when switching routes.
- State whether a counterexample refutes the target or only a stronger sufficient
  condition/proof mechanism. Do not turn failure of one route into a negative
  answer to the original problem, or expand scope without user direction.

## 4. Audit the result

- Build the focused module.
- Add the public result to `scripts/AxiomAudit.lean`.
- Run `scripts/check_axioms.sh`; update expected-type consumers in
  `scripts/CompletionContract.lean` for original-SPEC interfaces. A kernel-checked
  proof of a weaker type does not establish the requested theorem.
- Run the exact or high-precision independent check appropriate to the claim.
- Inspect assumptions for circularity or a hidden restatement of the conclusion.

## 5. Decide

Return exactly one of:

- proved theorem with exact public identifier;
- weaker theorem with the lost strength stated explicitly;
- certified counterexample or obstruction;
- `blocked (research)` with the precise missing lemma and evidence gathered.

Experiments can choose the next theorem; they cannot close the theorem.
