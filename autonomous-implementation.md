# Autonomous implementation workflow reference

> This file is a preserved workflow input reviewed during repository
> bootstrap. It is not the current project status, release definition of done,
> or claim ledger. Use `SPEC.md` for mathematical objectives, `README.md` and
> `FINAL_HANDOFF.md` for delivered status, and `AGENTS.md` plus the
> repository-scoped skills for active agent instructions. It remains excluded
> from the `kernelgecp` source distribution.

You are working in a Lean 4 / mathlib research repository.

Your task is to upgrade the repository's formalization workflow, project structure, CI, proof-discovery practices, computational verification, and mathlib-readiness. Perform the work directly and autonomously. Do not ask me for clarification or approval. Inspect the repository, infer the existing conventions, preserve all mathematical results, and make the best implementation decisions consistent with the instructions below.

The end state should make this repository easier to extend as a serious Lean research project, easier to audit, and closer to mathlib-quality engineering.

PRIMARY RULES

1. Do not weaken, delete, or restate existing mathematical results merely to make proofs easier.
2. Do not introduce `sorry`, `admit`, unsafe axioms, or new unproved assumptions.
3. Preserve existing public theorem names unless a rename is clearly required for correctness.
4. Keep the repository buildable throughout the work.
5. Prefer reusable lemmas over one-off tactical proofs.
6. Prefer existing mathlib APIs over reimplementing mathematics.
7. Do not add custom tactics or metaprogramming unless ordinary Lean/mathlib techniques genuinely cannot express the needed workflow.
8. Keep dependencies minimal.
9. Pin or respect the repository's existing Lean and mathlib versions.
10. Update documentation whenever implementation changes affect the documented workflow.
11. Run the complete verification suite before finishing.
12. If a requested tool or command has changed since these instructions were written, consult its current official documentation and implement the equivalent modern setup rather than stopping.

==================================================
PHASE 1: AUDIT THE CURRENT REPOSITORY
==================================================

Before modifying anything:

1. Read:
   - README.md
   - SPEC.md
   - FINDINGS.md
   - RESEARCH.md
   - APPLICATION.md
   - MATHLIB.md
   - lakefile.lean or lakefile.toml
   - lean-toolchain
   - .github/workflows/*
   - scripts/*
   - all top-level Lean modules
   - tests and Python package metadata, if present.

2. Identify:
   - headline public theorems;
   - theorem dependency structure;
   - existing computational checkers;
   - existing CI guarantees;
   - existing `#print axioms` checks;
   - places where theorem-search or generic infrastructure appears to have been recreated locally;
   - reusable lemmas that might ultimately belong in mathlib.

3. Do not spend time rewriting stable proofs simply for cosmetic reasons.

==================================================
PHASE 2: ADD A FORMALIZATION BLUEPRINT
==================================================

Introduce LeanProject / leanblueprint-style project documentation using the current recommended tooling.

Goal:

mathematical specification
    ->
dependency graph
    ->
Lean declaration
    ->
machine-checked completion status

Use the existing SPEC.md as the source of truth for mathematical intent.

Create an appropriate blueprint structure, normally something like:

blueprint/
    src/
        content.tex
        macros.tex
        ...
    web/
    print/

or the structure required by the current leanblueprint release.

The blueprint should document the principal theorem dependency graph, not every helper lemma.

At minimum include nodes for:

1. structural definitions;
2. central algebraic reduction(s);
3. main positive theorem(s);
4. counterexample or negative theorem(s), when applicable;
5. minimality/sharpness results;
6. application bridge results;
7. any major open conjecture still under investigation.

Each completed mathematical node should reference the actual Lean declaration implementing it.

Configure declaration checking so stale theorem names or nonexistent declarations cause verification failure where supported.

Do not duplicate the entire README in the blueprint.

Update README.md with a short "Formalization blueprint" section explaining:
- what the blueprint is;
- where it lives;
- how to build/view it;
- how declaration completion is checked.

If the blueprint tooling requires Python or another lightweight development dependency, document it clearly and add it only to development configuration.

==================================================
PHASE 3: STANDARDIZE PROOF DISCOVERY
==================================================

Create contributor/agent instructions that establish a search-first proof workflow.

Add this to the repository's existing agent rules, Cursor rules, AGENTS.md, CONTRIBUTING.md, or another appropriate contributor file.

The required workflow is:

BEFORE writing a bespoke proof of a reusable mathematical fact:

1. Search mathlib by theorem/type shape.
2. Try relevant editor or Lean discovery tools such as:
   - Loogle;
   - `#check`;
   - `#find`;
   - `exact?`;
   - `apply?`;
   - `library_search`;
   - `simp?`;
   - `aesop?`;
   - `grind?`.
3. Search by mathematical type structure rather than English description alone.
4. Only introduce a new generic lemma after confirming no suitable mathlib result exists.

Example mindset:

Instead of searching for:
    "trace invariant under permutation"

search for a shape resembling:
    Matrix.trace (A.reindex e e) = Matrix.trace A

For completed proofs, prefer stable explicit proofs over leaving exploratory automation in place.

When `grind?`, `simp?`, `exact?`, or similar tools propose a minimal proof, use the resulting explicit theorem list where reasonable.

==================================================
PHASE 4: USE AUTOMATION DELIBERATELY
==================================================

Audit suitable local proofs for opportunities to simplify proof search with modern Lean automation, but do not churn working proofs unnecessarily.

Establish these conventions:

1. `simp` for canonical rewriting.
2. `norm_num` for explicit arithmetic.
3. `ring` / `ring_nf` for polynomial identities.
4. `linarith` / `nlinarith` for arithmetic consequences.
5. `aesop` for structural proof search.
6. `grind` for appropriate first-order/algebraic goals.
7. `native_decide` only for finite decidable computational certification.
8. Explicit mathematical lemmas for central conceptual steps.

During proof development prefer:

    grind?

then replace it, when practical, with a constrained form such as:

    grind only [lemma1, lemma2, ...]

or another transparent proof that records the important dependencies.

Do not replace readable conceptual arguments with huge opaque automation calls merely to reduce line count.

==================================================
PHASE 5: IMPROVE API DESIGN
==================================================

Review the formal API for definitions or lemmas that are too specialized.

Apply this rule:

    Prove reusable results one abstraction level above the theorem that first required them.

Examples of desirable abstractions include:

- reindex/permutation invariance;
- principal-submatrix identities;
- block matrix identities;
- positivity preservation;
- determinant identities;
- finite expectation identities;
- norm equivalences under structural assumptions;
- generic combinatorial lemmas.

If the repository already contains such reusable infrastructure, preserve it.

Do not generalize aggressively when doing so introduces unnecessary typeclass complexity.

Prefer:
- small stable definitions;
- mathematically natural predicates;
- lemmas attached to predicates through namespaces;
- theorem names consistent with mathlib conventions;
- `[simp]` attributes only on genuinely canonical simplifications.

Avoid definitions whose names claim more mathematically than they provide.

==================================================
PHASE 6: UPGRADE CI AND LOCAL VERIFICATION
==================================================

Upgrade CI so the repository checks more than "it builds."

Implement the strongest appropriate subset of:

    lake build
    lake test
    lake lint
    rg '\bsorry\b' --glob '*.lean'
    python -m pytest

Use the actual supported commands for the pinned Lean/mathlib version.

If `lake lint` or `lake test` requires project setup, add the necessary lightweight configuration.

Retain the explicit no-`sorry` check even if another tool indirectly detects it.

Add a single local verification script, preferably:

    scripts/verify.sh

that reproduces the important CI checks.

It should:
- use `set -euo pipefail`;
- fail immediately on errors;
- print concise phase labels;
- execute Lean build/lint/test as supported;
- reject `sorry`;
- run Python tests if a Python package exists;
- run blueprint declaration checks if available.

Document:

    ./scripts/verify.sh

as the canonical pre-commit/pre-PR verification command.

==================================================
PHASE 7: AXIOM AUDITING
==================================================

For each headline theorem, audit its axioms.

Create a Lean audit file if one does not already exist, such as:

    AxiomAudit.lean

or integrate the checks into the existing theorem entrypoint.

Use:

    #print axioms theorem_name

for the main public results.

The expected dependencies should be limited to standard Lean foundations such as those already accepted by mathlib, for example:

    propext
    Classical.choice
    Quot.sound

Do not assert in documentation that a theorem is axiom-clean unless this is actually checked.

Where practical, make the audit part of CI or the local verification process.

==================================================
PHASE 8: FORMALIZE THE COMPUTATIONAL RESEARCH LOOP
==================================================

Standardize the repository's method for turning computational observations into theorems.

Document and support this pipeline:

    floating-point search
        ->
    high-precision confirmation
        ->
    rational/algebraic reconstruction
        ->
    exact independent computation
        ->
    Lean kernel proof

Create or improve RESEARCH.md so this methodology is explicit.

Any computationally discovered conjecture must be labeled as one of:

    PROVED
    COMPUTATIONALLY CERTIFIED ON FINITE INSTANCES
    NUMERICALLY OBSERVED
    CONJECTURED
    DISPROVED

Never blur these categories.

If the repository contains numerical counterexample searches, preserve enough metadata to reproduce them:
- matrix/problem parameters;
- random seed;
- precision;
- objective value;
- witness;
- code revision where appropriate.

==================================================
PHASE 9: PROPERTY-BASED AND RANDOMIZED TESTING
==================================================

If the repository has a Python application layer, strengthen it with property-style tests.

Use deterministic seeds.

Test generic identities on many small random instances where mathematically appropriate.

Examples:

- formal identity vs direct dense computation;
- Schur update vs recomputation;
- invariant under permutation/sign switching;
- exact residual vs direct matrix norm;
- structural predicate on dense vs sparse representation;
- positive theorem on random valid instances;
- known counterexample remains a counterexample;
- approximate algorithms never claim exact guarantees.

Use Hypothesis only if it provides real value and does not add unnecessary complexity. Otherwise deterministic randomized pytest loops are sufficient.

Keep exact certified examples separate from floating-point tests.

==================================================
PHASE 10: EXACT SMALL-INSTANCE CHECKERS
==================================================

Whenever appropriate, preserve or introduce exact finite checkers over:

    ℚ

or another exact algebraic representation.

Use them for:
- minimal counterexamples;
- small exhaustive censuses;
- determinant identities;
- pivot sequences;
- finite combinatorial inequalities.

Prefer:
- exact rational arithmetic;
- Cramer's rule or exact determinant formulas for very small matrices;
- `native_decide` for finite decidable checks.

Do not use `native_decide` as a substitute for a general mathematical theorem.

Clearly distinguish:

    exact finite verification

from:

    general symbolic proof.

==================================================
PHASE 11: MATHLIB-READINESS
==================================================

Review all generic helper lemmas and update MATHLIB.md.

Classify potential upstream contributions into:

A. Good immediate mathlib candidates
   Generic results independent of the project.

B. Needs API cleanup first
   Useful results currently stated with project-specific definitions.

C. Project-specific
   Results that should remain in this repository.

For each candidate record:

- theorem name;
- mathematical statement;
- existing mathlib namespace where it likely belongs;
- dependencies;
- why it is generally reusable;
- whether an equivalent theorem may already exist.

Do not submit upstream PRs automatically unless repository policy explicitly instructs you to do so.

The purpose is to make future upstreaming easy.

==================================================
PHASE 12: DOCUMENT PROJECT MATURITY
==================================================

Update README.md so a new reader can immediately distinguish:

1. the original open problem;
2. prior mathematical results;
3. what this repository proves;
4. what is formally machine checked;
5. what remains conjectural;
6. what the Python/application layer does;
7. how to reproduce verification.

Include a compact status table if useful.

Avoid overstating novelty.

Where a result reproduces prior literature, call it a formalization or machine-checked proof of that prior result.

Where the repository proves something genuinely new, state it precisely.

==================================================
PHASE 13: REPOSITORY ORGANIZATION
==================================================

If documentation is scattered or confusing, rationalize it without deleting useful content.

Preferred responsibilities:

README.md
    Project overview, headline results, verification status.

SPEC.md
    Mathematical scope and formal success criteria.

FINDINGS.md
    Accessible explanation of proved results.

RESEARCH.md
    Conjectures, experiments, negative results, research log.

APPLICATION.md
    Executable/Python interface and algorithm semantics.

MATHLIB.md
    Reusable upstream candidates.

CONTRIBUTING.md or AGENTS.md
    Lean proof workflow and repository engineering conventions.

blueprint/
    Mathematical dependency graph connected to Lean declarations.

Do not move files merely to make the tree look cleaner if existing links or workflows rely on them.

==================================================
PHASE 14: AGENT-SPECIFIC WORKFLOW RULES
==================================================

Add repository instructions for future Codex/Cursor agents.

They should include:

Before proving a Lean fact:

    1. Inspect the exact goal.
    2. Search mathlib.
    3. Try `exact?`, `apply?`, `simp?`, `aesop?`, `grind?`.
    4. Look for an existing project abstraction.
    5. Only then write a bespoke lemma.

Before adding a definition:

    1. Search mathlib for an existing concept.
    2. Check whether a weaker/more reusable definition suffices.
    3. Avoid encoding theorem assumptions into definitions unnecessarily.

Before claiming a mathematical result:

    1. Identify the exact Lean theorem.
    2. Confirm it builds.
    3. Check its axioms if it is headline-worthy.
    4. Distinguish theorem from finite computation or numerical evidence.
    5. Update the blueprint and documentation.

Before finishing any task:

    ./scripts/verify.sh

==================================================
PHASE 15: DO NOT OVERUSE METAPROGRAMMING
==================================================

Do not build custom tactics merely to automate repetitive project proofs unless there is a demonstrated recurring pattern that cannot be handled cleanly with:

    simp
    aesop
    grind
    norm_num
    ring
    linarith
    native_decide

If genuine recurring proof patterns emerge, prefer:
- helper lemmas;
- carefully chosen `[simp]` lemmas;
- local Aesop rules;

before writing custom metaprograms.

==================================================
PHASE 16: FINAL VERIFICATION
==================================================

Before finishing:

1. Run the entire verification script.
2. Confirm no Lean source contains `sorry`.
3. Confirm all headline theorems still exist.
4. Confirm blueprint declaration references resolve.
5. Confirm Python tests pass.
6. Confirm lint/test/build status.
7. Check git diff for accidental generated files, caches, binaries, large artifacts, or unrelated edits.
8. Ensure generated blueprint/build output is ignored appropriately unless it is intentionally committed.
9. Review documentation for stale commands and theorem names.
10. Verify no existing guarantee or theorem was silently weakened.

==================================================
DELIVERABLE
==================================================

Implement all appropriate changes directly.

At completion provide a concise report containing:

1. files added;
2. files materially changed;
3. CI/workflow improvements;
4. blueprint status;
5. proof-discovery rules added;
6. tests added;
7. mathlib candidates identified;
8. full verification commands run and their results;
9. anything deliberately not changed and why.

Do not leave a list of tasks for me to finish manually if they can reasonably be completed in the repository.

If an optional tool cannot be integrated due to incompatibility with the pinned Lean version, preserve the rest of the implementation, document the incompatibility precisely, and use the closest supported equivalent.
