# Repository agent instructions

Use the repository skill at `.agents/skills/phase-cadence/SKILL.md` for every
milestone-sized change.

When an open theorem requires iterative proof search or exact obstruction
work, compose it with `.agents/skills/lean-research-loop/SKILL.md`.

## Durable commands

```bash
uv sync --all-extras
lake exe cache get
./scripts/verify.sh
```

Focused checks may be used during development, but a phase is not verified
until the root command passes.

## Formalization rules

- Do not commit `sorry`, `admit`, custom axioms, or theorem statements whose
  assumptions already contain the conclusion.
- Public results are re-exported by `GECPKernelStructure/Theorems.lean` and
  audited by `scripts/check_axioms.sh`.
- Search mathlib before introducing a local abstraction. Record credible
  upstream candidates in `MATHLIB.md`.
- Prefer explicit finite and natural-number bounds. Derive asymptotic wording
  only after the explicit theorem exists.

## Numerical rules

- Numerical cross approximation uses factorizations/solves, never an explicit
  inverse.
- Results expose convergence, stopping, conditioning, precision, tie, and
  certification metadata.
- Experiments are deterministic for a fixed configuration and revision.
- Floating-point output is evidence about behavior, not proof of a theorem.

## Documentation cadence

Every phase updates README status plus the affected claim documents in the
same logical change. `SPEC.md` is authoritative: append status annotations and
delivered identifiers, but do not rewrite its mathematical goals.
