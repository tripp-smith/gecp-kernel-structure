# GECP Kernel Structure

Lean 4 formalization and reproducible Python research package for continuous
Gaussian elimination with complete pivoting (GECP), its positive-definite
pivoted-Cholesky baseline, and the fermionic DLR kernel.

> **Current phase:** D — Constructive separated approximation<br>
> **Phase state:** verified<br>
> **Last verification:** `./scripts/verify.sh` passed 2026-08-15 (29 tests)<br>
> **Verification command:** `./scripts/verify.sh`<br>
> **Delivery:** [Phase D draft PR #7](https://github.com/tripp-smith/gecp-kernel-structure/pull/7)<br>
> **Claim level:** explicit fermionic low-rank theorem proved; no GECP rate is claimed<br>
> **Workflow:** [`$phase-cadence`](.agents/skills/phase-cadence/SKILL.md)

## Status

| Phase | SPEC milestone | State | Named outputs | Verification evidence | Delivery |
| --- | --- | --- | --- | --- | --- |
| 0 | Bootstrap | complete | two skills, pinned toolchains, CI, package roots | skill validation; focused builds | merged PR #1 |
| A | Exact GECP core | complete | interpolation; run-selected core determinant/nonsingularity | Lean build; axiom audit; exact two-pivot check | merged PR #4 |
| B | PSD baseline | complete | same-domain Schur PSD; diagonal domination; canonical recursive GECP/Cholesky traces; fill bound | Lean axiom audit; tie regression | merged PR #5 |
| C | Fermionic structure/census | complete | exact derivatives/bounds; 128-bit GECP; interval-certified adaptive pivots; endpoint-resolved 21-case census | axiom audit; 27 Python tests; byte-identical census repeat | merged PR #6 |
| D | Separated approximation | verified | `expFamily_separatedApprox`; `fermionicKernel_separatedApprox`; verified dyadic evaluator | axiom audit; 29 tests; full verification | draft PR #7 |
| E | Structural GECP | blocked (research) | 21 exact surrogate cases; variable pivot-order obstruction | exact rational SHA below; no rate theorem | baseline merged in PR #1 |
| F | Sparse-\(\rho\) application | complete | `greenError_le_kernelError_mul_l1`; OMP/refinement | Lean build; two-atom regression | merged PR #1 |
| R | Release readiness | blocked (research) | auditable partial handoff; wheel and sdist | full verification passed | baseline merged in PR #1 |

Allowed states are `planned`, `in progress`, `verified`, `complete`, and
`blocked (research)`. A phase becomes `complete` only after its verified change
is merged and the next phase records that fact.

## Quick start

Prerequisites are [uv](https://docs.astral.sh/uv/) and the Lean toolchain
manager `elan`.

```bash
uv sync --all-extras
lake exe cache get
./scripts/verify.sh
```

The Python package is named `kernelgecp` and the Lean library is
`GECPKernelStructure`. See [APPLICATION.md](APPLICATION.md) for API usage,
[FINDINGS.md](FINDINGS.md) for claim status, and [RESEARCH.md](RESEARCH.md) for
active theorem contracts.

## Claim discipline

Documentation uses four distinct labels: **proved**, **observed**,
**conjectured**, and **not claimed**. In particular, this repository does not
currently claim the target continuous fermionic GECP rate or a solution of
Simons Problem 4.2. Numerical evidence never substitutes for a Lean theorem or
a separately certified analytic argument.

## Delivered evidence

- Endpoint-resolved 128-bit finite-grid census: all 21 cutoff/tolerance cases
  converged, byte-reproducible, SHA-256
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
- Exact geometric surrogates: every square minor for sizes 2–8 at
  `q ∈ {1/2, 2/3, 3/4}`, SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The exact surrogate pivot order changes with `q`; a universal
  parameter-independent pivot-order invariant is therefore not viable in that
  form.

At tolerance `1e-10`, the sampled rank grows from 8 at cutoff 1 to 124 at
cutoff `1e6`. This is numerical finite-grid evidence, not a continuous GECP
rate theorem. Continuous interval-certified pivoting is a separate supported
mode and is tested on synthetic and fermionic cases; it is not substituted for
the canonical finite-grid census or a Lean convergence proof.
