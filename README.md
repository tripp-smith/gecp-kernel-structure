# GECP Kernel Structure

Lean 4 formalization and reproducible Python research package for continuous
Gaussian elimination with complete pivoting (GECP), its positive-definite
pivoted-Cholesky baseline, and the fermionic DLR kernel.

> **Current phase:** R — Release-readiness audit<br>
> **Phase state:** blocked (research)<br>
> **Last verification:** `./scripts/verify.sh` and `uv build` passed 2026-08-15<br>
> **Verification command:** `./scripts/verify.sh`<br>
> **Delivery:** [draft PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1)<br>
> **Claim level:** core/PSD/kernel/Green theorems proved; D and E remain research-blocked<br>
> **Workflow:** [`$phase-cadence`](.agents/skills/phase-cadence/SKILL.md)

## Status

| Phase | SPEC milestone | State | Named outputs | Verification evidence | Delivery |
| --- | --- | --- | --- | --- | --- |
| 0 | Bootstrap | verified | two skills, pinned toolchains, CI, package roots | skill validation; focused builds | draft PR #1 |
| A | Exact GECP core | verified | interpolation; LDU determinant/nonsingularity | Lean build; axiom audit passed | draft PR #1 |
| B | PSD baseline | verified | Schur PSD; diagonal domination; canonical one-step rule; fill bound | Lean/Python focused checks | draft PR #1 |
| C | Fermionic structure/census | verified | reflection, centered form, continuity; stable API; 21-case grid census | Lean/Python checks; census SHA below | draft PR #1 |
| D | Separated approximation | blocked (research) | arithmetic scaffold; composite interpolant | numerical tests only; explicit uniform Lean theorem missing | draft PR #1 |
| E | Structural GECP | blocked (research) | 21 exact surrogate cases; variable pivot-order obstruction | exact rational SHA below; no rate theorem | draft PR #1 |
| F | Sparse-\(\rho\) application | verified | `greenError_le_kernelError_mul_l1`; OMP/refinement | Lean build; two-atom regression | draft PR #1 |
| R | Release readiness | blocked (research) | auditable partial handoff; wheel and sdist | full verification passed | draft PR #1 |

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

- Finite-grid census: 21 cutoff/tolerance cases, byte-reproducible, SHA-256
  `3381573b39addf3d5fb64c58f35f3923735b859a1d41953edbba947f88959676`.
- Exact geometric surrogates: every square minor for sizes 2–8 at
  `q ∈ {1/2, 2/3, 3/4}`, SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The exact surrogate pivot order changes with `q`; a universal
  parameter-independent pivot-order invariant is therefore not viable in that
  form.

The census currently uses the float64 finite-grid backend. The requested
128-bit adaptive/certified census is not represented by this data and remains
open work.
