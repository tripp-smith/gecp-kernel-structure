# Final Handoff Report

## Project summary

`gecp-kernel-structure` v1.0.0 delivers a `sorry`-free Lean 4 library and the
`kernelgecp` Python research package. The strongest fermionic result is a
machine-checked separated approximation with rank
`16p(s+1)` and uniform error `2⁻ᵖ` on
`[0,1] × [-2ˢ,2ˢ]`, which is the required
`O(log Λ log(1/ε))` low-rank scale.

Milestone E closes through the SPEC-approved obstruction outcome: exact strict
minor signs do not determine complete-pivot locations, while the refined
`CrossRatioControl` hypothesis is proved sufficient for geometric residual
decay. The repository does not claim that the fermionic residuals satisfy that
hypothesis with a cutoff-uniform contraction factor.

## Phase status

| Phase | State | Delivery |
| --- | --- | --- |
| 0 — bootstrap | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| A — exact GECP core | complete | [PR #4](https://github.com/tripp-smith/gecp-kernel-structure/pull/4) |
| B — PSD baseline | complete | [PR #5](https://github.com/tripp-smith/gecp-kernel-structure/pull/5) |
| C — fermionic structure/census | complete | [PR #6](https://github.com/tripp-smith/gecp-kernel-structure/pull/6) |
| D — separated approximation | complete | [PR #7](https://github.com/tripp-smith/gecp-kernel-structure/pull/7) |
| E — structural GECP | complete | [PR #8](https://github.com/tripp-smith/gecp-kernel-structure/pull/8) |
| F — sparse spectral application | complete | [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1) |
| R — release readiness | verified | [PR #9](https://github.com/tripp-smith/gecp-kernel-structure/pull/9) |

## Completed tasks

- T-001–T-003: Codex-native skills, pinned toolchains, repository skeleton,
  CI, axiom audit, and root verification command.
- T-004–T-005: exact GECP residual/run core plus finite-grid and
  arbitrary-precision numerical GECP.
- T-006: positive-semidefinite residual preservation and canonical recursive
  GECP/pivoted-Cholesky equivalence.
- T-007–T-008: fermionic identities, derivatives, endpoint-resolved grids,
  interval pivot certification, and deterministic census.
- T-009–T-010: explicit exponential-family and fermionic separated
  approximation theorems with an independent high-precision evaluator.
- T-011–T-012: exact structural census, minimized sign-regular obstruction,
  and cross-ratio sufficient condition.
- T-013–T-015: discrete-to-continuous certificate bridge, Green-function
  error transfer, and sparse two-atom application.
- T-016: v1.0.0 metadata, release smoke configuration, artifacts, status
  audit, and final handoff.

## Deferred research

No v1 implementation task is blocked. Open research and optional extensions
are:

- prove a cutoff-uniform `CrossRatioControl` for fermionic GECP, or replace it
  by a more intrinsic determinant/near-max-volume condition;
- formalize the Gimbutas–Marshall–Rokhlin selected-exponential Chebyshev
  construction as an alternate low-rank theorem;
- formalize stronger published pivoted-Cholesky Lipschitz rates;
- add AAA, vector fitting, and multidimensional sparse recovery in a later
  version.

## Public Lean results

The public theorem surface includes:

- exact GECP: `gecp_interpolates_selected_rows`,
  `gecp_interpolates_selected_cols`, `gecp_core_nonsingular`, and
  `gecp_core_det_eq_prod_pivots`;
- PSD baseline: `posDef_gecp_residual_posSemidefinite`,
  `posDef_gecp_residual_abs_le_diag`,
  `posDef_gecp_residual_sup_eq_diag_sup`,
  `gecp_eq_pivotedCholesky_of_posDef`, and
  `posDef_gecp_error_le_fillDistance`;
- fermionic structure: `fermionicKernel_reflection`,
  `fermionicKernel_centered`, `fermionicKernel_continuous`, and the exact
  derivative/bound theorems;
- low rank: `expFamily_separatedApprox` and
  `fermionicKernel_separatedApprox`;
- structural outcome:
  `signRegular_not_sufficient_for_parameterIndependent_pivots`,
  `residualUpdate_le_of_crossRatioControl`, and
  `gecp_error_le_geometric_of_crossRatioControl`;
- application: `grid_sup_le_max_add_lipschitz`,
  `approxPivot_of_grid_certificate`, and
  `greenError_le_kernelError_mul_l1`.

The axiom audit reports only Lean defaults (`propext`, `Quot.sound`, and where
needed `Classical.choice`).

## Python and API results

The package exports stable `FermionicKernel` evaluation, `gecp`,
`gecp_matrix`, pivoted Cholesky, solve-based cross approximation, interval
pivot certificates, deterministic census tools, sparse spectral recovery,
piecewise Chebyshev interpolation, and the Lean-matched
`fermionic_dyadic_taylor_approximation`.

Normal non-convergence is represented by `converged` and `stop_reason`.
Certified-pivot budget exhaustion is never reported as certification, and
high-precision runs retain canonical decimal strings.

## Verification commands run

```bash
./scripts/verify.sh
uv build
uv run python -m kernelgecp.census \
  --config experiments/release-smoke.toml \
  --output /tmp/gecp_release_smoke.jsonl
git status --short
```

## Verification results

- Lean build and public axiom audit pass; no `sorry`, `admit`, or custom axiom
  is present.
- Ruff, formatting, strict mypy, and all 30 pytest cases pass.
- `pip-audit` reports no known vulnerabilities; the local editable project is
  correctly skipped because it is not a PyPI dependency.
- The canonical 21-case census is byte-reproducible with SHA-256
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
- The exact-surrogate dataset has SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The 128-bit release smoke run converges at tolerance `1e-6` with rank 6 for
  cutoff 1 and rank 82 for cutoff `10⁶`. Its two-record JSONL has SHA-256
  `d7b5f53d61eb2ea49cd1ebde944b47a06f931b9d82687c8fa4458ac4c09907d5`.

## Known limitations and non-claims

- No solution of Simons Problem 4.2 is claimed.
- The low-rank theorem is not a GECP theorem.
- Census curves and pivot certificates are evidence, not a continuous GECP
  convergence proof.
- Arbitrary PSD tie-breaking is not claimed equivalent to diagonal pivoting;
  the formal theorem uses the documented diagonal-preferring canonical rule.
- The formal separated construction is dyadic truncated Taylor, not the
  selected-exponential Chebyshev variant.

## Numerical stability and performance

- Kernel evaluation uses sign-dependent overflow-safe formulas.
- GECP supports float64, exact rational surrogates, and arbitrary precision;
  the canonical census uses 128-bit arithmetic.
- Numerical cross approximation uses solves/factorizations rather than an
  explicit inverse.
- Every GECP result records pivots, residual history, conditioning, ties,
  convergence, and stopping metadata.
- The 21-case census converges in every tracked case; sampled rank at
  tolerance `1e-10` grows from 8 at cutoff 1 to 124 at cutoff `10⁶`.

## Security and dependency notes

The Lean revision and Python dependency graph are pinned. Runtime algorithms
need no network access, and the repository stores no credentials, personal
data, or telemetry. Generated environments, caches, and package builds remain
ignored.

## Files and artifacts delivered

- Lean library: `GECPKernelStructure`.
- Python package: `kernelgecp` v1.0.0.
- Canonical data and plot under `experiments/data/`.
- Release artifacts: `dist/kernelgecp-1.0.0-py3-none-any.whl` and
  `dist/kernelgecp-1.0.0.tar.gz`.
- Claim and research documentation: `README.md`, `FINDINGS.md`, `RESEARCH.md`,
  `APPLICATION.md`, `MATHLIB.md`, and `SPEC.md`.

## How to run locally

```bash
uv sync --all-extras
lake exe cache get
uv run python -c "import kernelgecp; print(kernelgecp.__version__)"
```

See `APPLICATION.md` for GECP, certified-pivot, low-rank, and sparse-recovery
examples.

## How to test

```bash
./scripts/verify.sh
```

Focused commands are documented in `AGENTS.md`; a phase is considered verified
only after the root command passes.

## Recommended next research step

Use interval-certified fermionic residuals to search for a cutoff-uniform
cross-ratio factor, then formalize the strongest surviving bound as a lemma
implying `CrossRatioControl`. A failure should be minimized into a certified
counterexample before changing the sufficient condition.
