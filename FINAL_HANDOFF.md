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

The post-v1 Python application phase demonstrates how those components work
together on realistic synthetic data: one universal fermionic basis compresses
two qualitatively different continuous spectra, sparse routines represent
clean and noisy Green functions under explicit dictionary assumptions, and
the PSD baseline selects spatial covariance landmarks. These are application
observations with held-out errors, not new formal convergence or inverse-
problem claims.

Post-v1 Phase G adds a rigorous rate obstruction and the first strict
fermionic block theorem. Lean proves that the first two complete pivots are
the two reflected cutoff corners for every positive cutoff, that no universal
factor below one can contract the first step, and that for `0<Λ≤1` the actual
two-step residual is at most one half of the initial complete-pivot magnitude.
This verifies the base block but does not supply the dyadic scaling lemma
needed for Conjecture G1.

## Implementation run metadata

This is a provenance record for the Codex implementation goal, captured at
handoff on 2026-08-15 at 17:39 EDT (UTC-04:00).

| Field | Recorded value |
| --- | --- |
| Codex surface | Codex desktop app |
| Declared model identity | OpenAI Codex based on the GPT-5 family |
| Exact serving model/version | Not exposed to the task runtime; no model slug or snapshot can be verified |
| Collaboration mode | Default mode; primary/root agent |
| Reasoning effort | Not exposed to the task runtime |
| Workspace execution mode | Unrestricted local filesystem/network access; approval policy `never` |
| Goal state | Complete |
| Goal-accounting usage at completion | 1,801,132 aggregate tokens and 10,747 seconds (2:59:07) |
| Usage scope | The implementation goal through v1 completion; this post-completion metadata update is not included because its goal counter is no longer available |

The aggregate token counter does not split input, cached input, reasoning, or
output tokens and is not an OpenAI invoice. The runtime also does not expose
the exact serving model. For an order-of-magnitude comparison only, the
current [OpenAI API pricing](https://openai.com/api/pricing/) lists standard
GPT-5.6 Sol text rates of $5 per million input tokens and $30 per million
output tokens. Applying those rates to the 1.801132-million-token aggregate
gives:

- an illustrative 80% input / 20% output estimate of **$18.01**;
- an uncached standard-rate sensitivity envelope of **$9.01–$54.03** if the
  whole aggregate were priced at the input or output rate, respectively.

These figures are API-equivalent estimates, not the actual Codex subscription,
credit, or internal compute cost. They exclude cache discounts, long-context
premiums, service-tier adjustments, tool-call charges, and any usage that was
not represented by the completed goal counter.

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
| R — release readiness | complete | [PR #9](https://github.com/tripp-smith/gecp-kernel-structure/pull/9), [closure PR #10](https://github.com/tripp-smith/gecp-kernel-structure/pull/10) |
| S — realistic synthetic applications | complete | implementation `6404311`; canonical evidence delivered directly to `main` |
| G — fermionic GECP rate | verified partial result; G1 open | draft [PR #11](https://github.com/tripp-smith/gecp-kernel-structure/pull/11) |

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
- T-017: typed synthetic-application API, Hubbard-like and gapped Green-
  function compression, clean/noisy sparse recovery, PSD covariance landmark
  selection, deterministic JSON, summary plot, and application documentation.
- T-018: variable-factor/block contraction infrastructure, exact one-step
  obstruction, exact first-two-pivot localization for every positive cutoff,
  rank-two secant/interpolation analysis, and the proved two-step half
  contraction on `0<Λ≤1`.

## Deferred research

No v1 implementation task is blocked. Open research and optional extensions
are:

- transport the proved `Λ≤1` two-corner half contraction through dyadic
  frequency layers, together with selected-pivot sign coherence or another
  nonexpansiveness invariant, to obtain the `O(log Λ)` block length required
  by Conjecture G1;
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
- Phase G structure:
  `gecp_error_le_product_of_crossRatioControl`,
  `gecp_error_le_dyadic_of_crossRatioProduct`,
  `fermionicKernel_no_uniform_firstStep_contraction`,
  `fermionicKernel_firstPivot_abs_le_reflectedCorner`, and
  `fermionicKernel_twoCornerResidual_le_half_initial`;
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

The research API additionally exports exact-form decimal evaluators
`fermionic_first_residual` and `fermionic_two_corner_residual`, census block
analysis, and focused interval-certified trajectories. The two-corner
evaluator independently checks the Lean base theorem; it is not used as the
source of that theorem.

It also exports `SyntheticApplicationConfig`,
`run_synthetic_application_suite`, typed result records,
`synthetic_spectral_density`, `fermionic_green_from_density`, quadrature
weights, and deterministic sensor-point generation. The application runner
combines existing public algorithms; it does not bypass their convergence,
conditioning, or stop metadata.

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
uv run python experiments/synthetic_applications.py \
  --config experiments/synthetic_applications.toml \
  --output experiments/data/synthetic_applications.json \
  --plot experiments/data/synthetic_applications.png
git status --short
```

## Verification results

- Lean build and public axiom audit pass; no `sorry`, `admit`, or custom axiom
  is present.
- Ruff, formatting, strict mypy, and all 45 pytest cases pass.
- `pip-audit` reports no known vulnerabilities; the local editable project is
  correctly skipped because it is not a PyPI dependency.
- The canonical 21-case census is byte-reproducible with SHA-256
  `a0a0a58271000ddd1efcf8514d3ae404eac359a900c7ba4ebc5c92336bf38179`.
- The exact-surrogate dataset has SHA-256
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`.
- The 128-bit release smoke run converges at tolerance `1e-6` with rank 6 for
  cutoff 1 and rank 82 for cutoff `10⁶`. Its two-record JSONL has SHA-256
  `d7b5f53d61eb2ea49cd1ebde944b47a06f931b9d82687c8fa4458ac4c09907d5`.
- The synthetic-application JSON is byte-identical across repeated runs and
  has SHA-256
  `f1e989153aa18afa915e3f75630d64019acd7939f65802dd5763b7e6892a8ac2`;
  the reviewed plot has SHA-256
  `0529263affd2fe10773f71784f06007650bff833374fbbebbb16f2507b2ab169`.
- At cutoff 8 and tolerance `1e-8`, a 12-pivot universal GECP basis gives
  Green errors `4.16e-10` and `2.91e-10` for the Hubbard-like and gapped
  densities, below the common `5.40e-9` discrete transfer bound.
- The known-transition fixture is recovered with four atoms to `3.33e-16`
  held-out error. The noisy blind scan uses eight effective atoms and reaches
  `8.56e-9` held-out function error without claiming atom identification.
- GECP and pivoted Cholesky agree on 31 of 72 covariance landmarks; the
  selected cross has maximum error `7.78e-7` and relative Frobenius error
  `1.14e-7`.

## Known limitations and non-claims

- No solution of Simons Problem 4.2 is claimed.
- The low-rank theorem is not a GECP theorem.
- Census curves and pivot certificates are evidence, not a continuous GECP
  convergence proof.
- The proved `0<Λ≤1` two-step half contraction is not extrapolated to larger
  cutoffs; the dyadic block scaling required by Conjecture G1 remains open.
- Arbitrary PSD tie-breaking is not claimed equivalent to diagonal pivoting;
  the formal theorem uses the documented diagonal-preferring canonical rule.
- The formal separated construction is dyadic truncated Taylor, not the
  selected-exponential Chebyshev variant.
- The post-v1 spectral fixtures are synthetic analogues, not material-
  specific calculations, experimental fits, or uncertainty-quantified
  analytic continuation.
- Accurate held-out Green-function recovery in the noisy dense scan does not
  identify the four generating atoms uniquely; the eight recovered atoms are
  an effective representation.
- Matching covariance pivots do not prove optimal sensor placement.

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
- Synthetic application configuration, canonical JSON, and reviewed summary
  plot under `experiments/` and `experiments/data/`.
- Release artifacts: `dist/kernelgecp-1.0.0-py3-none-any.whl` and
  `dist/kernelgecp-1.0.0.tar.gz`.
- Claim and research documentation: `README.md`, `FINDINGS.md`, `RESEARCH.md`,
  `APPLICATION.md`, `MATHLIB.md`, and `SPEC.md`.
- Preserved workflow input: `autonomous-implementation.md`, explicitly labeled
  as a reference rather than current status and excluded from package builds.

## How to run locally

```bash
uv sync --all-extras
lake exe cache get
uv run python -c "import kernelgecp; print(kernelgecp.__version__)"
uv run python experiments/synthetic_applications.py
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

Prove a dyadic renormalization/localization lemma that carries the exact
two-corner `Λ≤1` base contraction through successive frequency layers while
controlling intervening Schur updates. The strongest current route is to pair
selected-pivot `PivotCrossProductSignCoherent` nonexpansiveness with one strict
half reduction per `O(s+1)` pivots for `Λ≤2ˢ`; failure of either property should
be minimized and certified before the condition is weakened.
