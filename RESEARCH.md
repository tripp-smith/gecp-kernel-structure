# Research contracts

## Phase 0 — bootstrap

State: complete; merged in [PR #1](https://github.com/tripp-smith/gecp-kernel-structure/pull/1).

Deliverables are the repository-scoped skills, pinned Lean and Python
toolchains, package roots, verification scripts, CI, and progress metadata.
This phase makes no mathematical claim.

## Milestone A — exact selected core closure

State: complete; merged in [PR #4](https://github.com/tripp-smith/gecp-kernel-structure/pull/4).

The phase closes the remaining gap between the dependent `GECP.Run` object and
the determinant/nonsingularity API. The exact public contract is:

- `Run.SelectedIndex` enumerates the successful pivots in their run order;
- `Run.selectedCore` evaluates the original kernel at those selected rows and
  columns;
- `Run.finSelectedCore` reindexes that same matrix by `Fin run.pivots.length`;
- `gecp_core_det_eq_prod_pivots` proves that the determinant of the finite
  selected core is the ordered product of the actual residual pivots carried
  by `run`;
- `gecp_core_nonsingular` proves the selected core determinant is nonzero.

The proof must derive its block elimination step from `residualUpdate` and the
nonzero pivot witness already stored in `Run.step`. It may not assume external
LDU data. This phase does not add a GECP convergence or pivot-selection claim.

## Milestone B — recursive PSD/Cholesky equivalence

State: complete; merged in [PR #5](https://github.com/tripp-smith/gecp-kernel-structure/pull/5).

The phase strengthens the existing one-step diagonal-domination result to the
following exact contract for finite real PSD matrices:

- `diagonalResidual` is the symmetric same-domain Schur update at a positive
  diagonal pivot;
- `posDef_gecp_residual_posSemidefinite` proves that update remains PSD;
- `IsDiagonalPivot` and `IsCompletePivot` state the pivoted-Cholesky and GECP
  maximizer predicates;
- `PivotedCholeskyTrace` and `GECPTrace` record the same ordered residual
  updates with their respective pivot predicates;
- `gecp_eq_pivotedCholesky_of_posDef` proves the two recursive trace relations
  equivalent for every pivot list.

The complete-pivot relation deliberately uses a diagonal pivot location and
the repository's diagonal-preferring canonical tie policy. The theorem does
not claim that every arbitrary off-diagonal tie choice produces the same run,
and it does not claim the stronger published Lipschitz convergence rate.

## Milestone C — high-precision certified census closure

State: in progress.

The kernel identities and stable float64 evaluator are already proved and
tested. This closure phase replaces the remaining census placeholders with:

- a genuine `mpmath` finite-grid GECP backend that executes at
  `GECPConfig.precision_bits` and preserves canonical quantities as decimal
  strings;
- an outward-rounded interval branch-and-bound pivot certificate for the
  fermionic residual on the continuous rectangle;
- adaptive GECP results that retain every pivot certificate and never report
  certification after a cell-budget exhaustion;
- a canonical census whose recorded precision and algorithm match the backend
  actually used, with repeat-run byte equality;
- the kernel derivative bounds and grid-to-continuous theorem/API constants
  needed to interpret approximate pivots.

The numerical certificates validate computations; they are not substituted
for Lean proofs of the Milestone D or E rate theorems.

## Milestone D — constructive separated approximation

State: blocked (research).

Delivered:

- exact natural-number scale/order count definitions;
- independently tested piecewise composite-Chebyshev interpolation in Python;
- boundary and dense-grid regression checks.

Missing theorem contract:

- `expFamily_separatedApprox` must construct the published selected
  exponentials and prove the explicit uniform remainder bound;
- `fermionicKernel_separatedApprox` must combine positive, reflected negative,
  and central bands with an explicit natural-number rank and error budget.

Those theorem identifiers are deliberately absent. A numerical interpolant is
not a substitute for the analytic construction.

## Milestone E — structural GECP

State: blocked (research).

Exact checks completed for all requested geometric surrogates. There were no
zero minors and no sign mismatches. However, the exact pivot sequences depend
on `q`; for size eight the three requested ratios produce different row and
column orders. This is a certified obstruction to the proposed invariant
“GECP has one parameter-independent ordered pivot path.” It is not yet the
stronger SPEC-approved obstruction/refined sufficient condition needed to
close Milestone E.

The continuous finite-grid census also completed, but its float64 plateau at
large cutoff is not suitable theorem evidence. The next theorem contract is:

1. formalize sign regularity for geometric/exponential collocation minors;
2. state a quantitative cross-ratio hypothesis on Schur residuals;
3. prove that hypothesis implies an explicit near-max-volume comparison;
4. test whether the fermionic residual satisfies the hypothesis using
   outward-rounded or exact bounds.

Milestone E may close only with one rigorous outcome allowed by `SPEC.md`: the
target rate, a weaker improved rate, a structural-class theorem containing the
fermionic kernel, or a certified obstruction with a refined sufficient
condition. Experiments alone cannot close the milestone.

## Numerical evidence identifiers

- `experiments/data/exact_surrogates.json`:
  `ccd5d4b948629a6e9642bd5fc18c69c228709752c31c80113036c0e1037f0e82`
- `experiments/data/gecp_census.jsonl`:
  `3381573b39addf3d5fb64c58f35f3923735b859a1d41953edbba947f88959676`
