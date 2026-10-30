# Proposed repository

**Repository name:** `gecp-kernel-structure`

**GitHub description:**
Lean 4 formalization and computational study of structure-aware GECP bounds for continuous kernels, with pivoted Cholesky baselines and the fermionic DLR kernel from Simons Problem 4.2.

I would make the **fermionic kernel the principal research target**, not continuous positive-definite kernels in general. The workshop already describes the positive-definite case as having a relatively satisfactory answer through pivoted Cholesky/P-greedy theory. The interesting unresolved gap is the specific nonsymmetric DLR kernel

\[
K(t,\omega)=\frac{e^{-t\omega}}{1+e^{-\omega}},
\qquad
t\in[0,1],\quad \omega\in[-\Lambda,\Lambda],
\]

where current generic GECP theory gives a much weaker dependence on \(\Lambda\) than both observed behavior and known low-rank representations.

# SPEC.md

## 1. Project objective

Simons Problem 4.2 asks for stronger theoretical bounds for Gaussian elimination with complete pivoting, GECP, when it is applied to a structured continuous kernel.

The motivating kernel is the fermionic Lehmann kernel

\[
K_\Lambda(t,\omega)
=
\frac{e^{-t\omega}}{1+e^{-\omega}},
\qquad
(t,\omega)\in[0,1]\times[-\Lambda,\Lambda].
\]

Given selected nodes

\[
T_k=\{t_1,\ldots,t_k\},
\qquad
\Omega_k=\{\omega_1,\ldots,\omega_k\},
\]

GECP constructs the cross approximation

\[
\widehat K_k(t,\omega)
=
K(t,\Omega_k)
K(T_k,\Omega_k)^{-1}
K(T_k,\omega).
\]

This gives, for

\[
G(t)=\int_{-\Lambda}^{\Lambda}
K(t,\omega)\rho(\omega)\,d\omega,
\]

the induced approximation

\[
\widehat G(t)
=
\int_{-\Lambda}^{\Lambda}
\widehat K_k(t,\omega)\rho(\omega)\,d\omega
\]

and the immediate error transfer

\[
\|G-\widehat G\|_\infty
\le
\|K-\widehat K_k\|_\infty
\|\rho\|_1.
\]

The workshop records a generic analytic-kernel GECP estimate of the form

\[
k=O\!\left(\Lambda+\log(1/\varepsilon)\right)
\Longrightarrow
\|K-\widehat K_k\|_\infty\le\varepsilon,
\]

yet a rank-\(k\) approximation exists at the much better scale

\[
k=
O\!\left(
\log\Lambda\,
\log(1/\varepsilon)
\right).
\]

Empirically, GECP behaves much closer to this latter scale.

The central objective of this repository is:

\[
\boxed{
\text{Explain and prove why GECP exploits the special structure of }K_\Lambda.
}
\]

The ideal result is

\[
\boxed{
k
=
O\!\left(
\log(1+\Lambda)
\log(1/\varepsilon)
\right)
\quad\Longrightarrow\quad
\|K_\Lambda-\widehat K_k\|_\infty
\le
C\varepsilon
}
\]

with \(C\) constant or at most polylogarithmic in \(\Lambda\) and \(1/\varepsilon\).

A weaker but still meaningful result is any theorem that replaces the existing linear dependence on \(\Lambda\) by polynomial or polylogarithmic dependence.

## 2. Why this problem is a good successor project

The Nyström project succeeded by separating three layers:

1. an exact algebraic reduction;
2. a structural theorem;
3. a small executable implementation checked against formally certified examples.

Problem 4.2 admits a similar structure.

For positive-definite kernels, GECP reduces to pivoted Cholesky, also called P-greedy in the kernel literature. Existing theory already relates convergence to smoothness and geometric fill distance. Recent work proves an \(O(k^{-1/d})\) uniform residual rate for Lipschitz positive-definite kernels on compact subsets of \(\mathbb R^d\), improving to \(O(k^{-2/d})\) with one additional Lipschitz derivative.  Earlier P-greedy theory gives near-optimal convergence for Sobolev kernels and asymptotically uniform point distributions.

That positive-definite setting should be the **baseline/control theorem**, not the final research target.

The real opportunity is to isolate comparable structure in the fermionic kernel that survives despite its being rectangular and nonsymmetric.

## 3. Two-track research program

The project should have two independent tracks.

### Track A: continuous positive-definite baseline

Formalize enough of the known pivoted-Cholesky/P-greedy theory to establish the structural connection

\[
\text{GECP}
\quad\Longrightarrow\quad
\text{diagonal pivoting}
\quad\Longrightarrow\quad
\text{P-greedy}
\]

for symmetric positive-definite kernels.

This provides:

- a mathematically understood case;
- reusable residual and pivot machinery;
- a control experiment for the general GECP implementation;
- the correct conceptual vocabulary for kernel geometry.

### Track B: fermionic DLR kernel

Study

\[
K_\Lambda(t,\omega)
=
\frac{e^{-t\omega}}{1+e^{-\omega}}
\]

directly.

The goal is either:

1. prove a near-optimal GECP convergence theorem;
2. identify a stronger structural hypothesis containing the DLR kernel for which such a theorem holds;
3. or rigorously identify the obstruction preventing the desired bound.

Track B is the actual Problem 4.2 contribution.

## 4. Structural identities of the fermionic kernel

The first phase should formalize the elementary identities that generic analytic-kernel bounds ignore.

### 4.1 Centered representation

Write

\[
x=2t-1\in[-1,1],
\qquad
y=\omega/2.
\]

Then

\[
K(t,\omega)
=
\frac{e^{(1/2-t)\omega}}
{2\cosh(\omega/2)}
=
\frac{e^{-xy}}{2\cosh y}.
\]

Thus

\[
\boxed{
K_\Lambda(t,\omega)
=
w(\omega)\,e^{-x\omega/2},
\qquad
w(\omega)=\frac1{2\cosh(\omega/2)}.
}
\]

This decomposes the kernel into:

- a positive one-dimensional column weight;
- an exponential interaction kernel.

Column scaling does not alter the essential location of low-rank structure and preserves the sign pattern of minors.

### 4.2 Reflection symmetry

The kernel satisfies exactly

\[
\boxed{
K(t,\omega)
=
K(1-t,-\omega).
}
\]

This implies that a symmetric pivot configuration has a reflected counterpart and should be explicitly exploited in both analysis and computation.

### 4.3 Exponential-kernel structure

After centering and positive scaling, the core interaction is

\[
e^{-xy}.
\]

Exponential kernels have strong sign-regularity / total-positivity structure under ordered arguments. The project should investigate whether the relevant discretizations of the fermionic kernel inherit enough total positivity to control complete pivoting.

This is a **research hypothesis**, not an initial theorem claim.

The key possible route is:

\[
\text{ordered exponential structure}
\Longrightarrow
\text{controlled minors}
\Longrightarrow
\text{controlled GECP pivots}
\Longrightarrow
\text{near-max-volume cross}
\Longrightarrow
\text{near-optimal residual}.
\]

## 5. DLR approximation baseline

The Discrete Lehmann Representation already proves that imaginary-time Green's functions admit an exponential representation whose basis size scales as

\[
r=
O\!\left(
\log(\beta\omega_{\max})
\log(1/\varepsilon)
\right).
\]

The basis consists of selected exponentials, and corresponding interpolation nodes can be constructed numerically.

This repository should first isolate the kernel-level approximation theorem underneath that observation.

Target theorem:

> **Theorem D1, explicit separated approximation.**
> For every \(\Lambda\ge1\) and \(0<\varepsilon<1\), construct
>
> \[
> K_r(t,\omega)
> =
> \sum_{j=1}^r
> f_j(t)g_j(\omega)
> \]
>
> satisfying
>
> \[
> \|K_\Lambda-K_r\|_\infty\le\varepsilon
> \]
>
> with
>
> \[
> r
> \le
> C
> \log(1+\Lambda)
> \log(C/\varepsilon).
> \]

This theorem does **not** yet prove anything about GECP.

Its purpose is to establish a formally verified target against which the greedy algorithm can be compared.

## 6. Candidate proof of the low-rank bound

A practical formal route is a dyadic decomposition of the frequency domain.

Partition

\[
[1,\Lambda]
\]

into geometrically increasing intervals

\[
[1,2], [2,4], [4,8],\ldots
\]

with \(O(\log\Lambda)\) blocks, together with the corresponding negative-frequency intervals and a bounded neighborhood of zero.

On each scale, apply a fixed-degree polynomial or Chebyshev approximation to the rescaled exponential interaction.

If each block needs

\[
O(\log(1/\varepsilon))
\]

terms, then summing over

\[
O(\log\Lambda)
\]

frequency scales gives

\[
O(\log\Lambda\log(1/\varepsilon))
\]

separated terms.

This route is attractive for Lean since it breaks the global result into:

- interval rescaling;
- elementary exponential approximation;
- uniform error estimates;
- geometric-series bookkeeping.

The initial formal proof need not reproduce the most optimized DLR construction.

## 7. Continuous GECP definition

Define a residual recursively.

Set

\[
R_0(t,\omega)=K(t,\omega).
\]

At iteration \(j\), choose

\[
(t_j,\omega_j)
\in
\arg\max_{(t,\omega)}
|R_{j-1}(t,\omega)|.
\]

Let

\[
p_j=R_{j-1}(t_j,\omega_j).
\]

Then update

\[
R_j(t,\omega)
=
R_{j-1}(t,\omega)
-
\frac{
R_{j-1}(t,\omega_j)
R_{j-1}(t_j,\omega)
}{p_j}.
\]

The cross approximation is

\[
\widehat K_j=K-R_j.
\]

For compact domains and continuous residuals, existence of a maximizing pivot follows from compactness.

The Lean definition should distinguish:

- exact continuous GECP;
- GECP on a finite candidate grid;
- approximate pivoting with a multiplicative or additive pivot tolerance.

The implementation will necessarily use the latter two.

## 8. GECP algebraic invariants

Before attacking rates, formalize the exact identities.

After \(k\) successful pivots:

\[
R_k(t_i,\omega)=0,
\qquad
R_k(t,\omega_i)=0
\]

for every selected row or column node.

The approximation interpolates the kernel on the cross.

The selected core matrix

\[
K(T_k,\Omega_k)
\]

is nonsingular whenever all pivots are nonzero.

Its determinant satisfies

\[
\boxed{
\det K(T_k,\Omega_k)
=
\prod_{j=1}^k p_j
}
\]

up to the chosen ordering/sign convention.

This is the continuous analogue of the standard Gaussian-elimination determinant identity and connects pivot growth directly to volume.

## 9. Positive-definite baseline theorem

For a continuous symmetric positive-definite kernel

\[
K:\Omega\times\Omega\to\mathbb R,
\]

the row and column pivot coincide.

After selecting \(x_j\),

\[
R_j(x,y)
\]

remains positive semidefinite, and

\[
|R_j(x,y)|^2
\le
R_j(x,x)R_j(y,y).
\]

Therefore

\[
\boxed{
\|R_j\|_\infty
=
\max_x R_j(x,x).
}
\]

Hence complete pivoting chooses

\[
x_{j+1}
\in
\arg\max_x R_j(x,x),
\]

which is exactly pivoted Cholesky/P-greedy.

This equivalence should be completely formalized.

It is a crisp theorem with the same role as the Schur-complement identity in the Nyström project.

## 10. Baseline Lipschitz result

The full modern Lipschitz convergence theorem does not have to be the first Lean milestone.

A useful staged target is:

> If \(K\) is Lipschitz and SPD, bound the residual at \(x\) in terms of the distance from \(x\) to its closest selected pivot.

Then derive

\[
\|R_k\|_\infty
\le
C_K h(X_k),
\]

where

\[
h(X_k)
=
\sup_{x\in\Omega}
\min_{x_j\in X_k}
\|x-x_j\|
\]

is the fill distance.

Jeong and Townsend prove this type of bound and obtain

\[
\|R_k\|_\infty=O(k^{-1/d})
\]

for complete pivoting on a compact \(d\)-dimensional domain, with

\[
O(k^{-2/d})
\]

under stronger differentiability assumptions.

For this repository, formalizing the exact equivalence plus a simplified fill-distance theorem is sufficient for the baseline phase.

## 11. Main GECP conjecture

Define

\[
e_k(\Lambda)
=
\|K_\Lambda-\widehat K_k^{\mathrm{GECP}}\|_\infty.
\]

The central conjecture should be recorded as:

> **Conjecture G1.** There exist universal constants \(C,c>0\) such that
>
> \[
> e_k(\Lambda)
> \le
> C
> \exp\left(
> -c\frac{k}{\log(1+\Lambda)}
> \right).
> \]

This is equivalent, up to constants, to

\[
k
=
O\!\left(
\log(1+\Lambda)
\log(1/\varepsilon)
\right)
\]

for achieving error \(\varepsilon\).

Do not state this as a theorem until proved.

## 12. Approximation-to-greedy reduction

A major theoretical route should be to prove an abstract theorem of the following kind.

Suppose a kernel \(K\) has best rank-\(r\) error

\[
\sigma_r(K)
=
\inf_{\operatorname{rank}F\le r}
\|K-F\|_\infty.
\]

Can GECP satisfy

\[
\boxed{
\|R_k\|_\infty
\le
C(k)
\sigma_{\alpha k}(K)
}
\]

with \(C(k)\) polynomial rather than exponential for the structural class containing \(K_\Lambda\)?

For completely general matrices/kernels, this is false at the desired strength.

Thus the goal is to identify which property of the fermionic kernel kills the exponential growth factor.

Candidate properties:

- total positivity;
- sign regularity;
- monotonicity of derivatives;
- variation-diminishing behavior;
- reflection symmetry;
- log-concavity of column weights;
- ordered pivot geometry;
- bounded cross ratios;
- scale-local analyticity.

This abstract structural theorem may be the most reusable mathematical result from the project.

## 13. Total positivity research track

The centered kernel

\[
\widetilde K(x,y)=e^{-xy}
\]

is the first object to investigate.

For ordered point sets

\[
x_1<\cdots<x_m,
\qquad
y_1<\cdots<y_m,
\]

determine the exact signs of

\[
\det[e^{-x_i y_j}]_{i,j=1}^m.
\]

Positive row/column scaling then transfers corresponding sign-regularity statements to \(K_\Lambda\).

Research questions:

1. Is every relevant minor nonzero?
2. Is its sign determined entirely by \(m\)?
3. Does Schur complementation preserve the same kernel sign structure?
4. Do GECP residuals inherit monotone or sign-regular sections?
5. Can the maximum absolute residual be localized to predictable boundary or scale-transition regions?
6. Does complete pivoting become equivalent to a simpler nested node-selection rule?

A strong positive answer here could transform Problem 4.2 from arbitrary two-dimensional pivot search into a one-dimensional geometric problem.

## 14. Pivot geometry experiments

Implement high-precision GECP for the exact kernel over increasingly dense adaptive grids.

For

\[
\Lambda\in
\{1,10,10^2,10^3,10^4,10^5,10^6\}
\]

and tolerances down to at least \(10^{-12}\), record:

- pivot \(t_j\);
- pivot \(\omega_j\);
- pivot magnitude \(p_j\);
- residual supremum;
- core determinant;
- smallest singular value of the core;
- distance to boundaries;
- dyadic scale of \(|\omega_j|\);
- reflection partner;
- empirical rank required for each tolerance.

The main plots should test whether

\[
k/\log(1+\Lambda)
\]

collapses the residual curves onto a common exponential decay law.

## 15. Exact small-instance layer

As with the Nyström repository, floating-point evidence should be backed by exact small examples where possible.

Use rational or algebraic surrogate kernels such as:

\[
K_q(i,j)=q^{ij},
\]

or finite exponential matrices

\[
E_{ij}=x_i^{y_j}
\]

with rational \(x_i\).

These retain the exponential/total-positive structure but permit exact determinant and pivot comparisons.

Lean should certify:

- complete pivot sequence for small grids;
- determinant-product identity;
- sign of every relevant minor;
- monotonicity patterns suggested by the numerical experiments.

## 16. Discrete-to-continuous bridge

Practical GECP runs on a finite grid.

Let

\[
\mathcal T_h\subset[0,1],
\qquad
\Omega_h\subset[-\Lambda,\Lambda].
\]

Suppose each residual \(R_k\) has a known Lipschitz bound

\[
|R_k(z)-R_k(z')|
\le
L_k\|z-z'\|.
\]

Then

\[
\sup_D |R_k|
\le
\max_{z\in D_h}|R_k(z)|
+
L_k h.
\]

This provides a rigorous route from finite-grid pivot searches to the continuous GECP algorithm.

One target is an adaptive certified maximizer:

```text
evaluate residual
bound local Lipschitz constant
subdivide cells whose upper bounds exceed current maximum
terminate when pivot is certified within tolerance
```

This would make the computational experiments much stronger than ordinary dense-grid sampling.

## 17. Approximate GECP theorem

Numerical implementations cannot select the exact global maximizer.

Define \(\eta\)-complete pivoting by

\[
|R_k(t_{k+1},\omega_{k+1})|
\ge
\eta\|R_k\|_\infty,
\qquad
0<\eta\le1.
\]

Every structural convergence theorem should, where feasible, be strengthened to this approximate-pivot setting.

The target bound would have constants depending explicitly on \(\eta\), rather than assuming mathematically exact maximization.

This is essential if the theorem is intended to explain the practical DLR construction.

## 18. Green's-function consequence

Once a kernel error theorem is established, expose the downstream result as a separate theorem:

> If
>
> \[
> \|K-\widehat K_k\|_\infty\le\varepsilon,
> \]
>
> then for every \(\rho\in L^1[-\Lambda,\Lambda]\),
>
> \[
> \left\|
> K\rho-\widehat K_k\rho
> \right\|_\infty
> \le
> \varepsilon\|\rho\|_1.
> \]

This simple inequality is explicitly part of the workshop motivation.

It should be a public theorem because it connects the kernel result directly to its physical application.

## 19. Second research direction: sparse \(\rho\)

The workshop proposes a distinct, problem-dependent formulation.

Given a particular

\[
G=K\rho_0,
\]

find another representation

\[
G=K\rho
\]

with

\[
\|\rho\|_0
\]

as small as possible. The desired support might be substantially smaller than a universal cross approximation designed to approximate every function in the range of \(K\).

This should be treated as **Phase II**, not mixed into the GECP theorem.

The finite sparse model is

\[
G(t)
\approx
\sum_{j=1}^{s}
g_jK(t,\omega_j).
\]

The key research objective becomes

\[
\min s
\quad\text{subject to}\quad
\left\|
G-\sum_{j=1}^s
g_jK(\cdot,\omega_j)
\right\|_\infty
\le\varepsilon.
\]

## 20. Rational / sum-of-exponentials bridge

The workshop notes that in one dimension this perspective can sometimes be transformed into rational-approximation or sum-of-exponentials problems and mentions AAA-type methods as useful tools, yet says comparable methods are not currently practical for multidimensional integration kernels.

The project should not initially attempt a formalization of AAA.

Instead, define an interface:

\[
\text{sparse exponential representation}
\Longrightarrow
\text{sparse spectral measure}.
\]

Then experimentally compare:

- GECP/DLR universal nodes;
- AAA-derived representations;
- vector fitting;
- Prony-type recovery;
- nonlinear least squares;
- greedy matching pursuit over exponential atoms.

The central metric is

\[
s_{\mathrm{specific}}(\varepsilon)
\quad\text{vs.}\quad
k_{\mathrm{universal}}(\varepsilon).
\]

## 21. Lean architecture

Suggested module layout:

```text
GECPKernelStructure/
  Definitions.lean
  CrossApproximation.lean
  GECP/
    Definitions.lean
    Residual.lean
    Interpolation.lean
    Determinant.lean
    ApproxPivot.lean
  PositiveDefinite/
    ResidualPSD.lean
    PivotedCholesky.lean
    PowerFunction.lean
    FillDistance.lean
  Fermionic/
    Kernel.lean
    Symmetry.lean
    Centering.lean
    Derivatives.lean
    DyadicPartition.lean
    SeparatedApprox.lean
    TotalPositivity.lean
    GECPBounds.lean
  Discretization/
    Grid.lean
    LipschitzCertificate.lean
  GreenFunction.lean
  Computable.lean
  SmallInstanceChecks.lean
  Counterexamples/
  Theorems.lean
  MathlibReady.lean
```

The high-risk theorem files should depend on a stable low-level core, not the reverse.

## 22. Headline theorem names

Target names should be explicit and conservative.

```lean
gecp_interpolates_selected_rows
gecp_interpolates_selected_cols
gecp_core_det_eq_prod_pivots

gecp_eq_pivotedCholesky_of_posDef
posDef_gecp_residual_sup_eq_diag_sup

fermionicKernel_reflection
fermionicKernel_centered
fermionicKernel_continuous
fermionicKernel_separatedApprox

greenError_le_kernelError_mul_l1
```

If the main conjecture is solved:

```lean
fermionicKernel_gecp_error_le_exp
```

with a statement equivalent to

\[
\|R_k\|_\infty
\le
C e^{-ck/\log(1+\Lambda)}.
\]

No theorem name should contain `optimal` until a matching lower bound or explicit approximation-number comparison exists.

## 23. Python application layer

Package name:

```text
kernelgecp
```

Proposed API:

```python
from kernelgecp import (
    FermionicKernel,
    gecp,
    pivoted_cholesky,
    cross_approximation,
    evaluate_residual,
    certified_pivot,
    dlr_rank_bound,
)
```

Example:

```python
K = FermionicKernel(cutoff=1e5)

result = gecp(
    K,
    tol=1e-10,
    pivot="adaptive",
)

result.t_nodes
result.omega_nodes
result.pivots
result.residual
```

For SPD kernels:

```python
result = pivoted_cholesky(
    kernel,
    domain,
    tol=1e-10,
)
```

## 24. Numerical precision

This project will be more sensitive to floating-point effects than the Nyström project.

Support at least:

- `float64`;
- arbitrary precision through `mpmath`;
- optional exact arithmetic for synthetic rational kernels.

Every experiment used to support a mathematical conjecture should be rerunnable at higher precision.

GECP pivot ties or near-ties should be recorded rather than silently resolved.

## 25. Verification strategy

CI should include:

```bash
lake build
python -m pytest
```

and reject `sorry` in the Lean target.

Python tests should cover:

- reflection symmetry;
- centered formula;
- GECP interpolation identities;
- determinant equals pivot product;
- direct cross formula vs iterative residual;
- GECP = pivoted Cholesky for SPD matrices;
- PSD residual preservation;
- finite-grid max residual equals max diagonal for SPD examples;
- high-precision agreement with float64 for well-conditioned examples;
- DLR empirical rank scaling;
- exact exponential-kernel small cases.

## 26. Research census

Maintain a machine-readable experiment dataset.

Suggested schema:

```json
{
  "kernel": "fermionic",
  "Lambda": 100000,
  "tolerance": 1e-10,
  "precision_bits": 128,
  "algorithm": "gecp",
  "rank": 42,
  "max_residual": 8.2e-11,
  "t_nodes": [],
  "omega_nodes": [],
  "pivots": []
}
```

The dataset should be deterministic given configuration and code revision.

## 27. Success criteria

### Milestone A: exact GECP core

Formalize continuous/discrete GECP and prove:

\[
R_k(t_i,\omega)=0,
\qquad
R_k(t,\omega_i)=0,
\]

plus

\[
\det K(T_k,\Omega_k)
=
\prod_{j=1}^k p_j.
\]

### Milestone B: positive-definite equivalence

Prove GECP equals pivoted Cholesky/P-greedy for SPD kernels.

This establishes the known baseline identified by the workshop.

### Milestone C: DLR kernel structure

Machine-check:

\[
K(t,\omega)=K(1-t,-\omega)
\]

and

\[
K(t,\omega)
=
\frac{e^{-(2t-1)\omega/2}}
{2\cosh(\omega/2)}.
\]

Develop derivative, monotonicity and minor-sign infrastructure.

### Milestone D: explicit near-optimal separated approximation

Prove

\[
r
=
O(
\log(1+\Lambda)
\log(1/\varepsilon)
)
\]

suffices for a uniform separated approximation.

This reproduces the target scale underlying DLR theory in a Lean-friendly form. The DLR literature establishes this logarithmic-logarithmic basis-size scaling for imaginary-time Green's functions.

### Milestone E: structural GECP theorem

Deliver at least one of:

1. the desired near-optimal GECP rate;
2. a weaker polylogarithmic-\(\Lambda\) GECP rate;
3. a theorem for a structural kernel class containing \(K_\Lambda\);
4. a rigorous counterexample to a natural stronger conjecture plus a refined sufficient condition.

### Milestone F: sparse-\(\rho\) application

Demonstrate on specified Green's functions that problem-specific sparse representations can use fewer atoms than the universal GECP basis, with explicit approximation error.

This may remain computational rather than formally proved.

## 28. Claim discipline

Before Milestone B:

> Lean formalization of GECP and cross approximation for continuous kernels.

After Milestone B:

> Formal connection between GECP, pivoted Cholesky and P-greedy for positive-definite kernels.

After Milestone D:

> Machine-checked near-optimal separated approximation for the fermionic Lehmann kernel.

Only after Milestone E should the project claim:

> New structure-aware GECP bounds for the fermionic kernel from Simons Problem 4.2.

Do not say:

> Solved Problem 4.2

unless the theorem genuinely closes the \(\Lambda\)-dependence gap posed in the workshop.

Known pivoted-Cholesky and P-greedy convergence results must be attributed as prior work.

## 29. Documentation

Use the same repository discipline as the Nyström project:

```text
README.md
SPEC.md
FINDINGS.md
RESEARCH.md
APPLICATION.md
MATHLIB.md
```

`FINDINGS.md` should clearly separate:

- known prior results;
- formally reproduced prior results;
- new lemmas;
- computational observations;
- conjectures;
- proven new results.

`RESEARCH.md` should make the total-positivity track and the sparse-\(\rho\) track explicit so exploratory work does not leak into headline theorem claims.

## 30. Minimum publishable outcome

Even without solving the main GECP conjecture, this project is worthwhile if it produces:

1. a reusable Lean formalization of continuous cross approximation and GECP;
2. the exact equivalence between SPD GECP and pivoted Cholesky;
3. formally verified algebraic structure of the fermionic kernel;
4. a constructive

\[
O(\log\Lambda\log(1/\varepsilon))
\]

separated approximation theorem;
5. a high-precision empirical census of GECP pivot geometry across six orders of magnitude in \(\Lambda\);
6. either a promising structural conjecture or a certified obstruction.

The strongest outcome would connect the last two layers:

\[
\boxed{
\text{total/sign-regular fermionic structure}
\Rightarrow
\text{polynomially controlled GECP}
\Rightarrow
\|R_k\|_\infty
\lesssim
e^{-ck/\log\Lambda}.
}
\]

That would directly explain the empirical phenomenon Problem 4.2 is asking about.

## Recommended first attack

I would start even more narrowly than with pivoted Cholesky.

**Phase 1 should formalize the fermionic kernel identities and a discrete GECP implementation, then run a high-precision pivot census.**

The first research question I would test is:

\[
\boxed{
\text{Does GECP preserve a total/sign-regular structure in its successive Schur residuals?}
}
\]

If the answer is yes, that is potentially the structural mechanism the workshop problem is looking for. The exponential kernel hidden inside

\[
K(t,\omega)
=
\frac{e^{-(2t-1)\omega/2}}
{2\cosh(\omega/2)}
\]

makes this route much more specific and potentially much stronger than treating \(K\) merely as an analytic function.

That would give this project a similar shape to `nystrom-submodularity`: start with exact small structural facts, use computation to identify the invariant, and only then generalize it into the headline theorem.

---

## Implementation status annotation — 2026-08-15

This annotation records delivered identifiers without changing the goals above.

- Bootstrap: implemented and merged in PR #1 with repository skills `phase-cadence` and
  `lean-research-loop`, pinned Lean/mathlib 4.33.0, Python 3.12, CI, and a root
  verification command.
- Milestone A: `gecp_interpolates_selected_rows`,
  `gecp_interpolates_selected_cols`, `gecp_core_det_eq_prod_pivots`, and
  `gecp_core_nonsingular` are implemented for the finite selected core of every
  successful dependent `Run`. The determinant proof derives its block
  elimination step directly from `residualUpdate` and the run's stored nonzero
  pivot witnesses.
- Milestone B: same-domain Schur residual PSD preservation, diagonal
  domination, bound-equivalence, canonical diagonal-maximizer, recursive GECP
  and pivoted-Cholesky trace equivalence, and the conditional fill-distance
  theorem are implemented.
- Milestone C: reflection, centered-form, continuity, exact coordinate
  derivatives and domain bounds, stable numerical kernel, true 128-bit GECP,
  interval-certified adaptive pivoting, exact surrogate protocol, and an
  endpoint-resolved byte-reproducible finite-grid census are implemented.
- Milestone D: complete in merged PR #7. `expFamily_separatedApprox` proves an explicit
  `8p(s+1)`-term dyadic approximation of `exp(-tω)` on the positive-frequency
  domain, and `fermionicKernel_separatedApprox` proves a `16p(s+1)`-term
  approximation on `[0,1] × [-2ˢ,2ˢ]`, both with uniform error `2⁻ᵖ`.
  The construction is dyadic truncated Taylor rather than the published
  selected-exponential Chebyshev variant; it establishes the required
  logarithmic-logarithmic separated-rank scale without making a GECP claim.
- Milestone E: complete in merged PR #8 under the approved obstruction outcome.
  `signRegular_not_sufficient_for_parameterIndependent_pivots` is a minimized
  exact counterexample showing that strict minor signs do not determine a
  universal complete pivot. `CrossRatioControl`,
  `residualUpdate_le_of_crossRatioControl`, and
  `gecp_error_le_geometric_of_crossRatioControl` give the refined sufficient
  condition and its geometric rate. No claim is made that the fermionic
  residual sequence satisfies this condition with a contraction factor below
  one.
- Milestone F: the continuous-measure Green error transfer and tested sparse
  two-atom application are implemented.
- Release readiness: delivered in merged PR #9 and marked complete by PR #10 with the root
  verification command, package artifacts, and a 128-bit endpoint census
  smoke run; all mathematical and application milestones have merged
  deliveries.
- Implementation-run provenance: the exposed model family, unavailable exact
  serving version, Default collaboration mode, completion-time token and
  elapsed-time snapshot, and API-equivalent cost methodology are recorded in
  [FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). This is a
  status annotation only and does not modify the mathematical objectives or
  claim level of this specification.
- Post-v1 application phase S: complete. Public Python APIs
  `SyntheticApplicationConfig` and `run_synthetic_application_suite` connect
  the delivered fermionic GECP, Green-error transfer, sparse recovery, and PSD
  baseline to continuous Hubbard-like and gapped spectra, clean/noisy atomic
  spectra, and covariance landmarks. The canonical JSON and plot are under
  `experiments/data/`. These synthetic observations add application evidence
  only; they do not modify any formal theorem, Milestone E outcome, or
  non-claim above.
- Post-v1 research phase G: complete as a bounded
  obstruction-plus-base-case result in merged PR #11. Exact theorems now
  rule out cutoff-independent one-step contraction, prove variable-factor and
  dyadic block-product consequences, and show nonexpansiveness under
  `PivotCrossProductSignCoherent`. The first two corner pivots are proved
  complete for every positive cutoff, and their actual two-step residual is
  proved to contract by one half on `0<Λ≤1`. The 128-bit census supports, but
  does not prove, the extension to `2(s+1)`-pivot blocks when `Λ ≤ 2ˢ`.
  Conjecture G1 and the Simons Problem 4.2 claim remain open.
- Post-v1 verification phase H: locally verified with green CI on 2026-10-07
  and delivered in merged PR #12. The root and
  CI test entry points now use the synchronized Python interpreter, and the
  development lock resolves `urllib3` 2.8.0 after the prior 2.7.0 lock gained
  audit findings. The complete root command passes with 45 tests and no known
  dependency vulnerabilities in both maintained and fresh environments. This
  maintenance phase changes no theorem, experiment, or Phase G claim.
- Post-v1 research phase I: complete with green CI on 2026-10-07 and delivered
  in merged PR #13. The new
  `Run.borderedCore_det_eq_selectedCore_det_mul_finalResidual` theorem
  identifies each selected-core border determinant with the core determinant
  times the exact final residual, and
  `pivotCrossProductSignCoherent_iff_borderedMinorSignCoherent` reduces the
  residual nonexpansiveness condition to a four-bordered-minor sign condition
  on the original kernel. Exact surrogate diagnostics retain a first negative
  witness and all 21 canonical geometric runs pass. The root verification
  command passes with 47 tests. This does not prove the bordered-minor
  condition for the continuous fermionic kernel, Conjecture G1, or the Simons
  Problem 4.2 rate.
- Post-v1 research phase J: complete with green CI on 2026-10-07 and delivered
  in merged PR #14.
  `StrictSignRegularAtOrder` and
  `StrictSignRegular` formalize the ordered-minor hypothesis, and
  `strictSignRegular_borderedMinorSignCoherent` transfers its all-orders form
  to the Phase I four-minor condition for every successful run, including
  duplicate proposed nodes. Positive column scaling transfers the hypothesis
  from `exp(-t*omega)` to the fermionic kernel. Orders one and two are proved;
  the all-orders exponential theorem, dyadic localization, Conjecture G1, and
  the Simons Problem 4.2 rate remain open. The root verification command passes
  with 48 tests and no known dependency vulnerabilities.
- Post-v1 research phase K: complete with green CI on 2026-10-07 and delivered
  in merged PR #15.
  `strictSignRegularAtOrder_rowScale` completes positive coordinate scaling,
  `vandermonde_det_pos_of_strictMono` proves positivity of ordinary ordered
  Vandermonde determinants, and `expTaylorPrincipalMatrix_det_pos` proves that
  the `0, ..., n-1` exponential Taylor feature block has positive determinant
  at every order. This is one strictly positive Cauchy--Binet term, not the
  all-orders exponential sign theorem. Generalized Vandermonde positivity, the
  exponential-series determinant limit, Conjecture G1, and the Simons Problem
  4.2 rate remain open. The root command passes with 49 tests.
- Post-v1 research phase L: complete with green CI on 2026-10-08 and delivered
  in merged PR #16.
  `signVariations_lt_card_support` and `positiveRoots_lt_card_support` derive a
  sparse-polynomial positive-root bound from mathlib's Descartes rule;
  `generalizedVandermonde_det_ne_zero` converts that bound to power-matrix
  nonsingularity; and `geometric_generalizedVandermonde_det_pos` plus
  `generalizedVandermonde_det_pos` prove strict positivity on positive strictly
  increasing nodes and strictly increasing natural exponents. The zero-node
  boundary case, finite exponential Cauchy--Binet assembly, determinant limit,
  Conjecture G1, and the Simons Problem 4.2 rate remain open.
  The root verification command passes with 50 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase M: complete with green CI on 2026-10-08 and delivered
  in merged PR #17.
  `generalizedVandermonde_det_nonneg` extends Phase L to nonnegative strictly
  increasing nodes by uniform positive shifts and determinant continuity. The
  finite exponential Cauchy--Binet assembly, determinant limit, Conjecture G1,
  and the Simons Problem 4.2 rate remain open. The root command passes with 51
  tests and no known dependency vulnerabilities.
- Post-v1 research phase N: complete with green CI on 2026-10-08 and delivered
  in merged PR #18.
  `Matrix.det_mul_rect_eq_sum_injective` and `Matrix.det_mul_rect` provide a
  reusable rectangular Cauchy--Binet expansion, while `expTaylorMatrix_apply`
  and `expTaylorMatrix_det_pos` prove that every sufficiently long finite
  exponential Taylor feature matrix has positive determinant on nonnegative
  strictly increasing row and column tuples. The exponential-series
  determinant limit, all-orders exponential and fermionic strict sign
  regularity, Conjecture G1, and the Simons Problem 4.2 rate remain open. The
  root command passes with 52 tests and no known dependency vulnerabilities.
- Post-v1 research phase O: complete with green CI on 2026-10-08 and delivered
  in merged PR #19.
  `expTaylorMatrix_det_ge_principal` retains a fixed positive principal lower
  bound, `expTaylorMatrix_tendsto_expMatrix` and its determinant corollary
  formalize Taylor convergence, and `expMatrix_det_pos` proves strict
  positivity of the full `exp(x*y)` determinant for nonnegative strictly
  increasing node tuples. Arbitrary-node scaling, frequency-reversal
  orientation, all-orders exponential and fermionic strict sign regularity,
  Conjecture G1, and the Simons Problem 4.2 rate remain open. The root command
  passes with 53 tests and no known dependency vulnerabilities.
- Post-v1 research phase P: complete with green CI on 2026-10-08 and delivered
  in merged PR #20.
  `expMatrix_det_pos_of_strictMono` removes node-sign restrictions by positive
  exponential scaling, `Fin.sign_revPerm` computes the reversal orientation,
  and `expKernel_strictSignRegular` plus
  `fermionicKernel_strictSignRegular` prove the predicted ordered-minor signs
  at every order. `fermionicKernel_pivotCrossProductSignCoherent` makes the
  selected-residual nonexpansiveness premise unconditional. Dyadic strict
  contraction, Conjecture G1, and the Simons Problem 4.2 rate remain open. The
  root command passes with 54 tests and no known dependency vulnerabilities.
- Post-v1 research phase Q: complete with green CI on 2026-10-08 and delivered
  in merged PR #21.
  `Run.append` and `Run.finalResidual_append` compose dependent successful
  runs, `GECP.CompletePivotOn` and `GECP.Run.CompleteOn` encode
  restricted-domain exact complete pivoting, and
  `strictSignRegular_gecp_error_nonincreasing`
  propagates Phase P sign coherence into a uniform bound for every realized
  residual prefix. The fermionic corollary
  `fermionicKernel_gecp_error_le_cutoffCorner` bounds every exact
  cutoff-rectangle residual by the initial cutoff-corner pivot. This is
  factor-one nonexpansiveness, not dyadic strict contraction; Conjecture G1 and
  Problem 4.2 remain open. The root command passes with 54 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase T: complete with green CI on 2026-10-08 and delivered
  in merged PR #22.
  `strictSignRegular_gecp_pivotMagnitude_antitone` proves antitonicity of
  selected pivot magnitudes in every realized strictly sign-regular exact
  complete-pivot sequence;
  `GECP.Run.abs_finSelectedCore_det_eq_prod_abs_pivots` identifies the
  absolute selected-core determinant with the product of absolute pivots; and
  `strictSignRegular_gecp_error_pow_le_selectedCore_det` bounds the `n`th
  power of every rank-`n` domain residual entry by that determinant. This
  formalizes the GECP geometric-mean reduction, not determinant decay. A
  quantitative determinant estimate from the dyadic separated approximation
  or a near-volume argument, Conjecture G1, and Problem 4.2 remain open.
  The root command passes with 54 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase U: complete with green CI on 2026-10-08 and delivered
  in merged PR #23. `Matrix.det_add_eq_sum_columnChoices` expands a determinant
  over approximation/error column choices; mixed terms below the exterior
  degree threshold vanish through `Matrix.det_mul_rect_eq_zero_of_card_lt`;
  and `Matrix.abs_det_le_two_pow_mul_of_factors_approx` gives the explicit
  entrywise low-rank perturbation bound
  `2^n n! epsilon^(n-r) C^r`. The fermionic specialization
  `fermionicKernel_sample_det_le_two_pow` converts the existing
  `16p(s+1)`-term, `2^-p` separated approximation into determinant decay for
  arbitrary samples on `[0,1] x [-2^s,2^s]`. Combining this estimate with
  Phase T, optimizing `p`, and absorbing its explicit constants into the
  cutoff-uniform form of Conjecture G1 remain open. The root command passes
  with 55 tests and no known dependency vulnerabilities.
- Post-v1 research phase V: complete with green CI on 2026-10-08 and delivered
  in merged PR #24. Recursive selected-coordinate membership connects
  `Run.CompleteOn` to physical-domain samples, and
  `fermionicKernel_gecp_error_pow_le_two_pow` composes Phase T's residual-
  power inequality with Phase U's determinant decay. At
  `n = 32(2m+1)(s+1)`,
  `fermionicKernel_gecp_error_le_oddBlock` proves the explicit bound
  `|R_n(x,y)| <= 2n 2^-m`. The elementary choice `m = 16(s+1)` then gives
  `|R_n(x,y)| <= 1/2` by
  `n = 32(32(s+1)+1)(s+1)`. This is the first general strict contraction rate
  derived from the formal separated approximation, but its block is quadratic
  in `s+1`; the conjectured `2(s+1)` block, sharp constants, Conjecture G1,
  and Problem 4.2 remain open. The root command passes with 56 tests and no
  known dependency vulnerabilities.
- Post-v1 research phase W: complete with green CI on 2026-10-09 and delivered
  in merged PR #25. `oddBlock_scale_logarithmic` proves the exact scale
  condition for `m = 16 + 2 * Nat.clog 2 (s+1)`, and
  `fermionicKernel_gecp_error_le_half_logarithmicBlock` yields
  `|R_n(x,y)| <= 1/2` by
  `n = 32(33 + 4 * Nat.clog 2 (s+1))(s+1)`. This improves the certified
  general contraction block from quadratic to
  `O((s+1) log(s+1))` without strengthening the complete-pivot hypotheses.
  The conjectured `2(s+1)` block, removal of logarithmic overhead, optimal
  constants, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 57 tests and no known dependency vulnerabilities.
- Post-v1 research phase X: complete with green CI on 2026-10-09 and delivered
  in merged PR #26. `oddBlock_factor_le_accuracy` and
  `oddBlock_scale_accuracy` prove the exact arithmetic needed for any dyadic
  target order `q`, while
  `fermionicKernel_gecp_error_le_dyadicAccuracyBlock` yields
  `|R_n(x,y)| <= 2^-q` by
  `n = 32(33 + 4 * Nat.clog 2 (s+1) + 4q)(s+1)`. Thus the certified rank is
  `O((s+1)(log(s+1)+q))`, equivalently
  `O(log(1+Lambda)(log log(1+Lambda)+log(1/epsilon)))` after choosing a dyadic
  target. The requested-accuracy dependence now matches Conjecture G1; the
  additive logarithmic scale startup, the conjectured `2(s+1)` block, optimal
  constants, Conjecture G1 itself, and Problem 4.2 remain open. The root
  command passes with 58 tests and no known dependency vulnerabilities.
- Post-v1 research phase Y: complete with green CI on 2026-10-09 and delivered
  in merged PR #27. `sylvesterHadamard_det_square` certifies the exact
  determinant growth of the power-of-two Sylvester--Hadamard family, and
  `entrywise_determinant_constant_base_obstruction` proves that no universal
  bound `|det A| <= C^n` follows from unit entry bounds alone. This closes the
  generic determinant-prefactor shortcut: removing Phase X's additive
  logarithmic startup must exploit fermionic total positivity,
  divided-difference structure, or GECP-selected residual geometry. The
  obstruction is not a counterexample under those additional hypotheses and
  does not prove Conjecture G1 or solve Problem 4.2. The root command passes
  with 59 tests and no known dependency vulnerabilities.
- Post-v1 research phase Z: complete with green CI on 2026-10-09 and delivered
  in merged PR #28. `fermionicKernel_cross_product` and
  `fermionicKernel_residualUpdate_factor` expose the kernel's exact
  multiplicative-Monge structure. The orientation-free
  `fermionicKernel_residualUpdate_le_crossArea` bounds a one-step residual by
  its dimensionless time-frequency cross area, and
  `fermionicKernel_residualUpdate_le_half_of_crossArea` gives half contraction
  below the exact `log 2` area threshold, including both reflected corner
  orientations. This is a local theorem for the original kernel, not a claim
  that later GECP residuals preserve the same cross ratio; a global dyadic
  covering or residual-dominance theorem, Conjecture G1, and Problem 4.2
  remain open. The root command passes with 60 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase AA: complete with green CI on 2026-10-09 and
  delivered in merged PR #29. `fermionicKernel_eq_rpow` changes to the
  positive coordinate `z = exp(-omega)`, and
  `fermionicKernel_twoCornerResidual_eq_powerSecantError` identifies the exact
  nested residual after prescribed pivots `(0,omegaHigh)` and `(1,omegaLow)`
  with the normalized concave-power secant error on every strictly ordered
  asymmetric band. `fermionicKernel_twoCornerResidual_nonneg` proves its sign,
  while `fermionicKernel_twoCornerResidual_mem_Icc_minCrossArea` gives the
  two-sided tent bound by the minimum of the upper and reflected lower corner
  cross areas. This does not identify those prescribed pivots as complete on
  every subband, control the actual third pivot, prove a dyadic covering,
  establish Conjecture G1, or solve Problem 4.2. The root command passes with
  61 tests and no known dependency vulnerabilities.
- Post-v1 research phase AB: complete with green CI on 2026-10-09 and
  delivered in merged PR #30. The exact rational-coordinate witness
  `a=24^-3`, `b=24^3`, `z=(5/4)^3`, and `t=2/3` is transported through the
  Phase AA secant identity to the actual symmetric fermionic residual at
  cutoff `log 13824`. `fermionicKernel_twoCornerResidual_counterexample_value`
  evaluates that residual as `4495348/8973531`, and
  `powerSecant_twoThirds_counterexample_margin` proves its excess over half
  the initial pivot is the positive rational `222676/224338275`.
  `fermionicKernel_twoCornerResidual_not_half_contraction` packages the
  in-domain witness as a formal refutation of uniform two-pivot half
  contraction. This obstruction does not refute scale-dependent longer pivot
  blocks, Conjecture G1, or Problem 4.2; it shows that intermediate localized
  pivots are necessary at large cutoff. The root command passes with 62 tests
  and no known dependency vulnerabilities.
- Post-v1 research phase AC: complete with green CI on 2026-10-09 and
  delivered in merged PR #31. The exact one-corner comparisons and reflection
  identity bound the symmetric two-corner residual by `1/4` throughout both
  outer half-frequency bands. At the fixed center,
  `fermionicKernel_symmetricTwoCornerResidual_center` computes the residual as
  `1/2 - 1/(2*cosh(Lambda/2))`, which exceeds `1/4` for
  `2*log 4 <= Lambda`. The resulting
  `fermionicKernel_thirdCompletePivot_frequency_lt_halfCutoff` theorem uses
  the repository's `CompletePivotOn` predicate to force every actual third
  complete pivot into `abs omega < Lambda/2`, independently of ties. This is
  one rigorous dyadic localization step, not an iteration theorem for later
  residuals; Conjecture G1 and Problem 4.2 remain open. The root command passes
  with 63 tests and no known dependency vulnerabilities.
- Post-v1 research phase AD: complete with green CI on 2026-10-09 and
  delivered in merged PR #32.
  `residualUpdate_le_stripBound_of_signCoherent` proves that a local
  column-strip bound survives one sign-coherent update when complete-pivot
  maximality controls the selected-column factor, even when the pivot lies
  outside the strip. `stripBound_preserved_of_signCoherentCompletePivot`
  iterates the estimate, and `fermionicKernel_outerHalfBound_preserved`
  specializes it to keep the Phase AC outer-half bound `|R_n| <= 1/4` through
  every later realized exact complete-pivot residual. Therefore
  `fermionicKernel_laterCompletePivot_frequency_lt_halfCutoff_of_quarter_lt`
  forces every later pivot above one quarter into `|omega| < Lambda/2`. This
  is a fixed-strip invariant, not a shrinking-band recurrence or contraction
  below one quarter; Conjecture G1 and Problem 4.2 remain open. The root
  command passes with 64 tests and no known dependency vulnerabilities.
- Post-v1 research phase AE: complete with green CI on 2026-10-09 and
  delivered in merged PR #33.
  `Matrix.abs_det_le_of_factors_approx_except` and its closed coarse corollary
  prove that finitely many exactly represented columns add precisely their
  cardinality to the effective low-rank dimension while all other columns
  retain the approximation-error power. The fermionic specializations
  `fermionicKernel_sample_det_le_separatedApprox_except` and
  `fermionicKernel_sample_det_le_two_pow_except` require the dyadic cutoff
  only for nonexceptional frequencies. This removes the determinant-level
  obstruction to treating the two original cutoff corners as exceptions and
  all Phase AD-localized continuation columns at the smaller central scale;
  composition with the actual GECP selected core, a shrinking-band
  recurrence, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 65 tests and no known dependency vulnerabilities.
- Post-v1 research phase AF: complete with green CI on 2026-10-09 and
  delivered in merged PR #34.
  `fermionicKernel_gecp_error_pow_le_separatedApprox_except` and
  `fermionicKernel_gecp_error_pow_le_two_pow_except` compose the Phase AE
  determinant bounds with the actual selected core of a realized fermionic
  GECP prefix. The complete-pivot domain and evaluated residual may use a
  larger dyadic scale than the separated approximation applied to every
  nonexceptional selected column, and the closed exponent charges the exact
  exceptional cardinality. Packaging the first two recursive corner indices,
  deriving all remaining selected-column bounds from Phase AD, a shrinking-
  band recurrence, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 66 tests and no known dependency vulnerabilities.
- Post-v1 research phase AG: complete with green CI on 2026-10-09 and
  delivered in merged PR #35. `symmetricTwoCornerRun` prepends the
  exact prescribed cutoff pivots to any continuation of the two-corner
  residual, while `symmetricTwoCornerExceptional` identifies their two
  recursive selected indices and has proved cardinality two. Every
  nonexceptional selected column inherits the continuation's domain, and
  `fermionicKernel_symmetricTwoCornerRun_sample_det_le_two_pow` therefore
  applies the smaller-band determinant estimate with the exact rank shift
  `+2`. Deriving smaller-band `Run.CompleteOn` from Phase AD's above-quarter
  pivot condition, a shrinking-band recurrence, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 67 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AH: complete with green CI on 2026-10-09 and
  delivered in merged PR #36.
  `Run.ofResidualSequenceFrom` turns an indexed exact residual recurrence into
  a finite dependent run and proves its final residual, consecutive pivots,
  and exact selected-index cardinality. The symmetric two-corner specialization
  realizes every continuation residual from the original fermionic kernel.
  `fermionicKernel_aboveQuarterPrefix_completeOn_centralHalf` lifts Phase AD's
  pointwise localization to recursive central-half completeness, and
  `fermionicKernel_aboveQuarterPrefix_sample_det_le_two_pow` composes that fact
  with Phase AG's two-exception determinant bound. An explicit stopping bound
  for the above-quarter prefix, a repeated shrinking-band recurrence,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 68
  tests and no known dependency vulnerabilities.
- Post-v1 research phase AI: complete with green CI on 2026-10-09 and
  delivered in merged PR #37. The
  exact pivot-product identity and strict lower bounds on the two prescribed
  corner pivots turn every all-above-quarter continuation into a strict
  selected-core determinant lower bound. A binary-size/exponent-budget lemma
  contradicts Phase AH's smaller-band determinant upper bound at explicit
  parameters `p=24(s+1)` and continuation length `2048(s+1)^2`. Consequently,
  when `Lambda/2=2^s`, some continuation pivot before that length is at most
  `1/4`, and complete-pivot maximality bounds its whole cutoff-domain residual
  by `1/4`. A threshold-parametric nested-band restart, removal of the
  logarithmic startup, Conjecture G1, and Problem 4.2 remain open. The root
  command passes with 69 tests and no known dependency vulnerabilities.
- Post-v1 research phase AJ: complete with green CI on 2026-10-10 and
  delivered in merged PR #38. `Run.augmentedCore_step_det` and
  `Run.augmentedCore_det_eq_selectedCore_det_mul_finalResidual_minor` extend
  the one-border Schur identity to arbitrary finite residual minors while
  preserving the dependent pivot order. Under injectivity of the combined
  selected and border coordinates,
  `GECP.strictSignRegular_finalResidual_minor_oriented_pos` transfers the
  exact strict sign through explicit row and column orientation factors, and
  the generic and fermionic corollaries prove every such residual minor is
  nonzero. This does not yet factor the orientation into pointwise residual
  gauges, produce a threshold-parametric nested-band restart, establish
  Conjecture G1, or solve Problem 4.2. The root command passes with 70 tests
  and no known dependency vulnerabilities.
- Post-v1 research phase AK: complete with green CI on 2026-10-10 and
  delivered in merged PR #39.
  `anchoredSignGaugeKernel_pos_of_crossProductSignCoherent` converts
  the residual's balanced cross-product signs into anchor-based row and column
  signs that make every nonzero entry positive, while
  `abs_anchoredSignGaugeKernel` proves the transformation preserves entry
  magnitudes and `anchoredSignGaugeKernel_minor_det` gives its exact action on
  every finite minor. `Run.FreshRow`, `Run.FreshColumn`, and the Phase AJ
  nonvanishing theorem discharge the nonzero conditions for strictly
  sign-regular runs, yielding the fermionic specialization
  `fermionicKernel_finalResidual_anchoredSignGauge_pos`. This is order-one
  checkerboard normalization, not yet a compound-minor positivity theorem or
  a sharper determinant bound; it does not prove Conjecture G1 or solve
  Problem 4.2. The root command passes with 71 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase AL: complete with green CI on 2026-10-10 and
  delivered in merged PR #40. `Run.FreshRows`, `Run.FreshColumns`, and
  `Run.finalResidualMinorKernel` represent every fixed-order family of fresh
  final-residual minors as a compound kernel. The theorem
  `strictSignRegular_finalResidualMinorKernel_crossProductSignCoherent`
  multiplies four Phase AJ orientation-aware inequalities so all tuple signs
  cancel in squares, while the anchor-gauge and absolute-value corollaries
  produce a positive compound kernel without changing minor magnitudes. The
  fermionic specializations hold for every successful run. This does not yet
  derive a quantitative condensation inequality, remove the determinant
  prefactor, prove Conjecture G1, or solve Problem 4.2. The root command passes
  with 72 tests and no known dependency vulnerabilities.
- Post-v1 research phase AM: complete with green CI on 2026-10-10 and
  delivered in merged PR #41.
  `Run.augmentedCore_fin_two_mul_selectedCore_det` proves the exact
  selected-core Desnanot--Jacobi identity for two arbitrary appended rows and
  columns, while
  `strictSignRegular_augmentedCore_fin_two_condensation` uses bordered-minor
  sign coherence to replace the generic triangle factor two by a factor-one
  maximum. The fermionic specialization holds for every successful run. This
  homogeneous relation still requires an independent one-border-minor or
  stable-interpolation magnitude estimate before it can improve determinant
  decay; Conjecture G1 and Problem 4.2 remain open. The root command passes
  with 73 tests and no known dependency vulnerabilities.
- Post-v1 research phase AN: complete with green CI on 2026-10-10 and
  delivered in merged PR #42.
  `strictSignRegular_continuation_signCoherent` transfers selected-cross sign
  coherence through any successful prefix and finite continuation, while
  `strictSignRegular_completeContinuation_pivotsBounded` preserves an initial
  uniform envelope under complete pivoting. Consequently,
  `strictSignRegular_completeContinuation_selectedCore_det_le_pow` bounds the
  continuation determinant by `B^m` without the generic factorial, and the
  fermionic specialization holds for every successful prefix. An exact
  sign-regular witness shows the corresponding interpolation core need not be
  1-dominant. Applying this mechanism to the dyadic approximation still
  requires a compatible sign-regular remainder theorem; Conjecture G1 and
  Problem 4.2 remain open. The root command passes with 75 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AO: complete with green CI on 2026-10-10 and
  delivered in merged PR #43.
  `dyadicTaylorError` instantiates the positive-frequency
  error of the existing hard-masked dyadic Taylor construction, and
  `dyadicTaylorError_not_strictSignRegularAtOrder_two` proves an exact
  obstruction at its first mask transition. For rows `(1/2, 3/5)` and columns
  `(2, 5/2)`, the ordered error minor is positive instead of having the
  negative sign `expKernelSignature 2`. Hence the Phase AN factorial-free
  route cannot be applied directly to this approximation remainder. This does
  not exclude smooth overlaps, unmasked local tails, another compatible
  decomposition, Conjecture G1, or Problem 4.2. The root command passes with
  76 tests and no known dependency vulnerabilities.
- Post-v1 research phase AP: complete with green CI on 2026-10-10 and
  delivered in merged PR #44.
  `eighthOrderTail_eq_betaMoment` gives an exact
  unit-interval beta-weighted Laplace representation of the unmasked
  eighth-order Taylor tail, `eighthOrderTail_pos` proves strict positivity,
  and `eighthOrderTail_fin_two_det_neg_iff_betaMoment` reduces its order-two
  sign question to the normalized beta-moment kernel without changing the
  sign. The remaining beta-moment inequality, all-orders compatibility, a
  global smooth decomposition, Conjecture G1, and Problem 4.2 remain open.
  The root command passes with 77 tests and no known dependency vulnerabilities.
- Post-v1 research phase AQ: complete with green CI on 2026-10-10 and
  delivered in merged PR #45.
  `eighthOrder_unweighted_moment_gap` proves the exact
  `1/120` polynomial integral,
  `eighthOrderBetaMoment_two_mul_two_lt_one` proves the strict tilted first/
  second-moment gap for all nonnegative arguments, and
  `eighthOrderBetaMoment_variance_numerator_lt` closes the corresponding
  logarithmic-elasticity numerator on `0 <= x <= 2`. Differentiation, the
  order-two beta-moment kernel theorem, all-orders compatibility, Conjecture
  G1, and Problem 4.2 remain open. The root command passes with 78 tests and no
  known dependency vulnerabilities.
- Post-v1 research phase AR: complete with green CI on 2026-10-10 and
  delivered in merged PR #46. The public derivative theorems
  `eighthOrderBetaMoment_hasDerivAt` and
  `eighthOrderBetaMomentOne_hasDerivAt` justify differentiation under the beta
  integral, while `eighthOrderBetaElasticity_hasDerivAt` exposes the exact
  variance numerator and `eighthOrderBetaElasticity_strictAntiOn` proves strict
  decrease on `Set.Icc 0 2`. The final two-variable beta-moment determinant
  theorem, all-orders compatibility, Conjecture G1, and Problem 4.2 remain
  open. The root command passes with 79 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase AS: complete with green CI on 2026-10-10 and
  delivered in merged PR #47. The scale-ratio theorem
  `eighthOrderBetaMomentScaleRatio_strictMonoOn` transfers the Phase AR
  elasticity order to positive product-bounded frequency intervals,
  `eighthOrderBetaMoment_fin_two_det_neg` proves the normalized beta-moment
  determinant is strictly negative, and `eighthOrderTail_fin_two_det_neg`
  transfers that sign to the unmasked eighth-order Taylor tail. All-orders
  compatibility, a global smooth remainder, Conjecture G1, and Problem 4.2
  remain open. The root command passes with 80 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase AT: complete with green CI on 2026-10-10 and
  delivered in merged PR #48.
  `eighthOrderBetaMoment_scaled_eq_truncatedPower` proves the exact scaled beta
  identity, `eighthOrderTruncatedPower_integral_commonCutoff` moves every
  column to a shared integration domain, and
  `eighthOrderTail_eq_truncatedPowerComposition_commonCutoff` expresses the
  unmasked tail as a positive row scaling of an exponential/truncated-power
  composition. Total positivity of the truncated-power factor, continuous
  Cauchy--Binet strictness, all-orders tail sign regularity, Conjecture G1, and
  Problem 4.2 remain open. The root command passes with 81 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AU: complete with green CI on 2026-10-10 and
  delivered in merged PR #49.
  `Matrix.signature_mul_det_mul_rect_nonneg` proves nonnegativity of a signed
  rectangular product determinant from the common orientation of every
  ordered Cauchy--Binet summand,
  `Matrix.signedOrderedCauchyBinetTerm_le_signature_mul_det_mul_rect` retains
  any selected summand as a lower bound, and
  `Matrix.signature_mul_det_mul_rect_pos` obtains strictness from one positive
  witness. Total nonnegativity of the truncated-power factor, continuous
  Cauchy--Binet, all-orders tail sign regularity, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 82 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AV: complete with green CI on 2026-10-10 and
  delivered in merged PR #50.
  `GECP.TotallyNonnegativeAtOrder` and `GECP.TotallyNonnegative` add the weak
  ordered-minor vocabulary needed by the truncated-power factor;
  `eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_one` and
  `eighthOrderTruncatedPowerKernel_totallyNonnegativeAtOrder_two` prove its
  first two orders; and
  `eighthOrderTruncatedPowerKernel_fin_two_det_pos_iff` gives the exact paired
  support criterion for strictness. The all-orders truncated-power theorem,
  continuous Cauchy--Binet, all-orders tail sign regularity, Conjecture G1,
  and Problem 4.2 remain open. The root command passes with 83 tests and no
  known dependency vulnerabilities.
- Post-v1 research phase AW: complete with green CI on 2026-10-10 and
  delivered in merged PR #51.
  `normalizedTruncatedPowerKernel_zero` fixes the elementary strict order
  kernel at degree zero,
  `normalizedTruncatedPowerKernel_succ_integral` proves the exact normalized
  successor recurrence on its natural interval, and
  `normalizedTruncatedPowerKernel_succ_integral_commonDomain` extends that
  recurrence to any common containing interval. The existing seventh-degree
  factor is `7!` times the normalized degree-seven member. All-orders total
  nonnegativity, continuous Cauchy--Binet, all-orders tail sign regularity,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 84
  tests and no known dependency vulnerabilities.
- Post-v1 research phase AX: complete with green CI on 2026-10-10 and
  delivered in merged PR #52.
  `strictOrderKernel_minor_det_eq_zero_or_one` proves that every ordered minor
  of the normalized degree-zero base is exactly zero or one, while
  `strictOrderKernel_totallyNonnegative` and
  `normalizedTruncatedPowerKernel_zero_totallyNonnegative` package the
  all-orders consequence. Higher-degree total nonnegativity, continuous
  Cauchy--Binet, all-orders tail sign regularity, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 85 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AY: complete with green CI on 2026-10-11 and
  delivered in merged PR #53.
  `Matrix.detUpdateRowCLM` exposes the determinant as a continuous linear
  functional of one updated row,
  `Matrix.det_updateRow_intervalIntegral` commutes it with a Bochner interval
  integral, and the `_nonneg` and `_pos` corollaries transfer pointwise weak
  and strict signs. Multi-row continuous composition, higher-degree total
  nonnegativity, all-orders tail sign regularity, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 86 tests and no known
  dependency vulnerabilities.
- Post-v1 research phase AZ: complete with green CI on 2026-10-11 and
  delivered in merged PR #54.
  `TotallyNonnegativeAtOrder.det_minorMatrix_nonneg_of_monotone` covers weakly
  ordered samples by repeated-row or repeated-column degeneracy,
  `det_rowIntervalIntegral_nonneg` iterates the one-row determinant integral
  theorem over a finite set of rows, and
  `det_rowIntervalIntegral_kernel_nonneg` transfers total nonnegativity through
  ordered disjoint row intervals. The adjacent-row common-upper reduction,
  higher-degree total nonnegativity, all-orders tail sign regularity,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 87
  tests and no known dependency vulnerabilities.
- Post-v1 research phase BA: complete with green CI on 2026-10-11 and
  delivered in merged PR #55.
  `Matrix.det_successiveRowDifference` and
  `Matrix.det_commonUpperRowIntegral_eq_adjacent` reduce nested common-upper
  row integrals to ordered adjacent intervals,
  `det_commonUpperRowIntegral_kernel_nonneg` transfers total nonnegativity,
  and `normalizedTruncatedPowerKernel_succ_totallyNonnegativeAtOrder` raises
  the normalized degree. Consequently
  `normalizedTruncatedPowerKernel_totallyNonnegative` proves every degree and
  `eighthOrderTruncatedPowerKernel_totallyNonnegative` proves the all-orders
  seventh-degree factor theorem. Exact strict support, all-orders tail sign
  regularity, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 88 tests and no known dependency vulnerabilities.
- Post-v1 research phase BB: complete with green CI on 2026-10-11 and
  delivered in merged PR #56.
  `Matrix.sum_compositionSampleIntegrand_comp_perm` proves that finite
  permutation symmetrization of the row-wise composition integrand equals the
  product of the two sampled determinants, while
  `eighthOrderComposition_symmetrized_nonneg` gives the prescribed weak sign
  on every ordered intermediate sample using exponential strict sign
  regularity and the all-orders truncated-power theorem. Product-measure
  integration, strict support, all-orders tail sign regularity, Conjecture G1,
  and Problem 4.2 remain open. The root command passes with 89 tests and no
  known dependency vulnerabilities.
- Post-v1 research phase BC: complete with green CI on 2026-10-11 and
  delivered in merged PR #57.
  `Matrix.integral_comp_perm_pi` proves that a finite product of one common
  measure is invariant under coordinate permutations, and
  `Matrix.factorial_mul_det_integral_eq_integral_det_mul_det` packages the
  resulting finite-order Andreief identity without division. The proof also
  identifies the row-wise composition integral with the determinant of the
  one-dimensional pair integrals. Arbitrary-sample sign transfer, strict
  support, all-orders tail sign regularity, Conjecture G1, and Problem 4.2
  remain open. The root command passes with 90 tests and no known dependency
  vulnerabilities.
- Post-v1 research phase BD: complete with green CI on 2026-10-11 and
  delivered in merged PR #58.
  `Matrix.signed_det_mul_det_nonneg_of_strictMono` sorts injective samples by
  one common permutation, whose two determinant signs cancel, and handles
  repeated samples by determinant vanishing. Combined with the Phase BC
  Andreief identity, `Matrix.signature_mul_det_integral_nonneg` transfers
  ordered-sample weak signs through finite integral-kernel composition.
  Specializing to the exponential and seventh-degree truncated-power factors
  proves `eighthOrderTail_minor_signed_nonneg`: every ordered eighth-order
  Taylor-tail minor on positive times and nonnegative frequencies has the
  prescribed weak sign. Strict support, strict all-orders tail sign
  regularity, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 91 tests and no known dependency vulnerabilities.
- Post-v1 research phase BE: verified locally on 2026-10-11.
  `interlacingSample` places each intermediate coordinate strictly between
  adjacent positive frequency knots. The resulting truncated-power matrix is
  upper triangular with positive diagonal at every order. A continuous
  compactly supported bump promotes this explicit strict point to positive
  product-measure support, so Andreief gives
  `eighthOrderComposition_interval_det_pos`. Positive `t^8/7!` row scaling
  then proves `eighthOrderTail_minor_signed_pos`: every positive-order tail
  minor on strictly ordered positive times and frequencies has the prescribed
  strict sign. Translating this new compatible remainder structure into a
  cutoff-uniform GECP block contraction, Conjecture G1, and Problem 4.2 remain
  open. The root command passes with 92 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #59.
- Post-v1 research phase BF: verified locally on 2026-10-11.
  `positiveEighthOrderTailKernel_strictSignRegular` and its denominator-scaled
  fermionic counterpart turn the Phase BE ordered-minor theorem into global
  strict-sign-regular kernels on positive-real subtypes. The exact identity
  `fermionicKernel_eq_positiveEighthOrderTaylorApprox_add_tail` splits the
  positive-frequency kernel into eight separated Taylor terms and the
  compatible tail. On the positive unit band the tail is bounded by `1 / 8!`,
  and `positiveEighthOrderFermionicTailKernel_completeRun_selectedCore_det_le`
  gives the factorial-free complete-run determinant envelope
  `(1 / 8!)^rank`. Mixed determinants of the sum, the order-`8p` and dyadic
  extensions, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 93 tests and no known dependency vulnerabilities. The change
  passed CI and merged as PR #60.
- Post-v1 research phase BG: verified locally on 2026-10-11.
  The direct negative-exponential split has an order-two mixed-coefficient
  sign obstruction. The equivalent reflected identity
  `fermionicKernel_one_sub_eq_positive_exp` instead exposes the positive
  series `exp(u*omega)/(1+exp(omega))`. Eight positive monomial features give
  the low-rank component, while
  `positiveSeriesEighthOrderTailKernel_strictSignRegular` proves the tail is
  strictly totally positive. Denominator scaling preserves that structure and
  the unit-band `1 / 8!` error and `(1 / 8!)^rank` complete-run determinant
  envelopes. A 140-digit sweep observes common-sign mixed coefficients
  through order six, but the coefficientwise finite Cauchy--Binet theorem,
  dyadic assembly, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 94 tests and no known dependency vulnerabilities. The change
  passed CI and merged as PR #61.
- Post-v1 research phase BH: verified locally on 2026-10-11.
  `signature_mul_det_weightedFeatureProduct_eq_sum` proves the exact finite
  block-feature Cauchy--Binet expansion, and
  `expTaylor_blockFeatureCoefficient_nonneg` proves every fixed high-feature
  coefficient nonnegative for positive exponential Taylor features on
  nonnegative ordered nodes. This closes the finite coefficientwise gap from
  Phase BG. Passing to the exact infinite tail, extracting a quantitative
  mixed-determinant bound, dyadic assembly, Conjecture G1, and Problem 4.2
  remain open. The root command passes with 95 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #62.
- Post-v1 research phase BI: verified locally on 2026-10-11.
  `blockFeaturePolynomial_eq_mixedColumnPolynomial` identifies the Phase BH
  feature-count expansion with determinant multilinearity in physical columns,
  and `expTaylor_mixedColumnCoefficient_nonneg` proves every finite positive-
  Taylor fixed-cardinality mixed-column sum nonnegative. The exact infinite-
  tail limit, quantitative mixed-determinant bound, dyadic assembly,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 95
  tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #63.
- Post-v1 research phase BJ: verified locally on 2026-10-11.
  `mixedColumnCoefficient_tendsto` supplies the finite determinant-continuity
  passage, the low/high feature products are identified with finite eighth-
  order Taylor tails, and
  `positiveEighthOrderTail_mixedColumnCoefficient_nonneg` proves every exact
  infinite-tail mixed coefficient sum nonnegative. A quantitative bound on
  these coefficients, dyadic assembly, Conjecture G1, and Problem 4.2 remain
  open. The root command passes with 96 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #64.
- Post-v1 research phase BK: verified locally on 2026-10-11.
  `mixedColumnPolynomial_rank_amplification` combines low-rank coefficient
  vanishing with nonnegative mixed coefficients to prove
  `lambda^(n-k) det(E+LR) <= det(lambda E+LR)` for `lambda >= 1`.
  `positiveSeriesFermionic_rank_amplification` supplies the exact rank-eight
  fermionic specialization. A fixed-base amplified determinant envelope,
  dyadic assembly, Conjecture G1, and Problem 4.2 remain open. The root command
  passes with 97 tests and no known dependency vulnerabilities. The change
  passed CI and merged as PR #65.
- Post-v1 research phase BL: verified locally on 2026-10-11.
  `mixedColumnPolynomial_det_mono` turns the exact coefficient signs into
  monotonicity of the tail-weight determinant pencil.
  `positiveSeriesAmplifiedFermionicKernel_strictSignRegular` proves that every
  weight at least one is strictly totally positive, while
  `positiveSeriesAmplifiedFermionicKernel_factorial_abs_le_two` and its
  complete-run corollary give a unit-band entry bound of `2` and selected-core
  bound `2^n` at weight `8!`. The arbitrary ordered-sample determinant bridge,
  GECP composition, dyadic assembly, Conjecture G1, and Problem 4.2 remain
  open. The root command passes with 98 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #66.
- Post-v1 research phase BM: verified locally on 2026-10-11.
  `Matrix.det_le_prod_diagonal_of_strictlyTotallyPositive` proves the generic
  diagonal-product envelope by northwest Schur induction, and
  `minorMatrix_det_le_pow_of_strictSignRegular_one` transfers it to arbitrary
  ordered samples of a positive-minor kernel. Consequently
  `positiveSeriesAmplifiedFermionicKernel_factorial_minor_det_le_two_pow`
  bounds every unit-band amplified sample determinant by `2^n`, while
  `positiveSeriesFermionic_det_le_two_pow_div_factorial_pow` gives
  `det K_n <= 2^n / (8!)^(n-8)` for the original positive-series matrix when
  `n >= 8`. GECP composition, dyadic band/reflection assembly, Conjecture G1,
  and Problem 4.2 remain open. The root command passes with 99 tests and no
  known dependency vulnerabilities. The change passed CI and merged as PR #67.
- Post-v1 research phase BN: verified locally on 2026-10-12.
  `positiveSeriesFermionicKernel_abs_minor_det_le` transfers the Phase BM
  determinant decay to injective samples in arbitrary algorithmic order, and
  `positiveSeriesFermionicKernel_gecp_error_pow_le` composes it with the exact
  selected-core residual-power theorem. The explicit corollary
  `positiveSeriesFermionicKernel_gecp_error_le_half_ten` proves that ten exact
  complete pivots reduce every residual entry to at most one half on the
  reflected positive unit band. Full signed-frequency assembly, endpoint
  closure, dyadic rescaling, Conjecture G1, and Problem 4.2 remain open. The
  root command passes with 100 tests and no known dependency vulnerabilities.
  The change passed CI and merged as PR #68.
- Post-v1 research phase BO: verified locally on 2026-10-12.
  The affine frequency shift `z = 1 - omega` converts the full closed signed
  unit band into one positive exponential product band. Tail amplification at
  weight `128`, strict total positivity, and the diagonal-product determinant
  envelope prove `|det K_n| <= 2^n / 128^(n-8)` for every injective signed-band
  sample. `fermionicKernel_signedUnit_gecp_error_le_half_twelve` then proves
  that twelve exact complete pivots contract the original fermionic residual
  uniformly by one half on `[0,1] x [-1,1]`. Dyadic rescaling/continuation,
  Conjecture G1, and Problem 4.2 remain open. The root command passes with 101
  tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #69.
- Post-v1 research phase BP: verified locally on 2026-10-12.
  The center-parametric coordinate `z = center + 1 - omega` extends Phase BO
  from `[-1,1]` to every translated width-two frequency band without changing
  its constants. The weight-`128` amplified kernel remains strictly totally
  positive and bounded by two independently of `center`, giving
  `|det K_n| <= 2^n / 128^(n-8)` for arbitrary injective samples and a uniform
  twelve-pivot half contraction on
  `[0,1] x [center-1,center+1]`. Restarting this theorem after an earlier
  selected core, the shrinking-band recurrence, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 102 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #70.
- Post-v1 research phase BQ: complete; merged as PR #71 with green CI.
  `mixedColumnPolynomial_rank_amplification_abs` proves the absolute
  determinant amplification law, while `finalResidual_rank_amplification` and
  `augmentedCore_rank_amplification` show that an exact compatible split of
  the final Schur residual keeps exponent `n-k` independent of the selected
  prefix. The exact rational theorem
  `frozenColumnTaylor_rankAmplification_fails` proves that amplifying only
  fresh-column tails can instead reverse determinant sign, so the residual-
  level formulation is necessary. Constructing the compatible fermionic
  residual split, dyadic iteration, Conjecture G1, and Problem 4.2 remain open.
  The root command passes with 103 tests and no known dependency
  vulnerabilities. The change passed CI and merged as PR #71.
- Post-v1 research phase BR: complete; merged as PR #72 with green CI.
  `GECP.Run.conditionRows` applies the coefficients of a fixed GECP prefix
  linearly to an arbitrary companion matrix, and its additivity, scalar, and
  right-multiplication laws prove the exact transport identity
  `finalResidual = conditionRows(E) + conditionRows(L)*R` whenever
  `K = E + L*R`. The positive-series fermionic specialization therefore
  preserves its rank-eight right feature block through every prefix and
  composes with Phase BQ under the remaining conditioned-tail sign premise.
  An exact three-pivot regression verifies the identity and low-rank minors,
  but also finds a negative degree-three mixed coefficient, disproving the
  hoped-for automatic coefficientwise nonnegativity. A signed or real-rooted
  determinant-pencil comparison, dyadic iteration, Conjecture G1, and Problem
  4.2 remain open. The root command passes with 104 tests and no known
  dependency vulnerabilities. The change passed CI and merged as PR #72.
- Post-v1 research phase BS: complete; merged as PR #73 with green CI.
  Rank-eight coefficient vanishing reduces every canonical conditioned
  determinant pencil to degree at most eight after its forced zero is removed.
  Exact Lagrange interpolation at weights `56, 65, ..., 128` therefore gives
  a sign-free norming theorem with constant `518173`. At fresh order eighteen,
  nine amplified determinant envelopes of size `2^18` imply the exact
  physical determinant bound `2^-18`, and
  `positiveSeriesFermionicKernel_conditioned_det_le_halfPow` connects this
  directly to the Phase BR split. Exact regressions disprove real-rootedness
  and domination at the single weight `128`, while verifying the finite-weight
  certificate. The nine companion envelopes, dyadic iteration, Conjecture G1,
  and Problem 4.2 remain open. The root command passes with 105 tests and no
  known dependency vulnerabilities. The change passed CI and merged as PR #73.
- Post-v1 research phase BT: complete; merged as PR #74 with green CI.
  `GECP.companionBiUpdate` applies both fixed oblique pivot
  projectors to a companion, recovers the ordinary Schur update under
  self-application, annihilates the selected row and column, and transports
  `L*R` as separately conditioned feature factors. Its iteration through
  `GECP.Run.conditionBoth` gives the exact two-sided rank-eight fermionic
  residual split and exposes the conditioned low determinant as an ordered
  sum of matched feature-minor products. An exact `Fraction` regression keeps
  the one-sided mixed-sign witness but finds common orientation for every
  nonzero two-sided coefficient in 2,760 exhaustively enumerated deterministic
  samples. This does not prove all-orders matched-minor sign coherence, the
  nine Phase BS companion envelopes, residual half contraction, the dyadic
  recurrence, Conjecture G1, or Problem 4.2. The root command passes a
  3,568-job Lean build with 106 tests and no known dependency vulnerabilities.
  The change passed CI and merged as PR #74.
- Post-v1 research phase BU: complete; merged as PR #75 with green CI.
  Left and right feature-augmented core determinants factor exactly
  as `selectedCore.det` times the corresponding conditioned feature minor.
  Their matched product therefore contains `selectedCore.det ^ 2`, and
  `positiveSeriesConditionedLow_signedTerm_nonneg_iff` reduces every
  conditioned rank-eight Cauchy--Binet term sign to a product of two original
  bordered feature determinants. Exact `Fraction` arithmetic verifies 2,382
  left, right, and matched-product identities across multi-pivot runs and
  feature subsets. The remaining determinants form confluent exponential
  systems, but their all-orders sign theorem, the nine Phase BS companion
  envelopes, residual half contraction, the dyadic recurrence, Conjecture G1,
  and Problem 4.2 remain open. The root command passes a 3,570-job Lean build
  with 106 tests and no known dependency vulnerabilities. The change passed
  CI and merged as PR #75.
- Post-v1 research phase BV: complete; merged as PR #76 with green CI.
  `combinedFeatureLeft` and `combinedFeatureRight` adjoin fixed pivot blocks
  to arbitrary feature factors, and `weightedCombinedFeatureProduct` proves
  that marking those coordinates produces `mu • (U*V) + A*B`. Theorems
  `fixedBlockPluckerPairing_eq_mixedColumnCoefficient` and
  `fixedBlockPluckerPairing_eq_sum` identify the full-fixed-degree signed
  Plücker pairing with one physical mixed-column determinant coefficient.
  An exact three-pivot finite-Taylor witness has two negative and 54 positive
  bordered-feature products while their signed total is positive and equals
  that coefficient, disproving the Phase BU termwise-sign continuation while
  preserving its invariant aggregate. Positivity of the aggregate for
  arbitrary fermionic prefixes, the nine Phase BS companion envelopes,
  residual half contraction, the dyadic recurrence, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,571-job Lean build with
  107 tests and no known dependency vulnerabilities. The change passed CI and
  merged as PR #76.
- Post-v1 research phase BW: complete; merged as PR #77 with green CI.
  `borderedRankUpdatePencil_det` linearizes
  `det(mu • (U*V) + B)` as the determinant of the block matrix with corners
  `1`, `V`, `-mu • U`, and `B`, without invertibility assumptions. The mixed-
  column and Phase BV block-feature polynomials inherit this representation at
  every weight. Exact rational arithmetic verifies the identity at weight one
  and all nine Phase BS weights on the three-pivot witness. A negative signed
  four-cycle is invariant under row/column sign changes and proves that the
  bordered matrix cannot be gauged into an entrywise-nonnegative matrix, so a
  direct total-nonnegativity proof is unavailable. The nine companion
  envelopes, residual half contraction, the dyadic recurrence, Conjecture G1,
  and Problem 4.2 remain open. The root command passes a 3,572-job Lean build
  with 107 tests and no known dependency vulnerabilities. The change passed
  CI and merged as PR #77.
- Post-v1 research phase BX: complete; merged as PR #78 with green CI.
  `twoSidedCompanionSchur` gives the closed four-term projection formula for
  separate kernel and companion cross blocks around an invertible core.
  `doubleBorderedCompanionMatrix_det` and its absolute-value form prove that
  the corresponding double-bordered determinant factors by `det(C)^2`.
  Exact solve-based rational checks show that three iterative two-sided
  companion pivots agree entrywise with the closed formula at weight one and
  all nine Phase BS weights, and independently verify the determinant
  factorization. Negative order-two minors persist, so total nonnegativity is
  not obtained. Formal identification with arbitrary recursive
  `Run.conditionBoth`, the normalized double-border envelope, residual half
  contraction, the dyadic recurrence, Conjecture G1, and Problem 4.2 remain
  open. The root command passes a 3,573-job Lean build with 107 tests and no
  known dependency vulnerabilities. The change passed CI and merged as PR #78.
- Post-v1 research phase BY: complete; merged as PR #79 with green CI.
  `conditionRows_eq_selectedCoreSchur` and
  `conditionColumns_eq_selectedCoreSchur` prove the closed one-sided
  selected-core formulas, while
  `conditionBoth_eq_twoSidedCompanionSchur` identifies every arbitrary
  recursive two-sided run with the exact four-term Phase BX expression.
  `doubleBorderedCompanionMatrix_det_eq_conditionBoth_sample` specializes the
  squared-core determinant identity to arbitrary finite samples. Exact
  solve-based rational checks cover 144 deterministic randomized combinations
  of prefix length, companion weight, and fresh block. The normalized
  double-border envelope, residual half contraction, dyadic recurrence,
  Conjecture G1, and Problem 4.2 remain open. The root command passes a
  3,574-job Lean build with 107 tests and no known dependency vulnerabilities;
  the change passed CI and merged as PR #79.
- Post-v1 research phase BZ: complete; merged as PR #80 with green CI.
  `upperTangentMatrix` identifies the Phase BX double border as a
  distinguished minor of `[[K,H],[0,K]]`, while
  `upperTangentMatrix_mul_factorization` exposes the shared feature
  factorization through `[[1,W],[0,1]]`.
  `positiveSeriesFermionicKernel_det_le_halfPow_of_upperTangent` proves that
  the nine normalized tangent-minor bounds imply the Phase BS order-eighteen
  half-power estimate. Exact rational tests cover 24 randomized
  factorizations and 360 complete-prefix coefficient polynomials; all tested
  coefficients obey the explicit selected/fresh cross-inversion orientation,
  even though ordinary order-two minors can be negative. The oriented-minor
  theorem, normalized envelope, residual half contraction, dyadic recurrence,
  Conjecture G1, and Problem 4.2 remain open. The root command passes a
  3,577-job Lean build with 108 tests and no known dependency vulnerabilities;
  the change passed CI and merged as PR #80.
- Post-v1 research phase CA: complete; merged as PR #81 with green CI.
  `signature_mul_det_mul_mul_rect_eq_sum` expands an arbitrary rectangular
  three-factor determinant over two ordered intermediate flags, and
  `signature_mul_det_mul_mul_rect_nonneg_of_factor_orientations` composes
  separate left, middle, and right minor orientations. The upper-tangent
  specialization canonically reindexes both feature layers and exposes every
  sampled tangent determinant through those three minor families. Exact
  rational checks cover 24 shuffled rectangular triple expansions while
  retaining the Phase BZ coefficient and negative-minor guards. The actual
  exponential flag-minor signs, normalized tangent envelope, residual half
  contraction, dyadic recurrence, Conjecture G1, and Problem 4.2 remain open.
  The root command passes a 3,579-job Lean build with 108 tests and no known
  dependency vulnerabilities; the change passed CI and merged as PR #81.
- Post-v1 research phase CB: complete; merged as PR #82 with green CI.
  `reindex_upperTangent_diagonal_eq_blockDiagonal` identifies the diagonal
  tangent middle factor, after feature interleaving, with independent
  `2 x 2` cells `[[1,w],[0,1]]`.
  `tangentCell_totallyNonnegative` proves every local cell minor is
  nonnegative, while `upperTangentFeatureMiddle_oriented_minor_nonneg`
  transfers all-orders TN of the sparse interleaved matrix into the exact
  layer-major tuple orientation required by Phase CA. Exact rational checks
  exhaust 198,626 feature-major and 198,626 oriented layer-major minors
  through five cells. All-orders TN of the cell direct sum, the two outer
  flag-minor signs, normalized tangent envelope, Conjecture G1, and Problem
  4.2 remain open. The root command passes a 3,580-job Lean build with 109
  tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #82.
- Post-v1 research phase CC: complete; merged as PR #83 with green CI.
  `totallyNonnegative_of_nonnegative_of_noncrossing` proves a reusable
  determinant criterion: entrywise nonnegative kernels with noncrossing
  support are totally nonnegative at every order because all nonidentity
  Leibniz terms vanish. `interleavedTangentMiddle_totallyNonnegative`
  specializes it to arbitrary nonnegative diagonal tangent weights, and
  `upperTangentFeatureMiddle_diagonal_oriented_minor_nonneg` removes the
  conditional middle-factor premise from Phase CA. Exact checks cover mixed
  zero/positive weights, 225 crossing pairs, and 12,403 nonidentity matching
  terms while retaining the negative ordinary layer-major witness. The two
  outer flag-minor signs, normalized tangent envelope, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,581-job Lean build with
  110 tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #83.
- Post-v1 research phase CD: complete; merged as PR #84 with green CI.
  `phaseCD_negative_triple_term` proves that a positive-monomial,
  positive-diagonal three-feature tangent factorization has an ordered
  Cauchy--Binet summand equal to `-1/192`, and
  `phaseCD_no_factorwise_positive_orientation` formally rules out the
  separate factor-orientation interface used as Phase CA's sufficient
  criterion. The complete distinguished determinant is nevertheless proved
  equal to `24929/27648 > 0`. Exact enumeration finds two positive nonzero
  terms at two features and 10 positive plus 2 negative terms at three
  features, so the next route must preserve global cancellation rather than
  establish termwise signs. The normalized tangent envelope, order-eighteen
  contraction, dyadic recurrence, Conjecture G1, and Problem 4.2 remain open.
  The root command passes a 3,582-job Lean build with 111 tests and no known
  dependency vulnerabilities; the change passed CI and merged as PR #84.
- Post-v1 research phase CE: complete; merged as PR #85 with green CI.
  `upperTangentMatrix_factorized_distinguishedMinor_det` collapses every
  distinguished tangent minor of `K=A*B`, `H=A*W*B` to the squared selected-
  core determinant times one determinant of the separately conditioned
  feature factors. `signature_mul_conditionedFeature_det_eq_doubleFlagSum`
  identifies the entire ordered double-flag expansion with the same global
  determinant, preserving the cancellations that Phase CD proved essential.
  The finite positive-Taylor specialization exposes arbitrary diagonal degree
  weights directly. Exact rational checks cover 5,000 generalized monomial
  cases with arbitrary nonnegative weights and a signed-weight failure
  `-4/169`. The nonnegative diagonal-action sign theorem, normalized tangent
  envelope, order-eighteen contraction, dyadic recurrence, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,584-job Lean build with
  112 tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #85.
- Post-v1 research phase CF: complete; merged as PR #86 with green CI.
  `sparseFeatureProjection` packages a one-hot feature inclusion and its
  transpose as a nonnegative diagonal action.
  `det_mul_sparseFeatureProjection_mul` proves that inserting this projection
  isolates exactly one matched Plucker coordinate, and
  `mul_sparseFeatureProjection_mul_det_neg` converts a negative coordinate
  into a determinant sign failure. In the exact Phase BV witness the supports
  `(0,1,5)` and `(0,2,5)` both fail, refuting the broad nonnegative diagonal-
  action target suggested by Phase CE's dense census. All 90 contiguous
  cutoff/Phase-BS-weight combinations remain positive in that witness, so the
  physical step-weight family is the narrowed next target. Its grouped
  coefficient signs, the normalized tangent envelope, order-eighteen
  contraction, dyadic recurrence, Conjecture G1, and Problem 4.2 remain open.
  The root command passes a 3,585-job Lean build with 112 tests and no known
  dependency vulnerabilities; the change passed CI and merged as PR #86.
- Post-v1 research phase CG: complete; merged as PR #87 with green CI.
  `selectedCore_sq_mul_conditionedBlockFeatureCoefficient` identifies every
  conditioned tail-count coefficient with a grouped sum of original bordered
  Plucker products, while
  `selectedCore_sq_mul_conditionedStep_det_eq_sum` evaluates the entire
  contiguous-step determinant by weighting those same groups.
  `conditionedBlockFeatureCoefficient_nonneg_of_grouped` exposes the precise
  grouped sign premise without requiring the false coordinatewise signs. In
  the exact Phase BV witness all groups are nonnegative at every cutoff and 90
  direct Phase-BS-node evaluations equal the grouped polynomial values. The
  finite grouped sign theorem, infinite-tail passage, normalized tangent
  envelope, order-eighteen contraction, dyadic recurrence, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,586-job Lean build with
  112 tests and no known dependency vulnerabilities; the change passed CI and
  merged as PR #87.
- Post-v1 research phase CH: complete; merged as PR #88 with green CI.
  `phaseCH_middle_borderedProduct_neg` gives the minimal observed one-pivot
  positive-Taylor coordinate cancellation at the exact nodes
  `x=(1/8,1/4)`, `y=(6,22)`, with bordered products
  `561/64`, `-7/256`, and `801/32`.
  `phaseCH_proper_contiguous_groups_pos` proves that all proper nonempty
  prefix and suffix groups remain positive, and `phaseCH_firstStep_pos` and
  `phaseCH_secondStep_pos` prove the two nontrivial contiguous-step
  determinants positive for every nonnegative scalar weight. An independent
  exact census checks 6,468 two-feature cases without a negative coordinate
  and finds 800 three-feature cancellations, beginning with the formal
  witness. This isolates adjacent-degree compensation but does not prove the
  grouped sign for arbitrary nodes, prefixes, or fresh order. The infinite-
  tail passage, normalized tangent envelope, order-eighteen contraction,
  dyadic recurrence, Conjecture G1, and Problem 4.2 remain open. The root
  command passes a 3,587-job Lean build with 113 tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #88.
- Post-v1 research phase CI: complete; merged as PR #89 with green CI.
  `positiveTaylorThree_prefix_factor` and
  `positiveTaylorThree_suffix_factor` express both adjacent one-pivot,
  three-feature groups as the product of the ordered row gap, ordered column
  gap, and an explicit positive-coefficient polynomial.
  `positiveTaylorThree_prefix_pos` and
  `positiveTaylorThree_suffix_pos` prove strict positivity for arbitrary
  ordered nonnegative nodes, and the two step theorems preserve strict
  positivity for every nonnegative weight under their sharp boundary
  hypotheses. Exact arithmetic checks all coordinate and group factorizations
  on 2,198 boundary, super-physical, and randomized rational cases. This is
  the first general grouped-sign result, but it covers only one pivot and
  three features. Arbitrary feature count, compound fresh order, the infinite-
  tail passage, normalized tangent envelope, order-eighteen contraction,
  dyadic recurrence, Conjecture G1, and Problem 4.2 remain open. The root
  command passes a 3,588-job Lean build with 114 tests and no known dependency
  vulnerabilities; the change passed CI and merged as PR #89.
- Post-v1 research phase CJ: complete; merged as PR #90 with green CI.
  `two_mul_weightedCovarianceNumerator` and
  `weightedMass_mul_weightedCenteredGroup` formalize the weighted Chebyshev
  and conditional-covariance identities suggested by Phase CI.
  `weightedCenteredLowerGroup_nonneg` and
  `weightedCenteredUpperGroup_nonneg` prove both ordered block signs without
  division. `onePivotTaylorProduct_eq_geometricCenteredProduct` then rewrites
  every finite-Taylor bordered coordinate into the centered model, and the
  prefix, suffix, and step theorems prove all finite one-pivot contiguous
  grouped signs for arbitrary feature count and nonnegative step weight.
  Exact arithmetic checks 4,644 exhaustive small weighted groups and 6,782
  Taylor groups through 12 features, including 218 negative individual
  coordinates. Higher fresh order, the infinite-tail passage, normalized
  tangent envelope, order-eighteen contraction, dyadic recurrence, Conjecture
  G1, and Problem 4.2 remain open. The root command passes a 3,592-job Lean
  build with 115 tests and no known dependency vulnerabilities; the change
  passed CI and merged as PR #90.
- Post-v1 research phase CK: complete; merged as PR #91 with green CI.
  `det_add_smul_eq_prod_of_common_diagonal` factors a matrix
  pencil with common left/right factors and diagonal middle coordinates, and
  its nonnegative and strict variants turn a nonnegative generalized spectrum
  into an all-nonnegative-weight determinant sign without an inverse. Exact
  Phase BV data refute absolute-value `MTP2` in 372 of 3,136 subset-lattice
  comparisons while preserving nonnegative, log-concave, real-rooted
  threshold polynomials. A deterministic exact audit retains those signatures
  for 1,481 oriented polynomials from 200 compound Taylor cases through fresh
  order three, despite `MTP2` failures in 147 cases. The arbitrary-prefix
  generalized-spectrum theorem, infinite-tail passage, normalized tangent
  envelope, order-eighteen contraction, dyadic recurrence, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,593-job Lean build and
  116 tests.
- Post-v1 research phase CL: complete; merged as PR #92 with green CI.
  `finTwoPencil_det_expansion` gives the exact quadratic
  determinant pencil at fresh order two,
  `finTwoPencilDiscriminant_nonneg_iff` isolates the Lorentzian inequality
  `4*c0*c2 ≤ c1^2`, and `nonnegative_quadratic_factorization` constructs its
  two nonnegative spectral offsets. The matrix specialization yields an
  all-nonnegative-weight sign certificate without an inverse. Exact arithmetic
  checks 3,868 oriented threshold pencils from 500 deterministic finite-Taylor
  runs; all coefficients and discriminants pass, while 2,838 pencils fail both
  naive single-cross-minor sign criteria. The fermionic Lorentzian inequality,
  higher fresh order, infinite-tail passage, normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, and Problem
  4.2 remain open. The root command passes a 3,594-job Lean build and 117
  tests.
- Post-v1 research phase CM: complete; merged as PR #93 with green CI.
  `weightedMass_mul_weightedMass_mul_crossIncrement` proves an exact
  division-free two-block covariance decomposition, and
  `four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq` turns
  nonnegative weights, within-block comonotonicity, and cross-block order into
  the sharp Lorentzian inequality `4*c0*c2 ≤ c1^2`.
  `geometricCrossIncrement_lorentzian` packages the degree-prefix/suffix
  ratio-coordinate specialization needed by the next conditioned-factor
  bridge. Exact arithmetic checks 1,764 weighted splits and 11,024 oriented
  conditioned two-feature products with no failure. The representation of
  arbitrary conditioned Taylor rank-two factors in this model, the full
  fermionic fresh-order-two theorem, higher fresh order, infinite-tail
  passage, normalized tangent envelope, order-eighteen contraction, dyadic
  recurrence, Conjecture G1, and Problem 4.2 remain open. The root command
  passes a 3,594-job Lean build and 117 tests.
- Post-v1 research phase CN: complete; merged as PR #94 with green CI.
  `rankTwoMatchedCoordinate_eq_affine` factors matched left/right rank-two
  coordinates through chart weights and ratio increments. The within- and
  cross-block coefficient identities connect this representation to the
  Phase CM covariance theorem, while `rankTwoBlockMatrix_det` and
  `finTwoPencilMixed_rankTwoBlockMatrix` identify those coefficients with the
  original `2 x 2` determinant pencil.
  `rankTwoBlockPencil_lorentzian_of_affine` is the resulting division-free
  Phase CL certificate. The direct chart is not globally degree-monotone in
  458 of 500 exact cases and has opposite mean shifts in 276 of 2,868 splits.
  Nevertheless all 500 left/right configurations share an unoriented cyclic
  order and all 2,868 degree prefixes are circular intervals. Proving that
  projective-order invariant for arbitrary conditioned Taylor factors, then
  passing to the infinite tail and higher fresh order, remains open, as do
  the normalized tangent envelope, order-eighteen contraction, dyadic
  recurrence, Conjecture G1, and Problem 4.2. The root command passes a
  3,595-job Lean build and 117 tests.
- Post-v1 research phase CO: complete; merged as PR #95 with green CI.
  `finTwoPencilMixed_mul_mul` proves that a two-sided basis change scales the
  mixed coefficient by `det(P)det(Q)`, and
  `finTwoPencilDiscriminant_mul_mul` shows that the full discriminant scales
  by its square. Lorentzianity transfers through arbitrary changes and is
  equivalent under invertible ones.
  `rankTwoBlockPencil_lorentzian_of_projective_affine` composes this with the
  Phase CN chart certificate, allowing a different invertible affine chart
  for every cutoff. Exact arithmetic checks 512 invertible rational changes,
  including 271 negative scale products, plus a singular degeneration. The
  cyclic-prefix and separating-chart theorem for arbitrary conditioned Taylor
  factors, infinite-tail passage, higher fresh order, normalized tangent
  envelope, order-eighteen contraction, dyadic recurrence, Conjecture G1, and
  Problem 4.2 remain open. The root command passes a 3,595-job Lean build and
  118 tests.
- Post-v1 research phase CP: complete; merged as PR #96 with green CI.
  `finTwoPoleBasis_det` proves that the explicit pole basis
  `[[-c,1],[1,0]]` is invertible with determinant `-1`, and
  `finTwoPoleBasis_mul_rankTwoBlockMatrix_mul` computes its exact two-sided
  action on a rank-two feature block. `poleChartBase_mul_poleChartRatio`
  identifies the reciprocal affine coordinate away from pole collisions.
  `rankTwoBlockPencil_lorentzian_of_poles` composes these identities with the
  Phase CO projective transfer and the Phase CN ratio certificate. Exact
  rational checks cover exterior and interior pole choices, all five
  nontrivial cutoffs of a six-feature example, and a rejected collision. This
  closes the projective construction conditional on ordered reciprocal ratios;
  proving the conditioned-Taylor cyclic-prefix/separating-pole invariant,
  infinite-tail passage, higher fresh order, normalized tangent envelope,
  order-eighteen contraction, dyadic recurrence, Conjecture G1, and Problem
  4.2 remains open. The root command passes a 3,595-job Lean build and 119
  tests.
- Post-v1 research phase CQ: complete; merged as PR #97 with green CI.
  `finTwoScaledOrientedPoleBasis_det` and
  `finTwoScaledOrientedPoleBasis_mul_rankTwoBlockMatrix_mul` separate the two
  chart choices needed by the conditioned pencils: projective orientation and
  common base-coordinate sign. The matched coordinate scales by
  `scale²*orientation`, and
  `scaledOrientedPoleChartRatio_monovaryOn_of_matched_nonneg` derives
  reciprocal comonotonicity from positive transformed weights and oriented
  projective minor signs.
  `rankTwoBlockPencil_lorentzian_of_scaled_oriented_poles` composes the result
  through the Phase CO basis invariance. Exact rational gap search constructs
  a valid certificate on the first candidate pair for all 2,868 nontrivial
  Phase CN cutoffs, including cases needing negative orientation and negative
  raw chart scale. Proving the corresponding arbitrary-pivot conditioned-
  Taylor minor signs and separator existence, then the infinite-tail passage,
  higher fresh order, normalized tangent envelope, order-eighteen contraction,
  dyadic recurrence, Conjecture G1, and Problem 4.2 remains open. The root
  command passes a 3,595-job Lean build and 119 tests.
- Post-v1 research phase CR: complete; merged as PR #98 with green CI.
  `selectedCore_sq_mul_conditionedRankTwoMatchedCoordinate` proves the generic
  selected-core-square identity transporting each conditioned rank-two matched
  coordinate to the product of its two original feature-augmented core
  determinants.
  `orientedConditionedRankTwoMatchedCoordinate_nonneg_iff` gives the resulting
  two-way scalar-oriented sign transport. Exact rational reconstruction checks
  both augmented cores and the identity for all 11,024 feature pairs in the
  500 arbitrary-pivot cases. This is a representation bridge, not the still-
  open augmented Taylor minor-sign or cyclic-separator theorem. The root
  command passes a 3,596-job Lean build and 119 tests.
- Post-v1 research phase CS: complete; merged as PR #99 with green CI.
  `mul_nonneg_of_common_sign_schedule` proves that a shared nonzero schedule
  cancels after two independent reorientations, and
  `orientedConditionedRankTwoMatchedCoordinate_nonneg_of_common_sign_schedule`
  composes the fact through the Phase CR bridge. Exact reconstruction confirms
  shared nonzero schedules for all 11,024 left/right augmented pairs, with both
  signs occurring, but refutes the stronger conjecture that the raw schedule is
  fixed by pivot depth and the degree pair. Proving arbitrary-run chirotope
  synchronization remains open. The root command passes a 3,596-job Lean build
  and 119 tests.
- Post-v1 research phase CT: complete; merged as PR #100 with green CI. The
  delivered generic nonvanishing continuous-path endpoint sign theorem and
  its determinant specialization lock endpoint signs whenever the path stays
  nonzero. Exact canonical normalization and midpoint checks pass for all
  11,024 augmented determinant pairs. This isolates but does not prove
  nonvanishing of the confluent exponential determinant throughout the
  ordered positive unit chamber. The root command passes a 3,597-job Lean
  build and 119 tests.
- Post-v1 research phase CU: complete; merged as PR #101 with green CI. The
  exact exponential-cubic endpoint determinants reduce, after positive column
  scaling, to `83/194059800` and `-1/194059800`. Lean proves that the straight
  path remains in the ordered open unit chamber and contains a determinant
  zero, refuting universal chamber nonvanishing. The obstruction deliberately
  leaves open the narrower paired-subtuple geometry of the actual GECP
  endpoint determinants. The root command passes a 3,598-job Lean build and
  120 tests.
- Post-v1 research phase CV: complete; merged as PR #102 with green CI. The
  exact rational exponential kernel `K(i,j)=(101/100)^((i+1)(j+21))` has
  certified complete pivots `(3,3)` and then `(0,0)`, but its crossed
  degree-`(0,3)` augmented determinants have opposite signs. This refutes
  arbitrary-run common-
  chirotope synchronization, including the paired-subtuple strengthening,
  without claiming failure of every aggregate or absolute determinant route.
  The root command passes a 3,600-job Lean build and 121 tests.
- Post-v1 research phase CW: complete; merged as PR #103 with green CI. The
  division-free
  theorem
  `four_mul_weightedCovarianceNumerator_mul_le_crossIncrement_sq_of_meanShiftProduct_nonneg`
  proves the fresh-order-two Lorentzian inequality from nonnegative
  within-block covariances and one aggregate product of block-mean shifts.
  `rankTwoBlockPencil_lorentzian_of_affine_meanShift` lifts it to the matrix
  pencil. Exact validation covers 42 Phase CV witness pencils and 16,424
  broader conditioned pencils, without claiming the infinite-tail mean-shift
  sign or the final companion envelope. The root command passes a 3,600-job
  Lean build and 122 tests.
- Post-v1 research phase CX: complete; merged as PR #104 with green CI. The exact exponential
  kernel `K(i,j)=(21/20)^((i+1)(j+1))` is realized on open-unit nodes and has
  certified first complete pivot `(3,3)`. Its four-degree normalized
  conditioned matched-coordinate Hessian has an explicit positive-definite
  two-plane, so global one-positive-direction or Hodge--Riemann signature is
  false. All three nontrivial contiguous prefix/suffix planes nevertheless
  satisfy the required reverse Cauchy--Schwarz inequality. This narrows the
  infinite-tail target to a prefix-cone variation-diminishing/interlacing
  theorem rather than global polynomial stability. The root command passes a
  3,601-job Lean build and 123 tests.
- Post-v1 research phase CY: complete; merged as PR #105 with green CI. The delivered
  `quadraticTransferNodeNumerator_invariant` identifies a shared negative-node
  value for adjacent degree cutoffs, and
  `quadraticTransferNodeNumerator_outerFeature_eq` factors its numerator into
  two cumulative projected cross determinants.
  `finTwoPencil_transfer_lorentzian_of_projectedCross` proves both adjacent
  matrix discriminant inequalities under explicit endpoint, transfer, and
  projected-cross sign hypotheses. Exact finite-Taylor checks cover 13,924
  transfers and 9,170 strict interior signs, with rank-one boundary failures
  retained. The arbitrary-run analytic sign, infinite-tail passage,
  higher-fresh-order envelope, dyadic recurrence, G1, and Problem 4.2 remain
  open.
  Local root verification passes a 3,602-job Lean build and 123 tests.
- Post-v1 research phase CZ: complete; merged as PR #106 with green CI. The four-node nine-term Taylor
  witness has certified endpoint complete pivots and a strictly positive
  canonical transfer-node numerator for moving degree three. Both adjacent
  pencils retain positive endpoint determinants and positive discriminants.
  The degree-two and degree-four nodes supply negative alternate
  certificates. This refutes the universal strict-interior sign premise
  without refuting adaptive coverage or the full exponential target.
  `experiments/transfer_nodes.py` provides deterministic exact audits over all
  remaining fresh pairs, recording fixed-node failures and finite adaptive
  coverage. The bounded 100-case audit checks 54,845 transfers and covers
  47,496 interior cutoffs with four failed fixed nodes.
  `phaseCZ_interior_transferNode_pos`,
  `phaseCZ_alternate_transferNodes_neg`, and
  `phaseCZ_neighbor_discriminants_pos` are re-exported and axiom-audited.
  Root verification passes a 3,603-job Lean build and 126 tests.
- Post-v1 research phase DA: complete; merged as PR #107 with green CI. The analytic one-pivot base case
  is delivered by `expMonomialMinor_pos_of_pos_degrees`,
  `expMonomialMinor_zero_one_pos`, `expMonomialMinor_neg_of_zero_left`,
  `expMonomialMinor_signed_pos`, and `expMonomialMinor_matched_mul_pos`.
  These theorems concern actual exponential/monomial augmented minors on the
  closed unit interval, at positive rates at most one. They do not establish
  arbitrary-prefix adaptive coverage, the higher-order companion envelope,
  dyadic recurrence, G1, or Problem 4.2. A separate exact conditioned-Taylor
  regression rejects localization of bad nodes at `pivot_depth+1`.
  Local root verification passes a 3,605-job Lean build and 129 tests, with
  only the accepted public theorem axioms. Push and PR CI both pass.

- Post-v1 research phase DB: complete; merged as PR #108 with green CI. `phaseDB_not_adaptive_coverage`
  refutes neighboring canonical-node coverage on a coherent thirteen-term
  Taylor kernel after seven finite-grid complete pivots. Both adjacent
  numerators and all transfer sums are positive; the middle pencil retains a
  positive discriminant. `phaseDB_coherent_Taylor_kernel`,
  `phaseDB_residual_recurrence`, `phaseDB_complete_pivots`,
  `phaseDB_conditioned_features_reconstruct`,
  `phaseDB_actual_endpoints_and_transfers_pos`, and
  `phaseDB_actual_cutoff_discriminant_pos` expose the checked construction.
  The integer data is independently reproducible by
  `scripts/generate_coverage_certificate.py --check`. This does not refute
  the full exponential Lorentzian inequality or prove any continuous-domain
  complete-pivot claim, higher-order companion envelope, dyadic recurrence,
  G1, or Problem 4.2.
  Local root verification passes a 3,608-job Lean build and 133 tests with
  only the accepted public theorem axioms. Push and PR CI both pass.

- Phase DC status annotation (2026-10-19): an explicitly weaker all-orders
  one-step bridge is verified and merged as PR #109 with green CI, not the
  original arbitrary-prefix companion
  target. `Matrix.det_le_diagonal_mul_principalComplement_of_strictlyTotallyPositive`
  proves nonsymmetric singleton Fischer;
  `Matrix.abs_det_pivotSchur_le_matchedComplement_of_columnMax` bounds the
  arbitrary-pivot Schur minor by the original minor deleting the row with the
  pivot column's sorted rank. `GECP.localMinor_schur_le_of_strictTotalPositive`
  and `Fermionic.positiveSeriesFermionicKernel_localMinor_schur_le` preserve
  the same local columns at every order. Exact tests cover 2,520 one-step
  sample checks and 1,376 prefix/gauge checks; only the former bridge is
  formally proved. `experiments/companion_envelope.py` provides reproducible
  actual-exponential stress evidence through depth 64 with independent
  320/640-digit evaluations. These observations do not prove the companion
  envelope. Common all-orders residual sign-gauge transport, initial dyadic
  local budgets, arbitrary-prefix contraction, G1, and Problem 4.2 remain
  open; no original mathematical objective is replaced by this bridge.
  The root gate passes a 3,610-job Lean build, accepted axiom audit, and
  140 tests with no known dependency vulnerabilities.

- Phase DD status annotation (2026-10-19): verified and merged as PR #110
  with green push and PR CI.
  `GECP.Run.gaugedFinalResidual_strictSignRegular` proves that one
  magnitude-preserving row/column sign gauge makes every ordered fresh minor
  positive through any successful real-coordinate positive-minor kernel run.
  `abs_gaugedFinalResidual_minor_det` preserves every absolute minor;
  `finalResidual_localMinor_schur_le` discharges the arbitrary-prefix structural
  premise of the Phase DC matched-rank bound for the actual next residual.
  `Matrix.det_vandermonde_append` provides the underlying exact block identity
  over any commutative ring. `Run.negateColumns`,
  `fermionicKernel_reflectedRun_gaugedFinalResidual_strictSignRegular`, and
  `fermionicKernel_reflectedRun_finalResidual` transport original fermionic
  runs without reselecting pivots; signed frequencies and boundary times are
  included. Exact regressions check 3,825 common-gauge minors and 350 actual
  next-step bounds, with separate 160/320-digit actual-exponential evidence.
  This closes the common all-orders gauge claim, not the initial local dyadic
  determinant budget, iterated budget/pivot counting, G1, or Problem 4.2.
  Local root verification passes a 3,614-job Lean build and 146 tests, with
  only the accepted public axioms and no known dependency vulnerabilities.

- Phase DE status annotation (2026-10-19): verified; merged as PR #111 with
  green push and PR CI. `positivePanelTerms_error` and
  `fermionicPanelTerms_error` give `8p` separated terms with actual-kernel
  error `2^(-p)` on each signed central/dyadic panel, independently of its
  scale. `fermionicKernel_panel_sample_det_le_two_pow` supplies the initial
  rank-sensitive determinant estimate, and
  `fermionicKernel_panel_minor_det_le_dyadicBudget` proves the explicit local
  envelope `(2^(-q))^n` at `n=16*(33+4q)` for all natural panel/accuracy
  indices. Physical time endpoints, both frequency signs, repeated tuples,
  and arbitrary sample order are included. Exact feature/exponent arithmetic
  and independent actual-exponential precision checks supplement the formal
  proof; no target-order numerical determinant is treated as certification.
  This closes the initial local budget, not its backward transport through
  prefixes, panel-pivot counting, the global rate, G1, or Problem 4.2.
  Local root verification passes a 3,615-job Lean build and 151 tests with
  only the accepted public axioms and no known dependency vulnerabilities.

- Phase DF status annotation (2026-10-19): root verification passes 3,619
  Lean jobs, accepted public axioms, and 155 tests, with no known dependency
  vulnerabilities; verified and merged as PR #112 with green push and PR CI.
  `GECP.markedPivotProduct_le_initialMinor`
  constructs actual initial-kernel minors dominating selected pivot products;
  `panel_pivot_count_lt_of_initialBudget` combines exactly the budget order
  of marks with a fixed-order initial envelope. `fermionicPanelDomain_cover`,
  `fermionicPanelDomain_unique`, and `fermionicPanel_count` prove the signed
  finite partition. `fermionicKernel_gecp_error_le_localPanelBlock` proves
  error at most `2^(-q)` after `32*(s+1)*(33+4*q)` exact complete pivots on
  `[0,1] x [-2^s,2^s]`; `fermionicKernel_gecp_error_le_of_localPanelBlock`
  extends the guarantee to every later rank. The actual rate assumes only
  initialization, exact successful updates, realized prefixes, and complete
  maxima, not a residual determinant/convergence certificate. Independent
  finite exact and high-precision checks supplement the proof. This closes
  the dyadic finite global rate, not yet the arbitrary-real-cutoff/accuracy
  corollaries or the final full-SPEC / G1 / Problem 4.2 completion audit.

- Phase DG status annotation (2026-10-20): root verification passes 3,622
  Lean jobs, the permitted public axioms, 159 tests, and no known dependency
  vulnerabilities; verified and merged as PR #113 with green push and PR CI.
  The actual-domain
  `fermionicKernel_gecp_error_le_realCutoffBlock` retains the finite DF count
  for any `Lambda<=2^s` without padding the complete-pivot domain.
  `fermionicKernel_gecp_error_le_logarithmic` proves the G1-shaped bound for
  `Lambda>=1` with universal `C=1024` and `c=(log 2)^2/256`.
  `fermionicKernel_gecp_error_le_allCutoff_logarithmic` proves a
  `log(2+Lambda)` version for all positive cutoffs, and
  `fermionicKernel_gecp_error_le_realTolerance` gives an explicit ceiling rank
  for every positive requested tolerance. The separate realized-run premise
  is derived from initialization and exact successful updates.
  `fermionicKernel_smallCutoff_logarithmic_obstruction` proves that literal
  G1 cannot hold uniformly for all `Lambda>0`: for any proposed `C,c>0`,
  `Lambda=min(1,c^2/(16*C))` and the actual first complete pivot contradict
  that bound at rank one. This annotation preserves the original conjecture
  and records its parameter boundary rather than rewriting its goals.
  The final audit must distinguish the proved large-cutoff theorem, this
  rigorous small-cutoff obstruction, and the remaining finite-prefix,
  termination, cross-approximation, and full-SPEC completion obligations.

- Phase DH status annotation (2026-10-20): root verification passes 3,625
  Lean jobs, the permitted public axioms, 163 tests, and no known dependency
  vulnerabilities; verified and merged as PR #114 with green push and PR CI.
  `fermionicKernel_finiteRun_error_le_exp` proves the explicit large-cutoff
  logarithmic rate directly for any finite actual `Run fermionicKernel`
  satisfying recursive `CompleteOn` on its rectangle. No successful or
  complete pivot after the stored run is required. Finite dyadic-block,
  all-positive-cutoff and real-tolerance versions retain DG's constants.
  `residual_eq_zero_of_complete_zero` and the stopped-or-budget corollary
  include exact termination without a zero-pivot update.
  `Run.finalResidual_eq_sub_crossApproximation` identifies the residual with
  the original selected-core cross error, and
  `fermionicKernel_finiteRun_greenError_le_exp` bounds the actual difference
  of weighted Green integrals, deriving integrability from continuity and
  the uniform bounds. Its `_interval` corollary discharges the support premise
  on the physical interval. Independent checks cover finite prefixes, exact stops,
  selected-core solves, signed Green sums and independent precision.
  This closes the future-success, termination and physical-error interface
  gap, not the separate full-SPEC audit, approximate-pivot robustness, sharp
  constants, or the literal small-cutoff G1 obstruction already proved in DG.

- Phase DI status annotation (2026-10-20): original sections 1–30 are audited
  in `COMPLETION_AUDIT.md`, with 33 supported and seven incomplete requirement
  groups after the bounded closures. The full original goal is not complete.
  `fermionicKernel_separatedApprox_real` supplies D1's real cutoff/tolerance
  quantifiers with explicit rank `128*log(1+Lambda)*log(128/epsilon)`.
  `fermionicKernel_finiteRun_exists_completePivot` derives a compact maximum
  for every finite residual; concrete exact small-grid completeness/minor/
  pivot certificates are exported. Axiom gates now enforce the allowed set,
  report coverage and failing Lean exit status. Packaging excludes local
  duplicate notes and has an isolated wheel smoke check. Root verification
  passes 3,627 Lean jobs, 838 enforced axiom reports, 177 tests, and no known
  dependency vulnerabilities; verified and merged as PR #115 with green push
  and PR CI. Required work includes the continuous SPD/original-
  kernel Lipschitz baseline (not just finite canonical traces or an assumed
  residual Lipschitz bound), full result metadata, `1e-12` census depth and
  geometry/normalized plots, and listed Phase-II method comparisons. These
  gaps are recorded rather than removed from the mathematical goals above.

- **Phase DJ status (2026-10-20):**
  `KernelPosSemidef` lifts finite sampled Gram positivity to arbitrary kernel
  domains (including repeated samples). An actual finite diagonal run derives
  PSD of the residual and removed kernel, original-energy comparison, and
  `diagonalRun_error_le_original_fillDistance`: `abs R(x,y) <= 2*L*h` from
  original-kernel coordinatewise Lipschitz regularity and selected-node cover,
  without assuming residual regularity. Compact continuous domains admit a
  diagonal complete maximum. This closes audit R09; strict-SPD necessity for
  every complete pivot through fresh residuals remains open in R08. No weak-
  PSD tie necessity, asymptotic P-greedy rate, or new numerical certificate is
  claimed. Root verification passes 3,630 Lean jobs, 857 enforced axiom reports,
  180 tests and no known dependency vulnerabilities; green push and PR CI;
  merged as PR #116. The original mathematical goals above are unchanged.

- **Phase DK status (2026-10-20):** `completeRun_diagonal_of_posDef`
  derives diagonal necessity for every successful exact complete pivot on an
  arbitrary domain from strict SPD of the original kernel. Fresh strictness
  is propagated through injective rectangular Schur congruence, not assumed.
  `posDef_completeRun_iff_pivotedCholesky` proves equivalence with independently
  defined recursive diagonal maximization on actual residual updates.
  Residual PSD/diagonal domination, fresh strictness, compact continuous
  maximum existence and original-kernel fill control are exported. This
  closes R08's staged baseline scope; full published asymptotic rates, weak-
  PSD diagonal necessity, unique ties and numerical exactness are not claimed.
  The original goals above are unchanged; numerical metadata/census and
  Phase-II comparisons remain required work. Root verification passes 3,633
  Lean jobs, 876 enforced axiom reports, 183 tests and no known dependency
  vulnerabilities; green push and PR CI; merged as PR #117.

- **Phase DL status (2026-10-22):** typed `NumericalMetadata`,
  `ConditioningDiagnostic` and `RefinementDiagnostic` expose actual/requested
  precision, conditioning status/history, tie policy/counts/tolerance and
  certification/error scopes for Cholesky/sparse results, a GECP evidence view,
  and explicit non-iterative approximation scopes. Schema-2 synthetic summaries
  preserve these records; historical schema-1 committed artifacts are unchanged.
  Convergence reflects the returned reconstruction/training error, candidate
  exhaustion cannot reuse atoms, and optimizer failure/budget/rejection is
  separate from fit convergence. R24 is supported. This adds no Lean theorem,
  arbitrary-precision Cholesky/sparse implementation or continuous numerical
  certificate. Census and comparison requirements above remain unchanged.
  Root verification passes 3,633 Lean jobs, 876 enforced axiom reports,
  232 tests and a clean dependency audit after targeted transitive lock fixes;
  isolated wheel smoke and green push/PR CI; merged as PR #118.

- **Phase DM status (2026-10-22, verified; merged PR #119):** schema-2 `CensusRecord`
  implementation adds complete residual/determinant histories and boundary,
  dyadic-panel and nearest-selected-reflection diagnostics. The canonical
  configuration now requests four tolerances through `1e-12`, and normalized
  residual plotting plus independent higher-precision selected-core/held-out
  checks are implemented. The full repeated 28-case primary artifact and both
  plots are present; three actual order-32 GECP cases also converge and match
  the corresponding order-24 ranks, with different paths. Thirty independent
  192-bit core checks pass, with an explicit non-certification scope. Full
  root verification passes 3,633 Lean jobs, 876 axiom reports and 430 tests;
  exact-head push/PR CI pass and PR #119 is merged. This supports R25–R27;
  original mathematical objectives above are unchanged. R31 comparisons and
  the final whole-SPEC audit remain required. Pushed-main closure gates repeat
  before final handoff; no continuous numerical certificate is inferred.

- **Phase DN status (2026-10-22, verified; merged PR #120):** implements the six listed
  numerical comparisons in §§19–20 using shared time observations, explicit
  Matsubara/Prony bridges, original-kernel weight solves and retained failure
  diagnostics. The contract is frozen in RESEARCH; all 192 canonical fits and
  the reviewed PNG repeat bytes. Full root verification passes 513 tests,
  3,633 Lean jobs and 876 axiom reports; exact-head push/PR CI and installed-wheel
  checks (including SciPy 1.15) pass. PR #120 is merged; exact pushed-main
  closure gates repeat before handoff. Original mathematical goals above are
  unchanged. No globally optimal sparse support or unique recovery
  is inferred from a smallest tested passing support.

- **Phase DO audit status (2026-10-23, verified; merged PR #121):** the independent
  `FINAL_AUDIT.md` checks all thirty original sections against actual theorem
  assumptions, sources, artifacts and tests. It supports the original-kernel
  exact finite-run rate, separated approximation, strict-SPD/fill baseline and
  cross/Green consequences, but reopens R25: fixed order-24 seven-cutoff data
  and three fixed order-32 reruns do not implement §14's increasingly dense
  residual-adaptive grids. DO repairs interval lower-bound rounding, sampled
  adaptive tie metadata, a large-finite float64 update failure and missing
  independent centered-expression coverage, plus generic import layering and
  stale current documentation. Final `f0d082d` passes the 569-test root,
  3,634 Lean jobs, 876 axiom reports, exact-head push/PR CI and installed-wheel
  checks with pinned SciPy 1.18 and minimum 1.15.0. PR #121 is merged;
  exact pushed-main closure gates repeat before handoff.
  The next required target is the genuine seven-cutoff adaptive experiment and
  subsequent final closure. Original goals above are unchanged; the prior
  forty supported labels are not accepted as a completion proof. No ideal-
  exact-run or separately reconstructed-cross certificate, arbitrary-η rate,
  sparse optimality or literature priority is inferred from numerical evidence.

- **Phase DP status (2026-10-23, verified; merged PR #122):** implements
  `kernelgecp.adaptive_grids` and `adaptive_grid.run_adaptive_grid`, explicit
  nested physical axes and residual/surplus-driven selective refinement for
  all seven §14 cutoffs and four tolerances through `1e-12`. Both clean-source
  executions at `3cee4186fddf990abdd33dbe24c46ff908b1923f` converge in all seven
  cases/28 prefixes and repeat `experiments/data/adaptive_census.json` bytes.
  Actual stage counts are `3,2,3,4,4,5,6`; strict chosen ranks are
  `9,17,36,64,97,124,151`. Untouched seventh-fraction samples are evaluated
  only after grid/rank decisions. The normalized/stage PNG is reviewed and
  byte-repeated; independent `experiments.check_adaptive_census` MP192 checks
  pass all 21 selected cores and 28 final crosses. The full root and package
  gates pass locally: 713 tests, 3,634 Lean jobs, 876 allowed axiom reports and
  inspected archives/installed-wheel smokes with pinned/minimum SciPy. Final
  artifact-head push/PR CI `40017052450`/`40017056920` pass at `a6488cc`;
  final `3c38183` also passes the root and CI `40017937659`/`40017945103`.
  PR #122 merges as `3a6c01585d2a1407b13e9ff55378f7d3a1f2fc8f` at
  2026-10-23 20:20:48 UTC; exact pushed-main closure repeats before handoff.
  Elimination is MP128,
  coordinates retain 53-bit repr-decimal semantics, and every maximum is
  finite-sample evidence, not a continuous certificate. Original mathematical
  goals are unchanged. Final whole-scope closure follows this numerical phase.

- **Phase DQ source-audit status (2026-10-24 UTC; verified, merged PR #123):**
  rereads original §§1–30 and checks every explicit outcome, diagnostic, test,
  comparison method and proposed import against actual source and theorem
  dependencies in `FINAL_AUDIT.md`. `scripts/CompletionContract.lean` checks
  22 expected public theorem interfaces without introducing a new theorem.
  Fresh MP192 adaptive checking and adaptive/spectral PNG byte comparisons
  pass. The audit reproduces and repairs a one-ULP conditional Lipschitz pivot
  upper-bound error and rounded-log Taylor parameters that miss their actual
  cutoff/tolerance; exact-ratio boundary tests include subnormal inputs.
  Standalone certificate metadata distinguishes point/enclosure/export
  precision and retains uncertified coordinate-resolution outcomes. A live
  dependency audit also requires a development-only filelock 4.2.1 update;
  historical data/lock manifests are preserved. The original required scope
  is accounted for at its explicit theorem/conditional/experimental boundaries,
  but full root, archive/install, exact-head CI, PR merge and pushed-main
  closure must pass before overall completion. Original mathematical goals
  above are unchanged. No arbitrary-inexact-pivot rate, ideal-exact-run
  numerical certificate, global sparse optimum or priority is inferred.

  DQ code head `f2ffd5e344818e2c06ff5152a7f182ef7a31e3a8` subsequently passes
  the exact full root: 757 tests in 150.01 seconds, 3,634 Lean jobs, 876
  permitted axiom reports, 22 public-import examples and clean style/type/
  security gates. Push/PR CI `40050118292`/`40050129457` pass all four checks.
  Clean rebuilt archives and both outside-checkout installed environments
  pass, including 45 exports, adaptive success/budget failure and 48 six-method
  fits per environment. Draft PR #123 is open. Final documentation-head and
  pushed-main closure gates remain; these checks do not alter the goals above.

  Final documentation branch `0a6d3635529279beac9fb198c2d9de147c5ee2e5`
  repeats the exact root (757 tests in 149.97 seconds), rebuilt archives and
  both installed environments. Push/PR CI `40050605960`/`40050610230` pass all
  four checks. PR #123 merges at 2026-10-24 01:41:41 UTC as
  `d9978607b88a45d98d917f44318bcaf99b4c3239`; local main fast-forwards and
  immutable generation `3cee4186` remains an ancestor. Exact pushed-main
  root/CI/archive/install and synchronization checks precede final handoff
  and tracker completion. Required original-SPEC scope is complete at those
  explicit boundaries; the 100% estimate is not a positive answer to every
  research question, a continuous numerical certificate, or an external
  priority claim. The original mathematical goals above are unchanged.

- **Phase DR documentation status (2026-10-23 EDT; verified, merged PR #124):** finalizes
  the documentation under §29 without changing any original mathematical goal,
  theorem, numerical algorithm, historical data or publication state. DQ's
  actual exact-main closure is `e8c35eda7662566639c5eac6748d3032df87075f`:
  757 tests in 150.31 seconds, 3,634 Lean jobs, 876 allowed axiom reports,
  22 public-import examples, clean style/type/security and package/install
  checks, and successful CI `40051138962`. The original goal was marked
  complete afterward, with terminal cumulative accounting of 24,738,403
  aggregate tokens and 116,095 seconds (32:14:55), excluding this subsequent
  documentation task. Canonical exposed/unavailable runtime fields remain in
  [FINAL_HANDOFF.md](FINAL_HANDOFF.md#implementation-run-metadata). Dated
  historical pending statements retain their original phase context; current
  summaries record the completed outcome. DR repeats root/CI/package and
  merge/push verification before its own final handoff.

  Exact DR documentation source `51fa6b3` passed 759 tests in 144.90 seconds,
  unchanged Lean/axiom/consumer and style/type/security gates, CI
  `40054146369`/`40054160333`, clean archives and both installed environments.
  PR #124 merged on 2026-10-23 at 22:30:40 EDT as
  `a421f7c8ed8b1b8754f38116939e26181636fe82`. The original mathematical text,
  proofs, algorithms and historical data/manifests remain unchanged. DR is
  documentation/hygiene delivery, not a new mathematical or publication claim.

- **Phase DS workflow status (2026-10-24 EDT; verification/delivery in
  [PR #125](https://github.com/tripp-smith/gecp-kernel-structure/pull/125)):** updates the two
  workspace skills under §29 with original-obligation contracts, proof-route
  checkpoints and conditional numerical/verification/delivery guidance. This
  post-completion maintenance does not alter the original mathematical goals,
  completed scope decision, formal proofs, numerical code, dependency lock or
  historical experiment provenance. The phase requires skill validation,
  realistic behavioral checks, the root command and exact-revision PR CI before
  authorized merge/push; detailed contract and behavioral evidence are in
  `RESEARCH.md`, with the revision-specific check and delivery receipt in the PR.

- **Phase DT application extension (2026-10-24 EDT; verification/delivery in
  [PR #126](https://github.com/tripp-smith/gecp-kernel-structure/pull/126)):** establishes
  matched physical forward Green-function benchmarks against optional `pydlr`
  and `sparse-ir` under §§27–29. The existing mathematical goals and completion
  decision are unchanged. Identical observations, a common cached SVD fitter,
  untouched sampled validation, explicit failure records and separate setup/
  reuse/process-memory measurements precede any Rust optimization. This is a
  pilot application extension, not a new theorem, performance claim, certified
  numerical run, public API redesign or package publication. Contract and
  acceptance tests are recorded in `RESEARCH.md`. The frozen nine-case pilot
  retains eight successful sampled cases and a float64 GECP arithmetic-range
  failure at cutoff 1,000. All 24 successful forward workloads pass independent
  100-digit checking; two runs match numerical fields, not machine-dependent
  measurements. Existing numerical behavior is not changed to hide the failure.

- **Phase DU application extension (2026-10-24 EDT; final verification/delivery
  in [PR #127](https://github.com/tripp-smith/gecp-kernel-structure/pull/127)):**
  adds installed `FermionicBasis`, `FermionicFit`, `LeastSquaresPlan`, immutable
  diagnostics, cached solve/SVD designs, arbitrary complex batch axes, explicit
  physical tau/ordinary-index Matsubara conventions, versioned pickle-free
  persistence, typing and tutorials under §§27–29. The original mathematical
  goals and completion decision remain unchanged. Explicit float64/reference/
  automatic policies disclose selection, application and export precision;
  actual reference kernel evaluations/core solves recover the cutoff-1,000
  smoke case without changing the research `gecp` defaults or hiding its
  original failure. Sampled tests are not continuous or rounding certificates.
  Wheel/sdist rebuild and outside-checkout install/typing CI are introduced for
  three platforms and two Python versions; their actual results remain a
  delivery gate. `PRODUCTIZATION.md` tracks the broader benchmark, solver,
  Rust, numerical-bridge, release/publication and adoption requirements. A
  pending PyPI publisher is configured after user confirmation, but no package
  upload or completed full application goal is claimed.

- **Phase DV application extension (final verification/delivery in
  [PR #128](https://github.com/tripp-smith/gecp-kernel-structure/pull/128)):** extends
  §§23–29 with physical Matsubara adapters, matched-accuracy calibration,
  separately measured fresh-worker construction/reuse, complex matrix workloads
  and a small noninteracting Bethe self-consistency loop. The failed low-frequency
  sampling pilot is retained; independent physical accuracy is distinct from
  observed convergence. Fresh validation is fixed after pilot decisions and
  before canonical generation. Original mathematical goals, proved theorems and
  historical evidence are unchanged. `PRODUCTIZATION.md` keeps Rust evaluation,
  the numerical error bridge, release/publication and upstream adoption open.

  Frozen cohort A preserves an independently confirmed calibration miss;
  cohort B strengthens calibration, reserves fresh validation and passes all
  12 method/case records and three physical solver loops. Two runs have identical
  numerical fields; independent 95-digit pole/root/integral checks confirm the
  reported outcomes. Constructor/reuse/RSS scopes and different native eps/ranks
  remain explicit. A three-provider example and small standalone Sparse-IR recipe
  execute actual installed libraries. The concrete upstream proposal has not
  been submitted or accepted. No new theorem, Rust speedup or publication is
  claimed by this phase.
- **Phase DW application follow-on (final acceptance in
  [PR #129](https://github.com/tripp-smith/gecp-kernel-structure/pull/129)):** adds actual finite-query
  coefficient/core/solve/output discrepancy budgets and true selected-cross
  pivot quality on an explicit finite grid. Public conditional theorems are
  `Discretization.approxPivot_of_uniform_discrepancy` and
  `Discretization.dot_error_le_of_component_bounds`. The original continuous
  exact-GECP goals and delivered claims are unchanged. Numerical enclosure of
  exact physical binary64-input problems is distinguished from physical truth,
  truncated-SVD targets and a continuous inexact-trajectory rate. The six-case
  frozen protocol preserves an interval-resolution failure and a valid but
  large-discrepancy ill-conditioned fit, alongside useful observation designs.
  Two complete frozen artifacts are byte-identical; all five completed
  coefficient/prediction enclosures pass independent QR/exact-reproduction
  checks. The inconclusive case and physical/model errors remain explicit.
- **Phase DX application follow-on (delivered in
  [PR #130](https://github.com/tripp-smith/gecp-kernel-structure/pull/130)):** implements an explicit
  optional 128-bit Rust time-evaluation backend, preserving research defaults,
  portable archives, mpmath selection/frequency and NumPy fitting. Native input,
  arithmetic-state and resource boundaries are tested. Universal/source and
  native distributions have separate archive/license/source-identity gates.
  The frozen six-process production cohort passes the predeclared retention
  gate with 50.15% median application and 47.76% process-wall reduction for
  changing time meshes, physical errors below `2.57e-11` and paired differences
  below `4.26e-15`. Exact-main root, source and pure/native platform/install
  acceptance pass at `67e4b14`; the original mathematical goals and exact-GECP
  claims are unchanged.
- **Phase DY application release follow-on (delivered in
  [PR #131](https://github.com/tripp-smith/gecp-kernel-structure/pull/131)):** publishes
  [kernelgecp 1.1.0 on PyPI](https://pypi.org/project/kernelgecp/1.1.0/) and the
  [GitHub release](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0)
  from immutable tag/source `9a679930db326bda2456bfa91b3256f37359e4ce`.
  Exact-main root/CI/producer acceptance precedes promotion of the original
  three tested wheels and sdist through the authorized Trusted Publisher.
  Complete corresponding native source, notices and dependency inventory
  were public before upload. Workflow `40244470229` attempt 2 passes all six
  clean public OS/Python wheel installs and separate public sdist installation;
  public file hashes equal the tested manifest. Actual upstream outreach is
  [publicly recorded](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382);
  acceptance remains pending. The
  [release receipt](RELEASE_RECEIPT.md) preserves exact source/file provenance
  and the first-attempt metadata-visibility failure; only failed verification
  jobs were rerun, without another upload. Controlled functional application
  requirements are delivered at their stated numerical/platform boundaries.
  No original mathematical goal is rewritten or weakened.
- **Phase DZ published-release handoff (final acceptance in
  [PR #132](https://github.com/tripp-smith/gecp-kernel-structure/pull/132)):** synchronizes
  current documentation with actual public release evidence and extends only
  the future verifier's bounded metadata-readiness waits. Released tag/assets,
  original mathematical goals, proved claims and frozen numerical evidence
  remain unchanged. Final exact-head root/CI and merged-main delivery evidence
  belong to the external PR receipt rather than a self-referential source hash.
- **Phase EA invited upstream contribution (final acceptance in
  [PR #133](https://github.com/tripp-smith/gecp-kernel-structure/pull/133)):** delivers
  the maintainer-invited [Sparse-IR notebook PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10),
  preserving existing cells/metadata and adding normalized complex-matrix
  physical checks without a kernelgecp dependency. Isolated/full execution and
  the documented full book build pass under the upstream lock. Submission-time
  fork CI requires maintainer approval and is not claimed successful; current
  review/merge status is linked at the upstream PR. Published 1.1.0 assets,
  original mathematical goals and frozen scientific evidence are unchanged.

- **Accepted production integration extension EB–EJ (implementation requested):**
  the [five-field contract](PRODUCTION_INTEGRATIONS.md) adds supported optional
  TRIQS/CT-HYB, molecular Green GF2, EPW, MaxEnt and SciPy/GPyTorch workflows,
  frequency fitting and archive compatibility, realistic provider evidence,
  and a separately verified 1.2.0 release. The current EB contract is in
  RESEARCH.md and the acceptance ledger extends PRODUCTIZATION.md. All domain
  gates remain pending until their actual computations pass. This status
  annotation does not rewrite the original mathematical goals, strengthen
  certificates, or replace frozen 1.1.0 evidence.
- **Phase EB delivered identifiers (final acceptance in PR #134):** isolated
  provider specifications/build recipes, actual numerical probes under
  scripts/providers, fail-closed version/source/domain-evidence validation,
  export_conda_lock, platform-specific actual-package SHA256 locks/inventories,
  observed pinned source notices and the production-provider CI workflow.
  Six actual platform/provider executions pass at `d5f58a1`; final locked
  revision/main acceptance is recorded in the PR. This closes no EC–EJ
  scientific application gate and changes no original theorem or 1.1.0 asset.

- **Phase EC implementation contract:** adds cached ordinary-index Matsubara
  fitting, observation-domain diagnostics, v2 archives with actual public-v1
  compatibility checks, and optional actual TRIQS 4.0.2 single-Gf/block recipes.
  Numerical fits preserve supplied observations without implicit symmetry and
  do not strengthen the exact-GECP theorem or time-fit certificate scope.
  Root/provider/exact-head/main delivery remains required after EB acceptance.
- **Phase ED implementation contract:** deterministic interacting half-filled
  Bethe IPT, independent atom/finite-bath/finite-pole references and dense
  native Fourier refinement. Public numerical solver interfaces, frozen
  common-trajectory replay and independent closed-loop comparisons retain
  model/representation/convergence/certification distinctions. Exact-source
  acceptance and delivery remain required; no original theorem goal changes.
- **Phase EC delivered acceptance:** source `d0b21ca` and merged main `64e9f34`
  pass root verification, actual installed-provider tests and all applicable
  exact-revision workflows; PR #135 owns canonical branch/main receipts.
  Frequency fits and v2 archives target 1.2.0; public 1.1.0 remains immutable.
- **Phase ED implemented identifiers:** `IPTConfig`, `PhysicalBasis`, `IPTPlan`,
  immutable update/result/diagnostic records, optional `DenseIPTPlan`, finite
  Fock/Lehmann/pole references and `experiments.interacting_dmft`. Independent
  full campaign, root verification and delivery remain pending.
- **Phase EG implementation contract:** source-backed anisotropic EPW output,
  explicit Nambu/frequency-fit/archive/reconstruction and native restart
  interchange, with physical error/conditioning gates and authentic coarse
  MgB2 evidence. Scientific/root/CI/merge/main acceptance remains pending;
  no EPW internal acceleration or original mathematical claim is added.


- **Phase EI implementation contract:** isolated scalable covariance work adds
  an immutable column-oracle KernelCholeskyFactor, out-of-sample triangular
  extension, SciPy operators/shifted solves and optional frozen CPU float64
  GPyTorch conversion while preserving dense pivoted_cholesky. The contract
  in RESEARCH.md retains all UCI real-data, accuracy, memory, provider and
  exact-revision delivery gates in PRODUCTION_INTEGRATIONS.md. Implementation
  and domain acceptance are in progress; no original mathematical goal or
  published 1.1.0 evidence is changed.
- **Phase EI implementation identifiers (acceptance pending):**
  `KernelCholeskyFactor`, immutable extension/solve records, shifted SciPy
  operators/preconditioners and optional CPU provider conversion;
  experiments/covariance_benchmark.py plus scripts/audit_covariance.py bind
  attributed real observations, calibration isolation, numerical artifacts
  and source Git blobs. All scientific commands use the same locked GP
  runtime. Full real-data, root/CI and delivery gates remain mandatory;
  no synthetic/pilot result closes this extension or changes an original goal.
- **Phase EI frozen numerical evidence:** generation `203e58a` passes the
  six-cohort, three-repeat attributed UCI campaign, actual Mac CPU provider
  comparison and independently recomputed Git/blob/array gates. Public APIs
  and the dense interface remain distinct. [COVARIANCE.md](COVARIANCE.md)
  records finite representation errors, measured resource/timing scopes and
  non-claims. Final native/root, exact-head CI and delivery are still required;
  no continuous certificate or original mathematical statement changes.

- **Phase EF implementation contract:** supported molecular Green GF2 HDF5
  import with bound manifests, immutable matrix axes/static terms and complete
  iteration validation; actual H2/STO-3G two-geometry/native-grid refinement,
  independent finite-Fock/PySCF references and same-solution compression checks.
  Preparation follows the bounded contract in RESEARCH.md; scientific/root/CI
  acceptance and delivery remain open. Original mathematical goals are unchanged.
- **Phase ED delivered acceptance:** merged main `7530ec7` passes full root
  verification, all four exact-main workflows and canonical native/scientific
  acceptance in [PR #136](https://github.com/tripp-smith/gecp-kernel-structure/pull/136#issuecomment-6445705017).
  Historical generation remains `a3a64ba`; public 1.1.0 assets and original
  theorem statements are unchanged.
- **Phase EF native pilot:** source `a31d048` passes four actual Green GF2
  calculations (two H2 geometries, two native grids) in run 40527992216.
  This is distinct from final held-out representation/model/audit and delivery
  acceptance, which remain open in draft PR #137.

- **Phase EF comparison status:** frozen `56c3211` passes the declared
  calibration and fresh sampled validation gates; the earlier failed control
  run remains recorded in RESEARCH.md. Independent audit passes; full-root, exact-source CI and delivery remain open.

- **Production delivery status update:** ED, EF and EI are delivered under
  PR #136, #137 and #139 with exact-main root and CI acceptance. These later
  receipts supersede their earlier pending-delivery status annotations.
- **Phase EG implemented identifiers:** `EPWAnisotropicData`,
  `CauchyFrequencyPlan`, `EPWFrequencyCrossFit`, `fit_epw_frequency_cross`,
  the frozen `epw_frequency` campaign, independent `audit_epw_frequency`,
  and native `epw_restart` comparator. The explicit real Nambu parity policy
  is separate from generic Matsubara fitting and the specialized EPW IR format.
  Actual native source fixtures, failed routes and calibration are retained;
  final scientific, root, exact-source CI and delivery gates remain open.

- **Phase EG scientific status update:** frozen `80ea4bd` passes final
  imaginary-axis Green/field accuracy and independent auditing, but actual
  native default-order Padé comparison fails unchanged real-axis budgets.
  The failure and a bounded native conditioning study are recorded in RESEARCH.
  Successful reads alone do not close the fitted-restart obligation. Root at
  `8b52a3d` passes 2,374 tests/35 skips; later diagnostic changes need their own
  root and exact-head CI. No mandatory gate or original goal is weakened.


- **Phase EG development status:** the fixed seven-order native study is
  inconclusive with all candidates failing; native baseline and compiled-helper
  rechecks pass. Explicit lossless archive, 20 K roundtrip/auditor and a
  truth-assisted sensitivity diagnostic are developmental. No native
  sensitivity result or fitted/scored 20/30 K fields exist; required pure
  fitted-restart, root, CI and delivery gates remain open. See
  [EPW evidence](EPW_INTEROPERABILITY.md#development-checkpoint-native-continuation-and-archival).

- **Phase EG 20 K checkpoint:** frozen `a9c2c39` producer cases pass fixed
  model holdout and exact archive replay checks. Their first independent audit
  fails on adjacent reconstructed binary64 coordinates. A verifier-only
  operation-order correction retains the independent MP oracle, unchanged
  conditioning cap and adjacent-point rejection. Corrected audit, unmodified
  native restart comparison and delivery remain required; no new mathematical
  goal or fitted-restart success is claimed.
