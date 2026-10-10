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
