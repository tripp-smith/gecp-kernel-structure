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
- Milestone D: verified on its phase branch. `expFamily_separatedApprox` proves an explicit
  `8p(s+1)`-term dyadic approximation of `exp(-tω)` on the positive-frequency
  domain, and `fermionicKernel_separatedApprox` proves a `16p(s+1)`-term
  approximation on `[0,1] × [-2ˢ,2ˢ]`, both with uniform error `2⁻ᵖ`.
  The construction is dyadic truncated Taylor rather than the published
  selected-exponential Chebyshev variant; it establishes the required
  logarithmic-logarithmic separated-rank scale without making a GECP claim.
- Milestone E: blocked (research). Exact sign evidence and a
  parameter-independent pivot-order obstruction are recorded, but no approved
  GECP rate/structural theorem or sufficient-condition obstruction is claimed.
- Milestone F: the continuous-measure Green error transfer and tested sparse
  two-atom application are implemented.
- Release readiness: not complete while Milestones D and E retain the statuses
  above.
