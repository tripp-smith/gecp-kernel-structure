# Submitted Sparse-IR notebook contribution

## Destination and submitted change

The physical-accuracy and complex-matrix sampling recipe is submitted to the
[Sparse-IR ecosystem documentation](https://github.com/SpM-lab/sparse-ir-doc).
Its repository describes itself as the ecosystem theory/notation/library hub;
the current Python API is documented separately in the
[basis reference](https://sparse-ir.readthedocs.io/en/latest/basis.html).
The concrete executable example was [submitted publicly on 2026-10-25](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
after checking the current issue, notebook and contribution guidance. The exact
posted text was fetched and verified. A maintainer then
[invited an executable notebook PR](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6405100346).
[Upstream PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10) implements that
requested placement after the conditioning example. Its live page records
review, CI and merge status; submission and local checks do not establish
upstream acceptance.

The independently usable candidate is
[`examples/sparse_ir_physical_recipe.py`](examples/sparse_ir_physical_recipe.py).
It requires only NumPy and sparse-ir, not kernelgecp or our benchmark harness.
It uses native `TauSampling.fit` and two real poles with positive Hermitian
complex matrix residues summing to identity. Strictly disjoint time grids and
both frequency signs are checked against analytic truth. Following upstream's
notation, ordinary `m` is converted to reduced `n=2*m+1`. Tests use the actual
installed provider, canonical normalization invariants and an independent
Gauss–Legendre Fourier check.

The submitted patch adds three notebook cells:

- State beta, energy units, the physical Green-function minus sign, Fourier
  phase and reduced-frequency convention together.
- Preserve arbitrary trailing matrix axes and complex linearity.
- Check a requested physical absolute error in both domains, rather than
  treating a native truncation eps as that same external error contract.
- Distinguish independent physical truth from an iteration defect and from
  a rigorous continuous/rounding certificate.

## Why this is useful independently of GECP

Frozen [cohort A](experiments/data/dlr_applications_cohort_a.json.gz) retains
a calibration miss at cutoff 10; its final frequency error is 1.21e-8 against
target 1e-8. The independently checked follow-up
[cohort B](experiments/data/dlr_applications_cohort_b.json.gz) uses stronger
calibration, different native eps between methods and new validation points.
This is not a claim that Sparse-IR's documented native truncation semantics are
incorrect. It is a reason to teach the distinction between native parameters
and an application's physical error budget.

The separate solver pilot also illustrates observability: removing essential
low-frequency observations can leave a tiny sampled defect and poor unseen
predictions. The final workload retains those observations, verifies their
physical truth after the loop, and checks disjoint time/frequency points. These
lessons apply without installing GECP and without proposing a new IR method.

## Evidence and limits

The [three-provider example](examples/dlr_interoperability.py) and immutable
cohorts use pydlr 1.0.1 and sparse-ir 2.1.6, with source/configuration hashes,
shared observations, preserved failures and independent 95-digit checking.
The benchmark uses the same cached SVD fitter for all providers: it is not a
comparison of their native optimized fitting algorithms. Its single-machine
timings do not establish a universal fastest method or rank minimality.

The contribution uses only the documentation's existing NumPy and sparse-ir
dependencies. It does not ask upstream to absorb the GECP package, replace
native SVE/DLR algorithms, add a Rust dependency or endorse a numerical
certificate. The initial comment offered this self-contained recipe after the
conditioning example, with a comparison-page link. The notebook PR follows the
maintainer's requested placement. Its base commit
`c94bf6ae3d3953ca9e2081a8be9998db348cdd24` already has released-IR
plots/tables; the recipe addresses complex matrix physical-accuracy checks
rather than duplicating them. Any upstream patch should follow the destination's
licensing/contribution process; this proposal does not change licenses here.

## Recorded submission and local validation

Actual outreach is recorded at
https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382 . The
posted code runs with Python 3.12.12, NumPy 2.5.4, sparse-ir 2.1.6 and
pylibsparseir 0.12.0: rank 14, maximum sampled time error `2.76169e-11` and
frequency error `4.94085e-11`, against `1e-8`. The updated standalone recipe
repeats those results and all nine actual-provider interoperability tests pass.
This closed the original outreach requirement in `PRODUCTIZATION.md`.

PR #10's submitted commit is
`bc4659e14c1a5f319fd6dcc7fa42796c7b29e2bd`. Its only changed file,
`src/sparse_sampling_py.ipynb`, has SHA-256
`c93d98bcea0d0e1ef0b695931a78eda1001f6284d9986db20e1b2ecd94f0e16f`:
87 inserted lines, no deletions, and three new cells. The original eight cells
and notebook metadata remain unchanged; source outputs and execution counts
remain clear.

Using the upstream lock with Python 3.13.2, NumPy 2.3.5, sparse-ir 2.1.4 and
pylibsparseir 0.9.0, isolated new-cell execution and full-notebook execution
agree exactly. Both give rank 14, maximum entrywise absolute time error
`2.7616964271004463e-11` on 73 disjoint rows and frequency error
`4.94080361244419e-11` on 64 indices spanning both signs, against the unchanged
`1e-8` ceiling. Spectral positivity/normalization, coefficient/matrix shapes,
endpoint normalization, matrix adjoint symmetry and independent Fourier
quadrature pass; the quadrature discrepancy is below `4.78e-15`.
`uv run jupyter book build --all .` exits zero, with one warning that the
unchanged `README.md` is outside the toctree. Read-only independent review
checks the preserved notebook content, conventions and sampled claim scope.

The initial disposable validation harness failed because `nbformat.from_dict`
left list-valued JSON cell sources unnormalized for nbclient. Reading with
`nbformat.reads` corrected the harness; no notebook source fix was needed.
This failure is retained in the local contribution receipt rather than
represented as a failed numerical example.

At submission, [upstream CI 40248110146](https://github.com/SpM-lab/sparse-ir-doc/actions/runs/40248110146)
concluded `action_required` with no jobs: a permission gate, not passed upstream
execution. Consult [the live PR](https://github.com/SpM-lab/sparse-ir-doc/pull/10)
for later checks, review and merge decisions. Local contribution preparation
and validation are recorded in [PR #133](https://github.com/tripp-smith/gecp-kernel-structure/pull/133);
they do not imply maintainer acceptance, an upstream merge or a new scientific
result. The published kernelgecp 1.1.0 tag and artifacts are unchanged.
