# A small upstream contribution before a new basis backend

## Proposed destination and change

Propose a short physical-accuracy and complex-matrix sampling recipe to the
[Sparse-IR ecosystem documentation](https://github.com/SpM-lab/sparse-ir-doc).
Its repository describes itself as the ecosystem theory/notation/library hub;
the current Python API is documented separately in the
[basis reference](https://sparse-ir.readthedocs.io/en/latest/basis.html).
The concrete executable example was [submitted publicly on 2026-10-25](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
after checking the current issue, notebook and contribution guidance. The exact
posted text was fetched and verified. Maintainer acceptance remains pending.

The independently usable candidate is
[`examples/sparse_ir_physical_recipe.py`](examples/sparse_ir_physical_recipe.py).
It requires only NumPy and sparse-ir, not kernelgecp or our benchmark harness.
It uses native `TauSampling.fit` and two real poles with positive Hermitian
complex matrix residues summing to identity. Strictly disjoint time grids and
both frequency signs are checked against analytic truth. Following upstream's
notation, ordinary `m` is converted to reduced `n=2*m+1`. Tests use the actual
installed provider, canonical normalization invariants and an independent
Gauss–Legendre Fourier check.

The proposed patch is one executable recipe and a small documentation note:

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

## Evidence and limits to send with the proposal

The [three-provider example](examples/dlr_interoperability.py) and immutable
cohorts use pydlr 1.0.1 and sparse-ir 2.1.6, with source/configuration hashes,
shared observations, preserved failures and independent 95-digit checking.
The benchmark uses the same cached SVD fitter for all providers: it is not a
comparison of their native optimized fitting algorithms. Its single-machine
timings do not establish a universal fastest method or rank minimality.

Do not initially ask upstream to absorb the full GECP package, replace native
SVE/DLR algorithms, add a Rust dependency, or endorse a numerical certificate.
The submitted comment offers this small self-contained recipe as a notebook
cell after the conditioning example, with a comparison-page link. Current
upstream commit `c94bf6ae3d3953ca9e2081a8be9998db348cdd24` already has released-IR
plots/tables; the recipe addresses complex matrix physical-accuracy checks
rather than duplicating them. Any upstream patch should follow the destination's
licensing/contribution process; this proposal does not change licenses here.

## Adoption gate

Actual outreach is recorded at
https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382 . The
posted code runs with Python 3.12.12, NumPy 2.5.4, sparse-ir 2.1.6 and
pylibsparseir 0.12.0: rank 14, maximum sampled time error `2.76169e-11` and
frequency error `4.94085e-11`, against `1e-8`. The updated standalone recipe
repeats those results and all nine actual-provider interoperability tests pass.
This closes the outreach requirement in `PRODUCTIZATION.md`. Acceptance and
external adoption remain pending unless observed; no upstream PR or merger is
claimed.
