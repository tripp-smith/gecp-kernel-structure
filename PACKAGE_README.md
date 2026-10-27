# kernelgecp

Reusable cross-approximation bases for fermionic Green functions, with explicit
precision and conditioning diagnostics, complex matrix batches, safe archives,
and an optional Rust backend for repeated reference time evaluation.

```sh
python -m pip install kernelgecp==1.1.0
```

Python 3.12 and 3.13 are tested. Native wheels cover Linux x86_64 with glibc
2.28 or newer and macOS arm64 11 or newer. Other platforms, including Windows,
use the universal Python wheel. Installing a wheel needs no Rust compiler.
NumPy, SciPy and mpmath are the runtime dependencies.

## Fit once, evaluate in both domains

```python
import numpy as np
from kernelgecp import FermionicBasis, FermionicFit

beta, energy = 10.0, 0.3
basis = FermionicBasis(beta=beta, wmax=1.0, tol=1e-9)
tau = basis.tau_nodes
values = -np.exp(-tau * energy) / (1 + np.exp(-beta * energy))
fit = basis.fit_tau(values, tau)

print(fit.evaluate_tau(np.array([0.2, 1.0, 8.0])))
print(fit.evaluate_matsubara(np.arange(-2, 3)))
print(basis.diagnostics)

fit.save("green.npz")
restored = FermionicFit.load("green.npz")
```

`tau` uses physical units in `[0,beta]`, and `wmax` bounds physical spectral
energies. The convention is `G(tau)=-integral K(tau/beta,beta*energy) dρ`.
Matsubara inputs are **ordinary integer indices** `n`, giving frequencies
`(2*n+1)*pi/beta`. The one-pole frequency answer above is
`1/(1j*(2*n+1)*pi/beta-energy)`. Complex values and trailing matrix axes are
supported; `axis=` selects the observation axis. The basis caches factorizations
and evaluation designs for reuse.

The development checkout additionally supports Matsubara observation fitting
and optional TRIQS 4.0.2 interoperability, targeting 1.2.0. These are not part of
the immutable 1.1.0 package installed above. See the
[development API and archive migration notes](https://github.com/tripp-smith/gecp-kernel-structure/blob/main/APPLICATION.md#frequency-fitting-and-triqs-interoperability-unreleased-ec-source).

## Optional reference acceleration

```python
reference_basis = FermionicBasis(
    beta=20.0,
    wmax=50.0,
    tol=1e-10,
    precision="reference",
    precision_bits=128,
    reference_backend="rust",
)
```

Python remains the default. Rust accelerates only 128-bit reference time-kernel
evaluation and cached triangular substitution. Selection and frequency
evaluation remain mpmath; fits and final products remain NumPy/SciPy. Diagnostics
report each stage. A missing native module or unsupported precision raises a
clear error. Loading an archive restores Python execution by default; native
restoration requires an explicit `reference_backend="rust"` request.

The frozen rank-55 benchmark, including construction and fitting followed by
80 changing 64-row time meshes, measured a 50.15% median application-time
reduction on one macOS arm64 machine. The eight-query prefix improved 8.92%,
and warm cache hits showed little benefit. This is not a universal speed claim
or a claim to outperform specialized DLR/IR libraries. See the
[methods and complete experiment records](https://github.com/tripp-smith/gecp-kernel-structure/blob/main/METHODS.md).

## Accuracy and certification

`tol` controls the numerical basis construction; it is not a promised physical
error for arbitrary data. Check time and frequency predictions against the
accuracy required by your application. Diagnostics expose effective rank,
conditioning, stopping and arithmetic choices. Floating-point results are
numerical evidence.

`certify_fit_precision` can enclose actual coefficient/solve/output discrepancies
against an explicitly defined exact-data least-squares reference at supplied
queries. It can return inconclusive or large-error results. It does not bound
spectral/model errors or establish a continuous inexact-pivot convergence rate.
The Lean exact-run results retain their separate mathematical assumptions.

The Python/source project is Apache-2.0. Native wheels carry the composite
license expression for their bundled dependencies, full notices, a sanitized
dependency inventory and the matching complete source-asset identity. Source
and rebuild instructions accompany the
[versioned release](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0).

Read the [API and tutorials](https://github.com/tripp-smith/gecp-kernel-structure/blob/main/APPLICATION.md),
[runnable examples](https://github.com/tripp-smith/gecp-kernel-structure/tree/main/examples),
[formal scope audit](https://github.com/tripp-smith/gecp-kernel-structure/blob/main/FINAL_AUDIT.md),
and [citation metadata](https://github.com/tripp-smith/gecp-kernel-structure/blob/main/CITATION.cff).
