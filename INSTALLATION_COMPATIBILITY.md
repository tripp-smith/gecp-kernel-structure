# Installation and tested compatibility

The public release is currently kernelgecp 1.1.0. The production integrations
in this checkout target 1.2.0 and remain subject to the
[acceptance ledger](PRODUCTIZATION.md#production-integrations-12).

```sh
python -m pip install kernelgecp==1.1.0
```

That installs the released base API. To work with the new integrations before
publication, use a verified source revision and its committed locks:

```sh
uv sync --all-extras
./scripts/verify.sh
```

Base dependencies are NumPy, SciPy and mpmath. Optional provider stacks remain
isolated. Installing development extras does not install TRIQS, native Green,
Quantum ESPRESSO/EPW or Torch. Missing providers fail their dedicated workflows;
they never silently select a different scientific method.

| Component | Tested runtime | Platform boundary | Setup and example |
| --- | --- | --- | --- |
| Base Python package | Python >=3.12; exact tested versions recorded by distribution CI | Portable Linux, macOS and Windows package lanes | [Quickstart](PACKAGE_README.md), `python -m examples.reusable_basis` |
| Optional Rust evaluator | Source-bound abi3 wheels; release CI records target/toolchain | Linux x86_64 and macOS arm64 wheels; source build documented separately | [Native build](native/SOURCE_BUILD.md) |
| TRIQS adapter / deterministic DMFT | TRIQS 4.0.2 | Linux x86_64 and macOS arm64 actual provider tests | [TRIQS setup](environments/TRIQS.md), `python -m examples.triqs_interoperability` |
| CT-HYB DMFT | CT-HYB 4.0.1 with TRIQS 4.0.2 | Canonical full scientific runs: Linux x86_64 | [CT-HYB guide](CTHYB_DMFT.md) |
| Molecular GF2 | Green 1.0.0; PySCF 2.14.0 references | Actual native solver: Linux x86_64 | [Molecular guide](MOLECULAR_GREEN.md) |
| Superconductivity interchange | qe-7.5 / EPW 6.0 | Actual native restart: Linux x86_64 | [EPW guide](EPW_INTEROPERABILITY.md) |
| Positive scalar/diagonal continuation | MaxEnt 4.0.1 | Linux x86_64 and macOS arm64 installed-provider cases; full EH acceptance pending | [Spectroscopy](SPECTROSCOPY.md) |
| CPU GP operator conversion | Torch 2.14.1, GPyTorch 1.15.2, linear_operator 0.6.1 | Linux x86_64 and macOS arm64 actual provider tests | [GP setup](environments/gp/README.md), [covariance guide](COVARIANCE.md) |

Historical TRIQS APIs, untested operating-system/provider combinations, GPU
training and differentiation through pivot selection are not promised. The
macOS GP environment's upstream wheel boundary is arm64/macOS >=14; it does
not redefine base-package support.

Use the platform-specific explicit conda locks for TRIQS, followed by the
hash-pinned MaxEnt and DLR/IR requirements. Use `uv sync --project
environments/gp --frozen` for CPU GP and `uv sync --project
environments/green-science --frozen` for molecular comparison. Native builders,
source revisions, MPI/HDF5 variants, package inventories and thread settings
are documented in [environments](environments/README.md). Lock resolution alone
is not execution evidence.

Canonical calculations use one MPI rank where applicable and one thread in
OMP, OpenBLAS, MKL, vecLib, NumExpr and BLIS. Hosted jobs have a 330-minute cap;
fixed-bath acquisition has a shorter 240-minute driver deadline and durable
segments. One memory-intensive local workload runs at a time. No GPU or paid
compute is required by the accepted program.

## Small application examples

In the corresponding provider environment, run:

```sh
python -m examples.interacting_dmft
python -m examples.triqs_interoperability
python -m examples.molecular_green input.h5 sim.h5 run-manifest.json --archive fit.npz
python -m examples.epw_restart_archive bundle.json restored.imag_aniso
python -m examples.spectral_continuation observations.npz --mu0 1e-8
python -m examples.covariance_operator
```

EPW consumes the archive's JSON manifest and its adjacent NPZ payload. Its
example writes native text plus a new hash-bound context sidecar and refuses
to overwrite existing files. The spectroscopy example consumes an existing
bound scalar/diagonal dataset and prints training diagnostics; it does not
replace held-out scientific acceptance. Add `--gpytorch` to the covariance
example only in the optional CPU GP runtime. These examples are small interface
checks, not the full scientific campaigns or new performance claims.
