# Molecular Green GF2 interoperability

This integration targets kernelgecp 1.2.0. Public 1.1.0 artifacts are unchanged.
EF is delivered in [PR #137](https://github.com/tripp-smith/gecp-kernel-structure/pull/137#issuecomment-6481863146);
the complete 1.2 release still requires EH and EJ acceptance.

`kernelgecp.integrations.green.read_green_tau(input_path, results_path,
run_manifest, iteration=None)` reads actual restricted molecular Green 1.0.0
GF2 output. Install the isolated [scientific environment](environments/green-science/README.md)
for PySCF references; the importer alone needs the optional `green` extra
(h5py). Neither native Green nor PySCF is a default package dependency.

## Data contract

Supply a `kernelgecp.green-molecular-run.v1` JSON manifest bound to the actual
input, results, native convergence log and provider receipt. The importer checks
file hashes, complete iterations, units, restricted GF2 identity and declared
symmetric Lowdin storage. Native output alone does not establish provenance.
Use `scripts/providers/bind_green_manifest.py` with the native provider receipt;
do not manufacture a converged manifest from an existing checkpoint.

Returned arrays are immutable owned copies. The axis order is
`(tau, spin, k, orbital_row, orbital_column)`; the molecular k axis is retained.
The record preserves overlap, one-body matrices, chemical potential, static
self-energy, selected/source iterations, convergence evidence and energy
components. `iteration=None` selects the complete converged final checkpoint.
An explicitly selected complete earlier checkpoint is identified separately.
The distinct time endpoints are retained. No spin averaging, orbital rotation,
conjugation, static-self-energy subtraction or tail repair is implicit.

From a checkout containing the integration:

```sh
uv sync --all-extras
uv run python -m examples.molecular_green input.h5 sim.h5 run-manifest.json \
  --archive h2-fit.npz --wmax 4
```

The command fits the full time-dependent Green matrix and prints its observed
training residual. Supply a physical cutoff validated for your calculation;
the example's 4-Hartree default does not certify a new molecule. Keep the native
record and run manifest alongside the fit archive: fit coefficients alone do
not carry every molecular observable or its provenance.

## Scientific evidence and reproduction

The accepted cases are restricted H2/STO-3G at 0.74 and 1.5 Angstrom, beta=20
inverse Hartree, with symmetric Lowdin orthonormalization. The protocol in
`experiments/molecular_green.toml` declares the native grid, physical cutoffs,
separate dynamic-self-energy cutoff, calibration and untouched final queries.
Independent finite-temperature Lehmann references use PySCF integrals and
sector energies checked against FCI.

The comparison separates same-GF2 compression (limit 1e-6), native-grid and
physical-cutoff refinement, independently evaluated Dyson defect, and GF2 versus
FCI model discrepancy. A good compression fit does not make GF2 exact. Read the
canonical phase receipt for the accepted source, actual native runs and measured
values. Failed first-campaign cutoffs remain in
`experiments/data/molecular-green-6050ddd.tar.gz`; the replacement evidence is
`experiments/data/molecular-green-56c3211.tar.gz`. Their original request paths
are provenance, not portable inputs.

The [molecular native workflow](.github/workflows/molecular-green.yml) rebuilds
and runs the unmodified pinned solver on Linux without a local Docker daemon.
The [environment guide](environments/green-science/README.md) supplies optional
local container commands. Release validation must rerun the actual converged
solver, not infer success from exit zero, imports or frozen PR fixtures.
The importer supports compression, evaluation and archival; it does not replace
Green's internal transforms or accelerate its GF2 solver.
