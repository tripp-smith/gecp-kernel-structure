# Pinned Green provider runtime

This isolated Linux x86_64 image builds Green/WeakCoupling 1.0.0 and runs an
actual interacting GF2 calculation for H2. It is a provider-readiness gate for
the later molecular interoperability work. It does not install kernelgecp,
replace Green's internal representation, benchmark an acceleration, or certify
the physical GF2 approximation.

Build and run from the repository root, on a machine with a Docker daemon:

```sh
docker build -f environments/green/Dockerfile -t kernelgecp-green:eb .
mkdir -p artifacts
docker run --rm \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -v "$PWD/artifacts:/artifacts" kernelgecp-green:eb \
  /opt/venv/bin/python /workspace/scripts/providers/green_smoke.py \
  --output /artifacts/green.json --source-revision "$(git rev-parse HEAD)"
```

The GitHub provider workflow supplies the source revision and preserves the
artifacts. An unavailable local Docker daemon does not establish a provider
pass; native execution must succeed in that Linux workflow. Each invocation
requires a fresh artifact subdirectory, so an unsuccessful run's evidence is
not overwritten. A failed run removes any prior receipt at the requested output
path and does not create a new `status="passed"` receipt.

## What is pinned

- The Ubuntu 24.04 OCI base manifest is fixed by digest in the Dockerfile.
- `sources.json` fixes all eleven Green repositories, including transitive
  `green-opt`, to full Git revisions of the 1.0.0 release, and fixes Green's
  `magic_enum` dependency to its v0.9.3 commit. The builder checks
  each actual checkout and passes explicit FetchContent source-directory
  overrides with downloads/updates disabled during CMake configuration.
- `requirements.in` specifies the provider dependencies; `requirements.lock`
  pins the complete Python resolution and hashes. Installation requires binary
  distributions and verifies every installed distribution's archive hash.
  Green-mbtools is the released PyPI wheel, not a wheel rebuilt from the source
  checkout; the lock distinguishes that artifact from release source identity.
- The build records actual system-package versions, installed Python packages,
  native binary SHA256 and the source lock in `build-provenance.json`. Apt
  repositories remain external inputs; this is a source- and artifact-pinned
  build with recorded toolchain versions, not a bitwise-hermetic OS rebuild.

CMake must be at least 3.27 because the pinned green-grids/green-opt components
require it. The Ubuntu 24.04 environment supplies this. Compilation uses two
jobs; the molecular run uses one MPI process and one numerical-library thread.
No GPU or unrelated upstream large-system test suite is required by this gate.

The first native build at source `be38fce` failed during CMake generation because
`green-params` declares `magic_enum` through a nested FetchContent call. With
network fetching deliberately disabled, omitting that source left its target
undefined. The source lock now includes it and the builder supplies the matching
`FETCHCONTENT_SOURCE_DIR_MAGIC_ENUM` override. This repair keeps the disconnected
configuration policy; it is not evidence of a successful native GF2 run.

## Native calculation and acceptance

The smoke run generates new molecular input with the actual installed
Green-mbtools/PySCF provider and executes the actual built `mbpt.exe`:

- H2 bond length 0.74 Angstrom, STO-3G, restricted spin, symmetric Lowdin
  orthonormal storage, and a 512 MB integral-chunk budget.
- Physical inverse temperature 20 inverse Hartree, Green's shipped IR
  Lambda=10,000 grid, GF2, explicit 100-iteration budget, self-energy mixing
  weight 0.5, and energy-change stopping threshold 1e-10 Hartree.
- Actual native convergence marker, no maximum-iteration marker, complete
  final/current HDF5 iteration pair, finite 5D molecular arrays, final summed
  energy change below 1e-10, and maximum Green-function change below 1e-6.
- Nonzero dynamic self-energy and correlation energy, plus endpoint
  normalization and two-electron count within 1e-5. These are bounded sanity
  checks, not an exact physical solution or application-error certificate.

Do not confuse input generation's optional `--beta` Gaussian auxiliary-basis
parameter with the executable's physical `--BETA`; the former is not set here.
Green can exit zero and print `Completed.` after exhausting its iteration
budget. Neither that message nor the HDF5 `/iter` value is accepted as
convergence evidence.

The successful receipt is `artifacts/green.json`, with schema
`kernelgecp.provider-smoke.v1`, provider `green`, exact provider versions,
true checks, measurements and hashes of the retained evidence. Raw logs,
input/output HDF5, source/build provenance and a run manifest are preserved
under `artifacts/green/`. The run manifest records model, units, actual command
arguments and acceptance thresholds missing from Green's HDF5 output.

The later production adapter must independently evaluate Dyson consistency,
native-grid/cutoff refinement and representation errors. This availability gate
does not claim those later acceptance criteria.

Upstream references: [Green 1.0 source](https://github.com/Green-Phys/green-mbpt/tree/v1.0.0),
[molecular example](https://green-phys.org/docs/getting-started/examples/n2/),
[output contract](https://green-phys.org/docs/components/data-formats/output-files/),
and [convergence loop](https://github.com/Green-Phys/green-sc/blob/v1.0.0/src/green/sc/sc_loop.h).
