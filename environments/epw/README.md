# Pinned QE 7.5 / EPW 6.0 native runtime

This Linux x86_64 provider image builds unmodified Quantum ESPRESSO `qe-7.5`
(commit `770a0b2d12928a67048e2f3da8d10d057e52179e`) and EPW 6.0. The native
EPW driver declares version 6.0; the upstream EPW README's 5.9 label is stale.
The image uses an immutable Ubuntu 24.04 image digest, two build jobs, system
BLAS/LAPACK/FFTW, one MPI rank per calculation, and one OpenMP/BLAS thread.
`/opt/venv/bin/python` is the image's isolated system-Python 3.12 interpreter.
There is no default dependency on QE, MPI or a Fortran compiler in kernelgecp.

`sources.lock.json` pins SHA256, byte size and immutable source commits for QE,
its Wannier90, MBD and deviceXlib submodules, and the two pseudopotentials.
DeviceXlib is required by this upstream release even for a CPU build. The
verified archives are extracted unchanged, and initialized submodule markers
prevent QE's makefiles from silently fetching different sources. Source files
and licenses remain available in `/opt/qe`; nothing is downloaded at smoke-run
time. OS build-package versions are recorded in `/opt/epw-system-packages.txt`;
the apt repository resolution is recorded build provenance, not an immutable
snapshot lock for every system package.

## Build and run

Run from the repository root on a Linux x86_64 Docker host:

```sh
docker build --platform linux/amd64 -f environments/epw/Dockerfile -t kernelgecp-epw:qe7.5 .
mkdir -p experiments/scratch/provider-smoke
docker run --rm --platform linux/amd64 \
  -v "$PWD:/workspace:ro" \
  -v "$PWD/experiments/scratch/provider-smoke:/artifacts" \
  kernelgecp-epw:qe7.5 \
  /opt/venv/bin/python /workspace/scripts/providers/epw_smoke.py \
  --output /artifacts/epw.json --source-revision "$(git rev-parse HEAD)"
```

Native builds and calculations are intentionally performed in the provider
container rather than the default `uv` environment. A cross-architecture Docker
emulation run on macOS arm64 is not the canonical performance environment.
The initial local host has no running Docker daemon, so local static checks
alone do not establish native acceptance. Exact-revision provider CI and the
emitted calculation receipt are the execution evidence.

## Actual calculation and acceptance

The smoke executes the official coarse MgB2 regression in
[`test-suite/epw_super`](https://github.com/QEF/q-e/tree/qe-7.5/test-suite/epw_super):

1. `pw.x -in scf.in`: the original 6x6x6 SCF calculation must converge.
2. `ph.x -in ph.in`: the original 3x3x3 DFPT phonon calculation must finish.
3. Collect its actual dynamical matrices and response potentials into `save/`.
4. `pw.x -in scf_epw.in`, then `pw.x -in nscf_epw.in`: the 3x3x3 electronic
   mesh and full 27-point NSCF calculation used by EPW.
5. `epw.x -in epw1.in`: the native Wannier/electron-phonon calculation and
   isotropic Eliashberg solve at 15 K, including the original convergence
   parameters and real-axis continuation settings.

No input physics parameters are modified. The work-tree layout preserves the
upstream `../../pseudo/` path. Source-manifest hashes bind the exact input
bytes and UPFs. The probe checks successful process exits and native completion,
SCF convergence, native version banners, all required PH output, native EPW
iteration convergence, a finite positive electron-phonon coupling, and a finite
positive-gap imaginary-frequency mesh. A help/version command, compiled binary,
existing fixture, skipped calculation or import by itself cannot pass.

The JSON receipt has schema `kernelgecp.provider-smoke.v1`, provider `epw`, and
`status="passed"` only after all numerical checks succeed. It records the
calling repository revision, source-lock/native source identity, executable,
input/log/output hashes, and native numerical results. Exceptions write a failed
receipt and exit nonzero. Reusing an existing native work directory fails;
choose a new artifact output directory for a fresh run.

Under the output directory, `epw-native/test-suite/epw_super/` retains the native
logs, `MgB2.imag_iso_015.00`, `.a2f`, and calculation state. The native source
manifest is copied to `epw-native/`. Image build/configuration logs are kept at
`/opt/epw-build.log` and `/opt/epw-configure.log`. CI should retain execution and
build logs on failure, and archive selected outputs rather than indiscriminately
uploading all wavefunctions. These files can seed later domain-specific fixture
selection, but do not close its numerical acceptance gates.

The smoke is a native runtime/installation regression. Its coarse meshes and
upstream loose material convergence settings do not establish a converged MgB2
prediction, a production adapter, a GECP speedup, or a superconductivity accuracy
certificate. EPW's GPL source and upstream pseudopotential notices remain in
the provider image; no provider source is copied into the kernelgecp wheel.
