# Production provider environments

These are optional developer/CI runtimes for the EB–EJ production integration
contract. The base `pip install kernelgecp==1.1.0` dependencies and immutable
release artifacts are unchanged. Runtime availability is distinct from the
later domain's scientific accuracy/performance acceptance.

| Runtime | Pinned provider | Execution evidence |
| --- | --- | --- |
| [TRIQS](TRIQS.md) | TRIQS 4.0.2, CT-HYB 4.0.1, MaxEnt 4.0.1 | actual Linux x86_64/macOS arm64 numerical probes |
| [Green](green/README.md) | Green 1.0.0, PySCF 2.14.0 | native Linux interacting H2 GF2 convergence and stored-delta checks |
| [EPW](epw/README.md) | qe-7.5, EPW 6.0 | native Linux official coarse MgB2 SCF/DFPT/Eliashberg calculation |
| [GP](gp/README.md) | GPyTorch 1.15.2, CPU Torch 2.14.1 | actual Linux/macOS solves, determinant and Gaussian probability against SciPy |

All six platform/provider jobs pass at source `d5f58a1` in
[workflow 40403723259](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40403723259).
That is historical execution evidence; final frozen-source/root/CI and merged
main acceptance are recorded canonically in
[PR #134](https://github.com/tripp-smith/gecp-kernel-structure/pull/134).
Neither a successful build nor an import-only/skipped job closes a runtime gate.

## Exact conda replay

CI consumes platform-specific `@EXPLICIT` locks with artifact SHA256 hashes,
using micromamba 2.9.0-0. It then installs the separately hash-pinned pure
MaxEnt wheel with `--no-deps --require-hashes`. `triqs.yml` describes the initial
resolution; the explicit locks and companion inventories define replay.

```sh
micromamba create -n kernelgecp-triqs \
  -f environments/locks/triqs-linux-64-conda-explicit.txt
micromamba run -n kernelgecp-triqs python -m pip install \
  --no-deps --require-hashes -r environments/locks/triqs-maxent-pip.txt
micromamba run -n kernelgecp-triqs python scripts/providers/triqs_smoke.py \
  --source-revision "$(git rev-parse HEAD)" --output artifacts/triqs.json
```

Use `triqs-osx-arm64-conda-explicit.txt` on macOS arm64. Each companion JSON
records the actual package versions/builds, hashes, dependency constraints,
license metadata and platform. The exporter reads actual installed conda-meta;
it does not derive successful execution from package metadata.

The Mac lock was replayed into a fresh isolated conda environment and
re-exported byte-identically; a native resolvent calculation then passed.
The official micromamba binary also parsed all 131 SHA256-locked artifacts
successfully. Linux locks come from the actual successful Linux CI runtime.
Native builds record source/binary identities and actual toolchain packages;
their apt package manifests are observations, not a bitwise-hermetic OS claim.

## Scope, notices and resource limits

[provider-licenses.json](provider-licenses.json) records observed license/notice
files at exact source revisions. It explicitly records devXlib's missing
standalone license file and observed source-header notice. It does not infer
a repository-wide license from that header or replace the later release's
redistribution inventory. Provider runtimes and native applications are not
bundled into kernelgecp wheels by this phase.

Green and EPW builds use two jobs, and calculations one MPI rank and one
BLAS/OpenMP thread. Hosted jobs are capped at 330 minutes. Native Linux runs
are available in the dedicated workflow without requiring a local Docker
daemon. Local source-only checks never substitute for those executions.

Setup failures remain recorded: the actual TRIQS namespace is `triqs.mesh`,
Green's nested `magic_enum` needs an explicit disconnected source lock, and
molecular GF2 must use the generated integral directories explicitly. The
initial micromamba listing was not a replayable lock; the checked exporter
now writes canonical hash-locked explicit files. None was counted as success.
